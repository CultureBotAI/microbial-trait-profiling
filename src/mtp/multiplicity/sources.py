"""Verified execution and independently checked analysis bases for sidecars (#853)."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from mtp.count_input import profile_count_reader
from mtp.datamodel.multiple_testing_pydantic import SourceArtifact
from mtp.datamodel.workflow_pydantic import AnalysisBasis, ArmRun, ComparisonRun, WorkflowRun
from mtp.paths import file_sha256
from mtp.workflow.arms.runtime import verify_arm
from mtp.workflow.comparison import verify_comparison_run
from mtp.workflow.consumers import CoreBundle
from mtp.workflow.coverage import feature_filter_id
from mtp.workflow.runtime import artifact_path
from mtp.workflow.verification import verify_core_run


def source(path: Path) -> SourceArtifact:
    """Use content-addressed product identity; machine-local paths stay outside families."""
    digest = file_sha256(path)
    return SourceArtifact(identity=f"{path.name}@sha256:{digest}", sha256=digest)


def recorded_arm_output(run: ArmRun, step_id: str, filename: str) -> Path:
    """Resolve a unique file from a recorded stage, including directory inventories."""
    steps = [step for step in run.steps or [] if step.id == step_id]
    if len(steps) != 1:
        raise ValueError(f"Expected one recorded arm stage: {step_id}")
    paths: set[Path] = set()
    for recorded in steps[0].outputs or []:
        root = Path(recorded.path)
        if recorded.storage_kind == "file":
            if root.name == filename:
                paths.add(root.resolve())
            continue
        for entry in recorded.inventory:
            member = root / entry.path
            if not member.resolve().is_relative_to(root.resolve()):
                raise ValueError("Recorded arm member escapes its directory")
            if member.name == filename:
                paths.add(member.resolve())
    if len(paths) != 1:
        raise ValueError(f"Expected one recorded {step_id}/{filename}; found {len(paths)}")
    return paths.pop()


def recorded_comparison_output(run: ComparisonRun, filename: str) -> Path:
    """Require the native source to be inside the verified comparison inventory."""
    root = Path(run.manifest_path).parent.resolve()
    paths = {
        (root / entry.path).resolve() for entry in run.outputs if Path(entry.path).name == filename
    }
    if len(paths) != 1 or not all(path.is_relative_to(root) for path in paths):
        raise ValueError(f"Expected one comparison-owned {filename}; external source unsupported")
    return paths.pop()


def retained_features(profile: Path, catalog: Path) -> tuple[str, ...]:
    """Check catalog membership as a set; retain the actual profile column order."""
    columns = tuple(
        column
        for column in profile_count_reader(profile).count_columns
        if column != "taxa_count" and "mediadive.medium" not in column.lower()
    )
    with catalog.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    catalog_ids = [row["column_id"] for row in rows]
    traits = [row["column_id"] for row in rows if row["feature_class"] == "trait_annotation"]
    if len(set(catalog_ids)) != len(catalog_ids) or len(set(columns)) != len(columns):
        raise ValueError("Duplicate feature inventory")
    if set(columns) != set(traits):
        raise ValueError("Profile retained features do not match the core catalog")
    return columns


def validate_basis(
    basis: AnalysisBasis, *, product: str, features: tuple[str, ...], catalog: Path, scope: str
) -> None:
    """Reject producer metadata that cannot describe the independently read core."""
    if basis.profile_unit != product:
        raise ValueError("AnalysisBasis profile unit mismatch")
    if basis.feature_catalog_sha256 != file_sha256(catalog):
        raise ValueError("AnalysisBasis catalog identity mismatch")
    if basis.feature_filter_id != feature_filter_id(features):
        raise ValueError("AnalysisBasis ordered retained feature identity mismatch (#860)")
    if basis.retained_feature_count != len(features):
        raise ValueError("AnalysisBasis retained feature count mismatch")
    if scope not in {"off", "phenotype"}:
        raise ValueError("Unsupported strain rollup scope")
    if basis.rollup_policy != f"ncbi_taxonomy_closure_rollup_{scope}":
        raise ValueError("AnalysisBasis strain rollup scope mismatch")


def bind_analysis_basis(
    bundle: CoreBundle, basis: AnalysisBasis
) -> tuple[AnalysisBasis, tuple[str, ...]]:
    """Validate an upstream basis against a verified bundle; return ordered features.

    The caller supplies the actual upstream AnalysisBasis, never synthesizes a
    replacement from current files. The CoreBundle must come from open_core_bundle.
    Historical runs without a strain_rollup parameter explicitly used off.
    """
    parameters = {item.name: item.value for item in bundle.run.parameters or []}
    scope = parameters.get("strain_rollup", "off")
    catalog = bundle.path("feature_catalog")
    features = retained_features(bundle.path(bundle.profile_product), catalog)
    validate_basis(
        basis, product=bundle.profile_product, features=features, catalog=catalog, scope=scope
    )
    return basis, features


@dataclass(frozen=True)
class VerifiedSources:
    """Source records after runtime checks and independent feature-basis checks."""

    comparison: ComparisonRun
    core: WorkflowRun
    arms: dict[str, ArmRun]
    bases: dict[str, AnalysisBasis]
    features: dict[str, tuple[str, ...]]
    scope: str

    def product(self, name: str) -> Path:
        return artifact_path(self.core, name)


def open_verified_sources(path: Path) -> VerifiedSources:
    """Verify one completed comparison and each core/arm source before adaptation."""
    comparison = verify_comparison_run(path)
    core = verify_core_run(Path(comparison.core_manifest_path))
    arms = {
        name: verify_arm(Path(getattr(comparison, f"{name}_manifest_path")))
        for name in ("noel", "jed", "ilyas")
    }
    if (
        comparison.contrast_population != core.contrast_population
        or comparison.direction_basis != core.direction_basis
    ):
        raise ValueError("Comparison population/direction differs from the verified core")
    for name, arm in arms.items():
        if arm.arm != name or arm.core_manifest_sha256 != comparison.core_manifest_sha256:
            raise ValueError("Comparison arm identity or core binding mismatch")
    parameters = {item.name: item.value for item in core.parameters or []}
    scope = parameters.get("strain_rollup", "off")
    catalog = artifact_path(core, "feature_catalog")
    bases: dict[str, AnalysisBasis] = {}
    features: dict[str, tuple[str, ...]] = {}
    for product in ("count_profiles", "ora_profiles"):
        matches = [basis for basis in comparison.analysis_bases if basis.profile_unit == product]
        if len(matches) != 1:
            raise ValueError(f"Expected exactly one AnalysisBasis for {product}")
        features[product] = retained_features(artifact_path(core, product), catalog)
        validate_basis(
            matches[0], product=product, features=features[product], catalog=catalog, scope=scope
        )
        bases[product] = matches[0]
    return VerifiedSources(comparison, core, arms, bases, features, scope)
