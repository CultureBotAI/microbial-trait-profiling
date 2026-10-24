"""Bind downstream analyses to verified core products, never filenames alone (#629)."""

from __future__ import annotations

import argparse
import json
import uuid
from collections.abc import Sequence
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Literal

import yaml

from mtp.corpus_identity import CorpusIdentity, KgInput, measure_matrix
from mtp.count_input import ProfileCoverage
from mtp.datamodel.core_products_pydantic import SourceIdentity
from mtp.datamodel.workflow_pydantic import (
    ConsumerRecord,
    FileInventoryEntry,
    LegacyImportRecord,
    WorkflowRun,
)
from mtp.human_host_corpus import CORPUS_SEQUENCING_TYPES
from mtp.paths import file_sha256
from mtp.workflow.locations import LocationResolver, locate
from mtp.workflow.runtime import artifact_path
from mtp.workflow.verification import verify_core_run

ProfileKind = Literal["count_profiles", "ora_profiles"]
PROFILE_UNITS = {
    "count_profiles": "rolled_up_source_taxon_count",
    "ora_profiles": "distinct_matched_reported_taxon_count",
}


@dataclass(frozen=True)
class CoreBundle:
    """A successfully verified execution and its declared scientific input products."""

    manifest: Path
    run: WorkflowRun
    profile_product: ProfileKind
    cohort: str
    resolver: LocationResolver | None = None

    def path(self, product_id: str) -> Path:
        """Resolve only an artifact declared by this execution."""
        return artifact_path(self.run, product_id, self.resolver)


def open_core_bundle(
    manifest: Path,
    profile_product: ProfileKind = "count_profiles",
    *,
    expected_cohort: str | None = None,
) -> CoreBundle:
    """Verify bytes, producer completion, units, coverage and catalog before consumption."""
    run = verify_core_run(manifest)
    products = {product.id: product for product in run.definition.products}
    artifacts = {artifact.product_id: artifact for artifact in run.artifacts or []}
    profile = products[profile_product]
    artifact = artifacts[profile_product]
    if profile.kind != profile_product or artifact.unit != PROFILE_UNITS[profile_product]:
        raise ValueError(f"{profile_product}: incompatible profile kind or counting unit")
    if artifact.feature_policy != "include_mediadive":
        raise ValueError("Core catalog must retain MediaDive; filtering belongs to the consumer")
    cohort = artifact.cohort or "unknown"
    expected = expected_cohort or ("all" if profile_product == "count_profiles" else "Homo sapiens")
    if cohort != expected:
        raise ValueError(f"{profile_product}: cohort {cohort!r} does not match {expected!r}")
    coverage_id = "count_coverage" if profile_product == "count_profiles" else "ora_coverage"
    for product_id in (coverage_id, "feature_catalog", "experiments", "taxa_matrix"):
        artifact_path(run, product_id)
    if artifacts[coverage_id].cohort != cohort:
        raise ValueError("Coverage and profile cohorts disagree")
    return CoreBundle(manifest.resolve(), run, profile_product, cohort)


def bind_count_arguments(args: argparse.Namespace) -> CoreBundle | None:
    """Bind count-analysis CLI paths, rejecting explicit conflicting input files."""
    if args.core_run is None:
        return None
    bundle = open_core_bundle(args.core_run)
    bindings = {
        "profile_input": "count_profiles",
        "experiments_input": "experiments",
        "full_dump_input": "raw_dump",
    }
    if hasattr(args, "column_labels"):
        bindings["column_labels"] = "feature_metadata"
    for argument, product_id in bindings.items():
        supplied = getattr(args, argument)
        resolved = bundle.path(product_id)
        if supplied is not None and supplied.resolve() != resolved.resolve():
            raise ValueError(f"--core-run conflicts with explicit --{argument.replace('_', '-')}")
        setattr(args, argument, resolved)
    return bundle


def core_corpus_identity(bundle: CoreBundle, *, bugsigdb_vintage: str | None) -> CorpusIdentity:
    """Use only the bound run's KG bytes; a verified BugSigDB pin is not a KG pin."""
    path = bundle.path("source_identity")
    source = SourceIdentity.model_validate_json(path.read_text(encoding="utf-8"))
    return CorpusIdentity(
        matrix=measure_matrix(bundle.path(bundle.profile_product)),
        bugsigdb_vintage=bugsigdb_vintage,
        kg_release_tag=None,
        kg_provenance_path=str(path),
        kg_provenance_source="verified core execution input bytes; KG release identity unverified",
        kg_inputs=[
            KgInput(
                name=Path(record.path).name,
                declared_md5=None,
                declared_size=None,
                fetched_at=None,
                drive_file_id=None,
                present=True,
                path=record.path,
                sha256=record.sha256,
                measured_size=record.size_bytes,
            )
            for name in ("kg_edges", "kg_nodes")
            for record in [source.consumed_files[name]]
        ],
    )


def preflight_consumer_outputs(
    bundle: CoreBundle,
    record_path: Path,
    outputs: Sequence[Path] = (),
    *,
    output_directory: Path | None = None,
) -> None:
    """Reserve sidecars and protect the core before analysis writes any output (#654).

    Directory-based consumers use a dedicated destination disjoint from their
    input files and the immutable run bundle. Existing aliases in that directory
    are also checked, so a symlink or hard link cannot redirect a result write.
    This preflight does not claim protection against concurrent filesystem edits.
    """
    paths = [*outputs, record_path, record_path.with_suffix(".coverage.tsv")]
    seen_paths: set[Path] = set()
    seen_files: set[tuple[int, int]] = set()
    for path in paths:
        resolved = path.resolve()
        identity = _existing_file_identity(path)
        if resolved in seen_paths or (identity is not None and identity in seen_files):
            raise ValueError(f"Consumer output paths collide, including reserved sidecars: {path}")
        seen_paths.add(resolved)
        if identity is not None:
            seen_files.add(identity)
        if path.exists() and not path.is_file():
            raise ValueError(f"Consumer output must be a regular file: {path}")

    protected = {bundle.manifest.parent.resolve()}
    protected_files = [
        bundle.manifest,
        *(
            bundle.manifest.parent / "steps" / step.step_id / "_step.json"
            for step in bundle.run.steps or []
        ),
    ]
    for artifact in bundle.run.artifacts or []:
        path = bundle.path(artifact.product_id)
        protected.add(path.resolve())
        if path.is_dir():
            protected_files.extend(path / member.path for member in artifact.inventory)
        else:
            protected_files.append(path)
    protected_identities = {
        identity
        for path in protected_files
        if (identity := _existing_file_identity(path)) is not None
    }

    candidates = list(paths)
    if output_directory is not None:
        directory = output_directory.resolve()
        if any(_paths_overlap(directory, path) for path in protected):
            raise ValueError(
                f"Consumer output directory overlaps verified core inputs or run bundle: "
                f"{output_directory}; choose a separate analysis directory"
            )
        if any(sidecar.is_symlink() for sidecar in paths[-2:]):
            raise ValueError("Consumer output sidecars cannot redirect through symlinks")
        declared_paths = {path.absolute() for path in paths}
        for child in output_directory.rglob("*"):
            if child.is_symlink() and not child.exists():
                raise ValueError(f"Consumer output directory contains a dangling symlink: {child}")
            if child.is_file() and child.absolute() not in declared_paths:
                identity = _existing_file_identity(child)
                if child.resolve() in seen_paths or identity in seen_files:
                    raise ValueError(f"Consumer output paths collide through an alias: {child}")
                seen_paths.add(child.resolve())
                if identity is not None:
                    seen_files.add(identity)
            candidates.append(child)
    for candidate in candidates:
        resolved = candidate.resolve()
        if any(_paths_overlap(resolved, path) for path in protected) or (
            _existing_file_identity(candidate) in protected_identities
        ):
            raise ValueError(f"Consumer output overlaps a verified core input or run: {candidate}")


def _existing_file_identity(path: Path) -> tuple[int, int] | None:
    if not path.is_file():
        return None
    status = path.stat()
    return status.st_dev, status.st_ino


def _paths_overlap(left: Path, right: Path) -> bool:
    return left.is_relative_to(right) or right.is_relative_to(left)


def write_consumer_record(
    bundle: CoreBundle,
    path: Path,
    *,
    consumer: str,
    outputs: Sequence[Path],
    feature_policy: str = "exclude_mediadive",
    additional_products: Sequence[str] = (),
) -> Path:
    """Bind outputs and optional extra scientific inputs to one checked core execution."""
    if not outputs:
        raise ValueError("Consumer outputs must be a nonempty set of distinct files")
    coverage_id = "count_coverage" if bundle.profile_product == "count_profiles" else "ora_coverage"
    product_ids = [
        bundle.profile_product,
        coverage_id,
        "feature_catalog",
        "feature_metadata",
        "experiments",
        "taxa_matrix",
        "raw_dump",
        "studies",
        "source_identity",
    ]
    artifacts = {artifact.product_id: artifact for artifact in bundle.run.artifacts or []}
    products = {product.id: product for product in bundle.run.definition.products}
    for product_id in additional_products:
        if product_id in product_ids:
            raise ValueError(f"Additional consumed product is duplicated: {product_id}")
        if product_id not in products or product_id not in artifacts:
            raise ValueError(f"Additional consumed product is unresolved: {product_id}")
        if products[product_id].kind in PROFILE_UNITS:
            raise ValueError("Additional profile products require a separate consumer record")
        bundle.path(product_id)
        product_ids.append(product_id)
    preflight_consumer_outputs(bundle, path, outputs)
    current = verify_core_run(bundle.manifest)
    if (
        current.run_id != bundle.run.run_id
        or current.execution_fingerprint != bundle.run.execution_fingerprint
        or current.artifacts != bundle.run.artifacts
    ):
        raise ValueError("Core inputs or execution identity changed during downstream analysis")
    path.parent.mkdir(parents=True, exist_ok=True)
    coverage_path = build_consumer_coverage(
        bundle, path.with_suffix(".coverage.tsv"), feature_policy
    )
    coverage_facts = validate_consumer_coverage(bundle, coverage_path, feature_policy)
    record = ConsumerRecord(
        id=f"consumer:{uuid.uuid4()}",
        run_id=bundle.run.run_id,
        run_manifest_path=str(bundle.manifest),
        run_manifest_sha256=file_sha256(bundle.manifest),
        consumer=consumer,
        product_ids=product_ids,
        artifact_ids=[artifacts[product].id for product in product_ids],
        feature_policy=feature_policy,
        cohort=bundle.cohort,
        created_at=datetime.now(UTC),
        coverage_path=str(coverage_path.resolve()),
        feature_filter_id=str(coverage_facts["feature_filter_id"]),
        coverage_scope=str(coverage_facts["coverage_scope"]),
        outputs=[
            FileInventoryEntry(
                path=str(output.resolve()),
                sha256=file_sha256(output),
                size_bytes=output.stat().st_size,
            )
            for output in [*outputs, coverage_path]
        ],
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(record.model_dump_json(indent=2) + "\n", encoding="utf-8")
    return path


def verify_consumer_record(
    path: Path, *, resolver: LocationResolver | None = None
) -> ConsumerRecord:
    """Check a deck/figure input's entire recorded chain before accepting its outputs."""
    record = ConsumerRecord.model_validate_json(locate(path, resolver).read_text(encoding="utf-8"))
    manifest = locate(Path(record.run_manifest_path), resolver)
    if file_sha256(manifest) != record.run_manifest_sha256:
        raise ValueError("Consumer's core manifest changed")
    run = verify_core_run(manifest, resolver=resolver)
    if run.run_id != record.run_id:
        raise ValueError("Consumer references a different execution")
    artifacts = {artifact.id: artifact for artifact in run.artifacts or []}
    if len(record.product_ids) != len(record.artifact_ids):
        raise ValueError("Consumer's product and artifact references differ in length")
    if len(set(record.product_ids)) != len(record.product_ids):
        raise ValueError("Consumer repeats a product reference")
    for product, artifact_id in zip(record.product_ids, record.artifact_ids, strict=True):
        if artifact_id not in artifacts or artifacts[artifact_id].product_id != product:
            raise ValueError("Consumer has an unresolved product/artifact reference")
    profile_artifacts = [
        artifacts[item]
        for item in record.artifact_ids
        if artifacts[item].product_id in PROFILE_UNITS
    ]
    if len(profile_artifacts) != 1 or profile_artifacts[0].cohort != record.cohort:
        raise ValueError("Consumer cohort does not match its profile artifact")
    profile = profile_artifacts[0]
    definitions = {product.id: product for product in run.definition.products}
    if (
        profile.unit != PROFILE_UNITS[profile.product_id]
        or definitions[profile.product_id].kind != profile.product_id
    ):
        raise ValueError("Consumer profile has incompatible kind or counting unit")
    coverage = "count_coverage" if profile.product_id == "count_profiles" else "ora_coverage"
    if not {
        coverage,
        "feature_catalog",
        "feature_metadata",
        "experiments",
        "taxa_matrix",
        "raw_dump",
        "studies",
        "source_identity",
    }.issubset(record.product_ids):
        raise ValueError("Consumer is missing required core references")
    if record.feature_policy not in {"include_mediadive", "exclude_mediadive"}:
        raise ValueError("Unknown consumer feature policy")
    if not record.outputs or len({item.path for item in record.outputs}) != len(record.outputs):
        raise ValueError("Consumer outputs must be a nonempty set of distinct files")
    for output in record.outputs:
        actual = locate(Path(output.path), resolver)
        if actual.stat().st_size != output.size_bytes or file_sha256(actual) != output.sha256:
            raise ValueError(f"Consumer output changed: {actual}")
    if record.coverage_path not in {output.path for output in record.outputs}:
        raise ValueError("Consumer coverage is not bound to a hashed output")
    profile_kind: ProfileKind = (
        "count_profiles" if profile.product_id == "count_profiles" else "ora_profiles"
    )
    facts = validate_consumer_coverage(
        CoreBundle(manifest, run, profile_kind, record.cohort, resolver),
        locate(Path(record.coverage_path), resolver),
        record.feature_policy,
    )
    if (
        facts["feature_filter_id"] != record.feature_filter_id
        or facts["coverage_scope"] != record.coverage_scope
    ):
        raise ValueError("Consumer coverage filter or scope does not match the recorded policy")
    return record


def import_legacy_provenance(source: Path, output: Path) -> Path:
    """Preserve historical claims verbatim as unverified, even if they pass an old schema."""
    claims = yaml.safe_load(source.read_text(encoding="utf-8"))
    record = LegacyImportRecord(
        id=f"legacy:{uuid.uuid4()}",
        source_path=str(source.resolve()),
        source_sha256=file_sha256(source),
        imported_at=datetime.now(UTC),
        verification_status="imported_unverified",
        claims_json=json.dumps(claims, sort_keys=True, default=str),
    )
    if output.exists():
        raise FileExistsError(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(record.model_dump_json(indent=2) + "\n", encoding="utf-8")
    return output


def _consumer_coverage(
    bundle: CoreBundle, feature_policy: str
) -> tuple[dict[tuple[str, str, str], ProfileCoverage], tuple[str, ...]]:
    """Calculate annotation coverage for the input cohort and retained vocabulary."""
    from mtp.study_trait_profiles import collect_profile_taxa, read_taxa_trait_matrix
    from mtp.workflow.coverage import build_profile_coverage, feature_class

    if feature_policy not in {"include_mediadive", "exclude_mediadive"}:
        raise ValueError(f"Unknown consumer feature policy: {feature_policy!r}")
    matrix = read_taxa_trait_matrix(bundle.path("taxa_matrix"))
    columns = tuple(
        column
        for column in matrix.trait_columns
        if feature_policy == "include_mediadive" or feature_class(column) != "culture_medium"
    )
    cohort = bundle.cohort.strip()
    if not cohort or cohort == "unknown":
        raise ValueError("Consumer annotation coverage requires an explicit input cohort")
    profiles, _ = collect_profile_taxa(
        bundle.path("studies"),
        host_species=None if cohort == "all" else cohort,
        sequencing_types=CORPUS_SEQUENCING_TYPES,
    )
    return build_profile_coverage(profiles, matrix, feature_columns=columns), columns


def build_consumer_coverage(
    bundle: CoreBundle, path: Path, feature_policy: str = "exclude_mediadive"
) -> Path:
    """Write retained-feature annotation coverage before analysis-specific exclusions.

    The rows describe the core input cohort. Pairing, human-only selection by a
    downstream analysis, body-site selection, and zero-total exclusions require
    their own analysis records and are not inferred from this annotation table.
    A MediaDive-only taxon becomes matched-without-retained-annotations when
    medium features are excluded; it never becomes lookup-unmatched.
    """
    from mtp.workflow.coverage import write_profile_coverage

    coverage, _ = _consumer_coverage(bundle, feature_policy)
    path.parent.mkdir(parents=True, exist_ok=True)
    write_profile_coverage(path, coverage)
    return path


def validate_consumer_coverage(
    bundle: CoreBundle, path: Path, feature_policy: str = "exclude_mediadive"
) -> dict[str, object]:
    """Recompute retained-feature coverage and reject stale or all-feature substitutes."""
    from mtp.count_input import read_profile_coverage
    from mtp.workflow.coverage import feature_filter_id

    expected, columns = _consumer_coverage(bundle, feature_policy)
    observed = read_profile_coverage(path)
    if observed != expected:
        raise ValueError("Consumer coverage differs from the input cohort and retained features")
    return {
        "feature_filter_id": feature_filter_id(columns),
        "features": len(columns),
        "profiles": len(observed),
        "cohort": bundle.cohort,
        "coverage_scope": "input_cohort_before_analysis_specific_exclusions",
        "feature_policy": feature_policy,
    }
