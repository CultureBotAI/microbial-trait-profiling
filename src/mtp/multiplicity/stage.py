"""Additive post-arm correction, migration reports and strict source replay (#853)."""

from __future__ import annotations

import csv
import json
import sqlite3
import tempfile
from contextlib import closing
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import scipy
from git import Repo

from mtp.datamodel.multiple_testing_pydantic import (
    ExternalFamilyCorrection,
    PostArmMultiplicityRun,
    SourceArtifact,
)
from mtp.multiple_testing import canonical_sha256, load_policy
from mtp.multiplicity.families import FamilyInput, jed_family, noel_family
from mtp.multiplicity.inventory import correct_inventory, json_text, read_corrected_inventory
from mtp.multiplicity.native import native_family
from mtp.multiplicity.ora import ora_family
from mtp.multiplicity.sources import bind_recorded_hashes, open_verified_sources, source
from mtp.paths import file_sha256


def code_sha256() -> str:
    """Bind tracked canonical modules; ignore untracked scratch and editor lock files.

    Entry/final disk identity does not attest bytecode imported before run entry.
    """
    package = Path(__file__).resolve().parents[1]
    with Repo(package, search_parent_directories=True) as repository:
        if repository.working_tree_dir is None:
            raise ValueError("Multiplicity execution requires a Git working tree")
        root = Path(repository.working_tree_dir).resolve()
        relative_package = package.relative_to(root).as_posix()
        prefix = "" if relative_package == "." else relative_package + "/"
        paths = []
        for name, stage in repository.index.entries:
            if not (str(name).startswith(prefix) and str(name).endswith(".py")):
                continue
            if stage != 0:
                raise ValueError(
                    f"Resolve the tracked Python merge conflict before execution: {name}"
                )
            paths.append(root / name)
        if Path(__file__).resolve() not in paths:
            raise ValueError("Tracked code inventory must include the executing stage.py module")
    return canonical_sha256(
        {path.relative_to(package).as_posix(): file_sha256(path) for path in sorted(paths)}
    )


def migration_report(family: FamilyInput, corrected: Path, destination: Path) -> dict[str, int]:
    """Join native rows on disk; separate family, method, threshold and effect filters."""
    counts = {
        "native_reject_005": 0,
        "joint_bh_reject_005": 0,
        "by_reject_005": 0,
        "by_reject_0001": 0,
        "native_reported_hits": 0,
    }
    with (
        tempfile.TemporaryDirectory(prefix="mtp-migration-") as scratch,
        closing(sqlite3.connect(Path(scratch) / "native.sqlite")) as db,
        db,
    ):
        db.execute("CREATE TABLE native (id TEXT PRIMARY KEY, content TEXT NOT NULL)")
        for native in family.observations():
            value = {
                "native_q": native.native_q,
                "native_method": native.native_method,
                "native_family": native.native_family,
                "native_threshold": native.native_threshold,
                "source_row_sha256": canonical_sha256(native.row),
                "source_p": native.raw.p_value,
                "source_status": native.raw.status,
                "native_reported_hit": native.row.get("is_significant", ""),
            }
            db.execute(
                "INSERT INTO native VALUES (?, ?)", (native.raw.hypothesis_id, json_text(value))
            )
        with destination.open("x", encoding="utf-8", newline="") as handle:
            fields = [
                "hypothesis_id",
                "source_row_sha256",
                "source_p",
                "source_status",
                "native_family",
                "native_method",
                "native_q",
                "native_threshold",
                "native_reported_hit",
                "joint_bh_q",
                "by_q",
                "reject_by_005",
                "reject_by_0001",
                "family_changed",
            ]
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for (identity, content), result in zip(
                db.execute("SELECT id, content FROM native ORDER BY id"),
                read_corrected_inventory(corrected),
                strict=True,
            ):
                native_row = json.loads(content)
                if (
                    identity != result.hypothesis_id
                    or native_row["source_p"] != result.p_value
                    or native_row["source_status"] != result.status
                ):
                    raise ValueError("Native/sidecar report identity or p-value mismatch")
                tested = result.status == "tested"
                native_q = native_row["native_q"]
                counts["native_reject_005"] += tested and native_q <= 0.05
                counts["joint_bh_reject_005"] += result.reject_bh_005 is True
                counts["by_reject_005"] += result.reject_by_005 is True
                counts["by_reject_0001"] += result.reject_by_0001 is True
                counts["native_reported_hits"] += (
                    native_row["native_reported_hit"].lower() == "true"
                )
                writer.writerow(
                    {
                        "hypothesis_id": identity,
                        **native_row,
                        "joint_bh_q": result.q_bh,
                        "by_q": result.q_by,
                        "reject_by_005": result.reject_by_005,
                        "reject_by_0001": result.reject_by_0001,
                        "family_changed": family.metadata.family_id == "ilyas-joint-ora",
                    }
                )
    return counts


def _protect_sources(output: Path, roots: list[Path]) -> None:
    resolved = output.resolve()
    if output.exists() or output.is_symlink():
        raise FileExistsError("Multiplicity output must be a new directory")
    for root in roots:
        root = root.resolve()
        if resolved.is_relative_to(root) or root.is_relative_to(resolved):
            raise ValueError("Multiplicity sidecars must be disjoint from verified source bundles")


def run_post_arm(
    comparison: Path, output: Path, *, recorded_at: datetime | None = None
) -> PostArmMultiplicityRun:
    """Import unchanged native p-values and emit independently verifiable sidecars."""
    implementation = code_sha256()
    policy = load_policy().model_copy(deep=True)
    comparison_manifest = comparison / "comparison_run.json" if comparison.is_dir() else comparison
    comparison_digest = file_sha256(comparison_manifest)
    sources = bind_recorded_hashes(open_verified_sources(comparison), comparison_digest)
    instant = recorded_at if recorded_at is not None else datetime.now(UTC)
    _protect_sources(
        output,
        [
            Path(sources.comparison.manifest_path).parent,
            Path(sources.core.manifest_path).parent,
            *(Path(arm.manifest_path).parent for arm in sources.arms.values()),
        ],
    )
    inputs = [
        noel_family(sources),
        jed_family(sources),
        native_family(sources),
        ora_family(sources),
    ]
    output.mkdir(parents=True)
    corrections: list[ExternalFamilyCorrection] = []
    reports: list[SourceArtifact] = []
    summaries = []
    consumed = {
        artifact.identity: artifact for family in inputs for artifact in family.metadata.sources
    }
    observed_files = {
        path: sources.expected_digest(path) for family in inputs for path in family.source_paths
    }
    for family in inputs:
        relative = family.metadata.family_id
        destination = output / relative
        correction = correct_inventory(
            family.metadata,
            family.declarations(),
            (observation.raw for observation in family.observations()),
            destination,
            recorded_at=instant,
        )
        report_path = destination / "migration.tsv"
        counts = migration_report(family, destination / correction.corrected.identity, report_path)
        summaries.append(
            {
                "family": relative,
                "retained_features": family.metadata.retained_feature_count,
                "declared_tests": family.metadata.expected_family_size,
                **counts,
            }
        )
        for artifact in (correction.declared, correction.raw, correction.corrected):
            artifact.identity = f"{relative}/{artifact.identity}"
        corrections.append(correction)
        reports.append(
            SourceArtifact(identity=f"{relative}/migration.tsv", sha256=file_sha256(report_path))
        )
    summary_path = output / "migration.json"
    summary_path.write_text(
        json.dumps(
            {
                "families": summaries,
                "family_change": (
                    "Only ORA changes membership: separate native tails become "
                    "one joint family. Compare native BH with joint BH at equal q."
                ),
                "method_change": "Compare joint BH and BY at q<=0.05 within each family.",
                "threshold_change": "Compare BY q<=0.05 with BY q<=0.001 within each family.",
                "effects": (
                    "Effect/fold display filters are separate from multiplicity rejection. "
                    "Native reported ORA hits also require their fold cutoff."
                ),
                "precision": (
                    "Noel p/q exports use six significant digits; tiny differences "
                    "from native q may reflect export rounding."
                ),
                "scope": (
                    "Additive diagnostics only; no cross-family/arm error-control claim, "
                    "native overwrite, focal q reuse or selection-bias correction."
                ),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    reports.append(SourceArtifact(identity=summary_path.name, sha256=file_sha256(summary_path)))
    for path, digest in observed_files.items():
        if file_sha256(path) != digest or source(path).identity not in consumed:
            raise ValueError("Source changed from recorded hash while sidecars were being written")
    if code_sha256() != implementation:
        raise ValueError("Multiplicity implementation code changed during execution")
    if load_policy() != policy:
        raise ValueError("Multiplicity policy changed during execution")
    record = PostArmMultiplicityRun(
        version="1.1.0",
        recorded_at=instant,
        comparison_manifest_path=sources.comparison.manifest_path,
        comparison_manifest_sha256=comparison_digest,
        core_run_id=sources.core.run_id,
        contrast_population=sources.comparison.contrast_population,
        direction_basis=sources.comparison.direction_basis,
        strain_rollup=sources.scope,
        policy=policy,
        code_sha256=implementation,
        scipy_version=scipy.__version__,
        numpy_version=np.__version__,
        sources=sorted(consumed.values(), key=lambda artifact: artifact.identity),
        families=corrections,
        reports=reports,
    )
    (output / "multiplicity_run.json").write_text(
        record.model_dump_json(indent=2) + "\n", encoding="utf-8"
    )
    return record


def verify_post_arm(path: Path) -> PostArmMultiplicityRun:
    """Rebuild from externally verified native sources and compare every record/hash."""
    manifest = path / "multiplicity_run.json" if path.is_dir() else path
    root = manifest.parent.resolve()
    record = PostArmMultiplicityRun.model_validate_json(manifest.read_text(encoding="utf-8"))
    if record.version != "1.1.0" or record.policy != load_policy():
        raise ValueError("Unsupported multiplicity policy or stage version")
    artifacts = [
        *record.reports,
        *(
            artifact
            for family in record.families
            for artifact in (family.declared, family.raw, family.corrected)
        ),
    ]
    expected = {"multiplicity_run.json", *(artifact.identity for artifact in artifacts)}
    if len(expected) != len(artifacts) + 1:
        raise ValueError("Duplicate multiplicity output identities")
    actual = {item.relative_to(root).as_posix() for item in root.rglob("*") if item.is_file()}
    if actual != expected or any(item.is_symlink() for item in root.rglob("*")):
        raise ValueError("Multiplicity output inventory mismatch")
    for artifact in artifacts:
        target = root / artifact.identity
        if not target.resolve().is_relative_to(root) or file_sha256(target) != artifact.sha256:
            raise ValueError("Multiplicity artifact hash mismatch")
    with tempfile.TemporaryDirectory(prefix="mtp-verify-multiplicity-") as scratch:
        replay = run_post_arm(
            Path(record.comparison_manifest_path),
            Path(scratch) / "replay",
            recorded_at=record.recorded_at,
        )
    if record != replay:
        raise ValueError(
            "Multiplicity replay differs: source, metadata, policy, code or result tampering"
        )
    return record
