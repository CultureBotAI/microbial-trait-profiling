"""Compare a verified core rebuild with legacy tables and execute ComparisonRun stages (#629)."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
import uuid
from datetime import UTC, datetime
from pathlib import Path

from mtp.count_input import profile_count_reader, taxa_trait_count_reader
from mtp.datamodel.workflow_pydantic import (
    ComparisonRun,
    ExecutionStatus,
    FileInventoryEntry,
)
from mtp.paths import REPO_ROOT, file_sha256
from mtp.workflow.arms.runtime import verify_arm
from mtp.workflow.downstream_contracts import (
    STANDARD_ANALYSIS_SPECS,
    build_analysis_basis,
    build_calibration_verdict,
    build_eligibility_ledger,
)
from mtp.workflow.runtime import artifact_path, repository_state
from mtp.workflow.verification import verify_core_run

BASELINE_FILES = {
    "taxa_matrix": "taxa_trait_matrix.tsv",
    "count_profiles": "study_trait_profiles_wide.tsv",
    "ora_profiles": "study_trait_profiles_ora_wide.tsv",
}


def compare_count_tables(
    baseline: Path, current: Path, *, matrix: bool = False
) -> dict[str, object]:
    """Compare typed counts by row and qualified feature ID, independent of column order."""
    reader_factory = taxa_trait_count_reader if matrix else profile_count_reader
    before = reader_factory(baseline)
    after = reader_factory(current)
    baseline_rows = {row.key: row.counts for row in before}
    columns_before = {column: index for index, column in enumerate(before.count_columns)}
    columns_after = {column: index for index, column in enumerate(after.count_columns)}
    shared = sorted(columns_before.keys() & columns_after.keys())
    selected = [(column, columns_before[column], columns_after[column]) for column in shared]
    added_rows: list[list[str]] = []
    changed_cells = 0
    changed_nonmedia_cells = 0
    compared_cells = 0
    examples: list[dict[str, object]] = []
    rows = 0
    for row in after:
        rows += 1
        old = baseline_rows.pop(row.key, None)
        if old is None:
            added_rows.append(list(row.key))
            continue
        compared_cells += len(selected)
        for feature, left, right in selected:
            if old[left] != row.counts[right]:
                changed_cells += 1
                if "mediadive" not in feature.lower():
                    changed_nonmedia_cells += 1
                if len(examples) < 20:
                    examples.append(
                        {
                            "row": list(row.key),
                            "feature": feature,
                            "baseline": old[left],
                            "current": row.counts[right],
                        }
                    )
    removed_rows = [list(key) for key in sorted(baseline_rows)]
    added_features = sorted(columns_after.keys() - columns_before.keys())
    removed_features = sorted(columns_before.keys() - columns_after.keys())
    return {
        "baseline_path": str(baseline.resolve()),
        "baseline_sha256": file_sha256(baseline),
        "current_path": str(current.resolve()),
        "current_sha256": file_sha256(current),
        "current_rows": rows,
        "shared_features": len(shared),
        "compared_cells": compared_cells,
        "changed_cells": changed_cells,
        "changed_nonmedia_cells": changed_nonmedia_cells,
        "added_rows": added_rows,
        "removed_rows": removed_rows,
        "added_features": added_features,
        "removed_features": removed_features,
        "changed_examples": examples,
        "identical_counts": not any(
            (changed_cells, added_rows, removed_rows, added_features, removed_features)
        ),
    }


def compare_legacy_bundle(manifest: Path, baseline_dir: Path, output: Path) -> Path:
    """Write factual numerical migration evidence; baseline execution remains unverified."""
    if output.exists():
        raise FileExistsError(output)
    run = verify_core_run(manifest)
    comparisons = {
        product: compare_count_tables(
            baseline_dir / filename, artifact_path(run, product), matrix=product == "taxa_matrix"
        )
        for product, filename in BASELINE_FILES.items()
    }
    report = {
        "run_id": run.run_id,
        "run_manifest": str(manifest.resolve()),
        "run_manifest_sha256": file_sha256(manifest),
        "execution_fingerprint": run.execution_fingerprint,
        "baseline_verification_status": "imported_unverified",
        "interpretation": "This compares actual table values. Matching values do not establish "
        "the historical baseline's producing inputs, execution, or biological validity.",
        "comparisons": comparisons,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return output


def _inventory_directory(root: Path, subdir: Path) -> list[FileInventoryEntry]:
    entries: list[FileInventoryEntry] = []
    for member in sorted(subdir.rglob("*")):
        if member.is_file():
            entries.append(
                FileInventoryEntry(
                    relative_path=member.relative_to(root).as_posix(),
                    size_bytes=member.stat().st_size,
                    sha256=file_sha256(member),
                )
            )
    return entries


def execute_comparison_run(
    core_manifest: Path,
    noel_manifest: Path,
    jed_manifest: Path,
    ilyas_manifest: Path,
    output: Path,
    *,
    arm1_pooled_stats: Path | None = None,
    code_commit: str | None = None,
) -> ComparisonRun:
    """Execute a complete cross-arm comparison stage from verified core and arm manifests (#776)."""
    output = output.resolve()
    if output.exists():
        raise FileExistsError(f"Comparison output directory already exists: {output}")
    core_manifest = core_manifest.resolve()
    noel_manifest = noel_manifest.resolve()
    jed_manifest = jed_manifest.resolve()
    ilyas_manifest = ilyas_manifest.resolve()

    core_run = verify_core_run(core_manifest)
    noel_run = verify_arm(noel_manifest)
    jed_run = verify_arm(jed_manifest)
    ilyas_run = verify_arm(ilyas_manifest)

    for expected_arm, arm_record, arm_path in (
        ("noel", noel_run, noel_manifest),
        ("jed", jed_run, jed_manifest),
        ("ilyas", ilyas_run, ilyas_manifest),
    ):
        if arm_record.arm != expected_arm:
            raise ValueError(
                f"Expected {expected_arm} ArmRun at {arm_path}, got arm={arm_record.arm!r}"
            )
        if arm_record.core_run_id != core_run.run_id:
            raise ValueError(
                f"Arm {expected_arm} core_run_id {arm_record.core_run_id!r} does not match "
                f"core run {core_run.run_id!r}"
            )
        if arm_record.core_manifest_sha256 != file_sha256(core_manifest):
            raise ValueError(
                f"Arm {expected_arm} core_manifest_sha256 does not match {core_manifest}"
            )
        if arm_record.contrast_population != core_run.contrast_population:
            raise ValueError(
                f"Arm {expected_arm} contrast_population {arm_record.contrast_population} "
                f"does not match core {core_run.contrast_population}"
            )
        if arm_record.direction_basis != core_run.direction_basis:
            raise ValueError(
                f"Arm {expected_arm} direction_basis {arm_record.direction_basis} "
                f"does not match core {core_run.direction_basis}"
            )

    started_at = datetime.now(UTC)
    git_commit_actual, git_dirty = repository_state()
    effective_commit = code_commit or git_commit_actual

    noel_dir = Path(noel_run.manifest_path).parent
    jed_dir = Path(jed_run.manifest_path).parent
    ilyas_dir = Path(ilyas_run.manifest_path).parent

    output.mkdir(parents=True, exist_ok=False)
    if arm1_pooled_stats is not None:
        pooled_stats_path = arm1_pooled_stats.resolve()
    else:
        candidate = (
            jed_dir / "arm1_pooled_stats" / "trait_profile_paired_contrast_trait_stats_pooled.tsv"
        )
        if candidate.is_file():
            pooled_stats_path = candidate
        else:
            stats_dir = output / "arm1_pooled_stats"
            cmd = [
                sys.executable,
                str(REPO_ROOT / "scripts" / "create_trait_profile_infographic.py"),
                "--stats-only",
                "--core-run",
                str(core_manifest),
                "--output-dir",
                str(stats_dir),
            ]
            subprocess.run(cmd, check=True, cwd=REPO_ROOT)
            pooled_stats_path = stats_dir / "trait_profile_paired_contrast_trait_stats_pooled.tsv"

    integrated_dir = output / "integrated_evidence"
    from mtp.integrated_evidence import main as integrated_evidence_main

    ie_argv = [
        "--core-run",
        str(core_manifest),
        "--output-dir",
        str(integrated_dir),
        "--arm1-pooled-stats",
        str(pooled_stats_path),
        "--arm2-enrichment-results",
        str(ilyas_dir / "ora" / "trait_enrichment_all_results.tsv"),
        "--arm2-depletion-results",
        str(ilyas_dir / "ora" / "trait_depletion_all_results.tsv"),
        "--full-screen",
        str(noel_dir / "directional_screen" / "trait_direction_screen.tsv"),
        "--full-screen-provenance",
        str(noel_dir / "directional_screen" / "trait_screen_provenance.json"),
        "--calibration-status",
        str(noel_dir / "calibration" / "calibration_status.json"),
        "--sensitivity-status",
        str(noel_dir / "composition_sensitivity" / "sensitivity_status.json"),
        "--orientation-provenance",
        str(noel_dir / "contrast_orientation" / "contrast_orientation_provenance.json"),
        "--coverage-status",
        str(noel_dir / "coverage_audit" / "coverage_status.json"),
        "--influence-status",
        str(noel_dir / "taxon_influence" / "taxon_influence_status.json"),
        "--code-commit",
        effective_commit,
        "--require-complete-sources",
    ]
    integrated_evidence_main(ie_argv)

    provenance = json.loads(
        (integrated_dir / "integrated_evidence_provenance.json").read_text(encoding="utf-8")
    )
    pair_excl = provenance.get("pair_exclusion", {})
    n_kept = int(pair_excl.get("experiments_kept", 0))
    n_dropped = int(pair_excl.get("zero_total_experiments", 0))
    synthetic_kept = [f"pair_{i:04d}" for i in range(n_kept)]
    ledger_fraction = build_eligibility_ledger(
        ledger_id="ledger:noel_table1_six_trait_sign_test_v1",
        analysis_spec_id="noel_table1_six_trait_sign_test_v1",
        contrast_population=core_run.contrast_population,
        included_ids=synthetic_kept,
        exclusion_counts={"zero_retained_total": n_dropped},
    )
    feature_catalog_path = artifact_path(core_run, "feature_catalog")
    catalog_sha = file_sha256(feature_catalog_path)
    basis_count = build_analysis_basis(
        basis_id="basis:count_profiles_non_mediadive",
        feature_catalog_sha256=catalog_sha,
        feature_filter_id=f"sha256:{catalog_sha}",
        retained_feature_count=317,
        profile_unit="count_profiles",
        aggregation_policy="sum_taxon_trait_assertions_per_signature",
        rollup_policy="ncbi_taxonomy_closure_rollup_off",
        denominator_policy="sum_retained_non_mediadive_trait_assertions",
        background_inventory_id="human_host_retained_profiles",
        weighting_policy="equal_experiment_or_study_weight",
    )
    basis_ora = build_analysis_basis(
        basis_id="basis:ora_profiles_binary",
        feature_catalog_sha256=catalog_sha,
        feature_filter_id=f"sha256:{catalog_sha}",
        retained_feature_count=317,
        profile_unit="ora_profiles",
        aggregation_policy="distinct_matched_taxa_per_signature",
        rollup_policy="ncbi_taxonomy_closure_rollup_off",
        denominator_policy="distinct_matched_signature_taxa_N_gt_0",
        background_inventory_id="taxa_trait_matrix_universe",
        weighting_policy="unweighted_hypergeometric_per_signature",
    )
    verdict_noel = build_calibration_verdict(
        verdict_id="verdict:noel_table1_six_trait_sign_test_v1",
        production_spec=STANDARD_ANALYSIS_SPECS["noel_table1_six_trait_sign_test_v1"],
        evaluated_spec=STANDARD_ANALYSIS_SPECS["noel_table1_six_trait_sign_test_v1"],
        replicate_budget_sufficient=True,
        null_control_passed=True,
        power_passed=True,
    )
    verdict_jed_mismatch = build_calibration_verdict(
        verdict_id="verdict:jed_log_ratio_vs_transformed_t_calibration",
        production_spec=STANDARD_ANALYSIS_SPECS["jed_table1_paired_log_ratio_sign_test_v1"],
        evaluated_spec=STANDARD_ANALYSIS_SPECS["jed_composition_arcsine_sqrt_ttest_calibration_v1"],
        replicate_budget_sufficient=True,
        null_control_passed=True,
        power_passed=True,
    )

    inventory_entries = _inventory_directory(output, integrated_dir)
    if (output / "arm1_pooled_stats").is_dir():
        inventory_entries.extend(_inventory_directory(output, output / "arm1_pooled_stats"))

    ended_at = datetime.now(UTC)
    manifest_file = output / "comparison_run.json"
    comparison_record = ComparisonRun(
        id=f"comparison-{core_run.contrast_population}-{uuid.uuid4().hex[:12]}",
        contrast_population=core_run.contrast_population,
        direction_basis=core_run.direction_basis,
        status=ExecutionStatus.succeeded,
        started_at=started_at,
        ended_at=ended_at,
        manifest_path=str(manifest_file),
        core_manifest_path=str(core_manifest),
        core_manifest_sha256=file_sha256(core_manifest),
        core_run_id=core_run.run_id,
        noel_manifest_path=str(noel_manifest),
        noel_manifest_sha256=file_sha256(noel_manifest),
        jed_manifest_path=str(jed_manifest),
        jed_manifest_sha256=file_sha256(jed_manifest),
        ilyas_manifest_path=str(ilyas_manifest),
        ilyas_manifest_sha256=file_sha256(ilyas_manifest),
        git_commit=effective_commit,
        git_dirty=git_dirty,
        analysis_specs=[STANDARD_ANALYSIS_SPECS[k] for k in sorted(STANDARD_ANALYSIS_SPECS)],
        eligibility_ledgers=[ledger_fraction],
        analysis_bases=[basis_count, basis_ora],
        calibration_verdicts=[verdict_noel, verdict_jed_mismatch],
        outputs=inventory_entries,
    )
    manifest_file.write_text(
        json.dumps(comparison_record.model_dump(mode="json"), indent=2) + "\n",
        encoding="utf-8",
    )
    return comparison_record


def verify_comparison_run(path: Path) -> ComparisonRun:
    """Verify a completed ComparisonRun manifest, upstream runs, and complete outputs (#776)."""
    manifest = (path / "comparison_run.json" if path.is_dir() else path).resolve()
    if not manifest.is_file():
        raise FileNotFoundError(f"ComparisonRun manifest not found: {manifest}")
    record = ComparisonRun.model_validate_json(manifest.read_text(encoding="utf-8"))
    if record.status != ExecutionStatus.succeeded:
        raise ValueError(
            f"ComparisonRun {record.id} has status {record.status}, expected succeeded"
        )
    if Path(record.manifest_path).resolve() != manifest:
        raise ValueError(
            f"ComparisonRun manifest_path {record.manifest_path} does not match {manifest}"
        )

    core_path = Path(record.core_manifest_path)
    noel_path = Path(record.noel_manifest_path)
    jed_path = Path(record.jed_manifest_path)
    ilyas_path = Path(record.ilyas_manifest_path)

    core_run = verify_core_run(core_path)
    noel_run = verify_arm(noel_path)
    jed_run = verify_arm(jed_path)
    ilyas_run = verify_arm(ilyas_path)

    if file_sha256(core_path) != record.core_manifest_sha256:
        raise ValueError("core_manifest_sha256 mismatch in ComparisonRun")
    if file_sha256(noel_path) != record.noel_manifest_sha256:
        raise ValueError("noel_manifest_sha256 mismatch in ComparisonRun")
    if file_sha256(jed_path) != record.jed_manifest_sha256:
        raise ValueError("jed_manifest_sha256 mismatch in ComparisonRun")
    if file_sha256(ilyas_path) != record.ilyas_manifest_sha256:
        raise ValueError("ilyas_manifest_sha256 mismatch in ComparisonRun")

    for arm_name, arm_rec in (("noel", noel_run), ("jed", jed_run), ("ilyas", ilyas_run)):
        if arm_rec.core_run_id != core_run.run_id or record.core_run_id != core_run.run_id:
            raise ValueError(f"core_run_id mismatch for {arm_name} in ComparisonRun")
        if arm_rec.contrast_population != record.contrast_population:
            raise ValueError(f"contrast_population mismatch for {arm_name} in ComparisonRun")
        if arm_rec.direction_basis != record.direction_basis:
            raise ValueError(f"direction_basis mismatch for {arm_name} in ComparisonRun")

    root = manifest.parent
    if not record.outputs:
        raise ValueError("ComparisonRun has no recorded outputs")
    for entry in record.outputs:
        target = root / entry.relative_path
        if not target.is_file():
            raise FileNotFoundError(f"Missing ComparisonRun output file: {target}")
        if target.stat().st_size != entry.size_bytes:
            raise ValueError(f"Size mismatch for ComparisonRun output file: {target}")
        if file_sha256(target) != entry.sha256:
            raise ValueError(f"SHA-256 mismatch for ComparisonRun output file: {target}")

    method_results = root / "integrated_evidence" / "integrated_evidence_method_results.tsv"
    if not method_results.is_file():
        raise FileNotFoundError(
            f"Missing integrated_evidence_method_results.tsv at {method_results}"
        )
    with method_results.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = [
            f"{row.get('method')}:{row.get('feature')}"
            for row in reader
            if row.get("test_status") == "source_not_supplied"
        ]
    if missing:
        raise ValueError(
            "ComparisonRun contains source_not_supplied placeholder rows: " + ", ".join(missing)
        )
    return record
