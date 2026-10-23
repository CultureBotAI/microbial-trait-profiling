"""Compare a verified core rebuild with legacy tables and execute ComparisonRun stages (#629)."""

from __future__ import annotations

import asyncio
import csv
import json
import sys
import uuid
from collections.abc import Mapping, Sequence
from datetime import UTC, datetime
from pathlib import Path

from mtp.count_input import profile_count_reader, taxa_trait_count_reader
from mtp.datamodel.workflow_pydantic import (
    AnalysisBasis,
    AnalysisSpec,
    ArmRun,
    CalibrationVerdict,
    ComparisonRun,
    ContrastPopulation,
    DirectionBasis,
    EligibilityLedger,
    ExecutionStatus,
    FileInventoryEntry,
    StorageKind,
    WorkflowRun,
)
from mtp.paths import REPO_ROOT, file_sha256
from mtp.trait_guilds import OXYGEN_GUILDS, TABLE1_TRAITS, drop_zero_total_pairs, load_profiles
from mtp.workflow.arms.runtime import verify_arm
from mtp.workflow.coverage import feature_class, feature_filter_id
from mtp.workflow.downstream_contracts import (
    STANDARD_ANALYSIS_SPECS,
    build_analysis_basis,
    build_calibration_verdict,
    build_eligibility_ledger,
)
from mtp.workflow.operations import run_strain_rollup
from mtp.workflow.runtime import artifact_path, repository_state
from mtp.workflow.verification import verify_core_run

BASELINE_FILES = {
    "taxa_matrix": "taxa_trait_matrix.tsv",
    "count_profiles": "study_trait_profiles_wide.tsv",
    "ora_profiles": "study_trait_profiles_ora_wide.tsv",
}

#: Each integrated_evidence input an arm produces: the option, then the arm, step and
#: file name that arm's arm_run.json records (#814). No arm produces the pooled Arm 1
#: statistics, so the stage computes them; nothing produces the #184 coverage status,
#: so it is left out and integrated_evidence records it as not supplied.
ARM_EVIDENCE_INPUTS: tuple[tuple[str, str, str, str], ...] = (
    ("--arm2-enrichment-results", "ilyas", "enrichment", "trait_enrichment_all_results.tsv"),
    ("--arm2-depletion-results", "ilyas", "depletion", "trait_depletion_all_results.tsv"),
    ("--full-screen", "noel", "directional_screen", "trait_direction_screen.tsv"),
    ("--full-screen-provenance", "noel", "directional_screen", "trait_screen_provenance.json"),
    ("--calibration-status", "noel", "empirical_calibration", "calibration_status.json"),
    (
        "--orientation-provenance",
        "noel",
        "contrast_orientation",
        "contrast_orientation_provenance.json",
    ),
    ("--influence-status", "noel", "native_taxon_influence", "native_manifest.json"),
    ("--sensitivity-status", "jed", "composition_sensitivity", "sensitivity_status.json"),
)

#: The arm, step and file each calibration verdict is read from (#815).
CALIBRATION_SOURCES: tuple[tuple[str, str, str], ...] = (
    ("noel", "empirical_calibration", "calibration_status.json"),
    ("jed", "composition_calibration", "composition_calibration_summary.json"),
    ("jed", "composition_calibration", "composition_calibration_replicates.tsv"),
)

#: The procedure in Noel's calibration harness that evaluated each production spec (#815).
NOEL_CALIBRATED_PROCEDURES = {"noel_table1_six_trait_sign_test_v1": "fraction_sign_test"}

#: The protocol's negative control. It is declared in advance to fail the acceptance
#: criteria, since a criterion it passes is too weak to be evidence about any procedure
#: (methods/calibration.md); its failure is the expected outcome.
NOEL_NEGATIVE_CONTROL = "screened_fraction_sign_test"

#: Feature selection on the analysed corpus that a production spec's calibration leaves
#: open, recorded beside its verdict rather than in its label (#855).
SELECTION_UNCERTAINTY = {
    "noel_table1_six_trait_sign_test_v1": (
        "The six Table 1 traits were selected on this corpus before testing, so their "
        "p-values are exploratory. The negative control selects and then tests, but its "
        "inflation does not bound the historical selection's cost in either direction "
        "(#363). This verdict covers the sign test on pre-specified traits and does not "
        "validate the post-hoc six-trait analysis."
    )
}

#: The one label rule for selection. The declared protocol leaves the selection out of
#: the label, so a production verdict is anti-conservative only for its own failed null
#: criterion; True would label every selected spec anti-conservative whenever the
#: negative control fails, an owner decision (#855).
SELECTION_LABELS_ANTI_CONSERVATIVE = False

#: Noel's criteria that judge null error. Its power criterion only checks that power
#: is reported precisely; the protocol sets no threshold on power itself.
NULL_CRITERIA = ("per_test_level", "family_error")


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
                    path=member.relative_to(root).as_posix(),
                    size_bytes=member.stat().st_size,
                    sha256=file_sha256(member),
                )
            )
    return entries


async def _exec_command(cmd: list[str], cwd: Path) -> None:
    process = await asyncio.create_subprocess_exec(*cmd, cwd=cwd)
    returncode = await process.wait()
    if returncode != 0:
        raise RuntimeError(f"Command failed with exit code {returncode}: {cmd}")


def recorded_output(run: ArmRun, step: str, name: str) -> Path:
    """The file named ``name`` that ``step`` of a verified arm run recorded (#814)."""
    for item in run.steps or []:
        if item.id != step:
            continue
        matches = [
            Path(output.path)
            for output in item.outputs or []
            if output.storage_kind == StorageKind.file and Path(output.path).name == name
        ]
        if len(matches) != 1:
            raise ValueError(
                f"{run.arm} step {step!r} recorded {len(matches)} files named {name!r}: "
                f"{run.manifest_path}"
            )
        return matches[0]
    raise ValueError(f"{run.arm} run records no {step!r} step: {run.manifest_path}")


def _noel_procedure(statuses: object, procedure: str, path: Path) -> tuple[str, dict[str, str]]:
    """One procedure's overall status and criteria, checked against Noel's own rule (#844)."""
    from mtp.calibration import calibration_status

    entry = statuses.get(procedure) if isinstance(statuses, dict) else None
    overall = entry.get("overall") if isinstance(entry, dict) else None
    criteria = entry.get("criteria") if isinstance(entry, dict) else None
    if not isinstance(overall, str) or not isinstance(criteria, dict) or not criteria:
        raise ValueError(f"{path} records no overall status and criteria for {procedure!r}")
    verdicts = {str(name): str(value) for name, value in criteria.items()}
    rows = [
        {"procedure": procedure, "criterion": name, "verdict": value}
        for name, value in verdicts.items()
    ]
    derived = calibration_status(rows)[procedure]["overall"]
    if derived != overall:
        raise ValueError(
            f"{path} records {procedure!r} as {overall!r}, but its criteria give {derived!r}"
        )
    return overall, verdicts


def noel_calibration_verdict(path: Path, spec: AnalysisSpec) -> CalibrationVerdict:
    """Noel's verdict for ``spec``'s procedure on pre-specified traits (#815, #844, #855).

    The production procedure's own null criteria decide its error control. A failed one
    is anti-conservative at any budget, since the harness keeps a short run's fail only
    when its estimate already exceeds the tolerance. They count as passed only where the
    negative control failed them as declared; where it passed them they are too weak, and
    certification is withheld. The control's outcome and the features' selection on this
    corpus are recorded in their own fields; whether that selection enters the label is
    the one rule in ``SELECTION_LABELS_ANTI_CONSERVATIVE``. Noel's overall status maps
    explicitly: not_certified is an unmet replicate budget, and a failed power criterion
    is failed power. The protocol sets no threshold on power itself (#850).
    """
    content = json.loads(path.read_text(encoding="utf-8"))
    statuses = content.get("calibration_status") if isinstance(content, dict) else None
    procedure = NOEL_CALIBRATED_PROCEDURES[spec.id]
    overall, criteria = _noel_procedure(statuses, procedure, path)
    control, control_criteria = _noel_procedure(statuses, NOEL_NEGATIVE_CONTROL, path)
    # A passing sanity check is not evidence of calibration (#360); a failing one is.
    judged = [
        c for c in NULL_CRITERIA if c in criteria and not criteria[c].endswith("(sanity check)")
    ]
    too_weak = [c for c in judged if control_criteria.get(c) == "pass"]
    as_expected: bool | None = False if too_weak else None
    if judged and all(control_criteria.get(c) == "fail" for c in judged):
        as_expected = True
    null_control_passed: bool | None = None
    if "fail" in (criteria.get(c) for c in NULL_CRITERIA):
        null_control_passed = False
    elif as_expected and all(criteria[c] == "pass" for c in judged):
        null_control_passed = True
    selection = SELECTION_UNCERTAINTY.get(spec.id)
    selection_labelled = (
        SELECTION_LABELS_ANTI_CONSERVATIVE and bool(selection) and control == "fail"
    )
    verdict = build_calibration_verdict(
        verdict_id=f"verdict:{spec.id}",
        production_spec=spec,
        evaluated_spec=spec,
        replicate_budget_sufficient=not any(
            value.startswith(("insufficient_budget", "not_evaluated"))
            for value in (*criteria.values(), *control_criteria.values())
        ),
        null_control_passed=null_control_passed,
        power_passed=False if criteria.get("power_reported") == "fail" else None,
        anti_conservative_detected=selection_labelled,
    )
    read = "; ".join(
        f"{name} {status} [{', '.join(f'{c}={v}' for c, v in sorted(found.items()))}]"
        for name, status, found in (
            (procedure, overall, criteria),
            (NOEL_NEGATIVE_CONTROL, control, control_criteria),
        )
    )
    summary = f"{verdict.summary_message} Read from {path.name}: {read}."
    if as_expected is False:
        summary += (
            f" The negative control passed {', '.join(too_weak)}, which the protocol declares "
            "it must fail, so those criteria are too weak to be evidence and certification "
            "is withheld."
        )
    elif as_expected:
        summary += (
            f" The negative control failed {', '.join(judged)} as the protocol declares: in "
            "this harness selecting and then testing exceeds the tolerance, while "
            f"{procedure}'s null result on pre-specified traits stands as recorded."
        )
    if selection_labelled:
        summary += " The label counts the selection of the features as anti-conservative."
    if selection:
        summary += f" {selection}"
    summary += " Power is reported there without a threshold (#850)."
    return verdict.model_copy(
        update={
            "summary_message": summary,
            "negative_control_procedure": NOEL_NEGATIVE_CONTROL,
            "negative_control_outcome": control,
            "negative_control_as_expected": as_expected,
            "selection_uncertainty": selection,
        }
    )


def _composition_replicates(path: Path, summary: Mapping[str, object]) -> int:
    """The complete replicates in the composition harness's replicate table (#844, #856).

    Each replicate counts only with exactly one row per declared correction and threshold,
    carrying the summary's seed, stream, unit and family size, and the replicate numbers
    must run from 0 without gaps. The summary's per-cell counts must match the rows. A
    truncated, duplicated or inconsistent table is an error, never a shorter count.
    """
    from mtp.composition_calibration.composition import CORRECTIONS, THRESHOLDS

    cells = {(correction, threshold) for correction in CORRECTIONS for threshold in THRESHOLDS}
    identity = tuple(str(summary.get(key)) for key in ("seed", "stream", "unit", "n_traits"))
    found: dict[int, dict[tuple[str, float], int]] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            try:
                replicate, threshold = int(row["replicate"]), float(row["threshold"])
                significant = int(row["significant"])
                producer = (row["seed"], row["stream"], row["unit"], row["n_tests"])
            except (KeyError, TypeError, ValueError) as error:
                raise ValueError(f"{path} has a malformed replicate row: {row}") from error
            cell = (row["correction"], threshold)
            if producer != identity or cell not in cells or cell in found.get(replicate, {}):
                raise ValueError(
                    f"{path} has an inconsistent, undeclared or duplicate row {cell} for "
                    f"replicate {replicate}"
                )
            found.setdefault(replicate, {})[cell] = significant
    complete = len(found)
    if sorted(found) != list(range(complete)) or any(set(row) != cells for row in found.values()):
        raise ValueError(f"{path} is truncated or numbers its replicates outside 0..{complete - 1}")
    groups = summary.get("summary")
    recorded = {
        (group.get("correction"), group.get("threshold")): group
        for group in (groups if isinstance(groups, list) else [])
        if isinstance(group, Mapping)
    }
    if set(recorded) != cells or any(
        recorded[cell].get("replicates") != complete
        or recorded[cell].get("replicates_with_any_call")
        != sum(row[cell] > 0 for row in found.values())
        for cell in cells
    ):
        raise ValueError(f"The summary's per-cell replicate counts disagree with {path.name}")
    return complete


def jed_calibration_verdict(
    summary: Path,
    replicates: Path,
    production_spec: AnalysisSpec,
    evaluated_spec: AnalysisSpec,
) -> CalibrationVerdict:
    """The composition harness's replicate budget; it judges neither error nor power (#815).

    The summary always reports every requested replicate as completed, so completion is
    read from the validated replicate table and must agree with it (#844, #856).
    """
    from mtp.composition_calibration.composition import DEFAULT_REPLICATES

    content = json.loads(summary.read_text(encoding="utf-8"))
    if not isinstance(content, dict):
        raise ValueError(f"{summary} is not a JSON object")
    requested = content.get("replicates_requested")
    completed = content.get("replicates_completed")
    recorded = _composition_replicates(replicates, content)
    if not isinstance(requested, int) or not isinstance(completed, int) or completed != recorded:
        raise ValueError(
            f"{summary} reports {requested!r} requested and {completed!r} completed replicates; "
            f"{replicates.name} records {recorded} complete"
        )
    verdict = build_calibration_verdict(
        verdict_id="verdict:jed_log_ratio_vs_transformed_t_calibration",
        production_spec=production_spec,
        evaluated_spec=evaluated_spec,
        replicate_budget_sufficient=recorded >= max(requested, DEFAULT_REPLICATES),
        null_control_passed=None,
        power_passed=None,
    )
    verdict.summary_message += (
        f" Read from {summary.name} and {replicates.name}: {recorded} of {requested} "
        f"replicates completed against a default budget of {DEFAULT_REPLICATES}; it reports "
        "null error rates without a tolerance and runs no injected-effect scenario."
    )
    return verdict


def shared_basis_ledger(
    spec: AnalysisSpec, core_run: WorkflowRun, contrast_population: ContrastPopulation
) -> EligibilityLedger:
    """The studies integrated evidence's shared basis keeps for ``spec`` (#815, #844).

    It loads the core's count profiles with integrated evidence's guilds and filters and
    applies the same whole-pair zero-total exclusion, so the ledger exists before the
    stage writes anything. ``reconcile_ledger`` then checks what integrated evidence wrote.
    """
    if spec.sampling_unit != "study":
        raise ValueError(f"No shared-basis ledger for {spec.sampling_unit!r} units: {spec.id}")
    profiles, _ = load_profiles(
        artifact_path(core_run, "count_profiles"),
        artifact_path(core_run, "experiments"),
        {**TABLE1_TRAITS, **OXYGEN_GUILDS},
        reject_missing_counts=True,
    )
    exclusion = drop_zero_total_pairs(profiles)
    return build_eligibility_ledger(
        ledger_id=f"ledger:{spec.id}",
        analysis_spec_id=spec.id,
        contrast_population=contrast_population,
        included_ids={profile.study_id for profile in exclusion.profiles},
        exclusion_counts={"zero_retained_total": exclusion.studies_lost},
    )


def reconcile_ledger(ledger: EligibilityLedger, provenance: Mapping[str, object]) -> None:
    """Require integrated evidence's recorded exclusion to be the ledger's (#844).

    A record without the exclusion block is an error, never zero exclusions.
    """
    exclusion = provenance.get("pair_exclusion")
    if not isinstance(exclusion, Mapping):
        raise ValueError("integrated_evidence provenance has no pair_exclusion block")
    before = exclusion.get("studies_before")
    after = exclusion.get("studies_after")
    lost = exclusion.get("studies_lost")
    if not (isinstance(before, int) and isinstance(after, int) and isinstance(lost, int)):
        raise ValueError("pair_exclusion lacks its study counts")
    if (
        before != after + lost
        or exclusion.get("study_ids_after") != ledger.included_ids
        or (after, lost) != (ledger.included_count, ledger.excluded_count)
    ):
        raise ValueError(
            f"pair_exclusion (studies {before} -> {after}, {lost} lost) does not reconcile "
            f"with {ledger.id} ({ledger.included_count} included, "
            f"{ledger.excluded_count} excluded)"
        )


def retained_profile_features(catalog: Path, profiles: Path, *, ora: bool) -> tuple[str, ...]:
    """A profile table's retained (non-MediaDive) trait features, in its column order.

    They must be the feature catalog's trait_annotation features in the catalog's order,
    so a basis identifies the family its profiles actually hold (#860).
    """
    columns = profile_count_reader(profiles).count_columns
    if ora:
        if columns[0] != "taxa_count":
            raise ValueError(f"{profiles} lacks the ORA taxa_count denominator column")
        columns = columns[1:]
    retained = tuple(column for column in columns if feature_class(column) == "trait_annotation")
    with catalog.open(encoding="utf-8", newline="") as handle:
        catalogued = tuple(
            row["column_id"]
            for row in csv.DictReader(handle, delimiter="\t")
            if row["feature_class"] == "trait_annotation"
        )
    if retained != catalogued:
        raise ValueError(
            f"The retained features of {profiles} are not the feature catalog's "
            "trait_annotation features in catalog order"
        )
    return retained


def core_analysis_bases(core_run: WorkflowRun) -> list[AnalysisBasis]:
    """The count and ORA analysis bases of a verified core run.

    Each basis records the core's feature catalog by its bytes, the identity and number of
    the retained (non-MediaDive) trait features its profile table holds, in column order
    (#860), and the strain rollup scope the core declared, so a comparison built on the
    phenotype primary or the off sensitivity core says which it is (#788).
    """
    catalog = artifact_path(core_run, "feature_catalog")
    catalog_sha = file_sha256(catalog)
    count = retained_profile_features(catalog, artifact_path(core_run, "count_profiles"), ora=False)
    ora = retained_profile_features(catalog, artifact_path(core_run, "ora_profiles"), ora=True)
    rollup_policy = f"ncbi_taxonomy_closure_rollup_{run_strain_rollup(core_run)}"
    return [
        build_analysis_basis(
            basis_id="basis:count_profiles_non_mediadive",
            feature_catalog_sha256=catalog_sha,
            feature_filter_id=feature_filter_id(count),
            retained_feature_count=len(count),
            profile_unit="count_profiles",
            aggregation_policy="sum_taxon_trait_assertions_per_signature",
            rollup_policy=rollup_policy,
            denominator_policy="sum_retained_non_mediadive_trait_assertions",
            background_inventory_id="human_host_retained_profiles",
            weighting_policy="equal_experiment_or_study_weight",
        ),
        build_analysis_basis(
            basis_id="basis:ora_profiles_binary",
            feature_catalog_sha256=catalog_sha,
            feature_filter_id=feature_filter_id(ora),
            retained_feature_count=len(ora),
            profile_unit="ora_profiles",
            aggregation_policy="distinct_matched_taxa_per_signature",
            rollup_policy=rollup_policy,
            denominator_policy="distinct_matched_signature_taxa_N_gt_0",
            background_inventory_id="taxa_trait_matrix_universe",
            weighting_policy="unweighted_hypergeometric_per_signature",
        ),
    ]


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
    git_commit_actual, git_dirty = repository_state(REPO_ROOT)
    effective_commit = code_commit or git_commit_actual
    contrast_population = core_run.contrast_population or ContrastPopulation.all_paired
    direction_basis = core_run.direction_basis or DirectionBasis.group_1_vs_group_0

    # Every arm input comes from what its producing step recorded, and every input the
    # verdicts, the ledger and the analysis bases read is parsed before the stage writes
    # anything, so a bad input leaves no partial directory to block a rerun (#814, #844).
    arm_runs = {"noel": noel_run, "jed": jed_run, "ilyas": ilyas_run}
    evidence_inputs = [
        part
        for option, arm, step, name in ARM_EVIDENCE_INPUTS
        for part in (option, str(recorded_output(arm_runs[arm], step, name)))
    ]
    calibration = {
        name: recorded_output(arm_runs[arm], step, name) for arm, step, name in CALIBRATION_SOURCES
    }
    table1_spec = STANDARD_ANALYSIS_SPECS["noel_table1_six_trait_sign_test_v1"]
    verdict_noel = noel_calibration_verdict(calibration["calibration_status.json"], table1_spec)
    verdict_jed_mismatch = jed_calibration_verdict(
        calibration["composition_calibration_summary.json"],
        calibration["composition_calibration_replicates.tsv"],
        STANDARD_ANALYSIS_SPECS["jed_table1_paired_log_ratio_sign_test_v1"],
        STANDARD_ANALYSIS_SPECS["jed_composition_arcsine_sqrt_ttest_calibration_v1"],
    )
    ledger_fraction = shared_basis_ledger(table1_spec, core_run, contrast_population)
    analysis_bases = core_analysis_bases(core_run)

    output.mkdir(parents=True, exist_ok=False)
    if arm1_pooled_stats is not None:
        pooled_stats_path = arm1_pooled_stats.resolve()
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
        asyncio.run(_exec_command(cmd, REPO_ROOT))
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
        *evidence_inputs,
        "--code-commit",
        effective_commit,
        "--require-complete-sources",
    ]
    integrated_evidence_main(ie_argv)

    provenance = json.loads(
        (integrated_dir / "integrated_evidence_provenance.json").read_text(encoding="utf-8")
    )
    reconcile_ledger(ledger_fraction, provenance)

    inventory_entries = _inventory_directory(output, integrated_dir)
    if (output / "arm1_pooled_stats").is_dir():
        inventory_entries.extend(_inventory_directory(output, output / "arm1_pooled_stats"))

    ended_at = datetime.now(UTC)
    manifest_file = output / "comparison_run.json"
    comparison_record = ComparisonRun(
        id=f"comparison-{contrast_population}-{uuid.uuid4().hex[:12]}",
        contrast_population=contrast_population,
        direction_basis=direction_basis,
        status=ExecutionStatus.completed,
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
        analysis_bases=analysis_bases,
        calibration_verdicts=[verdict_noel, verdict_jed_mismatch],
        outputs=inventory_entries,
    )
    manifest_file.write_text(
        json.dumps(comparison_record.model_dump(mode="json"), indent=2) + "\n",
        encoding="utf-8",
    )
    return comparison_record


def _require_derived_bases(
    recorded: Sequence[AnalysisBasis], derived: Sequence[AnalysisBasis]
) -> None:
    """Refuse recorded bases other than the derived ones, each once and in order (#941)."""
    recorded_ids = [basis.id for basis in recorded]
    derived_ids = [basis.id for basis in derived]
    if recorded_ids != derived_ids:
        raise ValueError(
            f"ComparisonRun records analysis bases {recorded_ids}, but its core derives "
            f"{derived_ids}"
        )
    differences = [
        f"{basis.id} {field}"
        for basis, expected in zip(recorded, derived, strict=True)
        for field, value in basis.model_dump().items()
        if value != expected.model_dump()[field]
    ]
    if differences:
        raise ValueError(
            "ComparisonRun analysis bases differ from those its core derives: "
            + ", ".join(differences)
        )


def verify_comparison_run(path: Path) -> ComparisonRun:
    """Verify a completed ComparisonRun manifest, upstream runs, its analysis bases and
    complete outputs (#776, #817)."""
    manifest = (path / "comparison_run.json" if path.is_dir() else path).resolve()
    if not manifest.is_file():
        raise FileNotFoundError(f"ComparisonRun manifest not found: {manifest}")
    record = ComparisonRun.model_validate_json(manifest.read_text(encoding="utf-8"))
    if record.status != ExecutionStatus.completed:
        raise ValueError(
            f"ComparisonRun {record.id} has status {record.status}, expected completed"
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

    # The bases are derived again from the verified core: its recorded scope, its catalog
    # and the retained features of its profile tables (#817, #860)
    _require_derived_bases(record.analysis_bases, core_analysis_bases(core_run))

    root = manifest.parent
    if not record.outputs:
        raise ValueError("ComparisonRun has no recorded outputs")
    for entry in record.outputs:
        target = root / entry.path
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
            f"{row.get('method')}:{row.get('trait') or row.get('feature')}"
            for row in reader
            if row.get("test_status") == "source_not_supplied"
        ]
    if missing:
        raise ValueError(
            "ComparisonRun contains source_not_supplied placeholder rows: " + ", ".join(missing)
        )
    return record
