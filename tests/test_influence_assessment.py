"""Native execution is not a procedure verdict; preserve the handoff (#671)."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest

from mtp.influence_assessment import ASSESSMENT_VERSION, NATIVE_RESULT_COUNTS
from mtp.integrated_evidence import (
    ASSESSMENTS,
    AssessmentStatus,
    apply_statuses,
    assessment_rows,
    placeholder_records,
    read_status,
    write_view,
)
from mtp.paths import file_sha256
from mtp.taxon_contributions import write_table


def _write(path: Path, content: object) -> Path:
    path.write_text(json.dumps(content, indent=2) + "\n", encoding="utf-8")
    return path


def _native(tmp_path: Path, *, completed: bool = True) -> Path:
    family = tmp_path / "native_ora_families.tsv"
    family.write_text(
        "scenario\tscope\tbackground\tprofiles\n"
        "baseline\tpooled\tfixed_baseline\t10\n"
        "baseline\tpooled\trecomputed_remaining\t10\n"
        "exclude_genus_NCBITaxon:561\tpooled\tfixed_baseline\t8\n"
        "exclude_genus_NCBITaxon:561\tpooled\trecomputed_remaining\t8\n",
        encoding="utf-8",
    )
    scenarios = tmp_path / "scenarios.tsv"
    scenarios.write_text("scenario\tstatus\nexclude_genus_NCBITaxon:561\trun\n", encoding="utf-8")
    content: dict[str, object] = {
        "completed": completed,
        "protocol": {
            "version": "190-native-v2",
            "composition": {"weighting": "equal experiments", "minimum_pairs": 3},
            "ora": {"interpretation": "count-model diagnostic; not full-pipeline calibration"},
        },
        "protocol_sha256": "a" * 64,
        "contribution_manifest_sha256": "b" * 64,
        "input_sha256": {
            "study_trait_profiles_wide.tsv": "c" * 64,
            "experiments.tsv": "d" * 64,
        },
        "code_commit": "e" * 40,
        "code_dirty": False,
        "code_sha256": {"native_taxon_influence.py": "f" * 64},
        "runtime": {"python": "3.11"},
    }
    if completed:
        outputs = {family.name: file_sha256(family)}
        inputs = {scenarios.name: file_sha256(scenarios)}
        for name, count in NATIVE_RESULT_COUNTS.items():
            if count == "ora_families":
                continue
            path = tmp_path / name
            write_table(
                path,
                ["trait", "status"],
                [{"trait": f"T{i}", "status": "tested"} for i in range(12)],
            )
            outputs[name] = file_sha256(path)
        summary = tmp_path / "native_summary.md"
        summary.write_text("Native results\n", encoding="utf-8")
        outputs[summary.name] = file_sha256(summary)
        for scenario in ("baseline", "exclude_genus_NCBITaxon:561"):
            directory = tmp_path / scenario.replace(":", "_")
            directory.mkdir()
            for name, group in (
                ("removed_source_counts.tsv.gz", inputs),
                ("profile_counts.tsv.gz", inputs),
                ("native_ora_profiles.tsv.gz", outputs),
            ):
                path = directory / name
                write_table(path, ["trait", "status"], [{"trait": "T", "status": "tested"}])
                group[path.relative_to(tmp_path).as_posix()] = file_sha256(path)
        content["results"] = {
            "scenarios": 2,
            "scopes": 1,
            "composition_rows": 12,
            "composition_comparisons": 12,
            "ora_comparisons": 12,
            "ora_families": 4,
            "input_sha256": inputs,
            "output_sha256": outputs,
        }
    return _write(tmp_path / "native_manifest.json", content)


def _review(native_path: Path) -> Path:
    root = native_path.parent / "review"
    root.mkdir()
    (root / "run_manifest.json").write_bytes(native_path.read_bytes())
    names = (
        "native_table1_stats.tsv.gz",
        "native_composition_comparison.tsv.gz",
        "native_ora_fold_comparison.tsv.gz",
        "native_ora_families.tsv",
        "scenarios.tsv",
        "interpretation_handoff.tsv",
    )
    tables: dict[str, object] = {}
    for name in names:
        path = root / name
        row_count = 1
        if (native_path.parent / name).exists():
            path.write_bytes((native_path.parent / name).read_bytes())
            row_count = len(path.read_text().splitlines()) - 1
        elif name == "interpretation_handoff.tsv":
            path.write_text(
                "finding_id\tscope\tscenario\tflag_status\tintegration_rule\n"
                "lineage_example\tpooled\texclude_genus_NCBITaxon:561\tlineage_sensitive\t"
                "interpretation flag only; preserve separate method rows\n",
                encoding="utf-8",
            )
        else:
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write("trait\tstatus\nexample\ttested\n")
                if "comparison" in name:
                    handle.write("example\ttested\n" * 11)
                    row_count = 12
        tables[name] = {"rows": row_count, "sha256": file_sha256(path)}
    return _write(
        root / "publication_manifest.json",
        {
            "native_manifest_sha256": file_sha256(native_path),
            "publisher_sha256": "a" * 64,
            "source_reconciliation": {
                "mismatches": 0,
                "source_contribution_rows": 5,
                "profile_reported_taxon_rows": 5,
                "profile_count_rows": 5,
                "sha256": {"source_contributions.tsv.gz": "b" * 64},
            },
            "interpretation_source_sha256": {"effects.tsv.gz": "c" * 64},
            "tables": tables,
        },
    )


def _read(path: Path | None) -> AssessmentStatus:
    return read_status(path, "taxon_influence_status", "#190")


@pytest.mark.parametrize("supplied", [False, True])
def test_absent_evidence_does_not_assert_analysis_never_ran(tmp_path: Path, supplied: bool) -> None:
    status = _read(tmp_path / "missing.json" if supplied else None)
    assert status.default == "pending #190: influence evidence not supplied"
    assert not status.verdicts


def test_completed_native_diagnostic_has_no_procedure_verdict(tmp_path: Path) -> None:
    path = _native(tmp_path)
    status = _read(path)
    assert status.default == (
        "native influence results reported; no procedure-level verdict supplied #190"
    )
    assert status.verdicts == {}
    assert ("source", "native_manifest", "sha256", file_sha256(path)) in status.detail
    assert ("execution", "native_influence", "completed", "true") in status.detail
    assert any(
        kind == "native_ora_family" and "NCBITaxon:561" in scope
        for kind, scope, _, _ in status.detail
    )
    assert any(
        kind == "protocol" and "equal experiments" in value for kind, _, _, value in status.detail
    )
    records = placeholder_records("fraction_sign_test", "study", "not supplied", ["trait"])
    statuses = {key: AssessmentStatus(issue, {}, "pending") for key, issue in ASSESSMENTS.items()}
    statuses["taxon_influence_status"] = status
    applied = apply_statuses(records, statuses)
    assert applied[0].taxon_influence_status == status.default
    exported = assessment_rows(applied, statuses)
    assert any(row["scope_kind"] == "native_ora_family" for row in exported)


def test_incomplete_native_run_is_distinct_from_no_input(tmp_path: Path) -> None:
    status = _read(_native(tmp_path, completed=False))
    assert status.default == "native influence execution incomplete #190"
    assert not status.verdicts
    assert ("execution", "native_influence", "completed", "false") in status.detail


@pytest.mark.parametrize("publication_input", [False, True])
def test_review_handoff_retains_scoped_interpretation_without_global_verdict(
    tmp_path: Path,
    publication_input: bool,
) -> None:
    native = _native(tmp_path)
    review = _review(native)
    status = _read(review if publication_input else native)
    assert status.default.startswith("native influence results reported;")
    assert not status.verdicts
    assert ("source", "publication", "sha256", file_sha256(review)) in status.detail
    assert (
        "interpretation_flag",
        "lineage_example",
        "flag_status",
        "lineage_sensitive",
    ) in status.detail
    assert ("interpretation_flag", "lineage_example", "scope", "pooled") in status.detail


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("completed", "true", "completed must be a boolean"),
        ("protocol", {"version": "unknown"}, "incompatible"),
        ("protocol", [], "protocol must be an object"),
        ("input_sha256", {}, "must not be empty"),
        ("code_commit", "unknown", "full Git revision"),
        ("code_dirty", "false", "code_dirty must be a boolean"),
        ("results", {}, "results.scenarios"),
    ],
)
def test_incompatible_or_malformed_native_evidence_fails_closed(
    tmp_path: Path,
    field: str,
    value: object,
    message: str,
) -> None:
    path = _native(tmp_path)
    content = json.loads(path.read_text())
    content[field] = value
    _write(path, content)
    status = _read(path)
    assert status.default.startswith("unreadable #190:")
    assert message in status.default
    assert not status.verdicts


@pytest.mark.parametrize("content", [[], {"completed": True}, {"taxon_influence_status": "pass"}])
def test_unknown_and_global_verdict_contracts_are_not_copied_to_every_row(
    tmp_path: Path,
    content: object,
) -> None:
    status = _read(_write(tmp_path / "bad.json", content))
    assert status.default.startswith("unreadable #190:")
    assert not status.verdicts


def test_changed_native_artifact_and_incomplete_inventory_are_rejected(tmp_path: Path) -> None:
    path = _native(tmp_path)
    family = tmp_path / "native_ora_families.tsv"
    family.write_text(family.read_text().replace("recomputed_remaining", "fixed_baseline"))
    assert "changed evidence artifact" in _read(path).default
    content = json.loads(path.read_text())
    content["results"]["output_sha256"][family.name] = file_sha256(family)
    _write(path, content)
    assert "duplicate family identities" in _read(path).default


def test_malformed_table_is_unreadable_without_discarding_supplied_source_hash(
    tmp_path: Path,
) -> None:
    path = _native(tmp_path)
    family = tmp_path / "native_ora_families.tsv"
    family.write_text(family.read_text() + "baseline\tpooled\tfixed_baseline\t1\textra\n")
    content = json.loads(path.read_text())
    content["results"]["output_sha256"][family.name] = file_sha256(family)
    _write(path, content)
    status = _read(path)
    assert "malformed table row" in status.default
    assert ("source", "supplied_influence", "sha256", file_sha256(path)) in status.detail


def test_artifact_path_cannot_escape_the_native_directory(tmp_path: Path) -> None:
    path = _native(tmp_path)
    content = json.loads(path.read_text())
    content["results"]["output_sha256"]["../escape.tsv"] = "a" * 64
    _write(path, content)
    assert "artifact escapes" in _read(path).default


@pytest.mark.parametrize(
    ("group", "name"),
    [
        ("input_sha256", "baseline/profile_counts.tsv.gz"),
        ("input_sha256", "exclude_genus_NCBITaxon_561/removed_source_counts.tsv.gz"),
        ("output_sha256", "native_composition.tsv.gz"),
        ("output_sha256", "native_composition_comparison.tsv"),
        ("output_sha256", "native_ora_fold_comparison.tsv"),
        ("output_sha256", "baseline/native_ora_profiles.tsv.gz"),
    ],
)
@pytest.mark.parametrize("publication_input", [False, True])
def test_missing_declared_native_artifact_is_not_completed_evidence(
    tmp_path: Path, group: str, name: str, publication_input: bool
) -> None:
    native = _native(tmp_path)
    content = json.loads(native.read_text())
    del content["results"][group][name]
    _write(native, content)
    (tmp_path / name).unlink()
    path = _review(native) if publication_input else native
    status = _read(path)
    assert "lacks required artifacts" in status.default
    assert name in status.default
    assert not status.verdicts


@pytest.mark.parametrize("name", list(NATIVE_RESULT_COUNTS)[:-1])
def test_native_result_counts_must_match_actual_tables(tmp_path: Path, name: str) -> None:
    native = _native(tmp_path)
    content = json.loads(native.read_text())
    content["results"][NATIVE_RESULT_COUNTS[name]] += 1
    _write(native, content)
    assert f"native table row count disagrees with manifest: {name}" in _read(native).default


@pytest.mark.parametrize(
    "name", ["native_composition_comparison.tsv.gz", "native_ora_fold_comparison.tsv.gz"]
)
def test_self_consistent_truncated_compact_table_disagrees_with_native_results(
    tmp_path: Path, name: str
) -> None:
    native = _native(tmp_path)
    review = _review(native)
    table = review.parent / name
    write_table(table, ["trait", "status"], [{"trait": "T", "status": "tested"}])
    content = json.loads(review.read_text())
    content["tables"][name] = {"rows": 1, "sha256": file_sha256(table)}
    _write(review, content)
    assert "publication table row count disagrees with native results" in _read(review).default


def test_family_scenarios_must_match_selected_scenarios(tmp_path: Path) -> None:
    native = _native(tmp_path)
    scenarios = tmp_path / "scenarios.tsv"
    scenarios.write_text("scenario\tstatus\nexclude_other\trun\n", encoding="utf-8")
    content = json.loads(native.read_text())
    content["results"]["input_sha256"][scenarios.name] = file_sha256(scenarios)
    _write(native, content)
    assert "disagree with the selected scenario inventory" in _read(native).default


@pytest.mark.parametrize("mutation", ["native_hash", "reconciliation", "table_hash", "row_count"])
def test_review_mismatch_is_never_an_execution_success(tmp_path: Path, mutation: str) -> None:
    native = _native(tmp_path)
    path = _review(native)
    content = json.loads(path.read_text())
    if mutation == "native_hash":
        content["native_manifest_sha256"] = "f" * 64
    elif mutation == "reconciliation":
        content["source_reconciliation"]["mismatches"] = 1
    elif mutation == "table_hash":
        content["tables"]["native_table1_stats.tsv.gz"]["sha256"] = "f" * 64
    else:
        content["tables"]["native_table1_stats.tsv.gz"]["rows"] = 999
    _write(path, content)
    assert _read(native).default.startswith("unreadable #190:")
    assert not _read(path).verdicts


def test_different_input_corpus_cannot_receive_native_execution_standing(tmp_path: Path) -> None:
    status = read_status(
        _native(tmp_path),
        "taxon_influence_status",
        "#190",
        expected_input_sha256={"experiments.tsv": "a" * 64},
    )
    assert "different input bytes: experiments.tsv" in status.default


def _explicit(path: Path) -> Path:
    return _write(
        path,
        {
            "schema_version": ASSESSMENT_VERSION,
            "input_sha256": {
                "study_trait_profiles_wide.tsv": "c" * 64,
                "experiments.tsv": "d" * 64,
            },
            "procedure_verdicts": {
                "fraction_sign_test": {
                    "verdict": "lineage_sensitive",
                    "applies_to": "all_rows_of_procedure",
                    "basis": "Declared procedure-wide assessment of the diagnostic evidence",
                },
                "log_ratio_sign_test": {
                    "verdict": "not_assessed",
                    "applies_to": "all_rows_of_procedure",
                    "basis": "Assessment explicitly excludes this procedure",
                },
            },
        },
    )


def test_explicit_procedure_verdicts_never_leak_to_other_procedures(tmp_path: Path) -> None:
    status = _read(_explicit(tmp_path / "assessment.json"))
    assert status.for_procedure("fraction_sign_test") == "lineage_sensitive #190"
    assert status.for_procedure("log_ratio_sign_test") == "not_assessed #190"
    assert (
        status.for_procedure("guild_mcnemar")
        == "no procedure-level influence verdict supplied #190"
    )
    assert any(
        kind == "procedure" and scope == "fraction_sign_test" and key == "basis"
        for kind, scope, key, _ in status.detail
    )


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("verdict", {"nested": "pass"}),
        ("basis", ""),
        ("applies_to", "pooled_only"),
    ],
)
def test_partial_or_malformed_procedure_assessment_is_not_promoted(
    tmp_path: Path,
    field: str,
    value: object,
) -> None:
    path = _explicit(tmp_path / "assessment.json")
    content = json.loads(path.read_text())
    content["procedure_verdicts"]["fraction_sign_test"][field] = value
    _write(path, content)
    status = _read(path)
    assert status.default.startswith("unreadable #190:")
    assert not status.verdicts


def test_generated_dictionary_does_not_equate_pending_with_not_run(tmp_path: Path) -> None:
    path = write_view([], [], tmp_path / "view.md")
    text = path.read_text()
    assert "assessment has not run" not in text
    assert "does not establish that execution has not occurred" in text
    assert "never a scientific pass" in text
