"""The integrated evidence table: states, families and what a join must not lose (#178)."""

from __future__ import annotations

import csv
from dataclasses import replace
from pathlib import Path

import pytest

from bsdbng.integrated_evidence import (
    ALL_TIED,
    NOT_APPLICABLE,
    NOT_SIGNIFICANT,
    SIGNIFICANT,
    SOURCE_MISSING,
    TESTED,
    UNDECIDED,
    ZERO_ANNOTATION,
    Record,
    apply_statuses,
    fraction_sign_records,
    log_ratio_records,
    placeholder_records,
    read_status,
    summarize,
    write_table,
    write_view,
)
from bsdbng.table1_trait_maps import compute_direction_shifts
from bsdbng.trait_guilds import TABLE1_TRAITS, Profile

TRAITS = tuple(TABLE1_TRAITS)
STATUSES = {
    "calibration_status": "fraction_sign_test: pass",
    "sensitivity_status": "8 scenarios reported",
    "coverage_status": "pending #184",
    "orientation_status": "1052 classified",
    "taxon_influence_status": "pending #190",
}


def _profiles(pattern: str = "signal") -> list[Profile]:
    """Twelve studies whose counts give a signal, a tie, or no annotation at all."""
    profiles: list[Profile] = []
    for study in range(12):
        for direction in ("inc", "dec"):
            counts = dict.fromkeys(TRAITS, 5)
            if pattern == "signal" and direction == "inc":
                counts[TRAITS[0]] = 30
            if pattern == "tied":
                counts = dict.fromkeys(TRAITS, 5)
            if pattern == "empty":
                counts = dict.fromkeys(TRAITS, 0)
            profiles.append(
                Profile(
                    study_id=f"s{study}",
                    experiment_id=f"s{study}-1",
                    inc_dec=direction,
                    condition="C",
                    body_site="Feces",
                    guild_counts=counts,
                    total=sum(counts.values()) + (0 if pattern == "empty" else 20),
                )
            )
    return profiles


def test_counts_reconcile_with_the_source_run() -> None:
    profiles = _profiles()
    records = {r.trait: r for r in fraction_sign_records(profiles, source="x", q_threshold=0.05)}
    for shift in compute_direction_shifts(profiles).shifts:
        record = records[shift.trait]
        assert record.testable_units == shift.studies_scored
        assert record.ties == shift.studies_tied
        assert record.effect_value == shift.share_up
        assert record.q_value == shift.q_value
        assert record.unique_studies == 12
        assert record.family_size == 6


def test_an_all_tied_trait_is_not_a_negative_result() -> None:
    records = fraction_sign_records(_profiles("tied"), source="x", q_threshold=0.05)
    for record in records:
        assert record.test_status == ALL_TIED
        assert record.decision == UNDECIDED
        assert record.effect_value is None
        assert record.q_value is None
        assert record.interval_type == "unavailable"


def test_a_pair_with_no_positive_counts_is_recorded_as_zero_annotation() -> None:
    records = log_ratio_records(_profiles("empty"), source="x", q_threshold=0.05)
    assert records
    for record in records:
        assert record.test_status == ZERO_ANNOTATION
        assert record.decision == UNDECIDED
        assert record.q_value is None


def test_an_unsupplied_method_is_not_a_nonsignificant_result() -> None:
    records = placeholder_records("enrichment", "profile", "not supplied", TRAITS)
    assert len(records) == len(TRAITS)
    for record in records:
        assert record.test_status == SOURCE_MISSING
        assert record.decision == UNDECIDED
        assert record.q_value is None
        assert record.notes == "not supplied"


def test_each_method_keeps_its_own_family_and_units() -> None:
    profiles = _profiles()
    records = fraction_sign_records(profiles, source="x", q_threshold=0.05) + log_ratio_records(
        profiles, source="x", q_threshold=0.05
    )
    families = {(r.method, r.family_id, r.correction_method) for r in records}
    assert (
        "directional_consistency_fraction_sign_test",
        "table1_six_fractions",
        "Benjamini-Hochberg",
    ) in families
    assert (
        "directional_consistency_log_ratio_sign_test",
        "table1_fifteen_ratios",
        "Benjamini-Yekutieli",
    ) in families
    assert len({r.effect_units for r in records}) == 2
    summary = {row["method"]: row for row in summarize(records)}
    assert summary["directional_consistency_fraction_sign_test"]["families"] == (
        "table1_six_fractions"
    )
    assert all(
        row["combined_guarantee"] == "none; each family is corrected on its own"
        for row in summary.values()
    )


def test_display_state_is_independent_of_the_decision() -> None:
    record = fraction_sign_records(_profiles(), source="x", q_threshold=0.05)[0]
    masked = replace(record, display_status="masked_below_min_studies")
    assert masked.decision == record.decision
    assert record.display_status == NOT_APPLICABLE
    assert masked.test_status == TESTED


def test_statuses_are_attached_without_changing_any_result() -> None:
    records = fraction_sign_records(_profiles(), source="x", q_threshold=0.05)
    annotated = apply_statuses(records, STATUSES)
    for before, after in zip(records, annotated, strict=True):
        assert (after.decision, after.test_status, after.q_value) == (
            before.decision,
            before.test_status,
            before.q_value,
        )
        assert after.calibration_status == "fraction_sign_test: pass"
        assert after.coverage_status == "pending #184"
    with pytest.raises(ValueError, match="no status supplied"):
        apply_statuses(records, {"calibration_status": "x"})


def test_absent_assessment_files_read_as_pending(tmp_path: Path) -> None:
    assert read_status(None, "calibration_status", "#182") == "pending #182"
    assert read_status(tmp_path / "missing.json", "coverage_status", "#184") == "pending #184"
    broken = tmp_path / "broken.json"
    broken.write_text("not json", encoding="utf-8")
    assert read_status(broken, "coverage_status", "#184") == "unreadable #184"
    written = tmp_path / "calibration.json"
    written.write_text('{"calibration_status": {"a": "pass"}}', encoding="utf-8")
    assert read_status(written, "calibration_status", "#182") == "a: pass"


def test_every_field_survives_the_export(tmp_path: Path) -> None:
    profiles = _profiles()
    records = apply_statuses(
        fraction_sign_records(profiles, source="x", q_threshold=0.05)
        + placeholder_records("enrichment", "profile", "not supplied", TRAITS[:1]),
        STATUSES,
    )
    path = write_table(records, tmp_path / "evidence.tsv")
    with path.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == len(records)
    assert {"test_status", "decision", "display_status", "family_id", "interval_type"} <= set(
        rows[0]
    )
    unsupplied = next(r for r in rows if r["method"] == "enrichment")
    assert unsupplied["q_value"] == ""
    assert unsupplied["test_status"] == SOURCE_MISSING
    assert unsupplied["decision"] == UNDECIDED
    signal = next(r for r in rows if r["trait"] == TRAITS[0] and r["method"] != "enrichment")
    assert signal["decision"] in (SIGNIFICANT, NOT_SIGNIFICANT)
    view = write_view(records, summarize(records), tmp_path / "evidence.md")
    text = view.read_text(encoding="utf-8")
    assert "no pooled p-value" in text
    assert "source_not_supplied" in text


def test_a_record_carries_the_selection_note() -> None:
    record = fraction_sign_records(_profiles(), source="x", q_threshold=0.05)[0]
    assert "selected on this corpus" in record.selection_note
    assert isinstance(record, Record)
