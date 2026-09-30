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
from bsdbng.trait_guilds import TABLE1_TRAITS, Profile, drop_zero_total_pairs

#: Both decision thresholds every record builder now takes (#294).
THRESHOLDS = {"q_threshold": 0.05, "strict_q_threshold": 0.001}

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
    built = fraction_sign_records(profiles, source="x", **THRESHOLDS)
    records = {r.trait: r for r in built}
    for shift in compute_direction_shifts(profiles).shifts:
        record = records[shift.trait]
        assert record.testable_units == shift.studies_scored
        assert record.ties == shift.studies_tied
        assert record.effect_value == shift.share_up
        assert record.q_value == shift.q_value
        assert record.unique_studies == 12
        assert record.family_size == 6


def test_an_all_tied_trait_is_not_a_negative_result() -> None:
    records = fraction_sign_records(_profiles("tied"), source="x", **THRESHOLDS)
    for record in records:
        assert record.test_status == ALL_TIED
        assert record.decision == UNDECIDED
        assert record.effect_value is None
        assert record.q_value is None
        assert record.interval_type == "unavailable"


def test_a_pair_with_no_positive_counts_is_recorded_as_zero_annotation() -> None:
    records = log_ratio_records(_profiles("empty"), source="x", **THRESHOLDS)
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
    records = fraction_sign_records(profiles, source="x", **THRESHOLDS) + log_ratio_records(
        profiles, source="x", **THRESHOLDS
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
    record = fraction_sign_records(_profiles(), source="x", **THRESHOLDS)[0]
    masked = replace(record, display_status="masked_below_min_studies")
    assert masked.decision == record.decision
    assert record.display_status == NOT_APPLICABLE
    assert masked.test_status == TESTED


def test_statuses_are_attached_without_changing_any_result() -> None:
    records = fraction_sign_records(_profiles(), source="x", **THRESHOLDS)
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
        fraction_sign_records(profiles, source="x", **THRESHOLDS)
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
    record = fraction_sign_records(_profiles(), source="x", **THRESHOLDS)[0]
    assert "selected on this corpus" in record.selection_note
    assert isinstance(record, Record)


def test_a_zero_total_pair_is_excluded_rather_than_counted_as_a_tie() -> None:
    """A pair with no retained annotation on a side has an undefined share (#294).

    ``Profile.fraction`` returns 0.0 for a zero total, so such a pair enters the
    contrast as a measured share of nothing, the contrast comes out exactly zero,
    and the sign test counts it as a tie. That merges all-tied studies with
    zero-annotation studies, which the shared comparison basis of #178 forbids.
    """
    scored = _profiles()
    blank = [
        Profile(
            study_id="bsdb:blank",
            experiment_id="bsdb:blank-1",
            inc_dec=direction,
            condition="c",
            body_site="s",
            guild_counts=dict.fromkeys(TRAITS, 0),
            total=0,
        )
        for direction in ("inc", "dec")
    ]

    exclusion = drop_zero_total_pairs(scored + blank)

    assert exclusion.zero_total_experiments == 1
    assert exclusion.studies_lost == 1
    assert exclusion.studies_after == len({p.study_id for p in scored})
    assert all(p.study_id != "bsdb:blank" for p in exclusion.profiles)

    # Kept as a tie, the blank study would inflate the tie count; excluded, the
    # tallies match the corpus that has a defined share.
    with_blank = fraction_sign_records(scored + blank, source="x", **THRESHOLDS)
    without = fraction_sign_records(exclusion.profiles, source="x", **THRESHOLDS)
    assert {r.trait: r.ties for r in with_blank} != {r.trait: r.ties for r in without}


def test_each_row_is_decided_at_both_thresholds_with_their_sources() -> None:
    """One cutoff cannot speak for three families, so both are reported (#294)."""
    records = fraction_sign_records(_profiles(), source="x", **THRESHOLDS)

    for record in records:
        assert record.q_threshold == 0.05
        assert record.strict_q_threshold == 0.001
        # Each threshold says where it comes from rather than appearing as a bare number.
        assert record.q_threshold_source
        assert record.strict_q_threshold_source
        assert record.q_threshold_source != record.strict_q_threshold_source


def test_the_two_thresholds_decide_independently() -> None:
    """The strict column must read its own threshold, not copy the loose decision.

    The fixture's q values sit below both default thresholds, so a run at the
    defaults cannot tell the two columns apart - which is how a strict column
    that silently reuses ``q_threshold`` would pass. The thresholds here bracket
    the fixture's q instead, so the two decisions are required to disagree (#294).
    """
    records = fraction_sign_records(
        _profiles(), source="x", q_threshold=0.05, strict_q_threshold=1e-9
    )

    assert records
    for record in records:
        assert record.q_value is not None
        assert 1e-9 < record.q_value <= 0.05, "fixture must sit between the two thresholds"
        assert record.decision == SIGNIFICANT
        assert record.strict_decision == NOT_SIGNIFICANT


def test_the_stricter_threshold_can_only_shrink_the_significant_set() -> None:
    """A row significant at the stricter cutoff is significant at the looser one."""
    records = fraction_sign_records(_profiles(), source="x", **THRESHOLDS)

    assert records
    for record in records:
        if record.strict_decision == SIGNIFICANT:
            assert record.decision == SIGNIFICANT
