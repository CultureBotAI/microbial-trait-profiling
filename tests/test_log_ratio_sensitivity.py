"""Log-ratio sensitivity scenarios (#183)."""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from bsdbng.log_ratio_sensitivity import (
    INSUFFICIENT,
    PROTOCOL_PATH,
    TESTED,
    Row,
    complete_pairs,
    correct,
    eligibility_tables,
    leave_one_study_out,
    log_ratio_contrast,
    read_reported_taxa,
    run_scenarios,
    summarize,
)
from bsdbng.table1_trait_maps import compute_direction_shifts
from bsdbng.trait_composition import compute_log_ratio_shifts
from bsdbng.trait_guilds import TABLE1_TRAITS, Profile

TRAITS = tuple(TABLE1_TRAITS)


def _counts(seed: int, side: int) -> dict[str, int]:
    """Deterministic counts with zeros scattered through them."""
    values = {}
    for index, trait in enumerate(TRAITS):
        raw = (seed * 7 + index * 13 + side * 5) % 11
        values[trait] = 0 if raw in (0, 4) else raw * (1 + side)
    return values


def _corpus(studies: int = 12) -> list[Profile]:
    profiles: list[Profile] = []
    for study in range(studies):
        for experiment in range(1 + study % 3):
            seed = study * 3 + experiment
            site = "Feces" if study % 2 else "Saliva"
            for direction, side in (("inc", 1), ("dec", 0)):
                # One decreased profile has no retained annotation at all.
                empty = seed == 4 and side == 0
                counts = dict.fromkeys(TRAITS, 0) if empty else _counts(seed, side)
                total = 0 if empty else sum(counts.values()) + 40
                profiles.append(
                    Profile(
                        study_id=f"s{study}",
                        experiment_id=f"s{study}-{experiment}",
                        inc_dec=direction,
                        condition=f"c{study % 4}",
                        body_site=site,
                        guild_counts=counts,
                        total=total,
                    )
                )
    return profiles


def _rows(scenario: str, rows: list[Row]) -> list[Row]:
    return [row for row in rows if row.scenario == scenario]


def test_protocol_declares_every_scenario_the_code_runs() -> None:
    declared = {s["id"] for s in json.loads(PROTOCOL_PATH.read_text(encoding="utf-8"))["scenarios"]}
    pairs = complete_pairs(_corpus(), TRAITS)
    ran = {row.scenario for row in run_scenarios(pairs, TRAITS, min_studies=1)}
    assert ran | {"S8_leave_one_study_out"} == declared


def test_historical_scenario_reproduces_the_published_fraction_directions() -> None:
    profiles = _corpus()
    rows = _rows(
        "S0_historical_fraction",
        run_scenarios(complete_pairs(profiles, TRAITS), TRAITS, min_studies=1),
    )
    shifts = {shift.trait: shift for shift in compute_direction_shifts(profiles, TRAITS).shifts}
    for row in rows:
        shift = shifts[row.numerator]
        assert (row.studies_up, row.studies_down, row.studies_tied) == (
            shift.studies_up,
            shift.studies_down,
            shift.studies_tied,
        )
        if row.q_value is not None:
            assert row.q_value == pytest.approx(shift.q_value)


def test_complete_case_scenario_reproduces_the_log_ratio_module() -> None:
    profiles = _corpus()
    rows = _rows(
        "S2_ratio_complete_case",
        run_scenarios(complete_pairs(profiles, TRAITS), TRAITS, min_studies=1),
    )
    summaries = {
        (s.numerator, s.reference): s for s in compute_log_ratio_shifts(profiles).summaries
    }
    assert len(rows) == 15
    for row in rows:
        summary = summaries[(row.numerator, row.reference)]
        assert row.studies_eligible == summary.studies_eligible
        assert row.pairs_eligible == summary.eligible_experiments
        assert row.studies_partially_observed == summary.studies_partially_observed
        assert (row.studies_up, row.studies_down) == (summary.studies_up, summary.studies_down)
        if summary.by_q is not None and row.q_value is not None:
            assert row.q_value == pytest.approx(summary.by_q)


def test_matched_fractions_use_exactly_the_ratio_pairs() -> None:
    pairs = complete_pairs(_corpus(), TRAITS)
    rows = run_scenarios(pairs, TRAITS, min_studies=1)
    ratio = next(
        r
        for r in _rows("S2_ratio_complete_case", rows)
        if r.numerator == TRAITS[0] and r.reference == TRAITS[1]
    )
    matched = [
        r
        for r in _rows("S3_ratio_matched_fractions", rows)
        if r.reference == f"on {TRAITS[0]}/{TRAITS[1]} pairs"
    ]
    assert [r.numerator for r in matched] == [TRAITS[0], TRAITS[1]]
    for row in matched:
        assert row.pairs_eligible == ratio.pairs_eligible
        assert row.studies_eligible == ratio.studies_eligible


def test_zero_total_pairs_are_scored_historically_and_excluded_at_baseline() -> None:
    pairs = complete_pairs(_corpus(), TRAITS)
    rows = run_scenarios(pairs, TRAITS, min_studies=1)
    historical = _rows("S0_historical_fraction", rows)[0]
    baseline = _rows("S1_fraction_baseline", rows)[0]
    zero_total = sum(1 for p in pairs if p.inc.total <= 0 or p.dec.total <= 0)
    assert zero_total > 0
    assert historical.pairs_eligible - baseline.pairs_eligible == zero_total


def test_common_support_scenarios_share_one_population() -> None:
    pairs = complete_pairs(_corpus(studies=40), TRAITS)
    rows = _rows("S4_six_trait_common_support", run_scenarios(pairs, TRAITS, min_studies=1))
    assert len({(r.pairs_total, r.pairs_eligible) for r in rows}) == 1


def test_small_populations_are_described_but_not_tested() -> None:
    pairs = complete_pairs(_corpus(), TRAITS)
    row = summarize(
        scenario="x",
        comparison="log_ratio",
        numerator=TRAITS[0],
        reference=TRAITS[1],
        population="all",
        pairs=pairs,
        contrast=log_ratio_contrast(TRAITS[0], TRAITS[1]),
        min_studies=1000,
    )
    assert row.status == INSUFFICIENT
    assert row.p_value is None
    assert row.share_up is None
    assert row.mean_effect is not None


def test_untested_rows_keep_their_family_place_without_a_q() -> None:
    pairs = complete_pairs(_corpus(), TRAITS)
    tested = summarize(
        scenario="x",
        comparison="log_ratio",
        numerator=TRAITS[0],
        reference=TRAITS[2],
        population="all",
        pairs=pairs,
        contrast=log_ratio_contrast(TRAITS[0], TRAITS[2]),
        min_studies=1,
    )
    untested = summarize(
        scenario="x",
        comparison="log_ratio",
        numerator=TRAITS[0],
        reference=TRAITS[1],
        population="all",
        pairs=pairs,
        contrast=log_ratio_contrast(TRAITS[0], TRAITS[1]),
        min_studies=1000,
    )
    assert tested.status == TESTED
    corrected = correct([tested, untested], "f", "by")
    assert [row.family_size for row in corrected] == [2, 2]
    assert corrected[1].q_value is None
    assert corrected[0].q_value is not None
    assert corrected[0].q_value >= (tested.p_value or 0.0)


def test_leave_one_study_out_ranges_bracket_the_full_share() -> None:
    pairs = complete_pairs(_corpus(studies=30), TRAITS)
    for item in leave_one_study_out(pairs, TRAITS):
        if item.full_share_up is None or item.min_share_up is None or item.max_share_up is None:
            continue
        assert item.min_share_up <= item.full_share_up + 1e-12
        assert item.max_share_up >= item.full_share_up - 1e-12
        assert item.removals_changing_q_status >= 0


def test_eligibility_partitions_every_study(tmp_path: Path) -> None:
    pairs = complete_pairs(_corpus(), TRAITS)
    taxa = tmp_path / "taxa.tsv"
    taxa.write_text(
        "study_id\texperiment_id\tsignature_id\tdirection\ttaxon_id\ttaxon_name\ttaxonomic_rank\n"
        "s1\ts1-0\tx\tincreased\tNCBITaxon:1\ta\tspecies\n"
        "s1\ts1-0\tx\tdecreased\tNCBITaxon:2\tb\tgenus\n",
        encoding="utf-8",
    )
    reported = read_reported_taxa(taxa)
    assert (reported[("s1", "s1-0")].taxa, reported[("s1", "s1-0")].species) == (2, 1)
    overall, strata = eligibility_tables(pairs, TRAITS, reported)
    for row in overall:
        eligible = int(str(row["studies_eligible"]))
        excluded = int(str(row["studies_entirely_excluded"]))
        assert eligible + excluded == int(str(row["studies_total"]))
        assert row["coverage_sidecar_status"] == "pending #184"
    sites = [
        r
        for r in strata
        if r["stratum_kind"] == "body_site"
        and r["numerator"] == TRAITS[0]
        and r["reference"] == TRAITS[1]
    ]
    assert sum(int(str(r["pairs_total"])) for r in sites) == len(pairs)
    assert not math.isnan(float(str(overall[0]["pairs_eligible"])))
