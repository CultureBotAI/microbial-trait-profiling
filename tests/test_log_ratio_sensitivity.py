"""Log-ratio sensitivity scenarios (#183)."""

from __future__ import annotations

import json
import math
from dataclasses import replace
from itertools import combinations
from pathlib import Path

import pytest

from bsdbng.log_ratio_sensitivity import (
    FECES_ANY,
    INSUFFICIENT,
    NO_GUT_LABEL,
    OTHER_GUT,
    PROTOCOL_PATH,
    TESTED,
    Row,
    _detection_change,
    _estimate_agreement,
    _held_family_q,
    _intervals_overlap,
    _stratum_of,
    complete_pairs,
    correct,
    eligibility_tables,
    label_accounting,
    leave_one_stratum_out,
    leave_one_study_out,
    log_ratio_contrast,
    matched_reference,
    population_comparison,
    read_reported_taxa,
    run_scenarios,
    sign_test,
    study_means,
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


# Share of studies whose increased side carries the larger count, per trait. The first
# two traits carry a strong and opposite signal so that corrected q values land well
# below 1 and a BY family can be told apart from a BH one; the rest sit near the null.
_SIGNAL = (0.78, 0.22, 0.58, 0.50, 0.45, 0.53)
_SITES = ("Feces", "Feces", "Saliva", "Feces,Blood", "Rectum", "Vagina")


def _signal_corpus(studies: int = 96) -> list[Profile]:
    """A corpus with real signal, scattered zeros, several sites and several conditions.

    Without a signal every corrected q collapses to 1, which makes a q assertion
    vacuous and lets a BY family pass as a BH one. The direction of a trait is a
    deterministic function of the study, so every experiment of a study agrees and the
    study-level sign is unambiguous.
    """
    profiles: list[Profile] = []
    for study in range(studies):
        for experiment in range(1 + study % 2):
            seed = study * 2 + experiment
            counts: dict[str, dict[str, int]] = {"inc": {}, "dec": {}}
            for index, trait in enumerate(TRAITS):
                positive = (study * 37 + index * 11) % 100 < _SIGNAL[index] * 100
                low = 4 + (seed * 3 + index * 7) % 6
                high = low + 4
                counts["inc"][trait] = high if positive else low
                counts["dec"][trait] = low if positive else high
                # Scattered zeros, at different rates per side, so each ratio keeps a
                # different population and some studies are only partially observed.
                if (seed + index * 5) % 11 == 0:
                    counts["inc"][trait] = 0
                if (seed + index * 7) % 13 == 0:
                    counts["dec"][trait] = 0
            for direction in ("inc", "dec"):
                # One decreased profile has no retained annotation at all.
                empty = seed == 7 and direction == "dec"
                values = dict.fromkeys(TRAITS, 0) if empty else counts[direction]
                profiles.append(
                    Profile(
                        study_id=f"s{study:03d}",
                        experiment_id=f"s{study:03d}-{experiment}",
                        inc_dec=direction,
                        condition=f"c{study % 5}",
                        body_site=_SITES[study % len(_SITES)],
                        guild_counts=values,
                        total=0 if empty else sum(values.values()) + 50,
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
        assert item.max_abs_share_shift == pytest.approx(
            max(
                abs(item.min_share_up - item.full_share_up),
                abs(item.max_share_up - item.full_share_up),
            ),
            abs=1e-12,
        )


def _ratio_index(ratio: tuple[str, str]) -> int:
    return list(combinations(TRAITS, 2)).index(ratio)


def _removal_states(up: int, down: int, count: int) -> tuple[tuple[int, int], ...]:
    """The up/down counts reachable by removing ``count`` studies from one sign group."""
    return ((up - count, down), (up, down - count))


def _flip_states(up: int, down: int, count: int) -> tuple[tuple[int, int], ...]:
    """The up/down counts reachable by reversing ``count`` study signs."""
    return ((up - count, down + count), (up + count, down - count))


def _crosses_at(
    base: list[float], position: int, up: int, down: int, report_q: float
) -> bool | None:
    """Whether these counts put the held family's q at or below the rule."""
    if up < 0 or down < 0 or up + down == 0:
        return None
    p_value = sign_test([1.0] * up + [-1.0] * down)[3]
    assert p_value is not None
    return _held_family_q(base, position, p_value) <= report_q


def test_removal_flip_counts_are_the_minima_that_cross_the_rule() -> None:
    """A removal or flip count below the reported minimum must not cross the rule.

    The decision threshold is placed inside the range of q values that single
    removals produce, so at least one ratio has to cross at one removal. That makes
    the flip count non-zero, so forcing it to zero cannot pass.
    """
    pairs = complete_pairs(_signal_corpus(), TRAITS)
    ratio_pairs = list(combinations(TRAITS, 2))
    means = {r: study_means(pairs, log_ratio_contrast(*r)) for r in ratio_pairs}
    base = [sign_test(v.values())[3] or 1.0 for v in means.values()]

    # A threshold strictly between the full-data q and the largest single-removal q
    # of one ratio, so that removing one study does change that ratio's status.
    target = ratio_pairs[0]
    index = _ratio_index(target)
    full_q = _held_family_q(base, index, base[index])
    single = []
    for study in sorted(means[target]):
        rest = [v for other, v in means[target].items() if other != study]
        up, down, _t, p = sign_test(rest)
        if p is not None:
            single.append(_held_family_q(base, index, p))
    assert max(single) > full_q, "fixture must have a removal that raises q"
    report_q = (full_q + max(single)) / 2

    items = leave_one_study_out(pairs, TRAITS, report_q=report_q)
    by_ratio = {(i.numerator, i.reference): i for i in items}
    crossing = by_ratio[target]

    assert crossing.decision_rule_q == report_q
    assert crossing.family_size == 15
    assert crossing.correction == "Benjamini-Yekutieli"
    assert crossing.removals_changing_q_status > 0
    assert crossing.min_removals_to_cross == 1

    for item in items:
        ratio = (item.numerator, item.reference)
        position = _ratio_index(ratio)
        up, down, _tied, _p = sign_test(means[ratio].values())
        significant = _held_family_q(base, position, base[position]) <= report_q
        # A non-zero leave-one-out flip count and a minimum of one are the same fact.
        assert (item.removals_changing_q_status > 0) == (item.min_removals_to_cross == 1)
        for count, states in (
            (item.min_removals_to_cross, _removal_states),
            (item.min_flips_to_cross, _flip_states),
        ):
            if count is None:
                continue
            assert count >= 1
            # Nothing below the reported minimum may cross.
            for smaller in range(1, count):
                for r_up, r_down in states(up, down, smaller):
                    crossed = _crosses_at(base, position, r_up, r_down, report_q)
                    assert crossed is None or crossed == significant, (
                        f"{ratio} crosses at {smaller} < {count}"
                    )
            # The reported minimum itself must cross on one of the two paths.
            assert any(
                _crosses_at(base, position, r_up, r_down, report_q) not in (None, significant)
                for r_up, r_down in states(up, down, count)
            ), f"{ratio} does not cross at its reported minimum {count}"


def test_leave_one_stratum_out_only_tests_strata_above_the_minimum() -> None:
    """Each reported stratum removal must drop exactly that stratum's studies."""
    pairs = complete_pairs(_signal_corpus(), TRAITS)
    rows, skipped = leave_one_stratum_out(pairs, TRAITS, min_studies=4)
    assert rows, "fixture must produce testable strata"

    # The minimum is a real filter: raising it past every stratum reports nothing and
    # accounts for all of them as skipped.
    empty, all_skipped = leave_one_stratum_out(pairs, TRAITS, min_studies=10_000)
    assert empty == []
    assert all_skipped == len(rows) + skipped > 0

    for item in rows:
        contrast = log_ratio_contrast(item.numerator, item.reference)
        eligible = [p for p in pairs if contrast(p) is not None]
        inside = {p.study_id for p in eligible if _stratum_of(p, item.stratum_kind) == item.stratum}
        assert len(inside) == item.studies_eligible_in_stratum
        assert item.studies_eligible_in_stratum >= 4
        kept = [p for p in pairs if _stratum_of(p, item.stratum_kind) != item.stratum]
        assert item.studies_eligible_remaining == len(study_means(kept, contrast))
        assert item.studies_eligible_remaining < len(study_means(pairs, contrast))
        assert item.family_size == 15
        assert item.decision_rule_q == pytest.approx(0.05)

    assert {item.stratum_kind for item in rows} <= {"body_site", "condition", "gut_label_set"}


def test_label_set_strata_put_a_multi_label_feces_pair_with_feces() -> None:
    """A cell naming Feces alongside another site is fecal; the exact-label split is not."""
    pairs = complete_pairs(_signal_corpus(), TRAITS)
    by_site = {pair.body_site: pair for pair in pairs}
    assert "Feces,Blood" in by_site, "fixture must carry a multi-label Feces cell"

    assert by_site["Feces"].gut_stratum() == FECES_ANY
    assert by_site["Feces,Blood"].gut_stratum() == FECES_ANY
    assert by_site["Feces,Blood"].labels() == {"Feces", "Blood"}
    assert by_site["Rectum"].gut_stratum() == OTHER_GUT
    assert by_site["Saliva"].gut_stratum() == NO_GUT_LABEL
    assert by_site["Vagina"].gut_stratum() == NO_GUT_LABEL

    # S7, matching the whole cell, keeps the multi-label fecal pairs and the other gut.
    counts = label_accounting(pairs)
    assert counts["feces_labelled_pairs_inside_s7"] > 0
    assert counts["other_gut_pairs_inside_s7"] > 0
    assert counts["multi_label_pairs"] > 0
    assert (
        counts["exact_label_feces_pairs"] + counts["pairs_outside_exact_label_feces"]
        == counts["pairs_total"]
        == len(pairs)
    )

    # The three label-set strata partition every pair exactly once.
    strata = [FECES_ANY, OTHER_GUT, NO_GUT_LABEL]
    sizes = {name: sum(pair.gut_stratum() == name for pair in pairs) for name in strata}
    assert sum(sizes.values()) == len(pairs)
    assert all(size > 0 for size in sizes.values())


def test_label_set_scenarios_run_on_their_own_stratum() -> None:
    """Each of S9, S10 and S11 is computed on its stratum, not on every pair."""
    pairs = complete_pairs(_signal_corpus(), TRAITS)
    rows = run_scenarios(pairs, TRAITS, min_studies=1)
    for scenario, stratum in (
        ("S9_feces_any_label", FECES_ANY),
        ("S10_other_gut_no_feces", OTHER_GUT),
        ("S11_no_gut_label", NO_GUT_LABEL),
    ):
        scenario_rows = _rows(scenario, rows)
        expected = sum(pair.gut_stratum() == stratum for pair in pairs)
        assert len(scenario_rows) == 15
        assert {row.pairs_total for row in scenario_rows} == {expected}
        assert expected < len(pairs)
        assert {row.family for row in scenario_rows} == {f"{scenario.split('_')[0]}:fifteen_ratios"}


def test_population_comparison_pairs_each_matched_row_with_its_baseline() -> None:
    """Every S3 row is compared against the S1 row for the same trait, not another."""
    pairs = complete_pairs(_corpus(studies=30), TRAITS)
    rows = run_scenarios(pairs, TRAITS, min_studies=1)
    baseline = {r.numerator: r for r in _rows("S1_fraction_baseline", rows)}
    matched = {(r.numerator, r.reference): r for r in _rows("S3_ratio_matched_fractions", rows)}
    comparisons = population_comparison(rows, TRAITS)

    assert len(comparisons) == len(matched) == 30
    for item in comparisons:
        base = baseline[item.trait]
        rest = matched[(item.trait, matched_reference(item.ratio_numerator, item.ratio_reference))]
        assert item.baseline_studies_eligible == base.studies_eligible
        assert item.restricted_studies_eligible == rest.studies_eligible
        assert item.baseline_share_up == base.share_up
        assert item.restricted_share_up == rest.share_up
        assert item.baseline_family == base.family
        assert item.restricted_family == rest.family
        assert item.baseline_family_size != item.restricted_family_size
        if base.share_up is not None and rest.share_up is not None:
            assert item.share_shift == pytest.approx(rest.share_up - base.share_up, abs=1e-12)


def test_detection_change_separates_the_correction_from_the_raw_test() -> None:
    """A row with a raw p under the rule but q over it is attributed to the correction."""
    base = Row(
        scenario="S1",
        comparison="fraction",
        numerator="t",
        reference="",
        population="all",
        effect_units="u",
        studies_total=100,
        studies_eligible=100,
        studies_partially_observed=0,
        pairs_total=100,
        pairs_eligible=100,
        studies_up=60,
        studies_down=40,
        studies_tied=0,
        mean_effect=0.1,
        median_effect=0.1,
        share_up=0.6,
        ci_low=0.5,
        ci_high=0.7,
        p_value=0.001,
        status=TESTED,
        q_value=0.004,
    )
    lost_at_correction = replace(base, scenario="S3", p_value=0.02, q_value=0.4)
    lost_before = replace(base, scenario="S3", p_value=0.3, q_value=0.9)
    kept = replace(base, scenario="S3", p_value=0.001, q_value=0.01)

    assert _detection_change(base, lost_at_correction) == "lost at correction"
    assert _detection_change(base, lost_before) == "lost before correction"
    assert _detection_change(base, kept) == "retained"
    assert _detection_change(lost_before, kept) == "gained"
    assert _detection_change(lost_before, lost_before) == "absent in both"


def test_estimate_agreement_reports_a_moved_interval_even_when_the_side_flips() -> None:
    """A row whose interval has moved is not recorded only as a change of side."""
    low = Row(
        scenario="S1",
        comparison="fraction",
        numerator="t",
        reference="",
        population="all",
        effect_units="u",
        studies_total=10,
        studies_eligible=10,
        studies_partially_observed=0,
        pairs_total=10,
        pairs_eligible=10,
        studies_up=4,
        studies_down=6,
        studies_tied=0,
        mean_effect=-0.1,
        median_effect=-0.1,
        share_up=0.40,
        ci_low=0.30,
        ci_high=0.48,
        p_value=0.2,
        status=TESTED,
    )
    high_overlapping = replace(low, share_up=0.55, ci_low=0.45, ci_high=0.65)
    high_disjoint = replace(low, share_up=0.80, ci_low=0.70, ci_high=0.90)
    same_side_disjoint = replace(low, share_up=0.20, ci_low=0.10, ci_high=0.25)

    assert _estimate_agreement(low, high_overlapping, True) == "side reversed, intervals overlap"
    assert _estimate_agreement(low, high_disjoint, False) == "shifted and side reversed"
    assert _estimate_agreement(low, same_side_disjoint, False) == "shifted"
    assert _estimate_agreement(low, replace(low), True) == "consistent"
    assert _intervals_overlap(low, high_overlapping) is True
    assert _intervals_overlap(low, high_disjoint) is False


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
