"""Calibration scenarios, resampling and acceptance criteria (#182)."""

from __future__ import annotations

import pytest

from bsdbng.calibration import (
    ALPHA,
    INJECTED_TRAIT,
    SCREEN_KEEP,
    SEED,
    Outcome,
    Resample,
    Scenario,
    check_criteria,
    flip_directions,
    fraction_sign_test,
    inject_effect,
    replicate_rng,
    run_replicates,
    screened_fraction_sign_test,
    summarize,
    summarize_injection,
)
from bsdbng.guild_maps import study_direction_means
from bsdbng.table1_trait_maps import compute_direction_shifts
from bsdbng.trait_guilds import TABLE1_TRAITS, Profile

TRAITS = tuple(TABLE1_TRAITS)


def _scenario(name: str = "s") -> Scenario:
    """A pass-through scenario, for summaries built from outcomes rather than a corpus."""
    return Scenario(name, "t", lambda profiles, _rng: Resample(list(profiles)))


def _corpus(studies: int = 40) -> list[Profile]:
    """Paired profiles with varying counts, some studies holding two experiments."""
    profiles: list[Profile] = []
    for study in range(studies):
        for experiment in range(1 + study % 2):
            for direction, offset in (("inc", 3), ("dec", 1)):
                counts = {
                    trait: (study * 3 + index * 5 + offset) % 17
                    for index, trait in enumerate(TRAITS)
                }
                profiles.append(
                    Profile(
                        study_id=f"s{study}",
                        experiment_id=f"s{study}-{experiment}",
                        inc_dec=direction,
                        condition="C",
                        body_site="Feces",
                        guild_counts=counts,
                        total=sum(counts.values()) + 400,
                    )
                )
    return profiles


def test_flipping_relabels_whole_units_and_keeps_the_corpus_paired() -> None:
    profiles = _corpus()
    flipped = flip_directions(profiles, replicate_rng(182, "unit-test", 1), unit="study")
    assert len(flipped) == len(profiles)
    # Counts are untouched; only labels move.
    assert sorted(p.total for p in flipped) == sorted(p.total for p in profiles)
    for study in {p.study_id for p in profiles}:
        before = {(p.experiment_id, p.inc_dec, p.total) for p in profiles if p.study_id == study}
        after = {(p.experiment_id, p.inc_dec, p.total) for p in flipped if p.study_id == study}
        swapped = {
            (experiment, "dec" if direction == "inc" else "inc", total)
            for experiment, direction, total in before
        }
        assert after in (before, swapped)
    experiments = {(p.study_id, p.experiment_id) for p in flipped}
    for key in experiments:
        sides = {p.inc_dec for p in flipped if (p.study_id, p.experiment_id) == key}
        assert sides == {"inc", "dec"}


def test_flipping_by_experiment_can_split_a_study() -> None:
    profiles = _corpus()
    by_experiment = flip_directions(profiles, replicate_rng(182, "s", 7), unit="experiment")
    by_study = flip_directions(profiles, replicate_rng(182, "s", 7), unit="study")
    assert by_experiment != by_study
    with pytest.raises(ValueError, match="flip unit"):
        flip_directions(profiles, replicate_rng(182, "s", 1), unit="condition")


def test_injection_leaves_every_other_trait_bit_identical() -> None:
    """The ground-truth test of #293: "untouched" must mean untouched.

    The injected trait's share must rise, and every other trait's study contrast -
    the quantity the sign test scores - must be unchanged to the last bit. Raising
    the profile total instead lowers all of them, which manufactures real effects on
    the traits the false-discovery proportion counts as null.
    """
    profiles = _corpus()
    injected, effect = inject_effect(
        profiles, replicate_rng(SEED, "inject", 3), multiplier=4.0, study_share=0.5
    )
    assert effect.changed_profiles > 0
    before = study_direction_means(list(profiles), TRAITS)
    after = study_direction_means(injected, TRAITS)
    raised = 0
    for study, per_key in before.items():
        for trait in TRAITS:
            was = per_key[f"{trait}_inc"] - per_key[f"{trait}_dec"]
            now = after[study][f"{trait}_inc"] - after[study][f"{trait}_dec"]
            if trait != INJECTED_TRAIT:
                assert now == was, (study, trait)
            else:
                assert now >= was, (study, trait)
                raised += now > was
    assert raised > 0
    shifts_before = {s.trait: s for s in compute_direction_shifts(list(profiles)).shifts}
    shifts_after = {s.trait: s for s in compute_direction_shifts(injected).shifts}
    for trait in TRAITS:
        if trait == INJECTED_TRAIT:
            continue
        assert shifts_after[trait].p_value == shifts_before[trait].p_value, trait
        assert shifts_after[trait].studies_up == shifts_before[trait].studies_up, trait


def test_injection_touches_only_the_increased_side_and_holds_the_total() -> None:
    profiles = _corpus()
    injected, effect = inject_effect(
        profiles, replicate_rng(SEED, "inject", 3), multiplier=2.0, study_share=0.5
    )
    changed = [
        (before, after) for before, after in zip(profiles, injected, strict=True) if before != after
    ]
    assert changed
    assert effect.eligible_profiles >= effect.changed_profiles == len(changed)
    for before, after in changed:
        assert after.inc_dec == "inc"
        assert after.guild_counts[INJECTED_TRAIT] > before.guild_counts[INJECTED_TRAIT]
        # The counts come from the complement, so the denominator does not move.
        assert after.total == before.total
        assert sum(after.guild_counts.values()) <= after.total
        others = {t: c for t, c in after.guild_counts.items() if t != INJECTED_TRAIT}
        assert others == {t: c for t, c in before.guild_counts.items() if t != INJECTED_TRAIT}


def test_the_effective_lift_is_reported_not_just_the_nominal_one() -> None:
    """Rounding leaves small counts alone, so 1.1x is not a 1.1x lift."""
    profiles = _corpus()
    _, small = inject_effect(
        profiles, replicate_rng(SEED, "inject", 0), multiplier=1.1, study_share=1.0
    )
    _, large = inject_effect(
        profiles, replicate_rng(SEED, "inject", 0), multiplier=4.0, study_share=1.0
    )
    assert small.unchanged_after_rounding > 0
    assert small.mean_realized_multiplier is not None
    assert small.mean_realized_multiplier < small.nominal_multiplier
    rows = summarize_injection([small, small], _scenario("injected_1.1x_study_flip"))
    assert rows[0]["share_unchanged_after_rounding"] == pytest.approx(
        small.unchanged_after_rounding / small.eligible_profiles
    )
    assert float(str(rows[0]["mean_effective_multiplier"])) < 1.1
    assert large.mean_realized_multiplier is not None
    assert large.mean_realized_multiplier > small.mean_realized_multiplier


def test_the_screened_procedure_tests_only_what_it_selected() -> None:
    outcomes = screened_fraction_sign_test(_corpus())
    selected = [o for o in outcomes if o.selected]
    assert len(selected) == SCREEN_KEEP
    assert all(o.p_value is None for o in outcomes if not o.selected)
    assert len(outcomes) == len(TRAITS)


def test_a_ratio_naming_the_injected_trait_counts_as_injected() -> None:
    runs = [
        [
            Outcome("p", INJECTED_TRAIT, 0.01, 0.01),
            Outcome("p", f"{INJECTED_TRAIT}/sphere_shaped", 0.01, 0.01),
            Outcome("p", "sphere_shaped/arsenate_reduction", 0.5, 0.9),
        ]
    ]
    rows = {(r["label"], r["role"]) for r in summarize(runs, _scenario(), INJECTED_TRAIT)}
    assert (INJECTED_TRAIT, "injected_effect") in rows
    assert (f"{INJECTED_TRAIT}/sphere_shaped", "injected_effect") in rows
    assert ("sphere_shaped/arsenate_reduction", "null") in rows


def test_summary_counts_rejections_and_reports_monte_carlo_bounds() -> None:
    runs = [[Outcome("p", "a", 0.01, 0.01)], [Outcome("p", "a", 0.5, 0.9)]]
    row = next(r for r in summarize(runs, _scenario(), None) if r["label"] == "a")
    assert row["rejections_p"] == 1
    assert row["rejection_rate_p"] == pytest.approx(0.5)
    assert (
        float(str(row["rejection_rate_p_ci_low"]))
        < 0.5
        < float(str(row["rejection_rate_p_ci_high"]))
    )
    family = next(r for r in summarize(runs, _scenario(), None) if r["label"] == "FAMILY")
    assert family["rejection_rate_q"] == pytest.approx(0.5)
    assert family["mean_false_discovery_proportion"] == pytest.approx(0.5)


def test_a_failing_criterion_is_reported_as_failed() -> None:
    """An anticonservative procedure must not pass; the checker is not decorative."""
    rows = [
        {
            "scenario": "complete_null_study_flip",
            "procedure": "fraction_sign_test",
            "label": "trait",
            "role": "null",
            "rejection_rate_p": 0.4,
            "rejection_rate_p_ci_low": 0.3,
            "rejection_rate_q": 0.4,
        },
        {
            "scenario": "complete_null_study_flip",
            "procedure": "fraction_sign_test",
            "label": "FAMILY",
            "role": "family",
            "rejection_rate_p": 0.4,
            "rejection_rate_p_ci_low": 0.3,
            "rejection_rate_q": 0.4,
        },
    ]
    verdicts = {
        (r["criterion"], r["procedure"]): r["verdict"]
        for r in check_criteria(rows)
        if r["procedure"] == "fraction_sign_test"
    }
    assert verdicts[("per_test_level", "fraction_sign_test")] == "fail"
    assert verdicts[("family_error", "fraction_sign_test")] == "fail"


def test_a_short_deterministic_run_reproduces_itself() -> None:
    """The CI-sized check: a fixed seed gives the same outcomes every time."""
    profiles = _corpus(studies=30)
    scenario = Scenario(
        "complete_null_study_flip",
        "no trait differs",
        lambda profiles, rng: Resample(flip_directions(profiles, rng, unit="study")),
    )
    procedures = {"fraction_sign_test": fraction_sign_test}
    first = run_replicates(profiles, scenario, 2, seed=182, procedures=procedures)
    second = run_replicates(profiles, scenario, 2, seed=182, procedures=procedures)
    assert first == second
    assert len(first.runs) == 2
    assert {o.procedure for run in first.runs for o in run} == {"fraction_sign_test"}
    rows = summarize(first.runs, scenario, None)
    assert all(0.0 <= float(str(r["rejection_rate_p"])) <= 1.0 for r in rows)
    assert ALPHA == 0.05
