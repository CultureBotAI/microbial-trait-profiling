"""Calibration scenarios, resampling and acceptance criteria (#182, #293)."""

from __future__ import annotations

import json
from collections.abc import Sequence

import numpy as np
import pytest

from bsdbng.calibration import (
    ALPHA,
    INJECTED_TRAIT,
    INJECTION_MULTIPLIERS,
    INJECTION_STUDY_SHARE,
    NULL_REPLICATES,
    POWER_REPLICATES,
    PROTOCOL_PATH,
    SCREEN_KEEP,
    SEED,
    Outcome,
    Resample,
    Scenario,
    attainable_p_at_alpha,
    check_criteria,
    criterion_pass_probability,
    flip_assignments,
    flip_directions,
    fraction_sign_test,
    inject_effect,
    load_protocol,
    pooled_rate,
    replicate_counts,
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
#: The trials the declared null budget gives the least-powered procedure: one test
#: per replicate, over both complete-null scenarios.
DECLARED_TRIALS = 2 * NULL_REPLICATES


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


def _null_ledger(
    rejections: int,
    trials_per_replicate: int,
    replicates: int,
    procedure: str = "fraction_sign_test",
    *,
    families: int = 0,
) -> list[dict[str, object]]:
    """A replicate ledger holding exactly ``rejections`` rejections, spread evenly.

    Spreading them one per replicate keeps the measured design effect near one, so
    the criterion sees the same uncertainty a single-test procedure would give it.
    """
    rows: list[dict[str, object]] = []
    for index in range(replicates):
        scenario = "complete_null_study_flip" if index % 2 else "complete_null_experiment_flip"
        rows.append(
            {
                "scenario": scenario,
                "procedure": procedure,
                "replicate": index,
                "tests_evaluated": trials_per_replicate,
                "null_tests_evaluated": trials_per_replicate,
                "null_tests_rejected_p": 1 if index < rejections else 0,
                "family_rejected_q": 1 if index < families else 0,
            }
        )
    return rows


def _verdicts(
    rows: Sequence[dict[str, object]], ledger: Sequence[dict[str, object]]
) -> dict[tuple[str, str], str]:
    return {
        (str(r["criterion"]), str(r["procedure"])): str(r["verdict"])
        for r in check_criteria(rows, ledger)
    }


def test_flipping_relabels_whole_units_and_keeps_the_corpus_paired() -> None:
    profiles = _corpus()
    flipped = flip_directions(profiles, replicate_rng(SEED, "unit-test", 1), unit="study")
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
    by_experiment = flip_directions(profiles, replicate_rng(SEED, "s", 7), unit="experiment")
    by_study = flip_directions(profiles, replicate_rng(SEED, "s", 7), unit="study")
    assert by_experiment != by_study
    with pytest.raises(ValueError, match="flip unit"):
        flip_directions(profiles, replicate_rng(SEED, "s", 1), unit="condition")


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
    # Power is reported at q, so q carries its own interval (#293).
    assert (
        float(str(row["rejection_rate_q_ci_low"]))
        < 0.5
        < float(str(row["rejection_rate_q_ci_high"]))
    )
    assert row["smallest_observed_p"] == pytest.approx(0.01)
    family = next(r for r in summarize(runs, _scenario(), None) if r["label"] == "FAMILY")
    assert family["rejection_rate_q"] == pytest.approx(0.5)
    assert family["mean_false_discovery_proportion"] == pytest.approx(0.5)


def test_the_smallest_observed_p_is_not_the_attainable_resolution() -> None:
    """Two different quantities, which v1 reported under one name (#293)."""
    runs = [[Outcome("p", "a", 0.02, 0.02, scored=40)]]
    row = next(r for r in summarize(runs, _scenario(), None) if r["label"] == "a")
    assert row["smallest_observed_p"] == pytest.approx(0.02)
    below, above = attainable_p_at_alpha(40)
    assert below is not None
    assert above is not None
    assert below <= ALPHA < above
    assert row["attainable_p_below_alpha"] == pytest.approx(below)
    assert row["attainable_p_above_alpha"] == pytest.approx(above)
    # The observed p need not be attainable at all for this corpus size.
    assert below != pytest.approx(0.02)
    assert attainable_p_at_alpha(0) == (None, None)


def test_the_pooled_rate_charges_for_within_replicate_dependence() -> None:
    """Six labels that always agree carry the information of one, not six."""
    independent = pooled_rate([(1, 6)] * 20 + [(0, 6)] * 380)
    lockstep = pooled_rate([(6, 6)] * 20 + [(0, 6)] * 380)
    assert independent is not None
    assert lockstep is not None
    assert independent.trials == lockstep.trials == 2400
    assert lockstep.rate > independent.rate
    assert lockstep.design_effect > 5.0
    assert independent.design_effect < 1.5
    assert pooled_rate([(0, 0)]) is None


def test_the_amended_criterion_fails_a_rate_the_old_rule_accepted() -> None:
    """The heart of #293: at the declared budget 0.06 to 0.08 must all fail.

    v1 passed a label whose exact two-sided 95% lower bound was at most 0.05, which
    accepted 16/200 = 0.08 outright. Both rules are checked at the budget that
    separates them, because at 2400 trials even the old rule rejects 0.08 - the v1
    rule was anticonservative at the v1 budget, and the two defects travel together.
    """
    cases = {0.05: "pass", 0.06: "fail", 0.07: "fail", 0.08: "fail"}
    for rate, expected in cases.items():
        ledger = _null_ledger(round(rate * DECLARED_TRIALS), 1, DECLARED_TRIALS)
        verdict = _verdicts([], ledger)[("per_test_level", "fraction_sign_test")]
        assert verdict == expected, f"true rate {rate} gave {verdict}, expected {expected}"


def test_the_old_budget_cannot_demonstrate_control_at_any_rate() -> None:
    """At 200 replicates the v1 rule passed 16/200; the amended rule passes nothing.

    This is the test that pins the rule rather than the budget: 0.08 over 200 trials
    has an exact 95% lower bound of 0.047, so v1 called it a pass. The amended rule
    fails it, and fails a true 0.05 as well - 200 replicates are not enough to
    demonstrate control, which is why the declared budget rose to 1200.
    """
    for rate in (0.05, 0.08):
        ledger = _null_ledger(round(rate * 200), 1, 200)
        verdict = _verdicts([], ledger)[("per_test_level", "fraction_sign_test")]
        assert verdict == "fail", f"true rate {rate} at 200 trials gave {verdict}"


def test_the_criterion_power_meets_the_declared_requirement() -> None:
    """The reported power calculation the protocol promises, asserted."""
    declared = load_protocol()["criterion_calibration"]
    tolerance = float(declared["per_test_tolerance"])
    required = declared["required_operating_characteristics"]
    assert criterion_pass_probability(0.05, DECLARED_TRIALS, tolerance) >= float(
        required["pass_probability_at_0.05"]
    )
    assert 1.0 - criterion_pass_probability(0.07, DECLARED_TRIALS, tolerance) >= float(
        required["fail_probability_at_0.07"]
    )
    assert 1.0 - criterion_pass_probability(0.08, DECLARED_TRIALS, tolerance) >= float(
        required["fail_probability_at_0.08"]
    )


def test_the_criterion_power_formula_matches_a_simulation() -> None:
    """The analytic operating characteristics, checked against sampled ledgers."""
    rng = np.random.default_rng([SEED, 293])
    tolerance = 0.065
    for rate in (0.05, 0.07):
        draws = rng.binomial(DECLARED_TRIALS, rate, size=400)
        failures = 0
        for successes in draws:
            pooled = pooled_rate(
                [(1, 1)] * int(successes) + [(0, 1)] * (DECLARED_TRIALS - int(successes))
            )
            assert pooled is not None
            failures += pooled.upper_bound > tolerance
        empirical = failures / len(draws)
        analytic = 1.0 - criterion_pass_probability(rate, DECLARED_TRIALS, tolerance)
        assert empirical == pytest.approx(analytic, abs=0.06)
        if rate == 0.07:
            assert empirical >= 0.9


def test_the_family_criterion_fails_an_inflated_rate_and_accounts_for_error() -> None:
    """v1 compared a point estimate to 0.10, where BH claims 0.05 (#293)."""
    for rate, expected in {0.05: "pass", 0.08: "fail", 0.10: "fail"}.items():
        ledger = _null_ledger(0, 1, DECLARED_TRIALS, families=round(rate * DECLARED_TRIALS))
        verdict = _verdicts([], ledger)[("family_error", "fraction_sign_test")]
        assert verdict == expected, f"family rate {rate} gave {verdict}"


def test_the_guild_and_screened_procedures_are_in_the_applicable_checks() -> None:
    """Both were missing from the family check; the screened one is the control."""
    declared = {c["id"]: c for c in load_protocol()["acceptance_criteria"]}
    for criterion in ("per_test_level", "family_error"):
        assert "guild_mcnemar" in declared[criterion]["applies_to"], criterion
        assert "screened_fraction_sign_test" in declared[criterion]["applies_to"], criterion
    expected = declared["per_test_level"]["expected_verdict"]
    assert expected["screened_fraction_sign_test"] == "fail"


def test_a_negative_control_that_passes_is_flagged_as_unmet() -> None:
    """The control exists to show the criterion has teeth; silence would hide it."""
    ledger = _null_ledger(
        round(0.02 * DECLARED_TRIALS),
        1,
        DECLARED_TRIALS,
        procedure="screened_fraction_sign_test",
    )
    row = next(
        r
        for r in check_criteria([], ledger)
        if r["criterion"] == "per_test_level" and r["procedure"] == "screened_fraction_sign_test"
    )
    assert row["verdict"] == "pass"
    assert row["expected_verdict"] == "fail"
    assert row["expectation_met"] is False


def test_power_is_not_reported_as_passing_without_an_interval() -> None:
    """v1 passed whenever any injected row existed (#293)."""
    scenarios = [f"injected_{m:g}x_study_flip" for m in INJECTION_MULTIPLIERS]
    complete = [
        {
            "scenario": name,
            "procedure": "fraction_sign_test",
            "label": INJECTED_TRAIT,
            "role": "injected_effect",
            "replicates": POWER_REPLICATES,
            "rejection_rate_q": 0.5,
            "rejection_rate_q_ci_low": 0.45,
            "rejection_rate_q_ci_high": 0.55,
        }
        for name in scenarios
    ]
    assert _verdicts(complete, [])[("power_reported", "fraction_sign_test")] == "pass"

    wide = [dict(row) for row in complete]
    wide[0]["rejection_rate_q_ci_low"] = 0.1
    wide[0]["rejection_rate_q_ci_high"] = 0.9
    assert _verdicts(wide, [])[("power_reported", "fraction_sign_test")] == "fail"

    missing = [dict(row) for row in complete]
    missing[1]["rejection_rate_q_ci_high"] = None
    assert _verdicts(missing, [])[("power_reported", "fraction_sign_test")] == "fail"

    assert _verdicts(complete[:1], [])[("power_reported", "fraction_sign_test")] == "fail"


def test_a_short_deterministic_run_reproduces_itself() -> None:
    """The CI-sized check: a fixed seed gives the same outcomes every time."""
    profiles = _corpus(studies=30)
    scenario = Scenario(
        "complete_null_study_flip",
        "no trait differs",
        lambda profiles, rng: Resample(flip_directions(profiles, rng, unit="study")),
    )
    procedures = {"fraction_sign_test": fraction_sign_test}
    first = run_replicates(profiles, scenario, 2, seed=SEED, procedures=procedures)
    second = run_replicates(profiles, scenario, 2, seed=SEED, procedures=procedures)
    assert first == second
    assert len(first.runs) == 2
    assert {o.procedure for run in first.runs for o in run} == {"fraction_sign_test"}
    rows = summarize(first.runs, scenario, None)
    assert all(0.0 <= float(str(r["rejection_rate_p"])) <= 1.0 for r in rows)
    ledger = replicate_counts(first.runs, scenario, None)
    assert len(ledger) == 2
    assert all(r["null_tests_evaluated"] == len(TRAITS) for r in ledger)
    assert ALPHA == 0.05


def test_a_replicate_stream_is_addressable_without_replaying_it() -> None:
    """#263 needs replicate 7 on its own, not the seventh draw of a chain."""
    direct = flip_directions(_corpus(), replicate_rng(SEED, "shared", 7), unit="study")
    again = flip_directions(_corpus(), replicate_rng(SEED, "shared", 7), unit="study")
    other = flip_directions(_corpus(), replicate_rng(SEED, "shared", 8), unit="study")
    assert direct == again
    assert direct != other


def test_the_shared_flip_artefact_matches_the_scenario_it_names() -> None:
    """Replicate i of the artefact must be replicate i of the null it is named for."""
    profiles = _corpus(studies=12)
    rows = flip_assignments(profiles, replicates=3, seed=SEED)
    studies = sorted({p.study_id for p in profiles})
    assert len(rows) == 3 * len(studies)
    assert {r["experiment_id"] for r in rows} == {None}
    for replicate in range(3):
        declared = {
            str(r["study_id"]) for r in rows if r["replicate"] == replicate and r["flipped"] == 1
        }
        flipped = flip_directions(
            profiles,
            replicate_rng(SEED, "complete_null_study_flip", replicate),
            unit="study",
        )
        moved = {
            before.study_id
            for before, after in zip(profiles, flipped, strict=True)
            if before.inc_dec != after.inc_dec
        }
        assert declared == moved


def test_module_constants_match_the_declared_protocol() -> None:
    """One declaration, not two: the code reads its budget from the protocol."""
    declared = json.loads(PROTOCOL_PATH.read_text(encoding="utf-8"))
    assert declared["replicates"]["null"] == NULL_REPLICATES
    assert declared["replicates"]["power"] == POWER_REPLICATES
    assert declared["replicates"]["seed"] == SEED
    assert list(INJECTION_MULTIPLIERS) == declared["signal_injection"]["multipliers"]
    assert declared["signal_injection"]["study_share"] == INJECTION_STUDY_SHARE
    assert declared["signal_injection"]["injected_trait"] == INJECTED_TRAIT
    assert "ci_check" not in declared["replicates"]
