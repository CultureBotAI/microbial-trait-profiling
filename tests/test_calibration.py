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
    _family_bound,
    _pooled_bound,
    attainable_p_at_alpha,
    calibration_status,
    check_criteria,
    criterion_pass_probability,
    exact_sign_test_size,
    flip_assignments,
    flip_directions,
    fraction_sign_test,
    inject_effect,
    largest_passing_count,
    load_protocol,
    pooled_rate,
    replicate_counts,
    replicate_rng,
    run_calibration,
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
    fails it, because its estimate already exceeds the tolerance, and withholds any
    verdict on a true 0.05 - 200 replicates are not enough to demonstrate control,
    which is why the declared budget rose to 1200.
    """
    cases = {0.05: "insufficient_budget", 0.08: "fail"}
    for rate, expected in cases.items():
        ledger = _null_ledger(round(rate * 200), 1, 200)
        verdict = _verdicts([], ledger)[("per_test_level", "fraction_sign_test")]
        assert verdict == expected, f"true rate {rate} at 200 trials gave {verdict}"


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


def test_zero_rejections_in_ten_trials_do_not_certify_a_zero_rate() -> None:
    """The #240 counterexample: guild_mcnemar, 0/10 null rejections, bound reported as 0.

    The clustered variance is zero when no replicate rejects, so the normal bound
    collapsed onto the point estimate. The exact bound at ten trials is 0.2589, the
    same number the family criterion reported for the same ten replicates.
    """
    pooled = pooled_rate([(0, 1)] * 10)
    assert pooled is not None
    assert pooled.rate == 0.0
    assert pooled.effective_trials == 10
    assert pooled.upper_bound == pytest.approx(0.258866, abs=1e-6)


def test_a_zero_count_with_correlated_labels_counts_each_replicate_once() -> None:
    """Six labels per replicate and no rejection: the design effect is unmeasurable."""
    pooled = pooled_rate([(0, 6)] * 10)
    assert pooled is not None
    assert pooled.trials == 60
    assert pooled.effective_trials == 10
    assert pooled.upper_bound == pytest.approx(0.258866, abs=1e-6)


def test_the_quick_run_ledger_cannot_export_a_pass() -> None:
    """The exact shape of the reviewed quick run: 5 replicates x 2 null scenarios."""
    ledger = _null_ledger(0, 1, 10, procedure="guild_mcnemar")
    criteria = check_criteria([], ledger)
    verdicts = {
        str(r["criterion"]): str(r["verdict"])
        for r in criteria
        if r["procedure"] == "guild_mcnemar"
    }
    assert verdicts == {
        "per_test_level": "insufficient_budget",
        "family_error": "insufficient_budget",
    }
    status = calibration_status(criteria)["guild_mcnemar"]
    assert status["overall"] == "not_certified"


def test_a_pass_below_the_declared_budget_is_withheld() -> None:
    """A clean bound from too few replicates is not a certificate."""
    short = _null_ledger(0, 1, 400)
    full = _null_ledger(round(0.04 * DECLARED_TRIALS), 1, DECLARED_TRIALS)
    short_rows = {
        r["criterion"]: r
        for r in check_criteria([], short)
        if r["procedure"] == "fraction_sign_test"
    }
    assert short_rows["per_test_level"]["verdict"] == "insufficient_budget"
    assert short_rows["per_test_level"]["replicates_recorded"] == 400
    assert short_rows["per_test_level"]["replicates_required"] == DECLARED_TRIALS
    assert _verdicts([], full)[("per_test_level", "fraction_sign_test")] == "pass"


def test_the_overall_status_keeps_every_failed_criterion_visible() -> None:
    """#241 used to read only per_test_level; a failed family error must surface."""
    criteria = [
        {"criterion": "per_test_level", "procedure": "p", "verdict": "pass"},
        {"criterion": "family_error", "procedure": "p", "verdict": "fail"},
        {"criterion": "per_test_level", "procedure": "q", "verdict": "pass"},
        {"criterion": "power_reported", "procedure": "q", "verdict": "insufficient_budget"},
        {"criterion": "per_test_level", "procedure": "r", "verdict": "pass"},
        {"criterion": "selection_inflation_measured", "procedure": "r", "verdict": "reported"},
    ]
    status = calibration_status(criteria)
    assert status["p"]["overall"] == "fail"
    assert status["q"]["overall"] == "not_certified"
    assert status["r"]["overall"] == "pass"
    assert status["p"]["criteria"] == {"per_test_level": "pass", "family_error": "fail"}


def test_a_quick_end_to_end_run_certifies_nothing() -> None:
    """The reviewed command's shape, end to end: 5 null and 5 power replicates, seed 1."""
    procedures = {"fraction_sign_test": fraction_sign_test}
    result = run_calibration(_corpus(studies=30), procedures, 5, 5, seed=1)
    assert result.completed["complete_null_study_flip"] == 5
    status = calibration_status(result.criteria)
    assert all(entry["overall"] != "pass" for entry in status.values()), status
    assert "pass" not in {str(row["verdict"]) for row in result.criteria}


# --- #360: the flip makes each exact sign test see its own null -------------


def test_the_exact_size_of_the_sign_test_is_known_in_closed_form() -> None:
    """Six units reject only at 0/6 and 6/6; 710 is the McNemar's discordant count."""
    assert exact_sign_test_size(6) == pytest.approx(2 / 64)
    assert exact_sign_test_size(5) == 0.0
    assert exact_sign_test_size(710) == pytest.approx(0.0466189, abs=1e-6)
    assert exact_sign_test_size(0) is None


def test_the_study_flip_reproduces_the_exact_size_it_cannot_avoid() -> None:
    """Every scored count is fixed under the flip, so the rate converges on the size."""
    profiles = _corpus(studies=40)
    scenario = Scenario(
        "complete_null_study_flip",
        "no trait differs",
        lambda corpus, rng: Resample(flip_directions(corpus, rng, unit="study")),
    )
    run = run_replicates(
        profiles, scenario, 400, seed=SEED, procedures={"fraction_sign_test": fraction_sign_test}
    )
    rows = [r for r in summarize(run.runs, scenario, None) if r["label"] != "FAMILY"]
    for row in rows:
        size = row["exact_size_at_median_scored"]
        assert size is not None
        assert float(str(row["rejection_rate_p"])) == pytest.approx(float(str(size)), abs=0.04)
    scored = {o.scored for r in run.runs for o in r}
    unflipped = {o.scored for o in fraction_sign_test(profiles)}
    assert scored == unflipped


def test_per_test_level_is_a_sanity_check_for_the_published_procedures() -> None:
    """It is reported with its analytic size and kept out of the overall status."""
    ledger = _null_ledger(round(0.04 * DECLARED_TRIALS), 1, DECLARED_TRIALS)
    rows = [
        {
            "scenario": "complete_null_study_flip",
            "procedure": "fraction_sign_test",
            "label": trait,
            "role": "null",
            "exact_size_at_median_scored": 0.048,
        }
        for trait in TRAITS
    ]
    criteria = check_criteria(rows, ledger)
    per_test = next(
        r
        for r in criteria
        if r["criterion"] == "per_test_level" and r["procedure"] == "fraction_sign_test"
    )
    assert per_test["role"] == "sanity_check"
    assert per_test["analytic_exact_size"] == pytest.approx(0.048)
    assert "sanity check only" in str(per_test["detail"])
    verdicts = calibration_status(criteria)["fraction_sign_test"]["criteria"]
    assert isinstance(verdicts, dict)
    assert verdicts["per_test_level"] == "pass (sanity check)"
    # The negative control's per-test row stays judged: selection breaks the construction.
    screened = next(
        r
        for r in criteria
        if r["criterion"] == "per_test_level" and r["procedure"] == "screened_fraction_sign_test"
    )
    assert screened["role"] == "judged"


def test_a_failing_sanity_check_still_fails_the_procedure() -> None:
    criteria = [
        {
            "criterion": "per_test_level",
            "procedure": "p",
            "verdict": "fail",
            "role": "sanity_check",
        },
        {"criterion": "family_error", "procedure": "p", "verdict": "pass", "role": "judged"},
        {
            "criterion": "per_test_level",
            "procedure": "q",
            "verdict": "pass",
            "role": "sanity_check",
        },
        {"criterion": "family_error", "procedure": "q", "verdict": "pass", "role": "judged"},
    ]
    status = calibration_status(criteria)
    assert status["p"]["overall"] == "fail"
    assert status["q"]["overall"] == "pass"


# --- #361: the operating characteristic of the rule actually used -----------


def test_the_declared_budget_meets_the_requirement_under_the_implemented_rules() -> None:
    """Exact enumeration, both rules: 2400 trials missed 0.95 (0.925); 2800 meets it."""
    declared = load_protocol()["criterion_calibration"]
    required = declared["required_operating_characteristics"]
    pinned = declared["declared_operating_characteristics_at_2800_trials"]
    assert DECLARED_TRIALS == 2800
    for bound in (_pooled_bound, _family_bound):
        passes = criterion_pass_probability(0.05, DECLARED_TRIALS, 0.065, bound=bound)
        assert passes >= float(required["pass_probability_at_0.05"])
        assert passes == pytest.approx(float(pinned["pass_probability_at_0.05"]), abs=1e-4)
        at_007 = 1.0 - criterion_pass_probability(0.07, DECLARED_TRIALS, 0.065, bound=bound)
        assert at_007 >= float(required["fail_probability_at_0.07"])
        assert criterion_pass_probability(0.05, 2400, 0.065, bound=bound) == pytest.approx(
            0.9247, abs=1e-4
        )


def test_the_operating_characteristic_is_the_rules_own_threshold() -> None:
    """The enumerated threshold is exactly where the implemented verdict flips."""
    passing = largest_passing_count(DECLARED_TRIALS, 0.065)
    at = _null_ledger(passing, 1, DECLARED_TRIALS, procedure="guild_mcnemar")
    over = _null_ledger(passing + 1, 1, DECLARED_TRIALS, procedure="guild_mcnemar")
    assert _verdicts([], at)[("per_test_level", "guild_mcnemar")] == "pass"
    assert _verdicts([], over)[("per_test_level", "guild_mcnemar")] == "fail"


# --- #362: pin the declared negative control and the surviving mutants --------

#: Studies shifting up out of ten, per Table 1 trait, in TRAITS order.
SCREEN_UPS = (10, 0, 8, 6, 5, 3)


def _screen_corpus() -> list[Profile]:
    """Ten studies; trait i shifts up in SCREEN_UPS[i] of them and down in the rest.

    Both sides carry the same total, so each share is exactly SCREEN_UPS[i] / 10.
    """
    profiles: list[Profile] = []
    for study in range(10):
        inc: dict[str, int] = {}
        dec: dict[str, int] = {}
        for trait, ups in zip(TRAITS, SCREEN_UPS, strict=True):
            up = study < ups
            inc[trait], dec[trait] = (6, 4) if up else (4, 6)
        for direction, counts in (("inc", inc), ("dec", dec)):
            profiles.append(
                Profile(f"s{study}", f"s{study}-1", direction, "C", "Feces", counts, 1000)
            )
    return profiles


def test_the_screen_keeps_the_most_extreme_traits_and_corrects_over_them_only() -> None:
    """Known selection and known q-values: a reversed sort or a raw p would both fail."""
    outcomes = {o.label: o for o in screened_fraction_sign_test(_screen_corpus())}
    selected = {label for label, o in outcomes.items() if o.selected}
    assert selected == {TRAITS[0], TRAITS[1], TRAITS[2]}
    extreme = 2 / 1024
    eight_of_ten = 0.109375
    assert outcomes[TRAITS[0]].p_value == pytest.approx(extreme)
    assert outcomes[TRAITS[2]].p_value == pytest.approx(eight_of_ten)
    # BH over the three selected p-values, not the raw p and not over all six.
    assert outcomes[TRAITS[0]].q_value == pytest.approx(extreme * 3 / 2)
    assert outcomes[TRAITS[1]].q_value == pytest.approx(extreme * 3 / 2)
    assert outcomes[TRAITS[2]].q_value == pytest.approx(eight_of_ten)
    for trait in TRAITS[3:]:
        assert not outcomes[trait].selected
        assert outcomes[trait].p_value is None
        assert outcomes[trait].q_value is None


def test_both_complete_nulls_are_pooled_into_one_verdict() -> None:
    """The ledger alternates the two null scenarios; every replicate must be counted."""
    ledger = _null_ledger(round(0.04 * DECLARED_TRIALS), 1, DECLARED_TRIALS)
    scenarios = {str(r["scenario"]) for r in ledger}
    assert scenarios == {"complete_null_study_flip", "complete_null_experiment_flip"}
    row = next(
        r
        for r in check_criteria([], ledger)
        if r["criterion"] == "per_test_level" and r["procedure"] == "fraction_sign_test"
    )
    assert row["trials"] == DECLARED_TRIALS
    assert row["clusters"] == DECLARED_TRIALS


def test_the_clustered_variance_carries_its_small_sample_factor() -> None:
    """Three clusters of two tests: residual 2, so variance = 3/2 * 2 / 36."""
    pooled = pooled_rate([(1, 2), (0, 2), (2, 2)])
    assert pooled is not None
    assert pooled.rate == 0.5
    assert pooled.standard_error == pytest.approx((3 / 2 * 2 / 36) ** 0.5)


def test_a_p_or_q_exactly_at_the_threshold_counts_as_a_rejection() -> None:
    run = [
        [
            Outcome("fraction_sign_test", "a", ALPHA, 0.05),
            Outcome("fraction_sign_test", "b", 0.0500001, 0.0500001),
        ]
    ]
    (row,) = replicate_counts(run, _scenario("complete_null_study_flip"), None)
    assert row["null_tests_rejected_p"] == 1
    assert row["family_rejected_q"] == 1
    by_label = {r["label"]: r for r in summarize(run, _scenario("complete_null_study_flip"), None)}
    assert (by_label["a"]["rejections_p"], by_label["a"]["rejections_q"]) == (1, 1)
    assert (by_label["b"]["rejections_p"], by_label["b"]["rejections_q"]) == (0, 0)
    assert by_label["FAMILY"]["rejections_q"] == 1
