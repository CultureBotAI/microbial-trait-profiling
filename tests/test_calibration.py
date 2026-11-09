"""Calibration scenarios, resampling and acceptance criteria (#182, #293)."""

from __future__ import annotations

import json
from collections.abc import Sequence
from dataclasses import replace
from itertools import combinations
from pathlib import Path

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
from bsdbng.trait_guilds import OXYGEN_GUILDS, TABLE1_TRAITS, Profile
from mtp.calibration import guild_mcnemar_with
from mtp.guild_maps import compute_correlation
from mtp.trait_guilds import drop_zero_total_pairs

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


def test_guild_calibration_uses_production_whole_pair_eligibility() -> None:
    """A zero-total side must not dilute its study's valid directional contrast."""
    profiles: list[Profile] = []
    for study in range(3):
        for direction in ("inc", "dec"):
            aerobe = 6 if (direction == "inc") == (study > 0) else 4
            profiles.append(
                Profile(
                    study_id=f"s{study}",
                    experiment_id=f"s{study}-valid",
                    inc_dec=direction,
                    condition="C",
                    body_site="Feces",
                    guild_counts={"aerobe": aerobe, "anaerobe": 10 - aerobe},
                    total=10,
                )
            )
    profiles.extend(
        [
            replace(
                profiles[0], experiment_id="s0-invalid", guild_counts={"aerobe": 10, "anaerobe": 0}
            ),
            replace(
                profiles[1],
                experiment_id="s0-invalid",
                guild_counts={"aerobe": 0, "anaerobe": 0},
                total=0,
            ),
        ]
    )
    exclusion = drop_zero_total_pairs(profiles)
    production = compute_correlation(exclusion.profiles)
    raw = compute_correlation(profiles)
    assert raw.discordant != production.discordant
    outcome = guild_mcnemar_with()(profiles)[0]
    assert outcome.p_value == production.mcnemar_p
    assert outcome.scored == sum(production.discordant)
    for unit in ("study", "experiment"):
        flipped = flip_directions(profiles, replicate_rng(SEED, "eligibility", 3), unit=unit)
        retained = drop_zero_total_pairs(flipped)
        assert retained.zero_total_experiments == exclusion.zero_total_experiments == 1
        assert {(p.study_id, p.experiment_id) for p in retained.profiles} == {
            (p.study_id, p.experiment_id) for p in exclusion.profiles
        }
        assert (
            guild_mcnemar_with()(flipped)[0].p_value
            == compute_correlation(retained.profiles).mcnemar_p
        )


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


def test_the_tracked_trait_groups_are_pairwise_disjoint() -> None:
    """The injection headroom is the untracked complement only if no term is in two groups.

    Calibration loads ``{**TABLE1_TRAITS, **OXYGEN_GUILDS}`` and takes the headroom as
    ``total - sum(guild_counts)``. Before #1237, METPO:1000608 was both a Table 1 trait
    and an aerobe-guild member, so its counts were subtracted twice (#1345).
    """
    groups = {**TABLE1_TRAITS, **OXYGEN_GUILDS}
    assert len(groups) == len(TABLE1_TRAITS) + len(OXYGEN_GUILDS)
    for (left, left_terms), (right, right_terms) in combinations(groups.items(), 2):
        assert set(left_terms).isdisjoint(right_terms), (left, right)


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


def test_calibration_summary_marks_p_value_columns_as_uncorrected_diagnostics(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """#1410, #1483, #1485, #1488: raw-p columns and status JSON carry nominal diagnostic labels."""
    import csv
    import re

    import mtp.calibration as cal

    # #1746: p <= 0.05 < q distinguishes nominal per-test rejections_p (1) from
    # corrected FAMILY rejections_p / rejections_q (0).
    runs = [[Outcome("p", "a", 0.02, 0.12, scored=40)]]
    rows = summarize(runs, _scenario(), None)
    row = next(r for r in rows if r["label"] == "a")
    family = next(r for r in rows if r["label"] == "FAMILY")
    assert "rejections_p" in row
    assert "rejection_rate_p" in row
    assert "rejection_rate_p_ci_low" in row
    assert "rejection_rate_p_ci_high" in row
    assert "smallest_observed_p" in row
    assert "attainable_p_below_alpha" in row
    assert "attainable_p_above_alpha" in row
    assert "exact_size_at_median_scored" in row
    assert row["rejections_p"] == 1
    assert row["rejection_rate_p"] == pytest.approx(1.0)
    assert row["rejections_q"] == 0
    assert row["rejection_rate_q"] == pytest.approx(0.0)
    assert row["p_value_diagnostics_role"] == "calibration_diagnostic"
    assert row["p_value_diagnostics_correction"] == "none (nominal)"
    assert "p_value_diagnostics_role" not in family
    assert "p_value_diagnostics_correction" not in family
    assert family["rejections_p"] == 0
    assert family["rejections_q"] == 0
    assert family["rejection_rate_p"] == pytest.approx(0.0)
    assert family["rejection_rate_q"] == pytest.approx(0.0)
    # The FAMILY interval counts replicates with any q <= 0.05, not nominal p (#1794).
    assert (family["rejection_rate_p_ci_low"], family["rejection_rate_p_ci_high"]) == (
        pytest.approx(cal._interval(0, 1))
    )
    assert family["rejections_p"] == family["rejections_q"]
    assert family["rejection_rate_p"] == family["rejection_rate_q"]

    (rep_row,) = replicate_counts(runs, _scenario(), None)
    assert rep_row["null_tests_rejected_p"] == 1
    assert rep_row["family_rejected_q"] == 0
    assert rep_row["p_value_diagnostics_role"] == "calibration_diagnostic"
    assert rep_row["p_value_diagnostics_correction"] == "none (nominal)"

    mixed_runs = [
        [Outcome("p", "a", 0.02, 0.12, scored=40)],
        [Outcome("p", "a", 0.02, 0.04, scored=40)],
    ]
    mixed_rows = summarize(mixed_runs, _scenario(), None)
    mixed_row = next(r for r in mixed_rows if r["label"] == "a")
    mixed_family = next(r for r in mixed_rows if r["label"] == "FAMILY")
    assert mixed_row["rejections_p"] == 2
    assert mixed_row["rejections_q"] == 1
    assert mixed_family["rejections_p"] == 1
    assert mixed_family["rejections_q"] == 1
    assert mixed_family["rejection_rate_p"] == pytest.approx(0.5)
    assert mixed_family["rejection_rate_q"] == pytest.approx(0.5)
    assert (
        mixed_family["rejection_rate_p_ci_low"],
        mixed_family["rejection_rate_p_ci_high"],
    ) == pytest.approx(cal._interval(1, 2))

    note = str(cal.P_VALUE_DIAGNOSTICS_NOTE)
    assert note == (
        "On per-test rows, rejections_p, rejection_rate_p, rejection_rate_p_ci_low, "
        "rejection_rate_p_ci_high, smallest_observed_p, attainable_p_below_alpha, "
        "attainable_p_above_alpha, exact_size_at_median_scored, and replicate-level "
        "null_tests_rejected_p (at nominal p <= 0.05) are uncorrected nominal "
        "p-value calibration diagnostics, not inferential results; on FAMILY rows, "
        "rejections_p, rejection_rate_p, rejection_rate_p_ci_low, and "
        "rejection_rate_p_ci_high count replicates with any q <= 0.05 "
        "(a corrected family-wise rate, matching rejections_q)"
    )
    assert "ALPHA" not in note
    assert "Q_THRESHOLD" not in note
    for composition_only_token in (
        "null_p_value",
        "observed_p",
        "null_type_i_alpha_05",
        "type_i_error_metrics",
    ):
        assert re.search(rf"\b{re.escape(composition_only_token)}\b", note) is None

    monkeypatch.setattr(
        cal,
        "load_profiles",
        lambda *_a, **_k: (_corpus(studies=10), Outcome("p", "a", 0.02, 0.04)),
    )
    monkeypatch.setattr(cal, "read_bugsigdb_vintage", lambda *_a: "2026-04-24")
    monkeypatch.setattr(
        cal, "write_stats_report", lambda *_a: (tmp_path / "s.md", tmp_path / "s.tsv")
    )
    monkeypatch.setattr(cal, "file_sha256", lambda *_a: "0" * 64)
    monkeypatch.setattr(cal, "SHARED_FLIP_REPLICATES", 1)
    monkeypatch.setattr(cal, "SCENARIOS", (cal.NULL_SCENARIOS[0],))
    monkeypatch.setattr(cal, "procedures_for", lambda: {"fraction_sign_test": fraction_sign_test})

    assert (
        cal.main(
            [
                "--output-dir",
                str(tmp_path),
                "--replicates",
                "1",
                "--power-replicates",
                "1",
                "--seed",
                "1",
            ]
        )
        == 0
    )
    status = json.loads((tmp_path / "calibration_status.json").read_text(encoding="utf-8"))
    assert status["p_value_diagnostics"] == cal.P_VALUE_DIAGNOSTICS_NOTE
    assert status["not_calibrated_here"] == list(cal.NOT_CALIBRATED_HERE)

    with (tmp_path / "calibration_results.tsv").open(encoding="utf-8") as handle:
        written_results = list(csv.DictReader(handle, delimiter="\t"))
    non_family = [r for r in written_results if r["label"] != "FAMILY"]
    family_rows = [r for r in written_results if r["label"] == "FAMILY"]
    assert non_family
    assert family_rows
    assert {r["p_value_diagnostics_role"] for r in non_family} == {"calibration_diagnostic"}
    assert {r["p_value_diagnostics_correction"] for r in non_family} == {"none (nominal)"}
    assert {r["p_value_diagnostics_role"] for r in family_rows} == {""}
    assert {r["p_value_diagnostics_correction"] for r in family_rows} == {""}

    with (tmp_path / "calibration_replicates.tsv").open(encoding="utf-8") as handle:
        written_replicates = list(csv.DictReader(handle, delimiter="\t"))
    assert written_replicates
    assert {r["p_value_diagnostics_role"] for r in written_replicates} == {"calibration_diagnostic"}
    assert {r["p_value_diagnostics_correction"] for r in written_replicates} == {"none (nominal)"}


def test_a_nonfinite_mcnemar_p_is_not_counted_as_a_tested_replicate(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A non-finite McNemar p must not be passed or counted as a tested replicate (#1808)."""
    from types import SimpleNamespace

    import mtp.calibration as cal
    import mtp.guild_maps as gm

    profiles: list[Profile] = []
    for study in range(4):
        for direction in ("inc", "dec"):
            aerobe = 6 if (direction == "inc") == (study > 0) else 4
            profiles.append(
                Profile(
                    study_id=f"s{study}",
                    experiment_id=f"s{study}-1",
                    inc_dec=direction,
                    condition="C",
                    body_site="Feces",
                    guild_counts={"aerobe": aerobe, "anaerobe": 10 - aerobe},
                    total=10,
                )
            )
    baseline = compute_correlation(drop_zero_total_pairs(profiles).profiles)
    assert baseline.status == "tested"
    assert baseline.discordant == (3, 1)
    scenario = _scenario("complete_null_study_flip")

    for bad_p in (float("nan"), float("inf"), float("-inf")):
        with monkeypatch.context() as patch:
            patch.setattr(gm, "binomtest", lambda *_a, _p=bad_p, **_k: SimpleNamespace(pvalue=_p))
            (outcome,) = guild_mcnemar_with()(profiles)
        assert outcome.p_value is None
        assert outcome.q_value is None
        assert outcome.scored == sum(baseline.discordant)

        rows = summarize([[outcome]], scenario, None)
        summary = next(r for r in rows if r["label"] == "aerobe_vs_anaerobe")
        family = next(r for r in rows if r["label"] == "FAMILY")
        assert summary["replicates_tested"] == 0
        assert summary["rejections_p"] == 0
        assert summary["rejection_rate_p"] is None
        assert summary["rejections_q"] == 0
        assert summary["rejection_rate_q"] is None
        assert summary["smallest_observed_p"] is None
        assert family["rejections_p"] == 0
        assert family["rejections_q"] == 0

        (replicate,) = replicate_counts([[outcome]], scenario, None)
        assert replicate["tests_evaluated"] == 0
        assert replicate["null_tests_evaluated"] == 0
        assert replicate["null_tests_rejected_p"] == 0
        assert replicate["family_rejected_q"] == 0

    for overridden in (
        replace(baseline, status="tested", mcnemar_p=float("nan")),
        replace(baseline, status="tested", mcnemar_p=float("inf")),
        replace(baseline, status="tested", mcnemar_p=None),
        replace(baseline, status="not_tested", mcnemar_p=0.25),
    ):
        with monkeypatch.context() as patch:
            patch.setattr(cal, "compute_correlation", lambda *_a, _c=overridden, **_k: _c)
            (outcome,) = guild_mcnemar_with()(profiles)
        assert outcome.p_value is None
        assert outcome.q_value is None

    for bad_p in (float("nan"), float("inf"), float("-inf")):
        direct = Outcome("guild_mcnemar", "aerobe_vs_anaerobe", bad_p, bad_p, scored=4)
        rows = summarize([[direct]], scenario, None)
        summary = next(r for r in rows if r["label"] == "aerobe_vs_anaerobe")
        family = next(r for r in rows if r["label"] == "FAMILY")
        assert summary["replicates_tested"] == 0
        assert summary["rejections_p"] == 0
        assert summary["rejections_q"] == 0
        assert summary["smallest_observed_p"] is None
        assert family["rejections_p"] == 0
        assert family["rejections_q"] == 0

        (replicate,) = replicate_counts([[direct]], scenario, None)
        assert replicate["tests_evaluated"] == 0
        assert replicate["null_tests_evaluated"] == 0
        assert replicate["null_tests_rejected_p"] == 0
        assert replicate["family_rejected_q"] == 0

        bad_q = Outcome("guild_mcnemar", "aerobe_vs_anaerobe", 0.01, bad_p, scored=4)
        rows_q = summarize([[bad_q]], scenario, None)
        summary_q = next(r for r in rows_q if r["label"] == "aerobe_vs_anaerobe")
        family_q = next(r for r in rows_q if r["label"] == "FAMILY")
        assert summary_q["replicates_tested"] == 1
        assert summary_q["rejections_p"] == 1
        assert summary_q["rejections_q"] == 0
        assert family_q["rejections_p"] == 0
        assert family_q["rejections_q"] == 0

        (replicate_q,) = replicate_counts([[bad_q]], scenario, None)
        assert replicate_q["tests_evaluated"] == 1
        assert replicate_q["null_tests_evaluated"] == 1
        assert replicate_q["null_tests_rejected_p"] == 1
        assert replicate_q["family_rejected_q"] == 0

        # #1883: A FAMILY q rejection also requires a finite p_value.
        bad_p_good_q = Outcome("guild_mcnemar", "aerobe_vs_anaerobe", bad_p, 0.01, scored=4)
        rows_pq = summarize([[bad_p_good_q]], scenario, None)
        summary_pq = next(r for r in rows_pq if r["label"] == "aerobe_vs_anaerobe")
        family_pq = next(r for r in rows_pq if r["label"] == "FAMILY")
        assert summary_pq["replicates_tested"] == 0
        assert summary_pq["rejections_p"] == 0
        assert summary_pq["rejections_q"] == 0
        assert summary_pq["rejection_rate_p"] is None
        assert summary_pq["rejection_rate_q"] is None
        assert family_pq["rejections_p"] == 0
        assert family_pq["rejections_q"] == 0
        assert family_pq["rejection_rate_p"] == 0.0
        assert family_pq["rejection_rate_q"] == 0.0

        (replicate_pq,) = replicate_counts([[bad_p_good_q]], scenario, None)
        assert replicate_pq["tests_evaluated"] == 0
        assert replicate_pq["null_tests_evaluated"] == 0
        assert replicate_pq["null_tests_rejected_p"] == 0
        assert replicate_pq["family_rejected_q"] == 0


def test_calibration_protocol_and_docs_declare_ora_urn_scope_amendment() -> None:
    """#568: protocol, calibration module, and methods docs distinguish the three ORA urns."""
    import mtp.calibration as cal

    repo_root = PROTOCOL_PATH.resolve().parents[1]
    declared = load_protocol()

    assert declared["protocol"] == "calibration-v3"
    assert declared["declared_by"] == (
        "calibration-v1 was drafted by a coding agent (Claude Code) before any "
        "calibration result was computed. calibration-v2 was amended after review "
        "(#293) with the v1 results in hand. calibration-v3 repairs the guild "
        "population mismatch identified in #730, and its scope and not_covered "
        "declarations were amended in #568 to record the three ORA background/urn "
        "definitions and the separate arm-calibration harnesses (#550, #551). The "
        "amendments array records what changed, why, and what was already known "
        "when each change was made. Nothing here should be read as predeclared "
        "except where an amendment says it is unchanged from v1."
    )

    scope = str(declared["scope"])
    assert (
        "their adapters belong to their owners (#182 assigns them), and the "
        "assertion-level urn of #131 is a separate model question."
    ) not in scope
    assert scope == (
        "procedures owned by this repository's directional-consistency framework "
        "(src/mtp/calibration.py). As amended in #568, the enrichment arm's "
        "production hypergeometric procedure draws distinct reported taxa against "
        "the all-host taxa_trait_matrix.tsv background (M = 5,829 BugSigDB-reported "
        "taxa, including 645 non-human-only and 862 non-paired-human taxa; "
        "src/mtp/three_test_table/ora.py, calibrated separately in "
        "src/mtp/ora_calibration/ under #551), whereas the historical #131 "
        "assertion-count urn (106,993,585 assertions across 4,517 testable profiles) "
        "is retained only in the Table 1 native influence ORA reruns "
        "(src/mtp/native_taxon_influence.py, methods/table1_native_influence.md) and "
        "the harmonized baseline in methods/three_arm_statistics.md defines a "
        "proposed paired-human assertion-weighted background; the composition arm's "
        "native transformed procedure is calibrated separately in "
        "src/mtp/composition_calibration/ (#550)."
    )

    amendments = list(declared["amendments"])
    assert len(amendments) == 8
    assert amendments[0] == {
        "date": "2026-11-08",
        "protocol": "calibration-v3",
        "change": (
            "Amended scope and not_covered[0] to distinguish the three ORA "
            "background/urn definitions and record the separate arm-calibration "
            "harnesses (#550, #551) without changing the directional acceptance "
            "rules, eligibility policy, or replicate budget."
        ),
        "reason": (
            "Issue #568 showed that scope and not_covered[0] still carried the "
            "initial #182 owner-adapter and #131 assertion-level urn wording after "
            "production ORA switched to the distinct-taxon urn over "
            "taxa_trait_matrix.tsv (M = 5,829 BugSigDB-reported taxa across all "
            "hosts, including 645 non-human-only and 862 non-paired-human taxa; "
            "src/mtp/three_test_table/ora.py, calibrated in src/mtp/ora_calibration/ "
            "under #551), while the historical #131 assertion-count urn (106,993,585 "
            "assertions across 4,517 testable profiles) is retained only in the "
            "Table 1 native influence ORA reruns (src/mtp/native_taxon_influence.py, "
            "methods/table1_native_influence.md), the harmonized baseline in "
            "methods/three_arm_statistics.md defines a proposed paired-human "
            "assertion-weighted background, and the composition arm is calibrated in "
            "src/mtp/composition_calibration/ (#550)."
        ),
        "known_when_amended": (
            "The completed calibration-v3 full-budget run (58b2ce84), the production "
            "distinct-taxon ORA implementation and calibration harness (#551, #558), "
            "the paired composition calibration harness (#550), and the Table 1 "
            "native influence ORA reruns (#190)."
        ),
        "unchanged": (
            "Every directional procedure, eligibility policy, seed 182, null and "
            "injection design, deterministic flip-assignment stream, replicate "
            "budget, acceptance criterion, and tolerance."
        ),
    }
    assert [entry["date"] for entry in amendments[1:]] == [
        "2026-10-19",
        "2026-09-16",
        "2026-10-02",
        "2026-10-02",
        "2026-10-02",
        "2026-10-08",
        "2026-10-09",
    ]

    not_covered = list(declared["not_covered"])
    assert (
        "the enrichment arm's hypergeometric urn (#131) and the composition arm's "
        "native transform (owners' adapters)"
    ) not in not_covered
    assert not_covered[0] == (
        "the enrichment arm's hypergeometric procedures — production distinct-taxon "
        "urn over taxa_trait_matrix.tsv (M = 5,829 all-host taxa, including 645 "
        "non-human-only and 862 non-paired-human taxa; calibrated separately in "
        "src/mtp/ora_calibration/, #551), historical #131 assertion-count urn "
        "retained in the Table 1 native influence reruns (106,993,585 assertions "
        "across 4,517 testable profiles; methods/table1_native_influence.md), and "
        "proposed harmonized-baseline paired-human assertion-weighted background "
        "(methods/three_arm_statistics.md) — and the composition arm's native "
        "transformed procedure (calibrated separately in "
        "src/mtp/composition_calibration/, #550)"
    )

    doc = " ".join((cal.__doc__ or "").split())
    assert "need adapters from their owners" not in doc
    assert (
        "Only the directional-consistency procedures in this module are evaluated "
        "here. As recorded in the v3.1 (#568) amendment to "
        "methods/calibration_protocol.json, the enrichment arm's production "
        "distinct-taxon urn over taxa_trait_matrix.tsv (M = 5,829 all-host taxa, "
        "calibrated in src/mtp/ora_calibration/, #551), the historical #131 "
        "assertion-count urn retained only in the Table 1 native influence ORA "
        "reruns (106,993,585 assertions across 4,517 testable profiles, "
        "src/mtp/native_taxon_influence.py), the proposed harmonized-baseline "
        "paired-human assertion-weighted ORA background "
        "(methods/three_arm_statistics.md), and the composition arm's native "
        "transformed test (calibrated in src/mtp/composition_calibration/, #550) "
        "are outside this directional protocol."
    ) in doc

    assert cal.NOT_CALIBRATED_HERE == (
        "the enrichment arm's hypergeometric procedures (production M = 5,829 "
        "all-host distinct-taxon urn in src/mtp/ora_calibration/ #551; historical "
        "#131 assertion-count urn in Table 1 native influence reruns; proposed "
        "harmonized-baseline paired-human assertion-weighted background #568)",
        "the composition arm's native transformed test (calibrated in "
        "src/mtp/composition_calibration/ #550)",
        "the aggregated profile-support rule (#177)",
    )
    bsdbng_shim_source = (repo_root / "src" / "bsdbng" / "calibration.py").read_text(
        encoding="utf-8"
    )
    assert '"NOT_CALIBRATED_HERE",' in bsdbng_shim_source

    calibration_md = " ".join(
        (repo_root / "methods" / "calibration.md").read_text(encoding="utf-8").split()
    )
    assert "its adapter belongs to its owner" not in calibration_md
    assert (
        "The enrichment arm's hypergeometric procedures across all three documented "
        "ORA background definitions (`methods/three_arm_statistics.md`, #568): "
        "production ORA (`src/mtp/three_test_table/ora.py`) draws distinct reported "
        "taxa against the all-host `taxa_trait_matrix.tsv` background (`M = 5,829` "
        "BugSigDB-reported taxa, including 645 non-human-only and 862 "
        "non-paired-human taxa) rather than trait-annotation assertions (#131), "
        "where nested taxa inherit the same annotations and are not independent "
        "draws, and is calibrated separately at the production thresholds in "
        "`src/mtp/ora_calibration/` (#551); the Table 1 native influence ORA reruns "
        "(`src/mtp/native_taxon_influence.py`, `methods/table1_native_influence.md`) "
        "retain the historical #131 assertion-count urn (4,517 testable profiles "
        "against 106,993,585 assertions across 556,516 enrichment tests, with "
        "depletion fold <= 0.67); and the harmonized baseline in "
        "`methods/three_arm_statistics.md` defines a proposed paired-human "
        "assertion-weighted background."
    ) in calibration_md
    assert (
        "The composition arm's native transformed test, which is calibrated "
        "separately against the shared study-flip null in "
        "`src/mtp/composition_calibration/` (`methods/calibration_composition.md`, "
        "#550). [`calibration_protocol.json`](calibration_protocol.json) records "
        "this scope amendment explicitly (`v3.1`, #568) rather than silently "
        "rewriting its earlier #131 assertion-level and #182 owner-adapter "
        "declaration."
    ) in calibration_md

    three_arm_md = " ".join(
        (repo_root / "methods" / "three_arm_statistics.md").read_text(encoding="utf-8").split()
    )
    assert (
        "[`calibration_protocol.json`](calibration_protocol.json) (`v3.1`, #568) "
        "and [`calibration.md`](calibration.md) explicitly record how these three "
        "ORA background definitions relate to the directional calibration scope and "
        "the separate ORA calibration harness (`src/mtp/ora_calibration/`, #551)."
    ) in three_arm_md

    native_influence_md = " ".join(
        (repo_root / "methods" / "table1_native_influence.md").read_text(encoding="utf-8").split()
    )
    assert (
        "[`calibration_protocol.json`](calibration_protocol.json) (`v3.1`, #568) "
        "and [`calibration.md`](calibration.md) record this distinction between the "
        "historical #131 assertion-count urn retained here, the production "
        "`M = 5,829` all-host distinct-taxon urn calibrated in "
        "`src/mtp/ora_calibration/` (#551), and the proposed harmonized-baseline "
        "paired-human assertion-weighted background."
    ) in native_influence_md
