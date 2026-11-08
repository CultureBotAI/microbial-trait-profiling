"""Nominal narrative labels describe displayed endpoints, not a new FDR family (#1859)."""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import pandas as pd
import pytest

from mtp.cellulose_decr import analysis


def narrative_inputs() -> tuple[pd.DataFrame, ...]:
    """Seven primary rows with independently specified secondary endpoint identities."""
    compartments = pd.DataFrame(
        [
            {
                "compartment": analysis.FECES,
                "contrast_population": analysis.ALL_CONTRASTS,
                "short_condition": "CRC",
                "condition": "Colorectal cancer",
                "paired_exp_down_dec_gt_inc": 4,
                "paired_exp_up_inc_gt_dec": 2,
                "paired_exp_tied": 1,
                "share_down_untied": 2 / 3,
                "binom_p_two_sided": 0.3,
                "wilcoxon_p_two_sided": 1.0,
                "study_down": 2,
                "study_up": 1,
                "study_tied": 1,
                "study_share_down_untied": 2 / 3,
                "study_binom_p_two_sided": 0.7,
                "study_wilcoxon_p_two_sided": math.inf,
                "mean_paired_delta_inc_minus_dec": -0.1,
                "study_mean_delta_inc_minus_dec": -0.2,
                "n_studies_dec": 3,
                "n_studies_inc": 4,
                "n_profiles_dec": 5,
                "n_profiles_inc": 6,
            }
        ]
    )
    scenario_direction = pd.DataFrame(
        [
            {
                "scenario": scenario,
                "contrast_population": analysis.ALL_CONTRASTS,
                "short_condition": "CRC",
                "unit": unit,
                "binom_p_two_sided": value,
            }
            for scenario, unit, value in (
                ("baseline", "experiment", 0.0),
                ("alternative", "experiment", 0.4),
                ("no_pairs", "experiment", math.nan),
                ("nonfinite", "experiment", math.inf),
                ("baseline", "experiment", 0.0),
                ("baseline", "study", 0.2),
                ("alternative", "study", 0.2),
            )
        ]
    )
    direction_sentinels = [
        {
            **scenario_direction.iloc[0].to_dict(),
            "contrast_population": analysis.ORIENTED,
            "binom_p_two_sided": 0.91,
        },
        {
            **scenario_direction.iloc[0].to_dict(),
            "short_condition": "HIV",
            "binom_p_two_sided": 0.92,
        },
    ]
    scenario_direction = pd.concat(
        [scenario_direction, pd.DataFrame(direction_sentinels)], ignore_index=True
    )
    common = {
        "comparison": "PD vs. CRC",
        "intersection": 2,
        "union": 8,
        "jaccard_distance": 0.75,
        "u1_expected_jaccard_distance": 0.9,
        "u5_emp_weighted_expected_jaccard_distance": 0.8,
        "u6_signature_perm_expected_jaccard_distance": 0.7,
        "u6_signature_perm_expected_intersection": 2.5,
        "u1_fold_enrichment": 2.0,
    }
    pairwise = pd.DataFrame(
        [
            {
                **common,
                "endpoint": "pairwise_jaccard_distance",
                "u1_hypergeom_p": 1.0,
                "u5_emp_weighted_p": math.nan,
                "u6_signature_perm_p": 0.0,
            },
            {
                **common,
                "endpoint": "mean_pairwise_jaccard_distance",
                "u4_shared_trait_p": math.nan,
                "u5_emp_weighted_p": math.inf,
                "u6_signature_perm_p": 0.5,
            },
            {
                **common,
                "endpoint": "three_way_core",
                "u4_shared_trait_p": 1.0,
                "u5_emp_weighted_p": 0.1,
                "u6_signature_perm_p": math.nan,
            },
        ]
    )
    interaction_rows: list[dict[str, Any]] = []
    for short, condition, p in (
        ("CRC", "Colorectal cancer", 0.04),
        (analysis.POOLED_STRATIFIED, "Pooled, stratified by condition", 0.2),
    ):
        interaction_rows.append(
            {
                "short_condition": short,
                "condition": condition,
                "contrast_population": analysis.ALL_CONTRASTS,
                "unit": "experiment",
                "compartment_b": analysis.NON_FECAL_GI,
                "studies_in_both_compartments": 1,
                "conditional_or": 2.0,
                "conditional_or_ci95_low": 0.5,
                "conditional_or_ci95_high": 4.0,
                "a_down": 4,
                "a_up": 2,
                "b_down": 2,
                "b_up": 4,
                "fisher_p_two_sided": p,
            }
        )
    interaction = pd.DataFrame(interaction_rows)
    scenario_interaction = pd.DataFrame(
        [
            {**interaction_rows[index], "scenario": scenario, "fisher_p_two_sided": value}
            for index, scenario, value in (
                (0, "baseline", 0.6),
                (0, "alternative", 0.6),
                (0, "no_pairs", math.nan),
                (0, "baseline", 0.6),
                (1, "baseline", 0.0),
                (1, "nonfinite", math.inf),
            )
        ]
    )
    interaction_sentinels = [
        {
            **scenario_interaction.iloc[0].to_dict(),
            "contrast_population": analysis.ORIENTED,
            "fisher_p_two_sided": 0.91,
        },
        {
            **scenario_interaction.iloc[0].to_dict(),
            "short_condition": "HIV",
            "fisher_p_two_sided": 0.92,
        },
        {
            **scenario_interaction.iloc[0].to_dict(),
            "compartment_b": "Oral",
            "fisher_p_two_sided": 0.93,
        },
        {
            **scenario_interaction.iloc[0].to_dict(),
            "unit": "study",
            "fisher_p_two_sided": 0.94,
        },
    ]
    scenario_interaction = pd.concat(
        [scenario_interaction, pd.DataFrame(interaction_sentinels)], ignore_index=True
    )
    return compartments, scenario_direction, pairwise, interaction, scenario_interaction


def test_direction_narratives_label_both_tests_and_count_distinct_definitions() -> None:
    compartments, scenarios, *_ = narrative_inputs()
    try:
        rows = analysis.direction_summary_rows(compartments, scenarios)
    except ValueError as exc:
        raise AssertionError(
            "Unrelated scenario scopes must not enter direction narratives"
        ) from exc
    assert len(rows) == 2
    experiment, study = [row["summary_interpretation"] for row in rows]
    assert "Wilcoxon 2-sided nominal p = 1." in experiment
    assert "nominal sign-test p-values range from 0 to 0.4" in experiment
    assert "Wilcoxon 2-sided nominal p on study means = unavailable." in study
    assert "nominal sign-test p-values range from 0.2 to 0.2" in study
    assert experiment.endswith(
        "Nominal narrative endpoints: 1 Wilcoxon and 4 condition-definition sign tests; "
        "3 finite; 2 unavailable; not a correction family."
    )
    assert study.endswith(
        "Nominal narrative endpoints: 1 Wilcoxon and 2 condition-definition sign tests; "
        "2 finite; 1 unavailable; not a correction family."
    )
    for text in (experiment, study):
        assert "all distinct condition definitions represented by the range" in text
        assert "this p" not in text
        assert "measured cellulose degradation" in text


@pytest.mark.parametrize("empty", [False, True], ids=["unavailable", "no-definitions"])
def test_direction_range_reports_unavailable_without_inventing_endpoints(empty: bool) -> None:
    compartments, scenarios, *_ = narrative_inputs()
    selected = scenarios.loc[scenarios["unit"].eq("experiment")].copy()
    selected["binom_p_two_sided"] = math.nan
    if empty:
        selected = selected.iloc[0:0]
    text = analysis.direction_summary_rows(compartments, selected)[0]["summary_interpretation"]
    expected_range = "unavailable (no condition definitions)" if empty else "unavailable"
    assert f"nominal sign-test p-values range: {expected_range}." in text
    expected_count = 0 if empty else 4
    assert text.endswith(
        "Nominal narrative endpoints: 1 Wilcoxon and "
        f"{expected_count} condition-definition sign tests; "
        f"1 finite; {expected_count} unavailable; "
        "not a correction family."
    )


def test_pairwise_narrative_counts_u6_once_and_preserves_zero_and_one() -> None:
    _, _, pairwise, *_ = narrative_inputs()
    row = analysis.convergence_summary_rows(pairwise)[0]
    text = row["summary_interpretation"]
    assert "U1 exact nominal p = 1; U5 nominal p = unavailable; " in text
    assert "U6 (whole-signature permutation) nominal p = 0." in text
    assert text.endswith(
        "Nominal narrative endpoints: 2 finite; 1 unavailable; not a correction family."
    )
    assert row["p_value"] == 0.0


@pytest.mark.parametrize(
    ("index", "expected", "scope"),
    [
        (
            1,
            "U4 unavailable; U5 unavailable; U6 0.5 (nominal p-values)",
            "1 finite; 2 unavailable",
        ),
        (2, "U4 1; U5 0.1; U6 unavailable (nominal p-values)", "2 finite; 1 unavailable"),
    ],
    ids=["mean-distance", "three-way-core"],
)
def test_aggregate_narrative_labels_all_three_null_endpoints(
    index: int, expected: str, scope: str
) -> None:
    _, _, pairwise, *_ = narrative_inputs()
    row = analysis.convergence_summary_rows(pairwise)[index]
    assert row["null_or_reference_expectation"] == expected
    assert row["summary_interpretation"].endswith(
        f"Nominal narrative endpoints: {scope}; not a correction family."
    )


def test_interaction_ranges_name_the_actual_exact_test_and_definition_scope() -> None:
    *_, interaction, scenarios = narrative_inputs()
    try:
        rows = analysis.interaction_summary_rows(interaction, scenarios)
    except ValueError as exc:
        raise AssertionError(
            "Unrelated scenario scopes must not enter interaction narratives"
        ) from exc
    ordinary, stratified = [row["summary_interpretation"] for row in rows]
    assert "nominal Fisher exact p-values range from 0.6 to 0.6" in ordinary
    assert "nominal exact conditional p-values range from 0 to 0" in stratified
    assert ordinary.endswith(
        "Nominal narrative endpoints: 2 finite; 1 unavailable; not a correction family."
    )
    assert stratified.endswith(
        "Nominal narrative endpoints: 1 finite; 1 unavailable; not a correction family."
    )
    for text in (ordinary, stratified):
        assert "all distinct condition definitions represented by the range" in text
        assert "this p" not in text
        assert "not a test of anatomical specificity" in text


def test_no_scenario_column_is_one_supplied_definition_per_selected_row_scope() -> None:
    *_, interaction, _ = narrative_inputs()
    rows = analysis.interaction_summary_rows(interaction, interaction.copy())
    assert len(rows) == 2
    for row in rows:
        assert row["summary_interpretation"].endswith(
            "Nominal narrative endpoints: 1 finite; 0 unavailable; not a correction family."
        )


@pytest.mark.parametrize("kind", ["direction", "interaction"])
def test_duplicate_nan_scenario_is_accepted_and_counted_once(kind: str) -> None:
    compartments, direction, _, interaction, scenarios = narrative_inputs()
    source = direction if kind == "direction" else scenarios
    unavailable = source.loc[source["scenario"].eq("no_pairs")].copy()
    assert len(unavailable) == 1
    duplicates = pd.concat([unavailable, unavailable], ignore_index=True)
    before = duplicates.copy(deep=True)
    try:
        if kind == "direction":
            row = analysis.direction_summary_rows(compartments, duplicates)[0]
        else:
            row = analysis.interaction_summary_rows(interaction, duplicates)[0]
    except ValueError as error:
        raise AssertionError(
            "Two NaN copies of one scenario must be accepted as one unavailable endpoint"
        ) from error
    text = row["summary_interpretation"]
    finite = 1 if kind == "direction" else 0
    assert text.endswith(f"{finite} finite; 1 unavailable; not a correction family.")
    if kind == "direction":
        assert (
            "Nominal narrative endpoints: 1 Wilcoxon and 1 condition-definition sign test; " in text
        )
    assert "p-values range: unavailable." in text
    assert row["p_value"] == (0.3 if kind == "direction" else 0.04)
    pd.testing.assert_frame_equal(before, duplicates)


@pytest.mark.parametrize("kind", ["direction", "interaction"])
def test_conflicting_duplicate_scenario_endpoint_is_refused(kind: str) -> None:
    compartments, direction, _, interaction, scenarios = narrative_inputs()
    message = "Conflicting nominal p-values for one scenario endpoint"
    if kind == "direction":
        direction.loc[4, "binom_p_two_sided"] = 0.9
        with pytest.raises(ValueError, match=message):
            analysis.direction_summary_rows(compartments, direction)
    else:
        scenarios.loc[3, "fisher_p_two_sided"] = 0.9
        with pytest.raises(ValueError, match=message):
            analysis.interaction_summary_rows(interaction, scenarios)


def test_narrative_changes_preserve_primary_probabilities_bh_and_row_order() -> None:
    inputs = narrative_inputs()
    before = [frame.copy(deep=True) for frame in inputs]
    summary = analysis.significant_results_summary_table(*inputs)
    assert list(summary["test_id"]) == [
        "feces_CRC_experiment_direction_all_contrasts",
        "feces_CRC_study_direction_all_contrasts",
        "jaccard_PD_vs_CRC",
        "mean_pairwise_jaccard_distance",
        "three_way_core",
        "fisher_feces_vs_gi_CRC_experiment_all_contrasts",
        "fisher_feces_vs_gi_pooled_stratified_experiment_all_contrasts",
    ]
    np.testing.assert_allclose(
        summary["p_value"], [0.3, 0.7, 0.0, 0.5, math.nan, 0.04, 0.2], equal_nan=True
    )
    # Six finite primary hypotheses; sorted BH q values are 0, .12, .4, .45, .6, .7.
    np.testing.assert_allclose(
        summary["bh_q_within_table"], [0.45, 0.7, 0.0, 0.6, math.nan, 0.12, 0.4], equal_nan=True
    )
    assert set(summary["bh_family_size"]) == {6}
    assert set(summary["bh_unfiltered_summary_rows"]) == {7}
    assert set(summary["bh_selection_rule"]) == {analysis.SUMMARY_BH_SELECTION_RULE}
    for original, actual in zip(before, inputs, strict=True):
        pd.testing.assert_frame_equal(original, actual)
