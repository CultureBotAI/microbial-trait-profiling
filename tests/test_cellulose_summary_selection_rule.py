"""Pin the cellulose summary's complete selection rule and registry anchor (#1585)."""

from __future__ import annotations

import itertools
import json
from pathlib import Path

import pandas as pd

from mtp.cellulose_decr.analysis import SUMMARY_BH_SELECTION_RULE, interaction_summary_rows


def test_summary_selection_rule_states_the_complete_independent_contract() -> None:
    """Changing either units clause or any selected endpoint must change this contract."""
    expected = (
        "Concatenate direction_summary_rows (Feces all_contrasts and "
        "disease_vs_healthy_oriented, experiment and study units), "
        "convergence_summary_rows (U6 pairwise_jaccard_distance rows "
        "jaccard_PD_vs_HIV, jaccard_PD_vs_CRC, jaccard_HIV_vs_CRC, plus "
        "mean_pairwise_jaccard_distance and three_way_core), and "
        "interaction_summary_rows (Feces vs Non-fecal GI (mucosa/lumen) for CRC and "
        "Pooled, stratified by condition; all_contrasts and disease_vs_healthy_oriented, "
        "experiment and study units); apply Benjamini-Hochberg over "
        "the non-missing p_value entries in this selected summary table. Correcting "
        "within this selected summary table does not repair selection from the larger "
        "exploratory tables."
    )
    assert expected == SUMMARY_BH_SELECTION_RULE


def test_interaction_selection_keeps_both_populations_and_both_units() -> None:
    """The real row producer retains the eight stated combinations and rejects controls."""
    records: list[dict[str, object]] = []
    for population, unit, (short_condition, condition) in itertools.product(
        ("all_contrasts", "disease_vs_healthy_oriented"),
        ("experiment", "study"),
        (
            ("CRC", "Colorectal cancer"),
            (
                "Pooled, stratified by condition",
                "Pooled (PD + HIV + CRC), stratified by condition",
            ),
        ),
    ):
        records.append(
            {
                "contrast_population": population,
                "unit": unit,
                "short_condition": short_condition,
                "condition": condition,
                "compartment_b": "Non-fecal GI (mucosa/lumen)",
                "studies_in_both_compartments": 0,
                "conditional_or": 1.0,
                "conditional_or_ci95_low": 0.5,
                "conditional_or_ci95_high": 2.0,
                "a_down": 2,
                "a_up": 2,
                "b_down": 2,
                "b_up": 2,
                "fisher_p_two_sided": 0.5,
            }
        )
    records.extend(
        (
            {**records[0], "short_condition": "HIV", "condition": "HIV infection"},
            {**records[0], "compartment_b": "Oral"},
        )
    )
    interaction = pd.DataFrame.from_records(records)
    rows = interaction_summary_rows(interaction, interaction.copy())

    assert len(rows) == 8
    expected = set(
        itertools.product(
            ("all_contrasts", "disease_vs_healthy_oriented"),
            ("paired experiment", "study (mean delta)"),
            (
                "Colorectal cancer: Feces vs Non-fecal GI (mucosa/lumen)",
                "Pooled (PD + HIV + CRC), stratified by condition: "
                "Feces vs Non-fecal GI (mucosa/lumen)",
            ),
        )
    )
    assert {
        (row["contrast_population"], row["unit_of_observation"], row["scope_or_comparison"])
        for row in rows
    } == expected
    assert {row["result_category"] for row in rows} == {"3_compartment_interaction"}


def test_cellulose_registry_selection_rule_names_its_actual_anchor() -> None:
    """Pin the existing registry text without changing the registry or its central guard."""
    registry_path = Path(__file__).resolve().parents[1] / "methods/multiple_testing_families.yaml"
    registry = json.loads(registry_path.read_text(encoding="utf-8"))
    families = [
        family
        for family in registry["families"]
        if family["family_id"] == "mt-supp-cellulose-summary-bh"
    ]
    assert len(families) == 1
    family = families[0]
    assert family["size_source"]["selection_rule"] == (
        "Concatenate direction_summary_rows, convergence_summary_rows and "
        "interaction_summary_rows in significant_results_summary_table; pass their "
        "p_value column to benjamini_hochberg, which removes missing p from the correction vector."
    )
    assert {
        "path": "src/mtp/cellulose_decr/analysis.py",
        "function": "significant_results_summary_table",
        "procedures": ["BH"],
    } in family["anchors"]
