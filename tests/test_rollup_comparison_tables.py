"""Pure off/phenotype comparison contracts for #1041 (2026-10-25).

Inputs are asymmetric synthetic completed evidence and study ledgers. These
tests neither transform profiles nor recalculate a p value or correction.
"""

from __future__ import annotations

from dataclasses import fields
from typing import Any

import numpy as np
import pandas as pd
import pytest

from mtp.rollup_comparison.tables import build_tables

FA = "METPO:1000608"
PHENOTYPE_ONLY = "METPO:1000605"
OFF_ONLY = "TEST:off_only"
RANK_COLUMNS = [
    f"{effect}_rank_{order}"
    for effect in ("experiment_mean_delta_z", "study_mean_delta_z")
    for order in ("ascending", "descending", "absolute_descending")
]
Inputs = tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]


def scope(effects: dict[str, float], *, phenotype: bool) -> tuple[pd.DataFrame, pd.DataFrame]:
    studies = pd.DataFrame(index=pd.Index([f"s{i:02}" for i in range(60)], name="study_id"))
    records: list[dict[str, Any]] = []
    for position, (trait, effect) in enumerate(effects.items()):
        support = 49 if trait == "TEST:p1" else 50 if trait == FA else 60
        value = -2.0 if trait == "TEST:p3" else effect / 10
        studies[trait] = [value] * support + [0.0] * (60 - support)
        defined = effect != 0
        records.append(
            {
                "trait": trait,
                "trait_label": f"{'phenotype' if phenotype else 'off'} label {trait}",
                "composition_family_size": len(effects),
                "composition_profile_support": position + 100,
                "composition_pairs": 90,
                "composition_studies": 60,
                "experiment_mean_delta_z": effect,
                "study_mean_delta_z": float(studies[trait].mean()),
                "study_mean_delta_z_ci95_low": float(studies[trait].mean()) - 0.1,
                "study_mean_delta_z_ci95_high": float(studies[trait].mean()) + 0.2,
                "study_t_statistic": position + 0.125 if defined else np.nan,
                "study_t_p": 0.0001 if defined else np.nan,
                "study_t_q_by": 0.2 if defined else np.nan,
                "study_t_q_bh": 0.01 if defined else np.nan,
                "study_t_status": "tested" if defined else "untestable",
                "noel_family_size": len(effects),
                "noel_family_id": "noel-full-screen",
                "noel_status": "tested" if defined else "all_tied",
                "noel_studies_up": 13 if defined else 0,
                "noel_studies_down": 27 if defined else 0,
                "noel_studies_tied": 20 if defined else 60,
                "noel_studies_scored": 40 if defined else 0,
                "noel_q_by_nominal": 0.2 if defined else np.nan,
                "noel_q_bh_nominal": 0.01 if defined else np.nan,
                "ora_enriched_inc": 900 + position,
                "ora_direction_by_q": 0.00001,
                "carrier_direct_reported_taxa": 1000 + position,
                "native_selected_family_size": len(effects) - 1,
            }
        )
    evidence = pd.DataFrame(records)
    for effect_column in ("experiment_mean_delta_z", "study_mean_delta_z"):
        for order, values, ascending in (
            ("ascending", evidence[effect_column], True),
            ("descending", evidence[effect_column], False),
            ("absolute_descending", evidence[effect_column].abs(), False),
        ):
            evidence[f"{effect_column}_rank_{order}"] = values.rank(
                method="min", ascending=ascending
            ).astype("Int64")
    evidence.loc[evidence.trait == FA, "noel_q_by_nominal"] = 0.051 if phenotype else 0.05
    evidence.loc[evidence.trait == FA, "study_t_q_by"] = 0.05 if phenotype else 0.051
    evidence.loc[evidence.trait == "TEST:p2", "noel_q_by_nominal"] = 0.04 if phenotype else 0.08
    evidence.loc[evidence.trait == "TEST:n1", "noel_q_by_nominal"] = 0.01
    evidence.loc[evidence.trait == "TEST:p3", "study_t_q_by"] = 0.001
    return evidence, studies


@pytest.fixture
def inputs() -> Inputs:
    effects = {
        "TEST:p1": 10.0,
        "TEST:p2": 9.0,
        "TEST:p3": 8.0,
        "TEST:p4": 7.0,
        FA: 6.0,
        "TEST:tie": 6.0,
        "TEST:p7": 5.0,
        "TEST:p8": 4.0,
        "TEST:p9": 3.0,
        "TEST:p10": 2.0,
        "TEST:p11": 1.0,
        "TEST:n1": -9.0,
        "TEST:n2": -8.0,
        OFF_ONLY: -0.5,
        "TEST:zero": 0.0,
    }
    off, off_studies = scope(effects, phenotype=False)
    other = {trait: effect for trait, effect in effects.items() if trait != OFF_ONLY}
    other.update({PHENOTYPE_ONLY: 20.0, "TEST:new_high": 19.0})
    phenotype, phenotype_studies = scope(other, phenotype=True)
    return off.iloc[::-1], phenotype, off_studies, phenotype_studies.iloc[::-1]


def transition(inputs: Inputs, view: str) -> pd.DataFrame:
    rows = build_tables(*inputs).rank_transitions
    return rows.loc[rows.view == view].set_index("trait")


def test_catalog_uses_qualified_identity_and_retains_scope_labels(inputs: Inputs) -> None:
    result = build_tables(*inputs)
    catalog = result.trait_catalog.set_index("trait")
    assert len(result.shared_traits) == 14
    assert len(catalog) == 17
    assert catalog.loc[FA, "membership"] == "shared"
    assert catalog.loc[PHENOTYPE_ONLY, "membership"] == "phenotype_only"
    assert catalog.loc[OFF_ONLY, "membership"] == "off_only"
    assert bool(catalog.loc[FA, "in_off"])
    assert bool(catalog.loc[FA, "in_phenotype"])
    assert not bool(catalog.loc[PHENOTYPE_ONLY, "in_off"])
    assert catalog.loc[FA, "off_trait_label"] != catalog.loc[FA, "phenotype_trait_label"]
    assert pd.isna(catalog.loc[PHENOTYPE_ONLY, "off_trait_label"])


def test_side_by_side_preserves_all_parent_values_and_native_full_family_ranks(
    inputs: Inputs,
) -> None:
    result = build_tables(*inputs).shared_traits.set_index("trait")
    for prefix, source in zip(("off", "phenotype"), inputs[:2], strict=True):
        original = source.set_index("trait").loc[result.index]
        for column in original:
            pd.testing.assert_series_equal(
                result[f"{prefix}_{column}"], original[column], check_names=False
            )
        assert result.loc[FA, f"{prefix}_comparison_delta_z_mean"] == 6
        assert result.loc[FA, f"{prefix}_study_mean_delta_z"] == pytest.approx(0.5)
        assert result.loc[FA, f"{prefix}_study_mean_delta_z_ci95_low"] == pytest.approx(0.4)
        assert result.loc[FA, f"{prefix}_study_mean_delta_z_ci95_high"] == pytest.approx(0.7)
    assert result.loc[FA, "off_experiment_mean_delta_z_rank_descending"] == 5
    assert result.loc[FA, "phenotype_experiment_mean_delta_z_rank_descending"] == 7
    assert result.loc[FA, "off_composition_family_size"] == 15
    assert result.loc[FA, "phenotype_composition_family_size"] == 16
    assert len([name for name in result if str(name).removeprefix("off_") in RANK_COLUMNS]) == 6


def test_views_cover_union_and_keep_missing_scope_unavailable(inputs: Inputs) -> None:
    result = build_tables(*inputs)
    rows = result.rank_transitions
    assert len(rows) == 170
    assert rows.groupby("view").size().eq(17).all()
    row = transition(inputs, "comparison_mean_positive_top5_all").loc[PHENOTYPE_ONLY]
    assert row.membership_state == "entered"
    assert row.off_status == "unavailable"
    assert not row.off_available
    assert pd.isna(row.off_selected)
    assert row.phenotype_status == "selected"
    row = transition(inputs, "comparison_mean_negative_top5_all").loc[OFF_ONLY]
    assert row.membership_state == "left"
    assert row.phenotype_status == "unavailable"


def test_eligibility_records_direction_support_testability_and_absent_scope(inputs: Inputs) -> None:
    active = transition(inputs, "comparison_mean_positive_top5_active50")
    assert bool(active.loc[FA, "off_eligible"])
    assert not bool(active.loc["TEST:p1", "off_eligible"])
    assert not bool(active.loc["TEST:n1", "off_eligible"])
    assert pd.isna(active.loc[PHENOTYPE_ONLY, "off_eligible"])
    noel = transition(inputs, "noel_by_0_05")
    assert bool(noel.loc[FA, "off_eligible"])
    assert not bool(noel.loc["TEST:zero", "off_eligible"])


def test_full_scope_competition_ties_cross_the_top_five_cutoff(inputs: Inputs) -> None:
    view = transition(inputs, "comparison_mean_positive_top5_all")
    assert view.off_selected.fillna(False).sum() == 6
    assert view.loc[FA, "off_view_rank"] == view.loc["TEST:tie", "off_view_rank"] == 5
    assert view.loc[FA, "phenotype_view_rank"] == 7
    assert view.loc[FA, "membership_state"] == "left"
    assert view.loc["TEST:tie", "membership_state"] == "left"
    assert view.loc["TEST:p3", "off_status"] == "selected"
    assert inputs[0].loc[inputs[0].trait == "TEST:p3", "study_mean_delta_z"].lt(0).all()
    assert view.loc["TEST:zero", "off_status"] == "ineligible"
    assert view.loc["TEST:n1", "off_status"] == "ineligible"


def test_top_ten_uses_comparison_effect_and_both_direction_views(inputs: Inputs) -> None:
    positive = transition(inputs, "comparison_mean_positive_top10_all")
    assert positive.loc["TEST:p10", "off_status"] == "selected"
    assert positive.loc["TEST:p10", "phenotype_status"] == "not_selected"
    negative = transition(inputs, "comparison_mean_negative_top10_all")
    assert negative.loc["TEST:n1", "off_view_rank"] == 1
    assert negative.loc["TEST:n2", "off_view_rank"] == 2
    assert negative.loc["TEST:zero", "membership_state"] == "neither"
    assert negative.loc["TEST:p3", "off_status"] == "ineligible"


def test_active_support_is_study_delta_z_and_includes_fifty(inputs: Inputs) -> None:
    result = build_tables(*inputs).shared_traits.set_index("trait")
    assert result.loc[FA, "off_sdz_active_studies"] == 50
    assert result.loc[FA, "off_sdz_up"] == 50
    assert result.loc[FA, "off_sdz_down"] == 0
    assert result.loc[FA, "off_sdz_tied"] == 10
    assert result.loc[FA, "off_noel_studies_up"] == 13
    assert result.loc["TEST:p1", "off_sdz_active_studies"] == 49
    active = transition(inputs, "comparison_mean_positive_top5_active50")
    assert active.loc["TEST:p1", "off_status"] == "ineligible"
    assert active.loc[FA, "off_view_rank"] == 4
    assert active.loc[FA, "off_status"] == "selected"
    assert active.loc[FA, "phenotype_view_rank"] == 6
    assert result.loc[FA, "off_experiment_mean_delta_z_rank_descending"] == 5


def test_activity_uses_strict_epsilon_and_preserves_native_votes(inputs: Inputs) -> None:
    off, phenotype, studies, other_studies = (frame.copy(deep=True) for frame in inputs)
    studies[FA] = [2e-12, -2e-12, 1e-12, -1e-12] + [0.0] * 56
    off.loc[off.trait == FA, "study_mean_delta_z"] = float(studies[FA].mean())
    result = build_tables(off, phenotype, studies, other_studies).shared_traits.set_index("trait")
    assert result.loc[FA, "off_sdz_up"] == result.loc[FA, "off_sdz_down"] == 1
    assert result.loc[FA, "off_sdz_active_studies"] == 2
    assert result.loc[FA, "off_sdz_tied"] == 58
    assert result.loc[FA, "off_noel_studies_down"] == 27


def test_by_views_remain_separate_inclusive_and_preserve_untestable_rows(inputs: Inputs) -> None:
    noel = transition(inputs, "noel_by_0_05")
    t_test = transition(inputs, "study_t_by_0_05")
    assert noel.loc[FA, "membership_state"] == "left"
    assert t_test.loc[FA, "membership_state"] == "entered"
    assert noel.loc["TEST:p2", "membership_state"] == "entered"
    assert noel.loc["TEST:n1", "membership_state"] == "retained"
    assert noel.loc["TEST:p3", "membership_state"] == "neither"
    assert t_test.loc["TEST:p3", "membership_state"] == "retained"
    assert noel.loc["TEST:zero", "off_status"] == "untestable"
    assert t_test.loc["TEST:zero", "off_status"] == "untestable"
    assert pd.isna(noel.loc[FA, "off_view_rank"])
    assert "ora" not in " ".join(build_tables(*inputs).rank_transitions.view.unique())


def test_string_by_probability_retains_value_above_the_selection_boundary(inputs: Inputs) -> None:
    off, phenotype, studies, other_studies = (frame.copy(deep=True) for frame in inputs)
    off["noel_q_by_nominal"] = off.noel_q_by_nominal.astype(str)
    off.loc[off.trait == FA, "noel_q_by_nominal"] = "0.05000000000000001"
    view = transition((off, phenotype, studies, other_studies), "noel_by_0_05")
    assert view.loc[FA, "off_status"] == "not_selected"
    assert view.loc[FA, "off_selection_value"] == float("0.05000000000000001")


def test_adjacent_string_effects_remain_distinct_in_comparison_rank_views(inputs: Inputs) -> None:
    off, phenotype, studies, other_studies = (frame.copy(deep=True) for frame in inputs)
    off["experiment_mean_delta_z"] = off.experiment_mean_delta_z.astype(str)
    off.loc[off.trait == FA, "experiment_mean_delta_z"] = "7.000000000000001"
    off.loc[off.trait == "TEST:p4", "experiment_mean_delta_z"] = "7.0"
    view = transition((off, phenotype, studies, other_studies), "comparison_mean_positive_top5_all")
    assert view.loc[FA, "off_view_rank"] == 4
    assert view.loc["TEST:p4", "off_view_rank"] == 5
    assert view.loc[FA, "off_selection_value"] == float("7.000000000000001")


def test_input_permutation_does_not_change_results_or_modify_inputs(inputs: Inputs) -> None:
    saved = tuple(frame.copy(deep=True) for frame in inputs)
    expected = build_tables(*inputs)
    permuted = tuple(frame.iloc[::-1] for frame in inputs)
    observed = build_tables(*permuted)
    for field in fields(expected):
        pd.testing.assert_frame_equal(getattr(expected, field.name), getattr(observed, field.name))
    for before, after in zip(saved, inputs, strict=True):
        pd.testing.assert_frame_equal(before, after)


@pytest.mark.parametrize(
    "change",
    [
        "duplicate_trait",
        "blank_trait",
        "missing_ledger_trait",
        "duplicate_study",
        "nonfinite_study",
        "wrong_mean",
        "wrong_family",
    ],
)
def test_invalid_identity_or_ledger_alignment_fails(inputs: Inputs, change: str) -> None:
    off, phenotype, studies, other_studies = (frame.copy(deep=True) for frame in inputs)
    if change == "duplicate_trait":
        off = pd.concat([off, off.iloc[[0]]], ignore_index=True)
    elif change == "blank_trait":
        off.loc[off.index[0], "trait"] = " "
    elif change == "missing_ledger_trait":
        studies = studies.drop(columns=FA)
    elif change == "duplicate_study":
        studies = pd.concat([studies, studies.iloc[[0]]])
    elif change == "nonfinite_study":
        studies.iloc[0, 0] = np.nan
    elif change == "wrong_mean":
        off.loc[off.trait == FA, "study_mean_delta_z"] += 1
    else:
        off["composition_family_size"] = 14
    with pytest.raises(ValueError, match=r"identity|catalog|finite|mean|family"):
        build_tables(off, phenotype, studies, other_studies)


@pytest.mark.parametrize("column", ["sdz_up", "comparison_delta_z_mean"])
def test_derived_fields_cannot_silently_overwrite_parent_metadata(
    inputs: Inputs, column: str
) -> None:
    off, phenotype, studies, other_studies = inputs
    off = off.assign(**{column: 987})
    with pytest.raises(ValueError, match="collides"):
        build_tables(off, phenotype, studies, other_studies)
