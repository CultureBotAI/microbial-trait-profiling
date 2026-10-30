"""Exact stool population and support-first condition selection for #1093."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import pytest

from mtp.stool_group_contrasts import (
    DELTA_FRACTION_COLUMN,
    DELTA_PERCENTAGE_POINTS_COLUMN,
    exact_stool_profiles,
    select_stool_group_contrasts,
)


def cells() -> pd.DataFrame:
    """More rows or more significance must not outweigh unique-study support."""
    records: list[dict[str, object]] = []
    groups = [("Unsupported", 90), ("Beta", 8), ("Alpha", 8), ("Gamma", 5)]
    for condition, studies in groups:
        for trait in ("METPO:1", "METPO:2", "METPO:3"):
            records.append(
                {
                    "view_key": "grouped_direction_contrast",
                    "body_site_group": "Feces",
                    "condition": condition,
                    "trait_id": trait,
                    "study_count": studies,
                    "passes_support": condition != "Unsupported" and trait == "METPO:2",
                    DELTA_FRACTION_COLUMN: 0.125 if trait == "METPO:1" else -0.375,
                    "fisher_p_value": 0.99 if condition == "Alpha" else 0.0001,
                    "fisher_q_value": 1.0 if condition == "Alpha" else 0.0002,
                    "native_note": f"{condition}/{trait}",
                }
            )
    return pd.DataFrame(records, index=[f"cell-{index}" for index in range(len(records))])


def test_exact_stool_requires_literal_metadata_body_site_and_preserves_index() -> None:
    metadata = pd.DataFrame(
        {
            "body_site": ["Feces", "Feces,Blood", "feces", " Feces", "Feces ", "stool", None],
            "body_site_group": ["Other", "Feces", "Feces", "Feces", "Feces", "Feces", "Feces"],
            "study_id": list("abcdefg"),
        },
        index=[f"profile-{index}" for index in range(7)],
    )
    before = metadata.copy(deep=True)
    result = exact_stool_profiles(metadata)
    pd.testing.assert_frame_equal(result, metadata.iloc[[0]])
    result.loc["profile-0", "study_id"] = "changed"
    pd.testing.assert_frame_equal(metadata, before)


def test_exact_stool_never_substitutes_grouped_site() -> None:
    with pytest.raises(ValueError, match="body_site"):
        exact_stool_profiles(pd.DataFrame({"body_site_group": ["Feces"]}))


def test_group_selection_uses_support_then_studies_then_name_without_inference_selection() -> None:
    source = cells()
    before = source.copy(deep=True)
    selected, selection = select_stool_group_contrasts(source, group_count=2)
    expected = source.loc[source["condition"].isin(["Alpha", "Beta"])]
    pd.testing.assert_frame_equal(selected.drop(columns=DELTA_PERCENTAGE_POINTS_COLUMN), expected)
    assert selected[DELTA_PERCENTAGE_POINTS_COLUMN].tolist() == [12.5, -37.5, -37.5] * 2
    assert selected["passes_support"].tolist() == [False, True, False] * 2
    assert selection["condition"].tolist() == ["Alpha", "Beta", "Gamma", "Unsupported"]
    assert selection["stool_study_count"].tolist() == [8, 8, 5, 90]
    assert selection["eligible"].tolist() == [True, True, True, False]
    assert selection["rank"].iloc[:3].tolist() == [1, 2, 3]
    assert pd.isna(selection["rank"].iloc[3])
    assert selection["chosen"].tolist() == [True, True, False, False]
    pd.testing.assert_frame_equal(source, before)


def test_default_six_is_top_six_eligible_groups_without_padding_or_row_count_weighting() -> None:
    source = cells()
    extra = []
    for index in range(4):
        frame = source.loc[source["condition"].eq("Gamma")].copy()
        frame["condition"] = f"Extra {index}"
        frame["study_count"] = index + 10
        extra.append(frame)
    source = pd.concat([source, *extra], ignore_index=True)
    selected, selection = select_stool_group_contrasts(source)
    assert selection.loc[selection["chosen"], "condition"].tolist() == [
        "Extra 3",
        "Extra 2",
        "Extra 1",
        "Extra 0",
        "Alpha",
        "Beta",
    ]
    assert set(selected["condition"]) == {"Alpha", "Beta", *(f"Extra {i}" for i in range(4))}
    smaller, small_selection = select_stool_group_contrasts(cells())
    assert set(smaller["condition"]) == {"Alpha", "Beta", "Gamma"}
    assert small_selection["chosen"].sum() == 3


def test_group_ranking_is_independent_of_input_order_and_effect_size() -> None:
    source = cells().sample(frac=1, random_state=1093)
    source[DELTA_FRACTION_COLUMN] = 0.0
    source["fisher_p_value"] = np.nan
    source["fisher_q_value"] = np.nan
    _, selection = select_stool_group_contrasts(source, group_count=2)
    assert selection.loc[selection["chosen"], "condition"].tolist() == ["Alpha", "Beta"]


def test_undefined_delta_is_retained_without_changing_support_or_significance() -> None:
    source = cells()
    source.loc["cell-4", DELTA_FRACTION_COLUMN] = np.nan
    selected, _ = select_stool_group_contrasts(source)
    assert pd.isna(selected.loc["cell-4", DELTA_PERCENTAGE_POINTS_COLUMN])
    assert bool(selected.loc["cell-4", "passes_support"])
    assert selected.loc["cell-4", "fisher_p_value"] == source.loc["cell-4", "fisher_p_value"]


def test_no_eligible_groups_and_empty_cells_have_stable_empty_outputs() -> None:
    source = cells()
    source["passes_support"] = False
    selected, selection = select_stool_group_contrasts(source)
    assert selected.empty
    assert len(selection) == 4
    assert not selection["eligible"].any()
    assert selection["rank"].isna().all()
    assert not selection["chosen"].any()
    empty, empty_selection = select_stool_group_contrasts(source.iloc[:0])
    assert empty.empty
    assert empty_selection.empty
    assert list(empty.columns) == [*source.columns, DELTA_PERCENTAGE_POINTS_COLUMN]
    assert list(empty_selection.columns) == list(selection.columns)


@pytest.mark.parametrize("group_count", [0, -1, True, False, 2.0, "6", None])
def test_group_count_must_be_a_positive_nonboolean_integer(group_count: Any) -> None:
    with pytest.raises(ValueError, match="group_count"):
        select_stool_group_contrasts(cells(), group_count=group_count)


@pytest.mark.parametrize("count", [0, -1, True, False, 1.5, "8", np.nan, np.inf, None])
def test_malformed_native_study_counts_are_rejected(count: Any) -> None:
    source = cells().astype({"study_count": object})
    source.loc["cell-0", "study_count"] = count
    with pytest.raises(ValueError, match="study_count"):
        select_stool_group_contrasts(source)


@pytest.mark.parametrize("value", ["True", "False", 1, 0, None, np.nan])
def test_support_is_not_coerced_from_truthy_strings_or_numbers(value: Any) -> None:
    source = cells().astype({"passes_support": object})
    source.loc["cell-0", "passes_support"] = value
    with pytest.raises(ValueError, match="passes_support"):
        select_stool_group_contrasts(source)


def test_inconsistent_counts_and_duplicate_native_cells_are_rejected() -> None:
    source = cells()
    source.loc["cell-1", "study_count"] = 91
    with pytest.raises(ValueError, match="study_count"):
        select_stool_group_contrasts(source)
    with pytest.raises(ValueError, match=r"[Dd]uplicate"):
        select_stool_group_contrasts(pd.concat([cells(), cells().iloc[[0]]]))


def test_dataframe_index_is_not_a_cell_identity_and_duplicate_index_is_preserved() -> None:
    source = cells()
    source.index = pd.Index([0] * len(source))
    selected, _ = select_stool_group_contrasts(source, group_count=2)
    assert len(selected) == 6
    assert selected.index.tolist() == [0] * 6
    assert selected[DELTA_PERCENTAGE_POINTS_COLUMN].tolist() == [12.5, -37.5, -37.5] * 2


@pytest.mark.parametrize(
    ("column", "value"),
    [
        ("view_key", "individual_conditions_direction_contrast"),
        ("body_site_group", "Feces,Blood"),
        ("condition", ""),
        ("trait_id", " "),
    ],
)
def test_non_native_or_ambiguous_cells_are_rejected(column: str, value: str) -> None:
    source = cells()
    source.loc["cell-0", column] = value
    with pytest.raises(ValueError, match=column):
        select_stool_group_contrasts(source)


@pytest.mark.parametrize("column", ["condition", "study_count", "passes_support", "view_key"])
def test_missing_required_columns_are_rejected(column: str) -> None:
    with pytest.raises(ValueError, match=column):
        select_stool_group_contrasts(cells().drop(columns=column))


def test_existing_percentage_point_column_is_not_silently_overwritten() -> None:
    source = cells()
    source[DELTA_PERCENTAGE_POINTS_COLUMN] = 123.0
    with pytest.raises(ValueError, match=DELTA_PERCENTAGE_POINTS_COLUMN):
        select_stool_group_contrasts(source)
