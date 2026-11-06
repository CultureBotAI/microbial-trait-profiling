"""Nominal labels/counts for Sankey selection exports (#1430).

Synthetic supplied statistics only; the per-trait writer's plot call is mocked.
The eight existing selection fields must stay unchanged.
"""

from __future__ import annotations

import importlib
import math
import sys
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest

with patch.object(sys, "path", [str(Path(__file__).resolve().parents[1] / "scripts"), *sys.path]):
    sankeys = importlib.import_module("create_six_shared_trait_sankeys")

OLD_COLUMNS = [
    "trait_id",
    "trait_name",
    "slug",
    "base_color",
    "paired_q_value",
    "paired_p_value",
    "paired_mean_increased_minus_decreased",
    "paired_direction",
]
NEW_COLUMNS = [
    "paired_p_value_type",
    "displayed_nominal_test_count",
    "bh_family_size",
    "paired_q_value_family",
]
MIXED_P = (0.0, 1.0, math.nan, math.inf, -math.inf, 0.125)


def supplied_stats(p_values: tuple[float, ...]) -> pd.DataFrame:
    """Include an extra supplied row that the fixed display does not export."""
    records = [
        {
            "trait": trait.trait_id,
            "q_value": 0.75,
            "p_value": p,
            "mean_increased_minus_decreased": (-1) ** index * (index + 0.25),
            "direction": "increased" if index % 2 == 0 else "decreased",
        }
        for index, (trait, p) in enumerate(zip(sankeys.SHARED_TRAITS, p_values, strict=True))
    ]
    records.append(
        {
            "trait": "not-a-displayed-trait",
            "q_value": 0.0,
            "p_value": 0.0,
            "mean_increased_minus_decreased": 99.0,
            "direction": "increased",
        }
    )
    return pd.DataFrame.from_records(records).set_index("trait")


def existing_selection_columns(q_stats: pd.DataFrame) -> pd.DataFrame:
    """Independent expected export of the unchanged selection contract."""
    return pd.DataFrame(
        [
            {
                "trait_id": trait.trait_id,
                "trait_name": trait.label,
                "slug": trait.slug,
                "base_color": trait.base_color,
                "paired_q_value": q_stats.loc[trait.trait_id, "q_value"],
                "paired_p_value": q_stats.loc[trait.trait_id, "p_value"],
                "paired_mean_increased_minus_decreased": q_stats.loc[
                    trait.trait_id, "mean_increased_minus_decreased"
                ],
                "paired_direction": q_stats.loc[trait.trait_id, "direction"],
            }
            for trait in sankeys.SHARED_TRAITS
        ],
        columns=OLD_COLUMNS,
    )


@pytest.mark.parametrize(
    ("p_values", "expected_count", "expected_types"),
    [
        (MIXED_P, 3, ["nominal", "nominal", "not_tested", "not_tested", "not_tested", "nominal"]),
        ((math.nan,) * 6, 0, ["not_tested"] * 6),
        (
            (0.0, math.nan, math.nan, math.nan, math.nan, math.nan),
            1,
            ["nominal", "not_tested", "not_tested", "not_tested", "not_tested", "not_tested"],
        ),
        ((0.0, 0.2, 0.4, 0.6, 0.8, 1.0), 6, ["nominal"] * 6),
    ],
    ids=["mixed-finite", "all-unavailable", "zero-is-tested", "all-finite"],
)
def test_selection_nominal_metadata(
    p_values: tuple[float, ...], expected_count: int, expected_types: list[str]
) -> None:
    supplied = supplied_stats(p_values)
    before = supplied.copy(deep=True)

    result = sankeys.build_trait_selection(supplied)

    assert list(result.columns) == OLD_COLUMNS + NEW_COLUMNS
    assert result["paired_p_value_type"].tolist() == expected_types
    assert result["displayed_nominal_test_count"].tolist() == [expected_count] * len(result)
    assert result["bh_family_size"].tolist() == [len(supplied)] * len(result)
    assert result["paired_q_value_family"].tolist() == [
        f"pooled pair-level BH over {len(supplied)} retained catalog traits"
    ] * len(result)
    pd.testing.assert_frame_equal(result[OLD_COLUMNS], existing_selection_columns(supplied))
    pd.testing.assert_frame_equal(supplied, before)


def test_selection_uses_explicit_correction_family_size_when_present() -> None:
    supplied = supplied_stats(MIXED_P).assign(correction_family_size=340)

    result = sankeys.build_trait_selection(supplied)

    assert result["bh_family_size"].tolist() == [340] * len(result)
    assert result["paired_q_value_family"].tolist() == [
        "pooled pair-level BH over 340 retained catalog traits"
    ] * len(result)


def test_load_q_significance_records_full_catalog_size_before_shared_trait_filter() -> None:
    raw_catalog = supplied_stats(MIXED_P).reset_index()
    numeric_defaults = {
        "positive_pair_fraction": 0.5,
        "negative_pair_fraction": 0.5,
        "pair_count": 10,
        "profile_count": 20,
        "n_profiles_where_annotated": 15,
    }
    raw_without_col = raw_catalog.assign(**numeric_defaults)
    with patch.object(sankeys, "read_tsv", return_value=raw_without_col):
        loaded_without_col = sankeys.load_q_significance()
    assert len(loaded_without_col) == len(sankeys.SHARED_TRAITS)
    assert loaded_without_col["correction_family_size"].tolist() == [len(raw_catalog)] * len(
        sankeys.SHARED_TRAITS
    )

    raw_with_col = raw_without_col.assign(correction_family_size="340")
    with patch.object(sankeys, "read_tsv", return_value=raw_with_col):
        loaded_with_col = sankeys.load_q_significance()
    assert loaded_with_col["correction_family_size"].tolist() == [340] * len(sankeys.SHARED_TRAITS)


def test_selection_count_follows_exported_traits(monkeypatch: pytest.MonkeyPatch) -> None:
    supplied = supplied_stats(MIXED_P)
    chosen = (sankeys.SHARED_TRAITS[0], sankeys.SHARED_TRAITS[2], sankeys.SHARED_TRAITS[5])
    monkeypatch.setattr(sankeys, "SHARED_TRAITS", chosen)

    result = sankeys.build_trait_selection(supplied)

    assert result["trait_id"].tolist() == [trait.trait_id for trait in chosen]
    assert result["paired_p_value_type"].tolist() == ["nominal", "not_tested", "nominal"]
    assert result["displayed_nominal_test_count"].tolist() == [2, 2, 2]
    assert result["bh_family_size"].tolist() == [len(supplied)] * len(chosen)
    assert result["paired_q_value_family"].tolist() == [
        f"pooled pair-level BH over {len(supplied)} retained catalog traits"
    ] * len(chosen)
    pd.testing.assert_frame_equal(result[OLD_COLUMNS], existing_selection_columns(supplied))


def test_individual_selection_export_recounts_shown_tests(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    supplied = supplied_stats(MIXED_P)
    selection = sankeys.build_trait_selection(supplied)
    before = selection.copy(deep=True)
    empty_edges = pd.DataFrame({"analysis": pd.Series(dtype=str)})
    monkeypatch.setattr(sankeys, "OUTPUT_DIR", tmp_path)
    monkeypatch.setattr(sankeys, "REPO_ROOT", tmp_path)

    with patch.object(sankeys, "plot_single_trait_sankey"):
        sankeys.write_individual_outputs(
            empty_edges, empty_edges, empty_edges, empty_edges, selection
        )

    for trait, p in zip(sankeys.SHARED_TRAITS, MIXED_P, strict=True):
        saved = pd.read_csv(tmp_path / f"{trait.slug}_trait_selection.tsv", sep="\t")
        assert list(saved.columns) == OLD_COLUMNS + NEW_COLUMNS
        assert saved["paired_p_value_type"].tolist() == [
            "nominal" if math.isfinite(p) else "not_tested"
        ]
        assert saved["displayed_nominal_test_count"].tolist() == [int(math.isfinite(p))]
        assert saved["bh_family_size"].tolist() == [len(supplied)]
        assert saved["paired_q_value_family"].tolist() == [
            f"pooled pair-level BH over {len(supplied)} retained catalog traits"
        ]
        expected = existing_selection_columns(supplied)
        expected = expected.loc[expected["trait_id"] == trait.trait_id].reset_index(drop=True)
        pd.testing.assert_frame_equal(saved[OLD_COLUMNS], expected)
    pd.testing.assert_frame_equal(selection, before)
