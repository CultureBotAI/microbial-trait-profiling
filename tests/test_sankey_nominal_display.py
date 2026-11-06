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
SUMMARY_OLD_COLUMNS = [
    "trait_id",
    "trait_label",
    "slug",
    "base_color",
    "paired_q_value",
    "paired_mean_increased_minus_decreased",
    "paired_direction",
    "profile_count",
    "increased_profile_count",
    "decreased_profile_count",
    "annotated_profile_count",
    "increased_annotated_profile_count",
    "decreased_annotated_profile_count",
    "positive_z_profile_count",
    "increased_positive_z_profile_count",
    "decreased_positive_z_profile_count",
]
SUMMARY_NEW_COLUMNS = [
    "bh_family_size",
    "paired_q_value_family",
]
EXPECTED_FAMILY_340 = (
    "pooled pair-level BH over 340 support/variance-selected traits"
)
MIXED_P = (0.0, 1.0, math.nan, math.inf, -math.inf, 0.125)


def supplied_stats(
    p_values: tuple[float, ...],
    *,
    correction_family_size: object = 340,
    include_family_column: bool = True,
) -> pd.DataFrame:
    """Include an extra supplied row that the fixed display does not export."""
    records = [
        {
            "trait": trait.trait_id,
            "q_value": 0.75,
            "p_value": p,
            "mean_increased_minus_decreased": (-1) ** index * (index + 0.25),
            "direction": "increased" if index % 2 == 0 else "decreased",
            **({"correction_family_size": correction_family_size} if include_family_column else {}),
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
            **({"correction_family_size": correction_family_size} if include_family_column else {}),
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
    assert result["bh_family_size"].tolist() == [340] * len(result)
    assert result["paired_q_value_family"].tolist() == [EXPECTED_FAMILY_340] * len(result)
    assert "retained catalog traits" not in " ".join(result["paired_q_value_family"].tolist())
    pd.testing.assert_frame_equal(result[OLD_COLUMNS], existing_selection_columns(supplied))
    pd.testing.assert_frame_equal(supplied, before)


def test_selection_uses_explicit_correction_family_size_when_present() -> None:
    supplied = supplied_stats(MIXED_P, correction_family_size=312)

    result = sankeys.build_trait_selection(supplied)

    assert result["bh_family_size"].tolist() == [312] * len(result)
    assert result["paired_q_value_family"].tolist() == [
        "pooled pair-level BH over 312 support/variance-selected traits"
    ] * len(result)


def test_load_q_significance_and_selection_do_not_fallback_to_row_count_when_missing_or_blank() -> (
    None
):
    raw_catalog = supplied_stats(MIXED_P, include_family_column=False).reset_index()
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
    assert loaded_without_col["correction_family_size"].isna().all()
    selection_without_col = sankeys.build_trait_selection(loaded_without_col)
    assert selection_without_col["bh_family_size"].isna().all()
    assert selection_without_col["paired_q_value_family"].tolist() == [
        "family size not recorded"
    ] * len(sankeys.SHARED_TRAITS)

    raw_blank = raw_without_col.assign(correction_family_size="")
    with patch.object(sankeys, "read_tsv", return_value=raw_blank):
        loaded_blank = sankeys.load_q_significance()
    assert loaded_blank["correction_family_size"].isna().all()
    selection_blank = sankeys.build_trait_selection(loaded_blank)
    assert selection_blank["bh_family_size"].isna().all()
    assert selection_blank["paired_q_value_family"].tolist() == ["family size not recorded"] * len(
        sankeys.SHARED_TRAITS
    )

    raw_with_col = raw_without_col.assign(correction_family_size="340")
    with patch.object(sankeys, "read_tsv", return_value=raw_with_col):
        loaded_with_col = sankeys.load_q_significance()
    assert loaded_with_col["correction_family_size"].tolist() == [340] * len(sankeys.SHARED_TRAITS)
    selection_with_col = sankeys.build_trait_selection(loaded_with_col)
    assert selection_with_col["bh_family_size"].tolist() == [340] * len(sankeys.SHARED_TRAITS)
    assert selection_with_col["paired_q_value_family"].tolist() == [EXPECTED_FAMILY_340] * len(
        sankeys.SHARED_TRAITS
    )

    slice_without = (
        selection_with_col.drop(columns=["bh_family_size", "paired_q_value_family"])
        .iloc[[0]]
        .copy()
    )
    slice_meta = sankeys.selection_display_metadata(slice_without)
    assert slice_meta["bh_family_size"].isna().all()
    assert slice_meta["paired_q_value_family"].tolist() == ["family size not recorded"]


def test_build_trait_summary_appends_family_columns_additively() -> None:
    metadata = pd.DataFrame({"inc_dec": ["inc", "dec", "inc"]}, index=["p1", "p2", "p3"])
    raw_counts = pd.DataFrame(
        {trait.trait_id: [1, 0, 2] for trait in sankeys.SHARED_TRAITS},
        index=metadata.index,
    )
    zscores = pd.DataFrame(
        {trait.trait_id: [0.5, -0.5, 1.2] for trait in sankeys.SHARED_TRAITS},
        index=metadata.index,
    )
    supplied = supplied_stats(MIXED_P, correction_family_size=340)
    summary = sankeys.build_trait_summary(metadata, raw_counts, zscores, supplied)
    assert list(summary.columns) == SUMMARY_OLD_COLUMNS + SUMMARY_NEW_COLUMNS
    assert summary["bh_family_size"].tolist() == [340] * len(sankeys.SHARED_TRAITS)
    assert summary["paired_q_value_family"].tolist() == [EXPECTED_FAMILY_340] * len(
        sankeys.SHARED_TRAITS
    )

    missing_summary = sankeys.build_trait_summary(
        metadata,
        raw_counts,
        zscores,
        supplied_stats(MIXED_P, include_family_column=False),
    )
    assert list(missing_summary.columns) == SUMMARY_OLD_COLUMNS + SUMMARY_NEW_COLUMNS
    assert missing_summary["bh_family_size"].isna().all()
    assert missing_summary["paired_q_value_family"].tolist() == ["family size not recorded"] * len(
        sankeys.SHARED_TRAITS
    )

    # A blank recorded size is not filled from the row count (#1743).
    blank_summary = sankeys.build_trait_summary(
        metadata, raw_counts, zscores, supplied.assign(correction_family_size="")
    )
    assert blank_summary["bh_family_size"].isna().all()
    assert set(blank_summary["paired_q_value_family"]) == {"family size not recorded"}


@pytest.mark.parametrize("raw", [0, "0", -3, 340.7, "1.5", "abc", float("inf")])
def test_family_size_parser_rejects_sizes_that_are_not_positive_whole_numbers(raw: object) -> None:
    """Zero, negative, fractional or non-numeric sizes read as not recorded (#1743)."""
    assert sankeys._extract_bh_family_size(raw) is None
    assert sankeys._format_paired_q_value_family(raw) == "family size not recorded"


@pytest.mark.parametrize(("raw", "size"), [(340, 340), ("340", 340), (340.0, 340), ("340.0", 340)])
def test_family_size_parser_accepts_whole_sizes(raw: object, size: int) -> None:
    assert sankeys._extract_bh_family_size(raw) == size


def test_sankey_methods_paragraph_names_the_selection_q_family() -> None:
    """The six-trait selection's q names its family (#1744)."""
    methods = (Path(__file__).resolve().parents[1] / "trait-visualization/figure_methods.md").read_text(
        encoding="utf-8"
    )
    assert (
        "using `q_value <= 0.001`, where `q_value` is the pooled pair-level BH q over the "
        "support/variance-selected traits of that table (its `correction_family_size`)"
    ) in " ".join(methods.split())


def test_selection_count_follows_exported_traits(monkeypatch: pytest.MonkeyPatch) -> None:
    supplied = supplied_stats(MIXED_P)
    chosen = (sankeys.SHARED_TRAITS[0], sankeys.SHARED_TRAITS[2], sankeys.SHARED_TRAITS[5])
    monkeypatch.setattr(sankeys, "SHARED_TRAITS", chosen)

    result = sankeys.build_trait_selection(supplied)

    assert result["trait_id"].tolist() == [trait.trait_id for trait in chosen]
    assert result["paired_p_value_type"].tolist() == ["nominal", "not_tested", "nominal"]
    assert result["displayed_nominal_test_count"].tolist() == [2, 2, 2]
    assert result["bh_family_size"].tolist() == [340] * len(chosen)
    assert result["paired_q_value_family"].tolist() == [EXPECTED_FAMILY_340] * len(chosen)
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
        assert saved["bh_family_size"].tolist() == [340]
        assert saved["paired_q_value_family"].tolist() == [EXPECTED_FAMILY_340]
        expected = existing_selection_columns(supplied)
        expected = expected.loc[expected["trait_id"] == trait.trait_id].reset_index(drop=True)
        pd.testing.assert_frame_equal(saved[OLD_COLUMNS], expected)
    pd.testing.assert_frame_equal(selection, before)
