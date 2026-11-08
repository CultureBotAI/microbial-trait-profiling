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
    "bh_unfiltered_trait_count",
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
    "bh_unfiltered_trait_count",
    "paired_q_value_family",
]
SELECTION_RULE_NOTE = (
    "select_stratum_traits support, positive-variance, and MediaDive-exclusion filter"
)
EXPECTED_FAMILY_340 = (
    f"pooled pair-level BH over 340 support/variance-selected traits "
    f"(of 344 retained; {SELECTION_RULE_NOTE})"
)
EXPECTED_FAMILY_340_UNRECORDED_RETAINED = (
    f"pooled pair-level BH over 340 support/variance-selected traits "
    f"(retained count not recorded; {SELECTION_RULE_NOTE})"
)
MIXED_P = (0.0, 1.0, math.nan, math.inf, -math.inf, 0.125)


def supplied_stats(
    p_values: tuple[float, ...],
    *,
    correction_family_size: object = 340,
    retained_trait_count: object = 344,
    include_family_column: bool = True,
    include_retained_column: bool = True,
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
            **({"retained_trait_count": retained_trait_count} if include_retained_column else {}),
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
            **({"retained_trait_count": retained_trait_count} if include_retained_column else {}),
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
    assert result["bh_family_size"].dropna().tolist() == [340] * len(result)
    assert result["bh_unfiltered_trait_count"].dropna().tolist() == [344] * len(result)
    assert result["paired_q_value_family"].tolist() == [EXPECTED_FAMILY_340] * len(result)
    assert "retained catalog traits" not in " ".join(result["paired_q_value_family"].tolist())
    pd.testing.assert_frame_equal(result[OLD_COLUMNS], existing_selection_columns(supplied))
    pd.testing.assert_frame_equal(supplied, before)


def test_selection_uses_explicit_correction_family_size_when_present() -> None:
    supplied = supplied_stats(MIXED_P, correction_family_size=312, retained_trait_count=317)

    result = sankeys.build_trait_selection(supplied)

    assert result["bh_family_size"].dropna().tolist() == [312] * len(result)
    assert result["bh_unfiltered_trait_count"].dropna().tolist() == [317] * len(result)
    assert result["paired_q_value_family"].tolist() == [
        f"pooled pair-level BH over 312 support/variance-selected traits "
        f"(of 317 retained; {SELECTION_RULE_NOTE})"
    ] * len(result)


def test_load_q_significance_and_selection_do_not_fallback_to_row_count_when_missing_or_blank() -> (
    None
):
    raw_catalog = supplied_stats(
        MIXED_P, include_family_column=False, include_retained_column=False
    ).reset_index()
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
    assert loaded_without_col["retained_trait_count"].isna().all()
    selection_without_col = sankeys.build_trait_selection(loaded_without_col)
    assert selection_without_col["bh_family_size"].isna().all()
    assert selection_without_col["bh_unfiltered_trait_count"].isna().all()
    assert selection_without_col["paired_q_value_family"].tolist() == [
        "family size not recorded"
    ] * len(sankeys.SHARED_TRAITS)

    raw_blank = raw_without_col.assign(correction_family_size="", retained_trait_count="")
    with patch.object(sankeys, "read_tsv", return_value=raw_blank):
        loaded_blank = sankeys.load_q_significance()
    assert loaded_blank["correction_family_size"].isna().all()
    assert loaded_blank["retained_trait_count"].isna().all()
    selection_blank = sankeys.build_trait_selection(loaded_blank)
    assert selection_blank["bh_family_size"].isna().all()
    assert selection_blank["bh_unfiltered_trait_count"].isna().all()
    assert selection_blank["paired_q_value_family"].tolist() == ["family size not recorded"] * len(
        sankeys.SHARED_TRAITS
    )

    raw_with_col = raw_without_col.assign(correction_family_size="340", retained_trait_count="344")
    with patch.object(sankeys, "read_tsv", return_value=raw_with_col):
        loaded_with_col = sankeys.load_q_significance()
    assert loaded_with_col["correction_family_size"].tolist() == [340] * len(sankeys.SHARED_TRAITS)
    assert loaded_with_col["retained_trait_count"].tolist() == [344] * len(sankeys.SHARED_TRAITS)
    selection_with_col = sankeys.build_trait_selection(loaded_with_col)
    assert selection_with_col["bh_family_size"].dropna().tolist() == [340] * len(
        sankeys.SHARED_TRAITS
    )
    assert selection_with_col["bh_unfiltered_trait_count"].dropna().tolist() == [344] * len(
        sankeys.SHARED_TRAITS
    )
    assert selection_with_col["paired_q_value_family"].tolist() == [EXPECTED_FAMILY_340] * len(
        sankeys.SHARED_TRAITS
    )

    raw_only_family = raw_without_col.assign(correction_family_size="340")
    with patch.object(sankeys, "read_tsv", return_value=raw_only_family):
        loaded_only_family = sankeys.load_q_significance()
    selection_only_family = sankeys.build_trait_selection(loaded_only_family)
    assert selection_only_family["bh_family_size"].dropna().tolist() == [340] * len(
        sankeys.SHARED_TRAITS
    )
    assert selection_only_family["bh_unfiltered_trait_count"].isna().all()
    assert selection_only_family["paired_q_value_family"].tolist() == [
        EXPECTED_FAMILY_340_UNRECORDED_RETAINED
    ] * len(sankeys.SHARED_TRAITS)

    slice_without = (
        selection_with_col.drop(
            columns=["bh_family_size", "bh_unfiltered_trait_count", "paired_q_value_family"]
        )
        .iloc[[0]]
        .copy()
    )
    slice_meta = sankeys.selection_display_metadata(slice_without)
    assert slice_meta["bh_family_size"].isna().all()
    assert slice_meta["bh_unfiltered_trait_count"].isna().all()
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
    supplied = supplied_stats(MIXED_P, correction_family_size=340, retained_trait_count=344)
    summary = sankeys.build_trait_summary(metadata, raw_counts, zscores, supplied)
    assert list(summary.columns) == SUMMARY_OLD_COLUMNS + SUMMARY_NEW_COLUMNS
    assert summary["bh_family_size"].dropna().tolist() == [340] * len(sankeys.SHARED_TRAITS)
    assert summary["bh_unfiltered_trait_count"].dropna().tolist() == [344] * len(
        sankeys.SHARED_TRAITS
    )
    assert summary["paired_q_value_family"].tolist() == [EXPECTED_FAMILY_340] * len(
        sankeys.SHARED_TRAITS
    )

    missing_summary = sankeys.build_trait_summary(
        metadata,
        raw_counts,
        zscores,
        supplied_stats(MIXED_P, include_family_column=False, include_retained_column=False),
    )
    assert list(missing_summary.columns) == SUMMARY_OLD_COLUMNS + SUMMARY_NEW_COLUMNS
    assert missing_summary["bh_family_size"].isna().all()
    assert missing_summary["bh_unfiltered_trait_count"].isna().all()
    assert missing_summary["paired_q_value_family"].tolist() == ["family size not recorded"] * len(
        sankeys.SHARED_TRAITS
    )

    # A blank recorded size is not filled from the row count (#1743).
    blank_summary = sankeys.build_trait_summary(
        metadata,
        raw_counts,
        zscores,
        supplied.assign(correction_family_size="", retained_trait_count=""),
    )
    assert blank_summary["bh_family_size"].isna().all()
    assert blank_summary["bh_unfiltered_trait_count"].isna().all()
    assert set(blank_summary["paired_q_value_family"]) == {"family size not recorded"}

    unrecorded_retained_summary = sankeys.build_trait_summary(
        metadata,
        raw_counts,
        zscores,
        supplied.assign(retained_trait_count=""),
    )
    assert unrecorded_retained_summary["bh_family_size"].tolist() == [340] * len(
        sankeys.SHARED_TRAITS
    )
    assert unrecorded_retained_summary["bh_unfiltered_trait_count"].isna().all()
    assert unrecorded_retained_summary["paired_q_value_family"].tolist() == [
        EXPECTED_FAMILY_340_UNRECORDED_RETAINED
    ] * len(sankeys.SHARED_TRAITS)


@pytest.mark.parametrize("raw", [0, "0", -3, 340.7, "1.5", "abc", float("inf")])
def test_family_size_parser_rejects_sizes_that_are_not_positive_whole_numbers(raw: object) -> None:
    """Zero, negative, fractional or non-numeric sizes read as not recorded (#1743)."""
    assert sankeys._extract_bh_family_size(raw) is None
    assert sankeys._format_paired_q_value_family(raw) == "family size not recorded"
    assert (
        sankeys._format_paired_q_value_family(340, raw) == EXPECTED_FAMILY_340_UNRECORDED_RETAINED
    )


@pytest.mark.parametrize(("raw", "size"), [(340, 340), ("340", 340), (340.0, 340), ("340.0", 340)])
def test_family_size_parser_accepts_whole_sizes(raw: object, size: int) -> None:
    assert sankeys._extract_bh_family_size(raw) == size


def test_retained_trait_count_below_family_size_is_refused(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Refuse retained count below family size with ValueError (#1812, #1872, #1873, #1874)."""
    with pytest.raises(ValueError, match=r"300.*340"):
        sankeys._format_paired_q_value_family(340, 300)
    with pytest.raises(ValueError, match=r"339.*340.*METPO:1000608"):
        sankeys._format_paired_q_value_family(340, 339, trait_id="METPO:1000608")
    expected_equal = (
        f"pooled pair-level BH over 340 support/variance-selected traits "
        f"(of 340 retained; {SELECTION_RULE_NOTE})"
    )
    equal_stats = supplied_stats(MIXED_P, correction_family_size=340, retained_trait_count=340)
    metadata = pd.DataFrame({"inc_dec": ["inc", "dec", "inc"]}, index=["p1", "p2", "p3"])
    raw_counts = pd.DataFrame(
        {trait.trait_id: [1, 0, 2] for trait in sankeys.SHARED_TRAITS},
        index=metadata.index,
    )
    zscores = pd.DataFrame(
        {trait.trait_id: [0.5, -0.5, 1.2] for trait in sankeys.SHARED_TRAITS},
        index=metadata.index,
    )
    try:
        assert sankeys._format_paired_q_value_family(340, 340) == expected_equal
        assert sankeys._extract_retained_trait_count(340, 340, trait_id="METPO:1000608") == 340
        assert sankeys._extract_retained_trait_count(344, 340, trait_id="METPO:1000608") == 344
        assert sankeys._extract_retained_trait_count(None, 340, trait_id="METPO:1000608") is None
        equal_selection = sankeys.build_trait_selection(equal_stats)
        equal_summary = sankeys.build_trait_summary(metadata, raw_counts, zscores, equal_stats)
    except ValueError as exc:
        pytest.fail(f"consistent counts refused: {exc}")

    assert equal_selection["bh_unfiltered_trait_count"].dropna().tolist() == [340] * len(
        sankeys.SHARED_TRAITS
    )
    assert equal_selection["paired_q_value_family"].tolist() == [expected_equal] * len(
        sankeys.SHARED_TRAITS
    )

    inconsistent_stats = supplied_stats(
        MIXED_P, correction_family_size=340, retained_trait_count=300
    )
    with pytest.raises(ValueError, match=r"300.*340.*METPO:1000608"):
        sankeys.build_trait_selection(inconsistent_stats)

    raw_inconsistent = inconsistent_stats.reset_index().assign(
        positive_pair_fraction=0.5,
        negative_pair_fraction=0.5,
        pair_count=10,
        profile_count=20,
        n_profiles_where_annotated=15,
    )
    with (
        patch.object(sankeys, "read_tsv", return_value=raw_inconsistent),
        pytest.raises(ValueError, match=r"300.*340.*METPO:1000608"),
    ):
        sankeys.load_q_significance()

    # #1874: Every displayed row must be checked, not only the first row.
    third_trait_id = sankeys.SHARED_TRAITS[2].trait_id
    single_row_inconsistent = (
        supplied_stats(MIXED_P, correction_family_size=340, retained_trait_count=344)
        .reset_index()
        .assign(
            positive_pair_fraction=0.5,
            negative_pair_fraction=0.5,
            pair_count=10,
            profile_count=20,
            n_profiles_where_annotated=15,
        )
    )
    single_row_inconsistent.loc[
        single_row_inconsistent["trait"] == third_trait_id, "retained_trait_count"
    ] = 300
    with (
        patch.object(sankeys, "read_tsv", return_value=single_row_inconsistent),
        pytest.raises(ValueError, match=rf"300.*340.*{third_trait_id}"),
    ):
        sankeys.load_q_significance()

    slice_without = (
        equal_selection.drop(
            columns=["bh_family_size", "bh_unfiltered_trait_count", "paired_q_value_family"]
        )
        .iloc[[0]]
        .copy()
    )
    with pytest.raises(ValueError, match=r"300.*340.*METPO:1000608"):
        sankeys.selection_display_metadata(
            slice_without,
            default_bh_family_size=340,
            default_bh_unfiltered_trait_count=300,
        )
    with pytest.raises(ValueError, match=r"300.*340.*METPO:1000608"):
        sankeys.selection_display_metadata(
            slice_without.assign(bh_family_size=340, bh_unfiltered_trait_count=300)
        )

    with pytest.raises(ValueError, match=r"300.*340.*METPO:1000608"):
        sankeys.build_trait_summary(metadata, raw_counts, zscores, inconsistent_stats)

    assert equal_summary["bh_unfiltered_trait_count"].dropna().tolist() == [340] * len(
        sankeys.SHARED_TRAITS
    )
    assert equal_summary["paired_q_value_family"].tolist() == [expected_equal] * len(
        sankeys.SHARED_TRAITS
    )

    # #1872: main() must validate inputs before reset_output_dir() empties OUTPUT_DIR.
    out_dir = tmp_path / "six-shared-traits"
    out_dir.mkdir(parents=True, exist_ok=True)
    preserved_files = {
        "README.md": "# Existing Sankey README\n",
        "six_shared_trait_selection.tsv": "trait_id\tpaired_q_value\nMETPO:1000608\t1e-6\n",
        "six_shared_trait_summary.tsv": "trait_id\tprofile_count\nMETPO:1000608\t20\n",
    }
    for filename, content in preserved_files.items():
        (out_dir / filename).write_text(content, encoding="utf-8")

    monkeypatch.setattr(sankeys, "INVESTIGATION_DIR", tmp_path)
    monkeypatch.setattr(sankeys, "OUTPUT_DIR", out_dir)
    monkeypatch.setattr(sankeys, "load_profile_metadata", lambda: metadata)
    monkeypatch.setattr(sankeys, "load_raw_trait_counts", lambda _traits, _meta: raw_counts)
    monkeypatch.setattr(sankeys, "build_zscore_matrix", lambda _counts, _meta: zscores)
    monkeypatch.setattr(sankeys, "read_tsv", lambda _path: single_row_inconsistent)

    with pytest.raises(ValueError, match=rf"300.*340.*{third_trait_id}"):
        sankeys.main()

    assert sorted(p.name for p in out_dir.iterdir()) == sorted(preserved_files)
    for filename, content in preserved_files.items():
        assert (out_dir / filename).read_text(encoding="utf-8") == content


def test_sankey_methods_paragraph_names_the_selection_q_family() -> None:
    """The six-trait selection's q names its family, retained count, and rule (#1742, #1811)."""
    methods = (
        Path(__file__).resolve().parents[1] / "trait-visualization/figure_methods.md"
    ).read_text(encoding="utf-8")
    normalized = " ".join(methods.split())
    assert "`retained_trait_count` after the" not in normalized
    assert (
        "using `q_value <= 0.001`, where `q_value` is the pooled pair-level BH q over the "
        "340 support/variance-selected traits that `select_stratum_traits` kept in that "
        "table (its `correction_family_size`, of 344 retained non-MediaDive traits in "
        "`retained_trait_count` before the `select_stratum_traits` support-bounds, "
        "positive-variance, and MediaDive-exclusion filter; the selection keeps traits "
        "annotated in at least 6 and at most 95% of profiles with positive composition "
        "variance)."
    ) in normalized


def test_selection_count_follows_exported_traits(monkeypatch: pytest.MonkeyPatch) -> None:
    supplied = supplied_stats(MIXED_P)
    chosen = (sankeys.SHARED_TRAITS[0], sankeys.SHARED_TRAITS[2], sankeys.SHARED_TRAITS[5])
    monkeypatch.setattr(sankeys, "SHARED_TRAITS", chosen)

    result = sankeys.build_trait_selection(supplied)

    assert result["trait_id"].tolist() == [trait.trait_id for trait in chosen]
    assert result["paired_p_value_type"].tolist() == ["nominal", "not_tested", "nominal"]
    assert result["displayed_nominal_test_count"].tolist() == [2, 2, 2]
    assert result["bh_family_size"].dropna().tolist() == [340] * len(chosen)
    assert result["bh_unfiltered_trait_count"].dropna().tolist() == [344] * len(chosen)
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
        assert saved["bh_unfiltered_trait_count"].tolist() == [344]
        assert saved["paired_q_value_family"].tolist() == [EXPECTED_FAMILY_340]
        expected = existing_selection_columns(supplied)
        expected = expected.loc[expected["trait_id"] == trait.trait_id].reset_index(drop=True)
        pd.testing.assert_frame_equal(saved[OLD_COLUMNS], expected)
    pd.testing.assert_frame_equal(selection, before)
