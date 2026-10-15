"""Unit tests for composition_denominators module (Refs #174, #202).

Verifies non-MediaDive denominator normalization, MediaDive column exclusion,
zero-signal handling, cache validation, and key alignment.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from mtp.composition_denominators import (
    DENOMINATOR_DEFINITION,
    load_retained_profiles,
    retained_trait_counts,
    trait_fractions,
    validate_cached_totals,
)


def test_denominator_definition_constant() -> None:
    assert "retained non-MediaDive" in DENOMINATOR_DEFINITION


def test_retained_trait_counts_filters_mediadive() -> None:
    counts = pd.DataFrame(
        {
            "METPO:1000615": [1, 2],
            "METPO:2000517|mediadive.medium:1": [10, 20],
            "METPO:2000517|MediaDive.Medium:2": [30, 40],
            "METPO:1000616": [5, 6],
        },
        index=["p1", "p2"],
    )
    retained = retained_trait_counts(counts)
    assert list(retained.columns) == ["METPO:1000615", "METPO:1000616"]
    pd.testing.assert_frame_equal(
        retained,
        pd.DataFrame(
            {"METPO:1000615": [1, 2], "METPO:1000616": [5, 6]},
            index=["p1", "p2"],
        ),
    )


def test_trait_fractions_normalization() -> None:
    counts = pd.DataFrame(
        {
            "A": [10, 0, 0, 20],
            "B": [90, 0, 0, 80],
            "mediadive.medium:1": [900, 0, 700, 5],
        },
        index=["p1", "p2", "p3", "p4"],
    )
    fractions = trait_fractions(counts)
    assert list(fractions.columns) == ["A", "B"]
    assert fractions["A"].tolist() == pytest.approx([0.1, 0.0, 0.0, 0.2])
    assert fractions["B"].tolist() == pytest.approx([0.9, 0.0, 0.0, 0.8])
    assert fractions.iloc[1].eq(0.0).all()
    assert fractions.iloc[2].eq(0.0).all()
    assert not fractions.isna().any().any()


def test_trait_fractions_all_mediadive_matrix() -> None:
    only_media = pd.DataFrame(
        {"mediadive.medium:1": [10, 0], "mediadive.medium:2": [5, 20]},
        index=["p1", "p2"],
    )
    assert retained_trait_counts(only_media).shape == (2, 0)
    fractions = trait_fractions(only_media)
    assert fractions.shape == (2, 0)


def test_validate_cached_totals_success() -> None:
    counts = pd.DataFrame(
        {"A": [5, 10], "B": [15, 20], "mediadive.medium:1": [100, 200]},
        index=["s1|e1|inc", "s2|e2|dec"],
    )
    metadata = pd.DataFrame(
        {"total_trait_signal": [20, 30]},
        index=["s1|e1|inc", "s2|e2|dec"],
    )
    # Should pass without raising
    validate_cached_totals(counts, metadata)


def test_validate_cached_totals_detects_stale_cache() -> None:
    counts = pd.DataFrame(
        {"A": [5, 10], "B": [15, 20], "mediadive.medium:1": [100, 200]},
        index=["s1|e1|inc", "s2|e2|dec"],
    )
    stale_metadata = pd.DataFrame(
        {"total_trait_signal": [120, 230]},  # Includes MediaDive, stale!
        index=["s1|e1|inc", "s2|e2|dec"],
    )
    with pytest.raises(ValueError, match="Stale or incompatible non-MediaDive denominator cache"):
        validate_cached_totals(counts, stale_metadata)


def test_validate_cached_totals_rejects_duplicate_index() -> None:
    dup_counts = pd.DataFrame({"A": [1, 2]}, index=["dup", "dup"])
    metadata = pd.DataFrame({"total_trait_signal": [1, 2]}, index=["dup", "dup"])
    with pytest.raises(ValueError, match="Profile keys must be unique"):
        validate_cached_totals(dup_counts, metadata)


def test_validate_cached_totals_rejects_missing_keys() -> None:
    counts = pd.DataFrame({"A": [1, 2]}, index=["k1", "k2"])
    metadata = pd.DataFrame({"total_trait_signal": [1, 2]}, index=["k1", "missing_k"])
    with pytest.raises(ValueError, match="Cached metadata contains profile keys missing"):
        validate_cached_totals(counts, metadata)


def test_load_retained_profiles(tmp_path: Path) -> None:
    data = {
        "study_id": ["s1", "s2", "s3"],
        "experiment_id": ["e1", "e2", "e3"],
        "inc_dec": ["inc", "dec", "inc"],
        "trait_x": [10, 20, 0],
        "trait_y": [30, 40, 0],
        "mediadive.medium:1": [999, 888, 777],
    }
    df = pd.DataFrame(data)
    tsv_path = tmp_path / "profiles.tsv"
    df.to_csv(tsv_path, sep="\t", index=False)

    metadata = pd.DataFrame(
        {"total_trait_signal": [40, 60]},
        index=["s1|e1|inc", "s2|e2|dec"],
    )

    loaded = load_retained_profiles(tsv_path, metadata)
    assert list(loaded.columns) == ["trait_x", "trait_y"]
    assert loaded.index.tolist() == ["s1|e1|inc", "s2|e2|dec"]
    assert loaded.loc["s1|e1|inc", "trait_x"] == 10
    assert loaded.loc["s2|e2|dec", "trait_y"] == 40
