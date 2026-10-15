"""Composition denominators shared with the non-MediaDive trait-fraction arm.

Normalizes trait profiles over all retained non-MediaDive trait columns before
pooled/stratified trait selection and transformation. Zero-signal and
MediaDive-only rows remain in the corpus with all-zero fractions.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from mtp.trait_guilds import MEDIADIVE_MARKER

DENOMINATOR_DEFINITION = "sum of all retained non-MediaDive trait columns before selection"

__all__ = [
    "DENOMINATOR_DEFINITION",
    "load_retained_profiles",
    "retained_trait_counts",
    "trait_fractions",
    "validate_cached_totals",
]


def retained_trait_counts(counts: pd.DataFrame) -> pd.DataFrame:
    """Keep the full biological trait universe, independent of the tested subset."""
    columns = [column for column in counts.columns if MEDIADIVE_MARKER not in str(column).lower()]
    return counts.loc[:, columns]


def trait_fractions(counts: pd.DataFrame) -> pd.DataFrame:
    """Normalize retained counts; a zero retained denominator yields all-zero fractions.

    Zero-signal and MediaDive-only rows stay in the corpus. Their log1p values
    are zero and are included in column standardization like other rows; their
    final z-scores need not be zero. No biological signal is inferred from media.
    """
    retained = retained_trait_counts(counts)
    totals = retained.sum(axis=1)
    return retained.div(totals.where(totals.ne(0)), axis=0).fillna(0.0)


def validate_cached_totals(counts: pd.DataFrame, metadata: pd.DataFrame) -> None:
    """Reject stale caches; counts must include every retained input trait."""
    if not counts.index.is_unique or not metadata.index.is_unique:
        raise ValueError("Profile keys must be unique before checking cached totals")
    if not metadata.index.isin(counts.index).all():
        raise ValueError("Cached metadata contains profile keys missing from the raw matrix")
    expected = retained_trait_counts(counts).sum(axis=1).reindex(metadata.index)
    cached = pd.to_numeric(metadata["total_trait_signal"], errors="raise")
    mismatch = expected.ne(cached) | cached.isna()
    if mismatch.any():
        examples = ", ".join(map(str, metadata.index[mismatch][:3]))
        raise ValueError(
            f"Stale or incompatible non-MediaDive denominator cache: {examples}. "
            "Regenerate profile metadata from the same full raw matrix."
        )


def load_retained_profiles(path: Path, metadata: pd.DataFrame) -> pd.DataFrame:
    """Read the full denominator universe and align by keys, never by row position."""
    frame = pd.read_csv(path, sep="\t")
    keys = ["study_id", "experiment_id", "inc_dec"]
    counts = frame.drop(columns=keys).apply(pd.to_numeric, errors="raise")
    counts.index = frame[keys].astype(str).agg("|".join, axis=1)
    counts = retained_trait_counts(counts.groupby(level=0, sort=False).sum())
    validate_cached_totals(counts, metadata)
    return counts.loc[metadata.index]
