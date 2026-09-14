"""Composition denominators shared with the non-MediaDive trait-fraction arm."""

from __future__ import annotations

import pandas as pd

from bsdbng.trait_guilds import MEDIADIVE_MARKER

DENOMINATOR_DEFINITION = "sum of all retained non-MediaDive trait columns before selection"


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
