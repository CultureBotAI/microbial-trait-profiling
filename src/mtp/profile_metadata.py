"""Profile-keyed metadata frames with a ``profile_id`` column (#2417).

``trait-visualization/study_trait_profile_study_metadata.tsv`` is read with its first
column as the index. Files written before pandas 3.1 have a blank first header (#2396);
files written under pandas 3.1 name it ``profile_id``. Copying such an index into a
``profile_id`` column leaves a frame whose index level and column share a name, which
pandas refuses in merge, sort_values, groupby, reset_index and join. Readers that need
the column call :func:`with_profile_id_column`, which works for both header forms.
:func:`with_index_column` does the same for other keys (``trait_id`` in the trait metadata).
"""

from __future__ import annotations

import pandas as pd

__all__ = ["with_index_column", "with_profile_id_column"]


def with_index_column(frame: pd.DataFrame, column: str) -> pd.DataFrame:
    """Return a copy with an unnamed index and that index copied into ``column``."""
    result = frame.rename_axis(index=None)
    result[column] = result.index
    return result


def with_profile_id_column(frame: pd.DataFrame) -> pd.DataFrame:
    """Return a copy with an unnamed index and that index as the ``profile_id`` column."""
    return with_index_column(frame, "profile_id")
