"""Declared exact-stool population and support-first group selection (#1093).

The caller computes native grouped contrast cells from exact-stool metadata.
Selection uses the existing support flags and unique-study counts, never p/q
values. Conversion to percentage points does not change the native inference.
"""

from __future__ import annotations

import math
from numbers import Integral, Real

import numpy as np
import pandas as pd

DELTA_FRACTION_COLUMN = "increased_minus_decreased_positive_fraction"
DELTA_PERCENTAGE_POINTS_COLUMN = "increased_minus_decreased_positive_percentage_points"
SELECTION_COLUMNS = ("condition", "stool_study_count", "eligible", "rank", "chosen")
CELL_COLUMNS = (
    "view_key",
    "body_site_group",
    "condition",
    "trait_id",
    "study_count",
    "passes_support",
    DELTA_FRACTION_COLUMN,
)


def _require_columns(frame: pd.DataFrame, required: tuple[str, ...]) -> None:
    if not frame.columns.is_unique:
        raise ValueError("Duplicate table columns are not allowed")
    missing = set(required) - set(frame.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")


def exact_stool_profiles(metadata: pd.DataFrame) -> pd.DataFrame:
    """Keep only literal ``body_site == 'Feces'`` profiles, preserving rows/index.

    Grouped anatomy labels, multi-site values, spelling variants and whitespace
    are not substitutes for the declared metadata population.
    """
    _require_columns(metadata, ("body_site",))
    return metadata.loc[metadata["body_site"].eq("Feces")].copy()


def _percentage_points(value: object) -> float:
    if value is None or value is pd.NA:
        return math.nan
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(
            f"{DELTA_FRACTION_COLUMN} must contain numeric fractions or missing values"
        )
    numeric = float(value)
    if not math.isnan(numeric) and (not math.isfinite(numeric) or abs(numeric) > 1):
        raise ValueError(f"{DELTA_FRACTION_COLUMN} must be in [-1, 1] or missing")
    return 100 * numeric


def select_stool_group_contrasts(
    cells: pd.DataFrame, *, group_count: int = 6
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Keep all cells of the top supported groups by unique stool-study count.

    Input must contain native grouped-direction cells computed from exact-stool
    profiles. A group is eligible when at least one cell passes native support.
    Rank eligible groups by study count descending, then exact condition name.
    The selection ledger includes ineligible groups with missing rank. Selected
    cells retain their original row order, index, support flags and inference.
    """
    if isinstance(group_count, (bool, np.bool_)) or not isinstance(group_count, Integral):
        raise ValueError("group_count must be a positive nonboolean integer")
    if group_count <= 0:
        raise ValueError("group_count must be a positive nonboolean integer")
    _require_columns(cells, CELL_COLUMNS)
    if DELTA_PERCENTAGE_POINTS_COLUMN in cells.columns:
        raise ValueError(f"Refusing to overwrite existing {DELTA_PERCENTAGE_POINTS_COLUMN}")
    for column, expected in (
        ("view_key", "grouped_direction_contrast"),
        ("body_site_group", "Feces"),
    ):
        if any(not isinstance(value, str) or value != expected for value in cells[column]):
            raise ValueError(f"{column} must contain only {expected!r}")
    for column in ("condition", "trait_id"):
        if any(not isinstance(value, str) or not value.strip() for value in cells[column]):
            raise ValueError(f"{column} must contain nonblank identifiers")
    if cells.duplicated(["condition", "trait_id"]).any():
        raise ValueError("Duplicate condition/trait_id cells are not allowed")
    for value in cells["study_count"]:
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral) or value <= 0:
            raise ValueError("study_count must contain positive nonboolean integers")
    if any(not isinstance(value, (bool, np.bool_)) for value in cells["passes_support"]):
        raise ValueError("passes_support must contain boolean flags")
    # Validate all supplied effects, including cells in groups not selected.
    percentages = cells[DELTA_FRACTION_COLUMN].map(_percentage_points)
    records: list[dict[str, object]] = []
    for condition, group in cells.groupby("condition", sort=False):
        if group["study_count"].nunique() != 1:
            raise ValueError(f"Inconsistent study_count across traits in group {condition!r}")
        records.append(
            {
                "condition": condition,
                "stool_study_count": int(group["study_count"].iloc[0]),
                "eligible": bool(group["passes_support"].any()),
            }
        )
    selection = pd.DataFrame.from_records(records, columns=SELECTION_COLUMNS[:3])
    selection = selection.sort_values(
        ["eligible", "stool_study_count", "condition"], ascending=[False, False, True]
    ).reset_index(drop=True)
    selection["eligible"] = selection["eligible"].astype(bool)
    selection["rank"] = pd.array(
        [index + 1 if eligible else None for index, eligible in enumerate(selection["eligible"])],
        dtype="Int64",
    )
    selection["chosen"] = selection["rank"].le(group_count).fillna(False).astype(bool)
    chosen = selection.loc[selection["chosen"], "condition"]
    selected_rows = cells["condition"].isin(chosen)
    selected = cells.loc[selected_rows].copy()
    selected[DELTA_PERCENTAGE_POINTS_COLUMN] = percentages.loc[selected_rows].to_numpy()
    return selected, selection
