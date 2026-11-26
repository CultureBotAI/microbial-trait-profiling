"""Guild trait z-scores from the full retained profile table (#2399, #2403).

A guild's declared traits are not guaranteed to be among the study heatmap's top
traits (METPO:1000603 "anaerobic" is above the heatmap's profile-support cap), so
guild scores are computed from the retained counts by trait id, never by looking
up heatmap matrix columns.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

import pandas as pd

from mtp.composition_denominators import profile_trait_z_scores

__all__ = ["guild_trait_z_scores", "read_guild_trait_selection"]


def read_guild_trait_selection(path: Path) -> dict[str, str]:
    """Return {trait_id: profile-matrix object name} from a guild_trait_selection.tsv."""
    selection = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if "trait_id" not in selection.columns:
        raise ValueError(f"{path} has no trait_id column")
    name_column = (
        "profile_matrix_object_name"
        if "profile_matrix_object_name" in selection.columns
        else "trait_name"
    )
    if name_column not in selection.columns:
        raise ValueError(f"{path} has neither profile_matrix_object_name nor trait_name")
    trait_ids = selection["trait_id"].str.strip()
    names = selection[name_column].str.strip()
    if trait_ids.eq("").any() or names.eq("").any():
        raise ValueError(f"{path} has a blank trait_id or trait name")
    if trait_ids.duplicated().any() or names.duplicated().any():
        raise ValueError(f"{path} lists a trait more than once")
    if selection.empty:
        raise ValueError(f"{path} selects no traits")
    return dict(zip(trait_ids, names, strict=True))


def guild_trait_z_scores(counts: pd.DataFrame, traits: Mapping[str, str]) -> pd.DataFrame:
    """Score each guild trait (trait_id -> column name) from the full profile table."""
    return profile_trait_z_scores(counts, list(traits)).rename(columns=dict(traits))
