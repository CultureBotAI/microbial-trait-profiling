"""Null generation scenarios and RNG utilities for composition calibration (#263)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from mtp.calibration import _stream_key as stream_key
from mtp.calibration import replicate_rng


def generate_study_flip_null_profiles(
    profiles_df: pd.DataFrame,
    rng: np.random.Generator,
    *,
    unit: str = "study",
) -> pd.DataFrame:
    """Exchange increased and decreased direction labels for a random half of units.

    unit is 'study' or 'experiment'. Every profile of a flipped unit is relabelled,
    so the corpus stays paired and each unit keeps its own counts.
    """
    if unit not in ("study", "experiment"):
        raise ValueError(f"Unknown flip unit '{unit}'; expected 'study' or 'experiment'")

    if profiles_df.empty:
        return profiles_df.copy()

    if unit == "study":
        keys = sorted(profiles_df["study_id"].astype(str).unique(), key=repr)
    else:
        keys = sorted(
            zip(
                profiles_df["study_id"].astype(str),
                profiles_df["experiment_id"].astype(str),
                strict=True,
            ),
            key=repr,
        )

    draws = rng.random(len(keys))
    flipped_keys = {key for key, d in zip(keys, draws, strict=True) if d < 0.5}

    out_df = profiles_df.copy()
    if unit == "study":
        mask = out_df["study_id"].astype(str).isin(flipped_keys)
    else:
        tuples = list(
            zip(
                out_df["study_id"].astype(str),
                out_df["experiment_id"].astype(str),
                strict=True,
            )
        )
        mask = pd.Series([t in flipped_keys for t in tuples], index=out_df.index)

    swap_map = {
        "increased": "decreased",
        "decreased": "increased",
        "inc": "dec",
        "dec": "inc",
    }
    out_df.loc[mask, "inc_dec"] = out_df.loc[mask, "inc_dec"].map(
        lambda x: swap_map.get(str(x), str(x))
    )
    if "_norm_dir" in out_df.columns:
        norm_swap = {"inc": "dec", "dec": "inc"}
        out_df.loc[mask, "_norm_dir"] = out_df.loc[mask, "_norm_dir"].map(
            lambda x: norm_swap.get(str(x), str(x))
        )

    return out_df


__all__ = [
    "generate_study_flip_null_profiles",
    "replicate_rng",
    "stream_key",
]
