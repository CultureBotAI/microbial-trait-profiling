"""Null generation scenarios and RNG utilities for composition calibration (#263)."""

from __future__ import annotations

import hashlib

import numpy as np
import pandas as pd


def stream_key(name: str) -> int:
    """A stable integer stream ID for a named scenario, independent of its position."""
    return int.from_bytes(hashlib.sha256(name.encode("utf-8")).digest()[:8], "big")


def replicate_rng(seed: int, stream: str, replicate: int) -> np.random.Generator:
    """The generator for one replicate of one named stream.

    numpy.random.default_rng takes the whole triple as seed entropy, so every
    (seed, stream, replicate) is independently reproducible across arms.
    """
    return np.random.default_rng([seed, stream_key(stream), replicate])


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
