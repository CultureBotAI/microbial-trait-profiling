"""Null generation scenarios and RNG utilities for composition calibration (#263)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from mtp.calibration import _stream_key as stream_key
from mtp.calibration import replicate_rng
from mtp.composition_calibration.adapters.jed_composition import (
    PROFILE_NONTRAIT_COLUMNS,
    filter_paired_profiles,
)
from mtp.composition_denominators import retained_trait_counts


def study_flip_draw_order(profiles_df: pd.DataFrame) -> list[str]:
    """Order study IDs so non-zero-total paired studies draw first (#543).

    ``mtp.calibration`` and ``mtp.screen_calibration`` call ``drop_zero_total_pairs()``
    before drawing study flips, leaving the studies that have at least one paired
    experiment with positive retained non-MediaDive trait totals on both sides
    (939 studies on the canonical corpus). ``JedCompositionAdapter`` retains all
    959 paired studies in its corpus. Drawing in ``repr``-sorted order over the
    non-zero-total paired studies first, followed by any zero-total-only paired
    studies and any unpaired studies, guarantees that every study shared with
    ``calibration_flip_assignments.tsv`` receives the identical flip assignment
    in every replicate.
    """
    if profiles_df.empty or "study_id" not in profiles_df.columns:
        return []

    all_studies = set(profiles_df["study_id"].astype(str).unique())
    if "experiment_id" not in profiles_df.columns or "inc_dec" not in profiles_df.columns:
        return sorted(all_studies, key=repr)

    paired_df = filter_paired_profiles(profiles_df)
    if paired_df.empty:
        return sorted(all_studies, key=repr)

    retained = retained_trait_counts(
        paired_df.drop(columns=list(PROFILE_NONTRAIT_COLUMNS), errors="ignore")
        .apply(pd.to_numeric, errors="coerce")
        .fillna(0.0)
    )
    row_totals = retained.sum(axis=1)
    exp_min = (
        pd.DataFrame(
            {
                "study_id": paired_df["study_id"].astype(str).to_numpy(),
                "experiment_id": paired_df["experiment_id"].astype(str).to_numpy(),
                "total": row_totals.to_numpy(dtype=float),
            }
        )
        .groupby(["study_id", "experiment_id"], sort=False, as_index=False)["total"]
        .min()
    )
    nonzero_paired_studies = sorted(
        set(exp_min.loc[exp_min["total"] > 0, "study_id"].astype(str)),
        key=repr,
    )
    paired_studies = set(paired_df["study_id"].astype(str).unique())
    zero_paired_studies = sorted(paired_studies - set(nonzero_paired_studies), key=repr)
    unpaired_studies = sorted(all_studies - paired_studies, key=repr)
    return nonzero_paired_studies + zero_paired_studies + unpaired_studies


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
        keys: list[str] | list[tuple[str, str]] = study_flip_draw_order(profiles_df)
    else:
        keys = sorted(
            set(
                zip(
                    profiles_df["study_id"].astype(str),
                    profiles_df["experiment_id"].astype(str),
                    strict=True,
                )
            ),
            key=repr,
        )

    flipped_keys = {key for key in keys if rng.random() < 0.5}

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
    "study_flip_draw_order",
]
