"""Jed's paired log-ratio composition adapter for the calibration harness (#263).

Implements the paired composition inference procedure evaluated across studies:
1. Normalizes trait counts to profile fractions: f_ij = X_ij / sum_k X_ik.
2. Applies log-scale transformation: y_ij = log1p(1000 * f_ij).
3. Standardizes columns to z-scores across profiles: z_ij = (y_ij - mean_j) / std_j.
4. Computes matched increased-minus-decreased pair differences: Delta_pj = z_inc,j - z_dec,j.
5. Aggregates by study: mean_Delta_sj = mean_{p in s} Delta_pj.
6. Performs one-sample Student's t-test against zero across studies for each trait.
7. Corrects for multiple testing across traits via Benjamini-Hochberg (BH) or BY.
"""

from __future__ import annotations

import logging
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any, Literal

import numpy as np
import pandas as pd
from scipy import stats

from mtp.composition_calibration.protocol import MethodInferenceResult
from mtp.composition_denominators import retained_trait_counts, trait_fractions
from mtp.enrichment_stats import benjamini_hochberg
from mtp.ora_calibration.scenarios import shared_trait_columns
from mtp.study_trait_profiles import direction_to_inc_dec

logger = logging.getLogger(__name__)

RESULT_COLUMNS: list[str] = [
    "profile_id",
    "trait",
    "effect_size",
    "t_statistic",
    "p_value",
    "q_value",
    "is_significant",
    "mean_contrast",
    "median_contrast",
    "std_contrast",
    "ci95_lower",
    "ci95_upper",
    "n_studies",
    "n_pairs",
]

DEFAULT_PSEUDOCOUNT_SCALE = 1000.0


@dataclass(frozen=True)
class CompositionContrastData:
    """Precomputed pair-level and study-level contrasts for composition inference."""

    pair_contrasts: pd.DataFrame
    study_contrasts: pd.DataFrame
    n_pairs: int
    n_studies: int
    trait_columns: list[str]


def filter_paired_profiles(profiles_df: pd.DataFrame) -> pd.DataFrame:
    """Filter to experiments having both increased and decreased profiles within a study."""
    if profiles_df.empty:
        return profiles_df.copy()

    df = profiles_df.copy()
    if "_norm_dir" not in df.columns:
        df["_norm_dir"] = df["inc_dec"].astype(str).map(direction_to_inc_dec)

    # Keep only profiles mapped to valid 'inc' or 'dec'
    valid_mask = df["_norm_dir"].isin(["inc", "dec"])
    df = df[valid_mask].copy()

    # Identify experiments with both directions
    dir_counts = (
        df.groupby(["study_id", "experiment_id"])["_norm_dir"].nunique().reset_index(name="n_dirs")
    )
    paired_keys = dir_counts[dir_counts["n_dirs"] >= 2][["study_id", "experiment_id"]]

    paired_merged = df.merge(paired_keys, on=["study_id", "experiment_id"], how="inner")
    return paired_merged


def compute_profile_zscores(
    profiles_df: pd.DataFrame,
    trait_columns: Sequence[str],
    pseudocount_scale: float = DEFAULT_PSEUDOCOUNT_SCALE,
) -> pd.DataFrame:
    """Normalize counts to profile fractions, apply log1p, and standardize columns."""
    present_traits = [c for c in trait_columns if c in profiles_df.columns]
    if profiles_df.empty or not present_traits:
        return pd.DataFrame(columns=list(trait_columns), index=profiles_df.index)

    retained = retained_trait_counts(
        profiles_df[present_traits].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    )
    traits = list(retained.columns)
    if not traits:
        return pd.DataFrame(columns=[], index=profiles_df.index)

    fractions = trait_fractions(retained).to_numpy(dtype=float)

    # Log1p transformation scaled by pseudocount multiplier
    log_fractions = np.log1p(pseudocount_scale * fractions)

    # Standardize columns to z-scores across all profiles
    means = np.mean(log_fractions, axis=0, keepdims=True)
    stds = np.std(log_fractions, axis=0, keepdims=True, ddof=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        z_array = np.where(stds > 1e-12, (log_fractions - means) / stds, 0.0)

    zscores = pd.DataFrame(z_array, columns=traits, index=profiles_df.index)
    return zscores


def compute_composition_contrasts(
    profiles_df: pd.DataFrame,
    trait_columns: Sequence[str],
    pseudocount_scale: float = DEFAULT_PSEUDOCOUNT_SCALE,
) -> CompositionContrastData:
    """Compute matched experiment pair differences and study mean contrasts."""
    paired_df = filter_paired_profiles(profiles_df)
    present_traits = [c for c in trait_columns if c in paired_df.columns]
    if paired_df.empty or not present_traits:
        empty_df = pd.DataFrame(columns=list(trait_columns))
        return CompositionContrastData(
            pair_contrasts=empty_df,
            study_contrasts=empty_df,
            n_pairs=0,
            n_studies=0,
            trait_columns=list(trait_columns),
        )

    zscores = compute_profile_zscores(paired_df, present_traits, pseudocount_scale)
    traits = list(zscores.columns)
    if not traits:
        empty_df = pd.DataFrame(columns=[])
        return CompositionContrastData(
            pair_contrasts=empty_df,
            study_contrasts=empty_df,
            n_pairs=0,
            n_studies=0,
            trait_columns=[],
        )

    norm_dir = (
        paired_df["_norm_dir"]
        if "_norm_dir" in paired_df.columns
        else paired_df["inc_dec"].astype(str).map(direction_to_inc_dec)
    )
    study_col = paired_df["study_id"].astype(str)
    exp_col = paired_df["experiment_id"].astype(str)

    zscores["_study_id"] = study_col.values
    zscores["_experiment_id"] = exp_col.values
    zscores["_norm_dir"] = norm_dir.values

    inc_mask = zscores["_norm_dir"] == "inc"
    dec_mask = zscores["_norm_dir"] == "dec"

    inc_pairs = zscores[inc_mask].groupby(["_study_id", "_experiment_id"])[traits].mean()
    dec_pairs = zscores[dec_mask].groupby(["_study_id", "_experiment_id"])[traits].mean()

    common_pairs = inc_pairs.index.intersection(dec_pairs.index)
    pair_diff = inc_pairs.loc[common_pairs] - dec_pairs.loc[common_pairs]
    pair_diff.index.names = ["study_id", "experiment_id"]

    study_diff = pair_diff.groupby(level="study_id").mean()

    return CompositionContrastData(
        pair_contrasts=pair_diff,
        study_contrasts=study_diff,
        n_pairs=len(pair_diff),
        n_studies=len(study_diff),
        trait_columns=traits,
    )


class JedCompositionAdapter:
    """Method adapter for paired log-ratio trait composition inference (#263)."""

    method_id: str = "jed_composition_ttest"
    owner: str = "Jed"
    target_null: str = (
        "Study-level mean matched profile log-fraction contrast is zero across studies"
    )

    def _empty_result(self, metadata: dict[str, Any]) -> MethodInferenceResult:
        return MethodInferenceResult(
            method_id=self.method_id,
            owner=self.owner,
            target_null=self.target_null,
            results_df=pd.DataFrame(columns=RESULT_COLUMNS),
            metadata=metadata,
        )

    def run_inference(
        self,
        profiles_df: pd.DataFrame,
        taxa_matrix_df: pd.DataFrame,
        config: dict[str, Any] | None = None,
    ) -> MethodInferenceResult:
        """Run paired composition inference procedure on provided profile data."""
        config = config or {}
        alpha = float(config.get("alpha", 0.05))
        correction: Literal["bh", "by"] = config.get("correction", "bh")
        min_fold_raw = config.get("min_fold", None)
        min_contrast = float(min_fold_raw) if min_fold_raw is not None else None
        pseudocount_scale = float(config.get("pseudocount_scale", DEFAULT_PSEUDOCOUNT_SCALE))
        use_study_aggregation = bool(config.get("study_aggregation", True))

        # Determine trait universe: either specified in config or shared traits
        if config.get("traits"):
            candidate_columns = [c for c in config["traits"] if c in profiles_df.columns]
        elif not taxa_matrix_df.empty:
            candidate_columns = shared_trait_columns(profiles_df, taxa_matrix_df)
        else:
            meta_cols = {
                "study_id",
                "experiment_id",
                "inc_dec",
                "condition",
                "body_site",
                "taxa_count",
                "signature_size",
                "_norm_dir",
            }
            candidate_columns = [c for c in profiles_df.columns if c not in meta_cols]
        trait_columns = (
            retained_trait_counts(profiles_df[candidate_columns]).columns.tolist()
            if candidate_columns and not profiles_df.empty
            else [c for c in candidate_columns if "mediadive.medium:" not in str(c).lower()]
        )

        metadata: dict[str, Any] = {
            "method_id": self.method_id,
            "owner": self.owner,
            "alpha": alpha,
            "correction": correction,
            "min_contrast": min_contrast,
            "pseudocount_scale": pseudocount_scale,
            "study_aggregation": use_study_aggregation,
            "n_trait_columns": len(trait_columns),
            "n_tests": 0,
            "n_studies": 0,
            "n_pairs": 0,
        }

        if not trait_columns or profiles_df.empty:
            return self._empty_result(metadata)

        contrast_data = compute_composition_contrasts(
            profiles_df, trait_columns, pseudocount_scale=pseudocount_scale
        )
        metadata["n_studies"] = contrast_data.n_studies
        metadata["n_pairs"] = contrast_data.n_pairs

        # Decide testing matrix: study-aggregated or pair-level
        test_matrix = (
            contrast_data.study_contrasts if use_study_aggregation else contrast_data.pair_contrasts
        )
        n_obs = len(test_matrix)
        if n_obs < 2:
            return self._empty_result(metadata)

        records: list[dict[str, Any]] = []
        p_values: list[float] = []

        for trait in trait_columns:
            vals = test_matrix[trait].dropna().to_numpy(dtype=float)
            k_obs = len(vals)
            mean_c = float(np.mean(vals)) if k_obs > 0 else 0.0
            med_c = float(np.median(vals)) if k_obs > 0 else 0.0
            std_c = float(np.std(vals, ddof=1)) if k_obs >= 2 else 0.0

            if k_obs >= 2 and std_c > 1e-12:
                t_res = stats.ttest_1samp(vals, 0.0)
                t_stat = float(t_res.statistic)
                p_val = float(t_res.pvalue)
                se = std_c / np.sqrt(k_obs)
                t_crit = float(stats.t.ppf(1.0 - alpha / 2.0, k_obs - 1))
                ci_low = mean_c - t_crit * se
                ci_high = mean_c + t_crit * se
            else:
                t_stat = 0.0
                p_val = 1.0
                ci_low = float("nan")
                ci_high = float("nan")

            p_values.append(p_val)
            records.append(
                {
                    "profile_id": "pooled",
                    "trait": trait,
                    "effect_size": mean_c,
                    "t_statistic": t_stat,
                    "p_value": p_val,
                    "mean_contrast": mean_c,
                    "median_contrast": med_c,
                    "std_contrast": std_c,
                    "ci95_lower": ci_low,
                    "ci95_upper": ci_high,
                    "n_studies": contrast_data.n_studies,
                    "n_pairs": contrast_data.n_pairs,
                }
            )

        p_array = np.asarray(p_values, dtype=float)
        if correction == "by":
            # Benjamini-Yekutieli correction
            m = len(p_array)
            harmonic_sum = np.sum(1.0 / np.arange(1, m + 1))
            q_array = benjamini_hochberg(p_array) * harmonic_sum
            q_array = np.clip(q_array, 0.0, 1.0)
        else:
            q_array = benjamini_hochberg(p_array)

        metadata["n_tests"] = len(p_array)

        is_sig = q_array <= alpha
        if min_contrast is not None:
            effects = np.array([r["effect_size"] for r in records])
            is_sig = is_sig & (np.abs(effects) >= min_contrast)

        for i, rec in enumerate(records):
            rec["q_value"] = float(q_array[i])
            rec["is_significant"] = bool(is_sig[i])

        res_df = pd.DataFrame(records)[RESULT_COLUMNS]
        metadata["n_significant"] = int(np.sum(is_sig))

        return MethodInferenceResult(
            method_id=self.method_id,
            owner=self.owner,
            target_null=self.target_null,
            results_df=res_df,
            metadata=metadata,
        )


__all__ = [
    "DEFAULT_PSEUDOCOUNT_SCALE",
    "RESULT_COLUMNS",
    "CompositionContrastData",
    "JedCompositionAdapter",
    "compute_composition_contrasts",
    "compute_profile_zscores",
    "filter_paired_profiles",
]
