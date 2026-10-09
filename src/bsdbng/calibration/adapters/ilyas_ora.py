"""Ilyas ORA & Depletion framework adapter for the calibration harness."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

from bsdbng.calibration.protocol import MethodInferenceResult
from bsdbng.calibration.scenarios import (
    profile_identifiers,
    resolve_signature_sizes,
    shared_trait_columns,
    trait_population_counts,
)
from bsdbng.enrichment_stats import benjamini_hochberg
from bsdbng.study_trait_profiles import direction_to_inc_dec

RESULT_COLUMNS: list[str] = [
    "profile_id",
    "trait",
    "k_observed",
    "k_pop",
    "m_universe",
    "n_draws",
    "effect_size",
    "p_value",
    "q_value",
    "is_significant",
]

# Which (profile, trait) pairs enter the BH family. Recorded in metadata so that
# calibration reports state the test scope explicitly (see #263).
TEST_SCOPE = "pairs_with_k_observed_gt_0_and_k_pop_gt_0_and_n_draws_gt_0"


class IlyasORAAdapter:
    """Method adapter for hypergeometric Over-Representation Analysis (ORA) & Depletion."""

    method_id: str = "ilyas_ora_hypergeometric"
    owner: str = "Ilyas"
    target_null: str = (
        "Taxa in signature drawn uniformly from background universe without replacement"
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
        """Run ORA hypergeometric enrichment or depletion inference."""
        config = config or {}
        test_type = str(config.get("test_type", "enrichment"))
        if test_type not in {"enrichment", "depletion"}:
            raise ValueError(f"Unsupported test_type: {test_type}")
        alpha = float(config.get("alpha", 0.05))
        min_fold_raw = config.get("min_fold", 1.5)
        min_fold = float(min_fold_raw) if min_fold_raw is not None else 1.5

        trait_columns = shared_trait_columns(profiles_df, taxa_matrix_df)
        m_universe = len(taxa_matrix_df)
        metadata: dict[str, Any] = {
            "test_type": test_type,
            "alpha": alpha,
            "min_fold": min_fold,
            "test_scope": (
                TEST_SCOPE
                if test_type == "enrichment"
                else "pairs_with_k_pop_gt_0_and_n_draws_gt_0"
            ),
            "m_universe": m_universe,
            "n_trait_columns": len(trait_columns),
            "n_tests": 0,
        }
        if not trait_columns or m_universe <= 0:
            return self._empty_result(metadata)

        # 1. Production parity: filter to paired experiments only (having both inc and dec)
        if {"study_id", "experiment_id", "inc_dec"}.issubset(profiles_df.columns):
            norm_dirs = profiles_df["inc_dec"].astype(str).map(direction_to_inc_dec)
            valid_mask = norm_dirs.isin({"inc", "dec"})
            if not valid_mask.all():
                profiles_df = profiles_df[valid_mask].copy()
                norm_dirs = norm_dirs[valid_mask]

            exp_groups = norm_dirs.groupby(
                [profiles_df["study_id"], profiles_df["experiment_id"]]
            ).unique()
            paired_exps = {
                idx for idx, dirs in exp_groups.items() if "inc" in dirs and "dec" in dirs
            }
            exp_keys = list(zip(profiles_df["study_id"], profiles_df["experiment_id"], strict=True))
            paired_mask = [k in paired_exps for k in exp_keys]
            profiles_df = profiles_df[paired_mask].copy()

        if profiles_df.empty:
            return self._empty_result(metadata)

        # 2. Extract sizes, IDs, population counts, and observed counts
        signature_sizes = resolve_signature_sizes(profiles_df)
        profile_ids = profile_identifiers(profiles_df)
        k_pop_arr = trait_population_counts(taxa_matrix_df, trait_columns)
        observed = (
            profiles_df[trait_columns].apply(pd.to_numeric, errors="raise").to_numpy(dtype=np.int64)
        )
        if np.any(observed < 0):
            raise ValueError("Observed trait counts cannot be negative")

        # 3. Validate support constraints: M >= max(K, N) and N >= max(k_obs) and k_obs <= K
        max_k_pop = int(k_pop_arr.max()) if len(k_pop_arr) > 0 else 0
        max_n = int(signature_sizes.max()) if len(signature_sizes) > 0 else 0
        if m_universe < max(max_k_pop, max_n):
            raise ValueError(
                f"Background universe M ({m_universe}) is smaller than "
                f"maximum trait count K ({max_k_pop}) or maximum draw N ({max_n})"
            )
        if observed.size > 0:
            max_k_obs_per_row = observed.max(axis=1)
            if np.any(signature_sizes < max_k_obs_per_row):
                raise ValueError("Observed trait count k exceeds profile draw count N")
            if np.any(observed > k_pop_arr[np.newaxis, :]):
                raise ValueError("Observed trait count k exceeds population trait count K")

        if test_type == "enrichment":
            test_mask = (observed > 0) & (k_pop_arr[np.newaxis, :] > 0)
        else:
            test_mask = k_pop_arr[np.newaxis, :] > 0
        test_mask &= signature_sizes[:, np.newaxis] > 0
        if not np.any(test_mask):
            return self._empty_result(metadata)

        p_idx, t_idx = np.nonzero(test_mask)
        k_obs = observed[p_idx, t_idx]
        k_pop = k_pop_arr[t_idx]
        n_draws = signature_sizes[p_idx]

        if test_type == "enrichment":
            p_vals = hypergeom.sf(k_obs - 1, m_universe, k_pop, n_draws)
        else:
            p_vals = hypergeom.cdf(k_obs, m_universe, k_pop, n_draws)
        p_vals = np.asarray(p_vals, dtype=float)

        expected = n_draws * (k_pop / m_universe)
        effect_sizes = np.divide(
            k_obs, expected, out=np.zeros_like(expected, dtype=float), where=expected > 0
        )
        q_values = benjamini_hochberg(p_vals)
        metadata["n_tests"] = len(p_vals)

        if test_type == "enrichment":
            is_significant = (q_values <= alpha) & (effect_sizes >= min_fold)
        else:
            is_significant = (
                (q_values <= alpha) & (effect_sizes <= (1.0 / min_fold))
                if min_fold > 0
                else (q_values <= alpha)
            )

        res_df = pd.DataFrame(
            {
                "profile_id": profile_ids[p_idx],
                "trait": np.asarray(trait_columns, dtype=object)[t_idx],
                "k_observed": k_obs,
                "k_pop": k_pop,
                "m_universe": m_universe,
                "n_draws": n_draws,
                "effect_size": effect_sizes,
                "p_value": p_vals,
                "q_value": q_values,
                "is_significant": is_significant,
            }
        )
        return MethodInferenceResult(
            method_id=self.method_id,
            owner=self.owner,
            target_null=self.target_null,
            results_df=res_df,
            metadata=metadata,
        )
