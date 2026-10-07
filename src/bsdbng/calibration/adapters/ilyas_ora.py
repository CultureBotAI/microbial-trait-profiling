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
from bsdbng.enrichment_stats import benjamini_hochberg as benjamini_hochberg

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

        trait_columns = shared_trait_columns(profiles_df, taxa_matrix_df)
        m_universe = len(taxa_matrix_df)
        metadata: dict[str, Any] = {
            "test_type": test_type,
            "alpha": alpha,
            "test_scope": TEST_SCOPE,
            "m_universe": m_universe,
            "n_trait_columns": len(trait_columns),
            "n_tests": 0,
        }
        if not trait_columns or m_universe <= 0:
            return self._empty_result(metadata)

        signature_sizes = resolve_signature_sizes(profiles_df)
        profile_ids = profile_identifiers(profiles_df)
        k_pop_arr = trait_population_counts(taxa_matrix_df, trait_columns)
        observed = (
            profiles_df[trait_columns].apply(pd.to_numeric, errors="raise").to_numpy(dtype=np.int64)
        )

        test_mask = (observed > 0) & (k_pop_arr[np.newaxis, :] > 0)
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
                "is_significant": q_values <= alpha,
            }
        )
        return MethodInferenceResult(
            method_id=self.method_id,
            owner=self.owner,
            target_null=self.target_null,
            results_df=res_df,
            metadata=metadata,
        )
