"""Ilyas ORA & Depletion framework adapter for the calibration harness."""

from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

from bsdbng.calibration.protocol import MethodInferenceResult


class IlyasORAAdapter:
    """Method adapter for hypergeometric Over-Representation Analysis (ORA) & Depletion."""

    method_id: str = "ilyas_ora_hypergeometric"
    owner: str = "Ilyas"
    target_null: str = (
        "Taxa in signature drawn uniformly from background universe without replacement"
    )

    def run_inference(
        self,
        profiles_df: pd.DataFrame,
        taxa_matrix_df: pd.DataFrame,
        config: dict[str, Any] | None = None,
    ) -> MethodInferenceResult:
        """Run ORA hypergeometric enrichment and depletion inference."""
        config = config or {}
        test_type = config.get("test_type", "enrichment")  # "enrichment" or "depletion"
        alpha = config.get("alpha", 0.05)

        non_trait_cols = {
            "study_id",
            "experiment_id",
            "inc_dec",
            "condition",
            "body_site",
            "signature_size",
        }
        trait_columns = [col for col in profiles_df.columns if col not in non_trait_cols]

        m_universe = len(taxa_matrix_df)
        k_totals = (taxa_matrix_df[trait_columns] > 0).sum(axis=0).to_dict()

        results = []

        for _, row in profiles_df.iterrows():
            study_id = row.get("study_id", "study_1")
            exp_id = row.get("experiment_id", "exp_1")
            inc_dec = row.get("inc_dec", "inc")
            profile_id = f"{study_id}:{exp_id}:{inc_dec}"

            n_draws = int(row.get("signature_size", 10))

            for trait in trait_columns:
                k_observed = int(row[trait])
                k_pop = k_totals.get(trait, 0)

                if m_universe <= 0 or n_draws <= 0 or k_pop <= 0:
                    continue

                if test_type == "enrichment":
                    p_val = float(hypergeom.sf(k_observed - 1, m_universe, k_pop, n_draws))
                    expected_k = n_draws * (k_pop / m_universe)
                    effect_size = (k_observed / expected_k) if expected_k > 0 else 0.0
                else:
                    p_val = float(hypergeom.cdf(k_observed, m_universe, k_pop, n_draws))
                    expected_k = n_draws * (k_pop / m_universe)
                    effect_size = (k_observed / expected_k) if expected_k > 0 else 0.0

                results.append(
                    {
                        "profile_id": profile_id,
                        "trait": trait,
                        "k_observed": k_observed,
                        "k_pop": k_pop,
                        "m_universe": m_universe,
                        "n_draws": n_draws,
                        "effect_size": effect_size,
                        "p_value": p_val,
                    }
                )

        res_df = pd.DataFrame(results)
        if not res_df.empty:
            p_vals = res_df["p_value"].to_numpy()
            n_tests = len(p_vals)
            sorted_indices = np.argsort(p_vals)
            sorted_p = p_vals[sorted_indices]

            q_values = np.zeros(n_tests)
            cumulative_min = 1.0
            for i in range(n_tests - 1, -1, -1):
                rank = i + 1
                q = (sorted_p[i] * n_tests) / rank
                cumulative_min = min(cumulative_min, q)
                q_values[sorted_indices[i]] = cumulative_min

            res_df["q_value"] = np.clip(q_values, 0.0, 1.0)
            res_df["is_significant"] = res_df["q_value"] <= alpha
        else:
            res_df = pd.DataFrame(
                columns=[
                    "profile_id",
                    "trait",
                    "effect_size",
                    "p_value",
                    "q_value",
                    "is_significant",
                ]
            )

        return MethodInferenceResult(
            method_id=self.method_id,
            owner=self.owner,
            target_null=self.target_null,
            results_df=res_df,
            metadata={"test_type": test_type, "alpha": alpha},
        )
