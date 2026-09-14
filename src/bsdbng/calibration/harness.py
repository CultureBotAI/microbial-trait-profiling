"""Monte Carlo calibration harness for evaluating statistical pipelines."""

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

from bsdbng.calibration.metrics import (
    FDRMetric,
    TypeIErrorMetric,
    compute_discrete_ecdf,
    compute_empirical_fdr,
    compute_type_i_error,
)
from bsdbng.calibration.protocol import MethodAdapter
from bsdbng.calibration.scenarios import (
    generate_complete_null_ora_profiles,
    generate_mixed_signal_ora_profiles,
)


@dataclass
class CalibrationBenchmarkReport:
    """Consolidated calibration evaluation results."""

    method_id: str
    replicates: int
    complete_null_type_i: TypeIErrorMetric
    complete_null_ecdf: pd.DataFrame
    mixed_signal_fdr: FDRMetric
    raw_p_values_null: list[float]


def run_ora_calibration_benchmark(
    adapter: MethodAdapter,
    empirical_profiles_df: pd.DataFrame,
    taxa_matrix_df: pd.DataFrame,
    spike_traits: list[str],
    replicates: int = 50,
    seed: int = 42,
    alpha: float = 0.05,
    config: dict[str, Any] | None = None,
) -> CalibrationBenchmarkReport:
    """Run Monte Carlo calibration benchmark on an ORA adapter."""
    rng = np.random.default_rng(seed)
    config = config or {}

    null_p_values: list[float] = []
    mixed_discoveries: list[set[tuple[str, str]]] = []
    all_true_signal_pairs: set[tuple[str, str]] = set()

    for _ in range(replicates):
        # 1. Complete Null Replicate
        null_sim_df = generate_complete_null_ora_profiles(
            empirical_profiles_df, taxa_matrix_df, rng
        )
        null_res = adapter.run_inference(null_sim_df, taxa_matrix_df, config)
        if not null_res.results_df.empty:
            null_p_values.extend(null_res.results_df["p_value"].tolist())

        # 2. Mixed Signal Replicate
        mixed_sim_df, true_pairs = generate_mixed_signal_ora_profiles(
            empirical_profiles_df,
            taxa_matrix_df,
            spike_traits=spike_traits,
            signal_strength=0.6,
            rng=rng,
        )
        all_true_signal_pairs.update(true_pairs)
        mixed_res = adapter.run_inference(mixed_sim_df, taxa_matrix_df, config)

        replicate_disc = set()
        if not mixed_res.results_df.empty:
            sig_rows = mixed_res.results_df[mixed_res.results_df["is_significant"]]
            for _, sig_row in sig_rows.iterrows():
                replicate_disc.add((sig_row["profile_id"], sig_row["trait"]))
        mixed_discoveries.append(replicate_disc)

    type_i = compute_type_i_error(null_p_values, alpha=alpha)
    ecdf_df = compute_discrete_ecdf(null_p_values)
    fdr_metric = compute_empirical_fdr(mixed_discoveries, all_true_signal_pairs, target_q=alpha)

    return CalibrationBenchmarkReport(
        method_id=adapter.method_id,
        replicates=replicates,
        complete_null_type_i=type_i,
        complete_null_ecdf=ecdf_df,
        mixed_signal_fdr=fdr_metric,
        raw_p_values_null=null_p_values,
    )
