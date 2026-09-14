"""Evaluation metrics for statistical calibration, Type I error, and FDR control."""

import math
from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass
class TypeIErrorMetric:
    """Type I error rate evaluation at given significance level."""

    alpha: float
    empirical_rate: float
    monte_carlo_se: float
    replicates: int
    is_conservative: bool  # empirical_rate <= alpha + 2.58 * se (approx 99% upper bound)


@dataclass
class FDRMetric:
    """FDR and False Discovery Proportion (FDP) metric."""

    target_q: float
    empirical_fdr: float  # Mean FDP across replicates
    fdp_std: float
    replicates: int
    mean_discoveries: float
    mean_false_discoveries: float


def compute_type_i_error(
    p_values: np.ndarray | list[float], alpha: float = 0.05
) -> TypeIErrorMetric:
    """Compute empirical Type I error and Monte Carlo standard error under the null."""
    p_arr = np.asarray(p_values, dtype=float)
    p_arr = p_arr[~np.isnan(p_arr)]
    n = len(p_arr)
    if n == 0:
        return TypeIErrorMetric(
            alpha=alpha,
            empirical_rate=0.0,
            monte_carlo_se=0.0,
            replicates=0,
            is_conservative=True,
        )

    rate = float(np.mean(p_arr <= alpha))
    se = math.sqrt((rate * (1.0 - rate)) / n) if n > 1 else 0.0
    # Conservative or calibrated if within normal sampling tolerance
    is_conservative = rate <= (alpha + 2.58 * se)

    return TypeIErrorMetric(
        alpha=alpha,
        empirical_rate=rate,
        monte_carlo_se=se,
        replicates=n,
        is_conservative=is_conservative,
    )


def compute_empirical_fdr(
    replicate_discoveries: list[set[tuple[str, str]]],
    true_signal_pairs: set[tuple[str, str]],
    target_q: float = 0.05,
) -> FDRMetric:
    """Compute empirical False Discovery Proportion (FDP) and FDR across replicates.

    Each discovery is represented as a tuple (profile_id/study_id, trait_id).
    """
    if not replicate_discoveries:
        return FDRMetric(
            target_q=target_q,
            empirical_fdr=0.0,
            fdp_std=0.0,
            replicates=0,
            mean_discoveries=0.0,
            mean_false_discoveries=0.0,
        )

    fdp_list = []
    disc_counts = []
    false_counts = []

    for discoveries in replicate_discoveries:
        r = len(discoveries)
        disc_counts.append(r)
        if r == 0:
            fdp_list.append(0.0)
            false_counts.append(0)
        else:
            false_positives = len(discoveries - true_signal_pairs)
            false_counts.append(false_positives)
            fdp_list.append(false_positives / r)

    return FDRMetric(
        target_q=target_q,
        empirical_fdr=float(np.mean(fdp_list)),
        fdp_std=float(np.std(fdp_list, ddof=1)) if len(fdp_list) > 1 else 0.0,
        replicates=len(replicate_discoveries),
        mean_discoveries=float(np.mean(disc_counts)),
        mean_false_discoveries=float(np.mean(false_counts)),
    )


def compute_discrete_ecdf(p_values: np.ndarray | list[float], bins: int = 20) -> pd.DataFrame:
    """Compute empirical cumulative distribution function compared to theoretical uniform."""
    p_arr = np.asarray(p_values, dtype=float)
    p_arr = p_arr[~np.isnan(p_arr)]
    if len(p_arr) == 0:
        return pd.DataFrame(columns=["threshold", "empirical_cdf", "uniform_cdf", "difference"])

    thresholds = np.linspace(0.0, 1.0, bins + 1)
    empirical = [float(np.mean(p_arr <= t)) for t in thresholds]

    return pd.DataFrame(
        {
            "threshold": thresholds,
            "empirical_cdf": empirical,
            "uniform_cdf": thresholds,
            "difference": [emp - unif for emp, unif in zip(empirical, thresholds, strict=True)],
        }
    )
