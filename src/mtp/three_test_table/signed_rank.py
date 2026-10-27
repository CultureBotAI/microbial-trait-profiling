"""Versioned secondary study signed-rank inference (#1073; 2026-10-26).

Only the private test vector is zero-filtered. Small untied samples use the exact
signed-rank distribution; small tied samples enumerate every sign assignment.
Neither operation changes the source contrasts or their descriptive summaries.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.stats import PermutationMethod, wilcoxon

from mtp.guild_maps import SIGN_TOLERANCE

LEGACY_WILCOXON_CONTRACT = "study-mean-delta-z-wilcoxon-v1"
WILCOXON_CONTRACT = "study-mean-delta-z-wilcoxon-v2"
ZERO_TOLERANCE = SIGN_TOLERANCE
EXACT_MAX_N = 50
PERMUTATION_MAX_N = 13


@dataclass(frozen=True)
class SignedRankResult:
    """Secondary inference and effective support; untestable internal p remains one."""

    statistic: float
    pvalue: float
    nonzero_studies: int
    zero_studies: int
    method: str
    reason: str


def study_wilcoxon(values: NDArray[np.float64]) -> SignedRankResult:
    """Test finite study means with explicit zero, tie and small-sample policies."""
    vector = np.asarray(values, dtype=np.float64)
    if vector.ndim != 1:
        raise ValueError("Study Wilcoxon observations must be one-dimensional")
    if not np.isfinite(vector).all():
        raise ValueError("Study Wilcoxon observations must be finite")
    active = vector[np.abs(vector) > ZERO_TOLERANCE]
    count = len(active)
    zero_count = len(vector) - count
    if not count:
        return SignedRankResult(
            statistic=float("nan"),
            pvalue=1.0,
            nonzero_studies=0,
            zero_studies=zero_count,
            method="not_tested",
            reason="insufficient_support" if not len(vector) else "all_tied",
        )
    tied = len(np.unique(np.abs(active))) != count
    scipy_method: str | PermutationMethod
    if not tied and count <= EXACT_MAX_N:
        method = "exact"
        scipy_method = "exact"
    elif tied and count <= PERMUTATION_MAX_N:
        method = "exhaustive_sign_permutation"
        scipy_method = PermutationMethod(n_resamples=np.inf)
    else:
        method = "asymptotic"
        scipy_method = "asymptotic"
    result = wilcoxon(
        active,
        alternative="two-sided",
        zero_method="wilcox",
        correction=False,
        method=scipy_method,
    )
    statistic, pvalue = float(result.statistic), float(result.pvalue)
    if not np.isfinite(statistic) or not np.isfinite(pvalue) or not 0 <= pvalue <= 1:
        raise ValueError("Study Wilcoxon returned invalid finite inference")
    return SignedRankResult(statistic, pvalue, count, zero_count, method, "")


__all__ = [
    "EXACT_MAX_N",
    "LEGACY_WILCOXON_CONTRACT",
    "PERMUTATION_MAX_N",
    "WILCOXON_CONTRACT",
    "ZERO_TOLERANCE",
    "SignedRankResult",
    "study_wilcoxon",
]
