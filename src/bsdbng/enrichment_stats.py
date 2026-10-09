"""Consolidated statistical functions for trait enrichment and depletion analyses.

This module provides standard, validated implementations of:
- Benjamini-Hochberg FDR correction (:func:`benjamini_hochberg`)
- Hypergeometric survival function (upper tail, enrichment)
  (:func:`hypergeom_sf`, :func:`hypergeom_sf_many`)
- Hypergeometric cumulative distribution function (lower tail, depletion)
  (:func:`hypergeom_cdf`, :func:`hypergeom_cdf_many`)
- Pure-Python log-gamma fallbacks for verification and environments without SciPy
  (:func:`hypergeom_sf_pure_python`, :func:`hypergeom_cdf_pure_python`)
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any, overload

import numpy as np
from scipy.stats import hypergeom as _scipy_hypergeom


def _validate_p_values(p_arr: np.ndarray) -> None:
    """Validate that p-values are finite numbers in [0, 1]."""
    if not np.all(np.isfinite(p_arr)):
        raise ValueError("p-values must be finite real numbers")
    if np.any(p_arr < 0.0) or np.any(p_arr > 1.0):
        raise ValueError("p-values must be between 0.0 and 1.0 inclusive")


@overload
def benjamini_hochberg(p_values: list[float]) -> list[float]: ...


@overload
def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray: ...


@overload
def benjamini_hochberg(p_values: Sequence[float]) -> list[float]: ...


def benjamini_hochberg(
    p_values: Sequence[float] | np.ndarray,
) -> list[float] | np.ndarray:
    """Apply Benjamini-Hochberg (1995) FDR control to adjust p-values into q-values.

    Given a list or 1D array of p-values in [0, 1], computes adjusted q-values:
        q_{(i)} = min(1, min_{j >= i} (p_{(j)} * n / j))

    Guarantees:
    - Returns the same container type as input (list -> list, ndarray -> ndarray).
    - Monotonicity: if p_i <= p_j then q_i <= q_j.
    - Ties: identical p-values receive identical adjusted q-values.
    - Domination and bounding: p_i <= q_i <= 1.0 for all valid p_i in [0, 1].
    - Edge cases: empty sequence returns empty container; single value returns itself.
    """
    is_numpy = isinstance(p_values, np.ndarray)
    p_arr = np.asarray(p_values, dtype=float)
    n = len(p_arr)
    if n == 0:
        return np.empty(0, dtype=float) if is_numpy else []

    _validate_p_values(p_arr)

    if n == 1:
        val = float(np.clip(p_arr[0], 0.0, 1.0))
        if is_numpy:
            q_single = np.empty(1, dtype=float)
            q_single[0] = val
            return q_single
        return [val]

    order = np.argsort(p_arr)
    ranked = p_arr[order] * n / np.arange(1, n + 1)
    cummin = np.minimum.accumulate(ranked[::-1])[::-1]
    q_values = np.empty(n, dtype=float)
    q_values[order] = np.clip(cummin, 0.0, 1.0)

    if is_numpy:
        return q_values
    return [float(v) for v in q_values]


def _log_comb(n: int, k: int) -> float:
    """Compute ln(n choose k) using math.lgamma for numerical stability."""
    if k < 0 or k > n:
        return -float("inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def _check_hypergeom_params(
    population_size: int,
    population_successes: int,
    draws: int,
) -> None:
    """Validate standard urn hypergeometric parameters (M, K, N)."""
    if population_size < 0 or population_successes < 0 or draws < 0:
        raise ValueError("Invalid hypergeometric population or draw parameters")
    if population_successes > population_size or draws > population_size:
        raise ValueError("Invalid hypergeometric population or draw parameters")


def hypergeom_sf_pure_python(
    x: int,
    population_size: int,
    population_successes: int,
    draws: int,
) -> float:
    """Pure-Python log-gamma implementation of P[X >= x] for Hypergeometric(M, K, N)."""
    _check_hypergeom_params(population_size, population_successes, draws)
    k_min = max(0, draws - (population_size - population_successes))
    k_max = min(population_successes, draws)

    if x <= k_min:
        return 1.0
    if x > k_max:
        return 0.0

    log_denom = _log_comb(population_size, draws)
    log_probs = [
        _log_comb(population_successes, k)
        + _log_comb(population_size - population_successes, draws - k)
        - log_denom
        for k in range(x, k_max + 1)
    ]
    max_lp = max(log_probs)
    if max_lp == -float("inf"):
        return 0.0

    sf = math.exp(max_lp) * sum(math.exp(lp - max_lp) for lp in log_probs)
    return min(1.0, max(0.0, sf))


def hypergeom_cdf_pure_python(
    x: int,
    population_size: int,
    population_successes: int,
    draws: int,
) -> float:
    """Pure-Python log-gamma implementation of P[X <= x] for Hypergeometric(M, K, N)."""
    _check_hypergeom_params(population_size, population_successes, draws)
    k_min = max(0, draws - (population_size - population_successes))
    k_max = min(population_successes, draws)

    if x < k_min:
        return 0.0
    if x >= k_max:
        return 1.0

    log_denom = _log_comb(population_size, draws)
    log_probs = [
        _log_comb(population_successes, k)
        + _log_comb(population_size - population_successes, draws - k)
        - log_denom
        for k in range(k_min, min(x, k_max) + 1)
    ]
    max_lp = max(log_probs)
    if max_lp == -float("inf"):
        return 0.0

    cdf = math.exp(max_lp) * sum(math.exp(lp - max_lp) for lp in log_probs)
    return min(1.0, max(0.0, cdf))


def hypergeom_sf(
    x: int,
    population_size: int,
    population_successes: int,
    draws: int,
    *,
    use_pure_python: bool = False,
    scipy_hypergeom: Any = None,
) -> float:
    """Return P[X >= x] using SciPy or pure-Python fallback."""
    _check_hypergeom_params(population_size, population_successes, draws)
    if use_pure_python:
        return hypergeom_sf_pure_python(x, population_size, population_successes, draws)

    hg = scipy_hypergeom if scipy_hypergeom is not None else _scipy_hypergeom
    value = float(hg.sf(x - 1, population_size, population_successes, draws))
    if not math.isfinite(value):
        raise ValueError("Invalid hypergeometric population or draw parameters")
    return value


def hypergeom_cdf(
    x: int,
    population_size: int,
    population_successes: int,
    draws: int,
    *,
    use_pure_python: bool = False,
    scipy_hypergeom: Any = None,
) -> float:
    """Return P[X <= x] using SciPy or pure-Python fallback."""
    _check_hypergeom_params(population_size, population_successes, draws)
    if use_pure_python:
        return hypergeom_cdf_pure_python(x, population_size, population_successes, draws)

    hg = scipy_hypergeom if scipy_hypergeom is not None else _scipy_hypergeom
    value = float(hg.cdf(x, population_size, population_successes, draws))
    if not math.isfinite(value):
        raise ValueError("Invalid hypergeometric population or draw parameters")
    return value


def hypergeom_sf_many(
    counts: list[int],
    population_size: int,
    population_successes: list[int],
    draws: int,
    *,
    use_pure_python: bool = False,
    scipy_hypergeom: Any = None,
) -> list[float]:
    """Compute tail probabilities P[X >= x] for a vector of trait counts."""
    if len(counts) != len(population_successes):
        raise ValueError("Counts and population successes must have the same length")
    if not counts:
        return []

    for k in population_successes:
        _check_hypergeom_params(population_size, k, draws)

    if use_pure_python:
        return [
            hypergeom_sf_pure_python(x, population_size, k, draws)
            for x, k in zip(counts, population_successes, strict=True)
        ]

    hg = scipy_hypergeom if scipy_hypergeom is not None else _scipy_hypergeom
    values = hg.sf([x - 1 for x in counts], population_size, population_successes, draws)
    results = [float(v) for v in values]
    if any(not math.isfinite(v) for v in results):
        raise ValueError("Invalid hypergeometric population or draw parameters")
    return results


def hypergeom_cdf_many(
    counts: list[int],
    population_size: int,
    population_successes: list[int],
    draws: int,
    *,
    use_pure_python: bool = False,
    scipy_hypergeom: Any = None,
) -> list[float]:
    """Compute tail probabilities P[X <= x] for a vector of trait counts."""
    if len(counts) != len(population_successes):
        raise ValueError("Counts and population successes must have the same length")
    if not counts:
        return []

    for k in population_successes:
        _check_hypergeom_params(population_size, k, draws)

    if use_pure_python:
        return [
            hypergeom_cdf_pure_python(x, population_size, k, draws)
            for x, k in zip(counts, population_successes, strict=True)
        ]

    hg = scipy_hypergeom if scipy_hypergeom is not None else _scipy_hypergeom
    values = hg.cdf(counts, population_size, population_successes, draws)
    results = [float(v) for v in values]
    if any(not math.isfinite(v) for v in results):
        raise ValueError("Invalid hypergeometric population or draw parameters")
    return results
