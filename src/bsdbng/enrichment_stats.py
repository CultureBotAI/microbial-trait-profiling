"""Consolidated statistical functions for trait enrichment and depletion analyses (#366, #355).

This module provides standard, validated implementations of:
- Benjamini-Hochberg FDR correction (:func:`benjamini_hochberg`), delegating to
  :func:`bsdbng.trait_screen.step_up`
- Hypergeometric test (:func:`hypergeom_test`) for enrichment and depletion
- Fold enrichment ratio calculation (:func:`compute_fold_enrichment`)
- Hypergeometric survival function (upper tail, enrichment)
  (:func:`hypergeom_sf`, :func:`hypergeom_sf_many`)
- Hypergeometric cumulative distribution function (lower tail, depletion)
  (:func:`hypergeom_cdf`, :func:`hypergeom_cdf_many`)
- Pure-Python log-gamma reference oracles
  (:func:`hypergeom_sf_pure_python`, :func:`hypergeom_cdf_pure_python`)
  for verification and numerical validation against SciPy.
"""

from __future__ import annotations

import math
import operator
from collections.abc import Sequence
from typing import Any, overload

import numpy as np
from scipy.stats import hypergeom as _scipy_hypergeom

from bsdbng.trait_screen import step_up


def _validate_int(val: Any, name: str) -> int:
    """Validate that val is an integer index (rejecting booleans, floats, and strings)."""
    if isinstance(val, bool):
        raise TypeError(f"{name} must be an integer, got bool")
    try:
        return operator.index(val)
    except TypeError:
        raise TypeError(f"{name} must be an integer, got {type(val).__name__}") from None


def _validate_p_values(p_values: Sequence[float] | np.ndarray) -> np.ndarray:
    """Validate that p-values are numeric and in [0, 1]."""
    for v in p_values:
        if isinstance(v, (str, bool)):
            raise TypeError("p-values must be numeric")
    p_arr = np.asarray(p_values, dtype=float)
    if not np.all(np.isfinite(p_arr)):
        raise ValueError("p-values must be finite real numbers")
    if np.any(p_arr < 0.0) or np.any(p_arr > 1.0):
        raise ValueError("p-values must be between 0.0 and 1.0 inclusive")
    return p_arr


@overload
def benjamini_hochberg(
    p_values: list[float],
    *,
    family_size: int | None = None,
    method: str = "bh",
) -> list[float]: ...


@overload
def benjamini_hochberg(
    p_values: np.ndarray,
    *,
    family_size: int | None = None,
    method: str = "bh",
) -> np.ndarray: ...


@overload
def benjamini_hochberg(
    p_values: Sequence[float],
    *,
    family_size: int | None = None,
    method: str = "bh",
) -> list[float]: ...


def benjamini_hochberg(
    p_values: Sequence[float] | np.ndarray,
    *,
    family_size: int | None = None,
    method: str = "bh",
) -> list[float] | np.ndarray:
    """Apply Benjamini-Hochberg (1995) FDR control to adjust p-values into q-values.

    Delegates to :func:`bsdbng.trait_screen.step_up` for monotone rank step-up adjustment.

    Guarantees:
    - Returns the same container type as input (list -> list, ndarray -> ndarray).
    - Monotonicity: if p_i <= p_j then q_i <= q_j.
    - Ties: identical p-values receive identical adjusted q-values.
    - Domination and bounding: p_i <= q_i <= 1.0 for all valid p_i in [0, 1].
    - Edge cases: empty sequence returns empty container; single value returns itself.
    """
    is_numpy = isinstance(p_values, np.ndarray)
    p_arr = _validate_p_values(p_values)
    n = len(p_arr)
    if n == 0:
        return np.empty(0, dtype=float) if is_numpy else []

    eff_family = n if family_size is None else _validate_int(family_size, "family_size")
    adjusted = step_up(list(p_arr), method=method, family_size=eff_family)

    if is_numpy:
        return np.asarray(adjusted, dtype=float)
    return adjusted


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
    if population_size <= 0:
        raise ValueError("Population size must be strictly positive (M >= 1)")
    if population_successes < 0 or draws < 0:
        raise ValueError("Invalid hypergeometric population or draw parameters")
    if population_successes > population_size or draws > population_size:
        raise ValueError("Invalid hypergeometric population or draw parameters")


def compute_fold_enrichment(
    x: int,
    population_size: int,
    population_successes: int,
    draws: int,
) -> float:
    """Compute the fold enrichment ratio: observed / expected.

    The expected number of successes under the hypergeometric null is:
        E[X] = draws * (population_successes / population_size)

    Returns
    -------
    float
        Fold enrichment ratio. If expected == 0 and x > 0, returns inf.
        If expected == 0 and x == 0, returns 0.0.
    """
    x_int = _validate_int(x, "x")
    m_int = _validate_int(population_size, "population_size")
    k_int = _validate_int(population_successes, "population_successes")
    n_int = _validate_int(draws, "draws")

    _check_hypergeom_params(m_int, k_int, n_int)
    if x_int < 0 or x_int > n_int:
        raise ValueError("Invalid observed successes x")

    expected = (k_int / m_int) * n_int
    if expected <= 0.0:
        return float("inf") if x_int > 0 else 0.0
    return float(x_int / expected)


def hypergeom_sf_pure_python(
    x: int,
    population_size: int,
    population_successes: int,
    draws: int,
) -> float:
    """Pure-Python log-gamma reference oracle of P[X >= x] for Hypergeometric(M, K, N)."""
    x_int = _validate_int(x, "x")
    m_int = _validate_int(population_size, "population_size")
    k_int = _validate_int(population_successes, "population_successes")
    n_int = _validate_int(draws, "draws")

    _check_hypergeom_params(m_int, k_int, n_int)
    k_min = max(0, n_int - (m_int - k_int))
    k_max = min(k_int, n_int)

    if x_int <= k_min:
        return 1.0
    if x_int > k_max:
        return 0.0

    log_denom = _log_comb(m_int, n_int)
    log_probs = [
        _log_comb(k_int, k) + _log_comb(m_int - k_int, n_int - k) - log_denom
        for k in range(x_int, k_max + 1)
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
    """Pure-Python log-gamma reference oracle of P[X <= x] for Hypergeometric(M, K, N)."""
    x_int = _validate_int(x, "x")
    m_int = _validate_int(population_size, "population_size")
    k_int = _validate_int(population_successes, "population_successes")
    n_int = _validate_int(draws, "draws")

    _check_hypergeom_params(m_int, k_int, n_int)
    k_min = max(0, n_int - (m_int - k_int))
    k_max = min(k_int, n_int)

    if x_int < k_min:
        return 0.0
    if x_int >= k_max:
        return 1.0

    log_denom = _log_comb(m_int, n_int)
    log_probs = [
        _log_comb(k_int, k) + _log_comb(m_int - k_int, n_int - k) - log_denom
        for k in range(k_min, min(x_int, k_max) + 1)
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
    """Return P[X >= x] using SciPy or pure-Python reference fallback."""
    x_int = _validate_int(x, "x")
    m_int = _validate_int(population_size, "population_size")
    k_int = _validate_int(population_successes, "population_successes")
    n_int = _validate_int(draws, "draws")

    _check_hypergeom_params(m_int, k_int, n_int)
    if use_pure_python:
        return hypergeom_sf_pure_python(x_int, m_int, k_int, n_int)

    hg = scipy_hypergeom if scipy_hypergeom is not None else _scipy_hypergeom
    value = float(hg.sf(x_int - 1, m_int, k_int, n_int))
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
    """Return P[X <= x] using SciPy or pure-Python reference fallback."""
    x_int = _validate_int(x, "x")
    m_int = _validate_int(population_size, "population_size")
    k_int = _validate_int(population_successes, "population_successes")
    n_int = _validate_int(draws, "draws")

    _check_hypergeom_params(m_int, k_int, n_int)
    if use_pure_python:
        return hypergeom_cdf_pure_python(x_int, m_int, k_int, n_int)

    hg = scipy_hypergeom if scipy_hypergeom is not None else _scipy_hypergeom
    value = float(hg.cdf(x_int, m_int, k_int, n_int))
    if not math.isfinite(value):
        raise ValueError("Invalid hypergeometric population or draw parameters")
    return value


def hypergeom_test(
    x: int,
    population_size: int,
    population_successes: int,
    draws: int,
    *,
    alternative: str = "greater",
    use_pure_python: bool = False,
    scipy_hypergeom: Any = None,
) -> float:
    """Exact hypergeometric test under the 2x2 urn model.

    Parameters
    ----------
    x : int
        Number of successes observed in draw (k).
    population_size : int
        Total population / universe size (M).
    population_successes : int
        Number of successes in population (K).
    draws : int
        Number of draws / sample size (N).
    alternative : {"greater", "less", "upper", "lower", "enrichment", "depletion"}
        Test direction. "greater", "upper", or "enrichment" tests for
        over-representation (P[X >= x]). "less", "lower", or "depletion" tests
        for under-representation (P[X <= x]).
    use_pure_python : bool
        If True, use the pure-Python log-gamma implementation.
    scipy_hypergeom : Any
        Optional injectable scipy.stats.hypergeom object.

    Returns
    -------
    float
        Exact tail probability p-value in [0, 1].
    """
    alt = alternative.lower()
    if alt in ("greater", "upper", "enrichment"):
        return hypergeom_sf(
            x,
            population_size,
            population_successes,
            draws,
            use_pure_python=use_pure_python,
            scipy_hypergeom=scipy_hypergeom,
        )
    elif alt in ("less", "lower", "depletion"):
        return hypergeom_cdf(
            x,
            population_size,
            population_successes,
            draws,
            use_pure_python=use_pure_python,
            scipy_hypergeom=scipy_hypergeom,
        )
    else:
        raise ValueError(
            f"Unknown alternative {alternative!r}. Expected 'greater'/'upper' or 'less'/'lower'."
        )


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

    m_int = _validate_int(population_size, "population_size")
    n_int = _validate_int(draws, "draws")
    parsed_counts = [_validate_int(x, "count") for x in counts]
    parsed_successes = [_validate_int(k, "population_success") for k in population_successes]

    for k in parsed_successes:
        _check_hypergeom_params(m_int, k, n_int)

    if use_pure_python:
        return [
            hypergeom_sf_pure_python(x, m_int, k, n_int)
            for x, k in zip(parsed_counts, parsed_successes, strict=True)
        ]

    hg = scipy_hypergeom if scipy_hypergeom is not None else _scipy_hypergeom
    values = hg.sf([x - 1 for x in parsed_counts], m_int, parsed_successes, n_int)
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

    m_int = _validate_int(population_size, "population_size")
    n_int = _validate_int(draws, "draws")
    parsed_counts = [_validate_int(x, "count") for x in counts]
    parsed_successes = [_validate_int(k, "population_success") for k in population_successes]

    for k in parsed_successes:
        _check_hypergeom_params(m_int, k, n_int)

    if use_pure_python:
        return [
            hypergeom_cdf_pure_python(x, m_int, k, n_int)
            for x, k in zip(parsed_counts, parsed_successes, strict=True)
        ]

    hg = scipy_hypergeom if scipy_hypergeom is not None else _scipy_hypergeom
    values = hg.cdf(parsed_counts, m_int, parsed_successes, n_int)
    results = [float(v) for v in values]
    if any(not math.isfinite(v) for v in results):
        raise ValueError("Invalid hypergeometric population or draw parameters")
    return results
