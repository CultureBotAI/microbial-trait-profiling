"""Numerical reference-value and consolidation tests for enrichment statistics (#366, #355)."""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest
from scipy.stats import hypergeom as scipy_hypergeom

from bsdbng.enrichment_stats import (
    benjamini_hochberg,
    compute_fold_enrichment,
    hypergeom_cdf,
    hypergeom_cdf_many,
    hypergeom_cdf_pure_python,
    hypergeom_sf,
    hypergeom_sf_many,
    hypergeom_sf_pure_python,
    hypergeom_test,
)

# --- Benjamini-Hochberg (BH) Numerical Reference Tests ------------------------


def test_bh_empty_input() -> None:
    """Empty list or array returns empty container of matching type."""
    assert benjamini_hochberg([]) == []
    res_arr = benjamini_hochberg(np.empty(0, dtype=float))
    assert isinstance(res_arr, np.ndarray)
    assert len(res_arr) == 0


def test_bh_single_element() -> None:
    """Single p-value returns unchanged value in matching container."""
    assert benjamini_hochberg([0.042]) == pytest.approx([0.042])
    res = benjamini_hochberg(np.array([0.042]))
    assert isinstance(res, np.ndarray)
    np.testing.assert_allclose(res, [0.042])


def test_bh_four_element_analytical_reference() -> None:
    """Analytical benchmark: p = [0.01, 0.04, 0.03, 0.20].

    Sorted p: 0.01, 0.03, 0.04, 0.20
    Ranks:    1,    2,    3,    4
    p*n/rank: 0.04, 0.06, 0.16/3 (~0.05333), 0.20
    Step-up cumulative min: [0.04, 0.16/3, 0.16/3, 0.20]
    Mapped to original positions: [0.04, 0.16/3, 0.16/3, 0.20].
    """
    p_vals = [0.01, 0.04, 0.03, 0.20]
    expected = [0.04, 0.16 / 3, 0.16 / 3, 0.20]
    q_list = benjamini_hochberg(p_vals)
    assert isinstance(q_list, list)
    assert q_list == pytest.approx(expected, rel=1e-12)

    q_arr = benjamini_hochberg(np.array(p_vals))
    assert isinstance(q_arr, np.ndarray)
    np.testing.assert_allclose(q_arr, expected, rtol=1e-12)


def test_bh_five_element_analytical_reference() -> None:
    """Analytical benchmark: p = [0.01, 0.02, 0.03, 0.04, 0.05].

    p_i * 5 / i = 0.05 for all i.
    All adjusted values must equal 0.05.
    """
    p_vals = [0.01, 0.02, 0.03, 0.04, 0.05]
    expected = [0.05, 0.05, 0.05, 0.05, 0.05]
    assert benjamini_hochberg(p_vals) == pytest.approx(expected, rel=1e-12)


def test_bh_six_element_analytical_reference() -> None:
    """Analytical benchmark: p = [0.005, 0.01, 0.02, 0.04, 0.05, 0.10].

    p*6/rank:
    i=1: 0.005 * 6 / 1 = 0.03
    i=2: 0.01 * 6 / 2 = 0.03
    i=3: 0.02 * 6 / 3 = 0.04
    i=4: 0.04 * 6 / 4 = 0.06
    i=5: 0.05 * 6 / 5 = 0.06
    i=6: 0.10 * 6 / 6 = 0.10
    Cumulative min: [0.03, 0.03, 0.04, 0.06, 0.06, 0.10].
    """
    p_vals = [0.005, 0.01, 0.02, 0.04, 0.05, 0.10]
    expected = [0.03, 0.03, 0.04, 0.06, 0.06, 0.10]
    assert benjamini_hochberg(p_vals) == pytest.approx(expected, rel=1e-12)


def test_bh_ties_receive_identical_q_values() -> None:
    """Tied p-values must receive identical adjusted q-values regardless of input order."""
    # p = [0.02, 0.02, 0.02, 0.10], n = 4.
    # rank 1: 0.02*4/1 = 0.08
    # rank 2: 0.02*4/2 = 0.04
    # rank 3: 0.02*4/3 = 0.08/3 (~0.02667)
    # rank 4: 0.10*4/4 = 0.10
    # Step-up cummin: [0.08/3, 0.08/3, 0.08/3, 0.10].
    p_vals = [0.02, 0.02, 0.02, 0.10]
    expected = [0.08 / 3, 0.08 / 3, 0.08 / 3, 0.10]
    assert benjamini_hochberg(p_vals) == pytest.approx(expected, rel=1e-12)

    # Scramble order: ties separated
    scrambled = [0.02, 0.10, 0.02, 0.02]
    scrambled_exp = [0.08 / 3, 0.10, 0.08 / 3, 0.08 / 3]
    assert benjamini_hochberg(scrambled) == pytest.approx(scrambled_exp, rel=1e-12)


def test_bh_boundary_conditions() -> None:
    """Check boundary conditions at 0.0, 1.0, and capping at 1.0."""
    # All zeros
    assert benjamini_hochberg([0.0, 0.0, 0.0]) == [0.0, 0.0, 0.0]

    # All ones
    assert benjamini_hochberg([1.0, 1.0, 1.0]) == [1.0, 1.0, 1.0]

    # Zero, mid, one
    # 0.0*3/1 = 0.0; 0.5*3/2 = 0.75; 1.0*3/3 = 1.0
    assert benjamini_hochberg([0.0, 0.5, 1.0]) == pytest.approx([0.0, 0.75, 1.0], rel=1e-12)

    # Values that would exceed 1.0 must be capped at 1.0
    # p = [0.8, 0.9], n = 2:
    # rank 1: min(1.0, 0.8 * 2 / 1) = 1.0
    # rank 2: min(1.0, 0.9 * 2 / 2) = 0.9
    # Step-up cummin enforces monotone order: [0.9, 0.9].
    assert benjamini_hochberg([0.8, 0.9]) == pytest.approx([0.9, 0.9], rel=1e-12)


def test_bh_mathematical_invariants() -> None:
    """Test monotonicity, domination, and bounds across random p-values."""
    rng = np.random.default_rng(42)
    for _ in range(20):
        n = rng.integers(5, 50)
        p = rng.uniform(0.0, 1.0, size=n)
        q = benjamini_hochberg(p)

        # Domination: q_i >= p_i
        assert np.all(q >= p - 1e-12)

        # Bounds: 0 <= q_i <= 1
        assert np.all(q >= 0.0)
        assert np.all(q <= 1.0)

        # Monotonicity: if p_i <= p_j then q_i <= q_j
        sort_idx = np.argsort(p)
        sorted_q = q[sort_idx]
        assert np.all(np.diff(sorted_q) >= -1e-12)


def test_bh_rejects_invalid_p_values() -> None:
    """Invalid p-values (negative, >1, nan, inf) must raise ValueError."""
    with pytest.raises(ValueError, match=r"between 0\.0 and 1\.0"):
        benjamini_hochberg([-0.01, 0.05])
    with pytest.raises(ValueError, match=r"between 0\.0 and 1\.0"):
        benjamini_hochberg([0.05, 1.05])
    with pytest.raises(ValueError, match="finite real numbers"):
        benjamini_hochberg([0.05, float("nan")])
    with pytest.raises(ValueError, match="finite real numbers"):
        benjamini_hochberg([0.05, float("inf")])


def test_bh_family_size() -> None:
    """benjamini_hochberg supports effective family_size larger than input length."""
    # p = [0.01], family_size = 10 -> q = 0.01 * 10 / 1 = 0.10
    assert benjamini_hochberg([0.01], family_size=10) == pytest.approx([0.10])
    # p = [0.01, 0.02], family_size = 100
    assert benjamini_hochberg([0.01, 0.02], family_size=100) == pytest.approx([1.0, 1.0])


def test_bh_rejects_non_numeric() -> None:
    """benjamini_hochberg rejects non-numeric entries (bool, str) with TypeError."""
    bad_bools: list[Any] = [True, False]
    with pytest.raises(TypeError, match="p-values must be numeric"):
        benjamini_hochberg(bad_bools)
    bad_strs: list[Any] = ["0.05", "0.10"]
    with pytest.raises(TypeError, match="p-values must be numeric"):
        benjamini_hochberg(bad_strs)


# --- Hypergeometric Exact Analytical Values -----------------------------------


def test_hypergeom_small_analytical_values() -> None:
    """Exact hand-calculated values for M=4, K=2, N=2.

    Possible values of k in {0, 1, 2}.
    comb(2,0)*comb(2,2)/comb(4,2) = 1*1/6 = 1/6
    comb(2,1)*comb(2,1)/comb(4,2) = 2*2/6 = 4/6
    comb(2,2)*comb(2,0)/comb(4,2) = 1*1/6 = 1/6

    P[X >= 0] = 1.0
    P[X >= 1] = 5/6
    P[X >= 2] = 1/6

    P[X <= 0] = 1/6
    P[X <= 1] = 5/6
    P[X <= 2] = 1.0
    """
    M, K, N = 4, 2, 2
    # Upper tail (survival function / enrichment)
    assert hypergeom_sf(0, M, K, N) == pytest.approx(1.0)
    assert hypergeom_sf(1, M, K, N) == pytest.approx(5 / 6)
    assert hypergeom_sf(2, M, K, N) == pytest.approx(1 / 6)

    # Lower tail (cumulative distribution function / depletion)
    assert hypergeom_cdf(0, M, K, N) == pytest.approx(1 / 6)
    assert hypergeom_cdf(1, M, K, N) == pytest.approx(5 / 6)
    assert hypergeom_cdf(2, M, K, N) == pytest.approx(1.0)

    # Pure Python fallback must yield the exact same values
    assert hypergeom_sf_pure_python(0, M, K, N) == pytest.approx(1.0)
    assert hypergeom_sf_pure_python(1, M, K, N) == pytest.approx(5 / 6)
    assert hypergeom_sf_pure_python(2, M, K, N) == pytest.approx(1 / 6)

    assert hypergeom_cdf_pure_python(0, M, K, N) == pytest.approx(1 / 6)
    assert hypergeom_cdf_pure_python(1, M, K, N) == pytest.approx(5 / 6)
    assert hypergeom_cdf_pure_python(2, M, K, N) == pytest.approx(1.0)


def test_hypergeom_asymmetric_analytical_values() -> None:
    """Exact hand-calculated values for M=4, K=2, N=3.

    Support: k in {1, 2}.
    k=1: comb(2,1)*comb(2,2)/comb(4,3) = 2*1/4 = 2/4 = 0.5
    k=2: comb(2,2)*comb(2,1)/comb(4,3) = 1*2/4 = 2/4 = 0.5

    P[X >= 0] = 1.0
    P[X >= 1] = 1.0
    P[X >= 2] = 0.5
    P[X >= 3] = 0.0

    P[X <= 0] = 0.0
    P[X <= 1] = 0.5
    P[X <= 2] = 1.0
    """
    M, K, N = 4, 2, 3
    assert hypergeom_sf(0, M, K, N) == pytest.approx(1.0)
    assert hypergeom_sf(1, M, K, N) == pytest.approx(1.0)
    assert hypergeom_sf(2, M, K, N) == pytest.approx(0.5)
    assert hypergeom_sf(3, M, K, N) == pytest.approx(0.0)

    assert hypergeom_cdf(0, M, K, N) == pytest.approx(0.0)
    assert hypergeom_cdf(1, M, K, N) == pytest.approx(0.5)
    assert hypergeom_cdf(2, M, K, N) == pytest.approx(1.0)

    # Pure-python fallback
    assert hypergeom_sf_pure_python(0, M, K, N) == pytest.approx(1.0)
    assert hypergeom_sf_pure_python(1, M, K, N) == pytest.approx(1.0)
    assert hypergeom_sf_pure_python(2, M, K, N) == pytest.approx(0.5)
    assert hypergeom_sf_pure_python(3, M, K, N) == pytest.approx(0.0)

    assert hypergeom_cdf_pure_python(0, M, K, N) == pytest.approx(0.0)
    assert hypergeom_cdf_pure_python(1, M, K, N) == pytest.approx(0.5)
    assert hypergeom_cdf_pure_python(2, M, K, N) == pytest.approx(1.0)


# --- Pure-Python Fallback vs SciPy Reference Numerical Comparison -------------


@pytest.mark.parametrize(
    ("M", "K", "N", "x_values"),
    [
        (4, 2, 2, [0, 1, 2, 3]),
        (10, 4, 3, [0, 1, 2, 3, 4]),
        (100, 20, 10, [0, 1, 2, 3, 5, 8, 10, 12]),
        (500, 5, 100, [0, 1, 2, 3, 4, 5]),
        (2654, 300, 50, [0, 2, 5, 10, 15, 20, 30]),
    ],
)
def test_pure_python_fallback_matches_scipy(M: int, K: int, N: int, x_values: list[int]) -> None:
    """Pure-Python fallback must match SciPy to within 1e-10 across all regimes."""
    for x in x_values:
        # Survival function (P[X >= x])
        py_sf = hypergeom_sf_pure_python(x, M, K, N)
        sp_sf = float(scipy_hypergeom.sf(x - 1, M, K, N))
        assert py_sf == pytest.approx(sp_sf, abs=1e-10)

        # Cumulative distribution function (P[X <= x])
        py_cdf = hypergeom_cdf_pure_python(x, M, K, N)
        sp_cdf = float(scipy_hypergeom.cdf(x, M, K, N))
        assert py_cdf == pytest.approx(sp_cdf, abs=1e-10)

        # Wrapper with use_pure_python=True vs default
        assert hypergeom_sf(x, M, K, N, use_pure_python=True) == pytest.approx(
            hypergeom_sf(x, M, K, N, use_pure_python=False), abs=1e-10
        )
        assert hypergeom_cdf(x, M, K, N, use_pure_python=True) == pytest.approx(
            hypergeom_cdf(x, M, K, N, use_pure_python=False), abs=1e-10
        )


def test_hypergeom_vector_pure_python_matches_scipy() -> None:
    """Vectorized functions with use_pure_python=True must match SciPy output."""
    counts = [0, 2, 5, 10, 15]
    successes = [10, 50, 100, 200, 300]
    M = 2654
    N = 50

    sf_py = hypergeom_sf_many(counts, M, successes, N, use_pure_python=True)
    sf_sp = hypergeom_sf_many(counts, M, successes, N, use_pure_python=False)
    assert sf_py == pytest.approx(sf_sp, abs=1e-10)

    cdf_py = hypergeom_cdf_many(counts, M, successes, N, use_pure_python=True)
    cdf_sp = hypergeom_cdf_many(counts, M, successes, N, use_pure_python=False)
    assert cdf_py == pytest.approx(cdf_sp, abs=1e-10)


def test_hypergeom_boundary_extremes() -> None:
    """Test boundary extremes: K=0, N=0, K=M, N=M, x < k_min, x > k_max."""
    # K = 0 (no successes in population)
    assert hypergeom_sf(0, 100, 0, 10) == pytest.approx(1.0)
    assert hypergeom_sf(1, 100, 0, 10) == pytest.approx(0.0)
    assert hypergeom_cdf(0, 100, 0, 10) == pytest.approx(1.0)

    assert hypergeom_sf_pure_python(0, 100, 0, 10) == pytest.approx(1.0)
    assert hypergeom_sf_pure_python(1, 100, 0, 10) == pytest.approx(0.0)
    assert hypergeom_cdf_pure_python(0, 100, 0, 10) == pytest.approx(1.0)

    # N = 0 (zero draws)
    assert hypergeom_sf(0, 100, 20, 0) == pytest.approx(1.0)
    assert hypergeom_sf(1, 100, 20, 0) == pytest.approx(0.0)
    assert hypergeom_cdf(0, 100, 20, 0) == pytest.approx(1.0)

    assert hypergeom_sf_pure_python(0, 100, 20, 0) == pytest.approx(1.0)
    assert hypergeom_sf_pure_python(1, 100, 20, 0) == pytest.approx(0.0)
    assert hypergeom_cdf_pure_python(0, 100, 20, 0) == pytest.approx(1.0)

    # K = M (all items are successes in population)
    assert hypergeom_sf(10, 50, 50, 10) == pytest.approx(1.0)
    assert hypergeom_sf(11, 50, 50, 10) == pytest.approx(0.0)
    assert hypergeom_sf(9, 50, 50, 10) == pytest.approx(1.0)
    assert hypergeom_cdf(10, 50, 50, 10) == pytest.approx(1.0)
    assert hypergeom_cdf(9, 50, 50, 10) == pytest.approx(0.0)

    assert hypergeom_sf_pure_python(10, 50, 50, 10) == pytest.approx(1.0)
    assert hypergeom_sf_pure_python(11, 50, 50, 10) == pytest.approx(0.0)
    assert hypergeom_cdf_pure_python(10, 50, 50, 10) == pytest.approx(1.0)
    assert hypergeom_cdf_pure_python(9, 50, 50, 10) == pytest.approx(0.0)

    # N = M (entire urn drawn)
    assert hypergeom_sf(20, 50, 20, 50) == pytest.approx(1.0)
    assert hypergeom_sf(21, 50, 20, 50) == pytest.approx(0.0)
    assert hypergeom_sf(19, 50, 20, 50) == pytest.approx(1.0)
    assert hypergeom_cdf(20, 50, 20, 50) == pytest.approx(1.0)
    assert hypergeom_cdf(19, 50, 20, 50) == pytest.approx(0.0)

    assert hypergeom_sf_pure_python(20, 50, 20, 50) == pytest.approx(1.0)
    assert hypergeom_sf_pure_python(21, 50, 20, 50) == pytest.approx(0.0)
    assert hypergeom_cdf_pure_python(20, 50, 20, 50) == pytest.approx(1.0)
    assert hypergeom_cdf_pure_python(19, 50, 20, 50) == pytest.approx(0.0)


def test_hypergeom_input_validation() -> None:
    """Invalid parameter checks must raise ValueError or TypeError."""
    # K > M
    with pytest.raises(ValueError, match="Invalid hypergeometric"):
        hypergeom_sf(1, 4, 5, 2)
    with pytest.raises(ValueError, match="Invalid hypergeometric"):
        hypergeom_sf_pure_python(1, 4, 5, 2)
    with pytest.raises(ValueError, match="Invalid hypergeometric"):
        hypergeom_cdf(1, 4, 5, 2)
    with pytest.raises(ValueError, match="Invalid hypergeometric"):
        hypergeom_cdf_pure_python(1, 4, 5, 2)

    # N > M
    with pytest.raises(ValueError, match="Invalid hypergeometric"):
        hypergeom_sf(1, 4, 2, 5)
    with pytest.raises(ValueError, match="Invalid hypergeometric"):
        hypergeom_sf_pure_python(1, 4, 2, 5)

    # Negative population size
    with pytest.raises(ValueError, match="Population size must be strictly positive"):
        hypergeom_sf(1, -4, 2, 2)
    # Negative population successes or draws
    with pytest.raises(ValueError, match="Invalid hypergeometric"):
        hypergeom_sf(1, 4, -2, 2)
    with pytest.raises(ValueError, match="Invalid hypergeometric"):
        hypergeom_sf(1, 4, 2, -2)

    # Population size <= 0
    with pytest.raises(ValueError, match="Population size must be strictly positive"):
        hypergeom_sf(1, 0, 0, 0)

    # Rejection of bool, float, str for integer parameters
    bad_val: Any = True
    with pytest.raises(TypeError, match="must be an integer"):
        hypergeom_sf(bad_val, 100, 20, 10)
    bad_val = 1.5
    with pytest.raises(TypeError, match="must be an integer"):
        hypergeom_sf(1, bad_val, 20, 10)
    bad_val = "10"
    with pytest.raises(TypeError, match="must be an integer"):
        hypergeom_sf(1, 100, bad_val, 10)

    # Vector length mismatch
    with pytest.raises(ValueError, match="same length"):
        hypergeom_sf_many([1, 2], 4, [2], 2)
    with pytest.raises(ValueError, match="same length"):
        hypergeom_cdf_many([1, 2], 4, [2], 2)


def test_compute_fold_enrichment() -> None:
    """Test fold enrichment calculation across standard, edge, and error cases."""
    # Standard enrichment: M=100, K=20, N=10, x=4 -> Expected = 2.0, Fold = 2.0
    assert compute_fold_enrichment(4, 100, 20, 10) == pytest.approx(2.0)

    # Standard depletion: M=100, K=50, N=10, x=1 -> Expected = 5.0, Fold = 0.2
    assert compute_fold_enrichment(1, 100, 50, 10) == pytest.approx(0.2)

    # Boundary: K=M (all successes) -> Expected = N, Fold = 1.0
    assert compute_fold_enrichment(10, 50, 50, 10) == pytest.approx(1.0)

    # Boundary: N=M (entire urn drawn) -> Expected = K, Fold = 1.0
    assert compute_fold_enrichment(20, 50, 20, 50) == pytest.approx(1.0)

    # Zero observed: x=0, expected > 0 -> Fold = 0.0
    assert compute_fold_enrichment(0, 100, 20, 10) == pytest.approx(0.0)

    # Expected zero (K=0) with x=0 -> Fold = 0.0
    assert compute_fold_enrichment(0, 100, 0, 10) == pytest.approx(0.0)

    # Expected zero (K=0) with x>0 -> Fold = inf
    assert compute_fold_enrichment(1, 100, 0, 10) == float("inf")

    # Invalid x (negative or > draws)
    with pytest.raises(ValueError, match="Invalid observed successes"):
        compute_fold_enrichment(-1, 100, 20, 10)
    with pytest.raises(ValueError, match="Invalid observed successes"):
        compute_fold_enrichment(15, 100, 20, 10)

    # Invalid types
    bad_x: Any = True
    with pytest.raises(TypeError, match="must be an integer"):
        compute_fold_enrichment(bad_x, 100, 20, 10)


def test_hypergeom_test() -> None:
    """Test hypergeom_test direction, aliases, and equivalence to sf/cdf."""
    M, K, N, x = 100, 20, 10, 4
    # "greater", "upper", "enrichment" must match hypergeom_sf
    expected_upper = hypergeom_sf(x, M, K, N)
    for alt in ("greater", "upper", "enrichment", "GREATER", "Upper"):
        assert hypergeom_test(x, M, K, N, alternative=alt) == pytest.approx(expected_upper)
        assert hypergeom_test(x, M, K, N, alternative=alt, use_pure_python=True) == pytest.approx(
            expected_upper
        )

    # "less", "lower", "depletion" must match hypergeom_cdf
    expected_lower = hypergeom_cdf(x, M, K, N)
    for alt in ("less", "lower", "depletion", "LESS", "Lower"):
        assert hypergeom_test(x, M, K, N, alternative=alt) == pytest.approx(expected_lower)
        assert hypergeom_test(x, M, K, N, alternative=alt, use_pure_python=True) == pytest.approx(
            expected_lower
        )

    # Invalid alternative raises ValueError
    with pytest.raises(ValueError, match="Unknown alternative"):
        hypergeom_test(x, M, K, N, alternative="two-sided")
