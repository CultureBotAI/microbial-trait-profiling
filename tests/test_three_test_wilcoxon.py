"""Exact small-n study signed-rank contracts (#1073; 2026-10-26)."""

from __future__ import annotations

from itertools import product

import numpy as np
import pytest
from scipy import stats

from mtp.guild_maps import SIGN_TOLERANCE
from mtp.three_test_table.composition import WILCOXON_CONTRACT as COMPOSITION_WILCOXON_CONTRACT
from mtp.three_test_table.signed_rank import (
    EXACT_MAX_N,
    LEGACY_WILCOXON_CONTRACT,
    PERMUTATION_MAX_N,
    WILCOXON_CONTRACT,
    ZERO_TOLERANCE,
    study_wilcoxon,
)


def enumerated_signed_rank(values: list[float]) -> tuple[float, float]:
    """Enumerate equiprobable signs using hand-assigned average ranks, without SciPy."""
    magnitudes = [abs(value) for value in values]
    ranks = [
        sum(other < value for other in magnitudes)
        + (sum(other == value for other in magnitudes) + 1) / 2
        for value in magnitudes
    ]
    total = sum(ranks)
    positive = sum(rank for rank, value in zip(ranks, values, strict=True) if value > 0)
    observed = min(positive, total - positive)
    at_least_as_extreme = 0
    for assignment in product((False, True), repeat=len(values)):
        positive = sum(rank for rank, sign in zip(ranks, assignment, strict=True) if sign)
        at_least_as_extreme += min(positive, total - positive) <= observed
    return observed, at_least_as_extreme / 2 ** len(values)


def test_version_and_policy_constants_are_explicit() -> None:
    assert LEGACY_WILCOXON_CONTRACT == "study-mean-delta-z-wilcoxon-v1"
    assert WILCOXON_CONTRACT == "study-mean-delta-z-wilcoxon-v2"
    assert COMPOSITION_WILCOXON_CONTRACT == WILCOXON_CONTRACT
    assert ZERO_TOLERANCE == SIGN_TOLERANCE == 1e-12
    assert EXACT_MAX_N == 50
    assert PERMUTATION_MAX_N == 13


@pytest.mark.parametrize(
    ("values", "expected_p"), [([1.0], 1.0), ([1.0, 2.0, 3.0, 4.0, 5.0], 0.0625)]
)
def test_small_untied_one_sign_vectors_have_exact_probabilities(
    values: list[float], expected_p: float
) -> None:
    result = study_wilcoxon(np.array(values, dtype=np.float64))
    assert result.statistic == 0
    assert result.pvalue == expected_p
    assert result.nonzero_studies == len(values)
    assert result.zero_studies == 0
    assert result.method == "exact"
    assert result.reason == ""


@pytest.mark.parametrize(
    ("values", "expected_statistic", "expected_p"),
    [
        ([1.0, 1.0, -1.0], 2.0, 1.0),
        ([1.0, 1.0, 1.0, 1.0, -1.0], 3.0, 0.375),
        ([1.0, -2.0, 2.0, 4.0], 2.5, 0.5),
        ([-1.0, -1.0, -1.0], 0.0, 0.25),
    ],
)
def test_tied_small_vectors_match_independent_sign_enumeration(
    values: list[float], expected_statistic: float, expected_p: float
) -> None:
    expected = enumerated_signed_rank(values)
    assert expected == (expected_statistic, expected_p)
    result = study_wilcoxon(np.array(values, dtype=np.float64))
    assert (result.statistic, result.pvalue) == expected
    assert result.method == "exhaustive_sign_permutation"
    assert result.reason == ""


def test_zero_padding_does_not_change_effective_n_or_small_sample_method() -> None:
    values = np.zeros(548, dtype=np.float64)
    values[[1, 101, 201, 301, 401]] = [1, 2, 3, 4, 5]
    original = values.copy()
    result = study_wilcoxon(values)
    assert result.statistic == 0
    assert result.pvalue == 0.0625
    assert result.method == "exact"
    assert result.nonzero_studies == 5
    assert result.zero_studies == 543
    np.testing.assert_array_equal(values, original)


def test_tolerance_boundary_is_inclusive_and_next_float_is_active() -> None:
    above = np.nextafter(ZERO_TOLERANCE, np.inf)
    values = np.array([0, ZERO_TOLERANCE, -ZERO_TOLERANCE, above, -above], dtype=np.float64)
    original = values.copy()
    result = study_wilcoxon(values)
    assert result.nonzero_studies == 2
    assert result.zero_studies == 3
    assert result.statistic == 1.5
    assert result.pvalue == 1
    assert result.method == "exhaustive_sign_permutation"
    np.testing.assert_array_equal(values, original)


def test_distinct_adjacent_nonzero_magnitudes_are_not_rounded_into_ties() -> None:
    values = np.array([1.0, np.nextafter(1.0, np.inf)], dtype=np.float64)
    result = study_wilcoxon(values)
    assert result.method == "exact"
    assert result.statistic == 0
    assert result.pvalue == 0.5


@pytest.mark.parametrize(("length", "method"), [(50, "exact"), (51, "asymptotic")])
def test_untied_cutoff_uses_effective_support(length: int, method: str) -> None:
    values = np.arange(1, length + 1, dtype=np.float64)
    result = study_wilcoxon(values)
    expected = stats.wilcoxon(
        values, method="exact" if length == 50 else "approx", correction=False
    )
    assert result.method == method
    assert result.statistic == expected.statistic
    assert result.pvalue == expected.pvalue


@pytest.mark.parametrize(
    ("length", "method"), [(13, "exhaustive_sign_permutation"), (14, "asymptotic")]
)
def test_tied_cutoff_uses_exhaustive_signs_only_through_thirteen(length: int, method: str) -> None:
    values = np.ones(length, dtype=np.float64)
    result = study_wilcoxon(values)
    expected_p = (
        2 / 2**length
        if length == 13
        else stats.wilcoxon(values, method="approx", correction=False).pvalue
    )
    assert result.method == method
    assert result.statistic == 0
    assert result.pvalue == expected_p


@pytest.mark.parametrize(
    ("values", "reason"),
    [
        ([], "insufficient_support"),
        ([0.0], "all_tied"),
        ([0.0, 1e-12, -1e-12, np.nextafter(1e-12, 0.0)], "all_tied"),
    ],
)
def test_untestable_vectors_retain_internal_one_and_declared_support(
    values: list[float], reason: str
) -> None:
    result = study_wilcoxon(np.array(values, dtype=np.float64))
    assert np.isnan(result.statistic)
    assert result.pvalue == 1
    assert result.nonzero_studies == 0
    assert result.zero_studies == len(values)
    assert result.method == "not_tested"
    assert result.reason == reason


@pytest.mark.parametrize("values", [[np.nan], [np.inf], [-np.inf], [[1.0, 2.0]]])
def test_nonfinite_or_nonscalar_study_observations_fail_closed(values: list[object]) -> None:
    with pytest.raises(ValueError, match=r"finite|one-dimensional"):
        study_wilcoxon(np.asarray(values, dtype=np.float64))
