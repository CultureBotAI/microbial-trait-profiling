"""Tests for trait correlation, density adjustment, and near-twins (#264, #352)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

from bsdbng.density_adjustment import (
    DEFAULT_COLUMNS_PATH,
    DEFAULT_PROFILES_PATH,
    TABLE_1_TRAITS,
    TEST_PAIR_A,
    TWIN_A,
    TWIN_B,
    compute_correlation_metrics,
    compute_effective_tests,
    evaluate_conditional_twin_test,
    evaluate_table1_paired_effects,
    get_non_media_trait_columns,
    load_non_media_profile_matrix,
    normalize_traits_by_profile_composition,
    regress_out_pc1,
    run_density_audit,
)


def test_effective_tests_identity_and_collinear() -> None:
    """Verify effective test calculations on orthogonal and perfectly collinear structures."""
    # Orthogonal / uncorrelated case: 5 traits with equal eigenvalues of 1.0
    evals_orth = np.array([1.0, 1.0, 1.0, 1.0, 1.0])
    m_li_ji, m_galwey, m_nyholt = compute_effective_tests(evals_orth)
    assert pytest.approx(m_li_ji) == 5.0
    assert pytest.approx(m_galwey) == 5.0
    assert pytest.approx(m_nyholt) == 5.0

    # Perfectly collinear case: 5 traits with all variance in one component (5.0, 0, 0, 0, 0)
    evals_collinear = np.array([5.0, 0.0, 0.0, 0.0, 0.0])
    m_li_ji_col, m_galwey_col, m_nyholt_col = compute_effective_tests(evals_collinear)
    assert pytest.approx(m_li_ji_col) == 1.0
    assert pytest.approx(m_galwey_col) == 1.0
    assert pytest.approx(m_nyholt_col) == 1.0

    # Empty array
    assert compute_effective_tests(np.array([])) == (0.0, 0.0, 0.0)


def test_effective_tests_roundoff_eigenvalues() -> None:
    """Verify Li & Ji formula handles negative/roundoff eigenvalues properly (#352, #377)."""
    # Round-off negative and near-zero eigenvalues should not artificially add ~1.0
    evals = np.array([2.0, 1.0, 1e-11, -1e-11])
    m_li_ji, _m_galwey, _m_nyholt = compute_effective_tests(evals)
    assert pytest.approx(m_li_ji, abs=1e-8) == 2.0

    # Upper cap at P (total eigenvalues)
    evals_equal = np.array([1.0, 1.0, 1.0])
    m_li_ji_eq, _, _ = compute_effective_tests(evals_equal)
    assert pytest.approx(m_li_ji_eq, abs=1e-8) == 3.0
    assert m_li_ji_eq <= 3.0


def test_effective_tests_constant_columns() -> None:
    """Verify constant columns are dropped and do not inflate Li & Ji effective tests (#378)."""
    rng = np.random.default_rng(42)
    # 3 informative columns + 2 constant columns
    data_var = rng.normal(size=(50, 3))
    data_const1 = np.zeros((50, 1))
    data_const2 = np.full((50, 1), 7.0)
    data = np.hstack([data_var, data_const1, data_const2])
    cols = ["v1", "v2", "v3", "c1", "c2"]

    metrics = compute_correlation_metrics(data, cols, method="pearson")
    # Li & Ji effective tests should not exceed the number of variable traits (3)
    assert metrics.m_eff_li_ji <= 3.0 + 1e-7


def test_get_non_media_trait_columns(tmp_path: Path) -> None:
    """Verify that MediaDive traits (METPO:1004005) are strictly excluded."""
    cols_file = tmp_path / "columns.tsv"
    cols_data = pd.DataFrame(
        {
            "column_id": ["trait_1", "trait_2", "media_1", "trait_3"],
            "object_category": [
                "METPO:1000666",
                "METPO:1000623",
                "METPO:1004005",
                "METPO:2000103",
            ],
        }
    )
    cols_data.to_csv(cols_file, sep="\t", index=False)

    non_media = get_non_media_trait_columns(cols_file)
    assert non_media == ["trait_1", "trait_2", "trait_3"]
    assert "media_1" not in non_media


def test_load_non_media_profile_matrix(tmp_path: Path) -> None:
    """Verify that wide profile matrix is loaded and restricted to non-media traits."""
    cols_file = tmp_path / "columns.tsv"
    pd.DataFrame(
        {
            "column_id": ["t1", "t2"],
            "object_category": ["cat1", "cat2"],
        }
    ).to_csv(cols_file, sep="\t", index=False)

    prof_file = tmp_path / "profiles.tsv"
    pd.DataFrame(
        {
            "study_id": ["s1", "s1"],
            "experiment_id": ["e1", "e1"],
            "inc_dec": ["inc", "dec"],
            "t1": [5.0, 1.0],
            "t2": [10.0, 2.0],
        }
    ).to_csv(prof_file, sep="\t", index=False)

    df, cols = load_non_media_profile_matrix(prof_file, cols_file)
    assert cols == ["t1", "t2"]
    assert len(df) == 2

    # Test error if a column is missing
    bad_cols_file = tmp_path / "bad_columns.tsv"
    pd.DataFrame(
        {
            "column_id": ["t1", "missing_trait"],
            "object_category": ["cat1", "cat2"],
        }
    ).to_csv(bad_cols_file, sep="\t", index=False)

    with pytest.raises(ValueError, match="missing 1 expected trait columns"):
        load_non_media_profile_matrix(prof_file, bad_cols_file)


def test_regress_out_pc1_orthogonality() -> None:
    """Verify that regress_out_pc1 produces residuals strictly orthogonal to PC1."""
    rng = np.random.default_rng(42)
    # Generate 100 samples with 1 latent density factor driving 4 traits
    density = rng.normal(size=100)
    trait_matrix = np.column_stack(
        [
            2.0 * density + rng.normal(scale=0.5, size=100),
            1.5 * density + rng.normal(scale=0.5, size=100),
            3.0 * density + rng.normal(scale=0.5, size=100),
            rng.normal(scale=1.0, size=100),
        ]
    )

    residuals, pc1_scores, pc1_var_pct = regress_out_pc1(trait_matrix)
    assert residuals.shape == trait_matrix.shape
    assert pc1_var_pct > 50.0

    # Orthogonality check: dot product between residuals and centered PC1 should be ~0
    pc1_c = pc1_scores - pc1_scores.mean()
    for col_idx in range(residuals.shape[1]):
        dot_prod = np.dot(residuals[:, col_idx], pc1_c)
        assert pytest.approx(dot_prod, abs=1e-8) == 0.0


def test_compute_correlation_metrics_synthetic() -> None:
    """Verify compute_correlation_metrics summary metrics and upper triangle logic."""
    rng = np.random.default_rng(123)
    x = rng.normal(size=(50, 4))
    cols = ["col1", "col2", "col3", "col4"]

    metrics = compute_correlation_metrics(x, cols, method="pearson")
    assert metrics.n_samples == 50
    assert metrics.n_traits == 4
    assert -1.0 <= metrics.median_pairwise_r <= 1.0
    assert -1.0 <= metrics.mean_pairwise_r <= 1.0
    assert 0.0 <= metrics.pct_pairs_gt_05 <= 100.0
    assert 0.0 <= metrics.pct_pairs_gt_09 <= 100.0
    assert 0.0 < metrics.pc1_variance_pct <= 100.0
    assert 0.0 < metrics.m_eff_li_ji <= 4.0 + 1e-7

    # Test spearman method
    metrics_sp = compute_correlation_metrics(x, cols, method="spearman")
    assert metrics_sp.n_samples == 50
    assert metrics_sp.n_traits == 4

    # Test invalid method
    bad_method: Any = "invalid"
    with pytest.raises(ValueError, match="Unsupported correlation method"):
        compute_correlation_metrics(x, cols, method=bad_method)


def test_normalize_traits_by_profile_composition() -> None:
    """Verify profile normalization by row sum composition (Method B / TSS)."""
    prof_df = pd.DataFrame(
        {
            "study_id": ["s1", "s2"],
            "experiment_id": ["e1", "e2"],
            "inc_dec": ["inc", "inc"],
            "t1": [2.0, 0.0],
            "t2": [2.0, 0.0],
        }
    )

    # Test with Sequence[str]: list, tuple, pd.Index (#379)
    for cols in [["t1", "t2"], ("t1", "t2"), pd.Index(["t1", "t2"])]:
        norm = normalize_traits_by_profile_composition(prof_df, cols)
        assert norm.shape == (2, 2)
        # Row 1 sum is 4.0, normalized row should sum to 1.0 (0.5, 0.5)
        assert pytest.approx(norm[0].sum()) == 1.0
        assert pytest.approx(norm[0, 0]) == 0.5
        assert pytest.approx(norm[0, 1]) == 0.5
        # Row 2 was all zeros, should remain 0
        assert pytest.approx(norm[1].sum()) == 0.0

    # Test error if cols is empty
    with pytest.raises(ValueError, match="cols must not be empty"):
        normalize_traits_by_profile_composition(prof_df, [])


def test_evaluate_conditional_twin_test_vif() -> None:
    """Verify VIF and conditional test behavior on highly correlated twins (#380)."""
    rng = np.random.default_rng(99)
    depth = rng.uniform(50.0, 100.0, size=200)
    base = rng.normal(size=200)
    # Create near-twins with r ~ 0.95 and correlated with depth
    t_a = 0.5 * depth + base + rng.normal(scale=0.2, size=200)
    t_b = 0.5 * depth + base + rng.normal(scale=0.2, size=200)
    t_other = depth + rng.normal(scale=1.0, size=200)

    # 100 paired experiments (inc/dec)
    study_ids = [f"s_{i // 2}" for i in range(200)]
    exp_ids = [f"e_{i // 2}" for i in range(200)]
    inc_dec = ["inc", "dec"] * 100

    df = pd.DataFrame(
        {
            "study_id": study_ids,
            "experiment_id": exp_ids,
            "inc_dec": inc_dec,
            TWIN_A: t_a,
            TWIN_B: t_b,
            "t_other": t_other,
        }
    )

    summary = evaluate_conditional_twin_test(
        df, TWIN_A, TWIN_B, columns_for_depth=[TWIN_A, TWIN_B, "t_other"]
    )
    assert summary.trait_a == TWIN_A
    assert summary.trait_b == TWIN_B
    assert summary.profile_correlation > 0.90
    # For r > 0.90, VIF = 1 / (1 - r^2) > 5.0
    assert summary.vif_profile > 5.0
    assert summary.n_paired_experiments == 100
    assert summary.profile_partial_r_depth > 0.80
    assert summary.paired_slope_pvalue < 0.05
    assert summary.paired_slope_a_on_b > 0.0


def test_evaluate_table1_paired_effects_synthetic() -> None:
    """Verify paired experiment contrast evaluation across Table 1 traits (#383)."""
    study_ids = [f"s{i}" for i in range(20) for _ in range(2)]
    exp_ids = [f"e{i}" for i in range(20) for _ in range(2)]
    inc_dec = ["inc", "dec"] * 20

    # Ensure all Table 1 traits exist in test dataframe
    cols = list(TABLE_1_TRAITS.keys())
    data = {c: np.ones(40, dtype=float) for c in cols}
    # For TEST_PAIR_A, make inc > dec with strong signal
    inc_vals = np.array([5.0, 1.0] * 20)
    data[TEST_PAIR_A] = inc_vals

    prof_df = pd.DataFrame(
        {
            "study_id": study_ids,
            "experiment_id": exp_ids,
            "inc_dec": inc_dec,
            **data,
        }
    )

    adj_matrix = prof_df[cols].values.copy()
    effects = evaluate_table1_paired_effects(prof_df, adj_matrix, cols)
    assert len(effects) == len(cols)
    test_pair_eff = next(e for e in effects if e.trait_id == TEST_PAIR_A)
    assert test_pair_eff.raw_mean_diff > 0
    assert test_pair_eff.raw_wilcoxon_p < 0.05
    assert "survived" in test_pair_eff.status


def test_evaluate_table1_paired_effects_nan_fallback() -> None:
    """Verify Wilcoxon p-value falls back to NaN when non-zero pairs <= 10 (#383)."""
    study_ids = [f"s{i}" for i in range(5) for _ in range(2)]
    exp_ids = [f"e{i}" for i in range(5) for _ in range(2)]
    inc_dec = ["inc", "dec"] * 5

    cols = list(TABLE_1_TRAITS.keys())
    data = {c: np.zeros(10, dtype=float) for c in cols}
    # Set only 5 pairs for TEST_PAIR_A (<= 10 non-zero pairs)
    data[TEST_PAIR_A] = np.array([2.0, 1.0] * 5)

    prof_df = pd.DataFrame(
        {
            "study_id": study_ids,
            "experiment_id": exp_ids,
            "inc_dec": inc_dec,
            **data,
        }
    )

    adj_matrix = prof_df[cols].values.copy()
    effects = evaluate_table1_paired_effects(prof_df, adj_matrix, cols)
    test_pair_eff = next(e for e in effects if e.trait_id == TEST_PAIR_A)
    assert np.isnan(test_pair_eff.raw_wilcoxon_p)
    assert np.isnan(test_pair_eff.adj_wilcoxon_p)
    assert test_pair_eff.status == "insufficient non-zero pairs"


@pytest.mark.skipif(
    not (
        (DEFAULT_PROFILES_PATH.exists() or Path(f"{DEFAULT_PROFILES_PATH}.gz").exists())
        and DEFAULT_COLUMNS_PATH.exists()
    ),
    reason="Derived data files not available in test environment",
)
def test_full_density_audit_integration(tmp_path: Path) -> None:
    """Verify full density audit run against real derived data and validate core findings (#375)."""
    prof_path = (
        DEFAULT_PROFILES_PATH
        if DEFAULT_PROFILES_PATH.exists()
        else Path(f"{DEFAULT_PROFILES_PATH}.gz")
    )
    results = run_density_audit(
        profiles_path=prof_path,
        columns_path=DEFAULT_COLUMNS_PATH,
        output_dir=tmp_path,
    )

    # 1. Output files exist and are populated
    summary_tsv = tmp_path / "density_adjustment_summary.tsv"
    table1_tsv = tmp_path / "table1_density_adjusted_effects.tsv"
    report_txt = tmp_path / "density_adjustment_report.txt"

    assert summary_tsv.exists()
    assert table1_tsv.exists()
    assert report_txt.exists()

    summary_df = pd.read_csv(summary_tsv, sep="\t")
    assert len(summary_df) == 4
    assert set(summary_df["metric_set"]) == {
        "raw_pearson",
        "raw_spearman",
        "adjusted_pc1_pearson",
        "adjusted_profile_composition_pearson",
    }

    # 2. Structural sanity checks on real data
    raw = results["raw_pearson"]
    assert 0.0 < raw.pc1_variance_pct <= 100.0
    assert -1.0 <= raw.median_pairwise_r <= 1.0
    assert raw.n_traits > 0
    assert raw.n_samples > 0

    adj = results["adj_pc1_pearson"]
    assert -1.0 <= adj.median_pairwise_r <= 1.0

    twin = results["twin_summary"]
    assert -1.0 <= twin.profile_correlation <= 1.0
    assert twin.vif_profile >= 1.0

    t1_effects = results["table1_effects"]
    assert len(t1_effects) == len(TABLE_1_TRAITS)
    for eff in t1_effects:
        assert eff.trait_id in TABLE_1_TRAITS
        assert eff.status in {
            "insufficient non-zero pairs",
            "collapsed (attenuated to nonsignificant)",
            "enhanced (significant, larger magnitude)",
            "survived (significant post-adjustment)",
            "attenuated (nonsignificant post-adjustment)",
        }


def test_compute_effective_tests_non_positive_evals() -> None:
    """Verify compute_effective_tests handles zero or all-negative eigenvalues safely."""
    # Zero eigenvalues: should return (0.0, 0.0, float(p)) without ZeroDivisionError
    res_zero = compute_effective_tests([0.0, 0.0])
    assert res_zero == (0.0, 0.0, 2.0)

    # Negative eigenvalues
    res_neg = compute_effective_tests([-1.0, -2.0])
    assert res_neg == (0.0, 0.0, 2.0)

    # Single non-positive eigenvalue
    res_single = compute_effective_tests([0.0])
    assert res_single == (0.0, 0.0, 1.0)


def test_regress_out_pc1_degenerate_shapes() -> None:
    """Verify regress_out_pc1 handles matrices with < 2 rows or 0 columns safely."""
    # Single sample
    one_row = np.array([[1.0, 2.0, 3.0]])
    res_one, scores_one, var_one = regress_out_pc1(one_row)
    assert res_one.shape == (1, 3)
    assert len(scores_one) == 1
    assert var_one == 0.0

    # Zero samples
    zero_row = np.empty((0, 3))
    res_zero, scores_zero, var_zero = regress_out_pc1(zero_row)
    assert res_zero.shape == (0, 3)
    assert len(scores_zero) == 0
    assert var_zero == 0.0

    # Constant matrix (zero variance across all samples)
    const_mat = np.full((10, 4), 5.0)
    res_const, scores_const, var_const = regress_out_pc1(const_mat)
    assert res_const.shape == (10, 4)
    assert len(scores_const) == 10
    assert var_const == 0.0
    assert np.allclose(res_const, const_mat)


def test_evaluate_table1_paired_effects_degenerate() -> None:
    """Verify evaluate_table1_paired_effects handles empty or single-direction DataFrames safely."""
    cols = list(TABLE_1_TRAITS.keys())

    # Empty DataFrame
    empty_df = pd.DataFrame(columns=["study_id", "experiment_id", "inc_dec", *cols])
    adj_empty = np.empty((0, len(cols)))
    effects_empty = evaluate_table1_paired_effects(empty_df, adj_empty, cols)
    assert len(effects_empty) == len(cols)
    for eff in effects_empty:
        assert eff.raw_mean_diff == 0.0
        assert np.isnan(eff.raw_wilcoxon_p)
        assert eff.status == "insufficient non-zero pairs"

    # Single-direction DataFrame (inc only, no dec)
    single_dir_df = pd.DataFrame(
        {
            "study_id": ["s1", "s2"],
            "experiment_id": ["e1", "e2"],
            "inc_dec": ["inc", "inc"],
            **{c: [1.0, 2.0] for c in cols},
        }
    )
    adj_single = single_dir_df[cols].to_numpy(dtype=float)
    effects_single = evaluate_table1_paired_effects(single_dir_df, adj_single, cols)
    assert len(effects_single) == len(cols)
    for eff in effects_single:
        assert eff.status == "insufficient non-zero pairs"


def test_evaluate_conditional_twin_test_constant_and_unpaired() -> None:
    """Verify evaluate_conditional_twin_test handles constant traits and unpaired data."""
    # Constant traits: correlation is NaN, VIF is NaN (not Inf)
    const_df = pd.DataFrame(
        {
            "study_id": ["s1", "s1"],
            "experiment_id": ["e1", "e1"],
            "inc_dec": ["inc", "dec"],
            TWIN_A: [5.0, 5.0],
            TWIN_B: [5.0, 5.0],
        }
    )
    summary_const = evaluate_conditional_twin_test(const_df, TWIN_A, TWIN_B)
    assert np.isnan(summary_const.profile_correlation)
    assert np.isnan(summary_const.vif_profile)
    assert np.isnan(summary_const.profile_partial_r_depth)
    assert summary_const.profile_slope_a_on_b == 0.0
    assert np.isnan(summary_const.paired_slope_pvalue)

    # Single direction data (no paired contrasts)
    unpaired_df = pd.DataFrame(
        {
            "study_id": ["s1", "s2"],
            "experiment_id": ["e1", "e2"],
            "inc_dec": ["inc", "inc"],
            TWIN_A: [1.0, 2.0],
            TWIN_B: [3.0, 4.0],
        }
    )
    summary_unpaired = evaluate_conditional_twin_test(unpaired_df, TWIN_A, TWIN_B)
    assert summary_unpaired.n_paired_experiments == 0
    assert np.isnan(summary_unpaired.paired_correlation)
    assert np.isnan(summary_unpaired.vif_paired)
    assert summary_unpaired.paired_slope_a_on_b == 0.0
    assert np.isnan(summary_unpaired.paired_slope_pvalue)


def test_evaluate_table1_paired_effects_missing_keys_and_mismatched_shape() -> None:
    """Verify evaluate_table1_paired_effects handles missing IDs or shape mismatches safely."""
    cols = list(TABLE_1_TRAITS.keys())

    # Missing study_id column
    df_no_study = pd.DataFrame(
        {
            "experiment_id": ["e1", "e2"],
            "inc_dec": ["inc", "dec"],
            **{c: [1.0, 2.0] for c in cols},
        }
    )
    adj_mat = np.ones((2, len(cols)))
    effects_no_study = evaluate_table1_paired_effects(df_no_study, adj_mat, cols)
    assert len(effects_no_study) == len(cols)
    assert all(e.status == "insufficient non-zero pairs" for e in effects_no_study)

    # Mismatched row count (e.g. 3 rows in adj_mat, 2 rows in df)
    df_valid = pd.DataFrame(
        {
            "study_id": ["s1", "s1"],
            "experiment_id": ["e1", "e1"],
            "inc_dec": ["inc", "dec"],
            **{c: [1.0, 2.0] for c in cols},
        }
    )
    adj_bad_rows = np.ones((3, len(cols)))
    effects_bad_rows = evaluate_table1_paired_effects(df_valid, adj_bad_rows, cols)
    assert len(effects_bad_rows) == len(cols)
    assert all(e.status == "insufficient non-zero pairs" for e in effects_bad_rows)


def test_evaluate_conditional_twin_test_unknown_depth_columns_and_missing_keys() -> None:
    """Verify evaluate_conditional_twin_test handles unknown depth columns and missing IDs."""
    # Unknown depth columns should safely fall back to vec_a + vec_b
    df = pd.DataFrame(
        {
            "study_id": ["s1", "s1"],
            "experiment_id": ["e1", "e1"],
            "inc_dec": ["inc", "dec"],
            TWIN_A: [2.0, 4.0],
            TWIN_B: [1.0, 3.0],
        }
    )
    summary = evaluate_conditional_twin_test(
        df, TWIN_A, TWIN_B, columns_for_depth=["non_existent_col_1", "non_existent_col_2"]
    )
    assert summary.n_paired_experiments == 1
    assert summary.trait_a == TWIN_A
    assert summary.trait_b == TWIN_B
    # Depth fell back to vec_a + vec_b, standard deviation is > 0, partial_r evaluated
    assert not np.isnan(summary.profile_correlation)

    # Missing study_id in df
    df_no_id = pd.DataFrame(
        {
            TWIN_A: [2.0, 4.0],
            TWIN_B: [1.0, 3.0],
            "inc_dec": ["inc", "dec"],
        }
    )
    summary_no_id = evaluate_conditional_twin_test(df_no_id, TWIN_A, TWIN_B)
    assert summary_no_id.n_paired_experiments == 0
    assert np.isnan(summary_no_id.paired_correlation)

