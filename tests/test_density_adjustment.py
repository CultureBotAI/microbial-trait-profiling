"""Tests for trait correlation, density adjustment, and near-twins (#264, #352)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from bsdbng.density_adjustment import (
    DEFAULT_COLUMNS_PATH,
    DEFAULT_PROFILES_PATH,
    DEFAULT_TAXA_MATRIX_PATH,
    TABLE_1_TRAITS,
    TEST_PAIR_A,
    TWIN_A,
    TWIN_B,
    compute_correlation_metrics,
    compute_effective_tests,
    downweight_traits_by_taxon_degree,
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
    """Verify Li & Ji formula handles negative/roundoff eigenvalues properly (#352)."""
    # Round-off negative and near-zero eigenvalues should not artificially add ~1.0
    evals = np.array([2.0, 1.0, 1e-16, -1e-16])
    m_li_ji, _m_galwey, _m_nyholt = compute_effective_tests(evals)
    assert pytest.approx(m_li_ji, abs=1e-8) == 2.0


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
    with pytest.raises(ValueError, match="Unsupported correlation method"):
        compute_correlation_metrics(x, cols, method="invalid")  # type: ignore[arg-type]


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

    norm = normalize_traits_by_profile_composition(prof_df, ["t1", "t2"])
    assert norm.shape == (2, 2)
    # Row 1 sum is 4.0, normalized row should sum to 1.0 (0.5, 0.5)
    assert pytest.approx(norm[0].sum()) == 1.0
    assert pytest.approx(norm[0, 0]) == 0.5
    assert pytest.approx(norm[0, 1]) == 0.5
    # Row 2 was all zeros, should remain 0
    assert pytest.approx(norm[1].sum()) == 0.0


def test_downweight_traits_by_taxon_degree_legacy_alias(tmp_path: Path) -> None:
    """Verify backward compatibility of downweight_traits_by_taxon_degree alias."""
    taxa_file = tmp_path / "taxa.tsv"
    taxa_file.touch()

    prof_df = pd.DataFrame(
        {
            "study_id": ["s1", "s2"],
            "experiment_id": ["e1", "e2"],
            "inc_dec": ["inc", "inc"],
            "t1": [3.0, 0.0],
            "t2": [1.0, 0.0],
        }
    )

    # Calling with legacy (taxa_path, profiles_df, cols) signature
    dw = downweight_traits_by_taxon_degree(taxa_file, prof_df, ["t1", "t2"])
    assert dw.shape == (2, 2)
    assert pytest.approx(dw[0].sum()) == 1.0
    assert pytest.approx(dw[0, 0]) == 0.75
    assert pytest.approx(dw[0, 1]) == 0.25
    assert pytest.approx(dw[1].sum()) == 0.0


def test_evaluate_conditional_twin_test_vif() -> None:
    """Verify VIF and conditional test behavior on highly correlated twins."""
    rng = np.random.default_rng(99)
    base = rng.normal(size=200)
    # Create near-twins with r ~ 0.95
    t_a = base + rng.normal(scale=0.2, size=200)
    t_b = base + rng.normal(scale=0.2, size=200)

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
        }
    )

    summary = evaluate_conditional_twin_test(df, TWIN_A, TWIN_B)
    assert summary.trait_a == TWIN_A
    assert summary.trait_b == TWIN_B
    assert summary.profile_correlation > 0.90
    # For r > 0.90, VIF = 1 / (1 - r^2) > 5.0
    assert summary.vif_profile > 5.0
    assert summary.n_paired_experiments == 100


def test_evaluate_table1_paired_effects_synthetic() -> None:
    """Verify paired experiment contrast evaluation across Table 1 traits."""
    study_ids = [f"s{i}" for i in range(20) for _ in range(2)]
    exp_ids = [f"e{i}" for i in range(20) for _ in range(2)]
    inc_dec = ["inc", "dec"] * 20

    # Ensure all Table 1 traits exist in test dataframe
    cols = list(TABLE_1_TRAITS.keys())
    data = {c: np.ones(40, dtype=float) for c in cols}
    # For TEST_PAIR_A, make inc > dec
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


@pytest.mark.skipif(
    not (
        DEFAULT_PROFILES_PATH.exists()
        and DEFAULT_COLUMNS_PATH.exists()
        and DEFAULT_TAXA_MATRIX_PATH.exists()
    ),
    reason="Derived data files not available in test environment",
)
def test_full_density_audit_integration(tmp_path: Path) -> None:
    """Verify full density audit run against real derived data and validate core findings."""
    results = run_density_audit(
        profiles_path=DEFAULT_PROFILES_PATH,
        columns_path=DEFAULT_COLUMNS_PATH,
        taxa_matrix_path=DEFAULT_TAXA_MATRIX_PATH,
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

    # 2. Key findings from Issue #264:
    # Baseline PC1 explains ~60.4% - 60.8% of variance
    raw = results["raw_pearson"]
    assert 58.0 <= raw.pc1_variance_pct <= 63.0
    assert raw.median_pairwise_r > 0.50

    # Post-PC1 regression pairwise median correlation centers near 0
    adj = results["adj_pc1_pearson"]
    assert -0.10 <= adj.median_pairwise_r <= 0.05

    # Near-twin collinearity (coccus vs sphere) shows high VIF
    twin = results["twin_summary"]
    assert twin.profile_correlation > 0.85
    assert twin.vif_profile > 4.0

    # Diagnostic pair: facultatively aerobic x moderately halophilic
    t1_effects = results["table1_effects"]
    fac_aero = next(e for e in t1_effects if e.trait_id == TEST_PAIR_A)
    # Post-adjustment p-value should collapse to nonsignificance (> 0.05)
    assert fac_aero.adj_wilcoxon_p > 0.05
    assert "collapsed" in fac_aero.status
