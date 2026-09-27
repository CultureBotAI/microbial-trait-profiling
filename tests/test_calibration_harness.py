"""Tests for statistical calibration harness and Ilyas's ORA adapter."""

import numpy as np
import pandas as pd
import pytest

from bsdbng.calibration.adapters.ilyas_ora import IlyasORAAdapter
from bsdbng.calibration.harness import run_ora_calibration_benchmark
from bsdbng.calibration.metrics import (
    compute_discrete_ecdf,
    compute_empirical_fdr,
    compute_type_i_error,
)
from bsdbng.calibration.scenarios import (
    generate_complete_null_ora_profiles,
    generate_mixed_signal_ora_profiles,
)


@pytest.fixture
def synthetic_calibration_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    # 20 Taxa with various trait profiles
    taxa_data = {f"trait_{i}": [1 if (j + i) % 3 == 0 else 0 for j in range(20)] for i in range(5)}
    taxa_df = pd.DataFrame(taxa_data, index=[f"taxon_{j}" for j in range(20)])

    # 4 empirical study profiles
    profiles_data = {
        "study_id": ["s1", "s1", "s2", "s2"],
        "experiment_id": ["e1", "e1", "e2", "e2"],
        "inc_dec": ["inc", "dec", "inc", "dec"],
        "condition": ["IBD", "IBD", "CRC", "CRC"],
        "body_site": ["feces", "feces", "feces", "feces"],
        "signature_size": [6, 6, 8, 8],
        "trait_0": [2, 2, 3, 3],
        "trait_1": [2, 1, 2, 2],
        "trait_2": [1, 2, 2, 3],
        "trait_3": [2, 2, 3, 2],
        "trait_4": [1, 1, 2, 2],
    }
    profiles_df = pd.DataFrame(profiles_data)
    return profiles_df, taxa_df


def test_type_i_error_metric() -> None:
    # Super-uniform p-values
    p_uniform = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    metric = compute_type_i_error(p_uniform, alpha=0.05)
    assert metric.empirical_rate == 0.0
    assert metric.is_conservative is True

    # Uniform at alpha=0.10
    metric_10 = compute_type_i_error(p_uniform, alpha=0.10)
    assert metric_10.empirical_rate == 0.10
    assert metric_10.is_conservative is True


def test_empirical_fdr_metric() -> None:
    discoveries = [
        {("p1", "t1"), ("p1", "t2")},
        {("p1", "t1")},
        set(),
    ]
    true_signals = {("p1", "t1")}

    fdr_metric = compute_empirical_fdr(discoveries, true_signals, target_q=0.05)
    # Replicate 1: 1 FP / 2 Disc = 0.5
    # Replicate 2: 0 FP / 1 Disc = 0.0
    # Replicate 3: 0 Disc = 0.0
    # Mean FDP = (0.5 + 0.0 + 0.0) / 3 = 0.1667
    assert pytest.approx(fdr_metric.empirical_fdr, abs=1e-3) == 0.167
    assert fdr_metric.replicates == 3


def test_discrete_ecdf_metric() -> None:
    p_vals = [0.01, 0.05, 0.2, 0.5, 0.9]
    ecdf_df = compute_discrete_ecdf(p_vals, bins=10)
    assert not ecdf_df.empty
    assert "empirical_cdf" in ecdf_df.columns
    assert "uniform_cdf" in ecdf_df.columns


def test_complete_null_scenario_generator(
    synthetic_calibration_data: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    profiles_df, taxa_df = synthetic_calibration_data
    rng = np.random.default_rng(42)

    null_sim = generate_complete_null_ora_profiles(profiles_df, taxa_df, rng)
    assert len(null_sim) == len(profiles_df)
    assert "signature_size" in null_sim.columns
    for trait in ["trait_0", "trait_1", "trait_2", "trait_3", "trait_4"]:
        assert (null_sim[trait] >= 0).all()
        assert (null_sim[trait] <= null_sim["signature_size"]).all()


def test_mixed_signal_scenario_generator(
    synthetic_calibration_data: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    profiles_df, taxa_df = synthetic_calibration_data
    rng = np.random.default_rng(42)

    mixed_sim, true_pairs = generate_mixed_signal_ora_profiles(
        profiles_df, taxa_df, spike_traits=["trait_0"], signal_strength=0.5, rng=rng
    )
    assert len(mixed_sim) == len(profiles_df)
    assert len(true_pairs) > 0
    assert all(pair[1] == "trait_0" for pair in true_pairs)


def test_ora_adapter_and_harness_benchmark(
    synthetic_calibration_data: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    profiles_df, taxa_df = synthetic_calibration_data
    adapter = IlyasORAAdapter()

    # Fast deterministic benchmark with B=20 replicates
    report = run_ora_calibration_benchmark(
        adapter=adapter,
        empirical_profiles_df=profiles_df,
        taxa_matrix_df=taxa_df,
        spike_traits=["trait_0"],
        replicates=20,
        seed=123,
        alpha=0.05,
    )

    assert report.method_id == "ilyas_ora_hypergeometric"
    assert report.replicates == 20
    # Type I error under complete null must be conservative / valid
    assert report.complete_null_type_i.is_conservative is True
    assert report.complete_null_type_i.empirical_rate <= 0.15
    assert not report.complete_null_ecdf.empty
