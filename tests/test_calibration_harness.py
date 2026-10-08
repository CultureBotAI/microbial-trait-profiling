"""Tests for the statistical calibration harness and Ilyas's ORA adapter."""

from __future__ import annotations

import json
import subprocess
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

from bsdbng.calibration.adapters.ilyas_ora import IlyasORAAdapter
from bsdbng.calibration.harness import (
    EXIT_FAILURE,
    EXIT_SUCCESS,
    NULL_PVALUES_FILENAME,
    REPLICATES_FILENAME,
    SUMMARY_FILENAME,
    run_ora_calibration_benchmark,
)
from bsdbng.calibration.harness import main as harness_main
from bsdbng.calibration.metrics import (
    compute_discrete_ecdf,
    compute_empirical_fdr,
    compute_type_i_error,
)
from bsdbng.calibration.protocol import MethodInferenceResult
from bsdbng.calibration.scenarios import (
    generate_complete_null_ora_profiles,
    generate_mixed_signal_ora_profiles,
    resolve_signature_sizes,
)
from bsdbng.enrichment_stats import benjamini_hochberg

TRAITS = [f"METPO:{1000 + i}" for i in range(5)]
N_TAXA = 30
CalibrationData = tuple[pd.DataFrame, pd.DataFrame]


@pytest.fixture
def synthetic_calibration_data() -> CalibrationData:
    """Profiles and taxa matrix in the same layout the pipeline writes."""
    trait_values = {
        trait: [1 if (j + i) % 3 == 0 else 0 for j in range(N_TAXA)]
        for i, trait in enumerate(TRAITS)
    }
    taxa_df = pd.DataFrame(
        {
            "taxon_name": [f"Taxon {j}" for j in range(N_TAXA)],
            "taxonomic_rank": ["species"] * N_TAXA,
            **trait_values,
        },
        index=pd.Index([f"NCBITaxon:{j}" for j in range(N_TAXA)], name="taxon_id"),
    )
    profiles_df = pd.DataFrame(
        {
            "study_id": ["s1", "s1", "s2", "s2", "s3", "s3"],
            "experiment_id": ["e1", "e1", "e2", "e2", "e3", "e3"],
            "inc_dec": ["increased", "decreased"] * 3,
            "taxa_count": [6, 6, 8, 8, 5, 0],
            **{trait: [2, 2, 3, 2, 1, 0] for trait in TRAITS},
        }
    )
    return profiles_df, taxa_df


def _write_inputs(tmp_path: Path, data: CalibrationData) -> tuple[Path, Path]:
    profiles_df, taxa_df = data
    profiles_path = tmp_path / "profiles.tsv"
    taxa_path = tmp_path / "taxa.tsv"
    profiles_df.to_csv(profiles_path, sep="\t", index=False)
    taxa_df.to_csv(taxa_path, sep="\t")
    return profiles_path, taxa_path


def _cli_args(profiles: Path, taxa: Path, out: Path, *extra: str) -> list[str]:
    return [
        "--profiles",
        str(profiles),
        "--taxa-matrix",
        str(taxa),
        "--output-dir",
        str(out),
        "--replicates",
        "4",
        "--seed",
        "7",
        "--alpha",
        "0.05",
        *extra,
    ]


def _read_summary(out: Path) -> dict[str, Any]:
    summary: dict[str, Any] = json.loads((out / SUMMARY_FILENAME).read_text(encoding="utf-8"))
    return summary


# --- metrics -----------------------------------------------------------------


def test_type_i_error_metric() -> None:
    p_uniform = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    assert compute_type_i_error(p_uniform, alpha=0.05).empirical_rate == 0.0
    assert compute_type_i_error(p_uniform, alpha=0.10).empirical_rate == pytest.approx(0.10)


def test_empirical_fdr_metric() -> None:
    discoveries = [{("p1", "t1"), ("p1", "t2")}, {("p1", "t1")}, set()]
    fdr_metric = compute_empirical_fdr(discoveries, {("p1", "t1")}, target_q=0.05)
    assert fdr_metric.empirical_fdr == pytest.approx(0.5 / 3)
    assert fdr_metric.replicates == 3


def test_discrete_ecdf_metric() -> None:
    ecdf_df = compute_discrete_ecdf([0.01, 0.05, 0.2, 0.5, 0.9], bins=10)
    assert {"empirical_cdf", "uniform_cdf"} <= set(ecdf_df.columns)


def test_benjamini_hochberg_matches_reference() -> None:
    # Sorted p = [0.01, 0.03, 0.04, 0.20]; raw p*m/rank = [0.04, 0.06, 0.0533, 0.20];
    # step-up cumulative minimum gives [0.04, 0.0533, 0.0533, 0.20].
    p_values = np.array([0.01, 0.04, 0.03, 0.20])
    np.testing.assert_allclose(
        benjamini_hochberg(p_values), [0.04, 0.16 / 3, 0.16 / 3, 0.20], rtol=1e-12
    )


# --- scenarios ---------------------------------------------------------------


def test_signature_size_has_no_silent_fallback(
    synthetic_calibration_data: CalibrationData,
) -> None:
    profiles_df, _ = synthetic_calibration_data
    with pytest.raises(ValueError, match="no signature size column"):
        resolve_signature_sizes(profiles_df.drop(columns="taxa_count"))
    bad = profiles_df.assign(taxa_count=["6", "x", "8", "8", "5", "0"])
    with pytest.raises(ValueError, match="non-numeric"):
        resolve_signature_sizes(bad)


def test_complete_null_scenario_respects_signature_sizes(
    synthetic_calibration_data: CalibrationData,
) -> None:
    profiles_df, taxa_df = synthetic_calibration_data
    null_sim = generate_complete_null_ora_profiles(profiles_df, taxa_df, np.random.default_rng(1))
    assert len(null_sim) == len(profiles_df)
    assert null_sim["taxa_count"].tolist() == profiles_df["taxa_count"].tolist()
    for trait in TRAITS:
        assert (null_sim[trait] >= 0).all()
        assert (null_sim[trait] <= null_sim["taxa_count"]).all()
    assert (null_sim.loc[null_sim["taxa_count"] == 0, TRAITS] == 0).all().all()


def test_mixed_signal_scenario_generator(synthetic_calibration_data: CalibrationData) -> None:
    profiles_df, taxa_df = synthetic_calibration_data
    _, true_pairs = generate_mixed_signal_ora_profiles(
        profiles_df,
        taxa_df,
        spike_traits=[TRAITS[0]],
        signal_strength=0.5,
        rng=np.random.default_rng(1),
    )
    assert true_pairs
    assert all(trait == TRAITS[0] for _, trait in true_pairs)
    # The zero-size profile cannot carry a signal.
    assert ("s3:e3:decreased", TRAITS[0]) not in true_pairs


def test_mixed_signal_rejects_unknown_spike_trait(
    synthetic_calibration_data: CalibrationData,
) -> None:
    profiles_df, taxa_df = synthetic_calibration_data
    with pytest.raises(ValueError, match="Spike traits not found"):
        generate_mixed_signal_ora_profiles(profiles_df, taxa_df, spike_traits=["METPO:nope"])


# --- benchmark ---------------------------------------------------------------


def test_benchmark_records_every_replicate(synthetic_calibration_data: CalibrationData) -> None:
    profiles_df, taxa_df = synthetic_calibration_data
    report = run_ora_calibration_benchmark(
        adapter=IlyasORAAdapter(),
        empirical_profiles_df=profiles_df,
        taxa_matrix_df=taxa_df,
        spike_traits=[TRAITS[0]],
        replicates=6,
        seed=123,
    )
    assert report.method_id == "ilyas_ora_hypergeometric"
    assert report.status == "SUCCESS"
    assert [rec.replicate_index for rec in report.replicate_records] == list(range(6))
    assert report.complete_null_type_i is not None
    assert report.complete_null_type_i.replicates == 6
    assert report.complete_null_type_i.n_tests == len(report.null_pvalues)
    assert len(report.null_pvalues) > 0
    assert report.mixed_signal_fdr is not None
    assert report.mixed_signal_fdr.replicates == 6
    assert report.adapter_metadata["test_scope"]


def test_benchmark_is_reproducible(synthetic_calibration_data: CalibrationData) -> None:
    profiles_df, taxa_df = synthetic_calibration_data
    first = run_ora_calibration_benchmark(
        IlyasORAAdapter(), profiles_df, taxa_df, spike_traits=[TRAITS[0]], replicates=3, seed=9
    )
    second = run_ora_calibration_benchmark(
        IlyasORAAdapter(), profiles_df, taxa_df, spike_traits=[TRAITS[0]], replicates=3, seed=9
    )
    pd.testing.assert_frame_equal(first.null_pvalues, second.null_pvalues)


class _FlakyAdapter(IlyasORAAdapter):
    """Adapter whose first call raises, to exercise failed-replicate accounting."""

    def __init__(self) -> None:
        self.calls = 0

    def run_inference(
        self,
        profiles_df: pd.DataFrame,
        taxa_matrix_df: pd.DataFrame,
        config: dict[str, Any] | None = None,
    ) -> MethodInferenceResult:
        self.calls += 1
        if self.calls == 1:
            raise RuntimeError("boom")
        return super().run_inference(profiles_df, taxa_matrix_df, config)


def test_benchmark_counts_failed_replicates(synthetic_calibration_data: CalibrationData) -> None:
    profiles_df, taxa_df = synthetic_calibration_data
    report = run_ora_calibration_benchmark(
        _FlakyAdapter(), profiles_df, taxa_df, spike_traits=[], replicates=3, scenario="null"
    )
    assert report.status == "PARTIAL"
    assert report.failed_replicates == 1
    assert report.completed_replicates == 2
    assert "boom" in (report.replicate_records[0].error_message or "")


def test_benchmark_with_no_testable_pairs_is_failed(
    synthetic_calibration_data: CalibrationData,
) -> None:
    profiles_df, taxa_df = synthetic_calibration_data
    no_traits = taxa_df.assign(**dict.fromkeys(TRAITS, 0))
    report = run_ora_calibration_benchmark(
        IlyasORAAdapter(), profiles_df, no_traits, spike_traits=[], replicates=2, scenario="null"
    )
    assert report.status == "FAILED"
    assert report.invalid_replicates == 2
    assert report.complete_null_type_i is None


# --- CLI ---------------------------------------------------------------------


def test_cli_writes_results_and_provenance(
    synthetic_calibration_data: CalibrationData, tmp_path: Path
) -> None:
    profiles, taxa = _write_inputs(tmp_path, synthetic_calibration_data)
    out = tmp_path / "out"
    exit_code = harness_main(_cli_args(profiles, taxa, out, "--spike-traits", TRAITS[0]))
    assert exit_code == EXIT_SUCCESS

    summary = _read_summary(out)
    assert summary["status"] == "SUCCESS"
    assert summary["requested_replicates"] == 4
    assert summary["completed_replicates"] == 4
    assert summary["failed_replicates"] == 0
    assert summary["invalid_replicates"] == 0
    assert len(summary["replicates"]) == 4
    assert summary["n_null_pvalues"] > 0
    assert summary["complete_null_type_i"]["replicates"] == 4
    assert summary["complete_null_type_i"]["n_tests"] == summary["n_null_pvalues"]
    assert summary["mixed_signal_fdr"]["replicates"] == 4
    prov = summary["provenance"]
    assert prov["command_args"]["seed"] == 7
    assert prov["command_args"]["spike_traits"] == [TRAITS[0]]
    assert len(prov["inputs"]["profiles"]["sha256"]) == 64
    assert prov["inputs"]["taxa_matrix"]["rows"] == N_TAXA

    reps = pd.read_csv(out / REPLICATES_FILENAME, sep="\t")
    assert reps["status"].tolist() == ["completed"] * 4
    pvals = pd.read_csv(out / NULL_PVALUES_FILENAME, sep="\t")
    assert len(pvals) == summary["n_null_pvalues"]
    assert pvals["null_p_value"].between(0, 1).all()
    assert set(pvals["replicate_index"]) == {0, 1, 2, 3}


def _missing_taxa_matrix(profiles: Path, taxa: Path) -> tuple[Path, Path]:
    return profiles, taxa.with_name("missing.tsv")


def _empty_profiles(profiles: Path, taxa: Path) -> tuple[Path, Path]:
    profiles.write_text("study_id\texperiment_id\tinc_dec\ttaxa_count\n", encoding="utf-8")
    return profiles, taxa


def _profiles_without_signature_size(profiles: Path, taxa: Path) -> tuple[Path, Path]:
    df = pd.read_csv(profiles, sep="\t").drop(columns="taxa_count")
    df.to_csv(profiles, sep="\t", index=False)
    return profiles, taxa


@pytest.mark.parametrize(
    "mutate",
    [_missing_taxa_matrix, _empty_profiles, _profiles_without_signature_size],
    ids=["missing-taxa-matrix", "empty-profiles", "no-signature-size"],
)
def test_cli_invalid_inputs_exit_nonzero_and_record_failure(
    synthetic_calibration_data: CalibrationData,
    tmp_path: Path,
    mutate: Callable[[Path, Path], tuple[Path, Path]],
) -> None:
    profiles, taxa = mutate(*_write_inputs(tmp_path, synthetic_calibration_data))
    out = tmp_path / "out"
    assert harness_main(_cli_args(profiles, taxa, out)) == EXIT_FAILURE
    summary = _read_summary(out)
    assert summary["status"] == "FAILED"
    assert summary["provenance"]["completed_replicates"] == 0
    assert summary["provenance"]["error_message"]
    assert not (out / NULL_PVALUES_FILENAME).exists()


@pytest.mark.parametrize("bad_args", [["--replicates", "0"], ["--alpha", "1.5"]])
def test_cli_rejects_bad_arguments(tmp_path: Path, bad_args: list[str]) -> None:
    with pytest.raises(SystemExit) as exc:
        harness_main([*bad_args, "--output-dir", str(tmp_path)])
    assert exc.value.code != 0


def test_cli_failed_run_removes_stale_results(
    synthetic_calibration_data: CalibrationData, tmp_path: Path
) -> None:
    profiles, taxa = _write_inputs(tmp_path, synthetic_calibration_data)
    out = tmp_path / "out"
    assert harness_main(_cli_args(profiles, taxa, out)) == EXIT_SUCCESS
    assert (out / NULL_PVALUES_FILENAME).exists()
    assert harness_main(_cli_args(profiles, tmp_path / "gone.tsv", out)) == EXIT_FAILURE
    assert not (out / NULL_PVALUES_FILENAME).exists()
    assert _read_summary(out)["status"] == "FAILED"


def test_standard_module_command_executes_benchmark(
    synthetic_calibration_data: CalibrationData, tmp_path: Path
) -> None:
    """``python -m bsdbng.calibration.harness`` (what ``just calibrate`` runs) must do work."""
    profiles, taxa = _write_inputs(tmp_path, synthetic_calibration_data)
    out = tmp_path / "out"
    module_cmd = [sys.executable, "-m", "bsdbng.calibration.harness"]

    proc = subprocess.run(
        [*module_cmd, *_cli_args(profiles, taxa, out)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == EXIT_SUCCESS, proc.stderr
    summary = _read_summary(out)
    assert summary["status"] == "SUCCESS"
    assert summary["completed_replicates"] == 4
    assert summary["n_null_pvalues"] > 0
    assert summary["complete_null_type_i"] is not None

    bad = subprocess.run(
        [*module_cmd, *_cli_args(profiles, tmp_path / "missing.tsv", out)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert bad.returncode == EXIT_FAILURE
    assert _read_summary(out)["status"] == "FAILED"


def test_ilyas_ora_adapter_parity_with_production(
    synthetic_calibration_data: CalibrationData, tmp_path: Path
) -> None:
    """Verify IlyasORAAdapter produces identical results to scripts/trait_enrichment.py."""
    profiles_df, taxa_df = synthetic_calibration_data
    adapter = IlyasORAAdapter()
    adapter_res = adapter.run_inference(
        profiles_df, taxa_df, config={"alpha": 0.05, "min_fold": 1.5, "test_type": "enrichment"}
    )
    adapter_sig = adapter_res.results_df[adapter_res.results_df["is_significant"]]

    profiles_path, taxa_path = _write_inputs(tmp_path, synthetic_calibration_data)
    rep_txt = tmp_path / "report.txt"
    out_tsv = tmp_path / "out.tsv"
    subprocess.run(
        [
            sys.executable,
            "scripts/trait_enrichment.py",
            "--input",
            str(profiles_path),
            "--taxa-matrix",
            str(taxa_path),
            "--output",
            str(rep_txt),
            "--tsv-output",
            str(out_tsv),
            "--p-threshold",
            "0.05",
            "--min-fold",
            "1.5",
        ],
        check=True,
    )
    prod_df = pd.read_csv(out_tsv, sep="\t", comment="#")
    assert len(adapter_sig) == len(prod_df)
    for _, row in prod_df.iterrows():
        pid = f"{row['study_id']}:{row['experiment_id']}:{row['inc_dec']}"
        trait = row["trait"]
        match = adapter_sig[(adapter_sig["profile_id"] == pid) & (adapter_sig["trait"] == trait)]
        assert len(match) == 1
        np.testing.assert_allclose(match["p_value"].iloc[0], row["p_value"], rtol=1e-10)
        np.testing.assert_allclose(match["q_value"].iloc[0], row["fdr_q"], rtol=1e-10)
        np.testing.assert_allclose(match["effect_size"].iloc[0], row["fold"], rtol=1e-10)


def test_support_constraint_violations_raise(
    synthetic_calibration_data: CalibrationData,
) -> None:
    profiles_df, taxa_df = synthetic_calibration_data
    adapter = IlyasORAAdapter()

    # M < max(N): background universe smaller than signature size
    tiny_taxa = taxa_df.head(2)
    with pytest.raises(ValueError, match=r"Background universe M .* is smaller"):
        adapter.run_inference(profiles_df, tiny_taxa)

    # N < max(k_obs): signature size smaller than observed count
    bad_draw = profiles_df.copy()
    bad_draw.loc[0, "taxa_count"] = 1
    with pytest.raises(ValueError, match="Observed trait count k exceeds profile draw count N"):
        adapter.run_inference(bad_draw, taxa_df)

    # k_obs > K: observed count exceeds population count
    bad_pop = taxa_df.copy()
    bad_pop[TRAITS[0]] = 0
    with pytest.raises(ValueError, match="Observed trait count k exceeds population trait count K"):
        adapter.run_inference(profiles_df, bad_pop)


def test_complete_null_preserves_coannotation_bundles() -> None:
    """Whole-taxon sampling must preserve identical counts for perfectly co-annotated traits."""
    m_universe = 50
    co_presence = [1 if i % 2 == 0 else 0 for i in range(m_universe)]
    taxa_df = pd.DataFrame(
        {
            "taxon_name": [f"Taxon {j}" for j in range(m_universe)],
            "taxonomic_rank": ["species"] * m_universe,
            "TRAIT_A": co_presence,
            "TRAIT_B": co_presence,
        },
        index=pd.Index([f"NCBITaxon:{j}" for j in range(m_universe)], name="taxon_id"),
    )
    profiles_df = pd.DataFrame(
        {
            "study_id": ["s1", "s1"],
            "experiment_id": ["e1", "e1"],
            "inc_dec": ["increased", "decreased"],
            "taxa_count": [15, 20],
            "TRAIT_A": [5, 7],
            "TRAIT_B": [5, 7],
        }
    )
    rng = np.random.default_rng(42)
    for _ in range(10):
        null_df = generate_complete_null_ora_profiles(profiles_df, taxa_df, rng)
        np.testing.assert_array_equal(null_df["TRAIT_A"].to_numpy(), null_df["TRAIT_B"].to_numpy())


def test_replicate_monte_carlo_se_invariance_under_within_replicate_duplication() -> None:
    """Duplicating tests within replicates must not reduce Monte Carlo SE across replicates."""
    rep_rates = [0.04, 0.06, 0.05, 0.07, 0.03]
    p_vals_rep = []
    for r in rep_rates:
        n_sig = round(r * 100)
        p_vals_rep.append([0.01] * n_sig + [0.50] * (100 - n_sig))

    pooled_pvals = [p for sub in p_vals_rep for p in sub]
    metric_single = compute_type_i_error(pooled_pvals, alpha=0.05, replicate_rates=rep_rates)

    duplicated_rates = rep_rates
    duplicated_pvals = pooled_pvals * 10
    metric_dup = compute_type_i_error(
        duplicated_pvals, alpha=0.05, replicate_rates=duplicated_rates
    )

    assert metric_single.replicates == 5
    assert metric_dup.replicates == 5
    assert metric_single.monte_carlo_se == pytest.approx(metric_dup.monte_carlo_se)
    assert metric_single.empirical_rate == pytest.approx(metric_dup.empirical_rate)


def test_ilyas_ora_adapter_drops_unpaired_experiments(
    synthetic_calibration_data: CalibrationData,
) -> None:
    profiles_df, taxa_df = synthetic_calibration_data
    unpaired_row = pd.DataFrame(
        {
            "study_id": ["s_unpaired"],
            "experiment_id": ["e_unpaired"],
            "inc_dec": ["increased"],
            "taxa_count": [6],
            **{trait: [2] for trait in TRAITS},
        }
    )
    with_unpaired = pd.concat([profiles_df, unpaired_row], ignore_index=True)
    adapter = IlyasORAAdapter()
    res = adapter.run_inference(with_unpaired, taxa_df)
    assert "s_unpaired:e_unpaired:increased" not in res.results_df["profile_id"].values
