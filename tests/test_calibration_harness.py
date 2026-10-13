"""Tests for the statistical calibration harness, adapters, and scenarios (#348, #351, #343)."""

from __future__ import annotations

import json
import subprocess
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any, cast

import numpy as np
import pandas as pd
import pytest

from bsdbng.calibration import (
    SEED,
    SHARED_FLIP_REPLICATES,
    flip_directions,
    replicate_rng,
)
from bsdbng.ora_calibration import (
    IlyasORAAdapter,
    compute_type_i_error,
    generate_complete_null_ora_profiles,
    generate_mixed_signal_ora_profiles,
    generate_study_flip_null_ora_profiles,
)
from bsdbng.ora_calibration.harness import (
    EXIT_FAILURE,
    EXIT_SUCCESS,
    NULL_PVALUES_FILENAME,
    SUMMARY_FILENAME,
    run_ora_calibration_benchmark,
)
from bsdbng.ora_calibration.harness import main as harness_main
from bsdbng.ora_calibration.scenarios import flip_profile_dataframe
from bsdbng.paths import DERIVED_DATA_DIR
from bsdbng.trait_guilds import Profile

TRAITS: list[str] = ["TRAIT_A", "TRAIT_B", "TRAIT_C", "TRAIT_D", "TRAIT_E", "TRAIT_F"]
CalibrationData = tuple[pd.DataFrame, pd.DataFrame]


@pytest.fixture
def synthetic_calibration_data() -> CalibrationData:
    """Generate minimal paired profile and taxa-matrix tables."""
    profiles_df = pd.DataFrame(
        {
            "study_id": ["s1", "s1", "s2", "s2"],
            "experiment_id": ["e1", "e1", "e2", "e2"],
            "inc_dec": ["inc", "dec", "inc", "dec"],
            "taxa_count": [10, 12, 8, 9],
            "TRAIT_A": [4, 1, 3, 2],
            "TRAIT_B": [2, 5, 1, 4],
            "TRAIT_C": [1, 2, 0, 1],
            "TRAIT_D": [3, 0, 2, 1],
            "TRAIT_E": [0, 1, 1, 0],
            "TRAIT_F": [0, 0, 1, 2],
        }
    )
    rng = np.random.default_rng(1234)
    m_taxa = 20
    taxa_ids = [f"NCBITaxon:{1000 + i}" for i in range(m_taxa)]
    taxa_traits = rng.binomial(1, 0.4, size=(m_taxa, len(TRAITS)))
    taxa_traits[:4, :] = 1

    taxa_df = pd.DataFrame(
        taxa_traits,
        index=pd.Index(taxa_ids, name="taxon_id"),
        columns=TRAITS,
    )
    taxa_df.insert(0, "taxon_name", [f"Taxon {i}" for i in range(m_taxa)])
    taxa_df.insert(1, "taxonomic_rank", ["species"] * m_taxa)
    return profiles_df, taxa_df


def _write_inputs(tmp_path: Path, data: CalibrationData) -> tuple[Path, Path]:
    profiles_path = tmp_path / "profiles.tsv"
    taxa_path = tmp_path / "taxa.tsv"
    data[0].to_csv(profiles_path, sep="\t", index=False)
    data[1].to_csv(taxa_path, sep="\t")
    return profiles_path, taxa_path


def _read_summary(out_dir: Path) -> dict[str, Any]:
    return cast(
        dict[str, Any],
        json.loads((out_dir / SUMMARY_FILENAME).read_text(encoding="utf-8")),
    )


def _cli_args(profiles: Path, taxa: Path, out: Path, replicates: int = 4) -> list[str]:
    return [
        "--profiles",
        str(profiles),
        "--taxa-matrix",
        str(taxa),
        "--output-dir",
        str(out),
        "--replicates",
        str(replicates),
        "--seed",
        "777",
        "--scenario",
        "all",
        "--signal-strength",
        "0.6",
    ]


def test_reexport_of_study_flip_symbols_for_screen_calibration() -> None:
    """bsdbng.calibration must re-export study flip symbols for PR #319 (#343)."""
    assert SEED == 182
    assert SHARED_FLIP_REPLICATES == 100
    rng = replicate_rng(SEED, "complete_null_study_flip", 0)
    assert isinstance(rng, np.random.Generator)

    dummy_profiles = [
        Profile("s1", "e1", "inc", "cond1", "gut", {}, 0),
        Profile("s1", "e1", "dec", "cond1", "gut", {}, 0),
    ]
    flipped = flip_directions(dummy_profiles, rng, unit="study")
    assert len(flipped) == 2


def test_calibration_scenarios_binarize_taxa_matrix_counts() -> None:
    """Scenario generators must binarize count-matrix cells to presence (#348)."""
    m_taxa = 200
    n_draw = 10
    rng = np.random.default_rng(42)

    # Taxa matrix with cells 0, 1, 2, 3 (sub-taxa counts)
    raw_counts = rng.integers(0, 4, size=(m_taxa, len(TRAITS)))
    taxa_df = pd.DataFrame(
        raw_counts,
        index=pd.Index([f"NCBITaxon:{i}" for i in range(m_taxa)], name="taxon_id"),
        columns=TRAITS,
    )
    taxa_df.insert(0, "taxon_name", [f"Taxon {i}" for i in range(m_taxa)])
    taxa_df.insert(1, "taxonomic_rank", ["species"] * m_taxa)

    profiles_df = pd.DataFrame(
        {
            "study_id": ["s1", "s1"],
            "experiment_id": ["e1", "e1"],
            "inc_dec": ["increased", "decreased"],
            "taxa_count": [n_draw, n_draw],
            "TRAIT_A": [3, 2],
            "TRAIT_B": [2, 1],
            "TRAIT_C": [1, 2],
            "TRAIT_D": [1, 0],
            "TRAIT_E": [0, 1],
            "TRAIT_F": [0, 0],
        }
    )

    # 1. Complete null generation must not exceed draw count N
    for _ in range(5):
        null_df = generate_complete_null_ora_profiles(profiles_df, taxa_df, rng)
        for trait in TRAITS:
            assert (null_df[trait] <= n_draw).all(), (
                f"Simulated {trait} count exceeded signature size {n_draw}"
            )

        # IlyasORAAdapter should run without raising 'Observed trait count k exceeds N'
        adapter = IlyasORAAdapter()
        res = adapter.run_inference(null_df, taxa_df)
        assert len(res.results_df) > 0

    # 2. Mixed signal generation must also binarize
    mixed_df, _ = generate_mixed_signal_ora_profiles(
        profiles_df, taxa_df, spike_traits=["TRAIT_A"], signal_strength=0.6, rng=rng
    )
    for trait in TRAITS:
        assert (mixed_df[trait] <= n_draw).all()


def test_study_flip_null_generator_flips_units_and_preserves_counts(
    synthetic_calibration_data: CalibrationData,
) -> None:
    """Study flip null must exchange inc/dec labels while preserving counts (#263, #343)."""
    profiles_df, _ = synthetic_calibration_data
    rng = np.random.default_rng(42)
    flipped_df = generate_study_flip_null_ora_profiles(profiles_df, rng)

    assert len(flipped_df) == len(profiles_df)
    assert (flipped_df["taxa_count"].to_numpy() == profiles_df["taxa_count"].to_numpy()).all()
    for trait in TRAITS:
        assert (flipped_df[trait].to_numpy() == profiles_df[trait].to_numpy()).all()
    assert (flipped_df["inc_dec"].to_numpy() != profiles_df["inc_dec"].to_numpy()).any(), (
        "Expected at least one study direction label to flip under random null (#396)"
    )


def test_ilyas_ora_adapter_parity_with_production(
    synthetic_calibration_data: CalibrationData, tmp_path: Path
) -> None:
    """Verify IlyasORAAdapter produces identical enrichment results
    to scripts/trait_enrichment.py across both significant hits and complete test table (#393).
    """
    profiles_df, taxa_df = synthetic_calibration_data
    adapter = IlyasORAAdapter()
    adapter_res = adapter.run_inference(
        profiles_df, taxa_df, config={"alpha": 0.05, "min_fold": 1.5, "test_type": "enrichment"}
    )
    adapter_sig = adapter_res.results_df[adapter_res.results_df["is_significant"]]

    profiles_path, taxa_path = _write_inputs(tmp_path, synthetic_calibration_data)
    rep_txt = tmp_path / "report.txt"
    out_tsv = tmp_path / "out.tsv"
    all_tsv = tmp_path / "all_out.tsv"
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
            "--all-results-output",
            str(all_tsv),
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

    prod_all_df = pd.read_csv(all_tsv, sep="\t", comment="#")
    assert len(adapter_res.results_df) == len(prod_all_df)
    for _, row in prod_all_df.iterrows():
        pid = f"{row['study_id']}:{row['experiment_id']}:{row['inc_dec']}"
        trait = row["trait"]
        match = adapter_res.results_df[
            (adapter_res.results_df["profile_id"] == pid)
            & (adapter_res.results_df["trait"] == trait)
        ]
        assert len(match) == 1
        np.testing.assert_allclose(match["p_value"].iloc[0], row["p_value"], rtol=1e-10)
        np.testing.assert_allclose(match["q_value"].iloc[0], row["fdr_q"], rtol=1e-10)
        np.testing.assert_allclose(match["effect_size"].iloc[0], row["fold"], rtol=1e-10)


def test_ilyas_ora_adapter_depletion_parity_with_production(
    synthetic_calibration_data: CalibrationData, tmp_path: Path
) -> None:
    """Verify IlyasORAAdapter produces identical depletion results
    to scripts/trait_enrichment.py across hits and full table (#351, #393).
    """
    profiles_df, taxa_df = synthetic_calibration_data
    adapter = IlyasORAAdapter()
    adapter_res = adapter.run_inference(
        profiles_df,
        taxa_df,
        config={"alpha": 0.5, "min_fold": 1.0, "max_fold": 0.8, "test_type": "depletion"},
    )
    adapter_sig = adapter_res.results_df[adapter_res.results_df["is_significant"]]

    profiles_path, taxa_path = _write_inputs(tmp_path, synthetic_calibration_data)
    rep_txt = tmp_path / "report_dep.txt"
    out_tsv = tmp_path / "out_dep.tsv"
    all_tsv = tmp_path / "all_out_dep.tsv"
    subprocess.run(
        [
            sys.executable,
            "scripts/trait_enrichment.py",
            "--tail",
            "lower",
            "--input",
            str(profiles_path),
            "--taxa-matrix",
            str(taxa_path),
            "--output",
            str(rep_txt),
            "--tsv-output",
            str(out_tsv),
            "--all-results-output",
            str(all_tsv),
            "--p-threshold",
            "0.5",
            "--max-fold",
            "0.8",
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

    prod_all_df = pd.read_csv(all_tsv, sep="\t", comment="#")
    assert len(adapter_res.results_df) == len(prod_all_df)
    for _, row in prod_all_df.iterrows():
        pid = f"{row['study_id']}:{row['experiment_id']}:{row['inc_dec']}"
        trait = row["trait"]
        match = adapter_res.results_df[
            (adapter_res.results_df["profile_id"] == pid)
            & (adapter_res.results_df["trait"] == trait)
        ]
        assert len(match) == 1
        np.testing.assert_allclose(match["p_value"].iloc[0], row["p_value"], rtol=1e-10)
        np.testing.assert_allclose(match["q_value"].iloc[0], row["fdr_q"], rtol=1e-10)
        np.testing.assert_allclose(match["effect_size"].iloc[0], row["fold"], rtol=1e-10)


def _drop_sig_col(p: Path, t: Path) -> tuple[Path, Path]:
    df = pd.read_csv(p, sep="\t")
    no_sig = p.parent / "no_sig.tsv"
    df.drop(columns=["taxa_count"]).to_csv(no_sig, sep="\t", index=False)
    return no_sig, t


@pytest.mark.parametrize(
    "mutate",
    [
        lambda p, t: (p, t.parent / "missing_taxa.tsv"),
        lambda p, t: (p.parent / "empty_profiles.tsv", t),
        _drop_sig_col,
    ],
    ids=["missing-taxa-matrix", "empty-profiles", "no-signature-size"],
)
def test_cli_invalid_inputs_exit_nonzero_and_record_failure(
    synthetic_calibration_data: CalibrationData,
    tmp_path: Path,
    mutate: Callable[[Path, Path], tuple[Path, Path]],
) -> None:
    base_profiles, base_taxa = _write_inputs(tmp_path, synthetic_calibration_data)
    if "empty_profiles" in str(mutate):
        (tmp_path / "empty_profiles.tsv").write_text("", encoding="utf-8")
    profiles, taxa = mutate(base_profiles, base_taxa)
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
    """``python -m bsdbng.ora_calibration.harness`` (what ``just calibrate`` runs) must do work."""
    profiles, taxa = _write_inputs(tmp_path, synthetic_calibration_data)
    out = tmp_path / "out"
    module_cmd = [sys.executable, "-m", "bsdbng.ora_calibration.harness"]

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


def test_flip_profile_dataframe_experiment_unit_unsorted_input() -> None:
    """Experiment-level flip mask must align correctly with rows on unsorted input (#395)."""
    rows = [
        {"study_id": "s1", "experiment_id": "e1", "inc_dec": "inc", "taxa_count": 10},
        {"study_id": "s2", "experiment_id": "e2", "inc_dec": "dec", "taxa_count": 12},
        {"study_id": "s1", "experiment_id": "e1", "inc_dec": "dec", "taxa_count": 8},
        {"study_id": "s2", "experiment_id": "e2", "inc_dec": "inc", "taxa_count": 9},
        {"study_id": "s3", "experiment_id": "e3", "inc_dec": "inc", "taxa_count": 11},
        {"study_id": "s3", "experiment_id": "e3", "inc_dec": "dec", "taxa_count": 15},
    ]
    df = pd.DataFrame(rows)
    rng = np.random.default_rng(182)

    for _ in range(10):
        flipped = flip_profile_dataframe(df, rng, unit="experiment")
        for (sid, eid), _ in flipped.groupby(["study_id", "experiment_id"]):
            orig_group = df[(df["study_id"] == sid) & (df["experiment_id"] == eid)]
            orig_inc_idx = orig_group[orig_group["inc_dec"] == "inc"].index[0]
            orig_dec_idx = orig_group[orig_group["inc_dec"] == "dec"].index[0]

            new_inc_label = flipped.loc[orig_inc_idx, "inc_dec"]
            new_dec_label = flipped.loc[orig_dec_idx, "inc_dec"]

            if new_inc_label == "inc":
                assert new_dec_label == "dec"
            else:
                assert new_inc_label == "dec"
                assert new_dec_label == "inc"


def test_mixed_signal_twin_traits_and_unspiked_null_profiles() -> None:
    """Twin traits must be treated as true discoveries and unspiked stay null (#390)."""
    m_universe = 30
    taxa_ids = [f"NCBITaxon:{2000 + i}" for i in range(m_universe)]
    rng = np.random.default_rng(999)
    presence_a = (rng.random(m_universe) < 0.4).astype(np.int64)
    presence_b = presence_a.copy()  # perfect twin
    presence_c = (rng.random(m_universe) < 0.4).astype(np.int64)

    taxa_df = pd.DataFrame(
        {
            "taxon_name": [f"Taxon {i}" for i in range(m_universe)],
            "taxonomic_rank": ["species"] * m_universe,
            "TRAIT_A": presence_a,
            "TRAIT_TWIN": presence_b,
            "TRAIT_OTHER": presence_c,
        },
        index=pd.Index(taxa_ids, name="taxon_id"),
    )
    profiles_df = pd.DataFrame(
        {
            "study_id": ["s1", "s1", "s2", "s2"],
            "experiment_id": ["e1", "e1", "e2", "e2"],
            "inc_dec": ["inc", "dec", "inc", "dec"],
            "taxa_count": [10, 10, 10, 10],
            "TRAIT_A": [3, 3, 3, 3],
            "TRAIT_TWIN": [3, 3, 3, 3],
            "TRAIT_OTHER": [2, 2, 2, 2],
        }
    )
    mixed_df, true_spikes = generate_mixed_signal_ora_profiles(
        profiles_df,
        taxa_df,
        spike_traits=["TRAIT_A"],
        signal_strength=0.7,
        spike_fraction=0.5,
        rng=rng,
    )
    spiked_pids = {pid for pid, _ in true_spikes}
    assert len(spiked_pids) == 2

    for pid in spiked_pids:
        assert (pid, "TRAIT_A") in true_spikes
        assert (pid, "TRAIT_TWIN") in true_spikes
        assert (pid, "TRAIT_OTHER") not in true_spikes

    all_pids = set(
        mixed_df["study_id"] + ":" + mixed_df["experiment_id"] + ":" + mixed_df["inc_dec"]
    )
    unspiked_pids = all_pids - spiked_pids
    for pid in unspiked_pids:
        assert not any(p == pid for p, _ in true_spikes)


def test_mixed_signal_depletion_samples_from_negative_taxa() -> None:
    """Mixed signal depletion must sample from negative indices, forcing k_obs to zero (#390)."""
    m_universe = 30
    taxa_ids = [f"NCBITaxon:{3000 + i}" for i in range(m_universe)]
    # Half positive, half negative
    presence_a = np.array([1] * 15 + [0] * 15, dtype=np.int64)

    taxa_df = pd.DataFrame(
        {
            "taxon_name": [f"Taxon {i}" for i in range(m_universe)],
            "taxonomic_rank": ["species"] * m_universe,
            "TRAIT_A": presence_a,
        },
        index=pd.Index(taxa_ids, name="taxon_id"),
    )
    profiles_df = pd.DataFrame(
        {
            "study_id": ["s1", "s1"],
            "experiment_id": ["e1", "e1"],
            "inc_dec": ["inc", "dec"],
            "taxa_count": [8, 8],
            "TRAIT_A": [4, 4],
        }
    )
    rng = np.random.default_rng(42)
    mixed_df, true_spikes = generate_mixed_signal_ora_profiles(
        profiles_df,
        taxa_df,
        spike_traits=["TRAIT_A"],
        signal_strength=0.6,
        spike_fraction=1.0,
        rng=rng,
        test_type="depletion",
    )
    assert len(true_spikes) == 2
    for count in mixed_df["TRAIT_A"]:
        assert count == 0


def test_independent_scenario_streams_reproducibility(
    synthetic_calibration_data: CalibrationData,
) -> None:
    """Uniform null streams must be identical whether alone or with mixed scenario (#391)."""
    profiles_df, taxa_df = synthetic_calibration_data
    adapter = IlyasORAAdapter()

    rep_null_only = run_ora_calibration_benchmark(
        adapter,
        profiles_df,
        taxa_df,
        spike_traits=["TRAIT_A"],
        replicates=3,
        seed=123,
        scenario="null",
    )
    rep_all = run_ora_calibration_benchmark(
        adapter,
        profiles_df,
        taxa_df,
        spike_traits=["TRAIT_A"],
        replicates=3,
        seed=123,
        scenario="all",
    )

    assert rep_null_only.complete_null_type_i is not None
    assert rep_all.complete_null_type_i is not None
    np.testing.assert_array_equal(
        rep_null_only.null_pvalues["p_value"].to_numpy(),
        rep_all.null_pvalues["p_value"].to_numpy(),
    )
    assert (
        rep_null_only.complete_null_type_i.empirical_rate
        == rep_all.complete_null_type_i.empirical_rate
    )
    assert (
        rep_null_only.complete_null_type_i.monte_carlo_se
        == rep_all.complete_null_type_i.monte_carlo_se
    )


def test_type_i_conservativeness_requires_at_least_two_replicates() -> None:
    """Type I conservativeness must require R >= 2 to avoid trivial SE=0 passing (#396)."""
    single_rep = compute_type_i_error(
        [0.01] * 4 + [0.5] * 96,
        alpha=0.05,
        replicate_rates=[0.04],
    )
    assert single_rep.replicates == 1
    assert not single_rep.is_conservative

    multi_rep = compute_type_i_error(
        [0.01] * 20 + [0.5] * 480,
        alpha=0.05,
        replicate_rates=[0.04, 0.05, 0.03, 0.04, 0.04],
    )
    assert multi_rep.replicates == 5
    assert multi_rep.is_conservative


def test_canary_real_matrix_and_ora_profiles(tmp_path: Path) -> None:
    """Canary test running the benchmark harness on real derived files (#392)."""
    real_profiles = DERIVED_DATA_DIR / "study_trait_profiles_ora_wide.tsv"
    real_taxa = DERIVED_DATA_DIR / "taxa_trait_matrix.tsv"
    if not (real_profiles.is_file() and real_taxa.is_file()):
        pytest.skip("Real derived profiles and taxa matrix not present")

    out = tmp_path / "canary_out"
    code = harness_main(
        [
            "--profiles",
            str(real_profiles),
            "--taxa-matrix",
            str(real_taxa),
            "--output-dir",
            str(out),
            "--replicates",
            "1",
            "--max-profiles",
            "6",
            "--scenario",
            "null",
        ]
    )
    assert code == EXIT_SUCCESS
    summary = _read_summary(out)
    assert summary["status"] == "SUCCESS"
    assert "harness_check_passed" in summary
    assert "is_calibrated" in summary
