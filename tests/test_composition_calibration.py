"""Tests for paired composition arm calibration and JedCompositionAdapter (#263)."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import bsdbng.composition_calibration as bsdbng_cc
from mtp.composition_calibration import (
    STREAM,
    JedCompositionAdapter,
    MethodAdapter,
    MethodInferenceResult,
    compute_composition_contrasts,
    compute_profile_zscores,
    filter_paired_profiles,
    generate_study_flip_null_profiles,
    replicate_rng,
    run_composition_calibration,
    stream_key,
)
from mtp.composition_calibration.composition import (
    load_composition_trait_columns,
)
from mtp.composition_calibration.composition import (
    main as composition_main,
)


@pytest.fixture
def synthetic_paired_profiles() -> pd.DataFrame:
    """Synthetic profile dataset with 4 studies, 6 paired experiments, and 4 traits."""
    return pd.DataFrame(
        {
            "study_id": ["s1", "s1", "s1", "s1", "s2", "s2", "s3", "s3", "s4", "s4"],
            "experiment_id": ["e1", "e1", "e2", "e2", "e3", "e3", "e4", "e4", "e5", "e5"],
            "inc_dec": [
                "increased",
                "decreased",
                "increased",
                "decreased",
                "increased",
                "decreased",
                "increased",
                "decreased",
                "increased",
                "increased",  # unpaired experiment e5 in s4
            ],
            "condition": ["disease"] * 10,
            "body_site": ["feces"] * 10,
            "taxa_count": [10] * 10,
            "trait_A": [5, 1, 6, 2, 8, 3, 4, 1, 5, 5],
            "trait_B": [1, 5, 2, 7, 1, 4, 2, 6, 3, 3],
            "trait_C": [4, 4, 2, 1, 1, 3, 4, 3, 2, 2],
            "trait_const": [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
        }
    )


def test_filter_paired_profiles(synthetic_paired_profiles: pd.DataFrame) -> None:
    """Ensure unpaired experiments are dropped while paired experiments are retained."""
    filtered = filter_paired_profiles(synthetic_paired_profiles)
    assert len(filtered) == 8
    assert "e5" not in filtered["experiment_id"].values
    assert set(filtered["study_id"].unique()) == {"s1", "s2", "s3"}

    # Test empty dataframe
    empty = filter_paired_profiles(pd.DataFrame())
    assert empty.empty


def test_compute_profile_zscores(synthetic_paired_profiles: pd.DataFrame) -> None:
    """Verify log-ratio fractions and z-score standardizations."""
    traits = ["trait_A", "trait_B", "trait_C", "trait_const"]
    zscores = compute_profile_zscores(synthetic_paired_profiles, traits)

    assert zscores.shape == (10, 4)
    # Constant trait has zero std; should be filled with 0.0 without NaN
    assert np.allclose(zscores["trait_const"], 0.0)
    # Non-constant traits should have mean close to 0 and std close to 1
    for t in ["trait_A", "trait_B", "trait_C"]:
        assert abs(zscores[t].mean()) < 1e-10
        assert abs(zscores[t].std(ddof=0) - 1.0) < 1e-10


@pytest.mark.parametrize("selected", [["trait_A"], ["trait_C", "trait_A"]])
def test_selected_traits_preserve_full_denominator(
    synthetic_paired_profiles: pd.DataFrame, selected: list[str]
) -> None:
    """Testing one trait must not turn every positive fraction into one (#534)."""
    all_traits = ["trait_A", "trait_B", "trait_C", "trait_const"]
    full = compute_profile_zscores(synthetic_paired_profiles, all_traits)
    subset = compute_profile_zscores(synthetic_paired_profiles, selected)
    pd.testing.assert_frame_equal(subset, full[selected])
    assert not subset["trait_A"].eq(0).all()

    full_contrasts = compute_composition_contrasts(synthetic_paired_profiles, all_traits)
    subset_contrasts = compute_composition_contrasts(synthetic_paired_profiles, selected)
    pd.testing.assert_frame_equal(
        subset_contrasts.pair_contrasts, full_contrasts.pair_contrasts[selected]
    )
    pd.testing.assert_frame_equal(
        subset_contrasts.study_contrasts, full_contrasts.study_contrasts[selected]
    )
    assert subset_contrasts.n_pairs == full_contrasts.n_pairs == 4
    assert subset_contrasts.n_studies == full_contrasts.n_studies == 3


def test_denominator_keeps_untested_counts_and_excludes_numeric_metadata() -> None:
    """Untested biological counts matter; signature size and media do not (#534)."""
    profiles = pd.DataFrame(
        {
            "trait_A": [10.0, 10.0, 10.0, 0.0],
            "untested_trait": [10.0, 90.0, 990.0, 0.0],
            "study_id": [1, 2, 3, 4],
            "experiment_id": [11, 12, 13, 14],
            "inc_dec": ["inc", "dec", "inc", "dec"],
            "n_taxa": [900, 1, 80, 500],
            "taxa_count": [90, 10, 800, 50],
            "signature_size": [9, 100, 8, 5],
            "profile_id": [100, 200, 300, 400],
            "host_species": [9606, 9606, 9606, 9606],
            "_norm_dir": ["inc", "dec", "inc", "dec"],
            "mediadive.medium:1": [10000.0, 0.0, 50.0, 800.0],
        }
    )
    observed = compute_profile_zscores(profiles, ["trait_A"])
    transformed = np.log1p(1000 * np.array([0.5, 0.1, 0.01, 0.0]))
    expected = (transformed - transformed.mean()) / transformed.std(ddof=0)
    assert observed["trait_A"].to_numpy() == pytest.approx(expected)

    biological_only = profiles[["trait_A", "untested_trait"]]
    pd.testing.assert_frame_equal(observed, compute_profile_zscores(biological_only, ["trait_A"]))


def test_adapter_subset_changes_test_family_without_changing_raw_inference(
    synthetic_paired_profiles: pd.DataFrame,
) -> None:
    """Multiplicity uses selected traits while their effects and p-values stay fixed."""
    adapter = JedCompositionAdapter()
    full = adapter.run_inference(synthetic_paired_profiles, pd.DataFrame())
    subset = adapter.run_inference(
        synthetic_paired_profiles, pd.DataFrame(), {"traits": ["trait_A"]}
    )
    assert full.metadata["n_tests"] == 4
    assert subset.metadata["n_tests"] == 1
    raw_columns = [
        "effect_size",
        "t_statistic",
        "p_value",
        "mean_contrast",
        "median_contrast",
        "std_contrast",
        "ci95_lower",
        "ci95_upper",
        "n_studies",
        "n_pairs",
    ]
    pd.testing.assert_frame_equal(
        subset.results_df.set_index("trait")[raw_columns],
        full.results_df.set_index("trait").loc[["trait_A"], raw_columns],
    )
    assert subset.results_df.iloc[0]["effect_size"] > 0


@pytest.mark.parametrize("corpus", ["paired", "empty", "unpaired"])
@pytest.mark.parametrize(
    ("selected", "expected"),
    [
        (["n_taxa", "profile_id", "_norm_dir"], []),
        (["MediaDive.Medium", "missing_trait"], []),
        (["profile_id", "trait_A", "MediaDive.Medium", "missing_trait"], ["trait_A"]),
    ],
)
def test_nontrait_requests_do_not_become_empty_or_phantom_contrasts(
    synthetic_paired_profiles: pd.DataFrame,
    corpus: str,
    selected: list[str],
    expected: list[str],
) -> None:
    """Empty cohorts and rejected trait requests retain one output schema (#534)."""
    profiles = synthetic_paired_profiles.assign(
        n_taxa=200, profile_id=123, **{"MediaDive.Medium": 999}
    )
    profiles["_norm_dir"] = profiles["inc_dec"].map({"increased": "inc", "decreased": "dec"})
    if corpus == "empty":
        profiles = profiles.iloc[:0]
    elif corpus == "unpaired":
        profiles = profiles.loc[profiles["_norm_dir"].eq("inc")]

    zscores = compute_profile_zscores(profiles, selected)
    assert zscores.columns.tolist() == expected
    assert zscores.index.equals(profiles.index)
    contrasts = compute_composition_contrasts(profiles, selected)
    assert contrasts.trait_columns == expected
    assert contrasts.pair_contrasts.columns.tolist() == expected
    assert contrasts.study_contrasts.columns.tolist() == expected

    result = JedCompositionAdapter().run_inference(profiles, pd.DataFrame(), {"traits": selected})
    assert result.metadata["n_trait_columns"] == len(expected)
    if not expected or corpus != "paired":
        assert contrasts.n_pairs == 0
        assert contrasts.n_studies == 0
        assert result.results_df.empty
        assert result.metadata["n_tests"] == 0


def test_compute_composition_contrasts(synthetic_paired_profiles: pd.DataFrame) -> None:
    """Verify matched experiment differences and study aggregation."""
    traits = ["trait_A", "trait_B", "trait_C"]
    data = compute_composition_contrasts(synthetic_paired_profiles, traits)

    assert data.n_pairs == 4  # e1, e2 in s1; e3 in s2; e4 in s3
    assert data.n_studies == 3  # s1, s2, s3
    assert data.study_contrasts.shape == (3, 3)

    # trait_A is always higher in increased than decreased, so contrasts must be positive
    assert (data.study_contrasts["trait_A"] > 0).all()
    # trait_B is always higher in decreased than increased, so contrasts must be negative
    assert (data.study_contrasts["trait_B"] < 0).all()


def test_direction_flip_invariance(synthetic_paired_profiles: pd.DataFrame) -> None:
    """Flipping direction labels in a study must exactly negate its contrasts."""
    traits = ["trait_A", "trait_B", "trait_C"]
    orig_data = compute_composition_contrasts(synthetic_paired_profiles, traits)

    # Flip study s2
    flipped_df = synthetic_paired_profiles.copy()
    swap = {"increased": "decreased", "decreased": "increased"}
    mask = flipped_df["study_id"] == "s2"
    flipped_df.loc[mask, "inc_dec"] = flipped_df.loc[mask, "inc_dec"].map(swap)

    flipped_data = compute_composition_contrasts(flipped_df, traits)

    # s1 and s3 unchanged
    assert np.allclose(orig_data.study_contrasts.loc["s1"], flipped_data.study_contrasts.loc["s1"])
    assert np.allclose(orig_data.study_contrasts.loc["s3"], flipped_data.study_contrasts.loc["s3"])

    # s2 exactly negated
    assert np.allclose(orig_data.study_contrasts.loc["s2"], -flipped_data.study_contrasts.loc["s2"])


def test_generate_study_flip_null_profiles(synthetic_paired_profiles: pd.DataFrame) -> None:
    """Verify generate_study_flip_null_profiles swaps direction labels properly."""
    rng = replicate_rng(182, "test_stream", 0)
    flipped = generate_study_flip_null_profiles(synthetic_paired_profiles, rng, unit="study")

    assert len(flipped) == len(synthetic_paired_profiles)
    assert set(flipped.columns) == set(synthetic_paired_profiles.columns)

    # Unflipped studies preserve exact inc_dec
    for s_id in synthetic_paired_profiles["study_id"].unique():
        orig_s = synthetic_paired_profiles[synthetic_paired_profiles["study_id"] == s_id]
        flip_s = flipped[flipped["study_id"] == s_id]
        is_same = (orig_s["inc_dec"].values == flip_s["inc_dec"].values).all()
        is_swapped = (
            orig_s["inc_dec"].map({"increased": "decreased", "decreased": "increased"}).values
            == flip_s["inc_dec"].values
        ).all()
        assert is_same or is_swapped


def test_replicate_rng_determinism() -> None:
    """Verify stream_key and replicate_rng reproducibility."""
    key = stream_key("complete_null_study_flip")
    assert isinstance(key, int)

    rng1 = replicate_rng(182, "complete_null_study_flip", 0)
    rng2 = replicate_rng(182, "complete_null_study_flip", 0)
    draws1 = [rng1.random() for _ in range(5)]
    draws2 = [rng2.random() for _ in range(5)]

    assert np.allclose(draws1, draws2)
    assert abs(draws1[0] - 0.42837155989674747) < 1e-12


def test_adapter_protocol_compliance(synthetic_paired_profiles: pd.DataFrame) -> None:
    """Verify JedCompositionAdapter conforms to MethodAdapter protocol."""
    adapter = JedCompositionAdapter()
    assert isinstance(adapter, MethodAdapter)
    assert adapter.method_id == "jed_composition_ttest"
    assert adapter.owner == "Jed"

    res = adapter.run_inference(
        synthetic_paired_profiles,
        pd.DataFrame(),
        config={"traits": ["trait_A", "trait_B", "trait_C"], "alpha": 0.05},
    )
    assert isinstance(res, MethodInferenceResult)
    assert len(res.results_df) == 3
    assert set(res.results_df.columns) >= {
        "effect_size",
        "is_significant",
        "n_pairs",
        "n_studies",
        "p_value",
        "q_value",
        "t_statistic",
        "trait",
    }
    assert res.metadata["n_studies"] == 3
    assert res.metadata["n_pairs"] == 4


def test_run_composition_calibration_smoke(synthetic_paired_profiles: pd.DataFrame) -> None:
    """Run smoke calibration with 3 replicates on synthetic data."""
    traits = ["trait_A", "trait_B", "trait_C"]
    summary_dict, summary_df, replicates_df, null_pvals_df, _table1_df = (
        run_composition_calibration(
            profiles_df=synthetic_paired_profiles,
            trait_columns=traits,
            replicates=3,
            seed=182,
            stream=STREAM,
        )
    )

    assert summary_dict["replicates_completed"] == 3
    assert summary_dict["n_studies"] == 3
    assert summary_dict["n_traits"] == 3
    assert len(summary_df) == 6  # 2 corrections (BH, BY) x 3 thresholds (0.05, 0.01, 0.001)
    assert len(replicates_df) == 3 * 6
    assert len(null_pvals_df) == 3 * 3


def test_composition_cli_smoke(tmp_path: Path, synthetic_paired_profiles: pd.DataFrame) -> None:
    """Verify composition calibration CLI writes all artifacts successfully."""
    profiles_path = tmp_path / "profiles.tsv"
    columns_path = tmp_path / "columns.tsv"
    out_dir = tmp_path / "out"

    synthetic_paired_profiles.to_csv(profiles_path, sep="\t", index=False)
    cols_df = pd.DataFrame(
        {
            "column_id": ["trait_A", "trait_B", "trait_C"],
            "object_category": ["category_1", "category_1", "category_1"],
        }
    )
    cols_df.to_csv(columns_path, sep="\t", index=False)

    exit_code = composition_main(
        [
            "--profiles",
            str(profiles_path),
            "--columns",
            str(columns_path),
            "--output-dir",
            str(out_dir),
            "--replicates",
            "2",
            "--seed",
            "182",
        ]
    )

    assert exit_code == 0
    assert (out_dir / "composition_calibration_summary.json").is_file()
    assert (out_dir / "composition_calibration_summary.tsv").is_file()
    assert (out_dir / "composition_calibration_replicates.tsv").is_file()
    assert (out_dir / "composition_null_pvalues.tsv.gz").is_file()
    assert (out_dir / "table1_composition_calibration.tsv").is_file()

    summary_json = json.loads((out_dir / "composition_calibration_summary.json").read_text())
    assert summary_json["replicates_completed"] == 2
    assert summary_json["n_studies"] == 3


def test_bsdbng_shim_exports() -> None:
    """Verify bsdbng compatibility shim exposes JedCompositionAdapter and functions."""
    assert hasattr(bsdbng_cc, "JedCompositionAdapter")
    assert bsdbng_cc.JedCompositionAdapter is JedCompositionAdapter
    assert hasattr(bsdbng_cc, "run_composition_calibration")
    assert bsdbng_cc.run_composition_calibration is run_composition_calibration


def test_missing_or_malformed_columns_file_is_an_error(tmp_path: Path) -> None:
    """The columns file sets m for BH/BY; never fall back silently (#469)."""
    with pytest.raises(FileNotFoundError):
        load_composition_trait_columns(tmp_path / "absent.tsv")
    bad = tmp_path / "bad.tsv"
    bad.write_text("column_id\nfoo\n")
    with pytest.raises(ValueError, match="object_category"):
        load_composition_trait_columns(bad)
