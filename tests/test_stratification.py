"""Unit and integration tests for ORA background stratification (#118).

Tests:
1. Anatomical site classification (feces vs. feces_excluded).
2. Profile partitioning across strata.
3. Background annotation density decile binning.
4. Support status assignment and effect heterogeneity computation.
5. End-to-end execution producing stable filenames, required columns, and
   mandatory precautionary / diagnostic notes.
6. Single site stratum analysis (supporting analyze-site recipe).
7. CLI wrapper argument parsing and execution.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from mtp.stratification import (
    ANNOTATION_DENSITY_NOTE,
    PRECAUTIONARY_NOTE,
    PRIMARY_STRATA,
    StratumOutputs,
    build_sensitivity_summary_table,
    classify_site_stratum,
    compute_annotation_density_deciles,
    generate_sensitivity_reports,
    is_feces_site,
    partition_profiles,
    run_stratification_pipeline,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
CLI_WRAPPER = REPO_ROOT / "scripts" / "run_stratified_enrichment.py"


def test_is_feces_site_and_classification() -> None:
    """Verify body site classification accurately recognizes feces/stool terms."""
    assert is_feces_site("feces")
    assert is_feces_site("Feces")
    assert is_feces_site("stool")
    assert is_feces_site("Stool")
    assert is_feces_site("feces, lumen")
    assert is_feces_site("rectum, stool")

    assert not is_feces_site("saliva")
    assert not is_feces_site("vagina")
    assert not is_feces_site("oral cavity")
    assert not is_feces_site("skin")
    assert not is_feces_site(None)
    assert not is_feces_site("")

    assert classify_site_stratum("feces") == "feces"
    assert classify_site_stratum("Saliva") == "feces_excluded"
    assert classify_site_stratum(None) == "feces_excluded"


def test_partition_profiles(tmp_path: Path) -> None:
    """Verify wide profiles table is partitioned correctly by anatomical site."""
    profiles_tsv = tmp_path / "profiles.tsv"
    profiles_tsv.write_text(
        "study_id\texperiment_id\tinc_dec\ttaxa_count\tTraitA\tTraitB\n"
        "s1\te1\tinc\t2\t1\t1\n"
        "s2\te2\tinc\t3\t2\t0\n"
        "s3\te3\tdec\t1\t0\t1\n"
        "s4\te4\tdec\t4\t2\t2\n",
        encoding="utf-8",
    )

    experiments_tsv = tmp_path / "experiments.tsv"
    experiments_tsv.write_text(
        "study_id\texperiment_id\tbody_site\n"
        "s1\te1\tfeces\n"
        "s2\te2\tsaliva\n"
        "s3\te3\tstool\n"
        "s4\te4\tvagina\n",
        encoding="utf-8",
    )

    out_dir = tmp_path / "strata_out"
    partitions = partition_profiles(
        profiles_path=profiles_tsv,
        experiments_path=experiments_tsv,
        output_dir=out_dir,
        strata=PRIMARY_STRATA,
    )

    assert set(partitions.keys()) == {"pooled", "feces", "feces_excluded"}

    pooled_df = pd.read_csv(partitions["pooled"], sep="\t")
    feces_df = pd.read_csv(partitions["feces"], sep="\t")
    feces_ex_df = pd.read_csv(partitions["feces_excluded"], sep="\t")

    assert len(pooled_df) == 4
    assert len(feces_df) == 2
    assert set(feces_df["study_id"]) == {"s1", "s3"}
    assert len(feces_ex_df) == 2
    assert set(feces_ex_df["study_id"]) == {"s2", "s4"}

    # body_site should not be present in output profile tables
    assert "body_site" not in pooled_df.columns
    assert "body_site" not in feces_df.columns
    assert "body_site" not in feces_ex_df.columns


def test_partition_profiles_empty_error(tmp_path: Path) -> None:
    """Verify partition_profiles raises an informative ValueError for empty strata."""
    profiles_tsv = tmp_path / "profiles.tsv"
    profiles_tsv.write_text(
        "study_id\texperiment_id\tinc_dec\ttaxa_count\tTraitA\ns1\te1\tinc\t2\t1\n",
        encoding="utf-8",
    )
    experiments_tsv = tmp_path / "experiments.tsv"
    experiments_tsv.write_text(
        "study_id\texperiment_id\tbody_site\ns1\te1\tfeces\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="Stratum 'saliva' has 0 matching profiles"):
        partition_profiles(
            profiles_path=profiles_tsv,
            experiments_path=experiments_tsv,
            output_dir=tmp_path / "out",
            strata=["saliva"],
        )


def test_compute_annotation_density_deciles() -> None:
    """Verify background annotation density decile binning."""
    df = pd.DataFrame(
        {
            "trait": [f"Trait_{i}" for i in range(20)],
            "trait_total": list(range(1, 21)),
            "best_qvalue": [0.01 if i % 2 == 0 else 0.10 for i in range(20)],
            "n_profiles_tested": [10] * 20,
            "hit_rows": [1 if i % 2 == 0 else 0 for i in range(20)],
            "test_status": [
                "tested_significant" if i % 2 == 0 else "tested_nonsignificant" for i in range(20)
            ],
        }
    )

    deciles = compute_annotation_density_deciles(df, n_deciles=4)
    assert len(deciles) == 4
    assert list(deciles.columns) == [
        "decile",
        "min_trait_total",
        "max_trait_total",
        "total_traits",
        "significant_traits",
        "pct_significant",
    ]
    assert deciles["total_traits"].sum() == 20
    assert deciles["significant_traits"].sum() == 10

    # Empty DataFrame edge case
    empty_res = compute_annotation_density_deciles(pd.DataFrame(), n_deciles=4)
    assert empty_res.empty


def _write_ranked_table(path: Path, n_testable: int, rows: list[dict[str, object]]) -> None:
    """Write a synthetic complete trait summary with its explicit discovery state."""
    df = pd.DataFrame(rows)
    df["hit_rows"] = df["best_qvalue"].le(0.05).astype(int)
    df["test_status"] = "tested_nonsignificant"
    df.loc[df["hit_rows"] > 0, "test_status"] = "tested_significant"
    df.loc[df["n_profiles_tested"] == 0, "test_status"] = "zero_eligible_tests"
    with path.open("w", encoding="utf-8", newline="") as f:
        f.write(f"# Total-Testable-Profiles: {n_testable}\n")
        f.write(f"# Ranked-Traits: {len(df)}\n")
        df.to_csv(f, sep="\t", index=False)


def test_support_status_and_sensitivity_reports(tmp_path: Path) -> None:
    """Verify cross-stratum support status classification and report generation."""
    out_dir = tmp_path / "reports_out"
    out_dir.mkdir(parents=True)

    pooled_tsv = out_dir / "trait_enrichment_traits_overall_ranked_pooled.tsv"
    feces_tsv = out_dir / "trait_enrichment_traits_overall_ranked_feces.tsv"
    feces_ex_tsv = out_dir / "trait_enrichment_traits_overall_ranked_feces_excluded.tsv"

    pooled_rows = [
        {
            "rank": 1,
            "trait": "TraitBoth",
            "best_qvalue": 0.001,
            "mean_fold": 3.0,
            "n_profiles_tested": 10,
            "trait_total": 50,
            "prevalence": 0.2,
        },
        {
            "rank": 2,
            "trait": "TraitFecesOnly",
            "best_qvalue": 0.001,
            "mean_fold": 2.5,
            "n_profiles_tested": 10,
            "trait_total": 30,
            "prevalence": 0.2,
        },
        {
            "rank": 3,
            "trait": "TraitNonFecesOnly",
            "best_qvalue": 0.001,
            "mean_fold": 2.0,
            "n_profiles_tested": 10,
            "trait_total": 20,
            "prevalence": 0.2,
        },
        {
            "rank": 4,
            "trait": "TraitLowSupport",
            "best_qvalue": 0.001,
            "mean_fold": 2.0,
            "n_profiles_tested": 10,
            "trait_total": 10,
            "prevalence": 0.1,
        },
    ]

    feces_rows = [
        {
            "rank": 1,
            "trait": "TraitBoth",
            "best_qvalue": 0.01,
            "mean_fold": 3.5,
            "n_profiles_tested": 5,
            "trait_total": 30,
            "prevalence": 0.4,
        },
        {
            "rank": 2,
            "trait": "TraitFecesOnly",
            "best_qvalue": 0.01,
            "mean_fold": 3.0,
            "n_profiles_tested": 5,
            "trait_total": 25,
            "prevalence": 0.4,
        },
        {
            "rank": 3,
            "trait": "TraitNonFecesOnly",
            "best_qvalue": 0.50,
            "mean_fold": 1.0,
            "n_profiles_tested": 5,
            "trait_total": 20,
            "prevalence": 0.1,
        },
        {
            "rank": 4,
            "trait": "TraitLowSupport",
            "best_qvalue": 0.01,
            "mean_fold": 2.0,
            "n_profiles_tested": 1,
            "trait_total": 5,
            "prevalence": 0.2,
        },
    ]

    feces_ex_rows = [
        {
            "rank": 4,
            "trait": "TraitLowSupport",
            "best_qvalue": float("nan"),
            "mean_fold": float("nan"),
            "n_profiles_tested": 0,
            "trait_total": float("nan"),
            "prevalence": float("nan"),
        },
        {
            "rank": 1,
            "trait": "TraitBoth",
            "best_qvalue": 0.02,
            "mean_fold": 2.5,
            "n_profiles_tested": 5,
            "trait_total": 20,
            "prevalence": 0.4,
        },
        {
            "rank": 2,
            "trait": "TraitNonFecesOnly",
            "best_qvalue": 0.02,
            "mean_fold": 2.2,
            "n_profiles_tested": 5,
            "trait_total": 15,
            "prevalence": 0.4,
        },
        {
            "rank": 3,
            "trait": "TraitFecesOnly",
            "best_qvalue": 0.50,
            "mean_fold": 1.1,
            "n_profiles_tested": 5,
            "trait_total": 10,
            "prevalence": 0.1,
        },
    ]

    _write_ranked_table(pooled_tsv, 10, pooled_rows)
    _write_ranked_table(feces_tsv, 5, feces_rows)
    _write_ranked_table(feces_ex_tsv, 5, feces_ex_rows)

    stratum_outs = {
        "pooled": StratumOutputs(
            stratum="pooled",
            profiles_tsv=out_dir / "p.tsv",
            results_tsv=out_dir / "r.tsv",
            report_txt=out_dir / "rep.txt",
            overall_ranked_tsv=pooled_tsv,
            tested_traits_tsv=pooled_tsv,
            pairs_ranked_tsv=out_dir / "pairs_pooled.tsv",
        ),
        "feces": StratumOutputs(
            stratum="feces",
            profiles_tsv=out_dir / "p.tsv",
            results_tsv=out_dir / "r.tsv",
            report_txt=out_dir / "rep.txt",
            overall_ranked_tsv=feces_tsv,
            tested_traits_tsv=feces_tsv,
            pairs_ranked_tsv=out_dir / "pairs_feces.tsv",
        ),
        "feces_excluded": StratumOutputs(
            stratum="feces_excluded",
            profiles_tsv=out_dir / "p.tsv",
            results_tsv=out_dir / "r.tsv",
            report_txt=out_dir / "rep.txt",
            overall_ranked_tsv=feces_ex_tsv,
            tested_traits_tsv=feces_ex_tsv,
            pairs_ranked_tsv=out_dir / "pairs_feces_excluded.tsv",
        ),
    }

    per_trait_df, strata_df = build_sensitivity_summary_table(stratum_outs)
    assert len(strata_df) == 3
    assert len(per_trait_df) == 4

    status_map = dict(zip(per_trait_df["trait"], per_trait_df["support_status"], strict=True))
    assert status_map["TraitBoth"] == "both_strata_significant"
    assert status_map["TraitFecesOnly"] == "feces_only_significant"
    assert status_map["TraitNonFecesOnly"] == "feces_excluded_only_significant"
    assert status_map["TraitLowSupport"] == "insufficient_support"

    tsv_path, txt_path = generate_sensitivity_reports(stratum_outs, out_dir)
    assert tsv_path.exists()
    assert txt_path.exists()

    tsv_content = tsv_path.read_text(encoding="utf-8")
    txt_content = txt_path.read_text(encoding="utf-8")

    # Verify mandatory precautionary and diagnostic statements
    assert PRECAUTIONARY_NOTE in tsv_content
    assert PRECAUTIONARY_NOTE in txt_content
    assert ANNOTATION_DENSITY_NOTE in tsv_content
    assert ANNOTATION_DENSITY_NOTE in txt_content

    # Verify presence of support status breakdown
    assert "both_strata_significant" in txt_content
    assert "feces_only_significant" in txt_content
    assert "insufficient_support" in txt_content


def test_run_stratification_pipeline_e2e(tmp_path: Path) -> None:
    """End-to-end integration test of the background stratification pipeline."""
    # Create minimal wide profiles with 4 profiles, 2 traits, 2 studies
    # Study 1: Feces (2 profiles: inc and dec)
    # Study 2: Saliva (2 profiles: inc and dec)
    matrix_tsv = tmp_path / "taxa_trait_matrix.tsv"
    matrix_tsv.write_text(
        "taxon_id\tTrait1\tTrait2\n"
        "NCBITaxon:1\t1\t1\n"
        "NCBITaxon:2\t1\t0\n"
        "NCBITaxon:3\t0\t1\n"
        "NCBITaxon:4\t1\t1\n",
        encoding="utf-8",
    )

    profiles_tsv = tmp_path / "study_trait_profiles_wide.tsv"
    profiles_tsv.write_text(
        "study_id\texperiment_id\tinc_dec\ttaxa_count\tTrait1\tTrait2\n"
        "s1\te1\tinc\t4\t3\t2\n"
        "s1\te1\tdec\t4\t1\t3\n"
        "s2\te2\tinc\t4\t3\t1\n"
        "s2\te2\tdec\t4\t2\t2\n",
        encoding="utf-8",
    )

    experiments_tsv = tmp_path / "experiments.tsv"
    experiments_tsv.write_text(
        "study_id\texperiment_id\tbody_site\ns1\te1\tfeces\ns2\te2\tsaliva\n",
        encoding="utf-8",
    )

    out_dir = tmp_path / "derived"

    _stratum_outs, tsv_path, txt_path = run_stratification_pipeline(
        profiles_path=profiles_tsv,
        experiments_path=experiments_tsv,
        taxa_matrix=matrix_tsv,
        output_dir=out_dir,
        p_threshold=1.0,  # capture all tests to verify rankings
        min_fold=0.0,
        fdr=False,
        strata=PRIMARY_STRATA,
    )

    # 1. Check stable filenames
    expected_overall_files = {
        "pooled": out_dir / "trait_enrichment_traits_overall_ranked_pooled.tsv",
        "feces": out_dir / "trait_enrichment_traits_overall_ranked_feces.tsv",
        "feces_excluded": out_dir / "trait_enrichment_traits_overall_ranked_feces_excluded.tsv",
    }
    expected_pair_files = {
        "pooled": out_dir / "trait_enrichment_trait_study_pairs_ranked_pooled.tsv",
        "feces": out_dir / "trait_enrichment_trait_study_pairs_ranked_feces.tsv",
        "feces_excluded": (
            out_dir / "trait_enrichment_trait_study_pairs_ranked_feces_excluded.tsv"
        ),
    }

    for stratum, fpath in expected_overall_files.items():
        assert fpath.exists(), f"Missing overall ranked file: {fpath}"
        header_text = fpath.read_text(encoding="utf-8")
        assert f"# Body-Site-Group: {stratum}" in header_text

        # Verify n_profiles_tested and trait_total in header
        df = pd.read_csv(fpath, sep="\t", comment="#")
        assert "n_profiles_tested" in df.columns
        assert "trait_total" in df.columns

    for stratum, fpath in expected_pair_files.items():
        assert fpath.exists(), f"Missing pair ranked file: {fpath}"
        header_text = fpath.read_text(encoding="utf-8")
        assert f"# Body-Site-Group: {stratum}" in header_text

    # 2. Check summary TSV and TXT
    assert tsv_path.exists()
    assert txt_path.exists()

    txt_content = txt_path.read_text(encoding="utf-8")
    assert PRECAUTIONARY_NOTE in txt_content
    assert ANNOTATION_DENSITY_NOTE in txt_content


def test_run_stratification_pipeline_single_site(tmp_path: Path) -> None:
    """Verify run_stratification_pipeline with a single specific site (analyze-site)."""
    matrix_tsv = tmp_path / "taxa_trait_matrix.tsv"
    matrix_tsv.write_text(
        "taxon_id\tTrait1\tTrait2\n"
        "NCBITaxon:1\t1\t1\n"
        "NCBITaxon:2\t1\t0\n"
        "NCBITaxon:3\t0\t1\n"
        "NCBITaxon:4\t1\t1\n",
        encoding="utf-8",
    )

    profiles_tsv = tmp_path / "study_trait_profiles_wide.tsv"
    profiles_tsv.write_text(
        "study_id\texperiment_id\tinc_dec\ttaxa_count\tTrait1\tTrait2\n"
        "s1\te1\tinc\t4\t3\t2\n"
        "s1\te1\tdec\t4\t1\t3\n"
        "s2\te2\tinc\t4\t3\t1\n"
        "s2\te2\tdec\t4\t2\t2\n",
        encoding="utf-8",
    )

    experiments_tsv = tmp_path / "experiments.tsv"
    experiments_tsv.write_text(
        "study_id\texperiment_id\tbody_site\ns1\te1\tfeces\ns2\te2\tsaliva\n",
        encoding="utf-8",
    )

    out_dir = tmp_path / "derived" / "sites" / "feces"

    stratum_outs, _tsv_path, _txt_path = run_stratification_pipeline(
        profiles_path=profiles_tsv,
        experiments_path=experiments_tsv,
        taxa_matrix=matrix_tsv,
        output_dir=out_dir,
        p_threshold=1.0,
        min_fold=0.0,
        fdr=False,
        strata=["feces"],
    )

    assert "feces" in stratum_outs
    feces_out = stratum_outs["feces"]
    assert feces_out.overall_ranked_tsv.exists()
    assert feces_out.pairs_ranked_tsv.exists()

    df = pd.read_csv(feces_out.overall_ranked_tsv, sep="\t", comment="#")
    assert "body_site_group" in df.columns
    assert (df["body_site_group"] == "feces").all()


def test_cli_wrapper() -> None:
    """Verify run_stratified_enrichment.py CLI wrapper responds to --help."""
    res = subprocess.run(
        [sys.executable, str(CLI_WRAPPER), "--help"],
        capture_output=True,
        text=True,
    )
    assert res.returncode == 0
    assert "Evaluate ORA background stratification" in res.stdout


def test_complete_test_universe_reaches_density_diagnostic(tmp_path: Path) -> None:
    """Actual ORA -> exporter -> diagnostic retains null and single-hit traits (#664)."""
    matrix = tmp_path / "matrix.tsv"
    traits = ["Signal", "Null", "Singleton", "ObservedZero", "Unsupported"]
    pd.DataFrame(
        {
            "taxon_id": [f"NCBITaxon:{i}" for i in range(100)],
            "Signal": [int(i < 10) for i in range(100)],
            "Null": [int(i < 50) for i in range(100)],
            "Singleton": [int(i < 10) for i in range(100)],
            "ObservedZero": [int(i == 0) for i in range(100)],
            "Unsupported": [0] * 100,
        }
    ).to_csv(matrix, sep="\t", index=False)
    profiles = tmp_path / "profiles.tsv"
    rows = []
    for study in range(4):
        for direction in ("inc", "dec"):
            rows.append(
                [
                    f"s{study}",
                    "e1",
                    direction,
                    10,
                    10 if direction == "inc" else 0,
                    5,
                    10 if study == 0 and direction == "inc" else 0,
                    0,
                    0,
                ]
            )
    # Neither zero-draw paired profiles nor unpaired profiles belong in the denominator.
    rows.extend(
        [
            ["s4", "e1", "inc", 0, 0, 0, 0, 0, 0],
            ["s4", "e1", "dec", 0, 0, 0, 0, 0, 0],
            ["s5", "e1", "inc", 10, 10, 5, 10, 0, 0],
        ]
    )
    pd.DataFrame(
        rows, columns=["study_id", "experiment_id", "inc_dec", "taxa_count", *traits]
    ).to_csv(profiles, sep="\t", index=False)
    experiments = tmp_path / "experiments.tsv"
    pd.DataFrame(
        {
            "study_id": [f"s{i}" for i in range(6)],
            "experiment_id": ["e1"] * 6,
            "body_site": ["feces", "feces", "saliva", "saliva", "feces", "saliva"],
        }
    ).to_csv(experiments, sep="\t", index=False)
    outputs, summary_path, _ = run_stratification_pipeline(
        profiles_path=profiles,
        experiments_path=experiments,
        taxa_matrix=matrix,
        output_dir=tmp_path / "out",
    )
    pooled = outputs["pooled"]
    assert pooled.tested_traits_tsv is not None
    tested = pd.read_csv(pooled.tested_traits_tsv, sep="\t", comment="#").set_index("trait")
    ranked = pd.read_csv(pooled.overall_ranked_tsv, sep="\t", comment="#").set_index("trait")
    assert set(ranked.index) == {"Signal"}
    assert ranked.loc["Signal", "prevalence_summary"] == "4 / 8"
    assert ranked.loc["Signal", "prevalence"] == 0.5
    assert ranked.loc["Signal", "trait_total"] == 10  # Reference K, not profile sum 40.
    assert tested.loc["Singleton", "hit_rows"] == 1
    for trait in ("Null", "ObservedZero"):
        assert tested.loc[trait, "test_status"] == "tested_nonsignificant"
        assert tested.loc[trait, "n_profiles_tested"] == 8
        assert tested.loc[trait, "prevalence"] == 0
        assert tested.loc[trait, "hit_rows"] == 0
    assert tested.loc["Unsupported", "test_status"] == "zero_eligible_tests"
    assert tested.loc["Unsupported", "n_profiles_tested"] == 0
    assert pd.isna(tested.loc["Unsupported", "prevalence"])
    assert pooled.all_results_tsv is not None
    all_results = pd.read_csv(pooled.all_results_tsv, sep="\t", comment="#")
    assert len(all_results) == tested["n_profiles_tested"].sum() == 32
    pairs = pd.read_csv(pooled.pairs_ranked_tsv, sep="\t", comment="#")
    assert len(pairs) == 5  # The public pair ranking remains discovery-only.
    deciles = pd.read_csv(tmp_path / "out" / "annotation_density_deciles.tsv", sep="\t")
    assert deciles["total_traits"].sum() == 4
    assert deciles["significant_traits"].sum() == 2
    assert (deciles["pct_significant"] == 0).any()
    sensitivity = pd.read_csv(summary_path, sep="\t", comment="#").set_index("trait")
    assert sensitivity.loc["Null", "support_status"] == "neither_significant"
    assert sensitivity.loc["Unsupported", "support_status"] == "insufficient_support"
    assert pd.isna(sensitivity.loc["Unsupported", "feces_significant"])

    # Missing analysis must not become evidence of a tested negative result.
    missing_path = outputs["feces_excluded"].tested_traits_tsv
    assert missing_path is not None
    missing_path.unlink()
    missing, _ = build_sensitivity_summary_table(outputs)
    null = missing.set_index("trait").loc["Null"]
    assert null["support_status"] == "unavailable_comparison"
    assert null["feces_excluded_test_status"] == "unavailable_analysis"
    assert pd.isna(null["feces_excluded_prevalence"])


def test_density_rejects_selected_rankings() -> None:
    """A missing statistical universe cannot silently become a 100% discovery rate."""
    with pytest.raises(ValueError, match="complete trait summary"):
        compute_annotation_density_deciles(
            pd.DataFrame({"trait_total": [10], "best_qvalue": [0.01]})
        )


def test_report_only_export_preserves_unavailable_denominator(tmp_path: Path) -> None:
    """Legacy text discoveries cannot establish how many eligible profiles were tested."""
    report = tmp_path / "legacy_report.txt"
    report.write_text(
        "FDR q <= 0.001\nMin Fold >= 1.5\n"
        "s1\te1\tinc\n- TraitA, 0.0001, obs=5 exp=1 fold=5\n"
        "s2\te1\tinc\n- TraitA, 0.0002, obs=4 exp=1 fold=4\n"
    )
    output = tmp_path / "ranked.tsv"
    result = subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / "scripts" / "export_trait_enrichment_rankings.py"),
            "--input",
            str(report),
            "--results",
            str(tmp_path / "missing.tsv"),
            "--overall-output",
            str(output),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    row = pd.read_csv(output, sep="\t", comment="#").iloc[0]
    assert row["hit_rows"] == 2
    assert row["prevalence_summary"] == "2 / NA"
    assert pd.isna(row["n_profiles_tested"])
    assert pd.isna(row["prevalence"])


@pytest.mark.parametrize("draw_count", [0, 2])
def test_zero_discoveries_are_distinct_from_zero_tests(tmp_path: Path, draw_count: int) -> None:
    """An empty selected ranking never erases the analyzer's eligible universe."""
    matrix = tmp_path / "matrix.tsv"
    matrix.write_text("taxon_id\tNull\nNCBITaxon:1\t1\nNCBITaxon:2\t0\n")
    profiles = tmp_path / "profiles.tsv"
    profiles.write_text(
        "study_id\texperiment_id\tinc_dec\ttaxa_count\tNull\n"
        f"s1\te1\tinc\t{draw_count}\t0\n"
        f"s1\te1\tdec\t{draw_count}\t0\n"
    )
    experiments = tmp_path / "experiments.tsv"
    experiments.write_text("study_id\texperiment_id\tbody_site\ns1\te1\tfeces\n")
    outputs, _, _ = run_stratification_pipeline(
        profiles_path=profiles,
        experiments_path=experiments,
        taxa_matrix=matrix,
        output_dir=tmp_path / "out",
        strata=["pooled"],
    )
    tested_path = outputs["pooled"].tested_traits_tsv
    assert tested_path is not None
    row = pd.read_csv(tested_path, sep="\t", comment="#").iloc[0]
    assert f"# Total-Testable-Profiles: {2 if draw_count else 0}\n" in tested_path.read_text()
    deciles = pd.read_csv(tmp_path / "out" / "annotation_density_deciles.tsv", sep="\t")
    assert pd.read_csv(outputs["pooled"].overall_ranked_tsv, sep="\t", comment="#").empty
    assert row["hit_rows"] == 0
    if draw_count:
        assert row["test_status"] == "tested_nonsignificant"
        assert row["n_profiles_tested"] == 2
        assert row["prevalence"] == 0
        assert deciles["total_traits"].sum() == 1
        assert deciles["pct_significant"].sum() == 0
    else:
        assert row["test_status"] == "zero_eligible_tests"
        assert row["n_profiles_tested"] == 0
        assert pd.isna(row["prevalence"])
        assert deciles.empty
