import subprocess
import sys
from pathlib import Path


def test_enrichment_pipeline_fdr_propagation(tmp_path: Path) -> None:
    # 1. Create a dummy trait profile TSV
    input_tsv = tmp_path / "input_profiles.tsv"
    input_tsv.write_text(
        "study_id\texperiment_id\tinc_dec\tTraitA\tTraitB\n"
        "s1\te1\tinc\t10\t1\n"
        "s1\te1\tdec\t1\t10\n"
        "s1\te2\tinc\t0\t5\n"
        "s1\te2\tdec\t5\t0\n"
        "s2\te1\tinc\t10\t0\n"
        "s2\te1\tdec\t0\t10\n"
        "s3\te1\tinc\t0\t0\n"
        "s3\te1\tdec\t0\t0\n",
        encoding="utf-8",
    )

    report_txt = tmp_path / "report.txt"
    results_tsv = tmp_path / "results.tsv"

    # 2. Run trait_enrichment.py with default args (FDR=True)
    subprocess.run(
        [
            sys.executable,
            "scripts/trait_enrichment.py",
            "--input",
            str(input_tsv),
            "--output",
            str(report_txt),
            "--tsv-output",
            str(results_tsv),
            "--p-threshold",
            "1.0",
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    # Verify report header has FDR
    report_content = report_txt.read_text(encoding="utf-8")
    assert "FDR q < 1" in report_content

    overall_out = tmp_path / "overall.tsv"
    pair_out = tmp_path / "pairs.tsv"

    subprocess.run(
        [
            sys.executable,
            "scripts/export_trait_enrichment_rankings.py",
            "--input",
            str(report_txt),
            "--results",
            str(results_tsv),
            "--overall-output",
            str(overall_out),
            "--pair-output",
            str(pair_out),
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    # 4. Verify output content
    with overall_out.open("r", encoding="utf-8") as f:
        lines = f.readlines()
        assert any("# Stat-Type: fdr_q" in line for line in lines)
        header = next(line for line in lines if not line.startswith("#")).strip().split("\t")
        assert "best_qvalue" in header

    with pair_out.open("r", encoding="utf-8") as f:
        lines = f.readlines()
        assert any("# Stat-Type: fdr_q" in line for line in lines)
        header = next(line for line in lines if not line.startswith("#")).strip().split("\t")
        assert "q_value" in header


def test_enrichment_pipeline_no_fdr_propagation(tmp_path: Path) -> None:
    input_tsv = tmp_path / "input_profiles.tsv"
    input_tsv.write_text(
        "study_id\texperiment_id\tinc_dec\tTraitA\tTraitB\n"
        "s1\te1\tinc\t10\t1\n"
        "s1\te1\tdec\t1\t10\n",
        encoding="utf-8",
    )

    report_txt = tmp_path / "report.txt"
    results_tsv = tmp_path / "results.tsv"

    # Run with --no-fdr
    subprocess.run(
        [
            sys.executable,
            "scripts/trait_enrichment.py",
            "--input",
            str(input_tsv),
            "--output",
            str(report_txt),
            "--tsv-output",
            str(results_tsv),
            "--no-fdr",
            "--p-threshold",
            "1.0",
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    overall_out = tmp_path / "overall.tsv"

    subprocess.run(
        [
            sys.executable,
            "scripts/export_trait_enrichment_rankings.py",
            "--input",
            str(report_txt),
            "--results",
            str(results_tsv),
            "--overall-output",
            str(overall_out),
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    with overall_out.open("r", encoding="utf-8") as f:
        lines = f.readlines()
        assert any("# Stat-Type: p_value" in line for line in lines)
        header = next(line for line in lines if not line.startswith("#")).strip().split("\t")
        assert "best_pvalue" in header


def test_enrichment_pipeline_censorship_detection(tmp_path: Path) -> None:
    input_tsv = tmp_path / "input_profiles.tsv"
    # Ensure TraitA and TraitB appear in at least 2 profiles to pass prevalence filter
    input_tsv.write_text(
        "study_id\texperiment_id\tinc_dec\tTraitA\tTraitB\tTraitC\n"
        "s1\te1\tinc\t100\t100\t100\n"
        "s1\te1\tdec\t1\t1\t1\n"
        "s2\te1\tinc\t100\t100\t100\n"
        "s2\te1\tdec\t1\t1\t1\n",
        encoding="utf-8",
    )

    report_txt = tmp_path / "report.txt"
    results_tsv = tmp_path / "results.tsv"

    # Run with truncation cap of 1 trait
    subprocess.run(
        [
            sys.executable,
            "scripts/trait_enrichment.py",
            "--input",
            str(input_tsv),
            "--output",
            str(report_txt),
            "--tsv-output",
            str(results_tsv),
            "--max-traits-per-row",
            "1",
            "--p-threshold",
            "1.0",
            "--min-fold",
            "0",
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    # Force using the report
    results_tsv.unlink()
    overall_out = tmp_path / "overall.tsv"

    subprocess.run(
        [
            sys.executable,
            "scripts/export_trait_enrichment_rankings.py",
            "--input",
            str(report_txt),
            "--results",
            str(tmp_path / "notfound.tsv"),
            "--overall-output",
            str(overall_out),
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    with overall_out.open("r", encoding="utf-8") as f:
        lines = f.readlines()
        assert "# Max-Traits-Per-Profile-Cap: 1" in "".join(lines)
        assert "WARNING: Totals are censored" in "".join(lines)
        header = next(line for line in lines if not line.startswith("#")).strip().split("\t")
        assert "n_censored_profiles" in header

        data_rows = [line.strip().split("\t") for line in lines if not line.startswith("#")][1:]
        censored_idx = header.index("n_censored_profiles")
        # All hits in s1|e1|inc and s2|e1|inc should be marked censored because they reached cap 1
        for row in data_rows:
            assert int(row[censored_idx]) >= 1


def test_enrichment_pipeline_profile_fdr_propagation(tmp_path: Path) -> None:
    input_tsv = tmp_path / "input_profiles.tsv"
    input_tsv.write_text(
        "study_id\texperiment_id\tinc_dec\tTraitA\tTraitB\n"
        "s1\te1\tinc\t10\t1\n"
        "s1\te1\tdec\t1\t10\n"
        "s1\te2\tinc\t10\t1\n"
        "s1\te2\tdec\t1\t10\n"
        "s1\te3\tinc\t10\t1\n"
        "s1\te3\tdec\t1\t10\n"
        "s1\te4\tinc\t10\t1\n"
        "s1\te4\tdec\t1\t10\n"
        "s1\te5\tinc\t10\t1\n"
        "s1\te5\tdec\t1\t10\n",
        encoding="utf-8",
    )

    report_txt = tmp_path / "report.txt"
    results_tsv = tmp_path / "results.tsv"

    subprocess.run(
        [
            sys.executable,
            "scripts/trait_enrichment.py",
            "--input",
            str(input_tsv),
            "--output",
            str(report_txt),
            "--tsv-output",
            str(results_tsv),
            "--fdr-scope",
            "profile",
            "--p-threshold",
            "1.0",
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    overall_out = tmp_path / "overall.tsv"

    subprocess.run(
        [
            sys.executable,
            "scripts/export_trait_enrichment_rankings.py",
            "--input",
            str(report_txt),
            "--results",
            str(results_tsv),
            "--overall-output",
            str(overall_out),
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    with overall_out.open("r", encoding="utf-8") as f:
        lines = f.readlines()
        assert "# FDR-Scope: profile" in "".join(lines)
        header = next(line for line in lines if not line.startswith("#")).strip().split("\t")
        assert "prevalence_summary" in header
        data_row = [line for line in lines if not line.startswith("#")][1].strip().split("\t")
        prev_idx = header.index("prevalence_summary")
        assert "/" in data_row[prev_idx]
