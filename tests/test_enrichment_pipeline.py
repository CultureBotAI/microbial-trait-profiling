import csv
import subprocess
import sys
from pathlib import Path

def test_enrichment_pipeline_fdr_propagation(tmp_path: Path):
    # 1. Create a dummy trait profile TSV
    input_tsv = tmp_path / "input_profiles.tsv"
    # study_id, experiment_id, inc_dec, TraitA, TraitB
    # We need enough data to trigger some "hits"
    # Hypergeometric test: pop_size, pop_successes, draws, count
    # Let's make TraitA very enriched in one profile.
    input_tsv.write_text(
        "study_id\texperiment_id\tinc_dec\tTraitA\tTraitB\n"
        "s1\te1\tinc\t10\t1\n"
        "s1\te2\tinc\t0\t5\n"
        "s2\te1\tinc\t10\t0\n"
        "s3\te1\tinc\t0\t0\n",
        encoding="utf-8"
    )

    report_txt = tmp_path / "report.txt"
    results_tsv = tmp_path / "results.tsv"

    # 2. Run trait_enrichment.py with default args (FDR=True)
    subprocess.run([
        sys.executable, "bsdbng/scripts/trait_enrichment.py",
        "--input", str(input_tsv),
        "--output", str(report_txt),
        "--tsv-output", str(results_tsv),
        "--p-threshold", "1.0", # Ensure everything is "significant" for the test
    ], check=True, capture_output=True, text=True)

    # Verify report header has FDR
    report_content = report_txt.read_text(encoding="utf-8")
    assert "FDR q < 1" in report_content

    # 3. Run export_trait_enrichment_rankings.py
    # We don't specify output names, so it should use dynamic ones in data/derived
    # But for testing in a tmp_path, we should probably specify them or check where they went.
    # The script uses Path("data/derived") by default.
    # We can override them.
    
    # Let's run it and let it use default paths, but we'll mock the data/derived if needed?
    # Actually, let's just pass the paths to be sure.
    overall_out = tmp_path / "overall.tsv"
    pair_out = tmp_path / "pairs.tsv"
    
    subprocess.run([
        sys.executable, "bsdbng/scripts/export_trait_enrichment_rankings.py",
        "--input", str(report_txt),
        "--results", str(results_tsv),
        "--overall-output", str(overall_out),
        "--pair-output", str(pair_out),
    ], check=True, capture_output=True, text=True)

    # 4. Verify output content
    with overall_out.open("r", encoding="utf-8") as f:
        lines = f.readlines()
        assert any("# Stat-Type: fdr_q" in l for l in lines)
        
        # Check column names (first non-comment line)
        header = [l for l in lines if not l.startswith("#")][0].strip().split("\t")
        assert "best_qvalue" in header
        assert "median_qvalue" in header
        assert "best_pvalue" in header
        assert "prevalence_summary" in header
        
    with pair_out.open("r", encoding="utf-8") as f:
        lines = f.readlines()
        assert any("# Stat-Type: fdr_q" in l for l in lines)
        header = [l for l in lines if not l.startswith("#")][0].strip().split("\t")
        assert "q_value" in header

def test_enrichment_pipeline_no_fdr_propagation(tmp_path: Path):
    input_tsv = tmp_path / "input_profiles.tsv"
    input_tsv.write_text(
        "study_id\texperiment_id\tinc_dec\tTraitA\tTraitB\n"
        "s1\te1\tinc\t10\t1\n",
        encoding="utf-8"
    )

    report_txt = tmp_path / "report.txt"
    results_tsv = tmp_path / "results.tsv"

    # Run with --no-fdr
    subprocess.run([
        sys.executable, "bsdbng/scripts/trait_enrichment.py",
        "--input", str(input_tsv),
        "--output", str(report_txt),
        "--tsv-output", str(results_tsv),
        "--no-fdr",
        "--p-threshold", "1.0",
    ], check=True, capture_output=True, text=True)

    report_content = report_txt.read_text(encoding="utf-8")
    assert "p < 1" in report_content

    overall_out = tmp_path / "overall.tsv"
    
    subprocess.run([
        sys.executable, "bsdbng/scripts/export_trait_enrichment_rankings.py",
        "--input", str(report_txt),
        "--results", str(results_tsv),
        "--overall-output", str(overall_out),
    ], check=True, capture_output=True, text=True)

    with overall_out.open("r", encoding="utf-8") as f:
        lines = f.readlines()
        assert any("# Stat-Type: p_value" in l for l in lines)
        header = [l for l in lines if not l.startswith("#")][0].strip().split("\t")
        assert "best_pvalue" in header
        assert "median_pvalue" in header
        assert "best_qvalue" in header # Should still be there but maybe 1.0

def test_enrichment_pipeline_censorship_detection(tmp_path: Path):
    # Create input that will result in more than 1 hit
    input_tsv = tmp_path / "input_profiles.tsv"
    input_tsv.write_text(
        "study_id\texperiment_id\tinc_dec\tTraitA\tTraitB\tTraitC\n"
        "s1\te1\tinc\t10\t10\t10\n",
        encoding="utf-8"
    )

    report_txt = tmp_path / "report.txt"
    results_tsv = tmp_path / "results.tsv"

    # 1. Run with truncation cap of 2 traits
    subprocess.run([
        sys.executable, "bsdbng/scripts/trait_enrichment.py",
        "--input", str(input_tsv),
        "--output", str(report_txt),
        "--tsv-output", str(results_tsv),
        "--max-traits-per-row", "2",
        "--p-threshold", "1.0",
    ], check=True, capture_output=True, text=True)

    # 2. Run exporter forcing it to use the truncated report (results_tsv deleted)
    results_tsv.unlink()
    overall_out = tmp_path / "overall.tsv"
    
    subprocess.run([
        sys.executable, "bsdbng/scripts/export_trait_enrichment_rankings.py",
        "--input", str(report_txt),
        "--overall-output", str(overall_out),
    ], check=True, capture_output=True, text=True)

    with overall_out.open("r", encoding="utf-8") as f:
        lines = f.readlines()
        assert "# Max-Traits-Per-Profile-Cap: 2" in "".join(lines)
        assert "WARNING: Totals are censored" in "".join(lines)
        
        header = [l for l in lines if not l.startswith("#")][0].strip().split("\t")
        assert "n_censored_profiles" in header
        
        # Check first data row
        data_row = [l for l in lines if not l.startswith("#")][1].strip().split("\t")
        censored_idx = header.index("n_censored_profiles")
        assert int(data_row[censored_idx]) == 1 # Profile s1|e1|inc reached the cap

def test_enrichment_pipeline_profile_fdr_propagation(tmp_path: Path):
    input_tsv = tmp_path / "input_profiles.tsv"
    input_tsv.write_text(
        "study_id\texperiment_id\tinc_dec\tTraitA\tTraitB\n"
        "s1\te1\tinc\t10\t1\n"
        "s1\te2\tinc\t10\t1\n"
        "s1\te3\tinc\t0\t10\n"
        "s1\te4\tinc\t0\t10\n"
        "s1\te5\tinc\t0\t10\n",
        encoding="utf-8"
    )

    report_txt = tmp_path / "report.txt"
    results_tsv = tmp_path / "results.tsv"

    # Run with --fdr-scope profile
    subprocess.run([
        sys.executable, "bsdbng/scripts/trait_enrichment.py",
        "--input", str(input_tsv),
        "--output", str(report_txt),
        "--tsv-output", str(results_tsv),
        "--fdr-scope", "profile",
        "--p-threshold", "1.0",
    ], check=True, capture_output=True, text=True)

    overall_out = tmp_path / "overall.tsv"
    
    subprocess.run([
        sys.executable, "bsdbng/scripts/export_trait_enrichment_rankings.py",
        "--input", str(report_txt),
        "--results", str(results_tsv),
        "--overall-output", str(overall_out),
    ], check=True, capture_output=True, text=True)

    with overall_out.open("r", encoding="utf-8") as f:
        lines = f.readlines()
        assert "# FDR-Scope: profile" in "".join(lines)
        
        header = [l for l in lines if not l.startswith("#")][0].strip().split("\t")
        assert "prevalence_summary" in header
        
        # Check data row for prevalence summary format "X / Y"
        data_row = [l for l in lines if not l.startswith("#")][1].strip().split("\t")
        prev_idx = header.index("prevalence_summary")
        assert "/" in data_row[prev_idx]
