import csv
import subprocess
import sys
from pathlib import Path

import pytest


def test_enrichment_pipeline_fdr_propagation(tmp_path: Path) -> None:
    # 1. Create a dummy trait profile TSV
    input_tsv = tmp_path / "input_profiles.tsv"
    input_tsv.write_text(
        "study_id\texperiment_id\tinc_dec\ttaxa_count\tTraitA\tTraitB\n"
        "s1\te1\tinc\t10\t10\t1\n"
        "s1\te1\tdec\t10\t1\t10\n"
        "s1\te2\tinc\t5\t0\t5\n"
        "s1\te2\tdec\t5\t5\t0\n"
        "s2\te1\tinc\t10\t10\t0\n"
        "s2\te1\tdec\t10\t0\t10\n"
        "s3\te1\tinc\t2\t0\t0\n"
        "s3\te1\tdec\t2\t0\t0\n",
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


def test_enrichment_hypergeometric_taxon_sampling_values(tmp_path: Path) -> None:
    # Test specific 2x2 contingency table values
    # TraitA: present in 1 taxon in profile s1/e1/inc out of 2 taxa (k=1, N=2)
    # Background universe across paired profiles: M=4 taxa total, TraitA present in K=2 taxa.
    # hypergeom_sf(1, 4, 2, 2)
    input_tsv = tmp_path / "input_profiles.tsv"
    input_tsv.write_text(
        "study_id\texperiment_id\tinc_dec\ttaxa_count\tTraitA\tTraitB\n"
        "s1\te1\tinc\t2\t1\t0\n"
        "s1\te1\tdec\t2\t1\t2\n",
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
            "--p-threshold",
            "1.0",
            "--min-fold",
            "0",
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    with results_tsv.open("r", encoding="utf-8") as f:
        lines = f.readlines()
        assert any("# Background-Taxa-Universe-M: 4" in line for line in lines)
        header = next(line for line in lines if not line.startswith("#")).strip().split("\t")
        p_idx = header.index("p_value")
        fold_idx = header.index("fold")
        # TraitA row in s1/e1/inc: k=1, N=2, K=2, M=4.
        # Hypergeometric(M=4, K=2, N=2):
        # P(X=0) = comb(2,0)*comb(2,2)/comb(4,2) = 1/6
        # P(X=1) = comb(2,1)*comb(2,1)/comb(4,2) = 4/6 = 2/3
        # P(X=2) = comb(2,2)*comb(2,0)/comb(4,2) = 1/6
        # P(X >= 1) = 4/6 + 1/6 = 5/6 = 0.83333333...
        data_rows = [line.strip().split("\t") for line in lines if not line.startswith("#")][1:]
        trait_a_row = next(
            r
            for r in data_rows
            if r[header.index("trait")] == "TraitA" and r[header.index("inc_dec")] == "inc"
        )
        p_val = float(trait_a_row[p_idx])
        assert abs(p_val - 5 / 6) < 1e-4
        # Expected = 2 * (2 / 4) = 1.0, observed = 1, fold = 1.0
        fold_val = float(trait_a_row[fold_idx])
        assert abs(fold_val - 1.0) < 1e-4


def test_depletion_hypergeometric_taxon_sampling_values(tmp_path: Path) -> None:
    # Test specific 2x2 contingency table values for depletion
    # TraitA: present in 0 taxa in profile s1/e1/inc out of 2 taxa (k=0, N=2)
    # Background universe: M=4 taxa total, TraitA present in K=2 taxa.
    # hypergeom.cdf(0, 4, 2, 2) == P(X <= 0) == comb(2,0)*comb(2,2)/comb(4,2) == 1/6 == 0.1666666...
    input_tsv = tmp_path / "input_profiles.tsv"
    input_tsv.write_text(
        "study_id\texperiment_id\tinc_dec\ttaxa_count\tTraitA\tTraitB\n"
        "s1\te1\tinc\t2\t0\t2\n"
        "s1\te1\tdec\t2\t2\t0\n",
        encoding="utf-8",
    )
    report_txt = tmp_path / "depletion_report.txt"
    results_tsv = tmp_path / "depletion_results.tsv"

    subprocess.run(
        [
            sys.executable,
            "scripts/trait_depletion.py",
            "--input",
            str(input_tsv),
            "--output",
            str(report_txt),
            "--tsv-output",
            str(results_tsv),
            "--p-threshold",
            "1.0",
            "--max-fold",
            "2.0",
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    with results_tsv.open("r", encoding="utf-8") as f:
        lines = f.readlines()
        assert any("# Background-Taxa-Universe-M: 4" in line for line in lines)
        header = next(line for line in lines if not line.startswith("#")).strip().split("\t")
        p_idx = header.index("p_value")
        data_rows = [line.strip().split("\t") for line in lines if not line.startswith("#")][1:]
        trait_a_row = next(
            r
            for r in data_rows
            if r[header.index("trait")] == "TraitA" and r[header.index("inc_dec")] == "inc"
        )
        p_val = float(trait_a_row[p_idx])
        assert abs(p_val - 1 / 6) < 1e-4


def test_enrichment_pipeline_no_fdr_propagation(tmp_path: Path) -> None:
    input_tsv = tmp_path / "input_profiles.tsv"
    input_tsv.write_text(
        "study_id\texperiment_id\tinc_dec\ttaxa_count\tTraitA\tTraitB\n"
        "s1\te1\tinc\t10\t10\t1\n"
        "s1\te1\tdec\t10\t1\t10\n",
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
        "study_id\texperiment_id\tinc_dec\ttaxa_count\tTraitA\tTraitB\tTraitC\n"
        "s1\te1\tinc\t100\t100\t100\t100\n"
        "s1\te1\tdec\t100\t1\t1\t1\n"
        "s2\te1\tinc\t100\t100\t100\t100\n"
        "s2\te1\tdec\t100\t1\t1\t1\n",
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
        "study_id\texperiment_id\tinc_dec\ttaxa_count\tTraitA\tTraitB\n"
        "s1\te1\tinc\t10\t10\t1\n"
        "s1\te1\tdec\t10\t1\t10\n"
        "s1\te2\tinc\t10\t10\t1\n"
        "s1\te2\tdec\t10\t1\t10\n"
        "s1\te3\tinc\t10\t10\t1\n"
        "s1\te3\tdec\t10\t1\t10\n"
        "s1\te4\tinc\t10\t10\t1\n"
        "s1\te4\tdec\t10\t1\t10\n"
        "s1\te5\tinc\t10\t10\t1\n"
        "s1\te5\tdec\t10\t1\t10\n",
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


def test_repeated_taxon_fixture_matrix_wiring(tmp_path: Path) -> None:
    """Verify issue #131: 4 unique taxa repeated across 2 paired experiments report M=4, p=1/6

    with the explicit unique-taxon matrix, vs M=8, p=0.2142857 without matrix.
    Also verifies that a nonexistent matrix path fails instead of falling back silently.
    """
    matrix_tsv = tmp_path / "taxa_trait_matrix.tsv"
    matrix_tsv.write_text(
        "taxon_id\tTraitA\nNCBITaxon:1\t1\nNCBITaxon:2\t1\nNCBITaxon:3\t0\nNCBITaxon:4\t0\n",
        encoding="utf-8",
    )

    profiles_tsv = tmp_path / "profiles.tsv"
    profiles_tsv.write_text(
        "study_id\texperiment_id\tinc_dec\ttaxa_count\tTraitA\n"
        "s1\te1\tinc\t2\t2\n"
        "s1\te1\tdec\t2\t0\n"
        "s2\te2\tinc\t2\t2\n"
        "s2\te2\tdec\t2\t0\n",
        encoding="utf-8",
    )

    report_matrix = tmp_path / "report_matrix.txt"
    results_matrix = tmp_path / "results_matrix.tsv"

    # 1. With explicit unique-taxon matrix: M=4, K=2, N=2, k=2
    # P[X>=2] = comb(2,2)*comb(2,0)/comb(4,2) = 1/6 ~ 0.1666667
    subprocess.run(
        [
            sys.executable,
            "scripts/trait_enrichment.py",
            "--input",
            str(profiles_tsv),
            "--taxa-matrix",
            str(matrix_tsv),
            "--output",
            str(report_matrix),
            "--tsv-output",
            str(results_matrix),
            "--no-fdr",
            "--p-threshold",
            "1.0",
            "--min-fold",
            "0.0",
        ],
        check=True,
    )

    content_matrix = results_matrix.read_text(encoding="utf-8")
    assert "# Background-Taxa-Universe-M: 4" in content_matrix
    assert "# Background-Mode: TaxaMatrix" in content_matrix
    # p_value should be 1/6 (approx 0.166667)
    with results_matrix.open("r", encoding="utf-8") as f:
        reader = csv.DictReader((line for line in f if not line.startswith("#")), delimiter="\t")
        rows_matrix = list(reader)
        p_val_matrix = float(rows_matrix[0]["p_value"])
        assert abs(p_val_matrix - 1 / 6) < 1e-5

    # 2. Without matrix (empirical sum fallback): M=8, K=4, N=2, k=2
    # P[X>=2] = comb(4,2)*comb(4,0)/comb(8,2) = 6/28 = 3/14 ~ 0.2142857
    report_nomatrix = tmp_path / "report_nomatrix.txt"
    results_nomatrix = tmp_path / "results_nomatrix.tsv"
    subprocess.run(
        [
            sys.executable,
            "scripts/trait_enrichment.py",
            "--input",
            str(profiles_tsv),
            "--output",
            str(report_nomatrix),
            "--tsv-output",
            str(results_nomatrix),
            "--no-fdr",
            "--p-threshold",
            "1.0",
            "--min-fold",
            "0.0",
        ],
        check=True,
    )
    content_nomatrix = results_nomatrix.read_text(encoding="utf-8")
    assert "# Background-Taxa-Universe-M: 8" in content_nomatrix
    assert "# Background-Mode: EmpiricalProfileSum" in content_nomatrix
    with results_nomatrix.open("r", encoding="utf-8") as f:
        reader = csv.DictReader((line for line in f if not line.startswith("#")), delimiter="\t")
        rows_nomatrix = list(reader)
        p_val_nomatrix = float(rows_nomatrix[0]["p_value"])
        assert abs(p_val_nomatrix - 3 / 14) < 1e-5

    # 3. Nonexistent matrix path must fail loudly and not fall back silently
    res_bad = subprocess.run(
        [
            sys.executable,
            "scripts/trait_enrichment.py",
            "--input",
            str(profiles_tsv),
            "--taxa-matrix",
            str(tmp_path / "nonexistent_matrix.tsv"),
            "--output",
            str(report_matrix),
        ],
        capture_output=True,
        text=True,
    )
    assert res_bad.returncode != 0
    assert "Specified taxa matrix file does not exist" in res_bad.stderr

    # 4. Same test on trait_depletion.py
    report_dep_matrix = tmp_path / "dep_report_matrix.txt"
    results_dep_matrix = tmp_path / "dep_results_matrix.tsv"
    subprocess.run(
        [
            sys.executable,
            "scripts/trait_depletion.py",
            "--input",
            str(profiles_tsv),
            "--taxa-matrix",
            str(matrix_tsv),
            "--output",
            str(report_dep_matrix),
            "--tsv-output",
            str(results_dep_matrix),
            "--no-fdr",
            "--p-threshold",
            "1.0",
            "--max-fold",
            "10.0",
        ],
        check=True,
    )
    content_dep_matrix = results_dep_matrix.read_text(encoding="utf-8")
    assert "# Background-Taxa-Universe-M: 4" in content_dep_matrix
    assert "# Background-Mode: TaxaMatrix" in content_dep_matrix

    res_dep_bad = subprocess.run(
        [
            sys.executable,
            "scripts/trait_depletion.py",
            "--input",
            str(profiles_tsv),
            "--taxa-matrix",
            str(tmp_path / "nonexistent_matrix.tsv"),
            "--output",
            str(report_dep_matrix),
        ],
        capture_output=True,
        text=True,
    )
    assert res_dep_bad.returncode != 0
    assert "Specified taxa matrix file does not exist" in res_dep_bad.stderr


def test_malformed_taxa_count_strict_rejection(tmp_path: Path) -> None:
    """Verify issues #131 / #184: malformed, missing, negative, or invalid
    taxa_count values fail loudly with file/line/column context instead of
    silently substituting max trait counts.
    With background M=4, K=2 (TraitA present in 2 of 4 taxa), and profile trait count k=2:
    - Valid N=3 produces p=0.5, fold=1.3333
    - Replacing N with invalid must fail, not silently substitute N=2 (p=1/6, fold=2.0).
    """
    matrix_tsv = tmp_path / "taxa_trait_matrix.tsv"
    matrix_tsv.write_text(
        "taxon_id\tTraitA\nNCBITaxon:1\t1\nNCBITaxon:2\t1\nNCBITaxon:3\t0\nNCBITaxon:4\t0\n",
        encoding="utf-8",
    )

    valid_profiles_tsv = tmp_path / "valid_profiles.tsv"
    valid_profiles_tsv.write_text(
        "study_id\texperiment_id\tinc_dec\ttaxa_count\tTraitA\n"
        "s1\te1\tinc\t3\t2\n"
        "s1\te1\tdec\t3\t0\n",
        encoding="utf-8",
    )

    report_valid = tmp_path / "report_valid.txt"
    results_valid = tmp_path / "results_valid.tsv"

    # Valid N=3: P[X>=2] = hypergeom.sf(1, 4, 2, 3) = (comb(2,2)*comb(2,1)) / comb(4,3)
    # = (1 * 2) / 4 = 2/4 = 0.5. Expected = 3 * (2/4) = 1.5. Fold = 2 / 1.5 = 1.333333
    subprocess.run(
        [
            sys.executable,
            "scripts/trait_enrichment.py",
            "--input",
            str(valid_profiles_tsv),
            "--taxa-matrix",
            str(matrix_tsv),
            "--output",
            str(report_valid),
            "--tsv-output",
            str(results_valid),
            "--no-fdr",
            "--p-threshold",
            "1.0",
            "--min-fold",
            "0.0",
        ],
        check=True,
    )
    with results_valid.open("r", encoding="utf-8") as f:
        reader = csv.DictReader((line for line in f if not line.startswith("#")), delimiter="\t")
        rows = list(reader)
        assert abs(float(rows[0]["p_value"]) - 0.5) < 1e-5
        assert abs(float(rows[0]["fold"]) - 1.333333) < 1e-4
        assert int(rows[0]["count"]) == 2

    # Case 1: Invalid non-integer taxa_count string 'invalid'
    invalid_profiles_tsv = tmp_path / "invalid_profiles.tsv"
    invalid_profiles_tsv.write_text(
        "study_id\texperiment_id\tinc_dec\ttaxa_count\tTraitA\n"
        "s1\te1\tinc\tinvalid\t2\n"
        "s1\te1\tdec\t3\t0\n",
        encoding="utf-8",
    )
    res1 = subprocess.run(
        [
            sys.executable,
            "scripts/trait_enrichment.py",
            "--input",
            str(invalid_profiles_tsv),
            "--taxa-matrix",
            str(matrix_tsv),
            "--output",
            str(tmp_path / "out.txt"),
        ],
        capture_output=True,
        text=True,
    )
    assert res1.returncode != 0
    assert "malformed 'taxa_count' value 'invalid'" in res1.stderr
    assert "line 2" in res1.stderr

    # Case 2: Missing / empty taxa_count value
    missing_val_tsv = tmp_path / "missing_val.tsv"
    missing_val_tsv.write_text(
        "study_id\texperiment_id\tinc_dec\ttaxa_count\tTraitA\n"
        "s1\te1\tinc\t\t2\n"
        "s1\te1\tdec\t3\t0\n",
        encoding="utf-8",
    )
    res2 = subprocess.run(
        [
            sys.executable,
            "scripts/trait_enrichment.py",
            "--input",
            str(missing_val_tsv),
            "--taxa-matrix",
            str(matrix_tsv),
            "--output",
            str(tmp_path / "out.txt"),
        ],
        capture_output=True,
        text=True,
    )
    assert res2.returncode != 0
    assert "missing required 'taxa_count' column or value" in res2.stderr
    assert "line 2" in res2.stderr

    # Case 3: taxa_count < observed trait count k (taxa_count=1, TraitA=2)
    less_taxa_tsv = tmp_path / "less_taxa.tsv"
    less_taxa_tsv.write_text(
        "study_id\texperiment_id\tinc_dec\ttaxa_count\tTraitA\n"
        "s1\te1\tinc\t1\t2\n"
        "s1\te1\tdec\t3\t0\n",
        encoding="utf-8",
    )
    res3 = subprocess.run(
        [
            sys.executable,
            "scripts/trait_enrichment.py",
            "--input",
            str(less_taxa_tsv),
            "--taxa-matrix",
            str(matrix_tsv),
            "--output",
            str(tmp_path / "out.txt"),
        ],
        capture_output=True,
        text=True,
    )
    assert res3.returncode != 0
    assert "is less than observed trait count 2" in res3.stderr
    assert "line 2" in res3.stderr

    # Case 4: taxa_count > universe M (taxa_count=5, universe M=4)
    more_taxa_tsv = tmp_path / "more_taxa.tsv"
    more_taxa_tsv.write_text(
        "study_id\texperiment_id\tinc_dec\ttaxa_count\tTraitA\n"
        "s1\te1\tinc\t5\t2\n"
        "s1\te1\tdec\t3\t0\n",
        encoding="utf-8",
    )
    res4 = subprocess.run(
        [
            sys.executable,
            "scripts/trait_enrichment.py",
            "--input",
            str(more_taxa_tsv),
            "--taxa-matrix",
            str(matrix_tsv),
            "--output",
            str(tmp_path / "out.txt"),
        ],
        capture_output=True,
        text=True,
    )
    assert res4.returncode != 0
    assert "smaller than" in res4.stderr
    assert "maximum profile draw N (5)" in res4.stderr

    # Case 5: Same validations hold for trait_depletion.py
    res5 = subprocess.run(
        [
            sys.executable,
            "scripts/trait_depletion.py",
            "--input",
            str(invalid_profiles_tsv),
            "--taxa-matrix",
            str(matrix_tsv),
            "--output",
            str(tmp_path / "out.txt"),
        ],
        capture_output=True,
        text=True,
    )
    assert res5.returncode != 0
    assert "malformed 'taxa_count' value 'invalid'" in res5.stderr


def test_complete_unfiltered_export_254(tmp_path: Path) -> None:
    """Verify that --all-results-output emits a complete unfiltered per-trait table (#254).

    Checks:
    - Row count in the all-results file equals Total-Tests-Performed.
    - Zero discoveries still write all tested rows with is_significant=False.
    - Hand-computed M, N, K, k, expected, and fold match exactly.
    - Significant-hits TSV is a strict subset of the all-results TSV.
    - Depletion analyzer also produces the complete all-results TSV.
    """
    matrix_tsv = tmp_path / "taxa_trait_matrix.tsv"
    matrix_tsv.write_text(
        "taxon_id\tTraitA\tTraitB\n"
        "NCBITaxon:1\t1\t0\n"
        "NCBITaxon:2\t1\t1\n"
        "NCBITaxon:3\t0\t1\n"
        "NCBITaxon:4\t0\t0\n",
        encoding="utf-8",
    )

    profiles_tsv = tmp_path / "profiles.tsv"
    profiles_tsv.write_text(
        "study_id\texperiment_id\tinc_dec\ttaxa_count\tTraitA\tTraitB\n"
        "s1\te1\tinc\t3\t2\t0\n"
        "s1\te1\tdec\t3\t0\t1\n",
        encoding="utf-8",
    )

    report_txt = tmp_path / "report.txt"
    tsv_out = tmp_path / "results.tsv"
    all_out = tmp_path / "all_results.tsv"

    # Step 1: Run enrichment with strict thresholds (0 discoveries)
    subprocess.run(
        [
            sys.executable,
            "scripts/trait_enrichment.py",
            "--input",
            str(profiles_tsv),
            "--taxa-matrix",
            str(matrix_tsv),
            "--output",
            str(report_txt),
            "--tsv-output",
            str(tsv_out),
            "--all-results-output",
            str(all_out),
            "--p-threshold",
            "0.001",
            "--min-fold",
            "1.5",
        ],
        check=True,
    )

    assert all_out.exists()
    all_lines = all_out.read_text(encoding="utf-8").strip().splitlines()

    # Verify metadata header
    meta = {}
    data_rows = []
    for line in all_lines:
        if line.startswith("#"):
            parts = line[1:].strip().split(":", 1)
            if len(parts) == 2:
                meta[parts[0].strip()] = parts[1].strip()
        else:
            data_rows.append(line.split("\t"))

    assert meta["Total-Tests-Performed"] == "2"
    assert meta["Significant-Hits"] == "0"
    assert meta["Result-Status"] == "zero_discoveries"

    header = data_rows[0]
    expected_cols = [
        "study_id",
        "experiment_id",
        "inc_dec",
        "trait",
        "p_value",
        "fdr_q",
        "fold",
        "count",
        "expected",
        "testable_profiles",
        "background_m",
        "background_k",
        "is_significant",
        "correction_family",
        "eligibility_status",
    ]
    assert header == expected_cols

    # 2 tests performed, so 2 data rows
    rows = data_rows[1:]
    assert len(rows) == 2

    # Parse rows into dicts
    row_dicts = [dict(zip(header, r, strict=True)) for r in rows]

    # Verify hand-computed values:
    # TraitA in s1|e1|inc: k=2, N=3, M=4, K=2 (from matrix: NCBITaxon:1, 2)
    # Expected = 3 * (2/4) = 1.5. Fold = 2 / 1.5 = 1.333333
    row_a = next(r for r in row_dicts if r["trait"] == "TraitA" and r["inc_dec"] == "inc")
    assert row_a["study_id"] == "s1"
    assert row_a["experiment_id"] == "e1"
    assert row_a["count"] == "2"
    assert row_a["background_m"] == "4"
    assert row_a["background_k"] == "2"
    assert float(row_a["expected"]) == pytest.approx(1.5, rel=1e-3)
    assert float(row_a["fold"]) == pytest.approx(1.3333, rel=1e-3)
    assert row_a["is_significant"] == "False"
    assert row_a["eligibility_status"] == "testable"

    # In strict mode, tsv_out has zero discovery data rows
    tsv_lines = tsv_out.read_text(encoding="utf-8").strip().splitlines()
    tsv_data = [line for line in tsv_lines if not line.startswith("#")]
    assert len(tsv_data) == 1  # Header only

    # Step 2: Run with loose thresholds to produce discoveries and verify subset
    tsv_loose = tmp_path / "results_loose.tsv"
    all_loose = tmp_path / "all_results_loose.tsv"
    subprocess.run(
        [
            sys.executable,
            "scripts/trait_enrichment.py",
            "--input",
            str(profiles_tsv),
            "--taxa-matrix",
            str(matrix_tsv),
            "--output",
            str(report_txt),
            "--tsv-output",
            str(tsv_loose),
            "--all-results-output",
            str(all_loose),
            "--p-threshold",
            "1.0",
            "--min-fold",
            "1.0",
        ],
        check=True,
    )

    all_loose_lines = all_loose.read_text(encoding="utf-8").strip().splitlines()
    all_loose_data = [line.split("\t") for line in all_loose_lines if not line.startswith("#")]
    all_loose_dicts = [dict(zip(all_loose_data[0], r, strict=True)) for r in all_loose_data[1:]]

    tsv_loose_lines = tsv_loose.read_text(encoding="utf-8").strip().splitlines()
    tsv_loose_data = [line.split("\t") for line in tsv_loose_lines if not line.startswith("#")]
    tsv_loose_dicts = [dict(zip(tsv_loose_data[0], r, strict=True)) for r in tsv_loose_data[1:]]

    # Every significant row must exist in all_results with is_significant=True
    assert len(tsv_loose_dicts) > 0
    for sig_row in tsv_loose_dicts:
        matching = [
            r
            for r in all_loose_dicts
            if r["study_id"] == sig_row["study_id"]
            and r["experiment_id"] == sig_row["experiment_id"]
            and r["inc_dec"] == sig_row["inc_dec"]
            and r["trait"] == sig_row["trait"]
        ]
        assert len(matching) == 1
        assert matching[0]["is_significant"] == "True"

    # Step 3: Run depletion analysis with --all-results-output
    dep_tsv = tmp_path / "dep_results.tsv"
    dep_all = tmp_path / "dep_all.tsv"
    subprocess.run(
        [
            sys.executable,
            "scripts/trait_depletion.py",
            "--input",
            str(profiles_tsv),
            "--taxa-matrix",
            str(matrix_tsv),
            "--output",
            str(report_txt),
            "--tsv-output",
            str(dep_tsv),
            "--all-results-output",
            str(dep_all),
            "--p-threshold",
            "1.0",
            "--max-fold",
            "10.0",
        ],
        check=True,
    )

    assert dep_all.exists()
    dep_lines = dep_all.read_text(encoding="utf-8").strip().splitlines()
    dep_data = [line.split("\t") for line in dep_lines if not line.startswith("#")]
    assert dep_data[0] == expected_cols
    assert len(dep_data[1:]) > 0


def test_enrichment_and_depletion_zero_eligibility_and_aligned_denominators(tmp_path: Path) -> None:
    """Validate zero-annotation eligibility and aligned testability denominators (#132).

    Four-profile fixture across two paired experiments:
    s1/e1/inc: TraitA=10, TraitB=0, TraitC=0, taxa_count=10
    s1/e1/dec: TraitA=0, TraitB=10, TraitC=0, taxa_count=10
    s2/e1/inc: TraitA=10, TraitB=0, TraitC=0, taxa_count=10
    s2/e1/dec: TraitA=0, TraitB=10, TraitC=0, taxa_count=10

    TraitA is observed in inc arms; its 0 in dec arms is an observed zero.
    TraitB is observed in dec arms; its 0 in inc arms is an observed zero.
    TraitC is 0 in all arms (structural zero / coverage state): never tested.

    Both enrichment and depletion must report:
    - TraitA testable_profiles = 4
    - TraitB testable_profiles = 4
    - TraitC not tested (0 tests, 0 hits)
    - TraitA enriched in inc (2 hits), depleted in dec (2 hits)
    - TraitB enriched in dec (2 hits), depleted in inc (2 hits)
    """
    profiles_tsv = tmp_path / "profiles.tsv"
    profiles_tsv.write_text(
        "study_id\texperiment_id\tinc_dec\ttaxa_count\tTraitA\tTraitB\tTraitC\n"
        "s1\te1\tinc\t10\t10\t0\t0\n"
        "s1\te1\tdec\t10\t0\t10\t0\n"
        "s2\te1\tinc\t10\t10\t0\t0\n"
        "s2\te1\tdec\t10\t0\t10\t0\n",
        encoding="utf-8",
    )

    enr_report = tmp_path / "enr_report.txt"
    enr_tsv = tmp_path / "enr_results.tsv"
    res_enr = subprocess.run(
        [
            sys.executable,
            "scripts/trait_enrichment.py",
            "--input",
            str(profiles_tsv),
            "--output",
            str(enr_report),
            "--tsv-output",
            str(enr_tsv),
            "--p-threshold",
            "0.05",
            "--min-fold",
            "1.5",
        ],
        capture_output=True,
        text=True,
    )
    assert res_enr.returncode == 0, res_enr.stderr

    dep_report = tmp_path / "dep_report.txt"
    dep_tsv = tmp_path / "dep_results.tsv"
    res_dep = subprocess.run(
        [
            sys.executable,
            "scripts/trait_depletion.py",
            "--input",
            str(profiles_tsv),
            "--output",
            str(dep_report),
            "--tsv-output",
            str(dep_tsv),
            "--p-threshold",
            "0.05",
            "--max-fold",
            "0.67",
        ],
        capture_output=True,
        text=True,
    )
    assert res_dep.returncode == 0, res_dep.stderr

    # Parse TSV outputs
    def _read(tsv: Path) -> tuple[dict[str, str], list[str], list[dict[str, str]]]:
        lines = tsv.read_text(encoding="utf-8").splitlines()
        meta: dict[str, str] = {}
        body: list[str] = []
        for line in lines:
            if line.startswith("#"):
                k, _, v = line.lstrip("#").strip().partition(": ")
                meta[k] = v
            else:
                body.append(line)
        header = body[0].split("\t")
        rows = [dict(zip(header, row.split("\t"), strict=True)) for row in body[1:] if row]
        return meta, header, rows

    enr_meta, _, enr_rows = _read(enr_tsv)
    dep_meta, _, dep_rows = _read(dep_tsv)

    # 1. Total tests and metadata
    assert enr_meta["Total-Tests-Performed"] == "4"
    assert enr_meta["Significant-Hits"] == "4"
    assert enr_meta["Result-Status"] == "discoveries"

    # Depletion tests 4 profile-trait pairs for TraitA and 4 for TraitB;
    # TraitC is skipped (structural zero)
    assert dep_meta["Total-Tests-Performed"] == "8"
    assert dep_meta["Significant-Hits"] == "4"
    assert dep_meta["Result-Status"] == "discoveries"

    # 2. Denominators (testable_profiles) are aligned: 4 for TraitA, 4 for TraitB
    for r in enr_rows:
        assert r["testable_profiles"] == "4"
        assert r["trait"] in {"TraitA", "TraitB"}

    for r in dep_rows:
        assert r["testable_profiles"] == "4"
        assert r["trait"] in {"TraitA", "TraitB"}

    # 3. Direction of discoveries
    # TraitA: enriched in inc, depleted in dec (observed zero with count 0)
    enr_a = [r for r in enr_rows if r["trait"] == "TraitA"]
    dep_a = [r for r in dep_rows if r["trait"] == "TraitA"]
    assert len(enr_a) == 2
    assert all(r["inc_dec"] == "inc" for r in enr_a)
    assert len(dep_a) == 2
    assert all(r["inc_dec"] == "dec" and r["count"] == "0" for r in dep_a)

    # TraitB: enriched in dec, depleted in inc (observed zero with count 0)
    enr_b = [r for r in enr_rows if r["trait"] == "TraitB"]
    dep_b = [r for r in dep_rows if r["trait"] == "TraitB"]
    assert len(enr_b) == 2
    assert all(r["inc_dec"] == "dec" for r in enr_b)
    assert len(dep_b) == 2
    assert all(r["inc_dec"] == "inc" and r["count"] == "0" for r in dep_b)

    # TraitC: never tested, never reported
    assert not any(r["trait"] == "TraitC" for r in enr_rows)
    assert not any(r["trait"] == "TraitC" for r in dep_rows)
