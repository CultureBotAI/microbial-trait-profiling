"""Regression tests for #296: zero-discovery handoff from analyzer to exporter to summary.

The analyzer (scripts/trait_enrichment.py), the ranking exporter
(scripts/export_trait_enrichment_rankings.py) and the enrichment summary
(scripts/generate_enrichment_summary.py) must:

* finish successfully on a valid result with zero discoveries, preserving the
  declared output headers and recording status/counts;
* distinguish zero discoveries from zero eligible tests;
* fail with an actionable message (not an empty successful report) on
  malformed input.

Analysis thresholds are the script defaults; none are changed here.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from bsdbng.enrichment_results import (
    RESULT_COLUMNS,
    STATUS_DISCOVERIES,
    STATUS_ZERO_DISCOVERIES,
    STATUS_ZERO_ELIGIBLE_TESTS,
    MalformedResultsError,
    classify_result_status,
    load_results_tsv,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
ANALYZER = REPO_ROOT / "scripts" / "trait_enrichment.py"
EXPORTER = REPO_ROOT / "scripts" / "export_trait_enrichment_rankings.py"
SUMMARY = REPO_ROOT / "scripts" / "generate_enrichment_summary.py"

# Fixture from issue #296: one paired experiment, two traits -> four tests, zero hits.
ISSUE_296_PROFILES = (
    "study_id\texperiment_id\tinc_dec\ttaxa_count\tTraitA\tTraitB\n"
    "s1\te1\tinc\t2\t1\t1\n"
    "s1\te1\tdec\t2\t1\t1\n"
)

EXPECTED_OVERALL_HEADER = [
    "rank",
    "trait",
    "net_direction",
    "inc_hits",
    "dec_hits",
    "prevalence_summary",
    "stat_type",
    "best_qvalue",
    "median_qvalue",
    "hit_rows",
    "unique_profiles",
    "mean_fold",
    "max_fold",
    "testable_profiles",
    "best_pvalue",
    "min_p_direction",
    "min_p_study_id",
    "min_p_experiment_id",
    "min_p_profile_id",
    "prevalence",
]
EXPECTED_PAIR_HEADER = [
    "rank",
    "study_id",
    "experiment_id",
    "inc_dec",
    "profile_id",
    "trait",
    "p_value",
    "q_value",
    "fold",
    "obs",
    "exp",
    "testable_profiles",
    "hit_count",
    "is_censored_hit",
]


def _run(script: Path, *args: str | Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(script), *(str(a) for a in args)],
        capture_output=True,
        text=True,
        check=False,
        cwd=REPO_ROOT,
    )


def _read_tsv(path: Path) -> tuple[dict[str, str], list[str], list[list[str]]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    meta: dict[str, str] = {}
    body: list[str] = []
    for line in lines:
        if line.startswith("#"):
            key, _, value = line.lstrip("#").strip().partition(": ")
            meta[key] = value
        else:
            body.append(line)
    header = body[0].split("\t")
    rows = [line.split("\t") for line in body[1:] if line]
    return meta, header, rows


def _analyze(tmp_path: Path, profiles: str, *extra: str) -> tuple[Path, Path]:
    profiles_tsv = tmp_path / "profiles.tsv"
    profiles_tsv.write_text(profiles, encoding="utf-8")
    report = tmp_path / "report.txt"
    results = tmp_path / "results.tsv"
    proc = _run(
        ANALYZER,
        "--input",
        profiles_tsv,
        "--output",
        report,
        "--tsv-output",
        results,
        *extra,
    )
    assert proc.returncode == 0, proc.stderr
    return report, results


def _export(tmp_path: Path, report: Path, results: Path) -> subprocess.CompletedProcess[str]:
    return _run(
        EXPORTER,
        "--input",
        report,
        "--results",
        results,
        "--overall-output",
        tmp_path / "ranked.tsv",
        "--pair-output",
        tmp_path / "pairs.tsv",
    )


def _summarize(tmp_path: Path, results: Path) -> subprocess.CompletedProcess[str]:
    return _run(
        SUMMARY,
        "--results",
        results,
        "--experiments",
        tmp_path / "absent_full_dump.csv",
        "--mapping",
        tmp_path / "absent_mapping.tsv",
        "--trait-meta",
        tmp_path / "absent_columns.tsv",
        "--output-text",
        tmp_path / "summary.txt",
        "--output-image",
        tmp_path / "summary.png",
    )


def test_issue_296_zero_discovery_handoff(tmp_path: Path) -> None:
    """Analyzer -> exporter -> summary succeeds on the #296 fixture at default thresholds."""
    report, results = _analyze(tmp_path, ISSUE_296_PROFILES)

    # Analyzer: four tests, zero selected rows, declared columns, explicit status.
    meta, header, rows = _read_tsv(results)
    assert header == list(RESULT_COLUMNS)
    assert rows == []
    assert meta["Total-Tests-Performed"] == "4"
    assert meta["Significant-Hits"] == "0"
    assert meta["Result-Status"] == STATUS_ZERO_DISCOVERIES
    assert meta["P-Threshold"] == "0.001"
    assert meta["Min-Fold"] == "1.5"
    assert "Result status: zero_discoveries" in report.read_text(encoding="utf-8")

    # Exporter: succeeds and writes both declared schemas with zero rows.
    proc = _export(tmp_path, report, results)
    assert proc.returncode == 0, proc.stderr
    assert "KeyError" not in proc.stderr

    overall_meta, overall_header, overall_rows = _read_tsv(tmp_path / "ranked.tsv")
    assert overall_header == EXPECTED_OVERALL_HEADER
    assert overall_rows == []
    assert overall_meta["Result-Status"] == STATUS_ZERO_DISCOVERIES
    assert overall_meta["Total-Tests-Performed"] == "4"
    assert overall_meta["Significant-Hits"] == "0"
    assert overall_meta["Ranked-Traits"] == "0"
    assert overall_meta["Stat-Threshold"] == "0.001"
    assert overall_meta["Min-Fold"] == "1.5"

    pair_meta, pair_header, pair_rows = _read_tsv(tmp_path / "pairs.tsv")
    assert pair_header == EXPECTED_PAIR_HEADER
    assert pair_rows == []
    assert pair_meta["Result-Status"] == STATUS_ZERO_DISCOVERIES
    assert pair_meta["Total-Tests-Performed"] == "4"
    assert pair_meta["Exported-Rows"] == "0"

    # Downstream summary accepts the valid empty result and records its status.
    proc = _summarize(tmp_path, results)
    assert proc.returncode == 0, proc.stderr
    summary = (tmp_path / "summary.txt").read_text(encoding="utf-8")
    assert "Result-Status: zero_discoveries" in summary
    assert "Total-Tests-Performed: 4" in summary
    assert "Selected-Rows: 0" in summary
    assert "Directionality Bias: none" in summary
    assert (tmp_path / "summary.png").stat().st_size > 0


def test_zero_eligible_tests_is_distinct_from_zero_discoveries(tmp_path: Path) -> None:
    """Profiles without an inc/dec pair yield zero tests: status zero_eligible_tests."""
    unpaired = (
        "study_id\texperiment_id\tinc_dec\ttaxa_count\tTraitA\tTraitB\n"
        "s1\te1\tinc\t2\t1\t1\n"
        "s2\te1\tdec\t2\t1\t1\n"
    )
    report, results = _analyze(tmp_path, unpaired)

    meta, header, rows = _read_tsv(results)
    assert header == list(RESULT_COLUMNS)
    assert rows == []
    assert meta["Total-Tests-Performed"] == "0"
    assert meta["Result-Status"] == STATUS_ZERO_ELIGIBLE_TESTS

    proc = _export(tmp_path, report, results)
    assert proc.returncode == 0, proc.stderr
    overall_meta, overall_header, _ = _read_tsv(tmp_path / "ranked.tsv")
    assert overall_header == EXPECTED_OVERALL_HEADER
    assert overall_meta["Result-Status"] == STATUS_ZERO_ELIGIBLE_TESTS
    assert overall_meta["Total-Tests-Performed"] == "0"

    proc = _summarize(tmp_path, results)
    assert proc.returncode == 0, proc.stderr
    assert "Result-Status: zero_eligible_tests" in (tmp_path / "summary.txt").read_text(
        encoding="utf-8"
    )


def test_discoveries_path_still_ranks(tmp_path: Path) -> None:
    """A result with selected rows keeps status discoveries and the same declared headers."""
    profiles = (
        "study_id\texperiment_id\tinc_dec\ttaxa_count\tTraitA\tTraitB\n"
        "s1\te1\tinc\t10\t10\t1\n"
        "s1\te1\tdec\t10\t1\t10\n"
        "s2\te1\tinc\t10\t10\t0\n"
        "s2\te1\tdec\t10\t0\t10\n"
    )
    # Permissive threshold only to guarantee rows in a 4-profile fixture.
    report, results = _analyze(tmp_path, profiles, "--p-threshold", "1.0")
    meta, _, rows = _read_tsv(results)
    assert meta["Result-Status"] == STATUS_DISCOVERIES
    assert int(meta["Significant-Hits"]) == len(rows) > 0

    proc = _export(tmp_path, report, results)
    assert proc.returncode == 0, proc.stderr
    overall_meta, overall_header, overall_rows = _read_tsv(tmp_path / "ranked.tsv")
    assert overall_header == EXPECTED_OVERALL_HEADER
    assert overall_meta["Result-Status"] == STATUS_DISCOVERIES
    assert int(overall_meta["Ranked-Traits"]) == len(overall_rows) > 0
    pair_meta, pair_header, pair_rows = _read_tsv(tmp_path / "pairs.tsv")
    assert pair_header == EXPECTED_PAIR_HEADER
    assert int(pair_meta["Exported-Rows"]) == len(pair_rows) == len(rows)


_HEADER = "\t".join(RESULT_COLUMNS)
_META = "# FDR: True\n# P-Threshold: 0.001\n# Min-Fold: 1.5\n"


@pytest.mark.parametrize(
    ("content", "message"),
    [
        pytest.param(f"{_META}{_HEADER}\n", "no '# Total-Tests-Performed' metadata", id="no-count"),
        pytest.param(
            f"{_META}# Total-Tests-Performed: 4\nstudy_id\ttrait\n",
            "missing required columns",
            id="missing-columns",
        ),
        pytest.param(
            f"{_META}# Total-Tests-Performed: 4\n# Significant-Hits: 2\n{_HEADER}\n"
            "s1\te1\tinc\tTraitA\t0.0001\t0.0001\t2.0\t2\t1.0\t2\n",
            "declares 2 significant hits but the file contains 1 rows",
            id="truncated",
        ),
        pytest.param(
            f"{_META}# Total-Tests-Performed: 4\n{_HEADER}\n"
            "s1\te1\tinc\tTraitA\tnot-a-number\t0.0001\t2.0\t2\t1.0\t2\n",
            "invalid value 'not-a-number' in column 'p_value'",
            id="non-numeric",
        ),
        pytest.param(
            f"{_META}# Total-Tests-Performed: 0\n{_HEADER}\n"
            "s1\te1\tinc\tTraitA\t0.0001\t0.0001\t2.0\t2\t1.0\t2\n",
            "exceed the 0 tests performed",
            id="rows-exceed-tests",
        ),
        pytest.param("", "no column header line found", id="empty-file"),
    ],
)
def test_malformed_results_fail_actionably(tmp_path: Path, content: str, message: str) -> None:
    """Malformed results make both exporter and summary exit 1 with a specific message."""
    results = tmp_path / "results.tsv"
    results.write_text(content, encoding="utf-8")

    with pytest.raises(MalformedResultsError, match=r"[Rr]egenerate|exceed|missing"):
        load_results_tsv(results)

    proc = _export(tmp_path, tmp_path / "absent_report.txt", results)
    assert proc.returncode == 1
    assert message in proc.stderr
    assert "Traceback" not in proc.stderr
    assert not (tmp_path / "ranked.tsv").exists()
    assert not (tmp_path / "pairs.tsv").exists()

    proc = _summarize(tmp_path, results)
    assert proc.returncode == 1
    assert message in proc.stderr
    assert not (tmp_path / "summary.txt").exists()


def test_missing_results_fail_actionably(tmp_path: Path) -> None:
    """Missing analyzer output is an error, not an empty successful report."""
    proc = _export(tmp_path, tmp_path / "absent_report.txt", tmp_path / "absent.tsv")
    assert proc.returncode == 1
    assert "Run `scripts/trait_enrichment.py`" in proc.stderr

    proc = _summarize(tmp_path, tmp_path / "absent.tsv")
    assert proc.returncode == 1
    assert "Enrichment results TSV not found" in proc.stderr


def test_analyzer_rejects_profile_file_without_rows(tmp_path: Path) -> None:
    """A header-only profile file is malformed input, not a zero-eligible result."""
    profiles = tmp_path / "profiles.tsv"
    profiles.write_text("study_id\texperiment_id\tinc_dec\ttaxa_count\tTraitA\n", encoding="utf-8")
    proc = _run(ANALYZER, "--input", profiles, "--output", tmp_path / "report.txt")
    assert proc.returncode != 0
    assert "has a header but no profile rows" in proc.stderr


@pytest.mark.parametrize(
    ("tests", "rows", "status"),
    [
        (0, 0, STATUS_ZERO_ELIGIBLE_TESTS),
        (4, 0, STATUS_ZERO_DISCOVERIES),
        (4, 1, STATUS_DISCOVERIES),
    ],
)
def test_classify_result_status(tests: int, rows: int, status: str) -> None:
    assert classify_result_status(tests, rows) == status


def test_classify_result_status_rejects_impossible_counts() -> None:
    with pytest.raises(MalformedResultsError):
        classify_result_status(1, 2)
    with pytest.raises(MalformedResultsError):
        classify_result_status(-1, 0)
