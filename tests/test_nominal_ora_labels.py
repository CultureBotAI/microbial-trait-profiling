"""Nominal labels qualify existing ORA values and calls without changing them (#1503)."""

from __future__ import annotations

from pathlib import Path
from types import ModuleType

import pytest
from scipy.stats import false_discovery_control, hypergeom

from test_nominal_ora_pipeline import _table
from test_nominal_ora_pipeline import inputs as inputs
from test_nominal_ora_pipeline import producers as producers


@pytest.mark.parametrize("tail", ["upper", "lower"])
@pytest.mark.parametrize("mode", ["none", "global", "profile"])
def test_producer_and_rankings_label_the_actual_threshold_and_hits(
    tmp_path: Path,
    inputs: tuple[Path, Path],
    producers: tuple[ModuleType, ModuleType],
    capsys: pytest.CaptureFixture[str],
    tail: str,
    mode: str,
) -> None:
    analyzer, rankings = producers
    profiles, matrix = inputs
    report, selected, complete = (tmp_path / name for name in ("report.txt", "hits.tsv", "all.tsv"))
    threshold = 0.25 if tail == "upper" else 0.5
    analyzer.main(
        [
            "--input",
            str(profiles),
            "--taxa-matrix",
            str(matrix),
            "--output",
            str(report),
            "--tsv-output",
            str(selected),
            "--all-results-output",
            str(complete),
            "--tail",
            tail,
            "--p-threshold",
            str(threshold),
            "--min-fold",
            "1.2",
            "--max-fold",
            "0.8",
            "--fdr-scope",
            mode if mode != "none" else "global",
            "--no-fdr" if mode == "none" else "--fdr",
        ]
    )
    stdout = capsys.readouterr().out
    metadata, rows = _table(complete)
    _, selected_rows = _table(selected)
    expected_p = []
    for row in rows:
        background = 3 if row["trait"] == "TraitA" else 5
        count = int(row["count"])
        expected_p.append(
            float(
                hypergeom.sf(count - 1, 10, background, 2)
                if tail == "upper"
                else hypergeom.cdf(count, 10, background, 2)
            )
        )
    assert [float(row["p_value"]) for row in rows] == pytest.approx(expected_p)
    expected_q = list(expected_p)
    if mode != "none":
        groups = (
            [list(range(len(rows)))]
            if mode == "global"
            else [
                [i for i, row in enumerate(rows) if row["inc_dec"] == direction]
                for direction in ("inc", "dec")
            ]
        )
        for group in groups:
            values = false_discovery_control([expected_p[i] for i in group], method="bh")
            for i, value in zip(group, values, strict=True):
                expected_q[i] = float(value)
        assert [float(row["fdr_q"]) for row in rows] == pytest.approx(expected_q)
    else:
        assert all(row["fdr_q"] == "" for row in rows)
    decisions = [
        probability <= threshold
        and (float(row["fold"]) >= 1.2 if tail == "upper" else float(row["fold"]) <= 0.8)
        for row, probability in zip(rows, expected_q, strict=True)
    ]
    assert [row["is_significant"] == "True" for row in rows] == decisions
    assert len(selected_rows) == sum(decisions) == int(metadata["Significant-Hits"])
    assert float(metadata["P-Threshold"]) == threshold
    assert metadata["FDR-Scope"] == mode
    label = "nominal p" if mode == "none" else "FDR q"
    hits_label = "Nominal hits" if mode == "none" else "Significant hits"
    for path in (selected, complete):
        header, _ = _table(path)
        if mode == "none":
            assert header["Threshold-Statistic"] == label
            assert header["Hit-Statistic"] == label
        else:
            assert "Threshold-Statistic" not in header
            assert "Hit-Statistic" not in header
    text = report.read_text(encoding="utf-8")
    assert text.startswith(f"{label} <= {threshold:g}")
    assert f"{hits_label}: {sum(decisions)}" in text
    assert f"{hits_label}: {sum(decisions)}" in stdout
    if mode == "none":
        assert "Significant hits:" not in text + stdout
    assert len(
        rankings.parse_report(
            report, stat_type="p_value" if mode == "none" else "fdr_q", allow_empty=True
        )
    ) == sum(decisions)

    overall, pairs, tested = (
        tmp_path / name for name in ("overall.tsv", "pairs.tsv", "tested.tsv")
    )
    assert (
        rankings.main(
            [
                "--input",
                str(report),
                "--results",
                str(complete),
                "--unfiltered",
                "--profiles",
                str(profiles),
                "--overall-output",
                str(overall),
                "--pair-output",
                str(pairs),
                "--tested-output",
                str(tested),
            ]
        )
        == 0
    )
    exported_stdout = capsys.readouterr().out
    expected_stdout = "Analyzer nominal hits" if mode == "none" else "Analyzer significant hits"
    assert f"{expected_stdout}: {sum(decisions)}" in exported_stdout
    if mode == "none":
        assert "significant hits" not in exported_stdout
    for path in (overall, pairs, tested):
        header, _ = _table(path)
        if mode == "none":
            assert header["Threshold-Statistic"] == label
            assert header["Hit-Statistic"] == label
        else:
            assert "Threshold-Statistic" not in header
            assert "Hit-Statistic" not in header
        assert float(header["Stat-Threshold"]) == threshold
        assert header["FDR-Scope"] == mode
    filtered = rankings.parse_results_tsv(
        selected,
        stat_type="p_value" if mode == "none" else "fdr_q",
        threshold=threshold,
        min_fold=1.2,
        max_fold=0.8,
        tail=tail,
    )
    filtered_stdout = capsys.readouterr().out
    assert len(filtered) == sum(decisions)
    filtered_label = "nominal hits" if mode == "none" else "significant hits"
    assert f"Filtered to {sum(decisions)} {filtered_label}" in filtered_stdout
    if mode == "none":
        assert "significant hits" not in filtered_stdout


@pytest.mark.parametrize("prefix", ["p", "nominal p"])
@pytest.mark.parametrize("operator", ["<", "<="])
@pytest.mark.parametrize("tail", ["upper", "lower"])
def test_report_metadata_accepts_legacy_and_explicit_nominal_labels(
    tmp_path: Path, producers: tuple[ModuleType, ModuleType], prefix: str, operator: str, tail: str
) -> None:
    _, rankings = producers
    report = tmp_path / "report.txt"
    suffix = " (Depletion)" if tail == "lower" else ""
    hit_label = "Nominal hits" if prefix == "nominal p" else "Significant hits"
    report.write_text(
        f"{prefix} {operator} 0.037{suffix}\n"
        f"Tests performed: 19\n{hit_label}: 3\nResult status: discoveries\n",
        encoding="utf-8",
    )
    metadata = rankings.detect_metadata(report)
    assert metadata["stat_type"] == "p_value"
    assert metadata["threshold"] == 0.037
    assert metadata["tail"] == tail
    assert metadata["tests_performed"] == 19
    assert metadata["significant_hits"] == 3
