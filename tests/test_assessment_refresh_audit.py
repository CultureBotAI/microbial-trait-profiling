"""Assessment refreshes preserve scientific values and guild identity (#729/#734)."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest


def fixture_outputs(root: Path, *, refreshed: bool) -> tuple[Path, Path]:
    calibration = root / "calibration"
    integration = root / "integrated"
    calibration.mkdir(parents=True)
    integration.mkdir()
    for name in (
        "calibration_results.tsv",
        "calibration_replicates.tsv",
        "calibration_criteria.tsv",
    ):
        (calibration / name).write_text(
            "scenario\tprocedure\tlabel\trole\treplicate\tcriterion\treplicates\t"
            "replicates_tested\tp_value\ttruth\trequirement\ttolerance\t"
            "replicates_required\texpected_verdict\n"
            "s\tfraction_sign_test\ttrait\tnull\t0\tfamily_error\t1400\t1400\t0.25\t"
            "null\tfixed requirement\t0.065\t1400\t\n"
            "s\tguild_mcnemar\tguild\tnull\t0\tfamily_error\t1400\t1400\t"
            f"{'0.5' if refreshed else '0.1'}\tnull\tfixed requirement\t0.065\t1400\t\n",
            encoding="utf-8",
        )
    for name in ("calibration_flip_assignments.tsv", "calibration_injection.tsv"):
        (calibration / name).write_text("same retained bytes\n", encoding="utf-8")
    for name in ("integrated_evidence.tsv", "integrated_evidence_method_results.tsv"):
        (integration / name).write_text(
            "method\tp_value\tcalibration_status\tsensitivity_status\tcoverage_status\t"
            "orientation_status\ttaxon_influence_status\n"
            f"fraction\t0.001\t{'pass' if refreshed else 'insufficient_budget'}\tpending\t"
            "pending\tno verdict\tno verdict\n",
            encoding="utf-8",
        )
    for name in ("integrated_evidence_features.tsv", "integrated_evidence_families.tsv"):
        (integration / name).write_text("same retained bytes\n", encoding="utf-8")
    return calibration, integration


def run_audit(
    before: tuple[Path, Path], after: tuple[Path, Path], output: Path
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            sys.executable,
            "scripts/adhoc/issue729_audit_assessment_refresh.py",
            "--before-calibration",
            str(before[0]),
            "--after-calibration",
            str(after[0]),
            "--before-integration",
            str(before[1]),
            "--after-integration",
            str(after[1]),
            "--output",
            str(output),
        ],
        check=False,
        capture_output=True,
        text=True,
    )


@pytest.mark.parametrize(
    "corruption",
    [
        None,
        "scientific_p",
        "other_assessment",
        "non_guild",
        "flip",
        "lost_row",
    ],
)
def test_audit_allows_assessment_updates_but_rejects_scientific_change(
    tmp_path: Path,
    corruption: str | None,
) -> None:
    before = fixture_outputs(tmp_path / "before", refreshed=False)
    after = fixture_outputs(tmp_path / "after", refreshed=True)
    if corruption == "scientific_p":
        path = after[1] / "integrated_evidence_method_results.tsv"
        path.write_text(path.read_text().replace("0.001", "0.002"), encoding="utf-8")
    elif corruption == "other_assessment":
        path = after[1] / "integrated_evidence_method_results.tsv"
        path.write_text(path.read_text().replace("no verdict", "pass"), encoding="utf-8")
    elif corruption == "non_guild":
        path = after[0] / "calibration_results.tsv"
        path.write_text(path.read_text().replace("0.25", "0.75"), encoding="utf-8")
    elif corruption == "flip":
        (after[0] / "calibration_flip_assignments.tsv").write_text("changed\n", encoding="utf-8")
    elif corruption == "lost_row":
        path = before[1] / "integrated_evidence.tsv"
        text = path.read_text()
        path.write_text(text + text.splitlines(keepends=True)[1], encoding="utf-8")
    output = tmp_path / "audit.json"
    result = run_audit(before, after, output)
    assert result.returncode == (0 if corruption is None else 1), result.stderr
    audit = json.loads(output.read_text())
    assert audit["status"] == ("ok" if corruption is None else "failed")
    if corruption is None:
        assert audit["table_comparisons"][-1]["changed_columns"] == {"calibration_status": 1}


@pytest.mark.parametrize("corruption", ["lost", "duplicated", "remapped"])
@pytest.mark.parametrize(
    "table",
    [
        "calibration_results.tsv",
        "calibration_replicates.tsv",
        "calibration_criteria.tsv",
    ],
)
def test_guild_statistical_changes_must_not_hide_lost_or_remapped_results(
    tmp_path: Path,
    table: str,
    corruption: str,
) -> None:
    before = fixture_outputs(tmp_path / "before", refreshed=False)
    after = fixture_outputs(tmp_path / "after", refreshed=True)
    path = after[0] / table
    lines = path.read_text().splitlines(keepends=True)
    if corruption == "lost":
        lines = lines[:-1]
    elif corruption == "duplicated":
        lines.append(lines[-1])
    elif table == "calibration_criteria.tsv":
        lines[-1] = lines[-1].replace("family_error", "different_criterion")
    else:
        lines[-1] = lines[-1].replace("s\tguild_mcnemar", "different\tguild_mcnemar")
    path.write_text("".join(lines), encoding="utf-8")
    output = tmp_path / "audit.json"
    assert run_audit(before, after, output).returncode == 1
    result = json.loads(output.read_text())
    row = next(row for row in result["table_comparisons"] if Path(row["after"]).name == table)
    assert row["guild_row_identity_equal"] is False
    assert row["guild_rows_before"] == 1
    assert row["guild_rows_after"] == {"lost": 0, "duplicated": 2, "remapped": 1}[corruption]


@pytest.mark.parametrize(
    ("table", "old", "new"),
    [
        ("calibration_results.tsv", "\t1400\t1400\t", "\t400\t1400\t"),
        ("calibration_criteria.tsv", "fixed requirement", "relaxed requirement"),
        ("calibration_criteria.tsv", "\t0.065\t", "\t0.1\t"),
    ],
)
def test_guild_policy_cannot_change_with_its_numerical_outcomes(
    tmp_path: Path,
    table: str,
    old: str,
    new: str,
) -> None:
    before = fixture_outputs(tmp_path / "before", refreshed=False)
    after = fixture_outputs(tmp_path / "after", refreshed=True)
    path = after[0] / table
    lines = path.read_text().splitlines(keepends=True)
    lines[-1] = lines[-1].replace(old, new)
    path.write_text("".join(lines), encoding="utf-8")
    output = tmp_path / "audit.json"
    assert run_audit(before, after, output).returncode == 1
    result = json.loads(output.read_text())
    assert any(row.get("guild_fixed_policy_equal") is False for row in result["table_comparisons"])


def test_audit_refuses_to_overwrite_retained_evidence(tmp_path: Path) -> None:
    before = fixture_outputs(tmp_path / "before", refreshed=False)
    after = fixture_outputs(tmp_path / "after", refreshed=True)
    target = before[0] / "calibration_results.tsv"
    original = target.read_bytes()
    result = run_audit(before, after, target)
    assert result.returncode == 2
    assert "must be a new file" in result.stderr
    assert target.read_bytes() == original
