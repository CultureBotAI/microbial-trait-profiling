"""Bind consolidated conclusions to reported statuses and explicit provenance (#1069)."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import inspect
import re
from dataclasses import replace
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

from mtp.density_adjustment import Table1TraitEffect, family_bh_q_values
from mtp.trait_lists import DeclaredTraitList
from test_density_consolidated_conclusions import synthetic_results

_ROOT = Path(__file__).resolve().parents[1]
_CONCLUSIONS = "6. CONCLUSION & METHODOLOGICAL IMPLICATIONS\n"
_SURVIVED = "survived (significant post-adjustment)"
_ENHANCED = "enhanced (significant, larger magnitude)"
_COLLAPSED = "collapsed (attenuated to nonsignificant)"
_ADJUSTED_ONLY = "significant only post-adjustment"
_REVERSED = "reversed (significant post-adjustment)"
_UNDEFINED = "undefined magnitude change (raw mean zero, significant post-adjustment)"


@pytest.fixture(scope="module")
def audit_1069() -> ModuleType:
    path = _ROOT / "scripts/audit_trait_density_correlation.py"
    spec = importlib.util.spec_from_file_location("density_report_1069_under_test", path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _results(statuses: tuple[str, ...]) -> tuple[dict[str, Any], dict[str, Any]]:
    density, ontology = synthetic_results()
    template = density["table1_effects"][0]
    values = {
        _SURVIVED: (0.05, 0.02, 0.01, 0.03),
        _ENHANCED: (0.05, 0.02, 0.10, 0.03),
        _COLLAPSED: (0.05, 0.02, 0.01, 0.50),
        _ADJUSTED_ONLY: (0.05, 0.50, 0.01, 0.03),
        _REVERSED: (0.05, 0.02, -0.01, 0.03),
        _UNDEFINED: (0.00, 0.50, 0.01, 0.03),
    }
    density["table1_effects"] = [
        replace(
            template,
            trait_id=f"CUSTOM:{number}",
            trait_label=f"custom annotation {number}",
            status=status,
            raw_mean_diff=values[status][0],
            raw_wilcoxon_bh_q=values[status][1] + number / 1000,
            adj_mean_diff=values[status][2],
            adj_wilcoxon_bh_q=values[status][3] + number / 1000,
            bh_family_size=len(statuses),
        )
        for number, status in enumerate(statuses, start=1)
    ]
    return density, ontology


def _declared(effects: list[Table1TraitEffect]) -> DeclaredTraitList:
    identifiers = tuple(effect.trait_id for effect in effects)
    digest = hashlib.sha256("\n".join(identifiers).encode()).hexdigest()
    return DeclaredTraitList(identifiers, "custom_traits.txt", digest)


def _render(
    audit: ModuleType,
    output: Path,
    density: dict[str, Any],
    ontology: dict[str, Any],
    declared: bool,
    **kwargs: str,
) -> str:
    traits = _declared(density["table1_effects"]) if declared else None
    path: Path = audit.write_consolidated_report(output, density, ontology, traits, **kwargs)
    assert path.is_file()
    return path.read_text(encoding="utf-8")


def _section3(content: str) -> str:
    match = re.search(r"(?ms)^3\. .*?(?=^4\. )", content)
    assert match is not None, "Report is missing its effect table"
    return match.group()


def _status_rows(content: str) -> list[tuple[str, str]]:
    return [
        (line.split("|", 1)[0].strip(), line.rsplit("|", 1)[1].strip())
        for line in _section3(content).splitlines()
        if "| Raw mean:" in line
    ]


def _status_and_q_rows(content: str) -> list[tuple[str, str, str, str]]:
    rows = []
    for line in _section3(content).splitlines():
        if "| Raw mean:" not in line:
            continue
        columns = [column.strip() for column in line.split("|")]
        q_values = re.findall(r", q=([^ )]+)\)", line)
        assert len(q_values) == 2, "Each reported row must have raw and adjusted q fields"
        rows.append((columns[0], columns[-1], q_values[0], q_values[1]))
    return rows


def _survivors(content: str) -> list[str]:
    assert content.count(_CONCLUSIONS) == 1
    conclusions = content.split(_CONCLUSIONS, 1)[1]
    match = re.search(
        r"(?m)^3\. Significant before and after PC1 regression in the same direction "
        r"\(BH q <= 0\.05; survived or enhanced\): (.*)\.$",
        conclusions,
    )
    assert match is not None, "Section 6 must summarize the actual section 3 survivor set"
    return [] if match[1] == "none" else match[1].split(", ")


@pytest.mark.parametrize("declared", [False, True], ids=["default", "declared"])
@pytest.mark.parametrize(
    "statuses",
    [
        (_COLLAPSED, _ADJUSTED_ONLY, _REVERSED),
        (_SURVIVED, _ENHANCED, _COLLAPSED, _ADJUSTED_ONLY, _REVERSED, _UNDEFINED),
        (_ENHANCED,),
    ],
    ids=["no-survivors", "partial-survivors", "custom-singleton"],
)
def test_section6_agrees_with_section3_for_every_reported_set(
    audit_1069: ModuleType,
    tmp_path: Path,
    declared: bool,
    statuses: tuple[str, ...],
) -> None:
    density, ontology = _results(statuses)
    content = _render(audit_1069, tmp_path, density, ontology, declared)
    rows = _status_rows(content)
    assert rows == [(effect.trait_label, effect.status) for effect in density["table1_effects"]]
    expected = [label for label, status in rows if status in {_SURVIVED, _ENHANCED}]
    assert _survivors(content) == expected
    conclusions = content.split(_CONCLUSIONS, 1)[1]
    assert "facultatively aerobic" not in conclusions
    assert "Not written for a declared trait list" not in conclusions
    status_summary = re.search(r"(?ms)^2\. .*?(?=^3\. )", conclusions)
    assert status_summary is not None
    for label, status in rows:
        assert label in status_summary.group()
        assert f"'{status}'" in status_summary.group()
    quantified_rows = _status_and_q_rows(content)
    assert quantified_rows == [
        (
            effect.trait_label,
            effect.status,
            f"{effect.raw_wilcoxon_bh_q:.2e}",
            f"{effect.adj_wilcoxon_bh_q:.2e}",
        )
        for effect in density["table1_effects"]
    ]
    summary_rows = re.findall(
        r"(?m)^   (.+?): '(.*?)' \(raw BH q=([^,]+), adjusted BH q=([^)]+)\)\.$",
        status_summary.group(),
    )
    assert summary_rows == quantified_rows, "Every label must retain its own status and q fields"


def test_explicit_code_commit_and_unbound_note_do_not_change_table_values(
    audit_1069: ModuleType, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    assert "code_commit" in inspect.signature(audit_1069.write_consolidated_report).parameters
    density, ontology = _results((_SURVIVED, _COLLAPSED))
    monkeypatch.chdir(tmp_path)
    assert not (tmp_path / ".git").exists()
    historical = _render(audit_1069, tmp_path, density, ontology, False)
    assert "Provenance note: code_commit not supplied; producer commit is unverified." in (
        historical
    )
    commit = "fixture-commit"
    bound = _render(audit_1069, tmp_path, density, ontology, False, code_commit=commit)
    assert f"code_commit: {commit}\n" in bound
    assert "Provenance note: code_commit not supplied" not in bound
    assert _section3(bound) == _section3(historical)


def test_main_passes_explicit_code_commit_to_the_consolidated_report(
    audit_1069: ModuleType, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    density, ontology = _results((_SURVIVED, _COLLAPSED))
    density.update(summary_tsv=Path("summary.tsv"), table1_tsv=Path("effects.tsv"))
    ontology.update(
        pairs_tsv=Path("pairs.tsv"),
        summary_tsv=Path("ontology.tsv"),
        comparison_tsv=Path("comparison.tsv"),
    )
    commit = "fixture-cli-commit"
    args = argparse.Namespace(
        profiles=Path("unused_profiles.tsv"),
        columns=Path("unused_columns.tsv"),
        output_dir=tmp_path,
        code_commit=commit,
        trait_list=None,
    )
    captured: dict[str, Any] = {}

    def density_run(**kwargs: Any) -> dict[str, Any]:
        captured.update(kwargs)
        return density

    monkeypatch.setattr(audit_1069, "parse_args", lambda: args)
    monkeypatch.setattr(audit_1069, "run_density_audit", density_run)
    monkeypatch.setattr(audit_1069, "run_ontology_redundancy_audit", lambda **_kwargs: ontology)
    assert audit_1069.main() == 0
    assert captured["code_commit"] == commit
    report = (tmp_path / "density_adjustment_report.txt").read_text(encoding="utf-8")
    assert f"code_commit: {commit}\n" in report


def _historical_effects(document: str) -> list[Table1TraitEffect]:
    section = document.split("## 5. Paired Contrasts", 1)[1].split("## 6.", 1)[0]
    rows = [
        [cell.strip().strip("*") for cell in re.split(r"(?<!\\)\|", line)[1:-1]]
        for line in section.splitlines()
        if line.startswith("| `METPO:")
    ]
    assert len(rows) == 6

    def number(cell: str) -> float:
        return float(cell.strip("$"))

    def probability(cell: str) -> float:
        match = re.fullmatch(r"\$(\d+\.\d+) \\times 10\^\{(-?\d+)\}\$", cell)
        assert match is not None
        return float(f"{match[1]}e{match[2]}")

    raw_q = family_bh_q_values([probability(row[3]) for row in rows], len(rows))
    adjusted_q = family_bh_q_values([probability(row[5]) for row in rows], len(rows))
    density, _ = synthetic_results()
    template = density["table1_effects"][0]
    return [
        replace(
            template,
            trait_id=row[0].strip("`").replace(r"\|", "|"),
            trait_label=row[1],
            raw_mean_diff=number(row[2]),
            raw_wilcoxon_p=probability(row[3]),
            adj_mean_diff=number(row[4]),
            adj_wilcoxon_p=probability(row[5]),
            status=row[8].casefold(),
            raw_wilcoxon_bh_q=raw_q[index],
            adj_wilcoxon_bh_q=adjusted_q[index],
            bh_family_size=len(rows),
        )
        for index, row in enumerate(rows)
    ]


def test_documented_report_contract_matches_historical_table_and_generator(
    audit_1069: ModuleType, tmp_path: Path
) -> None:
    document = (_ROOT / "docs/density_adjustment_audit.md").read_text(encoding="utf-8")
    assert "### Consolidated report contract (#1069)" in document
    contract = document.split("### Consolidated report contract (#1069)", 1)[1].split("\n## ", 1)[0]
    assert "report section 3 corresponds to this document's section 5" in contract
    assert "report section 6" in contract
    assert "code_commit" in contract
    assert "Provenance note: code_commit not supplied; producer commit is unverified." in (contract)
    assert "leaves out the conclusions of section 6" not in document
    density, ontology = synthetic_results()
    density["table1_effects"] = _historical_effects(document)
    report = _render(audit_1069, tmp_path, density, ontology, False)
    expected = _survivors(report)
    match = re.search(r"(?m)^Historical surviving labels: (.*)\.$", contract)
    assert match is not None
    assert match[1].split("; ") == expected
