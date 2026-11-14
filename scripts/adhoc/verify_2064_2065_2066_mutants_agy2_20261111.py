#!/usr/bin/env python3
"""Verify tests/test_taxa_coverage_2064_2065_2066.py kills mutations of #2064, #2065, and #2066.

Purpose:
    Exercise tests/test_taxa_coverage_2064_2065_2066.py against synthetic mutants
    of scripts/create_taxa_coverage_histogram.py and docs/workflow.md in an
    isolated scratch directory, confirming baseline PASS (exit code 0) and
    call-phase AssertionError kills (exit code 1) via both pytest hook receipts
    and JUnit XML while verifying canonical files remain unchanged.
Inputs:
    - scripts/create_taxa_coverage_histogram.py
    - docs/workflow.md
    - tests/test_taxa_coverage_2064_2065_2066.py
Outputs:
    - <scratch-dir>/<canary-or-full-run>/verdicts.json
    - <scratch-dir>/verdicts.json
Issue:
    Closes #2064, #2065, #2066
Date:
    2026-11-11
"""

from __future__ import annotations

import argparse
import hashlib
import html
import importlib
import json
import os
import re
import shutil
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_REL = Path("scripts/create_taxa_coverage_histogram.py")
WORKFLOW_REL = Path("docs/workflow.md")
TEST_REL = Path("tests/test_taxa_coverage_2064_2065_2066.py")
STAGED_RELS = (SCRIPT_REL, WORKFLOW_REL, TEST_REL)
TEST_MODULE = "test_taxa_coverage_2064_2065_2066"
LOADED_MODULE = "taxa_coverage_histogram_2064"
TARGET_NODE_MARKERS = (
    "test_obo_form_rank_classes_under_taxonomic_rank_skip_without_counting_as_ncbitaxon_rank",
    "test_excluded_records_note_pins_always_zero_and_ncbitaxon_name_rank_scope",
    "test_non_rank_ncbitaxon_subject_under_taxonomic_rank_raises_distinct_error",
    "test_readme_pins_always_zero_strict_edge_parser_wording",
    "test_workflow_doc_853_pins_non_rank_ncbitaxon_subjects_under_taxonomic_rank",
)

_TESTCASE_RE = re.compile(
    r'<testcase(?P<attrs>(?:\s+\w+="[^"]*")*\s*)(?:/>|>(?P<body>.*?)</testcase>)',
    re.DOTALL,
)
_FAILURE_RE = re.compile(
    r'<failure(?P<attrs>(?:\s+\w+="[^"]*")*\s*)(?:/>|>(?P<body>.*?)</failure>)',
    re.DOTALL,
)
_ERROR_RE = re.compile(r'<error(?P<attrs>(?:\s+\w+="[^"]*")*\s*)', re.DOTALL)
_ATTR_RE = re.compile(r'(\w+)="([^"]*)"')


@dataclass(frozen=True)
class MutantSpec:
    """Single-target mutation specification."""

    name: str
    rel_path: Path
    needle: str
    replacement: str


_RANK_ERR_LINE = (
    "                        "
    "f\"Rank-class edge {row['id']} has a non-rank NCBITaxon subject: {subject}\"\n"
)

MUTANTS: tuple[MutantSpec, ...] = (
    MutantSpec(
        name="rank-root-refuses-any-nonletter-subject",
        rel_path=SCRIPT_REL,
        needle=(
            '                if subject.startswith("NCBITaxon:"):\n'
            "                    raise ValueError(\n" + _RANK_ERR_LINE + "                    )"
        ),
        replacement=(
            '                if not re.fullmatch(r"NCBITaxon:[A-Za-z_][A-Za-z0-9_]*", subject):\n'
            "                    raise ValueError(\n" + _RANK_ERR_LINE + "                    )"
        ),
    ),
    MutantSpec(
        name="rank-root-counts-obo-form",
        rel_path=SCRIPT_REL,
        needle='                if re.fullmatch(r"NCBITaxon:[A-Za-z_][A-Za-z0-9_]*", subject):',
        replacement=(
            "                if re.fullmatch("
            'r"(?:NCBITaxon:|OBO:NCBITaxon#_)[A-Za-z_][A-Za-z0-9_]*", subject):'
        ),
    ),
    MutantSpec(
        name="note-drops-always-zero",
        rel_path=SCRIPT_REL,
        needle=(
            '                "excluded_records counts malformed KG edge records '
            'dropped by the strict "\n'
            '                "edge parser (always 0); '
            'taxonomy_rank_class_edge_records_skipped counts "'
        ),
        replacement=(
            '                "excluded_records counts malformed KG edge records '
            'dropped by the strict "\n'
            '                "edge parser; '
            'taxonomy_rank_class_edge_records_skipped counts "'
        ),
    ),
    MutantSpec(
        name="note-drops-ncbitaxon-name",
        rel_path=SCRIPT_REL,
        needle=(
            '                "well-formed NCBITaxon:<name> ontology rank-class '
            'subclass_of records "\n'
            '                "skipped during organism-ancestry graph projection, '
            'while OBO-form rank-class "\n'
            '                "subjects take the non-NCBITaxon skip (#1008, #1274)."'
        ),
        replacement=(
            '                "well-formed ontology rank-class subclass_of '
            'records skipped during "\n'
            '                "organism-ancestry graph projection (#1008, #1274)."'
        ),
    ),
    MutantSpec(
        name="revert-rank-edge-error-message",
        rel_path=SCRIPT_REL,
        needle="f\"Rank-class edge {row['id']} has a non-rank NCBITaxon subject: {subject}\"",
        replacement="f\"Invalid taxonomy identifier in edge {row['id']}: {subject}\"",
    ),
    MutantSpec(
        name="rank-root-accepts-numeric-ncbitaxon",
        rel_path=SCRIPT_REL,
        needle=(
            '                if subject.startswith("NCBITaxon:"):\n'
            "                    raise ValueError(\n" + _RANK_ERR_LINE + "                    )"
        ),
        replacement=(
            '                if subject.startswith("NCBITaxon:"):\n                    continue'
        ),
    ),
    MutantSpec(
        name="revert-readme-always-zero",
        rel_path=SCRIPT_REL,
        needle=(
            '            "excluded_records is always 0: the strict KG edge parser raises on a "\n'
            '            "malformed record instead of dropping it).",'
        ),
        replacement=(
            '            "excluded_records: 0 confirms the strict KG edge '
            'parser dropped no record).",'
        ),
    ),
    MutantSpec(
        name="revert-workflow-doc-853",
        rel_path=WORKFLOW_REL,
        needle=(
            "Malformed NCBITaxon identifiers in organism ancestry edges and non-rank "
            "NCBITaxon subjects under OBO:NCBITaxon#_taxonomic_rank still fail. An input"
        ),
        replacement=(
            "Malformed NCBITaxon identifiers in organism ancestry edges still fail. An input"
        ),
    ),
)


class Receipt:
    """Collect call-phase outcomes and exception types from in-process pytest runs."""

    def __init__(self) -> None:
        self.call_outcomes: dict[str, str] = {}
        self.call_exc_names: dict[str, str] = {}
        self.call_details: dict[str, str] = {}
        self.non_call_failures: list[str] = []

    def pytest_runtest_logreport(self, report: pytest.TestReport) -> None:
        if report.when == "call":
            self.call_outcomes[report.nodeid] = report.outcome
        elif report.failed:
            self.non_call_failures.append(f"{report.when}:{report.nodeid}")

    def pytest_exception_interact(
        self,
        node: pytest.Item | pytest.Collector,
        call: pytest.CallInfo[object],
        report: pytest.CollectReport | pytest.TestReport,
    ) -> None:
        _ = report
        if call.when == "call" and call.excinfo is not None:
            self.call_exc_names[node.nodeid] = call.excinfo.typename
            self.call_details[node.nodeid] = str(call.excinfo.value).splitlines()[0]


def _ParseAttrs(raw: str) -> dict[str, str]:
    return {match.group(1): html.unescape(match.group(2)) for match in _ATTR_RE.finditer(raw)}


def _CanonicalDigests() -> dict[str, str]:
    return {
        str(rel_path): hashlib.sha256((REPO_ROOT / rel_path).read_bytes()).hexdigest()
        for rel_path in STAGED_RELS
    }


def _StageScratchWorkspace(workspace: Path) -> None:
    for rel_path in STAGED_RELS:
        source = REPO_ROOT / rel_path
        dest = workspace / rel_path
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, dest)


def _InspectJunitAssertionFailures(junit_xml: Path) -> tuple[list[str], list[str]]:
    if not junit_xml.is_file() or junit_xml.stat().st_size == 0:
        return [], ["missing_or_empty_junit_xml"]
    text = junit_xml.read_text(encoding="utf-8")
    killed_nodes: list[str] = []
    junit_errors: list[str] = []
    for tc_match in _TESTCASE_RE.finditer(text):
        tc_attrs = _ParseAttrs(tc_match.group("attrs"))
        name = tc_attrs.get("name", "")
        body = tc_match.group("body") or ""
        for err_match in _ERROR_RE.finditer(body):
            err_attrs = _ParseAttrs(err_match.group("attrs"))
            junit_errors.append(f"{name}:{err_attrs.get('type', 'error')}")
        if not any(marker in name for marker in TARGET_NODE_MARKERS):
            continue
        if _ERROR_RE.search(body):
            continue
        for fail_match in _FAILURE_RE.finditer(body):
            fail_attrs = _ParseAttrs(fail_match.group("attrs"))
            failure_type = fail_attrs.get("type", "")
            failure_msg = fail_attrs.get("message", "")
            failure_body = html.unescape(fail_match.group("body") or "")
            if (
                failure_type == "AssertionError"
                or failure_type.endswith(".AssertionError")
                or failure_msg.startswith("AssertionError")
                or "\nAssertionError:" in failure_body
                or ": AssertionError" in failure_body
            ):
                killed_nodes.append(name)
    return killed_nodes, junit_errors


def _RunTargetPytest(workspace: Path, junit_xml: Path) -> tuple[int, Receipt, list[str], list[str]]:
    receipt = Receipt()
    test_file = workspace / TEST_REL
    saved_cwd = Path.cwd()
    saved_sys_path = list(sys.path)
    for mod_name in (TEST_MODULE, f"tests.{TEST_MODULE}", LOADED_MODULE):
        sys.modules.pop(mod_name, None)
    importlib.invalidate_caches()
    os.environ["MPLBACKEND"] = "Agg"
    try:
        os.chdir(workspace)
        sys.path.insert(0, str(workspace / "tests"))
        code = int(
            pytest.main(
                [
                    "-q",
                    "-o",
                    "addopts=",
                    f"--rootdir={workspace}",
                    f"--junitxml={junit_xml}",
                    str(test_file),
                ],
                plugins=[receipt],
            )
        )
    finally:
        for mod_name in (TEST_MODULE, f"tests.{TEST_MODULE}", LOADED_MODULE):
            sys.modules.pop(mod_name, None)
        importlib.invalidate_caches()
        sys.path[:] = saved_sys_path
        os.chdir(saved_cwd)
    junit_assertion_nodes, junit_errors = _InspectJunitAssertionFailures(junit_xml)
    return code, receipt, junit_assertion_nodes, junit_errors


def _ApplyMutant(workspace: Path, spec: MutantSpec) -> None:
    target = workspace / spec.rel_path
    original = target.read_text(encoding="utf-8")
    count = original.count(spec.needle)
    if count != 1:
        raise RuntimeError(
            f"Expected 1 occurrence of needle for {spec.name} in {spec.rel_path}, found {count}"
        )
    target.write_text(original.replace(spec.needle, spec.replacement, 1), encoding="utf-8")


def _ValidateScratchDir(scratch_dir: Path, repo_root: Path = REPO_ROOT) -> Path:
    scratch_base = scratch_dir.resolve()
    repo_resolved = repo_root.resolve()
    if (
        scratch_base == repo_resolved
        or scratch_base.is_relative_to(repo_resolved)
        or repo_resolved.is_relative_to(scratch_base)
    ):
        raise ValueError(f"--scratch-dir must be disjoint from {repo_resolved}, got {scratch_base}")
    return scratch_base


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--scratch-dir",
        type=Path,
        required=True,
        help="Scratch directory outside the git worktree where verdicts.json is persisted.",
    )
    parser.add_argument(
        "--canary",
        action="store_true",
        help="Run baseline + first mutant only as a canary before the full sweep.",
    )
    args = parser.parse_args(argv)

    scratch_base = _ValidateScratchDir(args.scratch_dir)
    scratch_base.mkdir(parents=True, exist_ok=True)
    prefix = "canary-" if args.canary else "full-"
    run_dir = Path(tempfile.mkdtemp(prefix=prefix, dir=scratch_base))

    before_digests = _CanonicalDigests()

    baseline_ws = run_dir / "baseline"
    _StageScratchWorkspace(baseline_ws)
    baseline_xml = run_dir / "baseline-junit.xml"
    base_code, base_receipt, _, base_junit_errors = _RunTargetPytest(baseline_ws, baseline_xml)
    passed_targets = [
        nodeid
        for nodeid, outcome in base_receipt.call_outcomes.items()
        if outcome == "passed" and any(marker in nodeid for marker in TARGET_NODE_MARKERS)
    ]
    if (
        base_code != int(pytest.ExitCode.OK)
        or base_receipt.non_call_failures
        or base_junit_errors
        or len(passed_targets) < len(TARGET_NODE_MARKERS)
    ):
        raise RuntimeError(
            f"Baseline failed (code={base_code}, passed_targets={passed_targets}, "
            f"non_call_failures={base_receipt.non_call_failures}, junit_errors={base_junit_errors})"
        )

    selected = MUTANTS[:1] if args.canary else MUTANTS
    verdicts: list[dict[str, object]] = []
    for spec in selected:
        mutant_ws = run_dir / spec.name
        _StageScratchWorkspace(mutant_ws)
        _ApplyMutant(mutant_ws, spec)
        junit_xml = run_dir / f"{spec.name}-junit.xml"
        code, receipt, junit_assertion_nodes, junit_errors = _RunTargetPytest(mutant_ws, junit_xml)
        hook_assertion_nodes = [
            nodeid
            for nodeid, outcome in receipt.call_outcomes.items()
            if outcome == "failed"
            and receipt.call_exc_names.get(nodeid) == "AssertionError"
            and any(marker in nodeid for marker in TARGET_NODE_MARKERS)
        ]
        details = [
            receipt.call_details[nodeid]
            for nodeid in hook_assertion_nodes
            if receipt.call_details.get(nodeid)
        ]
        killed = (
            code == int(pytest.ExitCode.TESTS_FAILED)
            and not receipt.non_call_failures
            and not junit_errors
            and bool(hook_assertion_nodes)
            and bool(junit_assertion_nodes)
            and bool(details)
        )
        verdict: dict[str, object] = {
            "name": spec.name,
            "target": str(spec.rel_path),
            "exit_code": code,
            "killed": killed,
            "hook_assertion_nodes": hook_assertion_nodes,
            "junit_assertion_nodes": junit_assertion_nodes,
            "junit_errors": junit_errors,
            "non_call_failures": receipt.non_call_failures,
            "detail": details[0] if details else "",
        }
        verdicts.append(verdict)
        print(
            f"{spec.name}: killed={killed} exit_code={code} "
            f"hook_nodes={len(hook_assertion_nodes)} junit_nodes={len(junit_assertion_nodes)}"
        )
        if not killed:
            raise RuntimeError(f"Mutant survived or failed outside call AssertionError: {verdict}")

    after_digests = _CanonicalDigests()
    canonical_unchanged = before_digests == after_digests
    if not canonical_unchanged:
        raise RuntimeError(f"Canonical files changed during mutant run: {after_digests}")

    summary = {
        "mode": "canary" if args.canary else "full",
        "run_dir": str(run_dir),
        "baseline_passed_nodes": passed_targets,
        "mutant_count": len(verdicts),
        "all_killed": all(bool(item["killed"]) for item in verdicts),
        "canonical_unchanged": canonical_unchanged,
        "canonical_sha256": after_digests,
        "verdicts": verdicts,
    }
    run_verdicts = run_dir / "verdicts.json"
    top_verdicts = scratch_base / "verdicts.json"
    payload = json.dumps(summary, indent=2) + "\n"
    run_verdicts.write_text(payload, encoding="utf-8")
    top_verdicts.write_text(payload, encoding="utf-8")
    print(f"Wrote verdicts to {run_verdicts} and {top_verdicts}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
