#!/usr/bin/env python3
"""Scratch-isolated mutation verifier for issue #1033 Item 2.

Purpose:
    Verify that tests/test_figure_methods_1033.py kills mutations to
    trait-visualization/figure_methods.md (:19-25) and
    docs/reviews/2026-10-23-final-rerun-record.md (:1) using pytest --junitxml
    call-phase AssertionError verification on the target test node and writing a
    persistent verdicts.json in --scratch-dir without modifying the live worktree.

Inputs:
    - trait-visualization/figure_methods.md
    - docs/reviews/2026-10-23-final-rerun-record.md
    - tests/test_figure_methods_1033.py

Outputs:
    - <scratch-dir>/verdicts.json
    - <scratch-dir>/summary.json
    - <scratch-dir>/<case>/junit.xml

Issue: #1033
Date: 2026-11-11
"""

from __future__ import annotations

import argparse
import hashlib
import html
import importlib
import io
import json
import re
import shutil
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
FIGURE_METHODS_REL = Path("trait-visualization/figure_methods.md")
RERUN_RECORD_REL = Path("docs/reviews/2026-10-23-final-rerun-record.md")
TEST_REL = Path("tests/test_figure_methods_1033.py")
TARGET_TEST_NODE = "test_figure_methods_clustered_heatmaps_cites_canonical_final_rerun_record"


@dataclass(frozen=True)
class MutationSpec:
    """Single scratch-isolated mutation specification for #1033 Item 2."""

    name: str
    rel_path: Path
    needle: str
    replacement: str
    test_node: str
    expected_detail_substring: str


@dataclass(frozen=True)
class CallPhaseOutcome:
    """Recorded pytest call-phase and JUnit XML outcome for a single test node."""

    exit_code: int
    junit_xml_path: str
    junit_tests: int
    junit_failures: int
    junit_errors: int
    junit_skipped: int
    junit_testcase_matched: bool
    junit_assertion_error: bool
    setup_passed: bool
    call_outcome: str
    teardown_passed: bool
    detail: str


@dataclass(frozen=True)
class MutationVerdict:
    """Persistent verdict entry written to <scratch-dir>/verdicts.json."""

    name: str
    target_file: str
    test_node: str
    verdict: str
    killed: bool
    exit_code: int
    original_sha256: str
    mutated_sha256: str
    junit_xml_path: str
    call_phase_assertion_error: bool
    expected_detail_substring: str
    matched_detail: bool
    detail: str


MUTATIONS: tuple[MutationSpec, ...] = (
    MutationSpec(
        name="revert_unnamed_final_rerun_record",
        rel_path=FIGURE_METHODS_REL,
        needle=(
            "from the canonical final rerun record in\n"
            "[`docs/reviews/2026-10-23-final-rerun-record.md`]"
            "(../docs/reviews/2026-10-23-final-rerun-record.md)\n"
            "(#1004, reconciling the historical pre-fix tree and the corrected 19-table rerun)"
        ),
        replacement="from the final rerun record",
        test_node=TARGET_TEST_NODE,
        expected_detail_substring=(
            "Unnamed 'from the final rerun record before reusing them' must be replaced"
        ),
    ),
    MutationSpec(
        name="drop_canonical_record_path",
        rel_path=FIGURE_METHODS_REL,
        needle=(
            "[`docs/reviews/2026-10-23-final-rerun-record.md`]"
            "(../docs/reviews/2026-10-23-final-rerun-record.md)"
        ),
        replacement="[final rerun record](../docs/reviews/2026-10-23-final-rerun-record.md)",
        test_node=TARGET_TEST_NODE,
        expected_detail_substring=(
            "Clustered Heatmaps intro must explicitly cite "
            "`docs/reviews/2026-10-23-final-rerun-record.md`"
        ),
    ),
    MutationSpec(
        name="drop_relative_link",
        rel_path=FIGURE_METHODS_REL,
        needle="(../docs/reviews/2026-10-23-final-rerun-record.md)",
        replacement="(docs/reviews/2026-10-23-final-rerun-record.md)",
        test_node=TARGET_TEST_NODE,
        expected_detail_substring=(
            "Clustered Heatmaps intro must include relative link "
            "(../docs/reviews/2026-10-23-final-rerun-record.md)"
        ),
    ),
    MutationSpec(
        name="drop_canonical_label",
        rel_path=FIGURE_METHODS_REL,
        needle="canonical final rerun record",
        replacement="final rerun note",
        test_node=TARGET_TEST_NODE,
        expected_detail_substring=(
            "Clustered Heatmaps intro must identify the document as the "
            "canonical final rerun record"
        ),
    ),
    MutationSpec(
        name="drop_reconciling_pre_fix_and_19_table_clause",
        rel_path=FIGURE_METHODS_REL,
        needle="reconciling the historical pre-fix tree and the corrected 19-table rerun",
        replacement="recording the rerun",
        test_node=TARGET_TEST_NODE,
        expected_detail_substring=(
            "Clustered Heatmaps intro must state that the final rerun record reconciles "
            "the historical pre-fix tree and the corrected 19-table rerun"
        ),
    ),
    MutationSpec(
        name="drop_commit_511c891_1004",
        rel_path=FIGURE_METHODS_REL,
        needle="frozen code commit `511c891` (#1004)",
        replacement="frozen code commit (#1004)",
        test_node=TARGET_TEST_NODE,
        expected_detail_substring=(
            "Clustered Heatmaps intro must preserve frozen code commit `511c891` (#1004)"
        ),
    ),
    MutationSpec(
        name="drop_962_subsequent_runs_clause",
        rel_path=FIGURE_METHODS_REL,
        needle="added in #962",
        replacement="added earlier",
        test_node=TARGET_TEST_NODE,
        expected_detail_substring=(
            "Clustered Heatmaps intro must preserve the #962 subsequent-runs scope sentence"
        ),
    ),
    MutationSpec(
        name="corrupt_rerun_record_heading",
        rel_path=RERUN_RECORD_REL,
        needle="# Final rerun record: commit 511c891 (2026-10-23)",
        replacement="# Draft rerun notes",
        test_node=TARGET_TEST_NODE,
        expected_detail_substring=(
            "Unexpected header in docs/reviews/2026-10-23-final-rerun-record.md"
        ),
    ),
)


class _CallPhaseRecorder:
    """In-process pytest plugin recording setup/call/teardown outcomes for the target node."""

    def __init__(self, target_node: str) -> None:
        self.target_node = target_node
        self.setup_passed = False
        self.call_outcome = "not_run"
        self.teardown_passed = False
        self.call_longrepr = ""

    def pytest_runtest_logreport(self, report: pytest.TestReport) -> None:
        if not report.nodeid.endswith(f"::{self.target_node}"):
            return
        if report.when == "setup":
            self.setup_passed = report.outcome == "passed"
        elif report.when == "call":
            self.call_outcome = report.outcome
            self.call_longrepr = report.longreprtext
        elif report.when == "teardown":
            self.teardown_passed = report.outcome == "passed"


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _copy_baseline_tree(scratch_case_dir: Path, originals: dict[Path, bytes]) -> None:
    if scratch_case_dir.exists():
        shutil.rmtree(scratch_case_dir)
    for rel, content in originals.items():
        dst = scratch_case_dir / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_bytes(content)


def _extract_int_attr(tag_text: str, attr_name: str) -> int:
    match = re.search(rf'\b{re.escape(attr_name)}="(\d+)"', tag_text)
    if match is None:
        raise RuntimeError(f"Missing integer attribute {attr_name!r} in JUnit XML tag: {tag_text}")
    return int(match.group(1))


def _parse_junit_xml(
    junit_xml_path: Path, expected_test_node: str
) -> tuple[int, int, int, int, bool, bool, str]:
    if not junit_xml_path.is_file():
        raise RuntimeError(f"Expected JUnit XML file at {junit_xml_path}")
    raw_xml = junit_xml_path.read_text(encoding="utf-8")
    suite_match = re.search(r"<testsuite\b[^>]*>", raw_xml)
    if suite_match is None:
        raise RuntimeError(f"Missing <testsuite> element in {junit_xml_path}")
    suite_tag = suite_match.group(0)
    tests_count = _extract_int_attr(suite_tag, "tests")
    failures_count = _extract_int_attr(suite_tag, "failures")
    errors_count = _extract_int_attr(suite_tag, "errors")
    skipped_count = _extract_int_attr(suite_tag, "skipped")

    testcase_match = re.search(
        rf'<testcase\b[^>]*\bname="{re.escape(expected_test_node)}"[^>]*>(.*?)</testcase>',
        raw_xml,
        flags=re.DOTALL,
    )
    self_closing_match = re.search(
        rf'<testcase\b[^>]*\bname="{re.escape(expected_test_node)}"[^>]*/\s*>',
        raw_xml,
    )
    testcase_matched = testcase_match is not None or self_closing_match is not None
    body = testcase_match.group(1) if testcase_match is not None else ""
    unescaped_body = html.unescape(body)
    has_error_tag = "<error" in body
    failure_match = re.search(
        r'<failure\b[^>]*\bmessage="([^"]*)"[^>]*>(.*?)</failure>',
        body,
        flags=re.DOTALL,
    )
    assertion_error = False
    failure_detail = ""
    if failure_match is not None and not has_error_tag:
        msg = html.unescape(failure_match.group(1))
        text = html.unescape(failure_match.group(2))
        failure_detail = f"{msg}\n{text}".strip()
        assertion_error = "AssertionError:" in msg or "AssertionError:" in text
    elif unescaped_body:
        failure_detail = unescaped_body.strip()
    return (
        tests_count,
        failures_count,
        errors_count,
        skipped_count,
        testcase_matched,
        assertion_error,
        failure_detail,
    )


def _run_scratch_test_node(scratch_case_dir: Path, test_node: str) -> CallPhaseOutcome:
    test_file = scratch_case_dir / TEST_REL
    junit_xml_path = scratch_case_dir / "junit.xml"
    sys.modules.pop("test_figure_methods_1033", None)
    importlib.invalidate_caches()
    recorder = _CallPhaseRecorder(test_node)
    stdout_buf = io.StringIO()
    stderr_buf = io.StringIO()
    old_stdout, old_stderr = sys.stdout, sys.stderr
    try:
        sys.stdout = stdout_buf
        sys.stderr = stderr_buf
        exit_code = int(
            pytest.main(
                [
                    "-q",
                    "--tb=short",
                    "--import-mode=importlib",
                    f"--junitxml={junit_xml_path}",
                    "-o",
                    f"cache_dir={scratch_case_dir / '.pytest_cache'}",
                    f"{test_file}::{test_node}",
                ],
                plugins=[recorder],
            )
        )
    finally:
        sys.stdout = old_stdout
        sys.stderr = old_stderr
        sys.modules.pop("test_figure_methods_1033", None)

    (
        tests_count,
        failures_count,
        errors_count,
        skipped_count,
        testcase_matched,
        junit_assertion_error,
        junit_detail,
    ) = _parse_junit_xml(junit_xml_path, test_node)

    combined_output = "\n".join(
        part
        for part in (
            junit_detail,
            recorder.call_longrepr.strip(),
            stdout_buf.getvalue().strip(),
            stderr_buf.getvalue().strip(),
        )
        if part
    ).strip()

    return CallPhaseOutcome(
        exit_code=exit_code,
        junit_xml_path=str(junit_xml_path),
        junit_tests=tests_count,
        junit_failures=failures_count,
        junit_errors=errors_count,
        junit_skipped=skipped_count,
        junit_testcase_matched=testcase_matched,
        junit_assertion_error=junit_assertion_error,
        setup_passed=recorder.setup_passed,
        call_outcome=recorder.call_outcome,
        teardown_passed=recorder.teardown_passed,
        detail=combined_output or "passed",
    )


def _evaluate_mutation(
    spec: MutationSpec, scratch_root: Path, originals: dict[Path, bytes]
) -> MutationVerdict:
    case_dir = scratch_root / spec.name
    _copy_baseline_tree(case_dir, originals)

    target_path = case_dir / spec.rel_path
    original_bytes = originals[spec.rel_path]
    original_text = original_bytes.decode("utf-8")
    count = original_text.count(spec.needle)
    if count != 1:
        raise RuntimeError(
            f"Expected needle for {spec.name} to appear once in {spec.rel_path}, found {count}"
        )

    mutated_text = original_text.replace(spec.needle, spec.replacement, 1)
    mutated_bytes = mutated_text.encode("utf-8")
    target_path.write_bytes(mutated_bytes)

    outcome = _run_scratch_test_node(case_dir, spec.test_node)
    call_phase_assertion_error = (
        outcome.exit_code == 1
        and outcome.junit_tests == 1
        and outcome.junit_failures == 1
        and outcome.junit_errors == 0
        and outcome.junit_skipped == 0
        and outcome.junit_testcase_matched
        and outcome.junit_assertion_error
        and outcome.setup_passed
        and outcome.call_outcome == "failed"
        and outcome.teardown_passed
    )
    matched_detail = spec.expected_detail_substring in outcome.detail
    killed = call_phase_assertion_error and matched_detail
    return MutationVerdict(
        name=spec.name,
        target_file=str(spec.rel_path),
        test_node=spec.test_node,
        verdict="killed by call-phase AssertionError" if killed else "survived",
        killed=killed,
        exit_code=outcome.exit_code,
        original_sha256=_sha256_bytes(original_bytes),
        mutated_sha256=_sha256_bytes(mutated_bytes),
        junit_xml_path=outcome.junit_xml_path,
        call_phase_assertion_error=call_phase_assertion_error,
        expected_detail_substring=spec.expected_detail_substring,
        matched_detail=matched_detail,
        detail=outcome.detail,
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--scratch-dir",
        "--output-dir",
        dest="scratch_dir",
        type=Path,
        required=True,
        help="Persistent scratch directory outside the worktree for mutant copies and verdicts.",
    )
    parser.add_argument(
        "--canary",
        action="store_true",
        help="Run only the first mutation specification as a single-unit canary.",
    )
    parser.add_argument(
        "--mutation",
        type=str,
        default=None,
        help="Run only the named mutation specification.",
    )
    args = parser.parse_args(argv)

    scratch_root: Path = args.scratch_dir.resolve()
    if scratch_root == REPO_ROOT or scratch_root.is_relative_to(REPO_ROOT):
        raise ValueError("Scratch directory must be outside the live worktree")
    scratch_root.mkdir(parents=True, exist_ok=True)

    tracked_inputs = (FIGURE_METHODS_REL, RERUN_RECORD_REL, TEST_REL)
    originals = {rel: (REPO_ROOT / rel).read_bytes() for rel in tracked_inputs}

    selected: list[MutationSpec]
    if args.mutation is not None:
        selected = [m for m in MUTATIONS if m.name == args.mutation]
        if not selected:
            raise ValueError(f"Unknown mutation name: {args.mutation}")
    elif args.canary:
        selected = [MUTATIONS[0]]
    else:
        selected = list(MUTATIONS)

    verdicts: list[MutationVerdict] = []
    try:
        baseline_dir = scratch_root / "_control"
        _copy_baseline_tree(baseline_dir, originals)
        baseline_outcome = _run_scratch_test_node(baseline_dir, TARGET_TEST_NODE)
        baseline_ok = (
            baseline_outcome.exit_code == 0
            and baseline_outcome.junit_tests == 1
            and baseline_outcome.junit_failures == 0
            and baseline_outcome.junit_errors == 0
            and baseline_outcome.junit_skipped == 0
            and baseline_outcome.junit_testcase_matched
            and baseline_outcome.setup_passed
            and baseline_outcome.call_outcome == "passed"
            and baseline_outcome.teardown_passed
        )
        if not baseline_ok:
            raise RuntimeError(
                f"Baseline control test failed before mutation run: {baseline_outcome.detail}"
            )

        for spec in selected:
            verdict = _evaluate_mutation(spec, scratch_root, originals)
            verdicts.append(verdict)
    finally:
        for rel, content in originals.items():
            if (REPO_ROOT / rel).read_bytes() != content:
                raise RuntimeError(f"Live worktree input changed during verification: {rel}")

    killed_count = sum(1 for v in verdicts if v.killed)
    survived = [v.name for v in verdicts if not v.killed]
    mode = "canary" if args.canary else ("single" if args.mutation else "full")

    verdicts_payload = {
        "issue": 1033,
        "item": 2,
        "mode": mode,
        "scratch_dir": str(scratch_root),
        "control": {
            "test_node": TARGET_TEST_NODE,
            "verdict": "passed",
            "exit_code": baseline_outcome.exit_code,
            "junit_xml_path": baseline_outcome.junit_xml_path,
            "detail": baseline_outcome.detail,
        },
        "verdicts": [asdict(v) for v in verdicts],
    }
    summary_payload = {
        "issue": 1033,
        "item": 2,
        "mode": mode,
        "scratch_dir": str(scratch_root),
        "control_passed": True,
        "total": len(verdicts),
        "killed": killed_count,
        "survived": len(survived),
        "survived_names": survived,
        "all_killed": len(survived) == 0,
        "verdict_summaries": [
            {
                "name": v.name,
                "target_file": v.target_file,
                "test_node": v.test_node,
                "verdict": v.verdict,
                "killed": v.killed,
                "call_phase_assertion_error": v.call_phase_assertion_error,
                "matched_detail": v.matched_detail,
                "junit_xml_path": v.junit_xml_path,
                "detail_preview": v.detail.splitlines()[0] if v.detail else "",
            }
            for v in verdicts
        ],
    }

    verdicts_path = scratch_root / "verdicts.json"
    summary_path = scratch_root / "summary.json"
    verdicts_path.write_text(json.dumps(verdicts_payload, indent=2) + "\n", encoding="utf-8")
    summary_path.write_text(json.dumps(summary_payload, indent=2) + "\n", encoding="utf-8")

    print(json.dumps(summary_payload, indent=2))
    return 0 if not survived else 1


if __name__ == "__main__":
    raise SystemExit(main())
