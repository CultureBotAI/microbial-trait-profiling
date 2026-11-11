#!/usr/bin/env python3
"""Scratch-isolated mutation verifier for issues #2035 and #2047.

Purpose:
    Verify that unit tests in tests/test_enrichment_infographic_followup.py kill
    targeted mutations in scripts/analyze_enrichment_infographic_followup.py for:
      - #2035: refusing unlisted/unknown row conditions in save_traits_vs_organisms_heatmap
        instead of falling back to unlabelled '#64748b'
      - #2047: validating condition counts against COMPANION_HEATMAP_PALETTE_COLORS at the
        start of main() before creating output_dir or loading/writing tables, while keeping
        the backstop in save_traits_vs_organisms_heatmap.
    Every test run executes against an isolated scratch copy outside the worktree and
    verifies via pytest hooks and --junitxml that the target test item ran to its call
    phase and failed with an AssertionError.

Inputs:
    - scripts/analyze_enrichment_infographic_followup.py
    - tests/test_enrichment_infographic_followup.py

Outputs:
    - <scratch-dir>/verdicts.json and <scratch-dir>/<run-dir>/verdicts.json

Issue: #2035, #2047
Date: 2026-11-11
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import importlib
import io
import json
import os
import shutil
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT_REL = Path("scripts/analyze_enrichment_infographic_followup.py")
TEST_REL = Path("tests/test_enrichment_infographic_followup.py")
CONTROL_NODEIDS: tuple[str, ...] = (
    "tests/test_enrichment_infographic_followup.py::test_companion_heatmap_refuses_unlisted_row_conditions",
    "tests/test_enrichment_infographic_followup.py::test_main_refuses_excess_conditions_before_mkdir_or_loading_or_writing_tables",
    "tests/test_enrichment_infographic_followup.py::test_companion_heatmap_refuses_more_conditions_than_distinct_colours",
)


@dataclass(frozen=True)
class MutantSpec:
    """Single source-text mutation specification."""

    name: str
    issue: str
    rel_path: Path
    old: str
    new: str
    target_nodeid: str


@dataclass(frozen=True)
class TargetedPytestResult:
    """Structured result from running a single pytest nodeid in scratch."""

    rc: int
    output: str
    origin_ok: bool
    collected_nodeids: tuple[str, ...]
    setup_outcome: str | None
    call_outcome: str | None
    teardown_outcome: str | None
    call_exc_type: str | None
    call_is_assertion_error: bool
    call_exc_value: str | None
    call_tb_path: str | None
    call_tb_lineno: int | None
    non_call_exceptions: tuple[str, ...]
    junit_xml_path: str
    killed_by_call_assertion: bool


MUTANTS: tuple[MutantSpec, ...] = (
    MutantSpec(
        name="M1_restore_unknown_grey_fallback_2035",
        issue="#2035",
        rel_path=SCRIPT_REL,
        old=(
            "    if unlisted_conditions:\n"
            "        unique_unlisted = list(dict.fromkeys(unlisted_conditions))\n"
            "        msg = (\n"
            '            f"Companion heatmap row condition(s) {unique_unlisted} are not in "\n'
            '            f"TARGET_CONDITIONS {list(TARGET_CONDITIONS)}."\n'
            "        )\n"
            "        raise ValueError(msg)\n"
            "    row_colors = row_conditions.map(lambda condition: condition_palette[condition])"
        ),
        new=(
            "    row_colors = row_conditions.map(\n"
            '        lambda condition: condition_palette.get(condition, "#64748b")\n'
            "    )"
        ),
        target_nodeid=(
            "tests/test_enrichment_infographic_followup.py::"
            "test_companion_heatmap_refuses_unlisted_row_conditions[display-unknown]"
        ),
    ),
    MutantSpec(
        name="M2_allow_unlisted_non_unknown_condition_2035",
        issue="#2035",
        rel_path=SCRIPT_REL,
        old=(
            "    unlisted_conditions = [\n"
            "        clean_text(condition) for condition in row_conditions "
            "if condition not in condition_palette\n"
            "    ]\n"
            "    if unlisted_conditions:\n"
            "        unique_unlisted = list(dict.fromkeys(unlisted_conditions))\n"
            "        msg = (\n"
            '            f"Companion heatmap row condition(s) {unique_unlisted} are not in "\n'
            '            f"TARGET_CONDITIONS {list(TARGET_CONDITIONS)}."\n'
            "        )\n"
            "        raise ValueError(msg)\n"
            "    row_colors = row_conditions.map(lambda condition: condition_palette[condition])"
        ),
        new=(
            "    unlisted_conditions = [\n"
            "        clean_text(condition) for condition in row_conditions "
            'if clean_text(condition) == "unknown"\n'
            "    ]\n"
            "    if unlisted_conditions:\n"
            "        unique_unlisted = list(dict.fromkeys(unlisted_conditions))\n"
            "        msg = (\n"
            '            f"Companion heatmap row condition(s) {unique_unlisted} are not in "\n'
            '            f"TARGET_CONDITIONS {list(TARGET_CONDITIONS)}."\n'
            "        )\n"
            "        raise ValueError(msg)\n"
            "    row_colors = row_conditions.map(\n"
            '        lambda condition: condition_palette.get(condition, "#64748b")\n'
            "    )"
        ),
        target_nodeid=(
            "tests/test_enrichment_infographic_followup.py::"
            "test_companion_heatmap_refuses_unlisted_row_conditions[display-unlisted]"
        ),
    ),
    MutantSpec(
        name="M3_skip_unlisted_check_on_dominant_condition_column_2035",
        issue="#2035",
        rel_path=SCRIPT_REL,
        old=(
            "    if unlisted_conditions:\n"
            "        unique_unlisted = list(dict.fromkeys(unlisted_conditions))\n"
            "        msg = (\n"
            '            f"Companion heatmap row condition(s) {unique_unlisted} are not in "\n'
            '            f"TARGET_CONDITIONS {list(TARGET_CONDITIONS)}."\n'
            "        )\n"
            "        raise ValueError(msg)\n"
            "    row_colors = row_conditions.map(lambda condition: condition_palette[condition])"
        ),
        new=(
            '    if condition_color_column == "display_condition" and unlisted_conditions:\n'
            "        unique_unlisted = list(dict.fromkeys(unlisted_conditions))\n"
            "        msg = (\n"
            '            f"Companion heatmap row condition(s) {unique_unlisted} are not in "\n'
            '            f"TARGET_CONDITIONS {list(TARGET_CONDITIONS)}."\n'
            "        )\n"
            "        raise ValueError(msg)\n"
            "    row_colors = row_conditions.map(\n"
            '        lambda condition: condition_palette.get(condition, "#64748b")\n'
            "    )"
        ),
        target_nodeid=(
            "tests/test_enrichment_infographic_followup.py::"
            "test_companion_heatmap_refuses_unlisted_row_conditions[dominant-unknown]"
        ),
    ),
    MutantSpec(
        name="M4_skip_main_paired_top_early_guard_2047",
        issue="#2047",
        rel_path=SCRIPT_REL,
        old=(
            '    if args.condition_source == "paired-top":\n'
            "        validate_companion_heatmap_condition_count(args.top_paired_conditions)"
        ),
        new=(
            '    if False and args.condition_source == "paired-top":\n'
            "        validate_companion_heatmap_condition_count(args.top_paired_conditions)"
        ),
        target_nodeid=(
            "tests/test_enrichment_infographic_followup.py::"
            "test_main_refuses_excess_conditions_before_mkdir_or_loading_or_writing_tables[paired-top]"
        ),
    ),
    MutantSpec(
        name="M5_skip_main_target_conditions_early_guard_2047",
        issue="#2047",
        rel_path=SCRIPT_REL,
        old=(
            "    validate_companion_heatmap_condition_count(len(TARGET_CONDITIONS))\n"
            '    if args.condition_source == "paired-top":'
        ),
        new='    if args.condition_source == "paired-top":',
        target_nodeid=(
            "tests/test_enrichment_infographic_followup.py::"
            "test_main_refuses_excess_conditions_before_mkdir_or_loading_or_writing_tables[infographic]"
        ),
    ),
    MutantSpec(
        name="M6_mkdir_before_condition_count_guard_in_main_2047",
        issue="#2047",
        rel_path=SCRIPT_REL,
        old=(
            "    validate_companion_heatmap_condition_count(len(TARGET_CONDITIONS))\n"
            '    if args.condition_source == "paired-top":\n'
            "        validate_companion_heatmap_condition_count(args.top_paired_conditions)\n"
            "    args.output_dir.mkdir(parents=True, exist_ok=True)"
        ),
        new=(
            "    args.output_dir.mkdir(parents=True, exist_ok=True)\n"
            "    validate_companion_heatmap_condition_count(len(TARGET_CONDITIONS))\n"
            '    if args.condition_source == "paired-top":\n'
            "        validate_companion_heatmap_condition_count(args.top_paired_conditions)"
        ),
        target_nodeid=(
            "tests/test_enrichment_infographic_followup.py::"
            "test_main_refuses_excess_conditions_before_mkdir_or_loading_or_writing_tables[paired-top]"
        ),
    ),
    MutantSpec(
        name="M7_widen_condition_count_limit_to_nine_2047",
        issue="#2047",
        rel_path=SCRIPT_REL,
        old="    if condition_count > max_conditions:",
        new="    if condition_count > max_conditions + 1:",
        target_nodeid=(
            "tests/test_enrichment_infographic_followup.py::"
            "test_main_refuses_excess_conditions_before_mkdir_or_loading_or_writing_tables[paired-top]"
        ),
    ),
    MutantSpec(
        name="M8_remove_heatmap_backstop_call_2047",
        issue="#2047",
        rel_path=SCRIPT_REL,
        old=(
            "    validate_companion_heatmap_condition_count(len(TARGET_CONDITIONS))\n"
            "    condition_palette = {\n"
            "        condition: COMPANION_HEATMAP_PALETTE_COLORS[index]\n"
            "        for index, condition in enumerate(TARGET_CONDITIONS)\n"
            "    }"
        ),
        new=(
            "    condition_palette = {\n"
            "        condition: COMPANION_HEATMAP_PALETTE_COLORS[\n"
            "            index % len(COMPANION_HEATMAP_PALETTE_COLORS)\n"
            "        ]\n"
            "        for index, condition in enumerate(TARGET_CONDITIONS)\n"
            "    }"
        ),
        target_nodeid=(
            "tests/test_enrichment_infographic_followup.py::"
            "test_direct_save_heatmap_backstop_refuses_ninth_condition_with_value_error"
        ),
    ),
)


class _CallPhaseInspector:
    """Pytest plugin recording item origin, phase outcomes, and exception classes."""

    def __init__(self, expected_test_file: Path, expected_script: Path, expected_name: str) -> None:
        self.expected_test_file = expected_test_file
        self.expected_script = expected_script
        self.expected_name = expected_name
        self.origin_ok = False
        self.collected_nodeids: list[str] = []
        self.phase_outcomes: dict[str, str] = {}
        self.call_exception: dict[str, object] | None = None
        self.non_call_exceptions: list[str] = []

    def pytest_collection_modifyitems(self, session: Any, config: Any, items: list[Any]) -> None:
        del session, config
        self.collected_nodeids = [str(item.nodeid) for item in items]
        if len(items) != 1:
            self.origin_ok = False
            return
        item = items[0]
        item_path = Path(getattr(item, "path", getattr(item, "fspath", ""))).resolve()
        mod = getattr(item, "module", None)
        mod_file = Path(getattr(mod, "__file__", "") or "").resolve() if mod is not None else None
        mod_script = (
            Path(getattr(mod, "MODULE_PATH", "") or "").resolve() if mod is not None else None
        )
        followup_mod = getattr(mod, "followup", None) if mod is not None else None
        followup_file = (
            Path(getattr(followup_mod, "__file__", "") or "").resolve()
            if followup_mod is not None
            else None
        )
        self.origin_ok = (
            item.name == self.expected_name
            and item_path == self.expected_test_file
            and mod_file == self.expected_test_file
            and mod_script == self.expected_script
            and followup_file == self.expected_script
        )

    def pytest_runtest_logreport(self, report: Any) -> None:
        self.phase_outcomes[str(report.when)] = str(report.outcome)

    def pytest_exception_interact(self, node: Any, call: Any, report: Any) -> None:
        del node, report
        excinfo = getattr(call, "excinfo", None)
        if excinfo is None:
            return
        when = str(getattr(call, "when", "unknown"))
        exc_type = str(excinfo.type.__name__)
        is_assertion = bool(issubclass(excinfo.type, AssertionError))
        exc_value = str(excinfo.value)
        tb_path = ""
        tb_lineno = 0
        if excinfo.traceback:
            last_entry = excinfo.traceback[-1]
            tb_path = str(Path(str(last_entry.path)).resolve())
            tb_lineno = int(last_entry.lineno) + 1
        if when == "call":
            self.call_exception = {
                "exc_type": exc_type,
                "is_assertion_error": is_assertion,
                "exc_value": exc_value,
                "tb_path": tb_path,
                "tb_lineno": tb_lineno,
            }
        else:
            self.non_call_exceptions.append(f"{when}:{exc_type}:{exc_value}")


def _populate_scratch(scratch_root: Path) -> None:
    script_dst = scratch_root / SCRIPT_REL
    script_dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(ROOT / SCRIPT_REL, script_dst)
    test_dst = scratch_root / TEST_REL
    test_dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(ROOT / TEST_REL, test_dst)


def _worker(
    scratch_root_str: str,
    nodeid: str,
    out_file_str: str,
    report_file_str: str,
    junit_xml_str: str,
) -> None:
    """Run pytest rooted at scratch_root with importlib isolation and record phase report."""
    scratch_root = Path(scratch_root_str).resolve()
    prev_cwd = Path.cwd()
    saved_sys_modules = dict(sys.modules)
    saved_sys_path = list(sys.path)
    saved_meta_path = list(sys.meta_path)
    os.environ["GIT_CONFIG_GLOBAL"] = os.devnull
    os.environ["GIT_CONFIG_NOSYSTEM"] = "1"
    rel_test_path, _, test_name = nodeid.partition("::")
    test_file = (scratch_root / rel_test_path).resolve()
    expected_script = (scratch_root / SCRIPT_REL).resolve()
    inspector = _CallPhaseInspector(test_file, expected_script, test_name)
    buf = io.StringIO()
    try:
        for mod_name in list(sys.modules):
            if (
                mod_name == "tests"
                or mod_name.startswith("tests.")
                or "test_enrichment_infographic_followup" in mod_name
                or mod_name == "enrichment_infographic_followup"
            ):
                sys.modules.pop(mod_name, None)
        excluded_paths = {ROOT.resolve(), (ROOT / "tests").resolve()}
        sys.path[:] = [
            str(scratch_root),
            str(ROOT / "src"),
            str(scratch_root / "tests"),
        ] + [p for p in saved_sys_path if not p or Path(p).resolve() not in excluded_paths]
        sys.meta_path[:] = [
            hook for hook in saved_meta_path if type(hook).__name__ != "AssertionRewritingHook"
        ]
        os.chdir(scratch_root)
        importlib.invalidate_caches()
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
            rc = int(
                pytest.main(
                    [
                        "-q",
                        "--import-mode=importlib",
                        f"--rootdir={scratch_root}",
                        "-o",
                        "addopts=",
                        "-o",
                        f"cache_dir={scratch_root / '.pytest_cache'}",
                        f"--junitxml={junit_xml_str}",
                        f"{test_file}::{test_name}",
                    ],
                    plugins=[inspector],
                )
            )
            if inspector.origin_ok:
                print("PYTEST_ORIGIN_OK", flush=True)
    finally:
        for mod_name in list(sys.modules):
            if (
                mod_name == "tests"
                or mod_name.startswith("tests.")
                or "test_enrichment_infographic_followup" in mod_name
                or mod_name == "enrichment_infographic_followup"
            ):
                sys.modules.pop(mod_name, None)
        for mod_name, mod_obj in saved_sys_modules.items():
            if (
                mod_name == "tests"
                or mod_name.startswith("tests.")
                or "test_enrichment_infographic_followup" in mod_name
                or mod_name == "enrichment_infographic_followup"
            ):
                sys.modules[mod_name] = mod_obj
        sys.meta_path[:] = saved_meta_path
        sys.path[:] = saved_sys_path
        os.chdir(prev_cwd)
    output_text = buf.getvalue()
    Path(out_file_str).write_text(f"RC={rc}\n{output_text}", encoding="utf-8")
    report_payload = {
        "rc": rc,
        "output": output_text,
        "origin_ok": inspector.origin_ok,
        "collected_nodeids": inspector.collected_nodeids,
        "phase_outcomes": inspector.phase_outcomes,
        "call_exception": inspector.call_exception,
        "non_call_exceptions": inspector.non_call_exceptions,
        "junit_xml_path": junit_xml_str,
    }
    Path(report_file_str).write_text(json.dumps(report_payload, indent=2) + "\n", encoding="utf-8")


def _run_targeted_pytest(scratch_root: Path, nodeid: str) -> TargetedPytestResult:
    """Execute nodeid in scratch_root and verify call-phase AssertionError semantics."""
    out_file = scratch_root / "pytest_out.txt"
    report_file = scratch_root / "pytest_report.json"
    junit_xml_file = scratch_root / "junit.xml"
    _worker(
        str(scratch_root),
        nodeid,
        str(out_file),
        str(report_file),
        str(junit_xml_file),
    )
    if not out_file.exists() or not report_file.exists():
        raise RuntimeError("pytest report artifacts were not created")
    data = json.loads(report_file.read_text(encoding="utf-8"))
    rc = int(data["rc"])
    output = str(data["output"])
    origin_ok = bool(data["origin_ok"]) and "PYTEST_ORIGIN_OK" in output
    collected_nodeids = tuple(str(x) for x in data.get("collected_nodeids", []))
    phase_outcomes: dict[str, str] = {
        str(k): str(v) for k, v in dict(data.get("phase_outcomes", {})).items()
    }
    setup_outcome = phase_outcomes.get("setup")
    call_outcome = phase_outcomes.get("call")
    teardown_outcome = phase_outcomes.get("teardown")
    call_exc = data.get("call_exception")
    call_exc_type = str(call_exc["exc_type"]) if isinstance(call_exc, dict) else None
    call_is_assertion_error = (
        bool(call_exc["is_assertion_error"]) if isinstance(call_exc, dict) else False
    )
    call_exc_value = str(call_exc["exc_value"]) if isinstance(call_exc, dict) else None
    call_tb_path = str(call_exc["tb_path"]) if isinstance(call_exc, dict) else None
    call_tb_lineno = int(call_exc["tb_lineno"]) if isinstance(call_exc, dict) else None
    non_call_exceptions = tuple(str(x) for x in data.get("non_call_exceptions", []))
    rel_test_path, _, _ = nodeid.partition("::")
    expected_test_file = str((scratch_root / rel_test_path).resolve())
    junit_text = junit_xml_file.read_text(encoding="utf-8") if junit_xml_file.is_file() else ""
    junit_call_failure_ok = (
        'errors="0"' in junit_text
        and 'failures="1"' in junit_text
        and "<failure " in junit_text
        and "<error " not in junit_text
    )
    killed_by_call_assertion = (
        rc == 1
        and origin_ok
        and len(collected_nodeids) == 1
        and setup_outcome == "passed"
        and call_outcome == "failed"
        and teardown_outcome == "passed"
        and not non_call_exceptions
        and call_is_assertion_error
        and call_exc_type == "AssertionError"
        and call_tb_path == expected_test_file
        and junit_call_failure_ok
    )
    return TargetedPytestResult(
        rc=rc,
        output=output,
        origin_ok=origin_ok,
        collected_nodeids=collected_nodeids,
        setup_outcome=setup_outcome,
        call_outcome=call_outcome,
        teardown_outcome=teardown_outcome,
        call_exc_type=call_exc_type,
        call_is_assertion_error=call_is_assertion_error,
        call_exc_value=call_exc_value,
        call_tb_path=call_tb_path,
        call_tb_lineno=call_tb_lineno,
        non_call_exceptions=non_call_exceptions,
        junit_xml_path=str(junit_xml_file),
        killed_by_call_assertion=killed_by_call_assertion,
    )


def _extract_failure_detail(result: TargetedPytestResult) -> str:
    first_line = (result.call_exc_value or "").splitlines()[0] if result.call_exc_value else ""
    location = f"{TEST_REL}:{result.call_tb_lineno}" if result.call_tb_lineno else str(TEST_REL)
    if first_line:
        return (
            f"KILLED by call-phase AssertionError at {location} with PYTEST_ORIGIN_OK: {first_line}"
        )
    return f"KILLED by call-phase AssertionError at {location} with PYTEST_ORIGIN_OK"


def _write_verdicts(run_dir: Path, scratch_parent: Path, verdicts: list[dict[str, object]]) -> Path:
    verdicts_payload = json.dumps(verdicts, indent=2) + "\n"
    run_verdicts = run_dir / "verdicts.json"
    parent_verdicts = scratch_parent / "verdicts.json"
    run_verdicts.write_text(verdicts_payload, encoding="utf-8")
    parent_verdicts.write_text(verdicts_payload, encoding="utf-8")
    return run_verdicts


def verify_mutants(scratch_dir: Path, *, canary: bool = False) -> tuple[int, Path]:
    """Run control and selected mutants in persistent subdirectories of scratch_dir."""
    scratch_parent = scratch_dir.resolve()
    if scratch_parent.is_relative_to(ROOT):
        raise RuntimeError("--scratch-dir must be outside the repository root")
    scratch_parent.mkdir(parents=True, exist_ok=True)

    originals = {rel: (ROOT / rel).read_bytes() for rel in (SCRIPT_REL, TEST_REL)}
    selected = MUTANTS[:1] if canary else MUTANTS
    run_dir = Path(
        tempfile.mkdtemp(
            prefix="verify_2035_2047_canary_" if canary else "verify_2035_2047_all_",
            dir=scratch_parent,
        )
    )
    verdicts: list[dict[str, object]] = []

    control_nodeid = selected[0].target_nodeid
    control_dir = run_dir / "control"
    control_dir.mkdir(parents=True, exist_ok=True)
    _populate_scratch(control_dir)
    ctrl = _run_targeted_pytest(control_dir, control_nodeid)
    control_passed = (
        ctrl.rc == 0
        and ctrl.origin_ok
        and len(ctrl.collected_nodeids) == 1
        and ctrl.setup_outcome == "passed"
        and ctrl.call_outcome == "passed"
        and ctrl.teardown_outcome == "passed"
        and not ctrl.non_call_exceptions
    )
    if not control_passed:
        verdicts.append(
            {
                "case": "control",
                "name": "control",
                "nodeid": control_nodeid,
                "rc": ctrl.rc,
                "killed": False,
                "verdict": "CONTROL FAILED",
                "matched_marker": "",
                "setup_outcome": ctrl.setup_outcome,
                "call_outcome": ctrl.call_outcome,
                "teardown_outcome": ctrl.teardown_outcome,
                "detail": (
                    f"CONTROL FAILED: rc={ctrl.rc} origin_ok={ctrl.origin_ok} "
                    f"setup={ctrl.setup_outcome} call={ctrl.call_outcome} "
                    f"teardown={ctrl.teardown_outcome}"
                ),
                "log_file": str(control_dir / "pytest_out.txt"),
                "junit_xml": ctrl.junit_xml_path,
            }
        )
        _write_verdicts(run_dir, scratch_parent, verdicts)
        print(f"CONTROL FAILED: rc={ctrl.rc}\n{ctrl.output}")
        return 1, run_dir

    verdicts.append(
        {
            "case": "control",
            "name": "control",
            "nodeid": control_nodeid,
            "rc": ctrl.rc,
            "killed": False,
            "verdict": "CONTROL OK",
            "matched_marker": "PYTEST_ORIGIN_OK",
            "setup_outcome": ctrl.setup_outcome,
            "call_outcome": ctrl.call_outcome,
            "teardown_outcome": ctrl.teardown_outcome,
            "detail": (
                "CONTROL OK: target test passed setup, call, and teardown on "
                "unmutated scratch copy with verified PYTEST_ORIGIN_OK"
            ),
            "log_file": str(control_dir / "pytest_out.txt"),
            "junit_xml": ctrl.junit_xml_path,
        }
    )
    print("CONTROL OK (target test passed on unmutated scratch copy)")

    for mutant in selected:
        case_dir = run_dir / mutant.name
        case_dir.mkdir(parents=True, exist_ok=True)
        _populate_scratch(case_dir)
        target = case_dir / mutant.rel_path
        src_text = target.read_text(encoding="utf-8")
        if src_text.count(mutant.old) != 1:
            raise RuntimeError(
                f"Mutant {mutant.name} target snippet count={src_text.count(mutant.old)} != 1"
            )
        target.write_text(src_text.replace(mutant.old, mutant.new, 1), encoding="utf-8")
        mut_res = _run_targeted_pytest(case_dir, mutant.target_nodeid)
        if not mut_res.killed_by_call_assertion:
            reject_detail = (
                f"REJECTED [{mutant.name}]: rc={mut_res.rc} origin_ok={mut_res.origin_ok} "
                f"setup={mut_res.setup_outcome} call={mut_res.call_outcome} "
                f"teardown={mut_res.teardown_outcome} call_exc_type={mut_res.call_exc_type} "
                f"non_call_exceptions={list(mut_res.non_call_exceptions)}"
            )
            verdicts.append(
                {
                    "case": mutant.name,
                    "name": mutant.name,
                    "issue": mutant.issue,
                    "nodeid": mutant.target_nodeid,
                    "rc": mut_res.rc,
                    "killed": False,
                    "verdict": "REJECTED",
                    "matched_marker": "",
                    "setup_outcome": mut_res.setup_outcome,
                    "call_outcome": mut_res.call_outcome,
                    "teardown_outcome": mut_res.teardown_outcome,
                    "call_exc_type": mut_res.call_exc_type,
                    "call_tb_lineno": mut_res.call_tb_lineno,
                    "non_call_exceptions": list(mut_res.non_call_exceptions),
                    "detail": reject_detail,
                    "log_file": str(case_dir / "pytest_out.txt"),
                    "junit_xml": mut_res.junit_xml_path,
                }
            )
            _write_verdicts(run_dir, scratch_parent, verdicts)
            print(f"{reject_detail}\n{mut_res.output}")
            return 1, run_dir
        detail_str = _extract_failure_detail(mut_res)
        if not detail_str.strip():
            raise RuntimeError(f"Empty detail for {mutant.name}")
        verdicts.append(
            {
                "case": mutant.name,
                "name": mutant.name,
                "issue": mutant.issue,
                "nodeid": mutant.target_nodeid,
                "rc": mut_res.rc,
                "killed": True,
                "verdict": "KILLED",
                "matched_marker": "AssertionError",
                "setup_outcome": mut_res.setup_outcome,
                "call_outcome": mut_res.call_outcome,
                "teardown_outcome": mut_res.teardown_outcome,
                "call_exc_type": mut_res.call_exc_type,
                "call_tb_lineno": mut_res.call_tb_lineno,
                "non_call_exceptions": list(mut_res.non_call_exceptions),
                "detail": detail_str,
                "log_file": str(case_dir / "pytest_out.txt"),
                "junit_xml": mut_res.junit_xml_path,
            }
        )
        print(
            f"KILLED {mutant.name} ({mutant.issue}): "
            f"rc={mut_res.rc} AssertionError PYTEST_ORIGIN_OK"
        )

    for rel, orig_bytes in originals.items():
        current = (ROOT / rel).read_bytes()
        if hashlib.sha256(current).hexdigest() != hashlib.sha256(orig_bytes).hexdigest():
            raise RuntimeError(f"Worktree file modified during mutant verification: {rel}")

    verdicts_file = _write_verdicts(run_dir, scratch_parent, verdicts)
    if not verdicts_file.is_file():
        raise RuntimeError(f"Expected verdicts file does not exist: {verdicts_file}")
    print(
        f"ALL {len(selected)}/{len(selected)} #2035/#2047 MUTANTS KILLED CLEANLY; "
        f"worktree targets byte-identical; evidence: {run_dir}"
    )
    return 0, run_dir


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--canary",
        action="store_true",
        help="Run only the first mutant as a canary.",
    )
    parser.add_argument(
        "--scratch-dir",
        type=Path,
        required=True,
        help="Persistent directory outside the repository root for scratch runs and verdicts.json.",
    )
    args = parser.parse_args(argv)
    rc, _ = verify_mutants(args.scratch_dir, canary=args.canary)
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
