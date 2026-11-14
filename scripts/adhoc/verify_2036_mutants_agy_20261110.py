#!/usr/bin/env python3
"""Scratch-copy mutation verifier for #1997, #2036, #2037, and #2038.

Purpose:
  Verify that tests/test_cellulose_tufte_readme_scope_2036.py passes on an
  unmutated scratch copy (CONTROL) and fails in its pytest call phase with an
  actual AssertionError when each of the 8 scope/disposition mutations in
  tests/fixtures/cellulose_tufte_historical_sources/README.md is injected into
  an isolated scratch copy outside the worktree. Collection, setup, and
  teardown failures or non-AssertionError exceptions are rejected even if their
  messages contain the word "AssertionError".

Inputs:
  - tests/fixtures/cellulose_tufte_historical_sources/README.md
  - tests/test_cellulose_tufte_readme_scope_2036.py

Outputs:
  - verdicts.json in the persistent scratch run directory and parent
    --scratch-dir with non-empty detail strings and call-phase metadata for
    CONTROL and every mutant.
  - Exit 0 if CONTROL passes, all selected mutants die in the call phase by
    AssertionError with verified scratch origin, and canonical worktree files
    remain byte-identical; non-zero otherwise.

Issues: #1997, #2036, #2037, #2038
Date: 2026-11-10
"""

from __future__ import annotations

import argparse
import contextlib
import importlib
import io
import json
import os
import shutil
import sys
import tempfile
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pytest

ROOT = Path(__file__).resolve().parents[2]
README_REL = "tests/fixtures/cellulose_tufte_historical_sources/README.md"
TEST_REL = "tests/test_cellulose_tufte_readme_scope_2036.py"
TARGET_NODEID = (
    "tests/test_cellulose_tufte_readme_scope_2036.py::"
    "test_historical_sources_readme_scope_and_unpinned_disposition"
)


@dataclass(frozen=True)
class Mutant:
    """Single source-text mutation against a relative file path."""

    name: str
    rel_path: str
    old: str
    new: str
    nodeid: str
    expected_marker: str


@dataclass(frozen=True)
class TargetedPytestResult:
    """Structured outcome from running a single pytest node in an isolated scratch copy."""

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


MUTANTS: tuple[Mutant, ...] = (
    Mutant(
        name="M1_restore_repeated_scope_sentence_2038",
        rel_path=README_REL,
        old=(
            "The table below lists what `test_cellulose_tufte_statistics_sources.py` binds to\n"
            "identified historical sources; the other test files are summarised after it.\n"
            "They do not pin every numeric token or every repeated occurrence,"
        ),
        new=(
            "The table below lists what `test_cellulose_tufte_statistics_sources.py` binds to\n"
            "identified historical sources; the other test files are summarised after it. These\n"
            "tests bind selected displayed values and labels to identified historical sources. "
            "They do not pin every numeric token or every repeated occurrence,"
        ),
        nodeid=TARGET_NODEID,
        expected_marker="AssertionError",
    ),
    Mutant(
        name="M2_revert_overlap_summary_to_stale_2036",
        rel_path=README_REL,
        old=(
            "`test_cellulose_tufte_overlap_sources.py` pins selected set sizes, direction "
            "labels, percentages, distances and Plate III's nominal p-values for Plates I to IV "
            "and figures 7 to 9; it is summarised here, not tabulated."
        ),
        new=(
            "`test_cellulose_tufte_overlap_sources.py` pins selected set sizes and direction "
            "labels\nfor Plates I to IV and figures 7 to 9; it is summarised here, not tabulated."
        ),
        nodeid=TARGET_NODEID,
        expected_marker="AssertionError",
    ),
    Mutant(
        name="M3_drop_plate_three_nominal_p_from_overlap_summary_2036",
        rel_path=README_REL,
        old="percentages, distances and Plate III's nominal p-values",
        new="percentages and distances",
        nodeid=TARGET_NODEID,
        expected_marker="AssertionError",
    ),
    Mutant(
        name="M4_revert_mutation_summary_to_pins_above_2037",
        rel_path=README_REL,
        old=(
            "`test_cellulose_tufte_source_pin_mutations.py` checks that every tabulated pin, "
            "and one overlap pin, fails on targeted mutations."
        ),
        new=(
            "`test_cellulose_tufte_source_pin_mutations.py` checks that the pins\n"
            "above fail on targeted mutations."
        ),
        nodeid=TARGET_NODEID,
        expected_marker="AssertionError",
    ),
    Mutant(
        name="M5_revert_unpinned_tail_to_broader_source_pins_1997",
        rel_path=README_REL,
        old=(
            "The ten sourced-but-unpinned tokens from #1596/#1997 ('10 PD\u2229CRC', '= 40 taxa', "
            "'46 PD', and '27 shared (57.4%)' in fig7; '17 &lt; 27 (PD\u2013CRC)' in fig9; "
            "'17 of 32' in Plate I; '45%\u201357%' in Plate III; '50.4% shared' in Plate IV; "
            "'n = 2,335 pairs' in Plate V; 'Both FA (58.2% up' in Plate VI) have explicit "
            "disposition: won't pin (#1997), because the plates were retired from default "
            "generation by #1693 (#1572)."
        ),
        new="Broader source pins remain a separate task.",
        nodeid=TARGET_NODEID,
        expected_marker="AssertionError",
    ),
    Mutant(
        name="M6_omit_fig7_token_from_unpinned_list_1997",
        rel_path=README_REL,
        old="'27 shared (57.4%)'",
        new="'omitted fig7 token'",
        nodeid=TARGET_NODEID,
        expected_marker="AssertionError",
    ),
    Mutant(
        name="M7_omit_fig9_and_plate_tokens_from_unpinned_list_1997",
        rel_path=README_REL,
        old=(
            "'17 &lt; 27 (PD\u2013CRC)' in fig9; '17 of 32' in Plate I; "
            "'45%\u201357%' in Plate III; "
            "'50.4% shared' in Plate IV; 'n = 2,335 pairs' in Plate V; "
            "'Both FA (58.2% up' in Plate VI"
        ),
        new="'omitted remaining tokens'",
        nodeid=TARGET_NODEID,
        expected_marker="AssertionError",
    ),
    Mutant(
        name="M8_omit_wont_pin_retirement_disposition_1997",
        rel_path=README_REL,
        old=(
            "won't pin (#1997), because the plates were retired from default generation "
            "by #1693 (#1572)."
        ),
        new="deferred to a follow-up issue.",
        nodeid=TARGET_NODEID,
        expected_marker="AssertionError",
    ),
)


class _CallPhaseInspector:
    """Pytest plugin recording item origin, phase outcomes, and exception classes."""

    def __init__(self, expected_test_file: Path, expected_readme: Path, expected_name: str) -> None:
        self.expected_test_file = expected_test_file
        self.expected_readme = expected_readme
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
        mod_readme = (
            Path(getattr(mod, "README_PATH", "") or "").resolve() if mod is not None else None
        )
        self.origin_ok = (
            item.name == self.expected_name
            and item_path == self.expected_test_file
            and mod_file == self.expected_test_file
            and mod_readme == self.expected_readme
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
    readme_dst = scratch_root / README_REL
    readme_dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(ROOT / README_REL, readme_dst)
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
    expected_readme = (scratch_root / README_REL).resolve()
    inspector = _CallPhaseInspector(test_file, expected_readme, test_name)
    buf = io.StringIO()
    try:
        for mod_name in list(sys.modules):
            if (
                mod_name == "tests"
                or mod_name.startswith("tests.")
                or "test_cellulose_tufte_readme_scope_2036" in mod_name
            ):
                sys.modules.pop(mod_name, None)
        excluded_paths = {ROOT.resolve(), (ROOT / "tests").resolve()}
        sys.path[:] = [str(scratch_root), str(scratch_root / "tests")] + [
            p for p in saved_sys_path if not p or Path(p).resolve() not in excluded_paths
        ]
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
        sys.meta_path[:] = saved_meta_path
        sys.path[:] = saved_sys_path
        sys.modules.clear()
        sys.modules.update(saved_sys_modules)
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
    location = f"{TEST_REL}:{result.call_tb_lineno}" if result.call_tb_lineno else TEST_REL
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

    originals = {rel: (ROOT / rel).read_bytes() for rel in (README_REL, TEST_REL)}
    selected = MUTANTS[:1] if canary else MUTANTS
    run_dir = Path(
        tempfile.mkdtemp(
            prefix="verify_2036_canary_" if canary else "verify_2036_all_",
            dir=scratch_parent,
        )
    )
    verdicts: list[dict[str, object]] = []

    try:
        control_scratch = run_dir / "control"
        control_scratch.mkdir(parents=True, exist_ok=True)
        _populate_scratch(control_scratch)
        for nodeid in dict.fromkeys(m.nodeid for m in selected):
            control_res = _run_targeted_pytest(control_scratch, nodeid)
            control_passed = (
                control_res.rc == 0
                and control_res.origin_ok
                and len(control_res.collected_nodeids) == 1
                and control_res.setup_outcome == "passed"
                and control_res.call_outcome == "passed"
                and control_res.teardown_outcome == "passed"
                and not control_res.non_call_exceptions
            )
            if not control_passed:
                verdicts.append(
                    {
                        "case": "control",
                        "name": "control",
                        "nodeid": nodeid,
                        "rc": control_res.rc,
                        "killed": False,
                        "verdict": "CONTROL FAILED",
                        "matched_marker": "",
                        "setup_outcome": control_res.setup_outcome,
                        "call_outcome": control_res.call_outcome,
                        "teardown_outcome": control_res.teardown_outcome,
                        "detail": (
                            f"CONTROL FAILED: rc={control_res.rc} "
                            f"setup={control_res.setup_outcome} "
                            f"call={control_res.call_outcome} "
                            f"teardown={control_res.teardown_outcome}"
                        ),
                        "log_file": str(control_scratch / "pytest_out.txt"),
                        "junit_xml": control_res.junit_xml_path,
                    }
                )
                _write_verdicts(run_dir, scratch_parent, verdicts)
                print(
                    f"CONTROL FAILED for {nodeid} (rc={control_res.rc}):\n{control_res.output}",
                    file=sys.stderr,
                )
                return 1, run_dir
        verdicts.append(
            {
                "case": "control",
                "name": "control",
                "nodeid": TARGET_NODEID,
                "rc": 0,
                "killed": False,
                "verdict": "CONTROL OK",
                "matched_marker": "PYTEST_ORIGIN_OK",
                "setup_outcome": "passed",
                "call_outcome": "passed",
                "teardown_outcome": "passed",
                "detail": (
                    "CONTROL OK: target test passed setup, call, and teardown on "
                    "unmutated scratch copy with verified PYTEST_ORIGIN_OK"
                ),
                "log_file": str(control_scratch / "pytest_out.txt"),
                "junit_xml": str(control_scratch / "junit.xml"),
            }
        )
        print("CONTROL OK (target test passed on unmutated scratch copy)")

        for mutant in selected:
            scratch = run_dir / mutant.name
            scratch.mkdir(parents=True, exist_ok=True)
            _populate_scratch(scratch)
            target = scratch / mutant.rel_path
            text = target.read_text(encoding="utf-8")
            count = text.count(mutant.old)
            if count != 1:
                print(
                    f"FAIL [{mutant.name}]: expected 1 occurrence in "
                    f"{mutant.rel_path}, found {count}",
                    file=sys.stderr,
                )
                _write_verdicts(run_dir, scratch_parent, verdicts)
                return 1, run_dir
            target.write_text(text.replace(mutant.old, mutant.new, 1), encoding="utf-8")
            res = _run_targeted_pytest(scratch, mutant.nodeid)
            if not res.killed_by_call_assertion or res.call_exc_type != mutant.expected_marker:
                rejection_detail = (
                    f"REJECTED [{mutant.name}]: rc={res.rc} origin_ok={res.origin_ok} "
                    f"setup={res.setup_outcome} call={res.call_outcome} "
                    f"teardown={res.teardown_outcome} call_exc_type={res.call_exc_type} "
                    f"non_call_exceptions={list(res.non_call_exceptions)}"
                )
                verdicts.append(
                    {
                        "case": mutant.name,
                        "name": mutant.name,
                        "nodeid": mutant.nodeid,
                        "rc": res.rc,
                        "killed": False,
                        "verdict": "REJECTED",
                        "matched_marker": "",
                        "setup_outcome": res.setup_outcome,
                        "call_outcome": res.call_outcome,
                        "teardown_outcome": res.teardown_outcome,
                        "call_exc_type": res.call_exc_type,
                        "call_tb_lineno": res.call_tb_lineno,
                        "non_call_exceptions": list(res.non_call_exceptions),
                        "detail": rejection_detail,
                        "log_file": str(scratch / "pytest_out.txt"),
                        "junit_xml": res.junit_xml_path,
                    }
                )
                _write_verdicts(run_dir, scratch_parent, verdicts)
                print(f"{rejection_detail}\n{res.output}", file=sys.stderr)
                return 1, run_dir
            detail = _extract_failure_detail(res)
            if not detail.strip():
                raise RuntimeError(f"Empty detail for mutant {mutant.name}")
            verdicts.append(
                {
                    "case": mutant.name,
                    "name": mutant.name,
                    "nodeid": mutant.nodeid,
                    "rc": res.rc,
                    "killed": True,
                    "verdict": "KILLED",
                    "matched_marker": mutant.expected_marker,
                    "setup_outcome": res.setup_outcome,
                    "call_outcome": res.call_outcome,
                    "teardown_outcome": res.teardown_outcome,
                    "call_exc_type": res.call_exc_type,
                    "call_tb_lineno": res.call_tb_lineno,
                    "non_call_exceptions": list(res.non_call_exceptions),
                    "detail": detail,
                    "log_file": str(scratch / "pytest_out.txt"),
                    "junit_xml": res.junit_xml_path,
                }
            )
            print(f"KILLED {mutant.name}: rc=1 {mutant.expected_marker} PYTEST_ORIGIN_OK")
        run_verdicts = _write_verdicts(run_dir, scratch_parent, verdicts)
        if not run_verdicts.is_file():
            raise RuntimeError(f"Expected verdicts file missing: {run_verdicts}")
    finally:
        for rel, content in originals.items():
            if (ROOT / rel).read_bytes() != content:
                raise RuntimeError(f"Canonical worktree file changed during verification: {rel}")

    print(
        f"ALL {len(selected)}/{len(selected)} #1997/#2036/#2037/#2038 MUTANTS KILLED CLEANLY; "
        f"worktree targets byte-identical; evidence: {run_dir}"
    )
    return 0, run_dir


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--canary",
        action="store_true",
        help="Run control and the first mutant only",
    )
    parser.add_argument(
        "--scratch-dir",
        type=Path,
        required=True,
        help="Parent directory outside the repository for persistent scratch copies and verdicts",
    )
    args = parser.parse_args(argv)
    rc, _ = verify_mutants(args.scratch_dir, canary=bool(args.canary))
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
