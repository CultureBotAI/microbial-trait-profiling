"""Scratch-isolated mutant verifier for Issue #997 should-fix items.

Purpose:
    Verify that `tests/test_study_split_recheck_and_untestable_997.py` kills mutations
    to the pre-publication provenance recheck in `src/mtp/cellulose_decr/study_split_run.py`
    and the untestable-family selection status handling in
    `src/mtp/cellulose_decr/study_split.py` via call-phase `AssertionError`.

Inputs:
    - Worktree copies of `src/mtp/`, `methods/`, `schema/`, `workflows/`,
      and `tests/test_study_split_recheck_and_untestable_997.py`.

Outputs:
    - Per-run JUnit XML files (`junit_<label>.xml`) and `verdicts.json` inside a
      persistent subdirectory under `--scratch-dir`.

Issue: #997
Date: 2026-11-11
"""

from __future__ import annotations

import argparse
import importlib
import io
import json
import os
import shutil
import sys
import tempfile
from contextlib import redirect_stderr, redirect_stdout
from dataclasses import dataclass
from pathlib import Path

import pytest


@dataclass(frozen=True)
class MutantSpec:
    mutant_id: str
    description: str
    target_relpath: str
    needle: str
    replacement: str
    pytest_nodeid: str


MUTANTS: tuple[MutantSpec, ...] = (
    MutantSpec(
        mutant_id="M1_remove_pre_publication_recheck",
        description=(
            "Remove `if _provenance(paths, protocol, verified) != provenance:` guard "
            "before publication in study_split_run.execute."
        ),
        target_relpath="src/mtp/cellulose_decr/study_split_run.py",
        needle=(
            "    documents = _documents(verified, paths, protocol, timestamp)\n"
            "    if _provenance(paths, protocol, verified) != provenance:\n"
            '        raise ValueError("Study-split inputs or implementation '
            'changed during execution")'
        ),
        replacement="    documents = _documents(verified, paths, protocol, timestamp)",
        pytest_nodeid=(
            "tests/test_study_split_recheck_and_untestable_997.py::"
            "test_execute_refuses_publication_when_code_sha256_drifts_during_execution"
        ),
    ),
    MutantSpec(
        mutant_id="M2_wrong_recheck_error_message",
        description=(
            "Alter the ValueError message raised by the pre-publication provenance recheck "
            "in study_split_run.execute."
        ),
        target_relpath="src/mtp/cellulose_decr/study_split_run.py",
        needle=(
            "    if _provenance(paths, protocol, verified) != provenance:\n"
            '        raise ValueError("Study-split inputs or implementation '
            'changed during execution")'
        ),
        replacement=(
            "    if _provenance(paths, protocol, verified) != provenance:\n"
            '        raise ValueError("unexpected recheck error message")'
        ),
        pytest_nodeid=(
            "tests/test_study_split_recheck_and_untestable_997.py::"
            "test_execute_refuses_publication_when_code_sha256_drifts_during_execution"
        ),
    ),
    MutantSpec(
        mutant_id="M3_provenance_omits_recorded_inputs",
        description=(
            "Omit `_recorded_inputs(paths, verified)` from `study_split_run._provenance` "
            "so mid-execution input file drift is ignored."
        ),
        target_relpath="src/mtp/cellulose_decr/study_split_run.py",
        needle=(
            "def _provenance(\n"
            "    paths: dict[str, Path], protocol: Protocol, "
            "verified: VerifiedSources | None = None\n"
            ") -> dict[str, Any]:\n"
            "    return {**_execution_identity(protocol), "
            '"inputs": _recorded_inputs(paths, verified)}'
        ),
        replacement=(
            "def _provenance(\n"
            "    paths: dict[str, Path], protocol: Protocol, "
            "verified: VerifiedSources | None = None\n"
            ") -> dict[str, Any]:\n"
            '    return {**_execution_identity(protocol), "inputs": {}}'
        ),
        pytest_nodeid=(
            "tests/test_study_split_recheck_and_untestable_997.py::"
            "test_execute_refuses_publication_when_consumed_input_drifts_during_execution"
        ),
    ),
    MutantSpec(
        mutant_id="M4_execution_identity_omits_implementations",
        description=(
            "Omit `IMPLEMENTATIONS` hashing in `study_split_run._execution_identity` "
            "so implementation file drift during `_documents` is ignored."
        ),
        target_relpath="src/mtp/cellulose_decr/study_split_run.py",
        needle=(
            '"implementation": {name: file_sha256(REPO_ROOT / name) for name in IMPLEMENTATIONS},'
        ),
        replacement='"implementation": {},',
        pytest_nodeid=(
            "tests/test_study_split_recheck_and_untestable_997.py::"
            "test_execute_refuses_publication_when_implementation_file_drifts_during_execution"
        ),
    ),
    MutantSpec(
        mutant_id="M5_remove_no_disjoint_triple_override",
        description=(
            "Remove `if status == scanner.SELECTED and not candidates:` override in "
            "`study_split._training_selection`."
        ),
        target_relpath="src/mtp/cellulose_decr/study_split.py",
        needle=(
            "    status = selection.status\n"
            "    if status == scanner.SELECTED and not candidates:\n"
            '        status = "no_testable_study_disjoint_triple"\n'
            "    return candidates, status"
        ),
        replacement=("    status = selection.status\n    return candidates, status"),
        pytest_nodeid=(
            "tests/test_study_split_recheck_and_untestable_997.py::"
            "test_training_selection_and_run_splits_untestable_reasons"
            "[no_testable_study_disjoint_triple]"
        ),
    ),
    MutantSpec(
        mutant_id="M6_invert_candidates_guard_in_training_selection",
        description=(
            "Invert `not candidates` to `candidates` in `study_split._training_selection` "
            "so `selected` splits are mislabelled and empty disjoint-triple splits stay `selected`."
        ),
        target_relpath="src/mtp/cellulose_decr/study_split.py",
        needle=(
            "    if status == scanner.SELECTED and not candidates:\n"
            '        status = "no_testable_study_disjoint_triple"'
        ),
        replacement=(
            "    if status == scanner.SELECTED and candidates:\n"
            '        status = "no_testable_study_disjoint_triple"'
        ),
        pytest_nodeid=(
            "tests/test_study_split_recheck_and_untestable_997.py::"
            "test_training_selection_and_run_splits_untestable_reasons[selected]"
        ),
    ),
    MutantSpec(
        mutant_id="M7_relabel_no_traits_as_fewer_than_three_conditions",
        description=(
            "Collapse `no_traits` into `fewer_than_three_eligible_conditions` in "
            "`study_split._training_selection`."
        ),
        target_relpath="src/mtp/cellulose_decr/study_split.py",
        needle="    status = selection.status",
        replacement=(
            "    status = (\n"
            '        "fewer_than_three_eligible_conditions"\n'
            '        if selection.status == "no_traits"\n'
            "        else selection.status\n"
            "    )"
        ),
        pytest_nodeid=(
            "tests/test_study_split_recheck_and_untestable_997.py::"
            "test_training_selection_and_run_splits_untestable_reasons[no_traits]"
        ),
    ),
)


class _CallPhaseInspector:
    """In-process pytest plugin recording call-phase outcome and exception type."""

    def __init__(self) -> None:
        self.call_reports: list[pytest.TestReport] = []
        self.non_call_failures: list[tuple[str, str]] = []
        self.call_exc_types: list[str] = []

    def pytest_runtest_logreport(self, report: pytest.TestReport) -> None:
        if report.when == "call":
            self.call_reports.append(report)

    def pytest_exception_interact(
        self,
        node: object,
        call: pytest.CallInfo[object],
        report: pytest.TestReport,
    ) -> None:
        _ = node
        exc_name = call.excinfo.type.__name__ if call.excinfo is not None else "UnknownError"
        if report.when == "call":
            self.call_exc_types.append(exc_name)
        elif report.failed:
            self.non_call_failures.append((report.when, exc_name))


def _evict_project_modules() -> None:
    exact_names = {"mtp", "tests"}
    for mod_name in list(sys.modules):
        if mod_name in exact_names or mod_name.startswith(("mtp.", "tests.", "test_")):
            sys.modules.pop(mod_name, None)
    importlib.invalidate_caches()


def _populate_scratch_repo(repo_root: Path, scratch_repo: Path) -> None:
    ignore = shutil.ignore_patterns("__pycache__", "*.pyc")
    for directory in ("src", "methods", "schema", "workflows"):
        shutil.copytree(repo_root / directory, scratch_repo / directory, ignore=ignore)
    for filename in ("pyproject.toml", "Justfile", "uv.lock"):
        src_file = repo_root / filename
        if src_file.is_file():
            shutil.copy2(src_file, scratch_repo / filename)
    dst = scratch_repo / "tests/test_study_split_recheck_and_untestable_997.py"
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(
        repo_root / "tests/test_study_split_recheck_and_untestable_997.py",
        dst,
    )


def _run_pytest_node(
    repo_root: Path,
    scratch_repo: Path,
    nodeid: str,
    junit_xml: Path,
) -> tuple[int, _CallPhaseInspector, str]:
    inspector = _CallPhaseInspector()
    stdout_buf = io.StringIO()
    stderr_buf = io.StringIO()
    old_cwd = Path.cwd()
    old_sys_path = list(sys.path)
    old_bytecode = sys.dont_write_bytecode
    old_mpl = os.environ.get("MPLBACKEND")
    os.environ["MPLBACKEND"] = "Agg"
    sys.dont_write_bytecode = True
    try:
        _evict_project_modules()
        filtered_path = [
            entry
            for entry in old_sys_path
            if Path(entry or ".").resolve()
            not in {repo_root.resolve(), (repo_root / "src").resolve()}
        ]
        sys.path[:] = [str(scratch_repo / "src"), str(scratch_repo), *filtered_path]
        os.chdir(scratch_repo)
        with redirect_stdout(stdout_buf), redirect_stderr(stderr_buf):
            exit_code = int(
                pytest.main(
                    [
                        "-q",
                        "-o",
                        "addopts=",
                        f"--rootdir={scratch_repo}",
                        f"--junitxml={junit_xml}",
                        nodeid,
                    ],
                    plugins=[inspector],
                )
            )
    finally:
        os.chdir(old_cwd)
        sys.path[:] = old_sys_path
        sys.dont_write_bytecode = old_bytecode
        _evict_project_modules()
        if old_mpl is None:
            os.environ.pop("MPLBACKEND", None)
        else:
            os.environ["MPLBACKEND"] = old_mpl

    combined_output = (stdout_buf.getvalue() + "\n" + stderr_buf.getvalue()).strip()
    return exit_code, inspector, combined_output


def verify_mutants(
    repo_root: Path,
    scratch_parent: Path,
    *,
    canary_only: bool = False,
) -> tuple[Path, list[dict[str, object]]]:
    scratch_parent.mkdir(parents=True, exist_ok=True)
    run_label = "canary" if canary_only else "all"
    scratch_dir = Path(
        tempfile.mkdtemp(prefix=f"verify_997_{run_label}_", dir=str(scratch_parent))
    ).resolve()
    scratch_repo = scratch_dir / "repo"
    _populate_scratch_repo(repo_root, scratch_repo)

    junit_base = scratch_dir / "junit_baseline.xml"
    code, inspector, out = _run_pytest_node(
        repo_root,
        scratch_repo,
        "tests/test_study_split_recheck_and_untestable_997.py",
        junit_base,
    )
    if code != 0 or inspector.non_call_failures or not junit_base.is_file():
        raise RuntimeError(f"Baseline test failed in scratch_repo (exit={code}):\n{out}")

    selected = MUTANTS[:1] if canary_only else MUTANTS
    verdicts: list[dict[str, object]] = []

    for spec in selected:
        target_path = scratch_repo / spec.target_relpath
        original_text = target_path.read_text(encoding="utf-8")
        occurrences = original_text.count(spec.needle)
        if occurrences != 1:
            raise RuntimeError(
                f"{spec.mutant_id}: expected exactly 1 occurrence of needle in "
                f"{spec.target_relpath}, found {occurrences}"
            )
        mutated_text = original_text.replace(spec.needle, spec.replacement, 1)
        target_path.write_text(mutated_text, encoding="utf-8")
        junit_xml = scratch_dir / f"junit_{spec.mutant_id}.xml"
        try:
            exit_code, inspector, output = _run_pytest_node(
                repo_root, scratch_repo, spec.pytest_nodeid, junit_xml
            )
        finally:
            target_path.write_text(original_text, encoding="utf-8")

        junit_exists = junit_xml.is_file() and junit_xml.stat().st_size > 0
        call_failed = (
            len(inspector.call_reports) == 1
            and inspector.call_reports[0].failed
            and not inspector.non_call_failures
        )
        killed_by_assertion = (
            exit_code != 0
            and junit_exists
            and call_failed
            and inspector.call_exc_types == ["AssertionError"]
        )
        tail_lines = [line for line in output.splitlines() if line.strip()][-6:]
        detail = (
            f"exit_code={exit_code}; call_exc_types={inspector.call_exc_types}; "
            f"non_call_failures={inspector.non_call_failures}; "
            f"junit_xml={junit_xml.name}; summary={' | '.join(tail_lines)}"
        )
        verdict: dict[str, object] = {
            "mutant_id": spec.mutant_id,
            "description": spec.description,
            "target_relpath": spec.target_relpath,
            "pytest_nodeid": spec.pytest_nodeid,
            "status": "KILLED" if killed_by_assertion else "SURVIVED",
            "killed": killed_by_assertion,
            "exit_code": exit_code,
            "call_exc_types": inspector.call_exc_types,
            "non_call_failures": inspector.non_call_failures,
            "junit_xml": str(junit_xml),
            "detail": detail,
        }
        verdicts.append(verdict)
        if not killed_by_assertion:
            verdicts_path = scratch_dir / "verdicts.json"
            verdicts_path.write_text(json.dumps(verdicts, indent=2) + "\n", encoding="utf-8")
            raise RuntimeError(
                f"Mutant {spec.mutant_id} was not killed by call-phase AssertionError: {detail}"
            )

    verdicts_path = scratch_dir / "verdicts.json"
    verdicts_path.write_text(json.dumps(verdicts, indent=2) + "\n", encoding="utf-8")
    return verdicts_path, verdicts


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--scratch-dir",
        type=Path,
        required=True,
        help="Parent directory under which a persistent scratch run directory will be created.",
    )
    parser.add_argument(
        "--canary",
        action="store_true",
        help="Run only the first mutant (M1) after baseline verification.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    repo_root = Path(__file__).resolve().parents[2]
    verdicts_path, verdicts = verify_mutants(
        repo_root,
        args.scratch_dir.resolve(),
        canary_only=bool(args.canary),
    )
    killed_count = sum(1 for item in verdicts if item.get("killed") is True)
    print(
        f"Verified {killed_count}/{len(verdicts)} mutants killed by call-phase AssertionError. "
        f"Verdicts written to {verdicts_path}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
