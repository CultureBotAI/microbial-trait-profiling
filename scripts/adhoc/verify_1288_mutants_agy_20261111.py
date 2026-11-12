"""Scratch-isolated mutant verifier for Issue #1288 (`uv.lock` revision & `upload-time` guard).

Purpose:
    Verify that `tests/test_uv_lock_revision_1288.py` fails by call-phase `AssertionError`
    (with exact `pytest.ExitCode.TESTS_FAILED == 1`) when `uv.lock` regresses to
    `revision = 1` or `revision = 4`, when an `sdist` or `wheel` entry in `uv.lock`
    loses its `upload-time` field, or when the validator's revision, `sdist`, or
    `wheel` checks are weakened.

Inputs:
    - Worktree copies of `uv.lock`, `tests/test_uv_lock_revision_1288.py`, and
      `scripts/adhoc/verify_1288_mutants_agy_20261111.py`.

Outputs:
    - Per-run JUnit XML files (`junit_<label>.xml`) and `verdicts.json` inside a
      persistent subdirectory under `--scratch-dir` (strictly outside `REPO_ROOT`).

Issue: #1288
Date: 2026-11-11
"""

from __future__ import annotations

import argparse
import hashlib
import html
import importlib
import io
import json
import os
import re
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
        mutant_id="M1_lock_revision_1",
        description=(
            "uv.lock regresses from revision = 5 to revision = 1; must fail by "
            "call-phase AssertionError in test_committed_uv_lock_has_revision_5_and_upload_times."
        ),
        target_relpath="uv.lock",
        needle="version = 1\nrevision = 5\n",
        replacement="version = 1\nrevision = 1\n",
        pytest_nodeid=(
            "tests/test_uv_lock_revision_1288.py::"
            "test_committed_uv_lock_has_revision_5_and_upload_times"
        ),
    ),
    MutantSpec(
        mutant_id="M2_lock_revision_4",
        description=(
            "uv.lock regresses from revision = 5 to revision = 4; must fail by "
            "call-phase AssertionError in test_committed_uv_lock_has_revision_5_and_upload_times."
        ),
        target_relpath="uv.lock",
        needle="version = 1\nrevision = 5\n",
        replacement="version = 1\nrevision = 4\n",
        pytest_nodeid=(
            "tests/test_uv_lock_revision_1288.py::"
            "test_committed_uv_lock_has_revision_5_and_upload_times"
        ),
    ),
    MutantSpec(
        mutant_id="M3_strip_alabaster_sdist_upload_time",
        description=(
            "uv.lock strips upload-time from alabaster sdist; must fail by "
            "call-phase AssertionError in test_committed_uv_lock_has_revision_5_and_upload_times."
        ),
        target_relpath="uv.lock",
        needle='size = 24210, upload-time = "2024-07-26T18:15:03.762Z" }',
        replacement="size = 24210 }",
        pytest_nodeid=(
            "tests/test_uv_lock_revision_1288.py::"
            "test_committed_uv_lock_has_revision_5_and_upload_times"
        ),
    ),
    MutantSpec(
        mutant_id="M4_strip_alabaster_wheel_upload_time",
        description=(
            "uv.lock strips upload-time from alabaster wheel; must fail by "
            "call-phase AssertionError in test_committed_uv_lock_has_revision_5_and_upload_times."
        ),
        target_relpath="uv.lock",
        needle='size = 13929, upload-time = "2024-07-26T18:15:02.05Z" },',
        replacement="size = 13929 },",
        pytest_nodeid=(
            "tests/test_uv_lock_revision_1288.py::"
            "test_committed_uv_lock_has_revision_5_and_upload_times"
        ),
    ),
    MutantSpec(
        mutant_id="M5_weaken_min_revision_check_in_validator",
        description=(
            "validate_uv_lock_format weakens revision < 5 to revision < 1; must fail by "
            "call-phase AssertionError in "
            "test_validate_uv_lock_format_rejects_regressed_or_missing_revision."
        ),
        target_relpath="tests/test_uv_lock_revision_1288.py",
        needle=(
            "    if isinstance(revision, bool) or not isinstance(revision, int) or revision < 5:"
        ),
        replacement=(
            "    if isinstance(revision, bool) or not isinstance(revision, int) or revision < 1:"
        ),
        pytest_nodeid=(
            "tests/test_uv_lock_revision_1288.py::"
            "test_validate_uv_lock_format_rejects_regressed_or_missing_revision"
        ),
    ),
    MutantSpec(
        mutant_id="M6_skip_sdist_upload_time_check_in_validator",
        description=(
            "validate_uv_lock_format skips sdist upload-time validation; must fail by "
            "call-phase AssertionError in "
            "test_validate_uv_lock_format_rejects_missing_or_blank_sdist_upload_time."
        ),
        target_relpath="tests/test_uv_lock_revision_1288.py",
        needle="            if not _has_nonempty_upload_time(sdist):",
        replacement="            if False:",
        pytest_nodeid=(
            "tests/test_uv_lock_revision_1288.py::"
            "test_validate_uv_lock_format_rejects_missing_or_blank_sdist_upload_time"
        ),
    ),
    MutantSpec(
        mutant_id="M7_skip_wheel_upload_time_check_in_validator",
        description=(
            "validate_uv_lock_format skips wheel upload-time validation; must fail by "
            "call-phase AssertionError in "
            "test_validate_uv_lock_format_rejects_missing_or_blank_wheel_upload_time."
        ),
        target_relpath="tests/test_uv_lock_revision_1288.py",
        needle="                if not _has_nonempty_upload_time(wheel):",
        replacement="                if False:",
        pytest_nodeid=(
            "tests/test_uv_lock_revision_1288.py::"
            "test_validate_uv_lock_format_rejects_missing_or_blank_wheel_upload_time"
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


def validate_scratch_dir(scratch_dir: Path, repo_root: Path) -> Path:
    """Reject any scratch_dir lexically or resolved inside repo_root before mkdir."""
    resolved_repo = Path(repo_root).expanduser().resolve()
    lexical_repo = Path(os.path.abspath(Path(repo_root).expanduser()))
    lexical_scratch = Path(os.path.abspath(Path(scratch_dir).expanduser()))
    resolved_scratch = Path(scratch_dir).expanduser().resolve(strict=False)
    if any(
        candidate == root or candidate.is_relative_to(root)
        for candidate in (resolved_scratch, lexical_scratch)
        for root in (resolved_repo, lexical_repo)
    ):
        raise ValueError(
            f"--scratch-dir must be outside repository root {resolved_repo}; got {scratch_dir}"
        )
    return lexical_scratch


def _is_project_module(mod_name: str) -> bool:
    return mod_name in {"mtp", "tests"} or mod_name.startswith(("mtp.", "tests.", "test_"))


def _evict_project_modules() -> None:
    for mod_name in list(sys.modules):
        if _is_project_module(mod_name):
            sys.modules.pop(mod_name, None)
    importlib.invalidate_caches()


def _populate_scratch_repo(repo_root: Path, scratch_repo: Path) -> None:
    scratch_repo.mkdir(parents=True, exist_ok=True)
    for filename in ("pyproject.toml", "uv.lock"):
        src_file = repo_root / filename
        if src_file.is_file():
            shutil.copy2(src_file, scratch_repo / filename)
    for relpath in (
        "tests/test_uv_lock_revision_1288.py",
        "scripts/adhoc/verify_1288_mutants_agy_20261111.py",
    ):
        dst = scratch_repo / relpath
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(repo_root / relpath, dst)


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
    saved_project_modules = {
        name: mod for name, mod in sys.modules.items() if _is_project_module(name)
    }
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
        sys.path[:] = [str(scratch_repo), *filtered_path]
        os.chdir(scratch_repo)
        with redirect_stdout(stdout_buf), redirect_stderr(stderr_buf):
            exit_code = int(
                pytest.main(
                    [
                        "-q",
                        "-o",
                        "addopts=",
                        "-p",
                        "no:cacheprovider",
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
        sys.modules.update(saved_project_modules)
        importlib.invalidate_caches()
        if old_mpl is None:
            os.environ.pop("MPLBACKEND", None)
        else:
            os.environ["MPLBACKEND"] = old_mpl

    combined_output = (stdout_buf.getvalue() + "\n" + stderr_buf.getvalue()).strip()
    return exit_code, inspector, combined_output


def evaluate_mutant_outcome(
    exit_code: int,
    inspector: _CallPhaseInspector,
    junit_xml: Path,
    output: str,
) -> tuple[bool, str]:
    """Validate that a mutant run failed with TESTS_FAILED (1) and call-phase AssertionError."""
    junit_ok = False
    failure_text = ""
    if junit_xml.is_file() and junit_xml.stat().st_size > 0:
        xml = junit_xml.read_text(encoding="utf-8")
        has_error_or_skip = bool(re.search(r"<(?:error|skipped)\b", xml))
        failures = re.findall(r"<failure\b[^>]*>(.*?)</failure>", xml, re.DOTALL)
        if not has_error_or_skip and len(failures) == 1:
            failure_text = html.unescape(failures[0]).strip()
            junit_ok = bool(failure_text)

    call_failed = (
        len(inspector.call_reports) == 1
        and inspector.call_reports[0].failed
        and not inspector.non_call_failures
    )
    killed_by_assertion = (
        exit_code == int(pytest.ExitCode.TESTS_FAILED)
        and junit_ok
        and call_failed
        and inspector.call_exc_types == ["AssertionError"]
    )
    tail_lines = [line for line in output.splitlines() if line.strip()][-6:]
    detail = (
        f"exit_code={exit_code}; call_exc_types={inspector.call_exc_types}; "
        f"non_call_failures={inspector.non_call_failures}; "
        f"junit_xml={junit_xml.name}; summary={' | '.join(tail_lines)}"
    )
    return killed_by_assertion, detail


def _run_targeted_pytest(
    scratch_repo: Path,
    nodeid: str,
    junit_xml: Path,
    *,
    repo_root: Path | None = None,
) -> tuple[bool, str]:
    """Run one pytest nodeid and return (killed_by_call_assertion, detail)."""
    effective_root = repo_root if repo_root is not None else Path(__file__).resolve().parents[2]
    exit_code, inspector, output = _run_pytest_node(effective_root, scratch_repo, nodeid, junit_xml)
    return evaluate_mutant_outcome(exit_code, inspector, junit_xml, output)


def verify_mutants(
    repo_root: Path,
    scratch_parent: Path,
    *,
    canary_only: bool = False,
) -> tuple[Path, list[dict[str, object]]]:
    validated_parent = validate_scratch_dir(scratch_parent, repo_root)
    validated_parent.mkdir(parents=True, exist_ok=True)
    run_label = "canary" if canary_only else "all"
    scratch_dir = Path(
        tempfile.mkdtemp(prefix=f"verify_1288_{run_label}_", dir=str(validated_parent))
    ).resolve()
    scratch_repo = scratch_dir / "repo"

    tracked_relpaths = (
        "uv.lock",
        "tests/test_uv_lock_revision_1288.py",
        "scripts/adhoc/verify_1288_mutants_agy_20261111.py",
    )
    before_hashes = {
        rel: hashlib.sha256((repo_root / rel).read_bytes()).hexdigest() for rel in tracked_relpaths
    }

    _populate_scratch_repo(repo_root, scratch_repo)

    junit_base = scratch_dir / "junit_baseline.xml"
    code, inspector, out = _run_pytest_node(
        repo_root,
        scratch_repo,
        "tests/test_uv_lock_revision_1288.py",
        junit_base,
    )
    if code != int(pytest.ExitCode.OK) or inspector.non_call_failures or not junit_base.is_file():
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
            killed_by_assertion, detail = _run_targeted_pytest(
                scratch_repo,
                spec.pytest_nodeid,
                junit_xml,
                repo_root=repo_root,
            )
        finally:
            target_path.write_text(original_text, encoding="utf-8")

        after_hashes = {
            rel: hashlib.sha256((repo_root / rel).read_bytes()).hexdigest()
            for rel in tracked_relpaths
        }
        canonical_unchanged = before_hashes == after_hashes
        verdict: dict[str, object] = {
            "mutant_id": spec.mutant_id,
            "description": spec.description,
            "target_relpath": spec.target_relpath,
            "pytest_nodeid": spec.pytest_nodeid,
            "status": "KILLED" if (killed_by_assertion and canonical_unchanged) else "SURVIVED",
            "killed": bool(killed_by_assertion and canonical_unchanged),
            "canonical_unchanged": canonical_unchanged,
            "junit_xml": str(junit_xml),
            "detail": detail,
        }
        verdicts.append(verdict)
        if not killed_by_assertion or not canonical_unchanged:
            verdicts_path = scratch_dir / "verdicts.json"
            verdicts_path.write_text(json.dumps(verdicts, indent=2) + "\n", encoding="utf-8")
            raise RuntimeError(
                f"Mutant {spec.mutant_id} was not cleanly killed: {detail} "
                f"(canonical_unchanged={canonical_unchanged})"
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
        help="Parent directory outside REPO_ROOT under which a scratch run dir will be created.",
    )
    mode_group = parser.add_mutually_exclusive_group()
    mode_group.add_argument(
        "--canary",
        action="store_true",
        help="Run only the first mutant (M1) after baseline verification.",
    )
    mode_group.add_argument(
        "--all",
        action="store_true",
        help="Run all mutants after baseline verification (default when --canary is omitted).",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    repo_root = Path(__file__).resolve().parents[2]
    verdicts_path, verdicts = verify_mutants(
        repo_root,
        args.scratch_dir,
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
