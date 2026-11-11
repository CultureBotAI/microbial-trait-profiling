"""Scratch-isolated mutant verifier for Issue #1992.

Purpose:
    Verify that the direct assertion checks in `tests/test_missing_trait_fraction_1992.py`
    fail by call-phase `AssertionError` when `study_direction_means` omits or weakens
    the missing-trait guard (`guild not in profile.guild_counts`), when `not in` is
    inverted to `in`, when only the first profile or first guild is checked, when
    the `"missing is not zero"` error message is altered, or when test assertions are
    mutated.

Inputs:
    - Worktree copies of `src/mtp/` and `tests/test_missing_trait_fraction_1992.py`.

Outputs:
    - Per-run JUnit XML files (`junit_<label>.xml`) and `verdicts.json` inside a
      persistent subdirectory under `--scratch-dir`.

Issue: #1992
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
        mutant_id="M1_remove_missing_trait_guard",
        description=(
            "study_direction_means omits the guild not in profile.guild_counts check; "
            "must fail by AssertionError in "
            "test_comparison_rows_rejects_absent_table1_trait_and_accepts_explicit_zero."
        ),
        target_relpath="src/mtp/guild_maps.py",
        needle=(
            "            if guild not in profile.guild_counts:\n"
            '                raise ValueError(f"Missing trait {guild} for {key_id}; '
            'missing is not zero")\n'
        ),
        replacement="",
        pytest_nodeid=(
            "tests/test_missing_trait_fraction_1992.py::"
            "test_comparison_rows_rejects_absent_table1_trait_and_accepts_explicit_zero"
        ),
    ),
    MutantSpec(
        mutant_id="M2_invert_not_in_check",
        description=(
            "study_direction_means replaces `guild not in profile.guild_counts` with "
            "`guild in profile.guild_counts`; must fail by AssertionError when complete "
            "guild_counts are passed in "
            "test_comparison_rows_rejects_absent_table1_trait_and_accepts_explicit_zero."
        ),
        target_relpath="src/mtp/guild_maps.py",
        needle="            if guild not in profile.guild_counts:",
        replacement="            if guild in profile.guild_counts:",
        pytest_nodeid=(
            "tests/test_missing_trait_fraction_1992.py::"
            "test_comparison_rows_rejects_absent_table1_trait_and_accepts_explicit_zero"
        ),
    ),
    MutantSpec(
        mutant_id="M3_check_only_first_profile",
        description=(
            "study_direction_means checks guild_counts only on the first profile in "
            "`profiles`; must fail by AssertionError in "
            "test_study_direction_means_rejects_missing_trait_in_later_profile_and_later_guild."
        ),
        target_relpath="src/mtp/guild_maps.py",
        needle=(
            "    for profile in profiles:\n"
            "        key_id = (profile.study_id, profile.experiment_id)\n"
            "        for guild in guild_order:\n"
            "            if guild not in profile.guild_counts:"
        ),
        replacement=(
            "    for profile_idx, profile in enumerate(profiles):\n"
            "        key_id = (profile.study_id, profile.experiment_id)\n"
            "        for guild in guild_order:\n"
            "            if profile_idx == 0 and guild not in profile.guild_counts:"
        ),
        pytest_nodeid=(
            "tests/test_missing_trait_fraction_1992.py::"
            "test_study_direction_means_rejects_missing_trait_in_later_profile_and_later_guild"
        ),
    ),
    MutantSpec(
        mutant_id="M4_check_only_first_guild",
        description=(
            "study_direction_means checks guild_counts only for guild_order[0]; "
            "must fail by AssertionError in "
            "test_study_direction_means_rejects_missing_trait_in_later_profile_and_later_guild."
        ),
        target_relpath="src/mtp/guild_maps.py",
        needle="            if guild not in profile.guild_counts:",
        replacement="            if guild == guild_order[0] and guild not in profile.guild_counts:",
        pytest_nodeid=(
            "tests/test_missing_trait_fraction_1992.py::"
            "test_study_direction_means_rejects_missing_trait_in_later_profile_and_later_guild"
        ),
    ),
    MutantSpec(
        mutant_id="M5_alter_missing_is_not_zero_message",
        description=(
            "study_direction_means omits `; missing is not zero` from the ValueError "
            "message; must fail by AssertionError in "
            "test_compute_direction_shifts_rejects_missing_trait."
        ),
        target_relpath="src/mtp/guild_maps.py",
        needle='f"Missing trait {guild} for {key_id}; missing is not zero"',
        replacement='f"Missing trait {guild} for {key_id}"',
        pytest_nodeid=(
            "tests/test_missing_trait_fraction_1992.py::"
            "test_compute_direction_shifts_rejects_missing_trait"
        ),
    ),
    MutantSpec(
        mutant_id="M6_test_assertion_mutant_study_direction_means",
        description=(
            "Mutating the expected error string in "
            "test_study_direction_means_rejects_missing_trait_in_later_profile_and_later_guild "
            "must fail by call-phase AssertionError."
        ),
        target_relpath="tests/test_missing_trait_fraction_1992.py",
        needle=(
            "    assert str(caught) == ("
            "\"Missing trait anaerobe for ('study_2', 'exp_2'); missing is not zero\""
            ")"
        ),
        replacement=(
            "    assert str(caught) == ("
            "\"Missing trait aerobe for ('study_1', 'exp_1'); missing is not zero\""
            ")"
        ),
        pytest_nodeid=(
            "tests/test_missing_trait_fraction_1992.py::"
            "test_study_direction_means_rejects_missing_trait_in_later_profile_and_later_guild"
        ),
    ),
    MutantSpec(
        mutant_id="M7_test_assertion_mutant_comparison_rows_tied",
        description=(
            "Mutating the expected tied-study count for explicit zero count in "
            "test_comparison_rows_rejects_absent_table1_trait_and_accepts_explicit_zero "
            "must fail by call-phase AssertionError."
        ),
        target_relpath="tests/test_missing_trait_fraction_1992.py",
        needle='    assert cellulose_row["studies_tied"] == 1',
        replacement='    assert cellulose_row["studies_tied"] == 0',
        pytest_nodeid=(
            "tests/test_missing_trait_fraction_1992.py::"
            "test_comparison_rows_rejects_absent_table1_trait_and_accepts_explicit_zero"
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


def _is_project_module(mod_name: str) -> bool:
    return mod_name in {"mtp", "tests"} or mod_name.startswith(("mtp.", "tests.", "test_"))


def _evict_project_modules() -> None:
    for mod_name in list(sys.modules):
        if _is_project_module(mod_name):
            sys.modules.pop(mod_name, None)
    importlib.invalidate_caches()


def _populate_scratch_repo(repo_root: Path, scratch_repo: Path) -> None:
    ignore = shutil.ignore_patterns("__pycache__", "*.pyc")
    shutil.copytree(repo_root / "src", scratch_repo / "src", ignore=ignore)
    for filename in ("pyproject.toml", "Justfile", "uv.lock"):
        src_file = repo_root / filename
        if src_file.is_file():
            shutil.copy2(src_file, scratch_repo / filename)
    dst = scratch_repo / "tests/test_missing_trait_fraction_1992.py"
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(repo_root / "tests/test_missing_trait_fraction_1992.py", dst)
    verifier_rel = "scripts/adhoc/verify_1992_mutants_agy_20261111.py"
    verifier_dst = scratch_repo / verifier_rel
    verifier_dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(repo_root / verifier_rel, verifier_dst)


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
        sys.modules.update(saved_project_modules)
        importlib.invalidate_caches()
        if old_mpl is None:
            os.environ.pop("MPLBACKEND", None)
        else:
            os.environ["MPLBACKEND"] = old_mpl

    combined_output = (stdout_buf.getvalue() + "\n" + stderr_buf.getvalue()).strip()
    return exit_code, inspector, combined_output


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
    return killed_by_assertion, detail


def verify_mutants(
    repo_root: Path,
    scratch_parent: Path,
    *,
    canary_only: bool = False,
) -> tuple[Path, list[dict[str, object]]]:
    scratch_parent.mkdir(parents=True, exist_ok=True)
    run_label = "canary" if canary_only else "all"
    scratch_dir = Path(
        tempfile.mkdtemp(prefix=f"verify_1992_{run_label}_", dir=str(scratch_parent))
    ).resolve()
    scratch_repo = scratch_dir / "repo"
    _populate_scratch_repo(repo_root, scratch_repo)

    junit_base = scratch_dir / "junit_baseline.xml"
    code, inspector, out = _run_pytest_node(
        repo_root,
        scratch_repo,
        "tests/test_missing_trait_fraction_1992.py",
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
            killed_by_assertion, detail = _run_targeted_pytest(
                scratch_repo,
                spec.pytest_nodeid,
                junit_xml,
                repo_root=repo_root,
            )
        finally:
            target_path.write_text(original_text, encoding="utf-8")

        verdict: dict[str, object] = {
            "mutant_id": spec.mutant_id,
            "description": spec.description,
            "target_relpath": spec.target_relpath,
            "pytest_nodeid": spec.pytest_nodeid,
            "status": "KILLED" if killed_by_assertion else "SURVIVED",
            "killed": killed_by_assertion,
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
