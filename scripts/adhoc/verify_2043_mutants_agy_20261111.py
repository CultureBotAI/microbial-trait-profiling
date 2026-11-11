"""Scratch-isolated mutant verifier for Issue #2043.

Purpose:
    Verify that the direct assertion checks in `tests/test_arm_trait_key_policy.py`
    and `tests/test_contribution_caller_policy.py` fail by call-phase `AssertionError`
    (rather than `ValueError` or silent pass) when Noel omits or miswires
    `--trait-key-policy`, when `_step_trait_key_policy` defaults to union instead of
    `None`, or when `taxon_contributions.main` omits or hardcodes `trait_key_policy`
    when calling `load_graph`.

Inputs:
    - Worktree copies of `src/mtp/`, `methods/`, `schema/`, `workflows/`,
      `scripts/create_pd_trait_contributors.py`,
      `tests/fixtures/trait_key_policy/`,
      `tests/test_arm_trait_key_policy.py`, and
      `tests/test_contribution_caller_policy.py`.

Outputs:
    - Per-run JUnit XML files (`junit_<label>.xml`) and `verdicts.json` inside a
      persistent subdirectory under `--scratch-dir`.

Issue: #2043
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
        mutant_id="M1_noel_omits_policy_arg",
        description=(
            "Noel build_commands omits --trait-key-policy; must fail by AssertionError "
            "in test_noel_ledger_names_the_resolved_core_policy rather than ValueError."
        ),
        target_relpath="src/mtp/workflow/arms/noel.py",
        needle=(
            "    trait_policy_args = (\n"
            "        () if context.trait_key_policy is None "
            'else ("--trait-key-policy", context.trait_key_policy)\n'
            "    )"
        ),
        replacement="    trait_policy_args: tuple[str, ...] = ()",
        pytest_nodeid=(
            "tests/test_arm_trait_key_policy.py::"
            "test_noel_ledger_names_the_resolved_core_policy[raw_object_v1]"
        ),
    ),
    MutantSpec(
        mutant_id="M2_noel_raw_core_carries_union_policy",
        description=(
            "Noel build_commands hardcodes facultative_union_v1 instead of "
            "context.trait_key_policy; must fail by AssertionError in "
            "test_raw_core_ledger_command_carries_raw_object_v1."
        ),
        target_relpath="src/mtp/workflow/arms/noel.py",
        needle=(
            "    trait_policy_args = (\n"
            "        () if context.trait_key_policy is None "
            'else ("--trait-key-policy", context.trait_key_policy)\n'
            "    )"
        ),
        replacement=(
            "    trait_policy_args = (\n"
            "        () if context.trait_key_policy is None "
            'else ("--trait-key-policy", "facultative_union_v1")\n'
            "    )"
        ),
        pytest_nodeid=(
            "tests/test_arm_trait_key_policy.py::test_raw_core_ledger_command_carries_raw_object_v1"
        ),
    ),
    MutantSpec(
        mutant_id="M3_historical_raw_core_resolver_returns_union",
        description=(
            "resolve_trait_policy returns facultative_union_v1 for historical raw cores; "
            "must fail by AssertionError in test_raw_core_ledger_command_carries_raw_object_v1."
        ),
        target_relpath="src/mtp/workflow/trait_policy.py",
        needle="        return ResolvedTraitPolicy(HISTORICAL_MATRIX_POLICIES[identity], False)",
        replacement='        return ResolvedTraitPolicy("facultative_union_v1", False)',
        pytest_nodeid=(
            "tests/test_arm_trait_key_policy.py::test_raw_core_ledger_command_carries_raw_object_v1"
        ),
    ),
    MutantSpec(
        mutant_id="M4_step_policy_defaults_to_union",
        description=(
            "runtime._step_trait_key_policy omits parser.set_defaults(trait_key_policy=None) "
            "so omitted flag defaults to union; must fail by AssertionError in "
            "test_step_policy_parser_returns_none_when_flag_is_omitted."
        ),
        target_relpath="src/mtp/workflow/arms/runtime.py",
        needle=(
            "    parser = SCOPED_MODULES[command.argv[2]](_StepArguments)\n"
            "    parser.set_defaults(trait_key_policy=None)"
        ),
        replacement="    parser = SCOPED_MODULES[command.argv[2]](_StepArguments)",
        pytest_nodeid=(
            "tests/test_arm_trait_key_policy.py::"
            "test_step_policy_parser_returns_none_when_flag_is_omitted"
        ),
    ),
    MutantSpec(
        mutant_id="M5_taxon_contributions_ignores_policy",
        description=(
            "taxon_contributions.main omits trait_key_policy=args.trait_key_policy in load_graph; "
            "must fail by AssertionError in "
            "test_taxon_contributions_main_forwards_explicit_trait_key_policy_to_load_graph."
        ),
        target_relpath="src/mtp/taxon_contributions.py",
        needle=(
            "        strain_rollup=args.strain_rollup,\n"
            "        trait_key_policy=args.trait_key_policy,\n"
            "    )"
        ),
        replacement=("        strain_rollup=args.strain_rollup,\n    )"),
        pytest_nodeid=(
            "tests/test_contribution_caller_policy.py::"
            "test_taxon_contributions_main_forwards_explicit_trait_key_policy_to_load_graph[raw_object_v1]"
        ),
    ),
    MutantSpec(
        mutant_id="M6_taxon_contributions_hardcodes_default_policy",
        description=(
            "taxon_contributions.main passes DEFAULT_TRAIT_KEY_POLICY instead of "
            "args.trait_key_policy to load_graph; must fail by AssertionError in "
            "test_taxon_contributions_main_forwards_explicit_trait_key_policy_to_load_graph."
        ),
        target_relpath="src/mtp/taxon_contributions.py",
        needle=(
            "        strain_rollup=args.strain_rollup,\n"
            "        trait_key_policy=args.trait_key_policy,\n"
            "    )"
        ),
        replacement=(
            "        strain_rollup=args.strain_rollup,\n"
            "        trait_key_policy=DEFAULT_TRAIT_KEY_POLICY,\n"
            "    )"
        ),
        pytest_nodeid=(
            "tests/test_contribution_caller_policy.py::"
            "test_taxon_contributions_main_forwards_explicit_trait_key_policy_to_load_graph[raw_object_v1]"
        ),
    ),
    MutantSpec(
        mutant_id="M7_pairs_helper_drops_last_pair",
        description=(
            "Test helper _pairs drops the final adjacent pair; must fail by AssertionError "
            "in test_raw_core_ledger_command_carries_raw_object_v1."
        ),
        target_relpath="tests/test_arm_trait_key_policy.py",
        needle=(
            "def _pairs(argv: Sequence[str]) -> tuple[tuple[str, str], ...]:\n"
            "    return tuple((argv[i], argv[i + 1]) for i in range(len(argv) - 1))"
        ),
        replacement=(
            "def _pairs(argv: Sequence[str]) -> tuple[tuple[str, str], ...]:\n"
            "    return tuple((argv[i], argv[i + 1]) for i in range(len(argv) - 2))"
        ),
        pytest_nodeid=(
            "tests/test_arm_trait_key_policy.py::test_raw_core_ledger_command_carries_raw_object_v1"
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
    exact_names = {"mtp", "tests", "policy_pd_contributors", "create_pd_cellulose_analysis"}
    for mod_name in list(sys.modules):
        if mod_name in exact_names or mod_name.startswith(("mtp.", "tests.", "test_")):
            sys.modules.pop(mod_name, None)
    importlib.invalidate_caches()


def _populate_scratch_repo(repo_root: Path, scratch_repo: Path) -> None:
    ignore = shutil.ignore_patterns("__pycache__", "*.pyc")
    for directory in ("src", "scripts", "methods", "schema", "workflows"):
        shutil.copytree(repo_root / directory, scratch_repo / directory, ignore=ignore)
    shutil.copytree(
        repo_root / "tests/fixtures/trait_key_policy",
        scratch_repo / "tests/fixtures/trait_key_policy",
        ignore=ignore,
    )
    for filename in ("pyproject.toml", "Justfile", "uv.lock"):
        src_file = repo_root / filename
        if src_file.is_file():
            shutil.copy2(src_file, scratch_repo / filename)
    for test_rel in (
        "tests/test_arm_trait_key_policy.py",
        "tests/test_contribution_caller_policy.py",
    ):
        dst = scratch_repo / test_rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(repo_root / test_rel, dst)


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
        tempfile.mkdtemp(prefix=f"verify_2043_{run_label}_", dir=str(scratch_parent))
    ).resolve()
    scratch_repo = scratch_dir / "repo"
    _populate_scratch_repo(repo_root, scratch_repo)

    baseline_nodes = (
        "tests/test_arm_trait_key_policy.py",
        "tests/test_contribution_caller_policy.py",
    )
    for idx, node in enumerate(baseline_nodes, start=1):
        junit_base = scratch_dir / f"junit_baseline_{idx}.xml"
        code, inspector, out = _run_pytest_node(repo_root, scratch_repo, node, junit_base)
        if code != 0 or inspector.non_call_failures or not junit_base.is_file():
            raise RuntimeError(f"Baseline test {node} failed in scratch_repo (exit={code}):\n{out}")

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
