"""Isolated scratch-copy mutation verifier for issue #1973.

Purpose:
    Verify that non-ASCII POSIX home account pattern regressions (both literal
    UTF-8 and default JSON \\uXXXX escapes) in src/mtp/local_path_patterns.py
    and src/mtp/provenance_paths.py are caught by AssertionError in
    tests/test_local_path_patterns.py without mutating the canonical
    repository worktree.
Inputs:
    Canonical repository files under src/mtp/ and tests/.
Outputs:
    JSON summary printed to stdout; exit code 0 iff all 9 mutants fail with
    call-phase AssertionError, all selected test nodes pass unmutated controls,
    the broken-control self-check is rejected, and canonical file digests
    remain unchanged.
Issue:
    #1973
Date:
    2026-11-10
"""

from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import shutil
import sys
import tempfile
from collections.abc import Sequence
from contextlib import ExitStack
from dataclasses import dataclass
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
WATCHED_FILES = (
    "src/mtp/local_path_patterns.py",
    "src/mtp/provenance_paths.py",
    "tests/test_local_path_patterns.py",
)


@dataclass(frozen=True)
class MutantSpec:
    """Single source-text mutation against a target relative path."""

    name: str
    rel_path: str
    target: str
    replacement: str
    test_node: str


MUTANTS: tuple[MutantSpec, ...] = (
    MutantSpec(
        name="local_path_patterns_ascii_only_posix_accounts",
        rel_path="src/mtp/local_path_patterns.py",
        target=r"/(?:Users|home)/(?:[\w.-]|\\u[0-9A-Fa-f]{4})+",
        replacement=r"/(?:Users|home)/[A-Za-z0-9._-]+",
        test_node="test_posix_non_ascii_accounts_are_matched_in_full_across_both_consumers",
    ),
    MutantSpec(
        name="local_path_patterns_drop_json_unicode_escape",
        rel_path="src/mtp/local_path_patterns.py",
        target=r"/(?:Users|home)/(?:[\w.-]|\\u[0-9A-Fa-f]{4})+",
        replacement=r"/(?:Users|home)/[\w.-]+",
        test_node="test_posix_non_ascii_accounts_are_matched_in_full_across_both_consumers",
    ),
    MutantSpec(
        name="local_path_patterns_single_char_posix_account",
        rel_path="src/mtp/local_path_patterns.py",
        target=r"/(?:Users|home)/(?:[\w.-]|\\u[0-9A-Fa-f]{4})+",
        replacement=r"/(?:Users|home)/(?:[\w.-]|\\u[0-9A-Fa-f]{4})",
        test_node="test_posix_non_ascii_accounts_are_matched_in_full_across_both_consumers",
    ),
    MutantSpec(
        name="local_path_patterns_drop_home_root",
        rel_path="src/mtp/local_path_patterns.py",
        target=r"/(?:Users|home)/(?:[\w.-]|\\u[0-9A-Fa-f]{4})+",
        replacement=r"/(?:Users)/(?:[\w.-]|\\u[0-9A-Fa-f]{4})+",
        test_node="test_posix_non_ascii_accounts_are_matched_in_full_across_both_consumers",
    ),
    MutantSpec(
        name="local_path_patterns_drop_dot_hyphen_in_posix_account",
        rel_path="src/mtp/local_path_patterns.py",
        target=r"/(?:Users|home)/(?:[\w.-]|\\u[0-9A-Fa-f]{4})+",
        replacement=r"/(?:Users|home)/(?:\w|\\u[0-9A-Fa-f]{4})+",
        test_node="test_posix_non_ascii_accounts_are_matched_in_full_across_both_consumers",
    ),
    MutantSpec(
        name="provenance_slashless_ascii_only_posix_accounts",
        rel_path="src/mtp/provenance_paths.py",
        target=r"(?<![\w])(?<!\w-)(?:(?:Users|home)/(?:[\w.-]|\\u[0-9A-Fa-f]{4})+",
        replacement=r"(?<![\w])(?<!\w-)(?:(?:Users|home)/[A-Za-z0-9._-]+",
        test_node="test_slashless_posix_non_ascii_accounts_are_refused_in_provenance",
    ),
    MutantSpec(
        name="provenance_slashless_drop_json_unicode_escape",
        rel_path="src/mtp/provenance_paths.py",
        target=r"(?<![\w])(?<!\w-)(?:(?:Users|home)/(?:[\w.-]|\\u[0-9A-Fa-f]{4})+",
        replacement=r"(?<![\w])(?<!\w-)(?:(?:Users|home)/[\w.-]+",
        test_node="test_slashless_posix_non_ascii_accounts_are_refused_in_provenance",
    ),
    MutantSpec(
        name="provenance_slashless_single_char_posix_account",
        rel_path="src/mtp/provenance_paths.py",
        target=r"(?<![\w])(?<!\w-)(?:(?:Users|home)/(?:[\w.-]|\\u[0-9A-Fa-f]{4})+",
        replacement=r"(?<![\w])(?<!\w-)(?:(?:Users|home)/(?:[\w.-]|\\u[0-9A-Fa-f]{4})",
        test_node="test_slashless_posix_non_ascii_accounts_are_refused_in_provenance",
    ),
    MutantSpec(
        name="provenance_slashless_drop_home_root",
        rel_path="src/mtp/provenance_paths.py",
        target=r"(?<![\w])(?<!\w-)(?:(?:Users|home)/(?:[\w.-]|\\u[0-9A-Fa-f]{4})+",
        replacement=r"(?<![\w])(?<!\w-)(?:(?:Users)/(?:[\w.-]|\\u[0-9A-Fa-f]{4})+",
        test_node="test_slashless_posix_non_ascii_accounts_are_refused_in_provenance",
    ),
)


class _AssertionCapturePlugin:
    """Capture failed test reports, phases, and exception types from in-process pytest."""

    def __init__(self) -> None:
        self.failures: list[tuple[str, str, str]] = []

    def pytest_runtest_logreport(self, report: pytest.TestReport) -> None:
        if report.failed:
            exc_name = "UnknownFailure"
            if report.longrepr is not None and hasattr(report.longrepr, "reprcrash"):
                reprcrash = getattr(report.longrepr, "reprcrash", None)
                if reprcrash is not None and hasattr(reprcrash, "message"):
                    msg = str(reprcrash.message)
                    exc_name = msg.split(":", 1)[0].strip()
            self.failures.append((report.nodeid, report.when, exc_name))

    def pytest_collectreport(self, report: pytest.CollectReport) -> None:
        if report.failed:
            self.failures.append((str(report.nodeid), "collect", "CollectFailure"))


def _sha256_map(root: Path) -> dict[str, str]:
    digests: dict[str, str] = {}
    for rel in WATCHED_FILES:
        digests[rel] = hashlib.sha256((root / rel).read_bytes()).hexdigest()
    return digests


def _prepare_scratch_tree(scratch_root: Path) -> None:
    shutil.copytree(REPO_ROOT / "src" / "mtp", scratch_root / "src" / "mtp")
    tests_dir = scratch_root / "tests"
    tests_dir.mkdir(parents=True, exist_ok=True)
    for test_rel in ("test_local_path_patterns.py", "test_no_absolute_local_paths.py"):
        shutil.copy2(REPO_ROOT / "tests" / test_rel, tests_dir / test_rel)


def _purge_loaded_modules() -> None:
    for mod_name in list(sys.modules):
        if (
            mod_name == "mtp"
            or mod_name.startswith("mtp.")
            or mod_name in {"test_no_absolute_local_paths", "test_local_path_patterns"}
        ):
            sys.modules.pop(mod_name, None)
    importlib.invalidate_caches()


def _run_pytest_in_scratch(
    scratch_root: Path, test_node: str
) -> tuple[int, list[tuple[str, str, str]]]:
    plugin = _AssertionCapturePlugin()
    old_sys_path = list(sys.path)
    sys.path.insert(0, str(scratch_root / "tests"))
    sys.path.insert(0, str(scratch_root / "src"))
    _purge_loaded_modules()
    try:
        test_target = f"{scratch_root / 'tests' / 'test_local_path_patterns.py'}::{test_node}"
        code = int(
            pytest.main(
                [
                    "-q",
                    "-o",
                    "addopts=",
                    test_target,
                ],
                plugins=[plugin],
            )
        )
    finally:
        _purge_loaded_modules()
        sys.path[:] = old_sys_path
    return code, plugin.failures


def _verify_unmutated_controls(control_dir: Path, selected_nodes: Sequence[str]) -> None:
    for node in selected_nodes:
        control_code, control_failures = _run_pytest_in_scratch(control_dir, node)
        if control_code != 0 or control_failures:
            raise RuntimeError(
                f"Unmutated control failed for {node}: code={control_code}, "
                f"failures={control_failures}"
            )


def _verify_broken_control_refused(broken_dir: Path) -> None:
    test_file = broken_dir / "tests" / "test_local_path_patterns.py"
    text = test_file.read_text(encoding="utf-8")
    anchor = (
        "def test_slashless_posix_non_ascii_accounts_are_refused_in_provenance(\n"
        "    root: str, account: str, ascii_escaped: bool\n"
        ") -> None:\n"
    )
    if anchor not in text:
        raise RuntimeError("Slashless test anchor not found for broken-control self-check")
    test_file.write_text(
        text.replace(anchor, anchor + '    raise AssertionError("unrelated control failure")\n', 1),
        encoding="utf-8",
    )
    try:
        _verify_unmutated_controls(
            broken_dir,
            ("test_slashless_posix_non_ascii_accounts_are_refused_in_provenance",),
        )
    except RuntimeError:
        return
    raise RuntimeError("Broken slashless control was not refused by _verify_unmutated_controls")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Verify #1973 mutants in isolated scratch copies.")
    parser.add_argument(
        "--scratch-dir",
        type=Path,
        default=None,
        help="Optional scratch directory root for isolated mutant copies.",
    )
    args = parser.parse_args(argv)

    before_digests = _sha256_map(REPO_ROOT)
    for spec in MUTANTS:
        source_text = (REPO_ROOT / spec.rel_path).read_text(encoding="utf-8")
        if spec.target not in source_text:
            raise RuntimeError(
                f"Preflight failed: mutant {spec.name} target missing in {spec.rel_path}"
            )

    selected_nodes = tuple(dict.fromkeys(spec.test_node for spec in MUTANTS))
    results: list[dict[str, object]] = []

    with ExitStack() as stack:
        if args.scratch_dir is not None:
            args.scratch_dir.mkdir(parents=True, exist_ok=True)
            base_dir = Path(
                stack.enter_context(
                    tempfile.TemporaryDirectory(
                        prefix="verify_1973_mutants_", dir=str(args.scratch_dir)
                    )
                )
            )
        else:
            base_dir = Path(
                stack.enter_context(tempfile.TemporaryDirectory(prefix="verify_1973_mutants_"))
            )

        control_dir = base_dir / "control"
        _prepare_scratch_tree(control_dir)
        _verify_unmutated_controls(control_dir, selected_nodes)

        broken_control_dir = base_dir / "broken_control"
        _prepare_scratch_tree(broken_control_dir)
        _verify_broken_control_refused(broken_control_dir)

        for spec in MUTANTS:
            mutant_dir = base_dir / spec.name
            _prepare_scratch_tree(mutant_dir)
            target_file = mutant_dir / spec.rel_path
            original_text = target_file.read_text(encoding="utf-8")
            if spec.target not in original_text:
                raise RuntimeError(f"Mutant {spec.name} target not found in {spec.rel_path}")
            mutated_text = original_text.replace(spec.target, spec.replacement, 1)
            if mutated_text == original_text:
                raise RuntimeError(f"Mutant {spec.name} produced no change in {spec.rel_path}")
            target_file.write_text(mutated_text, encoding="utf-8")

            code, failures = _run_pytest_in_scratch(mutant_dir, spec.test_node)
            assertion_killed = (
                code == int(pytest.ExitCode.TESTS_FAILED)
                and bool(failures)
                and all(
                    when == "call" and exc_type == "AssertionError"
                    for _, when, exc_type in failures
                )
            )
            if not assertion_killed:
                raise RuntimeError(
                    f"Mutant {spec.name} was not cleanly killed by call-phase AssertionError: "
                    f"code={code}, failures={failures}"
                )
            results.append(
                {
                    "mutant": spec.name,
                    "rel_path": spec.rel_path,
                    "test_node": spec.test_node,
                    "exit_code": code,
                    "failures": failures,
                    "killed_by_assertion": assertion_killed,
                }
            )

    after_digests = _sha256_map(REPO_ROOT)
    if before_digests != after_digests:
        raise RuntimeError(
            f"Canonical digests changed during mutation run: {before_digests} != {after_digests}"
        )

    summary = {
        "issue": 1973,
        "mutant_count": len(results),
        "selected_control_nodes": list(selected_nodes),
        "broken_control_refused": True,
        "all_killed_by_assertion": True,
        "canonical_digests_unchanged": True,
        "canonical_sha256": after_digests,
        "mutants": results,
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
