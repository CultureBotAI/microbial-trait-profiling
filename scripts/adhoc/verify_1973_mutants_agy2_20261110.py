"""Isolated scratch-copy mutation verifier for issue #1973.

Purpose:
    Verify that non-ASCII POSIX home account pattern regressions in
    src/mtp/local_path_patterns.py and src/mtp/provenance_paths.py are caught
    by AssertionError in tests/test_local_path_patterns.py without mutating
    the canonical repository worktree.
Inputs:
    Canonical repository files under src/mtp/ and tests/.
Outputs:
    JSON summary printed to stdout; exit code 0 iff all 7 mutants fail with
    AssertionError and canonical file digests remain unchanged.
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
        target=r"/(?:Users|home)/[\w.-]+",
        replacement=r"/(?:Users|home)/[A-Za-z0-9._-]+",
        test_node="test_posix_non_ascii_accounts_are_matched_in_full_across_both_consumers",
    ),
    MutantSpec(
        name="local_path_patterns_single_char_posix_account",
        rel_path="src/mtp/local_path_patterns.py",
        target=r"/(?:Users|home)/[\w.-]+",
        replacement=r"/(?:Users|home)/[\w.-]",
        test_node="test_posix_non_ascii_accounts_are_matched_in_full_across_both_consumers",
    ),
    MutantSpec(
        name="local_path_patterns_drop_home_root",
        rel_path="src/mtp/local_path_patterns.py",
        target=r"/(?:Users|home)/[\w.-]+",
        replacement=r"/(?:Users)/[\w.-]+",
        test_node="test_posix_non_ascii_accounts_are_matched_in_full_across_both_consumers",
    ),
    MutantSpec(
        name="local_path_patterns_drop_dot_hyphen_in_posix_account",
        rel_path="src/mtp/local_path_patterns.py",
        target=r"/(?:Users|home)/[\w.-]+",
        replacement=r"/(?:Users|home)/\w+",
        test_node="test_posix_non_ascii_accounts_are_matched_in_full_across_both_consumers",
    ),
    MutantSpec(
        name="provenance_slashless_ascii_only_posix_accounts",
        rel_path="src/mtp/provenance_paths.py",
        target=r"(?<![\w])(?<!\w-)(?:(?:Users|home)/[\w.-]+",
        replacement=r"(?<![\w])(?<!\w-)(?:(?:Users|home)/[A-Za-z0-9._-]+",
        test_node="test_slashless_posix_non_ascii_accounts_are_refused_in_provenance",
    ),
    MutantSpec(
        name="provenance_slashless_single_char_posix_account",
        rel_path="src/mtp/provenance_paths.py",
        target=r"(?<![\w])(?<!\w-)(?:(?:Users|home)/[\w.-]+",
        replacement=r"(?<![\w])(?<!\w-)(?:(?:Users|home)/[\w.-]",
        test_node="test_slashless_posix_non_ascii_accounts_are_refused_in_provenance",
    ),
    MutantSpec(
        name="provenance_slashless_drop_home_root",
        rel_path="src/mtp/provenance_paths.py",
        target=r"(?<![\w])(?<!\w-)(?:(?:Users|home)/[\w.-]+",
        replacement=r"(?<![\w])(?<!\w-)(?:(?:Users)/[\w.-]+",
        test_node="test_slashless_posix_non_ascii_accounts_are_refused_in_provenance",
    ),
)


class _AssertionCapturePlugin:
    """Capture failed test reports and exception types from in-process pytest."""

    def __init__(self) -> None:
        self.failures: list[tuple[str, str]] = []

    def pytest_runtest_logreport(self, report: pytest.TestReport) -> None:
        if report.when == "call" and report.failed:
            exc_name = "UnknownFailure"
            if report.longrepr is not None and hasattr(report.longrepr, "reprcrash"):
                reprcrash = getattr(report.longrepr, "reprcrash", None)
                if reprcrash is not None and hasattr(reprcrash, "message"):
                    msg = str(reprcrash.message)
                    exc_name = msg.split(":", 1)[0].strip()
            self.failures.append((report.nodeid, exc_name))


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


def _run_pytest_in_scratch(scratch_root: Path, test_node: str) -> tuple[int, list[tuple[str, str]]]:
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
        control_code, control_failures = _run_pytest_in_scratch(
            control_dir,
            "test_posix_non_ascii_accounts_are_matched_in_full_across_both_consumers",
        )
        if control_code != 0 or control_failures:
            raise RuntimeError(
                f"Unmutated control failed in scratch copy: code={control_code}, "
                f"failures={control_failures}"
            )

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
            assertion_killed = code != 0 and any(
                exc_type == "AssertionError" for _, exc_type in failures
            )
            if not assertion_killed:
                raise RuntimeError(
                    f"Mutant {spec.name} was not killed by AssertionError: "
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
        "all_killed_by_assertion": True,
        "canonical_digests_unchanged": True,
        "canonical_sha256": after_digests,
        "mutants": results,
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
