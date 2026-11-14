"""Isolated scratch-copy mutant verifier for issue #903 (agy, 2026-11-09).

Purpose:
    Verify that tests/test_human_host_heatmap_filters.py kills source mutants of
    the dump-absent separation and subset outside-corpus fail-closed checks in
    scripts/create_taxa_guild_heatmap.py without modifying files in the worktree.
Inputs:
    - scripts/create_taxa_guild_heatmap.py (copied to a temporary scratch directory)
    - tests/test_human_host_heatmap_filters.py (copied to a temporary scratch directory)
Outputs:
    - Per-mutant verification summary printed to stdout; non-zero exit if any
      mutant survives, fails outside AssertionError, or if worktree files change.
Issue: #903
Date: 2026-11-09
"""

from __future__ import annotations

import hashlib
import importlib.util
import os
import shutil
import sys
import tempfile
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
TARGET_REL = Path("scripts/create_taxa_guild_heatmap.py")
TEST_REL = Path("tests/test_human_host_heatmap_filters.py")

MUTANTS: tuple[tuple[str, str, str], ...] = (
    (
        "M01_conflate_missing_with_outside_corpus",
        "    outside_corpus = int((in_dump & ~in_corpus).sum())",
        "    outside_corpus = int((~in_corpus).sum())",
    ),
    (
        "M02_zero_missing_from_dump_count",
        "    missing_from_dump = int((~in_dump).sum())",
        "    missing_from_dump = 0",
    ),
    (
        "M03_omit_missing_from_dump_fail_closed",
        "    if missing_from_dump > 0:",
        "    if missing_from_dump > 1:",
    ),
    (
        "M04_omit_subset_outside_corpus_fail_closed",
        '    if direction_basis == "disease_vs_comparison" and outside_corpus > 0:',
        '    if False and direction_basis == "disease_vs_comparison" and outside_corpus > 0:',
    ),
    (
        "M05_require_more_than_one_outside_corpus_drop",
        '    if direction_basis == "disease_vs_comparison" and outside_corpus > 0:',
        '    if direction_basis == "disease_vs_comparison" and outside_corpus > 1:',
    ),
    (
        "M06_omit_summary_direction_basis_forwarding",
        (
            "    taxa = load_bugsigdb_taxa_records(\n"
            "        bugsigdb_taxa_path, full_dump_path, direction_basis=direction_basis\n"
            "    )"
        ),
        "    taxa = load_bugsigdb_taxa_records(bugsigdb_taxa_path, full_dump_path)",
    ),
    (
        "M07_omit_main_direction_basis_forwarding",
        "        args.full_dump_input,\n        direction_basis=args.direction_basis,",
        "        args.full_dump_input,",
    ),
    (
        "M08_omit_direction_basis_validation",
        "    if direction_basis not in DIRECTION_BASIS_LABELS:",
        "    if False and direction_basis not in DIRECTION_BASIS_LABELS:",
    ),
    (
        "M09_omit_audit_summary_attr",
        '    taxa.attrs["audit_summary"] = audit_summary',
        "",
    ),
)


class _OutcomeCapturePlugin:
    def __init__(self, expected_test_path: Path) -> None:
        self.expected_test_path = expected_test_path.resolve()
        self.origin_checked = False
        self.passed = 0
        self.assertion_failures = 0
        self.other_failures: list[str] = []

    def pytest_collection_modifyitems(self, items: list[pytest.Item]) -> None:
        if not items:
            raise RuntimeError("No pytest items collected")
        for item in items:
            item_path = Path(str(item.fspath)).resolve()
            if item_path != self.expected_test_path:
                raise RuntimeError(f"Unexpected test origin: {item_path}")
        self.origin_checked = True

    def pytest_runtest_logreport(self, report: pytest.TestReport) -> None:
        if report.when == "call" and report.passed:
            self.passed += 1
        elif report.failed:
            longrepr = getattr(report, "longrepr", None)
            crash = getattr(longrepr, "reprcrash", None)
            msg = getattr(crash, "message", str(longrepr))
            if report.when == "call" and (
                msg.startswith(("AssertionError", "Failed: DID NOT RAISE", "assert "))
                or "AssertionError" in msg
            ):
                self.assertion_failures += 1
            else:
                self.other_failures.append(f"{report.nodeid}:{report.when}:{msg}")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _stage_scratch(scratch: Path) -> None:
    for rel in (TARGET_REL, TEST_REL):
        dst = scratch / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(REPO_ROOT / rel, dst)


def _run_in_scratch(scratch: Path) -> tuple[int, _OutcomeCapturePlugin]:
    os.environ["GIT_CONFIG_GLOBAL"] = "/dev/null"
    os.environ["GIT_CONFIG_NOSYSTEM"] = "1"
    os.environ.setdefault("MPLCONFIGDIR", str(scratch / "mpl-cache"))

    for mod_name in list(sys.modules):
        if mod_name in {
            "create_taxa_guild_heatmap",
            "test_human_host_heatmap_filters",
        } or mod_name.endswith(".test_human_host_heatmap_filters"):
            sys.modules.pop(mod_name, None)

    src_dir = str((REPO_ROOT / "src").resolve())
    if src_dir not in sys.path:
        sys.path.insert(0, src_dir)

    test_path = (scratch / TEST_REL).resolve()
    plugin = _OutcomeCapturePlugin(test_path)
    rc = int(
        pytest.main(
            [
                "-q",
                "-o",
                "addopts=--import-mode=importlib",
                "-k",
                (
                    "test_taxa_guild_reports_disjoint_row_losses_and_name_fallback "
                    "or test_taxa_guild_reports_zero_fallback_and_nonzero_filter_loss "
                    "or test_taxa_guild_separates_and_rejects_dump_absent_signatures "
                    "or test_taxa_guild_subset_mode_fails_when_rows_are_dropped_outside_corpus "
                    "or test_taxa_guild_rejects_unknown_direction_basis"
                ),
                str(test_path),
            ],
            plugins=[plugin],
        )
    )
    spec = importlib.util.spec_from_file_location("scratch_test_mod", test_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Failed to load spec for {test_path}")
    test_mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(test_mod)
    loaded = test_mod.load_script_module("create_taxa_guild_heatmap.py")
    expected_mod = (scratch / TARGET_REL).resolve()
    if loaded.__file__ is None or Path(loaded.__file__).resolve() != expected_mod:
        raise RuntimeError(f"Module origin mismatch: {getattr(loaded, '__file__', None)}")
    return rc, plugin


def main() -> int:
    tracked = (REPO_ROOT / TARGET_REL, REPO_ROOT / TEST_REL)
    before = {path: _sha256(path) for path in tracked}
    base_source = (REPO_ROOT / TARGET_REL).read_text(encoding="utf-8")

    with tempfile.TemporaryDirectory(prefix="verify_903_mutants_") as tmp_dir:
        scratch = Path(tmp_dir)
        _stage_scratch(scratch)
        ctrl_rc, ctrl_plugin = _run_in_scratch(scratch)
        print(
            f"CONTROL: rc={ctrl_rc} origin_checked={ctrl_plugin.origin_checked} "
            f"passed={ctrl_plugin.passed} assertion_failures={ctrl_plugin.assertion_failures} "
            f"other_failures={len(ctrl_plugin.other_failures)}"
        )
        if (
            ctrl_rc != 0
            or not ctrl_plugin.origin_checked
            or ctrl_plugin.assertion_failures != 0
            or ctrl_plugin.other_failures
        ):
            return 1

        killed = 0
        for name, target, replacement in MUTANTS:
            if target not in base_source:
                print(f"ERROR: target not found for {name}")
                return 1
            mutated = base_source.replace(target, replacement, 1)
            (scratch / TARGET_REL).write_text(mutated, encoding="utf-8")
            rc, plugin = _run_in_scratch(scratch)
            is_killed = (
                rc != 0
                and plugin.origin_checked
                and plugin.assertion_failures > 0
                and not plugin.other_failures
            )
            status = "KILLED" if is_killed else "SURVIVED_OR_ERROR"
            print(
                f"{name}: {status} (rc={rc}) origin_checked={plugin.origin_checked} "
                f"passed={plugin.passed} assertion_failures={plugin.assertion_failures} "
                f"other_failures={len(plugin.other_failures)}"
            )
            if not is_killed:
                for item in plugin.other_failures:
                    print("OTHER_FAILURE:", item)
                return 1
            killed += 1

    after = {path: _sha256(path) for path in tracked}
    if before != after:
        print("ERROR: worktree file hashes changed during verification!")
        return 1

    print(f"ALL {killed}/{len(MUTANTS)} MUTANTS KILLED BY ASSERTIONERROR; WORKTREE UNTOUCHED.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
