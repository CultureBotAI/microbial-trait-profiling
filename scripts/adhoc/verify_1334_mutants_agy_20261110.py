"""Isolated scratch-copy mutant verifier for issue #1334 (2026-11-10).

Purpose:
    Verify in an isolated temporary copy (never mutating the live worktree) that
    the regression tests in tests/test_native_taxon_influence.py kill at least 9
    mutants of src/mtp/native_taxon_influence.py and methods/multiple_testing.md
    via AssertionError.
Inputs:
    Worktree files src/mtp/*.py, tests/test_native_taxon_influence.py, and
    methods/multiple_testing.md.
Outputs:
    Per-mutant status table printed to stdout; non-zero exit if control fails,
    any mutant survives or fails without AssertionError, or any worktree file changes.
Issue:
    #1334
Date:
    2026-11-10
"""

from __future__ import annotations

import hashlib
import importlib
import os
import shutil
import sys
import tempfile
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
TARGET_FILES = (
    Path("src/mtp/native_taxon_influence.py"),
    Path("tests/test_native_taxon_influence.py"),
    Path("methods/multiple_testing.md"),
)

TEST_NAMES = (
    "test_constant_paired_effect_is_not_a_significant_native_result",
    "test_native_exports_end_to_end_with_an_empty_stratum",
    "test_unselected_but_testable_traits_are_excluded_from_native_bh_family",
    "test_selected_untestable_traits_are_padded_at_one_in_native_bh_family",
    "test_multiple_testing_methods_doc_covers_native_influence_and_pd_network_bh_families",
)

MUTANTS: tuple[tuple[str, Path, str, str], ...] = (
    (
        "M01_bh_over_tested_only",
        Path("src/mtp/native_taxon_influence.py"),
        '    selected_rows = [r for r in records.values() if r["selected"]]',
        '    selected_rows = [r for r in records.values() if r["p"] != ""]',
    ),
    (
        "M02_bh_over_all_records",
        Path("src/mtp/native_taxon_influence.py"),
        '    selected_rows = [r for r in records.values() if r["selected"]]',
        "    selected_rows = list(records.values())",
    ),
    (
        "M03_pad_untestable_with_zero",
        Path("src/mtp/native_taxon_influence.py"),
        'adjust([r["p"] if r["p"] != "" else 1.0 for r in selected_rows])',
        'adjust([r["p"] if r["p"] != "" else 0.0 for r in selected_rows])',
    ),
    (
        "M04_pad_untestable_with_half",
        Path("src/mtp/native_taxon_influence.py"),
        'adjust([r["p"] if r["p"] != "" else 1.0 for r in selected_rows])',
        'adjust([r["p"] if r["p"] != "" else 0.5 for r in selected_rows])',
    ),
    (
        "M05_family_tests_uses_selected_len",
        Path("src/mtp/native_taxon_influence.py"),
        '        row["family_tests"] = len(tested)',
        '        row["family_tests"] = len(selected_rows)',
    ),
    (
        "M06_family_selected_uses_tested_len",
        Path("src/mtp/native_taxon_influence.py"),
        '        row["family_selected_traits"] = len(selected_rows)',
        '        row["family_selected_traits"] = len(tested)',
    ),
    (
        "M07_summary_says_outside_vector",
        Path("src/mtp/native_taxon_influence.py"),
        '"untestable traits at p = 1.0 (q = 1.0) inside the correction vector. "',
        '"untestable traits receive q = 1.0 outside the correction vector. "',
    ),
    (
        "M08_methods_doc_omits_pd_network_rule",
        Path("methods/multiple_testing.md"),
        "passes_min_taxon_studies",
        "omitted_taxon_support_rule",
    ),
    (
        "M09_methods_doc_omits_shared_sign_test",
        Path("methods/multiple_testing.md"),
        "correct_shared_comparisons",
        "omitted_shared_comparisons",
    ),
)


class _OutcomeCapturePlugin:
    def __init__(self, expected_test_file: Path, expected_module_file: Path) -> None:
        self.expected_test_file = expected_test_file.resolve()
        self.expected_module_file = expected_module_file.resolve()
        self.passed = 0
        self.assertion_failures = 0
        self.other_failures: list[str] = []
        self.origin_checked = False

    def pytest_collection_modifyitems(self, items: list[pytest.Item]) -> None:
        for item in items:
            item_path = Path(str(item.fspath)).resolve()
            if item_path != self.expected_test_file:
                raise RuntimeError(f"Unexpected test origin: {item_path}")
        mod = sys.modules.get("mtp.native_taxon_influence")
        if mod is None or mod.__file__ is None:
            raise RuntimeError("mtp.native_taxon_influence not loaded during collection")
        if Path(mod.__file__).resolve() != self.expected_module_file:
            raise RuntimeError(f"Unexpected module origin: {mod.__file__}")
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


def _snapshot_hashes(root: Path) -> dict[Path, str]:
    return {rel: _sha256(root / rel) for rel in TARGET_FILES}


def _populate_scratch(scratch: Path) -> None:
    shutil.copytree(REPO_ROOT / "src" / "mtp", scratch / "src" / "mtp")
    (scratch / "tests").mkdir(parents=True, exist_ok=True)
    shutil.copy2(
        REPO_ROOT / "tests" / "test_native_taxon_influence.py",
        scratch / "tests" / "test_native_taxon_influence.py",
    )
    (scratch / "methods").mkdir(parents=True, exist_ok=True)
    shutil.copy2(
        REPO_ROOT / "methods" / "multiple_testing.md",
        scratch / "methods" / "multiple_testing.md",
    )
    shutil.copy2(
        REPO_ROOT / "methods" / "table1_native_influence_protocol.json",
        scratch / "methods" / "table1_native_influence_protocol.json",
    )


def _run_in_scratch(scratch: Path) -> tuple[int, _OutcomeCapturePlugin]:
    os.environ["GIT_CONFIG_GLOBAL"] = "/dev/null"
    os.environ["GIT_CONFIG_NOSYSTEM"] = "1"
    scratch_src = str((scratch / "src").resolve())
    repo_src = str((REPO_ROOT / "src").resolve())
    sys.path[:] = [p for p in sys.path if Path(p).resolve() != Path(repo_src)]
    if scratch_src not in sys.path:
        sys.path.insert(0, scratch_src)

    for mod_name in list(sys.modules):
        if (
            mod_name == "mtp"
            or mod_name.startswith("mtp.")
            or mod_name == "test_native_taxon_influence"
            or mod_name.endswith(".test_native_taxon_influence")
        ):
            sys.modules.pop(mod_name, None)
    importlib.invalidate_caches()

    test_path = (scratch / "tests" / "test_native_taxon_influence.py").resolve()
    mod_path = (scratch / "src" / "mtp" / "native_taxon_influence.py").resolve()
    plugin = _OutcomeCapturePlugin(test_path, mod_path)
    rc = int(
        pytest.main(
            [
                "-q",
                "-o",
                "addopts=--import-mode=importlib",
                "-k",
                " or ".join(TEST_NAMES),
                str(test_path),
            ],
            plugins=[plugin],
        )
    )
    return rc, plugin


def main() -> int:
    before = _snapshot_hashes(REPO_ROOT)
    with tempfile.TemporaryDirectory(prefix="verify_1334_mutants_") as tmp_dir:
        scratch = Path(tmp_dir)
        _populate_scratch(scratch)
        rc, control = _run_in_scratch(scratch)
        if (
            rc != 0
            or not control.origin_checked
            or control.passed != len(TEST_NAMES)
            or control.assertion_failures != 0
            or control.other_failures
        ):
            print("CONTROL FAILED:", rc, control.passed, control.other_failures)
            return 1
        print(f"CONTROL: rc=0 PYTEST_ORIGIN_OK ({control.passed} passed)")

        killed = 0
        for name, rel_path, target, replacement in MUTANTS:
            file_path = scratch / rel_path
            original = file_path.read_text(encoding="utf-8")
            if original.count(target) != 1:
                print(f"{name}: target count {original.count(target)} != 1 in {rel_path}")
                return 1
            file_path.write_text(original.replace(target, replacement, 1), encoding="utf-8")
            try:
                m_rc, outcome = _run_in_scratch(scratch)
            finally:
                file_path.write_text(original, encoding="utf-8")
            if (
                m_rc != 0
                and outcome.origin_checked
                and outcome.assertion_failures >= 1
                and not outcome.other_failures
            ):
                killed += 1
                print(
                    f"{name}: KILLED (rc={m_rc}, assertion_failures={outcome.assertion_failures})"
                )
            else:
                print(
                    f"{name}: SURVIVED/OTHER (rc={m_rc}, "
                    f"assertions={outcome.assertion_failures}, other={outcome.other_failures})"
                )
                return 1

    after = _snapshot_hashes(REPO_ROOT)
    if before != after:
        print("WORKTREE HASH MISMATCH:", before, after)
        return 1
    print(f"SUMMARY: {killed}/{len(MUTANTS)} mutants killed by AssertionError; worktree untouched.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
