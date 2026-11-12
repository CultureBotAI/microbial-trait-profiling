"""Isolated scratch-copy mutant verifier for issue #833 (agy, 2026-11-09).

Purpose:
    Verify that tests/test_cellulose_significant_results_summary_833.py kills
    mutations of the pooled sign-test bound, study/experiment counts, sensitivity
    table values, and recount script logic in an isolated temporary scratch copy
    without modifying any worktree file.
Inputs:
    - data/analysis/cellulose_decr/* (copied to a temporary scratch directory)
    - scripts/adhoc/recount_cellulose_readme_numbers_799.py (copied to scratch)
    - tests/test_cellulose_significant_results_summary_833.py (copied to scratch)
Outputs:
    - Per-mutant verification summary printed to stdout; non-zero exit if any
      mutant survives, fails outside AssertionError, or if worktree files change.
Issue: #833
Date: 2026-11-09
"""

from __future__ import annotations

import hashlib
import importlib.util
import shutil
import sys
import tempfile
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
ANALYSIS_REL = Path("data/analysis/cellulose_decr")
SUMMARY_REL = ANALYSIS_REL / "SIGNIFICANT_RESULTS_SUMMARY.md"
README_REL = ANALYSIS_REL / "README.md"
SENSITIVITY_REL = ANALYSIS_REL / "condition_definition_sensitivity.tsv"
RECOUNT_REL = Path("scripts/adhoc/recount_cellulose_readme_numbers_799.py")
TEST_REL = Path("tests/test_cellulose_significant_results_summary_833.py")
LE = "\u2264"

MUTANTS: tuple[tuple[str, Path, str, str], ...] = (
    (
        "M01_summary_stale_oriented_only_bound",
        SUMMARY_REL,
        f"p {LE} 1.21e-4 under every label",
        f"p {LE} 8.8e-5 under every label",
    ),
    (
        "M02_summary_uncalibrated_loose_bound",
        SUMMARY_REL,
        f"p {LE} 1.21e-4 under every label",
        f"p {LE} 1.50e-4 under every label",
    ),
    (
        "M03_summary_stale_pooled_all_contrasts_counts",
        SUMMARY_REL,
        "56 of 77 untied studies (p = 8.2e-5)",
        "58 of 78 untied studies (p = 1.9e-5)",
    ),
    (
        "M04_summary_stale_pooled_oriented_counts",
        SUMMARY_REL,
        "49 of 63 among\n     disease-versus-healthy contrasts (p = 1.1e-5)",
        "49 of 64 among\n     disease-versus-healthy contrasts (p = 2.4e-5)",
    ),
    (
        "M05_summary_table_a_pooled_row_drift",
        SUMMARY_REL,
        "| All three | all | 107 : 37 : 30 | 4.4e-9 | 56 : 21 : 12 | 8.2e-5 |",
        "| All three | all | 107 : 37 : 30 | 4.4e-9 | 56 : 21 : 12 | 8.8e-5 |",
    ),
    (
        "M06_summary_table_a2_synonyms_row_drift",
        SUMMARY_REL,
        "| + synonym labels | 13 : 5 : 3 | 0.096 | 10 : 4 : 2 | 0.18 |",
        "| + synonym labels | 13 : 5 : 3 | 0.096 | 10 : 2 : 2 | 0.039 |",
    ),
    (
        "M07_readme_stale_pooled_bound",
        README_REL,
        f"stays at p {LE} 1.21e-4 in every scenario.",
        f"stays at p {LE} 8.8e-5 in every scenario.",
    ),
    (
        "M08_sensitivity_tsv_worst_p_drift",
        SENSITIVITY_REL,
        "55\t21\t12\t0.723684\t0.000120629",
        "55\t21\t12\t0.723684\t0.000250000",
    ),
    (
        "M09_recount_script_scope_restricted_to_oriented",
        RECOUNT_REL,
        'if k[2:] == ("Pooled", "study")',
        'if k[1:] == ("disease_vs_healthy_oriented", "Pooled", "study")',
    ),
)


class _OutcomeCapturePlugin:
    def __init__(self, expected_test_file: Path) -> None:
        self.expected_test_file = expected_test_file.resolve()
        self.passed = 0
        self.assertion_failures = 0
        self.other_failures: list[str] = []
        self.origin_checked = False

    def pytest_collection_modifyitems(self, items: list[pytest.Item]) -> None:
        if not items:
            raise RuntimeError("No pytest items collected")
        for item in items:
            item_path = Path(str(item.fspath)).resolve()
            if item_path != self.expected_test_file:
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
    shutil.copytree(REPO_ROOT / ANALYSIS_REL, scratch / ANALYSIS_REL)
    (scratch / RECOUNT_REL.parent).mkdir(parents=True, exist_ok=True)
    shutil.copy2(REPO_ROOT / RECOUNT_REL, scratch / RECOUNT_REL)
    (scratch / TEST_REL.parent).mkdir(parents=True, exist_ok=True)
    shutil.copy2(REPO_ROOT / TEST_REL, scratch / TEST_REL)


def _run_in_scratch(scratch: Path, k_expr: str) -> tuple[int, _OutcomeCapturePlugin]:
    for mod_name in list(sys.modules):
        if "test_cellulose_significant_results_summary_833" in mod_name or (
            "recount_cellulose_readme_numbers_799" in mod_name
        ):
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
                k_expr,
                str(test_path),
            ],
            plugins=[plugin],
        )
    )
    spec = importlib.util.spec_from_file_location("scratch_test_833", test_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Failed to load spec for {test_path}")
    test_mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(test_mod)
    if Path(test_mod.ANALYSIS_DIR).resolve() != (scratch / ANALYSIS_REL).resolve():
        raise RuntimeError(f"ANALYSIS_DIR origin mismatch: {test_mod.ANALYSIS_DIR}")
    if Path(test_mod.RECOUNT_SCRIPT).resolve() != (scratch / RECOUNT_REL).resolve():
        raise RuntimeError(f"RECOUNT_SCRIPT origin mismatch: {test_mod.RECOUNT_SCRIPT}")
    return rc, plugin


def main() -> int:
    tracked = (
        REPO_ROOT / SUMMARY_REL,
        REPO_ROOT / README_REL,
        REPO_ROOT / SENSITIVITY_REL,
        REPO_ROOT / RECOUNT_REL,
        REPO_ROOT / TEST_REL,
    )
    before = {path: _sha256(path) for path in tracked}
    originals = {
        rel: (REPO_ROOT / rel).read_text(encoding="utf-8")
        for rel in (SUMMARY_REL, README_REL, SENSITIVITY_REL, RECOUNT_REL)
    }

    with tempfile.TemporaryDirectory(prefix="verify_833_mutants_") as tmp_dir:
        scratch = Path(tmp_dir)
        _stage_scratch(scratch)
        ctrl_rc, ctrl_plugin = _run_in_scratch(
            scratch,
            "test_pooled_sign_test_bound_matches_worst_scenario_across_all_populations "
            "or test_significant_results_summary_validator_rejects_perturbed_bounds_and_counts",
        )
        print(
            f"CONTROL: rc={ctrl_rc} origin_checked={ctrl_plugin.origin_checked} "
            f"passed={ctrl_plugin.passed} assertion_failures={ctrl_plugin.assertion_failures} "
            f"other_failures={len(ctrl_plugin.other_failures)}"
        )
        if (
            ctrl_rc != 0
            or not ctrl_plugin.origin_checked
            or ctrl_plugin.passed != 8
            or ctrl_plugin.assertion_failures != 0
            or ctrl_plugin.other_failures
        ):
            return 1

        killed = 0
        for name, rel_path, target, replacement in MUTANTS:
            for rel, content in originals.items():
                (scratch / rel).write_text(content, encoding="utf-8")
            base_source = originals[rel_path]
            if target not in base_source:
                print(f"ERROR: target not found for {name} in {rel_path}")
                return 1
            mutated = base_source.replace(target, replacement, 1)
            (scratch / rel_path).write_text(mutated, encoding="utf-8")
            rc, plugin = _run_in_scratch(
                scratch,
                "test_pooled_sign_test_bound_matches_worst_scenario_across_all_populations",
            )
            is_killed = (
                rc != 0
                and plugin.origin_checked
                and plugin.assertion_failures == 1
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
