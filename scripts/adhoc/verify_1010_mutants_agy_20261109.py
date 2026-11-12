"""Isolated scratch-copy mutant verifier for issue #1010 (agy, 2026-11-09).

Purpose:
    Verify that tests/test_taxa_coverage_taxonomy.py kills source mutants of the
    NCBITaxon:12908 / NCBITaxon:2787823 unclassified-sequences placeholder lineage
    classification and diagnosis rules in scripts/create_taxa_coverage_histogram.py
    without ever modifying files in the live worktree.
Inputs:
    - scripts/create_taxa_coverage_histogram.py (copied to a temporary scratch directory)
    - tests/test_taxa_coverage_histogram.py (copied to a temporary scratch directory)
    - tests/test_taxa_coverage_taxonomy.py (copied to a temporary scratch directory)
Outputs:
    - Per-mutant verification summary printed to stdout; non-zero exit if any
      mutant survives, fails outside AssertionError, or if worktree files change.
Issue: #1010
Date: 2026-11-09
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
TARGET_REL = Path("scripts/create_taxa_coverage_histogram.py")
HISTOGRAM_TEST_REL = Path("tests/test_taxa_coverage_histogram.py")
TAXONOMY_TEST_REL = Path("tests/test_taxa_coverage_taxonomy.py")
FIXTURE_REL = Path("tests/fixtures/coverage/taxonomy_rank_class.tsv")

MUTANTS: tuple[tuple[str, str, str], ...] = (
    (
        "M01_omit_12908_root",
        '    "NCBITaxon:12908": "unclassified_sequences",\n',
        "",
    ),
    (
        "M02_omit_2787823_root",
        '    "NCBITaxon:2787823": "unclassified_sequences",\n',
        "",
    ),
    (
        "M03_truncate_multiple_same_group_roots",
        '        return TaxonomyState(TAXONOMY_ROOTS[roots[0]], "resolved_root", roots)',
        '        return TaxonomyState(TAXONOMY_ROOTS[roots[0]], "resolved_root", roots[:1])',
    ),
    (
        "M04_allow_conflicting_root_groups",
        "        if len(groups) > 1:",
        "        if False:",
    ),
    (
        "M05_omit_unclassified_sequences_diagnosis_branch",
        '    if taxonomy.group == "unclassified_sequences":',
        "    if False:",
    ),
    (
        "M06_restore_exact_name_list_fallback",
        '    if taxonomy.group == "unclassified_sequences":',
        (
            '    if taxonomy.group == "unclassified_sequences" or lower_name in {'
            '"uncultured organism", "unidentified", "unclassified", "unknown"}:'
        ),
    ),
    (
        "M07_restore_dead_gut_metagenome_regex",
        r'        or re.search(r"\b(cag:\d+|uba\d+|ucg-\d+|rf39)\b", lower_name)',
        r'        or re.search(r"\b(cag:\d+|uba\d+|ucg-\d+|rf39|gut metagenome)\b", lower_name)',
    ),
    (
        "M08_use_missing_lineage_remediation_for_placeholder",
        '            "an unresolved label does not establish a multi-taxon composite",',
        '            "missing lineage does not establish a multi-taxon composite",',
    ),
    (
        "M09_skip_no_recognized_root_path_in_unresolved",
        '    if taxonomy.group == "unresolved":',
        '    if taxonomy.status in {"not_in_taxonomy_edges", "non_ncbi_identifier"}:',
    ),
)


class _OutcomeCapturePlugin:
    def __init__(self, expected_test_path: Path, expected_mod_path: Path) -> None:
        self.expected_test_path = expected_test_path.resolve()
        self.expected_mod_path = expected_mod_path.resolve()
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
        hist_mod = sys.modules.get("test_taxa_coverage_histogram")
        loaded = getattr(hist_mod, "coverage", None) if hist_mod is not None else None
        if loaded is None or loaded.__file__ is None:
            raise RuntimeError("create_taxa_coverage_histogram not loaded before test execution")
        loaded_path = Path(loaded.__file__).resolve()
        if loaded_path != self.expected_mod_path:
            raise RuntimeError(f"Module origin mismatch: {loaded_path} != {self.expected_mod_path}")
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
    for rel in (TARGET_REL, HISTOGRAM_TEST_REL, TAXONOMY_TEST_REL, FIXTURE_REL):
        dst = scratch / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(REPO_ROOT / rel, dst)


def _run_in_scratch(scratch: Path) -> tuple[int, _OutcomeCapturePlugin]:
    os.environ["GIT_CONFIG_GLOBAL"] = "/dev/null"
    os.environ["GIT_CONFIG_NOSYSTEM"] = "1"
    os.environ.setdefault("MPLCONFIGDIR", str(scratch / "mpl-cache"))

    for mod_name in list(sys.modules):
        if mod_name in {
            "create_taxa_coverage_histogram",
            "test_taxa_coverage_histogram",
            "test_taxa_coverage_taxonomy",
        } or mod_name.endswith((".test_taxa_coverage_histogram", ".test_taxa_coverage_taxonomy")):
            sys.modules.pop(mod_name, None)

    tests_dir = str((scratch / "tests").resolve())
    src_dir = str((REPO_ROOT / "src").resolve())
    sys.path = [p for p in sys.path if p not in {tests_dir, src_dir}]
    sys.path.insert(0, src_dir)
    sys.path.insert(0, tests_dir)

    hist_test_mod = importlib.import_module("test_taxa_coverage_histogram")
    coverage_mod = hist_test_mod.coverage
    expected_mod_path = (scratch / TARGET_REL).resolve()
    actual_mod_path = Path(str(coverage_mod.__file__)).resolve()
    if actual_mod_path != expected_mod_path:
        raise RuntimeError(f"Expected {expected_mod_path}, got {actual_mod_path}")
    coverage_mod.REPO_ROOT = REPO_ROOT

    test_path = (scratch / TAXONOMY_TEST_REL).resolve()
    plugin = _OutcomeCapturePlugin(test_path, expected_mod_path)
    rc = int(
        pytest.main(
            [
                "-q",
                "-o",
                "addopts=--import-mode=importlib",
                str(test_path),
            ],
            plugins=[plugin],
        )
    )
    return rc, plugin


def main() -> int:
    tracked = (
        REPO_ROOT / TARGET_REL,
        REPO_ROOT / HISTOGRAM_TEST_REL,
        REPO_ROOT / TAXONOMY_TEST_REL,
    )
    before = {path: _sha256(path) for path in tracked}
    base_source = (REPO_ROOT / TARGET_REL).read_text(encoding="utf-8")

    with tempfile.TemporaryDirectory(prefix="verify_1010_mutants_") as tmp_dir:
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
