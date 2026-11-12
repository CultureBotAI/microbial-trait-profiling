#!/usr/bin/env python3
"""Isolated scratch-copy mutation verifier for issue #1682.

Purpose:
    Prove that unit tests in tests/test_no_absolute_local_paths.py pass on the
    baseline tree and fail with AssertionError for each of 7 distinct mutations
    to the #1682 excluded-path-family scope rationale and boundary pins.
Inputs:
    tests/test_no_absolute_local_paths.py
Outputs:
    Per-mutant pass/kill summary printed to stdout; non-zero exit if baseline fails
    or any mutant survives or fails without AssertionError.
Issue:
    #1682
Date:
    2026-11-10
"""

from __future__ import annotations

import shutil
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
REL_TEST = Path("tests/test_no_absolute_local_paths.py")
TARGET_FILTER = (
    "test_each_enforcing_pattern_matches_its_family_and_not_the_other or "
    "test_safe_and_out_of_scope_samples_match_no_pattern or "
    "test_excluded_path_families_pin_scope_rationale_and_boundary or "
    "test_user_home_after_a_temp_root_is_still_reported"
)


@dataclass(frozen=True)
class MutantSpec:
    """Single source replacement mutation in the local-path guard test module."""

    name: str
    target: str
    replacement: str


MUTANTS: tuple[MutantSpec, ...] = (
    MutantSpec(
        name="remove_excluded_path_families_reference_from_module_docstring",
        target="What it deliberately leaves out (#1682, pinned in ``EXCLUDED_PATH_FAMILIES``):",
        replacement="What it deliberately leaves out:",
    ),
    MutantSpec(
        name="drop_external_volume_from_excluded_path_families",
        target='    "external_volume": ExcludedPathFamily(',
        replacement='    "unregistered_volume": ExcludedPathFamily(',
    ),
    MutantSpec(
        name="drop_issue_1682_tag_from_bare_root_tmp_rationale",
        target=(
            '            "resolved temp paths under private/tmp and per-user var/folders '
            'remain forbidden "\n'
            '            "(#1682)."'
        ),
        replacement=(
            '            "resolved temp paths under private/tmp and per-user var/folders '
            'remain forbidden."'
        ),
    ),
    MutantSpec(
        name="desync_out_of_scope_samples_from_excluded_families",
        target=(
            "OUT_OF_SCOPE_SAMPLES: tuple[str, ...] = (\n"
            '    "/" + "tmp/mtp-run/out.json",\n'
            '    "~/" + ".cache/uv",\n'
            '    "/" + "Volumes/External/mtp/out.tsv",\n'
            ")"
        ),
        replacement=(
            "OUT_OF_SCOPE_SAMPLES: tuple[str, ...] = (\n"
            '    "/" + "tmp/mtp-run/out.json",\n'
            '    "~/" + "different/cache",\n'
            '    "/" + "Volumes/External/mtp/out.tsv",\n'
            ")"
        ),
    ),
    MutantSpec(
        name="contaminate_bare_root_tmp_representative_with_private_tmp",
        target='            "/" + "tmp/scratch/table.tsv",',
        replacement='            "/" + "private/tmp/scratch/table.tsv",',
    ),
    MutantSpec(
        name="weaken_tilde_relative_home_forbidden_counterpart_to_tilde",
        target='            (USER_HOME, "/" + "Users/alice/.cache/uv"),',
        replacement='            (USER_HOME, "~/" + "alice/.cache/uv"),',
    ),
    MutantSpec(
        name="misclassify_external_volume_temp_counterpart_as_user_home",
        target=(
            '            (TEMP_ROOT, "/" + "Volumes/External" + "/" + "private/tmp/mtp/out.tsv"),'
        ),
        replacement=(
            '            (USER_HOME, "/" + "Volumes/External" + "/" + "private/tmp/mtp/out.tsv"),'
        ),
    ),
)


class _ReportPlugin:
    """Collect pytest failure representations in-process."""

    def __init__(self) -> None:
        self.failure_reprs: list[str] = []

    def pytest_runtest_logreport(self, report: pytest.TestReport) -> None:
        if report.failed and report.longrepr is not None:
            self.failure_reprs.append(str(report.longrepr))


def _stage_scratch_tree(scratch_root: Path) -> None:
    dest = scratch_root / REL_TEST
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(REPO_ROOT / REL_TEST, dest)


def _clear_cached_modules() -> None:
    for mod_name in (
        "test_no_absolute_local_paths",
        "tests.test_no_absolute_local_paths",
        "conftest",
        "tests.conftest",
        "tests",
    ):
        sys.modules.pop(mod_name, None)


def _run_pytest_in_scratch(scratch_root: Path) -> tuple[int, list[str]]:
    plugin = _ReportPlugin()
    _clear_cached_modules()
    code = int(
        pytest.main(
            [
                "-q",
                "-o",
                "addopts=",
                "--tb=long",
                "-k",
                TARGET_FILTER,
                str(scratch_root / REL_TEST),
            ],
            plugins=[plugin],
        )
    )
    _clear_cached_modules()
    return code, plugin.failure_reprs


def main() -> int:
    src_dir = str(REPO_ROOT / "src")
    if src_dir not in sys.path:
        sys.path.insert(0, src_dir)

    with tempfile.TemporaryDirectory(prefix="mtp-1682-mutants-") as tmp_dir:
        tmp_root = Path(tmp_dir)
        baseline_root = tmp_root / "baseline"
        _stage_scratch_tree(baseline_root)
        base_code, base_failures = _run_pytest_in_scratch(baseline_root)
        if base_code != 0:
            raise RuntimeError(f"Baseline tests failed (exit {base_code}): {base_failures}")
        print("BASELINE: PASS")

        killed = 0
        for index, mutant in enumerate(MUTANTS, start=1):
            mutant_root = tmp_root / f"mutant_{index}_{mutant.name}"
            _stage_scratch_tree(mutant_root)
            test_path = mutant_root / REL_TEST
            source_text = test_path.read_text(encoding="utf-8")
            if mutant.target not in source_text:
                raise RuntimeError(f"Mutant target not found for {mutant.name}")
            mutated_text = source_text.replace(mutant.target, mutant.replacement, 1)
            if mutated_text == source_text:
                raise RuntimeError(f"Mutant {mutant.name} produced no source change")
            test_path.write_text(mutated_text, encoding="utf-8")

            code, failures = _run_pytest_in_scratch(mutant_root)
            combined = "\n".join(failures)
            if code == 0:
                raise RuntimeError(f"Mutant {mutant.name} SURVIVED")
            if "AssertionError" not in combined:
                raise RuntimeError(
                    f"Mutant {mutant.name} failed without AssertionError:\n{combined}"
                )
            killed += 1
            print(f"MUTANT {index} ({mutant.name}): KILLED by AssertionError")

        if killed < 6:
            raise RuntimeError(f"Expected >= 6 killed mutants, got {killed}")
        print(f"SUMMARY: {killed}/{len(MUTANTS)} mutants killed by AssertionError.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
