#!/usr/bin/env python3
"""Isolated scratch-copy mutation verifier for issue #748.

Purpose:
    Prove that unit tests in tests/test_profile_enrichment_guild_assignment_diagram.py
    pass on the baseline tree and fail with AssertionError for each of 7 distinct
    mutations to scripts/create_profile_enrichment_guild_assignment_diagram.py.
Inputs:
    scripts/create_profile_enrichment_guild_assignment_diagram.py
    tests/test_profile_enrichment_guild_assignment_diagram.py
Outputs:
    Per-mutant pass/kill summary printed to stdout; non-zero exit if baseline fails
    or any mutant survives or fails without AssertionError.
Issue:
    #748
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
REL_SCRIPT = Path("scripts/create_profile_enrichment_guild_assignment_diagram.py")
REL_TEST = Path("tests/test_profile_enrichment_guild_assignment_diagram.py")


@dataclass(frozen=True)
class MutantSpec:
    """Single source replacement mutation in the workflow diagram script."""

    name: str
    target: str
    replacement: str


MUTANTS: tuple[MutantSpec, ...] = (
    MutantSpec(
        name="revert_box2_to_before_human_host_plotting_filters",
        target='"trait columns before plotting trait filters."',
        replacement='"trait columns before human-host plotting filters."',
    ),
    MutantSpec(
        name="omit_human_host_species_in_box2",
        target=("f\"{count_label(counts, 'raw_profile_count')} canonical {HUMAN_HOST_SPECIES} \""),
        replacement="f\"{count_label(counts, 'raw_profile_count')} canonical \"",
    ),
    MutantSpec(
        name="omit_sequencing_types_in_box2",
        target='f"{SEQUENCING_TYPES_TEXT} profiles x "',
        replacement='"profiles x "',
    ),
    MutantSpec(
        name="swap_raw_profile_count_with_study_count_in_box2",
        target=("f\"{count_label(counts, 'raw_profile_count')} canonical {HUMAN_HOST_SPECIES} \""),
        replacement=("f\"{count_label(counts, 'study_count')} canonical {HUMAN_HOST_SPECIES} \""),
    ),
    MutantSpec(
        name="swap_raw_trait_column_count_with_selected_trait_count_in_box2",
        target=("f\"{count_label(counts, 'raw_trait_column_count')} non-empty KG-Microbe \""),
        replacement=("f\"{count_label(counts, 'selected_trait_count')} non-empty KG-Microbe \""),
    ),
    MutantSpec(
        name="omit_non_empty_trait_column_qualifier_in_box2",
        target=("f\"{count_label(counts, 'raw_trait_column_count')} non-empty KG-Microbe \""),
        replacement="f\"{count_label(counts, 'raw_trait_column_count')} KG-Microbe \"",
    ),
    MutantSpec(
        name="bypass_build_workflow_stages_in_create_workflow_diagram",
        target="    stages = build_workflow_stages(counts)",
        replacement=(
            "    stages = [\n"
            "        *build_workflow_stages(counts)[:1],\n"
            "        {\n"
            '            "number": "2",\n'
            '            "title": "Trait profiles",\n'
            '            "body": (\n'
            "                f\"{count_label(counts, 'raw_profile_count')} \"\n"
            '                "study/experiment/direction profiles x "\n'
            "                f\"{count_label(counts, 'raw_trait_column_count')} KG-Microbe \"\n"
            '                "trait columns before human-host plotting filters."\n'
            "            ),\n"
            '            "accent": GREEN_700,\n'
            '            "facecolor": GREEN_50,\n'
            '            "edgecolor": "#a7f3d0",\n'
            "        },\n"
            "        *build_workflow_stages(counts)[2:],\n"
            "    ]"
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
    for rel_path in (REL_SCRIPT, REL_TEST):
        dest = scratch_root / rel_path
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(REPO_ROOT / rel_path, dest)


def _clear_cached_modules() -> None:
    for mod_name in (
        "test_profile_enrichment_guild_assignment_diagram",
        "tests.test_profile_enrichment_guild_assignment_diagram",
        "create_profile_enrichment_guild_assignment_diagram_test",
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

    with tempfile.TemporaryDirectory(prefix="mtp-748-mutants-") as tmp_dir:
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
            script_path = mutant_root / REL_SCRIPT
            source_text = script_path.read_text(encoding="utf-8")
            if mutant.target not in source_text:
                raise RuntimeError(f"Mutant target not found for {mutant.name}")
            mutated_text = source_text.replace(mutant.target, mutant.replacement, 1)
            if mutated_text == source_text:
                raise RuntimeError(f"Mutant {mutant.name} produced no source change")
            script_path.write_text(mutated_text, encoding="utf-8")

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
