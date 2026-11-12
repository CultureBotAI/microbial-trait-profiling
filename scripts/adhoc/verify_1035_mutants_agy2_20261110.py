#!/usr/bin/env python3
"""Verify #1035 tests kill palette-collision and table-derivation mutants on scratch copies.
Inputs:
  - Repository checkout containing scripts/create_taxa_coverage_histogram.py,
    scripts/generate_tufte_plates.py, scripts/build_slide_deck.py,
    tests/test_taxa_coverage_histogram.py, and tests/test_deck_builders.py.
Outputs:
  - Exit 0 when all 6 mutants fail with AssertionError and canonical files remain untouched.
Issue: #1035
Date: 2026-11-10
"""

from __future__ import annotations

import contextlib
import hashlib
import io
import os
import shutil
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
WATCHED_FILES = (
    "scripts/create_taxa_coverage_histogram.py",
    "scripts/generate_tufte_plates.py",
    "scripts/build_slide_deck.py",
    "tests/test_taxa_coverage_histogram.py",
    "tests/test_deck_builders.py",
)


@dataclass(frozen=True)
class MutantSpec:
    mutant_id: str
    rel_path: str
    old_text: str
    new_text: str
    pytest_target: str
    pytest_k: str


MUTANTS: tuple[MutantSpec, ...] = (
    MutantSpec(
        mutant_id="M01_viral_color_collides_with_candidatus",
        rel_path="scripts/generate_tufte_plates.py",
        old_text='"Non-Bacterial / Viral Taxon": "#3f6b4e",',
        new_text='"Non-Bacterial / Viral Taxon": "#157a6e",',
        pytest_target="tests/test_taxa_coverage_histogram.py",
        pytest_k="test_plate6_reads_compact_remediation_and_nonempty_rank_rows",
    ),
    MutantSpec(
        mutant_id="M02_taxonomy_unresolved_near_duplicate_color",
        rel_path="scripts/generate_tufte_plates.py",
        old_text='"Taxonomy Unresolved": "#9a6e14",',
        new_text='"Taxonomy Unresolved": "#8d5b33",',
        pytest_target="tests/test_taxa_coverage_histogram.py",
        pytest_k="test_plate6_reads_compact_remediation_and_nonempty_rank_rows",
    ),
    MutantSpec(
        mutant_id="M03_rank_color_collision",
        rel_path="scripts/create_taxa_coverage_histogram.py",
        old_text='"genus": "#009e73",',
        new_text='"genus": "#0072b2",',
        pytest_target="tests/test_taxa_coverage_histogram.py",
        pytest_k="test_plate6_reads_compact_remediation_and_nonempty_rank_rows",
    ),
    MutantSpec(
        mutant_id="M04_plate6_subtitle_ignores_rank_table",
        rel_path="scripts/generate_tufte_plates.py",
        old_text=(
            'higher_clause = f"Why Genus/Family signatures achieve {gen_fam_occ} trait coverage"'
        ),
        new_text='higher_clause = "Why Genus/Family signatures achieve 0.0%/0.0% trait coverage"',
        pytest_target="tests/test_deck_builders.py",
        pytest_k="test_committed_phenotype_coverage_tables_drive_deck_and_plate6",
    ),
    MutantSpec(
        mutant_id="M05_methods_plate1_ignores_total_occurrences",
        rel_path="scripts/generate_tufte_plates.py",
        old_text='occurrences_title = f"{total_occurrences:,} All-Profile Taxon Occurrences"',
        new_text='occurrences_title = "0 All-Profile Taxon Occurrences"',
        pytest_target="tests/test_deck_builders.py",
        pytest_k="test_committed_phenotype_coverage_tables_drive_deck_and_plate6",
    ),
    MutantSpec(
        mutant_id="M06_deck_plate6_caveat_ignores_alias_table",
        rel_path="scripts/build_slide_deck.py",
        old_text=(
            'alias_items.append(f"{name} ({taxon_id} → {target}, {occ} occs across {st} studies)")'
        ),
        new_text='alias_items.append("stale_alias_entry")',
        pytest_target="tests/test_deck_builders.py",
        pytest_k="test_committed_phenotype_coverage_tables_drive_deck_and_plate6",
    ),
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _prepare_scratch_tree(scratch_root: Path) -> None:
    for top_entry in (
        ".git",
        "uv.lock",
        "pyproject.toml",
        "src",
        "data",
        "trait-visualization",
        "methods",
    ):
        src_entry = REPO_ROOT / top_entry
        if src_entry.exists():
            (scratch_root / top_entry).symlink_to(src_entry, target_is_directory=src_entry.is_dir())
    for rel_file in WATCHED_FILES:
        dst = scratch_root / rel_file
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(REPO_ROOT / rel_file, dst)


def _run_pytest_in_scratch(scratch: Path, pytest_target: str, pytest_k: str) -> tuple[int, str]:
    old_cwd = Path.cwd()
    out = io.StringIO()
    err = io.StringIO()
    os.environ["GIT_CONFIG_GLOBAL"] = "/dev/null"
    os.environ["GIT_CONFIG_NOSYSTEM"] = "1"
    os.environ["MPLBACKEND"] = "Agg"
    scratch_modules = (
        "create_taxa_coverage_histogram",
        "generate_tufte_plates",
        "build_slide_deck",
        "test_taxa_coverage_histogram",
        "test_deck_builders",
    )
    try:
        os.chdir(scratch)
        for mod_name in scratch_modules:
            sys.modules.pop(mod_name, None)
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            rc = int(
                pytest.main(
                    [
                        "-o",
                        "addopts=",
                        str(scratch / pytest_target),
                        "-k",
                        pytest_k,
                        "-q",
                        "--tb=short",
                    ]
                )
            )
    finally:
        os.chdir(old_cwd)
        for mod_name in scratch_modules:
            sys.modules.pop(mod_name, None)
    return rc, out.getvalue() + "\n" + err.getvalue()


def main() -> int:
    before_hashes = {rel: _sha256(REPO_ROOT / rel) for rel in WATCHED_FILES}

    for spec in MUTANTS:
        with tempfile.TemporaryDirectory(prefix=f"mut_{spec.mutant_id}_") as tmp_dir:
            scratch = Path(tmp_dir)
            _prepare_scratch_tree(scratch)
            target_file = scratch / spec.rel_path
            source_text = target_file.read_text(encoding="utf-8")
            if spec.old_text not in source_text:
                raise RuntimeError(f"{spec.mutant_id}: target snippet not found in {spec.rel_path}")
            target_file.write_text(
                source_text.replace(spec.old_text, spec.new_text, 1),
                encoding="utf-8",
            )
            rc, combined = _run_pytest_in_scratch(scratch, spec.pytest_target, spec.pytest_k)
            if rc == 0 or "AssertionError" not in combined:
                print(
                    f"SURVIVED or wrong failure for {spec.mutant_id} (rc={rc}):\n{combined}",
                    file=sys.stderr,
                )
                return 1
            print(f"KILLED {spec.mutant_id} (rc={rc}, AssertionError confirmed)")

    after_hashes = {rel: _sha256(REPO_ROOT / rel) for rel in WATCHED_FILES}
    if before_hashes != after_hashes:
        raise RuntimeError("Canonical repository files changed during mutant verification!")
    print(f"All {len(MUTANTS)} mutants killed; canonical worktree files verified unchanged.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
