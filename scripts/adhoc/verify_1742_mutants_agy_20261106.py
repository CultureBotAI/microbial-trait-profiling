#!/usr/bin/env python3
"""Verify source-level mutants for #1742 (second half) in an isolated scratch copy.

Purpose: Confirm that every #1742 change (recording bh_unfiltered_trait_count from
retained_trait_count and stating the unfiltered retained count and select_stratum_traits
selection rule in paired_q_value_family across selection and summary tables) is killed
by an AssertionError in pytest, with import origin verified inside pytest and zero
writes to the worktree.
Inputs: Worktree scripts/create_six_shared_trait_sankeys.py and
tests/test_sankey_nominal_display.py.
Outputs: Exit code 0 if all mutants are killed by assertions; non-zero otherwise.
Issue: #1742
Date: 2026-11-06
"""

from __future__ import annotations

import contextlib
import io
import multiprocessing
import os
import shutil
import sys
import tempfile
from pathlib import Path

import pytest

WORKTREE = Path(__file__).resolve().parents[2]
ASSERTION_MARKERS = ("AssertionError", "Failed: DID NOT RAISE", "Regex pattern did not match")
ERROR_MARKERS = (
    "ImportError",
    "ModuleNotFoundError",
    "SyntaxError",
    "IndentationError",
    "NameError",
    "TypeError",
    "AttributeError",
    "INTERNALERROR",
    "ERROR collecting",
)

MUTANTS: list[tuple[str, str, str, str]] = [
    (
        "S1_format_paired_q_value_family_drops_retained_and_rule",
        "scripts/create_six_shared_trait_sankeys.py",
        (
            "    return (\n"
            '        f"pooled pair-level BH over {size} support/variance-selected traits "\n'
            '        f"({retained_clause}; {SELECTION_RULE_NOTE})"\n'
            "    )"
        ),
        '    return f"pooled pair-level BH over {size} support/variance-selected traits"',
    ),
    (
        "S2_format_paired_q_value_family_ignores_raw_retained",
        "scripts/create_six_shared_trait_sankeys.py",
        "    retained = _extract_bh_family_size(raw_retained)",
        "    retained = None",
    ),
    (
        "S3_build_trait_selection_unfiltered_count_always_na",
        "scripts/create_six_shared_trait_sankeys.py",
        (
            '                "bh_unfiltered_trait_count": pd.NA '
            "if retained_count is None else retained_count,"
        ),
        '                "bh_unfiltered_trait_count": pd.NA,',
    ),
    (
        "S4_build_trait_summary_unfiltered_count_always_na",
        "scripts/create_six_shared_trait_sankeys.py",
        (
            '            "bh_family_size": pd.NA if family_size is None else family_size,\n'
            '            "bh_unfiltered_trait_count": pd.NA '
            "if retained_count is None else retained_count,"
        ),
        (
            '            "bh_family_size": pd.NA if family_size is None else family_size,\n'
            '            "bh_unfiltered_trait_count": pd.NA,'
        ),
    ),
    (
        "S5_build_trait_summary_passes_none_for_retained_in_family_label",
        "scripts/create_six_shared_trait_sankeys.py",
        (
            '            "paired_q_value_family": '
            "_format_paired_q_value_family(family_size, retained_count),"
        ),
        '            "paired_q_value_family": _format_paired_q_value_family(family_size, None),',
    ),
    (
        "S6_methods_attaches_after_filter_to_retained_count_1811",
        "trait-visualization/figure_methods.md",
        "of 344 retained non-MediaDive traits in\n`retained_trait_count` before the",
        "of 344 unfiltered retained traits in\n`retained_trait_count` after the",
    ),
]


def _make_scratch_tree(tmp_root: Path) -> None:
    shutil.copytree(WORKTREE / "src", tmp_root / "src")
    shutil.copytree(WORKTREE / "scripts", tmp_root / "scripts")
    (tmp_root / "trait-visualization").mkdir()
    shutil.copy2(
        WORKTREE / "trait-visualization" / "figure_methods.md",
        tmp_root / "trait-visualization" / "figure_methods.md",
    )
    (tmp_root / "tests").mkdir()
    shutil.copy2(
        WORKTREE / "tests" / "test_sankey_nominal_display.py",
        tmp_root / "tests" / "test_sankey_nominal_display.py",
    )
    conftest = """import sys
from pathlib import Path

def pytest_collection_finish(session):
    mod = sys.modules["create_six_shared_trait_sankeys"]
    actual = Path(mod.__file__).resolve()
    expected = (
        Path(__file__).resolve().parents[1] / "scripts" / "create_six_shared_trait_sankeys.py"
    ).resolve()
    if actual != expected:
        raise RuntimeError(f"Import origin mismatch: {actual} != {expected}")
    print(f"PYTEST_ORIGIN_OK:{actual}")
"""
    (tmp_root / "tests" / "conftest.py").write_text(conftest, encoding="utf-8")


def _pytest_worker(tmp_root_str: str, out_file_str: str) -> None:
    tmp_root = Path(tmp_root_str)
    os.chdir(tmp_root)
    sys.path.insert(0, str(tmp_root / "src"))
    sys.path.insert(0, str(tmp_root / "scripts"))
    for mod_name in [
        k
        for k in list(sys.modules)
        if k in {"mtp", "create_six_shared_trait_sankeys"} or k.startswith(("mtp.", "bsdbng"))
    ]:
        sys.modules.pop(mod_name, None)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        rc = int(
            pytest.main(
                [
                    "-s",
                    "--tb=long",
                    "--rootdir",
                    str(tmp_root),
                    str(tmp_root / "tests" / "test_sankey_nominal_display.py"),
                ]
            )
        )
    Path(out_file_str).write_text(f"RC={rc}\n{buf.getvalue()}", encoding="utf-8")


def main() -> int:
    ctx = multiprocessing.get_context("spawn")
    for name, rel_path, target, replacement in MUTANTS:
        with tempfile.TemporaryDirectory(prefix=f"mutant_{name}_") as tmp_dir:
            tmp_root = Path(tmp_dir)
            _make_scratch_tree(tmp_root)
            file_path = tmp_root / rel_path
            text = file_path.read_text(encoding="utf-8")
            count = text.count(target)
            if count != 1:
                print(f"FAIL [{name}]: expected 1 target occurrence, found {count}")
                return 1
            file_path.write_text(text.replace(target, replacement, 1), encoding="utf-8")

            out_file = tmp_root / "pytest_out.txt"
            proc = ctx.Process(
                target=_pytest_worker,
                args=(str(tmp_root), str(out_file)),
            )
            proc.start()
            proc.join()
            if proc.exitcode != 0 or not out_file.exists():
                print(f"FAIL [{name}]: worker process failed with exitcode={proc.exitcode}")
                return 1
            raw = out_file.read_text(encoding="utf-8")
            first_line, _, out = raw.partition("\n")
            if first_line != "RC=1":
                print(f"FAIL [{name}]: expected RC=1, got {first_line}\n{out}")
                return 1
            expected_origin = (tmp_root / "scripts/create_six_shared_trait_sankeys.py").resolve()
            if f"PYTEST_ORIGIN_OK:{expected_origin}" not in out:
                print(f"FAIL [{name}]: missing PYTEST_ORIGIN_OK in output\n{out}")
                return 1
            if " failed" not in out:
                print(f"FAIL [{name}]: missing 'failed' summary in output\n{out}")
                return 1
            if any(err in out for err in ERROR_MARKERS):
                print(f"FAIL [{name}]: non-assertion error marker found in output\n{out}")
                return 1
            if not any(marker in out for marker in ASSERTION_MARKERS):
                print(f"FAIL [{name}]: no assertion failure marker in output\n{out}")
                return 1
            print(f"KILLED [{name}]: rc=1 assertion failure verified from scratch copy")
    print(f"ALL {len(MUTANTS)} #1742/#1811 MUTANTS KILLED CLEANLY.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
