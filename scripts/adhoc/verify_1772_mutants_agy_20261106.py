#!/usr/bin/env python3
"""Verify source-level mutants for #1772 in an isolated scratch copy.

Purpose: Confirm that every #1772 change (deriving expected stratum and row labels in
scripts/adhoc/issue727_compare_core_strata.py independently of the producer helpers
_read_tested_summary and _extract_q_or_p in mtp.stratification) is killed by an
AssertionError in pytest, with import origin verified inside pytest and zero writes to
the worktree.
Inputs: Worktree scripts/adhoc/issue727_compare_core_strata.py and
tests/test_stratification_stat_types.py.
Outputs: Exit code 0 if all mutants are killed by assertions; non-zero otherwise.
Issue: #1772
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

MUTANTS: list[tuple[str, str, str]] = [
    (
        "C1_reintroduce_producer_helper_call",
        "        parsed = _split_comments_and_csv(tested_path)",
        (
            "        from mtp.stratification import _read_tested_summary\n"
            "        _read_tested_summary(tested_path)\n"
            "        parsed = _split_comments_and_csv(tested_path)"
        ),
    ),
    (
        "C2_independent_row_stat_type_always_fdr_q",
        '    if p_val is not None:\n        return "p_value"',
        '    if p_val is not None:\n        return "fdr_q"',
    ),
    (
        "C3_nominal_fdr_none_scope_left_global",
        '            stat_type: str | None = "p_value"\n            fdr_scope = "none"',
        (
            '            stat_type: str | None = "p_value"\n'
            '            fdr_scope = raw_scope or "none"'
        ),
    ),
    (
        "C4_untested_trait_sig_basis_not_unavailable",
        'hit_statistic if int(row[col_idx["n_profiles_tested"]]) > 0 else "unavailable"',
        "hit_statistic",
    ),
    (
        "C5_p_threshold_header_ignored",
        'stat_threshold = comments.get("Stat-Threshold") or comments.get("P-Threshold")',
        'stat_threshold = comments.get("Stat-Threshold")',
    ),
]


def _make_scratch_tree(tmp_root: Path) -> None:
    shutil.copytree(WORKTREE / "src", tmp_root / "src")
    shutil.copytree(WORKTREE / "scripts", tmp_root / "scripts")
    (tmp_root / "tests").mkdir()
    shutil.copy2(
        WORKTREE / "tests" / "test_stratification_stat_types.py",
        tmp_root / "tests" / "test_stratification_stat_types.py",
    )
    conftest = """import importlib.util
from pathlib import Path
import mtp.stratification

def pytest_sessionstart(session):
    root = Path(__file__).resolve().parents[1]
    actual_strata = Path(mtp.stratification.__file__).resolve()
    expected_strata = (root / "src" / "mtp" / "stratification.py").resolve()
    if actual_strata != expected_strata:
        raise RuntimeError(f"Import origin mismatch: {actual_strata} != {expected_strata}")
    script_path = (root / "scripts" / "adhoc" / "issue727_compare_core_strata.py").resolve()
    spec = importlib.util.spec_from_file_location("issue727_compare_core_strata", script_path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    actual_script = Path(mod.__file__).resolve()
    if actual_script != script_path:
        raise RuntimeError(f"Script origin mismatch: {actual_script} != {script_path}")
    print(f"PYTEST_ORIGIN_OK:{actual_script}")
"""
    (tmp_root / "tests" / "conftest.py").write_text(conftest, encoding="utf-8")


def _pytest_worker(tmp_root_str: str, out_file_str: str) -> None:
    tmp_root = Path(tmp_root_str)
    os.chdir(tmp_root)
    sys.path.insert(0, str(tmp_root / "src"))
    for mod_name in [
        k for k in list(sys.modules) if k == "mtp" or k.startswith(("mtp.", "bsdbng"))
    ]:
        sys.modules.pop(mod_name, None)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        rc = int(
            pytest.main(
                [
                    "-s",
                    "--tb=long",
                    "-k",
                    "core_strata_comparison",
                    "--rootdir",
                    str(tmp_root),
                    str(tmp_root / "tests" / "test_stratification_stat_types.py"),
                ]
            )
        )
    Path(out_file_str).write_text(f"RC={rc}\n{buf.getvalue()}", encoding="utf-8")


def main() -> int:
    ctx = multiprocessing.get_context("spawn")
    for name, target, replacement in MUTANTS:
        with tempfile.TemporaryDirectory(prefix=f"mutant_{name}_") as tmp_dir:
            tmp_root = Path(tmp_dir)
            _make_scratch_tree(tmp_root)
            file_path = tmp_root / "scripts" / "adhoc" / "issue727_compare_core_strata.py"
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
            expected_origin = (
                tmp_root / "scripts" / "adhoc" / "issue727_compare_core_strata.py"
            ).resolve()
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
    print("ALL 5 #1772 MUTANTS KILLED CLEANLY.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
