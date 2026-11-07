#!/usr/bin/env python3
"""Verify source-level mutants for #1714 in an isolated scratch copy.

Purpose: Confirm that every #1714 change (returning NaN for undefined paired regression
slopes when N < 2 or delta_b is constant, and formatting non-finite slopes as NA in both
report writers while preserving 0.0000 for genuine zero slopes) is killed by an
AssertionError in pytest, with import origin verified inside pytest and zero writes to
the worktree.
Inputs: Worktree src/mtp/density_adjustment.py, scripts/audit_trait_density_correlation.py,
and tests/test_density_*.py.
Outputs: Exit code 0 if all mutants are killed by assertions; non-zero otherwise.
Issue: #1714
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

MUTANTS: list[tuple[str, str, str, str, list[str]]] = [
    (
        "Z1_revert_slope_paired_nan_to_zero",
        "src/mtp/density_adjustment.py",
        '    else:\n        slope_paired = float("nan")\n        pval_paired_slope = float("nan")',
        '    else:\n        slope_paired = 0.0\n        pval_paired_slope = float("nan")',
        [
            "tests/test_density_adjustment.py",
            "-k",
            "test_evaluate_conditional_twin_test_constant_and_unpaired",
        ],
    ),
    (
        "Z2_format_paired_regression_slope_no_na_branch",
        "src/mtp/density_adjustment.py",
        '    return f"{value:.4f}" if math.isfinite(value) else "NA"',
        '    return f"{value:.4f}"',
        [
            "tests/test_density_nominal_reporting.py",
            "-k",
            "test_format_paired_regression_slope",
        ],
    ),
    (
        "Z3_format_paired_regression_slope_always_na",
        "src/mtp/density_adjustment.py",
        '    return f"{value:.4f}" if math.isfinite(value) else "NA"',
        '    return "NA"',
        [
            "tests/test_density_nominal_reporting.py",
            "-k",
            "test_unmocked_undefined_vs_zero_paired_slope_report_and_readers",
        ],
    ),
    (
        "Z4_module_report_raw_float_format",
        "src/mtp/density_adjustment.py",
        '            f"Paired regression slope: "\n'
        '            f"{format_paired_regression_slope(twin_summary.paired_slope_a_on_b)} "',
        '            f"Paired regression slope: {twin_summary.paired_slope_a_on_b:.4f} "',
        [
            "tests/test_density_nominal_reporting.py",
            "-k",
            "test_unmocked_undefined_vs_zero_paired_slope_report_and_readers and module",
        ],
    ),
    (
        "Z5_consolidated_report_raw_float_format",
        "scripts/audit_trait_density_correlation.py",
        '            f"Paired regression slope: "\n'
        '            f"{format_paired_regression_slope(twin.paired_slope_a_on_b)} "',
        '            f"Paired regression slope: {twin.paired_slope_a_on_b:.4f} "',
        [
            "tests/test_density_nominal_reporting.py",
            "-k",
            "test_unmocked_undefined_vs_zero_paired_slope_report_and_readers and consolidated",
        ],
    ),
]


def _make_scratch_tree(tmp_root: Path) -> None:
    shutil.copytree(WORKTREE / "src", tmp_root / "src")
    shutil.copytree(WORKTREE / "scripts", tmp_root / "scripts")
    (tmp_root / "tests").mkdir()
    for test_name in ("test_density_adjustment.py", "test_density_nominal_reporting.py"):
        shutil.copy2(WORKTREE / "tests" / test_name, tmp_root / "tests" / test_name)
    conftest = """from pathlib import Path
import mtp.density_adjustment as _m

def pytest_sessionstart(session):
    actual = Path(_m.__file__).resolve()
    expected = (
        Path(__file__).resolve().parents[1] / "src" / "mtp" / "density_adjustment.py"
    ).resolve()
    if actual != expected:
        raise RuntimeError(f"Import origin mismatch: {actual} != {expected}")
    print(f"PYTEST_ORIGIN_OK:{actual}")
"""
    (tmp_root / "tests" / "conftest.py").write_text(conftest, encoding="utf-8")


def _pytest_worker(
    tmp_root_str: str,
    pytest_args: list[str],
    out_file_str: str,
) -> None:
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
                    "--tb=short",
                    "--rootdir",
                    str(tmp_root),
                    *[
                        str(tmp_root / arg) if arg.startswith("tests/") else arg
                        for arg in pytest_args
                    ],
                ]
            )
        )
    Path(out_file_str).write_text(f"RC={rc}\n{buf.getvalue()}", encoding="utf-8")


def main() -> int:
    ctx = multiprocessing.get_context("spawn")
    for name, rel_path, target, replacement, pytest_args in MUTANTS:
        with tempfile.TemporaryDirectory(prefix=f"mutant_{name}_") as tmp_dir:
            tmp_root = Path(tmp_dir)
            _make_scratch_tree(tmp_root)
            file_path = tmp_root / rel_path
            text = file_path.read_text(encoding="utf-8")
            count = text.count(target)
            if count != 1:
                print(f"FAIL [{name}]: expected 1 target occurrence in {rel_path}, found {count}")
                return 1
            file_path.write_text(text.replace(target, replacement, 1), encoding="utf-8")

            out_file = tmp_root / "pytest_out.txt"
            proc = ctx.Process(
                target=_pytest_worker,
                args=(str(tmp_root), pytest_args, str(out_file)),
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
            expected_origin = (tmp_root / "src/mtp/density_adjustment.py").resolve()
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
    print("ALL 5 #1714 MUTANTS KILLED CLEANLY.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
