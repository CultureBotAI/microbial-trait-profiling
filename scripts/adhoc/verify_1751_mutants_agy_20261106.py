#!/usr/bin/env python3
"""Verify source-level mutants for #1751 in an isolated scratch copy.

Purpose: Confirm that every #1751 change (setting status='not_tested' when mcnemar_p
is non-finite in compute_correlation, write_outputs, and main while preserving 'tested'
for finite mcnemar_p) is killed by an AssertionError in pytest, with import origin
verified inside pytest and zero writes to the worktree.
Inputs: Worktree src/mtp/guild_maps.py and tests/test_trait_guilds.py.
Outputs: Exit code 0 if all mutants are killed by assertions; non-zero otherwise.
Issue: #1751
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
        "M1_compute_correlation_ignores_nonfinite_mcnemar_p",
        "    elif mcnemar_p is None or not math.isfinite(mcnemar_p):\n"
        '        status = "not_tested"',
        '    elif False:\n        status = "not_tested"',
    ),
    (
        "M2_write_outputs_uses_raw_correlation_status",
        'f"aerobe_vs_anaerobe_paired_signs ({_mcnemar_status(correlation)})"',
        'f"aerobe_vs_anaerobe_paired_signs ({correlation.status})"',
    ),
    (
        "M3_main_cli_uses_raw_correlation_status",
        'f"[{_mcnemar_status(correlation)}]"',
        'f"[{correlation.status}]"',
    ),
    (
        "M4_mcnemar_status_returns_tested_for_nonfinite",
        "def _mcnemar_status(correlation: CorrelationResult) -> str:\n"
        '    """Never report a non-finite McNemar p as \'tested\' (#1751)."""\n'
        '    if correlation.status == "tested" and (\n'
        "        correlation.mcnemar_p is None or not math.isfinite(correlation.mcnemar_p)\n"
        "    ):\n"
        '        return "not_tested"',
        "def _mcnemar_status(correlation: CorrelationResult) -> str:\n"
        '    """Never report a non-finite McNemar p as \'tested\' (#1751)."""\n'
        '    if correlation.status == "tested" and (\n'
        "        correlation.mcnemar_p is None or not math.isfinite(correlation.mcnemar_p)\n"
        "    ):\n"
        '        return "tested"',
    ),
    (
        "M5_mcnemar_status_always_not_tested",
        "    return correlation.status",
        '    return "not_tested"',
    ),
]


def _make_scratch_tree(tmp_root: Path) -> None:
    shutil.copytree(WORKTREE / "src", tmp_root / "src")
    (tmp_root / "methods").mkdir()
    shutil.copy2(
        WORKTREE / "methods" / "guild_condition_body_site_heatmap.md",
        tmp_root / "methods" / "guild_condition_body_site_heatmap.md",
    )
    (tmp_root / "tests").mkdir()
    shutil.copy2(
        WORKTREE / "tests" / "test_trait_guilds.py",
        tmp_root / "tests" / "test_trait_guilds.py",
    )
    conftest = """from pathlib import Path
import mtp.guild_maps as _m

def pytest_sessionstart(session):
    actual = Path(_m.__file__).resolve()
    expected = (Path(__file__).resolve().parents[1] / "src" / "mtp" / "guild_maps.py").resolve()
    if actual != expected:
        raise RuntimeError(f"Import origin mismatch: {actual} != {expected}")
    print(f"PYTEST_ORIGIN_OK:{actual}")
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
                    "--rootdir",
                    str(tmp_root),
                    str(tmp_root / "tests" / "test_trait_guilds.py"),
                    "-k",
                    "test_guild_correlation_table_and_cli_label_nominal_p_values_with_test_counts",
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
            file_path = tmp_root / "src" / "mtp" / "guild_maps.py"
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
            expected_origin = (tmp_root / "src/mtp/guild_maps.py").resolve()
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
    print("ALL 5 #1751 MUTANTS KILLED CLEANLY.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
