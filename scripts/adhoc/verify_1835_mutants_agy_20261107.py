#!/usr/bin/env python3
"""Verify source-level mutants for #1835 in an isolated scratch copy.

Purpose:
  Confirm that the #1835 two-way slashless relation test and character-class candidate
  generation in tests/test_local_path_patterns.py kill one-sided widenings and narrowings
  of both src/mtp/local_path_patterns.py and src/mtp/provenance_paths.py (including
  Volumes, people, var/tmp, srv, and '+' in the user-name class) by AssertionError (rc=1)
  in pytest, with import origins verified inside pytest and zero writes to the worktree.

Inputs:
  - Worktree src/mtp/local_path_patterns.py
  - Worktree src/mtp/provenance_paths.py
  - Worktree tests/test_local_path_patterns.py
  - Worktree tests/test_no_absolute_local_paths.py

Outputs:
  - Exit code 0 if the unmutated control passes (rc=0) and all 12 mutants are killed by
    assertion failure (rc=1) with verified scratch import origins; non-zero otherwise.

Issue: #1835 (from PR #1822 delta review)
Date: 2026-11-07
"""

from __future__ import annotations

import argparse
import contextlib
import io
import multiprocessing
import os
import shutil
import sys
import tempfile
from pathlib import Path

import pytest

DEFAULT_WORKTREE = Path(__file__).resolve().parents[2]
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

LPP = "src/mtp/local_path_patterns.py"
PROV = "src/mtp/provenance_paths.py"
SLASHLESS_HOME = "(?:(?:Users|home)/[A-Za-z0-9._-]+"

MUTANTS: list[tuple[str, str, str, str]] = [
    (
        "M1_slashless_home_adds_volumes",
        PROV,
        SLASHLESS_HOME,
        "(?:(?:Users|home|Volumes)/[A-Za-z0-9._-]+",
    ),
    (
        "M2_slashless_home_adds_people",
        PROV,
        SLASHLESS_HOME,
        "(?:(?:Users|home|people)/[A-Za-z0-9._-]+",
    ),
    (
        "M3_slashless_temp_adds_var_tmp",
        PROV,
        '|private/tmp(?=/))"',
        '|private/tmp(?=/)|var/tmp(?=/))"',
    ),
    (
        "M4_shared_home_adds_srv",
        LPP,
        'r"/(?:Users|home)/',
        'r"/(?:Users|home|srv)/',
    ),
    (
        "M5_shared_name_adds_plus",
        LPP,
        "[A-Za-z0-9._-]+|",
        "[A-Za-z0-9._+-]+|",
    ),
    (
        "M6_slashless_home_adds_srv",
        PROV,
        SLASHLESS_HOME,
        "(?:(?:Users|home|srv)/[A-Za-z0-9._-]+",
    ),
    (
        "M7_slashless_name_adds_plus",
        PROV,
        SLASHLESS_HOME,
        "(?:(?:Users|home)/[A-Za-z0-9._+-]+",
    ),
    (
        "M8_shared_home_adds_volumes",
        LPP,
        'r"/(?:Users|home)/',
        'r"/(?:Users|home|Volumes)/',
    ),
    (
        "M9_shared_home_adds_people",
        LPP,
        'r"/(?:Users|home)/',
        'r"/(?:Users|home|people)/',
    ),
    (
        "M10_shared_temp_adds_var_tmp",
        LPP,
        '|/private/(?:tmp)(?=/)"',
        '|/private/(?:tmp)(?=/)|/var/(?:tmp)(?=/)"',
    ),
    (
        "M11_slashless_name_drops_dot",
        PROV,
        SLASHLESS_HOME,
        "(?:(?:Users|home)/[A-Za-z0-9_-]+",
    ),
    (
        "M12_shared_name_drops_dot",
        LPP,
        "[A-Za-z0-9._-]+|",
        "[A-Za-z0-9_-]+|",
    ),
]


def _make_scratch_tree(worktree: Path, tmp_root: Path) -> None:
    shutil.copytree(worktree / "src", tmp_root / "src")
    shutil.copy2(worktree / "pyproject.toml", tmp_root / "pyproject.toml")
    (tmp_root / "tests").mkdir()
    for test_name in ("test_local_path_patterns.py", "test_no_absolute_local_paths.py"):
        shutil.copy2(worktree / "tests" / test_name, tmp_root / "tests" / test_name)
    conftest = """from pathlib import Path
import mtp.local_path_patterns as _lpp
import mtp.provenance_paths as _prov
import test_no_absolute_local_paths as _guard

def pytest_sessionstart(session):
    root = Path(__file__).resolve().parents[1]
    expected_lpp = (root / "src" / "mtp" / "local_path_patterns.py").resolve()
    expected_prov = (root / "src" / "mtp" / "provenance_paths.py").resolve()
    expected_guard = (root / "tests" / "test_no_absolute_local_paths.py").resolve()
    assert Path(_lpp.__file__).resolve() == expected_lpp
    assert Path(_prov.__file__).resolve() == expected_prov
    assert Path(_guard.__file__).resolve() == expected_guard
    print(f"PYTEST_ORIGIN_OK:{expected_lpp}:{expected_prov}")
"""
    (tmp_root / "tests" / "conftest.py").write_text(conftest, encoding="utf-8")


def _pytest_worker(tmp_root_str: str, out_file_str: str) -> None:
    tmp_root = Path(tmp_root_str)
    os.chdir(tmp_root)
    sys.path.insert(0, str(tmp_root / "tests"))
    sys.path.insert(0, str(tmp_root / "src"))
    for mod_name in [
        k
        for k in list(sys.modules)
        if k in {"mtp", "test_no_absolute_local_paths"} or k.startswith(("mtp.", "bsdbng"))
    ]:
        sys.modules.pop(mod_name, None)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        rc = int(
            pytest.main(
                [
                    "-s",
                    "-o",
                    "addopts=",
                    "--tb=long",
                    "--rootdir",
                    str(tmp_root),
                    str(tmp_root / "tests" / "test_local_path_patterns.py"),
                ]
            )
        )
    Path(out_file_str).write_text(f"RC={rc}\n{buf.getvalue()}", encoding="utf-8")


def _run_one(ctx: multiprocessing.context.SpawnContext, tmp_root: Path) -> tuple[str, str]:
    out_file = tmp_root / "pytest_out.txt"
    proc = ctx.Process(target=_pytest_worker, args=(str(tmp_root), str(out_file)))
    proc.start()
    proc.join()
    if proc.exitcode != 0 or not out_file.exists():
        raise RuntimeError(f"worker process failed with exitcode={proc.exitcode}")
    raw = out_file.read_text(encoding="utf-8")
    first_line, _, out = raw.partition("\n")
    expected_lpp = (tmp_root / LPP).resolve()
    expected_prov = (tmp_root / PROV).resolve()
    if f"PYTEST_ORIGIN_OK:{expected_lpp}:{expected_prov}" not in out:
        raise RuntimeError(f"missing PYTEST_ORIGIN_OK in output:\n{out}")
    return first_line, out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worktree", type=Path, default=DEFAULT_WORKTREE)
    args = parser.parse_args()
    worktree = args.worktree.resolve()
    ctx = multiprocessing.get_context("spawn")

    with tempfile.TemporaryDirectory(prefix="mutant_control_") as tmp_dir:
        tmp_root = Path(tmp_dir)
        _make_scratch_tree(worktree, tmp_root)
        first_line, out = _run_one(ctx, tmp_root)
        if first_line != "RC=0":
            print(f"FAIL [CONTROL]: expected RC=0, got {first_line}\n{out}")
            return 1
        print("PASS [CONTROL]: rc=0 verified from scratch copy")

    failures = 0
    for name, rel_path, target, replacement in MUTANTS:
        with tempfile.TemporaryDirectory(prefix=f"mutant_{name}_") as tmp_dir:
            tmp_root = Path(tmp_dir)
            _make_scratch_tree(worktree, tmp_root)
            file_path = tmp_root / rel_path
            text = file_path.read_text(encoding="utf-8")
            count = text.count(target)
            if count != 1:
                print(f"FAIL [{name}]: expected 1 target occurrence, found {count}")
                return 1
            file_path.write_text(text.replace(target, replacement, 1), encoding="utf-8")
            first_line, out = _run_one(ctx, tmp_root)
            if first_line != "RC=1":
                print(f"SURVIVED [{name}]: expected RC=1, got {first_line}")
                failures += 1
                continue
            if any(err in out for err in ERROR_MARKERS):
                print(f"FAIL [{name}]: non-assertion error marker found in output\n{out}")
                return 1
            if not any(marker in out for marker in ASSERTION_MARKERS):
                print(f"FAIL [{name}]: no assertion failure marker in output\n{out}")
                return 1
            print(f"KILLED [{name}]: rc=1 assertion failure verified from scratch copy")

    if failures:
        print(f"FAILED: {failures}/{len(MUTANTS)} mutants survived.")
        return 1
    print(f"ALL {len(MUTANTS)} #1835 MUTANTS KILLED CLEANLY.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
