#!/usr/bin/env python3
"""Scratch-copy mutation verifier for Queue 12 (#1247 + #1248).

Purpose:
  Verify that every target test in tests/test_publication_figure_scan.py passes on a
  scratch copy of the branch (CONTROL) and fails with AssertionError/Failed when each
  of the 6 mutations in src/mtp/publication_figures/scan.py is injected into an
  isolated scratch copy outside the worktree (never mutating tracked files in place).

Inputs:
  - src/mtp/publication_figures/scan.py
  - tests/test_publication_figure_scan.py

Outputs:
  - Exit 0 if CONTROL passes and all 6 mutants die by AssertionError/Failed with
    verified scratch origin; non-zero otherwise.

Issues: #1247, #1248
Date: 2026-11-08
"""

from __future__ import annotations

import contextlib
import io
import multiprocessing
import os
import shutil
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]


@dataclass(frozen=True)
class Mutant:
    """Single source-text mutation against a relative file path."""

    name: str
    rel_path: str
    old: str
    new: str
    nodeid: str


MUTANTS: tuple[Mutant, ...] = (
    Mutant(
        name="M1_omit_trait_id_on_scan_page_1247",
        rel_path="src/mtp/publication_figures/scan.py",
        old='_wrap(bundle.labels[trait_id], 95) + "\\n" + _wrap(trait_id, 99)',
        new="_wrap(bundle.labels[trait_id], 95)",
        nodeid=(
            "tests/test_publication_figure_scan.py::"
            "test_display_label_uses_object_name_with_native_column_label_check"
            "[cellulose catabolic process]"
        ),
    ),
    Mutant(
        name="M2_omit_display_label_on_scan_page_1247",
        rel_path="src/mtp/publication_figures/scan.py",
        old='_wrap(bundle.labels[trait_id], 95) + "\\n" + _wrap(trait_id, 99)',
        new="_wrap(trait_id, 99)",
        nodeid=(
            "tests/test_publication_figure_scan.py::"
            "test_display_label_uses_object_name_with_native_column_label_check"
            "[cellulose catabolic process]"
        ),
    ),
    Mutant(
        name="M3_omit_trait_id_in_page_details_1247",
        rel_path="src/mtp/publication_figures/scan.py",
        old='"trait_id": trait_id,',
        new='"trait_id": "",',
        nodeid=(
            "tests/test_publication_figure_scan.py::"
            "test_display_label_uses_object_name_with_native_column_label_check"
            "[cellulose catabolic process]"
        ),
    ),
    Mutant(
        name="M4_omit_trait_label_in_page_details_1247",
        rel_path="src/mtp/publication_figures/scan.py",
        old='"trait_label": bundle.labels[trait_id],',
        new='"trait_label": trait_id,',
        nodeid=(
            "tests/test_publication_figure_scan.py::"
            "test_display_label_uses_object_name_with_native_column_label_check"
            "[cellulose catabolic process]"
        ),
    ),
    Mutant(
        name="M5_silent_object_name_fallback_when_column_absent_1248",
        rel_path="src/mtp/publication_figures/scan.py",
        old=(
            "    required = (\n"
            '        {"column_id", "object_name", "column_label"} '
            'if display else {"column_id", "column_label"}\n'
            "    )\n"
            '    _require(required.issubset(metadata), "Scan label metadata columns missing")\n'
            '    _require(not metadata["column_id"].duplicated().any(), '
            '"Duplicate scan metadata trait")\n'
            "    result = {\n"
            '        str(row["column_id"]): str(\n'
            '            (row["object_name"] or row["column_label"]) '
            'if display else row["column_label"]\n'
            "        )"
        ),
        new=(
            "    _require(\n"
            '        {"column_id", "column_label"}.issubset(metadata), '
            '"Scan label metadata columns missing"\n'
            "    )\n"
            '    _require(not metadata["column_id"].duplicated().any(), '
            '"Duplicate scan metadata trait")\n'
            "    result = {\n"
            '        str(row["column_id"]): str(\n'
            '            (row.get("object_name") or row["column_label"]) '
            'if display else row["column_label"]\n'
            "        )"
        ),
        nodeid=(
            "tests/test_publication_figure_scan.py::"
            "test_missing_object_name_column_in_metadata_is_refused"
        ),
    ),
    Mutant(
        name="M6_lowercase_rank_1_in_scan_note_1248",
        rel_path="src/mtp/publication_figures/scan.py",
        old='"Rank 1 = smallest nominal p-value. "',
        new='"rank 1 = smallest nominal p-value. "',
        nodeid=(
            "tests/test_publication_figure_scan.py::"
            "test_scan_explains_rank_null_and_readable_recorded_role"
        ),
    ),
)


def _populate_scratch(scratch_root: Path) -> None:
    shutil.copytree(ROOT / "src" / "mtp", scratch_root / "src" / "mtp")
    shutil.copytree(ROOT / "methods", scratch_root / "methods")
    (scratch_root / "tests").mkdir(parents=True, exist_ok=True)
    shutil.copy2(
        ROOT / "tests" / "test_publication_figure_scan.py",
        scratch_root / "tests" / "test_publication_figure_scan.py",
    )


def _worker(scratch_root_str: str, nodeid: str, out_file_str: str) -> None:
    """Run pytest in a clean child process rooted at scratch_root."""
    import importlib.util

    scratch_root = Path(scratch_root_str).resolve()
    os.chdir(scratch_root)
    sys.path.insert(0, str(scratch_root / "src"))
    for name in [
        key
        for key in list(sys.modules)
        if key in {"mtp", "test_publication_figure_scan"} or key.startswith("mtp.")
    ]:
        sys.modules.pop(name, None)
    rel_test_path, _, test_name = nodeid.partition("::")
    test_file = (scratch_root / rel_test_path).resolve()
    from mtp.publication_figures import scan as scan_mod

    expected_scan = (scratch_root / "src/mtp/publication_figures/scan.py").resolve()
    if Path(scan_mod.__file__ or "").resolve() != expected_scan:
        raise RuntimeError(f"Unexpected scan module origin: {scan_mod.__file__}")
    spec = importlib.util.spec_from_file_location("scan_test_mod", str(test_file))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load test spec for {test_file}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    if Path(mod.__file__ or "").resolve() != test_file:
        raise RuntimeError(f"Unexpected test module origin: {mod.__file__}")
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        print("PYTEST_ORIGIN_OK", flush=True)
        rc = int(pytest.main(["-q", "-o", "addopts=", f"{test_file}::{test_name}"]))
    Path(out_file_str).write_text(f"RC={rc}\n{buf.getvalue()}", encoding="utf-8")


def _run_case(scratch_root: Path, nodeid: str) -> tuple[int, str]:
    out_file = scratch_root / "pytest_out.txt"
    ctx = multiprocessing.get_context("spawn")
    proc = ctx.Process(target=_worker, args=(str(scratch_root), nodeid, str(out_file)))
    proc.start()
    proc.join()
    if proc.exitcode != 0 or not out_file.exists():
        raise RuntimeError(f"worker process failed with exitcode={proc.exitcode}")
    raw = out_file.read_text(encoding="utf-8")
    first, _, rest = raw.partition("\n")
    return int(first.removeprefix("RC=")), rest


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="verify_1247_control_") as tmp:
        scratch = Path(tmp)
        _populate_scratch(scratch)
        for nodeid in dict.fromkeys(m.nodeid for m in MUTANTS):
            rc, out = _run_case(scratch, nodeid)
            if rc != 0 or "PYTEST_ORIGIN_OK" not in out:
                print(f"CONTROL FAILED for {nodeid} (rc={rc}):\n{out}", file=sys.stderr)
                return 1
        print("CONTROL OK (all target tests pass on unmutated scratch copy)")

    for mutant in MUTANTS:
        with tempfile.TemporaryDirectory(prefix=f"verify_1247_{mutant.name}_") as tmp:
            scratch = Path(tmp)
            _populate_scratch(scratch)
            target = scratch / mutant.rel_path
            text = target.read_text(encoding="utf-8")
            count = text.count(mutant.old)
            if count != 1:
                print(
                    f"FAIL [{mutant.name}]: expected 1 occurrence in "
                    f"{mutant.rel_path}, found {count}",
                    file=sys.stderr,
                )
                return 1
            target.write_text(text.replace(mutant.old, mutant.new, 1), encoding="utf-8")
            rc, out = _run_case(scratch, mutant.nodeid)
            killed = (
                rc == 1
                and "PYTEST_ORIGIN_OK" in out
                and ("AssertionError" in out or "Failed: DID NOT RAISE" in out)
            )
            if not killed:
                print(
                    f"SURVIVED [{mutant.name}] (rc={rc}):\n{out}",
                    file=sys.stderr,
                )
                return 1
            print(f"KILLED [{mutant.name}] (rc={rc}, AssertionError/Failed, PYTEST_ORIGIN_OK)")

    print("ALL 6 #1247/#1248 MUTANTS KILLED CLEANLY")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
