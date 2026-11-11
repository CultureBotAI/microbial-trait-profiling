#!/usr/bin/env python3
"""Scratch-copy mutation verifier for #2068 and #2069.

Purpose:
  Verify that every target test in tests/test_publication_figure_scan.py passes on an
  unmutated scratch copy of the branch (CONTROL) and fails with the expected assertion
  marker when each of the 8 mutations in src/mtp/publication_figures/scan.py is injected
  into an isolated scratch copy outside the worktree (never mutating tracked files in place).

Inputs:
  - src/mtp/publication_figures/scan.py
  - tests/test_publication_figure_scan.py

Outputs:
  - Exit 0 if CONTROL passes, all selected mutants die by AssertionError or
    Failed: DID NOT RAISE with verified scratch origin, and canonical worktree files
    remain byte-identical; non-zero otherwise.

Issues: #2068, #2069
Date: 2026-11-10
"""

from __future__ import annotations

import argparse
import contextlib
import io
import json
import multiprocessing
import os
import shutil
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SCAN_REL = "src/mtp/publication_figures/scan.py"
TEST_REL = "tests/test_publication_figure_scan.py"
WRAPPED_CASE = (
    "tests/test_publication_figure_scan.py::"
    "test_display_label_uses_object_name_with_native_column_label_check["
    "capable of: an intentionally long synthetic display label that exceeds the "
    "ninety-five character wrap width of the scan heading so it spans two lines [X]]"
)
MISSING_TRUE_CASE = (
    "tests/test_publication_figure_scan.py::"
    "test_missing_object_name_column_in_metadata_is_refused[True]"
)
MISSING_FALSE_CASE = (
    "tests/test_publication_figure_scan.py::"
    "test_missing_object_name_column_in_metadata_is_refused[False]"
)
ASSERTION_MARKERS: tuple[str, ...] = (
    "AssertionError",
    "Failed: DID NOT RAISE",
)


@dataclass(frozen=True)
class Mutant:
    """Single source-text mutation against a relative file path."""

    name: str
    rel_path: str
    old: str
    new: str
    nodeid: str
    expected_marker: str


MUTANTS: tuple[Mutant, ...] = (
    Mutant(
        name="M1_revert_label_block_to_baseline_at_0905_2068",
        rel_path=SCAN_REL,
        old=(
            "        figure.text(\n"
            "            0.12,\n"
            "            0.938,\n"
            '            _wrap(bundle.labels[trait_id], 95) + "\\n" + _wrap(trait_id, 99),\n'
            "            color=INK,\n"
            '            va="top",\n'
            "        )"
        ),
        new=(
            "        figure.text(\n"
            "            0.12,\n"
            "            0.905,\n"
            '            _wrap(bundle.labels[trait_id], 95) + "\\n" + _wrap(trait_id, 99),\n'
            "            color=INK,\n"
            "        )"
        ),
        nodeid=WRAPPED_CASE,
        expected_marker="AssertionError",
    ),
    Mutant(
        name="M2_omit_va_top_at_0938_2068",
        rel_path=SCAN_REL,
        old='            va="top",',
        new='            va="baseline",',
        nodeid=WRAPPED_CASE,
        expected_marker="AssertionError",
    ),
    Mutant(
        name="M3_label_block_collides_with_title_2068",
        rel_path=SCAN_REL,
        old="            0.938,",
        new="            0.975,",
        nodeid=WRAPPED_CASE,
        expected_marker="AssertionError",
    ),
    Mutant(
        name="M4_label_block_collides_with_counts_line_2068",
        rel_path=SCAN_REL,
        old="            0.938,",
        new="            0.905,",
        nodeid=WRAPPED_CASE,
        expected_marker="AssertionError",
    ),
    Mutant(
        name="M5_require_object_name_only_when_display_true_2069",
        rel_path=SCAN_REL,
        old='    missing = sorted({"column_id", "object_name", "column_label"} - set(metadata))',
        new=(
            "    missing = sorted(\n"
            '        ({"column_id", "object_name", "column_label"} '
            'if display else {"column_id", "column_label"})\n'
            "        - set(metadata)\n"
            "    )"
        ),
        nodeid=MISSING_FALSE_CASE,
        expected_marker="Failed: DID NOT RAISE",
    ),
    Mutant(
        name="M6_omit_missing_column_names_from_error_message_2069",
        rel_path=SCAN_REL,
        old=(
            '    _require(not missing, f"Scan label metadata columns missing: '
            "{', '.join(missing)}\")"
        ),
        new='    _require(not missing, "Scan label metadata columns missing")',
        nodeid=MISSING_TRUE_CASE,
        expected_marker="AssertionError",
    ),
    Mutant(
        name="M7_report_only_first_missing_column_2069",
        rel_path=SCAN_REL,
        old="f\"Scan label metadata columns missing: {', '.join(missing)}\"",
        new='f"Scan label metadata columns missing: {missing[0]}"',
        nodeid=MISSING_TRUE_CASE,
        expected_marker="AssertionError",
    ),
    Mutant(
        name="M8_omit_column_id_from_required_columns_2069",
        rel_path=SCAN_REL,
        old='    missing = sorted({"column_id", "object_name", "column_label"} - set(metadata))',
        new='    missing = sorted({"object_name", "column_label"} - set(metadata))',
        nodeid=MISSING_TRUE_CASE,
        expected_marker="AssertionError",
    ),
)


def _populate_scratch(scratch_root: Path) -> None:
    shutil.copytree(ROOT / "src" / "mtp", scratch_root / "src" / "mtp")
    shutil.copytree(ROOT / "methods", scratch_root / "methods")
    (scratch_root / "tests").mkdir(parents=True, exist_ok=True)
    shutil.copy2(
        ROOT / TEST_REL,
        scratch_root / TEST_REL,
    )


def _worker(scratch_root_str: str, nodeid: str, out_file_str: str) -> None:
    """Run pytest in a clean child process rooted at scratch_root."""
    import importlib.util

    scratch_root = Path(scratch_root_str).resolve()
    os.chdir(scratch_root)
    os.environ["GIT_CONFIG_GLOBAL"] = os.devnull
    os.environ["GIT_CONFIG_NOSYSTEM"] = "1"
    mpl_dir = scratch_root.parent / ".mplconfig"
    mpl_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_dir))
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

    expected_scan = (scratch_root / SCAN_REL).resolve()
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
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--canary",
        action="store_true",
        help="Run control and the first mutant only",
    )
    parser.add_argument(
        "--scratch-dir",
        type=Path,
        default=None,
        help="Parent directory outside the repository for temporary scratch copies",
    )
    args = parser.parse_args()
    scratch_parent = args.scratch_dir.resolve() if args.scratch_dir is not None else None
    if scratch_parent is not None:
        if scratch_parent.is_relative_to(ROOT):
            raise RuntimeError("--scratch-dir must be outside the repository root")
        scratch_parent.mkdir(parents=True, exist_ok=True)

    originals = {rel: (ROOT / rel).read_bytes() for rel in (SCAN_REL, TEST_REL)}
    selected = MUTANTS[:1] if args.canary else MUTANTS
    run_dir = Path(
        tempfile.mkdtemp(
            prefix="verify_2068_2069_canary_" if args.canary else "verify_2068_2069_all_",
            dir=scratch_parent,
        )
    )
    verdicts: list[dict[str, object]] = []

    try:
        control_scratch = run_dir / "control"
        control_scratch.mkdir(parents=True, exist_ok=True)
        _populate_scratch(control_scratch)
        for nodeid in dict.fromkeys(m.nodeid for m in selected):
            rc, out = _run_case(control_scratch, nodeid)
            if rc != 0 or "PYTEST_ORIGIN_OK" not in out:
                print(f"CONTROL FAILED for {nodeid} (rc={rc}):\n{out}", file=sys.stderr)
                return 1
        verdicts.append({"case": "control", "rc": 0, "verdict": "CONTROL OK"})
        print("CONTROL OK (all target tests pass on unmutated scratch copy)")

        for mutant in selected:
            scratch = run_dir / mutant.name
            scratch.mkdir(parents=True, exist_ok=True)
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
            if rc != 1 or "PYTEST_ORIGIN_OK" not in out:
                print(
                    f"SURVIVED [{mutant.name}] (rc={rc}):\n{out}",
                    file=sys.stderr,
                )
                return 1
            if not any(marker in out for marker in ASSERTION_MARKERS):
                print(
                    f"FAIL [{mutant.name}]: no assertion marker in output\n{out}",
                    file=sys.stderr,
                )
                return 1
            matched_marker = next(m for m in ASSERTION_MARKERS if m in out)
            if matched_marker != mutant.expected_marker:
                print(
                    f"FAIL [{mutant.name}]: expected {mutant.expected_marker}, "
                    f"got {matched_marker}\n{out}",
                    file=sys.stderr,
                )
                return 1
            verdicts.append(
                {
                    "case": mutant.name,
                    "nodeid": mutant.nodeid,
                    "rc": rc,
                    "matched_marker": matched_marker,
                    "log_file": str(scratch / "pytest_out.txt"),
                }
            )
            print(f"KILLED {mutant.name}: rc=1 {matched_marker} PYTEST_ORIGIN_OK")
        (run_dir / "verdicts.json").write_text(
            json.dumps(verdicts, indent=2) + "\n", encoding="utf-8"
        )
    finally:
        if scratch_parent is None and run_dir.exists():
            shutil.rmtree(run_dir, ignore_errors=True)
        for rel, content in originals.items():
            if (ROOT / rel).read_bytes() != content:
                raise RuntimeError(f"Canonical worktree file changed during verification: {rel}")

    print(
        f"ALL {len(selected)}/{len(selected)} #2068/#2069 MUTANTS KILLED CLEANLY; "
        f"worktree targets byte-identical; evidence: {run_dir}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
