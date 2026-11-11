#!/usr/bin/env python3
"""Scratch-copy mutation verifier for #1997, #2036, #2037, and #2038.

Purpose:
  Verify that tests/test_cellulose_tufte_readme_scope_2036.py passes on an
  unmutated scratch copy (CONTROL) and fails with AssertionError when each of
  the 8 scope/disposition mutations in
  tests/fixtures/cellulose_tufte_historical_sources/README.md is injected into
  an isolated scratch copy outside the worktree.

Inputs:
  - tests/fixtures/cellulose_tufte_historical_sources/README.md
  - tests/test_cellulose_tufte_readme_scope_2036.py

Outputs:
  - verdicts.json in the scratch run directory (and parent scratch directory if
    supplied) with non-empty detail strings for CONTROL and every mutant.
  - Exit 0 if CONTROL passes, all selected mutants die by AssertionError with
    verified scratch origin, and canonical worktree files remain byte-identical.

Issues: #1997, #2036, #2037, #2038
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
README_REL = "tests/fixtures/cellulose_tufte_historical_sources/README.md"
TEST_REL = "tests/test_cellulose_tufte_readme_scope_2036.py"
TARGET_NODEID = (
    "tests/test_cellulose_tufte_readme_scope_2036.py::"
    "test_historical_sources_readme_scope_and_unpinned_disposition"
)
ASSERTION_MARKERS: tuple[str, ...] = ("AssertionError",)


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
        name="M1_restore_repeated_scope_sentence_2038",
        rel_path=README_REL,
        old=(
            "The table below lists what `test_cellulose_tufte_statistics_sources.py` binds to\n"
            "identified historical sources; the other test files are summarised after it.\n"
            "They do not pin every numeric token or every repeated occurrence,"
        ),
        new=(
            "The table below lists what `test_cellulose_tufte_statistics_sources.py` binds to\n"
            "identified historical sources; the other test files are summarised after it. These\n"
            "tests bind selected displayed values and labels to identified historical sources. "
            "They do not pin every numeric token or every repeated occurrence,"
        ),
        nodeid=TARGET_NODEID,
        expected_marker="AssertionError",
    ),
    Mutant(
        name="M2_revert_overlap_summary_to_stale_2036",
        rel_path=README_REL,
        old=(
            "`test_cellulose_tufte_overlap_sources.py` pins selected set sizes, direction "
            "labels, percentages, distances and Plate III's nominal p-values for Plates I to IV "
            "and figures 7 to 9; it is summarised here, not tabulated."
        ),
        new=(
            "`test_cellulose_tufte_overlap_sources.py` pins selected set sizes and direction "
            "labels\nfor Plates I to IV and figures 7 to 9; it is summarised here, not tabulated."
        ),
        nodeid=TARGET_NODEID,
        expected_marker="AssertionError",
    ),
    Mutant(
        name="M3_drop_plate_three_nominal_p_from_overlap_summary_2036",
        rel_path=README_REL,
        old="percentages, distances and Plate III's nominal p-values",
        new="percentages and distances",
        nodeid=TARGET_NODEID,
        expected_marker="AssertionError",
    ),
    Mutant(
        name="M4_revert_mutation_summary_to_pins_above_2037",
        rel_path=README_REL,
        old=(
            "`test_cellulose_tufte_source_pin_mutations.py` checks that every tabulated pin, "
            "and one overlap pin, fails on targeted mutations."
        ),
        new=(
            "`test_cellulose_tufte_source_pin_mutations.py` checks that the pins\n"
            "above fail on targeted mutations."
        ),
        nodeid=TARGET_NODEID,
        expected_marker="AssertionError",
    ),
    Mutant(
        name="M5_revert_unpinned_tail_to_broader_source_pins_1997",
        rel_path=README_REL,
        old=(
            "The ten sourced-but-unpinned tokens from #1596/#1997 ('10 PD\u2229CRC', '= 40 taxa', "
            "'46 PD', and '27 shared (57.4%)' in fig7; '17 &lt; 27 (PD\u2013CRC)' in fig9; "
            "'17 of 32' in Plate I; '45%\u201357%' in Plate III; '50.4% shared' in Plate IV; "
            "'n = 2,335 pairs' in Plate V; 'Both FA (58.2% up' in Plate VI) have explicit "
            "disposition: won't pin (#1997), because the plates were retired from default "
            "generation by #1693 (#1572)."
        ),
        new="Broader source pins remain a separate task.",
        nodeid=TARGET_NODEID,
        expected_marker="AssertionError",
    ),
    Mutant(
        name="M6_omit_fig7_token_from_unpinned_list_1997",
        rel_path=README_REL,
        old="'27 shared (57.4%)'",
        new="'omitted fig7 token'",
        nodeid=TARGET_NODEID,
        expected_marker="AssertionError",
    ),
    Mutant(
        name="M7_omit_fig9_and_plate_tokens_from_unpinned_list_1997",
        rel_path=README_REL,
        old=(
            "'17 &lt; 27 (PD\u2013CRC)' in fig9; '17 of 32' in Plate I; "
            "'45%\u201357%' in Plate III; "
            "'50.4% shared' in Plate IV; 'n = 2,335 pairs' in Plate V; "
            "'Both FA (58.2% up' in Plate VI"
        ),
        new="'omitted remaining tokens'",
        nodeid=TARGET_NODEID,
        expected_marker="AssertionError",
    ),
    Mutant(
        name="M8_omit_wont_pin_retirement_disposition_1997",
        rel_path=README_REL,
        old=(
            "won't pin (#1997), because the plates were retired from default generation "
            "by #1693 (#1572)."
        ),
        new="deferred to a follow-up issue.",
        nodeid=TARGET_NODEID,
        expected_marker="AssertionError",
    ),
)


def _populate_scratch(scratch_root: Path) -> None:
    readme_dst = scratch_root / README_REL
    readme_dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(ROOT / README_REL, readme_dst)
    test_dst = scratch_root / TEST_REL
    test_dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(ROOT / TEST_REL, test_dst)


def _worker(scratch_root_str: str, nodeid: str, out_file_str: str) -> None:
    """Run pytest in a clean child process rooted at scratch_root."""
    import importlib.util

    scratch_root = Path(scratch_root_str).resolve()
    os.chdir(scratch_root)
    os.environ["GIT_CONFIG_GLOBAL"] = os.devnull
    os.environ["GIT_CONFIG_NOSYSTEM"] = "1"
    for name in [
        key
        for key in list(sys.modules)
        if key in {"test_cellulose_tufte_readme_scope_2036", "readme_scope_test_mod"}
    ]:
        sys.modules.pop(name, None)
    rel_test_path, _, test_name = nodeid.partition("::")
    test_file = (scratch_root / rel_test_path).resolve()
    expected_readme = (scratch_root / README_REL).resolve()
    spec = importlib.util.spec_from_file_location("readme_scope_test_mod", str(test_file))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load test spec for {test_file}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    if Path(mod.__file__ or "").resolve() != test_file:
        raise RuntimeError(f"Unexpected test module origin: {mod.__file__}")
    resolved_readme = Path(getattr(mod, "README_PATH", "")).resolve()
    if resolved_readme != expected_readme:
        raise RuntimeError(f"Unexpected README_PATH origin: {resolved_readme}")
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


def _extract_failure_detail(output: str, matched_marker: str) -> str:
    for line in output.splitlines():
        stripped = line.strip()
        if matched_marker in stripped:
            return f"KILLED by {matched_marker} with PYTEST_ORIGIN_OK: {stripped}"
    return f"KILLED by {matched_marker} with PYTEST_ORIGIN_OK"


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

    originals = {rel: (ROOT / rel).read_bytes() for rel in (README_REL, TEST_REL)}
    selected = MUTANTS[:1] if args.canary else MUTANTS
    run_dir = Path(
        tempfile.mkdtemp(
            prefix="verify_2036_canary_" if args.canary else "verify_2036_all_",
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
        verdicts.append(
            {
                "case": "control",
                "name": "control",
                "nodeid": TARGET_NODEID,
                "rc": 0,
                "verdict": "CONTROL OK",
                "matched_marker": "PYTEST_ORIGIN_OK",
                "detail": (
                    "CONTROL OK: target test passed on unmutated scratch copy "
                    "with verified PYTEST_ORIGIN_OK"
                ),
                "log_file": str(control_scratch / "pytest_out.txt"),
            }
        )
        print("CONTROL OK (target test passed on unmutated scratch copy)")

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
            detail = _extract_failure_detail(out, matched_marker)
            if not detail.strip():
                raise RuntimeError(f"Empty detail for mutant {mutant.name}")
            verdicts.append(
                {
                    "case": mutant.name,
                    "name": mutant.name,
                    "nodeid": mutant.nodeid,
                    "rc": rc,
                    "verdict": "KILLED",
                    "matched_marker": matched_marker,
                    "detail": detail,
                    "log_file": str(scratch / "pytest_out.txt"),
                }
            )
            print(f"KILLED {mutant.name}: rc=1 {matched_marker} PYTEST_ORIGIN_OK")
        verdicts_payload = json.dumps(verdicts, indent=2) + "\n"
        (run_dir / "verdicts.json").write_text(verdicts_payload, encoding="utf-8")
        if scratch_parent is not None:
            (scratch_parent / "verdicts.json").write_text(verdicts_payload, encoding="utf-8")
    finally:
        if scratch_parent is None and run_dir.exists():
            shutil.rmtree(run_dir, ignore_errors=True)
        for rel, content in originals.items():
            if (ROOT / rel).read_bytes() != content:
                raise RuntimeError(f"Canonical worktree file changed during verification: {rel}")

    print(
        f"ALL {len(selected)}/{len(selected)} #1997/#2036/#2037/#2038 MUTANTS KILLED CLEANLY; "
        f"worktree targets byte-identical; evidence: {run_dir}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
