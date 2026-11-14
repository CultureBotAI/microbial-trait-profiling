#!/usr/bin/env python3
"""Adversarial mutation verifier for issue #1991 (2026-11-10).

Purpose:
    Verify that tests/test_verify_1838_mutants_1991.py kills at least 6 distinct
    source mutations in isolated scratch copies of
    scripts/adhoc/verify_1838_mutants_agy_20261107.py.
Inputs:
    --scratch-dir: optional root directory for isolated scratch copies.
Outputs:
    Exits 0 when the baseline passes and every mutant is killed by AssertionError.
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
import tempfile
from collections.abc import Sequence
from pathlib import Path
from types import ModuleType

import pytest


def _load_module(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load module spec from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def _run_tests_against_verifier_copy(test_module: ModuleType, verifier_copy: Path) -> None:
    monkeypatch = pytest.MonkeyPatch()
    try:
        monkeypatch.setattr(test_module, "VERIFIER_PATH", verifier_copy)
        test_module.test_main_prints_matched_marker_from_assertion_markers(monkeypatch)
        test_module.test_main_fails_when_no_assertion_marker_matches(monkeypatch)
    finally:
        monkeypatch.undo()


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    parser.add_argument(
        "--scratch-dir",
        type=Path,
        default=None,
        help="Optional directory for isolated scratch copies",
    )
    args = parser.parse_args(argv)

    repo_root = Path(__file__).resolve().parents[2]
    verifier_src = repo_root / "scripts/adhoc/verify_1838_mutants_agy_20261107.py"
    test_src = repo_root / "tests/test_verify_1838_mutants_1991.py"

    test_module = _load_module("test_verify_1838_mutants_1991_mod", test_src)
    verifier_text = verifier_src.read_text(encoding="utf-8")

    # Baseline check
    _run_tests_against_verifier_copy(test_module, verifier_src)

    mutants: list[tuple[str, str]] = [
        (
            "mutant_1_revert_hardcoded_assertion_error",
            verifier_text.replace(
                'print(f"KILLED {mutant.name}: rc=1 {matched_marker} PYTEST_ORIGIN_OK")',
                'print(f"KILLED {mutant.name}: rc=1 AssertionError PYTEST_ORIGIN_OK")',
                1,
            ),
        ),
        (
            "mutant_2_always_pick_first_marker",
            verifier_text.replace(
                "matched_marker = next(m for m in ASSERTION_MARKERS if m in out)",
                "matched_marker = ASSERTION_MARKERS[0]",
                1,
            ),
        ),
        (
            "mutant_3_always_pick_last_marker",
            verifier_text.replace(
                "matched_marker = next(m for m in ASSERTION_MARKERS if m in out)",
                "matched_marker = ASSERTION_MARKERS[-1]",
                1,
            ),
        ),
        (
            "mutant_4_drop_did_not_raise_from_assertion_markers",
            verifier_text.replace(
                '    "Failed: DID NOT RAISE",\n',
                '    "Failed: UNUSED_MARKER",\n',
                1,
            ),
        ),
        (
            "mutant_5_return_zero_when_no_assertion_marker_matches",
            verifier_text.replace(
                'no assertion marker in output\\n{out}")\n                return 1',
                'no assertion marker in output\\n{out}")\n                return 0',
                1,
            ),
        ),
        (
            "mutant_6_omit_pytest_origin_ok_in_killed_line",
            verifier_text.replace(
                'print(f"KILLED {mutant.name}: rc=1 {matched_marker} PYTEST_ORIGIN_OK")',
                'print(f"KILLED {mutant.name}: rc=1 {matched_marker}")',
                1,
            ),
        ),
    ]

    if args.scratch_dir is not None:
        args.scratch_dir.mkdir(parents=True, exist_ok=True)
    scratch_root = str(args.scratch_dir) if args.scratch_dir is not None else None

    killed = 0
    for name, mutated_verifier in mutants:
        if mutated_verifier == verifier_text:
            raise RuntimeError(f"Mutant {name} did not alter verifier source text")
        with tempfile.TemporaryDirectory(prefix=f"verify_1991_{name}_", dir=scratch_root) as raw:
            verifier_copy = Path(raw) / "verify_1838_mutants_agy_20261107.py"
            verifier_copy.write_text(mutated_verifier, encoding="utf-8")
            try:
                _run_tests_against_verifier_copy(test_module, verifier_copy)
            except AssertionError:
                killed += 1
            else:
                raise RuntimeError(f"Mutant survived: {name}")

    if killed != len(mutants) or killed < 6:
        raise RuntimeError(f"Expected at least 6 killed mutants, got {killed}/{len(mutants)}")
    print(f"All {killed}/{len(mutants)} mutants killed by AssertionError.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
