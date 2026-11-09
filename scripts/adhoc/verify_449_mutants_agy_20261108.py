#!/usr/bin/env python3
"""Scratch-copy mutation verifier for #449 (tests/test_tool_pins.py).

Purpose:
  Verify that every change introduced for #449 in `tests/test_tool_pins.py`
  (preserving all hook revs per repository URL rather than collapsing into a
  last-wins dict, requiring a single hook entry per repo URL, and rejecting
  marker-gated or compound `==` pins in `pyproject.toml`) is killed by the
  test suite when mutated in a temporary scratch copy, without mutating the
  live worktree.

Inputs:
  - Worktree files in the repository root (`tests/test_tool_pins.py`,
    `pyproject.toml`, `.pre-commit-config.yaml`, `uv.lock`).

Outputs:
  - Exit status 0 if the unmutated control passes (`PYTEST_ORIGIN_OK`) and all
    mutants (`M1`-`M5`) are killed by `AssertionError`/`Failed`.

Issue / Date:
  - Refs #449 (2026-11-08)
"""

from __future__ import annotations

import contextlib
import importlib.util
import io
import shutil
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]

TRACKED_REL_PATHS = (
    "pyproject.toml",
    ".pre-commit-config.yaml",
    "uv.lock",
    "tests/test_tool_pins.py",
)


@dataclass(frozen=True)
class Mutant:
    name: str
    rel_path: str
    old: str
    new: str
    expected_msg: str


MUTANTS: tuple[Mutant, ...] = (
    Mutant(
        name="M1_collapse_hook_revs_to_last_wins_list",
        rel_path="tests/test_tool_pins.py",
        old='revs.setdefault(str(repo["repo"]), []).append(str(repo["rev"]))',
        new='revs[str(repo["repo"])] = [str(repo["rev"])]',
        expected_msg="DID NOT RAISE AssertionError",
    ),
    Mutant(
        name="M2_collapse_duplicate_same_rev_via_set",
        rel_path="tests/test_tool_pins.py",
        old="assert len(repo_revs) == 1, (",
        new="assert len(set(repo_revs)) == 1, (",
        expected_msg="DID NOT RAISE AssertionError",
    ),
    Mutant(
        name="M3_allow_trailing_marker_after_equals_pin",
        rel_path="tests/test_tool_pins.py",
        old='exact_match = re.fullmatch(r"==\\s*([^\\s;,]+)", rest)',
        new='exact_match = re.match(r"==\\s*([^\\s;,]+)", rest)',
        expected_msg="DID NOT RAISE AssertionError",
    ),
    Mutant(
        name="M4_drop_compound_equals_refusal",
        rel_path="tests/test_tool_pins.py",
        old='elif "==" in rest or (name in PINNED_HOOKS and ";" in rest):',
        new='elif name in PINNED_HOOKS and ";" in rest:',
        expected_msg="marker-gated or compound pin",
    ),
    Mutant(
        name="M5_drop_pinned_hook_marker_refusal",
        rel_path="tests/test_tool_pins.py",
        old='elif "==" in rest or (name in PINNED_HOOKS and ";" in rest):',
        new='elif "==" in rest:',
        expected_msg="DID NOT RAISE AssertionError",
    ),
)


def _populate_scratch(scratch_root: Path) -> None:
    for rel in TRACKED_REL_PATHS:
        src = REPO_ROOT / rel
        dst = scratch_root / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


def _run_scratch_pytest(scratch_root: Path) -> tuple[int, str]:
    test_file = (scratch_root / "tests/test_tool_pins.py").resolve()
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        spec = importlib.util.spec_from_file_location("scratch_test_tool_pins", test_file)
        if spec is None or spec.loader is None:
            raise RuntimeError(f"Could not load spec for {test_file}")
        mod = importlib.util.module_from_spec(spec)
        sys.modules["scratch_test_tool_pins"] = mod
        try:
            spec.loader.exec_module(mod)
            if Path(mod.__file__ or "").resolve() != test_file:
                raise RuntimeError("mod.__file__ mismatch")
            if Path(mod.ROOT).resolve() != scratch_root.resolve():
                raise RuntimeError("mod.ROOT mismatch")
            print("PYTEST_ORIGIN_OK")
            sys.modules.pop("test_tool_pins", None)
            sys.modules.pop("tests.test_tool_pins", None)
            rc = int(pytest.main(["-o", "addopts=", "-q", str(test_file)]))
        finally:
            sys.modules.pop("scratch_test_tool_pins", None)
            sys.modules.pop("test_tool_pins", None)
            sys.modules.pop("tests.test_tool_pins", None)
    return rc, buf.getvalue()


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="verify_449_") as tmp:
        control_root = Path(tmp) / "control"
        _populate_scratch(control_root)
        rc, output = _run_scratch_pytest(control_root)
        if rc != 0 or "PYTEST_ORIGIN_OK" not in output:
            print("CONTROL FAILED:", output, file=sys.stderr)
            return 1
        print("Baseline unmutated scratch run passed (PYTEST_ORIGIN_OK)")

        for mutant in MUTANTS:
            m_root = Path(tmp) / mutant.name
            _populate_scratch(m_root)
            target = m_root / mutant.rel_path
            content = target.read_text(encoding="utf-8")
            if content.count(mutant.old) != 1:
                print(
                    f"Mutant {mutant.name}: expected 1 occurrence of target in {mutant.rel_path}, "
                    f"found {content.count(mutant.old)}",
                    file=sys.stderr,
                )
                return 1
            target.write_text(content.replace(mutant.old, mutant.new, 1), encoding="utf-8")
            m_rc, m_out = _run_scratch_pytest(m_root)
            if m_rc == 0 or "PYTEST_ORIGIN_OK" not in m_out or mutant.expected_msg not in m_out:
                print(
                    f"Mutant {mutant.name} SURVIVED or missed expected failure "
                    f"({mutant.expected_msg!r}):\n{m_out}",
                    file=sys.stderr,
                )
                return 1
            print(f"Mutant {mutant.name}: KILLED (rc={m_rc}, PYTEST_ORIGIN_OK)")

    print(f"All {len(MUTANTS)} mutants killed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
