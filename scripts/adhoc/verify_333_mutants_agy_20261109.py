#!/usr/bin/env python3
"""Isolated scratch-copy mutant verifier for issue #333.

Purpose:
    Verify that tests/test_table1_maps_methods_doc_333.py kills mutations to the
    matrix-specific host-filter drop explanation in
    methods/table1_trait_condition_body_site_maps.md with AssertionError in an
    isolated temporary directory, without modifying the live worktree.
Inputs:
    methods/table1_trait_condition_body_site_maps.md and
    tests/test_table1_maps_methods_doc_333.py in the repository worktree.
Outputs:
    Stdout control and mutant kill report (exit 0 iff control passes and all
    mutants die by AssertionError).
Issue/PR:
    Issue #333
Date:
    2026-11-09
"""

from __future__ import annotations

import contextlib
import hashlib
import io
import shutil
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
DOC_REL = Path("methods/table1_trait_condition_body_site_maps.md")
TEST_REL = Path("tests/test_table1_maps_methods_doc_333.py")


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


@dataclass(frozen=True)
class Mutant:
    name: str
    old: str
    new: str


MUTANTS: tuple[Mutant, ...] = (
    Mutant(
        name="M01_drop_canonical_0_row_sentence",
        old=(
            "The host filter drops 0 rows on the canonical `count_input_16s_wms_v1` profile\n"
            "matrix identified below (`519a095b…` / `43b2e230…`, 5,629 profile rows), which\n"
            "is built upstream from human-host 16S and shotgun profiles only (#656, #742)."
        ),
        new="The host filter runs upstream (#656, #742).",
    ),
    Mutant(
        name="M02_drop_1500_historical_matrix_attribution",
        old=(
            "`6e70724e…` over 7,329\n"
            "profile rows, where 1,500 rows were dropped as non-human (#333)"
        ),
        new="`6e70724e…` over 7,329\nprofile rows (#333)",
    ),
    Mutant(
        name="M03_drop_1482_historical_matrix_attribution",
        old=(
            "`a0c28978…` over 7,111 profile rows, where 1,482 rows were dropped as\n"
            "non-human (#222)"
        ),
        new="`a0c28978…` over 7,111 profile rows (#222)",
    ),
    Mutant(
        name="M04_restore_unqualified_parenthetical_beside_1500",
        old="none was dropped as\nnon-human, because the matrix holds human profiles only (#656)",
        new=(
            "(the v1.3.1 wide matrix is already human-only upstream, so that filter drops 0 rows "
            "on this run). 1,500 rows were dropped as non-human"
        ),
    ),
    Mutant(
        name="M05_add_unattributed_1500_drop_sentence",
        old="same corpus as the composition arm's 2,257 pairs.",
        new=(
            "same corpus as the composition arm's 2,257 pairs. "
            "1,500 rows were dropped as non-human."
        ),
    ),
    Mutant(
        name="M06_corrupt_canonical_row_count_5629",
        old="Of 5,629 profile rows read under the\n`count_input_16s_wms_v1` profile policy",
        new="Of 7,329 profile rows read under the\n`count_input_16s_wms_v1` profile policy",
    ),
)


class _OriginCheckPlugin:
    def __init__(self, expected_root: Path) -> None:
        self.expected_root = expected_root.resolve()
        self.checked = False

    def pytest_runtest_call(self, item: pytest.Item) -> None:
        module = getattr(item, "module", None)
        if module is None:
            raise RuntimeError("Missing test module")
        actual = Path(module.METHODS_PATH).resolve()
        expected = (self.expected_root / DOC_REL).resolve()
        _require(actual == expected, f"Test read {actual} instead of scratch {expected}")
        self.checked = True


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _run_in_scratch(scratch_root: Path) -> tuple[int, str, bool]:
    plugin = _OriginCheckPlugin(scratch_root)
    buf = io.StringIO()
    for mod_name in [k for k in sys.modules if k.startswith("test_table1_maps_methods_doc_333")]:
        sys.modules.pop(mod_name, None)
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        code = int(
            pytest.main(
                [
                    "-q",
                    "-o",
                    "addopts=",
                    "--rootdir",
                    str(scratch_root),
                    str(scratch_root / TEST_REL),
                ],
                plugins=[plugin],
            )
        )
    return code, buf.getvalue(), plugin.checked


def main() -> int:
    live_doc = REPO_ROOT / DOC_REL
    live_test = REPO_ROOT / TEST_REL
    before_doc = _sha256(live_doc)
    before_test = _sha256(live_test)

    with tempfile.TemporaryDirectory(prefix="verify-333-mutants-") as tmp:
        scratch = Path(tmp) / "repo"
        (scratch / "methods").mkdir(parents=True, exist_ok=True)
        (scratch / "tests").mkdir(parents=True, exist_ok=True)
        shutil.copy2(REPO_ROOT / "pyproject.toml", scratch / "pyproject.toml")
        shutil.copy2(live_doc, scratch / DOC_REL)
        shutil.copy2(live_test, scratch / TEST_REL)

        rc, out, checked = _run_in_scratch(scratch)
        if rc != 0 or not checked:
            sys.stderr.write(f"CONTROL FAILED (rc={rc}, checked={checked}):\n{out}\n")
            return 1
        print("CONTROL: rc=0 PYTEST_ORIGIN_OK")

        scratch_doc = scratch / DOC_REL
        original_text = scratch_doc.read_text(encoding="utf-8")

        for mutant in MUTANTS:
            _require(mutant.old in original_text, f"Snippet not found for {mutant.name}")
            try:
                scratch_doc.write_text(
                    original_text.replace(mutant.old, mutant.new, 1), encoding="utf-8"
                )
                rc, out, checked = _run_in_scratch(scratch)
                _require(checked, f"{mutant.name} did not execute scratch origin check")
                _require(rc != 0, f"{mutant.name} SURVIVED")
                _require(
                    "AssertionError" in out,
                    f"{mutant.name} failed without AssertionError:\n{out}",
                )
                for bad in ("SyntaxError", "ImportError", "ModuleNotFoundError"):
                    _require(bad not in out, f"{mutant.name} hit {bad}:\n{out}")
                print(f"KILLED (AssertionError): {mutant.name}")
            finally:
                scratch_doc.write_text(original_text, encoding="utf-8")

    _require(_sha256(live_doc) == before_doc, "Live methods doc was modified!")
    _require(_sha256(live_test) == before_test, "Live test file was modified!")
    print(f"ALL {len(MUTANTS)}/{len(MUTANTS)} MUTANTS KILLED BY ASSERTIONERROR")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
