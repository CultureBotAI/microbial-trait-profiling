"""Isolated scratch-copy mutation verifier for Issue #1258 (generator audited pins).

Purpose:
  Verify that tests/test_issue1231_generator_pins.py kills 10 distinct mutants of
  scripts/adhoc/issue1231_generate_historical_union_core_fixture.py in an isolated
  scratch directory, and verify scratch origin binding.
Inputs:
  - scripts/adhoc/issue1231_generate_historical_union_core_fixture.py
  - tests/fixtures/trait_key_policy/historical_union_core.tar.gz
  - tests/test_issue1231_generator_pins.py
Outputs:
  - Stdout summary table of control + 10 mutant runs; exits 0 iff all 10 mutants
    fail with AssertionError.
Issue: #1258
Date: 2026-11-09
"""

from __future__ import annotations

import contextlib
import io
import shutil
import sys
import tempfile
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
REL_FILES: tuple[str, ...] = (
    "scripts/adhoc/issue1231_generate_historical_union_core_fixture.py",
    "tests/fixtures/trait_key_policy/historical_union_core.tar.gz",
    "tests/test_issue1231_generator_pins.py",
)

MUTANTS: tuple[tuple[str, str, str], ...] = (
    (
        "M01_reintroduce_chunk_joined_commit_literal",
        'EXPECTED_MAIN_REFERENCE = "44897c6"\n',
        'EXPECTED_MAIN_REFERENCE = "".join(("4489", "7c6f"))\n',
    ),
    (
        "M02_reintroduce_chunk_joined_sha256_literal",
        "EXPECTED_MAIN_COMMIT, EXPECTED_MAIN_MATRIX_SHA256, EXPECTED_MAIN_MATRIX_SIZE = (\n"
        "    load_audited_producer_pins(COMMITTED_FIXTURE)\n"
        ")\n",
        "EXPECTED_MAIN_COMMIT, EXPECTED_MAIN_MATRIX_SHA256, EXPECTED_MAIN_MATRIX_SIZE = (\n"
        "    load_audited_producer_pins(COMMITTED_FIXTURE)\n"
        ")\n"
        'CHUNKED_DIGEST = "".join(("0517", "bd05", "a844", "81b4"))\n',
    ),
    (
        "M03_change_expected_main_reference",
        'EXPECTED_MAIN_REFERENCE = "44897c6"\n',
        'EXPECTED_MAIN_REFERENCE = "511c891"\n',
    ),
    (
        "M04_drop_commit_prefix_check",
        "        or not producer_commit.startswith(EXPECTED_MAIN_REFERENCE)\n",
        "        or False\n",
    ),
    (
        "M05_drop_git_commit_equality_check",
        '    if run_record.get("git_commit") != producer_commit:\n',
        '    if False and run_record.get("git_commit") != producer_commit:\n',
    ),
    (
        "M06_drop_matrix_path_check",
        '    if not isinstance(matrix_inv, dict) or matrix_inv.get("path") != MATRIX_SOURCE:\n',
        "    if not isinstance(matrix_inv, dict):\n",
    ),
    (
        "M07_drop_code_inventory_sha256_check",
        '        matching_entries[0].get("sha256") != matrix_sha256\n',
        "        False\n",
    ),
    (
        "M08_drop_code_inventory_size_check",
        '        or matching_entries[0].get("size_bytes") != matrix_size\n',
        "        or False\n",
    ),
    (
        "M09_override_expected_main_matrix_size_constant",
        "    return producer_commit, matrix_sha256, matrix_size\n",
        "    return producer_commit, matrix_sha256, matrix_size + 1\n",
    ),
    (
        "M10_corrupt_returned_commit_constant",
        "EXPECTED_MAIN_COMMIT, EXPECTED_MAIN_MATRIX_SHA256, EXPECTED_MAIN_MATRIX_SIZE = (\n"
        "    load_audited_producer_pins(COMMITTED_FIXTURE)\n"
        ")\n",
        "EXPECTED_MAIN_COMMIT, EXPECTED_MAIN_MATRIX_SHA256, EXPECTED_MAIN_MATRIX_SIZE = (\n"
        "    load_audited_producer_pins(COMMITTED_FIXTURE)\n"
        ")\n"
        'EXPECTED_MAIN_COMMIT = "44897c6" + "0" * 33\n',
    ),
)


def _populate_scratch(scratch_root: Path) -> None:
    for rel in REL_FILES:
        src = REPO_ROOT / rel
        dst = scratch_root / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


def _run_pytest_in_scratch(scratch_root: Path) -> tuple[int, str]:
    for mod_name in list(sys.modules):
        if mod_name in ("test_issue1231_generator_pins", "issue1231_generator") or (
            mod_name.startswith("tests.test_issue1231_generator_pins")
        ):
            sys.modules.pop(mod_name, None)

    buf = io.StringIO()
    test_file = scratch_root / "tests" / "test_issue1231_generator_pins.py"
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        rc = int(
            pytest.main(
                [
                    "-q",
                    "--tb=native",
                    "-o",
                    "addopts=",
                    str(test_file),
                ]
            )
        )
    return rc, buf.getvalue()


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="verify_1258_mutants_") as tmp_dir:
        scratch_root = Path(tmp_dir)
        _populate_scratch(scratch_root)

        control_rc, control_out = _run_pytest_in_scratch(scratch_root)
        if control_rc != 0:
            print("CONTROL FAILED UNEXPECTEDLY:\n" + control_out)
            return 1

        test_path = scratch_root / "tests" / "test_issue1231_generator_pins.py"
        original_test_text = test_path.read_text(encoding="utf-8")
        test_path.write_text(
            original_test_text
            + "\n\ndef test_scratch_origin_binding() -> None:\n"
            + f"    assert REPO_ROOT == Path({str(scratch_root)!r}).resolve()\n",
            encoding="utf-8",
        )
        origin_rc, origin_out = _run_pytest_in_scratch(scratch_root)
        test_path.write_text(original_test_text, encoding="utf-8")
        if origin_rc != 0:
            print("SCRATCH ORIGIN BINDING FAILED:\n" + origin_out)
            return 1

        print(f"Control + scratch-origin binding passed in {scratch_root}")

        killed = 0
        rel_target = "scripts/adhoc/issue1231_generate_historical_union_core_fixture.py"
        for name, target, replacement in MUTANTS:
            _populate_scratch(scratch_root)
            file_path = scratch_root / rel_target
            original = file_path.read_text(encoding="utf-8")
            if target not in original:
                raise RuntimeError(f"{name}: target {target!r} not found in {rel_target}")
            mutated = original.replace(target, replacement, 1)
            if mutated == original:
                raise RuntimeError(f"{name}: replacement produced no change in {rel_target}")
            file_path.write_text(mutated, encoding="utf-8")

            rc, out = _run_pytest_in_scratch(scratch_root)
            has_assertion_error = (
                "AssertionError" in out
                or "Failed: DID NOT RAISE" in out
                or "ValueError: Audited" in out
            )
            has_bad_error = any(
                err in out for err in ("SyntaxError", "ImportError", "ModuleNotFoundError")
            )
            if rc == 1 and has_assertion_error and not has_bad_error:
                killed += 1
                print(f"  [KILLED] {name} (rc={rc})")
            else:
                print(f"  [SURVIVED/INVALID] {name} (rc={rc}):\n{out}")
                return 1

        print(f"Summary: {killed}/{len(MUTANTS)} mutants killed.")
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
