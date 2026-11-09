#!/usr/bin/env python3
"""Scratch-copy mutation verifier for #1852 (unselected testable traits excluded from native BH).

Purpose:
  Verify that tests/test_native_taxon_influence.py kills source-level mutations of
  the eligibility and valid-test gates in src/mtp/native_taxon_influence.py
  (including the PR #1836 reviewer's surviving mutant d-code-test-unselected)
  via AssertionError (rc=1) when executed against an isolated scratch copy with
  verified import origin inside pytest.

Inputs:
  - src/mtp/native_taxon_influence.py
  - tests/test_native_taxon_influence.py

Outputs:
  - stdout summary of CONTROL (rc=0) and 11 killed mutants (rc=1, AssertionError,
    PYTEST_ORIGIN_OK).

Issue: #1852, #1892, #1893
Date: 2026-11-07
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

REPO_ROOT = Path(__file__).resolve().parents[2]
NTI_REL = Path("src/mtp/native_taxon_influence.py")
TEST_NTI_REL = Path("tests/test_native_taxon_influence.py")

ASSERTION_MARKERS = (
    "AssertionError",
    "Failed: DID NOT RAISE",
    "Regex pattern did not match",
)
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


@dataclass(frozen=True)
class Mutant:
    """One source-level mutation in native_taxon_influence.py."""

    name: str
    old: str
    new: str


MUTANTS: tuple[Mutant, ...] = (
    Mutant(
        name="M01_d_code_test_unselected_drop_eligible_from_valid",
        old='        valid = eligible and len(items) >= policy["minimum_pairs"] and sd > 1e-8',
        new='        valid = len(items) >= policy["minimum_pairs"] and sd > 1e-8',
    ),
    Mutant(
        name="M02_eligible_ignores_minimum_profile_support",
        old=(
            "        eligible = (\n"
            '            support >= policy["minimum_profile_support"]\n'
            "            and support <= max("
            '2, int(len(keys) * policy["maximum_profile_support_fraction"]))\n'
        ),
        new=(
            "        eligible = (\n"
            "            support <= max("
            '2, int(len(keys) * policy["maximum_profile_support_fraction"]))\n'
        ),
    ),
    Mutant(
        name="M03_eligible_ignores_maximum_profile_support_fraction",
        old=(
            "        eligible = (\n"
            '            support >= policy["minimum_profile_support"]\n'
            "            and support <= max("
            '2, int(len(keys) * policy["maximum_profile_support_fraction"]))\n'
            "            and var > 0\n"
            "        )"
        ),
        new=(
            "        eligible = (\n"
            '            support >= policy["minimum_profile_support"]\n'
            "            and var > 0\n"
            "        )"
        ),
    ),
    Mutant(
        name="M04_status_labels_unselected_as_untestable",
        old=(
            '            "status": "tested" if valid '
            'else "not_selected" if not eligible else "untestable",'
        ),
        new='            "status": "tested" if valid else "untestable",',
    ),
    Mutant(
        name="M05_unselected_row_given_sentinel_q_one",
        old='        if row["selected"] and row["p"] == "":\n            row["q"] = 1.0',
        new='        if row["p"] == "":\n            row["q"] = 1.0',
    ),
    Mutant(
        name="M06_family_selected_traits_counts_all_traits",
        old='        row["family_selected_traits"] = sum(r["selected"] for r in records.values())',
        new='        row["family_selected_traits"] = len(records)',
    ),
    Mutant(
        name="M07_family_tests_counts_all_traits_instead_of_tested",
        old='        row["family_tests"] = len(tested)',
        new='        row["family_tests"] = len(records)',
    ),
    Mutant(
        name="M08_family_tests_counts_selected_instead_of_tested",
        old='        row["family_tests"] = len(tested)',
        new='        row["family_tests"] = sum(r["selected"] for r in records.values())',
    ),
    Mutant(
        name="M09_p_computed_for_ineligible_row_even_when_valid_false",
        old='            "p": float(ttest_1samp(items, 0.0).pvalue) if valid else "",',
        new=(
            '            "p": float(ttest_1samp(items, 0.0).pvalue) '
            'if (len(items) >= policy["minimum_pairs"] and sd > 1e-8) else "",'
        ),
    ),
    Mutant(
        name="M10_bh_plus_unselected_testable_1892",
        old='    for row, q in zip(tested, adjust([r["p"] for r in tested]), strict=True):',
        new=(
            "    for row, q in zip(\n"
            "        tested,\n"
            "        adjust(\n"
            '            [r["p"] for r in tested]\n'
            "            + [\n"
            "                float(ttest_1samp(list(contrasts[k].values()), 0.0).pvalue)\n"
            "                for k, r in records.items()\n"
            '                if not r["selected"]\n'
            '                and len(contrasts[k]) >= policy["minimum_pairs"]\n'
            "                and float(np.std(list(contrasts[k].values()), ddof=1)) > 1e-8\n"
            "            ]\n"
            "        )[: len(tested)],\n"
            "        strict=True,\n"
            "    ):"
        ),
    ),
    Mutant(
        name="M11_variance_nonstrict_1893",
        old="            and var > 0\n        )",
        new="            and var >= 0\n        )",
    ),
)


def _setup_scratch(tmp_dir: Path) -> None:
    """Copy minimal source and test trees into the scratch directory."""
    for rel_dir in ("src/mtp", "src/bsdbng", "methods"):
        src_dir = REPO_ROOT / rel_dir
        if src_dir.exists():
            shutil.copytree(src_dir, tmp_dir / rel_dir)
    (tmp_dir / "tests").mkdir(parents=True, exist_ok=True)
    shutil.copy2(REPO_ROOT / TEST_NTI_REL, tmp_dir / TEST_NTI_REL)
    conftest = tmp_dir / "tests/conftest.py"
    expected_nti = (tmp_dir / NTI_REL).resolve()
    conftest.write_text(
        "from pathlib import Path\n"
        "import mtp.native_taxon_influence as _nti\n\n"
        "def pytest_runtest_setup(item):\n"
        "    actual = Path(_nti.__file__).resolve()\n"
        f"    expected = Path({str(expected_nti)!r}).resolve()\n"
        "    assert actual == expected, f'wrong origin: {actual} != {expected}'\n"
        "    print('PYTEST_ORIGIN_OK', flush=True)\n",
        encoding="utf-8",
    )


def _worker(tmp_dir_str: str, out_file_str: str) -> None:
    """Run pytest in a clean child process rooted at tmp_dir."""
    tmp_dir = Path(tmp_dir_str)
    os.chdir(tmp_dir)
    sys.path.insert(0, str(tmp_dir / "src"))
    for name in [k for k in list(sys.modules) if k == "mtp" or k.startswith(("mtp.", "bsdbng"))]:
        sys.modules.pop(name, None)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        rc = int(
            pytest.main(
                [
                    "-q",
                    "-s",
                    "--tb=long",
                    "--rootdir",
                    str(tmp_dir),
                    str(tmp_dir / TEST_NTI_REL),
                ]
            )
        )
    Path(out_file_str).write_text(f"RC={rc}\n{buf.getvalue()}", encoding="utf-8")


def _run_case(tmp_dir: Path) -> tuple[int, str]:
    out_file = tmp_dir / "pytest_out.txt"
    ctx = multiprocessing.get_context("spawn")
    proc = ctx.Process(target=_worker, args=(str(tmp_dir), str(out_file)))
    proc.start()
    proc.join()
    if proc.exitcode != 0 or not out_file.exists():
        raise RuntimeError(f"worker process failed with exitcode={proc.exitcode}")
    raw = out_file.read_text(encoding="utf-8")
    first, _, rest = raw.partition("\n")
    return int(first.removeprefix("RC=")), rest


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="mtp_1852_control_") as tmp:
        tmp_dir = Path(tmp)
        _setup_scratch(tmp_dir)
        rc, out = _run_case(tmp_dir)
        if rc != 0 or "PYTEST_ORIGIN_OK" not in out:
            print(f"FAIL [CONTROL]: rc={rc}\n{out}")
            return 1
        print("CONTROL: rc=0 PYTEST_ORIGIN_OK")

    for mutant in MUTANTS:
        with tempfile.TemporaryDirectory(prefix=f"mtp_1852_{mutant.name}_") as tmp:
            tmp_dir = Path(tmp)
            _setup_scratch(tmp_dir)
            target_path = tmp_dir / NTI_REL
            source = target_path.read_text(encoding="utf-8")
            count = source.count(mutant.old)
            if count != 1:
                print(f"FAIL [{mutant.name}]: expected 1 occurrence, found {count}")
                return 1
            target_path.write_text(source.replace(mutant.old, mutant.new, 1), encoding="utf-8")
            rc, out = _run_case(tmp_dir)
            if rc != 1:
                print(f"FAIL [{mutant.name}]: expected rc=1, got rc={rc}\n{out}")
                return 1
            if "PYTEST_ORIGIN_OK" not in out:
                print(f"FAIL [{mutant.name}]: missing PYTEST_ORIGIN_OK\n{out}")
                return 1
            if any(err in out for err in ERROR_MARKERS):
                print(f"FAIL [{mutant.name}]: non-assertion error marker in output\n{out}")
                return 1
            if not any(marker in out for marker in ASSERTION_MARKERS):
                print(f"FAIL [{mutant.name}]: no assertion marker in output\n{out}")
                return 1
            print(f"KILLED {mutant.name}: rc=1 AssertionError PYTEST_ORIGIN_OK")

    print(f"ALL {len(MUTANTS)} #1852 MUTANTS KILLED CLEANLY.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
