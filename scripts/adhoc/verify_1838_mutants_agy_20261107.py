#!/usr/bin/env python3
"""Scratch-copy mutation verifier for #1838 and PR #1886/#1962 review fixes.

Purpose:
  Verify that the empty-population comparison_rows regression test in
  tests/test_contrast_orientation.py kills source-level mutations of
  comparison_rows in src/mtp/contrast_orientation.py (including revert mutants
  #1889, non-default ratio trait mismatches #1888, docstring regression #1890,
  unknown ratio_traits names #1898, non-Table 1 profile columns #1988, unknown
  ratio trait message #1989, and docstring Table 1 rule #1990) via assertion
  failure (rc=1) when executed against an isolated scratch copy with verified
  import origin inside pytest.

Inputs:
  - src/mtp/contrast_orientation.py
  - tests/test_contrast_orientation.py

Outputs:
  - stdout summary of CONTROL (rc=0) and 19 killed mutants (rc=1, matched
    assertion marker, PYTEST_ORIGIN_OK).

Issue: #1838, #1888, #1889, #1890, #1898, #1988, #1989, #1990, #1991
Date: 2026-11-10
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
CO_REL = Path("src/mtp/contrast_orientation.py")
TEST_CO_REL = Path("tests/test_contrast_orientation.py")
TEST_CO_NAME = (
    "test_comparison_rows_emit_untestable_rows_for_an_empty_population_without_changing_others"
)

ASSERTION_MARKERS = (
    "AssertionError",
    "Failed: DID NOT RAISE",
    "Failed: an empty population must not abort the table:",
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
    """One source-level mutation in contrast_orientation.py."""

    name: str
    old: str
    new: str


MUTANTS: tuple[Mutant, ...] = (
    Mutant(
        name="M01_empty_population_skips_log_ratio_rows_entirely",
        old="        if not profile_list:\n            for numerator, reference in ratio_pairs:",
        new="        if not profile_list:\n            for numerator, reference in ():",
    ),
    Mutant(
        name="M02_empty_population_uses_hardcoded_table1_instead_of_declared_ratio_traits",
        old="    ratio_pairs = list(combinations(declared_ratio_traits, 2))",
        new="    ratio_pairs = list(combinations(TABLE1_TRAITS, 2))",
    ),
    Mutant(
        name="M03_populated_branch_drops_declared_ratio_traits_argument",
        old="        ratios = compute_log_ratio_shifts(profile_list, declared_ratio_traits)",
        new="        ratios = compute_log_ratio_shifts(profile_list)",
    ),
    Mutant(
        name="M04_comparison_rows_accepts_duplicate_traits",
        old=(
            "    if len(declared_ratio_traits) < 2 or len(set(declared_ratio_traits)) != len(\n"
            "        declared_ratio_traits\n"
            "    ):"
        ),
        new="    if len(declared_ratio_traits) < 2:",
    ),
    Mutant(
        name="M05_comparison_rows_accepts_single_trait",
        old=(
            "    if len(declared_ratio_traits) < 2 or len(set(declared_ratio_traits)) != len(\n"
            "        declared_ratio_traits\n"
            "    ):"
        ),
        new="    if len(set(declared_ratio_traits)) != len(declared_ratio_traits):",
    ),
    Mutant(
        name="M06_empty_population_log_ratio_family_label_hardcoded_15",
        old=(
            '                        "family": '
            'f"Benjamini-Yekutieli over {len(ratio_pairs)} ratios",'
        ),
        new='                        "family": "Benjamini-Yekutieli over 15 ratios",',
    ),
    Mutant(
        name="M07_empty_population_log_ratio_p_value_sentinel_one_instead_of_none",
        old='                        "p_value": None,\n                        "q_value": None,',
        new='                        "p_value": 1.0,\n                        "q_value": None,',
    ),
    Mutant(
        name="M08_empty_population_log_ratio_q_value_sentinel_one_instead_of_none",
        old='                        "p_value": None,\n                        "q_value": None,',
        new='                        "p_value": None,\n                        "q_value": 1.0,',
    ),
    Mutant(
        name="M09_empty_population_log_ratio_studies_scored_nonzero",
        old='                        "studies_scored": 0,',
        new='                        "studies_scored": 1,',
    ),
    Mutant(
        name="M10_empty_population_returns_early_dropping_later_populations",
        old=(
            "            continue\n"
            "        ratios = compute_log_ratio_shifts(profile_list, declared_ratio_traits)"
        ),
        new=(
            "            return rows\n"
            "        ratios = compute_log_ratio_shifts(profile_list, declared_ratio_traits)"
        ),
    ),
    Mutant(
        name="M11_empty_population_swaps_numerator_and_reference",
        old=(
            '                        "numerator": numerator,\n'
            '                        "reference": reference,'
        ),
        new=(
            '                        "numerator": reference,\n'
            '                        "reference": numerator,'
        ),
    ),
    Mutant(
        name="M12_revert_empty_branch_1889",
        old="        if not profile_list:\n            for numerator, reference in ratio_pairs:",
        new="        if False:\n            for numerator, reference in ratio_pairs:",
    ),
    Mutant(
        name="M13_empty_check_on_shifts_1889",
        old="        if not profile_list:\n            for numerator, reference in ratio_pairs:",
        new="        if not shifts.shifts:\n            for numerator, reference in ratio_pairs:",
    ),
    Mutant(
        name="M14_docstring_reverts_to_no_paired_profiles_1890",
        old="An\n    empty population emits untestable rows across both families.",
        new=(
            "A\n    population with no paired profiles emits untestable rows across both families."
        ),
    ),
    Mutant(
        name="M15_unknown_ratio_trait_accepted_1898",
        old=(
            "    for trait in declared_ratio_traits:\n"
            "        if trait not in known_ratio_traits:\n"
            "            raise ValueError(\n"
            '                f"Unknown ratio trait {trait}: ratio_traits must be Table 1 traits "\n'
            "                f\"({', '.join(TABLE1_TRAITS)})\"\n"
            "            )"
        ),
        new=(
            "    for trait in ():\n"
            "        if trait not in known_ratio_traits:\n"
            "            raise ValueError(\n"
            '                f"Unknown ratio trait {trait}: ratio_traits must be Table 1 traits "\n'
            "                f\"({', '.join(TABLE1_TRAITS)})\"\n"
            "            )"
        ),
    ),
    Mutant(
        name="M16_docstring_drops_ratio_traits_1898",
        old=(
            "    over the declared fractions and Benjamini-Yekutieli over the declared ratios\n"
            "    (``ratio_traits``, defaulting to ``TABLE1_TRAITS``; every named trait must be a\n"
            "    known trait in ``TABLE1_TRAITS``)."
        ),
        new="    over the declared fractions and Benjamini-Yekutieli over the declared ratios.",
    ),
    Mutant(
        name="M17_profile_columns_also_accepted_1988",
        old="    known_ratio_traits = set(TABLE1_TRAITS)\n",
        new=(
            "    known_ratio_traits = set(TABLE1_TRAITS) | {\n"
            "        name for _, group in populations for p in group for name in p.guild_counts\n"
            "    }\n"
        ),
    ),
    Mutant(
        name="M18_check_only_when_no_population_has_profiles_1988",
        old="    for trait in declared_ratio_traits:\n        if trait not in known_ratio_traits:",
        new=(
            "    for trait in (() if any(len(g) for _, g in populations) "
            "else declared_ratio_traits):\n"
            "        if trait not in known_ratio_traits:"
        ),
    ),
    Mutant(
        name="M19_docstring_keeps_name_drops_table1_rule_1990",
        old=(
            "    (``ratio_traits``, defaulting to ``TABLE1_TRAITS``; every named trait must be a\n"
            "    known trait in ``TABLE1_TRAITS``)."
        ),
        new="    (``ratio_traits``, defaulting to ``TABLE1_TRAITS``).",
    ),
)


def _setup_scratch(tmp_dir: Path) -> None:
    """Copy minimal source and test trees into the scratch directory."""
    for rel_dir in ("src/mtp", "src/bsdbng", "methods"):
        src_dir = REPO_ROOT / rel_dir
        if src_dir.exists():
            shutil.copytree(src_dir, tmp_dir / rel_dir)
    (tmp_dir / "tests").mkdir(parents=True, exist_ok=True)
    shutil.copy2(REPO_ROOT / TEST_CO_REL, tmp_dir / TEST_CO_REL)
    conftest = tmp_dir / "tests/conftest.py"
    expected_co = (tmp_dir / CO_REL).resolve()
    conftest.write_text(
        "from pathlib import Path\n"
        "import mtp.contrast_orientation as _co\n\n"
        "def pytest_runtest_setup(item):\n"
        "    actual = Path(_co.__file__).resolve()\n"
        f"    expected = Path({str(expected_co)!r}).resolve()\n"
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
                    str(tmp_dir / TEST_CO_REL),
                    "-k",
                    TEST_CO_NAME,
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
    with tempfile.TemporaryDirectory(prefix="mtp_1838_control_") as tmp:
        tmp_dir = Path(tmp)
        _setup_scratch(tmp_dir)
        rc, out = _run_case(tmp_dir)
        if rc != 0 or "PYTEST_ORIGIN_OK" not in out:
            print(f"FAIL [CONTROL]: rc={rc}\n{out}")
            return 1
        print("CONTROL: rc=0 PYTEST_ORIGIN_OK")

    for mutant in MUTANTS:
        with tempfile.TemporaryDirectory(prefix=f"mtp_1838_{mutant.name}_") as tmp:
            tmp_dir = Path(tmp)
            _setup_scratch(tmp_dir)
            target_path = tmp_dir / CO_REL
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
            matched_marker = next(m for m in ASSERTION_MARKERS if m in out)
            print(f"KILLED {mutant.name}: rc=1 {matched_marker} PYTEST_ORIGIN_OK")

    print(f"ALL {len(MUTANTS)} #1838 MUTANTS KILLED CLEANLY.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
