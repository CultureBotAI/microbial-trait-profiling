#!/usr/bin/env python3
"""Scratch-copy mutation verifier for #1522 native-influence AnalysisSpec policy.

Purpose:
  Verify that tests/test_downstream_contracts.py and tests/test_manifest_mt_producers.py
  kill source-level mutations of src/mtp/workflow/downstream_contracts.py (#1522) via
  AssertionError (rc=1) when executed against an isolated scratch copy with verified
  import origin inside pytest.

Inputs:
  - src/mtp/workflow/downstream_contracts.py
  - tests/test_downstream_contracts.py
  - tests/test_manifest_mt_producers.py

Outputs:
  - stdout summary of CONTROL (rc=0) and 6 killed mutants (rc=1, AssertionError,
    PYTEST_ORIGIN_OK).

Issue: #1522
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

REPO_ROOT = Path(__file__).resolve().parents[2]
TARGET_REL = Path("src/mtp/workflow/downstream_contracts.py")
TESTS_REL = (
    Path("tests/test_downstream_contracts.py"),
    Path("tests/test_manifest_mt_producers.py"),
)

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
    "KeyError",
    "INTERNALERROR",
    "ERROR collecting",
)


@dataclass(frozen=True)
class Mutant:
    """One source-level mutation in downstream_contracts.py."""

    name: str
    old: str
    new: str


EXPECTED_POLICY_BLOCK = (
    "        multiple_testing_policy=(\n"
    '            "BH_q_in_mt-supp-native-influence-composition_selected_per_scenario_scope_and_"\n'
    '            "mt-supp-native-influence-ora-enrichment_"'
    '"mt-supp-native-influence-ora-depletion_"\n'
    '            "selected_per_scenario_scope_tail_plus_mt-supp-native-influence-shared-sign_"\n'
    '            "executed_focal_per_scenario_scope_reference_or_background_and_"\n'
    '            "before_after_before_matched_after_matched_matched_change_population"\n'
    "        ),"
)

MUTANTS: tuple[Mutant, ...] = (
    Mutant(
        name="M1_revert_to_pre_1522_selected_tail_policy",
        old=(
            "        multiple_testing_policy=(\n"
            '            "BH_q_in_mt-supp-native-influence-composition_'
            'selected_per_scenario_scope_and_"\n'
            '            "mt-supp-native-influence-ora-enrichment_'
            'mt-supp-native-influence-ora-depletion_"\n'
            '            "selected_per_scenario_scope_tail_plus_'
            'mt-supp-native-influence-shared-sign_"\n'
            '            "executed_focal_per_scenario_scope_reference_or_background_and_"\n'
            '            "before_after_before_matched_after_matched_matched_change_population"\n'
            "        ),"
        ),
        new=(
            "        multiple_testing_policy=(\n"
            '            "BH_q_within_selected_native_influence_families_'
            'per_scenario_scope_tail"\n'
            "        ),"
        ),
    ),
    Mutant(
        name="M2_omit_shared_sign_family_id",
        old="mt-supp-native-influence-shared-sign_",
        new="native-influence-shared-sign_",
    ),
    Mutant(
        name="M3_label_shared_sign_as_selected_instead_of_executed_focal",
        old="executed_focal_per_scenario_scope_reference_or_background_and_",
        new="selected_per_scenario_scope_reference_or_background_and_",
    ),
    Mutant(
        name="M4_omit_reference_or_background_grouping",
        old="executed_focal_per_scenario_scope_reference_or_background_and_",
        new="executed_focal_per_scenario_scope_tail_and_",
    ),
    Mutant(
        name="M5_omit_matched_change_population",
        old="before_after_before_matched_after_matched_matched_change_population",
        new="before_after_before_matched_after_matched_population",
    ),
    Mutant(
        name="M6_mutate_eligibility_rule",
        old='eligibility_rule="historical_experiment_paired_contrasts",',
        new='eligibility_rule="paired_human_16S_WMS_nonzero_directional_difference_studies",',
    ),
)


def _setup_scratch(tmp_dir: Path) -> None:
    """Copy minimal source, methods, and test trees into the scratch directory."""
    for rel_dir in ("src/mtp", "src/bsdbng", "methods"):
        src_dir = REPO_ROOT / rel_dir
        if src_dir.exists():
            shutil.copytree(src_dir, tmp_dir / rel_dir)
    for rel_file in (*TESTS_REL, Path("pyproject.toml")):
        dst = tmp_dir / rel_file
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(REPO_ROOT / rel_file, dst)
    conftest = tmp_dir / "tests/conftest.py"
    expected_root = tmp_dir.resolve()
    conftest.write_text(
        "from pathlib import Path\n"
        "import mtp.paths as _paths\n"
        "import mtp.workflow.downstream_contracts as _dc\n\n"
        "def pytest_runtest_setup(item):\n"
        "    actual = _paths.REPO_ROOT.resolve()\n"
        f"    expected = Path({str(expected_root)!r}).resolve()\n"
        "    if actual != expected:\n"
        "        raise RuntimeError(f'wrong origin: {actual} != {expected}')\n"
        "    dc_path = Path(_dc.__file__).resolve()\n"
        "    if not dc_path.is_relative_to(expected / 'src'):\n"
        "        raise RuntimeError(f'wrong dc origin: {dc_path}')\n"
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
                    "--tb=short",
                    "--rootdir",
                    str(tmp_dir),
                    *(str(tmp_dir / rel) for rel in TESTS_REL),
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
    with tempfile.TemporaryDirectory(prefix="mtp_1522_control_") as tmp:
        tmp_dir = Path(tmp)
        _setup_scratch(tmp_dir)
        rc, out = _run_case(tmp_dir)
        if rc != 0 or "PYTEST_ORIGIN_OK" not in out:
            print(f"FAIL [CONTROL]: rc={rc}\n{out}")
            return 1
        print("CONTROL: rc=0 PYTEST_ORIGIN_OK")

    for mutant in MUTANTS:
        with tempfile.TemporaryDirectory(prefix=f"mtp_1522_{mutant.name}_") as tmp:
            tmp_dir = Path(tmp)
            _setup_scratch(tmp_dir)
            target_path = tmp_dir / TARGET_REL
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

    print(f"ALL {len(MUTANTS)} #1522 MUTANTS KILLED CLEANLY.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
