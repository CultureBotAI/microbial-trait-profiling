"""Scratch-copy mutation verifier for #1808 (non-finite McNemar p in readers).

Purpose:
  Verify that the regression tests in tests/test_integrated_evidence.py and
  tests/test_calibration.py kill source-level mutations of the #1808 non-finite
  McNemar p guards in src/mtp/integrated_evidence.py and src/mtp/calibration.py
  via AssertionError (rc=1) when executed against an isolated scratch copy with
  verified import origin inside pytest.

Inputs:
  - src/mtp/integrated_evidence.py
  - src/mtp/calibration.py
  - tests/test_integrated_evidence.py
  - tests/test_calibration.py

Outputs:
  - stdout summary of CONTROL (rc=0) and 12 killed mutants (rc=1, AssertionError,
    PYTEST_ORIGIN_OK).

Issue: #1808
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
IE_REL = Path("src/mtp/integrated_evidence.py")
CAL_REL = Path("src/mtp/calibration.py")
TEST_IE_REL = Path("tests/test_integrated_evidence.py")
TEST_CAL_REL = Path("tests/test_calibration.py")

TEST_IE_NAME = "test_a_nonfinite_mcnemar_p_is_not_tested_and_leaves_p_blank"
TEST_CAL_NAME = "test_a_nonfinite_mcnemar_p_is_not_counted_as_a_tested_replicate"


@dataclass(frozen=True)
class Mutant:
    """One source-level mutation in a target module and the test that must kill it."""

    name: str
    target_rel: Path
    test_rel: Path
    test_filter: str
    old: str
    new: str


MUTANTS: tuple[Mutant, ...] = (
    Mutant(
        name="M01_ie_status_always_tested_when_discordant",
        target_rel=IE_REL,
        test_rel=TEST_IE_REL,
        test_filter=TEST_IE_NAME,
        old=(
            "    elif mcnemar_p is None:\n"
            '        status = "not_tested"\n'
            "    else:\n"
            "        status = TESTED"
        ),
        new="    else:\n        status = TESTED",
    ),
    Mutant(
        name="M02_ie_status_insufficient_instead_of_not_tested",
        target_rel=IE_REL,
        test_rel=TEST_IE_REL,
        test_filter=TEST_IE_NAME,
        old='    elif mcnemar_p is None:\n        status = "not_tested"',
        new="    elif mcnemar_p is None:\n        status = INSUFFICIENT",
    ),
    Mutant(
        name="M03_ie_drop_isfinite_check_on_mcnemar_p",
        target_rel=IE_REL,
        test_rel=TEST_IE_REL,
        test_filter=TEST_IE_NAME,
        old=(
            '        if correlation.status == "tested"\n'
            "        and correlation.mcnemar_p is not None\n"
            "        and math.isfinite(correlation.mcnemar_p)"
        ),
        new="        if correlation.mcnemar_p is not None",
    ),
    Mutant(
        name="M04_ie_drop_correlation_status_check",
        target_rel=IE_REL,
        test_rel=TEST_IE_REL,
        test_filter=TEST_IE_NAME,
        old=(
            '        if correlation.status == "tested"\n'
            "        and correlation.mcnemar_p is not None\n"
            "        and math.isfinite(correlation.mcnemar_p)"
        ),
        new="        if correlation.mcnemar_p is not None and math.isfinite(correlation.mcnemar_p)",
    ),
    Mutant(
        name="M05_ie_pass_raw_correlation_mcnemar_p_bypassing_untested",
        target_rel=IE_REL,
        test_rel=TEST_IE_REL,
        test_filter=TEST_IE_NAME,
        old=(
            "            selection_note=GUILD_SELECTION_NOTE,\n"
            "            source=source,\n"
            "        )\n"
            "    ]\n"
            "    return [untested(record) for record in records]"
        ),
        new=(
            "            selection_note=GUILD_SELECTION_NOTE,\n"
            "            source=source,\n"
            "        )\n"
            "    ]\n"
            "    return list(records)"
        ),
    ),
    Mutant(
        name="M06_cal_guild_mcnemar_drop_isfinite_check",
        target_rel=CAL_REL,
        test_rel=TEST_CAL_REL,
        test_filter=TEST_CAL_NAME,
        old=(
            '            if correlation.status == "tested"\n'
            "            and correlation.mcnemar_p is not None\n"
            "            and math.isfinite(correlation.mcnemar_p)"
        ),
        new="            if correlation.mcnemar_p is not None",
    ),
    Mutant(
        name="M07_cal_guild_mcnemar_drop_status_check",
        target_rel=CAL_REL,
        test_rel=TEST_CAL_REL,
        test_filter=TEST_CAL_NAME,
        old=(
            '            if correlation.status == "tested"\n'
            "            and correlation.mcnemar_p is not None\n"
            "            and math.isfinite(correlation.mcnemar_p)"
        ),
        new=(
            "            if correlation.mcnemar_p is not None "
            "and math.isfinite(correlation.mcnemar_p)"
        ),
    ),
    Mutant(
        name="M08_cal_summarize_present_drop_isfinite",
        target_rel=CAL_REL,
        test_rel=TEST_CAL_REL,
        test_filter=TEST_CAL_NAME,
        old=(
            "            present = [\n"
            "                o\n"
            "                for o in tested\n"
            "                if o is not None\n"
            "                and o.selected\n"
            "                and o.p_value is not None\n"
            "                and math.isfinite(o.p_value)\n"
            "            ]"
        ),
        new=(
            "            present = [\n"
            "                o\n"
            "                for o in tested\n"
            "                if o is not None\n"
            "                and o.selected\n"
            "                and o.p_value is not None\n"
            "            ]"
        ),
    ),
    Mutant(
        name="M09_cal_summarize_rejected_q_drop_isfinite",
        target_rel=CAL_REL,
        test_rel=TEST_CAL_REL,
        test_filter=TEST_CAL_NAME,
        old=(
            "                if o.q_value is not None "
            "and math.isfinite(o.q_value) and o.q_value <= Q_THRESHOLD"
        ),
        new="                if o.q_value is not None and o.q_value <= Q_THRESHOLD",
    ),
    Mutant(
        name="M10_cal_summarize_family_calls_drop_isfinite",
        target_rel=CAL_REL,
        test_rel=TEST_CAL_REL,
        test_filter=TEST_CAL_NAME,
        old=(
            "                and o.p_value is not None\n"
            "                and math.isfinite(o.p_value)\n"
            "                and o.q_value is not None\n"
            "                and math.isfinite(o.q_value)\n"
            "                and o.q_value <= Q_THRESHOLD"
        ),
        new=(
            "                and o.q_value is not None\n"
            "                and o.q_value <= Q_THRESHOLD"
        ),
    ),
    Mutant(
        name="M11_cal_replicate_counts_tested_drop_isfinite",
        target_rel=CAL_REL,
        test_rel=TEST_CAL_REL,
        test_filter=TEST_CAL_NAME,
        old=(
            "            tested = [o for o in mine "
            "if o.p_value is not None and math.isfinite(o.p_value)]"
        ),
        new="            tested = [o for o in mine if o.p_value is not None]",
    ),
    Mutant(
        name="M12_cal_replicate_counts_family_rejected_q_drop_isfinite",
        target_rel=CAL_REL,
        test_rel=TEST_CAL_REL,
        test_filter=TEST_CAL_NAME,
        old=(
            '                    "family_rejected_q": int(\n'
            "                        any(\n"
            "                            o.q_value is not None\n"
            "                            and math.isfinite(o.q_value)\n"
            "                            and o.q_value <= Q_THRESHOLD\n"
            "                            for o in tested\n"
            "                        )\n"
            "                    ),"
        ),
        new=(
            '                    "family_rejected_q": int(\n'
            "                        any(\n"
            "                            o.q_value is not None and o.q_value <= Q_THRESHOLD\n"
            "                            for o in mine\n"
            "                        )\n"
            "                    ),"
        ),
    ),
)


def _setup_scratch(tmp_dir: Path) -> None:
    """Copy minimal source and test trees into the scratch directory."""
    for rel_dir in ("src/mtp", "src/bsdbng", "methods"):
        src_dir = REPO_ROOT / rel_dir
        if src_dir.exists():
            shutil.copytree(src_dir, tmp_dir / rel_dir)
    (tmp_dir / "tests").mkdir(parents=True, exist_ok=True)
    for rel_file in (TEST_IE_REL, TEST_CAL_REL):
        shutil.copy2(REPO_ROOT / rel_file, tmp_dir / rel_file)
    fixtures_src = REPO_ROOT / "tests/fixtures"
    if fixtures_src.exists():
        shutil.copytree(fixtures_src, tmp_dir / "tests/fixtures")
    conftest = tmp_dir / "tests/conftest.py"
    expected_ie = (tmp_dir / IE_REL).resolve()
    expected_cal = (tmp_dir / CAL_REL).resolve()
    conftest.write_text(
        "from pathlib import Path\n"
        "import mtp.calibration as cal\n"
        "import mtp.integrated_evidence as ie\n"
        "def pytest_sessionstart(session):\n"
        "    actual_ie = Path(ie.__file__).resolve()\n"
        "    actual_cal = Path(cal.__file__).resolve()\n"
        f"    expected_ie = Path({str(expected_ie)!r})\n"
        f"    expected_cal = Path({str(expected_cal)!r})\n"
        "    assert actual_ie == expected_ie\n"
        "    assert actual_cal == expected_cal\n"
        "    print('PYTEST_ORIGIN_OK', actual_ie, actual_cal)\n",
        encoding="utf-8",
    )


def _pytest_worker(
    tmp_dir_str: str, test_rel_str: str, test_filter: str, out_file_str: str
) -> None:
    tmp_dir = Path(tmp_dir_str)
    os.chdir(tmp_dir)
    sys.path.insert(0, str(tmp_dir / "src"))
    for mod_name in [
        k for k in list(sys.modules) if k == "mtp" or k.startswith(("mtp.", "bsdbng"))
    ]:
        sys.modules.pop(mod_name, None)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        rc = int(
            pytest.main(
                [
                    "-q",
                    "-s",
                    "-o",
                    "addopts=",
                    "--rootdir",
                    str(tmp_dir),
                    str(tmp_dir / test_rel_str),
                    "-k",
                    test_filter,
                ]
            )
        )
    Path(out_file_str).write_text(f"RC={rc}\n{buf.getvalue()}", encoding="utf-8")


def _run_pytest(
    ctx: multiprocessing.context.SpawnContext,
    tmp_dir: Path,
    test_rel: Path,
    test_filter: str,
) -> tuple[int, str]:
    out_file = tmp_dir / "pytest_out.txt"
    proc = ctx.Process(
        target=_pytest_worker,
        args=(str(tmp_dir), str(test_rel), test_filter, str(out_file)),
    )
    proc.start()
    proc.join()
    if proc.exitcode != 0 or not out_file.exists():
        raise RuntimeError(f"worker process failed with exitcode={proc.exitcode}")
    raw = out_file.read_text(encoding="utf-8")
    first_line, _, combined = raw.partition("\n")
    rc = int(first_line.removeprefix("RC="))
    return rc, combined


def main() -> int:
    ctx = multiprocessing.get_context("spawn")
    originals = {
        IE_REL: (REPO_ROOT / IE_REL).read_text(encoding="utf-8"),
        CAL_REL: (REPO_ROOT / CAL_REL).read_text(encoding="utf-8"),
    }
    with tempfile.TemporaryDirectory(prefix="verify_1808_mutants_") as raw_tmp:
        tmp_dir = Path(raw_tmp)
        _setup_scratch(tmp_dir)

        for label, test_rel, test_filter in (
            ("CONTROL_IE", TEST_IE_REL, TEST_IE_NAME),
            ("CONTROL_CAL", TEST_CAL_REL, TEST_CAL_NAME),
        ):
            rc, combined = _run_pytest(ctx, tmp_dir, test_rel, test_filter)
            if rc != 0 or "PYTEST_ORIGIN_OK" not in combined:
                print(f"{label} FAILED (rc={rc}):\n{combined}")
                return 1
            print(f"{label}: rc={rc} PYTEST_ORIGIN_OK")

        for mutant in MUTANTS:
            original_source = originals[mutant.target_rel]
            if original_source.count(mutant.old) != 1:
                print(
                    f"MUTANT TARGET NOT UNIQUE: {mutant.name} "
                    f"(count={original_source.count(mutant.old)})"
                )
                return 1
            mutated = original_source.replace(mutant.old, mutant.new, 1)
            target_path = tmp_dir / mutant.target_rel
            target_path.write_text(mutated, encoding="utf-8")
            try:
                rc, combined = _run_pytest(ctx, tmp_dir, mutant.test_rel, mutant.test_filter)
            finally:
                target_path.write_text(original_source, encoding="utf-8")

            has_origin = "PYTEST_ORIGIN_OK" in combined
            has_assertion = "AssertionError" in combined
            bad_exc = any(
                err in combined for err in ("SyntaxError", "ImportError", "NameError", "TypeError")
            )
            if rc != 1 or not has_origin or not has_assertion or bad_exc:
                print(
                    f"MUTANT {mutant.name} NOT CLEANLY KILLED "
                    f"(rc={rc}, origin={has_origin}, "
                    f"assertion={has_assertion}, bad_exc={bad_exc}):\n{combined}"
                )
                return 1
            print(f"KILLED {mutant.name}: rc={rc} AssertionError PYTEST_ORIGIN_OK")

    print(f"ALL {len(MUTANTS)} MUTANTS KILLED CLEANLY")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
