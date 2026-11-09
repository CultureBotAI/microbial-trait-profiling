#!/usr/bin/env python3
"""Scratch-copy mutation verifier for #1290 and #1027 portable evidence helpers.

Purpose:
  Verify that tests/test_uv_lock_audit.py kills source-level mutations of
  scripts/adhoc/issue1286_compare_uv_locks.py (#1290) and
  scripts/adhoc/issue928_guard_isolation.py (#1027) via AssertionError (rc=1)
  when executed against an isolated scratch copy with verified import origin
  inside pytest.

Inputs:
  - scripts/adhoc/issue1286_compare_uv_locks.py
  - scripts/adhoc/issue928_guard_isolation.py
  - tests/test_uv_lock_audit.py

Outputs:
  - stdout summary of CONTROL (rc=0) and 12 killed mutants (rc=1, AssertionError,
    PYTEST_ORIGIN_OK).

Issue: #1290, #1027, #2018, #2019, #2020, #2021
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
LOCK_REL = Path("scripts/adhoc/issue1286_compare_uv_locks.py")
PROBE_REL = Path("scripts/adhoc/issue928_guard_isolation.py")
TEST_REL = Path("tests/test_uv_lock_audit.py")

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
    """One source-level mutation in an evidence helper."""

    name: str
    rel_path: Path
    old: str
    new: str


MUTANTS: tuple[Mutant, ...] = (
    Mutant(
        name="M1_compare_locks_emits_resolved_path_1290",
        rel_path=LOCK_REL,
        old='            "path": portable.encode(path),',
        new='            "path": f"/abs/{path.name}",',
    ),
    Mutant(
        name="M2_compare_locks_omits_path_basis_1290",
        rel_path=LOCK_REL,
        old='        "path_basis": portable.declaration,\n',
        new="",
    ),
    Mutant(
        name="M3_guard_probe_emits_canonical_origin_absolute_1027",
        rel_path=PROBE_REL,
        old='        "canonical_origin": repo_paths.encode(PRODUCER),',
        new='        "canonical_origin": f"/abs/{PRODUCER.name}",',
    ),
    Mutant(
        name="M4_guard_probe_emits_isolated_origin_absolute_1027",
        rel_path=PROBE_REL,
        old='        "isolated_origin": output_paths.encode(copied),',
        new='        "isolated_origin": f"/abs/{copied.name}",',
    ),
    Mutant(
        name="M5_guard_probe_skips_output_name_account_check_1027",
        rel_path=PROBE_REL,
        old="    assert_portable_path_text(output.name)\n",
        new="",
    ),
    Mutant(
        name="M6_guard_probe_skips_provenance_root_precheck_1027",
        rel_path=PROBE_REL,
        old=(
            "    repo_root = args.provenance_root.resolve() "
            "if args.provenance_root is not None else ROOT\n"
        ),
        new="    repo_root = ROOT\n",
    ),
    Mutant(
        name="M7_lock_record_ignores_output_1290",
        rel_path=LOCK_REL,
        old="        output=args.output,\n",
        new="        output=None,\n",
    ),
    Mutant(
        name="M8_lock_no_result_check_2020",
        rel_path=LOCK_REL,
        old="    assert_portable_text(json.dumps(result, ensure_ascii=False, allow_nan=False))\n",
        new="",
    ),
    Mutant(
        name="M9_lock_ignores_provenance_root_1290",
        rel_path=LOCK_REL,
        old="[before, after], root=provenance_root)",
        new="[before, after], root=None)",
    ),
    Mutant(
        name="M10_guard_probe_skips_log_scrubbing_2018",
        rel_path=PROBE_REL,
        old='    (output / "pytest.log").write_text(log_text, encoding="utf-8")\n',
        new='    (output / "pytest.log").write_text(captured.getvalue(), encoding="utf-8")\n',
    ),
    Mutant(
        name="M11_guard_probe_omits_isolated_path_basis_2019",
        rel_path=PROBE_REL,
        old='        "isolated_path_basis": output_paths.declaration,\n',
        new="",
    ),
    Mutant(
        name="M12_guard_probe_no_isolation_text_check_1027",
        rel_path=PROBE_REL,
        old="    assert_portable_text(isolation_text)\n",
        new="",
    ),
    Mutant(
        name="M13_guard_probe_checkout_outputs_ignores_root_2054",
        rel_path=PROBE_REL,
        old='        "checkout_outputs": repo_paths.encode(checkout_outputs),',
        new='        "checkout_outputs": checkout_outputs.relative_to(ROOT).as_posix(),',
    ),
    Mutant(
        name="M14_guard_probe_unbounded_relative_replace_2055",
        rel_path=PROBE_REL,
        old=(
            '        if target.startswith("/"):\n'
            "            scrubbed = scrubbed.replace(target, placeholder)\n"
        ),
        new="        if True:\n            scrubbed = scrubbed.replace(target, placeholder)\n",
    ),
)


def _setup_scratch(tmp_dir: Path) -> None:
    """Copy minimal source, scripts, and test trees into the scratch directory."""
    for rel_dir in ("src/mtp", "src/bsdbng", "methods"):
        src_dir = REPO_ROOT / rel_dir
        if src_dir.exists():
            shutil.copytree(src_dir, tmp_dir / rel_dir)
    for rel_file in (
        LOCK_REL,
        PROBE_REL,
        Path("scripts/create_trait_profile_infographic.py"),
        TEST_REL,
        Path("tests/test_trait_profile_infographic.py"),
    ):
        dst = tmp_dir / rel_file
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(REPO_ROOT / rel_file, dst)
    (tmp_dir / "data/derived/infographic_investigation").mkdir(parents=True, exist_ok=True)
    conftest = tmp_dir / "tests/conftest.py"
    expected_root = tmp_dir.resolve()
    conftest.write_text(
        "from pathlib import Path\n"
        "import mtp.paths as _paths\n\n"
        "def pytest_runtest_setup(item):\n"
        "    actual = _paths.REPO_ROOT.resolve()\n"
        f"    expected = Path({str(expected_root)!r}).resolve()\n"
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
                    "--tb=short",
                    "--rootdir",
                    str(tmp_dir),
                    str(tmp_dir / TEST_REL),
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
    with tempfile.TemporaryDirectory(prefix="mtp_1290_control_") as tmp:
        tmp_dir = Path(tmp)
        _setup_scratch(tmp_dir)
        rc, out = _run_case(tmp_dir)
        if rc != 0 or "PYTEST_ORIGIN_OK" not in out:
            print(f"FAIL [CONTROL]: rc={rc}\n{out}")
            return 1
        print("CONTROL: rc=0 PYTEST_ORIGIN_OK")

    for mutant in MUTANTS:
        with tempfile.TemporaryDirectory(prefix=f"mtp_1290_{mutant.name}_") as tmp:
            tmp_dir = Path(tmp)
            _setup_scratch(tmp_dir)
            target_path = tmp_dir / mutant.rel_path
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

    print(f"ALL {len(MUTANTS)} #1290/#1027 MUTANTS KILLED CLEANLY.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
