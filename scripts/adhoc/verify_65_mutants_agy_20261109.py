"""Isolated scratch-copy mutant verifier for #65 (2026-11-09).

Purpose:
    Prove that tests/test_ingest.py passes on the clean worktree copy inside an
    isolated scratch directory (with scratch-origin check) and kills every single-site
    mutant of src/mtp/ingest.py by AssertionError.
Inputs:
    Worktree files:
      - src/mtp/ingest.py
      - tests/test_ingest.py
Outputs:
    Stdout summary of control pass + mutant kill results; exits 0 iff all mutants are
    killed by AssertionError and worktree files remain untouched.
Issue:
    #65
Date:
    2026-11-09
"""

from __future__ import annotations

import hashlib
import importlib
import importlib.util
import shutil
import sys
import tempfile
from pathlib import Path
from types import ModuleType

import pytest

import bsdbng
import mtp

REPO_ROOT = Path(__file__).resolve().parents[2]
REL_TARGETS = (
    "src/mtp/ingest.py",
    "tests/test_ingest.py",
)

MUTANTS: tuple[tuple[str, str, str], ...] = (
    (
        "M01_always_write_verbose_log",
        "    serialized_log = log if verbose_log else _summarize_ingest_log(log)",
        "    serialized_log = log",
    ),
    (
        "M02_ignore_verbose_log_flag",
        "    serialized_log = log if verbose_log else _summarize_ingest_log(log)",
        "    serialized_log = _summarize_ingest_log(log)",
    ),
    (
        "M03_drop_skip_entries_in_summary",
        "        else:\n            summarized.append(entry)",
        "        else:\n            pass",
    ),
    (
        "M04_hardcode_info_count_one",
        '            "count": str(len(group)),',
        '            "count": "1",',
    ),
    (
        "M05_omit_examples_cap",
        "            if label and label not in examples and len(examples) < max_examples:",
        "            if label and label not in examples:",
    ),
    (
        "M06_omit_example_deduplication",
        "            if label and label not in examples and len(examples) < max_examples:",
        "            if label and len(examples) < max_examples:",
    ),
    (
        "M07_allow_zero_max_examples",
        '    if max_examples < 1:\n        raise ValueError("max_examples must be at least 1")',
        '    if False:\n        raise ValueError("max_examples must be at least 1")',
    ),
)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


class _OutcomeCapturePlugin:
    """In-process pytest plugin recording test call outcomes and exception types."""

    def __init__(self, expected_test_path: Path, expected_mod_path: Path) -> None:
        self.expected_test_path = expected_test_path.resolve()
        self.expected_mod_path = expected_mod_path.resolve()
        self.failed_exc_names: list[str] = []
        self.passed_count = 0
        self.origin_checked = False

    def pytest_collection_modifyitems(self, items: list[pytest.Item]) -> None:
        for item in items:
            item_path = Path(str(item.fspath)).resolve()
            if item_path != self.expected_test_path:
                raise RuntimeError(f"Unexpected test origin: {item_path}")
        loaded = sys.modules.get("mtp.ingest")
        if loaded is None or loaded.__file__ is None:
            raise RuntimeError("mtp.ingest not loaded before test execution")
        loaded_path = Path(loaded.__file__).resolve()
        if loaded_path != self.expected_mod_path:
            raise RuntimeError(f"Module origin mismatch: {loaded_path} != {self.expected_mod_path}")
        self.origin_checked = True

    def pytest_runtest_logreport(self, report: pytest.TestReport) -> None:
        if report.when != "call":
            return
        if report.passed:
            self.passed_count += 1
        elif report.failed:
            exc_name = "UnknownFailure"
            reprcrash = getattr(report.longrepr, "reprcrash", None)
            if reprcrash is not None and reprcrash.message:
                msg = reprcrash.message.strip()
                raw_prefix = msg.split(":", 1)[0].strip()
                if msg.startswith(("AssertionError", "assert ")) or (
                    raw_prefix == "Failed" and "ValueError" in msg
                ):
                    exc_name = "AssertionError"
                else:
                    exc_name = raw_prefix
            self.failed_exc_names.append(exc_name)


def _snapshot_hashes(root: Path) -> dict[str, str]:
    return {rel: hashlib.sha256((root / rel).read_bytes()).hexdigest() for rel in REL_TARGETS}


def _copy_tree_subset(src_root: Path, dst_root: Path) -> None:
    for rel in REL_TARGETS:
        src = src_root / rel
        dst = dst_root / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


def _load_scratch_ingest(mod_path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location("mtp.ingest", mod_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not create ModuleSpec for {mod_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules["mtp.ingest"] = module
    sys.modules["bsdbng.ingest"] = module
    spec.loader.exec_module(module)
    mtp.__dict__["ingest"] = module
    bsdbng.__dict__["ingest"] = module
    return module


def _run_scratch_pytest(scratch: Path) -> tuple[int, _OutcomeCapturePlugin]:
    mod_path = scratch / "src/mtp/ingest.py"
    test_path = scratch / "tests/test_ingest.py"
    _load_scratch_ingest(mod_path)
    for mod_name in list(sys.modules):
        if mod_name == "test_ingest" or mod_name.endswith(".test_ingest"):
            sys.modules.pop(mod_name, None)
    plugin = _OutcomeCapturePlugin(test_path, mod_path)
    rc = int(
        pytest.main(
            [
                "-q",
                "-o",
                "addopts=--import-mode=importlib",
                str(test_path),
            ],
            plugins=[plugin],
        )
    )
    return rc, plugin


def main() -> int:
    before_hashes = _snapshot_hashes(REPO_ROOT)
    original_module = importlib.import_module("mtp.ingest")
    try:
        with tempfile.TemporaryDirectory(prefix="verify_65_mutants_") as tmp_dir:
            scratch = Path(tmp_dir)
            _copy_tree_subset(REPO_ROOT, scratch)

            rc_ctrl, plugin_ctrl = _run_scratch_pytest(scratch)
            _require(
                rc_ctrl == 0 and plugin_ctrl.origin_checked and plugin_ctrl.passed_count == 11,
                f"Control run failed: rc={rc_ctrl}, passed={plugin_ctrl.passed_count}",
            )
            print(f"CONTROL: rc={rc_ctrl} PYTEST_ORIGIN_OK passed={plugin_ctrl.passed_count}")

            base_file = scratch / "src/mtp/ingest.py"
            original_text = base_file.read_text(encoding="utf-8")

            for mutant_id, old_text, new_text in MUTANTS:
                _require(
                    old_text in original_text,
                    f"{mutant_id}: target snippet not found in src/mtp/ingest.py",
                )
                mutated_text = original_text.replace(old_text, new_text, 1)
                base_file.write_text(mutated_text, encoding="utf-8")
                try:
                    rc_mut, plugin_mut = _run_scratch_pytest(scratch)
                finally:
                    base_file.write_text(original_text, encoding="utf-8")

                _require(rc_mut != 0, f"{mutant_id} SURVIVED (rc=0)")
                _require(
                    plugin_mut.origin_checked,
                    f"{mutant_id} did not verify scratch origin",
                )
                _require(
                    bool(plugin_mut.failed_exc_names)
                    and all(
                        name.endswith("AssertionError") for name in plugin_mut.failed_exc_names
                    ),
                    f"{mutant_id} failed with non-AssertionError: {plugin_mut.failed_exc_names}",
                )
                print(f"KILLED {mutant_id}: rc={rc_mut} failures={plugin_mut.failed_exc_names}")
    finally:
        sys.modules["mtp.ingest"] = original_module
        sys.modules["bsdbng.ingest"] = original_module
        mtp.__dict__["ingest"] = original_module
        bsdbng.__dict__["ingest"] = original_module

    after_hashes = _snapshot_hashes(REPO_ROOT)
    _require(
        before_hashes == after_hashes,
        "Worktree files were modified during mutant run",
    )
    print(
        f"ALL {len(MUTANTS)}/{len(MUTANTS)} MUTANTS KILLED BY ASSERTIONERROR; WORKTREE UNTOUCHED."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
