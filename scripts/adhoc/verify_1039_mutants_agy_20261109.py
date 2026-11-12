"""Isolated scratch-copy mutant verifier for #1039 (2026-11-09).

Purpose:
    Prove that tests/test_mtp_imports.py passes on the clean worktree copy inside an
    isolated scratch directory (with scratch-origin check) and kills every single-site
    mutant of src/mtp/export_audit.py by AssertionError.
Inputs:
    Worktree files:
      - src/mtp/export_audit.py
      - tests/test_mtp_imports.py
Outputs:
    Stdout summary of control pass + mutant kill results; exits 0 iff all mutants are
    killed by AssertionError and worktree files remain untouched.
Issue:
    #1039
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
    "src/mtp/export_audit.py",
    "src/bsdbng/export_audit.py",
    "scripts/adhoc/issue937_audit_exports.py",
    "tests/test_mtp_imports.py",
)

_PYTEST_K_EXPR = (
    "test_export_audit_module_and_legacy_entrypoints_agree or "
    "test_declared_exports_rejects_all_mutations_and_non_top_level_bindings or "
    "test_public_definitions_discovers_nested_module_bindings_and_ignores_guards or "
    "test_export_omissions_detects_nested_public_definitions_and_reexports or "
    "test_annotated_shim_declaration_rejects_malformed_values"
)

MUTANTS: tuple[tuple[str, str, str], ...] = (
    (
        "M01_allow_augassign_all",
        (
            "        elif isinstance(node, ast.AugAssign):\n"
            '            if "__all__" in bound_names(node.target) or '
            "_target_mutates_all(node.target):\n"
            "                return True"
        ),
        "        elif isinstance(node, ast.AugAssign):\n            pass",
    ),
    (
        "M02_allow_mutating_method_call_all",
        (
            '            _root_target_name(node.func.value) == "__all__"\n'
            "            and node.func.attr in _MUTATING_ALL_METHODS"
        ),
        ('            _root_target_name(node.func.value) == "__all__"\n            and False'),
    ),
    (
        "M03_allow_subscript_assign_all",
        (
            "            if any(_target_mutates_all(target) for target in targets):\n"
            "                return True"
        ),
        "            if False:\n                return True",
    ),
    (
        "M04_allow_delete_all",
        (
            "        elif isinstance(node, ast.Delete):\n"
            "            if any(\n"
            '                "__all__" in bound_names(target) or _target_mutates_all(target)\n'
            "                for target in node.targets\n"
            "            ):\n"
            "                return True"
        ),
        "        elif isinstance(node, ast.Delete):\n            pass",
    ),
    (
        "M05_allow_nested_all_binding",
        (
            "            if id(node) not in valid_declaration_ids and any(\n"
            '                "__all__" in bound_names(target) for target in targets\n'
            "            ):\n"
            "                return True"
        ),
        "            if False:\n                return True",
    ),
    (
        "M06_skip_nested_if_statements_in_public_definitions",
        (
            "                else:\n"
            "                    yield from _iter_runtime_module_statements(stmt.body)\n"
            "                    yield from _iter_runtime_module_statements(stmt.orelse)"
        ),
        "                else:\n                    pass",
    ),
    (
        "M07_include_type_checking_guard_as_runtime",
        (
            "            if _is_type_checking_expr(stmt.test):\n"
            "                yield from _iter_runtime_module_statements(stmt.orelse)"
        ),
        (
            "            if _is_type_checking_expr(stmt.test):\n"
            "                yield from _iter_runtime_module_statements(stmt.body)"
        ),
    ),
    (
        "M08_ignore_nested_reexports_in_export_omissions",
        (
            "        # Canonical re-exports count; imported third-party libraries do not.\n"
            "        for node in _iter_runtime_module_statements(tree.body):"
        ),
        (
            "        # Canonical re-exports count; imported third-party libraries do not.\n"
            "        for node in tree.body:"
        ),
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
        loaded = sys.modules.get("mtp.export_audit")
        if loaded is None or loaded.__file__ is None:
            raise RuntimeError("mtp.export_audit not loaded before test execution")
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
                raw_prefix = reprcrash.message.split(":", 1)[0].strip()
                exc_name = (
                    "AssertionError"
                    if raw_prefix == "Failed" and "AssertionError" in reprcrash.message
                    else raw_prefix
                )
            self.failed_exc_names.append(exc_name)


def _snapshot_hashes(root: Path) -> dict[str, str]:
    return {rel: hashlib.sha256((root / rel).read_bytes()).hexdigest() for rel in REL_TARGETS}


def _copy_tree_subset(src_root: Path, dst_root: Path) -> None:
    for rel in REL_TARGETS:
        src = src_root / rel
        dst = dst_root / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


def _load_scratch_export_audit(mod_path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location("mtp.export_audit", mod_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not create ModuleSpec for {mod_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules["mtp.export_audit"] = module
    sys.modules["bsdbng.export_audit"] = module
    spec.loader.exec_module(module)
    mtp.__dict__["export_audit"] = module
    bsdbng.__dict__["export_audit"] = module
    return module


def _run_scratch_pytest(scratch: Path) -> tuple[int, _OutcomeCapturePlugin]:
    mod_path = scratch / "src/mtp/export_audit.py"
    test_path = scratch / "tests/test_mtp_imports.py"
    _load_scratch_export_audit(mod_path)
    for mod_name in list(sys.modules):
        if mod_name == "test_mtp_imports" or mod_name.endswith(".test_mtp_imports"):
            sys.modules.pop(mod_name, None)
    plugin = _OutcomeCapturePlugin(test_path, mod_path)
    rc = int(
        pytest.main(
            [
                "-q",
                "-o",
                "addopts=--import-mode=importlib",
                "-k",
                _PYTEST_K_EXPR,
                str(test_path),
            ],
            plugins=[plugin],
        )
    )
    return rc, plugin


def main() -> int:
    before_hashes = _snapshot_hashes(REPO_ROOT)
    original_module = importlib.import_module("mtp.export_audit")
    try:
        with tempfile.TemporaryDirectory(prefix="verify_1039_mutants_") as tmp_dir:
            scratch = Path(tmp_dir)
            _copy_tree_subset(REPO_ROOT, scratch)

            rc_ctrl, plugin_ctrl = _run_scratch_pytest(scratch)
            _require(
                rc_ctrl == 0 and plugin_ctrl.origin_checked and plugin_ctrl.passed_count >= 30,
                f"Control run failed: rc={rc_ctrl}, passed={plugin_ctrl.passed_count}",
            )
            print(f"CONTROL: rc={rc_ctrl} PYTEST_ORIGIN_OK passed={plugin_ctrl.passed_count}")

            base_file = scratch / "src/mtp/export_audit.py"
            original_text = base_file.read_text(encoding="utf-8")

            for mutant_id, old_text, new_text in MUTANTS:
                _require(
                    old_text in original_text,
                    f"{mutant_id}: target snippet not found in src/mtp/export_audit.py",
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
        sys.modules["mtp.export_audit"] = original_module
        sys.modules["bsdbng.export_audit"] = original_module
        mtp.__dict__["export_audit"] = original_module
        bsdbng.__dict__["export_audit"] = original_module

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
