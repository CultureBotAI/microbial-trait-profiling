"""Isolated mutation verifier for Issue #1023 (condition-module parser sibling import).

Purpose:
    Verify in an isolated temporary copy of the repository checkout that
    `tests/test_arm_direction_parser.py` binds to the isolated checkout and
    kills 6 targeted mutations in `scripts/create_condition_trait_module_heatmap.py`
    via `AssertionError`, without modifying any tracked file in the live worktree.
Inputs:
    `scripts/create_condition_trait_module_heatmap.py` and
    `tests/test_arm_direction_parser.py` in the current worktree.
Outputs:
    Per-control and per-mutant status summary printed to stdout; exit code 0
    when all controls pass and all 6 mutants die by `AssertionError`.
Issue / PR:
    Issue #1023
Date:
    2026-11-09
"""

from __future__ import annotations

import hashlib
import io
import os
import shutil
import sys
import tempfile
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from typing import Any

import pytest

import mtp.paths as mtp_paths
from mtp.workflow.arms import runtime as arm_runtime

REPO_ROOT = Path(__file__).resolve().parents[2]

TRACKED_TARGETS: tuple[str, ...] = (
    "scripts/create_condition_trait_module_heatmap.py",
    "tests/test_arm_direction_parser.py",
)

MUTANTS: tuple[tuple[str, str, str, str], ...] = (
    (
        "M01_move_sibling_import_to_module_top_level",
        "from mtp.human_host_corpus import CORPUS_SEQUENCING_TYPES\n",
        "from mtp.human_host_corpus import CORPUS_SEQUENCING_TYPES\n"
        "from create_study_trait_profile_heatmap import load_raw_experiment_metadata\n",
        "test_condition_module_parser_loads_without_sibling_import_regardless_of_order",
    ),
    (
        "M02_prepend_scripts_to_syspath_and_import_at_module_top_level",
        "from mtp.human_host_corpus import CORPUS_SEQUENCING_TYPES\n",
        "from mtp.human_host_corpus import CORPUS_SEQUENCING_TYPES\n"
        "import sys as _sys\n"
        "_sys.path.insert(0, str(Path(__file__).resolve().parent))\n"
        "from create_study_trait_profile_heatmap import load_raw_experiment_metadata\n",
        "test_condition_module_parser_loads_without_sibling_import_regardless_of_order",
    ),
    (
        "M03_import_sibling_inside_build_parser",
        "    return parser\n\n\ndef main() -> None:\n",
        "    from create_study_trait_profile_heatmap import load_raw_experiment_metadata\n"
        "    _ = load_raw_experiment_metadata\n"
        "    return parser\n\n\ndef main() -> None:\n",
        "test_condition_module_uses_its_execution_parser",
    ),
    (
        "M04_import_sibling_before_include_media_traits_error_in_main",
        "    if args.include_media_traits:\n"
        '        parser.error("MediaDive is excluded before composition normalization (#174).")\n'
        "\n"
        "    # This sibling-script import belongs to execution. The shared parser must be\n"
        "    # loadable without a script directory on sys.path or any scientific data access.\n"
        "    from create_study_trait_profile_heatmap import load_raw_experiment_metadata\n",
        "    from create_study_trait_profile_heatmap import load_raw_experiment_metadata\n"
        "    if args.include_media_traits:\n"
        '        parser.error("MediaDive is excluded before composition normalization (#174).")\n',
        "test_condition_module_parser_loads_without_sibling_import_regardless_of_order",
    ),
    (
        "M05_omit_load_raw_experiment_metadata_call_in_main",
        "    experiment_metadata = load_raw_experiment_metadata(args.full_dump_input)\n",
        '    raise RuntimeError("skipped deferred sibling call")\n',
        "test_condition_module_parser_loads_without_sibling_import_regardless_of_order",
    ),
    (
        "M06_pass_wrong_input_path_to_load_raw_experiment_metadata",
        "    experiment_metadata = load_raw_experiment_metadata(args.full_dump_input)\n",
        "    experiment_metadata = load_raw_experiment_metadata(args.studies_input)\n",
        "test_condition_module_parser_loads_without_sibling_import_regardless_of_order",
    ),
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _stage_isolated_checkout(root: Path) -> Path:
    isolated = (root / "repo").resolve()
    shutil.copytree(REPO_ROOT / "scripts", isolated / "scripts")
    shutil.copy2(REPO_ROOT / "pyproject.toml", isolated / "pyproject.toml")
    (isolated / "tests").mkdir(parents=True, exist_ok=True)
    shutil.copy2(
        REPO_ROOT / "tests/test_arm_direction_parser.py",
        isolated / "tests/test_arm_direction_parser.py",
    )
    return isolated


class _OriginCheckPlugin:
    def __init__(self, expected_test_file: Path, expected_root: Path) -> None:
        self.expected_test_file = expected_test_file.resolve()
        self.expected_root = expected_root.resolve()
        self.verified = False

    def pytest_collection_modifyitems(self, items: list[pytest.Item]) -> None:
        for item in items:
            if isinstance(item, pytest.Function):
                mod_file = getattr(item.module, "__file__", None)
                if mod_file is None or Path(mod_file).resolve() != self.expected_test_file:
                    raise AssertionError(f"Unexpected test module origin: {mod_file}")
                if getattr(item.module, "REPO_ROOT", None) != self.expected_root:
                    raise AssertionError("Test module REPO_ROOT does not point to isolated root")
                self.verified = True


def _run_pytest_inprocess(isolated: Path, expr: str) -> tuple[int, str]:
    os.environ["GIT_CONFIG_GLOBAL"] = "/dev/null"
    os.environ["GIT_CONFIG_NOSYSTEM"] = "1"
    runtime_mod: Any = arm_runtime
    old_paths_root = mtp_paths.REPO_ROOT
    old_runtime_root = runtime_mod.REPO_ROOT
    test_file = (isolated / "tests/test_arm_direction_parser.py").resolve()
    plugin = _OriginCheckPlugin(test_file, isolated)
    buf = io.StringIO()
    try:
        mtp_paths.REPO_ROOT = isolated
        runtime_mod.REPO_ROOT = isolated
        sys.modules.pop("test_arm_direction_parser", None)
        sys.modules.pop("create_study_trait_profile_heatmap", None)
        with redirect_stdout(buf), redirect_stderr(buf):
            code = int(
                pytest.main(
                    [
                        "-o",
                        "addopts=",
                        str(test_file),
                        "-k",
                        expr,
                        "-q",
                        "--tb=short",
                    ],
                    plugins=[plugin],
                )
            )
    finally:
        sys.modules.pop("test_arm_direction_parser", None)
        sys.modules.pop("create_study_trait_profile_heatmap", None)
        mtp_paths.REPO_ROOT = old_paths_root
        runtime_mod.REPO_ROOT = old_runtime_root
    if not plugin.verified:
        raise AssertionError("Pytest collection did not verify isolated test module origin")
    return code, buf.getvalue()


def main() -> int:
    if (mtp_paths.REPO_ROOT / "src/mtp").resolve() != (REPO_ROOT / "src/mtp").resolve():
        raise AssertionError("mtp imported from a different checkout than REPO_ROOT")
    before_hashes = {rel: _sha256(REPO_ROOT / rel) for rel in TRACKED_TARGETS}
    with tempfile.TemporaryDirectory(prefix="verify_1023_mutants_") as tmp_str:
        isolated = _stage_isolated_checkout(Path(tmp_str))

        ctrl_rc, ctrl_out = _run_pytest_inprocess(isolated, "condition_module or figure_parser")
        if ctrl_rc != 0:
            print("Control failed:\n", ctrl_out)
            return 1
        print(f"Controls + isolated origin check passed in {isolated}")

        target_file = isolated / "scripts/create_condition_trait_module_heatmap.py"
        original_text = target_file.read_text(encoding="utf-8")
        killed = 0
        for name, target, replacement, test_expr in MUTANTS:
            if original_text.count(target) != 1:
                print(f"[ERROR] {name}: target count {original_text.count(target)} != 1")
                return 1
            try:
                target_file.write_text(
                    original_text.replace(target, replacement, 1), encoding="utf-8"
                )
                rc, output = _run_pytest_inprocess(isolated, test_expr)
                if rc != 0 and "AssertionError" in output:
                    killed += 1
                    print(f"  [KILLED] {name} (rc={rc}, AssertionError)")
                else:
                    print(f"  [SURVIVED/UNEXPECTED] {name} (rc={rc})\n{output}")
                    return 1
            finally:
                target_file.write_text(original_text, encoding="utf-8")

    after_hashes = {rel: _sha256(REPO_ROOT / rel) for rel in TRACKED_TARGETS}
    if before_hashes != after_hashes:
        print("[ERROR] Live worktree hashes changed!")
        return 1
    print(f"Summary: {killed}/{len(MUTANTS)} mutants killed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
