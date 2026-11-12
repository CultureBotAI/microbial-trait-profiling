"""Isolated scratch-copy mutant verifier for #486 (drop no-op bandit[toml] extra).

Purpose:
  Verify in an isolated temporary directory copy that tests/test_bandit_config_486.py
  passes on the current worktree files, runs from the scratch copy, and kills 7/7
  mutations across pyproject.toml, uv.lock, .pre-commit-config.yaml, and the validator
  itself with AssertionError while leaving the worktree byte-identical.
Inputs:
  pyproject.toml, .pre-commit-config.yaml, uv.lock, and tests/test_bandit_config_486.py
Outputs:
  Stdout summary of control pass + 7/7 killed mutants (M01-M07).
Issue/PR:
  Closes #486
Date:
  2026-11-09
"""

from __future__ import annotations

import hashlib
import shutil
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
WATCHED_FILES: tuple[str, ...] = (
    "pyproject.toml",
    ".pre-commit-config.yaml",
    "uv.lock",
    "tests/test_bandit_config_486.py",
)


@dataclass(frozen=True)
class MutantSpec:
    """Single file-mutation specification."""

    label: str
    rel_path: str
    target: str
    replacement: str
    expected_substring: str


MUTANTS: tuple[MutantSpec, ...] = (
    MutantSpec(
        label="M01_pyproject_restore_toml_extra",
        rel_path="pyproject.toml",
        target='  "bandit>=1.9.4",\n',
        replacement='  "bandit[toml]>=1.9.4",\n',
        expected_substring="must not specify extras on bandit",
    ),
    MutantSpec(
        label="M02_uv_lock_restore_toml_extra",
        rel_path="uv.lock",
        target='    { name = "bandit", specifier = ">=1.9.4" },\n',
        replacement='    { name = "bandit", extras = ["toml"], specifier = ">=1.9.4" },\n',
        expected_substring="must not retain bandit extras",
    ),
    MutantSpec(
        label="M03_pre_commit_restore_additional_deps",
        rel_path=".pre-commit-config.yaml",
        target="      - id: bandit\n        args: [-c, pyproject.toml]\n",
        replacement=(
            "      - id: bandit\n"
            "        args: [-c, pyproject.toml]\n"
            '        additional_dependencies: ["bandit[toml]"]\n'
        ),
        expected_substring="additional_dependencies",
    ),
    MutantSpec(
        label="M04_pre_commit_drop_config_args",
        rel_path=".pre-commit-config.yaml",
        target="      - id: bandit\n        args: [-c, pyproject.toml]\n",
        replacement="      - id: bandit\n        args: []\n",
        expected_substring="Bandit hook must pass",
    ),
    MutantSpec(
        label="M05_pre_commit_rev_drift",
        rel_path=".pre-commit-config.yaml",
        target="  - repo: https://github.com/PyCQA/bandit\n    rev: 1.9.4\n",
        replacement="  - repo: https://github.com/PyCQA/bandit\n    rev: 1.8.3\n",
        expected_substring="Bandit pre-commit rev",
    ),
    MutantSpec(
        label="M06_pyproject_version_drift",
        rel_path="pyproject.toml",
        target='  "bandit>=1.9.4",\n',
        replacement='  "bandit>=1.8.3",\n',
        expected_substring="does not match uv.lock version",
    ),
    MutantSpec(
        label="M07_validator_allow_brackets",
        rel_path="tests/test_bandit_config_486.py",
        target='    assert "[" not in bandit_spec, (\n',
        replacement="    assert bool(bandit_spec), (\n",
        expected_substring="Regex pattern did not match",
    ),
)


class _PytestPlugin:
    """In-process pytest plugin capturing call-phase outcomes and AssertionError details."""

    def __init__(self, expected_test_file: Path) -> None:
        self.expected_test_file = expected_test_file.resolve()
        self.origin_ok = False
        self.assertion_failures: list[str] = []
        self.other_failures: list[str] = []

    def pytest_runtest_makereport(
        self,
        item: pytest.Item,
        call: pytest.CallInfo[None],
    ) -> None:
        if call.when != "call":
            return
        module = getattr(item, "module", None)
        if module is None:
            raise RuntimeError("Missing test module on pytest item")
        module_file = Path(str(getattr(module, "__file__", ""))).resolve()
        if module_file == self.expected_test_file:
            self.origin_ok = True
        if call.excinfo is not None:
            exc = call.excinfo.value
            detail = f"{type(exc).__name__}: {exc}"
            if isinstance(exc, AssertionError):
                self.assertion_failures.append(detail)
            else:
                self.other_failures.append(detail)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _snapshot_hashes(root: Path) -> dict[str, str]:
    return {rel: _sha256(root / rel) for rel in WATCHED_FILES}


def _copy_watched(src_root: Path, dst_root: Path) -> None:
    for rel in WATCHED_FILES:
        src = src_root / rel
        dst = dst_root / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


def _run_in_scratch(scratch_root: Path) -> tuple[int, _PytestPlugin]:
    test_file = scratch_root / "tests/test_bandit_config_486.py"
    plugin = _PytestPlugin(test_file)
    for mod_name in list(sys.modules):
        if "test_bandit_config_486" in mod_name:
            sys.modules.pop(mod_name, None)
    rc = int(
        pytest.main(
            [
                "-q",
                "-o",
                "addopts=",
                "-p",
                "no:cacheprovider",
                str(test_file),
            ],
            plugins=[plugin],
        )
    )
    return rc, plugin


def main() -> int:
    canonical_before = _snapshot_hashes(REPO_ROOT)

    with tempfile.TemporaryDirectory(prefix="verify_486_mutants_") as tmp_dir:
        scratch_root = Path(tmp_dir)
        _copy_watched(REPO_ROOT, scratch_root)
        scratch_before = _snapshot_hashes(scratch_root)
        _require(
            scratch_before == canonical_before,
            "Scratch copy hashes do not match canonical worktree",
        )

        rc, control_plugin = _run_in_scratch(scratch_root)
        _require(
            rc == 0
            and control_plugin.origin_ok
            and not control_plugin.assertion_failures
            and not control_plugin.other_failures,
            f"Control run failed: rc={rc}, origin_ok={control_plugin.origin_ok}",
        )
        print("CONTROL: rc=0 PYTEST_ORIGIN_OK")

        killed = 0
        for spec in MUTANTS:
            target_path = scratch_root / spec.rel_path
            original_text = target_path.read_text(encoding="utf-8")
            _require(
                original_text.count(spec.target) == 1,
                f"Mutant {spec.label} target count != 1 in {spec.rel_path}",
            )
            mutated_text = original_text.replace(spec.target, spec.replacement, 1)
            target_path.write_text(mutated_text, encoding="utf-8")
            try:
                m_rc, m_plugin = _run_in_scratch(scratch_root)
            finally:
                target_path.write_text(original_text, encoding="utf-8")

            _require(
                _snapshot_hashes(scratch_root) == scratch_before,
                f"Scratch restore mismatch after {spec.label}",
            )
            _require(
                m_rc != 0
                and m_plugin.origin_ok
                and len(m_plugin.assertion_failures) >= 1
                and not m_plugin.other_failures
                and any(
                    spec.expected_substring in detail for detail in m_plugin.assertion_failures
                ),
                (
                    f"Mutant {spec.label} was not killed as expected: "
                    f"rc={m_rc}, assertions={m_plugin.assertion_failures}, "
                    f"other={m_plugin.other_failures}"
                ),
            )
            killed += 1
            print(f"KILLED {spec.label}: {m_plugin.assertion_failures[0]}")

    canonical_after = _snapshot_hashes(REPO_ROOT)
    _require(
        canonical_after == canonical_before,
        "Canonical worktree files changed during mutation verification",
    )
    print(f"ALL {killed}/{len(MUTANTS)} MUTANTS KILLED BY ASSERTIONERROR")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
