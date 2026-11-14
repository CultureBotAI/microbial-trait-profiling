"""Isolated scratch-copy mutant verifier for #1318 (Dependabot linkml + linkml-runtime group).

Purpose:
  Verify in an isolated temporary directory copy that tests/test_dependabot_config_1318.py
  passes on the current worktree files, runs from the scratch copy, and kills 8/8
  mutations across .github/dependabot.yml and uv.lock with AssertionError while leaving
  the worktree byte-identical.
Inputs:
  .github/dependabot.yml, uv.lock, and tests/test_dependabot_config_1318.py
Outputs:
  Stdout summary of control pass + 8/8 killed mutants (M01-M08).
Issue/PR:
  Refs #1318 (Closes #1318 Dependabot group; uv.lock already resolves pydantic 2.14.2 >= 2.13.5)
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
    ".github/dependabot.yml",
    "uv.lock",
    "tests/test_dependabot_config_1318.py",
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
        label="M01_drop_groups_block",
        rel_path=".github/dependabot.yml",
        target=(
            "    groups:\n"
            "      linkml:\n"
            "        patterns:\n"
            '          - "linkml"\n'
            '          - "linkml-runtime"\n'
        ),
        replacement="",
        expected_substring="Expected groups.linkml under the 'uv' package-ecosystem entry",
    ),
    MutantSpec(
        label="M02_drop_linkml_pattern",
        rel_path=".github/dependabot.yml",
        target='          - "linkml"\n          - "linkml-runtime"\n',
        replacement='          - "linkml-runtime"\n',
        expected_substring="Expected groups.linkml.patterns ==",
    ),
    MutantSpec(
        label="M03_drop_linkml_runtime_pattern",
        rel_path=".github/dependabot.yml",
        target='          - "linkml"\n          - "linkml-runtime"\n',
        replacement='          - "linkml"\n',
        expected_substring="Expected groups.linkml.patterns ==",
    ),
    MutantSpec(
        label="M04_move_groups_under_github_actions",
        rel_path=".github/dependabot.yml",
        target=(
            '  - package-ecosystem: "github-actions"\n'
            '    directory: "/"\n'
            "    schedule:\n"
            '      interval: "weekly"\n'
            '  # "uv" (not "pip") so each bump regenerates uv.lock as well as pyproject.toml;\n'
            "  # CI syncs with `uv sync --locked` and fails on a stale lock (#420).\n"
            '  - package-ecosystem: "uv"\n'
            '    directory: "/"\n'
            "    schedule:\n"
            '      interval: "weekly"\n'
            "    # Group linkml and linkml-runtime into a single PR (#1318, #1043).\n"
            "    groups:\n"
            "      linkml:\n"
            "        patterns:\n"
            '          - "linkml"\n'
            '          - "linkml-runtime"\n'
        ),
        replacement=(
            '  - package-ecosystem: "github-actions"\n'
            '    directory: "/"\n'
            "    schedule:\n"
            '      interval: "weekly"\n'
            "    groups:\n"
            "      linkml:\n"
            "        patterns:\n"
            '          - "linkml"\n'
            '          - "linkml-runtime"\n'
            '  - package-ecosystem: "uv"\n'
            '    directory: "/"\n'
            "    schedule:\n"
            '      interval: "weekly"\n'
        ),
        expected_substring="groups.linkml must be configured under package-ecosystem 'uv'",
    ),
    MutantSpec(
        label="M05_add_linkml_to_ignore",
        rel_path=".github/dependabot.yml",
        target='      - dependency-name: "uv"\n',
        replacement='      - dependency-name: "uv"\n      - dependency-name: "linkml"\n',
        expected_substring="must not appear in the uv ignore list",
    ),
    MutantSpec(
        label="M06_add_linkml_runtime_to_ignore",
        rel_path=".github/dependabot.yml",
        target='      - dependency-name: "detect-secrets"\n',
        replacement=(
            '      - dependency-name: "detect-secrets"\n      - dependency-name: "linkml-runtime"\n'
        ),
        expected_substring="must not appear in the uv ignore list",
    ),
    MutantSpec(
        label="M07_change_pattern_to_linkml_only",
        rel_path=".github/dependabot.yml",
        target='          - "linkml-runtime"\n',
        replacement='          - "linkml-only"\n',
        expected_substring="Expected groups.linkml.patterns ==",
    ),
    MutantSpec(
        label="M08_uv_lock_pydantic_below_floor",
        rel_path="uv.lock",
        target='name = "pydantic"\nversion = "2.14.2"\n',
        replacement='name = "pydantic"\nversion = "2.13.4"\n',
        expected_substring="Expected locked pydantic >=",
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
    test_file = scratch_root / "tests/test_dependabot_config_1318.py"
    plugin = _PytestPlugin(test_file)
    for mod_name in list(sys.modules):
        if "test_dependabot_config_1318" in mod_name:
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

    with tempfile.TemporaryDirectory(prefix="verify_1318_mutants_") as tmp_dir:
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
