"""Isolated scratch-copy mutant verifier for #329 (trait_screen_manifest.json self-attestation).

Purpose:
  Verify in an isolated temporary directory copy that tests/test_trait_screen_manifest_329.py
  passes on the current worktree files, runs from the scratch copy, and kills 8/8
  mutations on methods/trait_screen_manifest.json with AssertionError while leaving
  the worktree byte-identical.
Inputs:
  methods/trait_screen_manifest.json and tests/test_trait_screen_manifest_329.py
Outputs:
  Stdout summary of control pass + 8/8 killed mutants (M01-M08).
Issue/PR:
  Closes #329
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
    "methods/trait_screen_manifest.json",
    "tests/test_trait_screen_manifest_329.py",
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
        label="M01_declared_before_first_run_false",
        rel_path="methods/trait_screen_manifest.json",
        target='  "declared_before_first_run": true,\n',
        replacement='  "declared_before_first_run": false,\n',
        expected_substring="Expected declared_before_first_run to remain boolean True",
    ),
    MutantSpec(
        label="M02_declared_by_omit_self_attestation",
        rel_path="methods/trait_screen_manifest.json",
        target=(
            '  "declared_by": "drafted by a coding agent (Claude Code) before any screen result '
            "was computed; not revised after seeing results (self-attestation only: cannot be "
            "verified from Git history because methods/trait_screen_manifest.json landed in "
            "commit 1f3d36c alongside methods/trait_direction_screen.md, which already quoted "
            'Li & Ji effective family size 112.0; #329)",\n'
        ),
        replacement=(
            '  "declared_by": "drafted by a coding agent (Claude Code) before any screen result '
            'was computed; not revised after seeing results",\n'
        ),
        expected_substring=(
            "Expected declared_by to state that pre-run declaration rests on self-attestation"
        ),
    ),
    MutantSpec(
        label="M03_drop_verification_field",
        rel_path="methods/trait_screen_manifest.json",
        target=(
            '  "declared_before_first_run_verification": "self-attestation only (#329): pre-run '
            "timing cannot be verified from Git history alone because "
            "methods/trait_screen_manifest.json was first committed in 1f3d36c together with "
            "methods/trait_direction_screen.md, which already quoted screen results (including "
            "Li & Ji effective family size 112.0). Future pre-registration manifests must be "
            "committed on their own before any code or document that produces or quotes "
            'results.",\n'
        ),
        replacement="",
        expected_substring="Expected declared_before_first_run_verification string",
    ),
    MutantSpec(
        label="M04_verification_omit_cannot_be_verified",
        rel_path="methods/trait_screen_manifest.json",
        target="pre-run timing cannot be verified from Git history alone",
        replacement="pre-run timing is verified from Git history alone",
        expected_substring="cannot be verified from Git history",
    ),
    MutantSpec(
        label="M05_verification_omit_1f3d36c",
        rel_path="methods/trait_screen_manifest.json",
        target="first committed in 1f3d36c together with",
        replacement="first committed together with",
        expected_substring="Expected declared_before_first_run_verification to cite commit 1f3d36c",
    ),
    MutantSpec(
        label="M06_verification_omit_112_0",
        rel_path="methods/trait_screen_manifest.json",
        target="(including Li & Ji effective family size 112.0)",
        replacement="(including effective family size)",
        expected_substring="Expected declared_before_first_run_verification to cite Li & Ji",
    ),
    MutantSpec(
        label="M07_verification_omit_future_pre_registration_rule",
        rel_path="methods/trait_screen_manifest.json",
        target="Future pre-registration manifests must be committed on their own before any code",
        replacement="Future pre-registration manifests may be committed alongside any code",
        expected_substring="committed on their own before",
    ),
    MutantSpec(
        label="M08_amendment_omit_329",
        rel_path="methods/trait_screen_manifest.json",
        target=(
            "    {\n"
            '      "date": "2026-11-09",\n'
            '      "change": "declared_by and declared_before_first_run_verification disclose '
            "that declared_before_first_run rests on self-attestation only and cannot be "
            "verified from Git history alone, and record the standalone pre-commit rule for "
            'future pre-registration manifests.",\n'
            '      "reason": "#329: methods/trait_screen_manifest.json landed in the same '
            "commit (1f3d36c) as methods/trait_direction_screen.md, which already quoted "
            "screen results (Li & Ji 112.0). Future pre-registration manifests must be "
            "committed on their own before the code or documents that produce or quote "
            'results.",\n'
            '      "known_when_amended": "the full screen result and git commit history of '
            '1f3d36c.",\n'
            '      "unchanged": "declared_before_first_run, the candidate rule, the test, '
            'both corrections, both family sizes and every output column"\n'
            "    }\n"
        ),
        replacement=(
            "    {\n"
            '      "date": "2026-11-09",\n'
            '      "change": "other note",\n'
            '      "reason": "other",\n'
            '      "known_when_amended": "none",\n'
            '      "unchanged": "all"\n'
            "    }\n"
        ),
        expected_substring="Expected exactly one #329 amendment entry",
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
    test_file = scratch_root / "tests/test_trait_screen_manifest_329.py"
    plugin = _PytestPlugin(test_file)
    for mod_name in list(sys.modules):
        if "test_trait_screen_manifest_329" in mod_name:
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

    with tempfile.TemporaryDirectory(prefix="verify_329_mutants_") as tmp_dir:
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
