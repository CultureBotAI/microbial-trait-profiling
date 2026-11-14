"""Isolated scratch-copy mutation verifier for Issue #1137 (2026-11-10).

Purpose:
  Verify that unit and documentation contract tests in ``tests/test_trait_key_policy.py``
  fail with ``AssertionError`` (pytest exit code 1) across 7 isolated mutants of
  ``src/mtp/trait_key_policy.py`` and ``methods/trait_key_policy.md`` while leaving the
  canonical working tree byte-identical before and after execution.

Inputs:
  - ``src/mtp/trait_key_policy.py``
  - ``methods/trait_key_policy.md``
  - ``tests/test_trait_key_policy.py``

Outputs:
  - JSON verification report printed to stdout (and optional ``--output-json``).
"""

from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import shutil
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
TARGET_FILES = (
    "src/mtp/trait_key_policy.py",
    "src/bsdbng/trait_key_policy.py",
    "methods/trait_key_policy.md",
    "tests/test_trait_key_policy.py",
)


@dataclass(frozen=True)
class MutantSpec:
    """Single isolated text mutation against a relative repository path."""

    name: str
    rel_path: str
    target: str
    replacement: str
    test_node: str


MUTANTS: tuple[MutantSpec, ...] = (
    MutantSpec(
        name="bypass_version_table_in_canonical",
        rel_path="src/mtp/trait_key_policy.py",
        target="    return definition.canonical_trait_id if obj in definition.member_ids else obj",
        replacement="    return FACULTATIVE_TRAIT_ID if obj in FACULTATIVE_MEMBER_IDS else obj",
        test_node="test_helpers_read_version_table_rather_than_unversioned_globals",
    ),
    MutantSpec(
        name="bypass_version_table_in_members",
        rel_path="src/mtp/trait_key_policy.py",
        target=(
            "    if obj_id == definition.canonical_trait_id:\n        return definition.member_ids"
        ),
        replacement=(
            "    if obj_id == FACULTATIVE_TRAIT_ID:\n        return FACULTATIVE_MEMBER_IDS"
        ),
        test_node="test_helpers_read_version_table_rather_than_unversioned_globals",
    ),
    MutantSpec(
        name="bypass_version_table_in_display_label",
        rel_path="src/mtp/trait_key_policy.py",
        target=(
            "        return effective_label.replace(definition.raw_label, definition.display_label)"
        ),
        replacement=(
            "        return effective_label.replace(FACULTATIVE_RAW_LABEL, FACULTATIVE_UNION_LABEL)"
        ),
        test_node="test_helpers_read_version_table_rather_than_unversioned_globals",
    ),
    MutantSpec(
        name="bypass_version_table_in_policy_metadata",
        rel_path="src/mtp/trait_key_policy.py",
        target=(
            "        asdict(item) for item in definition.source_vocabularies "
            "if item.object_id in members"
        ),
        replacement=(
            "        asdict(item) for item in FACULTATIVE_SOURCE_VOCABULARIES "
            "if item.object_id in members"
        ),
        test_node="test_helpers_read_version_table_rather_than_unversioned_globals",
    ),
    MutantSpec(
        name="allow_in_place_v1_mutation",
        rel_path="src/mtp/trait_key_policy.py",
        target="        if candidate != frozen_def:",
        replacement="        if False and candidate != frozen_def:",
        test_node="test_validate_trait_key_policy_definitions_rejects_in_place_v1_mutation",
    ),
    MutantSpec(
        name="drop_union_member_from_v1_table",
        rel_path="src/mtp/trait_key_policy.py",
        target='        member_ids=("METPO:1000605", "METPO:1000608", "METPO:1000612"),',
        replacement='        member_ids=("METPO:1000605", "METPO:1000608"),',
        test_node="test_policy_definitions_table_is_frozen_and_complete",
    ),
    MutantSpec(
        name="omit_methods_version_freezing_contract",
        rel_path="methods/trait_key_policy.md",
        target=(
            "Any future change to member CURIEs, canonical\n"
            "mappings, labels, or source vocabularies requires introducing a new policy version\n"
            "identifier"
        ),
        replacement=(
            "Any future change to member CURIEs, canonical\n"
            "mappings, labels, or source vocabularies updates the existing definition in place"
        ),
        test_node="test_methods_trait_key_policy_documents_version_frozen_table",
    ),
)


class _FailureRecorder:
    """Capture call-phase exception types from in-process pytest runs."""

    def __init__(self) -> None:
        self.call_failures: list[tuple[str, str]] = []

    def pytest_runtest_makereport(self, item: Any, call: Any) -> None:
        if call.when == "call" and call.excinfo is not None:
            self.call_failures.append((item.nodeid, call.excinfo.typename))


def _sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _snapshot_digests(root: Path) -> dict[str, str]:
    return {rel: _sha256_file(root / rel) for rel in TARGET_FILES}


def _copy_isolated_tree(scratch_root: Path) -> None:
    shutil.copytree(REPO_ROOT / "src", scratch_root / "src")
    shutil.copytree(REPO_ROOT / "methods", scratch_root / "methods")
    (scratch_root / "tests").mkdir(parents=True, exist_ok=True)
    shutil.copy2(
        REPO_ROOT / "tests/test_trait_key_policy.py",
        scratch_root / "tests/test_trait_key_policy.py",
    )
    shutil.copy2(REPO_ROOT / "pyproject.toml", scratch_root / "pyproject.toml")


def _purge_mtp_modules() -> None:
    for mod_name in list(sys.modules):
        if mod_name in ("mtp", "bsdbng", "test_trait_key_policy") or mod_name.startswith(
            ("mtp.", "bsdbng.")
        ):
            sys.modules.pop(mod_name, None)
    importlib.invalidate_caches()


def _run_pytest_node(scratch_root: Path, test_node: str) -> tuple[int, list[tuple[str, str]]]:
    saved_path = list(sys.path)
    recorder = _FailureRecorder()
    try:
        sys.path.insert(0, str(scratch_root / "src"))
        _purge_mtp_modules()
        target = f"{scratch_root / 'tests/test_trait_key_policy.py'}::{test_node}"
        code = int(
            pytest.main(
                [
                    "-q",
                    "-o",
                    f"cache_dir={scratch_root / '.pytest_cache'}",
                    target,
                ],
                plugins=[recorder],
            )
        )
        return code, recorder.call_failures
    finally:
        _purge_mtp_modules()
        sys.path[:] = saved_path


def run_verification() -> dict[str, Any]:
    before_digests = _snapshot_digests(REPO_ROOT)
    results: list[dict[str, Any]] = []

    with tempfile.TemporaryDirectory(prefix="verify_1137_mutants_") as tmp_dir:
        base_dir = Path(tmp_dir)
        control_root = base_dir / "control"
        _copy_isolated_tree(control_root)
        for spec in MUTANTS:
            ctrl_code, ctrl_failures = _run_pytest_node(control_root, spec.test_node)
            if ctrl_code != 0 or ctrl_failures:
                raise RuntimeError(
                    f"Positive control failed for {spec.test_node}: exit={ctrl_code}, "
                    f"failures={ctrl_failures}"
                )

        for spec in MUTANTS:
            mutant_root = base_dir / spec.name
            _copy_isolated_tree(mutant_root)
            target_file = mutant_root / spec.rel_path
            original_text = target_file.read_text(encoding="utf-8")
            occurrences = original_text.count(spec.target)
            if occurrences != 1:
                raise RuntimeError(
                    f"Mutant {spec.name} target occurred {occurrences} times in {spec.rel_path}"
                )
            mutated_text = original_text.replace(spec.target, spec.replacement, 1)
            if mutated_text == original_text:
                raise RuntimeError(f"Mutant {spec.name} did not alter {spec.rel_path}")
            target_file.write_text(mutated_text, encoding="utf-8")

            exit_code, failures = _run_pytest_node(mutant_root, spec.test_node)
            killed_by_assertion = (
                exit_code == 1
                and bool(failures)
                and all(exc_type == "AssertionError" for _, exc_type in failures)
            )
            if not killed_by_assertion:
                raise RuntimeError(
                    f"Mutant {spec.name} was not killed by AssertionError: "
                    f"exit_code={exit_code}, failures={failures}"
                )
            results.append(
                {
                    "mutant": spec.name,
                    "rel_path": spec.rel_path,
                    "test_node": spec.test_node,
                    "exit_code": exit_code,
                    "failures": failures,
                    "killed_by_assertion": killed_by_assertion,
                }
            )

    after_digests = _snapshot_digests(REPO_ROOT)
    if before_digests != after_digests:
        raise RuntimeError("Canonical working tree digests changed during mutant verification")

    return {
        "issue": 1137,
        "mutant_count": len(results),
        "all_killed_by_assertion": all(item["killed_by_assertion"] for item in results),
        "canonical_digests_unchanged": True,
        "canonical_sha256": after_digests,
        "mutants": results,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-json", type=Path, default=None)
    args = parser.parse_args()
    summary = run_verification()
    rendered = json.dumps(summary, indent=2, sort_keys=True)
    if args.output_json is not None:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(rendered + "\n", encoding="utf-8")
    print(rendered)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
