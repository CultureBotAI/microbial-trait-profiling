"""Isolated scratch-copy mutation verifier for Issue #1261 (2026-11-10).

Purpose:
  Verify that tests in ``tests/test_historical_trait_policy_core.py`` fail with
  ``AssertionError`` (pytest exit code 1) across 6 isolated mutants of
  ``tests/fixtures/trait_key_policy/README.md``, ``methods/three_test_table.md``,
  ``methods/delta_z_derived.md``, and ``tests/test_historical_trait_policy_core.py``
  while leaving the canonical working tree byte-identical before and after execution.

Inputs:
  - ``tests/fixtures/trait_key_policy/README.md``
  - ``tests/test_historical_trait_policy_core.py``
  - ``methods/three_test_table.md``
  - ``methods/delta_z_derived.md``

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
    "tests/fixtures/trait_key_policy/README.md",
    "tests/test_historical_trait_policy_core.py",
    "methods/three_test_table.md",
    "methods/delta_z_derived.md",
)
TEST_NODE = "test_fixture_readme_audit_command_and_methods_member_wording_are_accurate"


@dataclass(frozen=True)
class MutantSpec:
    """Single isolated text mutation against a relative repository path."""

    name: str
    rel_path: str
    target: str
    replacement: str


MUTANTS: tuple[MutantSpec, ...] = (
    MutantSpec(
        name="readme_stale_ingest_bugsigdb_glob",
        rel_path="tests/fixtures/trait_key_policy/README.md",
        target='     "$EXTRACTED_FIXTURE"/core/steps/ingest/studies/*.yaml \\',
        replacement='     "$EXTRACTED_FIXTURE"/core/steps/ingest_bugsigdb/studies/*.yaml \\',
    ),
    MutantSpec(
        name="readme_omit_export_tables_glob",
        rel_path="tests/fixtures/trait_key_policy/README.md",
        target='     "$EXTRACTED_FIXTURE"/core/steps/export_tables/derived/*.tsv \\\n',
        replacement="",
    ),
    MutantSpec(
        name="readme_omit_feature_catalog_glob",
        rel_path="tests/fixtures/trait_key_policy/README.md",
        target='     "$EXTRACTED_FIXTURE"/core/steps/feature_catalog/derived/*.tsv \\\n',
        replacement="",
    ),
    MutantSpec(
        name="three_test_table_overstate_member_read",
        rel_path="methods/three_test_table.md",
        target=(
            "nor do `three_test_table` or `delta_z_derived`\n"
            "parse `trait_member_ids_json` during statistical computation, though downstream\n"
            "consumers can read `trait_member_ids_json` from the core's `feature_metadata`\n"
            "when resolving policy membership"
        ),
        replacement=(
            "and instead read `trait_member_ids_json` from the core's `feature_metadata`\n"
            "when resolving policy membership"
        ),
    ),
    MutantSpec(
        name="delta_z_derived_overstate_member_read",
        rel_path="methods/delta_z_derived.md",
        target=(
            "nor do `three_test_table` or `delta_z_derived` parse `trait_member_ids_json`\n"
            "during statistical computation, though downstream consumers can read\n"
            "`trait_member_ids_json` from the core's `feature_metadata` when resolving policy\n"
            "membership"
        ),
        replacement=(
            "and instead read `trait_member_ids_json` from the core's `feature_metadata` when\n"
            "resolving policy membership"
        ),
    ),
    MutantSpec(
        name="helper_sha256_accept_non_hex",
        rel_path="tests/test_historical_trait_policy_core.py",
        target=(
            '    return isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value) is not None'
        ),
        replacement="    return isinstance(value, str) and len(value) == 64",
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
    shutil.copytree(
        REPO_ROOT / "tests/fixtures/trait_key_policy",
        scratch_root / "tests/fixtures/trait_key_policy",
    )
    shutil.copy2(
        REPO_ROOT / "tests/test_historical_trait_policy_core.py",
        scratch_root / "tests/test_historical_trait_policy_core.py",
    )
    shutil.copy2(REPO_ROOT / "pyproject.toml", scratch_root / "pyproject.toml")


def _purge_modules() -> None:
    for mod_name in list(sys.modules):
        if mod_name in (
            "mtp",
            "bsdbng",
            "test_historical_trait_policy_core",
        ) or mod_name.startswith(("mtp.", "bsdbng.")):
            sys.modules.pop(mod_name, None)
    importlib.invalidate_caches()


def _run_pytest_node(scratch_root: Path, test_node: str) -> tuple[int, list[tuple[str, str]]]:
    saved_path = list(sys.path)
    recorder = _FailureRecorder()
    try:
        sys.path.insert(0, str(scratch_root / "src"))
        _purge_modules()
        target = f"{scratch_root / 'tests/test_historical_trait_policy_core.py'}::{test_node}"
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
        _purge_modules()
        sys.path[:] = saved_path


def run_verification() -> dict[str, Any]:
    before_digests = _snapshot_digests(REPO_ROOT)
    results: list[dict[str, Any]] = []

    with tempfile.TemporaryDirectory(prefix="verify_1261_mutants_") as tmp_dir:
        base_dir = Path(tmp_dir)
        control_root = base_dir / "control"
        _copy_isolated_tree(control_root)
        ctrl_code, ctrl_failures = _run_pytest_node(control_root, TEST_NODE)
        if ctrl_code != 0 or ctrl_failures:
            raise RuntimeError(
                f"Positive control failed for {TEST_NODE}: exit={ctrl_code}, "
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

            exit_code, failures = _run_pytest_node(mutant_root, TEST_NODE)
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
                    "test_node": TEST_NODE,
                    "exit_code": exit_code,
                    "failures": failures,
                    "killed_by_assertion": killed_by_assertion,
                }
            )

    after_digests = _snapshot_digests(REPO_ROOT)
    if before_digests != after_digests:
        raise RuntimeError("Canonical working tree digests changed during mutant verification")

    return {
        "issue": 1261,
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
