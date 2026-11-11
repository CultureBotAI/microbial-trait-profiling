"""Isolated scratch-copy mutation verifier for issue #1250.

Purpose:
    Verify that regressions in display_trait_label resolution and recorded-table
    KG label preservation across src/mtp/ontology_audit.py and
    src/mtp/coupling.py are caught by call-phase AssertionError in
    tests/test_ontology_audit.py and tests/test_trait_coupling.py without
    mutating the canonical worktree.
Inputs:
    Canonical repository files under src/mtp/ and tests/.
Outputs:
    JSON summary printed to stdout; exit code 0 iff all 9 mutants fail with
    call-phase AssertionError, all controls pass, and canonical file digests
    remain unchanged.
Issue:
    #1250
Date:
    2026-11-10
"""

from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import shutil
import sys
import tempfile
from collections.abc import Sequence
from contextlib import ExitStack
from dataclasses import dataclass
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
WATCHED_FILES = (
    "src/mtp/ontology_audit.py",
    "src/mtp/coupling.py",
    "tests/test_ontology_audit.py",
    "tests/test_trait_coupling.py",
)


@dataclass(frozen=True)
class MutantSpec:
    """Single source-text mutation against a target relative path."""

    name: str
    rel_path: str
    target: str
    replacement: str
    test_file: str
    test_node: str


MUTANTS: tuple[MutantSpec, ...] = (
    MutantSpec(
        name="ontology_audit_ignore_member_ids",
        rel_path="src/mtp/ontology_audit.py",
        target="        members = read_recorded_member_ids(columns_path)",
        replacement="        members: dict[str, tuple[str, ...]] = {}",
        test_file="tests/test_ontology_audit.py",
        test_node="test_build_redundant_pairs_table_applies_display_trait_label",
    ),
    MutantSpec(
        name="ontology_audit_skip_display_trait_label",
        rel_path="src/mtp/ontology_audit.py",
        target=(
            "            str(column_id): display_trait_label(\n"
            "                str(column_id),\n"
            "                str(column_label),\n"
            "                members.get(str(column_id)),\n"
            "            )"
        ),
        replacement="            str(column_id): str(column_label)",
        test_file="tests/test_ontology_audit.py",
        test_node="test_build_redundant_pairs_table_applies_display_trait_label",
    ),
    MutantSpec(
        name="ontology_audit_leak_display_label_when_disabled",
        rel_path="src/mtp/ontology_audit.py",
        target="    if use_display_labels:",
        replacement="    if True:",
        test_file="tests/test_ontology_audit.py",
        test_node="test_build_redundant_pairs_table_applies_display_trait_label",
    ),
    MutantSpec(
        name="coupling_resolve_ignore_columns_path",
        rel_path="src/mtp/coupling.py",
        target="if columns_path is None:",
        replacement="if True:",
        test_file="tests/test_trait_coupling.py",
        test_node="test_coupling_applies_display_trait_label_for_union_and_preserves_raw",
    ),
    MutantSpec(
        name="coupling_resolve_skip_display_trait_label",
        rel_path="src/mtp/coupling.py",
        target="display_trait_label(trait_id, raw_name, members.get(trait_id))",
        replacement="raw_name",
        test_file="tests/test_trait_coupling.py",
        test_node="test_coupling_applies_display_trait_label_for_union_and_preserves_raw",
    ),
    MutantSpec(
        name="coupling_composition_arm_use_raw_table1_traits",
        rel_path="src/mtp/coupling.py",
        target=(
            "            trait_id=t,\n"
            "            trait_name=resolved_names.get(t, t),\n"
            '            primary_stat=stats_list[idx]["effect"],'
        ),
        replacement=(
            "            trait_id=t,\n"
            "            trait_name=TABLE_1_TRAITS.get(t, t),\n"
            '            primary_stat=stats_list[idx]["effect"],'
        ),
        test_file="tests/test_trait_coupling.py",
        test_node="test_coupling_applies_display_trait_label_for_union_and_preserves_raw",
    ),
    MutantSpec(
        name="coupling_loto_dropped_name_use_raw_table1_traits",
        rel_path="src/mtp/coupling.py",
        target="dropped_name = resolved_names.get(dropped, dropped)",
        replacement="dropped_name = TABLE_1_TRAITS.get(dropped, dropped)",
        test_file="tests/test_trait_coupling.py",
        test_node="test_coupling_applies_display_trait_label_for_union_and_preserves_raw",
    ),
    MutantSpec(
        name="coupling_report_drop_columns_path_resolution",
        rel_path="src/mtp/coupling.py",
        target="        corpus_metrics=corpus_metrics,\n        columns_path=columns_path,",
        replacement="        corpus_metrics=corpus_metrics,\n        columns_path=None,",
        test_file="tests/test_trait_coupling.py",
        test_node="test_coupling_applies_display_trait_label_for_union_and_preserves_raw",
    ),
    MutantSpec(
        name="coupling_recorded_lco_leak_display_names",
        rel_path="src/mtp/coupling.py",
        target="        focal_traits=t1_keys,\n        cluster_to_drop=CANONICAL_LCO_CLUSTER,",
        replacement=(
            "        focal_traits=t1_keys,\n"
            "        cluster_to_drop=CANONICAL_LCO_CLUSTER,\n"
            "        columns_path=columns_path,"
        ),
        test_file="tests/test_trait_coupling.py",
        test_node="test_coupling_applies_display_trait_label_for_union_and_preserves_raw",
    ),
)


class _AssertionCapturePlugin:
    """Capture failed test reports, phases, and exception types from in-process pytest."""

    def __init__(self) -> None:
        self.failures: list[tuple[str, str, str]] = []

    def pytest_runtest_logreport(self, report: pytest.TestReport) -> None:
        if report.failed:
            exc_name = "UnknownFailure"
            if report.longrepr is not None and hasattr(report.longrepr, "reprcrash"):
                reprcrash = getattr(report.longrepr, "reprcrash", None)
                if reprcrash is not None and hasattr(reprcrash, "message"):
                    msg = str(reprcrash.message).strip()
                    if msg.startswith(("AssertionError", "assert ")):
                        exc_name = "AssertionError"
                    else:
                        exc_name = msg.split(":", 1)[0].strip()
            self.failures.append((report.nodeid, report.when, exc_name))

    def pytest_collectreport(self, report: pytest.CollectReport) -> None:
        if report.failed:
            self.failures.append((str(report.nodeid), "collect", "CollectFailure"))


def _sha256_map(root: Path) -> dict[str, str]:
    digests: dict[str, str] = {}
    for rel in WATCHED_FILES:
        digests[rel] = hashlib.sha256((root / rel).read_bytes()).hexdigest()
    return digests


def _prepare_scratch_tree(scratch_root: Path) -> None:
    shutil.copytree(REPO_ROOT / "src" / "mtp", scratch_root / "src" / "mtp")
    tests_dir = scratch_root / "tests"
    tests_dir.mkdir(parents=True, exist_ok=True)
    for test_rel in ("test_ontology_audit.py", "test_trait_coupling.py"):
        shutil.copy2(REPO_ROOT / "tests" / test_rel, tests_dir / test_rel)


def _purge_loaded_modules() -> None:
    for mod_name in list(sys.modules):
        if (
            mod_name == "mtp"
            or mod_name.startswith("mtp.")
            or mod_name in {"test_ontology_audit", "test_trait_coupling"}
        ):
            sys.modules.pop(mod_name, None)
    importlib.invalidate_caches()


def _run_pytest_in_scratch(
    scratch_root: Path, test_file: str, test_node: str
) -> tuple[int, list[tuple[str, str, str]]]:
    plugin = _AssertionCapturePlugin()
    old_sys_path = list(sys.path)
    sys.path.insert(0, str(scratch_root / "src"))
    _purge_loaded_modules()
    try:
        test_target = f"{scratch_root / test_file}::{test_node}"
        code = int(
            pytest.main(
                [
                    "-q",
                    "-o",
                    "addopts=",
                    test_target,
                ],
                plugins=[plugin],
            )
        )
    finally:
        _purge_loaded_modules()
        sys.path[:] = old_sys_path
    return code, plugin.failures


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Verify #1250 mutants in isolated scratch copies.")
    parser.add_argument(
        "--scratch-dir",
        type=Path,
        default=None,
        help="Optional scratch directory root for isolated mutant copies.",
    )
    args = parser.parse_args(argv)

    before_digests = _sha256_map(REPO_ROOT)
    for spec in MUTANTS:
        source_text = (REPO_ROOT / spec.rel_path).read_text(encoding="utf-8")
        occurrences = source_text.count(spec.target)
        if occurrences != 1:
            raise RuntimeError(
                f"Preflight failed: mutant {spec.name} target count={occurrences} "
                f"in {spec.rel_path}"
            )

    selected_controls = tuple(dict.fromkeys((spec.test_file, spec.test_node) for spec in MUTANTS))
    results: list[dict[str, object]] = []

    try:
        with ExitStack() as stack:
            if args.scratch_dir is not None:
                args.scratch_dir.mkdir(parents=True, exist_ok=True)
                base_dir = Path(
                    stack.enter_context(
                        tempfile.TemporaryDirectory(
                            prefix="verify_1250_mutants_", dir=str(args.scratch_dir)
                        )
                    )
                )
            else:
                base_dir = Path(
                    stack.enter_context(tempfile.TemporaryDirectory(prefix="verify_1250_mutants_"))
                )

            control_dir = base_dir / "control"
            _prepare_scratch_tree(control_dir)
            for test_file, test_node in selected_controls:
                control_code, control_failures = _run_pytest_in_scratch(
                    control_dir, test_file, test_node
                )
                if control_code != 0 or control_failures:
                    raise RuntimeError(
                        f"Unmutated control failed in scratch copy ({test_node}): "
                        f"code={control_code}, failures={control_failures}"
                    )

            for spec in MUTANTS:
                mutant_dir = base_dir / spec.name
                _prepare_scratch_tree(mutant_dir)
                target_file = mutant_dir / spec.rel_path
                original_text = target_file.read_text(encoding="utf-8")
                if original_text.count(spec.target) != 1:
                    raise RuntimeError(f"Mutant {spec.name} target not unique in {spec.rel_path}")
                mutated_text = original_text.replace(spec.target, spec.replacement, 1)
                if mutated_text == original_text:
                    raise RuntimeError(f"Mutant {spec.name} produced no change in {spec.rel_path}")
                target_file.write_text(mutated_text, encoding="utf-8")

                code, failures = _run_pytest_in_scratch(mutant_dir, spec.test_file, spec.test_node)
                assertion_killed = (
                    code == int(pytest.ExitCode.TESTS_FAILED)
                    and bool(failures)
                    and all(
                        when == "call" and exc_type == "AssertionError"
                        for _, when, exc_type in failures
                    )
                )
                if not assertion_killed:
                    raise RuntimeError(
                        f"Mutant {spec.name} was not cleanly killed by call-phase AssertionError: "
                        f"code={code}, failures={failures}"
                    )
                results.append(
                    {
                        "mutant": spec.name,
                        "rel_path": spec.rel_path,
                        "test_node": spec.test_node,
                        "exit_code": code,
                        "failures": failures,
                        "killed_by_assertion": assertion_killed,
                    }
                )
    finally:
        after_digests = _sha256_map(REPO_ROOT)
        if before_digests != after_digests:
            raise RuntimeError(
                f"Canonical digests changed during mutation run: "
                f"{before_digests} != {after_digests}"
            )

    summary = {
        "issue": 1250,
        "mutant_count": len(results),
        "all_killed_by_assertion": True,
        "canonical_digests_unchanged": True,
        "canonical_sha256": after_digests,
        "mutants": results,
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
