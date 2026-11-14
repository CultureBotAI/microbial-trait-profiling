"""Isolated scratch-copy mutation verifier for Issues #966 and #1184.

Purpose:
  Verify that tests/test_adhoc_mutation_probe_isolation.py kills 10 distinct mutants across:
  - scripts/adhoc/issue959_probe_coverage_mutants.py
  - scripts/adhoc/issue1182_probe_scan_bindings.py
  - scripts/adhoc/issue1170_quadrant_stool_mutant.py
  inside an isolated scratch directory without modifying the canonical repository worktree.
Inputs:
  - scripts/adhoc/issue959_probe_coverage_mutants.py
  - scripts/adhoc/issue1182_probe_scan_bindings.py
  - scripts/adhoc/issue1170_quadrant_stool_mutant.py
  - tests/test_adhoc_mutation_probe_isolation.py
Outputs:
  - Stdout summary of control + 10 mutant runs; exits 0 iff all 10 mutants fail with AssertionError.
Issues: #966, #1184
Date: 2026-11-09
"""

from __future__ import annotations

import contextlib
import io
import shutil
import sys
import tempfile
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
COPIED_FILES: tuple[str, ...] = (
    "scripts/adhoc/issue959_probe_coverage_mutants.py",
    "scripts/adhoc/issue1182_probe_scan_bindings.py",
    "scripts/adhoc/issue1170_quadrant_stool_mutant.py",
    "scripts/create_taxa_coverage_histogram.py",
    "tests/test_taxa_coverage_histogram.py",
    "tests/test_cellulose_decr_selection_scan.py",
    "tests/test_cellulose_decr_selection_scan_scopes.py",
    "tests/test_delta_z_derived_tables.py",
    "tests/test_adhoc_mutation_probe_isolation.py",
    "pyproject.toml",
    "uv.lock",
)
MTP_CHECK_LINE = (
    '    if mtp_file is None or Path(mtp_file).resolve().parent != (ROOT / "src/mtp").resolve():\n'
)

MUTANTS: tuple[tuple[str, str, str, str, str], ...] = (
    (
        "M01_959_reintroduce_inplace_source_write",
        "scripts/adhoc/issue959_probe_coverage_mutants.py",
        '    isolated_source.write_text(producer_text, encoding="utf-8")\n',
        '    isolated_source.write_text(producer_text, encoding="utf-8")\n'
        "    if False:\n"
        '        SOURCE.write_text(producer_text, encoding="utf-8")\n',
        "test_adhoc_mutation_probes_define_isolation_and_checkout_guards",
    ),
    (
        "M02_959_drop_mtp_checkout_check",
        "scripts/adhoc/issue959_probe_coverage_mutants.py",
        MTP_CHECK_LINE,
        "    if False:\n",
        "test_issue959_probe_rejects_foreign_import_origin_before_writing",
    ),
    (
        "M03_959_drop_repo_root_checkout_check",
        "scripts/adhoc/issue959_probe_coverage_mutants.py",
        '    if (REPO_ROOT / "src/mtp").resolve() != (ROOT / "src/mtp").resolve():\n',
        "    if False:\n",
        "test_issue959_probe_rejects_foreign_import_origin_before_writing",
    ),
    (
        "M04_959_allow_collection_errors_as_verified",
        "scripts/adhoc/issue959_probe_coverage_mutants.py",
        "    if collection_errors:\n        return False\n",
        "    if False:\n        return False\n",
        "test_issue959_probe_is_verified_case_and_failure_evidence",
    ),
    (
        "M05_959_allow_setup_phase_failure_as_verified",
        "scripts/adhoc/issue959_probe_coverage_mutants.py",
        '            item["phase"] == "call"\n',
        "            True\n",
        "test_issue959_probe_is_verified_case_and_failure_evidence",
    ),
    (
        "M06_1182_reintroduce_inplace_source_write",
        "scripts/adhoc/issue1182_probe_scan_bindings.py",
        "        isolated_source.write_bytes(mutated)\n",
        "        isolated_source.write_bytes(mutated)\n"
        "        if False:\n"
        "            (ROOT / SOURCE).write_bytes(mutated)\n",
        "test_adhoc_mutation_probes_define_isolation_and_checkout_guards",
    ),
    (
        "M07_1182_drop_checkout_identity_check",
        "scripts/adhoc/issue1182_probe_scan_bindings.py",
        MTP_CHECK_LINE,
        "    if False:\n",
        "test_issue1182_probe_isolates_mutations_without_touching_canonical_source",
    ),
    (
        "M08_1182_drop_module_restore_in_finally",
        "scripts/adhoc/issue1182_probe_scan_bindings.py",
        "        sys.modules[module_name] = canonical_module\n",
        "        pass\n",
        "test_issue1182_probe_isolates_mutations_without_touching_canonical_source",
    ),
    (
        "M09_1170_reintroduce_inplace_source_write",
        "scripts/adhoc/issue1170_quadrant_stool_mutant.py",
        "        isolated_source.write_bytes(mutated)\n",
        "        isolated_source.write_bytes(mutated)\n"
        "        if False:\n"
        "            SOURCE.write_bytes(mutated)\n",
        "test_adhoc_mutation_probes_define_isolation_and_checkout_guards",
    ),
    (
        "M10_1170_drop_checkout_identity_check",
        "scripts/adhoc/issue1170_quadrant_stool_mutant.py",
        MTP_CHECK_LINE,
        "    if False:\n",
        "test_issue1170_probe_runs_isolated_mutations_without_touching_canonical_source",
    ),
)


def _stage_scratch(scratch: Path) -> None:
    for rel in COPIED_FILES:
        src = REPO_ROOT / rel
        dst = scratch / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
    for link_name in ("src", ".git"):
        (scratch / link_name).symlink_to(REPO_ROOT / link_name)


def _run_pytest_node(scratch: Path, node: str) -> tuple[int, str]:
    test_file = scratch / "tests/test_adhoc_mutation_probe_isolation.py"
    for mod_name in list(sys.modules):
        if "test_adhoc_mutation_probe_isolation" in mod_name or mod_name.startswith("issue"):
            sys.modules.pop(mod_name, None)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        rc = int(
            pytest.main(
                [
                    "-q",
                    "-p",
                    "no:cacheprovider",
                    "-o",
                    "addopts=",
                    f"{test_file}::{node}",
                ]
            )
        )
    return rc, buf.getvalue()


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="verify_966_1184_mutants_") as tmp:
        scratch = Path(tmp).resolve()
        _stage_scratch(scratch)
        originals = {
            rel: (scratch / rel).read_text(encoding="utf-8")
            for rel in (
                "scripts/adhoc/issue959_probe_coverage_mutants.py",
                "scripts/adhoc/issue1182_probe_scan_bindings.py",
                "scripts/adhoc/issue1170_quadrant_stool_mutant.py",
            )
        }

        # Positive controls on the isolated scratch copy
        for control_node in (
            "test_adhoc_mutation_probes_define_isolation_and_checkout_guards",
            "test_issue959_probe_rejects_foreign_import_origin_before_writing",
            "test_issue959_probe_is_verified_case_and_failure_evidence",
            "test_issue1182_probe_isolates_mutations_without_touching_canonical_source",
            "test_issue1170_probe_runs_isolated_mutations_without_touching_canonical_source",
        ):
            rc, out = _run_pytest_node(scratch, control_node)
            if rc != 0:
                print(f"Control {control_node} failed (rc={rc}):\n{out}")
                return 1
        loaded = sys.modules.get("test_adhoc_mutation_probe_isolation")
        expected_test = (scratch / "tests/test_adhoc_mutation_probe_isolation.py").resolve()
        if loaded is None or Path(loaded.__file__ or "").resolve() != expected_test:
            print("Scratch origin check failed for test_adhoc_mutation_probe_isolation.py")
            return 1
        print(f"Controls + scratch-origin binding passed in {scratch}")

        killed = 0
        for name, rel_path, target, replacement, test_node in MUTANTS:
            script_file = scratch / rel_path
            base_text = originals[rel_path]
            if base_text.count(target) != 1:
                print(f"Target count != 1 for {name} in {rel_path}")
                return 1
            try:
                script_file.write_text(base_text.replace(target, replacement, 1), encoding="utf-8")
                rc, out = _run_pytest_node(scratch, test_node)
            finally:
                script_file.write_text(base_text, encoding="utf-8")
            if rc != 0 and "AssertionError" in out:
                killed += 1
                print(f"  [KILLED] {name} (rc={rc})")
            else:
                print(f"  [SURVIVED/INVALID] {name} (rc={rc}):\n{out}")
                return 1

        print(f"Summary: {killed}/{len(MUTANTS)} mutants killed.")
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
