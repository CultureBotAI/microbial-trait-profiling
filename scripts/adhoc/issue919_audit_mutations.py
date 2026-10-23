"""Check #919 survivor mutations against audit tests (2026-10-23).

Inputs: current label-support source and synthetic tests. Outputs: pytest logs
and JSON mutation outcomes. Temporarily mutates the source, restoring exact
original bytes after each bounded test run; run only in the issue worktree.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import sys
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

import pytest

from mtp.paths import REPO_ROOT


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--mutation", choices=("all_input_taxa", "expanded_focal_family"), required=True
    )
    args = parser.parse_args()
    if "mtp.label_support_audit" in sys.modules:
        raise ValueError("Run each mutation in a fresh interpreter before importing the audit")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    source = REPO_ROOT / "src/mtp/label_support_audit.py"
    original = source.read_bytes()
    text = original.decode("utf-8")
    bound = 'bound_taxa = set(occurrences.loc[pair_keys.isin(paired), "taxon_id"])'
    focal = "FOCAL_TRAITS = tuple(AUDIT_TRAITS)\n"
    if text.count(bound) != 1 or text.count(focal) != 1 or text.count("FOCAL_TAXA = {") != 1:
        raise ValueError("Mutation targets changed; inspect the current implementation")
    mutants = {
        "all_input_taxa": (
            text.replace(bound, 'bound_taxa = set(occurrences["taxon_id"])'),
            "test_all_support_rules_share_zero_total_filtered_pairs",
        ),
        "expanded_focal_family": (
            text.replace(focal, "").replace("FOCAL_TAXA = {", focal + "FOCAL_TAXA = {"),
            "test_sidecar_records_sources_reconciliations_and_study_means",
        ),
    }
    name = args.mutation
    mutated, test = mutants[name]
    output = io.StringIO()
    try:
        source.write_bytes(mutated.encode("utf-8"))
        with redirect_stdout(output), redirect_stderr(output):
            code = pytest.main(["-q", f"{REPO_ROOT}/tests/test_label_support_audit.py::{test}"])
    finally:
        source.write_bytes(original)
    (args.output_dir / f"{name}.log").write_text(output.getvalue(), encoding="utf-8")
    record = {
        "original_source_sha256": hashlib.sha256(original).hexdigest(),
        "restored_source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "mutation": name,
        "pytest_exit": int(code),
        "killed": code == pytest.ExitCode.TESTS_FAILED,
    }
    (args.output_dir / f"{name}.json").write_text(
        json.dumps(record, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(record, indent=2))
    return 0 if record["killed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
