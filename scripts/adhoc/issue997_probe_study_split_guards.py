"""Probe #997 regression guards in isolated module copies; 2026-10-23.

Inputs: canonical study-split modules and synthetic regression tests. Outputs:
pytest.log and evidence.json beneath --output-dir. Canonical source bytes remain
unchanged; each invocation checks a baseline or one weakening and exact failures.
"""

from __future__ import annotations

import argparse
import importlib
import importlib.util
import io
import json
import sys
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

import pytest

from mtp.paths import REPO_ROOT, file_sha256

ROOT = Path(__file__).resolve().parents[2]
CORE_TESTS = ROOT / "tests/test_study_split_core.py"
SELECTION_TESTS = ROOT / "tests/test_study_split_selection.py"
GUARD_TEST = "test_change_after_analysis_reads_blocks_publication"
SELECTION_TEST = "test_repeated_split_report_has_no_aggregate_p_value"
TRAINING_TEST = "test_training_selection_rejects_the_heldout_only_trait_and_triple"
FINAL_GUARD = (
    "    if _provenance(paths, protocol, verified) != provenance:\n"
    '        raise ValueError("Study-split inputs or implementation changed during execution")\n'
)
MUTATIONS = {
    "training_selection_unrestricted": (
        "study_split",
        "        inputs,\n        split.training,\n        protocol=scan_protocol,\n",
        "        inputs,\n        None,\n        protocol=scan_protocol,\n",
        TRAINING_TEST,
    ),
    "final_recheck": ("study_split_run", FINAL_GUARD, "", GUARD_TEST),
    "no_traits": (
        "study_split",
        "    status = selection.status\n",
        "    status = (scanner.TOO_FEW_CONDITIONS if selection.status == scanner.NO_TRAITS "
        "else selection.status)\n",
        SELECTION_TEST,
    ),
    "no_disjoint_triple": (
        "study_split",
        '        status = "no_testable_study_disjoint_triple"\n',
        '        status = "selected"\n',
        SELECTION_TEST,
    ),
}
GUARD_NODE = f"{CORE_TESTS.relative_to(ROOT)}::{GUARD_TEST}"
SELECTION_NODE = f"{SELECTION_TESTS.relative_to(ROOT)}::{SELECTION_TEST}"
EXPECTED_FAILURES = {
    "training_selection_unrestricted": [f"{SELECTION_TESTS.relative_to(ROOT)}::{TRAINING_TEST}"],
    "final_recheck": [f"{GUARD_NODE}[source]", f"{GUARD_NODE}[code]"],
    "no_traits": [f"{SELECTION_NODE}[no_traits]"],
    "no_disjoint_triple": [f"{SELECTION_NODE}[no_testable_study_disjoint_triple]"],
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mutation", choices=["baseline", *MUTATIONS])
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if REPO_ROOT.resolve() != ROOT:
        raise AssertionError("Imported mtp belongs to a different checkout")
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    canonical = [
        *sorted((ROOT / "src/mtp/cellulose_decr").glob("*.py")),
        CORE_TESTS,
        SELECTION_TESTS,
        Path(__file__).resolve(),
    ]
    before = {str(path.relative_to(ROOT)): file_sha256(path) for path in canonical}
    expression = f"{GUARD_TEST} or {SELECTION_TEST} or {TRAINING_TEST}"
    original = None
    module_name = ""
    name = ""
    origins = {}
    captured = io.StringIO()
    try:
        if args.mutation != "baseline":
            name, old, new, expression = MUTATIONS[args.mutation]
            module_name = f"mtp.cellulose_decr.{name}"
            original = importlib.import_module(module_name)
            source = ROOT / f"src/mtp/cellulose_decr/{name}.py"
            if original.__file__ is None or Path(original.__file__).resolve() != source:
                raise AssertionError("Mutation target is not canonical")
            text = source.read_text()
            if text.count(old) != 1:
                raise AssertionError("Expected exactly one mutation target")
            copied = output / f"{name}.py"
            copied.write_text(text.replace(old, new))
            spec = importlib.util.spec_from_file_location(module_name, copied)
            if spec is None or spec.loader is None:
                raise AssertionError("Cannot load isolated mutant")
            mutant = importlib.util.module_from_spec(spec)
            sys.modules[module_name] = mutant
            spec.loader.exec_module(mutant)
            setattr(importlib.import_module("mtp.cellulose_decr"), name, mutant)
            if mutant.__file__ is None or Path(mutant.__file__).resolve() != copied:
                raise AssertionError("Mutation did not load from the isolated copy")
            origins = {"canonical": str(source), "isolated": str(copied)}
        with redirect_stdout(captured), redirect_stderr(captured):
            code = int(
                pytest.main(
                    [
                        "-q",
                        "-p",
                        "no:cacheprovider",
                        "--basetemp",
                        str(output / "pytest-tmp"),
                        str(CORE_TESTS),
                        str(SELECTION_TESTS),
                        "-k",
                        expression,
                    ]
                )
            )
    finally:
        if original is not None:
            sys.modules[module_name] = original
            setattr(importlib.import_module("mtp.cellulose_decr"), name, original)
    log = captured.getvalue()
    (output / "pytest.log").write_text(log)
    after = {str(path.relative_to(ROOT)): file_sha256(path) for path in canonical}
    expected = 0 if args.mutation == "baseline" else 1
    failed_nodes = [line.split()[1] for line in log.splitlines() if line.startswith("FAILED ")]
    expected_nodes = EXPECTED_FAILURES.get(args.mutation, [])
    marker = "DID NOT RAISE" if args.mutation == "final_recheck" else "AssertionError"
    if args.mutation == "training_selection_unrestricted":
        marker = "Training selection must exclude the held-out-only trait and triple"
    verified = (
        code == expected
        and before == after
        and set(failed_nodes) == set(expected_nodes)
        and (expected == 0 or marker in log)
    )
    evidence = {
        "mutation": args.mutation,
        "exit_code": code,
        "expected_code": expected,
        "test_expression": expression,
        "failed_nodes": failed_nodes,
        "expected_failed_nodes": expected_nodes,
        "origins": origins,
        "canonical_before": before,
        "canonical_after": after,
        "verified": verified,
    }
    (output / "evidence.json").write_text(json.dumps(evidence, indent=2) + "\n")
    print(json.dumps(evidence, indent=2))
    return int(not verified)


if __name__ == "__main__":
    raise SystemExit(main())
