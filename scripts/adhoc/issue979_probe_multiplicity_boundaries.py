"""Probe issue #979 in isolated copies; 2026-10-23.

Inputs: canonical multiplicity modules and synthetic regression tests. Outputs:
pytest.log and evidence.json under --output-dir. No canonical bytes are edited.
Each fresh invocation checks one weakening and its intended failing assertion.
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
from types import ModuleType

import pytest

from mtp.paths import REPO_ROOT, file_sha256

ROOT = Path(__file__).resolve().parents[2]
GUARDS = ROOT / "tests/test_multiplicity_regression_guards.py"
ADAPTERS = ROOT / "tests/test_multiplicity_adapters.py"
FOLD = 'fold >= 1.5 if tail == "upper" else fold <= 0.5'
METADATA = (
    "sources",
    "if expected is not None and artifact.sha256 != expected:",
    "if False:",
)
FINAL = (
    "stage",
    "path: sources.expected_digest(path) for family in inputs for path in family.source_paths",
    "path: file_sha256(path) for family in inputs for path in family.source_paths",
)
MUTATIONS = {
    "fold_upper": (
        [("ora", FOLD, FOLD.replace(">=", ">"))],
        "fold_boundary",
        "q/effect thresholds",
    ),
    "fold_lower": (
        [("ora", FOLD, FOLD.replace("<=", "<"))],
        "fold_boundary",
        "q/effect thresholds",
    ),
    "fold_both": (
        [("ora", FOLD, FOLD.replace(">=", ">").replace("<=", "<"))],
        "fold_boundary",
        "q/effect thresholds",
    ),
    "metadata_hash": ([METADATA], "non_rows_recorded_hash_checks", "Regex pattern did not match"),
    "final_recorded_hash": ([FINAL], "non_rows_recorded_hash_checks", "DID NOT RAISE ValueError"),
    "both_recorded_hashes": (
        [METADATA, FINAL],
        "non_rows_recorded_hash_checks",
        "DID NOT RAISE ValueError",
    ),
    "missing_stage": (
        [("stage", "if Path(__file__).resolve() not in paths:", "if False:")],
        "requires_its_own_tracked_module",
        "DID NOT RAISE ValueError",
    ),
    "hardcoded_prefix": (
        [("stage", "str(name).startswith(prefix)", 'str(name).startswith("src/mtp/")')],
        "tracks_the_actual_nested_package_prefix",
        "Tracked code inventory must include",
    ),
}
BOUNDARY_NODE = f"{GUARDS.relative_to(ROOT)}::test_ora_native_hit_includes_exact_fold_boundary"
HASH_NODE = f"{GUARDS.relative_to(ROOT)}::test_non_rows_recorded_hash_checks"
INVENTORY_NODE = f"{ADAPTERS.relative_to(ROOT)}::test_code_identity_requires_its_own_tracked_module"
EXPECTED_FAILURES = {
    "fold_upper": [f"{BOUNDARY_NODE}[upper-1.5]"],
    "fold_lower": [f"{BOUNDARY_NODE}[lower-0.5]"],
    "fold_both": [f"{BOUNDARY_NODE}[upper-1.5]", f"{BOUNDARY_NODE}[lower-0.5]"],
    "metadata_hash": [f"{HASH_NODE}[metadata]"],
    "final_recorded_hash": [f"{HASH_NODE}[final]"],
    "both_recorded_hashes": [f"{HASH_NODE}[metadata]", f"{HASH_NODE}[final]"],
    "missing_stage": [f"{INVENTORY_NODE}[False]", f"{INVENTORY_NODE}[True]"],
    "hardcoded_prefix": [
        f"{ADAPTERS.relative_to(ROOT)}::test_code_identity_tracks_the_actual_nested_package_prefix"
    ],
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
        *sorted((ROOT / "src/mtp/multiplicity").glob("*.py")),
        Path(__file__).resolve(),
        GUARDS,
        ADAPTERS,
        ROOT / "tests/multiplicity_fixture.py",
    ]
    before = {str(path.relative_to(ROOT)): file_sha256(path) for path in canonical}
    originals: dict[str, ModuleType] = {}
    origins: dict[str, dict[str, str]] = {}
    expression = "fold_boundary or non_rows_recorded_hash_checks or code_identity"
    marker = ""
    captured = io.StringIO()
    try:
        if args.mutation != "baseline":
            mutations, expression, marker = MUTATIONS[args.mutation]
            for name, old, new in mutations:
                module_name = f"mtp.multiplicity.{name}"
                original = importlib.import_module(module_name)
                source = ROOT / f"src/mtp/multiplicity/{name}.py"
                if original.__file__ is None or Path(original.__file__).resolve() != source:
                    raise AssertionError("Mutation target is not canonical")
                originals[name] = original
                text = source.read_text()
                if text.count(old) != 1:
                    raise AssertionError(f"Expected one mutation target: {old}")
                copied = output / f"{name}.py"
                copied.write_text(text.replace(old, new))
                spec = importlib.util.spec_from_file_location(module_name, copied)
                if spec is None or spec.loader is None:
                    raise AssertionError("Cannot load isolated mutant")
                mutant = importlib.util.module_from_spec(spec)
                sys.modules[module_name] = mutant
                spec.loader.exec_module(mutant)
                if name == "stage" and args.mutation not in {"missing_stage", "hardcoded_prefix"}:
                    # Keep the unmutated inventory function's canonical globals;
                    # relocating the target must not cause an unrelated hash failure.
                    vars(mutant)["code_sha256"] = original.code_sha256
                setattr(importlib.import_module("mtp.multiplicity"), name, mutant)
                if mutant.__file__ is None or Path(mutant.__file__).resolve() != copied:
                    raise AssertionError("Mutation did not load from its isolated copy")
                origins[name] = {"canonical": str(source), "isolated": str(copied)}
        with redirect_stdout(captured), redirect_stderr(captured):
            code = int(
                pytest.main(
                    [
                        "-q",
                        "-p",
                        "no:cacheprovider",
                        str(GUARDS),
                        str(ADAPTERS),
                        "-k",
                        expression,
                    ]
                )
            )
    finally:
        package = importlib.import_module("mtp.multiplicity")
        for name, original in originals.items():
            sys.modules[f"mtp.multiplicity.{name}"] = original
            setattr(package, name, original)
    log = captured.getvalue()
    (output / "pytest.log").write_text(log)
    after = {str(path.relative_to(ROOT)): file_sha256(path) for path in canonical}
    expected = 0 if args.mutation == "baseline" else 1
    failed_nodes = [line.split()[1] for line in log.splitlines() if line.startswith("FAILED ")]
    expected_nodes = EXPECTED_FAILURES.get(args.mutation, [])
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
        "expected_failure_marker": marker,
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
