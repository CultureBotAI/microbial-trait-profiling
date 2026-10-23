"""Probe multiplicity guards using isolated module copies (#974; 2026-10-23).

Inputs: this checkout's canonical adapter sources and regression tests. Outputs:
pytest.log and evidence.json under --output-dir. No canonical source is rewritten;
each fresh process loads at most one weakened copy and checks its intended failure.
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
TEST = ROOT / "tests/test_multiplicity_regression_guards.py"
MUTATIONS = {
    "native_support": (
        "native",
        "6 <= support[feature]",
        "1 <= support[feature]",
        "test_native_selection_guards_are_independent",
        "Native support_floor admitted",
    ),
    "native_variance": (
        "native",
        "and variance[feature] > 0",
        "and variance[feature] >= 0",
        "test_native_selection_guards_are_independent",
        "Native fraction_variance admitted",
    ),
    "ora_background": (
        "ora",
        'feature not in candidate_set or int(row["background_k"]) != prevalence[feature]',
        "feature not in candidate_set",
        "test_ora_background_prevalence_is_checked",
        "Failed: DID NOT RAISE ValueError",
    ),
    "ora_fold": (
        "ora",
        'q <= 0.001 and (fold >= 1.5 if tail == "upper" else fold <= 0.5)',
        "q <= 0.001",
        "test_ora_native_hit_needs_fold_as_well_as_q",
        "Failed: DID NOT RAISE ValueError",
    ),
    "eligible_reconcile": (
        "families",
        "if eligible != up + down + tied:",
        "if False:",
        "test_sign_eligible_count_reconciles_with_tallies",
        "Failed: DID NOT RAISE ValueError",
    ),
    "final_source_recheck": (
        "stage",
        "if file_sha256(path) != digest or source(path).identity not in consumed:",
        "if False:",
        "test_final_source_recheck_blocks_published_manifest",
        "Failed: DID NOT RAISE ValueError",
    ),
    "migration_threshold": (
        "stage",
        "tested and native_q <= 0.05",
        'tested and native_q <= native_row["native_threshold"]',
        "test_migration_native_q005_is_not_native_reported_threshold",
        "Native q<=0.05 comparison changed its threshold",
    ),
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mutation", choices=["baseline", *MUTATIONS])
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--legacy-tests", action="store_true")
    parser.add_argument("--expect-survive", action="store_true")
    args = parser.parse_args()
    if REPO_ROOT.resolve() != ROOT:
        raise AssertionError("Imported mtp belongs to a different checkout")
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    canonical = [
        *sorted((ROOT / "src/mtp/multiplicity").glob("*.py")),
        Path(__file__).resolve(),
        TEST,
        ROOT / "tests/test_multiplicity_adapters.py",
        ROOT / "tests/multiplicity_fixture.py",
    ]
    before = {str(path.relative_to(ROOT)): file_sha256(path) for path in canonical}
    expected_test, marker = "", ""
    origins: dict[str, str] = {}
    if args.mutation != "baseline":
        name, old, new, expected_test, marker = MUTATIONS[args.mutation]
        module_name = f"mtp.multiplicity.{name}"
        original = importlib.import_module(module_name)
        source = ROOT / f"src/mtp/multiplicity/{name}.py"
        if original.__file__ is None or Path(original.__file__).resolve() != source.resolve():
            raise AssertionError("Mutation target is not imported from this checkout")
        text = source.read_text()
        if text.count(old) != 1:
            raise AssertionError(f"Expected exactly one mutation target: {old}")
        copied = output / f"{name}.py"
        copied.write_text(text.replace(old, new))
        spec = importlib.util.spec_from_file_location(module_name, copied)
        if spec is None or spec.loader is None:
            raise AssertionError("Cannot load isolated mutant module")
        mutated = importlib.util.module_from_spec(spec)
        sys.modules[module_name] = mutated
        spec.loader.exec_module(mutated)
        if name == "stage":
            # Inventory paths are relative to the canonical module's __file__.
            # Preserve this unmutated function and its canonical globals so a
            # relocated mutant is tested for the intended guard, not its path.
            vars(mutated)["code_sha256"] = original.code_sha256
        package = importlib.import_module("mtp.multiplicity")
        setattr(package, name, mutated)
        origins = {"canonical": str(source), "isolated": str(copied)}
        if mutated.__file__ is None or Path(mutated.__file__).resolve() != copied:
            raise AssertionError("Mutation did not load from the isolated copy")
    tests = ROOT / "tests/test_multiplicity_adapters.py" if args.legacy_tests else TEST
    arguments = ["-q", "-p", "no:cacheprovider", str(tests)]
    if expected_test and not args.legacy_tests:
        arguments += ["-k", expected_test]
    captured = io.StringIO()
    with redirect_stdout(captured), redirect_stderr(captured):
        code = int(pytest.main(arguments))
    log = captured.getvalue()
    (output / "pytest.log").write_text(log)
    after = {str(path.relative_to(ROOT)): file_sha256(path) for path in canonical}
    expected_code = 0 if args.mutation == "baseline" or args.expect_survive else 1
    intended = expected_code == 0 or (marker in log and expected_test in log)
    verified = code == expected_code and before == after and intended
    record = {
        "mutation": args.mutation,
        "legacy_tests": args.legacy_tests,
        "exit_code": code,
        "expected_code": expected_code,
        "expected_test": expected_test,
        "expected_failure_marker": marker,
        "module_origins": origins,
        "canonical_before": before,
        "canonical_after": after,
        "verified": verified,
    }
    (output / "evidence.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record, indent=2))
    return int(not verified)


if __name__ == "__main__":
    raise SystemExit(main())
