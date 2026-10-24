"""Probe #976 contracts with isolated source copies (2026-10-23).

Inputs: canonical sources and synthetic regression tests. Outputs: pytest.log
and evidence.json under --output-dir. Canonical source bytes are never edited.
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
TEST = ROOT / "tests/test_multiplicity_cross_product.py"
REASONS = ROOT / "tests/test_multiplicity_reason_vocabulary.py"
MUTATIONS = {
    "noel_core_tallies": (
        "mtp.multiplicity.families",
        'if row["trait"] not in tallies or actual != tallies[row["trait"]]:',
        "if False:",
        TEST,
        "test_cross_product_status_and_correction_tampering_is_rejected and noel",
    ),
    "jed_ledger": (
        "mtp.multiplicity.families",
        "expected_support = asdict(ledger_support[identity])",
        "expected_support = {}",
        TEST,
        "test_cross_product_status_and_correction_tampering_is_rejected and jed",
    ),
    "native_moments": (
        "mtp.multiplicity.native",
        'if (\n                row["study_test_status"] != expected.status',
        'if False and (\n                row["study_test_status"] != expected.status',
        TEST,
        "test_cross_product_status_and_correction_tampering_is_rejected and native",
    ),
    "ora_native_bh": (
        "mtp.multiplicity.ora",
        "if exported != expected:",
        "if False:",
        TEST,
        "test_cross_product_status_and_correction_tampering_is_rejected and ora",
    ),
    "modern_reason": (
        "mtp.multiple_testing",
        'row.reason_vocabulary == "analysis_reason_v1"',
        "False",
        REASONS,
        "test_modern_reasons_reject_unknown_whitespace_and_missing_values",
    ),
    "corrected_reader": (
        "mtp.multiplicity.inventory",
        "            validate_raw_hypothesis(\n",
        "            (lambda row: row)(\n",
        REASONS,
        "test_corrected_reader_validates_raw_semantics",
    ),
}


class Outcomes:
    def __init__(self) -> None:
        self.intended: list[str] = []
        self.unexpected: list[str] = []

    def pytest_runtest_logreport(self, report: pytest.TestReport) -> None:
        if report.failed:
            if report.when == "call" and "DID NOT RAISE" in str(report.longrepr):
                self.intended.append(report.nodeid)
            else:
                self.unexpected.append(report.nodeid)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mutation", choices=["baseline", *MUTATIONS])
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if REPO_ROOT.resolve() != ROOT:
        raise RuntimeError("Imported mtp belongs to a different checkout")
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    watched = [
        *sorted((ROOT / "src/mtp/multiplicity").glob("*.py")),
        ROOT / "src/mtp/multiple_testing.py",
        ROOT / "src/mtp/native_composition.py",
        ROOT / "src/mtp/trait_composition.py",
        ROOT / "src/mtp/trait_screen.py",
        TEST,
        REASONS,
        Path(__file__).resolve(),
    ]
    before = {str(path.relative_to(ROOT)): file_sha256(path) for path in watched}
    origins = {}
    arguments = ["-q", "-p", "no:cacheprovider", str(TEST), str(REASONS)]
    if args.mutation != "baseline":
        name, old, replacement, test, selection = MUTATIONS[args.mutation]
        original = importlib.import_module(name)
        source = ROOT / "src" / Path(*name.split(".")).with_suffix(".py")
        if original.__file__ is None or Path(original.__file__).resolve() != source:
            raise RuntimeError("Unexpected mutation source import")
        content = source.read_text()
        if content.count(old) != 1:
            raise RuntimeError(f"Expected one mutation target: {old}")
        copied = output / source.name
        copied.write_text(content.replace(old, replacement))
        spec = importlib.util.spec_from_file_location(name, copied)
        if spec is None or spec.loader is None:
            raise RuntimeError("Cannot load isolated source")
        mutated = importlib.util.module_from_spec(spec)
        sys.modules[name] = mutated
        spec.loader.exec_module(mutated)
        if name == "mtp.multiple_testing":
            vars(mutated)["POLICY_PATH"] = original.POLICY_PATH
            vars(mutated)["_code_sha256"] = original._code_sha256
        parent, attribute = name.rsplit(".", 1)
        setattr(importlib.import_module(parent), attribute, mutated)
        if mutated.__file__ is None or Path(mutated.__file__).resolve() != copied:
            raise RuntimeError("Isolated import identity mismatch")
        origins = {"canonical": str(source), "isolated": str(copied)}
        arguments = ["-q", "-p", "no:cacheprovider", str(test), "-k", selection]
    outcomes = Outcomes()
    captured = io.StringIO()
    with redirect_stdout(captured), redirect_stderr(captured):
        status = int(pytest.main(arguments, plugins=[outcomes]))
    (output / "pytest.log").write_text(captured.getvalue())
    after = {str(path.relative_to(ROOT)): file_sha256(path) for path in watched}
    expected = 0 if args.mutation == "baseline" else 1
    verified = (
        before == after
        and status == expected
        and not outcomes.unexpected
        and (args.mutation == "baseline" or bool(outcomes.intended))
    )
    evidence = {
        "mutation": args.mutation,
        "exit_code": status,
        "expected_exit": expected,
        "module_origins": origins,
        "intended_failures": outcomes.intended,
        "unexpected_failures": outcomes.unexpected,
        "canonical_before": before,
        "canonical_after": after,
        "verified": verified,
    }
    (output / "evidence.json").write_text(json.dumps(evidence, indent=2) + "\n")
    print(json.dumps({"mutation": args.mutation, "verified": verified, "exit_code": status}))
    return int(not verified)


if __name__ == "__main__":
    raise SystemExit(main())
