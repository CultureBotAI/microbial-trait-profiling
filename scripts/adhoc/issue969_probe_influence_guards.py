"""Probe report/accounting guards in isolated modules (#969; 2026-10-23).

Inputs: this checkout's influence producer and focused synthetic tests. Outputs:
pytest.log, an isolated mutant source and evidence.json under --output-dir.
Run one case per fresh process; canonical source is never rewritten.
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
JOINT = "tests/test_influence_extension_joint.py"
CORE = "tests/test_influence_extension_core.py"
REPORT_TEST = "test_summary_bh_by_cells_are_in_exact_columns"
MUTATIONS = {
    "loso_bh_by": (
        "f\"{row['baseline_q_bh_six']} / {row['baseline_q_by_six']} | \"",
        "f\"{row['baseline_q_by_six']} / {row['baseline_q_bh_six']} | \"",
        JOINT,
        REPORT_TEST,
        "AssertionError",
    ),
    "total_matched_percent": (
        'row["effect_reduction_percent"],\n            row["matched_effect_reduction_percent"],',
        'row["matched_effect_reduction_percent"],\n            row["effect_reduction_percent"],',
        JOINT,
        REPORT_TEST,
        "AssertionError",
    ),
    "fraction_as_percent": (
        'row["effect_reduction_percent"],',
        'row["effect_reduction_fraction"],',
        JOINT,
        REPORT_TEST,
        "AssertionError",
    ),
    "annotation_focal_mass": (
        'row["removed_annotation_mass"],\n            row["removed_trait_annotation_mass"],',
        'row["removed_trait_annotation_mass"],\n            row["removed_annotation_mass"],',
        JOINT,
        REPORT_TEST,
        "AssertionError",
    ),
    "exact_zero_only": (
        "abs(float(before)) <= TIE",
        "float(before) == 0",
        JOINT,
        "test_reduction_uses_protocol_tolerance_for_before_effect",
        "AssertionError",
    ),
    "post_protocol_check": (
        "if protocol_sha != file_sha256(PROTOCOL):\n"
        '        raise ValueError("Extension protocol changed during analysis")',
        'if False:\n        raise ValueError("Extension protocol changed during analysis")',
        CORE,
        "test_protocol_change_after_analysis_entry_is_rejected_before_publication",
        "DID NOT RAISE",
    ),
    "assumed_graph_presence": (
        "if lineage in present",
        "if True",
        JOINT,
        "test_joint_distinguishes_absent_lineage_from_present_without_scope_mass",
        "AssertionError",
    ),
    "summed_union_mass": (
        '"removed_annotation_mass": removed_total,',
        '"removed_annotation_mass": sum(sum(c.values()) for c in individual_mass.values()),',
        JOINT,
        "test_joint_reports_individual_lineage_mass_without_summing_overlap",
        "AssertionError",
    ),
}


class Results:
    """Distinguish intended call assertions from setup, collection and unrelated errors."""

    def __init__(self) -> None:
        self.failed: list[tuple[str, str, str]] = []
        self.passed: list[str] = []

    def pytest_runtest_logreport(self, report: pytest.TestReport) -> None:
        if report.failed:
            self.failed.append((report.nodeid, report.when, str(report.longrepr)))
        elif report.when == "call" and report.passed:
            self.passed.append(report.nodeid)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mutation", choices=["baseline", *MUTATIONS])
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if REPO_ROOT.resolve() != ROOT:
        raise ValueError("Imported mtp belongs to a different checkout")
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    source = ROOT / "src/mtp/influence_extension.py"
    paths = [source, ROOT / JOINT, ROOT / CORE, Path(__file__).resolve()]
    before = {str(path.relative_to(ROOT)): file_sha256(path) for path in paths}
    expected_test = marker = ""
    copied: Path | None = None
    if args.mutation == "baseline":
        targets = sorted({f"{ROOT / spec[2]}::{spec[3]}" for spec in MUTATIONS.values()})
    else:
        old, new, test_file, expected_test, marker = MUTATIONS[args.mutation]
        original = importlib.import_module("mtp.influence_extension")
        if original.__file__ is None or Path(original.__file__).resolve() != source:
            raise ValueError("Canonical influence producer is not from this checkout")
        content = source.read_text()
        if content.count(old) != 1:
            raise ValueError(f"Expected one mutation anchor: {old}")
        copied = output / "influence_extension.py"
        copied.write_text(content.replace(old, new))
        spec = importlib.util.spec_from_file_location("mtp.influence_extension", copied)
        if spec is None or spec.loader is None:
            raise ValueError("Cannot load isolated influence producer")
        module = importlib.util.module_from_spec(spec)
        sys.modules["mtp.influence_extension"] = module
        spec.loader.exec_module(module)
        # This unchanged function uses its canonical globals for source inventory
        # paths. Relocation must not mask the deliberately weakened logical guard.
        vars(module)["implementation_hashes"] = original.implementation_hashes
        package = importlib.import_module("mtp")
        vars(package)["influence_extension"] = module
        if module.__file__ is None or Path(module.__file__).resolve() != copied:
            raise ValueError("Mutation did not load from its isolated source")
        targets = [f"{ROOT / test_file}::{expected_test}"]
    results = Results()
    captured = io.StringIO()
    with redirect_stdout(captured), redirect_stderr(captured):
        exit_code = int(pytest.main(["-q", "-p", "no:cacheprovider", *targets], plugins=[results]))
    (output / "pytest.log").write_text(captured.getvalue())
    after = {str(path.relative_to(ROOT)): file_sha256(path) for path in paths}
    baseline = args.mutation == "baseline"
    intended = bool(results.passed) and not results.failed if baseline else bool(results.failed)
    if not baseline:
        intended = intended and all(
            expected_test in node and when == "call" and marker in detail
            for node, when, detail in results.failed
        )
    verified = exit_code == (0 if baseline else 1) and intended and before == after
    evidence = {
        "mutation": args.mutation,
        "pytest_exit": exit_code,
        "passed": results.passed,
        "failed": [node for node, _, _ in results.failed],
        "expected_test": expected_test,
        "expected_failure_marker": marker,
        "isolated_source": str(copied) if copied is not None else None,
        "isolated_sha256": file_sha256(copied) if copied is not None else None,
        "canonical_before": before,
        "canonical_after": after,
        "verified": verified,
    }
    (output / "evidence.json").write_text(json.dumps(evidence, indent=2) + "\n")
    print(json.dumps({"mutation": args.mutation, "verified": verified, "pytest_exit": exit_code}))
    return int(not verified)


if __name__ == "__main__":
    raise SystemExit(main())
