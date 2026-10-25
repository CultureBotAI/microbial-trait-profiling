"""Probe sidecar guards using isolated source copies (#1077, #1081; 2026-10-25).

Inputs: this checkout's sidecar modules and regression tests. Outputs: a mutant
module, pytest.log and evidence.json in a fresh --output-dir; pytest and SQLite
use the explicit durable --scratch-dir. Canonical files are never rewritten.
Run baseline as the positive control; each mutant must fail its intended test.
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

from mtp.paths import file_sha256
from mtp.three_test_table.stage import durable_scratch, require_durable

ROOT = Path(__file__).resolve().parents[2]
MUTATIONS = {
    "ora_list_size_direction": (
        "ora",
        '"ora_direction_pairs_n_inc_gt_dec": population.inc_larger_pairs if background_k else 0,\n'
        '        "ora_direction_pairs_n_inc_lt_dec": '
        "population.dec_larger_pairs if background_k else 0,",
        '"ora_direction_pairs_n_inc_gt_dec": population.dec_larger_pairs if background_k else 0,\n'
        '        "ora_direction_pairs_n_inc_lt_dec": '
        "population.inc_larger_pairs if background_k else 0,",
        "test_ora_reports_list_size_balance_for_only_direction_eligible_pairs",
        "test_ora_reports_list_size_balance_for_only_direction_eligible_pairs",
        "List-size balance must retain increased/decreased ordering",
        1,
    ),
    "ora_nullable_counts": (
        "ora",
        '.astype(dict.fromkeys(_NULLABLE_COUNTS, "Int64"))',
        "",
        "test_ora_nullable_counts_preserve_integer_text_and_missing_background",
        "test_ora_nullable_counts_preserve_integer_text_and_missing_background",
        "ORA call counts must use nullable integers",
        1,
    ),
    "native_nullable_counts": (
        "native",
        '.astype(dict.fromkeys(_NULLABLE_COUNTS, "Int64"))',
        "",
        "test_native_nullable_counts_preserve_producer_integers_and_unselected_blanks",
        "test_native_nullable_counts_preserve_producer_integers_and_unselected_blanks",
        "Native selected counts must use nullable integers",
        2,
    ),
    "carrier_nullable_counts": (
        "carriers",
        '.astype(dict.fromkeys(_NULLABLE_COUNTS, "Int64"))',
        "",
        "test_carrier_nullable_counts_preserve_unresolved_and_observed_values",
        "test_carrier_nullable_counts_preserve_unresolved_and_observed_values",
        "Carrier counts must use nullable integers",
        2,
    ),
    "carrier_evidence_count": (
        "carriers",
        "if count != expected[key] or bool(direct) != support[key][0]:",
        "if bool(direct) != support[key][0]:",
        "test_carrier_rejects_self_rehashed_incomplete_or_false_evidence and missing_evidence",
        "test_carrier_rejects_self_rehashed_incomplete_or_false_evidence",
        "Failed: DID NOT RAISE",
        1,
    ),
    "carrier_direct_identity": (
        "carriers",
        "if count != expected[key] or bool(direct) != support[key][0]:",
        "if count != expected[key]:",
        "test_carrier_rejects_count_preserving_direct_carrier_conflict",
        "test_carrier_rejects_count_preserving_direct_carrier_conflict",
        "Failed: DID NOT RAISE",
        2,
    ),
    "carrier_manifest_recheck": (
        "carriers",
        "if file_sha256(manifest) != manifest_digest:",
        "if False:",
        "test_carrier_rejects_provenance_changed_during_read",
        "test_carrier_rejects_provenance_changed_during_read",
        "Failed: DID NOT RAISE",
        1,
    ),
    "ora_untestable_placeholder": (
        "ora",
        'record["ora_direction_p"] if record["ora_direction_status"] == "tested" else 1.0',
        'record["ora_direction_p"] if record["ora_direction_status"] == "tested" else 0.5',
        "test_ora_correction_receives_literal_one_for_each_untestable_trait",
        "test_ora_correction_receives_literal_one_for_each_untestable_trait",
        "Untestable ORA hypotheses must enter corrections at p=1",
        1,
    ),
}


class FailedCalls:
    """Distinguish the intended assertion failures from setup or collection errors."""

    def __init__(self) -> None:
        self.nodes: list[str] = []

    def pytest_runtest_logreport(self, report: pytest.TestReport) -> None:
        if report.when == "call" and report.failed:
            self.nodes.append(report.nodeid)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mutation", choices=["baseline", *MUTATIONS])
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--scratch-dir", required=True, type=Path)
    args = parser.parse_args()
    output = require_durable(args.output_dir)
    scratch = require_durable(args.scratch_dir)
    if output == scratch or output.is_relative_to(scratch) or scratch.is_relative_to(output):
        raise ValueError("Probe output and scratch must be separate directories")
    output.mkdir(parents=True, exist_ok=False)
    scratch.mkdir(parents=True, exist_ok=False)
    modules = ("carriers", "ora", "native")
    tests = [ROOT / f"tests/test_three_test_{name}.py" for name in modules]
    canonical = [
        *tests,
        *(ROOT / f"src/mtp/three_test_table/{name}.py" for name in modules),
        Path(__file__).resolve(),
    ]
    before = {str(path.relative_to(ROOT)): file_sha256(path) for path in canonical}
    arguments = ["-q", "-p", "no:cacheprovider", "--basetemp", str(scratch / "pytest")]
    expected_test, marker, expected_failures = "", "", 0
    if args.mutation != "baseline":
        name, old, new, selection, expected_test, marker, expected_failures = MUTATIONS[
            args.mutation
        ]
        module_name = f"mtp.three_test_table.{name}"
        original = importlib.import_module(module_name)
        source = ROOT / f"src/mtp/three_test_table/{name}.py"
        if original.__file__ is None or Path(original.__file__).resolve() != source.resolve():
            raise AssertionError("Mutation source is not imported from this checkout")
        text = source.read_text(encoding="utf-8")
        if text.count(old) != 1:
            raise AssertionError(f"Expected exactly one mutation target: {old}")
        copied = output / f"{name}.py"
        copied.write_text(text.replace(old, new), encoding="utf-8")
        spec = importlib.util.spec_from_file_location(module_name, copied)
        if spec is None or spec.loader is None:
            raise AssertionError("Cannot load isolated mutation")
        mutated = importlib.util.module_from_spec(spec)
        sys.modules[module_name] = mutated
        spec.loader.exec_module(mutated)
        package = importlib.import_module("mtp.three_test_table")
        setattr(package, name, mutated)
        if mutated.__file__ is None or Path(mutated.__file__).resolve() != copied:
            raise AssertionError("Mutation did not load from its isolated copy")
        tests = [ROOT / f"tests/test_three_test_{name}.py"]
        arguments += ["-k", selection]
    arguments.extend(str(path) for path in tests)
    captured = io.StringIO()
    failures = FailedCalls()
    with durable_scratch(scratch), redirect_stdout(captured), redirect_stderr(captured):
        code = int(pytest.main(arguments, plugins=[failures]))
    log = captured.getvalue()
    (output / "pytest.log").write_text(log, encoding="utf-8")
    after = {str(path.relative_to(ROOT)): file_sha256(path) for path in canonical}
    expected_code = 0 if args.mutation == "baseline" else 1
    intended = not expected_failures or (
        len(failures.nodes) == expected_failures
        and all(node.split("::", 1)[1].split("[", 1)[0] == expected_test for node in failures.nodes)
        and marker in log
    )
    verified = code == expected_code and before == after and intended
    record = {
        "mutation": args.mutation,
        "pytest_exit_code": code,
        "expected_exit_code": expected_code,
        "expected_test": expected_test,
        "failed_call_nodes": failures.nodes,
        "expected_failure_marker": marker,
        "canonical_before": before,
        "canonical_after": after,
        "verified": verified,
    }
    (output / "evidence.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"mutation": args.mutation, "verified": verified, "pytest_exit_code": code}))
    return int(not verified)


if __name__ == "__main__":
    raise SystemExit(main())
