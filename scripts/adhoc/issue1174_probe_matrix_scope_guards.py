"""Prove matrix scope-guard regression strength for #1174 (2026-10-28).

Inputs: frozen production source and three synthetic scope-guard regression tests.
Outputs: durable pytest logs, exact source mutations and a hash-bound JSON report.
Each mutant runs in a fresh spawned process. Only the mutated source is written,
and its original bytes are restored in finally before inspecting the result.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing
import platform
from collections import Counter
from contextlib import redirect_stderr, redirect_stdout
from importlib.metadata import version
from pathlib import Path
from typing import Any

import pytest

ROOT = Path(__file__).resolve().parents[2]
TEST_FILE = "tests/test_matrix_overlap_fidelity.py"
COVERAGE_TEST = "test_inactive_direction_is_excluded_from_coverage"
PD_TEST = "test_inactive_direction_is_excluded_from_tier12_pd_overlap"
CONTRADICTION_TEST = "test_contradictory_direction_parameter_is_refused_independently"
TESTS = (COVERAGE_TEST, PD_TEST, CONTRADICTION_TEST)
MUTATIONS = {
    "experiment_only_profile_key": (
        "src/mtp/matrix_overlap_inputs.py",
        "selected = _keys(tier_taxa, PROFILE_KEYS).isin(_keys(active, PROFILE_KEYS))",
        "selected = _keys(tier_taxa, EXPERIMENT_KEYS).isin(_keys(active, EXPERIMENT_KEYS))",
    ),
    "disabled_direction_parameter_refusal": (
        "src/mtp/matrix_overlap_coverage.py",
        'if "direction_basis" in params_by_name '
        'and params_by_name["direction_basis"] != direction_basis:',
        "if False:",
    ),
}
BOUND_FILES = (
    "src/mtp/matrix_overlap_coverage.py",
    "src/mtp/matrix_overlap_inputs.py",
    "src/mtp/trait_lists.py",
    "src/mtp/human_host_corpus.py",
    "src/mtp/workflow/coverage.py",
    "src/mtp/workflow/output_roots.py",
    "tests/test_matrix_overlap_coverage.py",
    TEST_FILE,
    "scripts/adhoc/issue1174_probe_matrix_scope_guards.py",
    "pyproject.toml",
    "uv.lock",
)


def digest(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def input_hashes() -> dict[str, str]:
    return {name: digest((ROOT / name).read_bytes()) for name in BOUND_FILES}


class Reports:
    def __init__(self) -> None:
        self.counts: Counter[str] = Counter()
        self.failures: list[dict[str, str]] = []

    def pytest_runtest_logreport(self, report: pytest.TestReport) -> None:
        if report.when == "call":
            self.counts[report.outcome] += 1
        if report.failed:
            self.failures.append(
                {"nodeid": report.nodeid, "phase": report.when, "failure": report.longreprtext}
            )


def pytest_arguments() -> list[str]:
    return ["-q", "--tb=short", *(f"{ROOT / TEST_FILE}::{name}" for name in TESTS)]


def worker(output: Path) -> None:
    reports = Reports()
    with (
        (output / "pytest.log").open("w", encoding="utf-8") as log,
        redirect_stdout(log),
        redirect_stderr(log),
    ):
        arguments = pytest_arguments()
        code = int(pytest.main(arguments, plugins=[reports]))
        (output / "pytest.json").write_text(
            json.dumps(
                {
                    "entrypoint": "pytest.main",
                    "arguments": arguments,
                    "exit_code": code,
                    "counts": dict(reports.counts),
                    "failures": reports.failures,
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
    raise SystemExit(code)


def intended_failure(case: str, record: dict[str, Any]) -> bool:
    if case == "baseline":
        return (
            record.get("exit_code") == 0
            and record.get("counts") == {"passed": 3}
            and record.get("failures") == []
        )
    expected_tests = (
        {COVERAGE_TEST, PD_TEST} if case == "experiment_only_profile_key" else {CONTRADICTION_TEST}
    )
    failures = record.get("failures", [])
    if record.get("exit_code") != 1 or record.get("counts") != {
        "passed": 3 - len(expected_tests),
        "failed": len(expected_tests),
    }:
        return False
    if len(failures) != len(expected_tests) or any(
        failure["phase"] != "call" for failure in failures
    ):
        return False
    by_test = {failure["nodeid"].split("::")[-1]: failure["failure"] for failure in failures}
    if set(by_test) != expected_tests:
        return False
    if case == "experiment_only_profile_key":
        return (
            "tier12_excluded_occurrence_rows" in by_test[COVERAGE_TEST]
            and "assert 1 == 2" in by_test[COVERAGE_TEST]
            and "mean_jaccard" in by_test[PD_TEST]
            and "assert 0.3333333333333333 == 0.5" in by_test[PD_TEST]
        )
    return "DID NOT RAISE" in by_test[CONTRADICTION_TEST]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    # Mutated analysis modules are imported only in the fresh workers.
    from mtp.workflow.output_roots import require_durable_output_path

    output = require_durable_output_path(args.output, allowed=False)
    if output.is_relative_to(ROOT) or ROOT.is_relative_to(output):
        raise ValueError("Campaign output must be disjoint from the source working tree")
    originals = {name: (ROOT / name).read_bytes() for name, _, _ in MUTATIONS.values()}
    for case, (name, before, _) in MUTATIONS.items():
        if originals[name].count(before.encode()) != 1:
            raise ValueError(f"{case} must match exactly one production expression")
    initial = input_hashes()
    output.mkdir(parents=True, exist_ok=False)
    records = []
    context = multiprocessing.get_context("spawn")
    for case in ("baseline", *MUTATIONS):
        directory = output / case
        directory.mkdir()
        target: Path | None = None
        mutation = None
        if case != "baseline":
            name, before, after = MUTATIONS[case]
            target = ROOT / name
            mutated = originals[name].replace(before.encode(), after.encode(), 1)
            mutation = {
                "file": name,
                "before": before,
                "after": after,
                "original_sha256": digest(originals[name]),
                "mutated_sha256": digest(mutated),
            }
        process = context.Process(target=worker, args=(directory,))
        try:
            if target is not None:
                target.write_bytes(mutated)
            process.start()
            process.join(timeout=300)
            if process.is_alive():
                raise TimeoutError(f"{case}: pytest worker exceeded 300 seconds")
            result_path = directory / "pytest.json"
            result = json.loads(result_path.read_text()) if result_path.exists() else {}
            result.update(
                {
                    "case": case,
                    "mutation": mutation,
                    "process_exit_code": process.exitcode,
                    "intended_outcome": process.exitcode == result.get("exit_code")
                    and intended_failure(case, result),
                }
            )
        finally:
            if process.is_alive():
                process.terminate()
                process.join(timeout=5)
            if target is not None:
                target.write_bytes(originals[target.relative_to(ROOT).as_posix()])
        if input_hashes() != initial:
            raise AssertionError("Frozen source/test/helper bytes changed or were not restored")
        records.append(result)
        print(f"{case}: {result.get('counts')}; intended={result['intended_outcome']}", flush=True)
        if case == "baseline" and not result["intended_outcome"]:
            break
    report = {
        "issue": 1174,
        "input_sha256": initial,
        "source_restored": all(
            (ROOT / name).read_bytes() == raw for name, raw in originals.items()
        ),
        "cases": records,
        "environment": {
            "python": platform.python_version(),
            **{name: version(name) for name in ("numpy", "pandas", "scipy", "pytest")},
        },
    }
    (output / "summary.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return int(len(records) != 3 or not all(record["intended_outcome"] for record in records))


if __name__ == "__main__":
    raise SystemExit(main())
