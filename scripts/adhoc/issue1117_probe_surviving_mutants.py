"""Prove M19/M22/M25 regression strength for #1117 (2026-10-28).

Inputs: the frozen matrix module and three synthetic regression tests.
Outputs: durable per-case pytest logs and a hash-bound JSON campaign report.
Each exact source mutation runs in a fresh spawned process; a finally block
restores the original production bytes before the next case or any exit.
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
SOURCE = ROOT / "src/mtp/matrix_overlap_coverage.py"
TEST_FILE = "tests/test_matrix_overlap_fidelity.py"
TESTS = {
    "M19": "test_wrong_direction_basis_is_refused_before_publication",
    "M22": "test_pd_profile_keys_exclude_a_decoy_taxon_that_changes_the_overlap",
    "M25": "test_automatic_commit_is_consumer_and_unread_product_is_not_declared_consumed",
}
MUTATIONS = {
    "M19": ('if direction_basis != "disease_vs_comparison":', "if False:"),
    "M22": (
        "t12_exp, t12_taxa = population.active_experiments, population.active_taxa",
        "t12_exp, t12_taxa = tier12_experiments, tier12_taxa",
    ),
    "M25": ('"code_commit": resolved_code_commit,', '"code_commit": bundle.run.git_commit,'),
}
BOUND_FILES = (
    "src/mtp/matrix_overlap_coverage.py",
    "src/mtp/matrix_overlap_inputs.py",
    "src/mtp/human_host_corpus.py",
    "src/mtp/workflow/coverage.py",
    "src/mtp/workflow/output_roots.py",
    "tests/test_matrix_overlap_coverage.py",
    TEST_FILE,
    "scripts/adhoc/issue1117_probe_surviving_mutants.py",
    "uv.lock",
)


def digest(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def source_hashes() -> dict[str, str]:
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


def worker(output: Path) -> None:
    reports = Reports()
    with (
        (output / "pytest.log").open("w", encoding="utf-8") as log,
        redirect_stdout(log),
        redirect_stderr(log),
    ):
        code = int(
            pytest.main(
                ["-q", "--tb=short", *(f"{ROOT / TEST_FILE}::{name}" for name in TESTS.values())],
                plugins=[reports],
            )
        )
        (output / "pytest.json").write_text(
            json.dumps(
                {"exit_code": code, "counts": dict(reports.counts), "failures": reports.failures},
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
    raise SystemExit(code)


def intended_failure(case: str, record: dict[str, Any]) -> bool:
    if case == "baseline":
        return record.get("exit_code") == 0 and record.get("counts") == {"passed": 4}
    failures = record.get("failures", [])
    failed_count = 2 if case == "M25" else 1
    if record.get("exit_code") != 1 or record.get("counts") != {
        "passed": 4 - failed_count,
        "failed": failed_count,
    }:
        return False
    if len(failures) != failed_count or any(failure["phase"] != "call" for failure in failures):
        return False
    if any(failure["nodeid"].split("::")[-1].split("[")[0] != TESTS[case] for failure in failures):
        return False
    text = "\n".join(failure["failure"] for failure in failures)
    if case == "M19":
        return "DID NOT RAISE" in text
    if case == "M22":
        return "mean_jaccard" in text and ("comparison failed" in text or "AssertionError" in text)
    return "consumer-current" in text and "core-parent" in text and "AssertionError" in text


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    # Import only the storage guard; the mutated module is first imported in each worker.
    from mtp.workflow.output_roots import require_durable_output_path

    output = require_durable_output_path(args.output, allowed=False)
    output.mkdir(parents=True, exist_ok=False)
    original = SOURCE.read_bytes()
    initial = source_hashes()
    records = []
    context = multiprocessing.get_context("spawn")
    for case in ("baseline", *MUTATIONS):
        directory = output / case
        directory.mkdir()
        mutated = original
        replacement = None
        if case != "baseline":
            before, after = MUTATIONS[case]
            if original.count(before.encode()) != 1:
                raise ValueError(f"{case} must match exactly one production expression")
            mutated = original.replace(before.encode(), after.encode(), 1)
            replacement = {"before": before, "after": after}
        try:
            SOURCE.write_bytes(mutated)
            process = context.Process(target=worker, args=(directory,))
            process.start()
            process.join()
            result_path = directory / "pytest.json"
            result = json.loads(result_path.read_text()) if result_path.exists() else {}
            result.update(
                {
                    "case": case,
                    "mutation": replacement,
                    "mutated_source_sha256": digest(mutated),
                    "process_exit_code": process.exitcode,
                    "intended_outcome": process.exitcode == result.get("exit_code")
                    and intended_failure(case, result),
                }
            )
        finally:
            SOURCE.write_bytes(original)
        if source_hashes() != initial:
            raise AssertionError("Frozen source/test/helper bytes changed or were not restored")
        records.append(result)
        print(f"{case}: {result.get('counts')}; intended={result['intended_outcome']}", flush=True)
    report = {
        "issue": 1117,
        "source_restored": SOURCE.read_bytes() == original,
        "source_sha256": initial,
        "cases": records,
        "environment": {
            "python": platform.python_version(),
            **{name: version(name) for name in ("numpy", "pandas", "scipy", "pytest")},
        },
        "M22_expectation": (
            "baseline mean 1/3; mutant admits inactive study and decoy, yielding 1/12"
        ),
    }
    (output / "summary.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return int(not all(record["intended_outcome"] for record in records))


if __name__ == "__main__":
    raise SystemExit(main())
