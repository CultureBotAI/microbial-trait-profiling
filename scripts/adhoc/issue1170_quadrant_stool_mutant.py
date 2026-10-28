"""Prove exact raw-stool quadrant selection (#1170, 2026-10-28).

Inputs: the derived-table source and independent row-2/row-3 synthetic tests.
Outputs: durable pytest logs and a hash-bound JSON proof under --output.
Only the quadrant is_stool assignment is mutated, in a spawned test process;
the original production bytes are restored in finally before returning.
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
SOURCE = ROOT / "src/mtp/delta_z_derived/tables.py"
TEST_FILE = "tests/test_delta_z_derived_tables.py"
TESTS = (
    "test_row_two_quadrants_exclude_multisite_nonaxis_pairs",
    "test_row_three_counts_every_trait_cell_against_exact_raw_stool_metadata",
)
BEFORE = 'result["is_stool"] = result.body_site.eq("Feces")'
AFTER = 'result["is_stool"] = result.body_site.str.contains("Feces", regex=False)'
SOURCES = (
    "src/mtp/delta_z_derived/tables.py",
    TEST_FILE,
    "scripts/adhoc/issue1170_quadrant_stool_mutant.py",
    "uv.lock",
)


def digest(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def hashes() -> dict[str, str]:
    return {name: digest((ROOT / name).read_bytes()) for name in SOURCES}


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
                ["-q", "--tb=short", *(f"{ROOT / TEST_FILE}::{name}" for name in TESTS)],
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


def intended_result(case: str, record: dict[str, Any]) -> bool:
    if case == "baseline":
        return record.get("exit_code") == 0 and record.get("counts") == {"passed": 2}
    failures = record.get("failures", [])
    return bool(
        record.get("exit_code") == 1
        and record.get("counts") == {"passed": 1, "failed": 1}
        and len(failures) == 1
        and failures[0]["nodeid"].endswith(f"::{TESTS[0]}")
        and failures[0]["phase"] == "call"
        and "n_comparisons" in failures[0]["failure"]
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    original = SOURCE.read_bytes()
    if original.count(BEFORE.encode()) != 1:
        raise ValueError("Expected exactly one quadrant-specific stool assignment")
    args.output.mkdir(parents=True, exist_ok=False)
    initial = hashes()
    records = []
    context = multiprocessing.get_context("spawn")
    for case in ("baseline", "contains_quadrant_only"):
        output = args.output / case
        output.mkdir()
        mutated = (
            original if case == "baseline" else original.replace(BEFORE.encode(), AFTER.encode(), 1)
        )
        try:
            SOURCE.write_bytes(mutated)
            process = context.Process(target=worker, args=(output,))
            process.start()
            process.join()
            result_path = output / "pytest.json"
            record = json.loads(result_path.read_text()) if result_path.exists() else {}
            record.update(
                {
                    "case": case,
                    "mutated_source_sha256": digest(mutated),
                    "process_exit_code": process.exitcode,
                    "intended_outcome": process.exitcode == record.get("exit_code")
                    and intended_result(case, record),
                }
            )
        finally:
            SOURCE.write_bytes(original)
        if hashes() != initial:
            raise AssertionError("Source, tests or helper changed or were not restored")
        records.append(record)
        print(f"{case}: {record.get('counts')}; intended={record['intended_outcome']}", flush=True)
    report = {
        "issue": 1170,
        "mutation": {"before": BEFORE, "after": AFTER},
        "source_restored": SOURCE.read_bytes() == original,
        "source_sha256": initial,
        "cases": records,
        "environment": {
            "python": platform.python_version(),
            **{name: version(name) for name in ("numpy", "pandas", "scipy", "pytest")},
        },
    }
    (args.output / "summary.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return int(not all(record["intended_outcome"] for record in records))


if __name__ == "__main__":
    raise SystemExit(main())
