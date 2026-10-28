"""Verify Delta-z statistical regression strength (#1115/#1154, 2026-10-28).

Inputs: synthetic mixed-zero, exact-stool count and active-support tests.
Outputs: durable pytest logs and JSON outcomes under --output. Each mutation
changes an in-memory callable in a spawned process; no production file changes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing
import platform
from collections import Counter
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from scipy import __version__ as SCIPY_VERSION
from scipy import stats

from mtp.delta_z_derived import tables

ROOT = Path(__file__).resolve().parents[2]
TEST_FILE = "tests/test_delta_z_derived_tables.py"
TESTS = (
    "test_locked_environment_records_its_reviewed_wilcoxon_contract",
    "test_mixed_zero_wilcoxon_pins_wilcox_and_no_continuity_correction",
    "test_row_three_counts_every_trait_cell_against_exact_raw_stool_metadata",
    "test_active_study_support_boundary_is_exactly_fifty",
)
MUTANTS = (
    "baseline",
    "pratt",
    "zsplit",
    "correction_true",
    "casefold_stool",
    "contains_stool",
    "strictly_over_fifty",
)
SOURCES = (
    "src/mtp/delta_z_derived/tables.py",
    TEST_FILE,
    "scripts/adhoc/issue1115_statistical_mutants.py",
    "uv.lock",
)


class Reports:
    def __init__(self) -> None:
        self.counts: Counter[str] = Counter()
        self.failures: list[str] = []

    def pytest_runtest_logreport(self, report: pytest.TestReport) -> None:
        if report.when == "call":
            self.counts[report.outcome] += 1
        if report.failed:
            self.failures.append(report.nodeid)


def hashes() -> dict[str, str]:
    return {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in SOURCES}


def worker(mutant: str, output: Path) -> int:
    before = hashes()
    wilcoxon = stats.wilcoxon
    metadata = tables._metadata
    summary = tables.trait_summary

    def wrong_wilcoxon(*args: Any, **kwargs: Any) -> Any:
        if mutant in {"pratt", "zsplit"}:
            kwargs["zero_method"] = mutant
        else:
            kwargs["correction"] = True
        return wilcoxon(*args, **kwargs)

    def wrong_stool(experiments: pd.DataFrame, pairs: pd.DataFrame) -> pd.DataFrame:
        result = metadata(experiments, pairs).copy()
        mask = (
            result.body_site.str.casefold().eq("feces")
            if mutant == "casefold_stool"
            else result.body_site.str.contains("Feces", regex=False)
        )
        result.loc[mask, "body_site"] = "Feces"
        return result

    def wrong_boundary(studies: pd.DataFrame, evidence: pd.DataFrame) -> pd.DataFrame:
        result = summary(studies, evidence)
        result["sdz_active_at_least_50"] = result.sdz_active_studies.gt(50)
        return result

    reports = Reports()
    with pytest.MonkeyPatch.context() as patch:
        if mutant in {"pratt", "zsplit", "correction_true"}:
            patch.setattr(stats, "wilcoxon", wrong_wilcoxon)
        elif mutant in {"casefold_stool", "contains_stool"}:
            patch.setattr(tables, "_metadata", wrong_stool)
        elif mutant == "strictly_over_fifty":
            patch.setattr(tables, "trait_summary", wrong_boundary)
        code = int(
            pytest.main(
                ["-q", "--tb=short", *(f"{ROOT / TEST_FILE}::{name}" for name in TESTS)],
                plugins=[reports],
            )
        )
    if before != hashes():
        raise AssertionError("Production, test or runner bytes changed during the campaign")
    record = {
        "mutant": mutant,
        "pytest_exit": code,
        "counts": dict(reports.counts),
        "failures": reports.failures,
        "source_sha256": before,
        "runtime": {
            "python": platform.python_version(),
            "scipy": SCIPY_VERSION,
            "scipy_wilcoxon_source_sha256": hashlib.sha256(
                Path(stats.__file__).with_name("_wilcoxon.py").read_bytes()
            ).hexdigest(),
        },
    }
    (output / "result.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    return code


def worker_entrypoint(mutant: str, output: Path) -> None:
    with (
        (output / "pytest.log").open("w", encoding="utf-8") as log,
        redirect_stdout(log),
        redirect_stderr(log),
    ):
        raise SystemExit(worker(mutant, output))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    context = multiprocessing.get_context("spawn")
    records = []
    for mutant in MUTANTS:
        output = args.output / mutant
        output.mkdir()
        process = context.Process(target=worker_entrypoint, args=(mutant, output))
        process.start()
        process.join()
        result_path = output / "result.json"
        record = json.loads(result_path.read_text()) if result_path.exists() else {}
        expected = 0 if mutant == "baseline" else 1
        valid = process.exitcode == expected and (expected == 0 or bool(record.get("failures")))
        records.append({**record, "expected_outcome": valid, "worker_exit": process.exitcode})
        print(f"{mutant}: exit={process.exitcode}; {record.get('counts', {})}", flush=True)
    (args.output / "summary.json").write_text(
        json.dumps(records, indent=2) + "\n", encoding="utf-8"
    )
    return int(not all(record["expected_outcome"] for record in records))


if __name__ == "__main__":
    raise SystemExit(main())
