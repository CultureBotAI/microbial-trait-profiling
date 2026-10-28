"""Prove #1129 tests detect sole-source 1000612 loss (2026-10-27).

Inputs: synthetic pytest fixtures and test_trait_canonicalization_profiles.py.
Outputs: durable pytest logs and JSON kill/survival records under --output.
Each child replaces only the in-memory edge-reader reference, filtering private
copies of its fixture inputs. No tracked production file is edited.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import multiprocessing
from collections import Counter
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from typing import Any

import pytest

from mtp import taxa_trait_matrix as matrix

ROOT = Path(__file__).resolve().parents[2]
TEST = "tests/test_trait_canonicalization_profiles.py"
MUTANTS = ("baseline", "drop_direct_1000612", "drop_strain_1000612")


class Reports:
    def __init__(self) -> None:
        self.counts: Counter[str] = Counter()
        self.failures: list[str] = []

    def pytest_runtest_logreport(self, report: pytest.TestReport) -> None:
        if report.when == "call":
            self.counts[report.outcome] += 1
        if report.failed:
            self.failures.append(report.nodeid)


def worker(mutant: str, output: Path) -> int:
    original = matrix.read_kg_trait_edges
    source = ROOT / "src/mtp/taxa_trait_matrix.py"
    source_before = source.read_bytes()
    filtered = output / "filtered-fixtures"
    filtered.mkdir()

    def drop_edges(path: Path, *args: Any, **kwargs: Any) -> matrix.KgTraitEdges:
        identity = hashlib.sha256(str(path.resolve()).encode()).hexdigest()
        target = filtered / f"{identity}.tsv"
        with path.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.reader(handle, delimiter="\t"))
        prefix = "NCBITaxon:" if mutant == "drop_direct_1000612" else "kgmicrobe.strain:"
        with target.open("w", encoding="utf-8", newline="") as handle:
            csv.writer(handle, delimiter="\t", lineterminator="\n").writerows(
                row
                for row in rows
                if not (len(row) >= 4 and row[1].startswith(prefix) and row[3] == "METPO:1000612")
            )
        return original(target, *args, **kwargs)

    reports = Reports()
    with pytest.MonkeyPatch.context() as patch:
        if mutant != "baseline":
            patch.setattr(matrix, "read_kg_trait_edges", drop_edges)
        code = int(pytest.main(["-q", "--tb=short", str(ROOT / TEST)], plugins=[reports]))
    if source.read_bytes() != source_before:
        raise AssertionError("Production source changed during mutation test")
    record = {
        "mutant": mutant,
        "pytest_exit": code,
        "counts": dict(reports.counts),
        "failures": reports.failures,
        "source_sha256": hashlib.sha256(source_before).hexdigest(),
    }
    (output / "result.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    return code


def worker_entrypoint(mutant: str, output: Path) -> None:
    """Run one mutant in a fresh Python interpreter with a private test log."""
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
        valid = process.exitcode == expected and (
            mutant == "baseline" or bool(record.get("failures"))
        )
        records.append({**record, "expected_outcome": valid, "worker_exit": process.exitcode})
        print(f"{mutant}: exit={process.exitcode}; {record.get('counts', {})}", flush=True)
        if not valid:
            break
    (args.output / "summary.json").write_text(
        json.dumps(records, indent=2) + "\n", encoding="utf-8"
    )
    return int(len(records) != len(MUTANTS) or not all(row["expected_outcome"] for row in records))


if __name__ == "__main__":
    raise SystemExit(main())
