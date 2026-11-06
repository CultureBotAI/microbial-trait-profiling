"""Run four reviewed-row regressions against an explicit producer (2026-10-23, #957).

Inputs: an unchanged producer script and the current synthetic coverage tests.
Outputs: pytest.log and evidence.json beneath --output-dir; no scientific inputs.
Uses isolated module bindings, never rewrites a producer or canonical test file.
"""

from __future__ import annotations

import argparse
import contextlib
import importlib.util
import json
import sys
from pathlib import Path
from types import ModuleType

import pytest

from mtp.paths import REPO_ROOT, file_sha256


def load_module(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if Path(module.__file__ or "").resolve() != path.resolve():
        raise RuntimeError("Unexpected producer import origin")
    return module


class Outcomes:
    """Inspect actual call failures, excluding collection/setup errors."""

    def __init__(self) -> None:
        self.passed: list[str] = []
        self.failed: list[str] = []
        self.unexpected: list[str] = []

    def pytest_runtest_logreport(self, report: pytest.TestReport) -> None:
        if report.when == "call" and report.passed:
            self.passed.append(report.nodeid)
        if report.failed:
            if report.when == "call" and "assert row.diagnosis_category == expected" in str(
                report.longrepr
            ):
                self.failed.append(report.nodeid)
            else:
                self.unexpected.append(report.nodeid)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--producer", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--expected-failures", type=int, choices=(0, 4), required=True)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    fixture_path = REPO_ROOT / "tests/test_taxa_coverage_histogram.py"
    regression_path = REPO_ROOT / "tests/test_taxa_coverage_taxonomy.py"
    current_producer = REPO_ROOT / "scripts/create_taxa_coverage_histogram.py"
    observed = [args.producer, fixture_path, regression_path, current_producer, Path(__file__)]
    before = {str(path.resolve()): file_sha256(path) for path in observed}
    binding = "test_taxa_coverage_histogram"
    previous = sys.modules.get(binding)
    outcomes = Outcomes()
    try:
        fixture = load_module(binding, fixture_path)
        vars(fixture)["coverage"] = load_module("issue957_selected_producer", args.producer)
        sys.modules[binding] = fixture
        with (
            (output / "pytest.log").open("w", encoding="utf-8") as log,
            contextlib.redirect_stdout(log),
            contextlib.redirect_stderr(log),
        ):
            status = pytest.main(
                [
                    str(regression_path),
                    "-q",
                    "-k",
                    "four_reviewed_taxa",
                    "--basetemp",
                    str(output / "fixtures"),
                ],
                plugins=[outcomes],
            )
    finally:
        if previous is None:
            sys.modules.pop(binding, None)
        else:
            sys.modules[binding] = previous
    after = {str(path.resolve()): file_sha256(path) for path in observed}
    valid = (
        before == after
        and not outcomes.unexpected
        and len(outcomes.failed) == args.expected_failures
        and len(outcomes.passed) == 4 - args.expected_failures
        and int(status) == (1 if args.expected_failures else 0)
    )
    (output / "evidence.json").write_text(
        json.dumps(
            {
                "issue": 957,
                "producer": str(args.producer.resolve()),
                "expected_failures": args.expected_failures,
                "pytest_exit": int(status),
                "passed": outcomes.passed,
                "intended_failures": outcomes.failed,
                "unexpected_failures": outcomes.unexpected,
                "sha256_before": before,
                "sha256_after": after,
                "verified": valid,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return 0 if valid else 1


if __name__ == "__main__":
    raise SystemExit(main())
