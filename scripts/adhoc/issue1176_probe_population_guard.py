"""Prove the population contradiction test in isolation (#1176; 2026-10-28).

Inputs: this checkout's matrix consumer and one synthetic regression test.
Outputs: fresh durable pytest log, isolated mutant and source-bound evidence.
Canonical source is never rewritten. Run one case per fresh process.
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
from mtp.workflow.output_roots import require_durable_output_path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = "src/mtp/matrix_overlap_coverage.py"
TEST_FILE = "tests/test_matrix_overlap_fidelity.py"
TEST_NAME = "test_contradictory_population_parameter_is_refused_independently"
ANCHOR = (
    'if (\n        "contrast_population" in params_by_name\n'
    '        and params_by_name["contrast_population"] != contrast_population\n    ):'
)


class Reports:
    """Only the intended call failure establishes that the guard mutant was killed."""

    def __init__(self) -> None:
        self.passed: list[str] = []
        self.failed: list[dict[str, str]] = []

    def pytest_runtest_logreport(self, report: pytest.TestReport) -> None:
        if report.failed:
            self.failed.append(
                {"nodeid": report.nodeid, "phase": report.when, "detail": report.longreprtext}
            )
        elif report.when == "call" and report.passed:
            self.passed.append(report.nodeid)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=("mutant", "baseline"))
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    if REPO_ROOT.resolve() != ROOT:
        raise ValueError("Imported mtp belongs to a different checkout")
    output = require_durable_output_path(args.output_dir.expanduser(), allowed=False)
    if output.is_relative_to(ROOT) or ROOT.is_relative_to(output):
        raise ValueError("Evidence output must be disjoint from the checkout")
    if args.output_dir.is_symlink():
        raise ValueError("Evidence output must not be a symlink")
    output.mkdir(parents=True, exist_ok=False)
    bound = (SOURCE, TEST_FILE, str(Path(__file__).relative_to(ROOT)), "pyproject.toml", "uv.lock")
    before = {name: file_sha256(ROOT / name) for name in bound}
    isolated = None
    if args.case == "mutant":
        source = (ROOT / SOURCE).read_text()
        if source.count(ANCHOR) != 1:
            raise ValueError("Expected exactly one population contradiction guard")
        isolated = output / "matrix_overlap_coverage.py"
        isolated.write_text(source.replace(ANCHOR, "if False:", 1))
        spec = importlib.util.spec_from_file_location("mtp.matrix_overlap_coverage", isolated)
        if spec is None or spec.loader is None:
            raise ValueError("Cannot load isolated matrix consumer")
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        vars(importlib.import_module("mtp"))["matrix_overlap_coverage"] = module
        if module.__file__ is None or Path(module.__file__).resolve() != isolated:
            raise ValueError("Mutant did not load from the isolated source")
    reports = Reports()
    log = io.StringIO()
    with redirect_stdout(log), redirect_stderr(log):
        result = int(
            pytest.main(
                ["-q", "-p", "no:cacheprovider", f"{ROOT / TEST_FILE}::{TEST_NAME}"],
                plugins=[reports],
            )
        )
    (output / "pytest.log").write_text(log.getvalue())
    after = {name: file_sha256(ROOT / name) for name in bound}
    if args.case == "baseline":
        intended = result == 0 and len(reports.passed) == 1 and not reports.failed
    else:
        intended = (
            result == 1
            and not reports.passed
            and len(reports.failed) == 1
            and reports.failed[0]["nodeid"].endswith(f"::{TEST_NAME}")
            and reports.failed[0]["phase"] == "call"
            and "DID NOT RAISE" in reports.failed[0]["detail"]
        )
    verified = intended and before == after
    record = {
        "issue": 1176,
        "case": args.case,
        "pytest_exit": result,
        "passed": reports.passed,
        "failed": reports.failed,
        "isolated_source": str(isolated) if isolated else None,
        "isolated_sha256": file_sha256(isolated) if isolated else None,
        "canonical_before": before,
        "canonical_after": after,
        "verified": verified,
    }
    (output / "evidence.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps({"case": args.case, "pytest_exit": result, "verified": verified}))
    return int(not verified)


if __name__ == "__main__":
    raise SystemExit(main())
