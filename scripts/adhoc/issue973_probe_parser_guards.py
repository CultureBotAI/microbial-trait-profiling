"""Compare arm-parser guards with an isolated original runtime (#973; 2026-10-23).

Inputs: this checkout's tests and optional saved original runtime source. Outputs:
pytest.log and evidence.json in --output-dir, including module/source identities.
The canonical checkout is never rewritten; failures must be test-call assertions.
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
SELECTED = (
    "rejects_help_without_printing or condition_module_uses_its_execution_parser "
    "or rejects_missing_required_heatmap_direction"
)


class Outcomes:
    """Require assertion failures during test calls, not collection or setup errors."""

    def __init__(self) -> None:
        self.failed: list[str] = []
        self.non_assertion_failures: list[str] = []
        self.passed: list[str] = []

    def pytest_runtest_logreport(self, report: pytest.TestReport) -> None:
        if report.when == "call" and report.passed:
            self.passed.append(report.nodeid)
        if report.failed:
            text = report.longreprtext
            if report.when == "call" and ("AssertionError" in text or "DID NOT RAISE" in text):
                self.failed.append(report.nodeid)
            else:
                self.non_assertion_failures.append(report.nodeid)

    def pytest_collectreport(self, report: pytest.CollectReport) -> None:
        if report.failed:
            self.non_assertion_failures.append(report.nodeid)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runtime-source", type=Path)
    parser.add_argument("--expected-failures", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if REPO_ROOT.resolve() != ROOT:
        raise ValueError("Imported mtp belongs to another checkout")
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    canonical = [
        ROOT / "src/mtp/workflow/arms/runtime.py",
        ROOT / "scripts/create_condition_trait_module_heatmap.py",
        ROOT / "tests/test_arm_direction_parser.py",
        ROOT / "tests/test_workflow_subset_integration.py",
        Path(__file__).resolve(),
    ]
    before = {str(path.relative_to(ROOT)): file_sha256(path) for path in canonical}
    original = importlib.import_module("mtp.workflow.arms.runtime")
    if original.__file__ is None or Path(original.__file__).resolve() != canonical[0]:
        raise ValueError("Runtime import is not from this checkout")
    origins = {"canonical": str(canonical[0])}
    if args.runtime_source is not None:
        copied = output / "original_runtime.py"
        copied.write_bytes(args.runtime_source.read_bytes())
        spec = importlib.util.spec_from_file_location("mtp.workflow.arms.runtime", copied)
        if spec is None or spec.loader is None:
            raise ValueError("Cannot load isolated runtime source")
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        vars(importlib.import_module("mtp.workflow.arms"))["runtime"] = module
        if module.__file__ is None or Path(module.__file__).resolve() != copied:
            raise ValueError("Isolated runtime did not load")
        origins.update(isolated=str(copied), source_sha256=file_sha256(copied))
    outcomes = Outcomes()
    captured = io.StringIO()
    with redirect_stdout(captured), redirect_stderr(captured):
        result = int(
            pytest.main(
                [
                    "-q",
                    "-p",
                    "no:cacheprovider",
                    str(canonical[2]),
                    str(canonical[3]),
                    "-k",
                    SELECTED,
                ],
                plugins=[outcomes],
            )
        )
    (output / "pytest.log").write_text(captured.getvalue())
    after = {str(path.relative_to(ROOT)): file_sha256(path) for path in canonical}
    expected_exit = 1 if args.expected_failures else 0
    verified = (
        result == expected_exit
        and len(outcomes.failed) == args.expected_failures
        and len(outcomes.failed) + len(outcomes.passed) == 15
        and not outcomes.non_assertion_failures
        and before == after
    )
    evidence = {
        "module_origins": origins,
        "pytest_exit": result,
        "expected_failures": args.expected_failures,
        "failed_calls": outcomes.failed,
        "passed_calls": outcomes.passed,
        "non_assertion_failures": outcomes.non_assertion_failures,
        "canonical_before": before,
        "canonical_after": after,
        "verified": verified,
    }
    (output / "evidence.json").write_text(json.dumps(evidence, indent=2) + "\n")
    print(json.dumps({"verified": verified, "failed": len(outcomes.failed)}))
    return int(not verified)


if __name__ == "__main__":
    raise SystemExit(main())
