"""Probe matrix follow-up guards with isolated modules (#1175; 2026-10-28).

Inputs: current producer/helper source and focused synthetic regression tests.
Outputs: fresh durable pytest logs, copied mutant source and hash-bound evidence.
Never rewrites tracked source; each CLI invocation uses one fresh interpreter.
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
CASES = {
    "orphan": (
        "src/mtp/matrix_overlap_inputs.py",
        '_covered(active, tier12_experiments, "Active profile")',
        "pass",
        "tests/test_matrix_overlap_fidelity.py",
        "test_active_profile_needs_tier12_metadata_even_when_source_metadata_exists",
        2,
    ),
    "provenance": (
        "scripts/adhoc/issue1117_compare_archived_output.py",
        'if provenance.get("declared_traits") != declared:',
        "if False:",
        "tests/test_issue1117_compare_archived_output.py",
        "test_provenance_declared_traits_are_checked_independently",
        4,
    ),
}


def load_copy(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ValueError("Cannot load the isolated guard mutant")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    if module.__file__ is None or Path(module.__file__).resolve() != path:
        raise ValueError("Mutant did not load from its isolated source")
    return module


class Reports:
    """Record call failures separately from collection/setup errors."""

    def __init__(self, parity_copy: Path | None = None) -> None:
        self.parity_copy = parity_copy
        self.passed: list[str] = []
        self.failed: list[dict[str, str]] = []

    def pytest_collection_modifyitems(self, items: list[pytest.Item]) -> None:
        if self.parity_copy is None:
            return
        modules = {item.module for item in items if isinstance(item, pytest.Function)}
        if len(modules) != 1:
            raise ValueError("Expected one parity test module")
        tests = modules.pop()
        # The fixture commits SCRIPT's bytes into its synthetic producer Git tree.
        # Both the executed helper and that bound source must be the same copy.
        tests.SCRIPT = self.parity_copy
        tests.HELPER = load_copy("issue1175_parity_mutant", self.parity_copy)

    def pytest_runtest_logreport(self, report: pytest.TestReport) -> None:
        if report.failed:
            self.failed.append(
                {"nodeid": report.nodeid, "phase": report.when, "detail": report.longreprtext}
            )
        elif report.when == "call" and report.passed:
            self.passed.append(report.nodeid)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("guard", choices=tuple(CASES))
    parser.add_argument("mode", choices=("mutant", "baseline"))
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if REPO_ROOT.resolve() != ROOT:
        raise ValueError("Imported mtp belongs to a different checkout")
    output = require_durable_output_path(args.output_dir.expanduser(), allowed=False)
    if output.is_relative_to(ROOT) or ROOT.is_relative_to(output):
        raise ValueError("Evidence output must be disjoint from the checkout")
    if args.output_dir.is_symlink():
        raise ValueError("Evidence output must not be a symlink")
    output.mkdir(parents=True, exist_ok=False)
    source, old, new, test_file, test_name, expected_count = CASES[args.guard]
    bound = (
        source,
        test_file,
        "src/mtp/matrix_overlap_coverage.py",
        "src/mtp/workflow/output_roots.py",
        "pyproject.toml",
        "uv.lock",
        str(Path(__file__).relative_to(ROOT)),
    )
    before = {name: file_sha256(ROOT / name) for name in bound}
    isolated = None
    if args.mode == "mutant":
        text = (ROOT / source).read_text()
        if text.count(old) != 1:
            raise ValueError("Expected exactly one guard mutation anchor")
        isolated = output / Path(source).name
        isolated.write_text(text.replace(old, new, 1))
        if args.guard == "orphan":
            module = load_copy("mtp.matrix_overlap_inputs", isolated)
            vars(importlib.import_module("mtp"))["matrix_overlap_inputs"] = module
    reports = Reports(isolated if args.guard == "provenance" else None)
    captured = io.StringIO()
    with redirect_stdout(captured), redirect_stderr(captured):
        result = int(
            pytest.main(
                ["-q", "-p", "no:cacheprovider", f"{ROOT / test_file}::{test_name}"],
                plugins=[reports],
            )
        )
    (output / "pytest.log").write_text(captured.getvalue())
    after = {name: file_sha256(ROOT / name) for name in bound}
    if args.mode == "baseline":
        intended = result == 0 and len(reports.passed) == expected_count and not reports.failed
    else:
        intended = (
            result == 1
            and not reports.passed
            and len(reports.failed) == expected_count
            and all(
                test_name in failure["nodeid"]
                and failure["phase"] == "call"
                and "DID NOT RAISE" in failure["detail"]
                for failure in reports.failed
            )
        )
    verified = intended and before == after
    evidence = {
        "issue": 1175,
        "guard": args.guard,
        "mode": args.mode,
        "pytest_exit": result,
        "passed": reports.passed,
        "failed": reports.failed,
        "isolated_source": str(isolated) if isolated else None,
        "isolated_sha256": file_sha256(isolated) if isolated else None,
        "canonical_before": before,
        "canonical_after": after,
        "verified": verified,
    }
    (output / "evidence.json").write_text(json.dumps(evidence, indent=2) + "\n")
    print(
        json.dumps(
            {"guard": args.guard, "mode": args.mode, "pytest_exit": result, "verified": verified}
        )
    )
    return int(not verified)


if __name__ == "__main__":
    raise SystemExit(main())
