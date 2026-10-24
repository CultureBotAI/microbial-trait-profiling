"""Verify #928 guard-test isolation; hardened for #950 (2026-10-23).

Inputs: this checkout's producer, synthetic guard tests and existing output tree.
Outputs: pytest.log and isolation.json in --output-dir. An isolated source copy
loses the guard; canonical source and checkout outputs are never rewritten.
"""

from __future__ import annotations

import argparse
import importlib.util
import io
import json
import sys
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from types import ModuleType

import pytest

from mtp.paths import REPO_ROOT, file_sha256

ROOT = Path(__file__).resolve().parents[2]
PRODUCER = ROOT / "scripts/create_trait_profile_infographic.py"
TEST = ROOT / "tests/test_trait_profile_infographic.py"
GUARD_TEST = "test_min_stratified_pairs_is_rejected_outside_site_stats"


def inventory(root: Path) -> dict[str, str]:
    """Include ignored and untracked files as well as tracked checkout outputs."""
    return {
        path.relative_to(root).as_posix(): file_sha256(path)
        for path in sorted(root.rglob("*"))
        if path.is_file()
    }


def load(path: Path, name: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ValueError("Cannot load producer source")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


class Probe:
    """Substitute only the isolated producer and distinguish assertions from crashes."""

    def __init__(self, producer: ModuleType) -> None:
        self.producer = producer
        self.failed: list[str] = []
        self.unexpected: list[str] = []

    def pytest_collection_modifyitems(self, items: list[pytest.Item]) -> None:
        for item in items:
            if isinstance(item, pytest.Function) and Path(item.module.__file__).resolve() == TEST:
                vars(item.module)["infographic"] = self.producer

    def pytest_runtest_logreport(self, report: pytest.TestReport) -> None:
        if report.failed:
            if (
                report.when == "call"
                and GUARD_TEST in report.nodeid
                and "DID NOT RAISE" in report.longreprtext
            ):
                self.failed.append(report.nodeid)
            else:
                self.unexpected.append(report.nodeid)

    def pytest_collectreport(self, report: pytest.CollectReport) -> None:
        if report.failed:
            self.unexpected.append(report.nodeid)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if REPO_ROOT.resolve() != ROOT:
        raise ValueError(
            "Installed mtp belongs to another checkout; use this checkout's environment"
        )
    canonical = load(PRODUCER, "issue950_canonical_infographic")
    checkout_outputs = vars(canonical)["INVESTIGATION_DIR"]
    if not isinstance(checkout_outputs, Path):
        raise ValueError("The producer must declare a Path output directory")
    output = args.output_dir.resolve()
    if output.is_relative_to(checkout_outputs):
        raise ValueError("Probe evidence must be outside the checkout's investigation outputs")
    before = inventory(checkout_outputs)
    canonical_paths = [PRODUCER, TEST, Path(__file__).resolve()]
    source_before = {str(path.relative_to(ROOT)): file_sha256(path) for path in canonical_paths}
    guard = (
        "        if not args.stratified_stats_only:\n"
        '            parser.error("--min-stratified-pairs requires --stratified-stats-only")\n'
    )
    text = PRODUCER.read_text()
    if text.count(guard) != 1:
        raise ValueError("Guard target changed; inspect the producer before mutation")
    output.mkdir(parents=True, exist_ok=False)
    copied = output / "infographic_without_guard.py"
    copied.write_text(text.replace(guard, ""))
    mutated = load(copied, "issue950_isolated_infographic")
    # Relocation must not change default paths. The shared fixture redirects these
    # exact declared outputs, including future constants, before main() executes.
    for name, value in vars(canonical).items():
        if isinstance(value, Path):
            vars(mutated)[name] = value
    probe = Probe(mutated)
    captured = io.StringIO()
    with redirect_stdout(captured), redirect_stderr(captured):
        code = int(
            pytest.main(["-q", "-p", "no:cacheprovider", f"{TEST}::{GUARD_TEST}"], plugins=[probe])
        )
    after = inventory(checkout_outputs)
    source_after = {str(path.relative_to(ROOT)): file_sha256(path) for path in canonical_paths}
    verified = (
        code == 1
        and len(probe.failed) == 2
        and not probe.unexpected
        and before == after
        and source_before == source_after
    )
    record = {
        "pytest_exit": code,
        "intended_failed_calls": probe.failed,
        "unexpected_failures": probe.unexpected,
        "canonical_origin": str(PRODUCER),
        "isolated_origin": str(copied),
        "isolated_sha256": file_sha256(copied),
        "checkout_outputs": checkout_outputs.relative_to(ROOT).as_posix(),
        "before": before,
        "after": after,
        "source_before": source_before,
        "source_after": source_after,
        "verified": verified,
    }
    (output / "pytest.log").write_text(captured.getvalue(), encoding="utf-8")
    (output / "isolation.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"verified": verified, "intended_failed_calls": probe.failed}))
    return int(not verified)


if __name__ == "__main__":
    raise SystemExit(main())
