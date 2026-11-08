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
import os
import sys
from collections.abc import Sequence
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from types import ModuleType

import pytest

from mtp.paths import REPO_ROOT, file_sha256
from mtp.provenance_paths import PortablePaths, assert_portable_path_text, assert_portable_text

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
    previous = sys.dont_write_bytecode
    sys.dont_write_bytecode = True
    try:
        spec.loader.exec_module(module)
    finally:
        sys.dont_write_bytecode = previous
    return module


class Probe:
    """Substitute only the isolated producer and distinguish assertions from crashes."""

    def __init__(self, producer: ModuleType) -> None:
        self.producer = producer
        self.failed: list[str] = []
        self.unexpected: list[str] = []
        self.temp_roots: set[Path] = set()

    def _record_temp_roots(self, config: pytest.Config) -> None:
        factory = getattr(config, "_tmp_path_factory", None)
        if factory is not None:
            basetemp = Path(factory.getbasetemp()).resolve()
            self.temp_roots.update((basetemp, basetemp.parent))

    def pytest_configure(self, config: pytest.Config) -> None:
        self._record_temp_roots(config)

    def pytest_runtest_setup(self, item: pytest.Item) -> None:
        self._record_temp_roots(item.config)

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


def _path_spellings(path: Path) -> set[str]:
    spellings = {path.as_posix(), path.resolve().as_posix()}
    for item in list(spellings):
        if item.startswith("/private/"):
            spellings.add(item.removeprefix("/private"))
    return {item for item in spellings if item and item != "/"}


def scrub_pytest_log(
    log_text: str,
    *,
    output: Path,
    repo_root: Path,
    temp_roots: set[Path],
) -> str:
    """Replace output, repository, and pytest temporary paths before portability checks."""
    replacements: list[tuple[str, str]] = []
    for spelling in _path_spellings(output):
        replacements.append((spelling, "<output-dir>"))
    for spelling in _path_spellings(repo_root):
        replacements.append((spelling, "<repo-root>"))
    cwd = Path.cwd().resolve()
    if cwd != repo_root.resolve():
        rel_repo = os.path.relpath(repo_root.resolve(), cwd).replace(os.sep, "/")
        if rel_repo and rel_repo != ".":
            replacements.append((rel_repo, "<repo-root>"))
    for temp_root in temp_roots:
        for spelling in _path_spellings(temp_root):
            replacements.append((spelling, "<pytest-tempdir>"))
    ordered = sorted(set(replacements), key=lambda item: len(item[0]), reverse=True)
    scrubbed = log_text
    for target, placeholder in ordered:
        scrubbed = scrubbed.replace(target, placeholder)
    if "..." in scrubbed:
        for target, placeholder in ordered:
            for index in range(1, max(1, len(target) - 5)):
                fragment = f"...{target[index:]}"
                if fragment in scrubbed:
                    scrubbed = scrubbed.replace(fragment, f".../{placeholder}")
    assert_portable_text(scrubbed)
    remainder = (
        scrubbed.replace("<output-dir>", "")
        .replace("<repo-root>", "")
        .replace("<pytest-tempdir>", "")
    )
    assert_portable_path_text(remainder)
    return scrubbed


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--provenance-root", type=Path, default=None)
    args = parser.parse_args(argv)
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
    assert_portable_path_text(output.name)
    copied = output / "infographic_without_guard.py"
    canonical_paths = [PRODUCER, TEST, Path(__file__).resolve()]
    repo_root = args.provenance_root.resolve() if args.provenance_root is not None else ROOT
    repo_paths = PortablePaths.create(
        output / "isolation.json",
        [*canonical_paths, copied] if args.provenance_root is not None else canonical_paths,
        root=repo_root,
    )
    before = inventory(checkout_outputs)
    source_before = {repo_paths.encode(path): file_sha256(path) for path in canonical_paths}
    guard = (
        "        if not args.stratified_stats_only:\n"
        '            parser.error("--min-stratified-pairs requires --stratified-stats-only")\n'
    )
    text = PRODUCER.read_text()
    if text.count(guard) != 1:
        raise ValueError("Guard target changed; inspect the producer before mutation")
    output.mkdir(parents=True, exist_ok=False)
    copied.write_text(text.replace(guard, ""))
    output_root = repo_root if output.is_relative_to(repo_root) else output
    output_paths = PortablePaths.create(output / "isolation.json", [copied], root=output_root)
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
    source_after = {repo_paths.encode(path): file_sha256(path) for path in canonical_paths}
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
        "path_basis": repo_paths.declaration,
        "isolated_path_basis": output_paths.declaration,
        "canonical_origin": repo_paths.encode(PRODUCER),
        "isolated_origin": output_paths.encode(copied),
        "isolated_sha256": file_sha256(copied),
        "checkout_outputs": checkout_outputs.relative_to(ROOT).as_posix(),
        "before": before,
        "after": after,
        "source_before": source_before,
        "source_after": source_after,
        "verified": verified,
    }
    isolation_text = json.dumps(record, indent=2) + "\n"
    assert_portable_text(isolation_text)
    log_text = scrub_pytest_log(
        captured.getvalue(),
        output=args.output_dir,
        repo_root=ROOT,
        temp_roots=probe.temp_roots,
    )
    (output / "pytest.log").write_text(log_text, encoding="utf-8")
    (output / "isolation.json").write_text(isolation_text, encoding="utf-8")
    print(json.dumps({"verified": verified, "intended_failed_calls": probe.failed}))
    return int(not verified)


if __name__ == "__main__":
    raise SystemExit(main())
