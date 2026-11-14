"""Probe #959 coverage callout and scope regressions in isolated copies (#959, #966).

Inputs: this checkout's coverage producer and synthetic geometry/CLI tests.
Outputs: per-mutant pytest logs and canonical-immutability/expected-assertion evidence.
Each case loads an isolated producer copy; canonical tracked files are never modified.
"""

from __future__ import annotations

import argparse
import importlib
import io
import json
import sys
from collections import Counter
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from typing import Any

import pytest

import mtp
from mtp.paths import REPO_ROOT, file_sha256

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "scripts/create_taxa_coverage_histogram.py"
TEST = ROOT / "tests/test_taxa_coverage_histogram.py"
SCRIPT_ASSIGNMENT = (
    'SCRIPT = Path(__file__).resolve().parents[1] / "scripts/create_taxa_coverage_histogram.py"'
)
SELECTION = "narrow_segments or scope_is_explicitly"
MUTANTS = (
    (
        "dropped_label",
        "        text.remove()\n        callout = axis.text(",
        "        text.remove()\n        continue\n        callout = axis.text(",
        "Missing segment callout:",
    ),
    (
        "callout_on_bar",
        "            row + bar_height / 2 + 0.08,",
        "            row,",
        "Callout is on its bar:",
    ),
    (
        "no_leader",
        "            axis.annotate(\n",
        "            continue\n            axis.annotate(\n",
        "Missing segment leader:",
    ),
    (
        "scope_hardcoded_at_entrypoint",
        "    run_coverage(**vars(args))",
        '    arguments = vars(args)\n    arguments["strain_rollup"] = "off"\n'
        "    run_coverage(**arguments)",
        "AssertionError: assert 'off' == 'phenotype'",
    ),
)


def verify_checkout_identity() -> None:
    """Ensure imported mtp and target paths resolve to this helper's own checkout (#966)."""
    mtp_file = getattr(mtp, "__file__", None)
    if mtp_file is None or Path(mtp_file).resolve().parent != (ROOT / "src/mtp").resolve():
        raise AssertionError("Imported mtp belongs to a different checkout")
    if (REPO_ROOT / "src/mtp").resolve() != (ROOT / "src/mtp").resolve():
        raise AssertionError("Imported mtp.paths.REPO_ROOT belongs to a different checkout")
    if SOURCE.resolve().parent != (ROOT / "scripts").resolve():
        raise AssertionError("Coverage producer path is outside this checkout")
    if TEST.resolve().parent != (ROOT / "tests").resolve():
        raise AssertionError("Coverage test path is outside this checkout")


class _IsolatedCasePlugin:
    def __init__(self, isolated_source: Path, isolated_test: Path) -> None:
        self.isolated_source = isolated_source.resolve()
        self.isolated_test = isolated_test.resolve()
        self.loaded_module_origin: str | None = None
        self.counts: Counter[str] = Counter()
        self.failures: list[dict[str, str]] = []
        self.collection_errors: list[str] = []

    def pytest_collectreport(self, report: pytest.CollectReport) -> None:
        if report.failed:
            self.collection_errors.append(report.longreprtext)

    def pytest_collection_modifyitems(self, items: list[pytest.Item]) -> None:
        modules = {item.module for item in items if isinstance(item, pytest.Function)}
        if not modules:
            self.collection_errors.append("No test functions collected")
            return
        for module in modules:
            mod_file = getattr(module, "__file__", None)
            if mod_file is None or Path(mod_file).resolve() != self.isolated_test:
                raise AssertionError(f"Unexpected test module origin: {mod_file}")
            coverage = getattr(module, "coverage", None)
            cov_file = getattr(coverage, "__file__", None)
            if (
                coverage is None
                or cov_file is None
                or Path(cov_file).resolve() != self.isolated_source
            ):
                raise AssertionError(f"Unexpected coverage producer origin: {cov_file}")
            coverage.REPO_ROOT = ROOT
            self.loaded_module_origin = str(Path(cov_file).resolve())

    def pytest_runtest_logreport(self, report: pytest.TestReport) -> None:
        if report.when == "call":
            self.counts[report.outcome] += 1
        if report.failed:
            self.failures.append(
                {
                    "nodeid": report.nodeid,
                    "phase": report.when,
                    "failure": report.longreprtext,
                }
            )


def is_verified_case(
    *,
    exit_code: int,
    expected: str | None,
    log: str,
    failures: list[dict[str, str]],
    collection_errors: list[str],
) -> bool:
    """Require intended call-phase assertion and exit code; reject collection/other errors."""
    if collection_errors:
        return False
    if expected is None:
        return exit_code == 0 and not failures
    return (
        exit_code == 1
        and bool(failures)
        and expected in log
        and "AssertionError" in log
        and all(
            item["phase"] == "call"
            and expected in item["failure"]
            and "AssertionError" in item["failure"]
            for item in failures
        )
    )


def _clear_coverage_modules() -> None:
    importlib.invalidate_caches()
    for module_name, module in sys.modules.copy().items():
        filename = getattr(module, "__file__", None)
        if module_name in {"test_taxa_coverage_histogram", "taxa_coverage_histogram"} or (
            filename
            and Path(filename).name
            in {"test_taxa_coverage_histogram.py", "create_taxa_coverage_histogram.py"}
        ):
            sys.modules.pop(module_name, None)


def run_case(
    output: Path,
    name: str,
    expected: str | None,
    *,
    source_text: str | None = None,
) -> dict[str, object]:
    """Execute one probe case against an isolated producer copy without touching tracked files."""
    verify_checkout_identity()
    case_dir = (output / f"isolated-{name}").resolve()
    case_dir.mkdir(parents=True, exist_ok=False)
    isolated_source = case_dir / "create_taxa_coverage_histogram.py"
    isolated_test = case_dir / "test_taxa_coverage_histogram.py"
    producer_text = SOURCE.read_text(encoding="utf-8") if source_text is None else source_text
    isolated_source.write_text(producer_text, encoding="utf-8")

    test_text = TEST.read_text(encoding="utf-8")
    if test_text.count(SCRIPT_ASSIGNMENT) != 1:
        raise ValueError(
            "Expected exactly one SCRIPT assignment in test_taxa_coverage_histogram.py"
        )
    isolated_test.write_text(
        test_text.replace(SCRIPT_ASSIGNMENT, f"SCRIPT = Path({str(isolated_source)!r})"),
        encoding="utf-8",
    )

    _clear_coverage_modules()
    plugin = _IsolatedCasePlugin(isolated_source, isolated_test)
    captured = io.StringIO()
    try:
        with redirect_stdout(captured), redirect_stderr(captured):
            code = int(
                pytest.main(
                    [
                        "-q",
                        "-p",
                        "no:cacheprovider",
                        "--basetemp",
                        str(output / name),
                        str(isolated_test),
                        "-k",
                        SELECTION,
                    ],
                    plugins=[plugin],
                )
            )
    finally:
        _clear_coverage_modules()
    log = captured.getvalue()
    (output / f"{name}.log").write_text(log, encoding="utf-8")
    verified = plugin.loaded_module_origin == str(isolated_source.resolve()) and is_verified_case(
        exit_code=code,
        expected=expected,
        log=log,
        failures=plugin.failures,
        collection_errors=plugin.collection_errors,
    )
    return {
        "case": name,
        "exit_code": code,
        "expected_assertion": expected,
        "isolated_source_sha256": file_sha256(isolated_source),
        "loaded_module_origin": plugin.loaded_module_origin,
        "verified": verified,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    verify_checkout_identity()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    before_source = file_sha256(SOURCE)
    before_test = file_sha256(TEST)
    source = SOURCE.read_text(encoding="utf-8")
    results: list[dict[str, object]] = []
    evidence: dict[str, Any] = {}
    try:
        for name, old, _new, _expected in MUTANTS:
            if source.count(old) != 1:
                raise ValueError(f"Expected exactly one mutation target for {name}")
        results.append(run_case(output, "baseline", None, source_text=source))
        for name, old, new, expected in MUTANTS:
            mutated_text = source.replace(old, new, 1)
            results.append(run_case(output, name, expected, source_text=mutated_text))
    finally:
        after_source = file_sha256(SOURCE)
        after_test = file_sha256(TEST)
        canonical_unchanged = after_source == before_source and after_test == before_test
        evidence = {
            "results": results,
            "source_sha256_before": before_source,
            "source_sha256_after": after_source,
            "test_sha256_before": before_test,
            "test_sha256_after": after_test,
            "canonical_unchanged": canonical_unchanged,
        }
        (output / "evidence.json").write_text(
            json.dumps(evidence, indent=2) + "\n", encoding="utf-8"
        )
        if not canonical_unchanged:
            raise AssertionError("Canonical source or test changed during isolated probe")
    print(json.dumps(evidence, indent=2))
    return int(
        len(results) != len(MUTANTS) + 1
        or not all(bool(result["verified"]) for result in results)
        or not bool(evidence["canonical_unchanged"])
    )


if __name__ == "__main__":
    raise SystemExit(main())
