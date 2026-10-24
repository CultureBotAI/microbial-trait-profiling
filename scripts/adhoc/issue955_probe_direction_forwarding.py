"""Probe CLI direction propagation using isolated source copies (#955/#962; 2026-10-23).

Inputs: this checkout's producer, shared labels and real main()-level regression.
Outputs: fresh per-case source/test copies, pytest logs and identity evidence under
--output-dir. Each case runs in a separate process; canonical files are never edited.
No production inputs, outputs or images are read or regenerated.
"""

from __future__ import annotations

import argparse
import io
import json
import multiprocessing
import sys
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

import pytest

from mtp.paths import REPO_ROOT, file_sha256

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "scripts/create_study_trait_profile_heatmap.py"
TEST_PATH = ROOT / "tests/test_study_trait_profile_heatmap.py"
TEST_NAME = "test_main_preserves_direction_basis_in_figures_and_drawn_values"
FORWARDING = "direction_basis=args.direction_basis,"
EXPECTED = {
    "baseline": None,
    "profile": "CLI basis did not reach profile legend",
    "paired": "CLI basis did not reach paired colorbar",
    "metadata": "CLI basis did not reach paired metadata",
}


def validate_checkout() -> None:
    """Refuse another checkout before reading producer bytes or creating outputs."""
    if REPO_ROOT.resolve() != ROOT:
        raise ValueError("Imported mtp belongs to another checkout")
    paths_module = sys.modules[file_sha256.__module__]
    if paths_module.__file__ is None or Path(paths_module.__file__).resolve() != (
        ROOT / "src/mtp/paths.py"
    ):
        raise ValueError("Imported hash helper belongs to another checkout")


def canonical_hashes() -> dict[str, str]:
    paths = (
        SOURCE,
        TEST_PATH,
        ROOT / "src/mtp/direction_labels.py",
        ROOT / "scripts/adhoc/issue842_render_direction_fixtures.py",
        Path(__file__).resolve(),
    )
    return {str(path.relative_to(ROOT)): file_sha256(path) for path in paths}


class Outcomes:
    """Accept only the named regression's intended call assertion, never a setup error."""

    def __init__(self, expected: str | None) -> None:
        self.expected = expected
        self.passed: list[str] = []
        self.failed: list[str] = []
        self.unexpected: list[str] = []

    def pytest_runtest_logreport(self, report: pytest.TestReport) -> None:
        intended = report.nodeid.endswith(f"::{TEST_NAME}") and report.when == "call"
        if intended and report.passed:
            self.passed.append(report.nodeid)
        if report.failed:
            if (
                intended
                and self.expected
                and f"AssertionError: {self.expected}" in (report.longreprtext)
            ):
                self.failed.append(report.nodeid)
            else:
                self.unexpected.append(report.nodeid)

    def pytest_collectreport(self, report: pytest.CollectReport) -> None:
        if report.failed:
            self.unexpected.append(report.nodeid)


def isolated_source(case: str) -> str:
    source = SOURCE.read_text()
    if source.count(FORWARDING) != 2:
        raise ValueError("Expected exactly two main() direction-forwarding call sites")
    if case in {"profile", "paired"}:
        position = source.find(FORWARDING) if case == "profile" else source.rfind(FORWARDING)
        return (
            source[:position]
            + 'direction_basis="group_1_vs_group_0",'
            + source[position + len(FORWARDING) :]
        )
    if case == "metadata":
        anchor = '"direction_basis": args.direction_basis,'
        if source.count(anchor) != 1:
            raise ValueError("Expected exactly one metadata direction-forwarding site")
        return source.replace(anchor, '"direction_basis": "group_1_vs_group_0",')
    return source


def copied_test(producer: Path) -> str:
    """Change only import/fixture locations; retain the actual numerical assertions."""
    source = TEST_PATH.read_text()
    anchors = {
        'MODULE_PATH = (\n    Path(__file__).resolve().parents[1] / "scripts" / '
        '"create_study_trait_profile_heatmap.py"\n)': f"MODULE_PATH = Path({str(producer)!r})",
        'str(MODULE_PATH.parent / "adhoc/issue842_render_direction_fixtures.py")': (
            f"str(Path({str(ROOT / 'scripts/adhoc/issue842_render_direction_fixtures.py')!r}))"
        ),
        "SPEC.loader.exec_module(heatmap)": (
            "SPEC.loader.exec_module(heatmap)\n"
            "assert Path(heatmap.__file__).resolve() == MODULE_PATH.resolve()"
        ),
    }
    for anchor, replacement in anchors.items():
        if source.count(anchor) != 1:
            raise ValueError(f"Expected one test-loader anchor: {anchor}")
        source = source.replace(anchor, replacement)
    return source


def run_case(output: Path, case: str) -> dict[str, object]:
    producer = output / "profile_heatmap.py"
    producer.write_text(isolated_source(case), encoding="utf-8")
    test_path = output / "test_direction_forwarding.py"
    test_path.write_text(copied_test(producer), encoding="utf-8")
    expected = EXPECTED[case]
    outcomes = Outcomes(expected)
    captured = io.StringIO()
    with redirect_stdout(captured), redirect_stderr(captured):
        code = int(
            pytest.main(
                [
                    "-q",
                    "-p",
                    "no:cacheprovider",
                    "--basetemp",
                    str(output / "pytest"),
                    f"{test_path}::{TEST_NAME}",
                ],
                plugins=[outcomes],
            )
        )
    (output / "pytest.log").write_text(captured.getvalue(), encoding="utf-8")
    return {
        "case": case,
        "pytest_exit": code,
        "expected_assertion": expected,
        "passed_calls": outcomes.passed,
        "failed_calls": outcomes.failed,
        "unexpected_failures": outcomes.unexpected,
        "isolated_producer": {"path": str(producer), "sha256": file_sha256(producer)},
        "isolated_test": {"path": str(test_path), "sha256": file_sha256(test_path)},
        "verified": code == (0 if expected is None else 1)
        and len(outcomes.passed) == int(expected is None)
        and len(outcomes.failed) == int(expected is not None)
        and not outcomes.unexpected,
    }


def child_case(output: Path, case: str) -> None:
    """Run a fixed internal worker in a fresh Python process, without command execution."""
    validate_checkout()
    before = canonical_hashes()
    output.mkdir(parents=True, exist_ok=False)
    sys.dont_write_bytecode = True
    evidence = run_case(output, case)
    after = canonical_hashes()
    evidence.update(canonical_before=before, canonical_after=after)
    evidence["verified"] = bool(evidence["verified"]) and before == after
    (output / "evidence.json").write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8")
    raise SystemExit(int(not evidence["verified"]))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    validate_checkout()
    before = canonical_hashes()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    cases = {}
    for case in EXPECTED:
        process = multiprocessing.get_context("spawn").Process(
            target=child_case, args=(output / case, case), name=f"direction-probe-{case}"
        )
        process.start()
        process.join()
        cases[case] = {
            "process_exit": process.exitcode,
            "evidence": str(output / case / "evidence.json"),
        }
        process.close()
    evidence: dict[str, object] = {
        "cases": cases,
        "verified": all(v["process_exit"] == 0 for v in cases.values()),
    }
    after = canonical_hashes()
    evidence.update(canonical_before=before, canonical_after=after)
    evidence["verified"] = bool(evidence["verified"]) and before == after
    (output / "evidence.json").write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"verified": evidence["verified"], "output": str(output)}))
    return int(not evidence["verified"])


if __name__ == "__main__":
    raise SystemExit(main())
