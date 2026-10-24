"""Probe #839/#929/#958 layout regressions against real rendered tests (2026-10-23).

Inputs: this checkout's guild renderer and regression test, synthetic profiles only.
Outputs: baseline/mutant diagnostic PNGs and JSON verdicts in --output-dir.
Constants are changed in memory; source files and scientific outputs are untouched.
"""

from __future__ import annotations

import argparse
import io
import json
import sys
from collections.abc import Iterator
from contextlib import contextmanager, redirect_stderr, redirect_stdout
from pathlib import Path
from typing import TypedDict

import pytest

from mtp import guild_maps

TEST_FUNCTIONS = (
    "tests/test_trait_guilds.py::test_the_plotting_area_starts_below_the_last_subtitle_line",
    "tests/test_table1_trait_maps.py::test_the_site_names_sit_below_the_last_subtitle_line",
)
EXPECTED_NODES = tuple(f"{test}[{value}]" for test in TEST_FUNCTIONS for value in (False, True))
# The guild layout without site headers must survive a site-header-only mutation.
# The other weakenings affect all four layouts; pin each observed regression case.
EXPECTED_FAILURES = {
    "baseline": (),
    "TITLE_BLOCK_HEIGHT": EXPECTED_NODES,
    "SUBTITLE_LINE_HEIGHT": EXPECTED_NODES,
    "SITE_HEADER_HEIGHT": EXPECTED_NODES[1:],
    "entry_counting": EXPECTED_NODES,
}
GEOMETRY_MESSAGES = (
    "Title/header clearance below ",
    "A two-line note changed title/header clearance by more than 4 layout px:",
)


@contextmanager
def fresh_test_modules() -> Iterator[None]:
    """Reimport tests for each mutation, then restore the caller's module bindings."""
    names = ("test_trait_guilds", "test_table1_trait_maps")
    previous = {name: sys.modules[name] for name in names if name in sys.modules}
    for name in names:
        sys.modules.pop(name, None)
    try:
        yield
    finally:
        for name in names:
            sys.modules.pop(name, None)
        sys.modules.update(previous)


class CaseResult(TypedDict):
    """One actual parameterized test, including failures outside the call phase."""

    nodeid: str
    status: str
    phase: str
    message: str


class CaseRecorder:
    """Reject crashes and unrelated assertions rather than calling them mutation kills."""

    def __init__(self) -> None:
        self.cases: dict[str, CaseResult] = {}
        self.collection_errors: list[str] = []

    def pytest_collectreport(self, report: pytest.CollectReport) -> None:
        if report.failed:
            self.collection_errors.append(report.longreprtext)

    def pytest_runtest_makereport(self, item: pytest.Item, call: pytest.CallInfo[object]) -> None:
        previous = self.cases.get(item.nodeid)
        if call.excinfo is not None:
            error = call.excinfo.value
            message = str(error)
            intended = (
                call.when == "call"
                and isinstance(error, AssertionError)
                and message.startswith(GEOMETRY_MESSAGES)
            )
            # An earlier setup/call error must not be overwritten by a later phase.
            if previous is None or previous["status"] != "error":
                self.cases[item.nodeid] = {
                    "nodeid": item.nodeid,
                    "status": "assertion_failure" if intended else "error",
                    "phase": call.when,
                    "message": f"{type(error).__name__}: {message}",
                }
        elif call.when == "call" and previous is None:
            self.cases[item.nodeid] = {
                "nodeid": item.nodeid,
                "status": "passed",
                "phase": "call",
                "message": "",
            }

    def results(self) -> list[CaseResult]:
        return [
            self.cases.get(node, {"nodeid": node, "status": "not_run", "phase": "", "message": ""})
            for node in EXPECTED_NODES
        ]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    root = Path(__file__).resolve().parents[2]
    outcomes: list[dict[str, object]] = []
    for constant in (
        None,
        "TITLE_BLOCK_HEIGHT",
        "SUBTITLE_LINE_HEIGHT",
        "SITE_HEADER_HEIGHT",
        "entry_counting",
    ):
        name = constant or "baseline"
        destination = args.output_dir / name
        destination.mkdir(exist_ok=True)
        recorder = CaseRecorder()
        with pytest.MonkeyPatch.context() as patch:
            if constant == "entry_counting":
                patch.setattr(
                    guild_maps,
                    "_title_block_height",
                    lambda lines: (
                        guild_maps.TITLE_BLOCK_HEIGHT + guild_maps.SUBTITLE_LINE_HEIGHT * len(lines)
                    ),
                )
            elif constant:
                patch.setattr(guild_maps, constant, getattr(guild_maps, constant) // 2)
            # Import each test after mutation, just as a fresh checkout run would.
            # Retained imported constants would spuriously reject a formula-only test.
            captured = io.StringIO()
            with fresh_test_modules(), redirect_stdout(captured), redirect_stderr(captured):
                code = pytest.main(
                    [
                        "-q",
                        "-p",
                        "no:cacheprovider",
                        "--basetemp",
                        str(destination / "fixtures"),
                        *(str(root / test) for test in TEST_FUNCTIONS),
                    ],
                    plugins=[recorder],
                )
            destination.joinpath("pytest.log").write_text(captured.getvalue(), encoding="utf-8")
            cases = recorder.results()
            failures = {case["nodeid"] for case in cases if case["status"] == "assertion_failure"}
            matched_expected_failures = failures == set(EXPECTED_FAILURES[name])
            valid = (
                not recorder.collection_errors
                and set(recorder.cases) == set(EXPECTED_NODES)
                and all(case["status"] in {"passed", "assertion_failure"} for case in cases)
            )
            rejected = (
                valid
                and code == pytest.ExitCode.TESTS_FAILED
                and bool(failures)
                and matched_expected_failures
            )
            passed = (
                valid
                and code == pytest.ExitCode.OK
                and all(case["status"] == "passed" for case in cases)
            )
            outcomes.append(
                {
                    "case": name,
                    "rejected": rejected,
                    "passed": passed,
                    "pytest_exit": int(code),
                    "tests": cases,
                    "collection_errors": recorder.collection_errors,
                    "expected_failures": EXPECTED_FAILURES[name],
                    "matched_expected_failures": matched_expected_failures,
                }
            )
    args.output_dir.joinpath("verdicts.json").write_text(
        json.dumps(outcomes, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(outcomes, indent=2))
    return 0 if outcomes[0]["passed"] and all(row["rejected"] for row in outcomes[1:]) else 1


if __name__ == "__main__":
    raise SystemExit(main())
