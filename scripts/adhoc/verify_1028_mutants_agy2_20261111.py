"""Verify #1028 docs/workflow.md site-stats provenance schema 1.1/1.2 guards.

Purpose: Copy docs/workflow.md, stratified_stats_provenance.json, and
tests/test_workflow_site_stats_provenance_1028.py into isolated scratch
directories under --scratch-dir and verify that each of 8 distinct mutations
to docs/workflow.md is killed by a call-phase AssertionError on the target
pytest node, verified independently via JUnit XML and a pytest hook receipt.
Inputs: docs/workflow.md,
trait-visualization/trait-profile-investigation/stratified_stats_provenance.json,
tests/test_workflow_site_stats_provenance_1028.py.
Outputs: verdicts.json, receipt.json, junit.xml, and pytest.log under --scratch-dir.
Refs: #1028
Date: 2026-11-11
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import html
import importlib
import io
import json
import os
import re
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Any

import pytest

ROOT = Path(__file__).resolve().parents[2]
DOC = "docs/workflow.md"
PROVENANCE = "trait-visualization/trait-profile-investigation/stratified_stats_provenance.json"
TEST = "tests/test_workflow_site_stats_provenance_1028.py"
TEST_MODULE = "test_workflow_site_stats_provenance_1028"
NODE = "test_workflow_documents_committed_schema_1_1_and_schema_1_2_site_stats_provenance"

MUTANTS: tuple[tuple[str, str, str, str], ...] = (
    (
        "committed_schema_1_2",
        "written by the final rerun (#1004) is schema version 1.1,",
        "written by the final rerun (#1004) is schema version 1.2,",
        NODE,
    ),
    (
        "wrong_entry_count",
        "with all eight of its\n`not_regenerated` entries",
        "with all seven of its\n`not_regenerated` entries",
        NODE,
    ),
    (
        "wrong_legacy_reason",
        '`"not refreshed by pooled/site statistics mode"`',
        '`"existing_stratum_table_not_refreshed"`',
        NODE,
    ),
    (
        "missing_rerun_issue",
        "written by the final rerun (#1004) is schema version 1.1,",
        "written by the final rerun is schema version 1.1,",
        NODE,
    ),
    (
        "missing_rewrite_clause",
        (
            "the next\n`just trait-profile-site-stats` run rewrites it compatibly "
            "under schema version 1.2,\nand nothing reads it."
        ),
        "nothing reads it.",
        NODE,
    ),
    (
        "wrong_rewrite_version",
        "rewrites it compatibly under schema version 1.2,",
        "rewrites it compatibly under schema version 1.1,",
        NODE,
    ),
    (
        "missing_existing_stratum_reason",
        "distinguishes `existing_stratum_table_not_refreshed` from",
        "distinguishes `legacy_stratum_table_not_refreshed` from",
        NODE,
    ),
    (
        "missing_not_produced_reason",
        "`not_produced_by_site_stats_mode` (the summaries and figure).",
        "`omitted_by_site_stats_mode` (the summaries and figure).",
        NODE,
    ),
)


class Receipt:
    """Persist pytest's phase, outcome, and exception site independently of XML."""

    def __init__(self, path: Path) -> None:
        self.path = path
        self.records: list[dict[str, Any]] = []

    @pytest.hookimpl(hookwrapper=True)
    def pytest_runtest_makereport(self, item: Any, call: Any) -> Any:
        outcome = yield
        report = outcome.get_result()
        error = call.excinfo
        self.records.append(
            {
                "nodeid": item.nodeid,
                "phase": call.when,
                "outcome": report.outcome,
                "assertion": error is not None and error.type is AssertionError,
                "site": str(error.traceback[-1].path) if error is not None else "",
            }
        )

    def pytest_sessionfinish(self) -> None:
        self.path.write_text(json.dumps(self.records, indent=2) + "\n", encoding="utf-8")


def _execute_scratch_pytest(scratch: Path, node: str) -> int:
    os.environ["MPLBACKEND"] = "Agg"
    module_keys = (TEST_MODULE, f"tests.{TEST_MODULE}")
    saved_modules = {k: sys.modules.pop(k) for k in module_keys if k in sys.modules}
    importlib.invalidate_caches()
    output = io.StringIO()
    try:
        with (
            contextlib.chdir(scratch),
            contextlib.redirect_stdout(output),
            contextlib.redirect_stderr(output),
        ):
            rc = int(
                pytest.main(
                    [
                        "-q",
                        "-o",
                        "addopts=",
                        "-p",
                        "no:cacheprovider",
                        "--import-mode=importlib",
                        f"--junitxml={scratch / 'junit.xml'}",
                        f"{TEST}::{node}",
                    ],
                    plugins=[Receipt(scratch / "receipt.json")],
                )
            )
    finally:
        for k in module_keys:
            sys.modules.pop(k, None)
        sys.modules.update(saved_modules)
    (scratch / "pytest.log").write_text(output.getvalue(), encoding="utf-8")
    (scratch / "exitcode.txt").write_text(str(rc), encoding="utf-8")
    return rc


def run_case(
    run: Path,
    name: str,
    node: str,
    content: str,
    mutant: bool,
    test_override: str | None = None,
) -> str:
    scratch = run / name
    (scratch / "docs").mkdir(parents=True)
    (scratch / "tests").mkdir()
    (scratch / PROVENANCE).parent.mkdir(parents=True)
    shutil.copy2(ROOT / PROVENANCE, scratch / PROVENANCE)
    if test_override is None:
        shutil.copy2(ROOT / TEST, scratch / TEST)
    else:
        (scratch / TEST).write_text(test_override, encoding="utf-8")
    (scratch / DOC).write_text(content, encoding="utf-8")
    rc = _execute_scratch_pytest(scratch, node)
    # Read only pytest's small, local JUnit dialect via stdlib regex/html; no XML
    # parser, entity expansion, or external resource resolution is used.
    xml = (scratch / "junit.xml").read_text(encoding="utf-8")
    cases = re.findall(r"<testcase\b([^>]*?)(?:/\s*>|>(.*?)</testcase>)", xml, re.DOTALL)
    if len(cases) != 1 or len(re.findall(r"<testcase\b", xml)) != 1:
        raise RuntimeError(f"Wrong JUnit case count: {name}")
    attributes = {
        key: html.unescape(value)
        for key, _quote, value in re.findall(r"""([\w:.-]+)\s*=\s*(["'])(.*?)\2""", cases[0][0])
    }
    if attributes.get("name") != node or attributes.get("classname") != TEST[:-3].replace("/", "."):
        raise RuntimeError(f"Wrong JUnit test identity: {name}")
    if re.search(r"<(?:error|skipped)\b", xml):
        raise RuntimeError(f"Setup/collection/teardown error or skip: {name}")
    failures = re.findall(r"<failure\b[^>]*>(.*?)</failure>", cases[0][1], re.DOTALL)
    if len(re.findall(r"<failure\b", xml)) != len(failures):
        raise RuntimeError(f"Malformed or misplaced JUnit failure: {name}")
    records = json.loads((scratch / "receipt.json").read_text(encoding="utf-8"))
    if (
        not isinstance(records, list)
        or len(records) != 3
        or [record.get("phase") for record in records] != ["setup", "call", "teardown"]
        or any(record.get("nodeid") != f"{TEST}::{node}" for record in records)
        or any(records[index].get("outcome") != "passed" for index in (0, 2))
    ):
        raise RuntimeError(f"Invalid pytest phase receipt: {name}")
    call = records[1]
    if not mutant:
        if rc != 0 or failures or call.get("outcome") != "passed" or call.get("assertion"):
            raise RuntimeError(f"Control did not pass: {name}")
        return "Target node passed in call phase (JUnit and JSON verified)"
    detail = html.unescape(failures[0]).strip() if len(failures) == 1 else ""
    if (
        rc != 1
        or not detail
        or call.get("outcome") != "failed"
        or call.get("assertion") is not True
        or Path(call.get("site") or "").resolve() != (scratch / TEST).resolve()
    ):
        raise RuntimeError(f"Not a target call-phase AssertionError: {name}")
    return detail


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scratch-dir", type=Path, required=True)
    parser.add_argument("--canary", action="store_true")
    args = parser.parse_args()
    parent = args.scratch_dir.resolve()
    if parent == ROOT or parent.is_relative_to(ROOT):
        raise ValueError("Scratch must be outside worktree")
    parent.mkdir(parents=True, exist_ok=True)
    run = Path(tempfile.mkdtemp(prefix="canary-" if args.canary else "full-", dir=parent))
    paths = [DOC, PROVENANCE, TEST, str(Path(__file__).resolve().relative_to(ROOT))]
    before = {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in paths}
    original = (ROOT / DOC).read_text(encoding="utf-8")
    selected = MUTANTS[:1] if args.canary else MUTANTS
    for name, old, _new, _node in selected:
        if original.count(old) != 1:
            raise RuntimeError(f"Nonunique mutation anchor: {name}")
    verdicts: list[dict[str, Any]] = []
    try:
        for name, old, new, node in selected:
            control = run_case(run, name + "-control", node, original, False)
            detail = run_case(run, name, node, original.replace(old, new, 1), True)
            verdicts.append(
                {
                    "mutant": name,
                    "node": f"{TEST}::{node}",
                    "killed": True,
                    "control": control,
                    "detail": detail,
                }
            )
            print(f"KILLED {name}: JUnit target call-phase AssertionError", flush=True)
    finally:
        after = {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in paths}
        (run / "verdicts.json").write_text(
            json.dumps(
                {
                    "verdicts": verdicts,
                    "canonical_unchanged": before == after,
                    "sha256": before,
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        if before != after:
            raise RuntimeError("Canonical files changed")
    print(f"PASS {len(verdicts)}/{len(selected)}; evidence: {run}")


if __name__ == "__main__":
    main()
