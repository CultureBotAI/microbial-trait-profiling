"""Verify #1275 taxonomy guards using scratch pytest and JUnit/JSON receipts.

Inputs: coverage source, src packages, recorded fixture and dedicated tests.
Outputs: persistent verdicts.json,
JUnit XML and pytest logs under --scratch-dir. Never edit the live worktree.
Date: 2026-11-11. Run --canary before the full set.
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import html
import io
import json
import multiprocessing
import os
import re
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Any

import pytest

ROOT = Path(__file__).resolve().parents[2]
DOC = "scripts/create_taxa_coverage_histogram.py"
TEST = "tests/test_rank_class_only_taxonomy_1275.py"
FIXTURE = "tests/fixtures/coverage/taxonomy_rank_class.tsv"
REJECT = "test_rank_class_only_is_rejected[species]"
MIXED = "test_rank_classes_are_counted_and_excluded_with_organism_edge"
MUTANTS = (
    (
        "allow_rank_only",
        "        if not accepted:\n",
        "        if not accepted and not rank_class_records:\n",
        REJECT,
    ),
    (
        "wrong_error",
        "Taxonomy input has no NCBITaxon subclass_of edges",
        "Taxonomy contains no usable edges",
        "test_rank_class_only_is_rejected[fixture]",
    ),
    (
        "omit_skip_count",
        "                    rank_class_records += 1\n",
        "                    rank_class_records += 0\n",
        MIXED,
    ),
    (
        "rank_is_accepted",
        "                    rank_class_records += 1\n",
        "                    rank_class_records += 1\n                    accepted += 1\n",
        REJECT,
    ),
    (
        "rank_is_not_skipped",
        "                    rank_class_records += 1\n                    continue\n",
        "                    rank_class_records += 1\n",
        MIXED,
    ),
    (
        "rank_pollutes_nodes",
        "                    rank_class_records += 1\n",
        "                    rank_class_records += 1\n"
        "                    self.nodes.add(subject)\n",
        MIXED,
    ),
)


class Receipt:
    """Persist pytest's actual phase, outcome and exception site independently of XML."""

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
        self.path.write_text(json.dumps(self.records, indent=2) + "\n")


def worker(root: str, node: str) -> None:
    scratch = Path(root)
    os.chdir(scratch)
    os.environ["MPLBACKEND"] = "Agg"
    os.environ["PYTHONPATH"] = str(scratch / "src")
    sys.path.insert(0, str(scratch))
    sys.path.insert(0, str(scratch / "src"))
    for module_name in list(sys.modules):
        if module_name in {"mtp", "bsdbng"} or module_name.startswith(("mtp.", "bsdbng.")):
            del sys.modules[module_name]
    output = io.StringIO()
    with contextlib.redirect_stdout(output), contextlib.redirect_stderr(output):
        rc = int(
            pytest.main(
                [
                    "-q",
                    "-o",
                    "addopts=",
                    "-p",
                    "no:cacheprovider",
                    f"--junitxml={scratch / 'junit.xml'}",
                    f"{TEST}::{node}",
                ],
                plugins=[Receipt(scratch / "receipt.json")],
            )
        )
    (scratch / "pytest.log").write_text(output.getvalue())
    (scratch / "exitcode.txt").write_text(str(rc))


def run_case(run: Path, name: str, node: str, content: str, mutant: bool) -> str:
    scratch = run / name
    (scratch / "scripts").mkdir(parents=True)
    (scratch / "tests").mkdir()
    shutil.copy2(ROOT / TEST, scratch / TEST)
    (scratch / FIXTURE).parent.mkdir(parents=True)
    shutil.copy2(ROOT / FIXTURE, scratch / FIXTURE)
    shutil.copytree(ROOT / "src", scratch / "src", ignore=shutil.ignore_patterns("__pycache__"))
    (scratch / DOC).write_text(content)
    proc = multiprocessing.get_context("spawn").Process(target=worker, args=(str(scratch), node))
    proc.start()
    proc.join(60)
    if proc.is_alive():
        proc.terminate()
        proc.join()
        raise RuntimeError(f"Timeout: {name}")
    if proc.exitcode != 0:
        raise RuntimeError(f"Worker failed: {name}")
    rc = int((scratch / "exitcode.txt").read_text())
    # This deliberately reads only pytest's small, local JUnit dialect; no XML
    # parser, entity expansion or external resource resolution is involved.
    xml = (scratch / "junit.xml").read_text()
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
    records = json.loads((scratch / "receipt.json").read_text())
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
    paths = [DOC, TEST, FIXTURE, str(Path(__file__).resolve().relative_to(ROOT))]
    before = {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in paths}
    original = (ROOT / DOC).read_text()
    selected = MUTANTS[:1] if args.canary else MUTANTS
    for name, old, _new, _node in selected:
        if original.count(old) != 1:
            raise RuntimeError(f"Nonunique mutation anchor: {name}")
    verdicts = []
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
                {"verdicts": verdicts, "canonical_unchanged": before == after, "sha256": before},
                indent=2,
            )
            + "\n"
        )
        if before != after:
            raise RuntimeError("Canonical files changed")
    print(f"PASS {len(verdicts)}/{len(selected)}; evidence: {run}")


if __name__ == "__main__":
    main()
