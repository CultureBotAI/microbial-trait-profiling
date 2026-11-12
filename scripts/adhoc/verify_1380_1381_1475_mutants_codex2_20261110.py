"""Verify #1380/#1381/#1475 in isolated copies; never edit canonical files.

Inputs: agent_docs.py and the dedicated regression module. Outputs: scratch
pytest logs and verdicts.json with native call-phase assertion details.
Date: 2026-11-10. Run --canary before the full run with --scratch-dir.
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import io
import json
import multiprocessing
import os
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Any

import pytest

ROOT = Path(__file__).resolve().parents[2]
SOURCE = "src/mtp/agent_docs.py"
TEST = "tests/test_agent_docs_1380_1381_1475.py"
CONTROL_GUARD = """        if not angle:
            candidate = text[position + 1] if char == "\\\\" and position + 1 < limit else char
            if (ord(candidate) < 0x20 and candidate not in " \\t\\r\\n") or ord(candidate) == 0x7F:
                return None
"""
MUTANTS = (
    (
        "cr_reference",
        'rf"(?:^|(?<=[\\r\\n]))[ \\t]*{_REF_LABEL}"',
        'rf"(?m)^[ \\t]*{_REF_LABEL}"',
        "test_reference_line_endings",
    ),
    ("bare_controls", CONTROL_GUARD, "", "test_bare_destination_controls"),
    (
        "title_limit",
        "title_end = _title_end(text, position, limit)",
        "title_end = _title_end(text, position)",
        "test_reference_title_limit",
    ),
    (
        "destination_limit",
        "dest = _destination(content, dest_start, limit)",
        "dest = _destination(content, dest_start)",
        "test_outside_code_destination_limit",
    ),
)


class Events:
    """Record actual pytest call exceptions, rather than searching printed text."""

    def __init__(self) -> None:
        self.records: list[dict[str, Any]] = []

    @pytest.hookimpl(hookwrapper=True)
    def pytest_runtest_makereport(self, item: Any, call: Any) -> Any:
        outcome = yield
        report = outcome.get_result()
        error = call.excinfo
        self.records.append(
            {
                "nodeid": item.nodeid,
                "phase": report.when,
                "outcome": report.outcome,
                "assertion": error is not None and error.type is AssertionError,
                "detail": str(error) if error is not None else "passed",
                "site": str(error.traceback[-1].path) if error is not None else "",
            }
        )


def worker(root: str, test_name: str) -> None:
    scratch = Path(root)
    os.chdir(scratch)
    os.environ["PYTHONPATH"] = str(scratch / "src")
    sys.path.insert(0, str(scratch / "src"))
    for name in list(sys.modules):
        if name == "mtp" or name.startswith("mtp."):
            del sys.modules[name]
    from mtp import agent_docs

    if Path(agent_docs.__file__).resolve() != scratch / SOURCE:
        raise RuntimeError("Wrong source origin")
    events = Events()
    output = io.StringIO()
    with contextlib.redirect_stdout(output), contextlib.redirect_stderr(output):
        rc = int(
            pytest.main(
                ["-q", "-o", "addopts=", "-p", "no:cacheprovider", f"{TEST}::{test_name}"],
                plugins=[events],
            )
        )
    (scratch / "pytest.log").write_text(output.getvalue())
    (scratch / "result.json").write_text(json.dumps({"rc": rc, "events": events.records}, indent=2))


def run_case(run: Path, name: str, node: str, source: str) -> dict[str, Any]:
    scratch = run / name
    (scratch / "src/mtp").mkdir(parents=True)
    (scratch / "tests").mkdir()
    shutil.copy2(ROOT / "src/mtp/__init__.py", scratch / "src/mtp/__init__.py")
    shutil.copy2(ROOT / TEST, scratch / TEST)
    (scratch / SOURCE).write_text(source)
    proc = multiprocessing.get_context("spawn").Process(target=worker, args=(str(scratch), node))
    proc.start()
    proc.join(60)
    if proc.is_alive():
        proc.terminate()
        proc.join()
        raise RuntimeError(f"Timeout: {name}")
    if proc.exitcode != 0:
        raise RuntimeError(f"Worker failed: {name}, {proc.exitcode}")
    return dict(json.loads((scratch / "result.json").read_text()))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scratch-dir", type=Path, required=True)
    parser.add_argument("--canary", action="store_true")
    args = parser.parse_args()
    parent = args.scratch_dir.resolve()
    if parent == ROOT or parent.is_relative_to(ROOT):
        raise ValueError("Scratch must be outside the worktree")
    parent.mkdir(parents=True, exist_ok=True)
    run = Path(tempfile.mkdtemp(prefix="canary-" if args.canary else "full-", dir=parent))
    files = [SOURCE, TEST, str(Path(__file__).resolve().relative_to(ROOT))]
    before = {path: hashlib.sha256((ROOT / path).read_bytes()).hexdigest() for path in files}
    source = (ROOT / SOURCE).read_text()
    selected = MUTANTS[:1] if args.canary else MUTANTS
    for name, old, _new, _node in selected:
        if source.count(old) != 1:
            raise RuntimeError(f"Nonunique mutation anchor: {name}")
    verdicts = []
    try:
        for name, old, new, node in selected:
            control = run_case(run, name + "-control", node, source)
            if control["rc"] != 0 or not any(e["phase"] == "call" for e in control["events"]):
                raise RuntimeError(f"Control failed: {name}")
            result = run_case(run, name, node, source.replace(old, new, 1))
            failures = [e for e in result["events"] if e["outcome"] == "failed"]
            if (
                result["rc"] != 1
                or not failures
                or not all(
                    e["phase"] == "call"
                    and e["assertion"]
                    and e["nodeid"].startswith(f"{TEST}::{node}")
                    and Path(e["site"]).resolve() == (run / name / TEST).resolve()
                    and e["detail"].strip()
                    for e in failures
                )
            ):
                raise RuntimeError(f"Not an intended AssertionError kill: {name}: {result}")
            detail = "\n".join(e["detail"] for e in failures)
            verdicts.append(
                {
                    "mutant": name,
                    "killed": True,
                    "detail": detail,
                    "control_rc": control["rc"],
                    "rc": result["rc"],
                }
            )
            print(f"KILLED {name}: call-phase AssertionError ({len(failures)} cases)", flush=True)
    finally:
        after = {path: hashlib.sha256((ROOT / path).read_bytes()).hexdigest() for path in files}
        if before != after:
            raise RuntimeError("Canonical files changed")
        (run / "verdicts.json").write_text(
            json.dumps(
                {"verdicts": verdicts, "canonical_unchanged": before == after, "sha256": before},
                indent=2,
            )
            + "\n"
        )
    print(f"PASS {len(verdicts)}/{len(selected)}; evidence: {run}")


if __name__ == "__main__":
    main()
