"""Prove buffered comparison guards using private module copies (#1215, 2026-10-29).

Inputs: frozen comparison/adapter modules and fourteen synthetic regression cases.
Outputs: fresh durable mutated copies, pytest logs and hash-bound JSON evidence.
Each case runs in a fresh process; shared production files are never rewritten.
The missing-inventory cases remain schema controls when the output loop is bypassed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing
import os
import platform
import sys
from collections import Counter
from contextlib import redirect_stderr, redirect_stdout
from importlib.machinery import SourceFileLoader
from importlib.metadata import version
from pathlib import Path
from typing import Any

import pytest
from git import Repo

ROOT = Path(__file__).resolve().parents[2]
COMPARISON = "src/mtp/workflow/comparison.py"
SOURCES = "src/mtp/multiplicity/sources.py"
TEST_FILE = "tests/test_comparison_run.py"
PATH_TEST = "test_buffered_comparison_rejects_wrong_manifest_path"
OUTPUT_TEST = "test_buffered_comparison_retains_output_inventory_checks"
ADAPTER_TEST = "test_buffered_adapter_rejects_one_core_population_or_direction_mismatch"
PATH_CASES = {f"{PATH_TEST}[{entry}]" for entry in ("comparison", "sources")}
OUTPUT_CASES = {
    f"{OUTPUT_TEST}[{mutation}-{entry}]"
    for mutation in ("empty_inventory", "missing_file", "hash", "size")
    for entry in ("comparison", "sources")
}
SCHEMA_CASES = {f"{OUTPUT_TEST}[missing_inventory-{entry}]" for entry in ("comparison", "sources")}
ADAPTER_CASES = {f"{ADAPTER_TEST}[{field}]" for field in ("population", "direction")}
TESTS = tuple(sorted(PATH_CASES | OUTPUT_CASES | SCHEMA_CASES | ADAPTER_CASES))
OUTPUT_BLOCK = """    if not record.outputs:
        raise ValueError("ComparisonRun has no recorded outputs")
    for entry in record.outputs:
        target = root / entry.path
        if not target.is_file():
            raise FileNotFoundError(f"Missing ComparisonRun output file: {target}")
        if target.stat().st_size != entry.size_bytes:
            raise ValueError(f"Size mismatch for ComparisonRun output file: {target}")
        if file_sha256(target) != entry.sha256:
            raise ValueError(f"SHA-256 mismatch for ComparisonRun output file: {target}")
"""
MUTATIONS = {
    "skip_buffered_manifest_path": (
        COMPARISON,
        "    if Path(record.manifest_path).resolve() != manifest:",
        "    if manifest_bytes is None and Path(record.manifest_path).resolve() != manifest:",
        PATH_CASES,
    ),
    "skip_buffered_output_inventory": (
        COMPARISON,
        OUTPUT_BLOCK,
        "    if manifest_bytes is None:\n"
        + "".join("    " + line for line in OUTPUT_BLOCK.splitlines(keepends=True)),
        OUTPUT_CASES,
    ),
    "skip_buffered_population_direction": (
        SOURCES,
        "    if (\n        comparison.contrast_population != core.contrast_population\n",
        "    if manifest_bytes is None and (\n"
        "        comparison.contrast_population != core.contrast_population\n",
        ADAPTER_CASES,
    ),
}
BOUND_FILES = (
    COMPARISON,
    SOURCES,
    TEST_FILE,
    "scripts/adhoc/issue1215_probe_buffer_guards.py",
    "scripts/adhoc/issue844_influence_assessment_fixture.py",
    "pyproject.toml",
    "uv.lock",
)


def digest(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


class Reports:
    def __init__(self) -> None:
        self.counts: Counter[str] = Counter()
        self.failures: list[dict[str, str]] = []

    def pytest_runtest_logreport(self, report: pytest.TestReport) -> None:
        if report.when == "call":
            self.counts[report.outcome] += 1
        if report.failed:
            self.failures.append(
                {"nodeid": report.nodeid, "phase": report.when, "failure": report.longreprtext}
            )


def worker(directory: Path, changed_file: str | None) -> None:
    reports = Reports()
    arguments = [
        "-q",
        "--tb=short",
        "--basetemp",
        str(directory / "pytest-tmp"),
        *(f"{ROOT / TEST_FILE}::{name}" for name in TESTS),
    ]
    with (
        (directory / "pytest.log").open("x", encoding="utf-8") as log,
        redirect_stdout(log),
        redirect_stderr(log),
    ):
        from mtp.multiplicity import sources
        from mtp.workflow import comparison

        loaded: dict[str, str] | None = None
        if changed_file is not None:
            module = comparison if changed_file == COMPARISON else sources
            copy = directory / "copy" / changed_file
            # Re-execute in the original namespace so test monkeypatches still reach it.
            SourceFileLoader(module.__name__, str(copy)).exec_module(module)
            if changed_file == COMPARISON:
                vars(sources)["verify_comparison_sources"] = comparison.verify_comparison_sources
            function = (
                comparison.verify_comparison_sources
                if changed_file == COMPARISON
                else sources.open_verified_sources
            )
            if Path(function.__code__.co_filename) != copy:
                raise RuntimeError("Worker did not load the isolated mutated module")
            loaded = {
                "module": module.__name__,
                "executed_source": str(copy),
                "sha256": digest(copy.read_bytes()),
            }
        code = int(pytest.main(arguments, plugins=[reports]))
        (directory / "pytest.json").write_text(
            json.dumps(
                {
                    "arguments": arguments,
                    "exit_code": code,
                    "counts": dict(reports.counts),
                    "failures": reports.failures,
                    "loaded_mutant": loaded,
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
    raise SystemExit(code)


def intended(case: str, result: dict[str, Any]) -> bool:
    if case == "baseline":
        return (
            result.get("exit_code") == 0
            and result.get("counts") == {"passed": 14}
            and result.get("failures") == []
        )
    expected = MUTATIONS[case][3]
    failures = result.get("failures", [])
    if result.get("exit_code") != 1 or result.get("counts") != {
        "passed": 14 - len(expected),
        "failed": len(expected),
    }:
        return False
    if len(failures) != len(expected) or any(item["phase"] != "call" for item in failures):
        return False
    found = {item["nodeid"].split("::")[-1]: item["failure"] for item in failures}
    return set(found) == expected and all("DID NOT RAISE" in text for text in found.values())


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--timeout-seconds", type=int, default=900)
    args = parser.parse_args()
    if args.timeout_seconds < 300:
        raise ValueError("Fixture watchdog must allow at least 300 seconds per worker")
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError("Choose a fresh mutant evidence directory")
    from mtp.workflow.output_roots import require_durable_output_path

    output = require_durable_output_path(args.output, allowed=False)
    if output.is_relative_to(ROOT) or ROOT.is_relative_to(output):
        raise ValueError("Mutant evidence must be disjoint from the source working tree")
    originals = {name: (ROOT / name).read_bytes() for name in BOUND_FILES}
    for case, (name, before, _after, _expected) in MUTATIONS.items():
        if originals[name].count(before.encode()) != 1:
            raise ValueError(f"{case}: mutation must match exactly one source block")
    repository = Repo(ROOT)
    producer = {
        "head": repository.head.commit.hexsha,
        "tree": repository.head.commit.tree.hexsha,
        "status": repository.git.status("--porcelain", "--untracked-files=all"),
    }
    output.mkdir(parents=True, exist_ok=False)
    records: list[dict[str, Any]] = []
    context = multiprocessing.get_context("spawn")
    for case in ("baseline", *MUTATIONS):
        directory = output / case
        directory.mkdir()
        changed_file, mutation = None, None
        if case != "baseline":
            changed_file, before, after, _expected = MUTATIONS[case]
            raw = originals[changed_file].replace(before.encode(), after.encode(), 1)
            copy = directory / "copy" / changed_file
            copy.parent.mkdir(parents=True)
            copy.write_bytes(raw)
            mutation = {
                "file": changed_file,
                "before": before,
                "after": after,
                "original_sha256": digest(originals[changed_file]),
                "mutant_sha256": digest(raw),
            }
        process = context.Process(target=worker, args=(directory, changed_file))
        try:
            process.start()
            process.join(timeout=args.timeout_seconds)
            timed_out = process.is_alive()
        finally:
            if process.is_alive():
                process.terminate()
                process.join(timeout=5)
                if process.is_alive():
                    process.kill()
                    process.join(timeout=5)
        result_path = directory / "pytest.json"
        evidence_saved = all(
            path.is_file() and path.stat().st_size > 0
            for path in (result_path, directory / "pytest.log")
        )
        result = json.loads(result_path.read_bytes()) if evidence_saved else {}
        result.update(
            case=case, mutation=mutation, process_exit_code=process.exitcode, timed_out=timed_out
        )
        result["evidence_saved"] = evidence_saved
        unchanged = all((ROOT / name).read_bytes() == raw for name, raw in originals.items())
        result["shared_files_unchanged"] = unchanged
        result["intended_outcome"] = (
            not timed_out
            and evidence_saved
            and unchanged
            and process.exitcode == result.get("exit_code")
            and intended(case, result)
        )
        records.append(result)
        print(f"{case}: {result.get('counts')}; intended={result['intended_outcome']}", flush=True)
        if not result["intended_outcome"]:
            break
    final = {name: digest((ROOT / name).read_bytes()) for name in BOUND_FILES}
    report = {
        "issue": 1215,
        "producer": producer,
        "watchdog_seconds": args.timeout_seconds,
        "input_sha256": {name: digest(raw) for name, raw in originals.items()},
        "final_sha256": final,
        "shared_files_unchanged": all(
            final[name] == digest(raw) for name, raw in originals.items()
        ),
        "no_shared_source_writes": True,
        "cases": records,
        "environment": {
            "python": platform.python_version(),
            "executable": sys.executable,
            "TMPDIR": os.environ.get("TMPDIR"),
            **{name: version(name) for name in ("numpy", "pandas", "scipy", "pytest")},
        },
    }
    (output / "summary.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return int(len(records) != 4 or not all(item["intended_outcome"] for item in records))


if __name__ == "__main__":
    raise SystemExit(main())
