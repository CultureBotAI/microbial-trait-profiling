"""Prove relocated scan byte and path binding guards for #1182 (2026-10-28).

Inputs: frozen production source and seven synthetic archive/byte regression cases.
Outputs: durable pytest logs, exact mutation descriptions and a hash-bound JSON report.
Each case runs in a fresh spawned process. Production bytes are restored in finally;
no historical F data, pipeline rerun, or real-data qualification is invoked.
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
from importlib.metadata import version
from pathlib import Path
from typing import Any

import pytest
from git import Repo

ROOT = Path(__file__).resolve().parents[2]
SOURCE = "src/mtp/cellulose_decr/selection_scan_scopes.py"
TEST_FILE = "tests/test_cellulose_decr_selection_scan_scopes.py"
ROUTE_TEST = "test_load_run_resolves_every_output_and_bound_core_without_live_reads"
BYTE_TEST = "test_loaded_bytes_must_still_match_the_verified_output_inventory"
ROUTE_CASES = tuple(f"{ROUTE_TEST}[{scope}]" for scope in ("off", "phenotype"))
BYTE_CASES = tuple(f"{BYTE_TEST}[{name}]" for name in ("summary", "provenance", "traits"))
EQUALITY_TEST = "test_loaded_consumer_record_must_equal_the_verified_record"
VERSION_TEST = "test_verified_snapshots_must_still_declare_v2_after_preflight"
TESTS = (*ROUTE_CASES, *BYTE_CASES, EQUALITY_TEST, VERSION_TEST)
MUTATIONS = {
    "bypass_bound_core_resolver": (
        "manifest = locate(Path(record.run_manifest_path), resolver).read_bytes()",
        "manifest = Path(record.run_manifest_path).read_bytes()",
    ),
    "drop_output_sha256_check": (
        "if digest != output.sha256 or len(raw) != output.size_bytes:",
        "if len(raw) != output.size_bytes:",
    ),
    "drop_verified_consumer_equality": (
        "if ConsumerRecord.model_validate_json(record_bytes) != record:",
        "if False:",
    ),
    "drop_final_v2_check": (
        "    _require_v2(directory, summary, provenance)",
        "    pass",
    ),
}
EXPECTED_FAILURES = {
    "bypass_bound_core_resolver": set(ROUTE_CASES),
    "drop_output_sha256_check": set(BYTE_CASES),
    "drop_verified_consumer_equality": {EQUALITY_TEST},
    "drop_final_v2_check": {VERSION_TEST},
}
BOUND_FILES = (
    SOURCE,
    "src/mtp/cellulose_decr/selection_scan.py",
    "src/mtp/workflow/consumers.py",
    "src/mtp/workflow/archive.py",
    "src/mtp/workflow/verification.py",
    "tests/test_cellulose_decr_selection_scan.py",
    TEST_FILE,
    "scripts/adhoc/issue1182_probe_scan_bindings.py",
    "pyproject.toml",
    "uv.lock",
)


def digest(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def input_hashes() -> dict[str, str]:
    return {name: digest((ROOT / name).read_bytes()) for name in BOUND_FILES}


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


def worker(output: Path) -> None:
    reports = Reports()
    arguments = [
        "-q",
        "--tb=short",
        "--basetemp",
        str(output / "pytest-tmp"),
        *(f"{ROOT / TEST_FILE}::{name}" for name in TESTS),
    ]
    with (
        (output / "pytest.log").open("w", encoding="utf-8") as log,
        redirect_stdout(log),
        redirect_stderr(log),
    ):
        code = int(pytest.main(arguments, plugins=[reports]))
        (output / "pytest.json").write_text(
            json.dumps(
                {
                    "entrypoint": "pytest.main",
                    "arguments": arguments,
                    "exit_code": code,
                    "counts": dict(reports.counts),
                    "failures": reports.failures,
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
    raise SystemExit(code)


def intended_failure(case: str, record: dict[str, Any]) -> bool:
    if case == "baseline":
        return (
            record.get("exit_code") == 0
            and record.get("counts") == {"passed": len(TESTS)}
            and record.get("failures") == []
        )
    expected = EXPECTED_FAILURES[case]
    failures = record.get("failures", [])
    if record.get("exit_code") != 1 or record.get("counts") != {
        "passed": len(TESTS) - len(expected),
        "failed": len(expected),
    }:
        return False
    if len(failures) != len(expected) or any(item["phase"] != "call" for item in failures):
        return False
    by_case = {item["nodeid"].split("::")[-1]: item["failure"] for item in failures}
    if set(by_case) != expected:
        return False
    if case == "bypass_bound_core_resolver":
        return all(
            "A relocated scan read an original path:" in failure and "run.json" in failure
            for failure in by_case.values()
        )
    return all("DID NOT RAISE" in failure for failure in by_case.values())


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    # Import the shared main-branch guard, not the source module being mutated.
    from mtp.workflow.output_roots import require_durable_output

    output = args.output.expanduser().resolve()
    require_durable_output(parser, output, allowed=False)
    if output.is_relative_to(ROOT) or ROOT.is_relative_to(output):
        raise ValueError("Campaign output must be disjoint from the source working tree")
    original = (ROOT / SOURCE).read_bytes()
    for case, (before, _) in MUTATIONS.items():
        if original.count(before.encode()) != 1:
            raise ValueError(f"{case} must match exactly one production expression")
    initial = input_hashes()
    repository = Repo(ROOT)
    producer = {
        "head": repository.head.commit.hexsha,
        "tree": repository.head.commit.tree.hexsha,
        "worktree_status": repository.git.status("--porcelain", "--untracked-files=all"),
    }
    output.mkdir(parents=True, exist_ok=False)
    records: list[dict[str, Any]] = []
    context = multiprocessing.get_context("spawn")
    for case in ("baseline", *MUTATIONS):
        directory = output / case
        directory.mkdir()
        mutated = original
        mutation = None
        if case != "baseline":
            before, after = MUTATIONS[case]
            mutated = original.replace(before.encode(), after.encode(), 1)
            mutation = {
                "file": SOURCE,
                "before": before,
                "after": after,
                "original_sha256": digest(original),
                "mutated_sha256": digest(mutated),
            }
        process = context.Process(target=worker, args=(directory,))
        result: dict[str, Any] = {}
        try:
            if case != "baseline":
                (ROOT / SOURCE).write_bytes(mutated)
            process.start()
            process.join(timeout=300)
            if process.is_alive():
                raise TimeoutError(f"{case}: pytest worker exceeded 300 seconds")
            result_path = directory / "pytest.json"
            result = json.loads(result_path.read_text()) if result_path.exists() else {}
            result.update(
                {
                    "case": case,
                    "mutation": mutation,
                    "process_exit_code": process.exitcode,
                    "intended_outcome": process.exitcode == result.get("exit_code")
                    and intended_failure(case, result),
                }
            )
        finally:
            if process.is_alive():
                process.terminate()
                process.join(timeout=5)
                if process.is_alive():
                    process.kill()
                    process.join(timeout=5)
            if case != "baseline":
                (ROOT / SOURCE).write_bytes(original)
        result["restored_source_sha256"] = digest((ROOT / SOURCE).read_bytes())
        if input_hashes() != initial:
            raise AssertionError("Frozen source/test/helper bytes changed or were not restored")
        records.append(result)
        print(f"{case}: {result.get('counts')}; intended={result['intended_outcome']}", flush=True)
        if case == "baseline" and not result["intended_outcome"]:
            break
    final = input_hashes()
    report = {
        "issue": 1182,
        "producer": producer,
        "input_sha256": initial,
        "final_sha256": final,
        "source_restored": (ROOT / SOURCE).read_bytes() == original,
        "cases": records,
        "environment": {
            "python": platform.python_version(),
            "executable": sys.executable,
            "platform": platform.platform(),
            "TMPDIR": os.environ.get("TMPDIR"),
            **{name: version(name) for name in ("numpy", "pandas", "scipy", "pytest")},
        },
    }
    (output / "summary.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return int(
        len(records) != 1 + len(MUTATIONS)
        or not all(record["intended_outcome"] for record in records)
    )


if __name__ == "__main__":
    raise SystemExit(main())
