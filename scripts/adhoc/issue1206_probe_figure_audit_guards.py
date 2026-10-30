"""Prove figure byte/coverage guards with isolated source copies (#1206, 2026-10-28).

Inputs: current input reader, auditor and synthetic regression/control cases.
Outputs: fresh durable copy-only mutants, pytest logs and hash-bound result JSON.
Shared production and tests are never changed. No real F/P files or figure run.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
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

from mtp.workflow.output_roots import require_durable_output_path

ROOT = Path(__file__).resolve().parents[2]
READER = "src/mtp/publication_figures/inputs.py"
AUDITOR = "scripts/adhoc/issue1046_check_figure_bundle.py"
INPUT_TEST = "tests/test_publication_figure_inputs.py"
AUDIT_TEST = "tests/test_publication_figure_audit.py"
COPY_FILES = (READER, AUDITOR, INPUT_TEST, AUDIT_TEST)
BOUND_FILES = (*COPY_FILES, "scripts/adhoc/issue1206_probe_figure_audit_guards.py", "uv.lock")
READ_CASE = "test_parses_checked_bytes_not_a_second_path_read"
AUDIT_CASE = "test_provenance_hash_binds_parsed_bytes"
TABLE_CASES = tuple(
    f"test_table_coverage_refused[{name}]" for name in ("thin", "thin_redeclare", "substitute")
)
QUADRANT_CASES = tuple(f"test_quadrant_coverage_refused[{name}]" for name in ("marks", "swap"))
SOURCE_CASES = tuple(
    f"test_table_source_parent_output_refused[{name}]"
    for name in ("direct", "logical_alias", "resolved_alias")
)
CASES = (
    (INPUT_TEST, READ_CASE),
    (AUDIT_TEST, AUDIT_CASE),
    (AUDIT_TEST, "test_table_covers_bound_parent_traits_across_pages"),
    *((AUDIT_TEST, name) for name in TABLE_CASES),
    (AUDIT_TEST, "test_quadrant_includes_every_population_comparison_and_axis"),
    *((AUDIT_TEST, name) for name in QUADRANT_CASES),
    *((AUDIT_TEST, name) for name in SOURCE_CASES),
)
MUTANTS = {
    "reader_second_read": (
        READER,
        "    return raw\n",
        "    return path.read_bytes()\n",
        {READ_CASE},
    ),
    "auditor_second_read": (
        AUDITOR,
        "    return raw\n",
        "    return path.read_bytes()\n",
        {AUDIT_CASE},
    ),
    "drop_table_coverage": (
        AUDITOR,
        "            table_coverage(pages, traits)",
        "            pass",
        set(TABLE_CASES),
    ),
    "drop_quadrant_coverage": (
        AUDITOR,
        "            quadrant_coverage(pages)",
        "            pass",
        set(QUADRANT_CASES),
    ),
    "drop_source_output_guard": (
        AUDITOR,
        "    protect_table_sources(output, aliases)\n",
        "",
        set(SOURCE_CASES),
    ),
}


def sha(raw: bytes) -> str:
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


def worker(directory: Path) -> None:
    reports = Reports()
    copied = directory / "copy"
    arguments = [
        "-q",
        "--tb=short",
        "--basetemp",
        str(directory / "pytest-tmp"),
        *(f"{copied / filename}::{case}" for filename, case in CASES),
    ]
    with (
        (directory / "pytest.log").open("x", encoding="utf-8") as log,
        redirect_stdout(log),
        redirect_stderr(log),
    ):
        # Only the copied reader is substituted; its dependencies retain the locked environment.
        from mtp import publication_figures

        name = "mtp.publication_figures.inputs"
        spec = importlib.util.spec_from_file_location(name, copied / READER)
        if spec is None or spec.loader is None:
            raise RuntimeError("Cannot load isolated input reader")
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
        vars(publication_figures)["inputs"] = module
        code = int(pytest.main(arguments, plugins=[reports]))
        (directory / "pytest.json").write_text(
            json.dumps(
                {
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


def intended(case: str, result: dict[str, Any]) -> bool:
    if case == "baseline":
        return (
            result.get("exit_code") == 0
            and result.get("counts") == {"passed": len(CASES)}
            and not result.get("failures")
        )
    expected = MUTANTS[case][3]
    failures = result.get("failures", [])
    if result.get("exit_code") != 1 or result.get("counts") != {
        "passed": len(CASES) - len(expected),
        "failed": len(expected),
    }:
        return False
    if len(failures) != len(expected) or any(item["phase"] != "call" for item in failures):
        return False
    by_case = {item["nodeid"].split("::")[-1]: item["failure"] for item in failures}
    if set(by_case) != expected:
        return False
    if case == "reader_second_read":
        return "Parsed value differs from the hashed first read" in by_case[READ_CASE]
    if case == "auditor_second_read":
        return "Parsed producer differs from the hashed first read" in by_case[AUDIT_CASE]
    return all("DID NOT RAISE" in failure for failure in by_case.values())


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = require_durable_output_path(args.output, allowed=False)
    if output.is_relative_to(ROOT) or ROOT.is_relative_to(output):
        raise ValueError("Probe output must be disjoint from the source working tree")
    if output.exists() or output.is_symlink():
        raise FileExistsError("Choose a fresh mutation-evidence directory")
    originals = {name: (ROOT / name).read_bytes() for name in BOUND_FILES}
    for filename, before, _after, _expected in MUTANTS.values():
        if originals[filename].count(before.encode()) != 1:
            raise ValueError("Each isolated mutation must match exactly one source expression")
    repository = Repo(ROOT)
    producer = {
        "head": repository.head.commit.hexsha,
        "tree": repository.head.commit.tree.hexsha,
        "status": repository.git.status("--porcelain", "--untracked-files=all"),
    }
    output.mkdir(parents=True, exist_ok=False)
    results: list[dict[str, Any]] = []
    context = multiprocessing.get_context("spawn")
    for case in ("baseline", *MUTANTS):
        directory = output / case
        directory.mkdir()
        mutation = None
        for filename in COPY_FILES:
            raw = originals[filename]
            if case != "baseline" and filename == MUTANTS[case][0]:
                _file, before, after, _expected = MUTANTS[case]
                raw = raw.replace(before.encode(), after.encode(), 1)
                mutation = {
                    "file": filename,
                    "before": before,
                    "after": after,
                    "original_sha256": sha(originals[filename]),
                    "mutant_sha256": sha(raw),
                }
            target = directory / "copy" / filename
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(raw)
        process = context.Process(target=worker, args=(directory,))
        try:
            process.start()
            process.join(timeout=180)
            if process.is_alive():
                raise TimeoutError(f"{case}: isolated pytest exceeded 180 seconds")
            result_path = directory / "pytest.json"
            result = json.loads(result_path.read_bytes()) if result_path.exists() else {}
            result.update(case=case, mutation=mutation, process_exit_code=process.exitcode)
            result["intended_outcome"] = process.exitcode == result.get("exit_code") and intended(
                case, result
            )
            results.append(result)
            print(
                f"{case}: {result.get('counts')}; intended={result['intended_outcome']}", flush=True
            )
        finally:
            if process.is_alive():
                process.terminate()
                process.join(timeout=5)
                if process.is_alive():
                    process.kill()
                    process.join(timeout=5)
        if any((ROOT / name).read_bytes() != raw for name, raw in originals.items()):
            raise AssertionError("Shared source/test/helper changed during the isolated probe")
        if not result["intended_outcome"]:
            break
    final = {name: sha((ROOT / name).read_bytes()) for name in BOUND_FILES}
    report = {
        "issue": 1206,
        "producer": producer,
        "input_sha256": {name: sha(raw) for name, raw in originals.items()},
        "final_sha256": final,
        "shared_sources_unchanged": all(final[name] == sha(raw) for name, raw in originals.items()),
        "no_shared_source_writes": True,
        "cases": results,
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
        len(results) != 1 + len(MUTANTS) or not all(item["intended_outcome"] for item in results)
    )


if __name__ == "__main__":
    raise SystemExit(main())
