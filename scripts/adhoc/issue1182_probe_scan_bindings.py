"""Prove relocated scan byte and path binding guards for #1182 in isolated copies (#1182, #1184).

Inputs: frozen production source and seven synthetic archive/byte regression cases.
Outputs: durable pytest logs, exact mutation descriptions and a hash-bound JSON report.
Each case runs in a fresh spawned process against an isolated source copy under --output;
canonical production bytes are never rewritten in place.
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
from collections.abc import Iterator
from contextlib import contextmanager, redirect_stderr, redirect_stdout
from importlib.metadata import version
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
from git import Repo

import mtp
from mtp.workflow.output_roots import add_temporary_output_option, require_durable_output

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


def verify_checkout_identity() -> None:
    """Ensure imported mtp and target files resolve to this helper's own checkout (#1184)."""
    mtp_file = getattr(mtp, "__file__", None)
    if mtp_file is None or Path(mtp_file).resolve().parent != (ROOT / "src/mtp").resolve():
        raise AssertionError("Imported mtp belongs to a different checkout")
    if (ROOT / SOURCE).resolve().parent != (ROOT / "src/mtp/cellulose_decr").resolve():
        raise AssertionError("Target source path is outside this checkout")
    if (ROOT / TEST_FILE).resolve().parent != (ROOT / "tests").resolve():
        raise AssertionError("Target test path is outside this checkout")


@contextmanager
def load_isolated_module(isolated_source: Path) -> Iterator[ModuleType]:
    """Load mtp.cellulose_decr.selection_scan_scopes from an isolated file copy (#1184)."""
    verify_checkout_identity()
    import mtp.cellulose_decr as parent_pkg
    from mtp.cellulose_decr import selection_scan_scopes as canonical_module

    module_name = "mtp.cellulose_decr.selection_scan_scopes"
    resolved_isolated = isolated_source.resolve()
    if resolved_isolated == (ROOT / SOURCE).resolve():
        raise AssertionError("Isolated module path must not be the canonical tracked source")
    spec = importlib.util.spec_from_file_location(module_name, resolved_isolated)
    if spec is None or spec.loader is None:
        raise AssertionError(f"Cannot load isolated module from {resolved_isolated}")
    module = importlib.util.module_from_spec(spec)
    test_path = (ROOT / TEST_FILE).resolve()
    try:
        sys.modules[module_name] = module
        parent_pkg.selection_scan_scopes = module
        spec.loader.exec_module(module)
        if Path(getattr(module, "__file__", "")).resolve() != resolved_isolated:
            raise AssertionError("Loaded module origin does not match isolated source copy")
        for loaded_name, loaded_mod in sys.modules.copy().items():
            filename = getattr(loaded_mod, "__file__", None)
            if filename and Path(filename).resolve() == test_path:
                sys.modules.pop(loaded_name, None)
        yield module
    finally:
        sys.modules[module_name] = canonical_module
        parent_pkg.selection_scan_scopes = canonical_module
        for loaded_name, loaded_mod in sys.modules.copy().items():
            filename = getattr(loaded_mod, "__file__", None)
            if filename and Path(filename).resolve() == test_path:
                sys.modules.pop(loaded_name, None)


class Reports:
    def __init__(self, expected_origin: Path | None = None) -> None:
        self.expected_origin = expected_origin.resolve() if expected_origin is not None else None
        self.loaded_module_origin: str | None = None
        self.counts: Counter[str] = Counter()
        self.failures: list[dict[str, str]] = []

    def pytest_collection_modifyitems(self, items: list[pytest.Item]) -> None:
        modules = {item.module for item in items if isinstance(item, pytest.Function)}
        for module in modules:
            scopes = getattr(module, "scopes", None)
            origin = getattr(scopes, "__file__", None)
            if origin is None:
                raise AssertionError("Test module scopes import has no __file__")
            resolved_origin = Path(origin).resolve()
            if self.expected_origin is not None and resolved_origin != self.expected_origin:
                raise AssertionError(
                    f"Expected isolated scopes origin {self.expected_origin}, got {resolved_origin}"
                )
            self.loaded_module_origin = str(resolved_origin)

    def pytest_runtest_logreport(self, report: pytest.TestReport) -> None:
        if report.when == "call":
            self.counts[report.outcome] += 1
        if report.failed:
            self.failures.append(
                {"nodeid": report.nodeid, "phase": report.when, "failure": report.longreprtext}
            )


def worker(output: Path, isolated_source: Path | None = None) -> None:
    target_source = (
        isolated_source
        if isolated_source is not None
        else output / "isolated_selection_scan_scopes.py"
    ).resolve()
    reports = Reports(expected_origin=target_source)
    arguments = [
        "-q",
        "--tb=short",
        "--basetemp",
        str(output / "pytest-tmp"),
        *(f"{ROOT / TEST_FILE}::{name}" for name in TESTS),
    ]
    with (
        load_isolated_module(target_source) as loaded_mod,
        (output / "pytest.log").open("w", encoding="utf-8") as log,
        redirect_stdout(log),
        redirect_stderr(log),
    ):
        code = int(pytest.main(arguments, plugins=[reports]))
        origin = reports.loaded_module_origin or str(Path(loaded_mod.__file__ or "").resolve())
        (output / "pytest.json").write_text(
            json.dumps(
                {
                    "entrypoint": "pytest.main",
                    "arguments": arguments,
                    "exit_code": code,
                    "loaded_module_origin": origin,
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
    add_temporary_output_option(parser)
    args = parser.parse_args()
    verify_checkout_identity()

    output = args.output.expanduser().resolve()
    require_durable_output(parser, output, allowed=args.allow_temporary_output)
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
        isolated_source = (directory / "isolated_selection_scan_scopes.py").resolve()
        isolated_source.write_bytes(mutated)
        process = context.Process(target=worker, args=(directory, isolated_source))
        result: dict[str, Any] = {}
        try:
            process.start()
            process.join(timeout=300)
            if process.is_alive():
                raise TimeoutError(f"{case}: pytest worker exceeded 300 seconds")
            result_path = directory / "pytest.json"
            result = json.loads(result_path.read_text()) if result_path.exists() else {}
            origin_ok = (
                result.get("loaded_module_origin") is not None
                and Path(str(result["loaded_module_origin"])).resolve() == isolated_source
            )
            result.update(
                {
                    "case": case,
                    "mutation": mutation,
                    "isolated_source_path": str(isolated_source),
                    "isolated_source_sha256": digest(isolated_source.read_bytes()),
                    "process_exit_code": process.exitcode,
                    "intended_outcome": process.exitcode == result.get("exit_code")
                    and origin_ok
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
        result["restored_source_sha256"] = digest((ROOT / SOURCE).read_bytes())
        if input_hashes() != initial:
            raise AssertionError("Frozen source/test/helper bytes changed or were not restored")
        records.append(result)
        print(f"{case}: {result.get('counts')}; intended={result['intended_outcome']}", flush=True)
        if case == "baseline" and not result["intended_outcome"]:
            break
    final = input_hashes()
    source_unchanged = (ROOT / SOURCE).read_bytes() == original and final == initial
    report = {
        "issue": 1182,
        "producer": producer,
        "input_sha256": initial,
        "final_sha256": final,
        "source_restored": source_unchanged,
        "source_unchanged": source_unchanged,
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
        or not source_unchanged
    )


if __name__ == "__main__":
    raise SystemExit(main())
