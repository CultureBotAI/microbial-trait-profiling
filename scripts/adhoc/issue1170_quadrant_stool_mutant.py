"""Prove exact raw-stool quadrant selection in isolated copies (#1170, #1184).

Inputs: the derived-table source and independent row-2/row-3 synthetic tests.
Outputs: durable pytest logs and a hash-bound JSON proof under --output.
Only an isolated copy of the quadrant is_stool assignment is mutated, in a spawned
test process; canonical production bytes are never rewritten in place.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import multiprocessing
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

import mtp

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "src/mtp/delta_z_derived/tables.py"
TEST_FILE = "tests/test_delta_z_derived_tables.py"
TESTS = (
    "test_row_two_quadrants_exclude_multisite_nonaxis_pairs",
    "test_row_three_counts_every_trait_cell_against_exact_raw_stool_metadata",
)
BEFORE = 'result["is_stool"] = result.body_site.map(is_stool_body_site)'
AFTER = 'result["is_stool"] = result.body_site.str.contains("Feces", regex=False)'
SOURCES = (
    "src/mtp/delta_z_derived/tables.py",
    TEST_FILE,
    "scripts/adhoc/issue1170_quadrant_stool_mutant.py",
    "uv.lock",
)


def digest(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def hashes() -> dict[str, str]:
    return {name: digest((ROOT / name).read_bytes()) for name in SOURCES}


def verify_checkout_identity() -> None:
    """Ensure imported mtp and target files resolve to this helper's own checkout (#1184)."""
    mtp_file = getattr(mtp, "__file__", None)
    if mtp_file is None or Path(mtp_file).resolve().parent != (ROOT / "src/mtp").resolve():
        raise AssertionError("Imported mtp belongs to a different checkout")
    if SOURCE.resolve().parent != (ROOT / "src/mtp/delta_z_derived").resolve():
        raise AssertionError("Target source path is outside this checkout")
    if (ROOT / TEST_FILE).resolve().parent != (ROOT / "tests").resolve():
        raise AssertionError("Target test path is outside this checkout")


@contextmanager
def load_isolated_module(isolated_source: Path) -> Iterator[ModuleType]:
    """Load mtp.delta_z_derived.tables from an isolated file copy (#1184)."""
    verify_checkout_identity()
    import mtp.delta_z_derived as parent_pkg
    from mtp.delta_z_derived import tables as canonical_module

    module_name = "mtp.delta_z_derived.tables"
    resolved_isolated = isolated_source.resolve()
    if resolved_isolated == SOURCE.resolve():
        raise AssertionError("Isolated module path must not be the canonical tracked source")
    spec = importlib.util.spec_from_file_location(module_name, resolved_isolated)
    if spec is None or spec.loader is None:
        raise AssertionError(f"Cannot load isolated module from {resolved_isolated}")
    module = importlib.util.module_from_spec(spec)
    test_path = (ROOT / TEST_FILE).resolve()
    try:
        sys.modules[module_name] = module
        parent_pkg.tables = module
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
        parent_pkg.tables = canonical_module
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
            tables_mod = getattr(module, "tables_module", None)
            origin = getattr(tables_mod, "__file__", None)
            if origin is None:
                raise AssertionError("Test module tables_module import has no __file__")
            resolved_origin = Path(origin).resolve()
            if self.expected_origin is not None and resolved_origin != self.expected_origin:
                raise AssertionError(
                    f"Expected isolated tables origin {self.expected_origin}, got {resolved_origin}"
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
        isolated_source if isolated_source is not None else output / "isolated_tables.py"
    ).resolve()
    reports = Reports(expected_origin=target_source)
    with (
        load_isolated_module(target_source) as loaded_mod,
        (output / "pytest.log").open("w", encoding="utf-8") as log,
        redirect_stdout(log),
        redirect_stderr(log),
    ):
        code = int(
            pytest.main(
                [
                    "-q",
                    "--tb=short",
                    "--basetemp",
                    str(output / "pytest-tmp"),
                    *(f"{ROOT / TEST_FILE}::{name}" for name in TESTS),
                ],
                plugins=[reports],
            )
        )
        origin = reports.loaded_module_origin or str(Path(loaded_mod.__file__ or "").resolve())
        (output / "pytest.json").write_text(
            json.dumps(
                {
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


def intended_result(case: str, record: dict[str, Any]) -> bool:
    if case == "baseline":
        return record.get("exit_code") == 0 and record.get("counts") == {"passed": 2}
    failures = record.get("failures", [])
    return bool(
        record.get("exit_code") == 1
        and record.get("counts") == {"passed": 1, "failed": 1}
        and len(failures) == 1
        and failures[0]["nodeid"].endswith(f"::{TESTS[0]}")
        and failures[0]["phase"] == "call"
        and "n_comparisons" in failures[0]["failure"]
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    verify_checkout_identity()
    output_root = args.output.expanduser().resolve()
    if output_root.is_relative_to(ROOT) or ROOT.is_relative_to(output_root):
        raise ValueError("Probe output must be disjoint from the source working tree")
    original = SOURCE.read_bytes()
    if original.count(BEFORE.encode()) != 1:
        raise ValueError("Expected exactly one quadrant-specific stool assignment")
    output_root.mkdir(parents=True, exist_ok=False)
    initial = hashes()
    records = []
    context = multiprocessing.get_context("spawn")
    for case in ("baseline", "contains_quadrant_only"):
        output = output_root / case
        output.mkdir()
        mutated = (
            original if case == "baseline" else original.replace(BEFORE.encode(), AFTER.encode(), 1)
        )
        isolated_source = (output / "isolated_tables.py").resolve()
        isolated_source.write_bytes(mutated)
        process = context.Process(target=worker, args=(output, isolated_source))
        record: dict[str, Any] = {}
        try:
            process.start()
            process.join(timeout=300)
            if process.is_alive():
                raise TimeoutError(f"{case}: pytest worker exceeded 300 seconds")
            result_path = output / "pytest.json"
            record = json.loads(result_path.read_text()) if result_path.exists() else {}
            origin_ok = (
                record.get("loaded_module_origin") is not None
                and Path(str(record["loaded_module_origin"])).resolve() == isolated_source
            )
            record.update(
                {
                    "case": case,
                    "isolated_source_path": str(isolated_source),
                    "mutated_source_sha256": digest(isolated_source.read_bytes()),
                    "process_exit_code": process.exitcode,
                    "intended_outcome": process.exitcode == record.get("exit_code")
                    and origin_ok
                    and intended_result(case, record),
                }
            )
        finally:
            if process.is_alive():
                process.terminate()
                process.join(timeout=5)
                if process.is_alive():
                    process.kill()
                    process.join(timeout=5)
        if hashes() != initial:
            raise AssertionError("Source, tests or helper changed or were not restored")
        records.append(record)
        print(f"{case}: {record.get('counts')}; intended={record['intended_outcome']}", flush=True)
    final = hashes()
    source_unchanged = SOURCE.read_bytes() == original and final == initial
    report = {
        "issue": 1170,
        "mutation": {"before": BEFORE, "after": AFTER},
        "source_restored": source_unchanged,
        "source_unchanged": source_unchanged,
        "source_sha256": initial,
        "final_sha256": final,
        "cases": records,
        "environment": {
            "python": platform.python_version(),
            **{name: version(name) for name in ("numpy", "pandas", "scipy", "pytest")},
        },
    }
    (output_root / "summary.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return int(not all(record["intended_outcome"] for record in records) or not source_unchanged)


if __name__ == "__main__":
    raise SystemExit(main())
