"""Verify influence guards in isolated source copies (#2078/#2079; 2026-11-08).

Inputs: this checkout's producers, synthetic tests and configuration. Outputs:
an exclusive directory containing the copied source, diff, pytest log and bound
evidence.json. Run exactly one case per fresh process, first with --control,
then with --control-evidence pointing to that control's evidence.json.

The caller supplies any heavy-job scheduling, RSS cap and timeout. A terminated
process has no successful receipt. This probes logical guards, not producer
inventory tampering: an extension copy retains the canonical inventory function.
Canonical source is never rewritten. No network or production data is needed.
"""

from __future__ import annotations

import argparse
import ast
import difflib
import hashlib
import importlib
import importlib.util
import json
import os
import stat
import sys
from contextlib import redirect_stderr, redirect_stdout
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

from mtp.paths import REPO_ROOT
from mtp.workflow.output_roots import add_temporary_output_option, require_durable_output

ROOT = Path(__file__).resolve().parents[2]
MAX_FILE_BYTES = 4 * 1024 * 1024
MAX_TOTAL_BYTES = 16 * 1024 * 1024
EXTENSION = "src/mtp/influence_extension.py"
NOEL = "src/mtp/influence_noel.py"
STATUS = "tests/test_influence_extension_status.py"
CHECK = "tests/test_influence_extension_noel_check.py"
CORE = "tests/test_influence_extension_noel_core.py"


@dataclass(frozen=True)
class Case:
    """One exact logical mutation and its public test witness."""

    name: str
    module: str
    source: str
    old: str
    new: str
    test: str


def _cases() -> dict[str, Case]:
    records: list[Case] = []

    def add(number: int, source: str, old: str, new: str, test: str) -> None:
        name = f"M{number:02}"
        records.append(Case(name, "mtp." + Path(source).stem, source, old, new, test))

    add(
        1,
        EXTENSION,
        '    if missing:\n        return "absent_lineages" if missing == set(lineages) '
        'else "partially_absent_lineage"\n',
        "",
        STATUS + "::test_lineage_presence_status_preserves_effect_and_directional_mass",
    )
    add(
        2,
        EXTENSION,
        "        for parents in graph.parents.values():\n",
        "        for parents in ():\n",
        STATUS
        + "::test_lineage_presence_status_preserves_effect_and_directional_mass[all_present]",
    )
    add(
        3,
        EXTENSION,
        '                            "removal_status": '
        "_removal_status(lineages, present, removed_mass),\n",
        '                            "removal_status": "removed_annotation_mass" '
        'if removed_mass else "no_removed_annotation_mass",\n',
        STATUS
        + "::test_absence_status_takes_precedence_even_when_one_direction_has_zero_removed_mass",
    )
    check = (
        "                _match(sidecar_values[field], value, "
        "f\"{scope}/{row['trait_key']}/{field}\")\n"
    )
    for number, field in enumerate(("p", "q_bh_six", "q_by_six"), 4):
        add(
            number,
            NOEL,
            check,
            f'                if field != "{field}":\n    ' + check,
            CHECK + f"::test_refuses_sidecar_mismatch_in_each_compared_quantity[{field}]",
        )
    add(
        7,
        NOEL,
        "ABSOLUTE_TOLERANCE = 0.0\n",
        "ABSOLUTE_TOLERANCE = 1e-12\n",
        CHECK + "::test_absolute_tolerance_does_not_hide_tiny_p_mismatches",
    )
    add(
        8,
        NOEL,
        "        by = false_discovery_control([1.0 if p is None else p for p in pvalues], "
        'method="by")\n',
        "        tested_by = iter(false_discovery_control([p for p in pvalues if p is not None], "
        'method="by"))\n'
        "        by = [1.0 if p is None else float(next(tested_by)) for p in pvalues]\n",
        CHECK + "::test_refuses_by_values_computed_after_dropping_untestable_slot",
    )
    add(
        9,
        NOEL,
        "                if key not in expected or key in selected:\n",
        "                if key not in expected:\n",
        CHECK + "::test_refuses_incomplete_or_ambiguous_native_family[duplicate]",
    )
    add(
        10,
        NOEL,
        "    if set(selected) != expected:\n        raise ValueError("
        'f"Missing Noel baseline family rows: {sorted(expected - set(selected))}")\n',
        "    for scope, trait in expected - set(selected):\n"
        "        template = next(iter(selected.values()))\n"
        '        selected[(scope, trait)] = {**template, "scope": scope, "trait_key": trait, '
        '"before_p": "", "before_q": ""}\n',
        CHECK + "::test_refuses_incomplete_or_ambiguous_native_family[missing]",
    )
    add(
        11,
        NOEL,
        "        if hashed.digest.hexdigest() != expected_sha256:\n",
        "        if hashlib.sha256(path.read_bytes()).hexdigest() != expected_sha256:\n",
        CHECK + "::test_consumed_gzip_hash_cannot_be_replaced_by_a_fresh_disk_hash",
    )
    add(
        12,
        EXTENSION,
        "    if noel_manifest is not None:\n        noel_check = check_noel_baseline(\n"
        '            inputs["noel_effects"], noel_sources["noel_effects"][1], baseline\n'
        "        )\n",
        "",
        CORE + "::test_verified_but_mismatched_noel_probability_refuses_before_publication",
    )
    add(
        13,
        EXTENSION,
        "    run = _verified_noel(core_manifest, noel_manifest)\n    result = {}\n",
        "    run = verify_arm(noel_manifest)\n    result = {}\n",
        CORE + "::test_verified_noel_from_different_core_refuses_before_publication",
    )
    add(
        14,
        NOEL,
        '                    "p_status": "matched",\n'
        '                    "bh_status": "matched",\n'
        '                    "by_status": "recomputed_from_noel_p",\n',
        '                    "status": "matched",\n',
        CHECK + "::test_checks_every_scope_family_and_labels_by_as_derived",
    )
    return {record.name: record for record in records}


CASES = _cases()
WATCHED = (
    EXTENSION,
    NOEL,
    STATUS,
    CHECK,
    CORE,
    "tests/test_influence_extension_joint.py",
    "tests/test_influence_extension_core.py",
    "tests/test_issue2078_influence_probe.py",
    "tests/conftest.py",
    "src/mtp/multiplicity/sources.py",
    "src/mtp/contribution_inputs.py",
    "src/mtp/taxon_contributions.py",
    "src/mtp/taxon_influence.py",
    "src/mtp/taxa_trait_matrix.py",
    "src/mtp/trait_guilds.py",
    "src/mtp/trait_key_policy.py",
    "src/mtp/workflow/operations.py",
    "src/mtp/workflow/consumers.py",
    "src/mtp/workflow/trait_policy.py",
    "src/mtp/workflow/arms/runtime.py",
    "scripts/adhoc/issue1030_noel_fixture.py",
    "scripts/adhoc/issue2078_probe_influence_guards.py",
    "methods/table1_influence_extension.md",
    "methods/table1_influence_extension_protocol.json",
    "methods/table1_taxon_influence_protocol.json",
    "methods/multiple_testing_families.yaml",
    "pyproject.toml",
    "uv.lock",
    "Justfile",
    ".secrets.baseline",
)


def read_bound(path: Path) -> tuple[bytes, dict[str, Any]]:
    """Read a bounded regular file and bind its bytes and filesystem identity."""
    before = path.lstat()
    if not stat.S_ISREG(before.st_mode) or before.st_size > MAX_FILE_BYTES:
        raise ValueError(f"Expected bounded regular input: {path.name}")
    with path.open("rb") as stream:
        raw = stream.read(MAX_FILE_BYTES + 1)
    after = path.lstat()
    fields = ("st_dev", "st_ino", "st_mode", "st_size", "st_mtime_ns")
    if (
        len(raw) > MAX_FILE_BYTES
        or len(raw) != before.st_size
        or any(getattr(before, key) != getattr(after, key) for key in fields)
    ):
        raise ValueError(f"Canonical input changed while reading: {path.name}")
    return raw, {
        "sha256": hashlib.sha256(raw).hexdigest(),
        **{key: getattr(after, key) for key in fields},
    }


def snapshot(root: Path, paths: tuple[str, ...]) -> dict[str, dict[str, Any]]:
    """Bind the explicit source/test/config inventory without scanning a checkout."""
    records = {name: read_bound(root / name)[1] for name in paths}
    if sum(record["st_size"] for record in records.values()) > MAX_TOTAL_BYTES:
        raise ValueError("Canonical input inventory exceeds its size bound")
    return records


def require_snapshot(root: Path, before: dict[str, dict[str, Any]]) -> None:
    if snapshot(root, tuple(before)) != before:
        raise ValueError("Canonical input changed during probe")


def mutated_source(case: Case, source: bytes, *, control: bool) -> bytes:
    """Even controls refuse a stale anchor instead of validating an obsolete case."""
    old = case.old.encode()
    if source.count(old) != 1:
        raise ValueError(f"Expected exactly one mutation anchor for {case.name}")
    changed = source.replace(old, case.new.encode(), 1)
    ast.parse(changed, filename=case.source)
    return source if control else changed


def require_module_origin(module: ModuleType, canonical: Path) -> None:
    file = getattr(module, "__file__", None)
    spec = getattr(module, "__spec__", None)
    origin = getattr(spec, "origin", None)
    if (
        not file
        or not origin
        or Path(file).resolve() != canonical.resolve()
        or Path(origin).resolve() != canonical.resolve()
    ):
        raise ValueError(f"Unexpected module origin: {module.__name__}")


class Results:
    """Use actual call exception types, never diagnostic-text assertion lookalikes."""

    def __init__(self) -> None:
        self.reports: list[dict[str, str]] = []
        self.exceptions: list[dict[str, Any]] = []
        self.collection_errors: list[str] = []

    def pytest_runtest_logreport(self, report: Any) -> None:
        self.reports.append(
            {"nodeid": report.nodeid, "when": report.when, "outcome": report.outcome}
        )

    def pytest_collectreport(self, report: Any) -> None:
        if report.failed:
            self.collection_errors.append(report.nodeid)

    def pytest_exception_interact(self, call: Any, report: Any) -> None:
        if call.excinfo is not None:
            error = call.excinfo.value
            detail = str(error)[:16000]
            assertion = isinstance(error, AssertionError) or (
                isinstance(error, pytest.fail.Exception) and detail.startswith("DID NOT RAISE")
            )
            self.exceptions.append(
                {
                    "nodeid": report.nodeid,
                    "when": getattr(report, "when", "collect"),
                    "type": type(error).__name__,
                    "detail": detail,
                    "assertion": assertion,
                }
            )

    def accepted(self, exit_code: int, target: str, *, control: bool) -> bool:
        if self.collection_errors or any(
            report["outcome"] == "skipped"
            or (report["when"] != "call" and report["outcome"] != "passed")
            for report in self.reports
        ):
            return False
        calls = [report for report in self.reports if report["when"] == "call"]
        if (
            not calls
            or len({report["nodeid"] for report in calls}) != len(calls)
            or any(
                not (report["nodeid"] == target or report["nodeid"].startswith(target + "["))
                for report in calls
            )
        ):
            return False
        failures = [report for report in calls if report["outcome"] == "failed"]
        if control:
            return (
                exit_code == 0
                and not self.exceptions
                and all(report["outcome"] == "passed" for report in calls)
            )
        return (
            exit_code == 1
            and bool(failures)
            and len(self.exceptions) == len(failures)
            and all(
                error["when"] == "call"
                and error["assertion"]
                and any(report["nodeid"] == error["nodeid"] for report in failures)
                for error in self.exceptions
            )
        )


def require_control_receipt(path: Path, case: Case, before: dict[str, dict[str, Any]]) -> str:
    raw, metadata = read_bound(path)
    record = json.loads(raw)
    expected = {
        "case": case.name,
        "mode": "control",
        "verified": True,
        "pytest_exit": 0,
        "target": case.test,
        "canonical_before": before,
        "canonical_after": before,
    }
    if not isinstance(record, dict) or any(
        record.get(key) != value for key, value in expected.items()
    ):
        raise ValueError("Mutation requires a matching passing control receipt")
    return str(metadata["sha256"])


def _install_copy(case: Case, copied: Path) -> None:
    original = importlib.import_module(case.module)
    require_module_origin(original, ROOT / case.source)
    spec = importlib.util.spec_from_file_location(case.module, copied)
    if spec is None or spec.loader is None:
        raise ValueError("Cannot load isolated influence source")
    module = importlib.util.module_from_spec(spec)
    sys.modules[case.module] = module
    spec.loader.exec_module(module)
    if case.source == EXTENSION:
        # Preserve canonical inventory paths; this is a logical-guard experiment.
        vars(module)["implementation_hashes"] = original.implementation_hashes
    setattr(importlib.import_module("mtp"), case.module.rsplit(".", 1)[1], module)
    require_module_origin(module, copied)
    if case.source == NOEL:
        # Install the checker before its consumer binds the imported function.
        extension = importlib.import_module("mtp.influence_extension")
        require_module_origin(extension, ROOT / EXTENSION)
        if extension.check_noel_baseline is not module.check_noel_baseline:
            raise ValueError("Extension did not bind the isolated Noel checker")


def run_probe(case: Case, output: Path, *, control: bool, control_evidence: Path | None) -> bool:
    if REPO_ROOT.resolve() != ROOT:
        raise ValueError("Imported mtp belongs to a different checkout")
    package = importlib.import_module("mtp")
    require_module_origin(package, ROOT / "src/mtp/__init__.py")
    if any(name in sys.modules for name in ("mtp.influence_extension", "mtp.influence_noel")):
        raise ValueError("Run one influence probe per fresh process")
    before = snapshot(ROOT, WATCHED)
    original, _ = read_bound(ROOT / case.source)
    copied_bytes = mutated_source(case, original, control=control)
    if control and control_evidence is not None:
        raise ValueError("A control cannot consume another control receipt")
    if not control and control_evidence is None:
        raise ValueError("Mutation requires a matching passing control receipt")
    control_sha = (
        None
        if control_evidence is None
        else require_control_receipt(control_evidence, case, before)
    )
    output.mkdir(parents=True, exist_ok=False)
    copied = output / Path(case.source).name
    with copied.open("xb") as stream:
        stream.write(copied_bytes)
    with (output / "mutation.diff").open("x", encoding="utf-8") as stream:
        stream.writelines(
            difflib.unified_diff(
                original.decode().splitlines(True),
                copied_bytes.decode().splitlines(True),
                fromfile=case.source,
                tofile=copied.name,
            )
        )
    _install_copy(case, copied)
    require_snapshot(ROOT, before)
    results = Results()
    arguments = [
        "-q",
        "--tb=short",
        "-o",
        "addopts=",
        "-p",
        "no:cacheprovider",
        f"--basetemp={output / 'pytest'}",
        case.test,
    ]
    previous_cwd = Path.cwd()
    old_environment = {
        key: os.environ.get(key)
        for key in ("PYTEST_ADDOPTS", "PYTEST_DISABLE_PLUGIN_AUTOLOAD", "PYTHONDONTWRITEBYTECODE")
    }
    old_bytecode = sys.dont_write_bytecode
    try:
        os.chdir(ROOT)
        os.environ.update(
            PYTEST_ADDOPTS="", PYTEST_DISABLE_PLUGIN_AUTOLOAD="1", PYTHONDONTWRITEBYTECODE="1"
        )
        sys.dont_write_bytecode = True
        with (
            (output / "pytest.log").open("x", encoding="utf-8") as log,
            redirect_stdout(log),
            redirect_stderr(log),
        ):
            exit_code = int(pytest.main(arguments, plugins=[results]))
    finally:
        os.chdir(previous_cwd)
        sys.dont_write_bytecode = old_bytecode
        for key, value in old_environment.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value
    after = snapshot(ROOT, WATCHED)
    isolated, _ = read_bound(copied)
    verified = (
        results.accepted(exit_code, case.test, control=control)
        and before == after
        and isolated == copied_bytes
    )
    evidence = {
        "case": case.name,
        "mode": "control" if control else "mutation",
        "target": case.test,
        "pytest_exit": exit_code,
        "reports": results.reports,
        "exceptions": results.exceptions,
        "collection_errors": results.collection_errors,
        "canonical_before": before,
        "canonical_after": after,
        "isolated_source": copied.name,
        "isolated_sha256": hashlib.sha256(isolated).hexdigest(),
        "control_evidence_sha256": control_sha,
        "inventory_scope": "logical_guards_with_canonical_producer_inventory",
        "verified": verified,
    }
    with (output / "evidence.json").open("x", encoding="utf-8") as stream:
        json.dump(evidence, stream, indent=2)
        stream.write("\n")
    return verified


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=CASES)
    parser.add_argument("--control", action="store_true")
    parser.add_argument("--control-evidence", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    add_temporary_output_option(parser)
    args = parser.parse_args(argv)
    output = args.output_dir.expanduser()
    require_durable_output(parser, output, allowed=args.allow_temporary_output)
    verified = run_probe(
        CASES[args.case],
        output.resolve(),
        control=args.control,
        control_evidence=args.control_evidence,
    )
    print(json.dumps({"case": args.case, "control": args.control, "verified": verified}))
    return int(not verified)


if __name__ == "__main__":
    raise SystemExit(main())
