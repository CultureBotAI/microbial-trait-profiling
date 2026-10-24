"""Probe #964 sensitivity regressions using isolated source copies (2026-10-23).

Inputs: one mutation and this checkout's synthetic sensitivity tests.
Outputs: a pytest log and JSON evidence under the requested output directory.
No tracked source is edited and no production data are read. The copied module's
__file__ is restored to its canonical identity after loading so controlled
inventory-provider fixtures can test the guard itself, not a foreign-path error.
Run each case in a fresh locked uv process.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from contextlib import ExitStack, redirect_stderr, redirect_stdout
from pathlib import Path
from unittest.mock import patch

import pytest

import mtp
import mtp.presence_depth_families as families
import mtp.presence_depth_sensitivity as sensitivity
from mtp.paths import file_sha256

ROOT = Path(__file__).resolve().parents[2]
TEST = ROOT / "tests/test_presence_depth_sensitivity.py"
SNAPSHOT = """    initially_loaded = implementation_sources()
    # Guard even not-yet-imported canonical modules against changes before a lazy
    # import. Only modules actually loaded are included in published provenance.
    source_snapshot = {
        path.relative_to(REPO_ROOT).as_posix(): file_sha256(path)
        for path in (REPO_ROOT / "src/mtp").rglob("*.py")
    }
    source_snapshot.update(initially_loaded)
"""
TALLY = """if (positive, negative, len(y) - positive - negative) != (
                        native_row.studies_up,
                        native_row.studies_down,
                        native_row.studies_tied,
                    ):"""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "case", choices=("baseline", "count_only", "drop_policy", "late_snapshot", "fixed_scope")
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    original_module = families if args.case == "fixed_scope" else sensitivity
    name = original_module.__name__
    filename = original_module.__file__
    if filename is None:
        raise AssertionError("Imported module has no source path")
    source = Path(filename).resolve()
    if source.parent != ROOT / "src/mtp":
        raise AssertionError("Imported module is outside this checkout")
    original = source.read_bytes()
    before = file_sha256(source)
    changed = original.decode("utf-8")
    if args.case == "count_only":
        if changed.count(TALLY) != 1:
            raise ValueError("Expected one sign-tally guard")
        changed = changed.replace(
            TALLY,
            "if len(y) != (native_row.studies_up + native_row.studies_down "
            "+ native_row.studies_tied):",
        )
    elif args.case == "drop_policy":
        targets = [
            line
            for line in changed.splitlines()
            if line.strip().startswith("paths.add(")
            and ("POLICY_PATH" in line or "multiple_testing_policy_v1.yaml" in line)
        ]
        if len(targets) != 1:
            raise ValueError("Expected one policy inventory entry")
        changed = changed.replace(targets[0], "    # Mutant omits the policy inventory entry.")
    elif args.case == "late_snapshot":
        target = "    columns, profiles, metadata = load_core_evidence(bundle)\n"
        if changed.count(SNAPSHOT) != 1 or changed.count(target) != 1:
            raise ValueError("Expected one snapshot and one core load")
        changed = changed.replace(SNAPSHOT, "").replace(target, target + SNAPSHOT)
    elif args.case == "fixed_scope":
        target = "            strain_rollup=strain_rollup,"
        if changed.count(target) != 1:
            raise ValueError("Expected one scope declaration")
        changed = changed.replace(target, '            strain_rollup="phenotype",')
    copied = output / source.name
    copied.write_text(changed, encoding="utf-8")
    spec = importlib.util.spec_from_file_location(name, copied)
    if spec is None or spec.loader is None:
        raise AssertionError("Cannot load isolated mutation")
    module = importlib.util.module_from_spec(spec)
    selection = {
        "baseline": "native_sign_tally or source_inventory or changed_execution_identity "
        "or verified_core_run",
        "count_only": "native_sign_tally",
        "drop_policy": "source_inventory or changed_execution_identity or verified_core_run",
        "late_snapshot": "changed_execution_identity",
        "fixed_scope": "verified_core_run",
    }[args.case]
    expected = {
        "baseline": None,
        "count_only": "Failed: DID NOT RAISE ValueError",
        "drop_policy": "KeyError: 'methods/multiple_testing_policy_v1.yaml'",
        "late_snapshot": "Failed: DID NOT RAISE ValueError",
        "fixed_scope": "AssertionError: assert {'phenotype'} == {'off'}",
    }[args.case]
    try:
        sys.modules[name] = module
        setattr(mtp, name.rsplit(".", 1)[1], module)
        spec.loader.exec_module(module)
        module.__file__ = str(source)
        with (
            ExitStack() as stack,
            (output / "pytest.log").open("w", encoding="utf-8") as log,
            redirect_stdout(log),
            redirect_stderr(log),
        ):
            if args.case == "fixed_scope":
                stack.enter_context(
                    patch.object(sensitivity, "declare_and_correct", module.declare_and_correct)
                )
            code = int(
                pytest.main(
                    [
                        "-q",
                        "-p",
                        "no:cacheprovider",
                        "--basetemp",
                        str(output / "fixtures"),
                        str(TEST),
                        "-k",
                        selection,
                    ]
                )
            )
    finally:
        sys.modules[name] = original_module
        setattr(mtp, name.rsplit(".", 1)[1], original_module)
    evidence = {
        "case": args.case,
        "pytest_exit": code,
        "expected_assertion": expected,
        "verified": code == (0 if expected is None else 1)
        and (expected is None or expected in (output / "pytest.log").read_text()),
        "source_sha256_before": before,
        "source_sha256_after": file_sha256(source),
        "source_unchanged": source.read_bytes() == original,
        "mutant_source_sha256": file_sha256(copied),
    }
    (output / "evidence.json").write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(evidence, indent=2))
    if not evidence["source_unchanged"]:
        raise AssertionError("Tracked source changed during isolated probe")
    return int(not evidence["verified"])


if __name__ == "__main__":
    raise SystemExit(main())
