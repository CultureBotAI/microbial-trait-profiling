"""Qualify section 6 rendering against the previous source for #1328 (2026-11-04).

Inputs: --baseline-ref, current report source, and the existing synthetic test
fixtures. Outputs: retained sources, identical input snapshots, before/after report
bytes, diffs, JSON equality checks, and a directly invoked regression test's proof.
Requires a fresh, explicit --output-dir in durable caller data outside this checkout.
Does not run either audit pipeline, read production inputs, or modify tracked files.
"""

from __future__ import annotations

import argparse
import copy
import difflib
import hashlib
import importlib.util
import json
import math
import sys
import traceback
from dataclasses import asdict, is_dataclass
from pathlib import Path
from types import ModuleType
from typing import Any

from git import Repo

import mtp
from mtp.trait_lists import DeclaredTraitList
from mtp.workflow.output_roots import require_durable_output_path

ROOT = Path(__file__).resolve().parents[2]
REPORT_SOURCE = "scripts/audit_trait_density_correlation.py"
FIXTURE_SOURCE = "tests/test_density_consolidated_conclusions.py"
SECTION6 = b"6. CONCLUSION & METHODOLOGICAL IMPLICATIONS\n"


def sha256(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def json_value(value: Any) -> Any:
    """Keep undefined/infinite synthetic values explicit in standards-compliant JSON."""
    if is_dataclass(value) and not isinstance(value, type):
        return json_value(asdict(value))
    if isinstance(value, dict):
        return {key: json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [json_value(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return {"nonfinite_float": str(value)}
    return value


def json_bytes(value: Any) -> bytes:
    text = json.dumps(json_value(value), indent=2, sort_keys=True, allow_nan=False)
    return (text + "\n").encode()


def load_module(name: str, source: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, source)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import retained source {source}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def compare_scenario(
    output: Path,
    before: ModuleType,
    after: ModuleType,
    density: dict[str, Any],
    ontology: dict[str, Any],
    trait_list: DeclaredTraitList | None,
) -> dict[str, Any]:
    """Render separate copies of one input and retain byte comparisons and all reports."""
    output.mkdir()
    inputs = (density, ontology, trait_list)
    original = json_bytes(inputs)
    (output / "input.json").write_bytes(original)
    reports = {}
    checks = {}
    for label, renderer in (("before", before), ("after", after)):
        destination = output / label
        destination.mkdir()
        supplied = copy.deepcopy(inputs)
        supplied_bytes = json_bytes(supplied)
        (destination / "input.json").write_bytes(supplied_bytes)
        checks[f"{label}_received_identical_input"] = supplied_bytes == original
        report = renderer.write_consolidated_report(destination, *supplied)
        content = report.read_bytes()
        if not content or content.count(SECTION6) != 1:
            raise ValueError(f"Expected one section 6 header in nonempty {report}")
        checks[f"{label}_input_unchanged"] = json_bytes(supplied) == original
        reports[label] = content
    before_prefix, _, before_section6 = reports["before"].partition(SECTION6)
    after_prefix, _, after_section6 = reports["after"].partition(SECTION6)
    delta = b"".join(
        difflib.diff_bytes(
            difflib.unified_diff,
            reports["before"].splitlines(keepends=True),
            reports["after"].splitlines(keepends=True),
            fromfile=b"before/density_adjustment_report.txt",
            tofile=b"after/density_adjustment_report.txt",
        )
    )
    (output / "report.diff").write_bytes(delta)
    checks["sections1_5_byte_identical"] = before_prefix == after_prefix
    if trait_list is None:
        checks["section6_changed"] = before_section6 != after_section6
        checks["nonempty_report_diff"] = bool(delta)
    else:
        checks["declared_trait_list_whole_report_identical"] = reports["before"] == reports["after"]
        checks["declared_trait_list_empty_diff"] = not delta
    result = {
        "input_sha256": sha256(original),
        "declared_trait_list": trait_list is not None,
        "report_bytes": {name: len(content) for name, content in reports.items()},
        "report_sha256": {name: sha256(content) for name, content in reports.items()},
        "sections1_5_sha256": {
            "before": sha256(before_prefix),
            "after": sha256(after_prefix),
        },
        "section6_sha256": {
            "before": sha256(before_section6),
            "after": sha256(after_section6),
        },
        "diff_bytes": len(delta),
        "checks": checks,
    }
    (output / "checks.json").write_bytes(json_bytes(result))
    return result


def regression_proof(
    output: Path, fixtures: ModuleType, before: ModuleType, after: ModuleType
) -> dict[str, Any]:
    """Invoke one unchanged test function against each source; this is not a pytest run."""
    output.mkdir()
    test_name = "test_conclusions_report_profile_and_paired_diagnostics_without_severity_claims"
    test = getattr(fixtures, test_name)
    parameters = {"profile_r": 0.0, "profile_vif": 1.0, "paired_r": 0.2, "paired_vif": 1.04}
    results = {}
    for label, renderer in (("before", before), ("after", after)):
        destination = output / label
        destination.mkdir()
        try:
            test(tmp_path=destination, audit_script=renderer, **parameters)
        except AssertionError:
            failure = traceback.format_exc()
            (destination / "assertion_failure.txt").write_text(failure, encoding="utf-8")
            results[label] = {"passed": False, "failure": "assertion_failure.txt"}
        else:
            results[label] = {"passed": True, "failure": None}
    result = {
        "mechanism": "Direct invocation of one existing test function; not the pytest runner",
        "test": f"{FIXTURE_SOURCE}::{test_name}",
        "parameters": parameters,
        "results": results,
        "checks": {
            "historical_regression_reproduced": not results["before"]["passed"],
            "current_regression_passed": results["after"]["passed"],
        },
    }
    (output / "proof.json").write_bytes(json_bytes(result))
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-ref", required=True, help="Immutable pre-change Git revision")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(f"mtp.__file__: {mtp.__file__}")
    expected_package = ROOT / "src/mtp/__init__.py"
    if mtp.__file__ is None or Path(mtp.__file__).resolve() != expected_package.resolve():
        raise RuntimeError(f"Expected mtp to import from {expected_package}")
    output = require_durable_output_path(args.output_dir.expanduser(), allowed=False)
    if output.is_relative_to(ROOT):
        parser.error("--output-dir must be outside the source checkout in durable caller data")
    if output.exists():
        parser.error(
            "--output-dir must be a fresh directory; existing evidence is never overwritten"
        )
    repository = Repo(ROOT)
    baseline_commit = repository.commit(args.baseline_ref)
    baseline_bytes = (baseline_commit.tree / REPORT_SOURCE).data_stream.read()
    watched = [
        REPORT_SOURCE,
        FIXTURE_SOURCE,
        str(Path(__file__).resolve().relative_to(ROOT)),
        "src/mtp/density_adjustment.py",
        "src/mtp/ontology_audit.py",
        "src/mtp/trait_lists.py",
    ]
    current_sources = {relative: (ROOT / relative).read_bytes() for relative in watched}
    output.mkdir(parents=True, exist_ok=False)
    source_dir = output / "sources"
    source_dir.mkdir()
    retained_before = source_dir / "baseline_report.py"
    retained_after = source_dir / "current_report.py"
    retained_fixtures = source_dir / "synthetic_fixtures.py"
    retained_before.write_bytes(baseline_bytes)
    retained_after.write_bytes(current_sources[REPORT_SOURCE])
    retained_fixtures.write_bytes(current_sources[FIXTURE_SOURCE])
    before = load_module("issue1328_baseline_report", retained_before)
    after = load_module("issue1328_current_report", retained_after)
    fixtures = load_module("issue1328_synthetic_fixtures", retained_fixtures)
    scenarios = {}
    for case in ("negative_correlation", "low_vif", "undefined"):
        density, ontology = fixtures.synthetic_results(case)
        scenarios[case] = compare_scenario(output / case, before, after, density, ontology, None)
    density, ontology = fixtures.synthetic_results()
    declared = DeclaredTraitList(
        tuple(effect.trait_id for effect in density["table1_effects"]), "fixture.txt", "a" * 64
    )
    scenarios["declared_trait_list"] = compare_scenario(
        output / "declared_trait_list", before, after, density, ontology, declared
    )
    regression = regression_proof(output / "regression", fixtures, before, after)
    checks = {
        "all_scenario_checks_passed": all(
            all(scenario["checks"].values()) for scenario in scenarios.values()
        ),
        "regression_proof_passed": all(regression["checks"].values()),
        "watched_sources_unchanged": all(
            (ROOT / relative).read_bytes() == content
            for relative, content in current_sources.items()
        ),
    }
    evidence = {
        "issue": 1328,
        "scope": "Synthetic report rendering only; no audit pipeline or production input",
        "baseline_ref_requested": args.baseline_ref,
        "baseline_commit": baseline_commit.hexsha,
        "imported_mtp_source": str(Path(mtp.__file__).resolve()),
        "baseline_report_source_sha256": sha256(baseline_bytes),
        "current_source_sha256": {
            relative: sha256(content) for relative, content in current_sources.items()
        },
        "scenarios": scenarios,
        "regression": regression,
        "checks": checks,
        "passed": all(checks.values()),
    }
    result_path = output / "qualification.json"
    result_path.write_bytes(json_bytes(evidence))
    if not result_path.is_file() or not result_path.stat().st_size:
        raise RuntimeError("Qualification evidence was not persisted")
    print(f"{'PASS' if evidence['passed'] else 'FAIL'}: {result_path}")
    return 0 if evidence["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
