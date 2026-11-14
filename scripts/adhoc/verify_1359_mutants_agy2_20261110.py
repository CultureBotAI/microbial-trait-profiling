"""Mutation-testing verifier for issue #1359 (contribution trait-key policy propagation).

Purpose:
    Preserve the seven original policy mutants and check provenance regressions
    using only the new provenance tests. Require AssertionError for isolated copies of
    src/mtp/contribution_inputs.py, src/mtp/taxon_contributions.py, and
    scripts/create_pd_trait_contributors.py.
Inputs:
    src/mtp/contribution_inputs.py
    src/mtp/taxon_contributions.py
    scripts/create_pd_trait_contributors.py
    tests/test_contribution_trait_policy.py
    tests/test_contribution_policy_provenance.py
    tests/test_taxon_policy_provenance.py
    tests/test_pd_policy_provenance.py
Outputs:
    Exit code 0 if baseline passes and all mutants are killed by AssertionError;
    optional exclusive JSON verdicts (--output); raises RuntimeError otherwise.
    --canary runs all controls and the first new mutant through the same harness.
    Temporary module copies and tiny synthetic fixtures use TMPDIR; no tree copy.
Issues: #1359, #2306, #2307
Date: 2026-11-11
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import importlib.util
import inspect
import io
import json
import sys
import tempfile
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType, SimpleNamespace
from typing import cast

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
CONTRIBUTION_INPUTS_PATH = REPO_ROOT / "src/mtp/contribution_inputs.py"
TAXON_CONTRIBUTIONS_PATH = REPO_ROOT / "src/mtp/taxon_contributions.py"
PD_CONTRIBUTORS_PATH = REPO_ROOT / "scripts/create_pd_trait_contributors.py"
sys.dont_write_bytecode = True

SOURCE_PATHS = {
    "contribution_inputs": CONTRIBUTION_INPUTS_PATH,
    "taxon_contributions": TAXON_CONTRIBUTIONS_PATH,
    "pd_contributors": PD_CONTRIBUTORS_PATH,
}
TEST_SCENARIOS: dict[str, tuple[str, tuple[str, ...]]] = {
    "resolve_policy": (
        "test_contribution_trait_policy",
        ("test_resolve_inputs_trait_key_policy_inherits_and_validates_core_and_column_metadata",),
    ),
    "taxon_contributions_cli": (
        "test_contribution_trait_policy",
        ("test_taxon_contributions_cli_inherits_core_manifest_policy_and_rejects_mismatch",),
    ),
    "pd_contributors": (
        "test_contribution_trait_policy",
        ("test_pd_contributors_inherits_core_policy_and_rejects_policy_mismatch",),
    ),
    "resolution_provenance": (
        "test_contribution_policy_provenance",
        (
            "test_policy_resolution_records_basis_and_every_consumed_source",
            "test_unrelated_manifest_is_parsed_metadata_not_core_verification",
            "test_policy_resolution_record_and_export_are_independent_immutable_values",
        ),
    ),
    "resolution_buffers": (
        "test_contribution_policy_provenance",
        ("test_policy_source_hashes_bind_the_same_buffers_that_were_parsed",),
    ),
    "resolution_guards": (
        "test_contribution_policy_provenance",
        (
            "test_policy_provenance_still_refuses_conflicting_declarations",
            "test_buffered_column_metadata_preserves_strict_table_validation",
        ),
    ),
    "taxon_provenance": (
        "test_taxon_policy_provenance",
        (
            "test_taxon_manifest_records_policy_basis_and_consumed_sources",
            "test_taxon_manifest_hashes_the_policy_buffers_used_for_its_counts",
            "test_taxon_cli_describes_policy_metadata_and_the_callers_verification_duty",
        ),
    ),
    "pd_provenance": (
        "test_pd_policy_provenance",
        (
            "test_pd_manifest_records_policy_resolution_sources",
            "test_pd_main_uses_one_captured_policy_resolution",
            "test_pd_loader_refuses_policy_that_contradicts_captured_resolution",
            "test_pd_policy_help_and_loader_do_not_claim_core_verification",
        ),
    ),
}


@dataclass(frozen=True)
class MutantSpec:
    name: str
    module_kind: str
    target: str
    replacement: str
    scenario: str


MUTANTS: tuple[MutantSpec, ...] = (
    MutantSpec(
        name="M1_contribution_inputs_ignore_core_manifest",
        module_kind="contribution_inputs",
        target="if manifest_path is not None:",
        replacement="if False:",
        scenario="resolve_policy",
    ),
    MutantSpec(
        name="M2_contribution_inputs_ignore_column_metadata",
        module_kind="contribution_inputs",
        target="if columns_path.is_file():",
        replacement="if False:",
        scenario="resolve_policy",
    ),
    MutantSpec(
        name="M3_contribution_inputs_skip_mismatch_guard",
        module_kind="contribution_inputs",
        target=(
            "if recorded_policy is not None and explicit is not None "
            "and explicit != recorded_policy:"
        ),
        replacement="if False:",
        scenario="resolve_policy",
    ),
    MutantSpec(
        name="M4_taxon_contributions_graph_forces_union",
        module_kind="taxon_contributions",
        target="trait_key_policy=policy,",
        replacement='trait_key_policy="facultative_union_v1",',
        scenario="taxon_contributions_cli",
    ),
    MutantSpec(
        name="M5_taxon_contributions_manifest_forces_union",
        module_kind="taxon_contributions",
        target='"trait_key_policy": policy,',
        replacement='"trait_key_policy": "facultative_union_v1",',
        scenario="taxon_contributions_cli",
    ),
    MutantSpec(
        name="M6_pd_contributors_graph_forces_union",
        module_kind="pd_contributors",
        target="trait_key_policy=resolved_policy,",
        replacement='trait_key_policy="facultative_union_v1",',
        scenario="pd_contributors",
    ),
    MutantSpec(
        name="M7_pd_contributors_manifest_forces_union",
        module_kind="pd_contributors",
        target='"trait_key_policy": policy,',
        replacement='"trait_key_policy": "facultative_union_v1",',
        scenario="pd_contributors",
    ),
    MutantSpec(
        "M8_default_is_misreported_as_explicit",
        "contribution_inputs",
        'InputTraitPolicyResolution(DEFAULT_TRAIT_KEY_POLICY, "default", source_sha256)',
        'InputTraitPolicyResolution(DEFAULT_TRAIT_KEY_POLICY, "explicit", source_sha256)',
        "resolution_provenance",
    ),
    MutantSpec(
        "M9_explicit_is_misreported_as_default",
        "contribution_inputs",
        'InputTraitPolicyResolution(explicit, "explicit", source_sha256)',
        'InputTraitPolicyResolution(explicit, "default", source_sha256)',
        "resolution_provenance",
    ),
    MutantSpec(
        "M10_manifest_is_misreported_as_columns",
        "contribution_inputs",
        'InputTraitPolicyResolution(manifest_policy, "core_manifest", source_sha256)',
        'InputTraitPolicyResolution(manifest_policy, "column_metadata", source_sha256)',
        "resolution_provenance",
    ),
    MutantSpec(
        "M11_columns_are_misreported_as_default",
        "contribution_inputs",
        'InputTraitPolicyResolution(column_policy, "column_metadata", source_sha256)',
        'InputTraitPolicyResolution(column_policy, "default", source_sha256)',
        "resolution_provenance",
    ),
    MutantSpec(
        "M12_manifest_hash_is_omitted",
        "contribution_inputs",
        'source_sha256["core_manifest"] = hashlib.sha256(manifest_data).hexdigest()',
        'source_sha256.pop("core_manifest", None)',
        "resolution_provenance",
    ),
    MutantSpec(
        "M13_columns_hash_is_omitted",
        "contribution_inputs",
        'source_sha256["taxa_trait_matrix_columns.tsv"] = hashlib.sha256(columns_data).hexdigest()',
        'source_sha256.pop("taxa_trait_matrix_columns.tsv", None)',
        "resolution_provenance",
    ),
    MutantSpec(
        "M14_manifest_hash_reads_later_bytes",
        "contribution_inputs",
        "hashlib.sha256(manifest_data).hexdigest()",
        "hashlib.sha256(manifest_path.read_bytes()).hexdigest()",
        "resolution_buffers",
    ),
    MutantSpec(
        "M15_columns_hash_reads_later_bytes",
        "contribution_inputs",
        "hashlib.sha256(columns_data).hexdigest()",
        "hashlib.sha256(columns_path.read_bytes()).hexdigest()",
        "resolution_buffers",
    ),
    MutantSpec(
        "M16_columns_parse_reopens_the_file",
        "contribution_inputs",
        "for row in _column_policy_rows(columns_data, columns_path)",
        "for row in rows(columns_path)",
        "resolution_buffers",
    ),
    MutantSpec(
        "M17_explicit_metadata_conflict_is_accepted",
        "contribution_inputs",
        "if recorded_policy is not None and explicit is not None and explicit != recorded_policy:",
        "if False:",
        "resolution_guards",
    ),
    MutantSpec(
        "M18_two_metadata_policies_are_accepted",
        "contribution_inputs",
        "        and manifest_policy != column_policy\n",
        "        and False\n",
        "resolution_guards",
    ),
    MutantSpec(
        "M19_taxon_resolution_provenance_is_omitted",
        "taxon_contributions",
        '"trait_key_policy_resolution": policy_resolution.as_dict(),',
        '"trait_key_policy_resolution": {},',
        "taxon_provenance",
    ),
    MutantSpec(
        "M20_taxon_manifest_resolves_later_metadata",
        "taxon_contributions",
        '"trait_key_policy_resolution": policy_resolution.as_dict(),',
        '"trait_key_policy_resolution": resolve_inputs_trait_key_policy_record('
        "args.inputs, explicit_policy, core_manifest=args.core_manifest).as_dict(),",
        "taxon_provenance",
    ),
    MutantSpec(
        "M21_pd_resolution_provenance_is_omitted",
        "pd_contributors",
        '"trait_key_policy_resolution": resolution.as_dict(),',
        '"trait_key_policy_resolution": {},',
        "pd_provenance",
    ),
    MutantSpec(
        "M22_pd_ledger_resolves_later_metadata",
        "pd_contributors",
        "        policy_resolution=resolution,\n",
        "",
        "pd_provenance",
    ),
    MutantSpec(
        "M23_pd_accepts_a_contradictory_captured_policy",
        "pd_contributors",
        "if trait_key_policy is not None and trait_key_policy != policy_resolution.policy:",
        "if False:",
        "pd_provenance",
    ),
    MutantSpec(
        "M24_taxon_help_claims_verification",
        "taxon_contributions",
        '"Optional core run.json to parse and validate trait-key policy metadata. "',
        '"Optional core run.json manifest to verify and inherit trait-key policy from. "',
        "taxon_provenance",
    ),
    MutantSpec(
        "M25_pd_help_claims_verification",
        "pd_contributors",
        '"Parse and validate policy metadata from a core run.json manifest; "',
        '"Verify and inherit trait-key policy from a core run.json manifest; "',
        "pd_provenance",
    ),
    MutantSpec(
        "M26_pd_loader_doc_claims_verification",
        "pd_contributors",
        "analysis here. Policy metadata is parsed and validated; the caller must verify the core",
        "analysis here. The policy is inherited from a verified core supplied by the caller",
        "pd_provenance",
    ),
)


def _load_module_from_text(module_name: str, source_text: str, filename: str) -> ModuleType:
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_file = Path(tmp_dir) / filename
        tmp_file.write_text(source_text, encoding="utf-8")
        spec = importlib.util.spec_from_file_location(module_name, tmp_file)
        if spec is None or spec.loader is None:
            raise RuntimeError(f"Could not build module spec for {module_name}")
        module = importlib.util.module_from_spec(spec)
        sys.modules[module_name] = module
        spec.loader.exec_module(module)
        return module


def _parameter_cases(function: Callable[..., None]) -> list[dict[str, object]]:
    """Expand only the declared, small parametrize cases from these shipped tests."""
    cases: list[dict[str, object]] = [{}]
    for marker in getattr(function, "pytestmark", []):
        if marker.name != "parametrize":
            continue
        raw_names, values = marker.args[:2]
        names = raw_names.split(",") if isinstance(raw_names, str) else raw_names
        names = [name.strip() for name in names]
        expanded = []
        for case in cases:
            for value in values:
                supplied = (value,) if len(names) == 1 else value
                expanded.append({**case, **dict(zip(names, supplied, strict=True))})
        cases = expanded
    return cases


@dataclass
class _Capture:
    stdout: io.StringIO
    stderr: io.StringIO

    def readouterr(self) -> SimpleNamespace:
        captured = SimpleNamespace(out=self.stdout.getvalue(), err=self.stderr.getvalue())
        for stream in (self.stdout, self.stderr):
            stream.seek(0)
            stream.truncate(0)
        return captured


def _invoke_test(function: Callable[..., None], directory: Path, pd_mod: ModuleType) -> int:
    cases = _parameter_cases(function)
    for index, parameters in enumerate(cases):
        case_path = directory / str(index)
        case_path.mkdir(parents=True)
        with pytest.MonkeyPatch.context() as monkeypatch:
            capture = _Capture(io.StringIO(), io.StringIO())
            fixtures = {
                "tmp_path": case_path,
                "monkeypatch": monkeypatch,
                "pd_contributors_mod": pd_mod,
                "capsys": capture,
            }
            arguments = dict(parameters)
            for name in inspect.signature(function).parameters:
                if name not in arguments:
                    if name not in fixtures:
                        raise RuntimeError(f"Unsupported fixture {name} in {function.__name__}")
                    arguments[name] = fixtures[name]
            try:
                with (
                    contextlib.redirect_stdout(capture.stdout),
                    contextlib.redirect_stderr(capture.stderr),
                ):
                    function(**arguments)
            except AssertionError as error:
                detail = str(error) or "assertion rejected the mutated behavior"
                raise AssertionError(f"{function.__name__}[{index}]: {detail}") from error
    return len(cases)


def _run_scenario(scenario: str, sources: dict[str, str], test_sources: dict[str, str]) -> int:
    module_names = (
        "mtp.contribution_inputs",
        "mtp.taxon_contributions",
        "mutant_pd_contributors_1359",
        *test_sources,
    )
    originals = {name: sys.modules.get(name) for name in module_names}
    total = 0
    try:
        with pytest.MonkeyPatch.context() as imports:
            for directory in (REPO_ROOT / "scripts", REPO_ROOT / "tests"):
                imports.syspath_prepend(str(directory))
            package = importlib.import_module("mtp")
            ci_mod = _load_module_from_text(
                "mtp.contribution_inputs", sources["contribution_inputs"], "contribution_inputs.py"
            )
            imports.setattr(package, "contribution_inputs", ci_mod, raising=False)
            tc_mod = _load_module_from_text(
                "mtp.taxon_contributions", sources["taxon_contributions"], "taxon_contributions.py"
            )
            tc_mod.__dict__["__file__"] = str(TAXON_CONTRIBUTIONS_PATH)
            imports.setattr(package, "taxon_contributions", tc_mod, raising=False)
            pd_mod = _load_module_from_text(
                "mutant_pd_contributors_1359",
                sources["pd_contributors"],
                "create_pd_trait_contributors.py",
            )
            modules = {}
            # Load the original fixture provider before tests that import its small builders.
            for name, source in test_sources.items():
                module = _load_module_from_text(name, source, f"{name}.py")
                module.__dict__["__file__"] = str(REPO_ROOT / "tests" / f"{name}.py")
                modules[name] = module
            scenarios = (
                TEST_SCENARIOS if scenario == "all" else {scenario: TEST_SCENARIOS[scenario]}
            )
            with tempfile.TemporaryDirectory(prefix="verify-1359-cases-") as temporary:
                for group, (module_name, functions) in scenarios.items():
                    for name in functions:
                        function = cast(Callable[..., None], getattr(modules[module_name], name))
                        total += _invoke_test(function, Path(temporary) / group / name, pd_mod)
    finally:
        for name, original in originals.items():
            if original is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = original
    return total


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--canary", action="store_true", help="Run controls and one new mutant.")
    parser.add_argument("--output", type=Path, help="Write exclusive JSON verdicts to this path.")
    arguments = parser.parse_args()
    sources = {kind: path.read_text(encoding="utf-8") for kind, path in SOURCE_PATHS.items()}
    test_names = dict.fromkeys(name for name, _functions in TEST_SCENARIOS.values())
    test_sources = {
        name: (REPO_ROOT / "tests" / f"{name}.py").read_text(encoding="utf-8")
        for name in test_names
    }
    bound_paths = [*SOURCE_PATHS.values(), *(REPO_ROOT / "tests" / f"{n}.py" for n in test_names)]
    before = {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in bound_paths}
    count = _run_scenario("all", sources, test_sources)
    verdicts = [{"name": "control", "status": "passed", "detail": f"{count} test cases passed"}]
    print(f"Baseline passed {count} #1359/#2306/#2307 cases.")
    selected = MUTANTS[7:8] if arguments.canary else MUTANTS
    if not selected:
        raise RuntimeError("No mutants selected")
    for mutant in selected:
        changed = dict(sources)
        source = changed[mutant.module_kind]
        if source.count(mutant.target) != 1:
            raise RuntimeError(f"Expected exactly one target for {mutant.name}: {mutant.target}")
        changed[mutant.module_kind] = source.replace(mutant.target, mutant.replacement, 1)
        try:
            _run_scenario(mutant.scenario, changed, test_sources)
        except AssertionError as error:
            detail = str(error) or "assertion rejected the mutated behavior"
            verdicts.append(
                {
                    "name": mutant.name,
                    "status": "killed",
                    "detail": detail,
                    "exception": "AssertionError",
                    "tests": mutant.scenario,
                }
            )
            print(f"[KILLED by AssertionError] {mutant.name}: {detail}")
        except Exception as error:
            raise RuntimeError(
                f"Mutant {mutant.name} raised unexpected {type(error).__name__}: {error}"
            ) from error
        else:
            raise RuntimeError(f"Mutant {mutant.name} survived")
    after = {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in bound_paths}
    if before != after:
        raise RuntimeError("Live source/test files changed during verification")
    report = {
        "verdicts": verdicts,
        "source_unchanged": True,
        "source_sha256": {
            path.relative_to(REPO_ROOT).as_posix(): value for path, value in before.items()
        },
    }
    if arguments.output is not None:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        with arguments.output.open("x", encoding="utf-8") as handle:
            json.dump(report, handle, indent=2)
            handle.write("\n")
    print(f"All {len(selected)}/{len(selected)} mutants killed by AssertionError.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
