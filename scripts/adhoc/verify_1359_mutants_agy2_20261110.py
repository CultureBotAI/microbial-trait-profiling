"""Mutation-testing verifier for issue #1359 (contribution trait-key policy propagation).

Purpose:
    Verify that the tests in tests/test_contribution_trait_policy.py fail with
    AssertionError when targeted mutations are applied to isolated copies of
    src/mtp/contribution_inputs.py, src/mtp/taxon_contributions.py, and
    scripts/create_pd_trait_contributors.py.
Inputs:
    src/mtp/contribution_inputs.py
    src/mtp/taxon_contributions.py
    scripts/create_pd_trait_contributors.py
    tests/test_contribution_trait_policy.py
Outputs:
    Exit code 0 if baseline passes and all 7 mutants are killed by AssertionError;
    raises RuntimeError otherwise.
Issue: #1359
Date: 2026-11-10
"""

from __future__ import annotations

import importlib.util
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
CONTRIBUTION_INPUTS_PATH = REPO_ROOT / "src/mtp/contribution_inputs.py"
TAXON_CONTRIBUTIONS_PATH = REPO_ROOT / "src/mtp/taxon_contributions.py"
PD_CONTRIBUTORS_PATH = REPO_ROOT / "scripts/create_pd_trait_contributors.py"
TEST_PATH = REPO_ROOT / "tests/test_contribution_trait_policy.py"


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
        target="if explicit is not None and explicit != recorded_policy:",
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


def _run_scenario(
    scenario: str,
    ci_source: str,
    tc_source: str,
    pd_source: str,
    test_mod: ModuleType,
) -> None:
    scripts_dir = str(REPO_ROOT / "scripts")
    added_scripts = False
    if scripts_dir not in sys.path:
        sys.path.insert(0, scripts_dir)
        added_scripts = True
    orig_ci = sys.modules.get("mtp.contribution_inputs")
    orig_tc = sys.modules.get("mtp.taxon_contributions")
    try:
        ci_mod = _load_module_from_text(
            "mtp.contribution_inputs", ci_source, "contribution_inputs.py"
        )
        tc_mod = _load_module_from_text(
            "mtp.taxon_contributions", tc_source, "taxon_contributions.py"
        )
        tc_mod.__dict__["__file__"] = str(TAXON_CONTRIBUTIONS_PATH)
        pd_mod = _load_module_from_text(
            "mutant_pd_contributors_1359", pd_source, "create_pd_trait_contributors.py"
        )
        test_mod.__dict__["resolve_inputs_trait_key_policy"] = (
            ci_mod.resolve_inputs_trait_key_policy
        )
        test_mod.__dict__["taxon_contributions"] = tc_mod

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            if scenario in {"all", "resolve_policy"}:
                fn_resolve = test_mod.__dict__[
                    "test_resolve_inputs_trait_key_policy_inherits_and_validates_core_and_column_metadata"
                ]
                fn_resolve(tmp_path / "raw", "raw_object_v1")
                fn_resolve(tmp_path / "union", "facultative_union_v1")
            if scenario in {"all", "taxon_contributions_cli"}:
                fn_cli = test_mod.__dict__[
                    "test_taxon_contributions_cli_inherits_core_manifest_policy_and_rejects_mismatch"
                ]
                mp = pytest.MonkeyPatch()
                try:
                    fn_cli(tmp_path / "cli", mp)
                finally:
                    mp.undo()
            if scenario in {"all", "pd_contributors"}:
                fn_pd = test_mod.__dict__[
                    "test_pd_contributors_inherits_core_policy_and_rejects_policy_mismatch"
                ]
                mp = pytest.MonkeyPatch()
                try:
                    fn_pd(tmp_path / "pd", mp, pd_mod)
                finally:
                    mp.undo()
    finally:
        if orig_ci is not None:
            sys.modules["mtp.contribution_inputs"] = orig_ci
        if orig_tc is not None:
            sys.modules["mtp.taxon_contributions"] = orig_tc
        sys.modules.pop("mutant_pd_contributors_1359", None)
        if added_scripts and scripts_dir in sys.path:
            sys.path.remove(scripts_dir)


def main() -> int:
    ci_original = CONTRIBUTION_INPUTS_PATH.read_text(encoding="utf-8")
    tc_original = TAXON_CONTRIBUTIONS_PATH.read_text(encoding="utf-8")
    pd_original = PD_CONTRIBUTORS_PATH.read_text(encoding="utf-8")
    test_source = TEST_PATH.read_text(encoding="utf-8")
    test_mod = _load_module_from_text(
        "test_contribution_trait_policy_runner_1359",
        test_source,
        "test_contribution_trait_policy.py",
    )
    test_mod.__dict__["__file__"] = str(TEST_PATH)

    _run_scenario("all", ci_original, tc_original, pd_original, test_mod)
    print("Baseline passed all #1359 scenarios.")

    killed = 0
    for mutant in MUTANTS:
        ci_text = ci_original
        tc_text = tc_original
        pd_text = pd_original
        if mutant.module_kind == "contribution_inputs":
            if mutant.target not in ci_text:
                raise RuntimeError(f"Missing target for {mutant.name}: {mutant.target}")
            ci_text = ci_text.replace(mutant.target, mutant.replacement, 1)
        elif mutant.module_kind == "taxon_contributions":
            if mutant.target not in tc_text:
                raise RuntimeError(f"Missing target for {mutant.name}: {mutant.target}")
            tc_text = tc_text.replace(mutant.target, mutant.replacement, 1)
        elif mutant.module_kind == "pd_contributors":
            if mutant.target not in pd_text:
                raise RuntimeError(f"Missing target for {mutant.name}: {mutant.target}")
            pd_text = pd_text.replace(mutant.target, mutant.replacement, 1)
        else:
            raise RuntimeError(f"Unknown module_kind: {mutant.module_kind}")

        try:
            _run_scenario(mutant.scenario, ci_text, tc_text, pd_text, test_mod)
        except AssertionError:
            killed += 1
            print(f"[KILLED] {mutant.name}")
        except Exception as exc:
            raise RuntimeError(
                f"Mutant {mutant.name} raised unexpected {type(exc).__name__}: {exc}"
            ) from exc
        else:
            raise RuntimeError(f"Mutant {mutant.name} survived!")

    if killed != len(MUTANTS):
        raise RuntimeError(f"Expected {len(MUTANTS)} killed mutants, got {killed}")
    print(f"All {killed}/{len(MUTANTS)} mutants killed by AssertionError.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
