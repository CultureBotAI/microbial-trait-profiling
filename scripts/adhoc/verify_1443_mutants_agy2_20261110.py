"""Mutation-testing verifier for issue #1443 (structural MT correction detection & scope).

Purpose:
    Verify that the structural correction detection, alias/partial/getattr/bare-script
    resolution, and policy scope tests in tests/test_multiple_testing_registry.py fail
    with AssertionError under targeted code/policy mutations.
Inputs:
    tests/test_multiple_testing_registry.py
    methods/multiple_testing_policy.md
Outputs:
    Exit code 0 if baseline passes and all 8 mutants are killed by AssertionError;
    raises RuntimeError otherwise.
Issue: #1443
Date: 2026-11-10
"""

from __future__ import annotations

import importlib.util
import sys
import tempfile
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType

REPO_ROOT = Path(__file__).resolve().parents[2]
TEST_PATH = REPO_ROOT / "tests/test_multiple_testing_registry.py"


@dataclass(frozen=True)
class MutantSpec:
    name: str
    target: str
    replacement: str
    test_names: tuple[str, ...]


MUTANTS: tuple[MutantSpec, ...] = (
    MutantSpec(
        name="M1_disable_accumulate_detection",
        target='dotted(node.func, limits).endswith("minimum.accumulate")',
        replacement="False",
        test_names=("test_structural_cumulative_minimum_and_step_up_loop_are_detected",),
    ),
    MutantSpec(
        name="M2_disable_step_up_formula",
        target="if _is_step_up_formula(node, limits):",
        replacement="if False:",
        test_names=("test_structural_cumulative_minimum_and_step_up_loop_are_detected",),
    ),
    MutantSpec(
        name="M3_disable_clamped_bonferroni",
        target="if _is_clamped_bonferroni(node, limits):",
        replacement="if False:",
        test_names=("test_inline_bonferroni_and_sidak_corrections_are_detected",),
    ),
    MutantSpec(
        name="M4_disable_sidak_formula",
        target="if _is_sidak_formula(node, limits):",
        replacement="if False:",
        test_names=("test_inline_bonferroni_and_sidak_corrections_are_detected",),
    ),
    MutantSpec(
        name="M5_disable_partial_resolution",
        target='return callee in {"partial", "functools.partial"}',
        replacement="return False",
        test_names=("test_reexported_alias_partial_getattr_and_bare_script_imports_are_detected",),
    ),
    MutantSpec(
        name="M6_disable_getattr_unwrapping",
        target='and dotted(node.func, limits) in {"getattr", "builtins.getattr"}',
        replacement="and False",
        test_names=("test_reexported_alias_partial_getattr_and_bare_script_imports_are_detected",),
    ),
    MutantSpec(
        name="M7_disable_bare_script_forwarding",
        target='if qualified.startswith("scripts."):',
        replacement="if False:",
        test_names=("test_reexported_alias_partial_getattr_and_bare_script_imports_are_detected",),
    ),
    MutantSpec(
        name="M8_regress_object_method_false_positive",
        target='return defined or "." not in name or name in NAMED_WRAPPERS',
        replacement="return True",
        test_names=("test_reexported_alias_partial_getattr_and_bare_script_imports_are_detected",),
    ),
)


def _load_module_from_source(module_name: str, source_text: str) -> ModuleType:
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_file = Path(tmp_dir) / "test_multiple_testing_registry_mutant.py"
        tmp_file.write_text(source_text, encoding="utf-8")
        spec = importlib.util.spec_from_file_location(module_name, tmp_file)
        if spec is None or spec.loader is None:
            raise RuntimeError(f"Could not create module spec for {module_name}")
        module = importlib.util.module_from_spec(spec)
        sys.modules[module_name] = module
        try:
            spec.loader.exec_module(module)
        finally:
            sys.modules.pop(module_name, None)
        module.__dict__["ROOT"] = REPO_ROOT
        module.__dict__["REGISTRY"] = REPO_ROOT / "methods/multiple_testing_families.yaml"
        return module


def _run_named_tests(module: ModuleType, test_names: tuple[str, ...]) -> None:
    for test_name in test_names:
        fn: Callable[[], None] = getattr(module, test_name)
        fn()


def main() -> int:
    original_source = TEST_PATH.read_text(encoding="utf-8")
    baseline_mod = _load_module_from_source("mt_registry_baseline_1443", original_source)
    all_new_tests = (
        "test_structural_cumulative_minimum_and_step_up_loop_are_detected",
        "test_inline_bonferroni_and_sidak_corrections_are_detected",
        "test_reexported_alias_partial_getattr_and_bare_script_imports_are_detected",
        "test_policy_enforcement_scope_and_out_of_scope_list_match_guard",
    )
    _run_named_tests(baseline_mod, all_new_tests)
    print("Baseline passed all 4 new #1443 test functions.")

    killed = 0
    for index, mutant in enumerate(MUTANTS, 1):
        if mutant.target not in original_source:
            raise RuntimeError(f"Target string missing for {mutant.name}: {mutant.target}")
        mutated_source = original_source.replace(mutant.target, mutant.replacement, 1)
        mutant_mod = _load_module_from_source(f"mt_registry_mutant_{index}", mutated_source)
        try:
            _run_named_tests(mutant_mod, mutant.test_names)
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
