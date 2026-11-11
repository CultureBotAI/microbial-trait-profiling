#!/usr/bin/env python3
"""Adversarial mutation verifier for issue #1635 (2026-11-10).

Purpose:
    Verify that tests/test_table1_native_results_readme_1635.py kills at least
    6 distinct mutations in isolated scratch copies of
    methods/table1-native-results/README.md.
Inputs:
    Repository root (defaults to two parent directories above this script).
Outputs:
    Exits 0 when the baseline passes and every mutant is killed by AssertionError.
"""

from __future__ import annotations

import importlib.util
import tempfile
from pathlib import Path
from types import ModuleType


def _load_module(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load module spec from {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> int:
    repo_root = Path(__file__).resolve().parents[2]
    readme_src = repo_root / "methods/table1-native-results/README.md"
    test_src = repo_root / "tests/test_table1_native_results_readme_1635.py"

    test_module = _load_module("test_table1_native_results_readme_1635", test_src)
    readme_text = readme_src.read_text(encoding="utf-8")

    # Verify baseline passes
    test_module.check_readme_text(readme_text)

    mutants: list[tuple[str, str]] = [
        (
            "mutant_1_revert_to_stale_all_supported_variable_traits",
            readme_text.replace(
                "Native p-values test the mean; native BH families contain finite p-values among "
                "selected supported variable traits, not just Table 1. Selected untestable traits "
                "receive q = 1.0 outside the correction vector.",
                "Native p-values test the mean; native BH families contain all supported variable "
                "traits, not just Table 1.",
                1,
            ),
        ),
        (
            "mutant_2_drop_untestable_traits_sentence",
            readme_text.replace(
                " Selected untestable traits receive q = 1.0 outside the correction vector.",
                "",
                1,
            ),
        ),
        (
            "mutant_3_wrong_untestable_q_value",
            readme_text.replace(
                "receive q = 1.0 outside the correction vector.",
                "receive q = 0.0 outside the correction vector.",
                1,
            ),
        ),
        (
            "mutant_4_inside_instead_of_outside_correction_vector",
            readme_text.replace(
                "outside the correction vector.",
                "inside the correction vector.",
                1,
            ),
        ),
        (
            "mutant_5_drop_finite_p_values_qualification",
            readme_text.replace(
                "contain finite p-values among selected supported variable traits",
                "contain selected supported variable traits",
                1,
            ),
        ),
        (
            "mutant_6_uncompressed_comparison_tsv_filename",
            readme_text.replace(
                "native_composition_comparison.tsv.gz",
                "native_composition_comparison.tsv",
                1,
            ),
        ),
        (
            "mutant_7_missing_composition_heading",
            readme_text.replace("## Composition\n\n", "## Summary\n\n", 1),
        ),
    ]

    killed = 0
    for name, mutated_readme in mutants:
        if mutated_readme == readme_text:
            raise RuntimeError(f"Mutant {name} did not alter README text")
        with tempfile.TemporaryDirectory(prefix=f"verify_1635_{name}_") as raw_dir:
            scratch_readme = Path(raw_dir) / "README.md"
            scratch_readme.write_text(mutated_readme, encoding="utf-8")
            try:
                test_module.check_readme_text(scratch_readme.read_text(encoding="utf-8"))
            except AssertionError:
                killed += 1
            else:
                raise RuntimeError(f"Mutant survived: {name}")

    if killed != len(mutants) or killed < 6:
        raise RuntimeError(f"Expected at least 6 killed mutants, got {killed}/{len(mutants)}")
    print(f"All {killed}/{len(mutants)} mutants killed by AssertionError.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
