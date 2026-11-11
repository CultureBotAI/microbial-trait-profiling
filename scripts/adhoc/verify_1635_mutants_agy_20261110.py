#!/usr/bin/env python3
"""Adversarial mutation verifier for issue #1635 (2026-11-10).

Purpose:
    Verify that tests/test_table1_native_results_readme_1635.py kills at least
    6 distinct README mutations plus a generator-only BH-policy drift mutation
    in isolated scratch copies.
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

import pytest


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
    gen_src = repo_root / "src/mtp/native_taxon_influence.py"

    test_module = _load_module("test_table1_native_results_readme_1635", test_src)
    readme_text = readme_src.read_text(encoding="utf-8")
    gen_text = gen_src.read_text(encoding="utf-8")

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

    # Mutant 8: generator-only drift where native_summary claims padded inside-vector BH
    # while transformed() still computes finite-only outside-vector BH.
    old_gen_lines = (
        (
            '        "Native p-values test the mean; native BH families contain "\n'
            '        "finite p-values among "\n'
        ).replace('contain "\n        "finite', "contain finite")
        + (
            '        "selected supported variable traits, not just Table 1. "\n'
            '        "Selected untestable traits "\n'
        ).replace('Table 1. "\n        "Selected', "Table 1. Selected")
        + ('        "receive q = 1.0 outside the correction vector. "')
    )
    new_gen_lines = (
        '        "Native p-values test the mean; native BH families contain all "\n'
        '        "selected supported "\n'
    ).replace('all "\n        "selected', "all selected") + (
        '        "variable traits (family_selected_traits), not just Table 1, padding selected "\n'
        '        "untestable traits at p = 1.0 (q = 1.0) inside the correction vector. "'
    )
    mutated_gen_text = gen_text.replace(old_gen_lines, new_gen_lines, 1)
    if mutated_gen_text == gen_text:
        raise RuntimeError("Generator drift mutant did not alter native_taxon_influence.py text")
    with tempfile.TemporaryDirectory(prefix="verify_1635_mutant_8_gen_drift_") as raw_dir:
        scratch_gen = Path(raw_dir) / "native_taxon_influence.py"
        scratch_gen.write_text(mutated_gen_text, encoding="utf-8")
        mutated_gen_module = _load_module("mutant_8_native_taxon_influence", scratch_gen)
        monkeypatch = pytest.MonkeyPatch()
        try:
            monkeypatch.setattr(test_module, "native_summary", mutated_gen_module.native_summary)
            try:
                test_module.check_readme_text(readme_text)
            except AssertionError:
                killed += 1
            else:
                raise RuntimeError("Mutant survived: mutant_8_generator_only_bh_wording_drift")
        finally:
            monkeypatch.undo()

    total_expected = len(mutants) + 1
    if killed != total_expected or killed < 8:
        raise RuntimeError(f"Expected {total_expected} killed mutants, got {killed}")
    print(f"All {killed}/{total_expected} mutants killed by AssertionError.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
