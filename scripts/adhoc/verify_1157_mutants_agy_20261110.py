"""Mutant verifier for Issue #1157 unified literal stool predicate (2026-11-10).

Purpose:
    Verify that the unit tests in tests/test_cellulose_decr_analysis.py and
    tests/test_matrix_overlap_coverage.py kill >= 6 non-trivial production
    mutants of the shared stool predicate and its call sites by raising
    AssertionError.
Inputs:
    In-memory monkeypatched mutations of mtp.delta_z_derived.tables,
    mtp.cellulose_decr.analysis, mtp.matrix_overlap_coverage, and
    mtp.density_adjustment.
Outputs:
    Per-mutant kill status printed to stdout; non-zero exit if any mutant
    survives.
Issue: #1157
Date: 2026-11-10
"""

from __future__ import annotations

import importlib.util
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from pathlib import Path
from types import ModuleType

import mtp.cellulose_decr.analysis as cellulose_analysis
import mtp.delta_z_derived.tables as delta_z_tables
import mtp.density_adjustment as density_adjustment
import mtp.matrix_overlap_coverage as matrix_overlap_coverage

REPO_ROOT = Path(__file__).resolve().parents[2]


def _load_test_module(name: str, relative_path: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / relative_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load {relative_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


test_cellulose_mod = _load_test_module(
    "test_cellulose_decr_analysis", "tests/test_cellulose_decr_analysis.py"
)
test_overlap_mod = _load_test_module(
    "test_matrix_overlap_coverage", "tests/test_matrix_overlap_coverage.py"
)


@contextmanager
def _patch_attr(target: object, attr: str, replacement: object) -> Iterator[None]:
    original = getattr(target, attr)
    setattr(target, attr, replacement)
    try:
        yield
    finally:
        setattr(target, attr, original)


def _run_target_tests() -> None:
    test_cellulose_mod.test_shared_stool_predicate_reexported_identically_across_modules()
    for val, expected in [
        ("Feces", True),
        ("feces", False),
        ("FECES", False),
        (" Feces ", False),
        ("Feces,Blood", False),
    ]:
        test_cellulose_mod.test_is_stool_body_site_accepts_only_exact_literal_feces(val, expected)
    for val, expected in [
        ("Feces", True),
        ("Feces,Blood", False),
        ("Oral", False),
    ]:
        test_cellulose_mod.test_validate_stool_body_site_accepts_literal_and_non_stool_sites(
            val, expected
        )
    for variant in ("feces", "FECES", " Feces ", "\tFeces\n"):
        test_cellulose_mod.test_validate_stool_body_site_refuses_case_and_whitespace_variants(
            variant
        )
        test_cellulose_mod.test_compartment_of_refuses_case_and_whitespace_feces_variants(variant)
    for site, compartment in [
        ("Feces", cellulose_analysis.FECES),
        ("Feces,Blood", cellulose_analysis.OTHER),
        ("Colon", cellulose_analysis.NON_FECAL_GI),
    ]:
        test_cellulose_mod.test_compartment_of_accepts_literal_feces_and_non_stool_sites(
            site, compartment
        )
    test_cellulose_mod.test_delta_z_tables_use_shared_is_stool_body_site()
    test_overlap_mod.test_parkinsons_stool_predicate_accepts_literal_feces_and_refuses_variants()
    test_overlap_mod.test_cascade_and_pd_overlap_refuse_case_or_whitespace_feces_variants()


def _expect_killed(name: str, mutate_and_run: Callable[[], None]) -> None:
    try:
        mutate_and_run()
    except AssertionError:
        print(f"KILLED: {name}")
        return
    raise RuntimeError(f"SURVIVED: {name}")


def main() -> int:
    # Baseline check: unmutated code passes all target tests
    _run_target_tests()

    # Mutant 1: is_stool_body_site folds case ("feces" -> True)
    def _m1() -> None:
        def _folded(body_site: object) -> bool:
            return isinstance(body_site, str) and body_site.lower() == "feces"

        with (
            _patch_attr(delta_z_tables, "is_stool_body_site", _folded),
            _patch_attr(cellulose_analysis, "is_stool_body_site", _folded),
            _patch_attr(matrix_overlap_coverage, "is_stool_body_site", _folded),
            _patch_attr(density_adjustment, "is_stool_body_site", _folded),
            _patch_attr(test_cellulose_mod, "is_stool_body_site", _folded),
            _patch_attr(test_overlap_mod, "is_stool_body_site", _folded),
        ):
            _run_target_tests()

    _expect_killed("M1_is_stool_case_folded", _m1)

    # Mutant 2: validate_stool_body_site silently returns False on lowercase "feces"
    def _m2() -> None:
        def _silent_false(body_site: object) -> bool:
            return isinstance(body_site, str) and body_site == "Feces"

        with (
            _patch_attr(delta_z_tables, "validate_stool_body_site", _silent_false),
            _patch_attr(cellulose_analysis, "validate_stool_body_site", _silent_false),
            _patch_attr(matrix_overlap_coverage, "validate_stool_body_site", _silent_false),
            _patch_attr(density_adjustment, "validate_stool_body_site", _silent_false),
            _patch_attr(test_cellulose_mod, "validate_stool_body_site", _silent_false),
            _patch_attr(test_overlap_mod, "validate_stool_body_site", _silent_false),
        ):
            _run_target_tests()

    _expect_killed("M2_validate_stool_silent_false_on_lowercase", _m2)

    # Mutant 3: validate_stool_body_site strips whitespace and accepts " Feces "
    def _m3() -> None:
        def _strip_accept(body_site: object) -> bool:
            if isinstance(body_site, str) and body_site.strip() == "Feces":
                return True
            if isinstance(body_site, str) and body_site.strip().lower() == "feces":
                raise ValueError("Refusing lowercase feces (#1157)")
            return False

        with (
            _patch_attr(delta_z_tables, "validate_stool_body_site", _strip_accept),
            _patch_attr(cellulose_analysis, "validate_stool_body_site", _strip_accept),
            _patch_attr(matrix_overlap_coverage, "validate_stool_body_site", _strip_accept),
            _patch_attr(density_adjustment, "validate_stool_body_site", _strip_accept),
            _patch_attr(test_cellulose_mod, "validate_stool_body_site", _strip_accept),
            _patch_attr(test_overlap_mod, "validate_stool_body_site", _strip_accept),
        ):
            _run_target_tests()

    _expect_killed("M3_validate_stool_strips_whitespace", _m3)

    # Mutant 4: compartment_of folds "feces" to FECES instead of calling validate_stool_body_site
    def _m4() -> None:
        def _legacy_compartment_of(body_site: str) -> str:
            if body_site.strip().lower() == "feces":
                return cellulose_analysis.FECES
            return cellulose_analysis.OTHER

        with (
            _patch_attr(cellulose_analysis, "compartment_of", _legacy_compartment_of),
            _patch_attr(test_cellulose_mod, "compartment_of", _legacy_compartment_of),
        ):
            _run_target_tests()

    _expect_killed("M4_compartment_of_case_folds_feces", _m4)

    # Mutant 5: _is_parkinsons_stool_row uses non-raising equality
    def _m5() -> None:
        def _silent_pd_stool(
            condition: object, body_site: object, *, exact_condition: bool = False
        ) -> bool:
            _ = exact_condition
            has_pd = "parkinson's disease" in str(condition).lower()
            return has_pd and isinstance(body_site, str) and body_site == "Feces"

        with (
            _patch_attr(matrix_overlap_coverage, "_is_parkinsons_stool_row", _silent_pd_stool),
            _patch_attr(test_overlap_mod, "_is_parkinsons_stool_row", _silent_pd_stool),
        ):
            _run_target_tests()

    _expect_killed("M5_parkinsons_stool_row_drops_lowercase_without_raising", _m5)

    # Mutant 6: validate_stool_body_site accepts compound "Feces,Blood" substring
    def _m6() -> None:
        def _substring_stool(body_site: object) -> bool:
            if isinstance(body_site, str) and "Feces" in body_site:
                return True
            if isinstance(body_site, str) and body_site.strip().lower() == "feces":
                raise ValueError("Refusing non-canonical Feces (#1157)")
            return False

        with (
            _patch_attr(delta_z_tables, "validate_stool_body_site", _substring_stool),
            _patch_attr(cellulose_analysis, "validate_stool_body_site", _substring_stool),
            _patch_attr(matrix_overlap_coverage, "validate_stool_body_site", _substring_stool),
            _patch_attr(density_adjustment, "validate_stool_body_site", _substring_stool),
            _patch_attr(test_cellulose_mod, "validate_stool_body_site", _substring_stool),
            _patch_attr(test_overlap_mod, "validate_stool_body_site", _substring_stool),
        ):
            _run_target_tests()

    _expect_killed("M6_validate_stool_accepts_compound_feces_substring", _m6)

    # Mutant 7: density_adjustment defines a separate predicate instead of re-exporting shared one
    def _m7() -> None:
        def _separate(body_site: object) -> bool:
            return isinstance(body_site, str) and body_site == "Feces"

        with _patch_attr(density_adjustment, "validate_stool_body_site", _separate):
            _run_target_tests()

    _expect_killed("M7_density_adjustment_separate_predicate_identity", _m7)

    print("All 7 mutants killed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
