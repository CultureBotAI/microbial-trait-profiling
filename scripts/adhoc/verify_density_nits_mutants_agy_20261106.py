"""Scratch-copy source-level mutant verifier for density-report review nits.

Purpose:
  Verify that source-level mutants for #1396, #1398, #1399, #1400, #1401, #1402,
  #1405, #1406, #1407, #1514, #1791, and #1792 are killed by assertion failure
  (rc=1) in an isolated temporary copy without touching the live worktree files,
  and that pytest imports the scratch copy (PYTEST_ORIGIN_OK).

Inputs:
  - src/mtp/density_adjustment.py
  - src/mtp/ontology_audit.py
  - scripts/audit_trait_density_correlation.py
  - tests/test_density_adjustment.py
  - tests/test_density_adjustment_q_values.py
  - tests/test_density_consolidated_conclusions.py
  - tests/test_ontology_audit.py

Outputs:
  - Exit 0 if all 11 mutants are killed by assertion failure (rc=1) with
    verified scratch import origin; non-zero otherwise.

Date: 2026-11-06
"""

from __future__ import annotations

import contextlib
import io
import multiprocessing as mp
import os
import shutil
import sys
import tempfile
from pathlib import Path

import pytest

MUTANTS: list[tuple[str, str, str, str, str]] = [
    (
        "N1_list_provenance_unlink_omitted",
        "src/mtp/density_adjustment.py",
        "(output_dir / DENSITY_PROVENANCE_FILENAME).unlink(missing_ok=True)",
        "pass",
        "tests/test_density_adjustment.py::test_declared_trait_list_run_drops_a_stray_default_provenance_record",
    ),
    (
        "N2_default_pc1_fit_traits_mentions_list",
        "src/mtp/density_adjustment.py",
        '"pc1_fit_traits": (\n'
        '            "every non-MediaDive trait column"\n'
        "            if trait_list is None\n"
        '            else "every non-MediaDive trait column, independent of the list"\n'
        "        ),",
        '"pc1_fit_traits": "every non-MediaDive trait column, independent of the list",',
        "tests/test_density_adjustment_q_values.py::test_default_run_records_machine_readable_provenance_matching_trait_list_schema",
    ),
    (
        "N3_mismatched_trait_ids_ignored",
        "src/mtp/density_adjustment.py",
        "if trait_list is not None and ids != list(trait_list.trait_ids):",
        "if False and trait_list is not None and ids != list(trait_list.trait_ids):",
        "tests/test_density_adjustment.py::test_density_provenance_rejects_mismatched_trait_ids_for_declared_list",
    ),
    (
        "N4_cli_provenance_print_omitted",
        "src/mtp/density_adjustment.py",
        "print(f\"Provenance JSON: {results['provenance']}\")",
        "pass",
        "tests/test_density_adjustment.py::test_density_cli_accepts_a_trait_list",
    ),
    (
        "N5_consolidated_fixed_claim_restored",
        "scripts/audit_trait_density_correlation.py",
        '"2. DENSITY AXIS ADJUSTMENT AND CORRELATION CENTERING",',
        '"2. DENSITY AXIS ADJUSTMENT: SEPARATING BIOLOGY FROM CURATION EFFORT",',
        "tests/test_density_consolidated_conclusions.py::test_conclusions_report_profile_and_paired_diagnostics_without_severity_claims",
    ),
    (
        "N6_ontology_clusters_constant_len",
        "src/mtp/ontology_audit.py",
        'int(df_pairs["cluster_id"].nunique())\n'
        '        if "cluster_id" in df_pairs.columns and not df_pairs.empty\n'
        "        else 0",
        "len(clusters)",
        "tests/test_ontology_audit.py::test_ontology_redundancy_audit_counts_observed_clusters_only",
    ),
    (
        "N7_module_report_fixed_claim_restored",
        "src/mtp/density_adjustment.py",
        '"2. DENSITY AXIS ADJUSTMENT AND CORRELATION CENTERING",',
        '"2. DENSITY AXIS ADJUSTMENT: SEPARATING BIOLOGY FROM CURATION EFFORT",',
        "tests/test_density_adjustment.py::test_module_report_omits_fixed_claims_and_formats_unavailable_sentinels_as_na",
    ),
    (
        "N8_sentinel_check_disabled",
        "src/mtp/density_adjustment.py",
        "if metrics.n_traits - metrics.constant_traits_excluded < 2:",
        "if False and metrics.n_traits - metrics.constant_traits_excluded < 2:",
        "tests/test_density_consolidated_conclusions.py::test_conclusions_do_not_treat_unavailable_zero_sentinels_as_observations",
    ),
    (
        "N9_plural_helper_unconditional",
        "scripts/audit_trait_density_correlation.py",
        "return singular if count == 1 else plural",
        "return plural",
        "tests/test_density_consolidated_conclusions.py::test_consolidated_report_pluralizes_singular_and_plural_counts",
    ),
    (
        "N10_profile_slope_and_partial_r_zero_sentinel_restored",
        "src/mtp/density_adjustment.py",
        'slope_prof = float("nan")',
        "slope_prof = 0.0",
        "tests/test_density_adjustment.py::test_evaluate_conditional_twin_test_constant_and_unpaired",
    ),
    (
        "N11_delta_a_rounding_noise_guard_disabled",
        "src/mtp/density_adjustment.py",
        "delta_a_is_constant = std_delta_a <= _EXACT_FIT_RESIDUAL_SHARE * max_abs_delta_a",
        "delta_a_is_constant = False",
        "tests/test_density_adjustment.py::test_evaluate_conditional_twin_test_constant_and_unpaired",
    ),
    (
        "N12_threshold_strict_less_T01",
        "src/mtp/density_adjustment.py",
        "delta_a_is_constant = std_delta_a <= _EXACT_FIT_RESIDUAL_SHARE * max_abs_delta_a",
        "delta_a_is_constant = std_delta_a < _EXACT_FIT_RESIDUAL_SHARE * max_abs_delta_a",
        "tests/test_density_adjustment.py::test_evaluate_conditional_twin_test_constant_and_unpaired",
    ),
    (
        "N13_threshold_too_generous_T02",
        "src/mtp/density_adjustment.py",
        "delta_a_is_constant = std_delta_a <= _EXACT_FIT_RESIDUAL_SHARE * max_abs_delta_a",
        "delta_a_is_constant = std_delta_a <= 1e-6 * max_abs_delta_a",
        "tests/test_density_adjustment.py::test_evaluate_conditional_twin_test_constant_and_unpaired",
    ),
    (
        "N14_threshold_absolute_not_relative_T03",
        "src/mtp/density_adjustment.py",
        "delta_a_is_constant = std_delta_a <= _EXACT_FIT_RESIDUAL_SHARE * max_abs_delta_a",
        "delta_a_is_constant = std_delta_a <= _EXACT_FIT_RESIDUAL_SHARE",
        "tests/test_density_adjustment.py::test_evaluate_conditional_twin_test_constant_and_unpaired",
    ),
    (
        "N15_threshold_scaled_by_delta_b_T04",
        "src/mtp/density_adjustment.py",
        "delta_a_is_constant = std_delta_a <= _EXACT_FIT_RESIDUAL_SHARE * max_abs_delta_a",
        (
            "delta_a_is_constant = std_delta_a <= _EXACT_FIT_RESIDUAL_SHARE * "
            "(float(np.max(np.abs(delta_b))) if n_paired else 0.0)"
        ),
        "tests/test_density_adjustment.py::test_evaluate_conditional_twin_test_constant_and_unpaired",
    ),
    (
        "N16_threshold_scaled_by_delta_a_without_operands_1824",
        "src/mtp/density_adjustment.py",
        "float(max(np.max(np.abs(delta_a)), np.max(np.abs(inc_a)), np.max(np.abs(dec_a))))",
        "float(np.max(np.abs(delta_a)))",
        "tests/test_density_adjustment.py::test_evaluate_conditional_twin_test_constant_and_unpaired",
    ),
    (
        "N17_refuse_only_on_length_P01",
        "src/mtp/density_adjustment.py",
        "if trait_list is not None and ids != list(trait_list.trait_ids):",
        "if trait_list is not None and len(ids) != len(trait_list.trait_ids):",
        "tests/test_density_adjustment.py::test_density_provenance_rejects_mismatched_trait_ids_for_declared_list",
    ),
    (
        "N18_refuse_only_on_set_P02",
        "src/mtp/density_adjustment.py",
        "if trait_list is not None and ids != list(trait_list.trait_ids):",
        "if trait_list is not None and set(ids) != set(trait_list.trait_ids):",
        "tests/test_density_adjustment.py::test_density_provenance_rejects_mismatched_trait_ids_for_declared_list",
    ),
    (
        "N19_clusters_equal_pair_count_K02",
        "src/mtp/ontology_audit.py",
        'int(df_pairs["cluster_id"].nunique())\n'
        '        if "cluster_id" in df_pairs.columns and not df_pairs.empty\n'
        "        else 0",
        "len(df_pairs)",
        "tests/test_ontology_audit.py::test_ontology_redundancy_audit_counts_observed_clusters_only",
    ),
    (
        "N20_positive_inf_prints_na_F02",
        "src/mtp/density_adjustment.py",
        '    if allow_positive_infinity and value == math.inf:\n        return "inf"\n',
        "",
        "tests/test_density_consolidated_conclusions.py::test_conclusions_distinguish_undefined_from_infinite_diagnostics",
    ),
    (
        "N21_metric_diff_nonfinite_unchecked_F05",
        "src/mtp/density_adjustment.py",
        (
            "    if not (math.isfinite(after_val) and math.isfinite(before_val)):\n"
            '        return "NA"\n'
        ),
        "",
        "tests/test_density_consolidated_conclusions.py::test_conclusions_distinguish_undefined_from_infinite_diagnostics",
    ),
    (
        "N22_module_section4_vif_swapped_F07",
        "src/mtp/density_adjustment.py",
        'f"Variance Inflation Factor (VIF): {twin_vif_prof}",',
        'f"Variance Inflation Factor (VIF): {twin_vif_pair}",',
        "tests/test_density_adjustment.py::test_density_audit_reports_a_declared_trait_list_without_changing_shared_rows",
    ),
    (
        "N23_orphan_collapsed_into_zero_clusters_1826",
        "scripts/audit_trait_density_correlation.py",
        "if n_collapsed > 0 and n_clusters == 0",
        "if False and n_collapsed > 0 and n_clusters == 0",
        "tests/test_density_consolidated_conclusions.py::test_consolidated_report_pluralizes_singular_and_plural_counts",
    ),
    (
        "N24_section5_missing_pc1_label_1827",
        "scripts/audit_trait_density_correlation.py",
        '"  collapsed + PC1 regression PC1 variance explained = "',
        '"  collapsed + PC1 regression variance explained = "',
        "tests/test_density_consolidated_conclusions.py::test_conclusions_report_observed_correlations_and_pc1_changes",
    ),
]


def _copy_scratch_tree(repo_root: Path, scratch_root: Path) -> None:
    shutil.copytree(repo_root / "src", scratch_root / "src")
    for rel in (
        "pyproject.toml",
        "scripts/audit_trait_density_correlation.py",
        "tests/test_density_adjustment.py",
        "tests/test_density_adjustment_q_values.py",
        "tests/test_density_consolidated_conclusions.py",
        "tests/test_ontology_audit.py",
        "methods/biological_interpretation_audit.md",
        "methods/trait_coupling_sensitivity.md",
        "docs/density_adjustment_audit.md",
    ):
        dst = scratch_root / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(repo_root / rel, dst)

    conftest = scratch_root / "tests" / "conftest.py"
    conftest.write_text(
        "from pathlib import Path\n"
        "import mtp.density_adjustment as da\n"
        "import mtp.ontology_audit as oa\n\n"
        "def pytest_sessionstart(session):\n"
        "    root = Path(__file__).resolve().parents[1]\n"
        "    da_expected = (root / 'src/mtp/density_adjustment.py').resolve()\n"
        "    oa_expected = (root / 'src/mtp/ontology_audit.py').resolve()\n"
        "    assert Path(da.__file__).resolve() == da_expected\n"
        "    assert Path(oa.__file__).resolve() == oa_expected\n"
        "    print('PYTEST_ORIGIN_OK')\n",
        encoding="utf-8",
    )


def _pytest_worker(scratch_root_str: str, pytest_target: str, out_file_str: str) -> None:
    scratch_root = Path(scratch_root_str)
    os.chdir(scratch_root)
    sys.path.insert(0, str(scratch_root / "src"))
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        rc = int(pytest.main(["-s", pytest_target]))
    Path(out_file_str).write_text(f"{rc}\n{buf.getvalue()}", encoding="utf-8")


def main() -> int:
    repo_root = Path(__file__).resolve().parents[2]
    ctx = mp.get_context("spawn")
    with tempfile.TemporaryDirectory(prefix="density_nits_mutants_") as tmp_dir:
        scratch_root = Path(tmp_dir)
        _copy_scratch_tree(repo_root, scratch_root)
        out_file = scratch_root / "pytest_out.txt"

        for name, rel_path, target, replacement, pytest_target in MUTANTS:
            target_file = scratch_root / rel_path
            original = target_file.read_text(encoding="utf-8")
            if target not in original:
                print(f"ERROR: target snippet not found for {name} in {rel_path}")
                return 2
            mutated = original.replace(target, replacement, 1)
            target_file.write_text(mutated, encoding="utf-8")
            try:
                proc = ctx.Process(
                    target=_pytest_worker,
                    args=(str(scratch_root), pytest_target, str(out_file)),
                )
                proc.start()
                proc.join()
            finally:
                target_file.write_text(original, encoding="utf-8")

            if proc.exitcode != 0 or not out_file.is_file():
                print(f"ERROR: worker crashed for {name} (exitcode={proc.exitcode})")
                return 1
            raw_out = out_file.read_text(encoding="utf-8")
            rc_str, _, combined = raw_out.partition("\n")
            rc = int(rc_str)
            has_origin = "PYTEST_ORIGIN_OK" in combined
            has_failed_summary = " failed" in combined
            has_no_error = (
                "=== ERRORS ===" not in combined
                and " error in " not in combined.lower()
                and " errors in " not in combined.lower()
                and "ImportError" not in combined
            )
            has_assertion = "AssertionError" in combined or "Failed: DID NOT RAISE" in combined
            if rc != 1 or not (
                has_origin and has_failed_summary and has_no_error and has_assertion
            ):
                print(f"SURVIVED OR INVALID KILL: {name} (rc={rc})")
                print(combined)
                return 1
            print(f"KILLED: {name} (rc={rc}, PYTEST_ORIGIN_OK)")

    print(f"All {len(MUTANTS)}/{len(MUTANTS)} density-nit mutants killed by assertion failure.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
