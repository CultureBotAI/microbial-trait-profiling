"""Compatibility shim forwarding bsdbng.density_adjustment to mtp.density_adjustment."""

import sys
from typing import Any

from mtp.density_adjustment import (
    DEFAULT_COLUMNS_PATH,
    DEFAULT_OUTPUT_DIR,
    DEFAULT_PROFILES_PATH,
    DENSITY_PROVENANCE_FILENAME,
    SIGNIFICANT_POST_ADJUSTMENT_STATUSES,
    TABLE_1_TRAITS,
    TEST_PAIR_A,
    TEST_PAIR_B,
    TRAIT_LIST_PROVENANCE_FILENAME,
    TWIN_A,
    TWIN_B,
    ConditionalTwinSummary,
    CorrelationMetrics,
    Table1TraitEffect,
    add_trait_list_argument,
    adjustment_status,
    compute_correlation_metrics,
    compute_effective_tests,
    declared_density_traits,
    density_provenance,
    evaluate_conditional_twin_test,
    evaluate_table1_paired_effects,
    family_bh_q_values,
    format_correlation_metric,
    format_metric_diff,
    format_p_and_q,
    format_paired_regression_slope,
    format_paired_regression_test,
    format_profile_regression_slope,
    format_q_value,
    format_report_number,
    get_non_media_trait_columns,
    load_non_media_profile_matrix,
    main,
    nominal_test_count_statement,
    normalize_traits_by_profile_composition,
    q_value_provenance,
    q_value_statement,
    read_non_media_trait_columns,
    regress_out_pc1,
    reported_traits_heading,
    run_density_audit,
    table1_traits_for_columns,
    write_density_provenance,
    write_trait_list_provenance,
)

__all__ = [
    "DEFAULT_COLUMNS_PATH",
    "DEFAULT_OUTPUT_DIR",
    "DEFAULT_PROFILES_PATH",
    "DENSITY_PROVENANCE_FILENAME",
    "SIGNIFICANT_POST_ADJUSTMENT_STATUSES",
    "TABLE_1_TRAITS",
    "TEST_PAIR_A",
    "TEST_PAIR_B",
    "TRAIT_LIST_PROVENANCE_FILENAME",
    "TWIN_A",
    "TWIN_B",
    "ConditionalTwinSummary",
    "CorrelationMetrics",
    "Table1TraitEffect",
    "add_trait_list_argument",
    "adjustment_status",
    "compute_correlation_metrics",
    "compute_effective_tests",
    "declared_density_traits",
    "density_provenance",
    "evaluate_conditional_twin_test",
    "evaluate_table1_paired_effects",
    "family_bh_q_values",
    "format_correlation_metric",
    "format_metric_diff",
    "format_p_and_q",
    "format_paired_regression_slope",
    "format_paired_regression_test",
    "format_profile_regression_slope",
    "format_q_value",
    "format_report_number",
    "get_non_media_trait_columns",
    "load_non_media_profile_matrix",
    "main",
    "nominal_test_count_statement",
    "normalize_traits_by_profile_composition",
    "q_value_provenance",
    "q_value_statement",
    "read_non_media_trait_columns",
    "regress_out_pc1",
    "reported_traits_heading",
    "run_density_audit",
    "table1_traits_for_columns",
    "write_density_provenance",
    "write_trait_list_provenance",
]

import mtp.density_adjustment as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.density_adjustment", run_name="__main__")
