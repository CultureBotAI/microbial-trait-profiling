"""Compatibility shim forwarding bsdbng.density_adjustment to mtp.density_adjustment."""

import sys
from typing import Any

from mtp.density_adjustment import (
    DEFAULT_COLUMNS_PATH,
    DEFAULT_OUTPUT_DIR,
    DEFAULT_PROFILES_PATH,
    TABLE_1_TRAITS,
    TEST_PAIR_A,
    TEST_PAIR_B,
    TWIN_A,
    TWIN_B,
    ConditionalTwinSummary,
    CorrelationMetrics,
    Table1TraitEffect,
    compute_correlation_metrics,
    compute_effective_tests,
    evaluate_conditional_twin_test,
    evaluate_table1_paired_effects,
    get_non_media_trait_columns,
    load_non_media_profile_matrix,
    main,
    normalize_traits_by_profile_composition,
    regress_out_pc1,
    run_density_audit,
)

__all__ = [
    "DEFAULT_COLUMNS_PATH",
    "DEFAULT_OUTPUT_DIR",
    "DEFAULT_PROFILES_PATH",
    "TABLE_1_TRAITS",
    "TEST_PAIR_A",
    "TEST_PAIR_B",
    "TWIN_A",
    "TWIN_B",
    "ConditionalTwinSummary",
    "CorrelationMetrics",
    "Table1TraitEffect",
    "compute_correlation_metrics",
    "compute_effective_tests",
    "evaluate_conditional_twin_test",
    "evaluate_table1_paired_effects",
    "get_non_media_trait_columns",
    "load_non_media_profile_matrix",
    "main",
    "normalize_traits_by_profile_composition",
    "regress_out_pc1",
    "run_density_audit",
]

import mtp.density_adjustment as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.density_adjustment", run_name="__main__")
