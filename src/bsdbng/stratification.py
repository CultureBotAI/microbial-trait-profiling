"""Compatibility shim forwarding bsdbng.stratification to mtp.stratification."""

import sys
from typing import Any

from mtp.stratification import (
    ANNOTATION_DENSITY_NOTE,
    DEFAULT_EXPERIMENTS_PATH,
    DEFAULT_OUTPUT_DIR,
    DEFAULT_PROFILES_PATH,
    DEFAULT_TAXA_MATRIX_PATH,
    PRECAUTIONARY_NOTE,
    PRIMARY_STRATA,
    StratumOutputs,
    build_sensitivity_summary_table,
    classify_site_stratum,
    compute_annotation_density_deciles,
    generate_sensitivity_reports,
    is_feces_site,
    main,
    parse_args,
    partition_profiles,
    run_stratification_pipeline,
    run_stratum_enrichment,
)

__all__ = [
    "ANNOTATION_DENSITY_NOTE",
    "DEFAULT_EXPERIMENTS_PATH",
    "DEFAULT_OUTPUT_DIR",
    "DEFAULT_PROFILES_PATH",
    "DEFAULT_TAXA_MATRIX_PATH",
    "PRECAUTIONARY_NOTE",
    "PRIMARY_STRATA",
    "StratumOutputs",
    "build_sensitivity_summary_table",
    "classify_site_stratum",
    "compute_annotation_density_deciles",
    "generate_sensitivity_reports",
    "is_feces_site",
    "main",
    "parse_args",
    "partition_profiles",
    "run_stratification_pipeline",
    "run_stratum_enrichment",
]

import mtp.stratification as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.stratification", run_name="__main__")
