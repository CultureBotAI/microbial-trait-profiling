"""Compatibility shim forwarding bsdbng.ora_calibration to mtp.ora_calibration."""

import sys
from typing import Any

from mtp.ora_calibration import (
    DEFAULT_MAX_FOLD,
    DEFAULT_MIN_FOLD,
    TEST_SCOPE_DEPLETION,
    TEST_SCOPE_ENRICHMENT,
    CalibrationBenchmarkReport,
    FDRMetric,
    IlyasORAAdapter,
    MethodAdapter,
    MethodInferenceResult,
    TypeIErrorMetric,
    compute_discrete_ecdf,
    compute_empirical_fdr,
    compute_type_i_error,
    flip_profile_dataframe,
    generate_complete_null_ora_profiles,
    generate_mixed_signal_ora_profiles,
    generate_study_flip_null_ora_profiles,
    run_ora_calibration_benchmark,
    save_calibration_results,
)

__all__ = [
    "DEFAULT_MAX_FOLD",
    "DEFAULT_MIN_FOLD",
    "TEST_SCOPE_DEPLETION",
    "TEST_SCOPE_ENRICHMENT",
    "CalibrationBenchmarkReport",
    "FDRMetric",
    "IlyasORAAdapter",
    "MethodAdapter",
    "MethodInferenceResult",
    "TypeIErrorMetric",
    "compute_discrete_ecdf",
    "compute_empirical_fdr",
    "compute_type_i_error",
    "flip_profile_dataframe",
    "generate_complete_null_ora_profiles",
    "generate_mixed_signal_ora_profiles",
    "generate_study_flip_null_ora_profiles",
    "run_ora_calibration_benchmark",
    "save_calibration_results",
]

import mtp.ora_calibration as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.ora_calibration", run_name="__main__")
