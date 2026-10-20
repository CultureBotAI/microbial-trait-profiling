"""Statistical calibration benchmark harness for BugSigDB ORA methods (#343, #348, #351, #389).

Provides Monte Carlo evaluation metrics, synthetic null scenarios, and method adapters
for ORA hypothesis testing and ranking exports.
"""

from importlib import import_module
from typing import TYPE_CHECKING, Any

from mtp.ora_calibration.adapters.ilyas_ora import (
    DEFAULT_MAX_FOLD,
    DEFAULT_MIN_FOLD,
    TEST_SCOPE_DEPLETION,
    TEST_SCOPE_ENRICHMENT,
    IlyasORAAdapter,
)
from mtp.ora_calibration.metrics import (
    FDRMetric,
    PowerMetric,
    TypeIErrorMetric,
    compute_discrete_ecdf,
    compute_empirical_fdr,
    compute_empirical_power,
    compute_type_i_error,
)
from mtp.ora_calibration.protocol import MethodAdapter, MethodInferenceResult
from mtp.ora_calibration.scenarios import (
    flip_profile_dataframe,
    generate_complete_null_ora_profiles,
    generate_mixed_signal_ora_profiles,
    generate_study_flip_null_ora_profiles,
)

if TYPE_CHECKING:
    from mtp.ora_calibration.harness import (
        CalibrationBenchmarkReport,
        run_ora_calibration_benchmark,
        save_calibration_results,
    )


def __getattr__(name: str) -> Any:
    """Preserve package exports without preloading the executable harness module."""
    if name in {
        "CalibrationBenchmarkReport",
        "run_ora_calibration_benchmark",
        "save_calibration_results",
    }:
        return getattr(import_module("mtp.ora_calibration.harness"), name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


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
    "PowerMetric",
    "TypeIErrorMetric",
    "compute_discrete_ecdf",
    "compute_empirical_fdr",
    "compute_empirical_power",
    "compute_type_i_error",
    "flip_profile_dataframe",
    "generate_complete_null_ora_profiles",
    "generate_mixed_signal_ora_profiles",
    "generate_study_flip_null_ora_profiles",
    "run_ora_calibration_benchmark",
    "save_calibration_results",
]
