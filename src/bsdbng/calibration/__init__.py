"""Statistical calibration suite for BugSigDB trait inference methods."""

from bsdbng.calibration.adapters.ilyas_ora import IlyasORAAdapter
from bsdbng.calibration.harness import CalibrationBenchmarkReport, run_ora_calibration_benchmark
from bsdbng.calibration.metrics import (
    FDRMetric,
    TypeIErrorMetric,
    compute_discrete_ecdf,
    compute_empirical_fdr,
    compute_type_i_error,
)
from bsdbng.calibration.protocol import MethodAdapter, MethodInferenceResult
from bsdbng.calibration.scenarios import (
    generate_complete_null_ora_profiles,
    generate_mixed_signal_ora_profiles,
)

__all__ = [
    "CalibrationBenchmarkReport",
    "FDRMetric",
    "IlyasORAAdapter",
    "MethodAdapter",
    "MethodInferenceResult",
    "TypeIErrorMetric",
    "compute_discrete_ecdf",
    "compute_empirical_fdr",
    "compute_type_i_error",
    "generate_complete_null_ora_profiles",
    "generate_mixed_signal_ora_profiles",
    "run_ora_calibration_benchmark",
]
