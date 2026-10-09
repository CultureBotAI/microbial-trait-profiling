"""Statistical calibration suite for BugSigDB trait inference methods (#343, #348, #351).

Re-exports shared study-level direction flip assignments and replicate RNG (#263, #343)
alongside method adapters, Monte Carlo evaluation metrics, and scenarios.
"""

from bsdbng.calibration.adapters.ilyas_ora import (
    TEST_SCOPE_DEPLETION,
    TEST_SCOPE_ENRICHMENT,
    IlyasORAAdapter,
)
from bsdbng.calibration.harness import run_ora_calibration_benchmark, save_calibration_results
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
    generate_study_flip_null_ora_profiles,
)
from bsdbng.calibration.study_flip import (
    SEED,
    SHARED_FLIP_REPLICATES,
    flip_assignments,
    flip_directions,
    flip_profile_dataframe,
    flipped_units,
    replicate_rng,
    unit_keys,
)

__all__ = [
    "SEED",
    "SHARED_FLIP_REPLICATES",
    "TEST_SCOPE_DEPLETION",
    "TEST_SCOPE_ENRICHMENT",
    "FDRMetric",
    "IlyasORAAdapter",
    "MethodAdapter",
    "MethodInferenceResult",
    "TypeIErrorMetric",
    "compute_discrete_ecdf",
    "compute_empirical_fdr",
    "compute_type_i_error",
    "flip_assignments",
    "flip_directions",
    "flip_profile_dataframe",
    "flipped_units",
    "generate_complete_null_ora_profiles",
    "generate_mixed_signal_ora_profiles",
    "generate_study_flip_null_ora_profiles",
    "replicate_rng",
    "run_ora_calibration_benchmark",
    "save_calibration_results",
    "unit_keys",
]
