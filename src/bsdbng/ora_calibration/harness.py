"""Compatibility shim.

Forwards bsdbng.ora_calibration.harness
to mtp.ora_calibration.harness.
"""

import sys
from typing import Any

from mtp.ora_calibration.harness import (
    ECDF_FILENAME,
    EXIT_FAILURE,
    EXIT_PARTIAL,
    EXIT_SUCCESS,
    NULL_PVALUES_FILENAME,
    REPLICATES_FILENAME,
    SCENARIOS,
    STATUS_FAILED,
    STATUS_PARTIAL,
    STATUS_SUCCESS,
    SUMMARY_FILENAME,
    CalibrationBenchmarkReport,
    CalibrationProvenance,
    DiscoveryThresholdMetrics,
    ReplicateRecord,
    _clear_previous_outputs,
    _load_inputs,
    _positive_int,
    _probability,
    _sha256,
    _software_versions,
    build_parser,
    logger,
    main,
    parse_args,
    run_ora_calibration_benchmark,
    save_calibration_results,
    select_spike_traits,
)

__all__ = [
    "ECDF_FILENAME",
    "EXIT_FAILURE",
    "EXIT_PARTIAL",
    "EXIT_SUCCESS",
    "NULL_PVALUES_FILENAME",
    "REPLICATES_FILENAME",
    "SCENARIOS",
    "STATUS_FAILED",
    "STATUS_PARTIAL",
    "STATUS_SUCCESS",
    "SUMMARY_FILENAME",
    "CalibrationBenchmarkReport",
    "CalibrationProvenance",
    "DiscoveryThresholdMetrics",
    "ReplicateRecord",
    "_clear_previous_outputs",
    "_load_inputs",
    "_positive_int",
    "_probability",
    "_sha256",
    "_software_versions",
    "build_parser",
    "logger",
    "main",
    "parse_args",
    "run_ora_calibration_benchmark",
    "save_calibration_results",
    "select_spike_traits",
]

import mtp.ora_calibration.harness as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.ora_calibration.harness", run_name="__main__")
