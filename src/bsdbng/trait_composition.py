"""Compatibility shim forwarding bsdbng.trait_composition to mtp.trait_composition."""

import sys
from typing import Any

from mtp.trait_composition import (
    ANALYSIS_VERSION,
    CONFIDENCE_LEVEL,
    LOG_RATIO_TIE_TOLERANCE,
    LogRatioResult,
    LogRatioSummary,
    StudyLogRatio,
    _display,
    _matched_profiles,
    _summarize,
    compute_log_ratio_shifts,
    file_sha256,
    main,
    paired_log2_ratio_change,
    write_log_ratio_outputs,
)

__all__ = [
    "ANALYSIS_VERSION",
    "CONFIDENCE_LEVEL",
    "LOG_RATIO_TIE_TOLERANCE",
    "LogRatioResult",
    "LogRatioSummary",
    "StudyLogRatio",
    "_display",
    "_matched_profiles",
    "_summarize",
    "compute_log_ratio_shifts",
    "file_sha256",
    "main",
    "paired_log2_ratio_change",
    "write_log_ratio_outputs",
]

import mtp.trait_composition as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.trait_composition", run_name="__main__")
