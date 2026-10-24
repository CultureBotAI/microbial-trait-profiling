"""Compatibility shim.

Forwards bsdbng.ora_calibration.metrics
to mtp.ora_calibration.metrics.
"""

import sys
from typing import Any

from mtp.ora_calibration.metrics import (
    FDRMetric,
    PowerMetric,
    TypeIErrorMetric,
    compute_discrete_ecdf,
    compute_empirical_fdr,
    compute_empirical_power,
    compute_type_i_error,
)

__all__ = [
    "FDRMetric",
    "PowerMetric",
    "TypeIErrorMetric",
    "compute_discrete_ecdf",
    "compute_empirical_fdr",
    "compute_empirical_power",
    "compute_type_i_error",
]

import mtp.ora_calibration.metrics as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.ora_calibration.metrics", run_name="__main__")
