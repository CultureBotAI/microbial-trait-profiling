"""Compatibility shim.

Forwards bsdbng.ora_calibration.adapters.ilyas_ora
to mtp.ora_calibration.adapters.ilyas_ora.
"""

import sys
from typing import Any

from mtp.ora_calibration.adapters.ilyas_ora import (
    DEFAULT_MAX_FOLD,
    DEFAULT_MIN_FOLD,
    RESULT_COLUMNS,
    TEST_SCOPE_DEPLETION,
    TEST_SCOPE_ENRICHMENT,
    TEST_SCOPE_PRE_SPECIFIED,
    IlyasORAAdapter,
)

__all__ = [
    "DEFAULT_MAX_FOLD",
    "DEFAULT_MIN_FOLD",
    "RESULT_COLUMNS",
    "TEST_SCOPE_DEPLETION",
    "TEST_SCOPE_ENRICHMENT",
    "TEST_SCOPE_PRE_SPECIFIED",
    "IlyasORAAdapter",
]

import mtp.ora_calibration.adapters.ilyas_ora as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.ora_calibration.adapters.ilyas_ora", run_name="__main__")
