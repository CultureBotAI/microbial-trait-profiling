"""Compatibility shim.

Forwards bsdbng.ora_calibration.adapters
to mtp.ora_calibration.adapters.
"""

from typing import Any

from mtp.ora_calibration.adapters import (
    DEFAULT_MAX_FOLD,
    DEFAULT_MIN_FOLD,
    TEST_SCOPE_DEPLETION,
    TEST_SCOPE_ENRICHMENT,
    IlyasORAAdapter,
)

__all__ = [
    "DEFAULT_MAX_FOLD",
    "DEFAULT_MIN_FOLD",
    "TEST_SCOPE_DEPLETION",
    "TEST_SCOPE_ENRICHMENT",
    "IlyasORAAdapter",
]

import mtp.ora_calibration.adapters as _target

# Preserve the shim package path so its child module aliases the canonical one.


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.ora_calibration.adapters", run_name="__main__")
