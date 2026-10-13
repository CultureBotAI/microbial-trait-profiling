"""Compatibility shim.

Forwards bsdbng.ora_calibration.protocol
to mtp.ora_calibration.protocol.
"""

import sys
from typing import Any

from mtp.ora_calibration.protocol import (
    MethodAdapter,
    MethodInferenceResult,
)

__all__ = [
    "MethodAdapter",
    "MethodInferenceResult",
]

import mtp.ora_calibration.protocol as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.ora_calibration.protocol", run_name="__main__")
