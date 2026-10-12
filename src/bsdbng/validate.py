"""Compatibility shim forwarding bsdbng.validate to mtp.validate."""

import sys
from typing import Any

from mtp.validate import (
    validate_all_studies,
)

__all__ = [
    "validate_all_studies",
]

import mtp.validate as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.validate", run_name="__main__")
