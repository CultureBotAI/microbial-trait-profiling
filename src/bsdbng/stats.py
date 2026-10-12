"""Compatibility shim forwarding bsdbng.stats to mtp.stats."""

import sys
from typing import Any

from mtp.stats import (
    print_stats,
)

__all__ = [
    "print_stats",
]

import mtp.stats as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.stats", run_name="__main__")
