"""Compatibility shim forwarding bsdbng to mtp."""

from typing import Any

from mtp import (
    __version__,
)

__all__ = [
    "__version__",
]


def __getattr__(name: str) -> Any:
    import mtp as _target

    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp", run_name="__main__")
