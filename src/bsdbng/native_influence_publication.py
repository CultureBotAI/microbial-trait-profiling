"""Compatibility shim forwarding bsdbng.native_influence_publication
to mtp.native_influence_publication.
"""

import sys
from typing import Any

from mtp.native_influence_publication import (
    main,
    publish,
    reconcile_exported_source_counts,
)

__all__ = [
    "main",
    "publish",
    "reconcile_exported_source_counts",
]

import mtp.native_influence_publication as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.native_influence_publication", run_name="__main__")
