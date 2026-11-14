"""Compatibility shim forwarding bsdbng.export_audit to mtp.export_audit."""

import sys
from typing import Any

from mtp.export_audit import (
    bound_names,
    cross_module_uses,
    declared_exports,
    export_omissions,
    function_arguments,
    import_origin,
    main,
    module_name,
    public_definitions,
    source_inventory,
)

__all__ = [
    "bound_names",
    "cross_module_uses",
    "declared_exports",
    "export_omissions",
    "function_arguments",
    "import_origin",
    "main",
    "module_name",
    "public_definitions",
    "source_inventory",
]

import mtp.export_audit as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.export_audit", run_name="__main__")
