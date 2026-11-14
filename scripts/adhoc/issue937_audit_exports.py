"""CLI wrapper forwarding canonical export audit to mtp.export_audit (#937, #1039).

Inputs: repository Python sources under src/mtp and scripts, including ignored files.
Output: JSON with source hashes and explicit cross-module public-use omissions.
"""

from __future__ import annotations

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


if __name__ == "__main__":
    main()
