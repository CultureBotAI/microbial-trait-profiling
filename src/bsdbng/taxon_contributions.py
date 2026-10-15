"""Compatibility shim forwarding bsdbng.taxon_contributions to mtp.taxon_contributions."""

import sys
from typing import Any

from mtp.taxon_contributions import (
    PROFILE_FIELDS,
    Ledger,
    Record,
    assertion_fields,
    export_ledger,
    main,
    profile_record,
    reconcile,
    write_table,
)

__all__ = [
    "PROFILE_FIELDS",
    "Ledger",
    "Record",
    "assertion_fields",
    "export_ledger",
    "main",
    "profile_record",
    "reconcile",
    "write_table",
]

import mtp.taxon_contributions as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.taxon_contributions", run_name="__main__")
