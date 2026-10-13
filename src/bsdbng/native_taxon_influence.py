"""Compatibility shim forwarding bsdbng.native_taxon_influence to mtp.native_taxon_influence."""

import sys
from typing import Any

from mtp.native_taxon_influence import (
    INPUT_NAMES,
    composition_q_sensitivity_table,
    compressed_rows,
    correct_shared_comparisons,
    export_native,
    fold_pairs,
    hypergeom_probabilities,
    main,
    native_summary,
    null_count_interval,
    ora_cells,
    reconstruct,
    run,
    transformed,
)

__all__ = [
    "INPUT_NAMES",
    "composition_q_sensitivity_table",
    "compressed_rows",
    "correct_shared_comparisons",
    "export_native",
    "fold_pairs",
    "hypergeom_probabilities",
    "main",
    "native_summary",
    "null_count_interval",
    "ora_cells",
    "reconstruct",
    "run",
    "transformed",
]

import mtp.native_taxon_influence as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.native_taxon_influence", run_name="__main__")
