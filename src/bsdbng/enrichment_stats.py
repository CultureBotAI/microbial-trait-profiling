"""Compatibility shim forwarding bsdbng.enrichment_stats to mtp.enrichment_stats."""

import sys
from typing import Any

from mtp.enrichment_stats import (
    _check_hypergeom_params,
    _log_comb,
    _validate_int,
    _validate_p_values,
    benjamini_hochberg,
    compute_fold_enrichment,
    hypergeom_cdf,
    hypergeom_cdf_many,
    hypergeom_cdf_pure_python,
    hypergeom_sf,
    hypergeom_sf_many,
    hypergeom_sf_pure_python,
    hypergeom_test,
)

__all__ = [
    "_check_hypergeom_params",
    "_log_comb",
    "_validate_int",
    "_validate_p_values",
    "benjamini_hochberg",
    "compute_fold_enrichment",
    "hypergeom_cdf",
    "hypergeom_cdf_many",
    "hypergeom_cdf_pure_python",
    "hypergeom_sf",
    "hypergeom_sf_many",
    "hypergeom_sf_pure_python",
    "hypergeom_test",
]

import mtp.enrichment_stats as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.enrichment_stats", run_name="__main__")
