"""Compatibility shim.

Forwards bsdbng.ora_calibration.scenarios
to mtp.ora_calibration.scenarios.
"""

import sys
from typing import Any

from mtp.ora_calibration.scenarios import (
    PROFILE_ID_COLUMNS,
    PROFILE_METADATA_COLUMNS,
    TAXA_METADATA_COLUMNS,
    _parse_strict_count,
    flip_profile_dataframe,
    generate_complete_null_ora_profiles,
    generate_mixed_signal_ora_profiles,
    generate_study_flip_null_ora_profiles,
    profile_identifiers,
    resolve_signature_sizes,
    shared_trait_columns,
    signature_size_column,
    trait_population_counts,
)

__all__ = [
    "PROFILE_ID_COLUMNS",
    "PROFILE_METADATA_COLUMNS",
    "TAXA_METADATA_COLUMNS",
    "_parse_strict_count",
    "flip_profile_dataframe",
    "generate_complete_null_ora_profiles",
    "generate_mixed_signal_ora_profiles",
    "generate_study_flip_null_ora_profiles",
    "profile_identifiers",
    "resolve_signature_sizes",
    "shared_trait_columns",
    "signature_size_column",
    "trait_population_counts",
]

import mtp.ora_calibration.scenarios as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.ora_calibration.scenarios", run_name="__main__")
