"""Framework adapters for statistical calibration (#348, #351)."""

from bsdbng.ora_calibration.adapters.ilyas_ora import (
    DEFAULT_MAX_FOLD,
    DEFAULT_MIN_FOLD,
    TEST_SCOPE_DEPLETION,
    TEST_SCOPE_ENRICHMENT,
    IlyasORAAdapter,
)

__all__ = [
    "DEFAULT_MAX_FOLD",
    "DEFAULT_MIN_FOLD",
    "TEST_SCOPE_DEPLETION",
    "TEST_SCOPE_ENRICHMENT",
    "IlyasORAAdapter",
]
