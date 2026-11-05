"""The PD profile paired test keeps zero contrasts, as scipy and main do (#1537)."""

from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import pandas as pd
import pytest
from scipy import stats

MODULE_PATH = (
    Path(__file__).resolve().parents[1] / "scripts/analyze_parkinsons_disease_trait_profiles.py"
)
SPEC = importlib.util.spec_from_file_location("pd_profile_analysis", MODULE_PATH)
assert SPEC is not None
assert SPEC.loader is not None
profiles = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = profiles
SPEC.loader.exec_module(profiles)


def test_safe_wilcoxon_passes_zero_contrasts_to_scipy() -> None:
    values = pd.Series([0.0, 0.0, 0.0, *[0.1 * (i + 1) for i in range(27)]])
    expected = float(stats.wilcoxon(values, alternative="two-sided").pvalue)
    without_zeros = float(stats.wilcoxon(values[values != 0], alternative="two-sided").pvalue)
    assert expected != pytest.approx(without_zeros)
    assert profiles.safe_wilcoxon(values) == pytest.approx(expected)


def test_safe_wilcoxon_is_undefined_for_all_zero_or_empty_contrasts() -> None:
    assert math.isnan(profiles.safe_wilcoxon(pd.Series([0.0, 0.0, 0.0])))
    assert math.isnan(profiles.safe_wilcoxon(pd.Series([float("nan")])))
