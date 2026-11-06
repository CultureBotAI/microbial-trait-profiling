"""Keep the three derived-BY reasons aligned with the anchor guard (#1580).

Only the registry and two named source files are read. No analysis module is
imported, and the observed-output limitation is distinguished from a visible
calculation elsewhere in the shared function.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pytest
import yaml

from test_multiple_testing_registry import discover

ROOT = Path(__file__).resolve().parents[1]
ADAPTER = "src/mtp/composition_calibration/adapters/jed_composition.py"
COMPOSITION = "src/mtp/composition_calibration/composition.py"
ANCHORS = {
    "mt-core-three-test-catalog-t": (ADAPTER, "JedCompositionAdapter.run_inference"),
    "mt-core-calibration-composition-observed": (COMPOSITION, "run_composition_calibration"),
    "mt-core-calibration-composition-null": (COMPOSITION, "run_composition_calibration"),
}


@dataclass(frozen=True)
class DerivedReason:
    path: str
    function: str
    text: str


@pytest.fixture(scope="module")
def reasons() -> dict[str, DerivedReason]:
    registry: dict[str, Any] = yaml.safe_load(
        (ROOT / "methods/multiple_testing_families.yaml").read_text(encoding="utf-8")
    )
    families = {family["family_id"]: family for family in registry["families"]}
    result = {}
    for family_id, (path, function) in ANCHORS.items():
        matches = [
            anchor
            for anchor in families[family_id]["anchors"]
            if anchor["path"] == path and anchor["function"] == function
        ]
        assert len(matches) == 1
        text = matches[0]["derived_procedures"]["BY"]
        assert isinstance(text, str)
        result[family_id] = DerivedReason(path, function, text)
    return result


@pytest.mark.parametrize("family_id", list(ANCHORS))
def test_derived_reason_names_the_visible_formula_check_and_source_review_limit(
    reasons: dict[str, DerivedReason], family_id: str
) -> None:
    reason = reasons[family_id]
    sites, _ = discover({reason.path: (ROOT / reason.path).read_text(encoding="utf-8")})
    assert "BY" in sites[(reason.path, reason.function)].get("derived_procedures", set())
    prose = " ".join(reason.text.casefold().split())
    assert "the guard checks the declaration, not the formula" not in prose
    assert "the guard checks a visible bh * h_m -> clip[0, 1] calculation" in prose
    assert "anchor" in prose
    assert "per-output link is source-reviewed" in prose


def test_observed_reason_names_the_null_witness_and_unchecked_observed_arithmetic(
    reasons: dict[str, DerivedReason],
) -> None:
    reason = reasons["mt-core-calibration-composition-observed"]
    prose = " ".join(reason.text.casefold().split())
    assert "shared anchor" in prose
    assert "null-replicate calculation" in prose
    assert "observed arithmetic is not checked" in prose
    assert "adapter's q_value column" in prose


def test_shared_anchor_witness_does_not_independently_check_observed_output() -> None:
    """Pin the actual limitation described by the observed family's reason."""
    source = (ROOT / COMPOSITION).read_text(encoding="utf-8")
    observed = "obs_by = np.clip(obs_bh * harmonic_sum, 0.0, 1.0)"
    null = "by_r = np.clip(bh_r * harmonic_sum, 0.0, 1.0)"
    assert source.count(observed) == source.count(null) == 1
    without_observed, _ = discover({COMPOSITION: source.replace(observed, "obs_by = obs_bh")})
    without_null, _ = discover({COMPOSITION: source.replace(null, "by_r = bh_r")})
    key = (COMPOSITION, "run_composition_calibration")
    assert "BY" in without_observed[key].get("derived_procedures", set())
    assert "BY" not in without_null[key].get("derived_procedures", set())
