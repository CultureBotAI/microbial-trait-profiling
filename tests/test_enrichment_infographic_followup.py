"""The enrichment follow-up names target conditions it cannot summarize (#586)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd

MODULE_PATH = (
    Path(__file__).resolve().parents[1] / "scripts" / "analyze_enrichment_infographic_followup.py"
)
SPEC = importlib.util.spec_from_file_location("enrichment_infographic_followup", MODULE_PATH)
assert SPEC is not None
assert SPEC.loader is not None
followup = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = followup
SPEC.loader.exec_module(followup)


def test_targets_below_the_profile_minimum_are_reported_with_their_counts() -> None:
    present, rare, absent = followup.TARGET_CONDITIONS
    metadata = pd.DataFrame(
        {
            "profile_id": ["p1", "p2", "p3"],
            "condition": [present, present, rare],
            "study_id": ["s1", "s1", "s2"],
            "experiment_id": ["e1", "e2", "e3"],
            "profile_trait_count": [5, 3, 4],
            "total_trait_signal": [10, 6, 8],
        }
    )

    depth = followup.summarize_condition_depth(metadata, min_profiles=2)

    assert depth["condition"].tolist() == [present]
    assert followup.dropped_target_conditions(metadata, depth) == {rare: 1, absent: 0}
