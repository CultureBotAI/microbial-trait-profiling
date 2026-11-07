"""Keep the four stratified ORA entries explicit about matched taxa counts (#1839)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest


@pytest.mark.parametrize("stratum", ["pooled", "feces", "feces-excluded", "custom-site"])
def test_stratified_ora_selection_defines_positive_taxa_count(stratum: str) -> None:
    """Matched taxa with no retained annotation still belong to the profile draw."""
    path = Path(__file__).resolve().parents[1] / "methods/multiple_testing_families.yaml"
    registry = json.loads(path.read_text(encoding="utf-8"))
    family_id = f"mt-core-stratified-ora-{stratum}"
    matches = [family for family in registry["families"] if family["family_id"] == family_id]
    assert len(matches) == 1, family_id
    selection = matches[0]["size_source"]["selection_rule"]
    assert (
        "positive taxa_count (the number of distinct matrix-matched reported taxa, "
        "including matched taxa with no retained annotation)"
    ) in selection, family_id
