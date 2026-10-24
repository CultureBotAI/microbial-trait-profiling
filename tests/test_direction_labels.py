"""Shared direction labels and safe regression-probe ownership (#962)."""

import importlib.util
import sys
from pathlib import Path

import pytest

from mtp.datamodel.workflow_pydantic import DirectionBasis
from mtp.direction_labels import DIRECTION_BASIS_LABELS

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize(
    "script", ["create_taxa_guild_heatmap.py", "create_study_trait_profile_heatmap.py"]
)
def test_figures_share_complete_schema_direction_labels(script: str) -> None:
    """Each figure consumes the same label vocabulary for every workflow basis."""
    spec = importlib.util.spec_from_file_location(
        "direction_label_figure", ROOT / "scripts" / script
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert module.DIRECTION_BASIS_LABELS is DIRECTION_BASIS_LABELS
    assert set(DIRECTION_BASIS_LABELS) == {basis.value for basis in DirectionBasis}
    assert dict(DIRECTION_BASIS_LABELS) == {
        "group_1_vs_group_0": "Group 1 vs Group 0",
        "disease_vs_comparison": "disease vs comparison",
    }


def test_probe_rejects_another_imported_checkout_before_writing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = ROOT / "scripts/adhoc/issue955_probe_direction_forwarding.py"
    spec = importlib.util.spec_from_file_location("direction_probe", path)
    assert spec is not None
    assert spec.loader is not None
    probe = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(probe)
    output = tmp_path / "outputs"
    monkeypatch.setattr(probe, "REPO_ROOT", tmp_path / "foreign-checkout")
    monkeypatch.setattr(sys, "argv", [str(path), "--output-dir", str(output)])
    with pytest.raises(ValueError, match="Imported mtp belongs to another checkout"):
        probe.main()
    assert not output.exists()
