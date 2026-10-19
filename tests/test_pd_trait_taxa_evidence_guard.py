"""The PD trait-taxa recipes refuse evidence built from other labels or another lollipop (#644)."""

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import pandas as pd
import pytest

SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))


def load_script(name: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(
        f"evidence_guard_{name}", SCRIPTS_DIR / f"{name}.py"
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


heatmap = load_script("create_pd_six_trait_taxon_heatmap")
exploratory = load_script("create_pd_exploratory_taxa_influence_heatmap")
contributors = load_script("create_pd_trait_contributors")

TRAITS = {"halo": 0.5, "aero": 0.3, "coccus": -0.08, "cellulose": -0.7}


def metadata() -> pd.DataFrame:
    """Two PD fecal profiles and one PD oral profile."""
    frame = pd.DataFrame(
        {
            "profile_id": ["s1|e1|inc", "s1|e1|dec", "s2|e2|inc"],
            "condition": [heatmap.PD_CONDITION] * 3,
            "body_site_group": [heatmap.FECES, heatmap.FECES, "Oral cavity"],
        }
    )
    return frame.set_index("profile_id", drop=False)


def evidence(profile_ids: list[str] | None = None) -> pd.DataFrame:
    profiles = profile_ids or ["s1|e1|inc", "s1|e1|dec"]
    rows = [(trait, profile) for trait in TRAITS for profile in profiles]
    return pd.DataFrame(rows, columns=["trait_id", "profile_id"])


def carriers(positive: set[str] | None = None) -> pd.DataFrame:
    positive = (
        {trait for trait, mean_z in TRAITS.items() if mean_z > 0} if positive is None else positive
    )
    return pd.DataFrame(
        {"trait_id": list(TRAITS), "positive_side_in_lollipop": [t in positive for t in TRAITS]}
    )


def summary() -> pd.DataFrame:
    return pd.DataFrame(
        {"trait_id": list(TRAITS), "increased_minus_decreased_mean_z": list(TRAITS.values())}
    )


def test_matching_evidence_passes() -> None:
    heatmap.check_trait_taxa_evidence(evidence(), carriers(), summary(), metadata())


def test_evidence_profile_outside_current_pd_feces_is_refused() -> None:
    stale = evidence(["s1|e1|inc", "s2|e2|inc"])

    with pytest.raises(ValueError, match=r"1 evidence profiles .*s2\|e2\|inc.*#584"):
        heatmap.check_trait_taxa_evidence(stale, carriers(), summary(), metadata())


def test_changed_positive_side_is_refused() -> None:
    four_positive = carriers({"halo", "aero", "coccus"})

    with pytest.raises(ValueError, match=r"3 positive-side traits.*has 2.*#584"):
        heatmap.check_trait_taxa_evidence(evidence(), four_positive, summary(), metadata())


def test_changed_trait_set_is_refused() -> None:
    sphere = pd.DataFrame({"trait_id": ["sphere"], "increased_minus_decreased_mean_z": [-0.1]})
    five_traits = pd.concat([summary(), sphere], ignore_index=True)

    with pytest.raises(ValueError, match=r"traits differ.*#584"):
        heatmap.check_trait_taxa_evidence(evidence(), carriers(), five_traits, metadata())


@pytest.mark.parametrize("module", [heatmap, exploratory, contributors], ids=lambda m: m.STEM)
def test_recipe_stops_before_writing_when_evidence_is_stale(
    module: ModuleType, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    evidence(["s1|e1|inc", "s2|e2|inc"]).to_csv(
        tmp_path / "pd_six_trait_taxon_profile_evidence.tsv", sep="\t", index=False
    )
    carriers().to_csv(tmp_path / "pd_six_trait_taxon_carriers.tsv", sep="\t", index=False)
    summary().to_csv(tmp_path / "pd_six_trait_direction_summary.tsv", sep="\t", index=False)
    before = sorted(path.name for path in tmp_path.iterdir())
    monkeypatch.setattr(module, "OUTPUT_DIR", tmp_path)
    monkeypatch.setattr(module, "load_profile_metadata", metadata)

    with pytest.raises(ValueError, match="#584"):
        module.main()

    assert sorted(path.name for path in tmp_path.iterdir()) == before
