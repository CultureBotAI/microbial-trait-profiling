"""Anaerobe v2 guild scores come from the full profile table, not the heatmap selection (#2399).

`anaerobic` (METPO:1000603) is a declared guild trait but is not among the heatmap's top
traits, so reading its z-score from the heatmap matrix raised KeyError: 'z_anaerobic'.
"""

import importlib.util
import math
import sys
from pathlib import Path
from types import ModuleType

import pandas as pd
import pytest

from mtp.composition_denominators import DENOMINATOR_DEFINITION

SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))


def load_script(name: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(f"guild_2399_{name}", SCRIPTS_DIR / f"{name}.py")
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


anaerobe = load_script("investigate_gut_anaerobe_depletion")

ANAEROBIC = "METPO:1000603"
OBLIGATE = "METPO:1000607"
OTHER = "METPO:1000999"
COUNTS = {
    # profile key: (anaerobic, obligately anaerobic, other)
    "s1|e1|inc": (1, 2, 7),
    "s1|e1|dec": (6, 3, 1),
    "s2|e2|inc": (0, 1, 9),
    "s2|e2|dec": (4, 4, 2),
}


def expected_z(trait_index: int) -> list[float]:
    """Hand-computed heatmap transform: z of log1p(1000 * fraction), population SD."""
    scaled = [math.log1p(1000 * values[trait_index] / sum(values)) for values in COUNTS.values()]
    mean = sum(scaled) / len(scaled)
    sd = math.sqrt(sum((value - mean) ** 2 for value in scaled) / len(scaled))
    return [(value - mean) / sd for value in scaled]


@pytest.fixture
def heatmap_without_anaerobic(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Heatmap selection holds only the obligate trait; the full table holds all three."""
    keys = list(COUNTS)
    wide = pd.DataFrame(
        {
            "study_id": [key.split("|")[0] for key in keys],
            "experiment_id": [key.split("|")[1] for key in keys],
            "inc_dec": [key.split("|")[2] for key in keys],
            ANAEROBIC: [values[0] for values in COUNTS.values()],
            OBLIGATE: [values[1] for values in COUNTS.values()],
            OTHER: [values[2] for values in COUNTS.values()],
        }
    )
    wide_path = tmp_path / "study_trait_profiles_wide.tsv"
    wide.to_csv(wide_path, sep="\t", index=False)

    vis_dir = tmp_path / "vis"
    vis_dir.mkdir()
    pd.DataFrame({"obligate_label": expected_z(1)}, index=keys).to_csv(
        vis_dir / "study_trait_profile_clustered_heatmap_matrix.tsv", sep="\t"
    )
    pd.DataFrame(
        {
            "study_id": wide["study_id"].to_numpy(),
            "experiment_id": wide["experiment_id"].to_numpy(),
            "inc_dec": wide["inc_dec"].to_numpy(),
            "total_trait_signal": [float(sum(values)) for values in COUNTS.values()],
            "denominator_definition": DENOMINATOR_DEFINITION,
        },
        index=keys,
    ).to_csv(vis_dir / "study_trait_profile_study_metadata.tsv", sep="\t")
    pd.DataFrame(
        {
            "trait_id": [OBLIGATE],
            "display_label": ["obligate_label"],
            "object_name": ["obligately anaerobic"],
        },
        index=[OBLIGATE],
    ).to_csv(vis_dir / "study_trait_profile_trait_metadata.tsv", sep="\t")
    monkeypatch.setattr(anaerobe, "VIS_DIR", vis_dir)
    monkeypatch.setattr(anaerobe, "PROFILE_MATRIX_PATH", wide_path)


@pytest.mark.usefixtures("heatmap_without_anaerobic")
def test_guild_trait_outside_heatmap_selection_is_scored_from_full_table() -> None:
    matrix, metadata, _, counts = anaerobe.load_study_heatmap()
    assert "anaerobic" not in matrix.columns

    scored = anaerobe.add_guild_scores(anaerobe.guild_trait_z_scores(counts), metadata)

    assert list(scored["z_anaerobic"]) == pytest.approx(expected_z(0))
    assert list(scored["z_obligately anaerobic"]) == pytest.approx(expected_z(1))
    # A guild trait inside the heatmap selection keeps exactly its heatmap z-score.
    assert list(scored["z_obligately anaerobic"]) == pytest.approx(
        list(matrix["obligately anaerobic"]), abs=1e-12
    )
    both = [(a + b) / 2 for a, b in zip(expected_z(0), expected_z(1), strict=True)]
    assert list(scored["guild_z_mean"]) == pytest.approx(both)

    summary = anaerobe.summarize_direction(scored).set_index("trait")
    assert list(summary.index) == ["anaerobic", "obligately anaerobic", "guild_z_mean"]
    inc = [z for key, z in zip(COUNTS, expected_z(0), strict=True) if key.endswith("inc")]
    dec = [z for key, z in zip(COUNTS, expected_z(0), strict=True) if key.endswith("dec")]
    assert summary.loc["anaerobic", "decreased_minus_increased_z"] == pytest.approx(
        sum(dec) / len(dec) - sum(inc) / len(inc)
    )


def test_missing_guild_trait_in_full_table_is_an_error() -> None:
    counts = pd.DataFrame({OBLIGATE: [1, 2], OTHER: [3, 4]}, index=["a|e|inc", "a|e|dec"])
    with pytest.raises(ValueError, match="METPO:1000603"):
        anaerobe.guild_trait_z_scores(counts)


def test_guild_scores_refuse_misaligned_profiles() -> None:
    guild_z = pd.DataFrame({"anaerobic": [0.0], "obligately anaerobic": [0.0]}, index=["p1"])
    metadata = pd.DataFrame({"inc_dec": ["inc"]}, index=["p2"])
    with pytest.raises(ValueError, match="not aligned"):
        anaerobe.add_guild_scores(guild_z, metadata)
