"""Anaerobe v2 guild scores come from the full profile table, not the heatmap selection (#2399).

`anaerobic` (METPO:1000603) is a declared guild trait but is not among the heatmap's top
traits, so reading its z-score from the heatmap matrix raised KeyError: 'z_anaerobic', and
the per-trait outputs (trait-enriched counts, inferred evidence, trait definitions and the
selection table) silently covered only the obligate trait.
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
MEDIA = "mediadive.medium:1"
COUNTS = {
    # profile key: (anaerobic, obligately anaerobic, other)
    "s1|e1|inc": (1, 2, 7),
    "s1|e1|dec": (6, 3, 1),
    "s2|e2|inc": (0, 1, 9),
    "s2|e2|dec": (4, 4, 2),
    "s3|e3|inc": (0, 0, 0),  # zero retained signal: all-zero fractions, still standardized
}
# MediaDive counts are outside the retained denominator and must not change any value.
MEDIA_COUNTS = [50, 0, 7, 900, 30]


def fraction(values: tuple[int, int, int], trait_index: int) -> float:
    total = sum(values)
    return values[trait_index] / total if total else 0.0


def expected_z(trait_index: int) -> list[float]:
    """Hand-computed heatmap transform: z of log1p(1000 * fraction), population SD."""
    scaled = [math.log1p(1000 * fraction(values, trait_index)) for values in COUNTS.values()]
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
            MEDIA: MEDIA_COUNTS,
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
            "condition": ["Crohn's disease"] * len(keys),
            "total_trait_signal": [float(sum(values)) for values in COUNTS.values()],
            "denominator_definition": DENOMINATOR_DEFINITION,
        },
        index=keys,
    ).to_csv(vis_dir / "study_trait_profile_study_metadata.tsv", sep="\t")
    pd.DataFrame(
        {
            "trait_id": [OBLIGATE],
            "profile_support": [4],
            "total_signal": [10.0],
            "composition_variance": [0.01],
            "display_label": ["obligate_label"],
            "column_label": ["obligately anaerobic"],
            "object_name": ["obligately anaerobic"],
            "object_category": ["biolink:PhenotypicQuality"],
            "display_category": ["PhenotypicQuality"],
            "selection_score": [0.02],
        },
        index=[OBLIGATE],
    ).to_csv(vis_dir / "study_trait_profile_trait_metadata.tsv", sep="\t")

    derived = tmp_path / "data" / "derived"
    derived.mkdir(parents=True)
    columns_path = derived / "taxa_trait_matrix_columns.tsv"
    pd.DataFrame(
        {
            "column_id": [ANAEROBIC, OBLIGATE, OTHER],
            "column_label": ["anaerobic", "obligately anaerobic", "other"],
            "object_name": ["anaerobic", "obligately anaerobic", "other"],
            "object_category": ["biolink:PhenotypicQuality"] * 3,
        }
    ).to_csv(columns_path, sep="\t", index=False)
    pd.DataFrame(
        {
            "taxon_id": ["NCBITaxon:1", "NCBITaxon:2", "NCBITaxon:3"],
            "taxon_name": ["a", "b", "c"],
            "taxonomic_rank": ["species"] * 3,
            ANAEROBIC: [1, 1, 0],
            OBLIGATE: [1, 0, 0],
            OTHER: [0, 1, 1],
        }
    ).to_csv(derived / "taxa_trait_matrix.tsv", sep="\t", index=False)

    monkeypatch.setattr(anaerobe, "VIS_DIR", vis_dir)
    monkeypatch.setattr(anaerobe, "PROFILE_MATRIX_PATH", wide_path)
    monkeypatch.setattr(anaerobe, "TRAIT_COLUMNS_PATH", columns_path)
    monkeypatch.setattr(anaerobe, "REPO_ROOT", tmp_path)


@pytest.mark.usefixtures("heatmap_without_anaerobic")
def test_guild_trait_outside_heatmap_selection_is_scored_from_full_table() -> None:
    matrix, metadata, _, counts = anaerobe.load_study_heatmap()
    assert "anaerobic" not in matrix.columns

    scored = anaerobe.add_guild_scores(anaerobe.guild_trait_z_scores(counts), metadata)

    assert list(scored["z_anaerobic"]) == pytest.approx(expected_z(0))
    assert list(scored["z_obligately anaerobic"]) == pytest.approx(expected_z(1))
    both = [(a + b) / 2 for a, b in zip(expected_z(0), expected_z(1), strict=True)]
    assert list(scored["guild_z_mean"]) == pytest.approx(both)
    either_positive = [a > 0 or b > 0 for a, b in zip(expected_z(0), expected_z(1), strict=True)]
    assert list(scored["has_enriched_guild_trait"]) == either_positive

    summary = anaerobe.summarize_direction(scored).set_index("trait")
    assert list(summary.index) == ["anaerobic", "obligately anaerobic", "guild_z_mean"]
    inc = [z for key, z in zip(COUNTS, expected_z(0), strict=True) if key.endswith("inc")]
    dec = [z for key, z in zip(COUNTS, expected_z(0), strict=True) if key.endswith("dec")]
    assert summary.loc["anaerobic", "decreased_minus_increased_z"] == pytest.approx(
        sum(dec) / len(dec) - sum(inc) / len(inc)
    )


@pytest.mark.usefixtures("heatmap_without_anaerobic")
def test_every_per_trait_output_lists_both_guild_traits() -> None:
    """Trait-enriched counts, inference metadata and the selection table keep "anaerobic"."""
    matrix, metadata, trait_metadata, counts = anaerobe.load_study_heatmap()
    guild_z = anaerobe.guild_trait_z_scores(counts)
    analysis_matrix, analysis_metadata = anaerobe.extend_with_guild_traits(
        matrix, trait_metadata, counts, guild_z
    )
    assert list(analysis_matrix["anaerobic"]) == pytest.approx(expected_z(0))
    assert list(analysis_matrix["obligately anaerobic"]) == list(matrix["obligately anaerobic"])
    assert dict(analysis_metadata["heatmap_selected"]) == {OBLIGATE: True, ANAEROBIC: False}

    row = analysis_metadata.loc[ANAEROBIC]
    anaerobic_counts = [values[0] for values in COUNTS.values()]
    fractions = [fraction(values, 0) for values in COUNTS.values()]
    mean = sum(fractions) / len(fractions)
    sample_variance = sum((value - mean) ** 2 for value in fractions) / (len(fractions) - 1)
    assert row["profile_support"] == sum(count > 0 for count in anaerobic_counts)
    assert row["total_signal"] == float(sum(anaerobic_counts))
    assert row["composition_variance"] == pytest.approx(sample_variance)
    assert row["object_name"] == "anaerobic"

    scored = anaerobe.add_guild_scores(guild_z, metadata)
    enriched = anaerobe.summarize_trait_enriched_profile_counts(
        analysis_matrix, scored, analysis_metadata
    )
    assert set(enriched.loc[enriched["is_guild_trait"], "trait_id"]) == {ANAEROBIC, OBLIGATE}

    inference = anaerobe.trait_metadata_for_inference(analysis_matrix, analysis_metadata)
    assert set(inference["trait_id"]) >= {ANAEROBIC, OBLIGATE}

    selection = anaerobe.build_guild_trait_selection(analysis_metadata).set_index("trait_id")
    assert selection.loc[ANAEROBIC, "profile_support"] == row["profile_support"]
    assert selection.loc[ANAEROBIC, "kg_microbe_taxa_with_trait_count"] == 2
    stats = selection[["profile_support", "total_signal", "composition_variance"]]
    assert not stats.isna().to_numpy().any()


@pytest.mark.usefixtures("heatmap_without_anaerobic")
def test_extension_is_a_no_op_when_every_guild_trait_is_in_the_heatmap() -> None:
    matrix, _, trait_metadata, counts = anaerobe.load_study_heatmap()
    guild_z = anaerobe.guild_trait_z_scores(counts)
    with_both = trait_metadata.copy()
    with_both.loc[ANAEROBIC] = with_both.loc[OBLIGATE]
    with_both.loc[ANAEROBIC, "trait_id"] = ANAEROBIC
    analysis_matrix, analysis_metadata = anaerobe.extend_with_guild_traits(
        matrix, with_both, counts, guild_z
    )
    assert analysis_matrix.equals(matrix)
    assert analysis_metadata["heatmap_selected"].all()


def test_missing_guild_trait_in_full_table_is_an_error() -> None:
    counts = pd.DataFrame({OBLIGATE: [1, 2], OTHER: [3, 4]}, index=["a|e|inc", "a|e|dec"])
    with pytest.raises(ValueError, match="METPO:1000603"):
        anaerobe.guild_trait_z_scores(counts)


def test_guild_scores_refuse_misaligned_profiles() -> None:
    guild_z = pd.DataFrame({"anaerobic": [0.0], "obligately anaerobic": [0.0]}, index=["p1"])
    metadata = pd.DataFrame({"inc_dec": ["inc"]}, index=["p2"])
    with pytest.raises(ValueError, match="not aligned"):
        anaerobe.add_guild_scores(guild_z, metadata)
