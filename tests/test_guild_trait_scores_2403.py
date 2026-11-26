"""Every guild score is computed from all of the guild's declared traits (#2403).

The maps, the Parkinson's disease analysis and the guild-assignment grid used to keep only
the guild traits found among the study heatmap's top-trait columns. "anaerobic"
(METPO:1000603) is not among them, so the anaerobe guild score silently became an
obligate-only score. These tests build a profile table whose heatmap holds neither guild
trait and check that each consumer scores both traits from the full table.
"""

import dataclasses
import importlib.util
import math
import sys
from pathlib import Path
from types import ModuleType

import pandas as pd
import pytest

from mtp.composition_denominators import DENOMINATOR_DEFINITION, profile_trait_z_scores
from mtp.guild_trait_scores import guild_trait_z_scores, read_guild_trait_selection

SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))


def load_script(name: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(f"guild_2403_{name}", SCRIPTS_DIR / f"{name}.py")
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


guild_maps = load_script("create_condition_body_site_guild_maps")
parkinsons = load_script("analyze_parkinsons_disease_trait_profiles")
assignment = load_script("create_profile_enrichment_guild_assignment_diagram")

ANAEROBIC = "METPO:1000603"
OBLIGATE = "METPO:1000607"
AEROBIC = "METPO:1000602"
MICROAEROPHILIC = "METPO:1000604"
OTHER = "METPO:1000999"
MEDIA = "mediadive.medium:1"
TRAITS = (ANAEROBIC, OBLIGATE, AEROBIC, MICROAEROPHILIC, OTHER)
COUNTS = {
    # profile key: counts in TRAITS order
    "s1|e1|inc": (1, 2, 5, 1, 7),
    "s1|e1|dec": (6, 3, 0, 2, 1),
    "s2|e2|inc": (0, 1, 4, 3, 9),
    "s2|e2|dec": (4, 4, 1, 0, 2),
    "s3|e3|inc": (0, 0, 0, 0, 0),  # zero retained signal
    "s3|e3|dec": (2, 0, 3, 1, 5),
}
# MediaDive counts are outside the retained denominator and must not change any value.
MEDIA_COUNTS = [50, 0, 7, 900, 30, 2]
# A profile in the full table but not in the heatmap's metadata: it must not enter the
# standardization population (#2415).
OUTSIDE_PROFILE = ("s9|e9|inc", (90, 1, 0, 0, 1), 0)
CONDITION = "Crohn's disease"


def expected_z(trait: str) -> list[float]:
    """Hand-computed heatmap transform: z of log1p(1000 * fraction), population SD."""
    index = TRAITS.index(trait)
    fractions = [values[index] / sum(values) if sum(values) else 0.0 for values in COUNTS.values()]
    scaled = [math.log1p(1000 * value) for value in fractions]
    mean = sum(scaled) / len(scaled)
    sd = math.sqrt(sum((value - mean) ** 2 for value in scaled) / len(scaled))
    return [(value - mean) / sd for value in scaled]


def expected_guild_mean(*traits: str) -> list[float]:
    columns = [expected_z(trait) for trait in traits]
    return [sum(values) / len(values) for values in zip(*columns, strict=True)]


def write_selection(path: Path, traits: dict[str, str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        {
            "guild": ["test guild"] * len(traits),
            "trait_id": list(traits),
            "trait_name": list(traits.values()),
            "profile_matrix_object_name": list(traits.values()),
        }
    ).to_csv(path, sep="\t", index=False)
    return path


@pytest.fixture
def profiles(tmp_path: Path) -> tuple[Path, Path]:
    """Write the full profile table and the heatmap's profile metadata.

    The full table holds one extra profile and lists its rows in reverse order, so
    scores must be standardized over the metadata's profiles and aligned by key.
    """
    keys = list(COUNTS)
    outside_key, outside_counts, outside_media = OUTSIDE_PROFILE
    wide_rows = [*zip(keys, COUNTS.values(), MEDIA_COUNTS, strict=True)]
    wide_rows = [*reversed(wide_rows), (outside_key, outside_counts, outside_media)]
    wide_path = tmp_path / "study_trait_profiles_wide.tsv"
    pd.DataFrame(
        {
            "study_id": [key.split("|")[0] for key, _, _ in wide_rows],
            "experiment_id": [key.split("|")[1] for key, _, _ in wide_rows],
            "inc_dec": [key.split("|")[2] for key, _, _ in wide_rows],
            **{
                trait: [values[index] for _, values, _ in wide_rows]
                for index, trait in enumerate(TRAITS)
            },
            MEDIA: [media for _, _, media in wide_rows],
        }
    ).to_csv(wide_path, sep="\t", index=False)

    vis_dir = tmp_path / "vis"
    vis_dir.mkdir()
    metadata_path = vis_dir / "study_trait_profile_study_metadata.tsv"
    pd.DataFrame(
        {
            "study_id": [key.split("|")[0] for key in keys],
            "experiment_id": [key.split("|")[1] for key in keys],
            "inc_dec": [key.split("|")[2] for key in keys],
            "condition": [CONDITION] * len(keys),
            "body_site": ["feces"] * len(keys),
            "body_site_group": ["Feces"] * len(keys),
            "condition_group": ["Digestive"] * len(keys),
            "total_trait_signal": [float(sum(values)) for values in COUNTS.values()],
            "denominator_definition": DENOMINATOR_DEFINITION,
        },
        index=keys,
    ).to_csv(metadata_path, sep="\t")
    return wide_path, metadata_path


def test_read_guild_trait_selection_returns_ids_and_names(tmp_path: Path) -> None:
    path = write_selection(
        tmp_path / "sel.tsv", {ANAEROBIC: "anaerobic", OBLIGATE: "obligately anaerobic"}
    )
    assert read_guild_trait_selection(path) == {
        ANAEROBIC: "anaerobic",
        OBLIGATE: "obligately anaerobic",
    }


@pytest.mark.parametrize(
    ("frame", "message"),
    [
        (pd.DataFrame({"trait_name": ["anaerobic"]}), "no trait_id column"),
        (pd.DataFrame({"trait_id": [ANAEROBIC], "trait_name": [""]}), "blank"),
        (
            pd.DataFrame({"trait_id": [ANAEROBIC, ANAEROBIC], "trait_name": ["a", "b"]}),
            "more than once",
        ),
        (pd.DataFrame({"trait_id": [], "trait_name": []}), "selects no traits"),
    ],
)
def test_read_guild_trait_selection_rejects_malformed_tables(
    tmp_path: Path, frame: pd.DataFrame, message: str
) -> None:
    path = tmp_path / "sel.tsv"
    frame.to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError, match=message):
        read_guild_trait_selection(path)


def test_shared_transform_matches_hand_computed_heatmap_z(profiles: tuple[Path, Path]) -> None:
    wide_path, metadata_path = profiles
    metadata = pd.read_csv(metadata_path, sep="\t", index_col=0)
    counts = guild_maps.load_retained_profiles(wide_path, metadata)
    assert list(counts.index) == list(COUNTS)
    z_scores = profile_trait_z_scores(counts, [ANAEROBIC, OBLIGATE])
    assert list(z_scores[ANAEROBIC]) == pytest.approx(expected_z(ANAEROBIC))
    assert list(z_scores[OBLIGATE]) == pytest.approx(expected_z(OBLIGATE))

    named = guild_trait_z_scores(counts, {ANAEROBIC: "anaerobic"})
    assert list(named.columns) == ["anaerobic"]
    with pytest.raises(ValueError, match="METPO:1000605"):
        profile_trait_z_scores(counts, ["METPO:1000605"])


def test_condition_body_site_maps_score_every_selected_trait(
    tmp_path: Path, profiles: tuple[Path, Path]
) -> None:
    wide_path, metadata_path = profiles
    metadata = pd.read_csv(metadata_path, sep="\t", index_col=0)
    counts = guild_maps.load_retained_profiles(wide_path, metadata)
    output_dir = tmp_path / "anaerobe-depletion-v2"
    write_selection(
        output_dir / "guild_trait_selection.tsv",
        {ANAEROBIC: "anaerobic", OBLIGATE: "obligately anaerobic"},
    )
    config = dataclasses.replace(guild_maps.ANAEROBE_CONFIG, output_dir=output_dir)

    scored = guild_maps.add_guild_scores(counts, metadata, config)

    assert list(scored["guild_z_mean"]) == pytest.approx(expected_guild_mean(ANAEROBIC, OBLIGATE))
    positives = [
        sum(z > 0 for z in pair)
        for pair in zip(expected_z(ANAEROBIC), expected_z(OBLIGATE), strict=True)
    ]
    assert list(scored["guild_positive_trait_count"]) == positives


def test_parkinsons_scores_record_every_selected_trait(
    tmp_path: Path, profiles: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    wide_path, metadata_path = profiles
    monkeypatch.setattr(parkinsons, "METADATA_PATH", metadata_path)
    monkeypatch.setattr(parkinsons, "PROFILE_MATRIX_PATH", wide_path)
    aerobe = dataclasses.replace(
        parkinsons.AEROBE_CONFIG,
        trait_selection_path=write_selection(
            tmp_path / "aerobe.tsv", {AEROBIC: "aerobic", MICROAEROPHILIC: "microaerophilic"}
        ),
    )
    anaerobe = dataclasses.replace(
        parkinsons.ANAEROBE_CONFIG,
        trait_selection_path=write_selection(
            tmp_path / "anaerobe.tsv", {ANAEROBIC: "anaerobic", OBLIGATE: "obligately anaerobic"}
        ),
    )
    monkeypatch.setattr(parkinsons, "ANALYSIS_CONFIGS", (aerobe, anaerobe))

    metadata, counts = parkinsons.load_profiles()
    scored = parkinsons.add_trait_scores(counts, metadata)

    assert list(scored["anaerobe_v2_score"]) == pytest.approx(
        expected_guild_mean(ANAEROBIC, OBLIGATE)
    )
    assert list(scored["aerobe_v2_score"]) == pytest.approx(
        expected_guild_mean(AEROBIC, MICROAEROPHILIC)
    )
    assert set(scored["anaerobe_depletion_v2_selected_trait_count"]) == {2}
    assert set(scored["aerobe_enrichment_v2_selected_trait_count"]) == {2}


def test_parkinsons_profiles_require_the_denominator_marker(
    tmp_path: Path, profiles: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    wide_path, metadata_path = profiles
    stale = pd.read_csv(metadata_path, sep="\t", index_col=0).drop(columns="denominator_definition")
    stale_path = tmp_path / "stale_metadata.tsv"
    stale.to_csv(stale_path, sep="\t")
    monkeypatch.setattr(parkinsons, "METADATA_PATH", stale_path)
    monkeypatch.setattr(parkinsons, "PROFILE_MATRIX_PATH", wide_path)
    with pytest.raises(ValueError, match="denominator marker"):
        parkinsons.load_profiles()


def test_guild_assignment_grid_scores_every_guild_trait(
    tmp_path: Path, profiles: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    wide_path, metadata_path = profiles
    grid_path = tmp_path / "grid.tsv"
    monkeypatch.setattr(assignment, "VIS_DIR", metadata_path.parent)
    monkeypatch.setattr(assignment, "STUDY_TRAIT_PROFILE_WIDE_TSV", wide_path)
    monkeypatch.setattr(assignment, "GRID_OUTPUT_TSV", grid_path)

    assignment.write_direction_grid_table()

    grid = pd.read_csv(grid_path, sep="\t").set_index(["guild", "profile_direction"])
    keys = list(COUNTS)
    for guild, traits in (
        ("anaerobe_v2_depletion", (ANAEROBIC, OBLIGATE)),
        ("aerobe_v2", (AEROBIC, MICROAEROPHILIC)),
    ):
        means = dict(zip(keys, expected_guild_mean(*traits), strict=True))
        inc = [means[key] for key in keys if key.endswith("inc")]
        dec = [means[key] for key in keys if key.endswith("dec")]
        assert grid.loc[(guild, "decreased"), "mean_guild_z"] == pytest.approx(sum(dec) / len(dec))
        assert grid.loc[(guild, "increased"), "mean_guild_z"] == pytest.approx(sum(inc) / len(inc))
