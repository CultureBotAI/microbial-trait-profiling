"""Profile metadata readers accept a blank or a named first header (#2417).

Metadata written before pandas 3.1 has a blank first header (#2396); metadata written by
`just trait-heatmaps` under pandas 3.1 names it ``profile_id``. Readers copy the index into
a ``profile_id`` column, and a frame whose index level and column share a name is refused
by pandas in merge, sort_values, groupby, reset_index and join. Every test here runs on
both header forms.
"""

import dataclasses
import importlib.util
import re
import sys
from pathlib import Path
from types import ModuleType

import pandas as pd
import pytest

from mtp.composition_denominators import DENOMINATOR_DEFINITION
from mtp.profile_metadata import with_index_column, with_profile_id_column

SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))


def load_script(name: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(f"index_2417_{name}", SCRIPTS_DIR / f"{name}.py")
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


oxygen_switches = load_script("investigate_taxa_oxygen_switches")
parkinsons = load_script("analyze_parkinsons_disease_trait_profiles")

HEADER_FORMS = pytest.mark.parametrize("index_name", [None, "profile_id"], ids=["blank", "named"])
TRAIT_A = "METPO:1000603"
TRAIT_B = "METPO:1000607"
PROFILES = {
    # profile key: (trait A count, trait B count)
    "1|e1|inc": (1, 3),
    "1|e1|dec": (4, 1),
    "2|e2|inc": (0, 2),
    "2|e2|dec": (3, 3),
}


def write_metadata(path: Path, frame: pd.DataFrame, index_name: str | None) -> Path:
    """Write metadata the way the heatmap producer does, with the given first header."""
    frame.rename_axis(index=index_name).to_csv(path, sep="\t")
    return path


def profile_metadata(condition: str = "Parkinson's disease") -> pd.DataFrame:
    keys = list(PROFILES)
    return pd.DataFrame(
        {
            "study_id": [int(key.split("|")[0]) for key in keys],
            "experiment_id": [key.split("|")[1] for key in keys],
            "inc_dec": [key.split("|")[2] for key in keys],
            "direction_label": ["increased" if k.endswith("inc") else "decreased" for k in keys],
            "condition": [condition] * len(keys),
            "body_site": ["feces"] * len(keys),
            "body_site_group": ["Feces"] * len(keys),
            "condition_group": ["Nervous system"] * len(keys),
            "host_species": ["Homo sapiens"] * len(keys),
            "title": [f"study {key}" for key in keys],
            "publication_year": [2020] * len(keys),
            "signature_count": [1] * len(keys),
            "taxon_count": [3] * len(keys),
            "total_trait_signal": [float(sum(values)) for values in PROFILES.values()],
            "denominator_definition": DENOMINATOR_DEFINITION,
        },
        index=keys,
    )


@HEADER_FORMS
def test_helper_leaves_no_ambiguous_profile_id(tmp_path: Path, index_name: str | None) -> None:
    path = write_metadata(tmp_path / "meta.tsv", profile_metadata(), index_name)
    read = pd.read_csv(path, sep="\t", index_col=0)
    assert read.index.name == index_name

    metadata = with_profile_id_column(read)

    assert metadata.index.name is None
    assert list(metadata["profile_id"]) == list(PROFILES)
    assert read.index.name == index_name, "the input frame is not modified"
    other = pd.DataFrame({"profile_id": list(PROFILES), "value": range(len(PROFILES))})
    assert len(other.merge(metadata, on="profile_id")) == len(PROFILES)
    assert list(metadata.sort_values(["study_id", "profile_id"])["profile_id"]) == sorted(PROFILES)
    assert metadata.groupby("profile_id").size().eq(1).all()
    assert "profile_id" in metadata.reset_index(drop=True).columns
    assert len(other.join(metadata[["study_id"]], on="profile_id")) == len(PROFILES)


def test_with_index_column_names_any_key() -> None:
    frame = pd.DataFrame({"label": ["a"]}, index=pd.Index(["METPO:1"], name="trait_id"))
    result = with_index_column(frame, "trait_id")
    assert result.index.name is None
    assert list(result["trait_id"]) == ["METPO:1"]
    assert list(result.sort_values("trait_id")["label"]) == ["a"]


@HEADER_FORMS
def test_taxa_oxygen_switches_merge_profiles_with_either_header(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, index_name: str | None
) -> None:
    metadata_path = write_metadata(tmp_path / "meta.tsv", profile_metadata(), index_name)
    taxa_path = tmp_path / "bugsigdb_taxa.tsv"
    pd.DataFrame(
        {
            "study_id": [1, 1, 2],
            "experiment_id": ["e1", "e1", "e2"],
            "signature_id": ["s1", "s2", "s3"],
            "direction": ["increased", "decreased", "decreased"],
            "taxon_id": ["NCBITaxon:1", "NCBITaxon:2", "NCBITaxon:1"],
            "taxon_name": ["Alpha", "Beta", "Alpha"],
            "taxonomic_rank": ["genus", "species", "genus"],
        }
    ).to_csv(taxa_path, sep="\t", index=False)
    annotations = pd.DataFrame(
        {
            "taxon_id": ["NCBITaxon:1", "NCBITaxon:2"],
            "taxon_name": ["Alpha", "Beta"],
            "taxonomic_rank": ["genus", "species"],
            "oxygen_class": ["anaerobic_only", "aerobic_only"],
            "aerobic_traits": ["", "aerobic"],
            "anaerobic_traits": ["anaerobic", ""],
            "aerobic_trait_signal": [0, 1],
            "anaerobic_trait_signal": [1, 0],
        }
    )
    monkeypatch.setattr(oxygen_switches, "PROFILE_METADATA_PATH", metadata_path)
    monkeypatch.setattr(oxygen_switches, "TAXA_PATH", taxa_path)

    metadata = oxygen_switches.load_profile_metadata()
    occurrences = oxygen_switches.load_taxa_occurrences(metadata, annotations)

    assert sorted(occurrences["profile_id"]) == ["1|e1|dec", "1|e1|inc", "2|e2|dec"]
    assert set(occurrences["oxygen_class"]) == {"anaerobic_only", "aerobic_only"}


@HEADER_FORMS
def test_parkinsons_feces_scores_sort_with_either_header(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, index_name: str | None
) -> None:
    keys = list(PROFILES)
    wide_path = tmp_path / "study_trait_profiles_wide.tsv"
    pd.DataFrame(
        {
            "study_id": [key.split("|")[0] for key in keys],
            "experiment_id": [key.split("|")[1] for key in keys],
            "inc_dec": [key.split("|")[2] for key in keys],
            TRAIT_A: [values[0] for values in PROFILES.values()],
            TRAIT_B: [values[1] for values in PROFILES.values()],
        }
    ).to_csv(wide_path, sep="\t", index=False)
    metadata_path = write_metadata(tmp_path / "meta.tsv", profile_metadata(), index_name)
    selection = tmp_path / "selection.tsv"
    pd.DataFrame(
        {"trait_id": [TRAIT_A, TRAIT_B], "trait_name": ["anaerobic", "obligately anaerobic"]}
    ).to_csv(selection, sep="\t", index=False)
    configs = tuple(
        dataclasses.replace(config, trait_selection_path=selection)
        for config in parkinsons.ANALYSIS_CONFIGS
    )
    monkeypatch.setattr(parkinsons, "METADATA_PATH", metadata_path)
    monkeypatch.setattr(parkinsons, "PROFILE_MATRIX_PATH", wide_path)
    monkeypatch.setattr(parkinsons, "ANALYSIS_CONFIGS", configs)

    metadata, counts = parkinsons.load_profiles()
    scores = parkinsons.parkinsons_feces_scores(parkinsons.add_trait_scores(counts, metadata))

    assert list(scores["profile_id"]) == ["1|e1|dec", "1|e1|inc", "2|e2|dec", "2|e2|inc"]
    assert scores.index.name is None


def test_helper_refuses_multiindex_and_conflicting_column() -> None:
    """Both edge cases fail loudly instead of silently (#2422)."""
    multi = pd.DataFrame(
        {"value": [1]}, index=pd.MultiIndex.from_tuples([("a", "b")], names=["x", "y"])
    )
    with pytest.raises(ValueError, match="MultiIndex"):
        with_profile_id_column(multi)
    conflicting = pd.DataFrame({"profile_id": ["other"]}, index=["1|e1|inc"])
    with pytest.raises(ValueError, match="disagrees"):
        with_profile_id_column(conflicting)
    consistent = pd.DataFrame({"profile_id": ["1|e1|inc"]}, index=["1|e1|inc"])
    assert list(with_profile_id_column(consistent)["profile_id"]) == ["1|e1|inc"]


INDEX_COPY = re.compile(r"""\[\s*["'](?:profile_id|trait_id|pair_id)["']\s*\]\s*=\s*\w+\.index\b""")


def test_no_reader_copies_a_metadata_index_into_a_key_column_by_hand() -> None:
    """Index-to-column copies of metadata keys go through mtp.profile_metadata (#2421)."""
    repo = SCRIPTS_DIR.parent
    offenders = [
        f"{path.relative_to(repo)}:{number}"
        for root in (repo / "scripts", repo / "src")
        for path in sorted(root.rglob("*.py"))
        if "adhoc" not in path.parts
        for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1)
        if INDEX_COPY.search(line)
    ]
    assert offenders == [], offenders
