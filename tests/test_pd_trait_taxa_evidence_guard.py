"""The PD trait-taxa recipes refuse evidence built from other profiles or another lollipop.

The evidence must cover the six lollipop traits with their current directions (#584) and the
exact current set of PD fecal profiles, in both directions (#644, #649).
"""

import importlib.util
import json
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


lookup = load_script("lookup_pd_trait_taxa")
heatmap = load_script("create_pd_six_trait_taxon_heatmap")
exploratory = load_script("create_pd_exploratory_taxa_influence_heatmap")
contributors = load_script("create_pd_trait_contributors")

TRAITS = {"halo": 0.5, "aero": 0.3, "coccus": -0.08, "cellulose": -0.7}


def metadata(pd_feces: tuple[str, ...] = ("s1|e1|inc", "s1|e1|dec")) -> pd.DataFrame:
    """The given PD fecal profiles, plus one PD oral and one other fecal profile."""
    ids = [*pd_feces, "s2|e2|inc", "s3|e3|dec"]
    frame = pd.DataFrame(
        {
            "profile_id": ids,
            "study_id": [profile.split("|")[0] for profile in ids],
            "experiment_id": [profile.split("|")[1] for profile in ids],
            "inc_dec": [profile.split("|")[2] for profile in ids],
            "condition": [heatmap.PD_CONDITION] * (len(pd_feces) + 1) + ["Colorectal cancer"],
            "body_site_group": [heatmap.FECES] * len(pd_feces) + ["Oral cavity", heatmap.FECES],
        }
    )
    return frame.set_index("profile_id", drop=False)


def summary(values: dict[str, float] = TRAITS) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "trait_id": list(values),
            "trait": [f"{trait_id} trait" for trait_id in values],
            "increased_minus_decreased_mean_z": list(values.values()),
        }
    )


def directions(values: dict[str, float] = TRAITS) -> dict[str, str]:
    return {trait_id: lookup.trait_direction(value) for trait_id, value in values.items()}


def evidence(
    profile_ids: tuple[str, ...] = ("s1|e1|inc", "s1|e1|dec"),
    values: dict[str, float] = TRAITS,
) -> pd.DataFrame:
    rows = [
        (trait_id, direction, profile)
        for trait_id, direction in directions(values).items()
        for profile in profile_ids
    ]
    return pd.DataFrame(rows, columns=["trait_id", "trait_direction", "profile_id"])


def carriers(values: dict[str, float] = TRAITS) -> pd.DataFrame:
    by_trait = directions(values)
    return pd.DataFrame({"trait_id": list(by_trait), "trait_direction": list(by_trait.values())})


def provenance(
    profiles: pd.DataFrame | None = None, values: dict[str, float] = TRAITS
) -> dict[str, object]:
    profile_ids = lookup.pd_feces_profile_ids(metadata() if profiles is None else profiles)
    return {
        **lookup.profile_set_fingerprint(profile_ids),
        "traits": lookup.trait_directions(summary(values)).to_dict(orient="records"),
    }


def check(**changes: object) -> None:
    """Run the guard on matching inputs, with any of them replaced."""
    inputs: dict[str, object] = {
        "evidence": evidence(),
        "carriers": carriers(),
        "provenance": provenance(),
        "trait_summary": summary(),
        "metadata": metadata(),
    }
    inputs.update(changes)
    heatmap.check_trait_taxa_evidence(**inputs)


def test_matching_evidence_passes() -> None:
    check()


def test_three_positive_traits_are_accepted() -> None:
    """The four-positive assumption is gone: any mix of directions is data (#584)."""
    values = {"halo": 0.5, "aero": 0.3, "arsenate": 0.2, "coccus": -0.08, "cellulose": -0.7}

    check(
        evidence=evidence(values=values),
        carriers=carriers(values),
        provenance=provenance(values=values),
        trait_summary=summary(values),
    )


def test_evidence_profile_outside_current_pd_feces_is_refused() -> None:
    with pytest.raises(ValueError, match=r"1 evidence profiles .*\(first: s2\|e2\|inc\)"):
        check(evidence=evidence(("s1|e1|inc", "s2|e2|inc")))


def test_pd_feces_profile_added_since_the_evidence_is_refused() -> None:
    """The check works in both directions: a new PD fecal profile is a change too (#649)."""
    grown = metadata(("s1|e1|inc", "s1|e1|dec", "s4|e4|inc"))

    with pytest.raises(
        ValueError, match=r"built from 2 PD fecal profiles .*current metadata has 3"
    ):
        check(metadata=grown)


def test_pd_feces_profile_swapped_for_another_is_refused() -> None:
    swapped = metadata(("s1|e1|inc", "s4|e4|dec"))

    with pytest.raises(
        ValueError, match=r"built from 2 PD fecal profiles .*current metadata has 2"
    ):
        check(metadata=swapped, evidence=evidence(("s1|e1|inc",)))


def test_flipped_direction_is_refused_and_named() -> None:
    flipped = {**TRAITS, "coccus": 0.08}

    with pytest.raises(
        ValueError, match=r"evidence trait directions differ .*coccus decreased, now increased"
    ):
        check(trait_summary=summary(flipped), provenance=provenance(values=flipped))


def test_flipped_direction_in_the_carrier_table_is_refused() -> None:
    flipped = {**TRAITS, "coccus": 0.08}

    with pytest.raises(ValueError, match=r"carrier table trait directions differ .*coccus"):
        check(
            evidence=evidence(values=flipped),
            trait_summary=summary(flipped),
            provenance=provenance(values=flipped),
        )


def test_changed_lollipop_value_is_refused() -> None:
    """Same directions, but the z-scores the evidence used have moved."""
    moved = {**TRAITS, "coccus": -0.09}

    with pytest.raises(ValueError, match=r"mean-z differences differ .* for coccus"):
        check(trait_summary=summary(moved))


def test_missing_and_extra_traits_are_listed() -> None:
    with_sphere = {**TRAITS, "sphere": -0.1}

    with pytest.raises(ValueError, match=r"evidence traits differ .*missing: sphere; extra: none"):
        check(trait_summary=summary(with_sphere), provenance=provenance(values=with_sphere))


def test_evidence_without_directions_is_refused() -> None:
    old_format = evidence().drop(columns="trait_direction")

    with pytest.raises(ValueError, match="no trait_direction column"):
        check(evidence=old_format)


def test_evidence_without_provenance_is_refused(tmp_path: Path) -> None:
    evidence().to_csv(tmp_path / lookup.EVIDENCE_NAME, sep="\t", index=False)
    carriers().to_csv(tmp_path / lookup.CARRIERS_NAME, sep="\t", index=False)
    summary().to_csv(tmp_path / lookup.SUMMARY_NAME, sep="\t", index=False)

    with pytest.raises(ValueError, match=r"has no pd_six_trait_taxon_evidence_provenance\.json"):
        heatmap.require_current_trait_taxa_evidence(metadata(), tmp_path)


def test_saved_matching_evidence_passes(tmp_path: Path) -> None:
    evidence().to_csv(tmp_path / lookup.EVIDENCE_NAME, sep="\t", index=False)
    carriers().to_csv(tmp_path / lookup.CARRIERS_NAME, sep="\t", index=False)
    summary().to_csv(tmp_path / lookup.SUMMARY_NAME, sep="\t", index=False)
    (tmp_path / lookup.PROVENANCE_NAME).write_text(json.dumps(provenance()), encoding="utf-8")

    heatmap.require_current_trait_taxa_evidence(metadata(), tmp_path)


def write_contributor_inputs(directory: Path) -> None:
    """Write what the contributors recipe reads before its first write.

    With these in place, a guard moved below that write lets the write happen, and the
    snapshot comparison sees it.
    """
    for name in (
        "pd_trait_contributors_cells.tsv",
        "pd_trait_contributors_rank_coverage.tsv",
        "pd_species_six_traits_influence_heatmap_traits.tsv",
    ):
        pd.DataFrame({"unused": [0]}).to_csv(directory / name, sep="\t", index=False)
    pd.DataFrame({"taxon": ["NCBITaxon:853"]}).to_csv(
        directory / "pd_trait_contributors_selected_taxa.tsv", sep="\t", index=False
    )
    pd.DataFrame(
        {
            "study_id": ["s1", "s3"],
            "experiment_id": ["e1", "e3"],
            "direction": ["decreased", "decreased"],
            "taxon_id": ["NCBITaxon:853", "NCBITaxon:853"],
            "taxon_name": ["Faecalibacterium prausnitzii"] * 2,
            "taxonomic_rank": ["species"] * 2,
        }
    ).to_csv(directory / "bugsigdb_taxa.tsv", sep="\t", index=False)


def snapshot(directory: Path) -> dict[str, bytes]:
    return {path.name: path.read_bytes() for path in directory.iterdir()}


def file_stamps(directory: Path) -> dict[str, tuple[int, int]]:
    """Size and modification time of each file, enough to notice a write."""
    if not directory.exists():
        return {}
    return {
        path.name: (path.stat().st_size, path.stat().st_mtime_ns) for path in directory.iterdir()
    }


@pytest.mark.parametrize("module", [heatmap, exploratory, contributors], ids=lambda m: m.STEM)
def test_recipe_stops_before_writing_when_evidence_is_stale(
    module: ModuleType, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    evidence(("s1|e1|inc", "s2|e2|inc")).to_csv(
        tmp_path / lookup.EVIDENCE_NAME, sep="\t", index=False
    )
    carriers().to_csv(tmp_path / lookup.CARRIERS_NAME, sep="\t", index=False)
    summary().to_csv(tmp_path / lookup.SUMMARY_NAME, sep="\t", index=False)
    (tmp_path / lookup.PROVENANCE_NAME).write_text(json.dumps(provenance()), encoding="utf-8")
    write_contributor_inputs(tmp_path)
    before = snapshot(tmp_path)
    # A path bound to the committed output folder at import would escape the redirect.
    committed_dir = module.OUTPUT_DIR
    committed_before = file_stamps(committed_dir)
    monkeypatch.setattr(module, "OUTPUT_DIR", tmp_path)
    monkeypatch.setattr(module, "load_profile_metadata", metadata)
    # load_taxon_specificity reads TAXA_PATH from the module that defines it.
    taxa_module = sys.modules[contributors.load_taxon_specificity.__module__]
    monkeypatch.setattr(taxa_module, "TAXA_PATH", tmp_path / "bugsigdb_taxa.tsv")

    with pytest.raises(ValueError, match="not current PD fecal profiles"):
        module.main()

    assert snapshot(tmp_path) == before
    assert file_stamps(committed_dir) == committed_before
