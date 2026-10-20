"""The PD trait-taxa recipes refuse evidence built from other profiles, inputs or lollipop.

The evidence must hold recorded lollipop traits with their directions (#584, #694), the exact
current set of PD fecal profiles, in both directions (#644, #649), and the inputs and summary
it was built from, unchanged (#693). The evidence and carrier tables must be the ones the
lookup wrote (#703).
"""

import importlib.util
import json
import re
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
INPUT_NAMES = [
    "bugsigdb_taxa.tsv",
    "taxa_trait_matrix.tsv",
    "study_trait_profiles_wide.tsv",
    "study_trait_profile_study_metadata.tsv",
    "pd_six_trait_direction_summary.tsv",
]
INPUTS = {name: [f"{index:08x}"] * 8 for index, name in enumerate(INPUT_NAMES)}
TABLES = {lookup.CARRIERS_NAME: ["0000000c"] * 8, lookup.EVIDENCE_NAME: ["0000000e"] * 8}


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
    profiles: pd.DataFrame | None = None,
    values: dict[str, float] = TRAITS,
    inputs: dict[str, list[str]] = INPUTS,
    tables: dict[str, list[str]] = TABLES,
) -> dict[str, object]:
    profile_ids = lookup.pd_feces_profile_ids(metadata() if profiles is None else profiles)
    return {
        **lookup.profile_set_fingerprint(profile_ids),
        "traits": lookup.trait_directions(summary(values)).to_dict(orient="records"),
        "input_sha256_chunks": inputs,
        "output_sha256_chunks": tables,
    }


def check(**changes: object) -> None:
    """Run the guard on matching inputs, with any of them replaced."""
    inputs: dict[str, object] = {
        "evidence": evidence(),
        "carriers": carriers(),
        "provenance": provenance(),
        "trait_summary": summary(),
        "metadata": metadata(),
        "current_inputs": INPUTS,
        "current_tables": TABLES,
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


def test_changed_input_file_is_refused_and_named() -> None:
    """A data/derived rebuild that keeps the PD profile ids is still a change (#693)."""
    rebuilt = {**INPUTS, "taxa_trait_matrix.tsv": ["ffffffff"] * 8}

    with pytest.raises(ValueError, match=r"taxa_trait_matrix\.tsv changed since it was built"):
        check(current_inputs=rebuilt)


def test_missing_input_file_is_refused() -> None:
    missing = {**INPUTS, "bugsigdb_taxa.tsv": []}

    with pytest.raises(ValueError, match=r"bugsigdb_taxa\.tsv changed since it was built"):
        check(current_inputs=missing)


def test_input_the_evidence_never_recorded_is_refused() -> None:
    """A newly tracked input counts as a change, so the check runs in both directions."""
    grown = {**INPUTS, "new_input.tsv": ["0000000a"] * 8}

    with pytest.raises(ValueError, match=r"new_input\.tsv changed since it was built"):
        check(current_inputs=grown)


def test_provenance_without_input_hashes_is_refused() -> None:
    old_format = {key: value for key, value in provenance().items() if key != "input_sha256_chunks"}

    with pytest.raises(ValueError, match="records no input hashes"):
        check(provenance=old_format)


@pytest.mark.parametrize("table", [lookup.EVIDENCE_NAME, lookup.CARRIERS_NAME])
def test_table_other_than_the_one_the_lookup_wrote_is_refused(table: str) -> None:
    """Without the hash, a trait whose rows were lost looks like one without carriers (#703)."""
    edited = {**TABLES, table: ["ffffffff"] * 8}

    with pytest.raises(
        ValueError,
        match=rf"the tables changed after `just pd-trait-taxa` wrote them: {re.escape(table)}",
    ):
        check(current_tables=edited)


def test_provenance_without_table_hashes_is_refused() -> None:
    old_format = {
        key: value for key, value in provenance().items() if key != "output_sha256_chunks"
    }

    with pytest.raises(ValueError, match="records no hashes of the tables it wrote"):
        check(provenance=old_format)


def test_flipped_lollipop_direction_is_refused_and_named() -> None:
    flipped = {**TRAITS, "coccus": 0.08}

    with pytest.raises(
        ValueError, match=r"its trait directions differ .*coccus decreased, now increased"
    ):
        check(trait_summary=summary(flipped))


def test_evidence_direction_other_than_its_provenance_is_refused() -> None:
    flipped = {**TRAITS, "coccus": 0.08}

    with pytest.raises(
        ValueError,
        match=r"the evidence trait directions differ from its provenance: "
        r"coccus decreased, recorded increased",
    ):
        check(
            provenance=provenance(values=flipped),
            trait_summary=summary(flipped),
            carriers=carriers(flipped),
        )


def test_carrier_direction_other_than_its_provenance_is_refused() -> None:
    flipped = {**TRAITS, "coccus": 0.08}

    with pytest.raises(ValueError, match=r"the carrier table trait directions differ .*coccus"):
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


def test_trait_without_pd_carriers_is_accepted() -> None:
    """No PD-reported taxon carries sphere shaped, so no rows name it (#694)."""
    with_sphere = {**TRAITS, "sphere": -0.1}

    check(trait_summary=summary(with_sphere), provenance=provenance(values=with_sphere))


def test_evidence_trait_missing_from_its_provenance_is_refused() -> None:
    with_sphere = {**TRAITS, "sphere": -0.1}

    with pytest.raises(ValueError, match="the evidence has traits its provenance does not record"):
        check(evidence=evidence(values=with_sphere))


def test_provenance_missing_a_lollipop_trait_is_refused() -> None:
    with_sphere = {**TRAITS, "sphere": -0.1}

    with pytest.raises(
        ValueError, match=r"its recorded traits differ .*missing: sphere; extra: none"
    ):
        check(trait_summary=summary(with_sphere))


def test_evidence_without_directions_is_refused() -> None:
    old_format = evidence().drop(columns="trait_direction")

    with pytest.raises(ValueError, match="no trait_direction column"):
        check(evidence=old_format)


def write_saved_evidence(directory: Path, profile_ids: tuple[str, ...]) -> None:
    """Write the tables, then a provenance that records their hashes, as the lookup does."""
    evidence(profile_ids).to_csv(directory / lookup.EVIDENCE_NAME, sep="\t", index=False)
    carriers().to_csv(directory / lookup.CARRIERS_NAME, sep="\t", index=False)
    summary().to_csv(directory / lookup.SUMMARY_NAME, sep="\t", index=False)
    written = provenance(tables=lookup.table_fingerprints(directory))
    (directory / lookup.PROVENANCE_NAME).write_text(json.dumps(written), encoding="utf-8")


def test_evidence_without_provenance_is_refused(tmp_path: Path) -> None:
    write_saved_evidence(tmp_path, ("s1|e1|inc", "s1|e1|dec"))
    (tmp_path / lookup.PROVENANCE_NAME).unlink()

    with pytest.raises(ValueError, match=r"has no pd_six_trait_taxon_evidence_provenance\.json"):
        heatmap.require_current_trait_taxa_evidence(metadata(), tmp_path)


def fingerprint_inputs_as(monkeypatch: pytest.MonkeyPatch, guard: ModuleType) -> None:
    """Make the guard see INPUTS as the current input hashes."""
    globals_ = guard.require_current_trait_taxa_evidence.__globals__
    monkeypatch.setitem(globals_, "input_fingerprints", lambda _inputs: INPUTS)


def test_saved_matching_evidence_passes(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    write_saved_evidence(tmp_path, ("s1|e1|inc", "s1|e1|dec"))
    fingerprint_inputs_as(monkeypatch, heatmap)

    heatmap.require_current_trait_taxa_evidence(metadata(), tmp_path)


def test_saved_tables_missing_a_traits_rows_are_refused(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Rows lost for a whole trait show only in the tables' hashes (#703).

    With every coccus row dropped from both tables, the traits left still match the provenance.
    """
    write_saved_evidence(tmp_path, ("s1|e1|inc", "s1|e1|dec"))
    for name in (lookup.EVIDENCE_NAME, lookup.CARRIERS_NAME):
        table = pd.read_csv(tmp_path / name, sep="\t")
        kept = table[table["trait_id"].ne("coccus")]
        kept.to_csv(tmp_path / name, sep="\t", index=False)
    fingerprint_inputs_as(monkeypatch, heatmap)

    with pytest.raises(
        ValueError,
        match=r"changed after `just pd-trait-taxa` wrote them: "
        r"pd_six_trait_taxon_carriers\.tsv, pd_six_trait_taxon_profile_evidence\.tsv",
    ):
        heatmap.require_current_trait_taxa_evidence(metadata(), tmp_path)


def test_saved_evidence_is_checked_against_the_files_it_read(tmp_path: Path) -> None:
    """Without the stand-in hashes, the guard hashes the real inputs, which differ from INPUTS."""
    write_saved_evidence(tmp_path, ("s1|e1|inc", "s1|e1|dec"))

    with pytest.raises(
        ValueError, match=r"pd_six_trait_direction_summary\.tsv.* changed since it was built"
    ):
        heatmap.require_current_trait_taxa_evidence(metadata(), tmp_path)


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
    write_saved_evidence(tmp_path, ("s1|e1|inc", "s2|e2|inc"))
    # The contributors recipe reads the species heatmap's trait table before it writes.
    lookup.trait_directions(summary()).to_csv(
        tmp_path / "pd_species_six_traits_influence_heatmap_traits.tsv", sep="\t", index=False
    )
    before = snapshot(tmp_path)
    # A path bound to the committed output folder at import would escape the redirect.
    committed_dir = module.OUTPUT_DIR
    committed_before = file_stamps(committed_dir)
    monkeypatch.setattr(module, "OUTPUT_DIR", tmp_path)
    monkeypatch.setattr(module, "load_profile_metadata", metadata)
    fingerprint_inputs_as(monkeypatch, module)

    with pytest.raises(ValueError, match=r"^[^;]*not current PD fecal profiles"):
        module.main()

    assert snapshot(tmp_path) == before
    assert file_stamps(committed_dir) == committed_before


def test_contributors_stop_when_the_heatmap_trait_table_is_stale(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The evidence is current, but the trait table the note quotes is not (#694)."""
    write_saved_evidence(tmp_path, ("s1|e1|inc", "s1|e1|dec"))
    stale = lookup.trait_directions(summary({**TRAITS, "coccus": 0.08}))
    stale.to_csv(
        tmp_path / "pd_species_six_traits_influence_heatmap_traits.tsv", sep="\t", index=False
    )
    before = snapshot(tmp_path)
    monkeypatch.setattr(contributors, "OUTPUT_DIR", tmp_path)
    monkeypatch.setattr(contributors, "load_profile_metadata", metadata)
    fingerprint_inputs_as(monkeypatch, contributors)

    with pytest.raises(ValueError, match=r"heatmap_traits\.tsv does not match .*coccus"):
        contributors.main()

    assert snapshot(tmp_path) == before
