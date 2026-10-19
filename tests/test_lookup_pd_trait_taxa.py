"""The PD trait-taxa lookup traces all six traits and carries their direction as data (#584)."""

import hashlib
import importlib.util
import math
import sys
from pathlib import Path

import pandas as pd
import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "lookup_pd_trait_taxa.py"
sys.path.insert(0, str(SCRIPT.parent))
SPEC = importlib.util.spec_from_file_location("lookup_pd_trait_taxa_under_test", SCRIPT)
assert SPEC is not None
assert SPEC.loader is not None
lookup = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(lookup)

# Three traits higher in increased profiles and three higher in decreased ones, as after #555.
LOLLIPOP = {
    "halo": 0.494,
    "arsenate": 0.420,
    "aerobic": 0.349,
    "coccus": -0.081,
    "sphere": -0.133,
    "cellulose": -0.701,
}


def summary(values: dict[str, float] = LOLLIPOP) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "trait_id": list(values),
            "trait": [f"{trait_id} trait" for trait_id in values],
            "increased_minus_decreased_mean_z": list(values.values()),
        }
    )


@pytest.mark.parametrize(
    ("difference", "direction"),
    [(0.081, "increased"), (-0.081, "decreased"), (0.0, "neither")],
)
def test_direction_is_the_sign_of_the_mean_z_difference(difference: float, direction: str) -> None:
    assert lookup.trait_direction(difference) == direction


def test_missing_mean_z_difference_has_no_direction() -> None:
    with pytest.raises(ValueError, match="no increased-minus-decreased mean z"):
        lookup.trait_direction(math.nan)


def test_all_six_traits_are_kept_whatever_their_sign() -> None:
    directions = lookup.trait_directions(summary())

    assert directions["trait_id"].tolist() == list(LOLLIPOP)
    assert directions["trait_direction"].tolist() == ["increased"] * 3 + ["decreased"] * 3


def test_profile_fingerprint_ignores_order_and_sees_membership() -> None:
    first = lookup.profile_set_fingerprint(["b", "a"])

    assert first == lookup.profile_set_fingerprint(["a", "b"])
    assert first["pd_feces_profiles"] == 2
    assert first != lookup.profile_set_fingerprint(["a", "c"])
    assert first != lookup.profile_set_fingerprint(["a", "b", "c"])


def test_pd_feces_profiles_exclude_other_conditions_and_sites() -> None:
    metadata = pd.DataFrame(
        {
            "profile_id": ["p2", "p1", "oral", "crc"],
            "condition": ["Parkinson's disease"] * 3 + ["Colorectal cancer"],
            "body_site_group": ["Feces", "Feces", "Oral cavity", "Feces"],
        }
    )

    assert lookup.pd_feces_profile_ids(metadata) == ["p1", "p2"]


def observed_occurrences() -> pd.DataFrame:
    """Taxon A in two increased profiles; taxon B in three decreased and one increased."""
    rows = [
        ("A", "s1", "s1|e1|inc", "increased", 2, 0),
        ("A", "s2", "s2|e1|inc", "increased", 2, 0),
        ("B", "s1", "s1|e1|dec", "decreased", 0, 3),
        ("B", "s2", "s2|e1|dec", "decreased", 0, 3),
        ("B", "s3", "s3|e1|dec", "decreased", 0, 3),
        ("B", "s3", "s3|e1|inc", "increased", 0, 3),
    ]
    frame = pd.DataFrame(
        rows,
        columns=["taxon_id", "study_id", "profile_id", "direction", "halo", "cellulose"],
    )
    frame["taxon_name"] = "Taxon " + frame["taxon_id"]
    frame["taxonomic_rank"] = "genus"
    frame["experiment_id"] = frame["profile_id"].str.split("|").str[1]
    return frame


def zscores() -> pd.DataFrame:
    profiles = ["s1|e1|inc", "s2|e1|inc", "s1|e1|dec", "s2|e1|dec", "s3|e1|dec", "s3|e1|inc"]
    return pd.DataFrame({"halo": [1.0] * 6, "cellulose": [1.0] * 6}, index=profiles)


def test_tables_carry_each_traits_direction() -> None:
    directions = lookup.trait_directions(summary({"halo": 0.494, "cellulose": -0.701}))

    carriers, evidence = lookup.build_carriers_and_evidence(
        observed_occurrences(), zscores(), directions
    )

    by_trait = evidence.groupby("trait_id")["trait_direction"].unique()
    assert by_trait["halo"].tolist() == ["increased"]
    assert by_trait["cellulose"].tolist() == ["decreased"]
    assert carriers.set_index("trait_id")["trait_direction"].to_dict() == {
        "halo": "increased",
        "cellulose": "decreased",
    }
    cellulose = carriers[carriers["trait_id"].eq("cellulose")].iloc[0]
    assert cellulose["trait_increased_minus_decreased_mean_z"] == -0.701
    assert (cellulose["increased_profiles"], cellulose["decreased_profiles"]) == (1, 3)


def test_carriers_are_ranked_on_their_traits_observed_side() -> None:
    """A decreased-direction trait lists the taxon with more decreased profiles first."""
    occurrences = observed_occurrences()
    occurrences["halo"] = 1
    occurrences["cellulose"] = 1
    directions = lookup.trait_directions(summary({"halo": 0.494, "cellulose": -0.701}))

    carriers, _ = lookup.build_carriers_and_evidence(occurrences, zscores(), directions)

    order = carriers.groupby("trait_id", sort=False)["taxon_id"].apply(list).to_dict()
    assert list(order) == ["halo", "cellulose"]
    assert order["halo"] == ["A", "B"]
    assert order["cellulose"] == ["B", "A"]


def test_note_names_every_trait_with_its_direction() -> None:
    directions = lookup.trait_directions(summary())
    carriers = pd.DataFrame(
        {
            "trait_id": list(LOLLIPOP),
            "taxon_name": "Escherichia",
            "taxonomic_rank": "genus",
            "increased_profiles": 9,
            "decreased_profiles": 2,
            "increased_positive_z_profiles": 9,
            "decreased_positive_z_profiles": 2,
        }
    )
    profiles = pd.DataFrame({"inc_dec": ["inc"] * 128 + ["dec"] * 119})

    note = lookup.lookup_note(directions, carriers, profiles)

    assert "128 increased and 119 decreased PD fecal profiles" in note
    assert "coccus trait decreased (-0.081)" in note
    assert "halo trait increased (+0.494)" in note
    assert note.count("Observed direction: decreased") == 3
    assert note.count("in positive-z decreased profiles, the trait's observed direction") == 3
    assert "four" not in note


def test_stale_summary_names_itself_and_how_to_rebuild_it() -> None:
    """A rebuilt metadata table moves the profile counts; the saved lollipop is then stale."""
    saved = summary().assign(increased_profiles=128)
    recomputed = saved.assign(increased_profiles=127)

    lookup.check_saved_summary(saved.copy(), saved)
    with pytest.raises(
        ValueError, match=r"pd_six_trait_direction_summary\.tsv no longer matches"
    ) as error:
        lookup.check_saved_summary(recomputed, saved)
    assert "scripts/create_pd_cellulose_analysis.py" in str(error.value)
    assert "increased_profiles" in str(error.value)


def test_tracked_inputs_are_the_files_the_evidence_reads(tmp_path: Path) -> None:
    inputs = lookup.tracked_inputs(tmp_path)

    assert sorted(inputs) == [
        "bugsigdb_taxa.tsv",
        "pd_six_trait_direction_summary.tsv",
        "study_trait_profile_study_metadata.tsv",
        "study_trait_profiles_wide.tsv",
        "taxa_trait_matrix.tsv",
    ]
    assert inputs["pd_six_trait_direction_summary.tsv"] == tmp_path / lookup.SUMMARY_NAME


def test_input_fingerprints_are_sha256_chunks_and_empty_when_missing(tmp_path: Path) -> None:
    present = tmp_path / "present.tsv"
    present.write_bytes(b"trait_id\tvalue\n")

    fingerprints = lookup.input_fingerprints(
        {"present.tsv": present, "missing.tsv": tmp_path / "missing.tsv"}
    )

    digest = hashlib.sha256(b"trait_id\tvalue\n").hexdigest()
    assert "".join(fingerprints["present.tsv"]) == digest
    assert [len(chunk) for chunk in fingerprints["present.tsv"]] == [8] * 8
    assert fingerprints["missing.tsv"] == []


def test_provenance_records_the_profiles_traits_and_input_hashes(tmp_path: Path) -> None:
    (tmp_path / lookup.SUMMARY_NAME).write_bytes(b"trait_id\n")
    directions = lookup.trait_directions(summary())

    provenance = lookup.evidence_provenance(["p2", "p1"], directions, tmp_path)

    assert provenance["pd_feces_profiles"] == 2
    assert [trait["trait_id"] for trait in provenance["traits"]] == list(LOLLIPOP)
    hashes = provenance["input_sha256_chunks"]
    assert sorted(hashes) == sorted(lookup.tracked_inputs(tmp_path))
    summary_digest = hashlib.sha256(b"trait_id\n").hexdigest()
    assert "".join(hashes[lookup.SUMMARY_NAME]) == summary_digest
