import pandas as pd
import pytest

from mtp.human_host_corpus import (
    UNRECORDED_HOST,
    excluded_host_note,
    excluded_host_summary,
    join_human_host_profiles,
)


def occurrences() -> pd.DataFrame:
    """Target-taxon occurrences of one human, one mouse and one unrecorded-host study."""
    return pd.DataFrame(
        [
            ("s1", "s1-1", "increased", "NCBITaxon:562", "s1|s1-1|inc"),
            ("s1", "s1-1", "decreased", "NCBITaxon:851", "s1|s1-1|dec"),
            ("s2", "s2-1", "increased", "NCBITaxon:562", "s2|s2-1|inc"),
            ("s2", "s2-1", "increased", "NCBITaxon:851", "s2|s2-1|inc"),
            ("s3", "s3-1", "decreased", "NCBITaxon:562", "s3|s3-1|dec"),
        ],
        columns=["study_id", "experiment_id", "direction", "taxon_id", "profile_id"],
    )


def profiles(profile_ids: list[str]) -> pd.DataFrame:
    """Human-host profile metadata indexed by profile ID."""
    return pd.DataFrame(
        {
            "study_id": [profile_id.split("|")[0] for profile_id in profile_ids],
            "experiment_id": [profile_id.split("|")[1] for profile_id in profile_ids],
            "condition": "Colorectal cancer",
        },
        index=pd.Index(profile_ids, name="profile_id"),
    )


def experiment_hosts() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "study_id": ["s1", "s2", "s3"],
            "experiment_id": ["s1-1", "s2-1", "s3-1"],
            "host_species": ["Homo sapiens", "Mus musculus", None],
        }
    )


def test_occurrences_outside_the_corpus_are_set_aside_with_their_host() -> None:
    joined, excluded = join_human_host_profiles(
        occurrences(), profiles(["s1|s1-1|inc", "s1|s1-1|dec"]), experiment_hosts()
    )
    assert joined["profile_id"].tolist() == ["s1|s1-1|inc", "s1|s1-1|dec"]
    assert joined["condition"].notna().all()
    assert {"study_id_taxa", "study_id"} <= set(joined.columns)
    assert excluded["profile_id"].tolist() == ["s2|s2-1|inc", "s2|s2-1|inc", "s3|s3-1|dec"]
    assert excluded["host_species"].tolist() == ["Mus musculus", "Mus musculus", UNRECORDED_HOST]


@pytest.mark.parametrize("unrecorded", [None, "", "unknown", "NA"])
def test_an_unrecorded_host_is_set_aside_as_not_recorded(unrecorded: str | None) -> None:
    """The corpus parser writes "unknown" for an empty host field (#713)."""
    hosts = experiment_hosts()
    hosts.loc[2, "host_species"] = unrecorded
    _, excluded = join_human_host_profiles(
        occurrences(), profiles(["s1|s1-1|inc", "s1|s1-1|dec"]), hosts
    )
    assert excluded.loc[excluded["study_id"] == "s3", "host_species"].tolist() == [UNRECORDED_HOST]


def test_duplicate_profile_metadata_raises() -> None:
    duplicated = profiles(["s1|s1-1|inc", "s1|s1-1|dec", "s1|s1-1|dec"])
    with pytest.raises(pd.errors.MergeError):
        join_human_host_profiles(occurrences(), duplicated, experiment_hosts())


def test_a_human_host_profile_missing_from_the_corpus_raises() -> None:
    with pytest.raises(ValueError, match="human-host profile"):
        join_human_host_profiles(occurrences(), profiles(["s1|s1-1|inc"]), experiment_hosts())


def test_an_experiment_without_a_host_record_raises() -> None:
    hosts = experiment_hosts().iloc[:2]
    with pytest.raises(ValueError, match="no experiment host record"):
        join_human_host_profiles(occurrences(), profiles(["s1|s1-1|inc", "s1|s1-1|dec"]), hosts)


def test_excluded_host_note_counts_one_direction() -> None:
    _, excluded = join_human_host_profiles(
        occurrences(), profiles(["s1|s1-1|inc", "s1|s1-1|dec"]), experiment_hosts()
    )
    assert excluded_host_note(excluded, "increased") == (
        "Human-host profiles only: 2 increased-profile occurrences (1 profile) "
        "from non-human hosts are excluded."
    )
    assert excluded_host_note(excluded, "decreased") == (
        "Human-host profiles only: 1 decreased-profile occurrence (1 profile) "
        "from non-human or unrecorded hosts is excluded."
    )
    assert excluded_host_note(excluded.iloc[0:0], "increased") == "Human-host profiles only."


def test_excluded_host_summary_agrees_in_number() -> None:
    _, excluded = join_human_host_profiles(
        occurrences().iloc[[2]], profiles(["s1|s1-1|inc"]), experiment_hosts()
    )
    assert excluded_host_summary(excluded) == (
        "1 occurrence (1 increased) in 1 profile from 1 study comes from "
        "experiments whose host is not recorded as Homo sapiens: Mus musculus 1."
    )


def test_excluded_host_summary_counts_directions_and_hosts() -> None:
    _, excluded = join_human_host_profiles(
        occurrences(), profiles(["s1|s1-1|inc", "s1|s1-1|dec"]), experiment_hosts()
    )
    assert excluded_host_summary(excluded) == (
        "3 occurrences (1 decreased, 2 increased) in 2 profiles from 2 studies come from "
        "experiments whose host is not recorded as Homo sapiens: Mus musculus 2, not recorded 1."
    )
