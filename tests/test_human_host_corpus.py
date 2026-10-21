import pandas as pd
import pytest

from mtp.human_host_corpus import (
    EXCLUDED_BY_HOST,
    EXCLUDED_BY_SEQUENCING_TYPE,
    UNRECORDED_HOST,
    UNRECORDED_SEQUENCING_TYPE,
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
            "sequencing_type": ["16S", "16S", "WMS"],
        }
    )


def other_method_occurrences() -> pd.DataFrame:
    """The occurrences above plus two of a human ITS experiment and one of a human
    experiment with no recorded sequencing type, both outside the corpus (#742)."""
    extra = pd.DataFrame(
        [
            ("s4", "s4-1", "increased", "NCBITaxon:562", "s4|s4-1|inc"),
            ("s4", "s4-1", "decreased", "NCBITaxon:851", "s4|s4-1|dec"),
            ("s5", "s5-1", "decreased", "NCBITaxon:562", "s5|s5-1|dec"),
        ],
        columns=occurrences().columns,
    )
    return pd.concat([occurrences(), extra], ignore_index=True)


def other_method_hosts() -> pd.DataFrame:
    extra = pd.DataFrame(
        {
            "study_id": ["s4", "s5"],
            "experiment_id": ["s4-1", "s5-1"],
            "host_species": ["Homo sapiens", "Homo sapiens"],
            "sequencing_type": ["ITS / ITS2", "unknown"],
        }
    )
    return pd.concat([experiment_hosts(), extra], ignore_index=True)


def test_occurrences_outside_the_corpus_are_set_aside_with_their_host() -> None:
    joined, excluded = join_human_host_profiles(
        occurrences(), profiles(["s1|s1-1|inc", "s1|s1-1|dec"]), experiment_hosts()
    )
    assert joined["profile_id"].tolist() == ["s1|s1-1|inc", "s1|s1-1|dec"]
    assert joined["condition"].notna().all()
    assert {"study_id_taxa", "study_id"} <= set(joined.columns)
    assert excluded["profile_id"].tolist() == ["s2|s2-1|inc", "s2|s2-1|inc", "s3|s3-1|dec"]
    assert excluded["host_species"].tolist() == ["Mus musculus", "Mus musculus", UNRECORDED_HOST]
    assert excluded["exclusion_reason"].eq(EXCLUDED_BY_HOST).all()


def test_human_occurrences_of_other_sequencing_methods_are_set_aside() -> None:
    """A Homo sapiens experiment outside 16S and WMS is an exclusion, not a key
    mismatch (#742)."""
    _, excluded = join_human_host_profiles(
        other_method_occurrences(), profiles(["s1|s1-1|inc", "s1|s1-1|dec"]), other_method_hosts()
    )
    by_study = excluded.drop_duplicates("study_id").set_index("study_id")
    assert by_study["exclusion_reason"].to_dict() == {
        "s2": EXCLUDED_BY_HOST,
        "s3": EXCLUDED_BY_HOST,
        "s4": EXCLUDED_BY_SEQUENCING_TYPE,
        "s5": EXCLUDED_BY_SEQUENCING_TYPE,
    }
    assert by_study.loc[["s4", "s5"], "sequencing_type"].tolist() == [
        "ITS / ITS2",
        UNRECORDED_SEQUENCING_TYPE,
    ]


def test_a_human_16s_or_wms_profile_missing_from_the_corpus_still_raises() -> None:
    hosts = other_method_hosts()
    hosts.loc[hosts["study_id"] == "s4", "sequencing_type"] = "WMS"
    with pytest.raises(ValueError, match="human-host profile"):
        join_human_host_profiles(
            other_method_occurrences(), profiles(["s1|s1-1|inc", "s1|s1-1|dec"]), hosts
        )


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
        "Human-host 16S and shotgun profiles only: 2 increased-profile occurrences (1 profile) "
        "from non-human hosts are excluded."
    )
    assert excluded_host_note(excluded, "decreased") == (
        "Human-host 16S and shotgun profiles only: 1 decreased-profile occurrence (1 profile) "
        "from non-human or unrecorded hosts is excluded."
    )
    assert excluded_host_note(excluded.iloc[0:0], "increased") == (
        "Human-host 16S and shotgun profiles only."
    )


def test_excluded_host_note_names_both_reasons() -> None:
    _, excluded = join_human_host_profiles(
        other_method_occurrences(), profiles(["s1|s1-1|inc", "s1|s1-1|dec"]), other_method_hosts()
    )
    assert excluded_host_note(excluded, "increased") == (
        "Human-host 16S and shotgun profiles only: 3 increased-profile occurrences (2 profiles) "
        "from non-human hosts and other sequencing methods are excluded."
    )
    assert excluded_host_note(excluded, "decreased") == (
        "Human-host 16S and shotgun profiles only: 3 decreased-profile occurrences (3 profiles) "
        "from non-human or unrecorded hosts and other or unrecorded sequencing methods are "
        "excluded."
    )


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


def test_excluded_host_summary_counts_each_reason() -> None:
    _, excluded = join_human_host_profiles(
        other_method_occurrences(), profiles(["s1|s1-1|inc", "s1|s1-1|dec"]), other_method_hosts()
    )
    assert excluded_host_summary(excluded) == (
        "6 occurrences (3 decreased, 3 increased) in 5 profiles from 4 studies lie outside "
        "the corpus: 3 come from experiments whose host is not recorded as Homo sapiens "
        "(Mus musculus 2, not recorded 1), and 3 from Homo sapiens experiments sequenced by "
        "neither 16S nor WMS (ITS / ITS2 2, not recorded 1)."
    )
