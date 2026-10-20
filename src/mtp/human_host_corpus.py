"""Keep BugSigDB taxon occurrences inside the human-host profile corpus (#587).

`data/derived/bugsigdb_taxa.tsv` lists the reported taxa of every BugSigDB host.
The profile metadata of the trait-profile investigations,
`trait-visualization/study_trait_profile_study_metadata.tsv`, keeps the profiles
whose experiment host is recorded as Homo sapiens in `data/raw/full_dump.csv`
(`scripts/create_study_trait_profile_heatmap.py`). An occurrence outside that
corpus has no condition or body site, so it is set aside and counted instead of
being drawn as an `unknown` condition.
"""

from __future__ import annotations

import pandas as pd

__all__ = [
    "HUMAN_HOST_SPECIES",
    "UNRECORDED_HOST",
    "excluded_host_note",
    "excluded_host_summary",
    "join_human_host_profiles",
]

HUMAN_HOST_SPECIES = "Homo sapiens"
UNRECORDED_HOST = "not recorded"
# Values that mean the source has no host for an experiment; "unknown" is what the
# corpus's full-dump parser writes for an empty field.
UNRECORDED_HOST_VALUES = {"", "unknown", "na", "nan"}
EXPERIMENT_KEY = ["study_id", "experiment_id"]


def join_human_host_profiles(
    occurrences: pd.DataFrame,
    profiles: pd.DataFrame,
    experiment_hosts: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Join profile metadata to occurrences in the corpus; return the rest with their host.

    `occurrences` has one row per taxon and profile, with `study_id`, `experiment_id`
    and `profile_id`; `profiles` is the corpus metadata indexed by profile ID;
    `experiment_hosts` gives `host_species` per study and experiment, read from the
    file the corpus filter read (#713). Columns that both tables carry keep the
    occurrence values with a `_taxa` suffix.

    An occurrence outside the corpus must come from an experiment whose host is not
    recorded as Homo sapiens: another species, or no host at all, which is reported
    as "not recorded". An occurrence of a Homo sapiens experiment, or of an
    experiment with no host row, raises, so a key mismatch cannot pass as a host
    exclusion.
    """
    in_corpus = occurrences["profile_id"].isin(profiles.index)
    excluded = occurrences.loc[~in_corpus].merge(
        experiment_hosts[[*EXPERIMENT_KEY, "host_species"]],
        on=EXPERIMENT_KEY,
        how="left",
        validate="many_to_one",
        indicator=True,
    )
    unmatched = excluded.loc[excluded["_merge"] == "left_only", "profile_id"]
    if not unmatched.empty:
        msg = (
            f"{unmatched.nunique():,} profile(s) outside the human-host corpus have no "
            f"experiment host record, e.g. {sorted(unmatched.unique())[:3]}"
        )
        raise ValueError(msg)
    excluded = excluded.drop(columns="_merge")
    hosts = excluded["host_species"].fillna("").astype(str).str.strip()
    excluded["host_species"] = hosts.where(
        ~hosts.str.lower().isin(UNRECORDED_HOST_VALUES), UNRECORDED_HOST
    )
    human = excluded.loc[excluded["host_species"] == HUMAN_HOST_SPECIES, "profile_id"]
    if not human.empty:
        msg = (
            f"{human.nunique():,} human-host profile(s) are missing from the profile "
            f"metadata, e.g. {sorted(human.unique())[:3]}"
        )
        raise ValueError(msg)

    joined = occurrences.loc[in_corpus].merge(
        profiles,
        left_on="profile_id",
        right_index=True,
        how="left",
        suffixes=("_taxa", ""),
        validate="many_to_one",
    )
    return joined, excluded


def _counted(count: int, singular: str, plural: str) -> str:
    """Return `count` with the noun form that agrees with it."""
    return f"{count:,} {singular if count == 1 else plural}"


def excluded_host_summary(excluded: pd.DataFrame) -> str:
    """Return a sentence counting the set-aside occurrences by direction and host."""
    if excluded.empty:
        return "No occurrence lies outside the human-host corpus."
    directions = excluded["direction"].value_counts()
    direction_text = ", ".join(
        f"{int(directions[direction]):,} {direction}" for direction in sorted(directions.index)
    )
    hosts = excluded["host_species"].value_counts()
    host_text = ", ".join(
        f"{host} {int(count):,}"
        for host, count in sorted(hosts.items(), key=lambda item: (-item[1], item[0]))
    )
    verb = "comes" if len(excluded) == 1 else "come"
    return (
        f"{_counted(len(excluded), 'occurrence', 'occurrences')} ({direction_text}) in "
        f"{_counted(excluded['profile_id'].nunique(), 'profile', 'profiles')} from "
        f"{_counted(excluded['study_id'].nunique(), 'study', 'studies')} {verb} from "
        f"experiments whose host is not recorded as {HUMAN_HOST_SPECIES}: {host_text}."
    )


def excluded_host_note(excluded: pd.DataFrame, direction: str) -> str:
    """Return a figure note stating how many `direction` occurrences were set aside."""
    selected = excluded[excluded["direction"] == direction]
    if selected.empty:
        return "Human-host profiles only."
    count = len(selected)
    profile_count = selected["profile_id"].nunique()
    hosts = (
        "non-human or unrecorded hosts"
        if selected["host_species"].eq(UNRECORDED_HOST).any()
        else "non-human hosts"
    )
    occurrences = _counted(
        count, f"{direction}-profile occurrence", f"{direction}-profile occurrences"
    )
    verb = "is" if count == 1 else "are"
    return (
        f"Human-host profiles only: {occurrences} "
        f"({_counted(profile_count, 'profile', 'profiles')}) from {hosts} {verb} excluded."
    )
