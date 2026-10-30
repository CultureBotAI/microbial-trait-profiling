"""Keep BugSigDB taxon occurrences inside the human-host profile corpus (#587).

`data/derived/bugsigdb_taxa.tsv` lists the reported taxa of every BugSigDB host.
The profile metadata of the trait-profile investigations,
`trait-visualization/study_trait_profile_study_metadata.tsv`, keeps the profiles
whose experiment host is recorded as Homo sapiens in `data/raw/full_dump.csv`
(`scripts/create_study_trait_profile_heatmap.py`). An occurrence outside that
corpus has no condition or body site, so it is set aside and counted instead of
being drawn as an `unknown` condition.

The corpus also keeps only bacterial community profiles: experiments sequenced by
16S amplicon or shotgun (`WMS`) (#742).
"""

from __future__ import annotations

from collections.abc import Collection

import pandas as pd

__all__ = [
    "CORPUS_DESCRIPTION",
    "CORPUS_SEQUENCING_TYPES",
    "EXCLUDED_BY_HOST",
    "EXCLUDED_BY_SEQUENCING_TYPE",
    "HUMAN_HOST_SPECIES",
    "UNRECORDED_HOST",
    "UNRECORDED_SEQUENCING_TYPE",
    "excluded_host_note",
    "excluded_host_summary",
    "experiment_in_profile_population",
    "join_human_host_profiles",
]

HUMAN_HOST_SPECIES = "Homo sapiens"
# BugSigDB "Sequencing type" values the corpus keeps. ITS / ITS2 and 18S profile fungi
# and other eukaryotes, which KG-Microbe traits barely cover; PCR is a targeted assay,
# not a community profile; an experiment with no recorded method is neither (#742).
CORPUS_SEQUENCING_TYPES = ("16S", "WMS")
UNRECORDED_HOST = "not recorded"
UNRECORDED_SEQUENCING_TYPE = "not recorded"
# Values that mean the source has no host (or sequencing type) for an experiment;
# "unknown" is what the corpus's full-dump parser writes for an empty field.
UNRECORDED_HOST_VALUES = {"", "unknown", "na", "nan"}
EXPERIMENT_KEY = ["study_id", "experiment_id"]
# Why an occurrence lies outside the corpus, in the `exclusion_reason` column
EXCLUDED_BY_HOST = "host"
EXCLUDED_BY_SEQUENCING_TYPE = "sequencing_type"
CORPUS_DESCRIPTION = "Human-host 16S and shotgun profiles only"


def experiment_in_profile_population(
    host: object,
    sequencing_type: object,
    *,
    host_species: str | None,
    sequencing_types: Collection[str] | None,
) -> bool:
    """Shared experiment selection for the canonical builder and coverage audit (#865)."""
    recorded_host = "" if host is None else str(host).strip()
    recorded_method = "" if sequencing_type is None else str(sequencing_type).strip()
    return (host_species is None or recorded_host == host_species.strip()) and (
        sequencing_types is None or recorded_method in {value.strip() for value in sequencing_types}
    )


def _recorded(values: pd.Series, unrecorded: str) -> pd.Series:
    """Return `values` stripped, with every empty or unknown value as `unrecorded`."""
    text = values.fillna("").astype(str).str.strip()
    return text.where(~text.str.lower().isin(UNRECORDED_HOST_VALUES), unrecorded)


def join_human_host_profiles(
    occurrences: pd.DataFrame,
    profiles: pd.DataFrame,
    experiment_hosts: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Join profile metadata to occurrences in the corpus; return the rest with their reason.

    `occurrences` has one row per taxon and profile, with `study_id`, `experiment_id`
    and `profile_id`; `profiles` is the corpus metadata indexed by profile ID;
    `experiment_hosts` gives `host_species` and `sequencing_type` per study and
    experiment, read from the file the corpus filter read (#713). Columns that both
    tables carry keep the occurrence values with a `_taxa` suffix.

    An occurrence outside the corpus must come from an experiment whose host is not
    recorded as Homo sapiens (another species, or no host at all, reported as "not
    recorded"), or from a Homo sapiens experiment whose sequencing type is not one of
    `CORPUS_SEQUENCING_TYPES` (#742). `exclusion_reason` says which. An occurrence of a
    Homo sapiens 16S or WMS experiment, or of an experiment with no host row, raises,
    so a key mismatch cannot pass as an exclusion.
    """
    in_corpus = occurrences["profile_id"].isin(profiles.index)
    excluded = occurrences.loc[~in_corpus].merge(
        experiment_hosts[[*EXPERIMENT_KEY, "host_species", "sequencing_type"]],
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
    excluded["host_species"] = _recorded(excluded["host_species"], UNRECORDED_HOST)
    excluded["sequencing_type"] = _recorded(excluded["sequencing_type"], UNRECORDED_SEQUENCING_TYPE)
    human = excluded["host_species"] == HUMAN_HOST_SPECIES
    in_corpus_method = excluded["sequencing_type"].isin(CORPUS_SEQUENCING_TYPES)
    missing = excluded.loc[human & in_corpus_method, "profile_id"]
    if not missing.empty:
        msg = (
            f"{missing.nunique():,} human-host profile(s) sequenced by "
            f"{' or '.join(CORPUS_SEQUENCING_TYPES)} are missing from the profile "
            f"metadata, e.g. {sorted(missing.unique())[:3]}"
        )
        raise ValueError(msg)
    excluded["exclusion_reason"] = EXCLUDED_BY_HOST
    excluded.loc[human, "exclusion_reason"] = EXCLUDED_BY_SEQUENCING_TYPE

    # The kept occurrences are checked too: profile metadata written before #742 would
    # otherwise count other methods under a "16S and shotgun only" note (#752)
    kept = occurrences.loc[in_corpus, [*EXPERIMENT_KEY, "profile_id"]].merge(
        experiment_hosts[[*EXPERIMENT_KEY, "host_species", "sequencing_type"]],
        on=EXPERIMENT_KEY,
        how="left",
        validate="many_to_one",
    )
    outside = kept.loc[
        (_recorded(kept["host_species"], UNRECORDED_HOST) != HUMAN_HOST_SPECIES)
        | ~_recorded(kept["sequencing_type"], UNRECORDED_SEQUENCING_TYPE).isin(
            CORPUS_SEQUENCING_TYPES
        ),
        "profile_id",
    ]
    if not outside.empty:
        msg = (
            f"{outside.nunique():,} profile(s) in the profile metadata are not "
            f"{HUMAN_HOST_SPECIES} {' or '.join(CORPUS_SEQUENCING_TYPES)} experiments, e.g. "
            f"{sorted(outside.unique())[:3]}; regenerate the metadata with "
            "`just trait-heatmaps` (#742)"
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


def _value_counts_text(values: pd.Series) -> str:
    """Return each value with its count, the commonest first."""
    counts = values.value_counts()
    return ", ".join(
        f"{value} {int(count):,}"
        for value, count in sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    )


def excluded_host_summary(excluded: pd.DataFrame) -> str:
    """Return a sentence counting the set-aside occurrences by direction and reason."""
    if excluded.empty:
        return "No occurrence lies outside the human-host corpus."
    directions = excluded["direction"].value_counts()
    direction_text = ", ".join(
        f"{int(directions[direction]):,} {direction}" for direction in sorted(directions.index)
    )
    by_host = excluded[excluded["exclusion_reason"] == EXCLUDED_BY_HOST]
    by_method = excluded[excluded["exclusion_reason"] == EXCLUDED_BY_SEQUENCING_TYPE]
    reasons = []
    if not by_host.empty:
        reasons.append(
            (
                len(by_host),
                f"experiments whose host is not recorded as {HUMAN_HOST_SPECIES}",
                _value_counts_text(by_host["host_species"]),
            )
        )
    if not by_method.empty:
        reasons.append(
            (
                len(by_method),
                f"{HUMAN_HOST_SPECIES} experiments with a sequencing method other than "
                f"{' or '.join(CORPUS_SEQUENCING_TYPES)}, or none recorded",
                _value_counts_text(by_method["sequencing_type"]),
            )
        )
    lead = (
        f"{_counted(len(excluded), 'occurrence', 'occurrences')} ({direction_text}) in "
        f"{_counted(excluded['profile_id'].nunique(), 'profile', 'profiles')} from "
        f"{_counted(excluded['study_id'].nunique(), 'study', 'studies')}"
    )
    if len(reasons) == 1:
        _, source, counts = reasons[0]
        verb = "comes" if len(excluded) == 1 else "come"
        return f"{lead} {verb} from {source}: {counts}."
    parts = [
        f"{count:,} {'comes' if count == 1 else 'come'} from {source} ({counts})"
        if index == 0
        else f"{count:,} from {source} ({counts})"
        for index, (count, source, counts) in enumerate(reasons)
    ]
    return f"{lead} lie outside the corpus: {', and '.join(parts)}."


def excluded_host_note(excluded: pd.DataFrame, direction: str) -> str:
    """Return a figure note stating how many `direction` occurrences were set aside."""
    selected = excluded[excluded["direction"] == direction]
    if selected.empty:
        return f"{CORPUS_DESCRIPTION}."
    count = len(selected)
    profile_count = selected["profile_id"].nunique()
    sources = []
    by_host = selected[selected["exclusion_reason"] == EXCLUDED_BY_HOST]
    if not by_host.empty:
        sources.append(
            "non-human or unrecorded hosts"
            if by_host["host_species"].eq(UNRECORDED_HOST).any()
            else "non-human hosts"
        )
    by_method = selected[selected["exclusion_reason"] == EXCLUDED_BY_SEQUENCING_TYPE]
    if not by_method.empty:
        sources.append(
            "other or unrecorded sequencing methods"
            if by_method["sequencing_type"].eq(UNRECORDED_SEQUENCING_TYPE).any()
            else "other sequencing methods"
        )
    occurrences = _counted(
        count, f"{direction}-profile occurrence", f"{direction}-profile occurrences"
    )
    verb = "is" if count == 1 else "are"
    return (
        f"{CORPUS_DESCRIPTION}: {occurrences} "
        f"({_counted(profile_count, 'profile', 'profiles')}) from {' and '.join(sources)} "
        f"{verb} excluded."
    )
