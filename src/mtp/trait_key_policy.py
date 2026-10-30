"""Versioned local trait-key semantics and declared source vocabularies (#1129).

The facultative union is an analysis mapping, not an upstream ontology merger.
Source labels document the METPO template mappings; they do not assert that a
particular matrix cell contains evidence from every listed source. Original KG
assertions, node labels, and source attribution remain authoritative evidence.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from typing import Literal

TraitKeyPolicy = Literal["raw_object_v1", "facultative_union_v1"]
TRAIT_KEY_POLICIES: tuple[TraitKeyPolicy, ...] = ("raw_object_v1", "facultative_union_v1")
DEFAULT_TRAIT_KEY_POLICY: TraitKeyPolicy = "facultative_union_v1"
FACULTATIVE_TRAIT_ID = "METPO:1000608"
FACULTATIVE_UNION_LABEL = "facultative (aerobe or anaerobe)"
FACULTATIVE_MEMBER_IDS = ("METPO:1000605", "METPO:1000608", "METPO:1000612")
TRAIT_POLICY_METADATA_COLUMNS = (
    "trait_key_policy",
    "trait_member_ids_json",
    "trait_source_vocabularies_json",
)


@dataclass(frozen=True)
class SourceVocabulary:
    """A declared source-label mapping, not an observed contribution count."""

    object_id: str
    source_id: str
    source_label: str


# METPO source template c6f4947, rows 80/83/87; OMP:0000087 has
# "facultative aerobe" as an exact synonym of facultative anaerobe (2024-03-25).
# BacDive distinguishes its source labels; BactoTraits already combines them.
FACULTATIVE_SOURCE_VOCABULARIES = (
    SourceVocabulary("METPO:1000605", "infores:bacdive", "facultative anaerobe"),
    SourceVocabulary("METPO:1000605", "infores:madin_etal", "facultative"),
    SourceVocabulary("METPO:1000608", "infores:bacdive", "facultative aerobe"),
    SourceVocabulary("METPO:1000608", "infores:madin_etal", "facultative"),
    SourceVocabulary("METPO:1000612", "infores:bactotraits", "Ox_facultative_aerobe_anaerobe"),
)


def parse_trait_key_policy(value: object) -> TraitKeyPolicy:
    """Refuse unknown policies instead of silently applying current semantics."""
    for policy in TRAIT_KEY_POLICIES:
        if value == policy:
            return policy
    raise ValueError(f"trait_key_policy must be one of {TRAIT_KEY_POLICIES!r}: {value!r}")


def canonical_object(obj: str, *, trait_key_policy: TraitKeyPolicy) -> str:
    """Map only declared union members; raw policy preserves the exact object ID."""
    if trait_key_policy == "raw_object_v1":
        return obj
    if trait_key_policy != "facultative_union_v1":
        raise ValueError(f"Unknown trait_key_policy: {trait_key_policy!r}")
    return FACULTATIVE_TRAIT_ID if obj in FACULTATIVE_MEMBER_IDS else obj


def canonical_object_name(obj: str, raw_name: str, *, trait_key_policy: TraitKeyPolicy) -> str:
    """Label the facultative union under facultative_union_v1; keep raw name under raw_object_v1."""
    policy = parse_trait_key_policy(trait_key_policy)
    if policy == "facultative_union_v1" and obj == FACULTATIVE_TRAIT_ID and raw_name:
        return FACULTATIVE_UNION_LABEL
    return raw_name


def policy_metadata(obj: str, *, trait_key_policy: TraitKeyPolicy) -> tuple[str, str, str]:
    """Return policy, member-object IDs, and declared source mappings as stable TSV cells.

    Member IDs exclude the predicate: a qualified or negated feature retains its
    original predicate separately. Unknown source vocabularies are an empty list,
    never invented from the canonical node label. Both JSON arrays are declarations,
    not claims that every member or source actually contributed to a given run.
    """
    policy = parse_trait_key_policy(trait_key_policy)
    if canonical_object(obj, trait_key_policy=policy) != obj:
        raise ValueError(f"Feature metadata object {obj!r} is not canonical under {policy}")
    members = (
        FACULTATIVE_MEMBER_IDS
        if policy == "facultative_union_v1" and obj == FACULTATIVE_TRAIT_ID
        else (obj,)
    )
    mappings = [
        asdict(item) for item in FACULTATIVE_SOURCE_VOCABULARIES if item.object_id in members
    ]
    return (
        policy,
        json.dumps(members, separators=(",", ":")),
        json.dumps(mappings, sort_keys=True, separators=(",", ":")),
    )
