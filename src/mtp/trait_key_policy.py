"""Versioned local trait-key semantics and declared source vocabularies (#1129).

The facultative union is an analysis mapping, not an upstream ontology merger.
Source labels document the METPO template mappings; they do not assert that a
particular matrix cell contains evidence from every listed source. Original KG
assertions, node labels, and source attribution remain authoritative evidence.
"""

from __future__ import annotations

import csv
import json
from collections.abc import Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Literal

TraitKeyPolicy = Literal["raw_object_v1", "facultative_union_v1"]
TRAIT_KEY_POLICIES: tuple[TraitKeyPolicy, ...] = ("raw_object_v1", "facultative_union_v1")
DEFAULT_TRAIT_KEY_POLICY: TraitKeyPolicy = "facultative_union_v1"
FACULTATIVE_TRAIT_ID = "METPO:1000608"
FACULTATIVE_RAW_LABEL = "facultatively aerobic"
FACULTATIVE_UNION_LABEL = "facultative (aerobe or anaerobe)"
FACULTATIVE_MEMBER_IDS = ("METPO:1000605", "METPO:1000608", "METPO:1000612")
FACULTATIVE_EXTRA_MEMBER_IDS = frozenset({"METPO:1000605", "METPO:1000612"})
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


def _parse_member_ids(member_ids: Sequence[str] | str | None) -> tuple[str, ...]:
    """Parse member IDs from a sequence or a recorded trait_member_ids_json cell."""
    if member_ids is None:
        return ()
    if isinstance(member_ids, str):
        text = member_ids.strip()
        if not text:
            return ()
        try:
            parsed = json.loads(text)
        except json.JSONDecodeError:
            return ()
        if not isinstance(parsed, list) or not all(isinstance(item, str) for item in parsed):
            return ()
        return tuple(item.strip() for item in parsed if item.strip())
    return tuple(str(item).strip() for item in member_ids if str(item).strip())


def display_trait_label(
    trait_id: str,
    recorded_label: str,
    member_ids: Sequence[str] | str | None = None,
) -> str:
    """Return the render-time display label for a trait from its recorded member IDs.

    Recorded column metadata (`object_name`, `column_label`) always preserves the raw
    KG-Microbe node name so historical core verification and column sort order remain
    stable (#1231, #1233). Downstream consumers display ``facultative (aerobe or anaerobe)``
    only when ``trait_id`` (or its object suffix after ``|``) is ``METPO:1000608`` and
    recorded ``member_ids`` include ``METPO:1000605`` or ``METPO:1000612`` (#1138, #1232).
    """
    cleaned_id = trait_id.strip()
    obj_id = cleaned_id.split("|", 1)[1].strip() if "|" in cleaned_id else cleaned_id
    if obj_id != FACULTATIVE_TRAIT_ID or not recorded_label:
        return recorded_label
    parsed_members = _parse_member_ids(member_ids)
    if not FACULTATIVE_EXTRA_MEMBER_IDS.intersection(parsed_members):
        return recorded_label
    if FACULTATIVE_UNION_LABEL in recorded_label:
        return recorded_label
    if FACULTATIVE_RAW_LABEL in recorded_label:
        return recorded_label.replace(FACULTATIVE_RAW_LABEL, FACULTATIVE_UNION_LABEL)
    capitalized_raw = FACULTATIVE_RAW_LABEL.capitalize()
    if capitalized_raw in recorded_label:
        return recorded_label.replace(capitalized_raw, FACULTATIVE_UNION_LABEL)
    return FACULTATIVE_UNION_LABEL


def read_recorded_member_ids(column_metadata_path: Path | None) -> dict[str, tuple[str, ...]]:
    """Read column_id and object_id -> member_ids from a column metadata TSV if present."""
    if column_metadata_path is None or not column_metadata_path.exists():
        return {}
    member_ids_by_trait: dict[str, tuple[str, ...]] = {}
    with column_metadata_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "trait_member_ids_json" not in reader.fieldnames:
            return {}
        for row in reader:
            members = _parse_member_ids(row.get("trait_member_ids_json"))
            if not members:
                continue
            column_id = (row.get("column_id") or "").strip()
            object_id = (row.get("object_id") or "").strip()
            if column_id:
                member_ids_by_trait[column_id] = members
            if object_id and object_id not in member_ids_by_trait:
                member_ids_by_trait[object_id] = members
    return member_ids_by_trait


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
