"""Versioned local trait-key semantics and declared source vocabularies (#1129, #1137).

The facultative union is an analysis mapping, not an upstream ontology merger.
Source labels document the METPO template mappings; they do not assert that a
particular matrix cell contains evidence from every listed source. Original KG
assertions, node labels, and source attribution remain authoritative evidence.

Each policy version in ``TRAIT_KEY_POLICY_DEFINITIONS`` freezes its canonical trait ID,
member CURIEs, raw/display labels, and declared source vocabularies. Any future change
to membership or source vocabularies requires adding a new policy version identifier
rather than mutating an existing version in place (#1137).
"""

from __future__ import annotations

import csv
import gzip
import json
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Literal, cast

TraitKeyPolicy = Literal["raw_object_v1", "facultative_union_v1"]
DEFAULT_TRAIT_KEY_POLICY: TraitKeyPolicy = "facultative_union_v1"
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


@dataclass(frozen=True)
class TraitKeyPolicyDefinition:
    """Immutable definition for a specific trait-key policy version (#1137)."""

    policy: TraitKeyPolicy
    canonical_trait_id: str
    member_ids: tuple[str, ...]
    raw_label: str
    display_label: str
    source_vocabularies: tuple[SourceVocabulary, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.policy, str) or not self.policy.strip():
            raise ValueError("Policy version identifier must be a non-empty string")
        if not isinstance(self.canonical_trait_id, str) or not self.canonical_trait_id.strip():
            raise ValueError("canonical_trait_id must be a non-empty string")
        if (
            not isinstance(self.member_ids, tuple)
            or not self.member_ids
            or not all(isinstance(item, str) and item.strip() for item in self.member_ids)
            or len(set(self.member_ids)) != len(self.member_ids)
        ):
            raise ValueError("member_ids must be a non-empty tuple of unique non-empty CURIEs")
        if self.canonical_trait_id not in self.member_ids:
            raise ValueError("canonical_trait_id must be included in member_ids")
        if not isinstance(self.raw_label, str) or not self.raw_label.strip():
            raise ValueError("raw_label must be a non-empty string")
        if not isinstance(self.display_label, str) or not self.display_label.strip():
            raise ValueError("display_label must be a non-empty string")
        if not isinstance(self.source_vocabularies, tuple) or not all(
            isinstance(item, SourceVocabulary) for item in self.source_vocabularies
        ):
            raise ValueError("source_vocabularies must be a tuple of SourceVocabulary records")

    @property
    def version(self) -> TraitKeyPolicy:
        """Return the frozen policy version identifier."""
        return self.policy

    @property
    def extra_member_ids(self) -> frozenset[str]:
        """Return non-canonical union member CURIEs for this policy version."""
        return frozenset(member for member in self.member_ids if member != self.canonical_trait_id)


# METPO source template c6f4947, rows 80/83/87; OMP:0000087 has
# "facultative aerobe" as an exact synonym of facultative anaerobe (2024-03-25).
# BacDive distinguishes its source labels; BactoTraits already combines them.
_V1_SOURCE_VOCABULARIES: tuple[SourceVocabulary, ...] = (
    SourceVocabulary("METPO:1000605", "infores:bacdive", "facultative anaerobe"),
    SourceVocabulary("METPO:1000605", "infores:madin_etal", "facultative"),
    SourceVocabulary("METPO:1000608", "infores:bacdive", "facultative aerobe"),
    SourceVocabulary("METPO:1000608", "infores:madin_etal", "facultative"),
    SourceVocabulary("METPO:1000612", "infores:bactotraits", "Ox_facultative_aerobe_anaerobe"),
)

_FROZEN_V1_DEFINITIONS: dict[TraitKeyPolicy, TraitKeyPolicyDefinition] = {
    "raw_object_v1": TraitKeyPolicyDefinition(
        policy="raw_object_v1",
        canonical_trait_id="METPO:1000608",
        member_ids=("METPO:1000608",),
        raw_label="facultatively aerobic",
        display_label="facultatively aerobic",
        source_vocabularies=_V1_SOURCE_VOCABULARIES,
    ),
    "facultative_union_v1": TraitKeyPolicyDefinition(
        policy="facultative_union_v1",
        canonical_trait_id="METPO:1000608",
        member_ids=("METPO:1000605", "METPO:1000608", "METPO:1000612"),
        raw_label="facultatively aerobic",
        display_label="facultative (aerobe or anaerobe)",
        source_vocabularies=_V1_SOURCE_VOCABULARIES,
    ),
}

TRAIT_KEY_POLICY_DEFINITIONS: Mapping[TraitKeyPolicy, TraitKeyPolicyDefinition] = MappingProxyType(
    dict(_FROZEN_V1_DEFINITIONS)
)
TRAIT_KEY_POLICIES: tuple[TraitKeyPolicy, ...] = tuple(TRAIT_KEY_POLICY_DEFINITIONS.keys())

_UNION_V1_DEFINITION = TRAIT_KEY_POLICY_DEFINITIONS["facultative_union_v1"]
FACULTATIVE_TRAIT_ID = _UNION_V1_DEFINITION.canonical_trait_id
FACULTATIVE_RAW_LABEL = _UNION_V1_DEFINITION.raw_label
FACULTATIVE_UNION_LABEL = _UNION_V1_DEFINITION.display_label
FACULTATIVE_MEMBER_IDS = _UNION_V1_DEFINITION.member_ids
FACULTATIVE_EXTRA_MEMBER_IDS = _UNION_V1_DEFINITION.extra_member_ids
FACULTATIVE_SOURCE_VOCABULARIES = _UNION_V1_DEFINITION.source_vocabularies


def validate_trait_key_policy_definitions(
    definitions: (
        Mapping[TraitKeyPolicy, TraitKeyPolicyDefinition] | Mapping[str, TraitKeyPolicyDefinition]
    ) = TRAIT_KEY_POLICY_DEFINITIONS,
    *,
    allow_additional_versions: bool = False,
) -> None:
    """Verify that frozen v1 policy definitions are intact and new definitions use new versions."""
    typed_definitions = cast(Mapping[str, TraitKeyPolicyDefinition], definitions)
    for required_version, frozen_def in _FROZEN_V1_DEFINITIONS.items():
        if required_version not in typed_definitions:
            raise ValueError(f"Missing frozen trait_key_policy version {required_version!r}")
        candidate = typed_definitions[required_version]
        if candidate != frozen_def:
            raise ValueError(
                f"Frozen trait_key_policy definition {required_version!r} was modified in place; "
                "introduce a new policy version identifier instead"
            )
    for key, definition in typed_definitions.items():
        if not isinstance(definition, TraitKeyPolicyDefinition):
            raise ValueError(f"Invalid policy definition for {key!r}: {definition!r}")
        if definition.policy != key:
            raise ValueError(
                f"Policy key {key!r} does not match definition.policy {definition.policy!r}"
            )
        if not allow_additional_versions and key not in _FROZEN_V1_DEFINITIONS:
            raise ValueError(
                f"Unregistered trait_key_policy version {key!r}; update TRAIT_KEY_POLICIES "
                "and TraitKeyPolicy explicitly when adding a new version"
            )


validate_trait_key_policy_definitions(TRAIT_KEY_POLICY_DEFINITIONS)


def parse_trait_key_policy(value: object) -> TraitKeyPolicy:
    """Refuse unknown policies instead of silently applying current semantics."""
    if isinstance(value, str) and value in TRAIT_KEY_POLICY_DEFINITIONS:
        return cast(TraitKeyPolicy, value)
    supported = tuple(TRAIT_KEY_POLICY_DEFINITIONS.keys())
    raise ValueError(f"trait_key_policy must be one of {supported!r}: {value!r}")


def get_trait_key_policy_definition(trait_key_policy: object) -> TraitKeyPolicyDefinition:
    """Return the version-frozen definition for ``trait_key_policy`` (#1137)."""
    policy = parse_trait_key_policy(trait_key_policy)
    return TRAIT_KEY_POLICY_DEFINITIONS[policy]


def canonical_trait_object_id(
    obj: str, *, trait_key_policy: TraitKeyPolicy = DEFAULT_TRAIT_KEY_POLICY
) -> str:
    """Map declared union members using the version-frozen definition for ``trait_key_policy``."""
    definition = get_trait_key_policy_definition(trait_key_policy)
    return definition.canonical_trait_id if obj in definition.member_ids else obj


def canonical_object(obj: str, *, trait_key_policy: TraitKeyPolicy) -> str:
    """Map only declared union members; raw policy preserves the exact object ID."""
    return canonical_trait_object_id(obj, trait_key_policy=trait_key_policy)


def policy_member_ids_for_trait(
    trait_id: str, *, trait_key_policy: TraitKeyPolicy = DEFAULT_TRAIT_KEY_POLICY
) -> tuple[str, ...]:
    """Return the version-frozen member object CURIEs for a canonical trait or column ID (#1137)."""
    definition = get_trait_key_policy_definition(trait_key_policy)
    cleaned_id = trait_id.strip()
    obj_id = cleaned_id.split("|", 1)[1].strip() if "|" in cleaned_id else cleaned_id
    if canonical_trait_object_id(obj_id, trait_key_policy=definition.policy) != obj_id:
        raise ValueError(
            f"Feature metadata object {trait_id!r} is not canonical under {definition.policy}"
        )
    if obj_id == definition.canonical_trait_id:
        return definition.member_ids
    return (obj_id,)


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
        except json.JSONDecodeError as exc:
            raise ValueError(f"Invalid trait_member_ids_json: {text!r}") from exc
        if not isinstance(parsed, list) or not all(isinstance(item, str) for item in parsed):
            raise ValueError(f"Invalid trait_member_ids_json: {text!r}")
        return tuple(item.strip() for item in parsed if item.strip())
    return tuple(str(item).strip() for item in member_ids if str(item).strip())


def expected_trait_display_label(
    trait_id: str,
    recorded_label: str | None = None,
    *,
    trait_key_policy: TraitKeyPolicy = DEFAULT_TRAIT_KEY_POLICY,
    member_ids: Sequence[str] | str | None = None,
) -> str:
    """Return the version-defined display label for ``trait_id`` under ``trait_key_policy``."""
    definition = get_trait_key_policy_definition(trait_key_policy)
    cleaned_id = trait_id.strip()
    obj_id = cleaned_id.split("|", 1)[1].strip() if "|" in cleaned_id else cleaned_id
    if recorded_label is None:
        if obj_id != definition.canonical_trait_id:
            return trait_id
        effective_label = definition.raw_label
    else:
        effective_label = recorded_label
    if obj_id != definition.canonical_trait_id or not effective_label:
        return effective_label
    effective_members = (
        definition.member_ids if member_ids is None else _parse_member_ids(member_ids)
    )
    if not definition.extra_member_ids.intersection(effective_members):
        return effective_label
    if definition.display_label in effective_label:
        return effective_label
    if definition.raw_label in effective_label:
        return effective_label.replace(definition.raw_label, definition.display_label)
    capitalized_raw = definition.raw_label.capitalize()
    if capitalized_raw in effective_label:
        return effective_label.replace(capitalized_raw, definition.display_label)
    return definition.display_label


def display_trait_label(
    trait_id: str,
    recorded_label: str,
    member_ids: Sequence[str] | str | None = None,
    *,
    trait_key_policy: TraitKeyPolicy = DEFAULT_TRAIT_KEY_POLICY,
) -> str:
    """Return the render-time display label for a trait from its recorded member IDs.

    Recorded column metadata (`object_name`, `column_label`) always preserves the raw
    KG-Microbe node name so historical core verification and column sort order remain
    stable (#1231, #1233). Downstream consumers display ``facultative (aerobe or anaerobe)``
    only when ``trait_id`` (or its object suffix after ``|``) is ``METPO:1000608`` and
    recorded ``member_ids`` include ``METPO:1000605`` or ``METPO:1000612`` (#1138, #1232).
    """
    parsed_members = _parse_member_ids(member_ids)
    return expected_trait_display_label(
        trait_id,
        recorded_label,
        trait_key_policy=trait_key_policy,
        member_ids=parsed_members,
    )


@dataclass(frozen=True)
class RecordedMemberEvidence:
    """Parsed member aliases and the number of usable source rows, not changed labels."""

    member_ids_by_trait: dict[str, tuple[str, ...]]
    member_id_rows: int


def read_recorded_member_evidence(column_metadata_path: Path | None) -> RecordedMemberEvidence:
    """Read declared members and count rows with nonempty members and an identifier.

    A column/object alias pair counts once; duplicate rows still count separately.
    Legacy metadata without member declarations retains the empty-mapping fallback.
    """
    if column_metadata_path is None:
        return RecordedMemberEvidence({}, 0)
    if not column_metadata_path.exists():
        raise FileNotFoundError(f"Column metadata TSV not found: {column_metadata_path}")
    member_ids_by_trait: dict[str, tuple[str, ...]] = {}
    member_id_rows = 0
    open_fn = gzip.open if column_metadata_path.suffix == ".gz" else open
    with open_fn(column_metadata_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "trait_member_ids_json" not in reader.fieldnames:
            return RecordedMemberEvidence({}, 0)
        for row in reader:
            members = _parse_member_ids(row.get("trait_member_ids_json"))
            if not members:
                continue
            column_id = (row.get("column_id") or "").strip()
            object_id = (row.get("object_id") or "").strip()
            if column_id or object_id:
                member_id_rows += 1
            if column_id:
                member_ids_by_trait[column_id] = members
            if object_id and object_id not in member_ids_by_trait:
                member_ids_by_trait[object_id] = members
    return RecordedMemberEvidence(member_ids_by_trait, member_id_rows)


def read_recorded_member_ids(column_metadata_path: Path | None) -> dict[str, tuple[str, ...]]:
    """Read column_id and object_id -> member_ids, preserving legacy fallback semantics."""
    return read_recorded_member_evidence(column_metadata_path).member_ids_by_trait


def policy_metadata(obj: str, *, trait_key_policy: TraitKeyPolicy) -> tuple[str, str, str]:
    """Return policy, member-object IDs, and declared source mappings as stable TSV cells.

    Member IDs exclude the predicate: a qualified or negated feature retains its
    original predicate separately. Unknown source vocabularies are an empty list,
    never invented from the canonical node label. Both JSON arrays are declarations,
    not claims that every member or source actually contributed to a given run.
    """
    definition = get_trait_key_policy_definition(trait_key_policy)
    if "|" in obj or canonical_trait_object_id(obj, trait_key_policy=definition.policy) != obj:
        raise ValueError(
            f"Feature metadata object {obj!r} is not canonical under {definition.policy}"
        )
    members = policy_member_ids_for_trait(obj, trait_key_policy=definition.policy)
    mappings = [
        asdict(item) for item in definition.source_vocabularies if item.object_id in members
    ]
    return (
        definition.policy,
        json.dumps(members, separators=(",", ":")),
        json.dumps(mappings, sort_keys=True, separators=(",", ":")),
    )


__all__ = [
    "DEFAULT_TRAIT_KEY_POLICY",
    "FACULTATIVE_EXTRA_MEMBER_IDS",
    "FACULTATIVE_MEMBER_IDS",
    "FACULTATIVE_RAW_LABEL",
    "FACULTATIVE_SOURCE_VOCABULARIES",
    "FACULTATIVE_TRAIT_ID",
    "FACULTATIVE_UNION_LABEL",
    "TRAIT_KEY_POLICIES",
    "TRAIT_KEY_POLICY_DEFINITIONS",
    "TRAIT_POLICY_METADATA_COLUMNS",
    "RecordedMemberEvidence",
    "SourceVocabulary",
    "TraitKeyPolicy",
    "TraitKeyPolicyDefinition",
    "canonical_object",
    "canonical_trait_object_id",
    "display_trait_label",
    "expected_trait_display_label",
    "get_trait_key_policy_definition",
    "parse_trait_key_policy",
    "policy_member_ids_for_trait",
    "policy_metadata",
    "read_recorded_member_evidence",
    "read_recorded_member_ids",
    "validate_trait_key_policy_definitions",
]
