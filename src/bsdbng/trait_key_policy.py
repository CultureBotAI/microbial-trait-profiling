"""Compatibility shim forwarding bsdbng.trait_key_policy to mtp.trait_key_policy."""

import sys
from typing import Any

from mtp.trait_key_policy import (
    DEFAULT_TRAIT_KEY_POLICY,
    FACULTATIVE_EXTRA_MEMBER_IDS,
    FACULTATIVE_MEMBER_IDS,
    FACULTATIVE_RAW_LABEL,
    FACULTATIVE_SOURCE_VOCABULARIES,
    FACULTATIVE_TRAIT_ID,
    FACULTATIVE_UNION_LABEL,
    TRAIT_KEY_POLICIES,
    TRAIT_KEY_POLICY_DEFINITIONS,
    TRAIT_POLICY_METADATA_COLUMNS,
    RecordedMemberEvidence,
    SourceVocabulary,
    TraitKeyPolicy,
    TraitKeyPolicyDefinition,
    canonical_object,
    canonical_trait_object_id,
    display_trait_label,
    expected_trait_display_label,
    get_trait_key_policy_definition,
    parse_trait_key_policy,
    policy_member_ids_for_trait,
    policy_metadata,
    read_recorded_member_evidence,
    read_recorded_member_ids,
    validate_trait_key_policy_definitions,
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

import mtp.trait_key_policy as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.trait_key_policy", run_name="__main__")
