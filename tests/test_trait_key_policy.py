"""Tests for version-frozen trait-key policy definitions and lookups (#1137)."""

from __future__ import annotations

import json
from collections.abc import Mapping
from dataclasses import FrozenInstanceError, replace
from pathlib import Path
from types import MappingProxyType
from typing import cast

import pytest

import mtp.trait_key_policy as policy_mod
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
    SourceVocabulary,
    TraitKeyPolicy,
    TraitKeyPolicyDefinition,
    canonical_object,
    canonical_trait_object_id,
    display_trait_label,
    expected_trait_display_label,
    get_trait_key_policy_definition,
    policy_member_ids_for_trait,
    policy_metadata,
    validate_trait_key_policy_definitions,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
RAW: TraitKeyPolicy = "raw_object_v1"
UNION: TraitKeyPolicy = "facultative_union_v1"
EXPECTED_V1_VOCABULARIES = (
    SourceVocabulary("METPO:1000605", "infores:bacdive", "facultative anaerobe"),
    SourceVocabulary("METPO:1000605", "infores:madin_etal", "facultative"),
    SourceVocabulary("METPO:1000608", "infores:bacdive", "facultative aerobe"),
    SourceVocabulary("METPO:1000608", "infores:madin_etal", "facultative"),
    SourceVocabulary("METPO:1000612", "infores:bactotraits", "Ox_facultative_aerobe_anaerobe"),
)


def test_policy_definitions_table_is_frozen_and_complete() -> None:
    """TRAIT_KEY_POLICY_DEFINITIONS is an immutable mapping covering every TraitKeyPolicy."""
    assert isinstance(TRAIT_KEY_POLICY_DEFINITIONS, MappingProxyType)
    assert tuple(TRAIT_KEY_POLICY_DEFINITIONS.keys()) == TRAIT_KEY_POLICIES
    with pytest.raises(TypeError):
        cast(dict[str, object], TRAIT_KEY_POLICY_DEFINITIONS)[RAW] = TRAIT_KEY_POLICY_DEFINITIONS[
            RAW
        ]

    raw_def = TRAIT_KEY_POLICY_DEFINITIONS[RAW]
    union_def = TRAIT_KEY_POLICY_DEFINITIONS[UNION]
    assert raw_def == get_trait_key_policy_definition(RAW)
    assert union_def == get_trait_key_policy_definition(UNION)

    assert raw_def.policy == RAW
    assert raw_def.version == RAW
    assert raw_def.canonical_trait_id == "METPO:1000608"
    assert raw_def.member_ids == ("METPO:1000608",)
    assert raw_def.extra_member_ids == frozenset()
    assert raw_def.raw_label == "facultatively aerobic"
    assert raw_def.display_label == "facultatively aerobic"
    assert raw_def.source_vocabularies == EXPECTED_V1_VOCABULARIES

    assert union_def.policy == UNION
    assert union_def.version == UNION
    assert union_def.canonical_trait_id == "METPO:1000608"
    assert union_def.member_ids == ("METPO:1000605", "METPO:1000608", "METPO:1000612")
    assert union_def.extra_member_ids == frozenset({"METPO:1000605", "METPO:1000612"})
    assert union_def.raw_label == "facultatively aerobic"
    assert union_def.display_label == "facultative (aerobe or anaerobe)"
    assert union_def.source_vocabularies == EXPECTED_V1_VOCABULARIES

    field_name = "member_ids"
    with pytest.raises(FrozenInstanceError):
        setattr(union_def, field_name, ("METPO:1000608",))


def test_backwards_compatible_constants_derive_from_union_v1_definition() -> None:
    """Module-level constants remain equal to the frozen facultative_union_v1 definition."""
    union_def = TRAIT_KEY_POLICY_DEFINITIONS[UNION]
    assert union_def.policy == DEFAULT_TRAIT_KEY_POLICY
    assert union_def.canonical_trait_id == FACULTATIVE_TRAIT_ID
    assert union_def.raw_label == FACULTATIVE_RAW_LABEL
    assert union_def.display_label == FACULTATIVE_UNION_LABEL
    assert union_def.member_ids == FACULTATIVE_MEMBER_IDS
    assert union_def.extra_member_ids == FACULTATIVE_EXTRA_MEMBER_IDS
    assert union_def.source_vocabularies == FACULTATIVE_SOURCE_VOCABULARIES


@pytest.mark.parametrize("bad_policy", ["", "facultative_union_v2", None, 42])
def test_unknown_policy_version_is_rejected_everywhere(bad_policy: object) -> None:
    """Lookup helpers refuse unknown or unversioned policy values."""
    with pytest.raises(ValueError, match="trait_key_policy"):
        get_trait_key_policy_definition(bad_policy)
    with pytest.raises(ValueError, match="trait_key_policy"):
        canonical_trait_object_id(
            "METPO:1000605", trait_key_policy=cast(TraitKeyPolicy, bad_policy)
        )
    with pytest.raises(ValueError, match="trait_key_policy"):
        policy_member_ids_for_trait(
            "METPO:1000608", trait_key_policy=cast(TraitKeyPolicy, bad_policy)
        )
    with pytest.raises(ValueError, match="trait_key_policy"):
        expected_trait_display_label(
            "METPO:1000608",
            "facultatively aerobic",
            trait_key_policy=cast(TraitKeyPolicy, bad_policy),
        )


def test_canonical_member_and_label_helpers_resolve_by_policy_version() -> None:
    """Versioned helpers resolve canonical IDs, member lists, and display labels per policy."""
    for member in ("METPO:1000605", "METPO:1000608", "METPO:1000612"):
        assert canonical_trait_object_id(member, trait_key_policy=RAW) == member
        assert canonical_object(member, trait_key_policy=RAW) == member
        assert canonical_trait_object_id(member, trait_key_policy=UNION) == "METPO:1000608"
        assert canonical_object(member, trait_key_policy=UNION) == "METPO:1000608"

    assert canonical_trait_object_id("GO:0030245", trait_key_policy=RAW) == "GO:0030245"
    assert canonical_trait_object_id("GO:0030245", trait_key_policy=UNION) == "GO:0030245"

    assert policy_member_ids_for_trait("METPO:1000608", trait_key_policy=RAW) == ("METPO:1000608",)
    assert policy_member_ids_for_trait("METPO:1000605", trait_key_policy=RAW) == ("METPO:1000605",)
    assert policy_member_ids_for_trait("METPO:1000608", trait_key_policy=UNION) == (
        "METPO:1000605",
        "METPO:1000608",
        "METPO:1000612",
    )
    assert policy_member_ids_for_trait("METPO:2000021|METPO:1000608", trait_key_policy=UNION) == (
        "METPO:1000605",
        "METPO:1000608",
        "METPO:1000612",
    )
    assert policy_member_ids_for_trait("GO:0030245", trait_key_policy=UNION) == ("GO:0030245",)

    for non_canonical in ("METPO:1000605", "METPO:1000612"):
        with pytest.raises(ValueError, match="not canonical"):
            policy_member_ids_for_trait(non_canonical, trait_key_policy=UNION)

    assert expected_trait_display_label("METPO:1000608", trait_key_policy=RAW) == (
        "facultatively aerobic"
    )
    assert expected_trait_display_label("METPO:1000608", trait_key_policy=UNION) == (
        "facultative (aerobe or anaerobe)"
    )
    assert (
        expected_trait_display_label(
            "METPO:1000608",
            "facultatively aerobic METPO:1000608",
            trait_key_policy=RAW,
        )
        == "facultatively aerobic METPO:1000608"
    )
    assert (
        expected_trait_display_label(
            "METPO:1000608",
            "facultatively aerobic METPO:1000608",
            trait_key_policy=UNION,
        )
        == "facultative (aerobe or anaerobe) METPO:1000608"
    )
    assert (
        expected_trait_display_label(
            "METPO:1000608",
            "Facultatively aerobic",
            trait_key_policy=UNION,
        )
        == "facultative (aerobe or anaerobe)"
    )
    assert (
        expected_trait_display_label(
            "METPO:1000608",
            "facultatively aerobic",
            trait_key_policy=UNION,
            member_ids=("METPO:1000608",),
        )
        == "facultatively aerobic"
    )
    assert (
        display_trait_label(
            "METPO:1000608",
            "facultatively aerobic",
            ("METPO:1000605", "METPO:1000608", "METPO:1000612"),
            trait_key_policy=RAW,
        )
        == "facultatively aerobic"
    )


def test_helpers_read_version_table_rather_than_unversioned_globals(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Mutating legacy module constants cannot alter lookup helpers; table entries drive them."""
    monkeypatch.setattr(policy_mod, "FACULTATIVE_TRAIT_ID", "METPO:9999999")
    monkeypatch.setattr(policy_mod, "FACULTATIVE_RAW_LABEL", "mutated raw")
    monkeypatch.setattr(policy_mod, "FACULTATIVE_UNION_LABEL", "mutated union")
    monkeypatch.setattr(policy_mod, "FACULTATIVE_MEMBER_IDS", ("METPO:9999999",))
    monkeypatch.setattr(policy_mod, "FACULTATIVE_EXTRA_MEMBER_IDS", frozenset())
    monkeypatch.setattr(policy_mod, "FACULTATIVE_SOURCE_VOCABULARIES", ())

    assert canonical_trait_object_id("METPO:1000605", trait_key_policy=UNION) == "METPO:1000608"
    assert policy_member_ids_for_trait("METPO:1000608", trait_key_policy=UNION) == (
        "METPO:1000605",
        "METPO:1000608",
        "METPO:1000612",
    )
    assert (
        expected_trait_display_label(
            "METPO:1000608", "facultatively aerobic", trait_key_policy=UNION
        )
        == "facultative (aerobe or anaerobe)"
    )
    _, members_json, vocab_json = policy_metadata("METPO:1000608", trait_key_policy=UNION)
    assert json.loads(members_json) == ["METPO:1000605", "METPO:1000608", "METPO:1000612"]
    assert len(json.loads(vocab_json)) == 5

    # Replace the version table with a synthetic definition to prove table lookup drives helpers.
    synthetic_vocab = (SourceVocabulary("METPO:1000777", "infores:synthetic", "synthetic_member"),)
    synthetic_union = TraitKeyPolicyDefinition(
        policy=UNION,
        canonical_trait_id="METPO:1000700",
        member_ids=("METPO:1000700", "METPO:1000777"),
        raw_label="synthetic raw trait",
        display_label="synthetic union (a or b)",
        source_vocabularies=synthetic_vocab,
    )
    monkeypatch.setattr(
        policy_mod,
        "TRAIT_KEY_POLICY_DEFINITIONS",
        MappingProxyType({RAW: TRAIT_KEY_POLICY_DEFINITIONS[RAW], UNION: synthetic_union}),
    )
    assert canonical_trait_object_id("METPO:1000777", trait_key_policy=UNION) == "METPO:1000700"
    assert canonical_trait_object_id("METPO:1000605", trait_key_policy=UNION) == "METPO:1000605"
    assert policy_member_ids_for_trait("METPO:1000700", trait_key_policy=UNION) == (
        "METPO:1000700",
        "METPO:1000777",
    )
    assert (
        expected_trait_display_label("METPO:1000700", "synthetic raw trait", trait_key_policy=UNION)
        == "synthetic union (a or b)"
    )
    assert (
        display_trait_label(
            "METPO:1000700",
            "synthetic raw trait",
            ("METPO:1000700", "METPO:1000777"),
            trait_key_policy=UNION,
        )
        == "synthetic union (a or b)"
    )
    _, syn_members, syn_vocabs = policy_metadata("METPO:1000700", trait_key_policy=UNION)
    assert json.loads(syn_members) == ["METPO:1000700", "METPO:1000777"]
    assert json.loads(syn_vocabs) == [
        {
            "object_id": "METPO:1000777",
            "source_id": "infores:synthetic",
            "source_label": "synthetic_member",
        }
    ]


def test_validate_trait_key_policy_definitions_rejects_in_place_v1_mutation() -> None:
    """In-place modification of a frozen v1 definition fails validation without a version bump."""
    validate_trait_key_policy_definitions(TRAIT_KEY_POLICY_DEFINITIONS)

    union_def = TRAIT_KEY_POLICY_DEFINITIONS[UNION]
    raw_def = TRAIT_KEY_POLICY_DEFINITIONS[RAW]

    mutated_cases: list[Mapping[str, TraitKeyPolicyDefinition]] = [
        {
            RAW: raw_def,
            UNION: replace(union_def, member_ids=("METPO:1000605", "METPO:1000608")),
        },
        {
            RAW: raw_def,
            UNION: replace(union_def, source_vocabularies=union_def.source_vocabularies[:-1]),
        },
        {
            RAW: raw_def,
            UNION: replace(union_def, display_label="edited display label"),
        },
        {
            RAW: replace(raw_def, member_ids=("METPO:1000605", "METPO:1000608")),
            UNION: union_def,
        },
        {
            RAW: raw_def,
        },
    ]
    for mutated in mutated_cases:
        raised = False
        message = ""
        try:
            validate_trait_key_policy_definitions(mutated)
        except ValueError as exc:
            raised = True
            message = str(exc).lower()
        assert raised
        assert any(token in message for token in ("frozen", "policy", "definition", "version"))

    # A new version key preserves v1 definitions while defining new semantics under v2.
    v2_def = TraitKeyPolicyDefinition(
        policy=cast(TraitKeyPolicy, "facultative_union_v2"),
        canonical_trait_id="METPO:1000608",
        member_ids=("METPO:1000605", "METPO:1000608"),
        raw_label="facultatively aerobic",
        display_label="facultative (v2 subset)",
        source_vocabularies=EXPECTED_V1_VOCABULARIES[:4],
    )
    extended = {
        RAW: raw_def,
        UNION: union_def,
        "facultative_union_v2": v2_def,
    }
    validate_trait_key_policy_definitions(extended, allow_additional_versions=True)
    with pytest.raises(ValueError, match="version"):
        validate_trait_key_policy_definitions(extended, allow_additional_versions=False)


def test_methods_trait_key_policy_documents_version_frozen_table() -> None:
    """methods/trait_key_policy.md documents TRAIT_KEY_POLICY_DEFINITIONS (#1137)."""
    content = (REPO_ROOT / "methods/trait_key_policy.md").read_text(encoding="utf-8")
    assert "TRAIT_KEY_POLICY_DEFINITIONS" in content
    assert "TraitKeyPolicyDefinition" in content
    assert "#1137" in content
    assert "new policy version" in content
