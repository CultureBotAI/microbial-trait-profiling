"""Versioned feature semantics preserve historical matrices and expose unions (#1129)."""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path
from typing import Literal

import pytest

from mtp.taxa_trait_matrix import (
    NodeMetadata,
    column_metadata_rows,
    create_taxa_trait_matrix,
    make_trait_key,
    read_kg_edges,
)
from mtp.workflow.contracts import read_table, validate_product
from mtp.workflow.coverage import write_feature_catalog
from mtp.workflow.matrix_contracts import validate_matrix_sources

Policy = Literal["raw_object_v1", "facultative_union_v1"]
RAW: Policy = "raw_object_v1"
UNION: Policy = "facultative_union_v1"
MEMBERS = ["METPO:1000605", "METPO:1000608", "METPO:1000612"]
CANONICAL = MEMBERS[1]
POLICY_COLUMNS = (
    "trait_key_policy",
    "trait_member_ids_json",
    "trait_source_vocabularies_json",
)
OLD_COLUMNS = (
    "column_id",
    "column_label",
    "predicate_id",
    "predicate_name",
    "object_id",
    "object_name",
    "object_category",
)


def _write(path: Path, header: list[str], rows: list[list[str]]) -> Path:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    return path


@pytest.fixture
def sources(tmp_path: Path) -> dict[str, Path]:
    edges = [
        ["t1", "NCBITaxon:1", "biolink:subclass_of", "NCBITaxon:10", "fixture:taxonomy"],
        ["t2", "NCBITaxon:2", "biolink:subclass_of", "NCBITaxon:10", "fixture:taxonomy"],
        ["t3", "NCBITaxon:3", "biolink:subclass_of", "NCBITaxon:10", "fixture:taxonomy"],
        ["a1", "NCBITaxon:1", "biolink:has_phenotype", MEMBERS[0], "infores:bacdive"],
        ["a2", "NCBITaxon:1", "biolink:has_phenotype", MEMBERS[1], "infores:madin_etal"],
        ["a3", "NCBITaxon:1", "biolink:has_phenotype", MEMBERS[2], "infores:bactotraits"],
        ["a4", "NCBITaxon:1", "biolink:has_phenotype", MEMBERS[0], "fixture:duplicate"],
        ["a5", "NCBITaxon:2", "biolink:has_phenotype", MEMBERS[2], "infores:bactotraits"],
        ["s2", "kgmicrobe.strain:two", "biolink:subclass_of", "NCBITaxon:2", "infores:bacdive"],
        ["s3", "kgmicrobe.strain:three", "biolink:subclass_of", "NCBITaxon:3", "infores:bacdive"],
        ["s4", "kgmicrobe.strain:two", "biolink:has_phenotype", MEMBERS[0], "infores:bacdive"],
        ["s5", "kgmicrobe.strain:three", "biolink:has_phenotype", MEMBERS[0], "infores:bacdive"],
    ]
    taxa = [
        ["increased", f"NCBITaxon:{i}", f"Named taxon {i}", "genus" if i == 10 else "species"]
        for i in (10, 1, 2, 3)
    ]
    nodes = [
        [member, "biolink:PhenotypicQuality", name]
        for member, name in zip(
            MEMBERS,
            ("facultatively anaerobic", "facultatively aerobic", "facultative oxygen preference"),
            strict=True,
        )
    ]
    return {
        "edges": _write(
            tmp_path / "edges.tsv",
            ["id", "subject", "predicate", "object", "primary_knowledge_source"],
            edges,
        ),
        "taxa": _write(
            tmp_path / "taxa.tsv", ["direction", "taxon_id", "taxon_name", "taxonomic_rank"], taxa
        ),
        "nodes": _write(tmp_path / "nodes.tsv", ["id", "category", "name"], nodes),
    }


def test_policy_is_explicit_and_cannot_leak_between_calls(sources: dict[str, Path]) -> None:
    for policy, expected in ((RAW, MEMBERS[0]), (UNION, CANONICAL), (RAW, MEMBERS[0])):
        assert (
            make_trait_key("biolink:has_phenotype", MEMBERS[0], trait_key_policy=policy) == expected
        )
        _, traits = read_kg_edges(sources["edges"], trait_key_policy=policy)
        assert traits["NCBITaxon:1"] == (set(MEMBERS) if policy == RAW else {CANONICAL})
    negated = "METPO:2000021"
    assert make_trait_key(negated, MEMBERS[0], trait_key_policy=RAW) == f"{negated}|{MEMBERS[0]}"
    assert make_trait_key(negated, MEMBERS[0], trait_key_policy=UNION) == f"{negated}|{CANONICAL}"
    assert (
        make_trait_key("METPO:2000001", "GO:0030245", trait_key_policy=UNION)
        == "METPO:2000001|GO:0030245"
    )


@pytest.mark.parametrize("policy", [RAW, UNION])
@pytest.mark.parametrize("scope", ["off", "phenotype"])
@pytest.mark.parametrize("modern", [False, True])
def test_matrix_and_independent_source_check_use_declared_policy(
    sources: dict[str, Path],
    tmp_path: Path,
    policy: Policy,
    scope: Literal["off", "phenotype"],
    modern: bool,
) -> None:
    matrix, metadata = tmp_path / "matrix.tsv", tmp_path / "metadata.tsv"
    create_taxa_trait_matrix(
        taxa_path=sources["taxa"],
        edges_path=sources["edges"],
        nodes_path=sources["nodes"],
        output_path=matrix,
        column_metadata_path=metadata,
        strain_rollup=scope,
        trait_key_policy=policy,
        include_policy_metadata=modern,
    )
    _, rows = read_table(matrix)
    genus = next(row for row in rows if row["taxon_id"] == "NCBITaxon:10")
    if policy == RAW:
        assert [int(genus[member]) for member in MEMBERS] == (
            [1, 1, 2] if scope == "off" else [3, 1, 2]
        )
    else:
        assert int(genus[CANONICAL]) == (2 if scope == "off" else 3)
        assert MEMBERS[0] not in genus
        assert MEMBERS[2] not in genus
    header, _ = read_table(metadata)
    assert header == ((*OLD_COLUMNS, *POLICY_COLUMNS) if modern else OLD_COLUMNS)
    validate_matrix_sources(
        matrix,
        sources["taxa"],
        sources["edges"],
        sources["nodes"],
        metadata,
        strain_rollup=scope,
        trait_key_policy=policy,
        include_policy_metadata=modern,
    )
    other = UNION if policy == RAW else RAW
    with pytest.raises(ValueError, match=r"(vocabulary|metadata|policy|count)"):
        validate_matrix_sources(
            matrix,
            sources["taxa"],
            sources["edges"],
            sources["nodes"],
            metadata,
            strain_rollup=scope,
            trait_key_policy=other,
            include_policy_metadata=modern,
        )


def test_union_metadata_discloses_declared_members_and_vocabularies_without_relabelling() -> None:
    nodes = {CANONICAL: NodeMetadata("facultatively aerobic", "biolink:PhenotypicQuality")}
    header, rows = column_metadata_rows({"taxon": Counter({CANONICAL: 1})}, nodes)
    row = dict(zip(header, rows[0], strict=True))
    assert row["object_name"] == "facultatively aerobic"
    assert row["column_label"] == f"facultatively aerobic {CANONICAL}"
    assert row["trait_key_policy"] == UNION
    assert json.loads(row["trait_member_ids_json"]) == MEMBERS
    mappings = json.loads(row["trait_source_vocabularies_json"])
    assert {
        tuple(mapping[key] for key in ("object_id", "source_id", "source_label"))
        for mapping in mappings
    } == {
        (MEMBERS[0], "infores:madin_etal", "facultative"),
        (MEMBERS[0], "infores:bacdive", "facultative anaerobe"),
        (MEMBERS[1], "infores:madin_etal", "facultative"),
        (MEMBERS[1], "infores:bacdive", "facultative aerobe"),
        (MEMBERS[2], "infores:bactotraits", "Ox_facultative_aerobe_anaerobe"),
    }
    raw_header, raw_rows = column_metadata_rows(
        {"taxon": Counter({CANONICAL: 1})},
        nodes,
        trait_key_policy=RAW,
    )
    raw = dict(zip(raw_header, raw_rows[0], strict=True))
    assert json.loads(raw["trait_member_ids_json"]) == [CANONICAL]
    assert len(json.loads(raw["trait_source_vocabularies_json"])) == 2


@pytest.mark.parametrize("modern", [False, True])
def test_catalog_preserves_matching_metadata_format(
    sources: dict[str, Path], tmp_path: Path, modern: bool
) -> None:
    metadata, catalog = tmp_path / "metadata.tsv", tmp_path / "catalog.tsv"
    create_taxa_trait_matrix(
        taxa_path=sources["taxa"],
        edges_path=sources["edges"],
        nodes_path=sources["nodes"],
        output_path=tmp_path / "matrix.tsv",
        column_metadata_path=metadata,
        trait_key_policy=RAW,
        include_policy_metadata=modern,
    )
    write_feature_catalog(metadata, catalog, trait_key_policy=RAW, include_policy_metadata=modern)
    header, rows = read_table(catalog)
    expected = (*OLD_COLUMNS, *POLICY_COLUMNS) if modern else OLD_COLUMNS
    assert header == (*expected, "feature_class", "annotation_evidence")
    facts = validate_product(
        catalog,
        "feature_catalog",
        context={"feature_metadata": metadata, "kg_nodes": sources["nodes"]},
        trait_key_policy=RAW,
        include_policy_metadata=modern,
    )
    assert facts["trait_key_policy"] == RAW
    assert rows[0]["object_name"]


@pytest.mark.parametrize(
    "damage", ["policy", "members", "vocabulary", "partial", "legacy_downgrade"]
)
def test_new_metadata_cannot_drop_or_forge_policy_disclosure(
    sources: dict[str, Path],
    tmp_path: Path,
    damage: str,
) -> None:
    metadata = tmp_path / "metadata.tsv"
    create_taxa_trait_matrix(
        taxa_path=sources["taxa"],
        edges_path=sources["edges"],
        nodes_path=sources["nodes"],
        output_path=tmp_path / "matrix.tsv",
        column_metadata_path=metadata,
    )
    header, rows = read_table(metadata)
    edited = list(header)
    if damage == "policy":
        rows[0]["trait_key_policy"] = RAW
    elif damage == "members":
        rows[0]["trait_member_ids_json"] = json.dumps([CANONICAL])
    elif damage == "vocabulary":
        rows[0]["trait_source_vocabularies_json"] = "[]"
    elif damage == "partial":
        edited.remove("trait_member_ids_json")
    else:
        edited = list(OLD_COLUMNS)
    _write(metadata, edited, [[row[column] for column in edited] for row in rows])
    with pytest.raises(ValueError, match=r"(metadata|policy|columns|vocabular|member)"):
        validate_product(metadata, "feature_metadata", trait_key_policy=UNION)


def test_unknown_policy_is_rejected() -> None:
    from mtp.trait_key_policy import parse_trait_key_policy

    for value in ("", "facultative_union_v2", None, 4):
        with pytest.raises(ValueError, match="trait_key_policy"):
            parse_trait_key_policy(value)
