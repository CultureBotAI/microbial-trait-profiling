"""Exercise standalone audits on independently counted raw/union fixtures (#1136)."""

from __future__ import annotations

import csv
import json
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Literal, cast

import pytest

from mtp import coassertion_audit, label_support_audit
from mtp.trait_key_policy import TraitKeyPolicy

FAC = "METPO:1000608"
BACTOTRAITS_FAC = "METPO:1000612"
BACDIVE_FAC = "METPO:1000605"
OBLIGATE = "METPO:1000606"
COCCUS = "METPO:1000668"
PAIR = "facultatively_anaerobic_and_obligately_aerobic"
Scope = Literal["off", "phenotype"]


def _write_table(path: Path, fields: list[str], rows: Sequence[Mapping[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _read_table(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


@pytest.fixture
def kg(tmp_path: Path) -> tuple[Path, Path]:
    """Three carriers differ by original object; the strain is linked to only one."""
    edges, nodes = tmp_path / "edges.tsv", tmp_path / "nodes.tsv"
    triples = [
        ("NCBITaxon:11", "biolink:subclass_of", "NCBITaxon:10", "infores:ncbi"),
        ("NCBITaxon:12", "biolink:subclass_of", "NCBITaxon:10", "infores:ncbi"),
        ("NCBITaxon:13", "biolink:subclass_of", "NCBITaxon:10", "infores:ncbi"),
        ("NCBITaxon:11", "biolink:has_phenotype", FAC, "infores:direct"),
        ("NCBITaxon:11", "biolink:has_phenotype", OBLIGATE, "infores:direct"),
        ("NCBITaxon:12", "biolink:has_phenotype", BACTOTRAITS_FAC, "infores:bactotraits"),
        ("NCBITaxon:12", "biolink:has_phenotype", OBLIGATE, "infores:direct"),
        ("NCBITaxon:13", "biolink:has_phenotype", OBLIGATE, "infores:direct"),
        ("kgmicrobe.strain:three", "biolink:subclass_of", "NCBITaxon:13", "infores:bacdive"),
        ("kgmicrobe.strain:three", "biolink:has_phenotype", BACDIVE_FAC, "infores:bacdive"),
        ("NCBITaxon:20", "biolink:has_phenotype", COCCUS, "infores:direct"),
    ]
    _write_table(
        edges,
        ["id", "subject", "predicate", "object", "primary_knowledge_source"],
        [
            {
                "id": f"edge-{index}",
                "subject": subject,
                "predicate": predicate,
                "object": obj,
                "primary_knowledge_source": source,
            }
            for index, (subject, predicate, obj, source) in enumerate(triples)
        ],
    )
    classes = [
        "METPO:1000602",
        "METPO:1000603",
        "METPO:1000604",
        BACDIVE_FAC,
        OBLIGATE,
        "METPO:1000607",
        FAC,
        "METPO:1000611",
        BACTOTRAITS_FAC,
        "METPO:1000614",
        "METPO:1000615",
        COCCUS,
    ]
    _write_table(
        nodes,
        ["id", "name"],
        [{"id": key, "name": f"Fixture class {key}"} for key in classes]
        + [
            {"id": f"NCBITaxon:{number}", "name": f"Fixture taxon {number}"}
            for number in (10, 11, 12, 13, 20)
        ],
    )
    return edges, nodes


def _label_inputs(tmp_path: Path, policy: TraitKeyPolicy, scope: Scope) -> dict[str, Path]:
    """Hand-counted matrix/profile products, independent of the matrix implementation."""
    paths = {name: tmp_path / f"{name}.tsv" for name in ("matrix", "profiles", "taxa")}
    raw = policy == "raw_object_v1"
    features = [FAC, OBLIGATE, COCCUS]
    if raw:
        features += [BACTOTRAITS_FAC, BACDIVE_FAC]
    genus = {FAC: 1 if raw else 2 if scope == "off" else 3, OBLIGATE: 3, COCCUS: 0}
    if raw:
        genus.update({BACTOTRAITS_FAC: 1, BACDIVE_FAC: 0 if scope == "off" else 1})
    reference = dict.fromkeys(features, 0)
    reference[COCCUS] = 1
    _write_table(
        paths["matrix"],
        ["taxon_id", "taxon_name", "taxonomic_rank", *features],
        [
            {
                "taxon_id": "NCBITaxon:10",
                "taxon_name": "Fixture genus",
                "taxonomic_rank": "genus",
                **genus,
            },
            {
                "taxon_id": "NCBITaxon:20",
                "taxon_name": "Fixture reference",
                "taxonomic_rank": "species",
                **reference,
            },
        ],
    )
    _write_table(
        paths["profiles"],
        ["study_id", "experiment_id", "inc_dec", *features],
        [
            {"study_id": "s1", "experiment_id": "e1", "inc_dec": "inc", **genus},
            {"study_id": "s1", "experiment_id": "e1", "inc_dec": "dec", **reference},
        ],
    )
    _write_table(
        paths["taxa"],
        ["study_id", "experiment_id", "direction", "taxon_id"],
        [
            {
                "study_id": "s1",
                "experiment_id": "e1",
                "direction": "increased",
                "taxon_id": "NCBITaxon:10",
            },
            {
                "study_id": "s1",
                "experiment_id": "e1",
                "direction": "decreased",
                "taxon_id": "NCBITaxon:20",
            },
        ],
    )
    return paths


@pytest.mark.parametrize("policy", ["raw_object_v1", "facultative_union_v1"])
def test_coassertion_public_audit_keeps_policy_specific_carriers(
    kg: tuple[Path, Path], tmp_path: Path, policy: TraitKeyPolicy
) -> None:
    output = tmp_path / "coassertion"
    record = json.loads(
        coassertion_audit.run_audit(*kg, output, trait_key_policy=policy).read_bytes()
    )
    summary = {
        row["scope"]: int(row["coasserted_taxon_nodes"])
        for row in _read_table(output / "coassertion_summary.tsv")
        if row["pair_id"] == PAIR
    }
    assert summary == (
        {"off": 1, "phenotype": 1, "full": 1}
        if policy == "raw_object_v1"
        else {"off": 2, "phenotype": 3, "full": 3}
    )
    evidence = [
        row
        for row in _read_table(output / "coassertion_evidence.tsv")
        if row["pair_id"] == PAIR and row["side"] == "left"
    ]
    assert {row["object"] for row in evidence} == (
        {FAC} if policy == "raw_object_v1" else {FAC, BACTOTRAITS_FAC, BACDIVE_FAC}
    )
    assert {row["trait_key"] for row in evidence} == {FAC}
    assert record["trait_key_policy"] == policy
    assert record["analysis_kind"] == "kg_annotation_census_not_a_core_or_arm_rerun"
    assert policy in (output / "README.md").read_text(encoding="utf-8")


@pytest.mark.parametrize("policy", ["raw_object_v1", "facultative_union_v1"])
@pytest.mark.parametrize("scope", ["off", "phenotype"])
def test_label_support_public_audit_reconciles_policy_specific_carriers(
    kg: tuple[Path, Path], tmp_path: Path, policy: TraitKeyPolicy, scope: Scope
) -> None:
    paths = _label_inputs(tmp_path, policy, scope)
    output = tmp_path / "support"
    record = label_support_audit.run_label_support_audit(
        taxa_path=paths["taxa"],
        matrix_path=paths["matrix"],
        profiles_path=paths["profiles"],
        edges_path=kg[0],
        nodes_path=kg[1],
        output_dir=output,
        strain_rollup=scope,
        trait_key_policy=policy,
    )
    support = next(
        row
        for row in _read_table(output / "taxon_trait_label_support.tsv")
        if row["taxon_id"] == "NCBITaxon:10" and row["trait_key"] == FAC
    )
    expected = 1 if policy == "raw_object_v1" else 2 if scope == "off" else 3
    assert int(support["carrier_count"]) == expected
    evidence = {
        row["carrier_taxon_id"]: row["primary_sources"]
        for row in _read_table(output / "carrier_source_evidence.tsv")
        if row["reported_taxon_id"] == "NCBITaxon:10" and row["trait_key"] == FAC
    }
    expected_evidence = {"NCBITaxon:11": "infores:direct"}
    if policy == "facultative_union_v1":
        expected_evidence["NCBITaxon:12"] = "infores:bactotraits"
        if scope == "phenotype":
            expected_evidence["NCBITaxon:13"] = "infores:bacdive"
    assert evidence == expected_evidence
    assert record["trait_key_policy"] == policy
    assert record["schema_version"] == "1.0"
    assert (
        record["verification_status"] == "hashes_and_counts_reconciled_not_core_execution_verified"
    )
    assert json.loads((output / "provenance.json").read_bytes())["trait_key_policy"] == policy


@pytest.mark.parametrize("policy", ["raw_object_v1", "facultative_union_v1"])
def test_coassertion_cli_forwards_policy(
    kg: tuple[Path, Path], tmp_path: Path, policy: TraitKeyPolicy
) -> None:
    output = tmp_path / "cli-coassertion"
    assert (
        coassertion_audit.main(
            [
                "--edges",
                str(kg[0]),
                "--nodes",
                str(kg[1]),
                "--output-dir",
                str(output),
                "--trait-key-policy",
                policy,
            ]
        )
        == 0
    )
    rows = _read_table(output / "coassertion_summary.tsv")
    row = next(row for row in rows if row["scope"] == "phenotype" and row["pair_id"] == PAIR)
    assert int(row["coasserted_taxon_nodes"]) == (1 if policy == "raw_object_v1" else 3)
    assert json.loads((output / "provenance.json").read_bytes())["trait_key_policy"] == policy


@pytest.mark.parametrize("policy", ["raw_object_v1", "facultative_union_v1"])
def test_label_support_cli_forwards_policy(
    kg: tuple[Path, Path], tmp_path: Path, policy: TraitKeyPolicy
) -> None:
    paths = _label_inputs(tmp_path, policy, "phenotype")
    output = tmp_path / "cli-support"
    assert (
        label_support_audit.main(
            [
                "--taxa",
                str(paths["taxa"]),
                "--matrix",
                str(paths["matrix"]),
                "--profiles",
                str(paths["profiles"]),
                "--edges",
                str(kg[0]),
                "--nodes",
                str(kg[1]),
                "--strain-rollup",
                "phenotype",
                "--output-dir",
                str(output),
                "--allow-temporary-output",
                "--trait-key-policy",
                policy,
            ]
        )
        == 0
    )
    rows = _read_table(output / "taxon_trait_label_support.tsv")
    row = next(row for row in rows if row["taxon_id"] == "NCBITaxon:10" and row["trait_key"] == FAC)
    assert int(row["carrier_count"]) == (1 if policy == "raw_object_v1" else 3)
    assert json.loads((output / "provenance.json").read_bytes())["trait_key_policy"] == policy


@pytest.mark.parametrize("scope", ["off", "phenotype"])
@pytest.mark.parametrize("matrix_policy", ["raw_object_v1", "facultative_union_v1"])
def test_label_support_refuses_a_matrix_from_the_other_policy(
    kg: tuple[Path, Path], tmp_path: Path, scope: Scope, matrix_policy: TraitKeyPolicy
) -> None:
    paths = _label_inputs(tmp_path, matrix_policy, scope)
    other_policy: TraitKeyPolicy = (
        "facultative_union_v1" if matrix_policy == "raw_object_v1" else "raw_object_v1"
    )
    output = tmp_path / "mismatched"
    with pytest.raises(ValueError, match="KG carrier count disagrees with matrix"):
        label_support_audit.run_label_support_audit(
            taxa_path=paths["taxa"],
            matrix_path=paths["matrix"],
            profiles_path=paths["profiles"],
            edges_path=kg[0],
            nodes_path=kg[1],
            output_dir=output,
            strain_rollup=scope,
            trait_key_policy=other_policy,
        )
    assert not output.exists()


def test_invalid_policy_refused_before_input_reads_or_output_creation(tmp_path: Path) -> None:
    missing, output = tmp_path / "missing.tsv", tmp_path / "invalid-output"
    invalid = cast(TraitKeyPolicy, "unrecognised_policy")
    with pytest.raises(ValueError, match="trait_key_policy"):
        coassertion_audit.run_audit(missing, missing, output, trait_key_policy=invalid)
    assert not output.exists()
    with pytest.raises(ValueError, match="trait_key_policy"):
        label_support_audit.run_label_support_audit(
            taxa_path=missing,
            matrix_path=missing,
            profiles_path=missing,
            edges_path=missing,
            nodes_path=missing,
            output_dir=output,
            strain_rollup="off",
            trait_key_policy=invalid,
        )
    assert not output.exists()
