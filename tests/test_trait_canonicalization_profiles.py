"""Asymmetric canonical trait carriers through matrix and profile operations (#1129)."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pytest
import yaml

from mtp.taxa_trait_matrix import StrainRollup, create_taxa_trait_matrix
from mtp.trait_key_policy import TraitKeyPolicy
from mtp.workflow.matrix_contracts import validate_matrix_sources
from mtp.workflow.operations import (
    count_profiles_operation,
    ora_profiles_operation,
    taxa_matrix_operation,
)

CANONICAL = "METPO:1000608"
OTHER = "METPO:1000603"
POLICY: TraitKeyPolicy = "facultative_union_v1"


@dataclass(frozen=True)
class Inputs:
    reported: Path
    edges: Path
    nodes: Path
    studies: Path

    def matrix_inputs(self) -> dict[str, Path]:
        return {"reported_taxa": self.reported, "kg_edges": self.edges, "kg_nodes": self.nodes}


def _table(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


@pytest.fixture
def inputs(tmp_path: Path) -> Inputs:
    """Three sole direct carriers, three sole strain carriers, and one shared carrier."""
    reported, edges, nodes = (tmp_path / name for name in ("taxa.tsv", "edges.tsv", "nodes.tsv"))
    reported.write_text(
        "direction\ttaxon_id\ttaxon_name\ttaxonomic_rank\n"
        "increased\tNCBITaxon:900\tFamily\tfamily\n"
        "increased\tNCBITaxon:100\tGenus\tgenus\n"
        + "".join(
            f"increased\tNCBITaxon:{taxon}\tSpecies {taxon}\tspecies\n" for taxon in range(101, 108)
        ),
        encoding="utf-8",
    )
    nodes.write_text(
        "id\tcategory\tname\n"
        "METPO:1000603\tbiolink:PhenotypicQuality\tanaerobic\n"
        "METPO:1000605\tbiolink:PhenotypicQuality\tfacultatively anaerobic\n"
        "METPO:1000608\tbiolink:PhenotypicQuality\tfacultatively aerobic\n"
        "METPO:1000612\tbiolink:PhenotypicQuality\tfacultative\n",
        encoding="utf-8",
    )
    records = [("NCBITaxon:100", "biolink:subclass_of", "NCBITaxon:900")]
    records.extend(
        (f"NCBITaxon:{taxon}", "biolink:subclass_of", "NCBITaxon:100") for taxon in range(101, 108)
    )
    for taxon, term in zip(
        range(101, 104), ("METPO:1000605", CANONICAL, "METPO:1000612"), strict=True
    ):
        records.append((f"NCBITaxon:{taxon}", "biolink:has_phenotype", term))
    for taxon, term in zip(
        range(104, 107), ("METPO:1000605", CANONICAL, "METPO:1000612"), strict=True
    ):
        records.extend(
            [
                (f"kgmicrobe.strain:{taxon}", "biolink:subclass_of", f"NCBITaxon:{taxon}"),
                (f"kgmicrobe.strain:{taxon}", "biolink:has_phenotype", term),
            ]
        )
    records.append(("NCBITaxon:101", "biolink:has_phenotype", OTHER))
    for subject in ("NCBITaxon:107", "kgmicrobe.strain:107a", "kgmicrobe.strain:107b"):
        if subject.startswith("kgmicrobe.strain:"):
            records.append((subject, "biolink:subclass_of", "NCBITaxon:107"))
        records.extend(
            (subject, "biolink:has_phenotype", term)
            for term in ("METPO:1000605", CANONICAL, "METPO:1000612")
        )
    edges.write_text(
        "id\tsubject\tpredicate\tobject\n"
        + "".join(
            f"e{i}\t{subject}\t{predicate}\t{obj}\n"
            for i, (subject, predicate, obj) in enumerate(records)
        ),
        encoding="utf-8",
    )
    experiments: list[dict[str, Any]] = []
    for number, increased, decreased in (
        (1, [100, 100], [104]),
        (2, [103], [106]),
        (3, [107, 107, 101], [102, 105]),
    ):
        experiment = f"bsdb:1129-{number}"
        experiments.append(
            {
                "id": experiment,
                "host_species": "Homo sapiens",
                "sequencing_type": "16S",
                "signatures": [
                    {
                        "id": f"{experiment}-{direction}",
                        "direction": direction,
                        "taxa": [
                            {"id": f"NCBITaxon:{taxon}", "taxon_name": f"Taxon {taxon}"}
                            for taxon in taxa
                        ],
                    }
                    for direction, taxa in (("increased", increased), ("decreased", decreased))
                ],
            }
        )
    studies = tmp_path / "studies"
    studies.mkdir()
    (studies / "bsdb_1129.yaml").write_text(
        yaml.safe_dump({"id": "bsdb:1129", "experiments": experiments}), encoding="utf-8"
    )
    return Inputs(reported, edges, nodes, studies)


@pytest.mark.parametrize("scope", ["off", "phenotype", "full"])
def test_native_matrix_preserves_each_sole_source_and_deduplicates_overlap(
    inputs: Inputs,
    tmp_path: Path,
    scope: StrainRollup,
) -> None:
    matrix, metadata = tmp_path / "matrix.tsv", tmp_path / "metadata.tsv"
    stats = create_taxa_trait_matrix(
        taxa_path=inputs.reported,
        edges_path=inputs.edges,
        nodes_path=inputs.nodes,
        output_path=matrix,
        column_metadata_path=metadata,
        strain_rollup=scope,
        trait_key_policy=POLICY,
    )
    rows = {row["taxon_id"]: row for row in _table(matrix)}
    assert stats.traits == 2
    assert set(next(iter(rows.values()))) == {
        "taxon_id",
        "taxon_name",
        "taxonomic_rank",
        CANONICAL,
        OTHER,
    }
    assert {taxon: int(row[CANONICAL]) for taxon, row in rows.items()} == {
        "NCBITaxon:900": 4 if scope == "off" else 7,
        "NCBITaxon:100": 4 if scope == "off" else 7,
        **{f"NCBITaxon:{taxon}": 1 for taxon in (101, 102, 103, 107)},
        **{f"NCBITaxon:{taxon}": int(scope != "off") for taxon in (104, 105, 106)},
    }
    assert {taxon: int(row[OTHER]) for taxon, row in rows.items()} == {
        f"NCBITaxon:{taxon}": int(taxon in {900, 100, 101})
        for taxon in (900, 100, *range(101, 108))
    }
    assert {row["column_id"] for row in _table(metadata)} == {CANONICAL, OTHER}
    assert {row["object_id"] for row in _table(metadata)} == {CANONICAL, OTHER}
    facts = validate_matrix_sources(
        matrix,
        inputs.reported,
        inputs.edges,
        inputs.nodes,
        metadata,
        strain_rollup=scope,
        trait_key_policy=POLICY,
    )
    assert facts["kg_source_aggregation_checked"] is True
    assert facts["kg_source_cells_checked"] == 18


@pytest.mark.parametrize("scope", ["off", "phenotype", "full"])
def test_core_count_and_ora_operations_preserve_asymmetric_canonical_carriers(
    inputs: Inputs,
    tmp_path: Path,
    scope: StrainRollup,
) -> None:
    matrix = {"taxa_matrix": tmp_path / "matrix.tsv", "feature_metadata": tmp_path / "metadata.tsv"}
    taxa_matrix_operation(
        inputs.matrix_inputs(),
        matrix,
        {"strain_rollup": scope, "trait_key_policy": POLICY},
    )
    validate_matrix_sources(
        matrix["taxa_matrix"],
        inputs.reported,
        inputs.edges,
        inputs.nodes,
        matrix["feature_metadata"],
        strain_rollup=scope,
        trait_key_policy=POLICY,
    )
    profile_inputs = {"studies": inputs.studies, "taxa_matrix": matrix["taxa_matrix"]}
    count_outputs = {
        "count_profiles": tmp_path / "counts.tsv",
        "count_coverage": tmp_path / "counts_coverage.tsv",
    }
    ora_outputs = {
        "ora_profiles": tmp_path / "ora.tsv",
        "ora_coverage": tmp_path / "ora_coverage.tsv",
    }
    count_profiles_operation(profile_inputs, count_outputs, {"count_host_species": "Homo sapiens"})
    ora_profiles_operation(profile_inputs, ora_outputs, {"ora_host_species": "Homo sapiens"})
    counts = {
        (row["experiment_id"], row["inc_dec"]): row
        for row in _table(count_outputs["count_profiles"])
    }
    ora = {
        (row["experiment_id"], row["inc_dec"]): row for row in _table(ora_outputs["ora_profiles"])
    }
    active = int(scope != "off")
    expected_counts = {
        ("bsdb:1129-1", "inc"): 4 if scope == "off" else 7,
        ("bsdb:1129-1", "dec"): active,
        ("bsdb:1129-2", "inc"): 1,
        ("bsdb:1129-2", "dec"): active,
        ("bsdb:1129-3", "inc"): 2,
        ("bsdb:1129-3", "dec"): 1 + active,
    }
    assert set(counts) == set(ora) == set(expected_counts)
    assert {key: int(row[CANONICAL]) for key, row in counts.items()} == expected_counts
    assert {key: int(row[CANONICAL]) for key, row in ora.items()} == {
        **expected_counts,
        ("bsdb:1129-1", "inc"): 1,
    }
    assert {key: int(row["taxa_count"]) for key, row in ora.items()} == {
        key: 2 if key[0] == "bsdb:1129-3" else 1 for key in expected_counts
    }
    for profiles in (counts, ora):
        assert {key: int(row[OTHER]) for key, row in profiles.items()} == {
            key: int(key in {("bsdb:1129-1", "inc"), ("bsdb:1129-3", "inc")})
            for key in expected_counts
        }
        assert all(
            "METPO:1000605" not in row and "METPO:1000612" not in row for row in profiles.values()
        )
