from pathlib import Path

from bsdbng.taxa_trait_matrix import create_taxa_trait_matrix


def test_create_taxa_trait_matrix_counts_self_and_descendants(tmp_path: Path) -> None:
    taxa_path = tmp_path / "taxa.tsv"
    nodes_path = tmp_path / "nodes.tsv"
    edges_path = tmp_path / "edges.tsv"
    output_path = tmp_path / "taxa_trait_matrix.tsv"
    column_metadata_path = tmp_path / "taxa_trait_matrix_columns.tsv"

    taxa_path.write_text(
        "\t".join(
            [
                "study_id",
                "experiment_id",
                "signature_id",
                "direction",
                "taxon_id",
                "taxon_name",
                "taxonomic_rank",
            ]
        )
        + "\n"
        + "\n".join(
            [
                "bsdb:1\tbsdb:1-1\tbsdb:1-1-1\tincreased\tNCBITaxon:900\tFamilyaceae\tfamily",
                "bsdb:1\tbsdb:1-1\tbsdb:1-1-1\tdecreased\tNCBITaxon:100\tGenus\tgenus",
                "bsdb:1\tbsdb:1-1\tbsdb:1-1-1\tincreased\tNCBITaxon:101\tSpecies a\tspecies",
                "bsdb:1\tbsdb:1-1\tbsdb:1-1-1\tunchanged\tNCBITaxon:102\tIgnored\tspecies",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    nodes_path.write_text(
        "\n".join(
            [
                "id\tcategory\tname",
                "METPO:1000615\tbiolink:PhenotypicQuality\tmesophilic",
                "METPO:1000603\tbiolink:PhenotypicQuality\tanaerobic",
                "METPO:1000636\tbiolink:PhenotypicQuality\tchemoheterotrophic",
                "METPO:2000006\tbiolink:OntologyClass\tuses as carbon source",
                "METPO:2000103\tbiolink:OntologyClass\tcapable of",
                "CHEBI:123\tbiolink:SmallMolecule\tacetate",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    edges_path.write_text(
        "\n".join(
            [
                "id\tsubject\tpredicate\tobject",
                "e1\tNCBITaxon:100\tbiolink:subclass_of\tNCBITaxon:900",
                "e2\tNCBITaxon:101\tbiolink:subclass_of\tNCBITaxon:100",
                "e3\tNCBITaxon:102\tbiolink:subclass_of\tNCBITaxon:100",
                "e4\tNCBITaxon:100\tbiolink:has_phenotype\tMETPO:1000615",
                "e5\tNCBITaxon:101\tbiolink:has_phenotype\tMETPO:1000603",
                "e6\tNCBITaxon:101\tbiolink:has_phenotype\tMETPO:1000603",
                "e7\tNCBITaxon:102\tbiolink:has_phenotype\tMETPO:1000615",
                "e8\tNCBITaxon:101\tMETPO:2000103\tMETPO:1000636",
                "e9\tNCBITaxon:101\tMETPO:2000103\tGO:0008150",
                "e10\tNCBITaxon:101\tMETPO:2000006\tCHEBI:123",
                "e11\tNCBITaxon:101\tMETPO:2000006\tCHEBI:123",
                "e12\tNCBITaxon:101\tbiolink:location_of\tNCBITaxon:999",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    stats = create_taxa_trait_matrix(
        taxa_path=taxa_path,
        kg_dir=tmp_path,
        output_path=output_path,
        column_metadata_path=column_metadata_path,
        edges_path=edges_path,
        nodes_path=nodes_path,
    )

    assert stats.reported_taxa == 3
    assert stats.traits == 5
    assert output_path.read_text(encoding="utf-8").splitlines() == [
        (
            "taxon_id\ttaxon_name\ttaxonomic_rank\tMETPO:1000603"
            "\tMETPO:2000103|GO:0008150\tMETPO:1000636\tMETPO:1000615"
            "\tMETPO:2000006|CHEBI:123"
        ),
        "NCBITaxon:900\tFamilyaceae\tfamily\t1\t1\t1\t2\t1",
        "NCBITaxon:100\tGenus\tgenus\t1\t1\t1\t2\t1",
        "NCBITaxon:101\tSpecies a\tspecies\t1\t1\t1\t0\t1",
    ]
    assert column_metadata_path.read_text(encoding="utf-8").splitlines() == [
        (
            "column_id\tcolumn_label\tpredicate_id\tpredicate_name\tobject_id"
            "\tobject_name\tobject_category"
        ),
        (
            "METPO:1000603\tanaerobic METPO:1000603\t\t\tMETPO:1000603"
            "\tanaerobic\tbiolink:PhenotypicQuality"
        ),
        (
            "METPO:2000103|GO:0008150\tcapable of: GO:0008150 [METPO:2000103]"
            "\tMETPO:2000103\tcapable of\tGO:0008150\t\t"
        ),
        (
            "METPO:1000636\tchemoheterotrophic METPO:1000636\t\t\tMETPO:1000636"
            "\tchemoheterotrophic\tbiolink:PhenotypicQuality"
        ),
        (
            "METPO:1000615\tmesophilic METPO:1000615\t\t\tMETPO:1000615"
            "\tmesophilic\tbiolink:PhenotypicQuality"
        ),
        (
            "METPO:2000006|CHEBI:123\tuses as carbon source: acetate CHEBI:123 "
            "[METPO:2000006]\tMETPO:2000006\tuses as carbon source\tCHEBI:123"
            "\tacetate\tbiolink:SmallMolecule"
        ),
    ]
