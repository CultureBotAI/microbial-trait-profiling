from pathlib import Path

import pytest

from bsdbng.taxa_trait_matrix import create_taxa_trait_matrix
from mtp import taxa_trait_matrix


def test_create_taxa_trait_matrix_counts_self_and_descendants(tmp_path: Path) -> None:
    taxa_path = tmp_path / "bugsigdb_taxa.tsv"
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


TAXONOMY_AND_TAXON_EDGES = [
    "e1\tNCBITaxon:100\tbiolink:subclass_of\tNCBITaxon:900",
    "e2\tNCBITaxon:101\tbiolink:subclass_of\tNCBITaxon:100",
    "e3\tNCBITaxon:102\tbiolink:subclass_of\tNCBITaxon:100",
    "e4\tNCBITaxon:101\tbiolink:has_phenotype\tMETPO:1000603",
    "e5\tNCBITaxon:102\tbiolink:has_phenotype\tMETPO:1000605",
]
# BacDive strain records (#774). Strains 1 and 2 belong to species 101, which itself carries
# METPO:1000603; strain 3 belongs to species 102; strain 4 has no taxon link and strain 5
# links only to a non-NCBITaxon genus node.
STRAIN_EDGES = [
    "s1\tkgmicrobe.strain:1\tbiolink:subclass_of\tNCBITaxon:101",
    "s2\tkgmicrobe.strain:1\tbiolink:has_phenotype\tMETPO:1000605",
    "s3\tkgmicrobe.strain:1\tbiolink:has_phenotype\tMETPO:1000603",
    "s4\tkgmicrobe.strain:2\tbiolink:subclass_of\tNCBITaxon:101",
    "s5\tkgmicrobe.strain:2\tbiolink:has_phenotype\tMETPO:1000605",
    "s6\tkgmicrobe.strain:3\tbiolink:subclass_of\tNCBITaxon:102",
    "s7\tkgmicrobe.strain:3\tMETPO:2000103\tGO:0030245",
    "s8\tkgmicrobe.strain:3\tbiolink:associated_with_sensitivity_to\tCHEBI:123",
    "s9\tkgmicrobe.strain:4\tbiolink:has_phenotype\tMETPO:1000608",
    "s10\tkgmicrobe.strain:5\tbiolink:subclass_of\tkgmicrobe.genus:Other",
    "s11\tkgmicrobe.strain:5\tbiolink:has_phenotype\tMETPO:1000668",
]


def _write_strain_fixture(directory: Path, edges: list[str]) -> tuple[Path, Path, Path]:
    taxa_path = directory / "bugsigdb_taxa.tsv"
    nodes_path = directory / "nodes.tsv"
    edges_path = directory / "edges.tsv"
    taxa_path.write_text(
        "direction\ttaxon_id\ttaxon_name\ttaxonomic_rank\n"
        "increased\tNCBITaxon:900\tFamilyaceae\tfamily\n"
        "decreased\tNCBITaxon:100\tGenus\tgenus\n"
        "increased\tNCBITaxon:101\tSpecies a\tspecies\n"
        "increased\tNCBITaxon:102\tSpecies b\tspecies\n",
        encoding="utf-8",
    )
    nodes_path.write_text(
        "id\tcategory\tname\n"
        "METPO:1000603\tbiolink:PhenotypicQuality\tanaerobic\n"
        "METPO:1000605\tbiolink:PhenotypicQuality\tfacultatively anaerobic\n"
        "METPO:1000608\tbiolink:PhenotypicQuality\tfacultatively aerobic\n"
        "METPO:1000668\tbiolink:PhenotypicQuality\tcoccus shaped\n"
        "METPO:2000103\tbiolink:OntologyClass\tcapable of\n"
        "GO:0030245\tbiolink:BiologicalProcess\tcellulose catabolic process\n",
        encoding="utf-8",
    )
    edges_path.write_text(
        "\n".join(["id\tsubject\tpredicate\tobject", *edges]) + "\n", encoding="utf-8"
    )
    return taxa_path, nodes_path, edges_path


def _matrix_lines(directory: Path, edges: list[str], **options: bool) -> list[str]:
    directory.mkdir()
    taxa_path, nodes_path, edges_path = _write_strain_fixture(directory, edges)
    output_path = directory / "taxa_trait_matrix.tsv"
    taxa_trait_matrix.create_taxa_trait_matrix(
        taxa_path=taxa_path,
        kg_dir=directory,
        output_path=output_path,
        column_metadata_path=directory / "taxa_trait_matrix_columns.tsv",
        edges_path=edges_path,
        nodes_path=nodes_path,
        **options,
    )
    return output_path.read_text(encoding="utf-8").splitlines()


def test_strain_edges_never_reach_the_default_matrix(tmp_path: Path) -> None:
    """Without the rollup, strain records leave the matrix byte-identical (#774)."""
    with_strains = _matrix_lines(tmp_path / "a", [*TAXONOMY_AND_TAXON_EDGES, *STRAIN_EDGES])
    without_strains = _matrix_lines(tmp_path / "b", TAXONOMY_AND_TAXON_EDGES)

    assert (
        with_strains
        == without_strains
        == [
            "taxon_id\ttaxon_name\ttaxonomic_rank\tMETPO:1000603\tMETPO:1000605",
            "NCBITaxon:900\tFamilyaceae\tfamily\t1\t1",
            "NCBITaxon:100\tGenus\tgenus\t1\t1",
            "NCBITaxon:101\tSpecies a\tspecies\t1\t0",
            "NCBITaxon:102\tSpecies b\tspecies\t0\t1",
        ]
    )


def test_strain_rollup_counts_strain_traits_at_the_linked_ncbitaxon_node(tmp_path: Path) -> None:
    """A strain's traits count once at its species, then roll up like the species' own.

    Species a carries anaerobic itself and through strain 1, and facultatively anaerobic
    through two strains: each counts once at the species. The genus and family count
    annotated NCBITaxon nodes, so facultatively anaerobic is 2 there (species a through its
    strains, species b directly). Strains without an NCBITaxon link add nothing (#774).
    """
    lines = _matrix_lines(
        tmp_path / "on", [*TAXONOMY_AND_TAXON_EDGES, *STRAIN_EDGES], strain_rollup=True
    )

    assert lines == [
        "taxon_id\ttaxon_name\ttaxonomic_rank\tMETPO:1000603"
        "\tMETPO:2000103|GO:0030245\tMETPO:1000605",
        "NCBITaxon:900\tFamilyaceae\tfamily\t1\t1\t2",
        "NCBITaxon:100\tGenus\tgenus\t1\t1\t2",
        "NCBITaxon:101\tSpecies a\tspecies\t1\t0\t1",
        "NCBITaxon:102\tSpecies b\tspecies\t0\t1\t1",
    ]


def test_strain_links_keep_ncbitaxon_parents_only(tmp_path: Path) -> None:
    """Strain 5's subclass_of a kgmicrobe.genus node is not a taxon link (#797)."""
    _, _, edges_path = _write_strain_fixture(tmp_path, [*TAXONOMY_AND_TAXON_EDGES, *STRAIN_EDGES])

    edges = taxa_trait_matrix.read_kg_trait_edges(edges_path, strains=True)

    assert edges.strain_taxa == {
        "kgmicrobe.strain:1": {"NCBITaxon:101"},
        "kgmicrobe.strain:2": {"NCBITaxon:101"},
        "kgmicrobe.strain:3": {"NCBITaxon:102"},
    }
    assert edges.strain_traits["kgmicrobe.strain:5"] == {"METPO:1000668"}
    default = taxa_trait_matrix.read_kg_trait_edges(edges_path)
    assert (default.strain_taxa, default.strain_traits) == ({}, {})


def test_roll_up_strain_traits_adds_to_every_linked_node_without_mutating() -> None:
    taxon_traits = {"NCBITaxon:101": {"METPO:1000603"}}
    strain_taxa = {"kgmicrobe.strain:1": {"NCBITaxon:101", "NCBITaxon:201"}}
    strain_traits = {
        "kgmicrobe.strain:1": {"METPO:1000603", "METPO:1000605"},
        "kgmicrobe.strain:2": {"METPO:1000608"},
    }

    rolled = taxa_trait_matrix.roll_up_strain_traits(taxon_traits, strain_taxa, strain_traits)

    assert rolled == {
        "NCBITaxon:101": {"METPO:1000603", "METPO:1000605"},
        "NCBITaxon:201": {"METPO:1000603", "METPO:1000605"},
    }
    assert taxon_traits == {"NCBITaxon:101": {"METPO:1000603"}}


def test_strain_rollup_is_opt_in_on_the_command_line() -> None:
    parser = taxa_trait_matrix.build_parser()
    assert parser.parse_args([]).strain_rollup is False
    assert parser.parse_args(["--strain-rollup"]).strain_rollup is True


@pytest.mark.parametrize("output", [[], ["--output", str(taxa_trait_matrix.DEFAULT_OUTPUT_PATH)]])
def test_strain_rollup_refuses_the_default_matrix_path(
    output: list[str], tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """The analyses read the default path, so a rolled-up matrix must go elsewhere (#797).

    The empty --kg-dir stops a build before it writes anything if the guard is missing.
    """
    with pytest.raises(SystemExit) as exited:
        taxa_trait_matrix.main(["--strain-rollup", "--kg-dir", str(tmp_path), *output])

    assert exited.value.code == 2
    assert "--strain-rollup needs an --output" in capsys.readouterr().err


def test_strain_rollup_writes_an_explicit_output(tmp_path: Path) -> None:
    taxa_path, nodes_path, edges_path = _write_strain_fixture(
        tmp_path, [*TAXONOMY_AND_TAXON_EDGES, *STRAIN_EDGES]
    )
    output_path = tmp_path / "rolled" / "taxa_trait_matrix.tsv"
    paths = ["--taxa", taxa_path, "--nodes", nodes_path, "--edges", edges_path]

    exit_code = taxa_trait_matrix.main(
        ["--strain-rollup", *map(str, paths), "--output", str(output_path)]
    )

    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8").splitlines()[0].endswith("\tMETPO:1000605")
