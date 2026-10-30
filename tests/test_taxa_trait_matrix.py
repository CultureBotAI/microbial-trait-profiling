import csv
import json
from pathlib import Path
from typing import TypedDict

import pytest

from bsdbng.taxa_trait_matrix import create_taxa_trait_matrix
from mtp import taxa_trait_matrix
from mtp.taxa_trait_matrix import StrainRollup


class MatrixOptions(TypedDict, total=False):
    strain_rollup: StrainRollup


def _metadata_rows(path: Path) -> list[list[str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    assert rows[0] == [
        "column_id",
        "column_label",
        "predicate_id",
        "predicate_name",
        "object_id",
        "object_name",
        "object_category",
        "trait_key_policy",
        "trait_member_ids_json",
        "trait_source_vocabularies_json",
    ]
    assert all(len(row) == 10 for row in rows)
    assert {row[7] for row in rows[1:]} == {"facultative_union_v1"}
    return rows


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
    metadata = _metadata_rows(column_metadata_path)
    assert ["\t".join(row[:7]) for row in metadata] == [
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
    assert [json.loads(row[8]) for row in metadata[1:]] == [
        ["METPO:1000603"],
        ["GO:0008150"],
        ["METPO:1000636"],
        ["METPO:1000615"],
        ["CHEBI:123"],
    ]
    assert [json.loads(row[9]) for row in metadata[1:]] == [[], [], [], [], []]


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


# Strain records for the scope rules (#788). Strain 1 also asserts a clinical class that the
# phenotype scope excludes; strains 1 and 2 disagree on glucose fermentation; strain 3 has a
# negated relation to a METPO object; species 102 asserts that it does not grow in a medium.
SCOPE_EDGES = [
    "n1\tkgmicrobe.strain:1\tbiolink:has_phenotype\tMETPO:1004004",
    "n2\tkgmicrobe.strain:1\tMETPO:2000011\tCHEBI:17234",
    "n3\tkgmicrobe.strain:2\tMETPO:2000037\tCHEBI:17234",
    "n4\tkgmicrobe.strain:2\tbiolink:has_phenotype\tMETPO:1000678",
    "n5\tkgmicrobe.strain:3\tMETPO:2000022\tMETPO:1000683",
    "n6\tNCBITaxon:102\tMETPO:2000518\tmediadive.medium:1",
]
ALL_EDGES = [*TAXONOMY_AND_TAXON_EDGES, *STRAIN_EDGES, *SCOPE_EDGES]


def _matrix_lines(
    directory: Path, edges: list[str], strain_rollup: StrainRollup | None = None
) -> list[str]:
    directory.mkdir()
    taxa_path, nodes_path, edges_path = _write_strain_fixture(directory, edges)
    output_path = directory / "taxa_trait_matrix.tsv"
    options: MatrixOptions = {} if strain_rollup is None else {"strain_rollup": strain_rollup}
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


def test_strain_edges_never_reach_the_off_matrix(tmp_path: Path) -> None:
    """Under off, strain records leave the matrix byte-identical (#774, #788)."""
    with_strains = _matrix_lines(tmp_path / "a", [*TAXONOMY_AND_TAXON_EDGES, *STRAIN_EDGES], "off")
    without_strains = _matrix_lines(tmp_path / "b", TAXONOMY_AND_TAXON_EDGES)

    assert (
        with_strains
        == without_strains
        == [
            "taxon_id\ttaxon_name\ttaxonomic_rank\tMETPO:1000603\tMETPO:1000608",
            "NCBITaxon:900\tFamilyaceae\tfamily\t1\t1",
            "NCBITaxon:100\tGenus\tgenus\t1\t1",
            "NCBITaxon:101\tSpecies a\tspecies\t1\t0",
            "NCBITaxon:102\tSpecies b\tspecies\t0\t1",
        ]
    )


def test_the_matrix_builder_defaults_to_the_primary_scope(tmp_path: Path) -> None:
    """The legacy path builds the owner's primary, as the core workflow does (#788)."""
    assert taxa_trait_matrix.PRIMARY_STRAIN_ROLLUP == "phenotype"
    default = _matrix_lines(tmp_path / "default", ALL_EDGES)

    assert default == _matrix_lines(tmp_path / "phenotype", ALL_EDGES, "phenotype")
    assert default != _matrix_lines(tmp_path / "off", ALL_EDGES, "off")


@pytest.mark.parametrize(
    ("strain_rollup", "cellulose"),
    [("phenotype", []), ("full", ["\tMETPO:2000103|GO:0030245"])],
)
def test_strain_rollup_counts_strain_traits_at_the_linked_ncbitaxon_node(
    tmp_path: Path, strain_rollup: StrainRollup, cellulose: list[str]
) -> None:
    """A strain's traits count once at its species, then roll up like the species' own.

    Species a carries anaerobic itself and through strain 1, and facultative oxygen
    preference through two strains: each counts once at the species. The genus and family
    count annotated NCBITaxon nodes, so facultative oxygen preference is 2 there (species a
    through its strains, species b directly). Strains without an NCBITaxon link add nothing
    (#774). The phenotype scope leaves out strain 3's capability edge, which has no bare
    METPO phenotype object (#788).
    """
    lines = _matrix_lines(
        tmp_path / "on", [*TAXONOMY_AND_TAXON_EDGES, *STRAIN_EDGES], strain_rollup
    )

    cell = ["\t1"] if cellulose else []
    no_cell = ["\t0"] if cellulose else []
    assert lines == [
        "taxon_id\ttaxon_name\ttaxonomic_rank\tMETPO:1000603"
        + "".join(cellulose)
        + "\tMETPO:1000608",
        "NCBITaxon:900\tFamilyaceae\tfamily\t1" + "".join(cell) + "\t2",
        "NCBITaxon:100\tGenus\tgenus\t1" + "".join(cell) + "\t2",
        "NCBITaxon:101\tSpecies a\tspecies\t1" + "".join(no_cell) + "\t1",
        "NCBITaxon:102\tSpecies b\tspecies\t0" + "".join(cell) + "\t1",
    ]


@pytest.mark.parametrize(
    ("strain_rollup", "added", "removed"),
    [
        ("off", set(), set()),
        ("phenotype", {"METPO:1000678"}, set()),
        (
            "full",
            {
                "METPO:1000678",
                "METPO:1004004",
                "METPO:2000011|CHEBI:17234",
                "METPO:2000103|GO:0030245",
            },
            {"METPO:2000518|mediadive.medium:1"},
        ),
    ],
)
def test_each_scope_has_its_own_vocabulary(
    tmp_path: Path, strain_rollup: StrainRollup, added: set[str], removed: set[str]
) -> None:
    """Off keeps the direct vocabulary, phenotype adds bare METPO phenotypes only, and
    full adds every strain key that does not assert an absence while dropping every one
    that does, the species' own "does not grow in" included (#788)."""
    direct = {"METPO:1000603", "METPO:1000608", "METPO:2000518|mediadive.medium:1"}

    header = _matrix_lines(tmp_path / strain_rollup, ALL_EDGES, strain_rollup)[0]

    assert set(header.split("\t")[3:]) == (direct | added) - removed


def test_phenotype_scope_excludes_handling_and_clinical_classes(tmp_path: Path) -> None:
    """Biosafety levels and pathogen labels are not ecophysiological traits (#788)."""
    assert set(taxa_trait_matrix.PHENOTYPE_EXCLUDED_TRAITS) == {
        "METPO:1001102",
        "METPO:1001103",
        "METPO:1001104",
        "METPO:1004002",
        "METPO:1004003",
        "METPO:1004004",
    }
    for trait in taxa_trait_matrix.PHENOTYPE_EXCLUDED_TRAITS:
        assert not taxa_trait_matrix.strain_trait_rolls_up(trait, "phenotype")
        assert taxa_trait_matrix.strain_trait_rolls_up(trait, "full")

    header = _matrix_lines(tmp_path / "phenotype", ALL_EDGES, "phenotype")[0]

    assert "METPO:1004004" not in header.split("\t")


def test_full_scope_never_carries_an_activity_and_its_absence(tmp_path: Path) -> None:
    """Strains 1 and 2 of species a disagree on glucose fermentation; under full only the
    activity survives, and no node keeps any negated key (#788)."""
    _, _, edges_path = _write_strain_fixture(tmp_path, ALL_EDGES)
    edges = taxa_trait_matrix.read_kg_trait_edges(edges_path, strains=True)

    traits = taxa_trait_matrix.scoped_taxon_traits(
        edges.taxon_traits, edges.strain_taxa, edges.strain_traits, "full"
    )

    assert "METPO:2000037|CHEBI:17234" in edges.strain_traits["kgmicrobe.strain:2"]
    assert "METPO:2000011|CHEBI:17234" in traits["NCBITaxon:101"]
    assert not any(taxa_trait_matrix.is_negated(key) for keys in traits.values() for key in keys)


def test_a_negated_relation_to_a_metpo_object_keeps_its_predicate(tmp_path: Path) -> None:
    """The absence of a phenotype is never keyed as the phenotype itself (#788)."""
    assert taxa_trait_matrix.make_trait_key("METPO:2000022", "METPO:1000683") == (
        "METPO:2000022|METPO:1000683"
    )
    assert taxa_trait_matrix.make_trait_key("biolink:has_phenotype", "METPO:1000683") == (
        "METPO:1000683"
    )
    _, _, edges_path = _write_strain_fixture(tmp_path, ALL_EDGES)

    keys = taxa_trait_matrix.read_kg_trait_edges(edges_path, strains=True).strain_traits

    assert "METPO:2000022|METPO:1000683" in keys["kgmicrobe.strain:3"]
    assert "METPO:1000683" not in keys["kgmicrobe.strain:3"]


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
    scoped = taxa_trait_matrix.read_kg_trait_edges(
        edges_path, strains=True, strain_rollup="phenotype"
    )
    assert scoped.strain_traits["kgmicrobe.strain:3"] == set()


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


def test_strain_rollup_is_the_primary_by_default_on_the_command_line(
    capsys: pytest.CaptureFixture[str],
) -> None:
    parser = taxa_trait_matrix.build_parser()
    assert parser.parse_args([]).strain_rollup == "phenotype"
    for scope in ("off", "phenotype", "full"):
        assert parser.parse_args(["--strain-rollup", scope]).strain_rollup == scope
    with pytest.raises(SystemExit):
        parser.parse_args(["--strain-rollup", "strains"])
    assert "invalid choice" in capsys.readouterr().err


@pytest.mark.parametrize("output", [[], ["--output", str(taxa_trait_matrix.DEFAULT_OUTPUT_PATH)]])
@pytest.mark.parametrize("scope", ["off", "full"])
def test_a_secondary_scope_refuses_the_default_matrix_path(
    output: list[str], scope: str, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """The analyses read the default path as the primary, so off and full go elsewhere.

    The empty --kg-dir stops a build before it writes anything if the guard is missing
    (#797, #788).
    """
    with pytest.raises(SystemExit) as exited:
        taxa_trait_matrix.main(["--strain-rollup", scope, "--kg-dir", str(tmp_path), *output])

    assert exited.value.code == 2
    assert f"--strain-rollup {scope} needs an --output" in capsys.readouterr().err


@pytest.mark.parametrize("scope", [[], ["--strain-rollup", "phenotype"]])
def test_the_primary_scope_may_build_the_default_matrix(scope: list[str], tmp_path: Path) -> None:
    """The guard lets the primary through; the empty --kg-dir then stops the build."""
    with pytest.raises(FileNotFoundError, match="KG file not found"):
        taxa_trait_matrix.main([*scope, "--kg-dir", str(tmp_path)])


def test_strain_rollup_writes_an_explicit_output(tmp_path: Path) -> None:
    taxa_path, nodes_path, edges_path = _write_strain_fixture(
        tmp_path, [*TAXONOMY_AND_TAXON_EDGES, *STRAIN_EDGES]
    )
    output_path = tmp_path / "rolled" / "taxa_trait_matrix.tsv"
    paths = ["--taxa", taxa_path, "--nodes", nodes_path, "--edges", edges_path]

    exit_code = taxa_trait_matrix.main(
        ["--strain-rollup", "phenotype", *map(str, paths), "--output", str(output_path)]
    )

    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8").splitlines()[0].endswith("\tMETPO:1000608")


def test_facultative_synonyms_unify_into_single_nonredundant_feature(tmp_path: Path) -> None:
    """METPO:1000605, METPO:1000608, and METPO:1000612 collapse to METPO:1000608 (#1123)."""
    for synonym in ("METPO:1000605", "METPO:1000608", "METPO:1000612"):
        assert (
            taxa_trait_matrix.make_trait_key("biolink:has_phenotype", synonym)
            == taxa_trait_matrix.FACULTATIVE_TRAIT_ID
            == "METPO:1000608"
        )
    assert (
        taxa_trait_matrix.make_trait_key("METPO:2000022", "METPO:1000605")
        == "METPO:2000022|METPO:1000608"
    )

    taxa_path = tmp_path / "bugsigdb_taxa.tsv"
    nodes_path = tmp_path / "nodes.tsv"
    edges_path = tmp_path / "edges.tsv"
    output_path = tmp_path / "taxa_trait_matrix.tsv"
    column_metadata_path = tmp_path / "taxa_trait_matrix_columns.tsv"

    taxa_path.write_text(
        "direction\ttaxon_id\ttaxon_name\ttaxonomic_rank\n"
        "increased\tNCBITaxon:100\tGenus\tgenus\n"
        "increased\tNCBITaxon:101\tSpecies a\tspecies\n",
        encoding="utf-8",
    )
    nodes_path.write_text(
        "id\tcategory\tname\n"
        "METPO:1000605\tbiolink:PhenotypicQuality\tfacultatively anaerobic\n"
        "METPO:1000608\tbiolink:PhenotypicQuality\tfacultatively aerobic\n"
        "METPO:1000612\tbiolink:PhenotypicQuality\tfacultative\n",
        encoding="utf-8",
    )
    edges_path.write_text(
        "id\tsubject\tpredicate\tobject\n"
        "e1\tNCBITaxon:101\tbiolink:subclass_of\tNCBITaxon:100\n"
        "e2\tNCBITaxon:101\tbiolink:has_phenotype\tMETPO:1000605\n"
        "e3\tNCBITaxon:101\tbiolink:has_phenotype\tMETPO:1000608\n"
        "e4\tNCBITaxon:101\tbiolink:has_phenotype\tMETPO:1000612\n"
        "s1\tkgmicrobe.strain:1\tbiolink:subclass_of\tNCBITaxon:101\n"
        "s2\tkgmicrobe.strain:1\tbiolink:has_phenotype\tMETPO:1000605\n"
        "s3\tkgmicrobe.strain:1\tbiolink:has_phenotype\tMETPO:1000612\n",
        encoding="utf-8",
    )

    stats = taxa_trait_matrix.create_taxa_trait_matrix(
        taxa_path=taxa_path,
        kg_dir=tmp_path,
        output_path=output_path,
        column_metadata_path=column_metadata_path,
        edges_path=edges_path,
        nodes_path=nodes_path,
        strain_rollup="phenotype",
    )

    assert stats.traits == 1
    assert output_path.read_text(encoding="utf-8").splitlines() == [
        "taxon_id\ttaxon_name\ttaxonomic_rank\tMETPO:1000608",
        "NCBITaxon:100\tGenus\tgenus\t1",
        "NCBITaxon:101\tSpecies a\tspecies\t1",
    ]
    metadata = _metadata_rows(column_metadata_path)
    assert ["\t".join(row[:7]) for row in metadata] == [
        (
            "column_id\tcolumn_label\tpredicate_id\tpredicate_name\tobject_id"
            "\tobject_name\tobject_category"
        ),
        (
            "METPO:1000608\tfacultative (aerobe or anaerobe) METPO:1000608\t\t\tMETPO:1000608"
            "\tfacultative (aerobe or anaerobe)\tbiolink:PhenotypicQuality"
        ),
    ]
    assert json.loads(metadata[1][8]) == ["METPO:1000605", "METPO:1000608", "METPO:1000612"]
    assert json.loads(metadata[1][9]) == [
        {
            "object_id": "METPO:1000605",
            "source_id": "infores:bacdive",
            "source_label": "facultative anaerobe",
        },
        {
            "object_id": "METPO:1000605",
            "source_id": "infores:madin_etal",
            "source_label": "facultative",
        },
        {
            "object_id": "METPO:1000608",
            "source_id": "infores:bacdive",
            "source_label": "facultative aerobe",
        },
        {
            "object_id": "METPO:1000608",
            "source_id": "infores:madin_etal",
            "source_label": "facultative",
        },
        {
            "object_id": "METPO:1000612",
            "source_id": "infores:bactotraits",
            "source_label": "Ox_facultative_aerobe_anaerobe",
        },
    ]
