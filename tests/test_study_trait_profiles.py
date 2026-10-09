import csv
from pathlib import Path

import yaml

from bsdbng.study_trait_profiles import create_study_trait_profiles
from bsdbng.taxa_trait_matrix import create_taxa_trait_matrix


def _write_study_yaml(study_dir: Path) -> None:
    study = {
        "id": "bsdb:99",
        "experiments": [
            {
                "id": "bsdb:99-1",
                "host_species": "Homo sapiens",
                "signatures": [
                    {
                        "id": "bsdb:99-1-1",
                        "direction": "increased",
                        "taxa": [
                            {"id": "NCBITaxon:1", "taxon_name": "Taxon one"},
                            {"id": "NCBITaxon:2", "taxon_name": "Taxon two"},
                        ],
                    },
                    {
                        "id": "bsdb:99-1-2",
                        "direction": "increased",
                        "taxa": [
                            {"id": "NCBITaxon:1", "taxon_name": "Taxon one"},
                        ],
                    },
                    {
                        "id": "bsdb:99-1-3",
                        "direction": "decreased",
                        "taxa": [
                            {"id": "NCBITaxon:3", "taxon_name": "Taxon three"},
                            {"id": "NCBITaxon:999", "taxon_name": "Missing taxon"},
                        ],
                    },
                ],
            }
        ],
    }
    (study_dir / "bsdb_99.yaml").write_text(
        yaml.dump(study, default_flow_style=False, sort_keys=False),
        encoding="utf-8",
    )


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_create_study_trait_profiles_sums_taxa_traits_by_direction(tmp_path: Path) -> None:
    study_dir = tmp_path / "studies"
    study_dir.mkdir()
    _write_study_yaml(study_dir)

    taxa_trait_matrix_path = tmp_path / "taxa_trait_matrix.tsv"
    taxa_trait_matrix_path.write_text(
        "\n".join(
            [
                "taxon_id\ttaxon_name\ttaxonomic_rank\tMETPO:1\tMETPO:2\tMETPO:3",
                "NCBITaxon:1\tTaxon one\tspecies\t2\t0\t1",
                "NCBITaxon:2\tTaxon two\tspecies\t0\t3\t1",
                "NCBITaxon:3\tTaxon three\tspecies\t1\t1\t0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    output_path = tmp_path / "study_trait_profiles_wide.tsv"

    stats = create_study_trait_profiles(
        study_dir=study_dir,
        taxa_trait_matrix_path=taxa_trait_matrix_path,
        output_path=output_path,
    )

    assert stats.profiles == 2
    assert stats.trait_columns == 3
    assert stats.taxon_observations == 5
    assert stats.profile_taxa == 4
    assert stats.matched_profile_taxa == 3
    assert stats.unmatched_profile_taxa == 1
    assert _read_tsv(output_path) == [
        {
            "study_id": "bsdb:99",
            "experiment_id": "bsdb:99-1",
            "inc_dec": "dec",
            "METPO:1": "1",
            "METPO:2": "1",
            "METPO:3": "0",
        },
        {
            "study_id": "bsdb:99",
            "experiment_id": "bsdb:99-1",
            "inc_dec": "inc",
            "METPO:1": "2",
            "METPO:2": "3",
            "METPO:3": "2",
        },
    ]


def test_create_study_trait_profiles_can_count_duplicate_taxa(tmp_path: Path) -> None:
    study_dir = tmp_path / "studies"
    study_dir.mkdir()
    _write_study_yaml(study_dir)

    taxa_trait_matrix_path = tmp_path / "taxa_trait_matrix.tsv"
    taxa_trait_matrix_path.write_text(
        "\n".join(
            [
                "taxon_id\ttaxon_name\ttaxonomic_rank\tMETPO:1\tMETPO:2",
                "NCBITaxon:1\tTaxon one\tspecies\t2\t1",
                "NCBITaxon:2\tTaxon two\tspecies\t0\t3",
                "NCBITaxon:3\tTaxon three\tspecies\t1\t1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    output_path = tmp_path / "study_trait_profiles_wide.tsv"

    stats = create_study_trait_profiles(
        study_dir=study_dir,
        taxa_trait_matrix_path=taxa_trait_matrix_path,
        output_path=output_path,
        deduplicate_taxa=False,
    )

    assert stats.profile_taxa == 5
    rows = _read_tsv(output_path)
    assert rows[1]["inc_dec"] == "inc"
    assert rows[1]["METPO:1"] == "4"
    assert rows[1]["METPO:2"] == "5"


def test_create_study_trait_profiles_writes_long_format(tmp_path: Path) -> None:
    study_dir = tmp_path / "studies"
    study_dir.mkdir()
    _write_study_yaml(study_dir)

    taxa_trait_matrix_path = tmp_path / "taxa_trait_matrix.tsv"
    taxa_trait_matrix_path.write_text(
        "\n".join(
            [
                "taxon_id\ttaxon_name\ttaxonomic_rank\tMETPO:1",
                "NCBITaxon:1\tTaxon one\tspecies\t2",
                "NCBITaxon:2\tTaxon two\tspecies\t3",
                "NCBITaxon:3\tTaxon three\tspecies\t1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    output_path = tmp_path / "study_trait_profiles_long.tsv"

    create_study_trait_profiles(
        study_dir=study_dir,
        taxa_trait_matrix_path=taxa_trait_matrix_path,
        output_path=output_path,
        output_format="long",
    )

    assert output_path.read_text(encoding="utf-8").splitlines() == [
        "trait\tbsdb:99_bsdb:99-1_dec\tbsdb:99_bsdb:99-1_inc",
        "METPO:1\t1\t5",
    ]


def test_create_study_trait_profiles_includes_all_mapped_kgmicrobe_traits(
    tmp_path: Path,
) -> None:
    study_dir = tmp_path / "studies"
    study_dir.mkdir()
    study = {
        "id": "bsdb:100",
        "experiments": [
            {
                "id": "bsdb:100-1",
                "signatures": [
                    {
                        "id": "bsdb:100-1-1",
                        "direction": "increased",
                        "taxa": [
                            {
                                "id": "NCBITaxon:900",
                                "taxon_name": "Familyaceae",
                                "taxonomic_rank": "family",
                            },
                            {
                                "id": "NCBITaxon:101",
                                "taxon_name": "Species a",
                                "taxonomic_rank": "species",
                            },
                        ],
                    },
                    {
                        "id": "bsdb:100-1-2",
                        "direction": "decreased",
                        "taxa": [
                            {
                                "id": "NCBITaxon:100",
                                "taxon_name": "Genus",
                                "taxonomic_rank": "genus",
                            },
                        ],
                    },
                ],
            }
        ],
    }
    (study_dir / "bsdb_100.yaml").write_text(
        yaml.dump(study, default_flow_style=False, sort_keys=False),
        encoding="utf-8",
    )

    taxa_path = tmp_path / "bugsigdb_taxa.tsv"
    taxa_path.write_text(
        "\n".join(
            [
                (
                    "study_id\texperiment_id\tsignature_id\tdirection\ttaxon_id"
                    "\ttaxon_name\ttaxonomic_rank"
                ),
                "bsdb:100\tbsdb:100-1\tbsdb:100-1-1\tincreased\tNCBITaxon:900\tFamilyaceae\tfamily",
                "bsdb:100\tbsdb:100-1\tbsdb:100-1-1\tincreased\tNCBITaxon:101\tSpecies a\tspecies",
                "bsdb:100\tbsdb:100-1\tbsdb:100-1-2\tdecreased\tNCBITaxon:100\tGenus\tgenus",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    nodes_path = tmp_path / "nodes.tsv"
    nodes_path.write_text(
        "\n".join(
            [
                "id\tcategory\tname",
                "METPO:1000001\tbiolink:PhenotypicQuality\tfamily trait",
                "METPO:1000002\tbiolink:PhenotypicQuality\tgenus trait",
                "METPO:1000003\tbiolink:PhenotypicQuality\tspecies trait",
                "METPO:1000004\tbiolink:PhenotypicQuality\tunreported descendant trait",
                "METPO:2000006\tbiolink:OntologyClass\tuses as carbon source",
                "CHEBI:123\tbiolink:SmallMolecule\tacetate",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    edges_path = tmp_path / "edges.tsv"
    edges_path.write_text(
        "\n".join(
            [
                "id\tsubject\tpredicate\tobject",
                "e1\tNCBITaxon:100\tbiolink:subclass_of\tNCBITaxon:900",
                "e2\tNCBITaxon:101\tbiolink:subclass_of\tNCBITaxon:100",
                "e3\tNCBITaxon:102\tbiolink:subclass_of\tNCBITaxon:100",
                "e4\tNCBITaxon:900\tbiolink:has_phenotype\tMETPO:1000001",
                "e5\tNCBITaxon:100\tbiolink:has_phenotype\tMETPO:1000002",
                "e6\tNCBITaxon:101\tbiolink:has_phenotype\tMETPO:1000003",
                "e7\tNCBITaxon:102\tbiolink:has_phenotype\tMETPO:1000004",
                "e8\tNCBITaxon:101\tMETPO:2000006\tCHEBI:123",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    taxa_trait_matrix_path = tmp_path / "taxa_trait_matrix.tsv"
    matrix_stats = create_taxa_trait_matrix(
        taxa_path=taxa_path,
        kg_dir=tmp_path,
        output_path=taxa_trait_matrix_path,
        column_metadata_path=tmp_path / "taxa_trait_matrix_columns.tsv",
        edges_path=edges_path,
        nodes_path=nodes_path,
    )
    output_path = tmp_path / "study_trait_profiles_wide.tsv"

    profile_stats = create_study_trait_profiles(
        study_dir=study_dir,
        taxa_trait_matrix_path=taxa_trait_matrix_path,
        output_path=output_path,
    )

    assert matrix_stats.reported_taxa == 3
    assert matrix_stats.traits == 5
    assert profile_stats.profile_taxa == 3
    assert profile_stats.matched_profile_taxa == 3
    assert profile_stats.unmatched_profile_taxa == 0

    matrix_rows = _read_tsv(taxa_trait_matrix_path)
    assert {row["taxon_id"] for row in matrix_rows} == {
        "NCBITaxon:100",
        "NCBITaxon:101",
        "NCBITaxon:900",
    }

    rows = {row["inc_dec"]: row for row in _read_tsv(output_path)}
    assert rows["inc"] == {
        "study_id": "bsdb:100",
        "experiment_id": "bsdb:100-1",
        "inc_dec": "inc",
        "METPO:1000001": "1",
        "METPO:1000002": "1",
        "METPO:1000003": "2",
        "METPO:1000004": "1",
        "METPO:2000006|CHEBI:123": "2",
    }
    assert rows["dec"] == {
        "study_id": "bsdb:100",
        "experiment_id": "bsdb:100-1",
        "inc_dec": "dec",
        "METPO:1000001": "0",
        "METPO:1000002": "1",
        "METPO:1000003": "1",
        "METPO:1000004": "1",
        "METPO:2000006|CHEBI:123": "1",
    }


def test_create_study_trait_profiles_binary_with_taxa_count_and_host_filter(
    tmp_path: Path,
) -> None:
    study_dir = tmp_path / "studies"
    study_dir.mkdir()
    study = {
        "id": "bsdb:99",
        "experiments": [
            {
                "id": "bsdb:99-1",
                "host_species": "Homo sapiens",
                "signatures": [
                    {
                        "id": "bsdb:99-1-1",
                        "direction": "increased",
                        "taxa": [
                            {"id": "NCBITaxon:1", "taxon_name": "Taxon one"},
                            {"id": "NCBITaxon:2", "taxon_name": "Taxon two"},
                        ],
                    },
                    {
                        "id": "bsdb:99-1-2",
                        "direction": "decreased",
                        "taxa": [
                            {"id": "NCBITaxon:3", "taxon_name": "Taxon three"},
                            {"id": "NCBITaxon:999", "taxon_name": "Missing taxon"},
                        ],
                    },
                ],
            },
            {
                "id": "bsdb:99-2",
                "host_species": "Mus musculus",
                "signatures": [
                    {
                        "id": "bsdb:99-2-1",
                        "direction": "increased",
                        "taxa": [
                            {"id": "NCBITaxon:1", "taxon_name": "Taxon one"},
                        ],
                    },
                ],
            },
        ],
    }
    (study_dir / "bsdb_99.yaml").write_text(
        yaml.dump(study, default_flow_style=False, sort_keys=False),
        encoding="utf-8",
    )

    taxa_trait_matrix_path = tmp_path / "taxa_trait_matrix.tsv"
    taxa_trait_matrix_path.write_text(
        "\n".join(
            [
                "taxon_id\ttaxon_name\ttaxonomic_rank\tMETPO:1\tMETPO:2\tMETPO:3",
                "NCBITaxon:1\tTaxon one\tspecies\t2\t0\t1",
                "NCBITaxon:2\tTaxon two\tspecies\t0\t3\t1",
                "NCBITaxon:3\tTaxon three\tspecies\t1\t1\t0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    # 1. Opt-in: binary=True, include_taxa_count=True, host_species="Homo sapiens"
    output_path = tmp_path / "study_trait_profiles_ora.tsv"
    stats = create_study_trait_profiles(
        study_dir=study_dir,
        taxa_trait_matrix_path=taxa_trait_matrix_path,
        output_path=output_path,
        binary=True,
        include_taxa_count=True,
        host_species="Homo sapiens",
    )

    assert stats.profiles == 2
    rows = _read_tsv(output_path)
    # Only Homo sapiens experiments present
    assert {r["experiment_id"] for r in rows} == {"bsdb:99-1"}

    # In dec: taxa are NCBITaxon:3 (matched) and NCBITaxon:999 (unmatched).
    # Matched taxa count N = 1 (satisfies #350)
    dec_row = next(r for r in rows if r["inc_dec"] == "dec")
    assert dec_row["taxa_count"] == "1"
    assert dec_row["METPO:1"] == "1"
    assert dec_row["METPO:2"] == "1"
    assert dec_row["METPO:3"] == "0"

    # In inc: taxa are NCBITaxon:1 and NCBITaxon:2.
    # Matched taxa count N = 2
    inc_row = next(r for r in rows if r["inc_dec"] == "inc")
    assert inc_row["taxa_count"] == "2"
    # Binary: even if NCBITaxon:1 has METPO:1=2, binary contribution is 1
    assert inc_row["METPO:1"] == "1"
    assert inc_row["METPO:2"] == "1"
    assert inc_row["METPO:3"] == "2"

    # 2. Default: host_species=None includes Mus musculus
    output_path_all = tmp_path / "study_trait_profiles_all_hosts.tsv"
    stats_all = create_study_trait_profiles(
        study_dir=study_dir,
        taxa_trait_matrix_path=taxa_trait_matrix_path,
        output_path=output_path_all,
    )
    assert stats_all.profiles == 3
    rows_all = _read_tsv(output_path_all)
    assert {r["experiment_id"] for r in rows_all} == {"bsdb:99-1", "bsdb:99-2"}
    assert "taxa_count" not in rows_all[0]
