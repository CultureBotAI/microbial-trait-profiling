import csv
from pathlib import Path

import yaml

from bsdbng.study_tsv import generate_tables


def _write_study_yaml(study_dir: Path) -> None:
    study = {
        "id": "bsdb:99",
        "source_record_id": "Study 99",
        "pmid": 12345678,
        "title": "Test Study Title",
        "publication_year": 2024,
        "doi": "10.1234/test",
        "url": None,
        "experiments": [
            {
                "id": "bsdb:99-1",
                "experiment_name": "Experiment 1",
                "group_0_name": "lean",
                "group_1_name": "obese",
                "signatures": [
                    {
                        "id": "bsdb:99-1-1",
                        "direction": "increased",
                        "taxa": [
                            {
                                "id": "NCBITaxon:1350",
                                "taxon_name": "Enterococcus",
                                "taxonomic_rank": "genus",
                            },
                            {
                                "id": "NCBITaxon:1485",
                                "taxon_name": "Clostridium",
                                "taxonomic_rank": "genus",
                            },
                        ],
                    },
                    {
                        "id": "bsdb:99-1-2",
                        "direction": "decreased",
                        "taxa": [
                            {
                                "id": "NCBITaxon:5678",
                                "taxon_name": "Bacteroides fragilis",
                                "taxonomic_rank": "species",
                            }
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


def test_generate_tables_writes_flattened_outputs(tmp_path: Path) -> None:
    study_dir = tmp_path / "studies"
    output_dir = tmp_path / "derived"
    study_dir.mkdir()
    _write_study_yaml(study_dir)

    tables = generate_tables(study_dir, output_dir)

    assert tables.taxa == output_dir / "taxa.tsv"

    study_rows = _read_tsv(tables.studies)
    assert study_rows == [
        {
            "study_id": "bsdb:99",
            "source_record_id": "Study 99",
            "pmid": "12345678",
            "title": "Test Study Title",
            "publication_year": "2024",
            "doi": "10.1234/test",
            "url": "",
            "experiment_count": "1",
            "signature_count": "2",
            "taxon_count": "3",
            "input_file": "bsdb_99.yaml",
        }
    ]

    experiment_rows = _read_tsv(tables.experiments)
    assert experiment_rows == [
        {
            "study_id": "bsdb:99",
            "experiment_id": "bsdb:99-1",
            "experiment_name": "Experiment 1",
            "group_0_name": "lean",
            "group_1_name": "obese",
            "signature_count": "2",
            "taxon_count": "3",
        }
    ]

    signature_rows = _read_tsv(tables.signatures)
    assert signature_rows == [
        {
            "study_id": "bsdb:99",
            "experiment_id": "bsdb:99-1",
            "signature_id": "bsdb:99-1-1",
            "direction": "increased",
            "taxon_count": "2",
        },
        {
            "study_id": "bsdb:99",
            "experiment_id": "bsdb:99-1",
            "signature_id": "bsdb:99-1-2",
            "direction": "decreased",
            "taxon_count": "1",
        },
    ]

    taxon_rows = _read_tsv(tables.taxa)
    assert len(taxon_rows) == 3
    assert taxon_rows[0]["taxon_name"] == "Enterococcus"
    assert taxon_rows[2]["taxonomic_rank"] == "species"


def test_generate_tables_writes_profile_matrices(tmp_path: Path) -> None:
    study_dir = tmp_path / "studies"
    output_dir = tmp_path / "derived"
    study_dir.mkdir()
    _write_study_yaml(study_dir)

    tables = generate_tables(study_dir, output_dir)

    wide_rows = _read_tsv(tables.taxon_profiles_wide)
    assert wide_rows == [
        {
            "study_id": "bsdb:99",
            "experiment_id": "bsdb:99-1",
            "inc_dec": "dec",
            "Bacteroides fragilis": "1",
            "Clostridium": "0",
            "Enterococcus": "0",
        },
        {
            "study_id": "bsdb:99",
            "experiment_id": "bsdb:99-1",
            "inc_dec": "inc",
            "Bacteroides fragilis": "0",
            "Clostridium": "1",
            "Enterococcus": "1",
        },
    ]

    long_rows = _read_tsv(tables.taxon_profiles_long)
    assert long_rows == [
        {
            "taxon": "Bacteroides fragilis",
            "bsdb:99_bsdb:99-1_dec": "1",
            "bsdb:99_bsdb:99-1_inc": "0",
        },
        {
            "taxon": "Clostridium",
            "bsdb:99_bsdb:99-1_dec": "0",
            "bsdb:99_bsdb:99-1_inc": "1",
        },
        {
            "taxon": "Enterococcus",
            "bsdb:99_bsdb:99-1_dec": "0",
            "bsdb:99_bsdb:99-1_inc": "1",
        },
    ]
