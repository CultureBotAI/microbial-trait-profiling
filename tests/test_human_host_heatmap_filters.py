import importlib.util
from pathlib import Path
from types import ModuleType


def load_script_module(script_name: str) -> ModuleType:
    module_path = Path(__file__).resolve().parents[1] / "scripts" / script_name
    spec = importlib.util.spec_from_file_location(script_name.removesuffix(".py"), module_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_taxa_guild_raw_records_keep_only_human_hosts(tmp_path: Path) -> None:
    module = load_script_module("create_taxa_guild_heatmap.py")
    full_dump_path = tmp_path / "full_dump.csv"
    full_dump_path.write_text(
        "\n".join(
            [
                "# comment",
                (
                    "BSDB ID,Abundance in Group 1,Host species,NCBI Taxonomy IDs,"
                    "MetaPhlAn taxon names,Taxon"
                ),
                "bsdb:study1/1/1,increased,Homo sapiens,562,s__Escherichia coli,",
                "bsdb:study2/1/1,increased,Mus musculus,573,s__Klebsiella pneumoniae,",
            ]
        ),
        encoding="utf-8",
    )

    records = module.load_bugsigdb_taxa_records(tmp_path / "unused.tsv", full_dump_path)

    assert records["study_id"].tolist() == ["bsdb:study1"]
    assert records["experiment_id"].tolist() == ["bsdb:study1-1"]
    assert records["taxon_id"].tolist() == ["NCBITaxon:562"]


def test_condition_taxon_pairs_keep_human_experiments(tmp_path: Path) -> None:
    module = load_script_module("create_condition_taxon_heatmap.py")
    full_dump_path = tmp_path / "full_dump.csv"
    full_dump_path.write_text(
        "\n".join(
            [
                "# comment",
                "BSDB ID,Condition,Host species",
                "bsdb:study1/1/1,Colorectal cancer,Homo sapiens",
                "bsdb:study1/2/1,Colorectal cancer,Mus musculus",
                "bsdb:study2/1/1,Obesity,Homo sapiens",
            ]
        ),
        encoding="utf-8",
    )

    pairs = module.load_study_condition_pairs(full_dump_path)

    assert pairs.to_dict("records") == [
        {
            "study_id": "bsdb:study1",
            "experiment_id": "bsdb:study1-1",
            "condition": "Colorectal cancer",
        },
        {
            "study_id": "bsdb:study2",
            "experiment_id": "bsdb:study2-1",
            "condition": "Obesity",
        },
    ]


def test_clustered_heatmap_matrix_keeps_human_experiments(tmp_path: Path) -> None:
    module = load_script_module("create_clustered_heatmap.py")
    full_dump_path = tmp_path / "full_dump.csv"
    full_dump_path.write_text(
        "\n".join(
            [
                "# comment",
                "BSDB ID,Host species",
                "bsdb:study1/1/1,Homo sapiens",
                "bsdb:study1/2/1,Mus musculus",
            ]
        ),
        encoding="utf-8",
    )
    taxa_path = tmp_path / "taxa.tsv"
    taxa_path.write_text(
        "\n".join(
            [
                "study_id\texperiment_id\tdirection\ttaxon_name\ttaxonomic_rank",
                "bsdb:study1\tbsdb:study1-1\tincreased\tEscherichia\tgenus",
                "bsdb:study1\tbsdb:study1-2\tincreased\tKlebsiella\tgenus",
            ]
        ),
        encoding="utf-8",
    )

    matrix, inc_dec = module.build_binary_matrix(
        taxa_path,
        top_taxa=10,
        max_profiles=10,
        full_dump_path=full_dump_path,
    )

    assert list(matrix.index) == ["bsdb:study1|bsdb:study1-1|inc"]
    assert list(matrix.columns) == ["Escherichia"]
    assert inc_dec.to_dict() == {"bsdb:study1|bsdb:study1-1|inc": "inc"}
