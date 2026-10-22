import csv
import importlib.util
from pathlib import Path
from types import ModuleType

import pytest

DUMP_FIELDS = [
    "BSDB ID",
    "Abundance in Group 1",
    "Condition",
    "Host species",
    "Sequencing type",
    "NCBI Taxonomy IDs",
    "MetaPhlAn taxon names",
]
TAXA_FIELDS = [
    "study_id",
    "experiment_id",
    "signature_id",
    "direction",
    "taxon_id",
    "taxon_name",
    "taxonomic_rank",
]
PARKINSONS = "Parkinson's disease"
SCLEROSIS = "Multiple sclerosis"
METAPHLAN_NAMES = {
    "562": "s__Escherichia coli",
    "573": "s__Klebsiella pneumoniae",
    "817": "s__Bacteroides fragilis",
    "841": "g__Roseburia",
    "572511": "g__Blautia",
}


def load_script_module(script_name: str) -> ModuleType:
    module_path = Path(__file__).resolve().parents[1] / "scripts" / script_name
    spec = importlib.util.spec_from_file_location(script_name.removesuffix(".py"), module_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def dump_row(
    bsdb_id: str, direction: str, condition: str, host: str, sequencing: str, *taxa: str
) -> tuple[str, ...]:
    """One full_dump.csv signature row, naming each NCBI taxon from METAPHLAN_NAMES."""
    names = ",".join(METAPHLAN_NAMES[taxon] for taxon in taxa)
    return bsdb_id, direction, condition, host, sequencing, ";".join(taxa), names


def write_full_dump(path: Path, rows: list[tuple[str, ...]]) -> Path:
    """Write a commented full_dump.csv whose rows follow DUMP_FIELDS."""
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("# BugSigDB export\n")
        writer = csv.writer(handle)
        writer.writerow(DUMP_FIELDS)
        writer.writerows(rows)
    return path


def write_taxa_table(path: Path, rows: list[tuple[str, str, str, str, str, str, str]]) -> Path:
    """Write a bugsigdb_taxa.tsv from (study, experiment, signature, direction, taxon, name,
    rank) rows, expanding the short IDs into the table's ID formats."""
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(TAXA_FIELDS)
        for study, experiment, signature, direction, taxon, name, rank in rows:
            writer.writerow(
                [
                    f"bsdb:{study}",
                    f"bsdb:{study}-{experiment}",
                    f"bsdb:{study}-{experiment}-{signature}",
                    direction,
                    f"NCBITaxon:{taxon}",
                    name,
                    rank,
                ]
            )
    return path


def write_subset_run(tmp_path: Path) -> tuple[Path, Path]:
    """Write a dump of eight studies and the taxa table of a subset that keeps three.

    The subset keeps the Parkinson's disease experiments of study1 and study2, with
    study1 reversed onto the disease-vs-comparison basis, and the multiple sclerosis
    experiment of study7. Study3's Parkinson's disease experiment and all four obesity
    studies stay outside it; obesity has the most studies in the dump.
    """
    human = "Homo sapiens"
    dump = write_full_dump(
        tmp_path / "full_dump.csv",
        [
            dump_row("bsdb:study1/1/1", "increased", PARKINSONS, human, "16S", "562", "841"),
            dump_row("bsdb:study1/1/2", "decreased", PARKINSONS, human, "16S", "817"),
            dump_row("bsdb:study2/1/1", "increased", PARKINSONS, human, "WMS", "562"),
            dump_row("bsdb:study2/1/2", "decreased", PARKINSONS, human, "WMS", "572511"),
            dump_row("bsdb:study3/1/1", "increased", PARKINSONS, human, "16S", "562", "841"),
            dump_row("bsdb:study3/1/2", "decreased", PARKINSONS, human, "16S", "817"),
            *(
                dump_row(f"bsdb:study{n}/1/1", "increased", "Obesity", human, "16S", "841")
                for n in (4, 5, 6, 8)
            ),
            dump_row("bsdb:study7/1/1", "increased", SCLEROSIS, human, "16S", "841"),
            dump_row("bsdb:study7/1/2", "decreased", SCLEROSIS, human, "16S", "817"),
        ],
    )
    taxa = write_taxa_table(
        tmp_path / "bugsigdb_taxa.tsv",
        [
            ("study1", "1", "1", "decreased", "562", "Escherichia coli", "species"),
            ("study1", "1", "1", "decreased", "841", "Roseburia", "genus"),
            ("study1", "1", "2", "increased", "817", "Bacteroides fragilis", "species"),
            ("study2", "1", "1", "increased", "562", "Escherichia coli", "species"),
            ("study2", "1", "2", "decreased", "572511", "Blautia", "genus"),
            ("study7", "1", "1", "increased", "841", "Roseburia", "genus"),
            ("study7", "1", "2", "decreased", "817", "Bacteroides fragilis", "species"),
        ],
    )
    return dump, taxa


def write_all_host_run(tmp_path: Path) -> tuple[Path, Path]:
    """Write a dump and the all-host taxa table that an all-paired run exports with it."""
    dump = write_full_dump(
        tmp_path / "full_dump.csv",
        [
            dump_row("bsdb:study1/1/1", "increased", "Obesity", "Homo sapiens", "WMS", "562"),
            dump_row("bsdb:study2/1/1", "increased", "Obesity", "Mus musculus", "16S", "573"),
            dump_row("bsdb:study3/1/1", "increased", "Obesity", "Homo sapiens", "PCR", "817"),
        ],
    )
    # An older export left a placeholder name and a guessed rank on study1's taxon
    taxa = write_taxa_table(
        tmp_path / "bugsigdb_taxa.tsv",
        [
            ("study1", "1", "1", "increased", "562", "taxon_562", "genus"),
            ("study2", "1", "1", "increased", "573", "Klebsiella pneumoniae", "species"),
            ("study3", "1", "1", "increased", "817", "Bacteroides fragilis", "species"),
        ],
    )
    return dump, taxa


def test_taxa_guild_raw_records_keep_only_human_16s_and_shotgun_hosts(tmp_path: Path) -> None:
    """A human PCR experiment is outside the corpus too (#742).

    An all-paired run's taxa table holds every host, so the dump's corpus rule is what
    keeps the figure on human 16S and shotgun signatures; the dump's names and ranks
    replace the table's.
    """
    module = load_script_module("create_taxa_guild_heatmap.py")
    full_dump_path, taxa_path = write_all_host_run(tmp_path)

    records = module.load_bugsigdb_taxa_records(taxa_path, full_dump_path)

    assert records["study_id"].tolist() == ["bsdb:study1"]
    assert records["experiment_id"].tolist() == ["bsdb:study1-1"]
    assert records["taxon_id"].tolist() == ["NCBITaxon:562"]
    assert records["taxon_name"].tolist() == ["Escherichia coli"]
    assert records["taxonomic_rank"].tolist() == ["species"]


def test_taxa_guild_counts_only_the_subset_profiles_in_their_own_direction(
    tmp_path: Path,
) -> None:
    """A subset run summarizes its own signatures, in its own direction basis (#768).

    Study3 and the obesity studies are in the dump but not in the subset; any of their
    profiles would raise a study, profile or direction count below.
    """
    module = load_script_module("create_taxa_guild_heatmap.py")
    full_dump_path, taxa_path = write_subset_run(tmp_path)

    records = module.load_bugsigdb_taxa_records(taxa_path, full_dump_path)
    summary, rank, studies, profiles = module.load_bugsigdb_taxa_summary(
        taxa_path, "species", [], full_dump_path
    )

    assert sorted(set(records["study_id"])) == ["bsdb:study1", "bsdb:study2", "bsdb:study7"]
    assert (rank, studies, profiles) == ("species", 3, 6)
    counts = ["bugsigdb_study_count", "bugsigdb_profile_count", "increased", "decreased"]
    # Study1's reversed experiment counts E. coli as decreased and B. fragilis as increased
    assert summary.loc["NCBITaxon:562", counts].tolist() == [2, 2, 1, 1]
    assert summary.loc["NCBITaxon:817", counts].tolist() == [2, 2, 1, 1]
    assert summary.loc["NCBITaxon:562", "taxon_name"] == "Escherichia coli"


def test_taxa_guild_fails_without_the_full_dump(tmp_path: Path) -> None:
    """Without the dump no corpus rule applies, so nothing falls back to every host (#747)."""
    module = load_script_module("create_taxa_guild_heatmap.py")
    _, taxa_path = write_all_host_run(tmp_path)

    with pytest.raises(FileNotFoundError, match="is missing; it decides which signatures"):
        module.load_bugsigdb_taxa_records(taxa_path, tmp_path / "missing_full_dump.csv")


def test_condition_taxon_conditions_and_denominators_come_from_the_subset(
    tmp_path: Path,
) -> None:
    """Condition ranking and study denominators count the subset's experiments only (#768).

    In the dump, obesity has four studies and Parkinson's disease three; the subset holds
    two Parkinson's disease studies and one multiple sclerosis study.
    """
    module = load_script_module("create_condition_taxon_heatmap.py")
    full_dump_path, taxa_path = write_subset_run(tmp_path)

    matrix, study_counts = module.build_condition_genus_matrix(
        taxa_path, full_dump_path, top_conditions=2, top_taxa=10
    )

    assert study_counts.to_dict() == {PARKINSONS: 2, SCLEROSIS: 1}
    assert sorted(matrix.columns) == ["Blautia", "Roseburia"]
    assert matrix.loc[PARKINSONS, ["Roseburia", "Blautia"]].tolist() == [0.5, 0.5]
    assert matrix.loc[SCLEROSIS, ["Roseburia", "Blautia"]].tolist() == [1.0, 0.0]


def test_condition_taxon_pairs_keep_human_16s_and_shotgun_experiments(tmp_path: Path) -> None:
    """ITS and untyped human experiments are outside the corpus (#742)."""
    module = load_script_module("create_condition_taxon_heatmap.py")
    full_dump_path = tmp_path / "full_dump.csv"
    full_dump_path.write_text(
        "\n".join(
            [
                "# comment",
                "BSDB ID,Condition,Host species,Sequencing type",
                "bsdb:study1/1/1,Colorectal cancer,Homo sapiens,16S",
                "bsdb:study1/2/1,Colorectal cancer,Mus musculus,16S",
                "bsdb:study2/1/1,Obesity,Homo sapiens,WMS",
                "bsdb:study2/2/1,Obesity,Homo sapiens,ITS / ITS2",
                "bsdb:study3/1/1,Obesity,Homo sapiens,NA",
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


def test_clustered_heatmap_matrix_keeps_human_16s_and_shotgun_experiments(
    tmp_path: Path,
) -> None:
    """A human 18S experiment is outside the corpus (#742)."""
    module = load_script_module("create_clustered_heatmap.py")
    full_dump_path = tmp_path / "full_dump.csv"
    full_dump_path.write_text(
        "\n".join(
            [
                "# comment",
                "BSDB ID,Host species,Sequencing type",
                "bsdb:study1/1/1,Homo sapiens,16S",
                "bsdb:study1/2/1,Mus musculus,16S",
                "bsdb:study1/3/1,Homo sapiens,18S",
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
                "bsdb:study1\tbsdb:study1-3\tincreased\tCandida\tgenus",
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


def test_an_export_without_sequencing_types_fails_instead_of_dropping_the_filters(
    tmp_path: Path,
) -> None:
    """With no corpus row, both heatmaps fail rather than keep every host (#753)."""
    full_dump_path = tmp_path / "full_dump.csv"
    full_dump_path.write_text(
        "\n".join(
            [
                "BSDB ID,Abundance in Group 1,Host species,NCBI Taxonomy IDs,MetaPhlAn taxon names",
                "bsdb:study1/1/1,increased,Homo sapiens,562,s__Escherichia coli",
                "bsdb:study1/2/1,increased,Mus musculus,573,s__Klebsiella pneumoniae",
            ]
        ),
        encoding="utf-8",
    )
    taxa_path = tmp_path / "taxa.tsv"
    taxa_path.write_text(
        "study_id\texperiment_id\tdirection\ttaxon_name\ttaxonomic_rank\n"
        "bsdb:study1\tbsdb:study1-2\tincreased\tKlebsiella\tgenus\n",
        encoding="utf-8",
    )
    clustered = load_script_module("create_clustered_heatmap.py")
    with pytest.raises(ValueError, match="no Homo sapiens 16S or WMS rows"):
        clustered.build_binary_matrix(
            taxa_path, top_taxa=10, max_profiles=10, full_dump_path=full_dump_path
        )
    guild = load_script_module("create_taxa_guild_heatmap.py")
    with pytest.raises(ValueError, match="no Homo sapiens 16S or WMS taxon rows"):
        guild.load_bugsigdb_taxa_records(tmp_path / "unused.tsv", full_dump_path)
