import importlib.util
from pathlib import Path

import pandas as pd
import pytest

MODULE_PATH = (
    Path(__file__).resolve().parents[1] / "scripts" / "create_study_trait_profile_heatmap.py"
)
SPEC = importlib.util.spec_from_file_location("create_study_trait_profile_heatmap", MODULE_PATH)
assert SPEC is not None
assert SPEC.loader is not None
heatmap = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(heatmap)


def test_filter_human_host_profiles_keeps_only_homo_sapiens() -> None:
    profile_matrix = pd.DataFrame(
        {"METPO:1": [1, 2, 3, 4]},
        index=[
            "study1|exp1|inc",
            "study1|exp1|dec",
            "study2|exp1|inc",
            "study3|exp1|inc",
        ],
    )
    profile_metadata = pd.DataFrame(
        {
            "study_id": ["study1", "study1", "study2", "study3"],
            "experiment_id": ["exp1", "exp1", "exp1", "exp1"],
            "inc_dec": ["inc", "dec", "inc", "inc"],
            "host_species": [
                "Homo sapiens",
                " Homo sapiens ",
                "Mus musculus",
                "Canis lupus familiaris",
            ],
        },
        index=profile_matrix.index,
    )

    filtered_matrix, filtered_metadata = heatmap.filter_human_host_profiles(
        profile_matrix,
        profile_metadata,
    )

    assert list(filtered_matrix.index) == ["study1|exp1|inc", "study1|exp1|dec"]
    assert list(filtered_metadata.index) == ["study1|exp1|inc", "study1|exp1|dec"]
    assert set(filtered_metadata["host_species"]) == {"Homo sapiens"}


def test_filter_human_host_profiles_fails_without_host_species() -> None:
    profile_matrix = pd.DataFrame({"METPO:1": [1]}, index=["study1|exp1|inc"])
    profile_metadata = pd.DataFrame({"study_id": ["study1"]}, index=profile_matrix.index)

    with pytest.raises(ValueError, match="missing host_species"):
        heatmap.filter_human_host_profiles(profile_matrix, profile_metadata)


def test_load_raw_experiment_metadata_maps_bsdb_id_to_experiment_id(tmp_path: Path) -> None:
    full_dump_path = tmp_path / "full_dump.csv"
    full_dump_path.write_text(
        "\n".join(
            [
                "# comment",
                "BSDB ID,Host species",
                "bsdb:study1/1/1,Homo sapiens",
                "bsdb:study1/2/1,Mus musculus",
                "bsdb:study2/3/1, Homo sapiens ",
            ]
        ),
        encoding="utf-8",
    )

    hosts = heatmap.load_raw_experiment_metadata(full_dump_path)

    assert hosts.loc[("bsdb:study1", "bsdb:study1-1"), "host_species"] == "Homo sapiens"
    assert hosts.loc[("bsdb:study1", "bsdb:study1-2"), "host_species"] == "Mus musculus"
    assert hosts.loc[("bsdb:study2", "bsdb:study2-3"), "host_species"] == "Homo sapiens"
    assert hosts.loc[("bsdb:study2", "bsdb:study2-3"), "condition"] == "unknown"


def write_two_experiment_study(tmp_path: Path) -> tuple[Path, Path]:
    """One study whose two experiments differ in condition and body site (#555)."""
    full_dump_path = tmp_path / "full_dump.csv"
    full_dump_path.write_text(
        "\n".join(
            [
                "# BugSigDB 2026-04-24_00:41_UTC",
                "BSDB ID,Condition,Body site,Location of subjects,Host species",
                "bsdb:study1/1/1,Crohn's disease,Feces,Sweden,Homo sapiens",
                "bsdb:study1/1/2,Crohn's disease,Feces,Sweden,Homo sapiens",
                "bsdb:study1/2/1,Ulcerative colitis,Colonic mucosa,Norway,Homo sapiens",
            ]
        ),
        encoding="utf-8",
    )
    studies_path = tmp_path / "studies.tsv"
    studies_path.write_text(
        "study_id\ttitle\tpublication_year\tsignature_count\ttaxon_count\n"
        "bsdb:study1\tInflammatory bowel disease cohort\t2020\t3\t12\n",
        encoding="utf-8",
    )
    return full_dump_path, studies_path


def test_join_profile_labels_takes_each_profiles_own_experiment_labels(tmp_path: Path) -> None:
    full_dump_path, studies_path = write_two_experiment_study(tmp_path)
    profile_ids = [
        "bsdb:study1|bsdb:study1-1|inc",
        "bsdb:study1|bsdb:study1-1|dec",
        "bsdb:study1|bsdb:study1-2|inc",
        "bsdb:study1|bsdb:study1-2|dec",
    ]
    profile_metadata = pd.DataFrame(
        {
            "study_id": ["bsdb:study1"] * 4,
            "experiment_id": ["bsdb:study1-1", "bsdb:study1-1", "bsdb:study1-2", "bsdb:study1-2"],
            "inc_dec": ["inc", "dec", "inc", "dec"],
        },
        index=profile_ids,
    )

    labelled = heatmap.join_profile_labels(profile_metadata, studies_path, full_dump_path)

    assert labelled["condition"].tolist() == [
        "Crohn's disease",
        "Crohn's disease",
        "Ulcerative colitis",
        "Ulcerative colitis",
    ]
    assert labelled["body_site"].tolist() == ["Feces", "Feces", "Colonic mucosa", "Colonic mucosa"]
    assert labelled["location"].tolist() == ["Sweden", "Sweden", "Norway", "Norway"]
    study_columns = [
        "study_majority_condition",
        "study_majority_body_site",
        "study_majority_location",
        "study_condition_count",
        "study_body_site_count",
    ]
    assert labelled[study_columns].drop_duplicates().to_dict("records") == [
        {
            "study_majority_condition": "Crohn's disease",
            "study_majority_body_site": "Feces",
            "study_majority_location": "Sweden",
            "study_condition_count": 2,
            "study_body_site_count": 2,
        }
    ]
    assert set(labelled["title"]) == {"Inflammatory bowel disease cohort"}


def test_load_raw_experiment_metadata_rejects_conflicting_experiment_labels(
    tmp_path: Path,
) -> None:
    full_dump_path = tmp_path / "full_dump.csv"
    full_dump_path.write_text(
        "\n".join(
            [
                "BSDB ID,Condition,Body site,Location of subjects,Host species",
                "bsdb:study1/1/1,Crohn's disease,Feces,Sweden,Homo sapiens",
                "bsdb:study1/1/2,Ulcerative colitis,Feces,Sweden,Homo sapiens",
            ]
        ),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="more than one condition label"):
        heatmap.load_raw_experiment_metadata(full_dump_path)
