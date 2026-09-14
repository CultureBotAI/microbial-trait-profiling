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


def test_load_raw_experiment_hosts_maps_bsdb_id_to_experiment_id(tmp_path: Path) -> None:
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

    hosts = heatmap.load_raw_experiment_hosts(full_dump_path)

    assert hosts.loc[("bsdb:study1", "bsdb:study1-1"), "host_species"] == "Homo sapiens"
    assert hosts.loc[("bsdb:study1", "bsdb:study1-2"), "host_species"] == "Mus musculus"
    assert hosts.loc[("bsdb:study2", "bsdb:study2-3"), "host_species"] == "Homo sapiens"
