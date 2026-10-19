import importlib.util
import sys
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest

MODULE_PATH = (
    Path(__file__).resolve().parents[1] / "scripts" / "create_condition_trait_module_heatmap.py"
)
SPEC = importlib.util.spec_from_file_location("create_condition_trait_module_heatmap", MODULE_PATH)
assert SPEC is not None
assert SPEC.loader is not None
heatmap = importlib.util.module_from_spec(SPEC)
with patch.object(sys, "path", [str(MODULE_PATH.parent), *sys.path]):
    SPEC.loader.exec_module(heatmap)


def test_load_aggregated_study_profiles_filters_to_human_experiments(tmp_path: Path) -> None:
    profile_path = tmp_path / "profiles.tsv"
    profile_path.write_text(
        "\n".join(
            [
                "study_id\texperiment_id\tinc_dec\tMETPO:1\tMETPO:2",
                "study1\tstudy1-1\tinc\t1\t0",
                "study1\tstudy1-1\tdec\t0\t2",
                "study2\tstudy2-1\tinc\t8\t8",
                "study3\tstudy3-1\tinc\t4\t0",
            ]
        ),
        encoding="utf-8",
    )

    study_matrix, profile_metadata = heatmap.load_aggregated_study_profiles(
        profile_path,
        human_experiment_ids={"study1-1", "study3-1"},
    )

    assert list(study_matrix.index) == ["study1", "study3"]
    assert study_matrix.loc["study1", "METPO:1"] == 1
    assert study_matrix.loc["study1", "METPO:2"] == 2
    assert profile_metadata.loc["study1", "profile_count"] == 2
    assert profile_metadata.loc["study3", "profile_count"] == 1

    # study_id is a column only, so the main() join and reset_index() work (#602).
    assert profile_metadata.index.name != "study_id"
    assert profile_metadata.reset_index()["study_id"].tolist() == ["study1", "study3"]
    titles = pd.DataFrame(
        {"title": ["first", "third"]}, index=pd.Index(["study1", "study3"], name="study_id")
    )
    assert profile_metadata.join(titles, on="study_id")["title"].tolist() == ["first", "third"]


PROFILE_HEADER = "study_id\texperiment_id\tinc_dec\tMETPO:1\tMETPO:2"
TWO_CONDITION_PROFILES = [
    "study1\tstudy1-1\tinc\t3\t0",
    "study1\tstudy1-1\tdec\t2\t1",
    "study1\tstudy1-2\tinc\t1\t2",
    "study1\tstudy1-3\tinc\t0\t5",
    "study2\tstudy2-1\tinc\t4\t4",
]


def two_condition_experiments() -> pd.DataFrame:
    """Study 1 has two Crohn's disease experiments and one ulcerative colitis experiment."""
    return pd.DataFrame(
        {
            "study_id": ["study1", "study1", "study1", "study2", "study3"],
            "experiment_id": ["study1-1", "study1-2", "study1-3", "study2-1", "study3-1"],
            "condition": [
                "Crohn's disease",
                "Crohn's disease",
                "Ulcerative colitis",
                "Ulcerative colitis",
                "Crohn's disease",
            ],
            "body_site": ["Feces", "Feces", "Colonic mucosa", "Colonic mucosa", "Feces"],
            "location": ["Sweden", "Sweden", "Sweden", "Norway", "Norway"],
            "host_species": ["Homo sapiens"] * 4 + ["Mus musculus"],
        }
    ).set_index(["study_id", "experiment_id"])


def assert_two_condition_modules(study_matrix: pd.DataFrame, study_metadata: pd.DataFrame) -> None:
    groups, group_metadata = heatmap.build_condition_modules(
        study_matrix, study_metadata, ["condition", "body_site"], 1, 1, 1, 1, False, 10
    )
    crohns = "Crohn's disease | Feces"
    colitis = "Ulcerative colitis | Colonic mucosa"
    counts = group_metadata[["study_count", "profile_count", "experiment_count"]]
    assert groups.loc[crohns].tolist() == [6, 3]
    assert groups.loc[colitis].tolist() == [4, 9]
    assert counts.loc[crohns].tolist() == [1, 3, 2]
    assert counts.loc[colitis].tolist() == [2, 2, 2]


def test_modules_take_each_experiments_own_condition_and_body_site(tmp_path: Path) -> None:
    """A two-condition study adds each experiment to its own module, not its majority's (#555)."""
    profile_path = tmp_path / "profiles.tsv"
    profile_path.write_text("\n".join([PROFILE_HEADER, *TWO_CONDITION_PROFILES]), encoding="utf-8")

    study_matrix, study_metadata = heatmap.load_aggregated_study_profiles(
        profile_path,
        experiment_metadata=two_condition_experiments(),
        group_by=["condition", "body_site"],
    )

    assert_two_condition_modules(study_matrix, study_metadata)


def test_modules_stay_aligned_when_the_human_host_filter_drops_rows(tmp_path: Path) -> None:
    """Dropped non-human rows must not shift labels onto other profiles' counts (#646)."""
    profile_path = tmp_path / "profiles.tsv"
    mouse_rows = ["study3\tstudy3-1\tinc\t50\t0", "study3\tstudy3-1\tdec\t0\t50"]
    rows = [*TWO_CONDITION_PROFILES[:2], *mouse_rows, *TWO_CONDITION_PROFILES[2:]]
    profile_path.write_text("\n".join([PROFILE_HEADER, *rows]), encoding="utf-8")

    study_matrix, study_metadata = heatmap.load_aggregated_study_profiles(
        profile_path,
        human_experiment_ids={"study1-1", "study1-2", "study1-3", "study2-1"},
        experiment_metadata=two_condition_experiments(),
        group_by=["condition", "body_site"],
    )

    assert "study3" not in set(study_metadata["study_id"])
    assert_two_condition_modules(study_matrix, study_metadata)


def test_experiment_metadata_without_group_by_raises(tmp_path: Path) -> None:
    """Without grouping fields every experiment would get its study's majority label (#646)."""
    profile_path = tmp_path / "profiles.tsv"
    profile_path.write_text("\n".join([PROFILE_HEADER, *TWO_CONDITION_PROFILES]), encoding="utf-8")

    no_group_by: list[list[str] | None] = [None, []]
    for group_by in no_group_by:
        with pytest.raises(ValueError, match="needs group_by"):
            heatmap.load_aggregated_study_profiles(
                profile_path,
                experiment_metadata=two_condition_experiments(),
                group_by=group_by,
            )
