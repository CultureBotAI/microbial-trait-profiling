import importlib.util
from pathlib import Path

MODULE_PATH = (
    Path(__file__).resolve().parents[1] / "scripts" / "create_condition_trait_module_heatmap.py"
)
SPEC = importlib.util.spec_from_file_location("create_condition_trait_module_heatmap", MODULE_PATH)
assert SPEC is not None
assert SPEC.loader is not None
heatmap = importlib.util.module_from_spec(SPEC)
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
