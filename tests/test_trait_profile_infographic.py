import importlib.util
from pathlib import Path

import pandas as pd

MODULE_PATH = (
    Path(__file__).resolve().parents[1] / "scripts" / "create_trait_profile_infographic.py"
)
SPEC = importlib.util.spec_from_file_location("create_trait_profile_infographic", MODULE_PATH)
assert SPEC is not None
assert SPEC.loader is not None
infographic = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(infographic)


def test_load_profile_matrix_reads_compressed_fallback(tmp_path: Path) -> None:
    matrix_path = tmp_path / "study_trait_profiles_wide.tsv"
    compressed_path = Path(f"{matrix_path}.gz")
    pd.DataFrame(
        {
            "study_id": ["study1", "study1"],
            "experiment_id": ["exp1", "exp1"],
            "inc_dec": ["inc", "dec"],
            "METPO:1": [2, 0],
        }
    ).to_csv(compressed_path, sep="\t", index=False)

    matrix = infographic.load_profile_matrix(matrix_path)

    assert list(matrix.index) == ["study1|exp1|inc", "study1|exp1|dec"]
    assert matrix.loc["study1|exp1|inc", "METPO:1"] == 2


def test_load_trait_labels_falls_back_to_profile_matrix_header(tmp_path: Path) -> None:
    matrix_path = tmp_path / "study_trait_profiles_wide.tsv"
    compressed_path = Path(f"{matrix_path}.gz")
    missing_trait_columns_path = tmp_path / "taxa_trait_matrix_columns.tsv"
    pd.DataFrame(
        {
            "study_id": ["study1"],
            "experiment_id": ["exp1"],
            "inc_dec": ["inc"],
            "METPO:1": [2],
            "pathways:trait": [1],
        }
    ).to_csv(compressed_path, sep="\t", index=False)

    labels = infographic.load_trait_labels(
        missing_trait_columns_path,
        profile_matrix_path=matrix_path,
    )

    assert list(labels.index) == ["METPO:1", "pathways:trait"]
    assert labels.loc["METPO:1", "trait_label"] == "METPO:1"
    assert labels.loc["pathways:trait", "display_category"] == "unknown"


def test_load_paired_inputs_maps_display_columns_to_stable_trait_ids(tmp_path: Path) -> None:
    matrix_path = tmp_path / "paired_matrix.tsv"
    metadata_path = tmp_path / "paired_metadata.tsv"
    trait_metadata_path = tmp_path / "trait_metadata.tsv"

    pd.DataFrame(
        {
            "oxygen\nMETPO:1": [1.0, 2.0],
            "nitrate reduction\npathways:nitrate_reduction": [0.0, 1.0],
        },
        index=["pair1", "pair2"],
    ).to_csv(matrix_path, sep="\t")
    pd.DataFrame(
        {
            "study_id": ["study1", "study2"],
            "increased_profile_ids": ["study1|exp1|inc", "study2|exp1|inc"],
            "decreased_profile_ids": ["study1|exp1|dec", "study2|exp1|dec"],
        },
        index=["pair1", "pair2"],
    ).to_csv(metadata_path, sep="\t")
    pd.DataFrame(
        {
            "trait_id": ["METPO:1", "METPO:2000103|pathways:nitrate_reduction"],
            "display_label": [
                "oxygen\nMETPO:1",
                "nitrate reduction\npathways:nitrate_reduction",
            ],
        },
        index=["METPO:1", "METPO:2000103|pathways:nitrate_reduction"],
    ).to_csv(trait_metadata_path, sep="\t")

    matrix, _metadata, trait_metadata = infographic.load_paired_inputs(
        matrix_path,
        metadata_path,
        trait_metadata_path,
    )

    assert list(matrix.columns) == ["METPO:1", "METPO:2000103|pathways:nitrate_reduction"]
    assert trait_metadata.loc["METPO:2000103|pathways:nitrate_reduction", "trait_label"] == (
        "nitrate reduction pathways:nitrate_reduction"
    )


def test_compute_trait_stats_writes_stable_ids_and_annotation_depth() -> None:
    matrix = pd.DataFrame(
        {
            "METPO:1": [1.0, 2.0, 3.0, 4.0],
            "METPO:2000103|pathways:nitrate_reduction": [0.0, 0.0, 0.0, 0.0],
        }
    )
    trait_metadata = pd.DataFrame(
        {
            "trait_id": ["METPO:1", "METPO:2000103|pathways:nitrate_reduction"],
            "trait_label": ["oxygen METPO:1", "nitrate reduction pathways:nitrate_reduction"],
            "selected_trait_rank": [1, 2],
            "profile_support": [4, 0],
            "total_signal": [10.0, 0.0],
            "composition_variance": [0.2, 0.0],
            "selection_score": [0.5, 0.0],
        },
        index=["METPO:1", "METPO:2000103|pathways:nitrate_reduction"],
    )
    annotation_counts = pd.Series({"METPO:1": 4, "METPO:2000103|pathways:nitrate_reduction": 0})

    stats = infographic.compute_trait_stats(
        matrix,
        trait_metadata,
        annotation_counts,
        stratum_kind="pooled",
        stratum="pooled",
        selected_profile_count=8,
    ).set_index("trait")

    assert "nitrate reduction\npathways:nitrate_reduction" not in stats.index
    assert stats.loc["METPO:1", "trait_label"] == "oxygen METPO:1"
    assert stats.loc["METPO:1", "n_profiles_where_annotated"] == 4
    assert stats.loc["METPO:1", "annotated_profile_fraction"] == 0.5
    assert stats.loc["METPO:2000103|pathways:nitrate_reduction", "q_value"] == 1.0


def test_write_stratified_trait_stats_emits_feces_excluded_and_body_site_files(
    tmp_path: Path,
) -> None:
    profile_matrix = pd.DataFrame(
        {
            "METPO:1": [10, 1, 9, 1, 0, 0, 1, 0],
            "METPO:2": [1, 10, 1, 9, 5, 1, 5, 1],
            "METPO:3": [0, 0, 5, 0, 10, 1, 9, 1],
        },
        index=[
            "study1|exp1|inc",
            "study1|exp1|dec",
            "study2|exp1|inc",
            "study2|exp1|dec",
            "study3|exp1|inc",
            "study3|exp1|dec",
            "study4|exp1|inc",
            "study4|exp1|dec",
        ],
    )
    pair_metadata = pd.DataFrame(
        {
            "study_id": ["study1", "study2", "study3", "study4"],
            "body_site_group": ["Feces", "Feces", "Oral", "Oral"],
            "increased_profile_ids": [
                "study1|exp1|inc",
                "study2|exp1|inc",
                "study3|exp1|inc",
                "study4|exp1|inc",
            ],
            "decreased_profile_ids": [
                "study1|exp1|dec",
                "study2|exp1|dec",
                "study3|exp1|dec",
                "study4|exp1|dec",
            ],
        },
        index=["pair1", "pair2", "pair3", "pair4"],
    )
    trait_labels = pd.DataFrame(
        {
            "trait_label": ["oxygen METPO:1", "nitrate METPO:2", "motility METPO:3"],
            "display_label": ["oxygen METPO:1", "nitrate METPO:2", "motility METPO:3"],
            "column_label": ["oxygen METPO:1", "nitrate METPO:2", "motility METPO:3"],
            "object_name": ["oxygen", "nitrate", "motility"],
            "object_category": ["biolink:OntologyClass"] * 3,
            "display_category": ["OntologyClass"] * 3,
        },
        index=["METPO:1", "METPO:2", "METPO:3"],
    )

    manifest = infographic.write_stratified_trait_stats(
        profile_matrix,
        pair_metadata,
        trait_labels,
        output_dir=tmp_path,
        min_pairs=2,
        max_traits=2,
        min_profile_support=1,
        max_profile_support_fraction=1.0,
    )

    assert set(manifest["filename_token"]) == {"feces", "feces_excluded", "oral"}
    for token in ["feces", "feces_excluded", "oral"]:
        output = tmp_path / f"trait_profile_paired_contrast_trait_stats_{token}.tsv"
        assert output.exists()
        stats = pd.read_csv(output, sep="\t")
        assert len(stats) == 2
        assert stats["trait"].str.startswith("METPO:").all()
        assert "selected_trait_rank" in stats.columns
        assert "n_profiles_where_annotated" in stats.columns

    assert (tmp_path / "trait_profile_paired_contrast_trait_stats_manifest.tsv").exists()
