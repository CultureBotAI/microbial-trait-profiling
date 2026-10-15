import importlib.util
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

MODULE_PATH = (
    Path(__file__).resolve().parents[1] / "scripts" / "create_trait_profile_infographic.py"
)
SPEC = importlib.util.spec_from_file_location("create_trait_profile_infographic", MODULE_PATH)
assert SPEC is not None
assert SPEC.loader is not None
infographic = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(infographic)


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
    assert stats.loc["METPO:1", "pair_count"] == 4
    assert stats.loc["METPO:1", "study_count"] == 4
    assert "ci95_lower" in stats.columns
    assert "ci95_upper" in stats.columns
    assert "study_ci95_lower" in stats.columns
    assert "study_ci95_upper" in stats.columns
    assert "study_p_value" in stats.columns
    assert "study_q_value" in stats.columns


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
        assert "study_count" in stats.columns
        assert "study_mean_increased_minus_decreased" in stats.columns
        assert "ci95_lower" in stats.columns
        assert "ci95_upper" in stats.columns
        assert "study_ci95_lower" in stats.columns
        assert "study_ci95_upper" in stats.columns
        assert "study_p_value" in stats.columns
        assert "study_q_value" in stats.columns

    assert (tmp_path / "trait_profile_paired_contrast_trait_stats_manifest.tsv").exists()

    uncapped_dir = tmp_path / "uncapped"
    uncapped_manifest = infographic.write_stratified_trait_stats(
        profile_matrix,
        pair_metadata,
        trait_labels,
        output_dir=uncapped_dir,
        min_pairs=2,
        min_profile_support=1,
        max_profile_support_fraction=1.0,
    )

    assert set(uncapped_manifest["selected_trait_count"]) == {3}
    for token in ["feces", "feces_excluded", "oral"]:
        output = uncapped_dir / f"trait_profile_paired_contrast_trait_stats_{token}.tsv"
        stats = pd.read_csv(output, sep="\t")
        assert len(stats) == 3


def test_compute_trait_stats_study_aggregation_and_confidence_intervals() -> None:
    # 4 pairs across 3 studies: study1 has 2 pairs, study2 has 1 pair, study3 has 1 pair
    matrix = pd.DataFrame(
        {
            "METPO:1": [1.0, 3.0, 2.0, 4.0],
            "METPO:2": [-2.0, -4.0, -1.0, -3.0],
        },
        index=["p1", "p2", "p3", "p4"],
    )
    study_ids = pd.Series(["s1", "s1", "s2", "s3"], index=["p1", "p2", "p3", "p4"])
    trait_metadata = pd.DataFrame(
        {
            "trait_id": ["METPO:1", "METPO:2"],
            "trait_label": ["oxygen METPO:1", "anaerobic METPO:2"],
            "selected_trait_rank": [1, 2],
            "profile_support": [4, 4],
            "total_signal": [10.0, 10.0],
            "composition_variance": [0.2, 0.2],
            "selection_score": [0.5, 0.5],
        },
        index=["METPO:1", "METPO:2"],
    )
    annotation_counts = pd.Series({"METPO:1": 4, "METPO:2": 4})

    stats = infographic.compute_trait_stats(
        matrix,
        trait_metadata,
        annotation_counts,
        stratum_kind="pooled",
        stratum="pooled",
        selected_profile_count=8,
        study_ids=study_ids,
    ).set_index("trait")

    # Pair-level:
    assert stats.loc["METPO:1", "pair_count"] == 4
    assert np.isclose(stats.loc["METPO:1", "mean_increased_minus_decreased"], 2.5)
    # Pair std = sqrt(5/3) = 1.290994, SE = 0.645497, df = 3, t_crit = 3.182446
    # 95% CI: [0.44575, 4.55425]
    assert np.isclose(stats.loc["METPO:1", "ci95_lower"], 2.5 - 3.182446 * 0.645497, atol=1e-3)
    assert np.isclose(stats.loc["METPO:1", "ci95_upper"], 2.5 + 3.182446 * 0.645497, atol=1e-3)

    # Study-level:
    # s1 mean: (1+3)/2 = 2.0; s2: 2.0; s3: 4.0 -> values [2.0, 2.0, 4.0]
    # study count = 3
    assert stats.loc["METPO:1", "study_count"] == 3
    assert np.isclose(stats.loc["METPO:1", "study_mean_increased_minus_decreased"], 8.0 / 3.0)
    # study variance = 4/3, std = sqrt(4/3), SE = 2/3, t = 4.0, df = 2
    assert np.isclose(stats.loc["METPO:1", "study_t_statistic"], 4.0)
    # t_crit(df=2) = 4.30265
    assert np.isclose(
        stats.loc["METPO:1", "study_ci95_lower"],
        8.0 / 3.0 - 4.30265 * (2.0 / 3.0),
        atol=1e-3,
    )
    assert np.isclose(
        stats.loc["METPO:1", "study_ci95_upper"],
        8.0 / 3.0 + 4.30265 * (2.0 / 3.0),
        atol=1e-3,
    )
    assert stats.loc["METPO:1", "study_p_value"] < 0.06
    assert np.isfinite(stats.loc["METPO:1", "study_q_value"])
    assert np.isclose(
        stats.loc["METPO:1", "study_minus_log10_p"],
        -np.log10(stats.loc["METPO:1", "study_p_value"]),
    )


def test_plot_top_traits_renders_with_study_level_stats() -> None:
    trait_stats = pd.DataFrame(
        {
            "trait": ["METPO:1", "METPO:2"],
            "trait_label": ["aerobic METPO:1", "anaerobic METPO:2"],
            "study_minus_log10_p": [8.5, 5.2],
            "study_mean_increased_minus_decreased": [0.21, -0.15],
            "study_ci95_lower": [0.14, -0.22],
            "study_ci95_upper": [0.28, -0.08],
            "minus_log10_p": [12.0, 8.0],
            "mean_increased_minus_decreased": [0.25, -0.18],
            "ci95_lower": [0.18, -0.25],
            "ci95_upper": [0.32, -0.11],
        }
    )
    fig, ax = plt.subplots()
    try:
        infographic.plot_top_traits(ax, trait_stats)
        assert "study-level p" in ax.get_xlabel()
        texts = [t.get_text() for t in ax.texts]
        assert any("[+0.14, +0.28]" in t for t in texts)
        assert any("[-0.22, -0.08]" in t for t in texts)
    finally:
        plt.close(fig)
