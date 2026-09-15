import importlib.util
import math
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pandas as pd

MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "create_six_shared_trait_cytoscape_heatmaps.py"
)
SPEC = importlib.util.spec_from_file_location(
    "create_six_shared_trait_cytoscape_heatmaps",
    MODULE_PATH,
)
assert SPEC is not None
assert SPEC.loader is not None
heatmaps = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = heatmaps
SPEC.loader.exec_module(heatmaps)


def test_prepare_pattern_cells_masks_by_fisher_p_value_and_support() -> None:
    view = SimpleNamespace(
        view_key="increased_grouped",
        view_label="Increased, grouped conditions",
        focus_direction="increased",
        condition_mode="grouped condition labels",
    )
    traits = pd.DataFrame(
        {
            "trait_id": ["trait1"],
            "trait_label": ["trait one"],
            "trait_rank": [0],
        }
    )
    edges = pd.DataFrame(
        [
            {
                "condition": "Condition A",
                "body_site_group": "Feces",
                "trait_id": "trait1",
                "trait_label": "trait one",
                "focus_direction": "increased",
                "study_count": 5,
                "focus_profile_count": 20,
                "focus_enriched_profile_count": 10,
                "focus_enriched_study_count": 4,
                "focus_no_enriched_trait_profile_count": 10,
                "comparator_profile_count": 20,
                "comparator_enriched_profile_count": 1,
                "focus_enriched_mean_guild_z": 1.2,
                "focus_enriched_minus_comparator_guild_z": 1.1,
                "rank_score": 1.0,
            },
            {
                "condition": "Condition B",
                "body_site_group": "Feces",
                "trait_id": "trait1",
                "trait_label": "trait one",
                "focus_direction": "increased",
                "study_count": 1,
                "focus_profile_count": 20,
                "focus_enriched_profile_count": 10,
                "focus_enriched_study_count": 1,
                "focus_no_enriched_trait_profile_count": 10,
                "comparator_profile_count": 20,
                "comparator_enriched_profile_count": 1,
                "focus_enriched_mean_guild_z": 1.2,
                "focus_enriched_minus_comparator_guild_z": 1.1,
                "rank_score": 1.0,
            },
        ]
    )

    cells = heatmaps.prepare_pattern_cells(edges, view, traits)

    assert bool(cells.loc[0, "displayed_in_heatmap"]) is True
    assert cells.loc[0, "fisher_p_value"] <= heatmaps.P_VALUE_CUTOFF
    assert bool(cells.loc[0, "passes_p_value_cutoff"]) is True
    assert bool(cells.loc[0, "selected_after_bh"]) is True
    assert cells.loc[0, "fisher_q_value"] <= heatmaps.Q_VALUE_CUTOFF
    assert math.isclose(
        cells.loc[0, "minus_log10_fisher_p_value"],
        -math.log10(cells.loc[0, "fisher_p_value"]),
    )
    assert bool(cells.loc[1, "displayed_in_heatmap"]) is False
    assert cells.loc[1, "blank_reason"] == "low positive-z profile or study support"


def test_prepare_contrast_cells_uses_signed_increased_minus_decreased_fraction() -> None:
    mode = heatmaps.ContrastMode(
        key="grouped",
        label="Grouped conditions",
        condition_mode="grouped condition labels",
        increased_view_key="increased_grouped",
        decreased_view_key="decreased_grouped",
    )
    traits = pd.DataFrame(
        {
            "trait_id": ["trait1"],
            "trait_label": ["trait one"],
            "trait_rank": [0],
        }
    )
    increased_edges = pd.DataFrame(
        [
            {
                "condition": "Condition A",
                "body_site_group": "Feces",
                "trait_id": "trait1",
                "trait_label": "trait one",
                "focus_direction": "increased",
                "study_count": 5,
                "focus_profile_count": 20,
                "focus_enriched_profile_count": 10,
                "comparator_profile_count": 20,
                "comparator_enriched_profile_count": 1,
            }
        ]
    )
    decreased_edges = pd.DataFrame(
        [
            {
                "condition": "Condition A",
                "body_site_group": "Feces",
                "trait_id": "trait1",
                "trait_label": "trait one",
                "focus_direction": "decreased",
                "study_count": 5,
                "focus_profile_count": 20,
                "focus_enriched_profile_count": 1,
                "comparator_profile_count": 20,
                "comparator_enriched_profile_count": 10,
            }
        ]
    )

    cells = heatmaps.prepare_contrast_cells(increased_edges, decreased_edges, mode, traits)

    assert bool(cells.loc[0, "displayed_in_heatmap"]) is True
    assert cells.loc[0, "increased_minus_decreased_positive_fraction"] == 0.45
    assert cells.loc[0, "fisher_p_value"] <= heatmaps.P_VALUE_CUTOFF
    assert cells.loc[0, "plot_signed_minus_log10_fisher_p_value"] > 0


def test_select_pattern_rows_prioritizes_significant_trait_recurrence() -> None:
    cells = pd.DataFrame(
        [
            {
                "condition_body_site_pair": "A | Feces",
                "displayed_in_heatmap": True,
                "selected_after_bh": False,
                "focus_positive_profile_count": 3,
                "plot_minus_log10_fisher_p_value": 2.1,
                "minus_log10_fisher_p_value": 2.1,
                "strength": 1.0,
            },
            {
                "condition_body_site_pair": "A | Feces",
                "displayed_in_heatmap": True,
                "selected_after_bh": False,
                "focus_positive_profile_count": 3,
                "plot_minus_log10_fisher_p_value": 2.0,
                "minus_log10_fisher_p_value": 2.0,
                "strength": 1.0,
            },
            {
                "condition_body_site_pair": "B | Feces",
                "displayed_in_heatmap": True,
                "selected_after_bh": True,
                "focus_positive_profile_count": 20,
                "plot_minus_log10_fisher_p_value": 8.0,
                "minus_log10_fisher_p_value": 8.0,
                "strength": 20.0,
            },
        ]
    )

    assert heatmaps.select_pattern_rows(cells, max_rows=1) == ["A | Feces"]


def test_contrasts_show_supported_effects_and_mark_strict_significance() -> None:
    cells = pd.DataFrame(
        {
            "condition_body_site_pair": ["A | Feces", "B | Feces", "C | Feces", "D | Feces"],
            "body_site_group": ["Feces"] * 4,
            "trait_id": ["trait1"] * 4,
            "fisher_p_value": [0.0099, 0.01, 0.0001, math.nan],
            "passes_support": [True, True, False, True],
            "increased_minus_decreased_positive_fraction": [0.5, 0.5, -0.5, math.nan],
            "increased_positive_profile_count": [10, 10, 1, 10],
            "decreased_positive_profile_count": [1, 1, 10, 0],
        }
    )
    finalized = heatmaps.finalize_inferential_cells(cells)
    assert heatmaps.select_contrast_rows(finalized) == ["A | Feces", "B | Feces"]
    assert heatmaps.select_contrast_rows(finalized.iloc[2:]) == []
    assert finalized["significant_in_heatmap"].tolist() == [True, False, False, False]
    assert finalized["displayed_in_heatmap"].tolist() == [True, True, False, False]
    assert finalized.loc[1, "blank_reason"] == ""
    annotations = heatmaps.annotation_matrix(
        finalized,
        cells["condition_body_site_pair"].tolist(),
        pd.DataFrame({"trait_id": ["trait1"]}),
        mode="contrast",
    )
    assert annotations["trait1"].tolist() == ["+0.50", "+0.50", "", ""]
    assert heatmaps.format_p_value(0.0099) == "<0.01"


def test_spearman_clustering_uses_pairwise_overlap_and_neutral_undefined_distance() -> None:
    matrix = pd.DataFrame(
        [
            [1.0, 2.0, 3.0, math.nan],
            [10.0, 20.0, 30.0, math.nan],
            [3.0, 2.0, 1.0, math.nan],
            [1.0, math.nan, math.nan, math.nan],
        ],
        index=["A | Feces", "B | Feces", "C | Feces", "D | Feces"],
        columns=["t1", "t2", "t3", "t4"],
    )
    cells = pd.DataFrame(
        {
            "condition_body_site_pair": matrix.index,
            "condition": ["A", "B", "C", "D"],
            "body_site_group": ["Feces"] * 4,
        }
    )

    clustering = heatmaps.build_spearman_clustering(
        matrix,
        cells,
        view_key="test",
        view_label="Test",
    )

    assert math.isclose(clustering.distances.loc["A | Feces", "B | Feces"], 0.0)
    assert math.isclose(clustering.distances.loc["A | Feces", "C | Feces"], 2.0)
    assert (
        clustering.distances.loc["A | Feces", "D | Feces"] == heatmaps.UNDEFINED_SPEARMAN_DISTANCE
    )
    assert bool(clustering.defined_correlations.loc["A | Feces", "D | Feces"]) is False
    assert sorted(clustering.clustered_row_labels) == sorted(matrix.index)
    low_information = clustering.tree_leaves.set_index("condition_body_site_pair")[
        "low_information_for_clustering"
    ]
    assert bool(low_information.loc["D | Feces"]) is True
    assert len(clustering.tree_nodes) == len(matrix) - 1

    membership = heatmaps.cut_grouped_dendrogram(clustering)
    assert membership["cluster_cut_fraction"].eq(heatmaps.GROUPED_CLUSTER_CUT_FRACTION).all()
    assert membership["cluster_id"].str.fullmatch(r"C\d+").all()
    top_cluster = membership.loc[membership["clustered_heatmap_order"].idxmax(), "cluster_id"]
    assert top_cluster == "C1"


def test_grouped_cluster_enrichment_labels_anatomy_and_condition() -> None:
    records = []
    for cluster_rank, (cluster_id, anatomy, condition) in enumerate(
        [("C1", "Feces", "Condition A"), ("C2", "Oral", "Condition B")],
        start=1,
    ):
        for offset in range(12):
            records.append(
                {
                    "view_key": "grouped_direction_contrast",
                    "cluster_id": cluster_id,
                    "cluster_rank": cluster_rank,
                    "cluster_size": 12,
                    "clustered_heatmap_order": (cluster_rank - 1) * 12 + offset + 1,
                    "condition_body_site_pair": f"{condition} {offset} | {anatomy}",
                    "condition": condition,
                    "body_site_group": anatomy,
                    "cluster_color": heatmaps.CLUSTER_COLORS[cluster_rank - 1],
                    "cluster_cut_fraction": heatmaps.GROUPED_CLUSTER_CUT_FRACTION,
                    "cluster_cut_distance": 0.7,
                }
            )
    membership = pd.DataFrame.from_records(records)

    enrichment = heatmaps.build_grouped_cluster_enrichment(membership)
    selected = enrichment[enrichment["selected_for_dendrogram_label"]]

    assert set(selected["label_type"]) == {"anatomy", "grouped_condition"}
    assert set(selected["cluster_id"]) == {"C1", "C2"}
    assert selected["fisher_q_value"].le(heatmaps.CLUSTER_LABEL_Q_CUTOFF).all()
    labels = heatmaps.dendrogram_cluster_labels(membership, enrichment)
    assert "Anatomy: Feces" in labels["C1"]
    assert "Condition: Condition A" in labels["C1"]


def test_cluster_comparison_flags_inverse_individual_pattern_as_outlier() -> None:
    grouped_matrix = pd.DataFrame(
        [[1.0, 2.0, 3.0]],
        index=["Group A | Feces"],
        columns=["t1", "t2", "t3"],
    )
    individual_matrix = pd.DataFrame(
        [
            [1.0, 2.0, 3.0],
            [3.0, 2.0, 1.0],
            [2.0, 4.0, 6.0],
            [4.0, 5.0, 8.0],
        ],
        index=["A1 | Feces", "A2 | Feces", "A3 | Feces", "A4 | Feces"],
        columns=["t1", "t2", "t3"],
    )
    grouped_cells = pd.DataFrame(
        {
            "condition_body_site_pair": grouped_matrix.index,
            "condition": ["Group A"],
            "body_site_group": ["Feces"],
        }
    )
    individual_cells = pd.DataFrame(
        {
            "condition_body_site_pair": individual_matrix.index,
            "condition": ["A1", "A2", "A3", "A4"],
            "body_site_group": ["Feces"] * 4,
        }
    )
    grouped = heatmaps.build_spearman_clustering(
        grouped_matrix,
        grouped_cells,
        view_key="grouped",
        view_label="Grouped",
    )
    individual = heatmaps.build_spearman_clustering(
        individual_matrix,
        individual_cells,
        view_key="individual",
        view_label="Individual",
    )
    rollup_map = pd.DataFrame(
        {
            "raw_condition": ["A1", "A2", "A3", "A4"],
            "body_site_group": ["Feces"] * 4,
            "grouped_condition": ["Group A"] * 4,
            "grouped_condition_body_site_pair": ["Group A | Feces"] * 4,
        }
    )

    comparison = heatmaps.compare_grouped_and_individual_clusters(
        grouped,
        individual,
        rollup_map,
    ).set_index("individual_condition")

    assert math.isclose(comparison.loc["A1", "individual_vs_grouped_spearman_rho"], 1.0)
    assert math.isclose(comparison.loc["A2", "individual_vs_grouped_spearman_rho"], -1.0)
    assert bool(comparison.loc["A2", "discordant_or_inverse_pattern"]) is True
    assert bool(comparison.loc["A2", "within_group_iqr_outlier"]) is True
    assert comparison.loc["A2", "global_difference_rank"] == 1


def test_readme_replacement_preserves_following_sections(tmp_path: Path) -> None:
    readme_path = tmp_path / "README.md"
    readme_path.write_text(
        "# Outputs\n\n"
        "## Condition-Body-Site Cytoscape Heatmaps\n\n"
        "old generated text\n\n"
        "## Following Section\n\n"
        "must remain\n",
        encoding="utf-8",
    )
    manifest = pd.DataFrame(
        [
            {
                "view_label": "Pattern",
                "heatmap_type": "trait_pattern",
                "png_path": "pattern.png",
                "clustered_png_path": pd.NA,
                "displayed_cell_count": 1,
                "bh_selected_cell_count": 0,
                "tested_cell_count": 2,
                "row_count": 1,
            }
        ]
    )

    heatmaps.append_readme_section(tmp_path, manifest)

    content = readme_path.read_text(encoding="utf-8")
    assert "old generated text" not in content
    assert "## Following Section\n\nmust remain" in content


def test_heatmap_writer_uses_full_edges_and_keeps_more_than_24_rows(tmp_path: Path) -> None:
    traits = pd.DataFrame(
        {
            "trait_id": ["trait1"],
            "trait_name": ["trait one"],
            "slug": ["trait_one"],
            "base_color": ["#ffffff"],
        }
    )
    traits.to_csv(tmp_path / heatmaps.TRAIT_SELECTION_FILENAME, sep="\t", index=False)
    views = []
    for mode in heatmaps.CONTRAST_MODES:
        for direction, view_key in [
            ("increased", mode.increased_view_key),
            ("decreased", mode.decreased_view_key),
        ]:
            display_filename = f"{view_key}_display_edges.tsv"
            views.append(
                {
                    "view_key": view_key,
                    "view_label": view_key,
                    "focus_direction": direction,
                    "condition_mode": mode.condition_mode,
                    "source_display_edges_path": display_filename,
                }
            )
            edges = pd.DataFrame(
                [
                    {
                        "condition": f"Condition {index}",
                        "body_site_group": "Feces",
                        "trait_id": "trait1",
                        "trait_label": "trait one",
                        "focus_direction": direction,
                        "study_count": 5,
                        "focus_profile_count": 20,
                        "focus_enriched_profile_count": 18 if direction == "increased" else 2,
                        "focus_enriched_study_count": 3,
                        "focus_no_enriched_trait_profile_count": (
                            2 if direction == "increased" else 18
                        ),
                        "comparator_profile_count": 20,
                        "comparator_enriched_profile_count": (
                            2 if direction == "increased" else 18
                        ),
                    }
                    for index in range(30)
                ]
            )
            edges.to_csv(tmp_path / f"{view_key}_edges.tsv", sep="\t", index=False)
            edges.iloc[:1].to_csv(tmp_path / display_filename, sep="\t", index=False)
    pd.DataFrame(views).to_csv(tmp_path / heatmaps.GRAPH_MANIFEST_FILENAME, sep="\t", index=False)

    with (
        patch.object(heatmaps, "plot_matrix_heatmap") as plot,
        patch.object(heatmaps, "plot_clustered_contrast_heatmap") as clustered_plot,
    ):
        manifest = heatmaps.write_six_shared_trait_cytoscape_heatmaps(tmp_path, tmp_path)

    contrasts = manifest[manifest["heatmap_type"] == "signed_direction_contrast"]
    assert contrasts["row_count"].tolist() == [30, 30]
    assert contrasts["displayed_cell_count"].tolist() == [30, 30]
    assert contrasts["tested_cell_count"].tolist() == [30, 30]
    assert contrasts["source_increased_edges_path"].tolist() == [
        "increased_grouped_edges.tsv",
        "increased_individual_conditions_edges.tsv",
    ]
    contrast_calls = [
        call
        for call in plot.call_args_list
        if call.kwargs["title"].startswith("Signed direction contrast:")
    ]
    assert [len(call.args[0]) for call in contrast_calls] == [30, 24, 6, 30, 24, 6]
    for call in contrast_calls:
        matrix = call.args[0]
        assert matrix.notna().all().all()
        assert call.kwargs["vmin"] == -1.0
        assert call.kwargs["vmax"] == 1.0
        assert call.kwargs["significance_borders"].to_numpy(dtype=bool).all()
        assert (np.char.find(call.args[1].to_numpy(dtype=str), "*") == -1).all()
    assert clustered_plot.call_count == 2
    grouped_cluster_call, individual_cluster_call = clustered_plot.call_args_list
    assert not grouped_cluster_call.kwargs["cluster_membership"].empty
    assert not grouped_cluster_call.kwargs["cluster_enrichment"].empty
    assert individual_cluster_call.kwargs["cluster_membership"].empty
    assert contrasts["clustered_png_path"].str.endswith("_clustered_heatmap.png").all()
    table = pd.read_csv(tmp_path / heatmaps.CONTRAST_TABLE_FILENAME, sep="\t")
    assert len(table) == 60
    assert not table.duplicated(["view_key", "condition_body_site_pair", "trait_id"]).any()
    assert table["displayed_in_heatmap"].all()
    pages = pd.read_csv(tmp_path / heatmaps.CONTRAST_PAGES_FILENAME, sep="\t")
    assert pages["row_count"].tolist() == [24, 6, 24, 6]
    assert pages["displayed_cell_count"].sum() == len(table)
    for _view, rows in table.groupby("view_key"):
        assert rows["heatmap_page"].value_counts().to_dict() == {1: 24, 2: 6}
        assert sorted(rows["clustered_heatmap_order"].tolist()) == list(range(1, 31))
    tree_nodes = pd.read_csv(tmp_path / heatmaps.CLUSTER_TREE_NODES_FILENAME, sep="\t")
    tree_leaves = pd.read_csv(tmp_path / heatmaps.CLUSTER_TREE_LEAVES_FILENAME, sep="\t")
    comparison = pd.read_csv(tmp_path / heatmaps.CLUSTER_COMPARISON_FILENAME, sep="\t")
    assert len(tree_nodes) == 58
    assert len(tree_leaves) == 60
    assert len(comparison) == 30
    membership = pd.read_csv(
        tmp_path / heatmaps.GROUPED_CLUSTER_MEMBERSHIP_FILENAME,
        sep="\t",
    )
    enrichment = pd.read_csv(
        tmp_path / heatmaps.GROUPED_CLUSTER_ENRICHMENT_FILENAME,
        sep="\t",
    )
    assert len(membership) == 30
    assert membership["cluster_id"].nunique() == 30
    assert set(enrichment["label_type"]) == {"anatomy", "grouped_condition"}
