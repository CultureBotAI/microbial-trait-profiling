import importlib.util
import math
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

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
    border_mask = heatmaps.flag_matrix(
        finalized,
        cells["condition_body_site_pair"].tolist(),
        pd.DataFrame({"trait_id": ["trait1"]}),
        flag_column="selected_after_bh",
    )
    assert border_mask["trait1"].tolist() == [True, True, False, False]
    assert heatmaps.format_p_value(0.0099) == "<0.01"


def test_signed_effect_cutoff_is_strictly_greater_than_point_zero_one() -> None:
    cells = pd.DataFrame(
        {
            "condition_body_site_pair": ["A | Feces", "B | Feces", "C | Feces"],
            "body_site_group": ["Feces"] * 3,
            "trait_id": ["trait1"] * 3,
            "fisher_p_value": [0.5, 0.5, 0.5],
            "passes_support": [True, True, True],
            "increased_minus_decreased_positive_fraction": [0.01, 0.010001, -0.02],
            "increased_positive_profile_count": [10, 10, 10],
            "decreased_positive_profile_count": [1, 1, 1],
        }
    )

    finalized = heatmaps.finalize_inferential_cells(cells)

    assert finalized["passes_absolute_effect_cutoff"].tolist() == [False, True, True]
    assert finalized["displayed_in_heatmap"].tolist() == [False, True, True]
    assert finalized.loc[0, "blank_reason"] == "absolute effect <=0.01"


def test_compact_dendrogram_label_keeps_each_supported_pattern_readable() -> None:
    label = (
        "Pattern: Increased-associated coccus shaped "
        "(23/23 rows; q=0.018; 1 BH cell); Increased-associated sphere shaped "
        "(24/27 rows; q=0.048; 2 BH cells)"
    )

    compact = heatmaps.compact_dendrogram_label(label, width=30)

    assert compact.splitlines() == [
        "Increased: coccus shaped",
        "23/23 | q=0.018 | 1 BH",
        "Increased: sphere shaped",
        "24/27 | q=0.048 | 2 BH",
    ]


def test_contrast_legend_reports_visible_scope_and_significance() -> None:
    mode = heatmaps.CONTRAST_MODES[0]
    cells = pd.DataFrame(
        {
            "condition_body_site_pair": ["A | Feces", "A | Feces", "B | Oral"],
            "condition": ["A", "A", "B"],
            "body_site_group": ["Feces", "Feces", "Oral"],
            "displayed_in_heatmap": [True, False, True],
            "selected_after_bh": [True, False, False],
            "nominal_p_value_selected": [True, False, False],
            "passes_support": [True, True, True],
            "fisher_p_value": [0.001, 0.5, 0.2],
            "passes_absolute_effect_cutoff": [True, False, True],
        }
    )

    summary = heatmaps.summarize_contrast_legend(
        cells,
        ["A | Feces", "B | Oral"],
        mode,
    )

    assert summary.condition_label_count == 2
    assert summary.anatomy_group_count == 2
    assert summary.row_count == 2
    assert summary.displayed_cell_count == 2
    assert summary.nominal_p_value_cell_count == 1
    assert summary.displayed_bh_selected_cell_count == 1
    assert summary.effect_cutoff_hidden_cell_count == 1


def test_dendrogram_labels_all_bh_significant_enrichments() -> None:
    row_labels = [f"F{index} | Feces" for index in range(5)] + [
        f"O{index} | Oral" for index in range(15)
    ]
    increasing = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    decreasing = list(reversed(increasing))
    matrix = pd.DataFrame(
        [increasing] * 5 + [decreasing] * 15,
        index=row_labels,
        columns=[f"trait{index}" for index in range(6)],
    )
    cells = pd.DataFrame(
        {
            "condition_body_site_pair": row_labels,
            "condition": [f"Condition {index % 2}" for index in range(20)],
            "body_site_group": ["Feces"] * 5 + ["Oral"] * 15,
            "trait_id": ["trait0"] * 20,
            "trait_label": ["trait one"] * 20,
            "displayed_in_heatmap": [True] * 20,
            "significant_in_heatmap": [False] * 20,
            "nominal_p_value_selected": [False] * 20,
            "selected_after_bh": [False] * 20,
            "increased_positive_profile_count": [5] * 20,
            "decreased_positive_profile_count": [1] * 20,
            "increased_minus_decreased_positive_fraction": [0.4] * 20,
            "plot_signed_minus_log10_fisher_p_value": [1.0] * 20,
        }
    )

    rows, summary, branch_tests, pattern_tests, linkage_matrix = heatmaps.clustered_contrast_rows(
        matrix, cells
    )

    assert linkage_matrix is not None
    assert not summary.empty
    assert summary["cluster_label_supported"].all()
    assert summary["minimum_branch_label_q_value"].le(heatmaps.Q_VALUE_CUTOFF).all()
    assert (
        summary["maximum_branch_label_purity"]
        .ge(heatmaps.MIN_SIGNIFICANT_BRANCH_LABEL_PURITY)
        .all()
    )
    assert summary["cluster_label"].str.contains("Anatomy:").all()
    assert not summary["cluster_label"].str.contains("Condition:").any()
    assert rows.loc[rows["spearman_cluster"].notna(), "cluster_label"].notna().all()
    assert rows.loc[rows["spearman_cluster"].isna(), "cluster_label"].isna().all()
    shown_tests = branch_tests["included_in_branch_label"]
    assert shown_tests.any()
    assert branch_tests.loc[shown_tests, "passes_display_criteria"].all()
    assert not pattern_tests["included_in_branch_label"].any()


def test_dendrogram_uses_descriptive_pattern_when_metadata_is_not_significant() -> None:
    row_labels = [f"R{index} | Feces" for index in range(20)]
    increasing = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    decreasing = list(reversed(increasing))
    matrix = pd.DataFrame(
        [increasing] * 10 + [decreasing] * 10,
        index=row_labels,
        columns=[f"trait{index}" for index in range(6)],
    )
    cells = pd.DataFrame(
        {
            "condition_body_site_pair": row_labels,
            "condition": ["Condition A"] * 20,
            "body_site_group": ["Feces"] * 20,
            "trait_id": ["trait0"] * 20,
            "trait_label": ["facultatively aerobic"] * 20,
            "displayed_in_heatmap": [True] * 20,
            "significant_in_heatmap": [True] + [False] * 9 + [True] + [False] * 9,
            "nominal_p_value_selected": [True] + [False] * 9 + [True] + [False] * 9,
            "selected_after_bh": [True] + [False] * 9 + [True] + [False] * 9,
            "increased_positive_profile_count": [5] * 20,
            "decreased_positive_profile_count": [1] * 20,
            "increased_minus_decreased_positive_fraction": [0.4] * 10 + [-0.4] * 10,
            "plot_signed_minus_log10_fisher_p_value": [1.0] * 20,
        }
    )

    _rows, summary, branch_tests, pattern_tests, linkage_matrix = heatmaps.clustered_contrast_rows(
        matrix, cells
    )

    assert linkage_matrix is not None
    assert not branch_tests["included_in_branch_label"].any()
    assert not summary.empty
    assert summary["cluster_label_source"].eq("bh_supported_trait_pattern").all()
    assert summary["cluster_label"].str.startswith("Pattern:").all()
    assert summary["cluster_label"].str.contains("q=").all()
    assert summary["cluster_label"].str.contains("BH cell").all()
    included_patterns = pattern_tests["included_in_branch_label"]
    assert included_patterns.any()
    assert pattern_tests.loc[included_patterns, "passes_bh"].all()
    assert pattern_tests.loc[included_patterns, "passes_display_criteria"].all()
    assert pattern_tests.loc[
        pattern_tests["passes_display_criteria"], "included_in_branch_label"
    ].all()


def test_heatmap_writer_uses_full_edges_and_keeps_more_than_24_rows(tmp_path: Path) -> None:
    (tmp_path / "README.md").write_text("# Fixture\n", encoding="utf-8")
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

    with patch.object(heatmaps, "plot_matrix_heatmap") as plot:
        manifest = heatmaps.write_six_shared_trait_cytoscape_heatmaps(tmp_path, tmp_path)

    contrasts = manifest[manifest["heatmap_type"] == "signed_direction_contrast"]
    assert contrasts["row_count"].tolist() == [30, 30]
    assert contrasts["condition_label_count"].tolist() == [30, 30]
    assert contrasts["anatomy_group_count"].tolist() == [1, 1]
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
    assert [len(call.args[0]) for call in contrast_calls] == [30, 24, 6]
    for call in contrast_calls:
        matrix = call.args[0]
        assert matrix.notna().all().all()
        assert call.kwargs["vmin"] == -1.0
        assert call.kwargs["vmax"] == 1.0
        assert call.kwargs["border_mask"].any().all()
    table = pd.read_csv(tmp_path / heatmaps.CONTRAST_TABLE_FILENAME, sep="\t")
    assert len(table) == 60
    assert not table.duplicated(["view_key", "condition_body_site_pair", "trait_id"]).any()
    assert table["displayed_in_heatmap"].all()
    grouped = table[table["view_key"].eq("grouped_direction_contrast")]
    assert grouped["spearman_cluster"].isna().all()
    assert (tmp_path / heatmaps.GROUPED_CONTRAST_CLUSTER_ROWS_FILENAME).exists()
    assert (tmp_path / heatmaps.GROUPED_CONTRAST_CLUSTER_SUMMARY_FILENAME).exists()
    assert (tmp_path / heatmaps.GROUPED_CONTRAST_BRANCH_TESTS_FILENAME).exists()
    assert (tmp_path / heatmaps.GROUPED_CONTRAST_PATTERN_TESTS_FILENAME).exists()
    pages = pd.read_csv(tmp_path / heatmaps.CONTRAST_PAGES_FILENAME, sep="\t")
    assert pages["view_key"].tolist() == [
        "individual_conditions_direction_contrast",
        "individual_conditions_direction_contrast",
    ]
    assert pages["row_count"].tolist() == [24, 6]
    assert pages["displayed_cell_count"].sum() == len(
        table[table["view_key"].eq("individual_conditions_direction_contrast")]
    )
    individual = table[table["view_key"].eq("individual_conditions_direction_contrast")]
    assert individual["heatmap_page"].value_counts().to_dict() == {1: 24, 2: 6}
    assert table[table["view_key"].eq("grouped_direction_contrast")]["heatmap_page"].isna().all()
    readme = (tmp_path / "README.md").read_text(encoding="utf-8")
    assert "Every internal dendrogram branch is tested" in readme
    assert "absolute magnitude is greater than 0.01" in readme
    assert "Issue #256 / PR #257" in readme
