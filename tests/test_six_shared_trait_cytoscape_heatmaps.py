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
        flag_column="significant_in_heatmap",
    )
    assert border_mask["trait1"].tolist() == [True, False, False, False]
    assert heatmaps.format_p_value(0.0099) == "<0.01"


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

    with patch.object(heatmaps, "plot_matrix_heatmap") as plot:
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
    assert grouped["spearman_cluster"].notna().all()
    assert grouped["spearman_cluster_label"].notna().all()
    assert (tmp_path / heatmaps.GROUPED_CONTRAST_CLUSTER_ROWS_FILENAME).exists()
    assert (tmp_path / heatmaps.GROUPED_CONTRAST_CLUSTER_SUMMARY_FILENAME).exists()
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
