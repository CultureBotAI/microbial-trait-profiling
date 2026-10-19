import importlib.util
import json
import math
import sys
from itertools import pairwise
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

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


def test_oral_selection_uses_inclusive_or_and_exact_disease_components() -> None:
    cells = pd.DataFrame(
        {
            "view_key": ["individual_conditions_direction_contrast"] * 7
            + ["grouped_direction_contrast"],
            "condition": [
                "Rheumatoid arthritis",
                "Periodontitis,Diabetes mellitus",
                "Root caries",
                "Colorectal cancer",
                "Head and neck carcinoma",
                "Mucositis",
                "Non-periodontitis comparison",
                "Oral / periodontal disease",
            ],
            "body_site_group": ["Oral cavity / dental", "Feces", "Feces"] + ["Feces"] * 5,
            "trait_id": ["trait"] * 8,
            "displayed_in_heatmap": [True] * 8,
            "selected_after_bh": [False, True] + [False] * 6,
            "fisher_p_value": [0.01] * 8,
            "fisher_q_value": [0.4, 0.03] + [0.4] * 6,
        }
    )
    cells["condition_body_site_pair"] = cells["condition"] + " | " + cells["body_site_group"]
    selected = heatmaps.select_oral_contrast_cells(cells)
    assert selected.index.tolist() == [0, 1, 2]
    assert selected["inclusion_reason"].tolist() == ["oral site", "oral disease", "oral disease"]
    pd.testing.assert_frame_equal(selected[cells.columns], cells.iloc[:3])
    with pytest.raises(ValueError, match="unique"):
        heatmaps.select_oral_contrast_cells(pd.concat([cells, cells.iloc[:1]]))
    cells["displayed_in_heatmap"] = "False"
    with pytest.raises(ValueError, match="Boolean"):
        heatmaps.select_oral_contrast_cells(cells)


@pytest.fixture
def oral_contrast_source(tmp_path: Path) -> pd.DataFrame:
    traits = pd.DataFrame(
        {
            "trait_id": [f"trait{i}" for i in range(6)],
            "trait_name": [f"Trait {i}" for i in range(6)],
            "slug": [f"trait{i}" for i in range(6)],
            "base_color": ["#000000"] * 6,
        }
    )
    traits.to_csv(tmp_path / heatmaps.TRAIT_SELECTION_FILENAME, sep="\t", index=False)
    records = []
    for condition, site, supported, q_value in [
        ("Gingivitis", "Oral cavity / dental", False, 1.0),
        ("Rheumatoid arthritis", "Oral cavity / dental", True, 0.4),
        ("Root caries", "Feces", True, 0.03),
        ("Colorectal cancer", "Feces", True, 0.001),
    ]:
        for trait_id in traits["trait_id"]:
            records.append(
                {
                    "view_key": "individual_conditions_direction_contrast",
                    "condition": condition,
                    "body_site_group": site,
                    "condition_body_site_pair": f"{condition} | {site}",
                    "trait_id": trait_id,
                    "study_count": 3 if supported else 1,
                    "passes_support": supported,
                    "passes_nominal_p_value_cutoff": True,
                    "passes_p_value_cutoff": True,
                    "displayed_in_heatmap": supported,
                    "nominal_p_value_selected": supported,
                    "significant_in_heatmap": supported,
                    "selected_after_bh": supported and q_value <= 0.05,
                    "fisher_p_value": 0.001,
                    "fisher_q_value": q_value,
                    "increased_minus_decreased_positive_fraction": 0.25,
                }
            )
    source = pd.DataFrame(records)
    source.loc[6, "study_count"] = 5
    source.to_csv(tmp_path / heatmaps.CONTRAST_TABLE_FILENAME, sep="\t", index=False)
    return source


def test_oral_writer_keeps_source_q_family_and_unsupported_diagnoses(
    tmp_path: Path, oral_contrast_source: pd.DataFrame
) -> None:
    source = oral_contrast_source
    with patch.object(heatmaps, "plot_matrix_heatmap") as plot:
        output = heatmaps.write_oral_trait_heatmaps(tmp_path, source_commit="source-sha")
    assert output.name == "oral_six_trait_heatmap.png"
    first = plot.call_args_list[0]
    matrix = first.args[0]
    assert len(matrix) == 3
    assert matrix.loc[matrix.index.str.startswith("Gingivitis")].isna().all().all()
    borders = first.kwargs["significance_borders"]
    assert borders.loc[borders.index.str.startswith("Root caries")].all().all()
    assert not borders.loc[borders.index.str.startswith("Rheumatoid")].any().any()
    assert "Rheumatoid arthritis | Oral | max n=5" in matrix.index
    subtitle = first.kwargs["subtitle"]
    assert "Dark border: BH-adjusted q <= 0.05 (6 oral cells)" in subtitle
    assert "12 oral cells met nominal Fisher p < 0.01" in subtitle
    assert "|effect|>0.01" in subtitle
    assert "at/below the absolute-effect cutoff" in subtitle
    assert "other rows with displayed cells" in subtitle
    rows = pd.read_csv(output.parent / "oral_six_trait_rows.tsv", sep="\t")
    rheumatoid = rows.loc[rows["condition"] == "Rheumatoid arthritis"]
    assert rheumatoid["max_trait_study_count"].item() == 5
    exported = pd.read_csv(output.parent / "oral_six_trait_cells.tsv", sep="\t")
    assert len(exported) == 18
    assert exported["fisher_q_value"].tolist() == source.iloc[:18]["fisher_q_value"].tolist()
    assert "nominal_p_value_selected" in exported.columns
    assert "significant_in_heatmap" in exported.columns
    assert exported["nominal_p_value_selected"].sum() == 12
    assert exported["significant_in_heatmap"].sum() == 12
    provenance = json.loads((output.parent / "oral_six_trait_provenance.json").read_text())
    assert provenance["source_individual_bh_family_size"] == 24
    assert provenance["source_artifact_commit"] == "source-sha"
    assert provenance["bh_significant_displayed_cells"] == 6
    assert provenance["bh_selected_cell_count"] == 6
    assert provenance["nominal_p_value_cell_count"] == 12
    assert provenance["significant_cell_count"] == 12
    assert provenance["oral_disease_only_rows"] == 1


@pytest.mark.parametrize("bad_pair", [None, "", "Rheumatoid arthritis | Oral cavity / dental"])
def test_oral_selection_rejects_missing_or_conflicting_pair_labels(
    oral_contrast_source: pd.DataFrame, bad_pair: str | None
) -> None:
    oral_contrast_source.loc[0, "condition_body_site_pair"] = bad_pair
    with pytest.raises(ValueError, match="canonical and unambiguous"):
        heatmaps.select_oral_contrast_cells(oral_contrast_source)


def test_oral_selection_rejects_ambiguous_pair_separator(
    oral_contrast_source: pd.DataFrame,
) -> None:
    cells = oral_contrast_source.iloc[:2].copy()
    cells["condition"] = ["A | B", "A"]
    cells["body_site_group"] = ["C", "B | C"]
    cells["condition_body_site_pair"] = "A | B | C"
    with pytest.raises(ValueError, match="canonical and unambiguous"):
        heatmaps.select_oral_contrast_cells(cells)


@pytest.mark.parametrize("bad_trait_id", [None, " ", "trait0"])
@pytest.mark.usefixtures("oral_contrast_source")
def test_oral_writer_rejects_incomplete_or_duplicate_trait_manifest(
    tmp_path: Path, bad_trait_id: str | None
) -> None:
    manifest_path = tmp_path / heatmaps.TRAIT_SELECTION_FILENAME
    manifest = pd.read_csv(manifest_path, sep="\t")
    manifest.loc[5, "trait_id"] = bad_trait_id
    manifest.to_csv(manifest_path, sep="\t", index=False)
    with pytest.raises(ValueError, match="six distinct, nonempty"):
        heatmaps.write_oral_trait_heatmaps(tmp_path)


def test_write_tsv_uses_stable_float_precision(tmp_path: Path) -> None:
    output_path = tmp_path / "floats.tsv"
    table = pd.DataFrame(
        {
            "merge_distance": [0.12345678901234567],
            "signed_minus_log10_fisher_p_value": [-1.2345678901234567],
        }
    )

    heatmaps.write_tsv(table, output_path)

    assert output_path.read_text(encoding="utf-8").splitlines()[1] == (
        "0.123456789012\t-1.23456789012"
    )


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
    assert finalized["nominal_p_value_selected"].tolist() == [True, False, False, False]
    assert finalized["displayed_in_heatmap"].tolist() == [True, True, False, False]
    assert finalized.loc[1, "blank_reason"] == ""
    borders = heatmaps.flag_matrix(
        finalized,
        cells["condition_body_site_pair"].tolist(),
        pd.DataFrame({"trait_id": ["trait1"]}),
        flag_column="selected_after_bh",
    )
    assert borders["trait1"].tolist() == [True, True, False, False]
    annotations = heatmaps.annotation_matrix(
        finalized,
        cells["condition_body_site_pair"].tolist(),
        pd.DataFrame({"trait_id": ["trait1"]}),
        mode="contrast",
    )
    assert annotations["trait1"].tolist() == ["+0.50", "+0.50", "", ""]
    assert heatmaps.format_p_value(0.0099) == "<0.01"


def test_contrast_absolute_effect_cutoff_is_strict_and_auditable() -> None:
    cells = pd.DataFrame(
        {
            "condition_body_site_pair": ["A | Feces", "B | Feces", "C | Saliva"],
            "condition": ["A", "B", "C"],
            "body_site_group": ["Feces", "Feces", "Saliva"],
            "trait_id": ["trait1"] * 3,
            "fisher_p_value": [0.5, 0.5, 0.5],
            "passes_support": [True, True, True],
            "increased_minus_decreased_positive_fraction": [0.01, 0.010001, -0.02],
            "increased_positive_profile_count": [4, 4, 4],
            "decreased_positive_profile_count": [3, 3, 3],
        }
    )

    finalized = heatmaps.finalize_inferential_cells(cells)

    assert finalized["passes_absolute_effect_cutoff"].tolist() == [False, True, True]
    assert finalized["displayed_in_heatmap"].tolist() == [False, True, True]
    assert finalized.loc[0, "blank_reason"] == "absolute effect <=0.01"

    mode = heatmaps.ContrastMode(
        key="grouped",
        label="Grouped conditions",
        condition_mode="grouped condition labels",
        increased_view_key="increased_grouped",
        decreased_view_key="decreased_grouped",
    )
    legend = heatmaps.summarize_contrast_legend(
        finalized,
        ["B | Feces", "C | Saliva"],
        mode,
    )
    assert legend.condition_label_kind == "grouped condition labels"
    assert legend.condition_label_count == 2
    assert legend.anatomy_group_count == 2
    assert legend.row_count == 2
    assert legend.displayed_cell_count == 2
    assert legend.effect_cutoff_hidden_cell_count == 0


def test_spearman_clustering_separates_rows_without_a_defined_neighbor() -> None:
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
    assert sorted(clustering.clustered_row_labels) == sorted(matrix.index[:3])
    assert clustering.insufficient_row_labels == ["D | Feces"]
    leaves = clustering.tree_leaves.set_index("condition_body_site_pair")
    assert bool(leaves.loc["D | Feces", "included_in_clustering"]) is False
    assert bool(leaves.loc["D | Feces", "low_information_for_clustering"]) is True
    assert leaves.loc["D | Feces", "clustering_exclusion_reason"] == (
        "no defined Spearman distance to any other row"
    )
    assert pd.isna(leaves.loc["D | Feces", "linkage_leaf_order"])
    assert leaves.loc["D | Feces", "clustered_heatmap_order"] == len(matrix)
    assert len(clustering.tree_nodes) == len(clustering.clustered_row_labels) - 1
    assert clustering.pair_count == 6
    assert clustering.undefined_pair_count == 3
    assert clustering.undefined_pair_fraction == 0.5
    assert clustering.distinct_distance_count == 3
    assert clustering.distance_tie_excess_count == 3
    assert clustering.distance_tie_fraction == 0.5
    assert clustering.linkage_pair_count == 3
    assert clustering.linkage_undefined_pair_count == 0
    assert clustering.linkage_distance_tie_excess_count == 1
    assert math.isclose(clustering.linkage_distance_tie_fraction, 1 / 3)


def test_spearman_clustering_is_invariant_to_input_row_permutation() -> None:
    matrix = pd.DataFrame(
        [
            [1.0, 2.0, 3.0, 4.0],
            [2.0, 4.0, 6.0, 8.0],
            [4.0, 3.0, 2.0, 1.0],
            [8.0, 6.0, 4.0, 2.0],
            [1.0, math.nan, math.nan, math.nan],
        ],
        index=["B | Feces", "A | Feces", "D | Feces", "C | Feces", "E | Feces"],
        columns=["t1", "t2", "t3", "t4"],
    )
    cells = pd.DataFrame(
        {
            "condition_body_site_pair": matrix.index,
            "condition": ["B", "A", "D", "C", "E"],
            "body_site_group": ["Feces"] * len(matrix),
        }
    )
    expected = heatmaps.build_spearman_clustering(
        matrix,
        cells,
        view_key="test",
        view_label="Test",
    )

    for seed in range(10):
        permuted_matrix = matrix.sample(frac=1.0, random_state=seed)
        permuted = heatmaps.build_spearman_clustering(
            permuted_matrix,
            cells,
            view_key="test",
            view_label="Test",
        )
        assert permuted.clustered_row_labels == expected.clustered_row_labels
        assert permuted.insufficient_row_labels == expected.insufficient_row_labels
        np.testing.assert_allclose(permuted.linkage_matrix, expected.linkage_matrix)


@pytest.mark.parametrize("preview", [False, True])
def test_clustered_heatmap_layout_preserves_row_order_and_readable_labels(
    tmp_path: Path, preview: bool
) -> None:
    conditions = [f"Condition with a long descriptive label {index:02}" for index in range(80)]
    labels = [f"{condition} | Oral cavity / dental" for condition in conditions]
    trait_order = pd.DataFrame(
        {
            "trait_id": [f"t{index}" for index in range(6)],
            "trait_label": [
                "facultatively aerobic",
                "cellulose catabolic process",
                "coccus shaped",
                "sphere shaped",
                "arsenate reduction",
                "moderately halophilic",
            ],
        }
    )
    matrix = pd.DataFrame(
        [np.roll(np.linspace(-1.0, 1.0, 6), index % 6) for index in range(len(labels))],
        index=labels,
        columns=trait_order["trait_id"],
    )
    cells = pd.DataFrame(
        {
            "condition_body_site_pair": labels,
            "condition": conditions,
            "body_site_group": ["Oral cavity / dental"] * len(labels),
        }
    )
    clustering = heatmaps.build_spearman_clustering(
        matrix, cells, view_key="grouped", view_label="Grouped"
    )
    borders = pd.DataFrame(False, index=matrix.index, columns=matrix.columns)
    borders.iloc[0, 0] = True
    annotations = matrix.map(lambda value: f"{value:.2f}")
    legend = heatmaps.ContrastLegendSummary(
        condition_label_kind="grouped condition labels",
        condition_label_count=len(labels),
        anatomy_group_count=1,
        row_count=len(labels),
        displayed_cell_count=matrix.size,
        displayed_bh_selected_cell_count=0,
        effect_cutoff_hidden_cell_count=0,
    )

    with (
        patch.object(heatmaps.Figure, "savefig", autospec=True) as savefig,
        patch.object(heatmaps.plt, "close"),
    ):
        heatmaps.plot_clustered_contrast_heatmap(
            clustering,
            borders,
            trait_order,
            tmp_path / "clustered.png",
            title="Clustered signed direction contrast",
            subtitle="Wrapped condition and body-site labels remain readable.",
            contrast_legend=legend,
            annotations=annotations,
            preview=preview,
        )
    figure = savefig.call_args.args[0]
    try:
        figure.canvas.draw()
        heatmap_axis = figure.axes[1]
        renderer = figure.canvas.get_renderer()
        legend_bounds = figure.legends[0].get_window_extent(renderer)
        colorbar_bounds = figure.axes[2].get_tightbbox(renderer)
        assert colorbar_bounds is not None
        assert legend_bounds.y1 < colorbar_bounds.y0
        assert colorbar_bounds.y1 < heatmap_axis.xaxis.label.get_window_extent(renderer).y0
        ordered_labels = clustering.clustered_row_labels + clustering.insufficient_row_labels
        np.testing.assert_allclose(
            heatmap_axis.images[0].get_array(), matrix.loc[ordered_labels].to_numpy()
        )
        assert len(heatmap_axis.patches) == 1
        tick_labels = heatmap_axis.get_yticklabels()
        if preview:
            assert not tick_labels
            assert not heatmap_axis.texts
            assert figure.get_figheight() == 8.5
        else:
            assert len(tick_labels) == len(labels)
            assert all(label.get_fontsize() >= 7.0 for label in tick_labels)
            assert all("\n" in label.get_text() for label in tick_labels)
            bounds = [label.get_window_extent() for label in tick_labels]
            assert all(left.y1 < right.y0 for left, right in pairwise(bounds))
            assert len(heatmap_axis.texts) == matrix.size
            assert [text.get_text() for text in heatmap_axis.texts] == (
                annotations.loc[ordered_labels].to_numpy().ravel().tolist()
            )
    finally:
        heatmaps.plt.close(figure)


@pytest.mark.parametrize("eligible_count", [3, 4])
@pytest.mark.parametrize("low_overlap_direction", [1.0, -1.0])
def test_cluster_comparison_uses_only_eligible_rows_for_outlier_threshold(
    eligible_count: int,
    low_overlap_direction: float,
) -> None:
    grouped_matrix = pd.DataFrame(
        [[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]],
        index=["Group A | Feces"],
        columns=["t1", "t2", "t3", "t4", "t5", "t6"],
    )
    individual_matrix = pd.DataFrame(
        [
            [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            [6.0, 5.0, 4.0, 3.0, 2.0, 1.0],
            [2.0, 4.0, 6.0, 8.0, 10.0, 12.0],
            [4.0, 5.0, 8.0, 9.0, 12.0, 13.0],
            [
                low_overlap_direction,
                2.0 * low_overlap_direction,
                3.0 * low_overlap_direction,
                math.nan,
                math.nan,
                math.nan,
            ],
        ],
        index=[
            "A1 | Feces",
            "A2 | Feces",
            "A3 | Feces",
            "A4 | Feces",
            "A5 | Feces",
        ],
        columns=["t1", "t2", "t3", "t4", "t5", "t6"],
    )
    if eligible_count == 3:
        individual_matrix = individual_matrix.drop(index="A4 | Feces")
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
            "condition": individual_matrix.index.str.removesuffix(" | Feces"),
            "body_site_group": ["Feces"] * len(individual_matrix),
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
    # The aggregate includes an undisplayed raw condition and repeated source
    # labels across trait edges. Literal commas belong to a single raw label.
    rollup_map = heatmaps.build_condition_rollup_map(
        pd.DataFrame(
            {
                "raw_conditions": ["A1; A2; A3; A4; A5; A6, with comma"] * 2,
                "body_site_group": ["Feces"] * 2,
                "condition": ["Group A"] * 2,
            }
        )
    )

    comparison = heatmaps.compare_grouped_and_individual_clusters(
        grouped,
        individual,
        rollup_map,
    ).set_index("individual_condition")

    assert math.isclose(comparison.loc["A1", "individual_vs_grouped_spearman_rho"], 1.0)
    assert math.isclose(comparison.loc["A2", "individual_vs_grouped_spearman_rho"], -1.0)
    assert bool(comparison.loc["A2", "pattern_diagnostic_eligible"]) is True
    assert bool(comparison.loc["A2", "discordant_or_inverse_pattern"]) is True
    assert bool(comparison.loc["A2", "within_group_iqr_outlier"]) is (eligible_count == 4)
    if eligible_count == 4:
        assert math.isclose(comparison.loc["A2", "within_group_outlier_threshold"], 1.25)
    else:
        assert comparison["within_group_outlier_threshold"].isna().all()
    assert comparison.loc["A2", "global_difference_rank"] == 1
    assert comparison["grouped_contributing_individual_condition_count"].eq(6).all()
    assert comparison["grouped_displayed_individual_condition_count"].eq(eligible_count + 1).all()
    assert comparison.loc["A5", "shared_supported_trait_count"] == 3
    assert bool(comparison.loc["A5", "pattern_diagnostic_eligible"]) is False
    assert bool(comparison.loc["A5", "discordant_or_inverse_pattern"]) is False
    assert pd.isna(comparison.loc["A5", "global_difference_rank"])
    assert pd.isna(comparison.loc["A5", "within_group_difference_rank"])
    assert bool(comparison.loc["A5", "within_group_iqr_outlier"]) is False


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
                "view_key": "pattern",
                "view_label": "Pattern",
                "heatmap_type": "trait_pattern",
                "png_path": "pattern.png",
                "clustered_png_path": pd.NA,
                "displayed_cell_count": 1,
                "significant_cell_count": 1,
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
    assert "`low positive-z profile or study support`" in content
    assert "absolute effect <=0.01" in content
    assert "neither direction emitted a candidate edge is absent" in content
    assert "fixed upstream labels" in content
    assert "Issue #256 / PR #257" in content


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
    for filename in heatmaps.OBSOLETE_CLUSTER_OUTPUTS:
        (tmp_path / filename).write_text("obsolete\n", encoding="utf-8")

    with (
        patch.object(heatmaps, "plot_matrix_heatmap") as plot,
        patch.object(heatmaps, "plot_clustered_contrast_heatmap") as clustered_plot,
    ):
        manifest = heatmaps.write_six_shared_trait_cytoscape_heatmaps(tmp_path, tmp_path)

    contrasts = manifest[manifest["heatmap_type"] == "signed_direction_contrast"]
    assert contrasts["row_count"].tolist() == [30, 30]
    assert contrasts["displayed_cell_count"].tolist() == [30, 30]
    assert contrasts["condition_label_count"].tolist() == [30, 30]
    assert contrasts["anatomy_group_count"].tolist() == [1, 1]
    assert contrasts["absolute_effect_cutoff"].tolist() == [
        heatmaps.SIGNED_EFFECT_ABS_CUTOFF,
        heatmaps.SIGNED_EFFECT_ABS_CUTOFF,
    ]
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
    assert [len(call.args[0]) for call in contrast_calls] == [24, 6, 30, 24, 6]
    for call in contrast_calls:
        matrix = call.args[0]
        assert matrix.notna().all().all()
        assert call.kwargs["vmin"] == -1.0
        assert call.kwargs["vmax"] == 1.0
        assert call.kwargs["significance_borders"].to_numpy(dtype=bool).all()
        assert call.kwargs["contrast_legend"].condition_label_count == len(matrix)
        assert (np.char.find(call.args[1].to_numpy(dtype=str), "*") == -1).all()
    assert clustered_plot.call_count == 4
    assert contrasts["clustered_png_path"].str.endswith("_clustered_heatmap.png").all()
    assert contrasts["primary_row_order"].tolist() == ["clustered", "ranked"]
    grouped_manifest = contrasts.iloc[0]
    assert grouped_manifest["preview_png_path"].endswith("_heatmap_preview.png")
    grouped_primary = clustered_plot.call_args_list[1]
    grouped_preview = clustered_plot.call_args_list[2]
    assert grouped_primary.args[3] == tmp_path / grouped_manifest["png_path"]
    assert grouped_preview.args[3] == tmp_path / grouped_manifest["preview_png_path"]
    assert not grouped_primary.kwargs["preview"]
    assert grouped_preview.kwargs["preview"]
    assert grouped_primary.kwargs["annotations"].notna().all().all()
    grouped_page_labels = [label for call in contrast_calls[:2] for label in call.args[0].index]
    assert grouped_page_labels == grouped_primary.args[0].insufficient_row_labels
    assert contrasts["clustered_row_count"].tolist() == [0, 0]
    assert contrasts["insufficient_overlap_row_count"].tolist() == [30, 30]
    assert contrasts["undefined_spearman_pair_count"].tolist() == [435, 435]
    assert contrasts["undefined_spearman_pair_fraction"].tolist() == [1.0, 1.0]
    assert contrasts["distance_tie_excess_count"].tolist() == [434, 434]
    assert contrasts["distance_tie_fraction"].tolist() == [434 / 435, 434 / 435]
    assert contrasts["linkage_pair_count"].tolist() == [0, 0]
    assert contrasts["linkage_distance_tie_fraction"].tolist() == [0.0, 0.0]
    assert not contrasts["dendrogram_cut_applied"].any()
    table = pd.read_csv(tmp_path / heatmaps.CONTRAST_TABLE_FILENAME, sep="\t")
    assert len(table) == 60
    assert not table.duplicated(["view_key", "condition_body_site_pair", "trait_id"]).any()
    assert table["displayed_in_heatmap"].all()
    pages = pd.read_csv(tmp_path / heatmaps.CONTRAST_PAGES_FILENAME, sep="\t")
    assert pages["row_count"].tolist() == [24, 6, 24, 6]
    assert pages["condition_label_count"].tolist() == [24, 6, 24, 6]
    assert pages["anatomy_group_count"].tolist() == [1, 1, 1, 1]
    assert pages["displayed_cell_count"].sum() == len(table)
    for _view, rows in table.groupby("view_key"):
        assert rows["heatmap_page"].value_counts().to_dict() == {1: 24, 2: 6}
        assert sorted(rows["clustered_heatmap_order"].tolist()) == list(range(1, 31))
    tree_nodes = pd.read_csv(tmp_path / heatmaps.CLUSTER_TREE_NODES_FILENAME, sep="\t")
    tree_leaves = pd.read_csv(tmp_path / heatmaps.CLUSTER_TREE_LEAVES_FILENAME, sep="\t")
    comparison = pd.read_csv(tmp_path / heatmaps.CLUSTER_COMPARISON_FILENAME, sep="\t")
    assert len(tree_nodes) == 0
    assert len(tree_leaves) == 60
    assert not tree_leaves["included_in_clustering"].any()
    assert len(comparison) == 30
    for filename in heatmaps.OBSOLETE_CLUSTER_OUTPUTS:
        assert not (tmp_path / filename).exists()


@pytest.mark.parametrize("supported_modes", [set(), {"grouped"}, {"individual_conditions"}])
def test_heatmap_writer_preserves_audit_schemas_when_a_view_has_no_displayed_rows(
    tmp_path: Path, supported_modes: set[str]
) -> None:
    traits = pd.DataFrame(
        {
            "trait_id": [f"trait{index}" for index in range(6)],
            "trait_name": [f"Trait {index}" for index in range(6)],
            "slug": [f"trait_{index}" for index in range(6)],
            "base_color": ["#ffffff"] * 6,
        }
    )
    traits.to_csv(tmp_path / heatmaps.TRAIT_SELECTION_FILENAME, sep="\t", index=False)
    views = []
    for mode in heatmaps.CONTRAST_MODES:
        supported = mode.key in supported_modes
        increased_positive, decreased_positive = (18, 2) if supported else (2, 1)
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
            focus_positive, comparator_positive = (
                (increased_positive, decreased_positive)
                if direction == "increased"
                else (decreased_positive, increased_positive)
            )
            edges = pd.DataFrame(
                [
                    {
                        "condition": "Group A" if mode.key == "grouped" else "Individual A",
                        "raw_conditions": "Individual A",
                        "body_site_group": "Feces",
                        "trait_id": trait.trait_id,
                        "trait_label": trait.trait_name,
                        "focus_direction": direction,
                        "study_count": 3,
                        "focus_profile_count": 20,
                        "focus_enriched_profile_count": focus_positive,
                        "focus_enriched_study_count": 3 if supported else 1,
                        "focus_no_enriched_trait_profile_count": 20 - focus_positive,
                        "comparator_profile_count": 20,
                        "comparator_enriched_profile_count": comparator_positive,
                    }
                    for trait in traits.itertuples(index=False)
                ]
            )
            edges.to_csv(tmp_path / f"{view_key}_edges.tsv", sep="\t", index=False)
            edges.to_csv(tmp_path / display_filename, sep="\t", index=False)
    pd.DataFrame(views).to_csv(tmp_path / heatmaps.GRAPH_MANIFEST_FILENAME, sep="\t", index=False)

    with (
        patch.object(heatmaps, "plot_matrix_heatmap"),
        patch.object(heatmaps, "plot_clustered_contrast_heatmap"),
    ):
        manifest = heatmaps.write_six_shared_trait_cytoscape_heatmaps(tmp_path, tmp_path)

    contrast_manifest = manifest.loc[manifest["heatmap_type"].eq("signed_direction_contrast")]
    assert contrast_manifest["row_count"].tolist() == [
        int(mode.key in supported_modes) for mode in heatmaps.CONTRAST_MODES
    ]
    contrasts = pd.read_csv(tmp_path / heatmaps.CONTRAST_TABLE_FILENAME, sep="\t")
    assert len(contrasts) == 12
    assert contrasts["displayed_in_heatmap"].sum() == 6 * len(supported_modes)
    leaves = pd.read_csv(tmp_path / heatmaps.CLUSTER_TREE_LEAVES_FILENAME, sep="\t")
    nodes = pd.read_csv(tmp_path / heatmaps.CLUSTER_TREE_NODES_FILENAME, sep="\t")
    comparison = pd.read_csv(tmp_path / heatmaps.CLUSTER_COMPARISON_FILENAME, sep="\t")
    assert len(leaves) == len(supported_modes)
    assert set(leaves["view_key"]) == {f"{mode}_direction_contrast" for mode in supported_modes}
    assert not leaves["included_in_clustering"].any()
    assert nodes.empty
    assert {"node_id", "left_child_id", "right_child_id"} <= set(nodes.columns)
    assert len(comparison) == int("individual_conditions" in supported_modes)
    assert comparison["individual_vs_grouped_spearman_rho"].isna().all()
    assert not comparison["pattern_diagnostic_eligible"].any()
    assert comparison["within_group_outlier_threshold"].isna().all()
    assert comparison["global_difference_rank"].isna().all()
    if "individual_conditions" in supported_modes:
        assert comparison["comparison_status"].tolist() == [
            "grouped row did not meet display support"
        ]


def test_oral_and_full_contrast_nominal_p_versus_bh_border_regression() -> None:
    source_path = heatmaps.OUTPUT_DIR / heatmaps.CONTRAST_TABLE_FILENAME
    assert source_path.exists(), f"Expected contrast table at {source_path}"
    contrasts = pd.read_csv(source_path, sep="\t")

    individual = contrasts[contrasts["view_key"] == "individual_conditions_direction_contrast"]
    assert len(individual) == 3144
    assert individual["passes_support"].sum() == 663
    assert individual["nominal_p_value_selected"].sum() == 29
    assert individual["significant_in_heatmap"].sum() == 29
    assert individual["selected_after_bh"].sum() == 1
    assert individual["displayed_in_heatmap"].sum() == 632
    assert (individual["nominal_p_value_selected"] == individual["significant_in_heatmap"]).all()

    # Oral cavity / dental subset within individual view
    oral_site = individual[individual["body_site_group"] == "Oral cavity / dental"]
    assert oral_site["passes_support"].sum() == 87
    assert oral_site["nominal_p_value_selected"].sum() == 2
    assert oral_site["significant_in_heatmap"].sum() == 2
    assert oral_site["selected_after_bh"].sum() == 0
    assert oral_site["displayed_in_heatmap"].sum() == 83

    # Oral export subset via select_oral_contrast_cells
    oral_export = heatmaps.select_oral_contrast_cells(contrasts)
    assert len(oral_export) == 522
    assert oral_export["condition_body_site_pair"].nunique() == 117
    assert oral_export["passes_support"].sum() == 87
    assert oral_export["nominal_p_value_selected"].sum() == 2
    assert oral_export["significant_in_heatmap"].sum() == 2
    assert oral_export["selected_after_bh"].sum() == 0
    assert oral_export["displayed_in_heatmap"].sum() == 83
    assert (oral_export["nominal_p_value_selected"] == oral_export["significant_in_heatmap"]).all()

    # Squamous cell carcinoma and oral squamous cell carcinoma stay separate oral rows.
    scc_pair = "Squamous cell carcinoma | Oral cavity / dental"
    oscc_pair = "Oral squamous cell carcinoma | Oral cavity / dental"
    for frame in (individual, oral_export):
        assert {scc_pair, oscc_pair} <= set(frame["condition_body_site_pair"])
    pair_conditions = individual.groupby("condition_body_site_pair")["condition"].unique()
    assert pair_conditions[scc_pair].tolist() == ["Squamous cell carcinoma"]
    assert pair_conditions[oscc_pair].tolist() == ["Oral squamous cell carcinoma"]

    # Verified regression cell: Colorectal cancer | Oral cavity / dental x coccus shaped
    # (METPO:1000668)
    crc_full = individual[
        (individual["condition"] == "Colorectal cancer")
        & (individual["body_site_group"] == "Oral cavity / dental")
        & (individual["trait_id"] == "METPO:1000668")
    ]
    assert len(crc_full) == 1
    cell = crc_full.iloc[0]
    assert math.isclose(cell["increased_minus_decreased_positive_fraction"], 5 / 7, rel_tol=1e-9)
    assert math.isclose(cell["fisher_p_value"], 0.00339366515837, rel_tol=1e-6)
    assert cell["fisher_p_value"] < heatmaps.P_VALUE_CUTOFF
    assert math.isclose(cell["fisher_q_value"], 0.587927267704, rel_tol=1e-6)
    assert cell["fisher_q_value"] > heatmaps.Q_VALUE_CUTOFF
    assert bool(cell["passes_support"]) is True
    assert bool(cell["nominal_p_value_selected"]) is True
    assert bool(cell["significant_in_heatmap"]) is True
    assert bool(cell["selected_after_bh"]) is False
    assert bool(cell["displayed_in_heatmap"]) is True

    # Verify that in flag_matrix for borders (selected_after_bh),
    # the nominal cell has no border in full view and oral export
    traits = pd.DataFrame({"trait_id": ["METPO:1000668"]})
    full_borders = heatmaps.flag_matrix(
        individual,
        ["Colorectal cancer | Oral cavity / dental"],
        traits,
        flag_column="selected_after_bh",
    )
    assert not full_borders.loc["Colorectal cancer | Oral cavity / dental", "METPO:1000668"]

    oral_borders = heatmaps.flag_matrix(
        oral_export,
        ["Colorectal cancer | Oral cavity / dental"],
        traits,
        flag_column="selected_after_bh",
    )
    assert not oral_borders.loc["Colorectal cancer | Oral cavity / dental", "METPO:1000668"]
