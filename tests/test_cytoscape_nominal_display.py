"""Nominal Fisher labels and distinct finite/displayed counts for issue #1430."""

import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import matplotlib.figure
import numpy as np
import pandas as pd
import pytest

SOURCE = (
    Path(__file__).resolve().parents[1] / "scripts/create_six_shared_trait_cytoscape_heatmaps.py"
)
SPEC = importlib.util.spec_from_file_location("cytoscape_nominal_display", SOURCE)
assert SPEC is not None
assert SPEC.loader is not None
heatmaps = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = heatmaps
SPEC.loader.exec_module(heatmaps)


@pytest.fixture
def display_inputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    traits = pd.DataFrame(
        {"trait_id": [f"t{i}" for i in range(6)], "trait_label": [f"Trait {i}" for i in range(6)]}
    )
    # Candidate=9, finite=7, displayed flag=7, nominal hits=2, visible finite=3.
    # The extra trait, capped row, missing p/value, and unsupported finite p each
    # distinguish source-view tests from entries represented by the actual matrix.
    records = [
        ("A", "t0", 0.001, True, 2.0),
        ("A", "t1", np.nan, True, 2.0),
        ("A", "t2", 0.4, False, 2.0),
        ("A", "t3", 0.002, True, np.nan),
        ("A", "outside", 0.003, True, 2.0),
        ("B", "t0", 0.004, True, 2.0),
        ("B", "t1", 0.006, True, 2.0),
        ("C", "t0", 0.005, True, 2.0),
        ("Z", "t0", np.inf, False, 2.0),
    ]
    cells = pd.DataFrame(
        records,
        columns=[
            "condition",
            "trait_id",
            "fisher_p_value",
            "displayed_in_heatmap",
            "plot_minus_log10_fisher_p_value",
        ],
    )
    cells["condition_body_site_pair"] = cells["condition"]
    cells["body_site_group"] = "Feces"
    cells["increased_minus_decreased_positive_fraction"] = (
        cells["plot_minus_log10_fisher_p_value"] / 10
    )
    cells["focus_positive_profile_count"] = 4
    cells["nominal_p_value_selected"] = cells.index.isin([0, 5])
    cells["significant_in_heatmap"] = cells["nominal_p_value_selected"]
    cells["selected_after_bh"] = cells.index == 0
    cells["passes_support"] = True
    cells["passes_absolute_effect_cutoff"] = True
    cells["fisher_q_value"] = 0.2
    cells["view_key"] = "fixture"
    cells["trait_rank"] = range(len(cells))
    return cells, traits


def test_nominal_blank_reason_and_annotation_keep_values() -> None:
    source = pd.DataFrame({"fisher_p_value": [0.2, np.nan, 0.001], "passes_support": [True] * 3})
    result = heatmaps.finalize_inferential_cells(source)
    assert result["blank_reason"].tolist() == ["nominal Fisher p>=0.01", "Fisher p=NA", ""]
    pd.testing.assert_series_equal(result["fisher_p_value"], source["fisher_p_value"])
    assert result["fisher_q_value"].tolist() == pytest.approx([0.3, 1.0, 0.003])
    assert result["displayed_in_heatmap"].tolist() == [False, False, True]
    cells = pd.DataFrame(
        {
            "condition_body_site_pair": ["A", "B"],
            "trait_id": ["t0"] * 2,
            "fisher_p_value": [0.001, 1e-20],
            "displayed_in_heatmap": [True] * 2,
            "focus_positive_profile_count": [4, 5],
        }
    )
    labels = heatmaps.annotation_matrix(
        cells, ["A", "B"], pd.DataFrame({"trait_id": ["t0"]}), mode="pattern"
    )
    assert labels.loc["A", "t0"].startswith("n=4\nnominal p=")
    assert labels.loc["B", "t0"].startswith("n=5\nnominal p<")


@pytest.mark.parametrize(
    "value_column",
    ["plot_minus_log10_fisher_p_value", "increased_minus_decreased_positive_fraction"],
)
def test_finite_counts_use_the_actual_row_trait_and_value_intersection(
    display_inputs: tuple[pd.DataFrame, pd.DataFrame], value_column: str
) -> None:
    cells, traits = display_inputs
    before = cells.copy(deep=True)
    assert heatmaps.heatmap_test_counts(cells, ["A", "B"], traits, value_column=value_column) == {
        "finite_test_count": 7,
        "displayed_finite_test_count": 3,
    }
    assert (
        heatmaps.heatmap_test_counts(cells, ["A"], traits, value_column=value_column)[
            "displayed_finite_test_count"
        ]
        == 1
    )
    assert (
        heatmaps.heatmap_test_counts(cells, [], traits, value_column=value_column)[
            "displayed_finite_test_count"
        ]
        == 0
    )
    pd.testing.assert_frame_equal(cells, before)


def test_duplicate_matrix_keys_count_source_tests_but_only_one_displayed_cell(
    display_inputs: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    cells, traits = display_inputs
    duplicated = pd.concat([cells, cells.iloc[[0]].copy()], ignore_index=True)
    assert heatmaps.heatmap_test_counts(
        duplicated, ["A", "B"], traits, value_column="plot_minus_log10_fisher_p_value"
    ) == {
        "finite_test_count": 8,
        "displayed_finite_test_count": 3,
    }


def test_pattern_caption_counts_capped_rows_without_relabeling_legacy_count(
    display_inputs: tuple[pd.DataFrame, pd.DataFrame], tmp_path: Path
) -> None:
    cells, traits = display_inputs
    view = SimpleNamespace(view_key="fixture", view_label="Fixture")
    with (
        patch.object(heatmaps, "select_pattern_rows", return_value=["A", "B"]),
        patch.object(heatmaps, "plot_matrix_heatmap") as plot,
    ):
        result = heatmaps.write_pattern_heatmap(cells, view, traits, tmp_path)
    assert result[1:] == (2, 7, 9)
    subtitle = plot.call_args.kwargs["subtitle"]
    assert "7 finite Fisher tests in the source view; 3 finite tested cells shown" in subtitle
    assert "nominal Fisher evidence" in subtitle
    assert "Exploratory cells require nominal Fisher p < 0.01" in subtitle
    assert "nominal" in plot.call_args.kwargs["colorbar_label"]


@pytest.mark.parametrize("mode_key", ["grouped", "individual_conditions"])
def test_contrast_captions_separate_page_counts_and_source_tests(
    display_inputs: tuple[pd.DataFrame, pd.DataFrame], tmp_path: Path, mode_key: str
) -> None:
    cells, traits = display_inputs
    mode = next(mode for mode in heatmaps.CONTRAST_MODES if mode.key == mode_key)
    with (
        patch.object(heatmaps, "select_contrast_rows", return_value=["A", "B"]),
        patch.object(heatmaps, "CONTRAST_ROWS_PER_PAGE", 1),
        patch.object(heatmaps, "plot_matrix_heatmap") as plot,
        patch.object(heatmaps, "plot_clustered_contrast_heatmap") as clustered,
    ):
        outputs = heatmaps.write_contrast_heatmap(cells, mode, traits, tmp_path)
    assert outputs.tested_cell_count == 9
    assert outputs.displayed_cell_count == 7
    for call in plot.call_args_list:
        expected_shown = 3 if len(call.args[0]) == 2 else {"A": 1, "B": 2}[call.args[0].index[0]]
        assert (
            f"7 finite Fisher tests in the source view; {expected_shown} finite tested cells shown"
            in call.kwargs["subtitle"]
        )
        assert "Nominal p and BH q" in call.kwargs["subtitle"]
        assert "source-view cells" in call.kwargs["subtitle"]
    for call in clustered.call_args_list:
        assert (
            "7 finite Fisher tests in the source view; 3 finite tested cells shown"
            in call.kwargs["subtitle"]
        )


def test_generated_manifest_and_pages_append_finite_counts(
    display_inputs: tuple[pd.DataFrame, pd.DataFrame], tmp_path: Path
) -> None:
    cells, traits = display_inputs
    views = pd.DataFrame(
        [
            {
                "view_key": key,
                "view_label": key,
                "focus_direction": direction,
                "condition_mode": mode.condition_mode,
                "source_display_edges_path": f"{key}_display_edges.tsv",
            }
            for mode in heatmaps.CONTRAST_MODES
            for key, direction in [
                (mode.increased_view_key, "increased"),
                (mode.decreased_view_key, "decreased"),
            ]
        ]
    )
    exports: dict[str, pd.DataFrame] = {}
    with (
        patch.object(heatmaps, "clean_legacy_outputs"),
        patch.object(heatmaps, "read_trait_order", return_value=traits),
        patch.object(
            heatmaps,
            "read_tsv",
            side_effect=lambda path: (
                views if path.name == heatmaps.GRAPH_MANIFEST_FILENAME else pd.DataFrame()
            ),
        ),
        patch.object(heatmaps, "prepare_pattern_cells", side_effect=lambda *_args: cells.copy()),
        patch.object(heatmaps, "prepare_contrast_cells", side_effect=lambda *_args: cells.copy()),
        patch.object(heatmaps, "select_pattern_rows", return_value=["A", "B"]),
        patch.object(heatmaps, "select_contrast_rows", return_value=["A", "B"]),
        patch.object(heatmaps, "CONTRAST_ROWS_PER_PAGE", 1),
        patch.object(heatmaps, "plot_matrix_heatmap"),
        patch.object(heatmaps, "plot_clustered_contrast_heatmap"),
        patch.object(heatmaps, "build_condition_rollup_map", return_value=pd.DataFrame()),
        patch.object(
            heatmaps, "compare_grouped_and_individual_clusters", return_value=pd.DataFrame()
        ),
        patch.object(
            heatmaps,
            "write_tsv",
            side_effect=lambda frame, path: exports.update({path.name: frame.copy()}),
        ),
    ):
        manifest = heatmaps.write_six_shared_trait_cytoscape_heatmaps(tmp_path, tmp_path)
    assert manifest["tested_cell_count"].tolist() == [9] * len(manifest)
    assert manifest["finite_test_count"].tolist() == [7] * len(manifest)
    assert manifest["displayed_finite_test_count"].tolist() == [3] * len(manifest)
    pages = exports[heatmaps.CONTRAST_PAGES_FILENAME]
    assert pages["finite_test_count"].tolist() == [7] * 4
    assert pages["displayed_finite_test_count"].tolist() == [1, 2, 1, 2]


def test_readme_names_candidates_and_finite_tests_separately(tmp_path: Path) -> None:
    path = tmp_path / "README.md"
    path.write_text("# Fixture\n", encoding="utf-8")
    manifest = pd.DataFrame(
        [
            {
                "view_label": "Fixture",
                "heatmap_type": "trait_pattern",
                "png_path": "figure.png",
                "displayed_cell_count": 7,
                "nominal_p_value_cell_count": 2,
                "significant_cell_count": 2,
                "bh_selected_cell_count": 1,
                "tested_cell_count": 9,
                "finite_test_count": 7,
                "displayed_finite_test_count": 3,
                "row_count": 2,
            }
        ]
    )
    heatmaps.append_readme_section(tmp_path, manifest)
    text = path.read_text(encoding="utf-8")
    assert "did not pass nominal Fisher p<0.01" in text
    assert "Cells require nominal Fisher p<0.01" in text
    assert "Cell color is nominal Fisher evidence" in text
    assert "smaller nominal p-values" in text
    header = next(line for line in text.splitlines() if line.startswith("| heatmap |"))
    row = next(line for line in text.splitlines() if line.startswith("| Fixture |"))
    fields = dict(
        zip(
            [v.strip() for v in header.split("|")[1:-1]],
            [v.strip() for v in row.split("|")[1:-1]],
            strict=True,
        )
    )
    assert fields["candidate cells"] == "9"
    assert fields["finite tests (source view)"] == "7"
    assert fields["finite tested cells shown"] == "3"
    assert (
        "The legacy `tested_cell_count` counts candidate cells, including untestable cells." in text
    )


def test_added_nominal_text_fits_synthetic_pattern_render(tmp_path: Path) -> None:
    traits = pd.DataFrame(
        {"trait_id": [f"t{i}" for i in range(6)], "trait_label": [f"Trait {i}" for i in range(6)]}
    )
    cells = pd.DataFrame(
        {
            "condition_body_site_pair": ["A"] * 6,
            "trait_id": traits["trait_id"],
            "fisher_p_value": [0.001] * 6,
            "displayed_in_heatmap": [True] * 6,
            "focus_positive_profile_count": [4] * 6,
            "plot_minus_log10_fisher_p_value": [3.0] * 6,
        }
    )
    captured: list[bool] = []

    def capture(fig: matplotlib.figure.Figure, *_args: object, **_kwargs: object) -> None:
        fig.canvas.draw()
        new_texts = [text for text in fig.axes[0].texts if "nominal p" in text.get_text()]
        assert len(new_texts) == 6
        for index, text in enumerate(new_texts):
            bounds = text.get_window_extent()
            left = fig.axes[0].transData.transform((index - 0.5, 0))[0]
            right = fig.axes[0].transData.transform((index + 0.5, 0))[0]
            assert left <= bounds.x0 < bounds.x1 <= right
        subtitle = next(text for text in fig.texts if "finite Fisher tests" in text.get_text())
        bounds = subtitle.get_window_extent()
        assert fig.bbox.contains(bounds.x0, bounds.y0)
        assert fig.bbox.contains(bounds.x1, bounds.y1)
        assert bounds.y0 > fig.axes[0].get_window_extent().y1
        captured.append(True)

    with (
        patch.object(heatmaps, "select_pattern_rows", return_value=["A"]),
        patch.object(matplotlib.figure.Figure, "savefig", capture),
    ):
        heatmaps.write_pattern_heatmap(
            cells, SimpleNamespace(view_key="fixture", view_label="Fixture"), traits, tmp_path
        )
    assert captured == [True]


def test_added_contrast_counts_fit_synthetic_overview_and_pages(
    display_inputs: tuple[pd.DataFrame, pd.DataFrame], tmp_path: Path
) -> None:
    cells, traits = display_inputs
    mode = next(mode for mode in heatmaps.CONTRAST_MODES if mode.key == "individual_conditions")
    captured: list[bool] = []

    def capture(fig: matplotlib.figure.Figure, *_args: object, **_kwargs: object) -> None:
        fig.canvas.draw()
        subtitle = next(text for text in fig.texts if "finite Fisher tests" in text.get_text())
        bounds = subtitle.get_window_extent()
        assert fig.bbox.contains(bounds.x0, bounds.y0)
        assert fig.bbox.contains(bounds.x1, bounds.y1)
        assert bounds.y0 > fig.axes[0].get_window_extent().y1
        captured.append(True)

    with (
        patch.object(heatmaps, "select_contrast_rows", return_value=["A", "B"]),
        patch.object(heatmaps, "CONTRAST_ROWS_PER_PAGE", 1),
        patch.object(heatmaps, "plot_clustered_contrast_heatmap"),
        patch.object(matplotlib.figure.Figure, "savefig", capture),
    ):
        heatmaps.write_contrast_heatmap(cells, mode, traits, tmp_path)
    assert captured == [True] * 3
