import importlib
import sys
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest
from matplotlib.figure import Figure
from matplotlib.patches import PathPatch

with patch.object(sys, "path", [str(Path(__file__).resolve().parents[1] / "scripts"), *sys.path]):
    sankeys = importlib.import_module("create_six_shared_trait_sankeys")


@pytest.fixture(autouse=True)
def headless_figures() -> None:
    sankeys.plt.switch_backend("Agg")


def direction_edges() -> pd.DataFrame:
    records = []
    for direction, count in [("increased", 40), ("decreased", 35)]:
        for rank, trait in enumerate(sankeys.SHARED_TRAITS):
            records.append(
                {
                    "condition": "Shared condition",
                    "body_site_group": "Feces",
                    "focus_direction": direction,
                    "trait_id": trait.trait_id,
                    "trait_rank": rank,
                    "focus_enriched_profile_count": count,
                    "focus_profile_count": count + 10,
                    "focus_enriched_study_ids": "study1;study2",
                    "background_no_enriched_profile_count": 10,
                    "focus_enriched_mean_guild_z": 1.0,
                    "focus_enriched_minus_comparator_guild_z": -0.5,
                    "rank_score": count,
                }
            )
    # These links must not become eligible just because two directions sum to 40.
    for direction in ["increased", "decreased"]:
        records.append(
            {
                **records[0],
                "condition": "Below threshold",
                "focus_direction": direction,
                "focus_enriched_profile_count": 20,
            }
        )
    return pd.DataFrame.from_records(records)


def test_bidirectional_selection_preserves_direction_counts(tmp_path: Path) -> None:
    edges = direction_edges()
    with (
        patch.object(sankeys, "OUTPUT_DIR", tmp_path),
        patch.object(sankeys, "REPO_ROOT", tmp_path),
        patch.object(sankeys, "plot_combined_direction_sankey"),
    ):
        manifest = sankeys.write_bidirectional_high_support_outputs(edges, edges)
    displayed = pd.read_csv(tmp_path / f"{sankeys.BIDIRECTIONAL_STEM}_display_edges.tsv", sep="\t")
    assert set(displayed["condition"]) == {"Shared condition"}
    assert len(displayed) == 12
    for direction, count, style in [("increased", 40, "solid"), ("decreased", 35, "dashed")]:
        subset = displayed[displayed["focus_direction"] == direction]
        assert len(subset) == 6
        assert subset["focus_enriched_profile_count"].eq(count).all()
        assert subset["line_style"].eq(style).all()
    assert manifest.iloc[0]["selected_pair_count"] == 1
    assert manifest.iloc[0]["display_edge_count"] == 12


def test_bidirectional_paths_have_distinct_offsets_and_correct_styles(tmp_path: Path) -> None:
    edges = direction_edges().query("condition == 'Shared condition'").copy()
    with patch.object(Figure, "savefig", autospec=True) as save:
        sankeys.plot_combined_direction_sankey(
            edges,
            "both",
            tmp_path / "test.png",
            color_scale_vmax=2.0,
            profile_node_scale=95,
        )
    figure = save.call_args.args[0]
    paths = [p for p in figure.axes[0].patches if isinstance(p, PathPatch)]
    assert len(paths) == 12  # No background underlays obscure the gaps.
    assert sum(p.get_linestyle() != "-" for p in paths) == 6
    assert sum(p.get_linestyle() == "-" for p in paths) == 6
    vertices = [np.asarray(p.get_path().vertices) for p in paths]
    assert len({tuple(v[1]) for v in vertices}) == 12
    assert len({tuple(v[0]) for v in vertices}) == 1
    assert len({tuple(v[-1]) for v in vertices}) == 1
    scales = sankeys.trait_color_scales(edges, shared_vmax=2.0)
    ordered = sankeys.add_parallel_offsets(edges, separate_directions=True).sort_values(
        "parallel_rank"
    )
    for (_, edge), curve in zip(ordered.iterrows(), paths, strict=True):
        cmap, norm = scales[edge.trait_id]
        # Shade reflects the mean positive z-score, even when the contrast is negative.
        assert curve.get_edgecolor()[:3] == pytest.approx(cmap(norm(1.0))[:3])
        expected_style = (
            sankeys.decreased_dash_style(curve.get_linewidth())
            if edge.focus_direction == "decreased"
            else "-"
        )
        assert curve.get_linestyle() == expected_style
    axis_points = figure.axes[0].get_position().height * figure.get_figheight() * 72
    midpoints = [(v[0, 1] + 3 * v[1, 1] + 3 * v[2, 1] + v[3, 1]) / 8 for v in vertices]
    for index in range(len(paths) - 1):
        separation = abs(midpoints[index] - midpoints[index + 1]) * axis_points
        half_widths = (paths[index].get_linewidth() + paths[index + 1].get_linewidth()) / 2
        assert separation - half_widths == pytest.approx(3.5)


def test_node_profiles_are_not_multiplied_by_traits() -> None:
    edges = direction_edges().query("condition == 'Shared condition'").copy()
    stats = sankeys.combined_node_stats(edges, "condition", count_profiles=True)
    assert stats.iloc[0]["profile_count"] == 95
    assert stats.iloc[0]["focus_enriched_profile_count"] == 6 * (40 + 35)
    small = sankeys.profile_node_radius(25, 100)
    large = sankeys.profile_node_radius(100, 100)
    assert small**2 / large**2 == pytest.approx(0.25)
    edges.loc[edges.index[0], "focus_profile_count"] = 99
    with pytest.raises(ValueError, match="disagree"):
        sankeys.combined_node_stats(edges, "condition", count_profiles=True)


def test_single_direction_curves_keep_solid_backgrounds(tmp_path: Path) -> None:
    edges = direction_edges().query("focus_direction == 'increased'")
    with patch.object(Figure, "savefig", autospec=True) as save:
        sankeys.plot_combined_direction_sankey(edges, "inc", tmp_path / "test.png")
    paths = [p for p in save.call_args.args[0].axes[0].patches if isinstance(p, PathPatch)]
    assert len(paths) == 2 * len(edges)
    assert all(p.get_linestyle() == "-" for p in paths)


def test_bidirectional_readme_note_preserves_other_sections(tmp_path: Path) -> None:
    path = tmp_path / "README.md"
    path.write_text("# Figures\n\n## Heatmaps\nExisting p/q notes.\n", encoding="utf-8")
    with patch.object(sankeys, "OUTPUT_DIR", tmp_path):
        sankeys.write_bidirectional_readme_note()
        first = path.read_text(encoding="utf-8")
        sankeys.write_bidirectional_readme_note()
    assert path.read_text(encoding="utf-8") == first
    assert "## Heatmaps\nExisting p/q notes." in first
