"""The enrichment follow-up names target conditions it cannot summarize (#586)."""

from __future__ import annotations

import importlib.util
import sys
from collections.abc import Callable
from pathlib import Path
from typing import cast

import numpy as np
import pandas as pd
import pytest
from matplotlib.backends.backend_agg import FigureCanvasAgg, RendererAgg
from matplotlib.collections import QuadMesh
from matplotlib.colors import rgb_to_hsv, to_hex, to_rgb
from matplotlib.figure import Figure
from matplotlib.text import Text

MODULE_PATH = (
    Path(__file__).resolve().parents[1] / "scripts" / "analyze_enrichment_infographic_followup.py"
)
SPEC = importlib.util.spec_from_file_location("enrichment_infographic_followup", MODULE_PATH)
assert SPEC is not None
assert SPEC.loader is not None
followup = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = followup
SPEC.loader.exec_module(followup)


@pytest.fixture(params=[False, True], ids=["supplied-order", "clustered-columns"])
def companion_heatmap_at_save(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, request: pytest.FixtureRequest
) -> tuple[Figure, pd.DataFrame]:
    """Draw the producer's actual artists at its save boundary (#1861)."""
    matrix = pd.DataFrame(
        [[0, 1, 1], [1, 0, 1]],
        index=["taxon-a", "taxon-b"],
        columns=["METPO:3", "METPO:1", "METPO:2"],
    )
    before = matrix.copy(deep=True)
    taxa = pd.DataFrame(
        {
            "taxon_name": ["Synthetic taxon A", "Synthetic taxon B"],
            "target_profile_count": [1, 1],
            "display_condition": [followup.TARGET_CONDITIONS[0]] * 2,
        },
        index=matrix.index,
    )
    traits = pd.DataFrame({"heatmap_label": ["Third", "First", "Second"]}, index=matrix.columns)
    captured: list[Figure] = []
    output = tmp_path / "heatmap.png"

    def at_save(figure: Figure, path: Path, **kwargs: object) -> None:
        assert path == output
        assert kwargs == {"dpi": 300, "bbox_inches": "tight"}
        draw = cast(Callable[[], None], FigureCanvasAgg(figure).draw)
        draw()
        captured.append(figure)

    monkeypatch.setattr(Figure, "savefig", at_save)
    displayed = followup.save_traits_vs_organisms_heatmap(
        matrix, taxa, traits, output, cluster_columns=request.param
    )
    pd.testing.assert_frame_equal(matrix, before)
    assert len(captured) == 1
    return captured[0], displayed


def test_companion_heatmap_mesh_strokes_have_a_half_point_floor(
    companion_heatmap_at_save: tuple[Figure, pd.DataFrame],
) -> None:
    figure, _ = companion_heatmap_at_save
    meshes = figure.findobj(QuadMesh)
    assert len(meshes) >= 2, "Both the binary matrix and condition strip must be inspected"
    assert all(np.min(mesh.get_linewidth()) >= 0.5 for mesh in meshes)


def test_companion_heatmap_colorbar_text_clears_the_condition_strip(
    companion_heatmap_at_save: tuple[Figure, pd.DataFrame],
) -> None:
    figure, _ = companion_heatmap_at_save
    get_renderer = cast(Callable[[], RendererAgg], FigureCanvasAgg(figure).get_renderer)
    renderer = get_renderer()
    colorbar = next(axis for axis in figure.axes if axis.get_ylabel() == "Trait present")
    other_mesh_axes = [
        axis for axis in figure.axes if axis is not colorbar and axis.findobj(QuadMesh)
    ]
    assert len(other_mesh_axes) == 2
    words = [
        text for text in colorbar.findobj(Text) if text.get_visible() and text.get_text().strip()
    ]
    assert words
    assert any(text.get_text().startswith("0") for text in words)
    for text in words:
        bounds = text.get_window_extent(renderer)
        for axis in other_mesh_axes:
            assert not bounds.overlaps(axis.get_window_extent(renderer)), text.get_text()


def test_companion_heatmap_layout_keeps_every_labelled_binary_cell(
    companion_heatmap_at_save: tuple[Figure, pd.DataFrame],
) -> None:
    figure, expected = companion_heatmap_at_save
    axis = next(
        axis
        for axis in figure.axes
        if axis.get_ylabel() == "Organism observed in target-condition profiles"
    )
    meshes = axis.findobj(QuadMesh)
    assert len(meshes) == 1
    rows = [label.get_text() for label in axis.get_yticklabels()]
    columns = [label.get_text() for label in axis.get_xticklabels()]
    actual = pd.DataFrame(
        np.asarray(meshes[0].get_array()).reshape(len(rows), len(columns)),
        index=rows,
        columns=columns,
    )
    pd.testing.assert_frame_equal(
        actual.sort_index().sort_index(axis=1),
        expected.sort_index().sort_index(axis=1),
        check_dtype=False,
    )


def _cie_lightness(rgb: tuple[float, float, float]) -> float:
    channels = np.asarray(rgb, dtype=float)
    linear = np.where(channels <= 0.04045, channels / 12.92, ((channels + 0.055) / 1.055) ** 2.4)
    luminance = float(0.2126 * linear[0] + 0.7152 * linear[1] + 0.0722 * linear[2])
    return 116.0 * (luminance ** (1.0 / 3.0)) - 16.0 if luminance > 0.008856 else 903.3 * luminance


def test_companion_heatmap_multi_condition_palette_and_legend_clear_dendrogram(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Multi-condition strip colours stay distinct and clear all axes (#1922, #2005, #2006)."""
    conditions = (
        *followup.TARGET_CONDITIONS,
        "Multiple sclerosis",
        "Stroke",
        "Constipation",
        "Colorectal cancer",
        "Parkinson's disease",
    )[:8]
    assert len(conditions) == 8
    monkeypatch.setattr(followup, "TARGET_CONDITIONS", conditions)
    matrix = pd.DataFrame(
        [[int((i >> j) & 1) for j in range(4)] for i in range(len(conditions))],
        index=[f"taxon-{i}" for i in range(len(conditions))],
        columns=["METPO:1", "METPO:2", "METPO:3", "METPO:4"],
    )
    taxa = pd.DataFrame(
        {
            "taxon_name": [f"Taxon {i}" for i in range(len(conditions))],
            "target_profile_count": list(range(1, len(conditions) + 1)),
            "display_condition": list(conditions),
        },
        index=matrix.index,
    )
    traits = pd.DataFrame(
        {"heatmap_label": ["First", "Second", "Third", "Fourth"]}, index=matrix.columns
    )
    captured: list[Figure] = []

    def at_save(figure: Figure, _path: Path, **_kwargs: object) -> None:
        draw = cast(Callable[[], None], FigureCanvasAgg(figure).draw)
        draw()
        captured.append(figure)

    monkeypatch.setattr(Figure, "savefig", at_save)
    followup.save_traits_vs_organisms_heatmap(
        matrix, taxa, traits, tmp_path / "multi.png", cluster_columns=True
    )
    assert len(captured) == 1
    figure = captured[0]
    get_renderer = cast(Callable[[], RendererAgg], FigureCanvasAgg(figure).get_renderer)
    renderer = get_renderer()
    dendrogram = next(axis for axis in figure.axes if axis.get_legend() is not None)
    legend = dendrogram.get_legend()
    assert legend is not None
    patches = legend.get_patches()
    assert len(patches) == len(conditions)
    legend_hexes = [to_hex(patch.get_facecolor()) for patch in patches]
    assert len(set(legend_hexes)) == len(conditions), legend_hexes
    trait_rgb = to_rgb("#155e75")
    trait_hue = float(rgb_to_hsv(trait_rgb)[0])
    palette_rgbs = [to_rgb(patch.get_facecolor()) for patch in patches]
    for color_hex, color_rgb in zip(legend_hexes, palette_rgbs, strict=True):
        hue_diff = abs(float(rgb_to_hsv(color_rgb)[0]) - trait_hue)
        assert min(hue_diff, 1.0 - hue_diff) >= 0.03, color_hex
    data_rgbs = [trait_rgb, *palette_rgbs]
    lightness = [_cie_lightness(rgb) for rgb in data_rgbs]
    for i in range(len(lightness)):
        for j in range(i + 1, len(lightness)):
            assert abs(lightness[i] - lightness[j]) >= 8.0, (i, j, lightness)

    heatmap_ax = next(axis for axis in figure.axes if axis.get_xlabel() == "Trait")
    strip_ax = next(
        axis
        for axis in figure.axes
        if axis is not heatmap_ax
        and axis.collections
        and len(axis.collections[0].get_facecolor()) == len(conditions)
    )
    strip_mesh = cast(QuadMesh, strip_ax.collections[0])
    strip_hexes = [
        to_hex(cast(tuple[float, float, float, float], color))
        for color in strip_mesh.get_facecolor()
    ]
    label_to_index = {
        followup.display_taxon_label(taxon_id, taxa): idx
        for idx, taxon_id in enumerate(matrix.index)
    }
    row_labels = [tick.get_text() for tick in heatmap_ax.get_yticklabels()]
    assert len(row_labels) == len(conditions)
    assert strip_hexes == [legend_hexes[label_to_index[label]] for label in row_labels]
    strip_in_condition_order = [""] * len(conditions)
    for row_pos, label in enumerate(row_labels):
        strip_in_condition_order[label_to_index[label]] = strip_hexes[row_pos]
    assert strip_in_condition_order == legend_hexes

    dendrogram_box = dendrogram.get_window_extent(renderer)
    legend_box = legend.get_window_extent(renderer)
    assert not legend_box.overlaps(dendrogram_box)
    for axis in figure.axes:
        if axis is dendrogram:
            continue
        assert not legend_box.overlaps(axis.get_window_extent(renderer))
        tight_box = axis.get_tightbbox(renderer)
        assert tight_box is not None
        assert not legend_box.overlaps(tight_box)
    legend_words = [
        text for text in legend.findobj(Text) if text.get_visible() and text.get_text().strip()
    ]
    assert legend_words
    for word in legend_words:
        assert not word.get_window_extent(renderer).overlaps(dendrogram_box), word.get_text()


def test_targets_below_the_profile_minimum_are_reported_with_their_counts() -> None:
    present, rare, absent = followup.TARGET_CONDITIONS
    metadata = pd.DataFrame(
        {
            "profile_id": ["p1", "p2", "p3"],
            "condition": [present, present, rare],
            "study_id": ["s1", "s1", "s2"],
            "experiment_id": ["e1", "e2", "e3"],
            "profile_trait_count": [5, 3, 4],
            "total_trait_signal": [10, 6, 8],
        }
    )

    depth = followup.summarize_condition_depth(metadata, min_profiles=2)

    assert depth["condition"].tolist() == [present]
    assert followup.dropped_target_conditions(metadata, depth) == {rare: 1, absent: 0}


def _profiles(rows: list[tuple[str, str, str, str]]) -> pd.DataFrame:
    """Profiles from (profile_id, study_id, experiment_id, condition) rows."""
    frame = pd.DataFrame(rows, columns=["profile_id", "study_id", "experiment_id", "condition"])
    frame["inc_dec"] = frame["profile_id"].str.rsplit("|", n=1).str[-1]
    frame["direction_label"] = frame["inc_dec"].map({"inc": "increased", "dec": "decreased"})
    frame["body_site"] = "Feces"
    frame["location"] = "Japan"
    frame["host_species"] = "Homo sapiens"
    frame["title"] = "A study"
    frame["experiment_label"] = "Experiment 1: control vs case"
    frame["profile_trait_count"] = 3
    frame["total_trait_signal"] = 6.0
    return frame.set_index("profile_id", drop=False)


def test_title_table_lists_corpus_experiments_and_marks_missing_titles(tmp_path: Path) -> None:
    target = followup.TARGET_CONDITIONS[0]
    metadata = _profiles(
        [
            ("s1|s1-1|inc", "s1", "s1-1", target),
            ("s1|s1-1|dec", "s1", "s1-1", target),
            ("s2|s2-1|inc", "s2", "s2-1", target),
        ]
    )
    experiments = pd.DataFrame(
        {
            "study_id": ["s1", "s2", "s3"],
            "experiment_id": ["s1-1", "s2-1", "s3-1"],
            "condition": [target] * 3,
            "efo_id": ["MONDO:1"] * 3,
            "experiment_label": ["Experiment 1: control vs case"] * 3,
            "group_1_definition": ["cases"] * 3,
            "sequencing_type": ["16S"] * 3,
        }
    )
    studies = pd.DataFrame(
        {"study_id": ["s1"], "title": ["Study one"], "publication_year": ["2020"]}
    ).assign(doi="", url="")

    table = followup.write_target_title_table(
        metadata, experiments, studies, tmp_path / "titles.tsv"
    ).set_index("experiment_id")

    # s3-1 carries the label in the dump but has no profile in this corpus.
    assert table.index.tolist() == ["s1-1", "s2-1"]
    assert table.loc["s1-1", "profile_count"] == 2
    assert table.loc["s1-1", "title"] == "Study one"
    assert table.loc["s2-1", "title"] == followup.MISSING_TITLE
    assert table.loc["s2-1", "publication_year"] == "unknown"


def test_cross_target_check_leaves_out_targets_it_never_queried(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    queried, unqueried = "Queried condition", "Rare condition"
    monkeypatch.setattr(followup, "TARGET_CONDITIONS", (queried, unqueried))
    metadata = _profiles(
        [
            ("s1|s1-1|inc", "s1", "s1-1", queried),
            ("s2|s2-1|inc", "s2", "s2-1", queried),
            ("s3|s3-1|inc", "s3", "s3-1", unqueried),
            ("s4|s4-1|inc", "s4", "s4-1", "Other condition"),
        ]
    )
    values = pd.DataFrame(
        [[3, 1, 0], [2, 2, 0], [0, 1, 4], [1, 0, 2]],
        index=metadata.index,
        columns=["METPO:1", "METPO:2", "METPO:3"],
    )

    neighbors, crosscheck, _ = followup.compute_profile_neighbors(
        values, metadata, top_n=5, min_condition_profiles=2
    )

    assert set(neighbors["query_condition"]) == {queried}
    assert crosscheck["query_condition"].tolist() == [queried]


def _paired_trait_fixture(tmp_path: Path) -> tuple[Path, pd.DataFrame]:
    """Deliberately disagreeing row, p-value and q-value orders."""
    path = tmp_path / "paired-traits.tsv"
    pd.DataFrame(
        {
            "trait": ["Third METPO:3", "First METPO:1", "Again METPO:3", "Second METPO:2"],
            "p_value": [0.8, 0.01, 0.2, 0.4],
            "q_value": [0.2, 0.7, 0.5, 0.01],
        }
    ).to_csv(path, sep="\t", index=False)
    labels = pd.DataFrame(
        {
            "display_name": ["First", "Second", "Third", "Extra"],
            "display_category": ["phenotype"] * 4,
        },
        index=["METPO:1", "METPO:2", "biolink:has_phenotype|METPO:3", "METPO:4"],
    )
    return path, labels


def test_paired_trait_table_order_ignores_p_and_q_values(tmp_path: Path) -> None:
    path, labels = _paired_trait_fixture(tmp_path)
    expected = ["biolink:has_phenotype|METPO:3", "METPO:1", "METPO:2"]
    assert followup.load_top_paired_traits(path, 3, labels) == expected

    changed = pd.read_csv(path, sep="\t")
    changed["p_value"] = [0.001, 0.9, 0.6, 0.3]
    changed["q_value"] = [0.8, 0.04, 0.1, 0.6]
    changed_path = tmp_path / "changed-p-and-q.tsv"
    changed.to_csv(changed_path, sep="\t", index=False)
    assert followup.load_top_paired_traits(changed_path, 3, labels) == expected


@pytest.mark.parametrize("fill", [False, True], ids=["no-fill", "fill"])
def test_main_heatmap_caption_matches_supplied_order_and_clustering(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, fill: bool
) -> None:
    path, labels = _paired_trait_fixture(tmp_path)
    metadata = _profiles(
        [
            ("s1|s1-1|inc", "s1", "s1-1", followup.TARGET_CONDITIONS[0]),
            ("s1|s1-1|dec", "s1", "s1-1", followup.TARGET_CONDITIONS[0]),
        ]
    )
    values = pd.DataFrame(
        [[3, 1, 0, 20], [0, 2, 4, 10]], index=metadata.index, columns=labels.index
    )
    binary = pd.DataFrame(
        [[1, 1, 0, 1], [0, 1, 1, 1]],
        index=["taxon-a", "taxon-b"],
        columns=labels.index,
    )
    output = tmp_path / "output"
    monkeypatch.setattr(followup, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "followup",
            "--output-dir",
            str(output),
            "--artifact-prefix",
            "fixture_",
            "--trait-source",
            "paired-top",
            "--paired-trait-stats-input",
            str(path),
            "--top-paired-traits",
            "3",
            "--max-heatmap-traits",
            "4",
            *(["--fill-heatmap-traits"] if fill else []),
        ],
    )
    for name in ("load_experiments", "load_studies", "load_profile_labels"):
        monkeypatch.setattr(followup, name, lambda *_args: pd.DataFrame())
    monkeypatch.setattr(followup, "load_trait_labels", lambda *_args: labels)
    monkeypatch.setattr(followup, "load_profiles", lambda *_args: (values, metadata, None))
    for name in (
        "summarize_condition_depth",
        "write_target_title_table",
        "summarize_autism_labels",
        "summarize_guild_connections",
    ):
        monkeypatch.setattr(followup, name, lambda *_args: pd.DataFrame())
    monkeypatch.setattr(followup, "dropped_target_conditions", lambda *_args: {})
    monkeypatch.setattr(
        followup,
        "compute_profile_neighbors",
        lambda *_args: (pd.DataFrame(), pd.DataFrame(), pd.DataFrame()),
    )
    monkeypatch.setattr(followup, "write_report", lambda *_args: None)
    expected_order = ["biolink:has_phenotype|METPO:3", "METPO:1", "METPO:2"]

    def build_matrix(
        actual_values: pd.DataFrame,
        actual_metadata: pd.DataFrame,
        actual_labels: pd.DataFrame,
        _full_dump_path: Path,
        _matrix_path: Path,
        max_traits: int,
        _max_taxa: int,
        preferred_traits: list[str],
        fill_with_high_signal_traits: bool,
        balance_taxa_by_condition: bool,
    ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        assert preferred_traits == expected_order
        assert max_traits == (4 if fill else 3)
        assert fill_with_high_signal_traits is fill
        assert balance_taxa_by_condition is False
        selected = followup.select_heatmap_traits(
            actual_values,
            actual_metadata,
            actual_labels,
            max_traits,
            preferred_traits,
            fill_with_high_signal_traits,
        )
        return (
            binary[selected].copy(),
            pd.DataFrame(index=binary.index),
            actual_labels.loc[selected],
        )

    captured: list[tuple[pd.DataFrame, Path, str, bool]] = []

    def capture_heatmap(
        matrix: pd.DataFrame,
        taxa_metadata: pd.DataFrame,
        trait_metadata: pd.DataFrame,
        output_path: Path,
        title: str,
        cluster_columns: bool,
    ) -> pd.DataFrame:
        assert taxa_metadata.index.tolist() == matrix.index.tolist()
        assert trait_metadata.index.tolist() == matrix.columns.tolist()
        captured.append((matrix.copy(), output_path, title, cluster_columns))
        return matrix.copy()

    monkeypatch.setattr(followup, "build_traits_vs_organisms_heatmap", build_matrix)
    monkeypatch.setattr(followup, "save_traits_vs_organisms_heatmap", capture_heatmap)
    followup.main()

    assert len(captured) == 1
    matrix, image_path, title, cluster_columns = captured[0]
    assert image_path == output / "fixture_traits_vs_organisms_heatmap.png"
    assert matrix.index.tolist() == ["taxon-a", "taxon-b"]
    assert matrix.columns.tolist() == [*expected_order, *(["METPO:4"] if fill else [])]
    assert matrix.to_numpy().tolist() == (
        [[0, 1, 1, 1], [1, 0, 1, 1]] if fill else [[0, 1, 1], [1, 0, 1]]
    )
    assert cluster_columns is fill
    if fill:
        assert title == (
            "Companion Heatmap: Relevant Traits vs Organisms in Infographic Outlier Conditions\n"
            "Trait columns include paired-profile top traits plus high-signal "
            "target-condition traits"
        )
        assert "supplied paired-trait table order" not in title
    else:
        assert title == (
            "Companion Heatmap: Our Top Paired-Profile Traits vs Organisms in "
            "Infographic Outlier Conditions\n"
            "Trait columns follow the supplied paired-trait table order"
        )
    assert "p-value" not in title
    assert "q-value" not in title
    saved_matrix = pd.read_csv(
        output / "fixture_traits_vs_organisms_heatmap_matrix.tsv", sep="\t", index_col=0
    )
    pd.testing.assert_frame_equal(saved_matrix, matrix)


def test_companion_heatmap_refuses_more_conditions_than_distinct_colours(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A ninth condition must not silently reuse a strip colour (#2005 delta)."""
    conditions = tuple(f"Condition {i}" for i in range(9))
    monkeypatch.setattr(followup, "TARGET_CONDITIONS", conditions)
    matrix = pd.DataFrame([[1, 0], [0, 1]], index=["taxon-a", "taxon-b"], columns=["T:1", "T:2"])
    taxa = pd.DataFrame(
        {
            "taxon_name": ["Taxon A", "Taxon B"],
            "target_profile_count": [1, 1],
            "display_condition": [conditions[0], conditions[8]],
        },
        index=matrix.index,
    )
    traits = pd.DataFrame({"heatmap_label": ["One", "Two"]}, index=matrix.columns)
    output = tmp_path / "nine.png"
    with pytest.raises(
        ValueError, match=r"^The companion heatmap has 8 distinct condition colours; 9"
    ):
        followup.save_traits_vs_organisms_heatmap(matrix, taxa, traits, output)
    assert not output.exists()
