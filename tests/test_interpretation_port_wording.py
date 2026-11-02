"""Guard the annotation wording ported from Jed's PR #193 (#1262)."""

from __future__ import annotations

import inspect
from pathlib import Path

import numpy as np
import pytest
from PIL import Image
from plotly.graph_objects import Figure

from mtp import guild_maps, table1_trait_maps


@pytest.fixture
def annotation_matrix() -> guild_maps.GuildMatrix:
    """Include a shaded annotation-zero cell beside an observed fraction."""
    return guild_maps.GuildMatrix(
        conditions=["condition"],
        body_sites=["Feces"],
        columns=["Feces|aerobe_increased", "Feces|aerobe_decreased"],
        values=[[None, 0.25]],
        study_counts=[[3, 3]],
        row_clusters=[1],
        pooled=[],
        dropped_empty_conditions=[],
        dropped_body_sites=[],
        sampled_no_signal=[[True, False]],
    )


@pytest.mark.parametrize("renderer", ["combined", "body_site"])
@pytest.mark.parametrize("mark_no_signal", [True, False])
@pytest.mark.parametrize("group_noun", ["guild", "trait"])
def test_heatmap_zero_caption_describes_annotation_evidence(
    annotation_matrix: guild_maps.GuildMatrix,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    renderer: str,
    mark_no_signal: bool,
    group_noun: str,
) -> None:
    figures: list[Figure] = []

    def capture(figure: Figure, *_args: object, **_kwargs: object) -> None:
        figures.append(figure)

    monkeypatch.setattr(Figure, "write_image", capture)
    output = tmp_path / "heatmap.png"
    if renderer == "combined":
        guild_maps.write_figure(
            annotation_matrix,
            output,
            "fixture",
            guild_order=("aerobe",),
            mark_no_signal=mark_no_signal,
            group_noun=group_noun,
        )
    else:
        guild_maps.write_body_site_panels(
            annotation_matrix,
            output,
            "fixture",
            sites=("Feces",),
            guild_order=("aerobe",),
            mark_no_signal=mark_no_signal,
            group_noun=group_noun,
        )

    assert len(figures) == 1
    title = str(figures[0].layout.title.text)
    assert "no taxon carried" not in title
    assert "an absent trait" not in title
    caption = f"Dark grey = sampled, with no recorded {group_noun} annotation signal."
    assert (caption in title) is mark_no_signal
    assert "summed trait-assertion counts" in title
    assert "Directions are curated signature labels." in title
    assert "Annotation shares do not measure patient traits or exposure." in title
    # The wording change must not rename the machine-readable no-signal state.
    assert guild_maps.CELL_NO_SIGNAL == "no_guild_signal"
    assert annotation_matrix.values == [[None, 0.25]]
    assert annotation_matrix.sampled_no_signal == [[True, False]]


def test_no_signal_documentation_does_not_claim_trait_absence() -> None:
    description = inspect.getdoc(guild_maps.write_figure)
    assert description is not None
    assert "no group annotations were recorded" in description
    assert "zero mapped signal" in description
    assert "an absent trait" not in description
    assert "no recorded annotation" in inspect.getsource(guild_maps.GuildMatrix)


def test_correlation_documentation_does_not_equate_measurements() -> None:
    description = inspect.getdoc(table1_trait_maps.build_shift_figure)
    assert description is not None
    assert "numerical overlap" in description
    assert "not semantic equivalence or independent evidence" in description
    assert "close to the same measurement" not in description


def test_standalone_shift_figure_names_the_annotation_contrast() -> None:
    shift = table1_trait_maps.TraitShift(
        trait="coccus_shaped",
        studies_up=2,
        studies_down=1,
        studies_tied=1,
        share_up=2 / 3,
        ci_low=0.1,
        ci_high=0.95,
        p_value=1.0,
        q_value=1.0,
    )
    result = table1_trait_maps.ShiftResult(shifts=[shift], spearman={}, studies=4)
    figure = table1_trait_maps.build_shift_figure(result, "fixture", ("coccus_shaped",))
    title = str(figure.layout.title.text)
    annotations = {str(annotation.text) for annotation in figure.layout.annotations}

    assert "paired annotation contrast" in title
    assert "paired directional shift" not in title
    assert "Share of studies with positive annotation contrast" in annotations
    assert "Share of studies shifted toward the increased side" not in annotations
    assert "annotation" in str(figure.layout.xaxis.title.text)
    assert "mean increased annotation share" in title
    assert "mean decreased annotation share" in title
    assert "Directions are curated signature labels." in title
    assert (
        "Shared annotations do not establish patient-level change or independent validation."
        in title
    )
    # Keep the existing inference and support displays and saved values intact.
    assert "Benjamini-Hochberg across the six traits" in title
    assert "2/3 studies, q=1" in annotations
    assert list(figure.data[0].x) == [2 / 3]
    assert list(figure.data[0].error_x.array) == [0.95 - 2 / 3]
    assert list(figure.data[0].error_x.arrayminus) == [2 / 3 - 0.1]


@pytest.mark.parametrize(
    "vintage",
    [
        "synthetic fixture, 2026-10-30",
        "synthetic fixture, " + "extended export provenance " * 8,
        None,
    ],
    ids=["synthetic", "long-vintage", "missing-vintage"],
)
def test_rendered_shift_subtitle_lines_fit_canvas_and_clear_panel_headers(
    tmp_path: Path, vintage: str | None
) -> None:
    """#1276: measure real rendered ink, including clipped subtitle edge pixels."""
    traits = table1_trait_maps.TRAIT_ORDER
    result = table1_trait_maps.ShiftResult(
        shifts=[
            table1_trait_maps.TraitShift(
                trait=trait,
                studies_up=500,
                studies_down=400,
                studies_tied=22,
                share_up=500 / 900,
                ci_low=0.52,
                ci_high=0.59,
                p_value=0.001,
                q_value=0.002,
            )
            for trait in traits
        ],
        spearman={
            f"{left}~{right}": (0.3, 0.01, 900)
            for index, left in enumerate(traits)
            for right in traits[index + 1 :]
        },
        studies=922,
    )
    figure = table1_trait_maps.build_shift_figure(result, vintage, traits)
    figure.layout.title.text = "Synthetic fixture: " + str(figure.layout.title.text)
    assert figure.layout.width == 1650
    # These diagnostic colors affect only paint, not text or layout geometry.
    figure.update_layout(title_font_color="#ff00ff")
    for annotation in figure.layout.annotations:
        annotation.font.color = "#00ffff"
    output = tmp_path / "shift-ink-bounds.png"
    figure.write_image(str(output), scale=1)
    with Image.open(output) as image:
        pixels = np.asarray(image.convert("RGB"))
    red, green, blue = (pixels[:, :, channel] for channel in range(3))
    title = (red > 200) & (blue > 200) & (green < 80)
    headers = (green > 200) & (blue > 200) & (red < 80)
    title_rows = np.flatnonzero(title.any(axis=1))
    header_rows = np.flatnonzero(headers.any(axis=1))
    assert title_rows.size, "Rendered title ink must be measurable"
    assert header_rows.size, "Rendered header ink must be measurable"
    # Test every visible horizontal ink band; no character-count/font-width estimate.
    bands = np.split(title_rows, np.flatnonzero(np.diff(title_rows) > 1) + 1)
    assert len(bands) >= 3, "The main title and subtitle must all be rendered"
    for band in bands:
        columns = np.flatnonzero(title[band].any(axis=0))
        assert 16 <= int(columns.min()) <= int(columns.max()) < 1650 - 16, (
            "Rendered subtitle line reaches a canvas edge: "
            f"rows={int(band.min())}:{int(band.max())}, "
            f"columns={int(columns.min())}:{int(columns.max())}"
        )
    clearance = int(header_rows.min()) - int(title_rows.max())
    assert clearance >= 20, f"Rendered subtitle overlaps panel headers: {clearance}px"
