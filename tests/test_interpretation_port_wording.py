"""Guard the annotation wording ported from Jed's PR #193 (#1262)."""

from __future__ import annotations

import inspect
from pathlib import Path

import pytest
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
def test_heatmap_zero_caption_describes_annotation_evidence(
    annotation_matrix: guild_maps.GuildMatrix,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    renderer: str,
    mark_no_signal: bool,
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
        )
    else:
        guild_maps.write_body_site_panels(
            annotation_matrix,
            output,
            "fixture",
            sites=("Feces",),
            guild_order=("aerobe",),
            mark_no_signal=mark_no_signal,
        )

    assert len(figures) == 1
    title = str(figures[0].layout.title.text)
    assert "no taxon carried" not in title
    assert "an absent trait" not in title
    caption = "Dark grey = sampled, with no recorded annotation signal."
    assert (caption in title) is mark_no_signal
    assert "summed trait-assertion counts" in title
    assert "Directions are curated signature labels." in title
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
    # Keep the existing inference and support displays and saved values intact.
    assert "Benjamini-Hochberg across the six traits" in title
    assert "2/3 studies, q=1" in annotations
    assert list(figure.data[0].x) == [2 / 3]
    assert list(figure.data[0].error_x.array) == [0.95 - 2 / 3]
    assert list(figure.data[0].error_x.arrayminus) == [2 / 3 - 0.1]
