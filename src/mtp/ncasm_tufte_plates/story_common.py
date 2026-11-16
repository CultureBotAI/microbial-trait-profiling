"""Canvas constants, colours and panel helpers for the ten data-story plots.

Ported from ``tufte_slides_common.py`` in the slide-wizard deck folder
``slides/NCASM-disease-microbiome-traits/`` (sha256 ``06880fae25fb``, 2026-10-21; #2264).
Only what the data-story plots use is kept: the 24-slide deck's own helpers
(``init_slide``, the table, bar and lollipop drawers, ``save_and_verify_slide``) and the
``tufte_slides/`` folder the module used to create on import are not part of the plates.
"""

from __future__ import annotations

from matplotlib.axes import Axes
from matplotlib.figure import Figure

from mtp.ncasm_tufte_plates.plates_common import PLATE_HEIGHT_IN, PLATE_WIDTH_IN
from mtp.ncasm_tufte_plates.tufte_style import (
    ACCENT_CRIMSON,
    ACCENT_OCHRE,
    ACCENT_SLATE,
    BAR_NEUTRAL,
    BG_PAPER,
    INK_BODY,
    INK_MAIN,
    RULE_LIGHT,
)

# The same canvas as plates I-XX; see ``plates_common.PLATE_WIDTH_IN`` for why the width
# is 13.33 in and not the 13.333 in the builders ask for.
SLIDE_W = PLATE_WIDTH_IN
SLIDE_H = PLATE_HEIGHT_IN
SLIDE_DPI = 200

SLIDE_SIZES = {
    "kicker": 12.0,
    "title": 21.5,
    "subtitle": 13.8,
    "panel_title": 15.8,
    "axis_label": 13.5,
    "tick_label": 13.0,
    "footer": 11.0,
}

INK = INK_MAIN
MUTED = INK_BODY
FAINT = BAR_NEUTRAL
PAPER = BG_PAPER
RULE = RULE_LIGHT
ACCENT = ACCENT_CRIMSON  # #9E2A2B crimson
SECONDARY = ACCENT_SLATE  # #2B4C6F slate blue
TERTIARY = "#386641"  # forest green
HIGHLIGHT = ACCENT_OCHRE  # #8C5A2B ochre


def despine(
    ax: Axes,
    top: bool = True,
    right: bool = True,
    left: bool = False,
    bottom: bool = False,
) -> None:
    """Hide the named spines and show the others."""
    for side, hide in (("top", top), ("right", right), ("left", left), ("bottom", bottom)):
        if side in ax.spines:
            ax.spines[side].set_visible(not hide)


def add_panel(
    fig: Figure,
    rect: tuple[float, float, float, float],
    heading: str | None = None,
    subnote: str | None = None,
) -> Axes:
    """Add an axes in figure coordinates (left, bottom, width, height) with an optional
    heading above it. The builders pass a ``subnote`` that was never drawn."""
    left, bottom, _width, height = rect
    _ = subnote  # Subnotes below panel titles intentionally omitted
    if heading:
        fig.text(
            left,
            bottom + height + 0.020,
            heading,
            fontsize=SLIDE_SIZES["panel_title"],
            fontweight="bold",
            color=INK,
            ha="left",
            va="bottom",
        )
    return fig.add_axes(rect, facecolor=PAPER)
