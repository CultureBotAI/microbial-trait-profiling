"""Tufte-plate helpers used by the NCASM plates: paper-and-ink tokens, a plate
scaffold, range frames, direct-label layout, and deterministic saving.

Copied from ``figure_wizard/tufte.py`` in CultureBotAI/figure-wizard (file revision
``fc1197b``, 2026-10-17; BSD 3-Clause License, Copyright (c) 2026, CultureBotAI), which
the original builders imported from a local checkout (#2264). Only the helpers the
30 plates use are kept. The drawing calls are unchanged; the differences are type
annotations and explicit keyword arguments.

BSD 3-Clause notice for this file: redistribution and use in source and binary forms,
with or without modification, are permitted provided that the above copyright notice,
this list of conditions and the disclaimer in figure-wizard's LICENSE are retained.
"""

from __future__ import annotations

import warnings
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Any, Literal

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.transforms import offset_copy
from matplotlib.typing import RcKeyType

__all__ = [
    "ACCENT_CRIMSON",
    "ACCENT_OCHRE",
    "ACCENT_SLATE",
    "BAR_NEUTRAL",
    "BG_CARD",
    "BG_PAPER",
    "INK_BODY",
    "INK_MAIN",
    "INK_MUTED",
    "RULE_LIGHT",
    "SERIF",
    "TINT_SLATE",
    "add_panel_header",
    "apply_range_frame",
    "draw_slopegraph_series",
    "draw_tufte_quartile_rule",
    "make_tufte_plate",
    "relax_vertical_positions",
    "save_plate",
    "setup_tufte_rcparams",
    "tufte_rcparams",
]

# Archival paper and ink
BG_PAPER = "#FCFBF7"
BG_CARD = "#F4F1EA"
INK_MAIN = "#181715"
INK_BODY = "#4A4742"
INK_MUTED = "#78736A"
RULE_LIGHT = "#D6D1C7"

# Semantic accents: one entity per accent, fixed across a document
ACCENT_CRIMSON = "#9E2A2B"
ACCENT_SLATE = "#2B4C6F"
TINT_SLATE = "#EBF1F5"
ACCENT_OCHRE = "#8C5A2B"
BAR_NEUTRAL = "#BFB8AC"

SERIF = ["Palatino", "Palatino Linotype", "Book Antiqua", "Georgia", "Baskerville", "DejaVu Serif"]


def tufte_rcparams() -> dict[RcKeyType, Any]:
    """The plate style as an rcParams dict; ``svg.fonttype`` and ``svg.hashsalt`` are left
    to ``save_plate``."""
    return {
        "figure.facecolor": BG_PAPER,
        "axes.facecolor": BG_PAPER,
        "savefig.facecolor": BG_PAPER,
        "font.family": "serif",
        "font.serif": SERIF,
        "text.color": INK_MAIN,
        "axes.labelcolor": INK_BODY,
        "xtick.color": INK_BODY,
        "ytick.color": INK_BODY,
        "axes.edgecolor": INK_MAIN,
        "axes.linewidth": 0.75,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }


def setup_tufte_rcparams() -> None:
    """Apply the plate style to matplotlib's rcParams.

    The build wraps every plate in ``matplotlib.rc_context``, so the change ends with
    that plate."""
    mpl.rcParams.update(tufte_rcparams())


def make_tufte_plate(
    kicker: str,
    title: str,
    subtitle: str,
    source_note: str,
    plate_id: str,
    figsize: tuple[float, float] = (13.0, 7.6),
) -> Figure:
    """A landscape plate with kicker, headline, italic subtitle, header rule and a
    provenance footer. Panels go in y = 0.08..0.85."""
    setup_tufte_rcparams()
    fig = plt.figure(figsize=figsize, facecolor=BG_PAPER)
    fig.text(0.04, 0.955, kicker.upper(), fontsize=8.8, fontweight="bold", color=ACCENT_CRIMSON)
    fig.text(0.04, 0.918, title, fontsize=15.8, fontweight="bold", color=INK_MAIN)
    fig.text(0.04, 0.886, subtitle, fontsize=9.6, fontstyle="italic", color=INK_BODY)
    fig.add_artist(
        Line2D(
            [0.04, 0.96],
            [0.868, 0.868],
            transform=fig.transFigure,
            color=INK_MAIN,
            linewidth=0.8,
        )
    )
    fig.add_artist(
        Line2D(
            [0.04, 0.96],
            [0.052, 0.052],
            transform=fig.transFigure,
            color=RULE_LIGHT,
            linewidth=0.65,
        )
    )
    fig.text(0.04, 0.026, f"Source: {source_note}", fontsize=7.8, color=INK_MUTED)
    fig.text(0.96, 0.026, plate_id, fontsize=7.8, color=INK_MUTED, ha="right")
    return fig


def add_panel_header(ax: Axes, title: str, subcaption: str, gap_pt: float = 4.0) -> None:
    """Bold panel title above an italic reading guide, both above the axes.

    Offsets are in points, so short axes do not make them collide."""
    sub = ax.annotate(
        subcaption,
        xy=(0, 1),
        xycoords="axes fraction",
        xytext=(0, gap_pt),
        textcoords="offset points",
        fontsize=8.6,
        fontstyle="italic",
        color=INK_BODY,
        va="bottom",
    )
    ax.annotate(
        title,
        xy=(0, 1),
        xycoords=sub,
        xytext=(0, gap_pt),
        textcoords="offset points",
        fontsize=10.8,
        fontweight="bold",
        color=INK_MAIN,
        va="bottom",
    )


def _frame_axis(
    ax: Axes,
    side: Literal["bottom", "left"],
    bounds: tuple[float, float] | None,
    outward_pt: float,
    keep_labels: bool,
) -> None:
    spine = ax.spines[side]
    if side == "bottom":
        if bounds is None:
            spine.set_visible(False)
            ax.tick_params(axis="x", which="both", length=0, labelbottom=keep_labels)
            return
        spine.set_visible(True)
        spine.set_bounds(bounds[0], bounds[1])
        spine.set_position(("outward", outward_pt))
        ax.tick_params(
            axis="x", which="major", length=mpl.rcParams["xtick.major.size"], labelbottom=True
        )
        ax.tick_params(axis="x", which="minor", length=mpl.rcParams["xtick.minor.size"])
        return
    if bounds is None:
        spine.set_visible(False)
        ax.tick_params(axis="y", which="both", length=0, labelleft=keep_labels)
        return
    spine.set_visible(True)
    spine.set_bounds(bounds[0], bounds[1])
    spine.set_position(("outward", outward_pt))
    ax.tick_params(axis="y", which="major", length=mpl.rcParams["ytick.major.size"], labelleft=True)
    ax.tick_params(axis="y", which="minor", length=mpl.rcParams["ytick.minor.size"])


def apply_range_frame(
    ax: Axes,
    x_bounds: tuple[float, float] | None = None,
    y_bounds: tuple[float, float] | None = None,
    outward_pt: float = 5.0,
    keep_labels: bool = False,
) -> None:
    """Clip each visible spine to the data's min and max (a Tufte range frame).

    An axis given no bounds loses its spine and tick marks, and also its tick labels
    unless ``keep_labels`` is set. An axis given bounds gets its spine, ticks and labels
    back, so the call can be repeated."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    _frame_axis(ax, "bottom", x_bounds, outward_pt, keep_labels)
    _frame_axis(ax, "left", y_bounds, outward_pt, keep_labels)


def relax_vertical_positions(
    y_targets: Sequence[float],
    min_gap: float,
    y_min: float | None = None,
    y_max: float | None = None,
) -> list[float]:
    """Label positions at least ``min_gap`` apart, as close as possible to the targets
    (least squares), in the original order and within [y_min, y_max].

    Solved exactly: with q_i = p_i - i*gap the constraint becomes q non-decreasing, which
    pool-adjacent-violators solves; uniform bounds on q are then a clip. If the labels
    cannot fit in the bounds, they are spread evenly over them with a warning."""
    targets = np.asarray(y_targets, dtype=float)
    if min_gap < 0:
        raise ValueError(f"min_gap must be >= 0, got {min_gap}")
    if y_min is not None and y_max is not None and y_min > y_max:
        raise ValueError(f"y_min ({y_min}) is above y_max ({y_max})")
    if np.isnan(targets).any():
        raise ValueError("label targets contain NaN")
    n = len(targets)
    low = -np.inf if y_min is None else y_min
    if n <= 1:
        high_single = np.inf if y_max is None else y_max
        return [float(np.clip(y, low, high_single)) for y in targets]
    order = np.argsort(targets, kind="stable")
    shifted = targets[order] - min_gap * np.arange(n)
    blocks: list[tuple[float, float]] = []  # (mean, weight)
    for value in shifted:
        blocks.append((float(value), 1.0))
        while len(blocks) > 1 and blocks[-2][0] > blocks[-1][0]:
            m2, w2 = blocks.pop()
            m1, w1 = blocks.pop()
            blocks.append(((m1 * w1 + m2 * w2) / (w1 + w2), w1 + w2))
    q = np.concatenate([np.full(int(w), m) for m, w in blocks])
    high = np.inf if y_max is None else y_max - min_gap * (n - 1)
    if y_min is not None and y_max is not None and low > high:
        warnings.warn(
            f"{n} labels need {min_gap * (n - 1):g} but only {y_max - y_min:g} is available; "
            "spreading them evenly",
            stacklevel=2,
        )
        placed = np.linspace(y_min, y_max, n)
    else:
        placed = np.clip(q, low, high) + min_gap * np.arange(n)
    result = np.empty(n)
    result[order] = placed
    return [float(value) for value in result]


def draw_slopegraph_series(
    ax: Axes,
    x_positions: Sequence[float],
    y_values: Sequence[float],
    color: str,
    right_label: str,
    right_y_relaxed: float | None = None,
    linewidth: float = 1.6,
    linestyle: str = "-",
    markersize: float = 5.0,
    is_focal: bool = False,
    label_pad_pt: float = 10.0,
) -> None:
    """One slopegraph trajectory with node dots and a direct right-hand label.

    A label moved by ``relax_vertical_positions`` gets a thin leader line."""
    ax.plot(
        x_positions,
        y_values,
        color=color,
        linewidth=linewidth,
        linestyle=linestyle,
        zorder=4 if is_focal else 3,
    )
    ax.scatter(
        x_positions,
        y_values,
        color=color,
        s=markersize**2,
        edgecolor=BG_PAPER,
        linewidth=0.8,
        zorder=5 if is_focal else 4,
    )
    x_end, y_end = x_positions[-1], y_values[-1]
    y_text = y_end if right_y_relaxed is None else right_y_relaxed
    displaced = abs(y_text - y_end) > 1e-9
    figure = ax.get_figure(root=True)
    if figure is None:
        raise ValueError("the slopegraph axes belong to no figure")
    ax.annotate(
        right_label,
        xy=(x_end, y_end),
        xycoords="data",
        xytext=(x_end, y_text),
        textcoords=offset_copy(ax.transData, fig=figure, x=label_pad_pt, units="points"),
        fontsize=8.4,
        fontweight="bold" if is_focal else "normal",
        color=color if is_focal else INK_MAIN,
        va="center",
        arrowprops=(
            {
                "arrowstyle": "-",
                "color": INK_MUTED,
                "lw": 0.6,
                "shrinkA": 1,
                "shrinkB": markersize / 2 + 1,
            }
            if displaced
            else None
        ),
    )


def draw_tufte_quartile_rule(
    ax: Axes,
    values: Iterable[float],
    y_rule: float,
    unit_suffix: str = "%",
    tick_pt: float = 6.0,
    label_pad_pt: float = 6.0,
    label_width_pt: float = 36.0,
) -> None:
    """A minimal quartile summary ([Min, Q1] · median · [Q3, Max]) along ``y_rule``.

    Ticks mark the exact values. Their labels sit below the rule and are spread sideways
    when they would overlap (``label_width_pt`` is the room one label needs). Constant
    data gets a single label."""
    data = np.asarray(list(values), dtype=float)
    if data.size == 0:
        raise ValueError("draw_tufte_quartile_rule needs at least one value")
    v_min, q1, med, q3, v_max = (float(v) for v in np.percentile(data, [0, 25, 50, 75, 100]))
    ax.plot([v_min, q1], [y_rule, y_rule], color=INK_MUTED, linewidth=0.9)
    ax.plot([q3, v_max], [y_rule, y_rule], color=INK_MUTED, linewidth=0.9)
    ax.plot([q1, q3], [y_rule, y_rule], color=INK_MAIN, linewidth=2.2)
    ax.scatter([med], [y_rule], color=INK_MAIN, s=28, edgecolor=BG_PAPER, linewidth=1.2, zorder=5)
    stats = [(v_min, "Min"), (q1, "Q1"), (med, "Median"), (q3, "Q3"), (v_max, "Max")]
    for val, _ in stats:
        ax.plot(
            [val],
            [y_rule],
            marker="|",
            markersize=2 * tick_pt,
            markeredgewidth=0.8,
            color=INK_MAIN,
            linestyle="none",
        )
    figure = ax.get_figure(root=True)
    if figure is None:
        raise ValueError("the quartile-rule axes belong to no figure")
    below = offset_copy(ax.transData, fig=figure, y=-(tick_pt + label_pad_pt), units="points")
    if v_max == v_min:
        ax.annotate(
            f"All values\n{v_min:.1f}{unit_suffix}",
            xy=(v_min, y_rule),
            xytext=(v_min, y_rule),
            textcoords=below,
            fontsize=7.6,
            color=INK_BODY,
            ha="center",
            va="top",
        )
        return
    x0, x1 = ax.get_xlim()
    width_pt = ax.get_window_extent().width * 72.0 / figure.dpi
    gap = label_width_pt / width_pt * max(x1 - x0, v_max - v_min)
    anchors = relax_vertical_positions([val for val, _ in stats], min_gap=gap)
    for (val, name), x_label in zip(stats, anchors, strict=True):
        ax.annotate(
            f"{name}\n{val:.1f}{unit_suffix}",
            xy=(val, y_rule),
            xytext=(x_label, y_rule),
            textcoords=below,
            fontsize=7.6,
            color=INK_BODY,
            ha="center",
            va="top",
        )


def save_plate(
    fig: Figure,
    stem: Path,
    *,
    formats: Sequence[str] = ("svg", "pdf", "png"),
    dpi: int = 300,
    svg_fonttype: str | None = None,
    salt: str | None = None,
    facecolor: str,
) -> list[Path]:
    """Save ``fig`` as <stem>.<fmt> for each format; returns the paths written.

    The stem keeps any dots, and missing parent directories are created. With ``salt``
    set, SVG ids are seeded from it and no dates are written, so SVG, PDF and PNG output
    does not change from run to run in one environment. PDF and PS embed TrueType fonts.
    rcParams are restored afterwards."""
    stem.parent.mkdir(parents=True, exist_ok=True)
    rc: dict[RcKeyType, Any] = {"pdf.fonttype": 42, "ps.fonttype": 42}
    if svg_fonttype is not None:
        rc["svg.fonttype"] = svg_fonttype
    if salt is not None:
        rc["svg.hashsalt"] = salt
    written: list[Path] = []
    with mpl.rc_context(rc):
        for fmt in (f.lower().lstrip(".") for f in formats):
            path = stem.with_name(f"{stem.name}.{fmt}")
            meta: dict[str, Any] = {}
            if salt is not None:
                if fmt == "svg":
                    meta = {"Date": None}
                elif fmt == "pdf":
                    meta = {"CreationDate": None}
            fig.savefig(
                path,
                dpi=dpi,
                format=fmt,
                metadata=meta or None,
                facecolor=facecolor,
            )
            written.append(path)
    return written
