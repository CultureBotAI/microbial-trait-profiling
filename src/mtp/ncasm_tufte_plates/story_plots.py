"""The ten "gray-dot universe" data-story plots, which the build writes as plates 21-30.

Ported from ``build_tufte_10_data_story_plots.py`` in the slide-wizard deck folder
``slides/NCASM-disease-microbiome-traits/`` (sha256 ``8d24c320dffb``, 2026-10-22; #2264).
The drawing calls are unchanged. The differences: each plot takes its output folder as an
argument (the module no longer creates a folder on import or edits ``sys.path``), SVG and
PDF files are written without dates and with a fixed id seed, and the standalone
ten-page PDF and audit JSON of the original ``main()`` are left to
``mtp.ncasm_tufte_plates.build``.

Every plot here draws its points from a seeded ``numpy`` random generator and types its
numbers as literals (#775, #776). Each therefore saves with ``random_marks=True``, which
stamps ``RANDOM_MARKS_BANNER`` across the plate, and is listed in ``RANDOM_MARK_PLOTS``,
which ``mtp.ncasm_tufte_plates.build`` refuses to draw unless asked. A plot rebuilt from
run tables drops both. ``tests/test_ncasm_tufte_plates.py`` fails when a function that
calls a random generator is missing from either.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import TypedDict

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse, FancyBboxPatch
from matplotlib.text import Text
from matplotlib.transforms import Bbox

from mtp.ncasm_tufte_plates.data_status import RANDOM_MARKS_BANNER
from mtp.ncasm_tufte_plates.story_common import (
    ACCENT,
    FAINT,
    HIGHLIGHT,
    INK,
    MUTED,
    PAPER,
    RULE,
    SECONDARY,
    SLIDE_DPI,
    SLIDE_H,
    SLIDE_SIZES,
    SLIDE_W,
    TERTIARY,
    add_panel,
    despine,
)
from mtp.ncasm_tufte_plates.tufte_style import setup_tufte_rcparams

TOTAL_PLOTS = 10


class StoryPlotInfo(TypedDict):
    """What one data-story plot wrote, with the layout check's findings."""

    stem: str
    png_path: str
    svg_path: str
    pdf_path: str
    out_of_bounds: list[str]
    overlaps: list[tuple[str, str]]


def init_story_plot(
    plot_num: int,
    archetype_tag: str,
    title: str,
    subtitle: str,
    footer_source: str,
) -> Figure:
    """Create a 13.333x7.5 in 16:9 Tufte data-story canvas."""
    setup_tufte_rcparams()
    mpl.rcParams.update(
        {
            "font.size": SLIDE_SIZES["tick_label"],
            "axes.titlesize": SLIDE_SIZES["panel_title"],
            "axes.labelsize": SLIDE_SIZES["axis_label"],
            "xtick.labelsize": SLIDE_SIZES["tick_label"],
            "ytick.labelsize": SLIDE_SIZES["tick_label"],
        }
    )

    fig = plt.figure(figsize=(SLIDE_W, SLIDE_H), dpi=SLIDE_DPI, facecolor=PAPER)

    kicker = f"DATA STORY PLOT {plot_num:02d} OF {TOTAL_PLOTS:02d}  \u00b7  {archetype_tag.upper()}"
    fig.text(
        0.04,
        0.945,
        kicker,
        fontsize=SLIDE_SIZES["kicker"],
        fontweight="bold",
        color=ACCENT,
        ha="left",
        va="top",
    )
    fig.text(
        0.04,
        0.908,
        title,
        fontsize=SLIDE_SIZES["title"],
        fontweight="bold",
        color=INK,
        ha="left",
        va="top",
    )
    fig.text(
        0.04,
        0.855,
        subtitle,
        fontsize=SLIDE_SIZES["subtitle"],
        fontstyle="italic",
        color=MUTED,
        ha="left",
        va="top",
    )

    fig.add_artist(
        Line2D(
            [0.04, 0.96],
            [0.822, 0.822],
            transform=fig.transFigure,
            color=INK,
            linewidth=1.0,
        )
    )
    fig.add_artist(
        Line2D(
            [0.04, 0.96],
            [0.058, 0.058],
            transform=fig.transFigure,
            color=RULE,
            linewidth=0.7,
        )
    )
    fig.text(
        0.04,
        0.030,
        footer_source,
        fontsize=SLIDE_SIZES["footer"],
        color=MUTED,
        ha="left",
        va="center",
    )
    fig.text(
        0.96,
        0.030,
        f"Tufte Data Story ({plot_num}/{TOTAL_PLOTS})",
        fontsize=SLIDE_SIZES["footer"],
        color=MUTED,
        ha="right",
        va="center",
    )
    return fig


def stamp_random_marks(fig: Figure) -> None:
    """Say across the top of the plate that its points are random draws (#775)."""
    fig.text(
        0.5,
        0.992,
        RANDOM_MARKS_BANNER,
        fontsize=11.5,
        fontweight="bold",
        color=PAPER,
        ha="center",
        va="top",
        bbox={"facecolor": ACCENT, "edgecolor": "none", "boxstyle": "square,pad=0.35"},
        gid="random-marks-banner",
        zorder=10,
    )


def save_and_verify_story_plot(
    fig: Figure, stem: str, out_dir: Path, *, random_marks: bool
) -> StoryPlotInfo:
    """Save PNG/SVG/PDF into ``out_dir`` and report overlapping and out-of-bounds text.

    ``random_marks`` is true for a plot whose points come from a random generator; the
    plate is then stamped with ``RANDOM_MARKS_BANNER``."""
    if random_marks:
        stamp_random_marks(fig)
    out_dir.mkdir(parents=True, exist_ok=True)
    png_path = out_dir / f"{stem}.png"
    svg_path = out_dir / f"{stem}.svg"
    pdf_path = out_dir / f"{stem}.pdf"

    fig.canvas.draw()
    fig_bbox = fig.bbox

    text_boxes: list[tuple[str, Bbox]] = []
    out_of_bounds: list[str] = []
    overlaps: list[tuple[str, str]] = []

    all_texts: list[Text] = list(fig.texts)
    for ax in fig.axes:
        all_texts.extend(ax.texts)

    for t in all_texts:
        if not t.get_visible():
            continue
        s = t.get_text().strip()
        if not s:
            continue
        bb = t.get_window_extent()
        if (
            bb.x0 < fig_bbox.x0 - 4
            or bb.y0 < fig_bbox.y0 - 4
            or bb.x1 > fig_bbox.x1 + 4
            or bb.y1 > fig_bbox.y1 + 4
        ):
            out_of_bounds.append(f"{s[:40]} ({bb.x0:.1f},{bb.y0:.1f},{bb.x1:.1f},{bb.y1:.1f})")
        text_boxes.append((s, bb))

    for i in range(len(text_boxes)):
        s1, b1 = text_boxes[i]
        b1_s = b1.expanded(0.90, 0.82)
        for j in range(i + 1, len(text_boxes)):
            s2, b2 = text_boxes[j]
            b2_s = b2.expanded(0.90, 0.82)
            if b1_s.overlaps(b2_s):
                overlaps.append((s1[:35], s2[:35]))

    # The original wrote the SVG and PDF with random ids and the current date. The id seed
    # and the absent dates make all three files repeatable; the drawing is the same.
    with mpl.rc_context({"svg.hashsalt": stem}):
        fig.savefig(png_path, dpi=SLIDE_DPI, facecolor=PAPER)
        fig.savefig(svg_path, facecolor=PAPER, metadata={"Date": None})
        fig.savefig(pdf_path, facecolor=PAPER, metadata={"CreationDate": None})
    plt.close(fig)

    return {
        "stem": stem,
        "png_path": str(png_path),
        "svg_path": str(svg_path),
        "pdf_path": str(pdf_path),
        "out_of_bounds": out_of_bounds,
        "overlaps": overlaps,
    }


def draw_story_callout(
    ax: Axes,
    x: float,
    y: float,
    w: float,
    h: float,
    big_num: str,
    label: str,
    sub: str = "",
    color: str = ACCENT,
) -> None:
    """Draw a callout card with stacked hero number, label and subnote."""
    box = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.012,rounding_size=0.02",
        facecolor="#F4F3EE",
        edgecolor=color,
        linewidth=1.3,
    )
    ax.add_patch(box)
    ax.text(
        x + 0.05 * w,
        y + 0.74 * h,
        big_num,
        fontsize=20.0,
        fontweight="bold",
        color=color,
        ha="left",
        va="center",
    )
    ax.text(
        x + 0.05 * w,
        y + 0.42 * h,
        label,
        fontsize=12.6,
        fontweight="bold",
        color=INK,
        ha="left",
        va="center",
    )
    if sub:
        ax.text(
            x + 0.05 * w,
            y + 0.15 * h,
            sub,
            fontsize=11.2,
            color=MUTED,
            ha="left",
            va="center",
        )


def draw_quartile_rule_h(
    ax: Axes,
    q_min: float,
    q1: float,
    med: float,
    q3: float,
    q_max: float,
    y: float,
    color: str = INK,
    lw_whisker: float = 1.2,
    lw_iqr: float = 4.0,
    dot_size: float = 65.0,
    zorder: int = 4,
) -> None:
    """Draw a horizontal quartile rule: whisker, thick [q1..q3] and a hollow median dot."""
    ax.plot([q_min, q_max], [y, y], color=color, lw=lw_whisker, zorder=zorder)
    ax.plot([q1, q3], [y, y], color=color, lw=lw_iqr, zorder=zorder + 1)
    ax.scatter([med], [y], color=PAPER, edgecolor=color, s=dot_size, lw=1.8, zorder=zorder + 2)


# =============================================================================
# PLOT 01: 2D Joint Effect-vs-Consistency Cloud with Marginal Histograms
# Inspired by: tufte_023.png, tufte_101.jpg, tufte_012.png
# =============================================================================
def plot_01_joint_effect_consistency(out_dir: Path) -> StoryPlotInfo:
    fig = init_story_plot(
        1,
        "Tufte Joint Scatter + Marginal Histograms (cf. tufte_023, tufte_101)",
        "Six Global Trait Signals Emerge From the Cloud of 239 Evaluable Traits",
        "Plotting all 239 evaluable traits across 935 contrasts (\u0394z) and 548 studies (% "
        "majority votes) reveals two opposing tails.",
        "Source: Arm 1 contrast shift (\u0394z, n = 935) vs. Arm 3 study majority vote (n = 548) "
        "across 239 evaluable traits (of 317 curated).",
    )

    # Central 2D scatter panel + top/right marginal histograms
    ax = fig.add_axes((0.08, 0.14, 0.55, 0.45), facecolor=PAPER)
    ax_top = fig.add_axes((0.08, 0.595, 0.55, 0.085), facecolor=PAPER)
    ax_right = fig.add_axes((0.635, 0.14, 0.075, 0.45), facecolor=PAPER)
    fig.text(
        0.08,
        0.735,
        "A. 239-Trait Joint Cloud: Contrast Shift (\u0394z) vs. Study Vote (%)",
        fontsize=SLIDE_SIZES["panel_title"],
        fontweight="bold",
        color=INK,
        ha="left",
        va="top",
    )
    fig.text(
        0.08,
        0.703,
        "Gray dots = 233 background traits; top/right step histograms show marginal distributions",
        fontsize=12.0,
        color=MUTED,
        ha="left",
        va="top",
    )

    rng = np.random.default_rng(101)
    # Generate realistic correlated background of 233 traits centered near (0, 50%)
    cov = [[0.0022, 0.095], [0.095, 7.2]]
    pts = rng.multivariate_normal([0.002, 50.1], cov, size=233)
    bg_x = np.clip(pts[:, 0], -0.11, 0.105)
    bg_y = np.clip(pts[:, 1], 42.5, 55.2)

    ax.axvline(0.0, color=RULE, lw=1.1, ls="--", zorder=1)
    ax.axhline(50.0, color=RULE, lw=1.1, ls="--", zorder=1)

    ax.scatter(
        bg_x, bg_y, color=FAINT, edgecolor="#8E887E", linewidth=0.4, s=46, alpha=0.72, zorder=2
    )

    # Rug ticks along bottom & left range frames (cf. tufte_012.png)
    ax.scatter(bg_x, np.full_like(bg_x, 36.5), marker="|", color=MUTED, s=28, alpha=0.45, zorder=2)
    ax.scatter(
        np.full_like(bg_y, -0.214), bg_y, marker="_", color=MUTED, s=28, alpha=0.45, zorder=2
    )

    # Focal 6 traits: (dz, study_pct_inc, label, color, tx, ty, ha)
    focal = [
        (
            -0.17,
            40.6,
            "#1 Cellulose catabolism\n(\u0394z = \u22120.25 \u00b7 40.6% Inc / 59.4% Dec \u00b7 q "
            "= 1.5e-7)",
            SECONDARY,
            -0.04,
            39.2,
            "left",
        ),
        (
            -0.15,
            46.4,
            "Obligate anaerobe\n(\u0394z = \u22120.19 \u00b7 46.4% Inc / 53.6% Dec)",
            HIGHLIGHT,
            -0.05,
            44.2,
            "left",
        ),
        (
            0.16,
            57.2,
            "Facultative anaerobe\n(\u0394z = +0.21 \u00b7 57.2% Inc \u00b7 q = 3.2e-4)",
            ACCENT,
            0.03,
            59.8,
            "right",
        ),
        (
            0.13,
            58.0,
            "Halophilic (\u0394z = +0.14 \u00b7 58.0% Inc \u00b7 q = 0.0022)",
            ACCENT,
            -0.02,
            57.6,
            "right",
        ),
        (
            0.14,
            56.2,
            "Coccus (\u0394z = +0.16 \u00b7 56.2% Inc \u00b7 q = 0.0047)",
            ACCENT,
            0.05,
            52.0,
            "right",
        ),
        (
            0.10,
            58.8,
            "Arsenate red. (\u0394z = +0.25 \u00b7 58.8% Inc)",
            TERTIARY,
            -0.01,
            54.6,
            "right",
        ),
    ]
    for fx, fy, flab, fcol, tx, ty, ha in focal:
        ax.scatter([fx], [fy], color=fcol, s=125, edgecolor=PAPER, lw=1.3, zorder=5)
        ax.plot([fx, tx], [fy, ty], color=fcol, lw=1.3, zorder=4)
        ax.text(tx, ty, flab, fontsize=11.8, fontweight="bold", color=fcol, ha=ha, va="center")

    ax.set_xlim(-0.22, 0.21)
    ax.set_ylim(36.0, 62.0)
    ax.set_xticks([-0.15, -0.08, 0.0, 0.08, 0.15])
    ax.set_xticklabels(["\u22120.15", "\u22120.08", "0.00", "+0.08", "+0.15"])
    ax.set_yticks([40, 45, 50, 55, 60])
    ax.set_yticklabels(["40%", "45%", "50%", "55%", "60%"])
    ax.set_xlabel(
        "Arm 1 Mean Within-Contrast Shift \u0394z (935 Contrasts)",
        fontsize=SLIDE_SIZES["axis_label"],
    )
    ax.set_ylabel("Arm 3 Study Majority Vote (% Increased)", fontsize=SLIDE_SIZES["axis_label"])
    despine(ax, top=True, right=True, left=False, bottom=False)

    # Top marginal histogram + Tufte quartile rule
    all_x = np.concatenate([bg_x, np.array([f[0] for f in focal])])
    bins_x = np.linspace(-0.20, 0.20, 25)
    ax_top.hist(
        all_x,
        bins=bins_x.tolist(),
        color=FAINT,
        edgecolor=INK,
        lw=0.8,
        histtype="stepfilled",
        alpha=0.65,
    )
    ax_top.axvline(0.0, color=RULE, lw=1.0, ls="--")
    ax_top.scatter(
        [-0.17, 0.16, 0.13, 0.14],
        [12, 12, 12, 12],
        color=[SECONDARY, ACCENT, ACCENT, ACCENT],
        s=45,
        zorder=5,
    )
    ax_top.set_xlim(-0.22, 0.21)
    ax_top.axis("off")

    # Right marginal histogram
    all_y = np.concatenate([bg_y, np.array([f[1] for f in focal])])
    bins_y = np.linspace(36.0, 61.0, 22)
    ax_right.hist(
        all_y,
        bins=bins_y.tolist(),
        orientation="horizontal",
        color=FAINT,
        edgecolor=INK,
        lw=0.8,
        histtype="stepfilled",
        alpha=0.65,
    )
    ax_right.axhline(50.0, color=RULE, lw=1.0, ls="--")
    ax_right.scatter(
        [14, 14, 14], [40.6, 58.0, 57.2], color=[SECONDARY, ACCENT, ACCENT], s=45, zorder=5
    )
    ax_right.set_ylim(36.0, 62.0)
    ax_right.axis("off")

    # Right-hand interpretation panel
    ax_info = add_panel(
        fig,
        (0.74, 0.15, 0.22, 0.58),
        "B. Traitome Summary",
        "Full denominator benchmark",
    )
    ax_info.set_xlim(0, 1)
    ax_info.set_ylim(0, 1)
    ax_info.axis("off")
    draw_story_callout(
        ax_info,
        0.0,
        0.70,
        1.0,
        0.25,
        "239 / 317",
        "Evaluable Traits",
        "233 stay in the gray null cloud",
        SECONDARY,
    )
    draw_story_callout(
        ax_info,
        0.0,
        0.37,
        1.0,
        0.25,
        "r = 0.84",
        "Cross-Arm Concordance",
        "Contrast \u0394z vs. Study Vote %",
        ACCENT,
    )
    draw_story_callout(
        ax_info,
        0.0,
        0.04,
        1.0,
        0.25,
        "q \u2264 0.05",
        "6 Focal FDR Traits",
        "Cellulose #1 (q = 1.5e-7)",
        TERTIARY,
    )

    return save_and_verify_story_plot(
        fig, "plot_01_joint_effect_consistency_cloud", out_dir, random_marks=True
    )


# =============================================================================
# PLOT 02: 2D Paired Contrast Parity Cloud (y = x) + Orthogonal Difference Strip
# Inspired by: tufte_006.jpg, tufte_058.jpg, tufte_012.png
# =============================================================================
def plot_02_paired_parity_variance(out_dir: Path) -> StoryPlotInfo:
    fig = init_story_plot(
        2,
        "Paired Parity Cloud + Orthogonal Residual Strip (cf. tufte_006, tufte_058)",
        "Why Within-Contrast Pairing Is Essential: 51\u00d7 Variance Collapse Across 935 Contrasts",
        "Plotting raw Increased vs. Decreased signature densities reveals huge study-level spread "
        "along y = x that cancels upon pairing.",
        "Source: 935 Tier 1+2 paired contrasts (689 non-tied for cellulose catabolism: 413 Dec > "
        "Inc vs. 276 Inc > Dec, p = 4.6\u00d710\u207b\u00b9\u2070).",
    )

    ax_a = add_panel(
        fig,
        (0.07, 0.15, 0.43, 0.56),
        "A. 935 Contrasts in 2D Signature Space (f_inc vs. f_dec)",
        "Points above the y = x parity line have higher trait density in Decreased taxa",
    )
    rng = np.random.default_rng(202)
    # 413 contrasts where Dec > Inc, 276 where Inc > Dec, 246 tied (mostly near 0)
    base_dec = rng.beta(1.8, 5.2, size=413) * 0.58 + 0.04
    inc_when_dec = base_dec * rng.uniform(0.15, 0.82, size=413)

    base_inc = rng.beta(1.6, 5.5, size=276) * 0.52 + 0.03
    dec_when_inc = base_inc * rng.uniform(0.18, 0.85, size=276)

    tied_v = rng.beta(0.8, 8.0, size=246) * 0.30

    # Plot parity diagonal
    ax_a.plot([0, 0.62], [0, 0.62], color=INK, lw=1.4, ls="--", zorder=3)
    ax_a.scatter(
        tied_v, tied_v + rng.normal(0, 0.006, size=246), color=FAINT, s=22, alpha=0.45, zorder=1
    )
    ax_a.scatter(base_inc, dec_when_inc, color="#B89B9B", s=28, alpha=0.55, zorder=2)
    ax_a.scatter(inc_when_dec, base_dec, color=SECONDARY, s=30, alpha=0.55, zorder=2)

    # Centroid callout
    mean_inc = 0.118
    mean_dec = 0.182
    ax_a.scatter([mean_inc], [mean_dec], color=ACCENT, s=150, edgecolor=PAPER, lw=1.5, zorder=6)
    ax_a.plot([mean_inc, 0.34], [mean_dec, 0.53], color=ACCENT, lw=1.4, zorder=5)
    ax_a.text(
        0.35,
        0.53,
        "Mean Contrast Centroid\n(f_dec = 18.2% vs f_inc = 11.8%)",
        fontsize=12.2,
        fontweight="bold",
        color=ACCENT,
        ha="left",
        va="center",
    )
    ax_a.text(
        0.04,
        0.56,
        "Above Parity (Dec > Inc):\n413 contrasts (59.9%)\np = 4.6\u00d710\u207b\u00b9\u2070, q = "
        "1.5\u00d710\u207b\u2077",
        fontsize=12.2,
        fontweight="bold",
        color=SECONDARY,
        ha="left",
        va="top",
    )
    ax_a.text(
        0.35,
        0.07,
        "Below Parity (Inc > Dec):\n276 contrasts (40.1%)",
        fontsize=12.0,
        color=MUTED,
        ha="left",
        va="bottom",
    )

    ax_a.set_xlim(-0.02, 0.64)
    ax_a.set_ylim(-0.02, 0.64)
    ax_a.set_xticks([0.0, 0.15, 0.30, 0.45, 0.60])
    ax_a.set_xticklabels(["0%", "15%", "30%", "45%", "60%"])
    ax_a.set_yticks([0.0, 0.15, 0.30, 0.45, 0.60])
    ax_a.set_yticklabels(["0%", "15%", "30%", "45%", "60%"])
    ax_a.set_xlabel("Increased Signature Trait Density (f_inc)", fontsize=SLIDE_SIZES["axis_label"])
    ax_a.set_ylabel("Decreased Signature Trait Density (f_dec)", fontsize=SLIDE_SIZES["axis_label"])
    despine(ax_a, top=True, right=True)

    # Panel B: Unpaired vs Paired Variance Strip Cloud
    ax_b = add_panel(
        fig,
        (0.56, 0.15, 0.40, 0.56),
        "B. Study Batch Variance Collapse Upon Pairing",
        "Gray dots = 548 study means before vs. after within-contrast differencing",
    )
    despine(ax_b, top=True, right=True, left=True, bottom=False)
    ax_b.yaxis.set_visible(False)
    ax_b.axvline(0.0, color=INK, lw=1.1, ls="--", zorder=1)

    unpaired_dev = rng.normal(0.04, 0.21, size=548)
    unpaired_dev = np.clip(unpaired_dev, -0.46, 0.48)
    y_unp = rng.uniform(0.62, 0.82, size=548)
    ax_b.scatter(unpaired_dev, y_unp, color=FAINT, s=26, alpha=0.65, zorder=2)
    draw_quartile_rule_h(ax_b, -0.42, -0.11, 0.03, 0.18, 0.45, y=0.55, color=MUTED)
    ax_b.text(
        -0.46,
        0.90,
        "1. Unpaired Signature Density (1,870 sigs)  \u00b7  Study R\u00b2 = 9.01%",
        fontsize=12.0,
        fontweight="bold",
        color=MUTED,
        ha="left",
        va="center",
    )

    paired_diff = np.concatenate([inc_when_dec - base_dec, base_inc - dec_when_inc])
    paired_diff = np.clip(paired_diff, -0.46, 0.42)
    y_prd = rng.uniform(0.14, 0.34, size=len(paired_diff))
    ax_b.scatter(paired_diff, y_prd, color=SECONDARY, s=24, alpha=0.45, zorder=2)
    draw_quartile_rule_h(ax_b, -0.38, -0.14, -0.064, 0.02, 0.28, y=0.06, color=ACCENT)
    ax_b.text(
        -0.46,
        0.42,
        "2. Paired Shift \u0394 = f_inc \u2212 f_dec  \u00b7  Study R\u00b2 = 0.18% (51\u00d7 "
        "drop)",
        fontsize=12.0,
        fontweight="bold",
        color=ACCENT,
        ha="left",
        va="center",
    )
    ax_b.text(
        -0.064,
        -0.05,
        "Median \u0394 = \u22126.4 pp\n(413 Dec : 276 Inc)",
        fontsize=11.8,
        fontweight="bold",
        color=ACCENT,
        ha="center",
        va="top",
    )

    ax_b.set_xlim(-0.50, 0.52)
    ax_b.set_ylim(-0.18, 0.98)
    ax_b.set_xticks([-0.4, -0.2, 0.0, 0.2, 0.4])
    ax_b.set_xticklabels(["\u221240 pp", "\u221220 pp", "0 pp", "+20 pp", "+40 pp"])
    ax_b.set_xlabel(
        "Centered Trait Density / Within-Contrast Shift \u0394 (pp)",
        fontsize=SLIDE_SIZES["axis_label"],
    )

    return save_and_verify_story_plot(
        fig, "plot_02_paired_parity_variance_cloud", out_dir, random_marks=True
    )


# =============================================================================
# PLOT 03: Taxonomic Rank Occurrence Landscape (Mapped Gray/Slate vs Unmapped Rings)
# Inspired by: tufte_058.jpg, tufte_063.png, tufte_071.png
# =============================================================================
def plot_03_rank_occurrence_landscape(out_dir: Path) -> StoryPlotInfo:
    fig = init_story_plot(
        3,
        "Multi-Stratum Occurrence Cloud + Quartile Rules (cf. tufte_063, tufte_071)",
        "Why Occurrence Coverage Reaches 91.3% Despite 39.8% Unmapped Singleton Species",
        "High-frequency genera sit on the right and are 97.4% annotated; unmapped taxa (hollow "
        "crimson rings) cluster at 1\u20132 occurrences.",
        "Source: 3,221 Tier 1+2 unique TaxIDs (19,388 occurrences) across 935 contrasts: 69.9% "
        "unique / 91.3% occurrence-weighted coverage.",
    )

    ax = add_panel(
        fig,
        (0.24, 0.15, 0.72, 0.56),
        "Contrast Occurrences per Taxon Across Six Ranks (Log10 Scale \u00b7 N = 3,221 TaxIDs)",
        "Solid gray/slate dots = KG-Microbe trait-annotated taxa (2,251); hollow crimson rings = "
        "unmapped taxa (mostly n = 1\u20132 species)",
    )
    despine(ax, top=True, right=True, left=True, bottom=False)
    ax.yaxis.set_visible(False)

    ranks = [
        ("Class (n = 50)", 84.0, 98.7, 0.2, 0.9, 1.45, 1.95, 2.55, 42, 8),
        ("Phylum (n = 45)", 86.7, 98.4, 0.3, 1.0, 1.55, 2.05, 2.62, 39, 6),
        ("Order (n = 99)", 79.8, 96.9, 0.0, 0.6, 1.15, 1.68, 2.38, 79, 20),
        ("Genus (n = 842)", 85.7, 97.4, 0.0, 0.48, 0.95, 1.52, 2.48, 260, 50),
        ("Family (n = 202)", 86.6, 96.9, 0.0, 0.45, 0.90, 1.48, 2.32, 160, 27),
        ("Species (n = 1,983)", 60.2, 77.7, 0.0, 0.0, 0.30, 0.78, 1.95, 280, 180),
    ]

    rng = np.random.default_rng(303)
    for idx, (rname, u_cov, o_cov, qmin, q1, med, q3, qmax, n_map_pts, n_unmap_pts) in enumerate(
        ranks
    ):
        y0 = float(len(ranks) - 1 - idx)
        # Mapped points shifted toward higher occurrences
        map_x = np.clip(rng.normal(med + 0.15, 0.48, size=n_map_pts), qmin, qmax)
        map_y = y0 + rng.uniform(-0.18, 0.18, size=n_map_pts)
        ax.scatter(
            map_x, map_y, color=FAINT, edgecolor=SECONDARY, lw=0.35, s=26, alpha=0.68, zorder=2
        )

        # Unmapped points clustered near 1-2 occurrences (log10 = 0.0 to 0.45)
        unmap_x = np.clip(rng.exponential(0.22, size=n_unmap_pts), 0.0, 0.95)
        unmap_y = y0 + rng.uniform(-0.20, 0.20, size=n_unmap_pts)
        ax.scatter(
            unmap_x, unmap_y, facecolor="none", edgecolor=ACCENT, lw=1.0, s=34, alpha=0.75, zorder=3
        )

        # Quartile rule underneath each strip
        draw_quartile_rule_h(
            ax, qmin, q1, med, q3, qmax, y=y0 - 0.28, color=INK, lw_iqr=3.2, dot_size=45
        )

        # Left label
        ax.text(
            -0.12,
            y0,
            f"{rname}\n{u_cov:.1f}% uniq \u00b7 {o_cov:.1f}% occ",
            fontsize=12.0,
            fontweight="bold" if "Genus" in rname or "Species" in rname else "normal",
            color=INK,
            ha="right",
            va="center",
            clip_on=False,
        )

    # Call out high-occurrence mapped genera on Genus row (y = 2.0)
    apex_genera = [
        (2.46, 2.0, "Streptococcus (192 exps)", SECONDARY, 2.22, 2.42),
        (2.36, 2.0, "Faecalibacterium (186)", ACCENT, 1.72, 2.42),
        (2.31, 2.0, "Prevotella (183)", HIGHLIGHT, 2.48, 1.62),
    ]
    for gx, gy, glab, gcol, tx, ty in apex_genera:
        ax.scatter([gx], [gy], color=gcol, s=95, edgecolor=PAPER, lw=1.2, zorder=6)
        ax.plot([gx, tx], [gy, ty], color=gcol, lw=1.2, zorder=5)
        ax.text(
            tx, ty, glab, fontsize=11.5, fontweight="bold", color=gcol, ha="center", va="center"
        )

    # Call out singleton unmapped species cluster on Species row (y = 0.0)
    ax.text(
        1.18,
        0.02,
        "789 unmapped species (39.8%)\ncluster at 1\u20132 occurrences (8.7% of total mass)",
        fontsize=11.8,
        fontweight="bold",
        color=ACCENT,
        ha="left",
        va="center",
    )

    ax.set_xlim(-0.08, 2.85)
    ax.set_ylim(-0.55, 5.55)
    ax.set_xticks([0.0, 0.477, 1.0, 1.477, 2.0, 2.477])
    ax.set_xticklabels(["1 occ", "3 occ", "10 occ", "30 occ", "100 occ", "300 occ"])
    ax.set_xlabel(
        "Contrast Occurrences per Unique Taxon Across 935 Contrasts (Log Scale)",
        fontsize=SLIDE_SIZES["axis_label"],
    )

    return save_and_verify_story_plot(
        fig, "plot_03_rank_occurrence_bubble_rug", out_dir, random_marks=True
    )


# =============================================================================
# PLOT 04: Multi-Stratum Gray-Dot Raincloud & Quartile Strips Across 7 Major Conditions
# Inspired by: tufte_021.png, tufte_059.jpg, tufte_071.png, tufte_130.webp
# =============================================================================
def plot_04_seven_conditions_logratio_strips(out_dir: Path) -> StoryPlotInfo:
    fig = init_story_plot(
        4,
        "Multi-Stratum Contrast Cloud + Tufte Quartile Strips (cf. tufte_059, tufte_071)",
        "Every Major Condition Shifts Negative on the Paired Anaerobe / Facultative Log-Ratio",
        "Plotting all 388 individual contrasts across the Seven Major Conditions shows consistent "
        "anaerobic depletion without bar-chart aggregation.",
        "Source: Within-contrast paired log-ratio ln(Obligate Anaerobe / Facultative) across 388 "
        "contrasts in the Seven Major Conditions (548 studies total).",
    )

    ax = add_panel(
        fig,
        (0.28, 0.15, 0.66, 0.56),
        "Within-Contrast Paired Log-Ratio ln(Obligate Anaerobe / Facultative) Across 388 Contrasts",
        "Each dot = 1 contrast; horizontal rule = Q1 \u2013 [hollow Median] \u2013 Q3; dashed "
        "line = 0.0 balanced parity",
    )
    despine(ax, top=True, right=True, left=True, bottom=False)
    ax.yaxis.set_visible(False)
    ax.axvline(0.0, color=INK, lw=1.2, ls="--", zorder=2)

    conditions = [
        ("COVID-19 (Acute viral \u00b7 n = 29)", 29, -0.34, -0.78, -0.04, 63.8, ACCENT),
        ("Parkinson's Disease (Neurological \u00b7 n = 82)", 82, -0.31, -0.72, 0.02, 61.9, ACCENT),
        ("Inflammatory Bowel Disease (GI \u00b7 n = 24)", 24, -0.46, -0.92, -0.10, 67.9, ACCENT),
        ("Colorectal Cancer (Neoplastic \u00b7 n = 65)", 65, -0.29, -0.68, 0.05, 60.7, SECONDARY),
        (
            "Obesity / Metabolic (Metabolic \u00b7 n = 110)",
            110,
            -0.24,
            -0.60,
            0.08,
            57.9,
            SECONDARY,
        ),
        (
            "HIV Infection (Immunodeficiency \u00b7 n = 30)",
            30,
            -0.38,
            -0.82,
            -0.02,
            65.5,
            SECONDARY,
        ),
        ("Irritable Bowel Syndrome (GI \u00b7 n = 48)", 48, -0.42, -0.88, -0.05, 64.7, HIGHLIGHT),
    ]

    rng = np.random.default_rng(404)
    for idx, (cname, n_c, med, q1, q3, pct_down, col) in enumerate(conditions):
        y0 = float(len(conditions) - 1 - idx)
        pts = rng.normal(med, 0.52, size=n_c)
        pts = np.clip(pts, -1.75, 1.25)
        y_jit = y0 + rng.uniform(-0.18, 0.18, size=n_c)

        # Plot gray dots for individual contrasts
        ax.scatter(
            pts, y_jit, color=FAINT, edgecolor="#8A847A", lw=0.35, s=30, alpha=0.72, zorder=2
        )
        # Overlay Tufte quartile rule
        draw_quartile_rule_h(
            ax,
            np.min(pts),
            q1,
            med,
            q3,
            np.max(pts),
            y=y0 - 0.25,
            color=col,
            lw_iqr=3.8,
            dot_size=55,
        )

        ax.text(
            -1.88,
            y0,
            cname,
            fontsize=12.2,
            fontweight="bold" if idx < 3 else "normal",
            color=INK,
            ha="right",
            va="center",
            clip_on=False,
        )
        ax.text(
            1.34,
            y0,
            f"Med {med:+.2f}  ({pct_down:.1f}% down)",
            fontsize=12.0,
            fontweight="bold",
            color=col,
            ha="left",
            va="center",
            clip_on=False,
        )

    ax.set_xlim(-1.85, 2.15)
    ax.set_ylim(-0.55, 6.55)
    ax.set_xticks([-1.5, -1.0, -0.5, 0.0, 0.5, 1.0])
    ax.set_xticklabels(
        [
            "\u22121.5 (Anaerobe loss)",
            "\u22121.0",
            "\u22120.5",
            "0.0 (Parity)",
            "+0.5",
            "+1.0 (Anaerobe gain)",
        ]
    )
    ax.set_xlabel(
        "Paired Contrast Log-Ratio ln(Obligate Anaerobe / Facultative)  \u00b7  Global 548-Study "
        "Sign: 60.7% Down (q = 9.6\u00d710\u207b\u2076)",
        fontsize=SLIDE_SIZES["axis_label"],
    )

    return save_and_verify_story_plot(
        fig, "plot_04_seven_conditions_logratio_strips", out_dir, random_marks=True
    )


# =============================================================================
# PLOT 05: 2D Ecological Redox Phase Space & Compartment Centroid Vectors
# Inspired by: tufte_058.jpg, tufte_042.jpg, tufte_094.png
# =============================================================================
def plot_05_redox_phase_space(out_dir: Path) -> StoryPlotInfo:
    fig = init_story_plot(
        5,
        "2D Quadrant Phase Space + Trajectory Vectors (cf. tufte_042, tufte_058)",
        "Coupled Oxygen Rewiring: Stool Contrasts Pull Into Aerobe-Up / Anaerobe-Down",
        "Across 918 active Tier 1+2 contrasts, 346 fall in Q4 (Aerobe UP & Anaerobe DOWN) vs. 236 "
        "in Q2 (OR = 1.47, p = 5.9\u00d710\u207b\u2076), driven by Feces (OR = 2.07).",
        "Source: 2\u00d72 McNemar discordance & compartment polarity across 918 active Tier 1+2 "
        "contrasts in 540 studies (Feces n = 592, Oral n = 142, GI mucosa n = 49).",
    )

    ax = add_panel(
        fig,
        (0.08, 0.15, 0.57, 0.56),
        "A. 918 Contrasts in (\u0394 Facultative Aerobe, \u0394 Obligate Anaerobe) Phase Space",
        "Gray dots = 918 contrasts; vectors & ellipses = compartment centroids",
    )
    rng = np.random.default_rng(505)
    # Generate 918 points matching McNemar quadrant counts (Q4=346, Q2=236, Q1+Q3=336)
    q4_x = np.abs(rng.normal(0.16, 0.14, size=346))
    q4_y = -np.abs(rng.normal(0.17, 0.14, size=346))
    q2_x = -np.abs(rng.normal(0.13, 0.12, size=236))
    q2_y = np.abs(rng.normal(0.13, 0.12, size=236))
    q13_x = rng.normal(0.0, 0.11, size=336)
    q13_y = np.sign(q13_x) * np.abs(rng.normal(0.09, 0.09, size=336))

    all_x = np.clip(np.concatenate([q4_x, q2_x, q13_x]), -0.48, 0.48)
    all_y = np.clip(np.concatenate([q4_y, q2_y, q13_y]), -0.48, 0.48)

    ax.axvline(0.0, color=INK, lw=1.1, ls="--", zorder=2)
    ax.axhline(0.0, color=INK, lw=1.1, ls="--", zorder=2)
    ax.scatter(all_x, all_y, color=FAINT, edgecolor="#8E887E", lw=0.3, s=26, alpha=0.58, zorder=1)

    # Highlight Q4 vs Q2 corners
    ax.text(
        0.46,
        -0.44,
        "Q4: Aerobe UP + Anaerobe DOWN\n346 contrasts (McNemar OR = 1.47, p = 5.9e-6)",
        fontsize=12.0,
        fontweight="bold",
        color=ACCENT,
        ha="right",
        va="bottom",
    )
    ax.text(
        -0.46,
        0.44,
        "Q2: Reverse (Aerobe DOWN + Anaerobe UP)\n236 contrasts",
        fontsize=11.8,
        fontweight="bold",
        color=SECONDARY,
        ha="left",
        va="top",
    )

    # Compartment trajectory vectors + 95% confidence ellipses
    compartments = [
        (
            "Feces / Stool (n = 592)\n263 vs 127 \u00b7 OR = 2.07, p = 5.0e-12",
            0.22,
            -0.24,
            0.26,
            0.22,
            -35,
            ACCENT,
            0.24,
            -0.10,
        ),
        (
            "Oral Cavity (n = 142)\n40 vs 54 \u00b7 Polarity \u22120.33",
            -0.07,
            0.06,
            0.16,
            0.14,
            25,
            SECONDARY,
            -0.24,
            0.20,
        ),
        (
            "Non-Fecal GI Mucosa (n = 49)\n11 vs 15 \u00b7 Polarity \u22120.36",
            -0.06,
            0.08,
            0.15,
            0.14,
            30,
            TERTIARY,
            0.16,
            0.34,
        ),
    ]
    for clab, cx, cy, ew, eh, ang, col, tx, ty in compartments:
        ell = Ellipse(
            (cx, cy),
            width=ew,
            height=eh,
            angle=ang,
            facecolor=col,
            alpha=0.14,
            edgecolor=col,
            lw=1.5,
            ls="-",
            zorder=3,
        )
        ax.add_patch(ell)
        ax.annotate(
            "",
            xy=(cx, cy),
            xytext=(0.0, 0.0),
            arrowprops={"arrowstyle": "-|>", "color": col, "lw": 2.6, "mutation_scale": 16},
            zorder=5,
        )
        ax.scatter([cx], [cy], color=col, s=90, edgecolor=PAPER, lw=1.2, zorder=6)
        ax.plot([cx, tx], [cy, ty], color=col, lw=1.1, zorder=4)
        ax.text(tx, ty, clab, fontsize=11.6, fontweight="bold", color=col, ha="center", va="center")

    ax.set_xlim(-0.50, 0.50)
    ax.set_ylim(-0.50, 0.50)
    ax.set_xticks([-0.4, -0.2, 0.0, 0.2, 0.4])
    ax.set_xticklabels(["\u22120.4", "\u22120.2", "0.0", "+0.2", "+0.4"])
    ax.set_yticks([-0.4, -0.2, 0.0, 0.2, 0.4])
    ax.set_yticklabels(["\u22120.4", "\u22120.2", "0.0", "+0.2", "+0.4"])
    ax.set_xlabel(
        "Within-Contrast Shift in Facultative Aerobe Fraction (\u0394z)",
        fontsize=SLIDE_SIZES["axis_label"],
    )
    ax.set_ylabel(
        "Within-Contrast Shift in Obligate Anaerobe Fraction (\u0394z)",
        fontsize=SLIDE_SIZES["axis_label"],
    )
    despine(ax, top=True, right=True)

    ax_b = add_panel(
        fig,
        (0.70, 0.15, 0.26, 0.56),
        "B. Compartment Rewiring",
        "Opposite baseline O\u2082 gradients",
    )
    ax_b.set_xlim(0, 1)
    ax_b.set_ylim(0, 1)
    ax_b.axis("off")
    draw_story_callout(
        ax_b,
        0.0,
        0.70,
        1.0,
        0.25,
        "OR = 2.07",
        "Stool Coupled Shift",
        "263 vs 127 (p = 5.0e-12)",
        ACCENT,
    )
    draw_story_callout(
        ax_b,
        0.0,
        0.37,
        1.0,
        0.25,
        "+0.76",
        "Stool Redox Polarity",
        "Anoxic lumen gains aerobes",
        SECONDARY,
    )
    draw_story_callout(
        ax_b,
        0.0,
        0.04,
        1.0,
        0.25,
        "\u22120.36 / \u22120.33",
        "Mucosa & Oral Polarity",
        "No aerobe gain at O\u2082-exposed sites",
        TERTIARY,
    )

    return save_and_verify_story_plot(
        fig, "plot_05_redox_phase_space_quadrants", out_dir, random_marks=True
    )


# =============================================================================
# PLOT 06: Organism-Level Morphology x Oxygen Co-Annotation Cloud Across 842 Genera
# Inspired by: tufte_006.jpg, tufte_095.png
# =============================================================================
def plot_06_organism_guild_cloud(out_dir: Path) -> StoryPlotInfo:
    fig = init_story_plot(
        6,
        "2D Guild Propensity Cloud + Richness Bubbles (cf. tufte_006, tufte_095)",
        "Why Cell Shape Tracks Oxygen Tolerance: Two Opposing Multi-Trait Bacterial Guilds",
        "Plotting all 842 Tier 1+2 genera shows that Coccus + Facultative pathobionts enrich in "
        "disease while Rod + Anaerobic fiber fermenters deplete.",
        "Source: Genus-level directional propensity (Inc \u2212 Dec)/(Inc + Dec) vs. total "
        "occurrences across 842 Tier 1+2 genera (935 contrasts).",
    )

    ax = add_panel(
        fig,
        (0.08, 0.15, 0.86, 0.56),
        "842 Tier 1+2 Genera: Net Disease Shift Propensity vs. Contrast Occurrences (Log10 Scale)",
        "Gray dots = 831 background genera; Slate = Rod + Obligate Anaerobe + Cellulose guild; "
        "Crimson = Coccus + Facultative + Halophilic guild",
    )
    rng = np.random.default_rng(606)
    n_bg = 831
    bg_occ = np.clip(rng.exponential(0.45, size=n_bg), 0.0, 1.85)
    # Funnel narrowing at higher occurrence
    bg_prop = rng.normal(0.0, 0.42 / (1.0 + 0.8 * bg_occ), size=n_bg)
    bg_prop = np.clip(bg_prop, -0.88, 0.88)
    bg_sizes = 18 + 35 * (bg_occ / 2.0)

    ax.axvline(0.0, color=INK, lw=1.1, ls="--", zorder=2)
    ax.scatter(
        bg_prop, bg_occ, s=bg_sizes, color=FAINT, edgecolor="#8A847A", lw=0.35, alpha=0.62, zorder=1
    )

    # Left guild: Rod / Bacillus + Obligate Anaerobe + Cellulose degraders (Depleted in disease)
    rod_anaerobes = [
        (-0.18, 2.46, "Bacteroides (288 occ \u00b7 Rod + Anaerobe + Cell.)", TERTIARY, -0.24, 2.50),
        (
            -0.46,
            2.36,
            "Faecalibacterium (229 occ \u00b7 Rod + Anaerobe, 0 Cell.)",
            SECONDARY,
            -0.50,
            2.36,
        ),
        (-0.34, 2.26, "Ruminococcus (182 occ \u00b7 Anaerobe + Cell.)", SECONDARY, -0.50, 2.10),
        (-0.52, 2.18, "Roseburia (151 occ \u00b7 Rod + Anaerobe + Cell.)", SECONDARY, -0.54, 1.84),
        (-0.41, 1.90, "Eubacterium (79 occ \u00b7 Rod + Anaerobe + Cell.)", SECONDARY, -0.52, 1.56),
        (-0.44, 1.82, "Lachnospira (66 occ \u00b7 Rod + Anaerobe + Cell.)", SECONDARY, -0.52, 1.28),
    ]
    for px, py, plab, col, tx, ty in rod_anaerobes:
        ax.scatter([px], [py], s=155, color=col, edgecolor=PAPER, lw=1.3, zorder=5)
        ax.plot([px, tx + 0.01], [py, ty], color=col, lw=1.2, zorder=4)
        ax.text(tx, ty, plab, fontsize=11.4, fontweight="bold", color=col, ha="right", va="center")

    # Right guild: Coccus + Facultative Anaerobe + Halophilic (Enriched in disease)
    coccus_fac = [
        (0.42, 2.31, "Streptococcus (204 occ \u00b7 Coccus + Facultative)", ACCENT, 0.48, 2.36),
        (
            0.39,
            2.08,
            "Escherichia (120 occ \u00b7 Facultative + Arsenate red.)",
            HIGHLIGHT,
            0.52,
            2.10,
        ),
        (0.48, 2.02, "Enterococcus (105 occ \u00b7 Coccus + Fac. + Halo.)", ACCENT, 0.54, 1.84),
        (0.51, 1.65, "Staphylococcus (45 occ \u00b7 Coccus + Fac. + Halo.)", ACCENT, 0.56, 1.58),
        (0.36, 1.54, "Gemella (35 occ \u00b7 Coccus + Facultative)", HIGHLIGHT, 0.46, 1.32),
    ]
    for px, py, plab, col, tx, ty in coccus_fac:
        ax.scatter([px], [py], s=155, color=col, edgecolor=PAPER, lw=1.3, zorder=5)
        ax.plot([px, tx - 0.01], [py, ty], color=col, lw=1.2, zorder=4)
        ax.text(tx, ty, plab, fontsize=11.4, fontweight="bold", color=col, ha="left", va="center")

    ax.set_xlim(-1.18, 1.06)
    ax.set_ylim(-0.08, 2.68)
    ax.set_xticks([-0.6, -0.3, 0.0, 0.3, 0.6])
    ax.set_xticklabels(
        [
            "\u22120.60 (Depleted in Disease)",
            "\u22120.30",
            "0.00 (Balanced)",
            "+0.30",
            "+0.60 (Enriched in Disease)",
        ]
    )
    ax.set_yticks([0.0, 0.699, 1.301, 2.0, 2.477])
    ax.set_yticklabels(["1 occ", "5 occ", "20 occ", "100 occ", "300 occ"])
    ax.set_xlabel(
        "Genus Net Disease Shift Propensity (Increased \u2212 Decreased) / Total Occurrences",
        fontsize=SLIDE_SIZES["axis_label"],
    )
    ax.set_ylabel("Contrast Occurrences (Log10)", fontsize=SLIDE_SIZES["axis_label"])
    despine(ax, top=True, right=True)

    return save_and_verify_story_plot(
        fig, "plot_06_organism_guild_mechanism_cloud", out_dir, random_marks=True
    )


# =============================================================================
# PLOT 07: Unified 201-Trait Gray-Dot Strip + Five 10,000-Permutation Null Barcode Strips
# Inspired by: slide_17_90_pathway_benchmark.png + tufte_070.png + tufte_035.png
# =============================================================================
def plot_07_pathway_universe_and_nulls(out_dir: Path) -> StoryPlotInfo:
    fig = init_story_plot(
        7,
        "Unified Gray-Dot Empirical Strip + Permutation Null Barcodes (cf. slide_17, tufte_070)",
        "Cellulose Catabolism (44.4%) Beats All Five Permutation Nulls and All 200 Other Traits",
        "Placing 50,000 permutation draws (U\u2081\u2013U\u2085) and all 201 qualifying real "
        "pathways on one axis proves both statistical and biological extremity.",
        "Source: 10,000-permutation null distributions (U1\u2013U5) and 201 qualifying real "
        "traits across Tier 1+2 PD (71), CRC (44), and HIV (20) stool contrasts.",
    )

    ax = add_panel(
        fig,
        (0.28, 0.15, 0.66, 0.56),
        "Mean Pairwise Cross-Disease Jaccard Overlap Across PD, CRC, and HIV Stool (%)",
        "Top 5 rows = 10,000-draw permutation nulls (U\u2081\u2013U\u2085); Bottom row = 201 real "
        "traits (gray dots + quartile rule)",
    )
    despine(ax, top=True, right=True, left=True, bottom=False)
    ax.yaxis.set_visible(False)

    # Vertical benchmark line at Observed Cellulose = 44.4%
    ax.axvline(44.4, color=ACCENT, lw=1.8, ls="--", zorder=3)

    rng = np.random.default_rng(707)
    null_specs = [
        ("U\u2081: Graph-Wide Taxa Null", 8.2, 2.1, 4.2, 13.5, "Mean 8.2%  (Obs = 5.4\u00d7)"),
        (
            "U\u2082: Tier 1+2 Human Taxa Null (2,251)",
            11.4,
            2.4,
            6.5,
            17.2,
            "Mean 11.4%  (Obs = 3.9\u00d7)",
        ),
        (
            "U\u2083: Stool Cellulose Urn Null",
            14.8,
            2.8,
            8.8,
            21.8,
            "Mean 14.8%  (Obs = 3.0\u00d7)",
        ),
        (
            "U\u2085: Genus+Species Depleted Null (23)",
            20.0,
            3.4,
            12.8,
            28.4,
            "Mean 20.0%  (G+S Obs 42.9% = 2.1\u00d7, p = 7.0e-4)",
        ),
        (
            "U\u2084: All Depleted Stool Taxa Null (402)",
            21.8,
            3.2,
            14.0,
            29.8,
            "Mean 21.8%  (Obs 44.4% = 2.0\u00d7, p = 3.4e-6)",
        ),
    ]

    for idx, (nname, mu, sd, lo, hi, rtxt) in enumerate(null_specs):
        y0 = 5.4 - idx * 0.85
        draws = np.clip(rng.normal(mu, sd, size=180), lo, hi)
        # Barcode ticks + gray cloud
        ax.scatter(
            draws, np.full_like(draws, y0), marker="|", color=FAINT, s=65, alpha=0.55, zorder=2
        )
        draw_quartile_rule_h(
            ax,
            lo,
            mu - 0.67 * sd,
            mu,
            mu + 0.67 * sd,
            hi,
            y=y0 - 0.18,
            color=SECONDARY if idx >= 3 else MUTED,
            lw_iqr=3.2,
            dot_size=42,
        )
        ax.text(
            -0.8,
            y0,
            nname,
            fontsize=11.8,
            fontweight="bold" if idx >= 3 else "normal",
            color=INK,
            ha="right",
            va="center",
            clip_on=False,
        )
        ax.text(
            hi + 0.8,
            y0,
            rtxt,
            fontsize=11.4,
            fontweight="bold" if idx >= 3 else "normal",
            color=SECONDARY if idx >= 3 else MUTED,
            ha="left",
            va="center",
        )

    # Separator rule above empirical 201-trait strip
    ax.axhline(1.38, color=RULE, lw=0.9, ls=":")

    # Bottom strip: The 201 real biological pathways (Gray-dot empirical universe)
    bg_x = np.concatenate(
        [
            rng.uniform(7.2, 27.4, 50),
            rng.uniform(27.4, 32.0, 50),
            rng.uniform(32.0, 36.2, 50),
            rng.uniform(36.2, 41.0, 48),
        ]
    )
    bg_y = rng.uniform(0.35, 0.92, size=len(bg_x))
    ax.scatter(bg_x, bg_y, color=FAINT, edgecolor="#8A847A", lw=0.35, s=42, alpha=0.80, zorder=2)
    draw_quartile_rule_h(
        ax, 7.2, 27.4, 32.0, 36.2, 44.4, y=0.08, color=INK, lw_iqr=4.2, dot_size=65
    )
    ax.text(
        -0.8,
        0.62,
        "Empirical Universe:\n201 Qualifying Traits\n(Median 32.0%, IQR 27.4\u201336.2%)",
        fontsize=12.0,
        fontweight="bold",
        color=INK,
        ha="right",
        va="center",
        clip_on=False,
    )

    # Top 3 GO BP callouts (well separated horizontally so no overlap)
    top3 = [
        (44.4, 0.78, "#1 Cellulose catabolism (44.4% \u00b7 44 taxa)", ACCENT, 36.5, 1.18),
        (43.5, 0.58, "#2 Nitrogen fixation (43.5% \u00b7 33 taxa)", SECONDARY, 16.5, 1.18),
        (41.3, 0.40, "#3 Xylan catabolism (41.3% \u00b7 55 taxa)", TERTIARY, 38.2, -0.22),
    ]
    for xv, yv, tlab, col, tx, ty in top3:
        ax.scatter([xv], [yv], color=col, s=115, edgecolor=PAPER, lw=1.2, zorder=6)
        ax.plot([xv, tx], [yv, ty], color=col, lw=1.2, zorder=5)
        ax.text(tx, ty, tlab, fontsize=11.4, fontweight="bold", color=col, ha="center", va="center")

    ax.set_xlim(0.0, 48.5)
    ax.set_ylim(-0.45, 5.90)
    ax.set_xticks([0, 10, 21.8, 32.0, 44.4])
    ax.set_xticklabels(
        ["0%", "10%", "21.8% (U\u2084 Null)", "32.0% (201-Trait Med.)", "44.4% (Cellulose #1)"]
    )
    ax.set_xlabel(
        "Mean Pairwise Jaccard Overlap Across PD, CRC & HIV Stool Contrasts (%)",
        fontsize=SLIDE_SIZES["axis_label"],
    )

    return save_and_verify_story_plot(
        fig, "plot_07_pathway_universe_and_five_nulls", out_dir, random_marks=True
    )


# =============================================================================
# PLOT 08: Barycentric 3-Disease Guild Affinity Simplex of All 402 Depleted Stool Taxa
# Inspired by: tufte_094.png, tufte_012.png, tufte_095.png
# =============================================================================
def plot_08_barycentric_disease_simplex(out_dir: Path) -> StoryPlotInfo:
    fig = init_story_plot(
        8,
        "Barycentric 3-Disease Simplex Cloud (cf. tufte_094, tufte_095)",
        "How the 44 Cellulose-Degrading Taxa Pull Into the 3-Disease Convergence Zone",
        "Plotting all 402 depleted stool taxa in a PD\u2013CRC\u2013HIV simplex shows 358 "
        "background taxa at the corners vs. 44 guild taxa converging.",
        "Source: 402 unique depleted stool taxa across Tier 1+2 PD (n = 71), CRC (n = 44), and "
        "HIV (n = 20) contrasts (44 cellulose guild taxa: 12 in all 3, 15 in 2, 17 in 1).",
    )

    ax = add_panel(
        fig,
        (0.06, 0.14, 0.58, 0.57),
        "A. Barycentric Disease-Sharing Simplex (N = 402 Depleted Stool Taxa)",
        "Gray = 358 background taxa; Crimson = 12 all-3 core; Slate = 15 two-disease; Ochre = 17 "
        "one-disease",
    )
    ax.axis("off")

    # Equilateral triangle vertices: Top = PD (0.5, 0.88), Bottom-Left = CRC (0.12, 0.14),
    # Bottom-Right = HIV (0.88, 0.14)
    v_pd = np.array([0.50, 0.86])
    v_crc = np.array([0.14, 0.14])
    v_hiv = np.array([0.86, 0.14])

    tri_x = [v_pd[0], v_crc[0], v_hiv[0], v_pd[0]]
    tri_y = [v_pd[1], v_crc[1], v_hiv[1], v_pd[1]]
    ax.plot(tri_x, tri_y, color=INK, lw=1.4, zorder=2)

    # Internal dashed barycentric medians meeting at centroid
    cent = (v_pd + v_crc + v_hiv) / 3.0
    for v in (v_pd, v_crc, v_hiv):
        ax.plot([v[0], cent[0]], [v[1], cent[1]], color=RULE, lw=0.9, ls="--", zorder=1)

    rng = np.random.default_rng(808)

    def bary_to_xy(w_pd: float, w_crc: float, w_hiv: float) -> np.ndarray:
        s = w_pd + w_crc + w_hiv
        return (w_pd * v_pd + w_crc * v_crc + w_hiv * v_hiv) / s

    # 358 background non-cellulose taxa (mostly corners & edges)
    bg_pts = []
    for _ in range(228):  # single-disease corners
        corner = rng.choice([0, 1, 2], p=[0.46, 0.34, 0.20])
        w = np.array([0.04, 0.04, 0.04])
        w[corner] = rng.uniform(0.75, 0.95)
        w += rng.uniform(0.0, 0.08, size=3)
        bg_pts.append(bary_to_xy(*w))
    for _ in range(95):  # 2-disease edges
        edge = rng.choice([0, 1, 2])
        w = np.array([0.45, 0.45, 0.45])
        w[edge] = rng.uniform(0.02, 0.08)
        w += rng.normal(0, 0.06, size=3)
        w = np.clip(w, 0.02, 1.0)
        bg_pts.append(bary_to_xy(*w))
    for _ in range(35):  # 3-disease center background
        w = rng.uniform(0.25, 0.45, size=3)
        bg_pts.append(bary_to_xy(*w))
    bg_arr = np.array(bg_pts)
    ax.scatter(
        bg_arr[:, 0],
        bg_arr[:, 1],
        color=FAINT,
        edgecolor="#8E887E",
        lw=0.3,
        s=26,
        alpha=0.62,
        zorder=2,
    )

    # 17 single-disease cellulose taxa (9 PD, 6 CRC, 2 HIV) in Ochre
    single_pts = []
    for c_idx, cnt in [(0, 9), (1, 6), (2, 2)]:
        for _ in range(cnt):
            w = np.array([0.06, 0.06, 0.06])
            w[c_idx] = rng.uniform(0.72, 0.88)
            w += rng.uniform(0.0, 0.06, size=3)
            single_pts.append(bary_to_xy(*w))
    s_arr = np.array(single_pts)
    ax.scatter(s_arr[:, 0], s_arr[:, 1], color=HIGHLIGHT, edgecolor=PAPER, lw=0.9, s=68, zorder=4)

    # 15 two-disease cellulose taxa (6 PD+CRC, 6 PD+HIV, 3 HIV+CRC) in Slate
    two_pts = []
    for missing_idx, cnt in [(2, 6), (1, 6), (0, 3)]:
        for _ in range(cnt):
            w = np.array([0.46, 0.46, 0.46])
            w[missing_idx] = rng.uniform(0.04, 0.08)
            w += rng.uniform(-0.06, 0.06, size=3)
            two_pts.append(bary_to_xy(*w))
    t_arr = np.array(two_pts)
    ax.scatter(t_arr[:, 0], t_arr[:, 1], color=SECONDARY, edgecolor=PAPER, lw=1.0, s=82, zorder=5)

    # 12 all-three core cellulose taxa in Crimson near center
    core_pts = []
    for _ in range(12):
        w = rng.uniform(0.28, 0.42, size=3)
        core_pts.append(bary_to_xy(*w))
    c_arr = np.array(core_pts)
    ax.scatter(c_arr[:, 0], c_arr[:, 1], color=ACCENT, edgecolor=PAPER, lw=1.2, s=110, zorder=6)

    # Vertex labels
    ax.text(
        0.50,
        0.93,
        "Parkinson's Disease Stool (n = 71)\n33 depleted cellulose taxa",
        fontsize=12.2,
        fontweight="bold",
        color=INK,
        ha="center",
        va="bottom",
    )
    ax.text(
        0.02,
        0.06,
        "Colorectal Cancer Stool (n = 44)\n27 depleted cellulose taxa",
        fontsize=12.0,
        fontweight="bold",
        color=INK,
        ha="left",
        va="top",
    )
    ax.text(
        0.98,
        0.06,
        "HIV Infection Stool (n = 20)\n23 depleted cellulose taxa",
        fontsize=12.0,
        fontweight="bold",
        color=INK,
        ha="right",
        va="top",
    )

    # Callout to central 12-taxon core
    ax.plot([cent[0], 0.82], [cent[1], 0.68], color=ACCENT, lw=1.4, zorder=5)
    ax.text(
        0.83,
        0.68,
        "12 All-Three Core Taxa (27.3%)\nRoseburia, Ruminococcus,\nBacteroides, "
        "Clostridium,\nLachnospira, Anaerostipes",
        fontsize=11.8,
        fontweight="bold",
        color=ACCENT,
        ha="left",
        va="center",
    )
    ax.plot([0.29, 0.03], [0.52, 0.66], color=SECONDARY, lw=1.3, zorder=5)
    ax.text(
        0.02,
        0.66,
        "6 PD \u2229 CRC Only\n(Full PD \u2229 CRC = 18/42, 42.9%)",
        fontsize=11.6,
        fontweight="bold",
        color=SECONDARY,
        ha="left",
        va="center",
    )

    ax.set_xlim(0.0, 1.15)
    ax.set_ylim(0.0, 1.02)

    ax_b = add_panel(
        fig,
        (0.71, 0.15, 0.25, 0.56),
        "B. Guild Convergence",
        "44 cellulose taxa vs. 402 U\u2084",
    )
    ax_b.set_xlim(0, 1)
    ax_b.set_ylim(0, 1)
    ax_b.axis("off")
    draw_story_callout(
        ax_b,
        0.0,
        0.70,
        1.0,
        0.25,
        "61.4%",
        "Shared in \u22652 Diseases",
        "27 / 44 cellulose taxa (12 in all 3)",
        ACCENT,
    )
    draw_story_callout(
        ax_b,
        0.0,
        0.37,
        1.0,
        0.25,
        "44.4%",
        "Mean Pairwise Jaccard",
        "2.0\u00d7 U\u2084 baseline (21.8%, p = 3.4e-6)",
        SECONDARY,
    )
    draw_story_callout(
        ax_b,
        0.0,
        0.04,
        1.0,
        0.25,
        "42.9%",
        "Genus + Species Only",
        "23 taxa (2.1\u00d7 vs 20.0%, p = 7.0e-4)",
        TERTIARY,
    )

    return save_and_verify_story_plot(
        fig, "plot_08_barycentric_disease_simplex", out_dir, random_marks=True
    )


# =============================================================================
# PLOT 09: Anatomical Compartment Contrast Barcode & Quartile-Gap Strip Plot
# Inspired by: tufte_070.png, tufte_059.jpg, tufte_071.png
# =============================================================================
def plot_09_compartment_barcode_strips(out_dir: Path) -> StoryPlotInfo:
    fig = init_story_plot(
        9,
        "Anatomical Barcode Tick-Strips + Quartile Rules (cf. tufte_070, tufte_071)",
        "Cellulose Degrader Loss Is Strictly a Luminal Stool Signal, Absent in Mucosa & Oral Sites",
        "Plotting every contrast as a barcode tick shows strong depletion in stool "
        "(63.5%\u201375.4% down) vs. null/reversed shifts in mucosa and oral sites.",
        "Source: Compartment-stratified within-contrast cellulose shifts (\u0394z) across 935 "
        "Tier 1+2 contrasts (Stool vs. Oral Fisher OR = 5.10, p = 5.3\u00d710\u207b\u2076).",
    )

    ax = add_panel(
        fig,
        (0.28, 0.15, 0.66, 0.56),
        "Individual Contrast Cellulose Shifts (\u0394z) Across Stool, Mucosal, and Oral Sites",
        "Each vertical tick = 1 contrast; horizontal rule = Q1 \u2013 [hollow Median] \u2013 Q3; "
        "dashed line = 0.0 parity",
    )
    despine(ax, top=True, right=True, left=True, bottom=False)
    ax.yaxis.set_visible(False)
    ax.axvline(0.0, color=INK, lw=1.2, ls="--", zorder=2)

    strata = [
        (
            "Stool: Parkinson's (n = 61 non-tied)",
            61,
            -0.71,
            -1.15,
            -0.18,
            "75.4% Dec (46:15 \u00b7 p = 9.0e-6)",
            ACCENT,
        ),
        (
            "Stool: HIV Infection (n = 17 non-tied)",
            17,
            -0.47,
            -0.92,
            0.04,
            "70.6% Dec (12:5 \u00b7 p = 0.044)",
            ACCENT,
        ),
        (
            "Stool: Colorectal Cancer (n = 33 non-tied)",
            33,
            -0.43,
            -0.85,
            0.08,
            "69.7% Dec (23:10 \u00b7 p = 0.036)",
            ACCENT,
        ),
        (
            "All Tier 1+2 Stool / Feces (n = 496)",
            180,
            -0.25,
            -0.68,
            0.12,
            "63.5% Dec (315:181 \u00b7 p = 1.9e-9)",
            SECONDARY,
        ),
        (
            "Non-Fecal GI Mucosa Biopsies (n = 32)",
            32,
            -0.05,
            -0.42,
            0.34,
            "53.1% Dec (17:15 \u00b7 p = 0.86 \u00b7 Null)",
            MUTED,
        ),
        (
            "Oral Cavity Pooled (n = 52)",
            52,
            0.22,
            -0.18,
            0.64,
            "34.6% Dec (18:34 \u00b7 OR = 5.10 vs Stool)",
            HIGHLIGHT,
        ),
    ]

    rng = np.random.default_rng(909)
    for idx, (sname, n_pts, med, q1, q3, rtxt, col) in enumerate(strata):
        y0 = float(len(strata) - 1 - idx)
        pts = np.clip(rng.normal(med, 0.52, size=n_pts), -1.85, 1.35)
        # Barcode ticks + subtle gray dots
        ax.scatter(
            pts,
            np.full_like(pts, y0 + 0.06),
            marker="|",
            color=FAINT if col == MUTED else col,
            s=95,
            alpha=0.48,
            zorder=2,
        )
        ax.scatter(
            pts, y0 + rng.uniform(-0.10, 0.18, size=n_pts), color=FAINT, s=18, alpha=0.55, zorder=1
        )
        draw_quartile_rule_h(
            ax,
            np.min(pts),
            q1,
            med,
            q3,
            np.max(pts),
            y=y0 - 0.22,
            color=col,
            lw_iqr=3.8,
            dot_size=52,
        )

        ax.text(
            -1.95,
            y0,
            sname,
            fontsize=12.0,
            fontweight="bold" if idx < 4 else "normal",
            color=INK,
            ha="right",
            va="center",
            clip_on=False,
        )
        ax.text(
            1.42,
            y0,
            rtxt,
            fontsize=11.8,
            fontweight="bold" if col != MUTED else "normal",
            color=col,
            ha="left",
            va="center",
            clip_on=False,
        )

    ax.set_xlim(-1.92, 2.50)
    ax.set_ylim(-0.55, 5.55)
    ax.set_xticks([-1.5, -1.0, -0.5, 0.0, 0.5, 1.0])
    ax.set_xticklabels(
        [
            "\u22121.5 (Depleted)",
            "\u22121.0",
            "\u22120.5",
            "0.0 (Parity)",
            "+0.5",
            "+1.0 (Enriched)",
        ]
    )
    ax.set_xlabel(
        "Within-Contrast Cellulose Catabolism Shift (\u0394z)  \u00b7  Pooled PD+CRC+HIV Stool: "
        "73.0% Dec (81:30, p = 1.4\u00d710\u207b\u2076)",
        fontsize=SLIDE_SIZES["axis_label"],
    )

    return save_and_verify_story_plot(
        fig, "plot_09_compartment_barcode_strip_forest", out_dir, random_marks=True
    )


# =============================================================================
# PLOT 10: 2,251-Taxon Leave-One-Taxon-Out (LOTO) Jackknife Cloud & Platform Strips
# Inspired by: tufte_016.jpg, tufte_063.png, tufte_130.webp
# =============================================================================
def plot_10_loto_and_platform_strips(out_dir: Path) -> StoryPlotInfo:
    fig = init_story_plot(
        10,
        "2,251-Taxon LOTO Jackknife Cloud + Assay Strips (cf. tufte_016, tufte_130)",
        "Zero Sign Flips Across 2,251 Taxon Drops and Full Concordance Across 16S & WMS",
        "Dropping any single taxon out of 2,251 never flips a global trait direction; both 16S (n "
        "= 807) and Shotgun WMS (n = 128) confirm the shift.",
        "Source: Leave-One-Taxon-Out (LOTO) sensitivity across 2,251 Tier 1+2 taxa and sequencing "
        "platform stratification across 935 contrasts.",
    )

    ax_a = add_panel(
        fig,
        (0.21, 0.15, 0.39, 0.56),
        "A. 2,251 Leave-One-Taxon-Out Jackknife Cloud",
        "Gray dots = 2,251 recalculated shifts when each taxon is dropped; 0 sign flips",
    )
    despine(ax_a, top=True, right=True, left=True, bottom=False)
    ax_a.yaxis.set_visible(False)
    ax_a.axvline(50.0, color=INK, lw=1.2, ls="--", zorder=2)

    loto_traits = [
        ("Cellulose catabolism", 40.1, 41.9, "40.1% (max 41.9%, Ruminococcus)", SECONDARY),
        ("Obligate anaerobe", 48.0, 49.1, "48.0% (max 49.1%, Bacteroides)", HIGHLIGHT),
        ("Arsenate reduction", 57.1, 55.6, "57.1% (min 55.6%, Escherichia)", ACCENT),
        ("Halophilic", 57.1, 55.5, "57.1% (min 55.5%, Staphylococcus)", ACCENT),
        ("Facultative aerobe", 55.7, 54.2, "55.7% (min 54.2%, Streptococcus)", ACCENT),
        ("Coccus shaped", 55.2, 53.5, "55.2% (min 53.5%, Streptococcus)", TERTIARY),
    ]

    rng = np.random.default_rng(1010)
    for idx, (tname, full_pct, max_drop_pct, drop_lab, col) in enumerate(loto_traits):
        y0 = float(len(loto_traits) - 1 - idx)
        # 2,250 taxa barely move the percentage; a few high-occurrence genera shift it by 0.5-1.9 pp
        jitter_x = full_pct + rng.normal(0.0, 0.22, size=220)
        tail_x = np.linspace(full_pct, max_drop_pct, 18)
        all_loto = np.concatenate([jitter_x, tail_x])
        y_jit = y0 + rng.uniform(-0.16, 0.16, size=len(all_loto))

        ax_a.scatter(
            all_loto, y_jit, color=FAINT, edgecolor="#8A847A", lw=0.3, s=26, alpha=0.70, zorder=2
        )
        ax_a.scatter([full_pct], [y0], color=col, s=95, edgecolor=PAPER, lw=1.2, zorder=5)
        ax_a.scatter([max_drop_pct], [y0], facecolor=PAPER, edgecolor=col, s=65, lw=1.8, zorder=5)

        ax_a.text(
            35.2,
            y0,
            tname,
            fontsize=11.8,
            fontweight="bold",
            color=INK,
            ha="right",
            va="center",
            clip_on=False,
        )
        if full_pct < 50.0:
            ax_a.text(
                max_drop_pct + 0.6, y0, drop_lab, fontsize=10.8, color=col, ha="left", va="center"
            )
        else:
            ax_a.text(
                full_pct + 0.6, y0, drop_lab, fontsize=10.8, color=col, ha="left", va="center"
            )

    ax_a.set_xlim(35.8, 68.5)
    ax_a.set_ylim(-0.55, 5.55)
    ax_a.set_xticks([40, 45, 50, 55, 60])
    ax_a.set_xticklabels(["40% Inc", "45%", "50% (Null)", "55%", "60% Inc"])
    ax_a.set_xlabel(
        "Recalculated Contrast % Increased Across 2,251 Taxon Drops",
        fontsize=SLIDE_SIZES["axis_label"],
    )

    # Panel B: 16S (n=807) vs Shotgun WMS (n=128) Contrast Strips
    ax_b = add_panel(
        fig,
        (0.73, 0.15, 0.21, 0.56),
        "B. 16S vs. Shotgun WMS",
        "Assay-stratified contrast clouds",
    )
    despine(ax_b, top=True, right=True, left=False, bottom=False)
    ax_b.axvline(0.0, color=INK, lw=1.1, ls="--", zorder=2)

    plat_rows = [
        ("Cellulose\n16S (807)", -0.18, -0.52, 0.14, "62.8% Dec", SECONDARY),
        ("Cellulose\nWMS (128)", -0.25, -0.64, 0.08, "69.0% Dec", ACCENT),
        ("Facultative\n16S (807)", 0.16, -0.14, 0.48, "55.5% Inc", SECONDARY),
        ("Facultative\nWMS (128)", 0.21, -0.10, 0.56, "56.6% Inc", ACCENT),
    ]
    for idx, (_plab, med, q1, q3, ptxt, col) in enumerate(plat_rows):
        y0 = float(len(plat_rows) - 1 - idx)
        pts = np.clip(rng.normal(med, 0.42, size=85), -1.1, 1.1)
        ax_b.scatter(
            pts, y0 + rng.uniform(-0.15, 0.15, size=85), color=FAINT, s=20, alpha=0.65, zorder=1
        )
        draw_quartile_rule_h(
            ax_b, -0.95, q1, med, q3, 0.95, y=y0 - 0.22, color=col, lw_iqr=3.5, dot_size=48
        )
        ax_b.text(
            1.05,
            y0,
            ptxt,
            fontsize=11.5,
            fontweight="bold",
            color=col,
            ha="left",
            va="center",
            clip_on=False,
        )

    ax_b.set_yticks([3, 2, 1, 0])
    ax_b.set_yticklabels([r[0] for r in plat_rows], fontsize=11.2, color=INK)
    ax_b.set_xlim(-1.15, 1.85)
    ax_b.set_ylim(-0.55, 3.55)
    ax_b.set_xticks([-0.8, 0.0, 0.8])
    ax_b.set_xticklabels(["\u22120.8", "0.0", "+0.8"])
    ax_b.set_xlabel("Contrast Shift \u0394z", fontsize=SLIDE_SIZES["axis_label"])

    return save_and_verify_story_plot(
        fig, "plot_10_loto_jackknife_and_platform_strips", out_dir, random_marks=True
    )


PLOT_BUILDERS: list[Callable[[Path], StoryPlotInfo]] = [
    plot_01_joint_effect_consistency,
    plot_02_paired_parity_variance,
    plot_03_rank_occurrence_landscape,
    plot_04_seven_conditions_logratio_strips,
    plot_05_redox_phase_space,
    plot_06_organism_guild_cloud,
    plot_07_pathway_universe_and_nulls,
    plot_08_barycentric_disease_simplex,
    plot_09_compartment_barcode_strips,
    plot_10_loto_and_platform_strips,
]

# The plots whose points come from a random generator (#775). The build draws these only
# with --audit-random-marks; remove a plot here when it reads run tables instead.
RANDOM_MARK_PLOTS: frozenset[Callable[[Path], StoryPlotInfo]] = frozenset(PLOT_BUILDERS)

PLOT_SPEAKER_NOTES: dict[int, str] = {
    1: (
        "Data Story Plot 01 (Joint Effect-Size x Directional Consistency Cloud): Plotting all 317 "
        "non-MediaDive biological traits "
        "in the joint plane of paired effect size (Delta z) vs. contrast sign fraction (% "
        "Increased across 918 active Tier 1+2 contrasts in 540 studies) "
        "with marginal step-outline histograms and dot-dash rugs shows why 40 traits survive "
        "dual-FDR criteria (Arm 1 BY <= 0.05 & Arm 3 BH <= 0.05), led by Cellulose catabolism "
        "(Delta z = -0.25, 40.1% Inc / 59.9% Dec, q = 1.5e-7) and Arsenate reduction (Delta z = "
        "+0.25, 57.1% Inc, q = 1.3e-4)."
    ),
    2: (
        "Data Story Plot 02 (Paired Parity Cloud & Variance Stabilization): Plotting all 317 "
        "traits on the parity plane of mean "
        "trait fraction in Increased vs. Decreased sets across 935 Tier 1+2 contrasts (548 "
        "studies) reveals tight pairing along y = x that "
        "absorbs between-study baseline variance (9.5x variance collapse from unpaired SD 0.210 "
        "to paired SD 0.022; curation R^2 drops 51x to 0.18%)."
    ),
    3: (
        "Data Story Plot 03 (Taxonomic Rank x Occurrence Coverage Landscape): Plotting all 3,221 "
        "Tier 1+2 TaxIDs across 6 taxonomic "
        "ranks against log10 signature occurrences shows that high-frequency taxa (>= 10 "
        "contrasts) achieve 91.3% occurrence-weighted KG-Microbe coverage (97.4% at Genus), "
        "whereas the 789 unmapped species are overwhelmingly singletons."
    ),
    4: (
        "Data Story Plot 04 (Seven Major Conditions Log-Ratio Strip-Matrix): Across all 7 major "
        "conditions (387 contrasts), "
        "within-study log-ratio distributions ln(Anaerobe/Facultative) and "
        "ln(Cellulose/Aerotolerant) shift left of zero in 6 of 7 "
        "conditions, led by IBD (-0.62), Colorectal Cancer (-0.54), and Parkinson's Disease "
        "(-0.48)."
    ),
    5: (
        "Data Story Plot 05 (Coupled 2D Redox Phase Space): Plotting all 918 active Tier 1+2 "
        "contrasts in the 2D plane of (Delta Facultative "
        "Aerobe, Delta Obligate Anaerobe) shows 346 contrasts in the coupled Q4 quadrant (Aerobe "
        "Up + Anaerobe Down) vs. 236 in Q2 "
        "(McNemar OR = 1.47, p = 5.9e-6; Stool OR = 2.07, p = 5.0e-12)."
    ),
    6: (
        "Data Story Plot 06 (Organism-Level Guild Mechanism Cloud): Plotting all 842 Tier 1+2 "
        "genera by their KG-Microbe trait "
        "annotations and empirical contrast win-rate shows that the 44 cellulose-degrading guild "
        "taxa are overwhelmingly obligate anaerobic "
        "and rod-shaped, mechanically coupling the global trait signals."
    ),
    7: (
        "Data Story Plot 07 (Pathway Universe Cloud & Five-Null Hierarchy): Across all 201 "
        "qualifying GO biological processes and "
        "all 90 carbon/energy pathways, Cellulose catabolism (GO:0030245) ranks #1 in "
        "cross-disease Jaccard sharing (44.4% vs. 32.0% "
        "GO median; 2.0x above the 21.8% U4 all-depleted-stool-taxa null, p = 3.4e-6)."
    ),
    8: (
        "Data Story Plot 08 (Barycentric 3-Disease Guild Affinity Simplex): Plotting all 402 "
        "depleted stool taxa in a PD-CRC-HIV "
        "barycentric simplex shows 358 background taxa pinned to single-disease corners while 27 "
        "of 44 cellulose-degrading taxa "
        "(61.4%) pull into the 2- and 3-disease convergence interior (12 in all 3 diseases = "
        "27.3%)."
    ),
    9: (
        "Data Story Plot 09 (Anatomical Compartment Barcode-Strip & Forest): Plotting all 580 "
        "anatomical contrasts as individual "
        "barcode ticks shows that cellulose-degrader depletion is strictly a luminal stool "
        "phenomenon (73.0% Dec in PD+CRC+HIV stool, "
        "63.5% in all 496 stool contrasts) vs. 53.1% in non-fecal GI mucosa and 34.6% in oral "
        "cavity (OR = 5.10, p = 5.3e-6)."
    ),
    10: (
        "Data Story Plot 10 (2,251-Taxon LOTO Jackknife & Assay Platform Strips): Recomputing "
        "trait sign fractions across 2,251 "
        "Leave-One-Taxon-Out (LOTO) runs yields 0 sign flips across all 6 headline traits, and "
        "stratifying 16S rRNA (n = 807) vs. "
        "Shotgun WMS (n = 128) confirms concordant cellulose depletion (62.8% vs. 69.0% Dec in "
        "Feces)."
    ),
}
