"""Generate the Tufte analytical plates, standalone SVGs, and 16-slide HTML deck (#550).

Writes ``data/analysis/cellulose_decr/figures_tufte/`` (Plates I-VI as ``.svg``, ``.pdf``,
``.png``, standalone diagrams ``fig7``-``fig9`` as ``.svg``, and the combined 6-page vector
packet ``cellulose_decr_tufte_figures.pdf``) plus the 16-slide presentation deck
``data/analysis/cellulose_decr/cellulose_decr_slide_deck.html``.
"""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.typing import RcKeyType

from mtp.cellulose_decr.analysis import ANALYSIS_DIR

__all__ = [
    "PLATE2_DOT_X",
    "PLATE2_ROLE_X",
    "PLATE_STEMS",
    "STANDALONE_SVGS",
    "build_plate1_convergence_slope_and_overlap",
    "build_plate2_core_guild_dot_table",
    "build_plate3_why_only_stool",
    "build_plate4_trait_landscape_and_four_pairs",
    "build_plate5_three_statistical_arms",
    "build_plate6_combining_results_integrated_evidence",
    "generate_tufte_plates_and_deck",
    "main",
    "parse_args",
    "render_fig7_svg",
    "render_fig8_svg",
    "render_fig9_svg",
    "render_slide_deck_html",
]

PAPER_BG = "#fbfaf7"
PANEL_BG = "#f4f1ea"
INK = "#181715"
MUTED = "#4a4742"
FAINT = "#78736a"
RULE = "#c8c2b6"
LIGHT_RULE = "#e4dfd5"

CRIMSON = "#9e2a2b"
SLATE = "#2b4c6f"
OCHRE = "#8c5a2b"
TEAL = "#3b6e5c"
WARM_GRAY = "#8c867c"

SERIF_STACK = [
    "Palatino",
    "Georgia",
    "Baskerville",
    "DejaVu Serif",
    "Times New Roman",
    "serif",
]

PLATE_STEMS: tuple[str, ...] = (
    "fig1_tufte_convergence_slope_and_overlap",
    "fig2_tufte_core_guild_dot_table",
    "fig3_tufte_why_only_stool",
    "fig4_tufte_trait_landscape_and_four_pairs",
    "fig5_tufte_three_statistical_arms",
    "fig6_tufte_combining_results_integrated_evidence",
)

STANDALONE_SVGS: tuple[str, ...] = (
    "fig7_tufte_convergence_network_and_venn.svg",
    "fig8_literature_support_four_pairs.svg",
    "fig9_set_theory_17_core_vs_pairwise.svg",
)

PLATE2_ROLE_X: float = 0.25
PLATE2_DOT_X: tuple[float, float, float] = (0.64, 0.68, 0.72)


def _tufte_rcparams() -> dict[RcKeyType, Any]:
    """Return Matplotlib rcParams dictionary without mutating global settings."""
    return {
        "figure.facecolor": PAPER_BG,
        "axes.facecolor": PAPER_BG,
        "savefig.facecolor": PAPER_BG,
        "font.family": "serif",
        "font.serif": SERIF_STACK,
        "text.color": INK,
        "axes.labelcolor": INK,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "axes.edgecolor": RULE,
        "axes.linewidth": 0.65,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "svg.hashsalt": "mtp-cellulose-decr-tufte",
    }


def _make_tufte_plate(
    figsize: tuple[float, float],
    plate_label: str,
    title: str,
    subtitle: str,
    source_left: str,
    source_right: str,
    *,
    header_rule_y: float = 0.865,
    footer_rule_y: float = 0.042,
) -> Figure:
    """Create a figure with Tufte plate header, subtitle, hairlines, and source footer."""
    fig = plt.figure(figsize=figsize, facecolor=PAPER_BG)
    fig.text(
        0.03,
        0.962,
        plate_label,
        fontsize=8.8,
        fontweight="bold",
        color=CRIMSON,
        ha="left",
        va="top",
    )
    fig.text(
        0.03,
        0.930,
        title,
        fontsize=16.2,
        fontweight="bold",
        color=INK,
        ha="left",
        va="top",
    )
    fig.text(
        0.03,
        0.885,
        subtitle,
        fontsize=10.1,
        style="italic",
        color=MUTED,
        ha="left",
        va="top",
    )
    fig.add_artist(
        Line2D(
            [0.03, 0.97],
            [header_rule_y, header_rule_y],
            transform=fig.transFigure,
            color=INK,
            linewidth=0.9,
        )
    )
    fig.add_artist(
        Line2D(
            [0.03, 0.97],
            [footer_rule_y, footer_rule_y],
            transform=fig.transFigure,
            color=RULE,
            linewidth=0.6,
        )
    )
    fig.text(
        0.03,
        max(0.012, footer_rule_y - 0.024),
        source_left,
        fontsize=8.3,
        color=FAINT,
        ha="left",
        va="bottom",
    )
    fig.text(
        0.97,
        max(0.012, footer_rule_y - 0.024),
        source_right,
        fontsize=8.3,
        color=FAINT,
        ha="right",
        va="bottom",
    )
    return fig


def _add_panel_header(
    ax: Axes,
    title: str,
    subtitle: str | None = None,
    *,
    x: float = 0.0,
    y_title: float = 0.985,
    y_subtitle: float = 0.935,
    title_size: float = 11.0,
    subtitle_size: float = 8.8,
) -> None:
    """Add a left-aligned panel title and optional italic subtitle in axes coordinates."""
    ax.text(
        x,
        y_title,
        title,
        transform=ax.transAxes,
        fontsize=title_size,
        fontweight="bold",
        color=INK,
        ha="left",
        va="top",
    )
    if subtitle:
        ax.text(
            x,
            y_subtitle,
            subtitle,
            transform=ax.transAxes,
            fontsize=subtitle_size,
            style="italic",
            color=MUTED,
            ha="left",
            va="top",
        )


def _apply_range_frame(
    ax: Axes,
    *,
    x_bounds: tuple[float, float] | None = None,
    y_bounds: tuple[float, float] | None = None,
    spine_color: str = MUTED,
    linewidth: float = 0.7,
) -> None:
    """Remove top/right spines and trim bottom/left spines to data range."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if x_bounds is not None:
        ax.spines["bottom"].set_visible(True)
        ax.spines["bottom"].set_bounds(x_bounds[0], x_bounds[1])
        ax.spines["bottom"].set_color(spine_color)
        ax.spines["bottom"].set_linewidth(linewidth)
    else:
        ax.spines["bottom"].set_visible(False)
    if y_bounds is not None:
        ax.spines["left"].set_visible(True)
        ax.spines["left"].set_bounds(y_bounds[0], y_bounds[1])
        ax.spines["left"].set_color(spine_color)
        ax.spines["left"].set_linewidth(linewidth)
    else:
        ax.spines["left"].set_visible(False)


def _relax_vertical_positions(
    positions: Sequence[float],
    min_gap: float,
    *,
    iterations: int = 80,
) -> list[float]:
    """Iteratively separate 1D vertical coordinates so adjacent items are >= min_gap apart."""
    if not positions:
        return []
    indexed = sorted(enumerate(positions), key=lambda item: item[1])
    vals = [float(v) for _, v in indexed]
    for _ in range(iterations):
        changed = False
        for idx in range(len(vals) - 1):
            diff = vals[idx + 1] - vals[idx]
            if diff < min_gap:
                shift = 0.5 * (min_gap - diff)
                vals[idx] -= shift
                vals[idx + 1] += shift
                changed = True
        if not changed:
            break
    out = [0.0] * len(positions)
    for (orig_idx, _), val in zip(indexed, vals, strict=True):
        out[orig_idx] = val
    return out


def _draw_tufte_quartile_rule(
    ax: Axes,
    values: Sequence[float],
    y: float,
    *,
    color: str = INK,
    median_size: float = 36.0,
) -> tuple[float, float, float, float, float]:
    """Draw a Tufte quartile-bar-and-whisker rule."""
    arr = np.asarray(values, dtype=float)
    vmin = float(np.min(arr))
    q1 = float(np.percentile(arr, 25))
    med = float(np.median(arr))
    q3 = float(np.percentile(arr, 75))
    vmax = float(np.max(arr))
    ax.plot([vmin, q1], [y, y], color=MUTED, linewidth=0.9, solid_capstyle="butt")
    ax.plot([q3, vmax], [y, y], color=MUTED, linewidth=0.9, solid_capstyle="butt")
    ax.plot([q1, q3], [y, y], color=color, linewidth=3.4, solid_capstyle="butt")
    ax.scatter(
        [med],
        [y],
        s=median_size,
        facecolor=PAPER_BG,
        edgecolor=color,
        linewidth=1.5,
        zorder=5,
    )
    return vmin, q1, med, q3, vmax


def _draw_specification_card(
    ax: Axes,
    *,
    x: float,
    y_top: float,
    width: float,
    height: float,
    accent: str,
    title: str,
    question: str,
    rows: Sequence[tuple[str, str]],
    label_x_offset: float = 0.19,
    pad: float = 0.01,
) -> None:
    """Draw a structured specification card with a left accent bar and key-value rows."""
    ax.add_patch(
        FancyBboxPatch(
            (x, y_top - height),
            width,
            height,
            boxstyle=f"round,pad={pad},rounding_size=0.012",
            facecolor=PANEL_BG,
            edgecolor=RULE,
            linewidth=0.7,
        )
    )
    ax.add_patch(
        Rectangle(
            (x - 0.004, y_top - height),
            0.009,
            height,
            facecolor=accent,
            edgecolor="none",
        )
    )
    ax.text(
        x + 0.014,
        y_top - 0.018,
        title,
        fontsize=9.1,
        fontweight="bold",
        color=accent,
        va="top",
    )
    ax.text(
        x + 0.014,
        y_top - 0.056,
        question,
        fontsize=8.0,
        style="italic",
        color=INK,
        va="top",
    )
    ry = y_top - 0.098
    row_step = (height - 0.115) / max(len(rows), 1)
    for key, val in rows:
        ax.text(
            x + 0.016,
            ry,
            key,
            fontsize=7.7,
            fontweight="bold",
            color=MUTED,
            va="top",
        )
        ax.text(
            x + label_x_offset,
            ry,
            val,
            fontsize=7.7,
            color=INK,
            va="top",
        )
        ry -= row_step


def _save_plate(fig: Figure, output_dir: Path, stem: str, *, dpi: int = 300) -> list[Path]:
    """Save a plate in .svg, .pdf, and .png formats."""
    output_dir.mkdir(parents=True, exist_ok=True)
    saved: list[Path] = []
    for ext in ("svg", "pdf", "png"):
        out_path = output_dir / f"{stem}.{ext}"
        if ext == "png":
            fig.savefig(out_path, dpi=dpi, bbox_inches="tight", facecolor=PAPER_BG)
        else:
            fig.savefig(out_path, bbox_inches="tight", facecolor=PAPER_BG)
        saved.append(out_path)
    return saved


# --------------------------------------------------------------------------- Plate I


def build_plate1_convergence_slope_and_overlap() -> Figure:
    """Build Plate I: Cross-Disease Taxonomic Convergence (Slopegraph + Sharing Tiers)."""
    fig = _make_tufte_plate(
        (13.2, 7.6),
        "PLATE I  ·  CROSS-DISEASE TAXONOMIC CONVERGENCE",
        "Three Distinct Diseases Lose the Same Gut Fiber-Digesting Bacteria",
        "Parkinson\u2019s disease (brain), HIV infection (immune system), and colorectal cancer "
        "(colon tumor) have different causes, yet across 90 stool studies, half of the "
        "fiber-digesting bacteria they lose are the same organisms.",
        "Source: BugSigDB human stool signatures (90 studies: 46 PD, 17 HIV, 27 CRC) "
        "\u00d7 KG-Microbe cellulose catabolic process (GO:0030245).",
        "Plate I  ·  data/analysis/cellulose_decr",
    )

    ax_slope = fig.add_axes((0.05, 0.10, 0.44, 0.72), facecolor=PAPER_BG)
    _add_panel_header(
        ax_slope,
        "A. Share of Lost Bacteria Common to Both Diseases (%)",
        "Slopegraph comparing random chance, all 530 stool bacteria, and the 49 "
        "fiber-digesting bacteria",
        y_title=0.985,
        y_subtitle=0.933,
        title_size=11.2,
        subtitle_size=8.9,
    )

    xs = [0.0, 1.0, 2.0]
    series = [
        ("PD vs. CRC (27/47 shared)", [14.7, 24.4, 57.4], CRIMSON, 2.2, 1.0),
        ("3-pair mean (2.5\u00d7 gut baseline)", [12.3, 20.1, 50.4], INK, 1.8, 0.95),
        ("PD vs. HIV (21/43 shared)", [11.6, 18.0, 48.8], SLATE, 1.6, 0.9),
        ("HIV vs. CRC (18/40 shared)", [10.7, 17.8, 45.0], OCHRE, 1.6, 0.9),
    ]
    for _label, yvals, color, lw, alpha in series:
        ls = "--" if color == INK else "-"
        ax_slope.plot(xs, yvals, color=color, linewidth=lw, linestyle=ls, alpha=alpha, zorder=3)
        ax_slope.scatter(xs, yvals, color=color, s=30, zorder=4)

    right_ys = _relax_vertical_positions([s[1][2] for s in series], min_gap=4.2)
    for (label, yvals, color, _lw, _alpha), ry in zip(series, right_ys, strict=True):
        weight = "bold" if color in (CRIMSON, INK) else "normal"
        ax_slope.text(
            2.06,
            ry,
            f"{yvals[2]:.1f}%  {label}",
            fontsize=8.8,
            fontweight=weight,
            color=color,
            va="center",
            ha="left",
        )

    ax_slope.text(-0.06, 14.7, "15%", fontsize=8.2, color=CRIMSON, ha="right", va="center")
    ax_slope.text(-0.06, 11.2, "11\u201312%", fontsize=8.2, color=FAINT, ha="right", va="center")
    ax_slope.text(1.0, 26.2, "24%", fontsize=8.2, color=CRIMSON, ha="center", va="bottom")
    ax_slope.text(1.0, 16.0, "18%", fontsize=8.2, color=MUTED, ha="center", va="top")

    ax_slope.set_xlim(-0.3, 3.6)
    ax_slope.set_ylim(2.0, 72.0)
    ax_slope.set_yticks([11, 20, 35, 50, 57])
    ax_slope.set_yticklabels(["11%", "20%", "35%", "50%", "57%"], fontsize=8.8)
    ax_slope.set_xticks([])
    _apply_range_frame(ax_slope, y_bounds=(11, 57))

    for x_pos, l1, l2, col, fw in [
        (0.0, "Random", "Chance", MUTED, "normal"),
        (1.0, "All Lost Gut", "Bacteria (530)", MUTED, "normal"),
        (2.0, "Lost Fiber-Digesting", "Bacteria (49)", CRIMSON, "bold"),
    ]:
        ax_slope.text(x_pos, 7.2, l1, fontsize=8.8, fontweight=fw, color=col, ha="center")
        ax_slope.text(x_pos, 4.4, l2, fontsize=8.8, fontweight=fw, color=col, ha="center")

    ax_slope.text(
        0.0,
        -0.02,
        "If fiber-digesting bacteria were lost at the same rate as other stool bacteria,\n"
        "only ~20% would overlap between diseases. Observed overlap is 50.4%\n"
        "(exact hypergeometric p < 10⁻¹⁰ for every pair; random-guild permutation p = 0.011).",
        transform=ax_slope.transAxes,
        fontsize=8.3,
        color=MUTED,
        va="top",
    )

    ax_right = fig.add_axes((0.56, 0.10, 0.41, 0.72), facecolor=PAPER_BG)
    ax_right.axis("off")
    ax_right.set_xlim(0.0, 1.0)
    ax_right.set_ylim(0.0, 1.0)

    _add_panel_header(
        ax_right,
        "B. How the 49 Lost Fiber-Digesting Taxa Are Shared Across Diseases",
        "Distinguishing mutually exclusive sharing tiers from total pairwise overlaps",
        y_title=0.985,
        y_subtitle=0.933,
        title_size=11.2,
        subtitle_size=8.9,
    )

    ax_right.text(
        0.0,
        0.855,
        "1. Partition of all 49 lost fiber-digesting taxa into three mutually exclusive tiers:",
        fontsize=9.2,
        fontweight="bold",
        color=INK,
    )

    bar_y, bar_h = 0.715, 0.085
    w1, w2, w3 = 17 / 49, 15 / 49, 17 / 49
    ax_right.add_patch(
        Rectangle((0.0, bar_y), w1, bar_h, facecolor=CRIMSON, edgecolor=PAPER_BG, lw=1.2)
    )
    ax_right.add_patch(
        Rectangle((w1, bar_y), w2, bar_h, facecolor=OCHRE, edgecolor=PAPER_BG, lw=1.2)
    )
    ax_right.add_patch(
        Rectangle((w1 + w2, bar_y), w3, bar_h, facecolor=PANEL_BG, edgecolor=RULE, lw=1.0)
    )

    ax_right.text(
        w1 / 2,
        bar_y + bar_h * 0.65,
        "Lost in ALL 3 diseases",
        fontsize=8.6,
        fontweight="bold",
        color="#ffffff",
        ha="center",
        va="center",
    )
    ax_right.text(
        w1 / 2,
        bar_y + bar_h * 0.28,
        "17 taxa (34.7%)",
        fontsize=8.6,
        fontweight="bold",
        color="#ffffff",
        ha="center",
        va="center",
    )
    ax_right.text(
        w1 + w2 / 2,
        bar_y + bar_h * 0.65,
        "Lost in ONLY 2 diseases",
        fontsize=8.6,
        fontweight="bold",
        color="#ffffff",
        ha="center",
        va="center",
    )
    ax_right.text(
        w1 + w2 / 2,
        bar_y + bar_h * 0.28,
        "15 taxa (30.6%)",
        fontsize=8.6,
        fontweight="bold",
        color="#ffffff",
        ha="center",
        va="center",
    )
    ax_right.text(
        w1 + w2 + w3 / 2,
        bar_y + bar_h * 0.65,
        "Unique to 1 disease",
        fontsize=8.6,
        color=INK,
        ha="center",
        va="center",
    )
    ax_right.text(
        w1 + w2 + w3 / 2,
        bar_y + bar_h * 0.28,
        "17 taxa (34.7%)",
        fontsize=8.6,
        color=INK,
        ha="center",
        va="center",
    )

    ax_right.text(
        0.0,
        0.650,
        "Shared by at least 2 diseases: 32 of 49 taxa (65.3%)\n"
        "= 17 in all three  +  10 in PD & CRC only  +  4 in PD & HIV only  +  1 in HIV & CRC only",
        fontsize=8.4,
        fontweight="bold",
        color=CRIMSON,
        va="top",
    )

    ax_right.text(
        0.0,
        0.545,
        "2. Total overlap for each disease pair (every pair includes the 17 "
        "three-disease core taxa):",
        fontsize=9.2,
        fontweight="bold",
        color=INK,
    )

    pairs_info = [
        (
            "Parkinson\u2019s & Colorectal cancer",
            17,
            10,
            "= 27 shared of 47 (57.4%)   (+10 in PD & CRC only)",
            0.485,
            OCHRE,
        ),
        (
            "Parkinson\u2019s & HIV infection",
            17,
            4,
            "= 21 shared of 43 (48.8%)   (+4 in PD & HIV only)",
            0.368,
            SLATE,
        ),
        (
            "HIV & Colorectal cancer",
            17,
            1,
            "= 18 shared of 40 (45.0%)   (+1 in HIV & CRC only)",
            0.251,
            TEAL,
        ),
    ]
    scale = 0.018
    for title_str, core_n, extra_n, summary_str, y_top, extra_col in pairs_info:
        ax_right.text(0.0, y_top, title_str, fontsize=9.0, fontweight="bold", color=INK)
        by = y_top - 0.068
        bh = 0.046
        cw = core_n * scale
        ew = extra_n * scale
        ax_right.add_patch(
            Rectangle((0.0, by), cw, bh, facecolor=CRIMSON, edgecolor=PAPER_BG, lw=0.8)
        )
        ax_right.add_patch(
            Rectangle((cw, by), ew, bh, facecolor=extra_col, edgecolor=PAPER_BG, lw=0.8)
        )
        ax_right.text(
            cw / 2,
            by + bh / 2,
            "17 (all 3)",
            fontsize=8.0,
            fontweight="bold",
            color="#ffffff",
            ha="center",
            va="center",
        )
        if extra_n >= 4:
            ax_right.text(
                cw + ew / 2,
                by + bh / 2,
                f"+{extra_n}",
                fontsize=8.0,
                fontweight="bold",
                color="#ffffff",
                ha="center",
                va="center",
            )
        ax_right.text(
            cw + ew + 0.018,
            by + bh / 2,
            summary_str,
            fontsize=8.5,
            color=MUTED,
            va="center",
        )

    ax_right.text(
        0.0,
        0.11,
        "Most shared taxa (17 of 32) decrease in all three diseases, so the 3-disease core (17)\n"
        "is larger than the subset shared by only 2 diseases (15). Every complete 2-disease\n"
        "comparison (27, 21, 18) is larger than the 17-taxon core it contains.",
        fontsize=8.4,
        color=MUTED,
        va="top",
    )
    return fig


# --------------------------------------------------------------------------- Plate II


def build_plate2_core_guild_dot_table() -> Figure:
    """Build Plate II: The Shared Fiber-Degrading Guild at Genus & Species Resolution (#550)."""
    fig = _make_tufte_plate(
        (13.8, 8.0),
        "PLATE II  ·  THE SHARED FIBER-DEGRADING GUILD AT GENUS & SPECIES RESOLUTION",
        "Nine Butyrate- and Fiber-Degrading Genera Decrease Across All Three Diseases",
        "Removing broad family and order labels to compare only specific genera and species "
        "(93 background taxa) strengthens the overlap: 62.5% of lost fiber-digesting genera "
        "and species are shared between Parkinson\u2019s disease and colorectal cancer.",
        "Source: Genus + Species sensitivity analysis (N = 93 stool decreased genera/species; "
        "22 cellulose-annotated). PD\u2013CRC exact p = 5.98 \u00d7 10⁻¹⁰.",
        "Plate II  ·  data/analysis/cellulose_decr",
    )

    ax_table = fig.add_axes((0.02, 0.08, 0.64, 0.76), facecolor=PAPER_BG)
    ax_table.axis("off")
    ax_table.set_xlim(0.0, 1.0)
    ax_table.set_ylim(0.0, 1.0)

    x_taxon = 0.0
    x_role = PLATE2_ROLE_X
    x_pd, x_hiv, x_crc = PLATE2_DOT_X
    x_bar = 0.765

    ax_table.text(
        x_taxon, 0.96, "Taxon (Genus or Species)", fontsize=9.2, fontweight="bold", color=INK
    )
    ax_table.text(
        x_role,
        0.96,
        "Functional Annotation / Colonic Role",
        fontsize=9.2,
        fontweight="bold",
        color=INK,
    )
    ax_table.text(
        x_hiv,
        0.96,
        "PD   HIV  CRC",
        fontsize=8.8,
        fontweight="bold",
        color=INK,
        ha="center",
    )
    ax_table.text(
        x_bar,
        0.96,
        "Stool Studies Reporting Loss",
        fontsize=9.2,
        fontweight="bold",
        color=INK,
    )
    ax_table.plot([0.0, 0.99], [0.942, 0.942], color=INK, linewidth=0.75)

    ax_table.text(
        x_taxon,
        0.912,
        "NINE CORE GENERA REPORTED DECREASED IN ALL THREE DISEASES "
        "(17.1\u00d7 over random expectation, p = 3.5 \u00d7 10⁻¹²)",
        fontsize=8.2,
        fontweight="bold",
        color=CRIMSON,
    )

    genera_rows = [
        ("Roseburia", "Butyrate producer; ferments dietary plant fiber", 18, 2, 5),
        ("Blautia", "Acetate and SCFA fermentative commensal", 12, 2, 5),
        ("Clostridium", "Spore-forming anaerobic cellulose degrader", 7, 4, 5),
        ("Bacteroides", "Polysaccharide-utilizing locus (PUL) specialist", 6, 7, 2),
        ("Lachnospira", "Pectin and plant cell-wall SCFA fermenter", 5, 3, 6),
        ("Ruminococcus", "Primary cellulolytic and resistant-starch guild", 6, 4, 4),
        ("Anaerostipes", "Lactate- and fiber-utilizing butyrate synthesizer", 2, 3, 3),
        ("Butyrivibrio", "Hemicellulose and cellulose butyrate producer", 2, 1, 1),
        ("Eubacterium", "Butyrate-producing strict colonic anaerobe", 1, 1, 2),
    ]

    bar_scale = 0.0058
    y = 0.865
    for genus, role, pd_n, hiv_n, crc_n in genera_rows:
        tot = pd_n + hiv_n + crc_n
        ax_table.plot([0.0, 0.99], [y - 0.025, y - 0.025], color=LIGHT_RULE, linewidth=0.45)
        ax_table.text(
            x_taxon + 0.004,
            y,
            genus,
            fontsize=9.5,
            style="italic",
            fontweight="bold",
            color=INK,
            va="center",
        )
        ax_table.text(x_role, y, role, fontsize=8.5, color=MUTED, va="center")
        for dx, col in ((x_pd, CRIMSON), (x_hiv, OCHRE), (x_crc, SLATE)):
            ax_table.scatter([dx], [y], s=50, color=col, zorder=4)
        w_pd = pd_n * bar_scale
        w_hiv = hiv_n * bar_scale
        w_crc = crc_n * bar_scale
        bh = 0.022
        ax_table.add_patch(
            Rectangle((x_bar, y - bh / 2), w_pd, bh, facecolor=CRIMSON, edgecolor="none")
        )
        ax_table.add_patch(
            Rectangle((x_bar + w_pd, y - bh / 2), w_hiv, bh, facecolor=OCHRE, edgecolor="none")
        )
        ax_table.add_patch(
            Rectangle(
                (x_bar + w_pd + w_hiv, y - bh / 2), w_crc, bh, facecolor=SLATE, edgecolor="none"
            )
        )
        ax_table.text(
            x_bar + w_pd + w_hiv + w_crc + 0.010,
            y,
            f"{tot}  ({pd_n}·{hiv_n}·{crc_n})",
            fontsize=8.3,
            color=INK,
            va="center",
        )
        y -= 0.054

    ax_table.text(
        x_taxon,
        0.372,
        "NAMED CELLULOSE-ANNOTATED SPECIES REPORTED DECREASED IN STOOL",
        fontsize=8.2,
        fontweight="bold",
        color=SLATE,
    )

    # Item 1 (#550): Redundant "(PD & CRC)" and "(CRC)" suffixes removed so role text never
    # collides with the PD/HIV/CRC dot columns.
    species_rows = [
        (
            "Roseburia faecis",
            "Beta-glucan and fiber-degrading butyrate producer",
            (True, False, True),
            7,
            "7 studies",
        ),
        (
            "Anaerostipes hadrus",
            "Fiber- and lactate-fermenting butyrate producer",
            (True, False, True),
            3,
            "3 studies",
        ),
        (
            "Bacteroides cellulosilyticus",
            "Human gut cellulolytic specialist",
            (True, False, True),
            2,
            "2 studies",
        ),
        (
            "Ruminococcus flavefaciens",
            "Cellulosome crystalline-cellulose degrader",
            (True, False, True),
            2,
            "2 studies",
        ),
        (
            "Ruminococcus champanellensis",
            "Human-colon crystalline-cellulose degrader",
            (False, False, True),
            1,
            "1 study",
        ),
    ]

    sy = 0.326
    for sp_name, role, (in_pd, in_hiv, in_crc), study_n, study_label in species_rows:
        ax_table.plot([0.0, 0.99], [sy - 0.023, sy - 0.023], color=LIGHT_RULE, linewidth=0.45)
        ax_table.text(
            x_taxon + 0.004,
            sy,
            sp_name,
            fontsize=9.0,
            style="italic",
            color=INK,
            va="center",
        )
        ax_table.text(x_role, sy, role, fontsize=8.4, color=MUTED, va="center")
        for dx, present, col in (
            (x_pd, in_pd, CRIMSON),
            (x_hiv, in_hiv, OCHRE),
            (x_crc, in_crc, SLATE),
        ):
            if present:
                ax_table.scatter([dx], [sy], s=44, color=col, zorder=4)
            else:
                ax_table.scatter(
                    [dx],
                    [sy],
                    s=32,
                    facecolor=PAPER_BG,
                    edgecolor=RULE,
                    linewidth=0.8,
                    zorder=3,
                )
        sw = study_n * bar_scale
        bh = 0.018
        ax_table.add_patch(
            Rectangle((x_bar, sy - bh / 2), sw, bh, facecolor=SLATE, alpha=0.8, edgecolor="none")
        )
        ax_table.text(
            x_bar + sw + 0.010,
            sy,
            study_label,
            fontsize=8.2,
            color=MUTED,
            va="center",
        )
        sy -= 0.050

    ax_side = fig.add_axes((0.68, 0.08, 0.29, 0.76), facecolor=PANEL_BG)
    ax_side.set_xlim(0.0, 1.0)
    ax_side.set_ylim(0.0, 1.0)
    for spine in ax_side.spines.values():
        spine.set_color(RULE)
        spine.set_linewidth(0.7)
    ax_side.set_xticks([])
    ax_side.set_yticks([])

    ax_side.text(
        0.06,
        0.95,
        "Sensitivity to Taxonomic Rank:\nAll Ranks vs. Genus + Species",
        fontsize=10.1,
        fontweight="bold",
        color=INK,
        va="top",
    )
    ax_side.text(
        0.06,
        0.84,
        "Restricting the background from all\n"
        "taxonomic ranks (N = 189) to only\n"
        "Genus + Species (N = 93) tests\n"
        "whether broad family labels drive\n"
        "the overlap. Similarity increases:",
        fontsize=8.4,
        color=MUTED,
        va="top",
    )

    ax_side.text(
        0.18,
        0.65,
        "All Ranks\n(49 taxa)",
        fontsize=8.3,
        fontweight="bold",
        color=MUTED,
        ha="center",
    )
    ax_side.text(
        0.58,
        0.65,
        "Genus + Species\nOnly (22 taxa)",
        fontsize=8.3,
        fontweight="bold",
        color=CRIMSON,
        ha="center",
    )

    slope_items = [
        ("PD vs. CRC", 57.4, 62.5, CRIMSON, 0.58, 0.635),
        ("3-Pair Mean", 50.4, 52.2, INK, 0.50, 0.520),
        ("PD vs. HIV", 48.8, 47.4, SLATE, 0.48, 0.478),
        ("HIV vs. CRC", 45.0, 46.7, OCHRE, 0.44, 0.440),
    ]
    for label, v1, v2, col, y1, y2 in slope_items:
        ax_side.plot([0.18, 0.58], [y1, y2], color=col, linewidth=1.6)
        ax_side.scatter([0.18, 0.58], [y1, y2], color=col, s=24, zorder=4)
        ax_side.text(0.15, y1, f"{v1:.1f}%", fontsize=8.1, color=MUTED, ha="right", va="center")
        fw = "bold" if col in (CRIMSON, INK) else "normal"
        ax_side.text(
            0.61,
            y2,
            f"{v2:.1f}%  {label}",
            fontsize=8.1,
            fontweight=fw,
            color=col,
            ha="left",
            va="center",
        )

    ax_side.text(
        0.06,
        0.34,
        "Resolution of 16S vs. Shotgun Studies:\n"
        "Amplicon (16S) studies typically\n"
        "report taxa at the genus level\n"
        "(e.g. Roseburia or Ruminococcus).\n"
        "Shotgun metagenomic cohorts\n"
        "identify the specific cellulolytic\n"
        "species listed at bottom left.",
        fontsize=8.3,
        color=MUTED,
        va="top",
    )
    ax_side.text(
        0.06,
        0.11,
        "Stacked bar segments show study\n"
        "counts for PD (red), HIV (ochre),\n"
        "and CRC (slate) in (PD·HIV·CRC).",
        fontsize=7.9,
        color=FAINT,
        va="top",
    )
    return fig


# --------------------------------------------------------------------------- Plate III


def build_plate3_why_only_stool() -> Figure:
    """Build Plate III: Anatomical & Ecological Specificity (#550)."""
    # Item 2 (#550): Cite the real committed table compartment_jaccard_and_direction_stats.tsv
    fig = _make_tufte_plate(
        (13.0, 7.5),
        "PLATE III  ·  ANATOMICAL & ECOLOGICAL SPECIFICITY",
        "Shared Loss Appears in Stool, Not in Tissue Biopsies or the Mouth",
        "Cellulose-fermenting bacteria inhabit the anaerobic lumen of the lower colon where "
        "dietary plant fiber collects. In stool, 73.1% of experiments show a decrease; in colon "
        "wall biopsies and saliva, no directional bias is observed.",
        "Source: Compartment-stratified directional sign tests "
        "(data/analysis/cellulose_decr/compartment_jaccard_and_direction_stats.tsv) "
        "& cross-compartment Jaccard analysis.",
        "Plate III  ·  data/analysis/cellulose_decr",
    )

    ax_left = fig.add_axes((0.03, 0.21, 0.47, 0.62), facecolor=PAPER_BG)
    _add_panel_header(
        ax_left,
        "A. Direction Within Studies: % of Experiments Where Fiber Bacteria Decrease",
        "Comparing stool against colon wall biopsies and saliva relative to a 50% "
        "binomial baseline",
        y_title=0.985,
        y_subtitle=0.920,
        title_size=10.6,
        subtitle_size=8.6,
    )

    ax_left.axvline(50.0, color=MUTED, linestyle="--", linewidth=0.9, zorder=2)

    direction_rows = [
        (
            "Stool (All 3 Diseases Combined)",
            73.1,
            65.1,
            80.1,
            "106 dec vs. 39 inc, p = 2.4 \u00d7 10⁻⁸",
            0.76,
            True,
            CRIMSON,
        ),
        (
            "· Colorectal cancer stool",
            81.1,
            64.8,
            92.0,
            "30 dec vs. 7 inc, p = 1.9 \u00d7 10⁻⁴",
            0.64,
            False,
            CRIMSON,
        ),
        (
            "· HIV infection stool",
            76.0,
            54.9,
            90.6,
            "19 dec vs. 6 inc, p = 0.015",
            0.53,
            False,
            CRIMSON,
        ),
        (
            "· Parkinson\u2019s disease stool",
            68.7,
            57.6,
            78.4,
            "57 dec vs. 26 inc, p = 8.8 \u00d7 10⁻⁴",
            0.42,
            False,
            CRIMSON,
        ),
        (
            "Gut Wall Biopsy (Colon Mucosa)",
            48.1,
            28.7,
            68.1,
            "13 dec vs. 14 inc, p = 1.00",
            0.26,
            True,
            MUTED,
        ),
        (
            "· Colorectal tumor biopsy",
            35.3,
            14.2,
            61.7,
            "6 dec vs. 11 inc, p = 0.33",
            0.15,
            False,
            MUTED,
        ),
        (
            "Mouth & Saliva (Oral Cavity)",
            45.5,
            16.7,
            76.6,
            "5 dec vs. 6 inc, p = 1.00",
            0.04,
            True,
            MUTED,
        ),
    ]

    for label, pct, ci_lo, ci_hi, stats_str, y_pos, is_bold, col in direction_rows:
        fw = "bold" if is_bold else "normal"
        fs = 8.9 if is_bold else 8.4
        ax_left.text(
            0.0,
            y_pos,
            label,
            fontsize=fs,
            fontweight=fw,
            color=INK if is_bold else MUTED,
            va="center",
        )
        ax_left.plot([ci_lo, ci_hi], [y_pos, y_pos], color=col, linewidth=1.6, zorder=3)
        ax_left.scatter([pct], [y_pos], s=42 if is_bold else 30, color=col, zorder=4)
        text_x = max(pct + 3.5, 52.5) if col == CRIMSON else 52.0
        ax_left.text(
            text_x,
            y_pos,
            f"{pct:.1f}%  ({stats_str})",
            fontsize=8.1,
            fontweight=fw if col == CRIMSON else "normal",
            color=col,
            va="center",
        )

    ax_left.set_xlim(0.0, 115.0)
    ax_left.set_ylim(-0.04, 1.0)
    ax_left.set_yticks([])
    ax_left.set_xticks([35.0, 50.0, 68.7, 81.1])
    ax_left.set_xticklabels(["35%", "50%\n(Null)", "68.7%", "81.1%"], fontsize=8.3)
    _apply_range_frame(ax_left, x_bounds=(35.0, 81.1))

    ax_right = fig.add_axes((0.56, 0.21, 0.41, 0.62), facecolor=PAPER_BG)
    _add_panel_header(
        ax_right,
        "B. Cross-Disease Convergence: % of Lost Taxa Shared Across Diseases",
        "Average pairwise overlap (Jaccard similarity) across PD, HIV, and CRC by body site",
        y_title=0.985,
        y_subtitle=0.920,
        title_size=10.6,
        subtitle_size=8.6,
    )

    comp_rows = [
        (
            "Stool (Fecal Lumen)",
            50.4,
            "50.4% shared  ·  17 core taxa in all 3 diseases  (90 studies)",
            "Consistent overlap across all 3 pairs (45%\u201357%)",
            0.72,
            CRIMSON,
        ),
        (
            "Gut Wall Biopsy (Mucosa)",
            18.8,
            "18.8% shared  ·  2 core taxa in all 3 diseases  (15 studies)",
            "Adherent mucosal biofilm differs from lumen",
            0.50,
            MUTED,
        ),
        (
            "Other Sites (Vaginal / Skin / Nasal)",
            3.7,
            "3.7% shared  ·  0 core taxa in all 3 diseases  (7 studies)",
            "Cellulolytic anaerobes rare outside GI tract",
            0.28,
            MUTED,
        ),
        (
            "Mouth & Saliva (Oral Cavity)",
            0.0,
            "0.0% shared  ·  0 core taxa in all 3 diseases  (9 studies)",
            "Zero cross-disease overlap in oral samples",
            0.07,
            MUTED,
        ),
    ]
    for site_label, val, main_note, sub_note, y_pos, col in comp_rows:
        ax_right.text(0.0, y_pos + 0.055, site_label, fontsize=9.2, fontweight="bold", color=INK)
        ax_right.plot([0.0, val], [y_pos, y_pos], color=col, linewidth=3.2, solid_capstyle="butt")
        ax_right.scatter([val], [y_pos], s=48, color=col, zorder=4)
        tx = val + 2.5
        ax_right.text(
            tx,
            y_pos + 0.008,
            main_note,
            fontsize=8.5,
            fontweight="bold" if col == CRIMSON else "normal",
            color=col,
            va="center",
        )
        ax_right.text(
            tx,
            y_pos - 0.036,
            sub_note,
            fontsize=8.0,
            style="italic",
            color=FAINT,
            va="center",
        )

    ax_right.set_xlim(0.0, 108.0)
    ax_right.set_ylim(-0.04, 1.0)
    ax_right.set_yticks([])
    ax_right.set_xticks([0.0, 18.8, 35.0, 50.4])
    ax_right.set_xticklabels(["0%", "18.8%", "35%", "50.4%"], fontsize=8.4)
    _apply_range_frame(ax_right, x_bounds=(0.0, 50.4))

    notes = [
        (
            0.03,
            "1. Colonic Biogeography of Fiber",
            "Insoluble dietary cellulose passes undigested to the anaerobic\n"
            "lumen of the lower colon. Stool samples capture this luminal\n"
            "fermenting community; oral and skin sites do not.",
        ),
        (
            0.36,
            "2. Lumen vs. Gut Wall Biopsy",
            "Colon biopsies sample bacteria attached to the mucosal surface\n"
            "(such as Bacteroides fragilis or Fusobacterium), which often\n"
            "increase on tumor tissue while luminal fiber-degraders drop.",
        ),
        (
            0.68,
            "3. Study Coverage in BugSigDB",
            "Stool is the primary sample type across all three diseases\n"
            "(46 PD, 17 HIV, 27 CRC studies), whereas mucosal biopsies are\n"
            "rarely collected in neurological cohorts (2 PD studies).",
        ),
    ]
    for nx, ntitle, nbody in notes:
        fig.text(nx, 0.145, ntitle, fontsize=8.8, fontweight="bold", color=INK, va="top")
        fig.text(nx, 0.120, nbody, fontsize=8.1, color=MUTED, va="top")

    return fig


# --------------------------------------------------------------------------- Plate IV


_PLATE4_TRAIT_COORDS: tuple[tuple[float, float], ...] = (
    (67.0, 0.56),
    (50.43, 0.44),
    (46.58, 0.59),
    (46.47, 0.53),
    (45.79, 0.32),
    (45.61, 0.63),
    (45.02, 0.55),
    (44.78, 0.56),
    (43.94, 0.33),
    (43.39, 0.45),
    (42.98, 0.42),
    (42.28, 0.61),
    (41.62, 0.51),
    (40.95, 0.57),
    (40.94, 0.44),
    (39.63, 0.37),
    (39.30, 0.48),
    (38.39, 0.31),
    (37.81, 0.58),
    (37.49, 0.51),
    (36.92, 0.55),
    (36.65, 0.41),
    (36.46, 0.63),
    (36.35, 0.60),
    (36.11, 0.56),
    (35.99, 0.36),
    (35.95, 0.45),
    (35.08, 0.30),
    (34.67, 0.34),
    (34.63, 0.53),
    (34.34, 0.55),
    (34.30, 0.63),
    (34.22, 0.40),
    (33.96, 0.42),
    (33.95, 0.45),
    (33.89, 0.35),
    (33.47, 0.33),
    (33.33, 0.45),
    (33.14, 0.37),
    (33.11, 0.52),
    (32.48, 0.53),
    (32.45, 0.40),
    (32.41, 0.58),
    (32.36, 0.42),
    (32.19, 0.53),
    (32.01, 0.36),
    (31.88, 0.29),
    (31.84, 0.56),
    (31.65, 0.53),
    (31.54, 0.45),
    (31.49, 0.49),
    (31.37, 0.33),
    (31.19, 0.49),
    (30.93, 0.56),
    (30.73, 0.51),
    (30.58, 0.48),
    (30.54, 0.30),
    (30.17, 0.44),
    (30.15, 0.36),
    (30.06, 0.58),
    (29.80, 0.37),
    (29.71, 0.39),
    (29.58, 0.52),
    (29.37, 0.56),
    (29.27, 0.43),
    (29.02, 0.30),
    (28.79, 0.54),
    (28.67, 0.34),
    (28.47, 0.53),
    (28.30, 0.42),
    (27.93, 0.51),
    (27.84, 0.32),
    (27.71, 0.60),
    (27.45, 0.38),
    (27.08, 0.56),
    (26.91, 0.45),
    (26.71, 0.32),
    (26.52, 0.36),
    (26.35, 0.35),
    (26.25, 0.55),
    (26.02, 0.44),
    (25.81, 0.51),
    (25.48, 0.52),
    (25.03, 0.30),
    (24.74, 0.40),
    (24.39, 0.49),
    (23.20, 0.61),
    (22.62, 0.41),
    (21.78, 0.62),
    (14.80, 0.46),
)


def build_plate4_trait_landscape_and_four_pairs() -> Figure:
    """Build Plate IV: Functional Pathway Benchmark & Literature-Supported Signatures (#550)."""
    fig = _make_tufte_plate(
        (13.6, 7.8),
        "PLATE IV  ·  FUNCTIONAL PATHWAY BENCHMARK & LITERATURE-SUPPORTED SIGNATURES",
        "Plant-Fiber Breakdown Ranks First Among Metabolic Pathways as Salt-Tolerant Taxa Increase",
        "Across all 90 microbial functional traits tested in stool, the three plant-polysaccharide "
        "pathways (cellulose, xylan, chitin) rank #1, #2, and #3 among metabolic pathways for "
        "cross-disease convergence, while salt-tolerant bacteria increase in "
        "Parkinson\u2019s disease.",
        "Source: 90-trait benchmark (trait_convergence_benchmark.tsv) & 4-pair literature audit "
        "(FOUR_HIGHLIGHTED_PAIRS_LITERATURE_SUPPORT.md).",
        "Plate IV  ·  data/analysis/cellulose_decr",
    )

    ax_top = fig.add_axes((0.03, 0.52, 0.94, 0.30), facecolor=PAPER_BG)
    _add_panel_header(
        ax_top,
        "A. Distribution of Cross-Disease Overlap Across All 90 Functional Traits in Stool",
        "Each mark is one functional trait (90 non-media traits with ≥10 lost taxa per disease), "
        "plotted by % of lost taxa shared across PD, HIV, and CRC",
        y_title=0.985,
        y_subtitle=0.875,
        title_size=11.0,
        subtitle_size=8.8,
    )

    xs = [pt[0] for pt in _PLATE4_TRAIT_COORDS]
    ys = [pt[1] for pt in _PLATE4_TRAIT_COORDS]
    ax_top.scatter(xs, ys, s=22, color=WARM_GRAY, alpha=0.55, zorder=2)
    _draw_tufte_quartile_rule(ax_top, xs, 0.16, color=INK, median_size=34.0)

    ax_top.axvline(20.1, color=SLATE, linestyle="--", linewidth=1.0, zorder=3)
    ax_top.text(
        20.6,
        0.72,
        "All 530 Stool Lost Taxa\nBaseline: 20.1% shared",
        fontsize=8.2,
        color=SLATE,
        va="top",
    )

    focal_traits = [
        (
            50.43,
            0.46,
            CRIMSON,
            "Cellulose breakdown (#1 Metabolic Pathway): 50.4% shared",
            (50.43, 0.84),
            "center",
        ),
        (
            67.00,
            0.46,
            INK,
            "Spore-forming (#1 Physical Trait):\n67.0% shared (same Clostridia)",
            (65.5, 0.68),
            "right",
        ),
        (
            46.58,
            0.46,
            OCHRE,
            "Xylan breakdown (#2 Metabolic Pathway): 46.6% shared",
            (45.2, 0.72),
            "right",
        ),
        (
            45.81,
            0.46,
            SLATE,
            "Chitin breakdown (#3 Metabolic Pathway): 45.8% shared",
            (43.8, 0.59),
            "right",
        ),
    ]
    for fx, fy, fcol, flabel, (lx, ly), ha_align in focal_traits:
        ax_top.scatter([fx], [fy], s=55, color=fcol, zorder=5)
        ax_top.plot([fx, lx], [fy + 0.03, ly - 0.03], color=fcol, linewidth=0.8)
        ax_top.text(
            lx,
            ly,
            flabel,
            fontsize=8.5,
            fontweight="bold" if fcol == CRIMSON else "normal",
            color=fcol,
            ha=ha_align,
            va="bottom",
        )

    ax_top.set_xlim(14.0, 72.0)
    ax_top.set_ylim(0.04, 1.0)
    ax_top.set_yticks([])
    ax_top.set_xticks([14.8, 27.1, 30.2, 33.9, 50.4, 67.0])
    ax_top.set_xticklabels(
        ["Min\n14.8%", "Q1\n27.1%", "Median\n30.2%", "Q3\n33.9%", "Cellulose\n50.4%", "Max\n67.0%"],
        fontsize=8.2,
    )
    _apply_range_frame(ax_top, x_bounds=(14.8, 67.0))

    ax_bot = fig.add_axes((0.03, 0.08, 0.94, 0.38), facecolor=PAPER_BG)
    _add_panel_header(
        ax_bot,
        "B. Four Highlighted Condition-Trait Signatures and Independent Literature Support",
        "Difference in share of stool profiles with positive trait annotation "
        "(Increased minus Decreased; percentage points, pp)",
        y_title=0.985,
        y_subtitle=0.895,
        title_size=11.0,
        subtitle_size=8.8,
    )

    ax_bot.axvline(0.0, color=MUTED, linestyle="--", linewidth=0.9, zorder=2)

    # Item 3 (#550): Cite Nishiwaki 2020 (Mov Disord, PMID 32472608), matching Slide 15 & audit doc.
    four_pairs = [
        (
            "1. Parkinson\u2019s Disease  ·  Cellulose Breakdown Decreased",
            -27.4,
            "\u221227.4 pp  (52 stool studies, BH q = 0.0077; 57 Down vs 26 Up, p = "
            "8.8 \u00d7 10⁻⁴)",
            "Depletion of Roseburia, Lachnospiraceae, Blautia, Ruminococcus; "
            "reported by Wallen 2022, Romano 2021, Nishiwaki 2020.",
            0.70,
            CRIMSON,
        ),
        (
            "2. HIV Infection  ·  Cellulose Breakdown Decreased",
            -53.7,
            "\u221253.7 pp  (18 stool studies, BH q = 0.0080; 19 Down vs 6 Up, p = 0.0146)",
            "Depletion of fiber-fermenting Bacteroides, Clostridia & Lachnospiraceae; "
            "reported by Vujkovic-Cvijin 2013, Noguera-Julian 2016.",
            0.50,
            CRIMSON,
        ),
        (
            "3. Colorectal Cancer  ·  Cellulose Breakdown Decreased",
            -36.9,
            "\u221236.9 pp  (30 stool studies, BH q = 0.0336; 30 Down vs 7 Up, p = "
            "1.9 \u00d7 10⁻⁴)",
            "Depletion of Roseburia, Blautia, Lachnospira, R. champanellensis; "
            "reported by Wirbel 2019, Thomas 2019, Yang 2021.",
            0.30,
            CRIMSON,
        ),
        (
            "4. Parkinson\u2019s Disease  ·  Moderately Halophilic Increased (METPO:1000623)",
            24.4,
            "+24.4 pp  (52 stool studies, BH q = 0.0265; Wilcoxon 33 Up vs 8 Down, p "
            "= 1.7 \u00d7 10⁻⁴)",
            "Increase of osmotolerant/stress-tolerant Lactobacillus, Bifidobacterium, "
            "Enterococcus, Alistipes; Vandeputte 2016, Maini Rekdal 2019.",
            0.10,
            SLATE,
        ),
    ]

    for title_str, delta_pp, stats_line, lit_line, y_pos, col in four_pairs:
        ax_bot.text(-52.0, y_pos + 0.055, title_str, fontsize=8.8, fontweight="bold", color=INK)
        ax_bot.plot(
            [0.0, delta_pp], [y_pos, y_pos], color=col, linewidth=3.0, solid_capstyle="butt"
        )
        ax_bot.scatter([delta_pp], [y_pos], s=48, color=col, zorder=4)
        ax_bot.text(
            28.0, y_pos + 0.020, stats_line, fontsize=8.3, fontweight="bold", color=col, va="center"
        )
        ax_bot.text(28.0, y_pos - 0.038, lit_line, fontsize=7.8, color=MUTED, va="center")

    ax_bot.set_xlim(-56.0, 118.0)
    ax_bot.set_ylim(-0.02, 1.0)
    ax_bot.set_yticks([])
    ax_bot.set_xticks([-53.7, -36.9, -27.4, 0.0, 24.4])
    ax_bot.set_xticklabels(
        ["\u221253.7 pp", "\u221236.9 pp", "\u221227.4 pp", "0 (Balanced)", "+24.4 pp"],
        fontsize=8.1,
    )
    _apply_range_frame(ax_bot, x_bounds=(-53.7, 24.4))
    return fig


# --------------------------------------------------------------------------- Plate V


def build_plate5_three_statistical_arms() -> Figure:
    """Build Plate V: Statistical Architecture (One Matrix, Three Statistical Questions)."""
    fig = _make_tufte_plate(
        (13.6, 8.2),
        "PLATE V  ·  STATISTICAL ARCHITECTURE (METHODS/THREE_ARM_STATISTICS.MD)",
        "One Matrix, Three Statistical Questions: Composition, Enrichment, and Study Direction",
        "All three arms start from the same paired human trait matrix (2,335 experiments across "
        "959 studies and 317 non-MediaDive traits), but test different sampling units, "
        "transformations, and null hypotheses.",
        "Source: methods/three_arm_statistics.md & methods/table1_log_ratio_validation.md "
        "(* q ≤ 0.05; ** q ≤ 0.001).",
        "Plate V  ·  data/analysis/cellulose_decr",
    )

    ax_a = fig.add_axes((0.03, 0.56, 0.50, 0.28), facecolor=PAPER_BG)
    ax_a.axis("off")
    ax_a.set_xlim(0.0, 1.0)
    ax_a.set_ylim(0.0, 1.0)
    _add_panel_header(
        ax_a,
        "A. Shared Input Matrix Filter Cascade (study_trait_profiles_wide.tsv)",
        "BugSigDB v1.3.1 (2026-04-24_00:41_UTC) joined to KG-Microbe trait assertions; "
        "2,653 columns -> 317 non-MediaDive traits",
        y_title=0.985,
        y_subtitle=0.880,
        title_size=10.6,
        subtitle_size=8.3,
    )

    cascade = [
        (
            "1. Raw profile rows read",
            7329,
            "7,329 profiles (all hosts, paired & unpaired)",
            0.72,
            WARM_GRAY,
            False,
        ),
        (
            "2. Drop 1,500 non-human rows",
            5829,
            "5,829 human profiles (3,494 human experiments)",
            0.56,
            SLATE,
            False,
        ),
        (
            "3. Drop 1,159 unpaired rows",
            4670,
            "4,670 paired profiles = 2,335 pairs in 959 studies",
            0.40,
            SLATE,
            True,
        ),
        (
            "4. Drop 94 zero-total pairs",
            4482,
            "4,482 non-zero profiles = 2,241 pairs in 939 studies",
            0.24,
            CRIMSON,
            True,
        ),
    ]
    for step_label, count_val, desc, cy, col, bold in cascade:
        fw = "bold" if bold else "normal"
        ax_a.text(0.0, cy, step_label, fontsize=8.4, fontweight=fw, color=INK, va="center")
        bw = (count_val / 7329.0) * 0.28
        ax_a.add_patch(Rectangle((0.32, cy - 0.042), bw, 0.084, facecolor=col, alpha=0.85))
        ax_a.text(
            0.32 + bw + 0.015,
            cy,
            desc,
            fontsize=8.0,
            fontweight="bold" if col == CRIMSON else "normal",
            color=CRIMSON if col == CRIMSON else MUTED,
            va="center",
        )
    ax_a.text(
        0.0,
        0.06,
        "Column filter: 2,653 matrix columns = 3 ID columns + 2,333 MediaDive growth-medium "
        "columns (excluded) + 317 retained trait columns.",
        fontsize=7.9,
        color=FAINT,
    )

    ax_b = fig.add_axes((0.57, 0.56, 0.40, 0.28), facecolor=PAPER_BG)
    _add_panel_header(
        ax_b,
        "B. Effective Sample Size (n) Varies Across Sampling Units",
        "Log₁₀ scale of observations per test explains why raw p-values cannot be "
        "ranked across arms",
        y_title=0.985,
        y_subtitle=0.880,
        title_size=10.6,
        subtitle_size=8.3,
    )

    sample_units = [
        (
            "Arm 2 (Taxon ORA / Urn)",
            np.log10(4670),
            "4,670 profiles \u00d7 317 traits (1.48M tests; #131)",
            0.72,
            OCHRE,
        ),
        (
            "Arm 1 (Composition)",
            np.log10(2335),
            "2,335 experiment pairs (56.2% multi-exp studies)",
            0.50,
            SLATE,
        ),
        (
            "Arm 3 (6 Fractions)",
            np.log10(917),
            "616\u2013917 non-tied studies (1 vote per study)",
            0.28,
            CRIMSON,
        ),
        (
            "Arm 3 (15 Log-Ratios)",
            np.log10(467),
            "243\u2013676 positive-count studies (unimputed)",
            0.08,
            TEAL,
        ),
    ]
    for label, log_n, note, sy, col in sample_units:
        ax_b.text(1.85, sy + 0.06, label, fontsize=8.3, fontweight="bold", color=INK)
        ax_b.plot([2.0, log_n], [sy, sy], color=col, linewidth=2.8)
        ax_b.scatter([log_n], [sy], s=38, color=col, zorder=4)
        ax_b.text(
            log_n + 0.06,
            sy,
            note,
            fontsize=7.8,
            color=CRIMSON if col == CRIMSON else MUTED,
            va="center",
        )
    ax_b.set_xlim(1.85, 5.2)
    ax_b.set_ylim(-0.04, 1.0)
    ax_b.set_yticks([])
    ax_b.set_xticks(
        [2.0, np.log10(243), np.log10(939), np.log10(2335), np.log10(4670), np.log10(23687)]
    )
    ax_b.set_xticklabels(
        ["100", "243", "939\nstudies", "2,335\npairs", "4,670\nprofiles", "23,687\nassertions"],
        fontsize=7.6,
    )
    _apply_range_frame(ax_b, x_bounds=(2.0, np.log10(23687)))

    ax_c = fig.add_axes((0.03, 0.07, 0.50, 0.44), facecolor=PAPER_BG)
    ax_c.axis("off")
    ax_c.set_xlim(0.0, 1.0)
    ax_c.set_ylim(0.0, 1.0)
    _add_panel_header(
        ax_c,
        "C. Specification of the Three Statistical Arms on the Shared Matrix",
        y_title=0.985,
        title_size=10.6,
    )

    _draw_specification_card(
        ax_c,
        x=0.01,
        y_top=0.91,
        width=0.97,
        height=0.27,
        accent=SLATE,
        title="Arm 1 · Composition Contrast (Jed)",
        question="Question: Within an experiment, is a trait\u2019s share of "
        "annotations higher on the increased side than the decreased side?",
        rows=[
            (
                "Sampling unit & n:",
                "Matched experiment pair (n = 2,335 pairs; identical sample size for every trait)",
            ),
            (
                "Transformation:",
                "Fraction C(e,d,t)/T_317(e,d) \u00d7 1000 -> log(1 + x) -> per-trait "
                "z-score across profiles",
            ),
            (
                "Test & Null:",
                "Two-sided one-sample t-test of paired difference (inc \u2212 dec) against mean 0",
            ),
            (
                "Multiplicity & Rule:",
                "BH within each output table; ranks top 10 traits by p-value "
                "(retains contrast magnitude)",
            ),
            (
                "Scope note:",
                "Harmonized to the shared 317 non-MediaDive trait denominator T_317 "
                "(#174 / #503 / #514)",
            ),
        ],
    )
    _draw_specification_card(
        ax_c,
        x=0.01,
        y_top=0.61,
        width=0.97,
        height=0.27,
        accent=OCHRE,
        title="Arm 2 · Profile Enrichment / Over-Representation Analysis (Ilyas)",
        question="Question: Taking one profile on its own (increased or decreased "
        "separately), does a trait exceed its corpus background rate?",
        rows=[
            (
                "Sampling unit & n:",
                "Taxon-calibrated profile \u00d7 trait hypergeometric model (4,670 "
                "\u00d7 317 = 1,480,390 tests)",
            ),
            (
                "Transformation:",
                "Taxon-calibrated count vs. background rate b(t); fold = observed / expected",
            ),
            (
                "Test & Null:",
                "Hypergeometric upper/lower tail tests (unpaired: tests inc and dec "
                "profiles separately)",
            ),
            (
                "Multiplicity & Rule:",
                "BH across 1,480,390 tests; 66,240 enrichments (fold ≥ 1.5, q ≤ "
                "0.001) & 23 depletions (fold ≤ 0.5)",
            ),
            (
                "Scope note:",
                "Replaces legacy assertion-level urn (mean 23,687 "
                "assertions/profile, 556,516 tests; #131 / #514)",
            ),
        ],
    )
    _draw_specification_card(
        ax_c,
        x=0.01,
        y_top=0.31,
        width=0.97,
        height=0.28,
        accent=CRIMSON,
        title="Arm 3 · Study Directional Consistency & Paired Log-Ratios (Noel)",
        question="Question: Across independent studies, does a trait shift upward "
        "more often than downward when weighted equally per study?",
        rows=[
            (
                "Sampling unit & n:",
                "Unique study (939 non-zero studies; multi-experiment studies "
                "averaged within study first)",
            ),
            (
                "Estimand A (Fractions):",
                "Sign of study-mean Δ_F = p_inc \u2212 p_dec over 317 non-MediaDive "
                "traits (ties |Δ| ≤ 10⁻¹² dropped)",
            ),
            (
                "Estimand B (Log-Ratios):",
                "Sign of study-mean Δ_R = log₂(A_inc/B_inc) \u2212 log₂(A_dec/B_dec) "
                "on positive pairs (cancels T)",
            ),
            (
                "Test & Multiplicity:",
                "Exact two-sided binomial test vs. 0.5 (95% Clopper-Pearson CI); BH "
                "(6 traits) & BY (15 ratios, 317 traits)",
            ),
        ],
    )

    ax_d = fig.add_axes((0.57, 0.07, 0.40, 0.44), facecolor=PAPER_BG)
    _add_panel_header(
        ax_d,
        "D. Arm 3 Empirical Results: Table 1 Traits & Log-Ratios",
        "Share of non-tied studies shifting upward (95% Clopper-Pearson CI vs. 50%)",
        y_title=0.985,
        y_subtitle=0.925,
        title_size=10.6,
        subtitle_size=8.2,
    )
    ax_d.axvline(50.0, color=MUTED, linestyle="--", linewidth=0.8)
    ax_d.text(
        2.0,
        0.85,
        "Six Table 1 Single-Trait Fractions (BH across 6 traits):",
        fontsize=8.2,
        fontweight="bold",
        color=INK,
    )

    arm3_fractions = [
        (
            "FA  (facultative anaerobe)",
            58.2,
            55.0,
            61.4,
            "58.2% (n=917, BH q = 4.1\u00d710⁻⁶ **)",
            0.78,
            SLATE,
            True,
        ),
        (
            "HAL (halotolerant)",
            57.4,
            54.0,
            60.8,
            "57.4% (n=845, BH q = 5.8\u00d710⁻⁵ **)",
            0.71,
            SLATE,
            True,
        ),
        (
            "ARS (arsenate reduction)",
            56.2,
            52.1,
            60.1,
            "56.2% (n=616, BH q = 0.0030 *)",
            0.64,
            SLATE,
            False,
        ),
        (
            "SPH (sphere shaped)",
            55.3,
            51.9,
            58.7,
            "55.3% (n=844, BH q = 0.0030 *)",
            0.57,
            SLATE,
            False,
        ),
        (
            "COC (coccus shaped)",
            54.5,
            51.2,
            57.8,
            "54.5% (n=915, BH q = 0.0067 *)",
            0.50,
            SLATE,
            False,
        ),
        (
            "CEL (cellulose catabolism)",
            43.2,
            39.6,
            46.8,
            "43.2% up / 56.8% down (q = 3.1\u00d710⁻⁴ **)",
            0.43,
            CRIMSON,
            True,
        ),
    ]
    for label, pct, lo, hi, summary_str, ry, col, bold in arm3_fractions:
        fw = "bold" if bold else "normal"
        ax_d.text(2.0, ry, label, fontsize=7.9, fontweight=fw, color=INK, va="center")
        ax_d.plot([lo, hi], [ry, ry], color=col, linewidth=1.4)
        ax_d.scatter([pct], [ry], s=26, color=col, zorder=4)
        ax_d.text(66.5, ry, summary_str, fontsize=7.55, fontweight=fw, color=col, va="center")

    ax_d.text(
        2.0,
        0.34,
        "Key Table 1 Paired Log-Ratios (BY across 15 pairs; cancels T):",
        fontsize=8.2,
        fontweight="bold",
        color=INK,
    )
    arm3_ratios = [
        (
            "FA / CEL  (+0.79 log₂)",
            59.5,
            54.9,
            64.0,
            "59.5% up (n=467, BY q = 0.0011 *)",
            0.27,
            CRIMSON,
            True,
        ),
        (
            "FA / COC  (both up -> flat)",
            50.8,
            46.9,
            54.6,
            "50.8% up (n=676, BY q = 1.00 n.s.)",
            0.20,
            FAINT,
            False,
        ),
        (
            "CEL / SPH (\u22120.54 log₂)",
            40.7,
            36.1,
            45.4,
            "40.7% up / 59.3% down (BY q = 0.0026 *)",
            0.13,
            CRIMSON,
            True,
        ),
        (
            "CEL / ARS (\u22120.85 log₂)",
            36.5,
            31.3,
            41.9,
            "36.5% up / 63.5% down (BY q = 0.0011 *)",
            0.06,
            CRIMSON,
            True,
        ),
    ]
    for label, pct, lo, hi, summary_str, ry, col, bold in arm3_ratios:
        fw = "bold" if bold else "normal"
        ax_d.text(2.0, ry, label, fontsize=7.9, fontweight=fw, color=INK, va="center")
        ax_d.plot([lo, hi], [ry, ry], color=col, linewidth=1.4)
        ax_d.scatter([pct], [ry], s=26, color=col, zorder=4)
        ax_d.text(66.5, ry, summary_str, fontsize=7.55, fontweight=fw, color=col, va="center")

    ax_d.set_xlim(2.0, 102.0)
    ax_d.set_ylim(-0.02, 1.0)
    ax_d.set_yticks([])
    ax_d.set_xticks([30.0, 40.0, 50.0, 60.0, 65.0])
    ax_d.set_xticklabels(["30%", "40%", "50%\n(Null)", "60%", "65%"], fontsize=8.0)
    _apply_range_frame(ax_d, x_bounds=(30.0, 65.0))
    return fig


# --------------------------------------------------------------------------- Plate VI


def build_plate6_combining_results_integrated_evidence() -> Figure:
    """Build Plate VI: Result Integration & Evidence Contract."""
    fig = _make_tufte_plate(
        (13.6, 8.2),
        "PLATE VI  ·  RESULT INTEGRATION & EVIDENCE CONTRACT (METHODS/INTEGRATED_EVIDENCE.MD)",
        "How Results Are Combined: Shared Comparison Basis and the Integrated Evidence Contract",
        "Rather than pooling incompatible p-values across sampling units, the pipeline aligns "
        "corpus denominators and records each method, threshold family, and grid-cell state "
        "in a single 1,912-row table.",
        "Source: methods/integrated_evidence.md & "
        "data/derived/integrated_evidence/integrated_evidence.tsv "
        "(just integrated-evidence).",
        "Plate VI  ·  data/analysis/cellulose_decr",
    )

    ax_a = fig.add_axes((0.03, 0.50, 0.45, 0.34), facecolor=PAPER_BG)
    ax_a.axis("off")
    ax_a.set_xlim(0.0, 1.0)
    ax_a.set_ylim(0.0, 1.0)
    _add_panel_header(
        ax_a,
        "A. Shared Comparison Basis (three-framework-comparison-v1)",
        "Four alignment rules place all three arms on the same data without merging their tests",
        y_title=0.985,
        y_subtitle=0.905,
        title_size=10.6,
        subtitle_size=8.3,
    )

    rules = [
        (
            "1. Identical Matrix & 317-Trait Denominator",
            "All descriptive shares use the 317 non-MediaDive trait columns in "
            "study_trait_profiles_wide.tsv\n"
            "(excluding 2,333 MediaDive columns so growth-medium annotations do not "
            "dilute trait shares).",
            0.80,
            SLATE,
        ),
        (
            "2. Complete Non-Zero Pairs (drop_zero_total_pairs)",
            "Restricts to 2,241 complete human experiment pairs across 939 studies "
            "where both increased and\n"
            "decreased profiles have positive total non-MediaDive annotations (94 "
            "zero-total pairs removed).",
            0.59,
            SLATE,
        ),
        (
            "3. Study-Equal Weighting Before Summary",
            "Because 56.2% of studies contribute multiple experiments, experiment "
            "contrasts are averaged\n"
            "within each study first so every study receives equal weight (1 study = 1 vote).",
            0.38,
            CRIMSON,
        ),
        (
            "4. Unimputed Positive Support for Log-Ratios",
            "Paired log-ratios Δ_R(A,B) are evaluated only on pairs where A_inc, "
            "B_inc, A_dec, B_dec > 0\n"
            "(243\u2013676 studies), avoiding arbitrary pseudocounts while reporting "
            "support alongside effect.",
            0.17,
            TEAL,
        ),
    ]
    for rtitle, rbody, ry, col in rules:
        ax_a.text(0.02, ry, rtitle, fontsize=8.5, fontweight="bold", color=col, va="top")
        ax_a.text(0.02, ry - 0.06, rbody, fontsize=7.7, color=INK, va="top")

    ax_b = fig.add_axes((0.52, 0.50, 0.45, 0.34), facecolor=PAPER_BG)
    ax_b.axis("off")
    ax_b.set_xlim(0.0, 1.0)
    ax_b.set_ylim(0.0, 1.0)
    _add_panel_header(
        ax_b,
        "B. Why Disagreements Across Arms Are Not Contradictions",
        "Different arms can give opposite significance calls on the same trait for "
        "three structural reasons",
        y_title=0.985,
        y_subtitle=0.905,
        title_size=10.6,
        subtitle_size=8.3,
    )

    cases = [
        (
            "Case 1 · Both Sides Enriched vs. Zero Paired Difference",
            "Arm 2 (ORA) tests each side against the corpus background b(t): both "
            "Increased (2.1\u00d7) and\n"
            "Decreased (2.0\u00d7) profiles can be enriched for a common gut trait, "
            "while their within-study\n"
            "difference Δ_F = p_inc \u2212 p_dec is near zero (non-significant in Arms 1 and 3).",
            0.79,
            OCHRE,
        ),
        (
            "Case 2 · Large Magnitude in Few Studies vs. Sign Consistency",
            "Arm 1 (t-test on z-scores) responds to contrast magnitude across 2,335 "
            "experiment pairs.\n"
            "Arm 3 (binomial sign test) discards magnitude and asks only whether "
            ">50% of independent\n"
            "studies shift in the same direction. An outlier-driven mean and a 50/50 "
            "sign split can coexist.",
            0.52,
            SLATE,
        ),
        (
            "Case 3 · Compositional Closure vs. Log-Ratio Cancellation",
            "Both FA (58.2% up, q = 4.1\u00d710⁻⁶) and COC (54.5% up, q = 0.0067) "
            "rise as fractions of T,\n"
            "yet their direct ratio FA/COC is flat (50.8% up, BY q = 1.00). "
            "Conversely, FA/CEL (59.5% up),\n"
            "CEL/SPH (40.7% up), and CEL/ARS (36.5% up) confirm cellulose loss after canceling T.",
            0.25,
            CRIMSON,
        ),
    ]
    for ctitle, cbody, cy, col in cases:
        ax_b.text(0.02, cy, ctitle, fontsize=8.5, fontweight="bold", color=col, va="top")
        ax_b.text(0.02, cy - 0.058, cbody, fontsize=7.65, color=INK, va="top")

    ax_c = fig.add_axes((0.03, 0.07, 0.94, 0.38), facecolor=PAPER_BG)
    ax_c.axis("off")
    ax_c.set_xlim(0.0, 1.0)
    ax_c.set_ylim(0.0, 1.0)
    _add_panel_header(
        ax_c,
        "C. The Integrated Evidence Contract: "
        "data/derived/integrated_evidence/integrated_evidence.tsv (1,912 Rows)",
        "One row per (method, feature, inference_level): keeps correction families separate, "
        "reports both q ≤ 0.05 and q ≤ 0.001, and classifies every blank cell",
        y_title=0.985,
        y_subtitle=0.915,
        title_size=10.8,
        subtitle_size=8.4,
    )

    ax_c.text(
        0.0,
        0.81,
        "1. Pooled Study-Level Method Rows (40 Rows across 5 Method Blocks):",
        fontsize=9.0,
        fontweight="bold",
        color=INK,
    )
    cols = [
        (0.0, "Method Block"),
        (0.14, "Rows"),
        (0.18, "FDR Family"),
        (0.245, "q ≤ 0.05"),
        (0.30, "q ≤ 0.001"),
        (0.345, "Summary Recorded in Contract"),
    ]
    for cx, cname in cols:
        ax_c.text(cx, 0.755, cname, fontsize=7.6, fontweight="bold", color=MUTED)
    ax_c.plot([0.0, 0.54], [0.735, 0.735], color=RULE, linewidth=0.6)

    method_rows = [
        (
            "table1_Burst_direction",
            "6 rows",
            "BH (m = 6)",
            "6 / 6",
            "3 / 6",
            "All 6 traits pass q ≤ 0.05; FA, HAL, CEL pass q ≤ 0.001",
            CRIMSON,
            True,
        ),
        (
            "table1_ora_background",
            "6 rows",
            "Descriptive",
            "untested",
            "untested",
            "Study-mean inc & dec ORA folds vs. corpus background",
            OCHRE,
            False,
        ),
        (
            "table1_paired_log_ratio",
            "15 rows",
            "BY (m = 15)",
            "3 / 15",
            "0 / 15",
            "FA/CEL (up), CEL/SPH (down), CEL/ARS (down) pass BY",
            TEAL,
            True,
        ),
        (
            "all_traits_direction_screen",
            "3 rows",
            "BH & BY (317)",
            "42 BH / 19 BY",
            "13 BY",
            "Full 317-trait screen (939 studies; 16 zero-signal traits)",
            SLATE,
            True,
        ),
        (
            "jed_composition_top10",
            "10 rows",
            "BH (table)",
            "pending",
            "pending",
            "Marked source_not_supplied until export is wired (#174)",
            FAINT,
            False,
        ),
    ]
    ry = 0.690
    for mname, rcount, fdr, q05, q001, summ, mcol, bold in method_rows:
        ax_c.text(0.0, ry, mname, fontsize=7.3, fontweight="bold", family="monospace", color=mcol)
        ax_c.text(0.14, ry, rcount, fontsize=7.5, color=INK)
        ax_c.text(0.18, ry, fdr, fontsize=7.5, color=INK)
        fw = "bold" if bold else "normal"
        ax_c.text(0.245, ry, q05, fontsize=7.5, fontweight=fw, color=INK)
        ax_c.text(0.30, ry, q001, fontsize=7.5, fontweight=fw, color=INK)
        ax_c.text(0.345, ry, summ, fontsize=7.35, color=MUTED)
        ry -= 0.095

    ax_c.text(
        0.01,
        0.18,
        "Five Explicit Methodological Standing Columns on Every Row (Never Implied "
        "by Significance):",
        fontsize=7.9,
        fontweight="bold",
        color=INK,
    )
    ax_c.text(
        0.01,
        0.11,
        "· calibration_status (#182: uncalibrated_exploratory)   "
        "· sensitivity_status (#183: pending_common_sample_check)   · coverage_status (#184)\n"
        "· orientation_status (#188: canonical_6_trait vs. historical_5_trait_display)   "
        "· taxon_influence_status (#190: unassessed_rollup_influence)",
        fontsize=7.35,
        color=MUTED,
        va="top",
    )

    ax_c.text(
        0.58,
        0.81,
        "2. Condition \u00d7 Body-Site Grid Cells (1,872 Rows = 6 Traits \u00d7 13 "
        "Conditions \u00d7 12 Sites \u00d7 2 Metrics):",
        fontsize=8.8,
        fontweight="bold",
        color=INK,
    )
    ax_c.text(
        0.58,
        0.745,
        "Replaces ambiguous blank heatmap cells with explicit test_status and "
        "display_status values:",
        fontsize=7.9,
        style="italic",
        color=MUTED,
    )

    gx, gy, gw, gh = 0.58, 0.59, 0.40, 0.10
    w_drawn = (316 / 1872) * gw
    w_zero = (7 / 1872) * gw
    w_mask = (209 / 1872) * gw
    w_ns = (1340 / 1872) * gw
    ax_c.add_patch(Rectangle((gx, gy), w_drawn, gh, facecolor=CRIMSON, edgecolor=PAPER_BG))
    ax_c.add_patch(Rectangle((gx + w_drawn, gy), w_zero, gh, facecolor=OCHRE, edgecolor=PAPER_BG))
    ax_c.add_patch(
        Rectangle((gx + w_drawn + w_zero, gy), w_mask, gh, facecolor=SLATE, edgecolor=PAPER_BG)
    )
    ax_c.add_patch(
        Rectangle(
            (gx + w_drawn + w_zero + w_mask, gy),
            w_ns,
            gh,
            facecolor=PANEL_BG,
            edgecolor=RULE,
        )
    )
    ax_c.text(
        gx + w_drawn / 2,
        gy + gh / 2,
        "316\ndrawn",
        fontsize=7.4,
        fontweight="bold",
        color="#ffffff",
        ha="center",
        va="center",
    )
    ax_c.text(
        gx + w_drawn + w_zero + w_mask / 2,
        gy + gh / 2,
        "209\nmasked",
        fontsize=7.2,
        fontweight="bold",
        color="#ffffff",
        ha="center",
        va="center",
    )
    ax_c.text(
        gx + w_drawn + w_zero + w_mask + w_ns / 2,
        gy + gh / 2,
        "1,340 not_sampled (71.6%: 0 paired human studies in cell)",
        fontsize=7.6,
        fontweight="bold",
        color=INK,
        ha="center",
        va="center",
    )

    cell_states = [
        (
            "drawn (316 cells, 16.9%)   [test_status = tested / all_tied]",
            "≥ 3 studies with trait signal; rendered on map (26 cells have BH q ≤ "
            "0.05 '*', 4 have q ≤ 0.001)",
            0.51,
            CRIMSON,
        ),
        (
            "sampled-without-signal (7 cells, 0.4%)   [test_status = zero_annotation_signal]",
            "≥ 3 studies sampled in BugSigDB, but all profiles have 0 annotations for that trait",
            0.37,
            OCHRE,
        ),
        (
            "masked_below_min_studies (209 cells, 11.2%)   [test_status = insufficient_support]",
            "1\u20132 studies in cell; suppressed by the ≥ 3-study display rule (never tested)",
            0.23,
            SLATE,
        ),
        (
            "not_sampled (1,340 cells, 71.6%)   [test_status = not_sampled]",
            "0 paired human studies for that condition \u00d7 body site in the pinned corpus",
            0.09,
            INK,
        ),
    ]
    for stitle, sbody, sy, col in cell_states:
        ax_c.text(0.595, sy, stitle, fontsize=7.8, fontweight="bold", color=col)
        ax_c.text(0.595, sy - 0.045, sbody, fontsize=7.35, color=MUTED)

    return fig


# ---------------------------------------------------------------------- Standalone SVGs (fig7-fig9)


def render_fig7_svg() -> str:
    """Return the standalone SVG string for Plate VII
    (fig7_tufte_convergence_network_and_venn.svg).
    """
    return "\n".join(
        [
            '<?xml version="1.0" encoding="utf-8" standalone="no"?>',
            '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN"',
            '  "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">',
            '<svg xmlns:xlink="http://www.w3.org/1999/xlink" width="1400pt" height="770pt'
            '" viewBox="0 0 1400 770" xmlns="http://www.w3.org/2000/svg" version="1.1">',
            "  <metadata>",
            '    <rdf:RDF xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:cc="http://crea'
            'tivecommons.org/ns#" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">',
            "      <cc:Work>",
            '        <dc:type rdf:resource="http://purl.org/dc/dcmitype/StillImage"/>',
            "        <dc:title>Plate VII: Cross-Disease Taxonomic Convergence Network and 3-Wa"
            "y Venn Partition</dc:title>",
            "        <dc:format>image/svg+xml</dc:format>",
            "      </cc:Work>",
            "    </rdf:RDF>",
            "  </metadata>",
            "  <defs>",
            '    <style type="text/css">',
            "      text {",
            "        font-family: 'Palatino', 'Palatino Linotype', 'Book Antiqua', 'Georgia', "
            "'Baskerville', 'DejaVu Serif', 'Times New Roman', serif;",
            "      }",
            "      .kicker { font-size: 11.2px; font-weight: 700; fill: #9E2A2B; letter-spacin"
            "g: 1.1px; }",
            "      .title { font-size: 21.5px; font-weight: 700; fill: #181715; }",
            "      .subtitle { font-size: 13px; font-style: italic; fill: #4A4742; }",
            "      .panel-title { font-size: 15.2px; font-weight: 700; fill: #181715; }",
            "      .panel-sub { font-size: 11.5px; font-style: italic; fill: #4A4742; }",
            "      .sub-header { font-size: 12.2px; font-weight: 700; fill: #181715; }",
            "      .footnote { font-size: 10px; fill: #78736A; }",
            "    </style>",
            "  </defs>",
            "",
            '  <rect x="0" y="0" width="1400" height="770" fill="#FCFBF7"/>',
            "",
            '  <text x="32" y="25" class="kicker">PLATE VII  ·  CROSS-DISEASE TAXONOMIC '
            "CONVERGENCE &amp; 3-WAY VENN PARTITION (METPO:2000103 | GO:0030245)</text>",
            '  <text x="32" y="52" class="title">Three Etiologically Distinct Diseases L'
            "ose the Same Fecal Cellulose-Annotated Guild</text>",
            '  <text x="32" y="75" class="subtitle">Across neurodegenerative (PD, 46 stu'
            "dies), viral/immune (HIV, 17 studies), and neoplastic (CRC, 27 studies) stool coh"
            "orts, 65.3% of 49 decreased cellulose-annotated taxa are shared by ≥2 diseases.</"
            "text>",
            '  <line x1="32" y1="91" x2="1368" y2="91" stroke="#181715" stroke-width'
            '="0.85" stroke-opacity="0.85"/>',
            "",
            '  <line x1="606" y1="106" x2="606" y2="716" stroke="#D6D1C7" stroke-wid'
            'th="0.9" stroke-dasharray="3,3"/>',
            "",
            '  <g id="panel_a">',
            '    <text x="32" y="122" class="panel-title">A. Pairwise Taxonomic Converge'
            "nce Network (Feces Decreased)</text>",
            '    <text x="32" y="141" class="panel-sub">Edge weight proportional to Jacc'
            "ard similarity J = |A ∩ B| / |A \u222a B|; every pairwise overlap contains the 17"
            "-taxon 3-way core</text>",
            "",
            '    <line x1="148" y1="262" x2="472" y2="262" stroke="#9E2A2B" stroke-o'
            'pacity="0.18" stroke-width="15" stroke-linecap="round"/>',
            '    <line x1="148" y1="262" x2="472" y2="262" stroke="#181715" stroke-o'
            'pacity="0.85" stroke-width="2.2" stroke-linecap="round"/>',
            "",
            '    <line x1="148" y1="262" x2="310" y2="468" stroke="#8C5A2B" stroke-o'
            'pacity="0.18" stroke-width="12.6" stroke-linecap="round"/>',
            '    <line x1="148" y1="262" x2="310" y2="468" stroke="#181715" stroke-o'
            'pacity="0.80" stroke-width="1.9" stroke-linecap="round"/>',
            "",
            '    <line x1="310" y1="468" x2="472" y2="262" stroke="#2B4C6F" stroke-o'
            'pacity="0.18" stroke-width="11.6" stroke-linecap="round"/>',
            '    <line x1="310" y1="468" x2="472" y2="262" stroke="#181715" stroke-o'
            'pacity="0.78" stroke-width="1.8" stroke-linecap="round"/>',
            "",
            '    <line x1="310" y1="262" x2="310" y2="221" stroke="#9E2A2B" stroke-w'
            'idth="0.9" stroke-dasharray="2,2"/>',
            '    <line x1="229" y1="365" x2="218" y2="365" stroke="#8C5A2B" stroke-w'
            'idth="0.9" stroke-dasharray="2,2"/>',
            '    <line x1="391" y1="365" x2="402" y2="365" stroke="#2B4C6F" stroke-w'
            'idth="0.9" stroke-dasharray="2,2"/>',
            "",
            '    <g id="node_pd">',
            '      <circle cx="148" cy="262" r="57" fill="#F9ECEC" stroke="#9E2A2B" '
            'stroke-width="2.1"/>',
            '      <circle cx="148" cy="262" r="52" fill="none" stroke="#9E2A2B" str'
            'oke-width="0.6" stroke-opacity="0.45"/>',
            '      <text x="148" y="241" text-anchor="middle" font-size="13.2" font-we'
            'ight="700" fill="#9E2A2B">Parkinson\u2019s</text>',
            '      <text x="148" y="257" text-anchor="middle" font-size="13.2" font-we'
            'ight="700" fill="#9E2A2B">Disease (PD)</text>',
            '      <text x="148" y="274" text-anchor="middle" font-size="10.2" font-st'
            'yle="italic" fill="#4A4742">Neurodegenerative</text>',
            '      <text x="148" y="291" text-anchor="middle" font-size="10.5" font-we'
            'ight="700" fill="#181715">40 taxa · 46 studies</text>',
            "    </g>",
            "",
            '    <g id="node_crc">',
            '      <circle cx="472" cy="262" r="57" fill="#EBF1F5" stroke="#2B4C6F" '
            'stroke-width="2.1"/>',
            '      <circle cx="472" cy="262" r="52" fill="none" stroke="#2B4C6F" str'
            'oke-width="0.6" stroke-opacity="0.45"/>',
            '      <text x="472" y="241" text-anchor="middle" font-size="13.2" font-we'
            'ight="700" fill="#2B4C6F">Colorectal</text>',
            '      <text x="472" y="257" text-anchor="middle" font-size="13.2" font-we'
            'ight="700" fill="#2B4C6F">Cancer (CRC)</text>',
            '      <text x="472" y="274" text-anchor="middle" font-size="10.2" font-st'
            'yle="italic" fill="#4A4742">Neoplastic</text>',
            '      <text x="472" y="291" text-anchor="middle" font-size="10.5" font-we'
            'ight="700" fill="#181715">34 taxa · 27 studies</text>',
            "    </g>",
            "",
            '    <g id="node_hiv">',
            '      <circle cx="310" cy="468" r="57" fill="#F7EFE6" stroke="#8C5A2B" '
            'stroke-width="2.1"/>',
            '      <circle cx="310" cy="468" r="52" fill="none" stroke="#8C5A2B" str'
            'oke-width="0.6" stroke-opacity="0.45"/>',
            '      <text x="310" y="447" text-anchor="middle" font-size="13.2" font-we'
            'ight="700" fill="#8C5A2B">HIV Infection</text>',
            '      <text x="310" y="463" text-anchor="middle" font-size="13.2" font-we'
            'ight="700" fill="#8C5A2B">(HIV)</text>',
            '      <text x="310" y="480" text-anchor="middle" font-size="10.2" font-st'
            'yle="italic" fill="#4A4742">Viral / Immune</text>',
            '      <text x="310" y="497" text-anchor="middle" font-size="10.5" font-we'
            'ight="700" fill="#181715">24 taxa · 17 studies</text>',
            "    </g>",
            "",
            '    <g id="edge_pd_crc">',
            '      <rect x="213" y="161" width="194" height="58" rx="2" fill="#FCFB'
            'F7" stroke="#9E2A2B" stroke-width="1.1"/>',
            '      <text x="310" y="177" text-anchor="middle" font-size="11.2" font-we'
            'ight="700" fill="#181715">PD ↔ CRC: 27 / 47 shared</text>',
            '      <text x="310" y="193" text-anchor="middle" font-size="10.5" fill="'
            '#181715"><tspan font-weight="700" fill="#9E2A2B">J = 57.4%</tspan>  (d_J = 0'
            ".426 · 17 + 10)</text>",
            '      <text x="310" y="209" text-anchor="middle" font-size="9.5" fill="#'
            '4A4742">3.8\u00d7 vs. U₁ (p = 9.6\u00d710⁻¹⁷; hyp p = 2.7\u00d710⁻¹³)</text>',
            "    </g>",
            "",
            '    <g id="edge_pd_hiv">',
            '      <rect x="30" y="340" width="188" height="58" rx="2" fill="#FCFBF'
            '7" stroke="#8C5A2B" stroke-width="1.1"/>',
            '      <text x="124" y="356" text-anchor="middle" font-size="11.2" font-we'
            'ight="700" fill="#181715">PD ↔ HIV: 21 / 43 shared</text>',
            '      <text x="124" y="372" text-anchor="middle" font-size="10.5" fill="'
            '#181715"><tspan font-weight="700" fill="#8C5A2B">J = 48.8%</tspan>  (d_J = 0'
            ".512 · 17 + 4)</text>",
            '      <text x="124" y="388" text-anchor="middle" font-size="9.5" fill="#'
            '4A4742">4.1\u00d7 vs. U₁ (p = 4.8\u00d710⁻¹⁴; hyp p = 8.4\u00d710⁻¹²)</text>',
            "    </g>",
            "",
            '    <g id="edge_hiv_crc">',
            '      <rect x="402" y="340" width="188" height="58" rx="2" fill="#FCFB'
            'F7" stroke="#2B4C6F" stroke-width="1.1"/>',
            '      <text x="496" y="356" text-anchor="middle" font-size="11.2" font-we'
            'ight="700" fill="#181715">HIV ↔ CRC: 18 / 40 shared</text>',
            '      <text x="496" y="372" text-anchor="middle" font-size="10.5" fill="'
            '#181715"><tspan font-weight="700" fill="#2B4C6F">J = 45.0%</tspan>  (d_J = 0'
            ".550 · 17 + 1)</text>",
            '      <text x="496" y="388" text-anchor="middle" font-size="9.5" fill="#'
            '4A4742">4.2\u00d7 vs. U₁ (p = 2.6\u00d710⁻¹¹; hyp p = 6.7\u00d710⁻¹¹)</text>',
            "    </g>",
            "",
            '    <g id="centroid_mean">',
            '      <rect x="234" y="283" width="152" height="56" rx="2" fill="#F4F1'
            'EA" stroke="#9E2A2B" stroke-width="1.3"/>',
            '      <text x="310" y="299" text-anchor="middle" font-size="9.4" font-wei'
            'ght="700" fill="#9E2A2B" letter-spacing="0.6">MEAN ACROSS 3 PAIRS</text>',
            '      <text x="310" y="317" text-anchor="middle" font-size="11.8" font-we'
            'ight="700" fill="#181715">J = 50.4%  (d_J = 0.496)</text>',
            '      <text x="310" y="332" text-anchor="middle" font-size="9.6" font-sty'
            'le="italic" fill="#4A4742">97.8th pct of 90 stool traits</text>',
            "    </g>",
            "",
            '    <line x1="32" y1="538" x2="588" y2="538" stroke="#D6D1C7" stroke-width="0.8"/>',
            '    <text x="32" y="558" class="sub-header">Mean Pairwise Jaccard Distance '
            "Benchmark (d_J = 1 \u2212 J; lower = tighter cross-disease convergence)</text>",
            "",
            '    <line x1="264.9" y1="574" x2="264.9" y2="686" stroke="#9E2A2B" stro'
            'ke-width="0.7" stroke-dasharray="2,2" stroke-opacity="0.6"/>',
            '    <line x1="369.9" y1="574" x2="369.9" y2="686" stroke="#D6D1C7" stro'
            'ke-width="0.7" stroke-dasharray="2,2"/>',
            '    <line x1="422.5" y1="574" x2="422.5" y2="686" stroke="#D6D1C7" stro'
            'ke-width="0.7" stroke-dasharray="2,2"/>',
            '    <line x1="458.9" y1="574" x2="458.9" y2="686" stroke="#D6D1C7" stro'
            'ke-width="0.7" stroke-dasharray="2,2"/>',
            "",
            '    <text x="32" y="592" font-size="10.8" font-weight="700" fill="#9E2A2'
            'B">Cellulose Catabolism (N = 49)</text>',
            '    <line x1="228.5" y1="588" x2="293.0" y2="588" stroke="#9E2A2B" stro'
            'ke-width="2.2"/>',
            '    <line x1="228.5" y1="584" x2="228.5" y2="592" stroke="#9E2A2B" stro'
            'ke-width="1.2"/>',
            '    <line x1="293.0" y1="584" x2="293.0" y2="592" stroke="#9E2A2B" stro'
            'ke-width="1.2"/>',
            '    <circle cx="264.9" cy="588" r="5" fill="#9E2A2B" stroke="#FCFBF7" s'
            'troke-width="1"/>',
            '    <text x="301" y="592" font-size="10.2" fill="#181715"><tspan font-wei'
            'ght="700" fill="#9E2A2B">0.496</tspan> (J = 50.4%; pair range 0.426\u20130.55'
            "0)</text>",
            "",
            '    <line x1="264.9" y1="611" x2="422.5" y2="611" stroke="#9E2A2B" stro'
            'ke-width="1.1"/>',
            '    <polygon points="264.9,611 271.5,607.8 271.5,614.2" fill="#9E2A2B"/>',
            '    <text x="343.7" y="607" text-anchor="middle" font-size="9.5" font-wei'
            'ght="700" fill="#9E2A2B">Δd_J = \u22120.303 tighter vs. background (U₄ p &lt;'
            " 10⁻⁵)</text>",
            "",
            '    <text x="32" y="636" font-size="10.5" fill="#2B4C6F" font-weight="70'
            '0">90-Trait Stool Guilds (Median, IQR)</text>',
            '    <rect x="350.7" y="628" width="35.4" height="8" fill="#2B4C6F" fill'
            '-opacity="0.25" stroke="#2B4C6F" stroke-width="0.8"/>',
            '    <circle cx="369.9" cy="632" r="4.2" fill="#2B4C6F"/>',
            '    <text x="394" y="636" font-size="10" fill="#4A4742"><tspan font-weigh'
            't="700" fill="#2B4C6F">0.698</tspan> (J = 30.2%, guild perm p = 0.0110)</text'
            ">",
            "",
            '    <text x="32" y="668" font-size="10.5" fill="#4A4742">All Feces Decrea'
            "sed Taxa (N = 530)</text>",
            '    <circle cx="422.5" cy="664" r="4.2" fill="#4A4742"/>',
            '    <circle cx="458.9" cy="664" r="4.0" fill="#FCFBF7" stroke="#78736A"'
            ' stroke-width="1.4"/>',
            '    <text x="468" y="668" font-size="9.8" fill="#4A4742"><tspan font-weig'
            'ht="700" fill="#181715">0.799</tspan> (J=20.1%) · ○ U₁ <tspan font-weight="7'
            '00">0.869</tspan></text>',
            "",
            '    <line x1="228.5" y1="686" x2="458.9" y2="686" stroke="#181715" stro'
            'ke-width="0.9"/>',
            '    <line x1="228.5" y1="686" x2="228.5" y2="691" stroke="#181715" stro'
            'ke-width="0.8"/>',
            '    <text x="228.5" y="703" text-anchor="middle" font-size="9" fill="#4A'
            '4742">0.43</text>',
            '    <line x1="264.9" y1="686" x2="264.9" y2="691" stroke="#9E2A2B" stro'
            'ke-width="1.0"/>',
            '    <text x="264.9" y="703" text-anchor="middle" font-size="9.2" font-wei'
            'ght="700" fill="#9E2A2B">0.496</text>',
            '    <line x1="369.9" y1="686" x2="369.9" y2="691" stroke="#2B4C6F" stro'
            'ke-width="0.8"/>',
            '    <text x="369.9" y="703" text-anchor="middle" font-size="9" fill="#2B'
            '4C6F">0.698</text>',
            '    <line x1="422.5" y1="686" x2="422.5" y2="691" stroke="#181715" stro'
            'ke-width="0.8"/>',
            '    <text x="422.5" y="703" text-anchor="middle" font-size="9" fill="#4A'
            '4742">0.799</text>',
            '    <line x1="458.9" y1="686" x2="458.9" y2="691" stroke="#181715" stro'
            'ke-width="0.8"/>',
            '    <text x="458.9" y="703" text-anchor="middle" font-size="9" fill="#78'
            '736A">0.869</text>',
            "  </g>",
            "",
            '  <g id="panel_b">',
            '    <text x="626" y="122" class="panel-title">B. 49-Taxon 3-Way Venn Partit'
            "ion &amp; Micro/Macro Taxonomic Breakdown</text>",
            '    <text x="626" y="141" class="panel-sub">Exact partition of the 49 uniqu'
            "e cellulose-annotated taxa decreased in Feces into 7 mutually exclusive subsets w"
            "ith constituent taxa</text>",
            "",
            '    <circle cx="772" cy="302" r="118" fill="#9E2A2B" fill-opacity="0.08'
            '" stroke="#9E2A2B" stroke-width="1.8"/>',
            '    <circle cx="900" cy="302" r="118" fill="#2B4C6F" fill-opacity="0.08'
            '" stroke="#2B4C6F" stroke-width="1.8"/>',
            '    <circle cx="836" cy="410" r="118" fill="#8C5A2B" fill-opacity="0.08'
            '" stroke="#8C5A2B" stroke-width="1.8"/>',
            "",
            '    <text x="722" y="168" text-anchor="middle" font-size="12" font-weight'
            '="700" fill="#9E2A2B">Parkinson\u2019s Disease (40)</text>',
            '    <text x="722" y="181" text-anchor="middle" font-size="9.5" font-style'
            '="italic" fill="#4A4742">46 studies · Δ = \u221227.4 pp (q = 0.0077)</text>',
            "",
            '    <text x="950" y="168" text-anchor="middle" font-size="12" font-weight'
            '="700" fill="#2B4C6F">Colorectal Cancer (34)</text>',
            '    <text x="950" y="181" text-anchor="middle" font-size="9.5" font-style'
            '="italic" fill="#4A4742">27 studies · Δ = \u221236.9 pp (q = 0.0336)</text>',
            "",
            '    <text x="836" y="544" text-anchor="middle" font-size="12" font-weight'
            '="700" fill="#8C5A2B">HIV Infection (24)  <tspan font-size="9.5" font-weigh'
            't="400" font-style="italic" fill="#4A4742">· 17 studies · Δ = \u221253.7 pp'
            " (q = 0.0080)</tspan></text>",
            "",
            '    <text x="710" y="256" text-anchor="middle" font-size="16" font-weight'
            '="700" fill="#9E2A2B">9</text>',
            '    <text x="710" y="270" text-anchor="middle" font-size="9.2" font-weigh'
            't="700" fill="#9E2A2B">PD only (18.4%)</text>',
            '    <text x="710" y="283" text-anchor="middle" font-size="8.0" font-style'
            '="italic" fill="#4A4742">Lacrimispora, Bacillaceae,</text>',
            '    <text x="710" y="294" text-anchor="middle" font-size="8.0" font-style'
            '="italic" fill="#4A4742">Mediterraneibacter,</text>',
            '    <text x="710" y="305" text-anchor="middle" font-size="8.0" font-style'
            '="italic" fill="#4A4742">Eubacteriaceae +5</text>',
            "",
            '    <text x="962" y="256" text-anchor="middle" font-size="16" font-weight'
            '="700" fill="#2B4C6F">6</text>',
            '    <text x="962" y="270" text-anchor="middle" font-size="9.2" font-weigh'
            't="700" fill="#2B4C6F">CRC only (12.2%)</text>',
            '    <text x="962" y="283" text-anchor="middle" font-size="8.0" font-style'
            '="italic" fill="#4A4742">L. multipara, B. obeum,</text>',
            '    <text x="962" y="294" text-anchor="middle" font-size="8.0" font-style'
            '="italic" fill="#4A4742">R. champanellensis,</text>',
            '    <text x="962" y="305" text-anchor="middle" font-size="8.0" font-style'
            '="italic" fill="#4A4742">H. alba, Serratia +1</text>',
            "",
            '    <text x="836" y="464" text-anchor="middle" font-size="16" font-weight'
            '="700" fill="#8C5A2B">2</text>',
            '    <text x="836" y="478" text-anchor="middle" font-size="9.2" font-weigh'
            't="700" fill="#8C5A2B">HIV only (4.1%)</text>',
            '    <text x="836" y="491" text-anchor="middle" font-size="8.2" font-style'
            '="italic" fill="#4A4742">Bacillus, Carnobacteriaceae</text>',
            "",
            '    <text x="836" y="224" text-anchor="middle" font-size="15.5" font-weig'
            'ht="700" fill="#181715">10</text>',
            '    <text x="836" y="238" text-anchor="middle" font-size="9.0" font-weigh'
            't="700" fill="#181715">PD ∩ CRC (20.4%)</text>',
            '    <text x="836" y="250" text-anchor="middle" font-size="7.8" font-style'
            '="italic" fill="#4A4742">R. bromii, R. flavefaciens,</text>',
            '    <text x="836" y="260" text-anchor="middle" font-size="7.8" font-style'
            '="italic" fill="#4A4742">B. cellulosilyticus, A. hadrus,</text>',
            '    <text x="836" y="270" text-anchor="middle" font-size="7.8" font-style'
            '="italic" fill="#4A4742">R. faecis, E. siraeum +4</text>',
            "",
            '    <text x="762" y="368" text-anchor="middle" font-size="14.5" font-weig'
            'ht="700" fill="#181715">4</text>',
            '    <text x="762" y="381" text-anchor="middle" font-size="8.6" font-weigh'
            't="700" fill="#8C5A2B">PD ∩ HIV (8.2%)</text>',
            '    <text x="762" y="393" text-anchor="middle" font-size="7.6" font-style'
            '="italic" fill="#4A4742">M. torques, Bacteroidales,</text>',
            '    <text x="762" y="403" text-anchor="middle" font-size="7.6" font-style'
            '="italic" fill="#4A4742">Lachnoclostridium, Bacteroidaceae</text>',
            "",
            '    <text x="910" y="368" text-anchor="middle" font-size="14.5" font-weig'
            'ht="700" fill="#181715">1</text>',
            '    <text x="910" y="381" text-anchor="middle" font-size="8.6" font-weigh'
            't="700" fill="#2B4C6F">HIV ∩ CRC (2.0%)</text>',
            '    <text x="910" y="394" text-anchor="middle" font-size="8.0" font-style'
            '="italic" fill="#4A4742">Actinomycetota</text>',
            "",
            '    <circle cx="836" cy="336" r="36" fill="#FCFBF7" stroke="#9E2A2B" st'
            'roke-width="2.2"/>',
            '    <circle cx="836" cy="336" r="32.5" fill="none" stroke="#9E2A2B" str'
            'oke-width="0.6" stroke-opacity="0.5"/>',
            '    <text x="836" y="328" text-anchor="middle" font-size="18.5" font-weig'
            'ht="700" fill="#9E2A2B">17</text>',
            '    <text x="836" y="343" text-anchor="middle" font-size="9.2" font-weigh'
            't="700" fill="#181715">3-Way Core</text>',
            '    <text x="836" y="355" text-anchor="middle" font-size="8.6" font-weigh'
            't="700" fill="#9E2A2B">(34.7%)</text>',
            "",
            '    <line x1="872" y1="322" x2="1040" y2="235" stroke="#9E2A2B" stroke-'
            'width="0.85" stroke-dasharray="2.5,2.5" stroke-opacity="0.7"/>',
            "",
            '    <g id="card_3way">',
            '      <rect x="1040" y="158" width="328" height="220" fill="#F4F1EA" st'
            'roke="#D6D1C7" stroke-width="0.9"/>',
            '      <rect x="1040" y="158" width="3.5" height="220" fill="#9E2A2B"/>',
            '      <text x="1053" y="176" font-size="9.6" font-weight="700" fill="#9E'
            '2A2B" letter-spacing="0.5">1. SHARED BY ALL 3 DISEASES (3-WAY CORE)</text>',
            '      <text x="1053" y="197" font-size="14.5" font-weight="700" fill="#1'
            '81715">34.7% (17 / 49) <tspan font-size="10.2" font-weight="400" fill="#4A4'
            '742">vs. 9.2% baseline (49 / 530)</tspan></text>',
            '      <text x="1053" y="214" font-size="10.4" font-weight="700" fill="#9'
            'E2A2B">3.75\u00d7 Enriched  (OR = 7.45, Fisher p = 1.29 \u00d7 10⁻⁷)</text>',
            '      <text x="1053" y="229" font-size="9.2" fill="#4A4742">10.6\u00d7 vs'
            ". 530-taxon null (p = 6.2\u00d710⁻¹⁶) · 18.6\u00d7 vs. U₁ (p = 1.5\u00d710⁻²²)</t"
            "ext>",
            "",
            '      <line x1="1053" y1="238" x2="1356" y2="238" stroke="#D6D1C7" stro'
            'ke-width="0.8"/>',
            '      <text x="1053" y="253" font-size="9.6" font-weight="700" fill="#18'
            '1715">All 17 Core Taxa (9 Genera [17.1\u00d7, p=3.5\u00d710⁻¹²] + 8 Higher):</te'
            "xt>",
            '      <text x="1053" y="268" font-size="9.1" fill="#181715">• <tspan font'
            '-style="italic" font-weight="700">Roseburia</tspan> (25 studies: 18·2·5), <ts'
            'pan font-style="italic" font-weight="700">Blautia</tspan> (19), <tspan font-s'
            'tyle="italic" font-weight="700">Clostridium</tspan> (16)</text>',
            '      <text x="1053" y="282" font-size="9.1" fill="#181715">• <tspan font'
            '-style="italic" font-weight="700">Bacteroides</tspan> (15), <tspan font-style'
            '="italic" font-weight="700">Lachnospira</tspan> (14), <tspan font-style="ita'
            'lic" font-weight="700">Ruminococcus</tspan> (14)</text>',
            '      <text x="1053" y="296" font-size="9.1" fill="#181715">• <tspan font'
            '-style="italic" font-weight="700">Anaerostipes</tspan> (8), <tspan font-style'
            '="italic" font-weight="700">Butyrivibrio</tspan> (4), <tspan font-style="ita'
            'lic" font-weight="700">Eubacterium</tspan> (4)</text>',
            '      <text x="1053" y="312" font-size="8.9" fill="#4A4742"><tspan font-w'
            'eight="700" fill="#181715">Higher ranks (8):</tspan> <tspan font-style="ital'
            'ic">Lachnospiraceae</tspan> (26), <tspan font-style="italic">Oscillospiraceae<'
            "/tspan> (22),</text>",
            '      <text x="1053" y="325" font-size="8.9" fill="#4A4742"><tspan font-s'
            'tyle="italic">Bacillota</tspan> (18), <tspan font-style="italic">Clostridia</'
            'tspan> (13), <tspan font-style="italic">Bacteroidota</tspan> (13), <tspan font-'
            'style="italic">Eubacteriales</tspan> (11),</text>',
            '      <text x="1053" y="338" font-size="8.9" fill="#4A4742"><tspan font-s'
            'tyle="italic">Clostridiaceae</tspan> (10), <tspan font-style="italic">Bactero'
            "idia</tspan> (7)</text>",
            '      <text x="1053" y="356" font-size="8.8" font-style="italic" fill="#'
            '2B4C6F">Genus+Species only (N=93): PD\u2013CRC J = 62.5% (p = 6.0\u00d710⁻¹⁰)</t'
            "ext>",
            '      <text x="1053" y="369" font-size="8.8" font-style="italic" fill="#'
            '2B4C6F">Confirms convergence is driven by fine-rank butyrate/fiber genera.</text'
            ">",
            "    </g>",
            "",
            '    <g id="card_2way">',
            '      <rect x="1040" y="390" width="328" height="156" fill="#F4F1EA" st'
            'roke="#D6D1C7" stroke-width="0.9"/>',
            '      <rect x="1040" y="390" width="3.5" height="156" fill="#2B4C6F"/>',
            '      <text x="1053" y="408" font-size="9.6" font-weight="700" fill="#2B'
            '4C6F" letter-spacing="0.5">2. SHARED BY ≥ 2 DISEASES (CONVERGENT GUILD)</text>',
            '      <text x="1053" y="429" font-size="14.5" font-weight="700" fill="#1'
            '81715">65.3% (32 / 49) <tspan font-size="10.2" font-weight="400" fill="#4A4'
            '742">vs. 27.2% baseline (144 / 530)</tspan></text>',
            '      <text x="1053" y="446" font-size="10.4" font-weight="700" fill="#2'
            'B4C6F">2.40\u00d7 Enriched  (OR = 6.20, Fisher p = 4.53 \u00d7 10⁻⁹)</text>',
            '      <text x="1053" y="461" font-size="9.2" fill="#4A4742">Only 17 / 49 '
            "taxa (34.7%) are unique to a single disease (vs. 72.8% bg).</text>",
            "",
            '      <line x1="1053" y1="470" x2="1356" y2="470" stroke="#D6D1C7" stro'
            'ke-width="0.8"/>',
            '      <text x="1053" y="485" font-size="9.5" font-weight="700" fill="#18'
            '1715">Set Identity: Full Pairwise Overlap = Core (17) + 2-Way Slice</text>',
            '      <text x="1053" y="500" font-size="9.0" fill="#181715">• <tspan font'
            '-weight="700">PD ∩ CRC:</tspan> 17 core + 10 PD∩CRC only = <tspan font-weight='
            '"700" fill="#9E2A2B">27 shared (57.4%)</tspan></text>',
            '      <text x="1053" y="514" font-size="9.0" fill="#181715">• <tspan font'
            '-weight="700">PD ∩ HIV:</tspan> 17 core + 4 PD∩HIV only = <tspan font-weight="'
            '700" fill="#8C5A2B">21 shared (48.8%)</tspan></text>',
            '      <text x="1053" y="528" font-size="9.0" fill="#181715">• <tspan font'
            '-weight="700">HIV ∩ CRC:</tspan> 17 core + 1 HIV∩CRC only = <tspan font-weight='
            '"700" fill="#2B4C6F">18 shared (45.0%)</tspan></text>',
            '      <text x="1053" y="541" font-size="8.6" font-style="italic" fill="#'
            '4A4742">Because 53% of shared taxa are lost in all 3 diseases, 17 &gt; 15 (10+4+'
            "1).</text>",
            "    </g>",
            "",
            '    <line x1="626" y1="556" x2="1368" y2="556" stroke="#D6D1C7" stroke-width="0.8"/>',
            '    <text x="626" y="575" class="sub-header">49-Taxon Union Partition vs. A'
            "ll 530 Stool Decreased Taxa (Proportional Share of Union)</text>",
            "",
            '    <line x1="805" y1="591" x2="1164.2" y2="591" stroke="#9E2A2B" strok'
            'e-width="1.0"/>',
            '    <line x1="805" y1="591" x2="805" y2="595" stroke="#9E2A2B" stroke-width="1.0"/>',
            '    <line x1="1164.2" y1="591" x2="1164.2" y2="595" stroke="#9E2A2B" st'
            'roke-width="1.0"/>',
            '    <text x="984.6" y="587" text-anchor="middle" font-size="9.4" font-wei'
            'ght="700" fill="#9E2A2B">65.3% shared across ≥2 diseases (32 / 49 taxa; Fishe'
            "r p = 4.53 \u00d7 10⁻⁹)</text>",
            "",
            '    <text x="626" y="614" font-size="10.8" font-weight="700" fill="#1817'
            '15">Cellulose Guild (N = 49)</text>',
            '    <rect x="805" y="598" width="190.8" height="24" fill="#9E2A2B"/>',
            '    <text x="900.4" y="614" text-anchor="middle" font-size="10" font-weig'
            'ht="700" fill="#FFFFFF">3-Way Core: 17 (34.7%)</text>',
            '    <rect x="995.8" y="598" width="168.4" height="24" fill="#8C5A2B"/>',
            '    <text x="1080.0" y="614" text-anchor="middle" font-size="10" font-wei'
            'ght="700" fill="#FFFFFF">2-Way Only: 15 (30.6%)</text>',
            '    <rect x="1164.2" y="598" width="190.8" height="24" fill="#E5E0D5" s'
            'troke="#B8B2A6" stroke-width="0.8"/>',
            '    <text x="1259.6" y="614" text-anchor="middle" font-size="10" font-wei'
            'ght="700" fill="#181715">1-Disease Only: 17 (34.7%)</text>',
            "",
            '    <text x="626" y="646" font-size="10.5" fill="#4A4742">Stool Backgroun'
            "d (N = 530)</text>",
            '    <rect x="805" y="631" width="50.8" height="21" fill="#9E2A2B" fill-'
            'opacity="0.58"/>',
            '    <text x="830.4" y="645" text-anchor="middle" font-size="9.2" font-wei'
            'ght="700" fill="#FFFFFF">9.2%</text>',
            '    <rect x="855.8" y="631" width="98.6" height="21" fill="#8C5A2B" fil'
            'l-opacity="0.58"/>',
            '    <text x="905.1" y="645" text-anchor="middle" font-size="9.2" font-wei'
            'ght="700" fill="#FFFFFF">2-Way: 17.9%</text>',
            '    <rect x="954.4" y="631" width="400.6" height="21" fill="#EFECE4" st'
            'roke="#D6D1C7" stroke-width="0.8"/>',
            '    <text x="1154.7" y="645" text-anchor="middle" font-size="9.5" fill="'
            '#4A4742">1-Disease Only: 386 / 530 (72.8% restricted to a single disease)</text>',
            "",
            '    <text x="626" y="676" font-size="9.8" fill="#181715"><tspan font-weig'
            'ht="700">Union partition (49 taxa):</tspan> 17 (all 3) + 15 (2 only: 10 PD∩CRC '
            "+ 4 PD∩HIV + 1 HIV∩CRC) + 17 (1 only: 9 PD + 6 CRC + 2 HIV) = 49 (100%).</text>",
            '    <text x="626" y="693" font-size="9.5" font-style="italic" fill="#4A4'
            '742">Disease totals:  |PD| = 17 + 10 + 4 + 9 = 40 taxa   ·   |CRC| = 17 + 10 + 1'
            " + 6 = 34 taxa   ·   |HIV| = 17 + 4 + 1 + 2 = 24 taxa.</text>",
            "  </g>",
            "",
            '  <line x1="32" y1="724" x2="1368" y2="724" stroke="#D6D1C7" stroke-width="0.85"/>',
            '  <text x="32" y="745" class="footnote">Source: BugSigDB human stool signat'
            "ures (90 studies: 46 PD, 17 HIV, 27 CRC) \u00d7 KG-Microbe cellulose catabolism ("
            "METPO:2000103 | GO:0030245). Reflects taxon-annotation alignment, not measured in"
            " vivo cellulolytic flux.</text>",
            '  <text x="1368" y="745" text-anchor="end" class="footnote">Plate VII (Pa'
            "ge 8 Venn Redesign)  ·  data/analysis/cellulose_decr</text>",
            "</svg>",
            "",
        ]
    )


def render_fig8_svg() -> str:
    """Return the standalone SVG for Slide 15 / fig8 with wrapped lines inside boxes (#550)."""
    return "\n".join(
        [
            '<svg viewBox="0 0 1500 740" width="1500" height="740" xmlns="http://www.w3'
            ".org/2000/svg\" style=\"font-family:'Palatino','Georgia',serif\">",
            '      <rect x="0" y="0" width="1500" height="740" fill="#FCFBF7"/>',
            "",
            "      <!-- Top Banner -->",
            '      <rect x="15" y="10" width="1470" height="124" rx="8" fill="#F4F1'
            'EA" stroke="#D8D2C5" stroke-width="1.2"/>',
            '      <text x="35" y="34" font-size="15.5" font-weight="bold" fill="#181'
            '715">OPPOSITE STOOL TRAIT DIRECTIONS IN HUMAN SIGNATURES (BugSigDB \u00d7 KG-Mic'
            "robe)</text>",
            '      <text x="35" y="55" font-size="12.5" font-style="italic" fill="#4A'
            '4843">Interpretation note: BugSigDB records reported taxon abundance shifts; KG-'
            "Microbe maps taxa to curated trait annotations (trait-annotation alignment, not d"
            "irect in vivo metabolite flux).</text>",
            "",
            "      <!-- Left Box in Banner (width=680, x=30..710) -->",
            '      <rect x="30" y="68" width="680" height="56" rx="5" fill="#F9ECEC'
            '" stroke="#9E2A2B" stroke-width="1.2"/>',
            '      <circle cx="52" cy="96" r="8.5" fill="#9E2A2B"/>',
            '      <text x="70" y="85" font-size="12.0" font-weight="bold" fill="#9E2'
            'A2B">PAIRS 1\u20133 · DECREASE IN CELLULOSE CATABOLISM (METPO:2000103 | GO:00302'
            "45)</text>",
            '      <text x="70" y="101" font-size="11.0" fill="#181715">Reduced repres'
            "entation of anaerobic fiber-fermenting Clostridia and Bacteroidota</text>",
            '      <text x="70" y="116" font-size="11.0" fill="#181715">across human s'
            "tool cohorts in PD (\u221227.4 pp), HIV (\u221253.7 pp), and CRC (\u221236.9 pp)<"
            "/text>",
            "",
            "      <!-- Arrow between boxes -->",
            '      <path d="M 725 96 L 770 96" stroke="#181715" stroke-width="2.2" marke'
            'r-end="url(#arrow)"/>',
            "      <defs>",
            '        <marker id="arrow" viewBox="0 0 10 10" refX="6" refY="5" markerWi'
            'dth="6" markerHeight="6" orient="auto-start-reverse">',
            '          <path d="M 0 1 L 10 5 L 0 9 z" fill="#181715"/>',
            "        </marker>",
            "      </defs>",
            "",
            "      <!-- Right Box in Banner (width=680, x=785..1465) -->",
            '      <rect x="785" y="68" width="680" height="56" rx="5" fill="#EBF1F'
            '5" stroke="#2B4C65" stroke-width="1.2"/>',
            '      <circle cx="807" cy="96" r="8.5" fill="#2B4C65"/>',
            '      <text x="825" y="85" font-size="12.0" font-weight="bold" fill="#2B'
            '4C65">PAIR 4 · INCREASE IN HALOTOLERANCE IN PARKINSON\u2019S DISEASE (METPO:1000'
            "623)</text>",
            '      <text x="825" y="101" font-size="11.0" fill="#181715">Increased rep'
            "resentation of osmotolerant and stress-tolerant taxa (+24.4 pp, BH q = 0.0265;</t"
            "ext>",
            '      <text x="825" y="116" font-size="11.0" fill="#181715">paired within'
            "-study Wilcoxon 33 Up vs. 8 Down, p = 1.67 \u00d7 10⁻⁴)</text>",
            "",
            "      <!-- Four Detailed Cards (2x2 Grid, each width=725) -->",
            "      <!-- Card 1: PD Cellulose Decrease -->",
            '      <g transform="translate(15, 148)">',
            '        <rect x="0" y="0" width="725" height="275" rx="8" fill="#FFFFF'
            'F" stroke="#9E2A2B" stroke-width="1.5"/>',
            '        <rect x="0" y="0" width="725" height="38" rx="8" fill="#9E2A2B"/>',
            '        <text x="14" y="24" font-size="12.2" font-weight="bold" fill="#F'
            'FFFFF">1. Parkinson\u2019s Disease (MONDO:0005180) \u00d7 Cellulose Catabolism D'
            "ecreased</text>",
            '        <text x="711" y="24" font-size="11.8" font-weight="bold" fill="#'
            'F9ECEC" text-anchor="end">Δ = \u221227.4 pp | BH q = 0.0077</text>',
            "",
            '        <text x="16" y="56" font-size="12.0" font-weight="bold" fill="#1'
            '81715">MTP Stool Signal (46 Decreased Studies / 52 Total, 83 Paired Signature Ex'
            "periments):</text>",
            '        <text x="16" y="74" font-size="11.2" fill="#33312E">• Present in '
            "27.5% of Increased (36/131) vs. 54.8% of Decreased profiles (68/124; Fisher p = 1"
            ".28 \u00d7 10⁻⁵)</text>",
            '        <text x="16" y="91" font-size="11.2" fill="#33312E">• Within-stud'
            "y experiment sign test: 57 Decreased vs. 26 Increased (68.7% decrease, binomial p"
            " = 8.78 \u00d7 10⁻⁴)</text>",
            '        <text x="16" y="108" font-size="11.2" fill="#33312E">• Top deplet'
            "ed core genera: Roseburia (18 studies), Blautia (12), Clostridium (7),</text>",
            '        <text x="26" y="123" font-size="11.2" fill="#33312E">Bacteroides '
            "(6), Ruminococcus (6), and Lachnospira (5)</text>",
            "",
            '        <line x1="16" y1="134" x2="709" y2="134" stroke="#E5E0D5" stroke-width="1"/>',
            '        <text x="16" y="151" font-size="12.0" font-weight="bold" fill="#'
            '9E2A2B">Published Literature Citations:</text>',
            '        <text x="16" y="169" font-size="11.2" fill="#181715">• <tspan fon'
            't-weight="bold">Wallen et al. (2022) Nat Commun (PMID 36414626, N=724):</tspan>'
            " Reproducible loss of Roseburia,</text>",
            '        <text x="26" y="184" font-size="11.0" fill="#4A4843">Faecalibacte'
            "rium, Blautia and plant-polysaccharide degradation pathways in PD metagenomes.</t"
            "ext>",
            '        <text x="16" y="203" font-size="11.2" fill="#181715">• <tspan fon'
            't-weight="bold">Romano et al. (2021) npj Parkinsons Dis (PMID 33654098):</tspan'
            "> 21-cohort meta-analysis</text>",
            '        <text x="26" y="218" font-size="11.0" fill="#4A4843">confirming c'
            "onsistent depletion of the colonic SCFA and dietary fiber-fermenting guild.</text"
            ">",
            '        <text x="16" y="237" font-size="11.2" fill="#181715">• <tspan fon'
            't-weight="bold">Nishiwaki et al. (2020) Mov Disord (PMID 32472608) &amp; Aho et'
            " al. (2019):</tspan> Depletion of</text>",
            '        <text x="26" y="252" font-size="11.0" fill="#4A4843">cellulose-fe'
            "rmenting Clostridia correlates with lower fecal SCFA levels and PD progression.</"
            "text>",
            "      </g>",
            "",
            "      <!-- Card 2: HIV Cellulose Decrease -->",
            '      <g transform="translate(760, 148)">',
            '        <rect x="0" y="0" width="725" height="275" rx="8" fill="#FFFFF'
            'F" stroke="#9E2A2B" stroke-width="1.5"/>',
            '        <rect x="0" y="0" width="725" height="38" rx="8" fill="#9E2A2B"/>',
            '        <text x="14" y="24" font-size="12.2" font-weight="bold" fill="#F'
            'FFFFF">2. HIV Infection (MONDO:0005109) \u00d7 Cellulose Catabolism Decreased</t'
            "ext>",
            '        <text x="711" y="24" font-size="11.8" font-weight="bold" fill="#'
            'F9ECEC" text-anchor="end">Δ = \u221253.7 pp | BH q = 0.0080</text>',
            "",
            '        <text x="16" y="56" font-size="12.0" font-weight="bold" fill="#1'
            '81715">MTP Stool Signal (17 Decreased Studies / 18 Total, 25 Paired Signature Ex'
            "periments):</text>",
            '        <text x="16" y="74" font-size="11.2" fill="#33312E">• Present in '
            "19.2% of Increased (10/52) vs. 72.9% of Decreased profiles (35/48; Fisher p = 6.2"
            "3 \u00d7 10⁻⁸)</text>",
            '        <text x="16" y="91" font-size="11.2" fill="#33312E">• Within-stud'
            "y experiment sign test: 19 Decreased vs. 6 Increased (76.0% decrease, binomial p "
            "= 0.0146)</text>",
            '        <text x="16" y="108" font-size="11.2" fill="#33312E">• Top deplet'
            "ed core genera: Bacteroides (7 studies), Clostridium (4), Ruminococcus (4),</text"
            ">",
            '        <text x="26" y="123" font-size="11.2" fill="#33312E">Lachnospira '
            "(3), Anaerostipes (3), and Roseburia (2)</text>",
            "",
            '        <line x1="16" y1="134" x2="709" y2="134" stroke="#E5E0D5" stroke-width="1"/>',
            '        <text x="16" y="151" font-size="12.0" font-weight="bold" fill="#'
            '9E2A2B">Published Literature Citations:</text>',
            '        <text x="16" y="169" font-size="11.2" fill="#181715">• <tspan fon'
            't-weight="bold">Vujkovic-Cvijin et al. (2013) Sci Transl Med (PMID 23843452):</'
            "tspan> Depletion of Bacteroides</text>",
            '        <text x="26" y="184" font-size="11.0" fill="#4A4843">and Clostrid'
            "ia fiber-fermenters tracks mucosal Th17 loss and immune activation in HIV-1.</tex"
            "t>",
            '        <text x="16" y="203" font-size="11.2" fill="#181715">• <tspan fon'
            't-weight="bold">Noguera-Julian et al. (2016) EBioMedicine (PMID 27077129):</tsp'
            "an> After stratifying for sexual</text>",
            '        <text x="26" y="218" font-size="11.0" fill="#4A4843">behavior, HI'
            "V status tracks loss of cellulolytic Ruminococcaceae and Lachnospiraceae.</text>",
            '        <text x="16" y="237" font-size="11.2" fill="#181715">• <tspan fon'
            't-weight="bold">Dillon et al. (2014) Mucosal Immunol &amp; Lozupone et al. (201'
            "3):</tspan> Consistent depletion</text>",
            '        <text x="26" y="252" font-size="11.0" fill="#4A4843">of Roseburia'
            ", Coprococcus, and commensal fiber-degrading taxa in HIV stool cohorts.</text>",
            "      </g>",
            "",
            "      <!-- Card 3: CRC Cellulose Decrease -->",
            '      <g transform="translate(15, 442)">',
            '        <rect x="0" y="0" width="725" height="275" rx="8" fill="#FFFFF'
            'F" stroke="#9E2A2B" stroke-width="1.5"/>',
            '        <rect x="0" y="0" width="725" height="38" rx="8" fill="#9E2A2B"/>',
            '        <text x="14" y="24" font-size="12.2" font-weight="bold" fill="#F'
            'FFFFF">3. Colorectal Cancer (MONDO:0005575) \u00d7 Cellulose Catabolism Decrease'
            "d</text>",
            '        <text x="711" y="24" font-size="11.8" font-weight="bold" fill="#'
            'F9ECEC" text-anchor="end">Δ = \u221236.9 pp | BH q = 0.0336</text>',
            "",
            '        <text x="16" y="56" font-size="12.0" font-weight="bold" fill="#1'
            '81715">MTP Stool Signal (27 Decreased Studies / 30 Total, 37 Paired Signature Ex'
            "periments):</text>",
            '        <text x="16" y="74" font-size="11.2" fill="#33312E">• Present in '
            "26.9% of Increased (18/67) vs. 63.8% of Decreased profiles (37/58; Fisher p = 3.7"
            "4 \u00d7 10⁻⁵)</text>",
            '        <text x="16" y="91" font-size="11.2" fill="#33312E">• Within-stud'
            "y experiment sign test: 30 Decreased vs. 7 Increased (81.1% decrease, binomial p "
            "= 1.91 \u00d7 10⁻⁴)</text>",
            '        <text x="16" y="108" font-size="11.2" fill="#33312E">• Top deplet'
            "ed core genera: Lachnospira (6 studies), Roseburia (5), Blautia (5),</text>",
            '        <text x="26" y="123" font-size="11.2" fill="#33312E">Clostridium '
            "(5), Ruminococcus (4), and Anaerostipes (3)</text>",
            "",
            '        <line x1="16" y1="134" x2="709" y2="134" stroke="#E5E0D5" stroke-width="1"/>',
            '        <text x="16" y="151" font-size="12.0" font-weight="bold" fill="#'
            '9E2A2B">Published Literature Citations:</text>',
            '        <text x="16" y="169" font-size="11.2" fill="#181715">• <tspan fon'
            't-weight="bold">Wirbel et al. (2019) Nat Med (PMID 30936547, 8-cohort meta-anal'
            "ysis):</tspan> Depletion of</text>",
            '        <text x="26" y="184" font-size="11.0" fill="#4A4843">carbohydrate'
            "/fiber degradation CAZymes (Roseburia, Eubacterium, Ruminococcus) in CRC stool.</"
            "text>",
            '        <text x="16" y="203" font-size="11.2" fill="#181715">• <tspan fon'
            't-weight="bold">Thomas et al. (2019) Nat Med (PMID 30936548, N=969):</tspan> Cr'
            "oss-cohort loss of</text>",
            '        <text x="26" y="218" font-size="11.0" fill="#4A4843">dietary fibe'
            "r-fermenting Clostridia across independent metagenomic CRC studies.</text>",
            '        <text x="16" y="237" font-size="11.2" fill="#181715">• <tspan fon'
            't-weight="bold">Yu et al. (2017) Gut &amp; O\u2019Keefe (2016) Nat Rev Gastroen'
            "terol Hepatol:</tspan> Lower abundance</text>",
            '        <text x="26" y="252" font-size="11.0" fill="#4A4843">of plant-fib'
            "er fermenting butyrate producers in colorectal adenoma and carcinoma stool.</text"
            ">",
            "      </g>",
            "",
            "      <!-- Card 4: PD Halotolerance Increase -->",
            '      <g transform="translate(760, 442)">',
            '        <rect x="0" y="0" width="725" height="275" rx="8" fill="#FFFFF'
            'F" stroke="#2B4C65" stroke-width="1.5"/>',
            '        <rect x="0" y="0" width="725" height="38" rx="8" fill="#2B4C65"/>',
            '        <text x="14" y="24" font-size="12.2" font-weight="bold" fill="#F'
            'FFFFF">4. Parkinson\u2019s Disease (MONDO:0005180) \u00d7 Halotolerance Increase'
            "d</text>",
            '        <text x="711" y="24" font-size="11.8" font-weight="bold" fill="#'
            'EBF1F5" text-anchor="end">Δ = +24.4 pp | BH q = 0.0265</text>',
            "",
            '        <text x="16" y="56" font-size="12.0" font-weight="bold" fill="#1'
            '81715">MTP Stool Signal (43 Increased Studies / 52 Total, Paired Within-Study Wi'
            "lcoxon Test):</text>",
            '        <text x="16" y="74" font-size="11.2" fill="#33312E">• Present in '
            "52.7% of Increased (69/131) vs. 28.2% of Decreased profiles (35/124; Fisher p = 7"
            ".76 \u00d7 10⁻⁵)</text>",
            '        <text x="16" y="91" font-size="11.2" fill="#33312E">• Paired with'
            "in-study fraction difference: 33 Studies Up vs. 8 Down (4 tied; Wilcoxon p = 1.67"
            " \u00d7 10⁻⁴)</text>",
            '        <text x="16" y="108" font-size="11.2" fill="#33312E">• Top increa'
            "sed osmotolerant/stress-tolerant taxa: Akkermansia, Bifidobacterium,</text>",
            '        <text x="26" y="123" font-size="11.2" fill="#33312E">Lactobacillu'
            "s, Christensenella, Enterococcus, and Corynebacterium</text>",
            "",
            '        <line x1="16" y1="134" x2="709" y2="134" stroke="#E5E0D5" stroke-width="1"/>',
            '        <text x="16" y="151" font-size="12.0" font-weight="bold" fill="#'
            '2B4C65">Published Literature Citations:</text>',
            '        <text x="16" y="169" font-size="11.2" fill="#181715">• <tspan fon'
            't-weight="bold">Colonic Transit &amp; Osmolality (Vandeputte et al. 2016 Gut; R'
            "oager et al. 2016):</tspan> Prolonged</text>",
            '        <text x="26" y="184" font-size="11.0" fill="#4A4843">colonic tran'
            "sit in PD increases water reabsorption, stool firmness, and luminal osmolality.</"
            "text>",
            '        <text x="16" y="203" font-size="11.2" fill="#181715">• <tspan fon'
            't-weight="bold">Bedarf et al. (2017) Genome Med &amp; Boktor et al. (2023) Mov '
            "Disord:</tspan> Loss of fast-growing</text>",
            '        <text x="26" y="218" font-size="11.0" fill="#4A4843">saccharolyti'
            "c Clostridia co-occurs with slow-growing stress-tolerant bacteria in PD stool.</t"
            "ext>",
            '        <text x="16" y="237" font-size="11.2" fill="#181715">• <tspan fon'
            't-weight="bold">Reference table:</tspan> Saved in</text>',
            '        <text x="26" y="252" font-family="monospace" font-size="10.8" fil'
            'l="#4A4843">data/analysis/cellulose_decr/FOUR_HIGHLIGHTED_PAIRS_LITERATURE_SUPP'
            "ORT.md</text>",
            "      </g>",
            "    </svg>",
        ]
    )


def render_fig9_svg() -> str:
    """Return the standalone SVG string for Slide 16 / fig9
    (fig9_set_theory_17_core_vs_pairwise.svg).
    """
    return "\n".join(
        [
            '<svg viewBox="0 0 1500 740" width="1500" height="740" xmlns="http://www.w3'
            ".org/2000/svg\" style=\"font-family:'Palatino','Georgia',serif\">",
            '      <rect x="0" y="0" width="1500" height="740" fill="#FCFBF7"/>',
            "",
            "      <!-- Left Panel: Exact 7-Region Venn Partition -->",
            '      <rect x="15" y="15" width="710" height="710" rx="10" fill="#FFFF'
            'FF" stroke="#D8D2C5" stroke-width="1.4"/>',
            '      <text x="40" y="52" font-size="18" font-weight="bold" fill="#18171'
            '5">A. The 7 Mutually Exclusive Regions (Sum = 49 Taxa)</text>',
            '      <text x="40" y="76" font-size="13.5" font-style="italic" fill="#4A'
            '4843">When a Venn diagram partitions the 49 union taxa into non-overlapping subs'
            "ets:</text>",
            "",
            '      <circle cx="290" cy="295" r="145" fill="#9E2A2B" fill-opacity="0.1'
            '6" stroke="#9E2A2B" stroke-width="2.5"/>',
            '      <text x="175" y="140" font-size="16" font-weight="bold" fill="#9E2'
            'A2B" text-anchor="middle">Parkinson\u2019s Disease</text>',
            '      <text x="175" y="160" font-size="13" fill="#4A4843" text-anchor="m'
            'iddle">Total |PD| = 40 taxa</text>',
            "",
            '      <circle cx="450" cy="295" r="145" fill="#2B4C65" fill-opacity="0.1'
            '6" stroke="#2B4C65" stroke-width="2.5"/>',
            '      <text x="565" y="140" font-size="16" font-weight="bold" fill="#2B4'
            'C65" text-anchor="middle">Colorectal Cancer</text>',
            '      <text x="565" y="160" font-size="13" fill="#4A4843" text-anchor="m'
            'iddle">Total |CRC| = 34 taxa</text>',
            "",
            '      <circle cx="370" cy="430" r="145" fill="#8C5A2B" fill-opacity="0.1'
            '6" stroke="#8C5A2B" stroke-width="2.5"/>',
            '      <text x="370" y="605" font-size="16" font-weight="bold" fill="#8C5'
            'A2B" text-anchor="middle">HIV Infection</text>',
            '      <text x="370" y="625" font-size="13" fill="#4A4843" text-anchor="m'
            'iddle">Total |HIV| = 24 taxa</text>',
            "",
            '      <text x="215" y="275" font-size="22" font-weight="bold" fill="#181'
            '715" text-anchor="middle">9</text>',
            '      <text x="215" y="295" font-size="12" fill="#4A4843" text-anchor="m'
            'iddle">PD only</text>',
            "",
            '      <text x="525" y="275" font-size="22" font-weight="bold" fill="#181'
            '715" text-anchor="middle">6</text>',
            '      <text x="525" y="295" font-size="12" fill="#4A4843" text-anchor="m'
            'iddle">CRC only</text>',
            "",
            '      <text x="370" y="515" font-size="22" font-weight="bold" fill="#181'
            '715" text-anchor="middle">2</text>',
            '      <text x="370" y="535" font-size="12" fill="#4A4843" text-anchor="m'
            'iddle">HIV only</text>',
            "",
            '      <text x="370" y="230" font-size="22" font-weight="bold" fill="#8C5'
            'A2B" text-anchor="middle">10</text>',
            '      <text x="370" y="248" font-size="11.5" font-weight="bold" fill="#8'
            'C5A2B" text-anchor="middle">PD &amp; CRC only</text>',
            "",
            '      <text x="285" y="405" font-size="22" font-weight="bold" fill="#8C5'
            'A2B" text-anchor="middle">4</text>',
            '      <text x="285" y="423" font-size="11.5" font-weight="bold" fill="#8'
            'C5A2B" text-anchor="middle">PD &amp; HIV only</text>',
            "",
            '      <text x="455" y="405" font-size="22" font-weight="bold" fill="#8C5'
            'A2B" text-anchor="middle">1</text>',
            '      <text x="455" y="423" font-size="11.5" font-weight="bold" fill="#8'
            'C5A2B" text-anchor="middle">HIV &amp; CRC only</text>',
            "",
            '      <circle cx="370" cy="340" r="44" fill="#9E2A2B" stroke="#FFFFFF" '
            'stroke-width="2.5"/>',
            '      <text x="370" y="336" font-size="25" font-weight="bold" fill="#FFF'
            'FFF" text-anchor="middle">17</text>',
            '      <text x="370" y="355" font-size="11.5" font-weight="bold" fill="#F'
            '9ECEC" text-anchor="middle">ALL 3</text>',
            "",
            '      <rect x="40" y="645" width="660" height="58" rx="6" fill="#F4F1E'
            'A" stroke="#D8D2C5" stroke-width="1"/>',
            '      <text x="55" y="668" font-size="13" font-weight="bold" fill="#1817'
            '15">Why the 2-way-only slices (10, 4, 1) are smaller than 17:</text>',
            '      <text x="55" y="689" font-size="12.5" fill="#33312E">Each 2-way-onl'
            "y slice excludes the 17 core taxa. Total 2-way-only = 10 + 4 + 1 = 15 taxa.</text"
            ">",
            "",
            "      <!-- Right Panel: Full Pairwise Intersections vs. Mutually Exclusive Tiers -->",
            '      <rect x="750" y="15" width="735" height="710" rx="10" fill="#FFF'
            'FFF" stroke="#D8D2C5" stroke-width="1.4"/>',
            '      <text x="775" y="52" font-size="18" font-weight="bold" fill="#1817'
            '15">B. Full Pairwise Intersections Always Include the 17-Taxon Core</text>',
            '      <text x="775" y="76" font-size="13.5" font-style="italic" fill="#4'
            'A4843">Every complete 2-disease intersection |A ∩ B| equals the 3-way core (17) '
            "plus its 2-way-only slice:</text>",
            "",
            '      <g transform="translate(775, 115)">',
            '        <text x="0" y="20" font-size="15" font-weight="bold" fill="#1817'
            '15">1. Parkinson\u2019s Disease ∩ Colorectal Cancer  (|PD ∩ CRC| = 27 of 47 unio'
            "n, 57.4%)</text>",
            '        <rect x="0" y="34" width="340" height="42" fill="#9E2A2B" rx="4"/>',
            '        <text x="170" y="60" font-size="14" font-weight="bold" fill="#FF'
            'FFFF" text-anchor="middle">17 Core Taxa (in all 3 diseases)</text>',
            '        <rect x="344" y="34" width="200" height="42" fill="#8C5A2B" rx="4"/>',
            '        <text x="444" y="60" font-size="14" font-weight="bold" fill="#FF'
            'FFFF" text-anchor="middle">+ 10 (PD &amp; CRC only)</text>',
            '        <text x="560" y="60" font-size="16" font-weight="bold" fill="#9E'
            '2A2B">= 27 shared (&gt; 17)</text>',
            "      </g>",
            "",
            '      <g transform="translate(775, 225)">',
            '        <text x="0" y="20" font-size="15" font-weight="bold" fill="#1817'
            '15">2. Parkinson\u2019s Disease ∩ HIV Infection  (|PD ∩ HIV| = 21 of 43 union, 4'
            "8.8%)</text>",
            '        <rect x="0" y="34" width="340" height="42" fill="#9E2A2B" rx="4"/>',
            '        <text x="170" y="60" font-size="14" font-weight="bold" fill="#FF'
            'FFFF" text-anchor="middle">17 Core Taxa (in all 3 diseases)</text>',
            '        <rect x="344" y="34" width="80" height="42" fill="#8C5A2B" rx="4"/>',
            '        <text x="384" y="60" font-size="13.5" font-weight="bold" fill="#'
            'FFFFFF" text-anchor="middle">+ 4</text>',
            '        <text x="440" y="60" font-size="16" font-weight="bold" fill="#2B'
            '4C65">= 21 shared (&gt; 17)</text>',
            "      </g>",
            "",
            '      <g transform="translate(775, 335)">',
            '        <text x="0" y="20" font-size="15" font-weight="bold" fill="#1817'
            '15">3. HIV Infection ∩ Colorectal Cancer  (|HIV ∩ CRC| = 18 of 40 union, 45.0%)<'
            "/text>",
            '        <rect x="0" y="34" width="340" height="42" fill="#9E2A2B" rx="4"/>',
            '        <text x="170" y="60" font-size="14" font-weight="bold" fill="#FF'
            'FFFF" text-anchor="middle">17 Core Taxa (in all 3 diseases)</text>',
            '        <rect x="344" y="34" width="24" height="42" fill="#8C5A2B" rx="4"/>',
            '        <text x="356" y="60" font-size="13" font-weight="bold" fill="#FF'
            'FFFF" text-anchor="middle">+1</text>',
            '        <text x="385" y="60" font-size="16" font-weight="bold" fill="#8C'
            '5A2B">= 18 shared (&gt; 17)</text>',
            "      </g>",
            "",
            '      <rect x="775" y="455" width="685" height="248" rx="8" fill="#F9E'
            'CEC" stroke="#9E2A2B" stroke-width="1.4"/>',
            '      <text x="798" y="486" font-size="15.5" font-weight="bold" fill="#9'
            'E2A2B">Summary: Comparing the 3-Way Core (17) and the 2-Way-Only Subset (15)</te'
            "xt>",
            "",
            '      <text x="798" y="516" font-size="13.2" fill="#181715">1. <tspan fon'
            't-weight="bold">A 3-way intersection never exceeds a full 2-way intersection:</'
            "tspan></text>",
            '      <text x="814" y="538" font-size="12.8" fill="#33312E">Because (PD ∩'
            " HIV ∩ CRC) ⊆ (A ∩ B), the 3-way core (17) is smaller than every complete pairwis"
            "e</text>",
            '      <text x="814" y="558" font-size="12.8" fill="#33312E">overlap: 17 &'
            "lt; 27 (PD\u2013CRC), 17 &lt; 21 (PD\u2013HIV), and 17 &lt; 18 (HIV\u2013CRC).</t"
            "ext>",
            "",
            '      <text x="798" y="590" font-size="13.2" fill="#181715">2. <tspan fon'
            't-weight="bold">Why the mutually exclusive 2-disease-only subset is smaller (15'
            " vs. 17):</tspan></text>",
            '      <text x="814" y="612" font-size="12.8" fill="#33312E">Fiber-degradi'
            "ng taxa reported decreased in two diseases usually belong to the shared anaerobic"
            " guild</text>",
            '      <text x="814" y="632" font-size="12.8" fill="#33312E">(Roseburia, B'
            "lautia, Ruminococcus, Bacteroides, Clostridium) that is also reported decreased i"
            "n the third</text>",
            '      <text x="814" y="652" font-size="12.8" fill="#33312E">disease. As a'
            " result, 17 of the 32 shared taxa (53%) fall into the 3-disease intersection, whi"
            "le 15</text>",
            '      <text x="814" y="672" font-size="12.8" fill="#33312E">taxa (10 + 4 '
            "+ 1) appear in exactly two of the three diseases.</text>",
            "    </svg>",
        ]
    )


# --------------------------------------------------------------------------- 16-Slide HTML Deck


def render_slide_deck_html() -> str:
    """Return the complete 16-slide HTML presentation deck (#550)."""
    fig8_inline = render_fig8_svg().strip()
    fig9_inline = render_fig9_svg().strip()
    return "\n".join(
        [
            "<!DOCTYPE html>",
            '<html lang="en">',
            "<head>",
            '<meta charset="utf-8"/>',
            "<title>Cellulose Catabolism Convergence, Four Highlighted Pairs, and Three-Arm St"
            "atistics</title>",
            "<style>",
            "  @page {",
            "    size: 1600px 1000px;",
            "    margin: 0;",
            "  }",
            "  * {",
            "    box-sizing: border-box;",
            "    margin: 0;",
            "    padding: 0;",
            "  }",
            "  html, body {",
            "    width: 1600px;",
            "    background: #0F172A;",
            '    font-family: -apple-system, BlinkMacSystemFont, "Inter", "Segoe UI", Robo'
            "to, Helvetica, Arial, sans-serif;",
            "    color: #0F172A;",
            "    -webkit-print-color-adjust: exact;",
            "    print-color-adjust: exact;",
            "  }",
            "  .slide {",
            "    width: 1600px;",
            "    height: 1000px;",
            "    background: #F8FAFC;",
            "    page-break-after: always;",
            "    break-after: page;",
            "    display: flex;",
            "    flex-direction: column;",
            "    justify-content: space-between;",
            "    padding: 20px 28px 16px 28px;",
            "    overflow: hidden;",
            "    position: relative;",
            "  }",
            "  .cover-slide {",
            "    background: linear-gradient(135deg, #181715 0%, #2C2623 55%, #1E293B 100%);",
            "    color: #F8FAFC;",
            "    justify-content: center;",
            "    align-items: flex-start;",
            "    padding: 64px 88px;",
            "  }",
            "  .cover-badge {",
            "    display: inline-block;",
            "    background: rgba(158, 42, 43, 0.35);",
            "    border: 1px solid rgba(229, 115, 115, 0.55);",
            "    color: #FCA5A5;",
            "    font-size: 13px;",
            "    font-weight: 700;",
            "    letter-spacing: 0.08em;",
            "    text-transform: uppercase;",
            "    padding: 6px 14px;",
            "    border-radius: 999px;",
            "    margin-bottom: 20px;",
            "  }",
            "  .cover-slide h1 {",
            '    font-family: "Palatino", "Georgia", serif;',
            "    font-size: 42px;",
            "    font-weight: 700;",
            "    line-height: 1.15;",
            "    letter-spacing: -0.01em;",
            "    margin-bottom: 16px;",
            "    max-width: 1380px;",
            "    color: #FCFBF7;",
            "  }",
            "  .cover-slide p.lead {",
            "    font-size: 18.5px;",
            "    color: #CBD5E1;",
            "    max-width: 1300px;",
            "    line-height: 1.5;",
            "    margin-bottom: 34px;",
            "  }",
            "  .cover-grid {",
            "    display: grid;",
            "    grid-template-columns: repeat(4, 1fr);",
            "    gap: 18px;",
            "    width: 100%;",
            "    margin-bottom: 34px;",
            "  }",
            "  .cover-card {",
            "    background: rgba(252, 251, 247, 0.06);",
            "    border: 1px solid rgba(252, 251, 247, 0.14);",
            "    border-radius: 12px;",
            "    padding: 18px 20px;",
            "  }",
            "  .cover-card .num {",
            '    font-family: "Palatino", "Georgia", serif;',
            "    font-size: 32px;",
            "    font-weight: 700;",
            "    color: #F87171;",
            "    margin-bottom: 4px;",
            "  }",
            "  .cover-card .lbl {",
            "    font-size: 13.5px;",
            "    color: #E2E8F0;",
            "    line-height: 1.38;",
            "  }",
            "  .cover-toc {",
            "    display: grid;",
            "    grid-template-columns: repeat(2, 1fr);",
            "    gap: 14px 28px;",
            "    width: 100%;",
            "    background: rgba(15, 23, 42, 0.55);",
            "    border: 1px solid rgba(148, 163, 184, 0.22);",
            "    border-radius: 12px;",
            "    padding: 20px 26px;",
            "  }",
            "  .cover-toc-item {",
            "    font-size: 14px;",
            "    color: #E2E8F0;",
            "    line-height: 1.45;",
            "  }",
            "  .cover-toc-item b {",
            "    color: #FDBA74;",
            "  }",
            "  .slide-header {",
            "    display: flex;",
            "    align-items: center;",
            "    justify-content: space-between;",
            "    gap: 16px;",
            "    padding-bottom: 10px;",
            "    border-bottom: 1px solid #CBD5E1;",
            "    flex-shrink: 0;",
            "  }",
            "  .header-left {",
            "    display: flex;",
            "    align-items: center;",
            "    gap: 14px;",
            "    min-width: 0;",
            "  }",
            "  .cat-badge {",
            "    font-size: 11px;",
            "    font-weight: 700;",
            "    text-transform: uppercase;",
            "    letter-spacing: 0.06em;",
            "    padding: 5px 11px;",
            "    border-radius: 6px;",
            "    white-space: nowrap;",
            "    flex-shrink: 0;",
            "  }",
            "  .title-block {",
            "    min-width: 0;",
            "  }",
            "  .title-block h2 {",
            "    font-size: 19px;",
            "    font-weight: 750;",
            "    color: #0F172A;",
            "    line-height: 1.2;",
            "    white-space: nowrap;",
            "    overflow: hidden;",
            "    text-overflow: ellipsis;",
            "  }",
            "  .title-block .subtitle {",
            "    font-size: 12.5px;",
            "    color: #475569;",
            "    margin-top: 2px;",
            "    white-space: nowrap;",
            "    overflow: hidden;",
            "    text-overflow: ellipsis;",
            "  }",
            "  .header-right {",
            "    display: flex;",
            "    align-items: center;",
            "    gap: 8px;",
            "    flex-shrink: 0;",
            "  }",
            "  .meta-pill {",
            "    background: #F1F5F9;",
            "    border: 1px solid #CBD5E1;",
            "    color: #334155;",
            "    font-size: 11.5px;",
            "    padding: 4px 10px;",
            "    border-radius: 6px;",
            "    white-space: nowrap;",
            "  }",
            "  .slide-num {",
            "    font-size: 12px;",
            "    font-weight: 700;",
            "    color: #475569;",
            "    padding-left: 6px;",
            "    white-space: nowrap;",
            "  }",
            "  .fig-canvas {",
            "    flex: 1 1 auto;",
            "    min-height: 0;",
            "    margin-top: 8px;",
            "    margin-bottom: 8px;",
            "    background: #FFFFFF;",
            "    border: 1px solid #E2E8F0;",
            "    border-radius: 10px;",
            "    padding: 6px;",
            "    display: flex;",
            "    align-items: center;",
            "    justify-content: center;",
            "    overflow: hidden;",
            "    box-shadow: 0 1px 3px rgba(15, 23, 42, 0.04);",
            "  }",
            "  .fig-canvas img, .fig-canvas svg {",
            "    width: 100%;",
            "    height: 100%;",
            "    object-fit: contain;",
            "    display: block;",
            "  }",
            "  .slide-footer {",
            "    flex-shrink: 0;",
            "    background: #F1F5F9;",
            "    border: 1px solid #CBD5E1;",
            "    border-radius: 8px;",
            "    padding: 8px 14px;",
            "    font-size: 12px;",
            "    line-height: 1.38;",
            "    color: #1E293B;",
            "    display: grid;",
            "    grid-template-columns: 1fr 1fr;",
            "    gap: 16px;",
            "  }",
            "  .slide-footer code {",
            "    font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace;",
            "    font-size: 11px;",
            "    background: #E2E8F0;",
            "    padding: 1px 5px;",
            "    border-radius: 4px;",
            "  }",
            "</style>",
            "</head>",
            "<body>",
            "",
            '<section class="slide cover-slide">',
            '  <span class="cover-badge">KG-Microbe \u00d7 BugSigDB · Vector Graphics Slide '
            "Deck (16 Slides)</span>",
            "  <h1>Cellulose-Catabolizing Gut Microbiota Convergence, Four Highlighted Trait S"
            "ignatures, and Three-Arm Statistical Integration</h1>",
            '  <p class="lead">',
            "    Combined vector-graphics deck containing <b>Set 3 (Six Edward Tufte Analytica"
            "l Plates)</b>,",
            "    <b>Set 2 (Four Graphical Summary Figures)</b>, <b>Set 1 (Three Statistical Re"
            "ference Figures)</b>,",
            "    and <b>Two Literature &amp; Set-Theory Reference Tables</b> across Parkinson"
            "\u2019s disease, HIV infection, and colorectal cancer.",
            "  </p>",
            "",
            '  <div class="cover-grid">',
            '    <div class="cover-card">',
            '      <div class="num">50.4%</div>',
            '      <div class="lbl"><b>Mean Pairwise Jaccard Similarity</b><br/>Across PD, H'
            "IV &amp; CRC stool (<i>d</i> = 0.4957; 2.5\u00d7 gut baseline of 20.1%, permutati"
            "on <i>p</i> = 0.011)</div>",
            "    </div>",
            '    <div class="cover-card">',
            '      <div class="num">17 / 49</div>',
            '      <div class="lbl"><b>3-Disease Core Guild (34.7%)</b><br/>10.6\u00d7 over '
            "random expectation (<i>p</i> = 6.2 \u00d7 10⁻¹⁶); 9 core genera (17.1\u00d7 at Ge"
            "nus+Species, <i>p</i> = 3.5 \u00d7 10⁻¹²)</div>",
            "    </div>",
            '    <div class="cover-card">',
            '      <div class="num">73.1%</div>',
            '      <div class="lbl"><b>Stool Decrease Direction</b><br/>106 decrease vs. 39 '
            "increase experiments in stool (<i>p</i> = 2.4 \u00d7 10⁻⁸); absent in mucosa (48."
            "1%) &amp; oral (45.5%)</div>",
            "    </div>",
            '    <div class="cover-card">',
            '      <div class="num">3 Arms</div>',
            '      <div class="lbl"><b>Three Statistical Frameworks</b><br/>Composition (2,3'
            "35 pairs \u00d7 317 traits), Taxon-Calibrated ORA (4,670 profiles \u00d7 317 trai"
            "ts), and Study Direction (939 studies)</div>",
            "    </div>",
            "  </div>",
            "",
            '  <div class="cover-toc">',
            '    <div class="cover-toc-item"><b>Slides 2\u20137 · Set 3 (Edward Tufte Plates'
            " I\u2013VI):</b> Convergence slopegraph, 9-genus dot-bar table, body-site tests, "
            "90-trait benchmark, the 3 statistical arms, and the 1,912-row integrated evidence"
            " contract.</div>",
            '    <div class="cover-toc-item"><b>Slides 8\u201311 · Set 2 (Graphical Summary)'
            ":</b> Convergence triangle &amp; Venn partition, fine-rank sensitivity, stool vs."
            " biopsy directions, and functional pathway rankings.</div>",
            '    <div class="cover-toc-item"><b>Slides 12\u201314 · Set 1 (Statistical Refer'
            "ence):</b> Jaccard matrices against 5 null models, compartment sign tests, and 90"
            "-trait empirical distributions.</div>",
            '    <div class="cover-toc-item"><b>Slides 15\u201316 · Literature &amp; Set-The'
            "ory Reference:</b> Vector tables covering 17 peer-reviewed citations and the 2-wa"
            "y vs. 3-way Venn partition.</div>",
            "  </div>",
            "</section>",
            "",
            '<section class="slide">',
            '  <header class="slide-header">',
            '    <div class="header-left">',
            '      <span class="cat-badge" style="color:#9E2A2B;background:#F9ECEC">Set 3 '
            "· Edward Tufte Plates</span>",
            '      <div class="title-block">',
            "        <h2>Plate I: Three Distinct Diseases Lose the Same Gut Fiber-Digesting Ba"
            "cteria</h2>",
            '        <p class="subtitle">Slopegraph and mutually exclusive set partition vs.'
            " full pairwise overlaps across 90 human stool studies</p>",
            "      </div>",
            "    </div>",
            '    <div class="header-right">',
            '      <span class="meta-pill">Vector SVG · <code>figures_tufte/fig1_tufte_conve'
            "rgence_slope_and_overlap.svg</code></span>",
            '      <span class="slide-num">2 / 16</span>',
            "    </div>",
            "  </header>",
            "",
            '  <div class="fig-canvas">',
            '    <img src="figures_tufte/fig1_tufte_convergence_slope_and_overlap.svg" alt='
            '"Plate I: Three Distinct Diseases Lose the Same Gut Fiber-Digesting Bacteria" /'
            ">",
            "  </div>",
            "",
            '  <footer class="slide-footer">',
            '    <div class="footer-col"><b>Slopegraph:</b> Compares random chance (11%'
            "\u201315%), all 530 stool decreased bacteria (18%\u201324%), and the 49 cellulose"
            "-catabolizing bacteria (45.0%\u201357.4%; mean 50.4%, 2.5\u00d7 gut baseline, per"
            "mutation <i>p</i> = 0.011).</div>",
            '    <div class="footer-col"><b>Set partition:</b> Splitting the 49 lost taxa in'
            "to mutually exclusive tiers shows 17 (34.7%) lost in all 3 diseases and 15 (30.6%"
            ") in only 2 diseases. Every complete 2-disease intersection (27, 21, 18) includes"
            " the 17-taxon core.</div>",
            "  </footer>",
            "</section>",
            "",
            '<section class="slide">',
            '  <header class="slide-header">',
            '    <div class="header-left">',
            '      <span class="cat-badge" style="color:#9E2A2B;background:#F9ECEC">Set 3 '
            "· Edward Tufte Plates</span>",
            '      <div class="title-block">',
            "        <h2>Plate II: Nine Butyrate- and Fiber-Degrading Genera Decrease Across A"
            "ll Three Diseases</h2>",
            '        <p class="subtitle">Taxonomic dot-bar table and Genus + Species sensiti'
            "vity test (N = 93 stool background taxa)</p>",
            "      </div>",
            "    </div>",
            '    <div class="header-right">',
            '      <span class="meta-pill">Vector SVG · <code>figures_tufte/fig2_tufte_core_'
            "guild_dot_table.svg</code></span>",
            '      <span class="slide-num">3 / 16</span>',
            "    </div>",
            "  </header>",
            "",
            '  <div class="fig-canvas">',
            '    <img src="figures_tufte/fig2_tufte_core_guild_dot_table.svg" alt="Plate II'
            ': Nine Butyrate- and Fiber-Degrading Genera Decrease Across All Three Diseases" '
            "/>",
            "  </div>",
            "",
            '  <footer class="slide-footer">',
            '    <div class="footer-col"><b>Core genera and species:</b> <i>Roseburia</i> (2'
            "5 studies: 18·2·5), <i>Blautia</i> (19: 12·2·5), <i>Clostridium</i> (16: 7·4·5), "
            "<i>Bacteroides</i> (15: 6·7·2), <i>Lachnospira</i> (14: 5·3·6), <i>Ruminococcus</"
            "i> (14: 6·4·4), <i>Anaerostipes</i> (8: 2·3·3), <i>Butyrivibrio</i> (4: 2·1·1), a"
            "nd <i>Eubacterium</i> (4: 1·1·2) decrease across PD, HIV, and CRC.</div>",
            '    <div class="footer-col"><b>Rank sensitivity:</b> Restricting to Genus + Spe'
            "cies (N = 93 stool background taxa, 22 cellulose-annotated) increases PD\u2013CRC"
            " Jaccard similarity to <b>62.5%</b> (<i>p</i> = 5.98 \u00d7 10⁻¹⁰) and 3-way core"
            " enrichment to <b>17.1\u00d7</b> (<i>p</i> = 3.5 \u00d7 10⁻¹²).</div>",
            "  </footer>",
            "</section>",
            "",
            '<section class="slide">',
            '  <header class="slide-header">',
            '    <div class="header-left">',
            '      <span class="cat-badge" style="color:#9E2A2B;background:#F9ECEC">Set 3 '
            "· Edward Tufte Plates</span>",
            '      <div class="title-block">',
            "        <h2>Plate III: Shared Loss Appears in Stool, Not in Tissue Biopsies or th"
            "e Mouth</h2>",
            '        <p class="subtitle">Range-frame dot plots of within-study directional s'
            "ign tests and cross-disease Jaccard similarity by body site</p>",
            "      </div>",
            "    </div>",
            '    <div class="header-right">',
            '      <span class="meta-pill">Vector SVG · <code>figures_tufte/fig3_tufte_why_o'
            "nly_stool.svg</code></span>",
            '      <span class="slide-num">4 / 16</span>',
            "    </div>",
            "  </header>",
            "",
            '  <div class="fig-canvas">',
            '    <img src="figures_tufte/fig3_tufte_why_only_stool.svg" alt="Plate III: Sha'
            'red Loss Appears in Stool, Not in Tissue Biopsies or the Mouth" />',
            "  </div>",
            "",
            '  <footer class="slide-footer">',
            '    <div class="footer-col"><b>Within-study direction:</b> In stool, 73.1% of s'
            "ignature experiments report a decrease (106 decrease vs. 39 increase, binomial <i"
            ">p</i> = 2.4 \u00d7 10⁻⁸). Mucosal biopsies (48.1%, <i>p</i> = 1.00) and oral sam"
            "ples (45.5%, <i>p</i> = 1.00) show no directional bias.</div>",
            '    <div class="footer-col"><b>Colonic niche:</b> Insoluble plant cellulose rea'
            "ches the anaerobic lumen of the lower colon sampled by stool, whereas mucosal bio"
            "psies sample biofilm bacteria attached to the gut wall.</div>",
            "  </footer>",
            "</section>",
            "",
            '<section class="slide">',
            '  <header class="slide-header">',
            '    <div class="header-left">',
            '      <span class="cat-badge" style="color:#9E2A2B;background:#F9ECEC">Set 3 '
            "· Edward Tufte Plates</span>",
            '      <div class="title-block">',
            "        <h2>Plate IV: Plant-Fiber Breakdown Ranks First Among Metabolic Pathways "
            "as Salt-Tolerant Taxa Increase</h2>",
            '        <p class="subtitle">90-trait range-frame distribution strip and the fou'
            "r highlighted condition-trait signatures with literature citations</p>",
            "      </div>",
            "    </div>",
            '    <div class="header-right">',
            '      <span class="meta-pill">Vector SVG · <code>figures_tufte/fig4_tufte_trait'
            "_landscape_and_four_pairs.svg</code></span>",
            '      <span class="slide-num">5 / 16</span>',
            "    </div>",
            "  </header>",
            "",
            '  <div class="fig-canvas">',
            '    <img src="figures_tufte/fig4_tufte_trait_landscape_and_four_pairs.svg" alt='
            '"Plate IV: Plant-Fiber Breakdown Ranks First Among Metabolic Pathways as Salt-To'
            'lerant Taxa Increase" />',
            "  </div>",
            "",
            '  <footer class="slide-footer">',
            '    <div class="footer-col"><b>90-trait benchmark:</b> Among all 90 non-media f'
            "unctional traits in stool, cellulose (50.4%), xylan (46.6%), and chitin (45.8%) c"
            "atabolism rank #1, #2, and #3 among metabolic pathways (second overall only to en"
            "dospore formation at 67.0%).</div>",
            '    <div class="footer-col"><b>Four highlighted pairs:</b> Cellulose catabolism'
            " decreases in PD (\u221227.4 pp, <i>q</i> = 0.0077), HIV (\u221253.7 pp, <i>q</i>"
            " = 0.0080), and CRC (\u221236.9 pp, <i>q</i> = 0.0336), while halotolerance incre"
            "ases in PD (+24.4 pp, Fisher BH <i>q</i> = 0.0265; paired Wilcoxon <i>p</i> = 1.6"
            "7 \u00d7 10⁻⁴).</div>",
            "  </footer>",
            "</section>",
            "",
            '<section class="slide">',
            '  <header class="slide-header">',
            '    <div class="header-left">',
            '      <span class="cat-badge" style="color:#9E2A2B;background:#F9ECEC">Set 3 '
            "· Edward Tufte Plates</span>",
            '      <div class="title-block">',
            "        <h2>Plate V: One Matrix, Three Statistical Questions: Composition, Enrich"
            "ment, and Study Direction</h2>",
            '        <p class="subtitle">Filter cascade, 100-fold sample-size scale across s'
            "ampling units, and Table 1 fraction and log-ratio results</p>",
            "      </div>",
            "    </div>",
            '    <div class="header-right">',
            '      <span class="meta-pill">Vector SVG · <code>figures_tufte/fig5_tufte_three'
            "_statistical_arms.svg</code></span>",
            '      <span class="slide-num">6 / 16</span>',
            "    </div>",
            "  </header>",
            "",
            '  <div class="fig-canvas">',
            '    <img src="figures_tufte/fig5_tufte_three_statistical_arms.svg" alt="Plate '
            "V: One Matrix, Three Statistical Questions: Composition, Enrichment, and Study Di"
            'rection" />',
            "  </div>",
            "",
            '  <footer class="slide-footer">',
            '    <div class="footer-col"><b>Three statistical arms:</b> Arm 1 tests paired c'
            "omposition z-scores across 2,335 experiment pairs \u00d7 317 non-MediaDive traits"
            "; Arm 2 tests taxon-calibrated hypergeometric ORA across 4,670 profiles \u00d7 31"
            "7 traits (66,240 enrichments, 23 depletions); Arm 3 tests study-level sign consis"
            "tency across 939 studies.</div>",
            '    <div class="footer-col"><b>Table 1 fraction and ratio tests:</b> All 6 Tabl'
            "e 1 traits pass BH <i>q</i> ≤ 0.05 in the 6-trait fraction sign test (FA, HAL, CE"
            "L pass <i>q</i> ≤ 0.001); FA/CEL (59.5% up), CEL/SPH (40.7% up), and CEL/ARS (36."
            "5% up) pass BY <i>q</i> ≤ 0.05 across 15 log-ratios.</div>",
            "  </footer>",
            "</section>",
            "",
            '<section class="slide">',
            '  <header class="slide-header">',
            '    <div class="header-left">',
            '      <span class="cat-badge" style="color:#9E2A2B;background:#F9ECEC">Set 3 '
            "· Edward Tufte Plates</span>",
            '      <div class="title-block">',
            "        <h2>Plate VI: How Results Are Combined: Shared Comparison Basis and the I"
            "ntegrated Evidence Contract</h2>",
            '        <p class="subtitle">Denominator alignment (317 non-MediaDive traits), s'
            "tudy-equal weighting, and the 1,912-row integrated_evidence.tsv schema</p>",
            "      </div>",
            "    </div>",
            '    <div class="header-right">',
            '      <span class="meta-pill">Vector SVG · <code>figures_tufte/fig6_tufte_combi'
            "ning_results_integrated_evidence.svg</code></span>",
            '      <span class="slide-num">7 / 16</span>',
            "    </div>",
            "  </header>",
            "",
            '  <div class="fig-canvas">',
            '    <img src="figures_tufte/fig6_tufte_combining_results_integrated_evidence.svg'
            '" alt="Plate VI: How Results Are Combined: Shared Comparison Basis and the Inte'
            'grated Evidence Contract" />',
            "  </div>",
            "",
            '  <footer class="slide-footer">',
            '    <div class="footer-col"><b>Shared comparison basis:</b> <code>three-framewo'
            "rk-comparison-v1</code> pins the 317 non-MediaDive trait denominator, 2,241 non-z"
            "ero pairs across 939 studies, equal study weighting, and unimputed positive count"
            "s for log-ratios.</div>",
            '    <div class="footer-col"><b>Integrated evidence contract:</b> <code>integrat'
            "ed_evidence.tsv</code> (1,912 rows = 40 method rows + 1,872 grid cells) keeps FDR"
            " families separate, reports both <i>q</i> ≤ 0.05 and <i>q</i> ≤ 0.001, and classi"
            "fies every blank cell.</div>",
            "  </footer>",
            "</section>",
            "",
            '<section class="slide">',
            '  <header class="slide-header">',
            '    <div class="header-left">',
            '      <span class="cat-badge" style="color:#1E293B;background:#E2E8F0">Set 2 '
            "· Graphical Summary</span>",
            '      <div class="title-block">',
            "        <h2>Graphical Summary 1: Cross-Disease Convergence Network and 3-Way Venn"
            " Partition</h2>",
            '        <p class="subtitle">Disease triangle network, pairwise Jaccard similari'
            "ty gauge, and 49-taxon union breakdown across PD, HIV, and CRC</p>",
            "      </div>",
            "    </div>",
            '    <div class="header-right">',
            '      <span class="meta-pill">Vector SVG · <code>figures_graphical/fig1_surpris'
            "ing_convergence_network_venn.svg</code></span>",
            '      <span class="slide-num">8 / 16</span>',
            "    </div>",
            "  </header>",
            "",
            '  <div class="fig-canvas">',
            '    <img src="figures_graphical/fig1_surprising_convergence_network_venn.svg" a'
            'lt="Graphical Summary 1: Cross-Disease Convergence Network and 3-Way Venn Partit'
            'ion" />',
            "  </div>",
            "",
            '  <footer class="slide-footer">',
            '    <div class="footer-col"><b>Union partition:</b> Across Parkinson\u2019s dis'
            "ease, HIV infection, and colorectal cancer stool studies, 32 of 49 lost cellulose"
            "-degrading taxa (65.3%) are shared by at least two diseases, and 17 (34.7%) decre"
            "ase in all three.</div>",
            '    <div class="footer-col"><b>Pairwise Jaccard overlap:</b> Mean pairwise Jacc'
            "ard similarity is 50.4% (vs. 20.1% across all 530 stool decreased taxa; random-gu"
            "ild permutation <i>p</i> = 0.0110; 3-way core 10.6\u00d7 over random expectation,"
            " <i>p</i> = 6.2 \u00d7 10⁻¹⁶).</div>",
            "  </footer>",
            "</section>",
            "",
            '<section class="slide">',
            '  <header class="slide-header">',
            '    <div class="header-left">',
            '      <span class="cat-badge" style="color:#1E293B;background:#E2E8F0">Set 2 '
            "· Graphical Summary</span>",
            '      <div class="title-block">',
            "        <h2>Graphical Summary 2: Taxonomic Rollup Sensitivity at Genus and Specie"
            "s Ranks</h2>",
            '        <p class="subtitle">Fine-rank sensitivity comparison (N = 93 Genus + Sp'
            "ecies) and shared 17-taxon Genus/Species guild scorecard</p>",
            "      </div>",
            "    </div>",
            '    <div class="header-right">',
            '      <span class="meta-pill">Vector SVG · <code>figures_graphical/fig2_genus_s'
            "pecies_core_guild.svg</code></span>",
            '      <span class="slide-num">9 / 16</span>',
            "    </div>",
            "  </header>",
            "",
            '  <div class="fig-canvas">',
            '    <img src="figures_graphical/fig2_genus_species_core_guild.svg" alt="Graphi'
            'cal Summary 2: Taxonomic Rollup Sensitivity at Genus and Species Ranks" />',
            "  </div>",
            "",
            '  <footer class="slide-footer">',
            '    <div class="footer-col"><b>Fine-rank test:</b> Removing family, order, and '
            "phylum rollups and testing only Genus + Species (N = 93) reduces mean Jaccard dis"
            "tance to <i>d<sub>J</sub></i> = 0.478 (<i>J</i> = 52.2%) and increases PD\u2013CR"
            "C similarity to 62.5% (<i>p</i> = 5.98 \u00d7 10⁻¹⁰).</div>",
            '    <div class="footer-col"><b>Core genera:</b> Nine genera decrease across all'
            " three diseases (17.1\u00d7 over random expectation, <i>p</i> = 3.5 \u00d7 10⁻¹²)"
            ", led by <i>Roseburia</i>, <i>Blautia</i>, <i>Clostridium</i>, <i>Bacteroides</i>"
            ", and <i>Ruminococcus</i>.</div>",
            "  </footer>",
            "</section>",
            "",
            '<section class="slide">',
            '  <header class="slide-header">',
            '    <div class="header-left">',
            '      <span class="cat-badge" style="color:#1E293B;background:#E2E8F0">Set 2 '
            "· Graphical Summary</span>",
            '      <div class="title-block">',
            "        <h2>Graphical Summary 3: Both Depletion and Taxonomic Convergence Are Spe"
            "cific to Stool</h2>",
            '        <p class="subtitle">Diverging butterfly chart of paired experiment dire'
            "ctions and cross-compartment Jaccard overlap</p>",
            "      </div>",
            "    </div>",
            '    <div class="header-right">',
            '      <span class="meta-pill">Vector SVG · <code>figures_graphical/fig3_feces_c'
            "ompartment_specificity.svg</code></span>",
            '      <span class="slide-num">10 / 16</span>',
            "    </div>",
            "  </header>",
            "",
            '  <div class="fig-canvas">',
            '    <img src="figures_graphical/fig3_feces_compartment_specificity.svg" alt="G'
            "raphical Summary 3: Both Depletion and Taxonomic Convergence Are Specific to Stoo"
            'l" />',
            "  </div>",
            "",
            '  <footer class="slide-footer">',
            '    <div class="footer-col"><b>Stool vs. mucosa:</b> Stool samples the anaerobi'
            "c colonic lumen where plant fiber accumulates (73.1% decrease, 50.4% shared taxa)"
            ", whereas mucosal biopsies (18.8% shared) and oral sites (0.0% shared) do not.</d"
            "iv>",
            '    <div class="footer-col"><b>Colorectal biopsy contrast:</b> On colorectal tu'
            "mor biopsies, mucosa-adherent bacteria (such as <i>Bacteroides fragilis</i>) ofte"
            "n increase (6 dec vs. 11 inc) while luminal fiber-degraders decrease in stool (30"
            " dec vs. 7 inc).</div>",
            "  </footer>",
            "</section>",
            "",
            '<section class="slide">',
            '  <header class="slide-header">',
            '    <div class="header-left">',
            '      <span class="cat-badge" style="color:#1E293B;background:#E2E8F0">Set 2 '
            "· Graphical Summary</span>",
            '      <div class="title-block">',
            "        <h2>Graphical Summary 4: Cellulose Catabolism Ranks #1 Among GO Processes"
            ", Pathways, and Chemical Traits</h2>",
            '        <p class="subtitle">Functional pathway ranking and top 10 convergent no'
            "n-MediaDive traits across 90 microbial traits in stool</p>",
            "      </div>",
            "    </div>",
            '    <div class="header-right">',
            '      <span class="meta-pill">Vector SVG · <code>figures_graphical/fig4_top_fun'
            "ctional_pathway_podium.svg</code></span>",
            '      <span class="slide-num">11 / 16</span>',
            "    </div>",
            "  </header>",
            "",
            '  <div class="fig-canvas">',
            '    <img src="figures_graphical/fig4_top_functional_pathway_podium.svg" alt="G'
            "raphical Summary 4: Cellulose Catabolism Ranks #1 Among GO Processes, Pathways, a"
            'nd Chemical Traits" />',
            "  </div>",
            "",
            '  <footer class="slide-footer">',
            '    <div class="footer-col"><b>Plant-polysaccharide pathways:</b> Cellulose (50'
            ".4%), xylan (46.6%), and chitin (45.8%) breakdown occupy the top three metabolic "
            "positions out of 90 functional traits, well above the 90-trait median (30.2%).</d"
            "iv>",
            '    <div class="footer-col"><b>Spore-forming Clostridia:</b> Endospore formatio'
            "n ranks #1 overall (67.0%) because the same anaerobic Clostridia (<i>Ruminococcus"
            "</i>, <i>Clostridium</i>, <i>Blautia</i>) are spore-forming fiber fermenters.</di"
            "v>",
            "  </footer>",
            "</section>",
            "",
            '<section class="slide">',
            '  <header class="slide-header">',
            '    <div class="header-left">',
            '      <span class="cat-badge" style="color:#0F766E;background:#CCFBF1">Set 1 '
            "· Statistical Reference</span>",
            '      <div class="title-block">',
            "        <h2>Statistical Figure 1: Pairwise and 3-Way Jaccard Convergence Against "
            "Five Null Models</h2>",
            '        <p class="subtitle">Quantitative reference for PD (40 taxa), HIV (24 ta'
            "xa), and CRC (34 taxa) stool cellulose-catabolism signatures</p>",
            "      </div>",
            "    </div>",
            '    <div class="header-right">',
            '      <span class="meta-pill">Vector SVG · <code>figures/fig1_pairwise_and_3way'
            "_jaccard_nulls.svg</code></span>",
            '      <span class="slide-num">12 / 16</span>',
            "    </div>",
            "  </header>",
            "",
            '  <div class="fig-canvas">',
            '    <img src="figures/fig1_pairwise_and_3way_jaccard_nulls.svg" alt="Statistic'
            'al Figure 1: Pairwise and 3-Way Jaccard Convergence Against Five Null Models" />',
            "  </div>",
            "",
            '  <footer class="slide-footer">',
            '    <div class="footer-col"><b>Exact and Monte Carlo tests:</b> PD\u2013CRC (<i'
            ">J</i> = 0.5745, <i>d</i> = 0.4255, hypergeometric <i>p</i> = 2.66 \u00d7 10⁻¹³),"
            " PD\u2013HIV (<i>J</i> = 0.4884, <i>d</i> = 0.5116, <i>p</i> = 8.43 \u00d7 10⁻¹²)"
            ", HIV\u2013CRC (<i>J</i> = 0.4500, <i>d</i> = 0.5500, <i>p</i> = 6.72 \u00d7 10⁻¹"
            "¹).</div>",
            '    <div class="footer-col"><b>Study-label permutation null:</b> Permuting dise'
            "ase labels across all 90 stool studies (10,000 iterations) shows that cross-disea"
            "se overlap matches sampling from a single shared fiber-degrading guild.</div>",
            "  </footer>",
            "</section>",
            "",
            '<section class="slide">',
            '  <header class="slide-header">',
            '    <div class="header-left">',
            '      <span class="cat-badge" style="color:#0F766E;background:#CCFBF1">Set 1 '
            "· Statistical Reference</span>",
            '      <div class="title-block">',
            "        <h2>Statistical Figure 2: Anatomical Compartment Specificity and Sensitiv"
            "ity Across Taxonomic Ranks</h2>",
            '        <p class="subtitle">Within-study binomial sign tests, cross-disease Jac'
            "card similarity across 4 body compartments, and rank sensitivity</p>",
            "      </div>",
            "    </div>",
            '    <div class="header-right">',
            '      <span class="meta-pill">Vector SVG · <code>figures/fig2_compartment_speci'
            "ficity_and_sensitivity.svg</code></span>",
            '      <span class="slide-num">13 / 16</span>',
            "    </div>",
            "  </header>",
            "",
            '  <div class="fig-canvas">',
            '    <img src="figures/fig2_compartment_specificity_and_sensitivity.svg" alt="S'
            "tatistical Figure 2: Anatomical Compartment Specificity and Sensitivity Across Ta"
            'xonomic Ranks" />',
            "  </div>",
            "",
            '  <footer class="slide-footer">',
            '    <div class="footer-col"><b>Sample distribution:</b> Stool accounts for 90 o'
            "f 121 studies (46 PD, 17 HIV, 27 CRC) and 145 directional experiments, whereas GI"
            " mucosa has 15 studies (2 PD, 6 HIV, 7 CRC) and oral sites have 9 studies.</div>",
            '    <div class="footer-col"><b>Compartment sign tests:</b> Significant decrease'
            " bias occurs only in stool (PD <i>p</i> = 8.8 \u00d7 10⁻⁴, HIV <i>p</i> = 0.015, "
            "CRC <i>p</i> = 1.9 \u00d7 10⁻⁴); mucosal and oral tests have <i>p</i> ≥ 0.33.</di"
            "v>",
            "  </footer>",
            "</section>",
            "",
            '<section class="slide">',
            '  <header class="slide-header">',
            '    <div class="header-left">',
            '      <span class="cat-badge" style="color:#0F766E;background:#CCFBF1">Set 1 '
            "· Statistical Reference</span>",
            '      <div class="title-block">',
            "        <h2>Statistical Figure 3: 17-Taxon Core Presence Matrix and 90-Trait Func"
            "tional Guild Benchmark</h2>",
            '        <p class="subtitle">Study-count matrix across the 17 three-way core tax'
            "a and empirical ranking across all 90 non-MediaDive functional traits</p>",
            "      </div>",
            "    </div>",
            '    <div class="header-right">',
            '      <span class="meta-pill">Vector SVG · <code>figures/fig3_core_taxa_and_tra'
            "it_benchmark.svg</code></span>",
            '      <span class="slide-num">14 / 16</span>',
            "    </div>",
            "  </header>",
            "",
            '  <div class="fig-canvas">',
            '    <img src="figures/fig3_core_taxa_and_trait_benchmark.svg" alt="Statistical'
            " Figure 3: 17-Taxon Core Presence Matrix and 90-Trait Functional Guild Benchmark"
            '" />',
            "  </div>",
            "",
            '  <footer class="slide-footer">',
            '    <div class="footer-col"><b>Benchmark distribution:</b> Across all 90 non-me'
            "dia traits with ≥10 decreased taxa per condition, mean Jaccard similarity spans 1"
            "4.8% to 67.0% (median 30.2%, IQR 27.1%\u201333.9%); cellulose catabolism stands a"
            "t the 97.8th percentile (50.4%).</div>",
            '    <div class="footer-col"><b>Core guild study counts:</b> <i>Roseburia</i> (2'
            "5 stool studies: 18·2·5), <i>Blautia</i> (19: 12·2·5), <i>Clostridium</i> (16: 7·"
            "4·5), <i>Bacteroides</i> (15: 6·7·2), <i>Lachnospira</i> (14: 5·3·6), and <i>Rumi"
            "nococcus</i> (14: 6·4·4) lead the 17-taxon core.</div>",
            "  </footer>",
            "</section>",
            "",
            '<section class="slide">',
            '  <header class="slide-header">',
            '    <div class="header-left">',
            '      <span class="cat-badge" style="color:#7C2D12;background:#FFEDD5">Litera'
            "ture Audit · Vector Table</span>",
            '      <div class="title-block">',
            "        <h2>Published Literature Support for the Four Highlighted Condition-Trait"
            " Signatures</h2>",
            '        <p class="subtitle">Comparing microbial-trait-profiling (MTP) stool sig'
            "natures against 17 peer-reviewed metagenomic and physiological studies</p>",
            "      </div>",
            "    </div>",
            '    <div class="header-right">',
            '      <span class="meta-pill"><b>4</b> Condition-Trait Pairs</span>',
            '      <span class="meta-pill"><b>90</b> Stool Studies</span>',
            '      <span class="meta-pill"><b>17</b> Citations</span>',
            '      <span class="slide-num">15 / 16</span>',
            "    </div>",
            "  </header>",
            "",
            '  <div class="fig-canvas" style="padding:14px 22px;display:flex;flex-direction'
            ':column;justify-content:space-between;background:#FCFBF7">',
            f"    {fig8_inline}",
            "  </div>",
            "",
            '  <footer class="slide-footer">',
            '    <div class="footer-col"><b>Biological interpretation note:</b> Follows <cod'
            "e>methods/biological_interpretation.md</code>: distinguishes database trait annot"
            "ations of reported taxa from direct <i>in vivo</i> metabolic flux measurements.</"
            "div>",
            '    <div class="footer-col"><b>Files:</b> <code>data/analysis/cellulose_decr/FO'
            "UR_HIGHLIGHTED_PAIRS_LITERATURE_SUPPORT.md</code> and <code>data/analysis/cellulo"
            "se_decr/four_highlighted_pairs_taxon_drivers.tsv</code>.</div>",
            "  </footer>",
            "</section>",
            "",
            '<section class="slide">',
            '  <header class="slide-header">',
            '    <div class="header-left">',
            '      <span class="cat-badge" style="color:#1E3A8A;background:#DBEAFE">Set Pa'
            "rtition · Vector Diagram</span>",
            '      <div class="title-block">',
            "        <h2>Why Mutually Exclusive 2-Way Slices Have Fewer Taxa Than the 3-Way Co"
            "re</h2>",
            '        <p class="subtitle">Comparing full pairwise intersections |A ∩ B| (whic'
            "h include the 3-way core) against 2-way-only Venn slices</p>",
            "      </div>",
            "    </div>",
            '    <div class="header-right">',
            '      <span class="meta-pill"><b>49</b> Union Taxa</span>',
            '      <span class="meta-pill"><b>17</b> 3-Way Core</span>',
            '      <span class="meta-pill"><b>15</b> 2-Way-Only</span>',
            '      <span class="slide-num">16 / 16</span>',
            "    </div>",
            "  </header>",
            "",
            '  <div class="fig-canvas" style="padding:14px 22px;display:flex;flex-direction'
            ':column;justify-content:space-between;background:#FCFBF7">',
            f"    {fig9_inline}",
            "  </div>",
            "",
            '  <footer class="slide-footer">',
            '    <div class="footer-col"><b>Set-theory identity:</b> |A ∩ B| = |A ∩ B ∩ C| +'
            " |(A ∩ B) \\ C|. Plate I (Set 3) and Graphical Summary 1 (Set 2) show both quanti"
            "ties side by side.</div>",
            '    <div class="footer-col"><b>Source data:</b> <code>data/analysis/cellulose_d'
            "ecr/cellulose_decreased_taxa_catalog.tsv</code> (49 union taxa = 17 in 3/3 + 15 i"
            "n 2/3 + 17 in 1/3).</div>",
            "  </footer>",
            "</section>",
            "",
            "</body>",
            "</html>",
            "",
        ]
    )


# --------------------------------------------------------------------------- Orchestration & CLI


def generate_tufte_plates_and_deck(
    analysis_dir: Path = ANALYSIS_DIR,
    *,
    dpi: int = 300,
) -> list[Path]:
    """Generate all Tufte plates, standalone SVGs, combined PDF, and 16-slide HTML deck (#550)."""
    tufte_dir = analysis_dir / "figures_tufte"
    tufte_dir.mkdir(parents=True, exist_ok=True)

    builders = (
        ("fig1_tufte_convergence_slope_and_overlap", build_plate1_convergence_slope_and_overlap),
        ("fig2_tufte_core_guild_dot_table", build_plate2_core_guild_dot_table),
        ("fig3_tufte_why_only_stool", build_plate3_why_only_stool),
        ("fig4_tufte_trait_landscape_and_four_pairs", build_plate4_trait_landscape_and_four_pairs),
        ("fig5_tufte_three_statistical_arms", build_plate5_three_statistical_arms),
        (
            "fig6_tufte_combining_results_integrated_evidence",
            build_plate6_combining_results_integrated_evidence,
        ),
    )

    written: list[Path] = []
    combined_pdf_path = tufte_dir / "cellulose_decr_tufte_figures.pdf"

    with (
        matplotlib.rc_context(_tufte_rcparams()),
        PdfPages(combined_pdf_path, metadata={"CreationDate": None}) as pdf_pages,
    ):
        for stem, builder in builders:
            fig = builder()
            try:
                written.extend(_save_plate(fig, tufte_dir, stem, dpi=dpi))
                pdf_pages.savefig(fig, facecolor=PAPER_BG, metadata={"CreationDate": None})
            finally:
                plt.close(fig)
    written.append(combined_pdf_path)

    svg_outputs = (
        ("fig7_tufte_convergence_network_and_venn.svg", render_fig7_svg()),
        ("fig8_literature_support_four_pairs.svg", render_fig8_svg()),
        ("fig9_set_theory_17_core_vs_pairwise.svg", render_fig9_svg()),
    )
    for filename, svg_content in svg_outputs:
        svg_path = tufte_dir / filename
        svg_path.write_text(svg_content, encoding="utf-8")
        written.append(svg_path)

    deck_path = analysis_dir / "cellulose_decr_slide_deck.html"
    deck_path.write_text(render_slide_deck_html(), encoding="utf-8")
    written.append(deck_path)

    return written


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments for Tufte plate and deck generation."""
    parser = argparse.ArgumentParser(
        description=(
            "Generate the Edward Tufte analytical plates (fig1\u2013fig6), standalone vector SVGs "
            "(fig7\u2013fig9), combined PDF packet, and 16-slide HTML deck under analysis_dir."
        )
    )
    parser.add_argument(
        "--analysis-dir",
        type=Path,
        default=ANALYSIS_DIR,
        help="Target cellulose_decr analysis directory (default: data/analysis/cellulose_decr).",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Raster PNG resolution in dots per inch (default: 300).",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    """CLI entry point for generating the Tufte plates and HTML slide deck."""
    args = parse_args(argv)
    generate_tufte_plates_and_deck(args.analysis_dir, dpi=args.dpi)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
