#!/usr/bin/env python3
"""Generate a clean workflow diagram for trait enrichment and Sankey visualization."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath

DEFAULT_PNG_PATH = Path("data/derived/profile_enrichment_guild_assignment_diagram.png")
DEFAULT_SVG_PATH = Path("data/derived/profile_enrichment_guild_assignment_diagram.svg")


def draw_rounded_card(
    ax: plt.Axes,
    x: float,
    y: float,
    width: float,
    height: float,
    facecolor: str,
    edgecolor: str,
    corner_radius: float = 0.015,
    linewidth: float = 1.2,
    zorder: int = 2,
) -> None:
    """Draw a card rectangle with rounded corners."""
    box = patches.FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle=f"round,pad=0.0,rounding_size={corner_radius}",
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=linewidth,
        zorder=zorder,
    )
    ax.add_patch(box)


def draw_step_badge(
    ax: plt.Axes,
    cx: float,
    cy: float,
    num: int,
    color: str,
    radius: float = 0.011,
    zorder: int = 4,
) -> None:
    """Draw a numbered step badge circle."""
    circle = patches.Circle(
        (cx, cy),
        radius=radius,
        facecolor=color,
        edgecolor="none",
        zorder=zorder,
    )
    ax.add_patch(circle)
    ax.text(
        cx,
        cy,
        str(num),
        ha="center",
        va="center",
        fontsize=9.5,
        fontweight="bold",
        color="#ffffff",
        zorder=zorder + 1,
    )


def draw_bezier_ribbon(
    ax: plt.Axes,
    x0: float,
    y0_bot: float,
    y0_top: float,
    x1: float,
    y1_bot: float,
    y1_top: float,
    color: str,
    alpha: float = 0.7,
    zorder: int = 3,
) -> None:
    """Draw a smooth curved Sankey-like ribbon between two vertical segments."""
    dx = (x1 - x0) * 0.5
    path_data = [
        (MplPath.MOVETO, (x0, y0_bot)),
        (MplPath.CURVE4, (x0 + dx, y0_bot)),
        (MplPath.CURVE4, (x1 - dx, y1_bot)),
        (MplPath.CURVE4, (x1, y1_bot)),
        (MplPath.LINETO, (x1, y1_top)),
        (MplPath.CURVE4, (x1 - dx, y1_top)),
        (MplPath.CURVE4, (x0 + dx, y0_top)),
        (MplPath.CURVE4, (x0, y0_top)),
        (MplPath.CLOSEPOLY, (x0, y0_bot)),
    ]
    codes, verts = zip(*path_data, strict=True)
    path = MplPath(verts, codes)
    patch = patches.PathPatch(
        path,
        facecolor=color,
        edgecolor="none",
        alpha=alpha,
        zorder=zorder,
    )
    ax.add_patch(patch)


def draw_mini_sankey_glyph(
    ax: plt.Axes,
    x: float,
    y: float,
    width: float,
    height: float,
    theme: str = "guild",
) -> None:
    """Draw an embedded mini Sankey diagram glyph inside a branch card."""
    node_w = width * 0.08
    left_x = x
    right_x = x + width - node_w

    n_h = height * 0.65
    top_y = y + (height - n_h) / 2

    # Draw left and right node pillars
    draw_rounded_card(
        ax,
        left_x,
        top_y,
        node_w,
        n_h,
        facecolor="#ffffff",
        edgecolor="#cbd5e1",
        corner_radius=0.005,
        linewidth=1.0,
        zorder=4,
    )
    draw_rounded_card(
        ax,
        right_x,
        top_y,
        node_w,
        n_h,
        facecolor="#ffffff",
        edgecolor="#cbd5e1",
        corner_radius=0.005,
        linewidth=1.0,
        zorder=4,
    )

    # Subtitle above
    ax.text(
        x + width / 2,
        top_y + n_h + 0.018,
        "lines = condition -> body site edges\nwidth = profile support (n)",
        ha="center",
        va="bottom",
        fontsize=6.5,
        color="#64748b",
        linespacing=1.2,
        zorder=5,
    )

    # Labels underneath
    ax.text(
        left_x + node_w / 2,
        top_y - 0.012,
        "conditions",
        ha="center",
        va="top",
        fontsize=7,
        color="#64748b",
        zorder=5,
    )
    ax.text(
        right_x + node_w / 2,
        top_y - 0.012,
        "body sites",
        ha="center",
        va="top",
        fontsize=7,
        color="#64748b",
        zorder=5,
    )

    if theme == "guild":
        c_blue = "#60a5fa"
        c_gray = "#e2e8f0"
        # Background grey link
        draw_bezier_ribbon(
            ax,
            left_x + node_w,
            top_y + 0.01,
            top_y + 0.025,
            right_x,
            top_y + n_h - 0.025,
            top_y + n_h - 0.01,
            color=c_gray,
            alpha=0.8,
            zorder=2,
        )
        # Blue link 1 (top to middle)
        draw_bezier_ribbon(
            ax,
            left_x + node_w,
            top_y + n_h - 0.03,
            top_y + n_h - 0.01,
            right_x,
            top_y + n_h - 0.03,
            top_y + n_h - 0.01,
            color=c_blue,
            alpha=0.75,
            zorder=3,
        )
        # Blue link 2 (middle to bottom)
        draw_bezier_ribbon(
            ax,
            left_x + node_w,
            top_y + 0.02,
            top_y + 0.04,
            right_x,
            top_y + 0.02,
            top_y + 0.04,
            color=c_blue,
            alpha=0.75,
            zorder=3,
        )
        # Blue crossing link
        draw_bezier_ribbon(
            ax,
            left_x + node_w,
            top_y + n_h * 0.45,
            top_y + n_h * 0.58,
            right_x,
            top_y + 0.06,
            top_y + 0.08,
            color=c_blue,
            alpha=0.6,
            zorder=3,
        )
    else:
        c_slate = "#94a3b8"
        c_light = "#cbd5e1"
        draw_bezier_ribbon(
            ax,
            left_x + node_w,
            top_y + n_h - 0.035,
            top_y + n_h - 0.01,
            right_x,
            top_y + n_h - 0.035,
            top_y + n_h - 0.01,
            color=c_slate,
            alpha=0.6,
            zorder=3,
        )
        draw_bezier_ribbon(
            ax,
            left_x + node_w,
            top_y + 0.015,
            top_y + 0.035,
            right_x,
            top_y + 0.015,
            top_y + 0.035,
            color=c_slate,
            alpha=0.6,
            zorder=3,
        )
        draw_bezier_ribbon(
            ax,
            left_x + node_w,
            top_y + n_h * 0.45,
            top_y + n_h * 0.6,
            right_x,
            top_y + 0.04,
            top_y + 0.06,
            color=c_light,
            alpha=0.6,
            zorder=2,
        )
        draw_bezier_ribbon(
            ax,
            left_x + node_w,
            top_y + 0.04,
            top_y + 0.06,
            right_x,
            top_y + n_h * 0.45,
            top_y + n_h * 0.6,
            color=c_light,
            alpha=0.6,
            zorder=2,
        )


def plot_workflow_diagram(
    output_png: Path,
    output_svg: Path | None = None,
) -> None:
    """Generate the publication-ready workflow infographic diagram."""
    fig, ax = plt.subplots(figsize=(16, 9), dpi=300)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # Canvas Background
    bg = patches.Rectangle((0, 0), 1, 1, facecolor="#f8fafc", edgecolor="none", zorder=0)
    ax.add_patch(bg)

    # 1. Header Section
    ax.text(
        0.035,
        0.95,
        "BugSigDB data to trait-enrichment Sankey diagrams",
        fontsize=21,
        fontweight="bold",
        color="#0f172a",
        va="top",
        ha="left",
    )
    subtitle = (
        "End-to-end analytical pipeline transforming BugSigDB human microbiome signatures "
        "and KG-Microbe functional traits\ninto hypergeometric enrichment profiles, "
        "directional contrasts, and condition-to-anatomy Sankey maps."
    )
    ax.text(
        0.035,
        0.905,
        subtitle,
        fontsize=8.5,
        color="#475569",
        va="top",
        ha="left",
        linespacing=1.35,
    )

    # 2. Top Process Cards (Steps 1 - 5)
    step_defs = [
        {
            "num": 1,
            "title": "BugSigDB Ingestion",
            "body": (
                "• Ingest raw study signatures\n"
                "• Filter for Homo sapiens hosts\n"
                "• Standardize taxonomy & NCBI IDs\n"
                "• LinkML schema data validation"
            ),
            "bg": "#f8fafc",
            "border": "#cbd5e1",
            "badge": "#1e293b",
        },
        {
            "num": 2,
            "title": "Trait Mapping",
            "body": (
                "• KG-Microbe graph integration\n"
                "• Taxa mapped to traits & guilds\n"
                "• Build wide study presence matrix\n"
                "• Clean growth medium exclusions"
            ),
            "bg": "#f0fdf4",
            "border": "#86efac",
            "badge": "#16a34a",
        },
        {
            "num": 3,
            "title": "Statistical Testing",
            "body": (
                "• Hypergeometric over-representation\n"
                "• One-sided exact p-value tests\n"
                "• Filter: p < 0.001 & fold > 1.5\n"
                "• Benjamini-Hochberg FDR q <= 0.001"
            ),
            "bg": "#fefce8",
            "border": "#fde047",
            "badge": "#ca8a04",
        },
        {
            "num": 4,
            "title": "Rollup Aggregation",
            "body": (
                "• Condition & anatomy rollups\n"
                "• Multi-study edge consolidation\n"
                "• Contrast inc vs dec directions\n"
                "• Generate structured audit TSVs"
            ),
            "bg": "#eff6ff",
            "border": "#93c5fd",
            "badge": "#2563eb",
        },
        {
            "num": 5,
            "title": "Visual Diagrams",
            "body": (
                "• Curved Bézier Sankey ribbons\n"
                "• Ribbon width = profile count (n)\n"
                "• Grouped 6-trait & guild maps\n"
                "• Tripartite network diagrams"
            ),
            "bg": "#fff1f2",
            "border": "#fca5a5",
            "badge": "#e11d48",
        },
    ]

    card_y = 0.655
    card_h = 0.185
    card_w = 0.162
    spacing = 0.024
    start_x = 0.033

    for i, step in enumerate(step_defs):
        x = start_x + i * (card_w + spacing)
        draw_rounded_card(
            ax,
            x,
            card_y,
            card_w,
            card_h,
            facecolor=step["bg"],
            edgecolor=step["border"],
            corner_radius=0.012,
        )

        badge_x = x + 0.020
        badge_y = card_y + card_h - 0.028
        draw_step_badge(ax, badge_x, badge_y, step["num"], color=step["badge"], radius=0.010)

        # Title
        ax.text(
            badge_x + 0.014,
            badge_y,
            step["title"],
            fontsize=9.2,
            fontweight="bold",
            color="#0f172a",
            va="center",
            ha="left",
            zorder=4,
        )

        # Body
        ax.text(
            x + 0.010,
            badge_y - 0.022,
            step["body"],
            fontsize=7.6,
            color="#334155",
            va="top",
            ha="left",
            linespacing=1.35,
            zorder=4,
        )

        # Arrow connector to next step
        if i < len(step_defs) - 1:
            arrow_x = x + card_w + spacing * 0.45
            arrow_y = card_y + card_h * 0.5
            ax.annotate(
                "",
                xy=(arrow_x + 0.005, arrow_y),
                xytext=(arrow_x - 0.005, arrow_y),
                arrowprops={
                    "arrowstyle": "-|>",
                    "color": "#94a3b8",
                    "lw": 1.2,
                    "mutation_scale": 10,
                },
                zorder=3,
            )

    # 3. Midsection Text Header
    ax.text(
        0.035,
        0.598,
        "The workflow produces dual Sankey & Network visualization products from shared profiles.",
        fontsize=10.2,
        fontweight="bold",
        color="#0f172a",
        va="top",
        ha="left",
    )
    ax.text(
        0.035,
        0.575,
        (
            "Trait-enrichment guild & 6-trait Sankeys highlight key microbial phenotypes; "
            "tripartite networks provide unbiased condition-to-anatomy overviews."
        ),
        fontsize=7.8,
        color="#64748b",
        va="top",
        ha="left",
    )

    # 4. Two-Branch Cards
    branch_y = 0.295
    branch_h = 0.245
    branch_w = 0.448
    bx1 = 0.033
    bx2 = 0.518

    # Branch A: Trait-Enrichment & Guild Sankeys
    draw_rounded_card(
        ax,
        bx1,
        branch_y,
        branch_w,
        branch_h,
        facecolor="#f0f7ff",
        edgecolor="#93c5fd",
        corner_radius=0.012,
    )
    ax.text(
        bx1 + 0.016,
        branch_y + branch_h - 0.024,
        "Trait-Enrichment & Guild Sankeys",
        fontsize=11.0,
        fontweight="bold",
        color="#1d4ed8",
        va="top",
        ha="left",
        zorder=4,
    )
    ax.text(
        bx1 + 0.016,
        branch_y + branch_h - 0.052,
        "Each curved ribbon links a disease condition group to\n"
        "key microbial traits (aerobe, anaerobe, spore, biofilm)\n"
        "and sampled body-site targets. Ribbon width reflects\n"
        "study/profile count or sample support.",
        fontsize=7.8,
        color="#334155",
        va="top",
        ha="left",
        linespacing=1.35,
        zorder=4,
    )

    pill1_y = branch_y + 0.062
    pill2_y = branch_y + 0.020
    draw_rounded_card(
        ax,
        bx1 + 0.016,
        pill1_y,
        0.265,
        0.034,
        facecolor="#ffffff",
        edgecolor="#bfdbfe",
        corner_radius=0.008,
    )
    ax.text(
        bx1 + 0.024,
        pill1_y + 0.017,
        "Grouped 6-Trait & Guild Sankey Figures",
        fontsize=7.4,
        fontweight="bold",
        color="#1e40af",
        va="center",
        ha="left",
        zorder=5,
    )

    draw_rounded_card(
        ax,
        bx1 + 0.016,
        pill2_y,
        0.265,
        0.034,
        facecolor="#ffffff",
        edgecolor="#bfdbfe",
        corner_radius=0.008,
    )
    ax.text(
        bx1 + 0.024,
        pill2_y + 0.017,
        "Enrichment rank TSVs & FDR audit tables",
        fontsize=7.4,
        fontweight="bold",
        color="#1e40af",
        va="center",
        ha="left",
        zorder=5,
    )

    # Mini-Sankey for Branch A
    draw_mini_sankey_glyph(ax, bx1 + 0.295, branch_y + 0.025, 0.138, 0.155, theme="guild")

    # Branch B: Directional & Tripartite Networks
    draw_rounded_card(
        ax,
        bx2,
        branch_y,
        branch_w,
        branch_h,
        facecolor="#f8fafc",
        edgecolor="#cbd5e1",
        corner_radius=0.012,
    )
    ax.text(
        bx2 + 0.016,
        branch_y + branch_h - 0.024,
        "Directional & Tripartite Networks",
        fontsize=11.0,
        fontweight="bold",
        color="#0f172a",
        va="top",
        ha="left",
        zorder=4,
    )
    ax.text(
        bx2 + 0.016,
        branch_y + branch_h - 0.052,
        "Direct condition-to-body-site mappings across increased\n"
        "and decreased profiles, paired with Cytoscape-inspired\n"
        "tripartite networks (Condition -> Trait -> Anatomy)\n"
        "without guild-specific filtering.",
        fontsize=7.8,
        color="#334155",
        va="top",
        ha="left",
        linespacing=1.35,
        zorder=4,
    )

    draw_rounded_card(
        ax,
        bx2 + 0.016,
        pill1_y,
        0.265,
        0.034,
        facecolor="#ffffff",
        edgecolor="#e2e8f0",
        corner_radius=0.008,
    )
    ax.text(
        bx2 + 0.024,
        pill1_y + 0.017,
        "Increased / Decreased Tripartite Networks",
        fontsize=7.4,
        fontweight="bold",
        color="#1e293b",
        va="center",
        ha="left",
        zorder=5,
    )

    draw_rounded_card(
        ax,
        bx2 + 0.016,
        pill2_y,
        0.265,
        0.034,
        facecolor="#ffffff",
        edgecolor="#e2e8f0",
        corner_radius=0.008,
    )
    ax.text(
        bx2 + 0.024,
        pill2_y + 0.017,
        "Clustered Heatmaps & Frequency Histograms",
        fontsize=7.4,
        fontweight="bold",
        color="#1e293b",
        va="center",
        ha="left",
        zorder=5,
    )

    # Mini-Sankey for Branch B
    draw_mini_sankey_glyph(ax, bx2 + 0.295, branch_y + 0.025, 0.138, 0.155, theme="direction")

    # 5. Bottom Shared Infrastructure Cards (4 items)
    bot_defs = [
        {
            "title": "LinkML Data Architecture",
            "body": (
                "Schema-defined study models\n"
                "with Pydantic runtime checks\n"
                "and strict taxon verification\n"
                "across all study YAMLs."
            ),
        },
        {
            "title": "Taxa-Trait Mapping",
            "body": (
                "Knowledge graph integration\n"
                "linking NCBI microbial taxa\n"
                "to standardized phenotypic\n"
                "and metabolic traits."
            ),
        },
        {
            "title": "Hypergeometric Pipeline",
            "body": (
                "Exact one-sided Fisher/\n"
                "hypergeometric over-representation\n"
                "and depletion tests with FDR\n"
                "Benjamini-Hochberg correction."
            ),
        },
        {
            "title": "Visual Reporting Suite",
            "body": (
                "Multi-tiered reporting including\n"
                "grouped 6-trait Sankey maps,\n"
                "tripartite networks, clustered\n"
                "heatmaps, and ranking tables."
            ),
        },
    ]

    bot_y = 0.08
    bot_h = 0.165
    bot_w = 0.218
    bot_spacing = 0.020
    bot_start_x = 0.033

    for i, bitem in enumerate(bot_defs):
        bx = bot_start_x + i * (bot_w + bot_spacing)
        draw_rounded_card(
            ax,
            bx,
            bot_y,
            bot_w,
            bot_h,
            facecolor="#ffffff",
            edgecolor="#cbd5e1",
            corner_radius=0.012,
        )

        ax.text(
            bx + 0.012,
            bot_y + bot_h - 0.024,
            bitem["title"],
            fontsize=9.2,
            fontweight="bold",
            color="#0f172a",
            va="top",
            ha="left",
            zorder=4,
        )

        ax.text(
            bx + 0.012,
            bot_y + bot_h - 0.052,
            bitem["body"],
            fontsize=7.8,
            color="#475569",
            va="top",
            ha="left",
            linespacing=1.35,
            zorder=4,
        )

        if i < len(bot_defs) - 1:
            arrow_x = bx + bot_w + bot_spacing * 0.45
            arrow_y = bot_y + bot_h * 0.5
            ax.annotate(
                "",
                xy=(arrow_x + 0.005, arrow_y),
                xytext=(arrow_x - 0.005, arrow_y),
                arrowprops={
                    "arrowstyle": "-|>",
                    "color": "#94a3b8",
                    "lw": 1.2,
                    "mutation_scale": 10,
                },
                zorder=3,
            )

    # 6. Footer Section
    ax.text(
        0.5,
        0.035,
        (
            "Summary workflow of the trait-profile enrichment and visualization pipeline. "
            "Reproducible via Justfile recipe 'workflow-diagram'."
        ),
        fontsize=7.5,
        color="#64748b",
        va="center",
        ha="center",
    )

    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(output_png.resolve()), dpi=300, bbox_inches="tight")
    print(f"Saved workflow diagram PNG to {output_png}")

    if output_svg:
        output_svg.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(output_svg.resolve()), format="svg", bbox_inches="tight")
        print(f"Saved workflow diagram SVG to {output_svg}")

    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate workflow infographic diagram for BugSigDB trait enrichment pipeline."
        )
    )
    parser.add_argument(
        "--output-png",
        type=Path,
        default=DEFAULT_PNG_PATH,
        help="Path to save output PNG",
    )
    parser.add_argument(
        "--output-svg",
        type=Path,
        default=DEFAULT_SVG_PATH,
        help="Path to save output SVG",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    plot_workflow_diagram(output_png=args.output_png, output_svg=args.output_svg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
