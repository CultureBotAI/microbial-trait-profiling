"""Create static Cytoscape-inspired PNG graph views for the six shared traits."""

from __future__ import annotations

from pathlib import Path
from textwrap import wrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgb
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, PathPatch
from matplotlib.path import Path as MplPath

from create_six_shared_trait_graphs import (
    GRAPH_VIEWS,
    OUTPUT_DIR,
    REPO_ROOT,
    color_from_signal,
    joined_unique_count,
    read_display_edges,
    value_for_color,
)

MANIFEST_FILENAME = "six_shared_trait_static_graph_manifest.tsv"

SLATE_50 = "#f8fafc"
SLATE_200 = "#e2e8f0"
SLATE_400 = "#94a3b8"
SLATE_500 = "#64748b"
SLATE_700 = "#334155"
SLATE_900 = "#0f172a"
NODE_FILL = "#ffffff"


def wrapped_plot_text(value: str, width: int = 150) -> str:
    """Wrap long figure annotations so they do not widen the saved canvas."""
    return "\n".join(wrap(value, width=width))


def wrap_node_label(value: object, width: int = 22) -> str:
    """Wrap a graph-node label with a compact width."""
    text = str(value or "").strip() or "unknown"
    if len(text) > 42:
        text = f"{text[:39].rstrip()}..."
    return "\n".join(wrap(text, width=max(10, width))) or "unknown"


def weighted_signal(group: pd.DataFrame) -> float:
    """Return count-weighted color intensity for an aggregated graph edge."""
    values = group["color_value"].astype(float)
    weights = group["focus_enriched_profile_count"].astype(float)
    if float(weights.sum()) <= 0:
        return float(values.mean()) if not values.empty else 0.0
    return float(np.average(values, weights=weights))


def ordered_node_labels(stats: pd.DataFrame) -> list[str]:
    """Return node labels ordered by trait breadth and positive-z support."""
    return (
        stats.sort_values(
            ["traits", "n", "links", "rank_score"],
            ascending=[False, False, False, False],
        )
        .index.astype(str)
        .tolist()
    )


def evenly_spaced_positions(
    labels: list[str],
    *,
    x: float,
    top: float = 0.82,
    bottom: float = 0.16,
) -> dict[str, tuple[float, float]]:
    """Return node positions along a vertical column."""
    if not labels:
        return {}
    if len(labels) == 1:
        return {labels[0]: (x, 0.50)}
    return {
        label: (x, y)
        for label, y in zip(labels, np.linspace(top, bottom, len(labels)), strict=True)
    }


def scale_values(
    values: pd.Series,
    *,
    output_min: float,
    output_max: float,
) -> pd.Series:
    """Scale positive values with a square-root transform."""
    numeric = values.astype(float)
    positive = numeric[numeric > 0]
    if positive.empty:
        return pd.Series(output_min, index=values.index, dtype=float)
    min_count = float(positive.min())
    max_count = float(positive.max())
    if max_count <= min_count:
        return pd.Series((output_min + output_max) / 2, index=values.index, dtype=float)
    scaled = (np.sqrt(numeric.clip(lower=min_count)) - np.sqrt(min_count)) / (
        np.sqrt(max_count) - np.sqrt(min_count)
    )
    return output_min + scaled.clip(0.0, 1.0) * (output_max - output_min)


def node_stats(frame: pd.DataFrame, column: str) -> pd.DataFrame:
    """Summarize support for condition or body-site nodes."""
    records: list[dict[str, object]] = []
    for label, group in frame.groupby(column, sort=False):
        records.append(
            {
                column: label,
                "links": len(group),
                "traits": group["trait_id"].nunique(),
                "n": int(group["focus_enriched_profile_count"].sum()),
                "s": joined_unique_count(group["focus_enriched_study_ids"]),
                "rank_score": float(group["rank_score"].max()),
            }
        )
    return pd.DataFrame.from_records(records).set_index(column)


def trait_stats(frame: pd.DataFrame) -> pd.DataFrame:
    """Summarize support for the six central trait nodes."""
    records: list[dict[str, object]] = []
    for trait_id, group in frame.groupby("trait_id", sort=False):
        records.append(
            {
                "trait_id": trait_id,
                "trait_label": str(group["trait_label"].iloc[0]),
                "trait_slug": str(group["trait_slug"].iloc[0]),
                "trait_color": str(group["trait_color"].iloc[0]),
                "trait_rank": int(group["trait_rank"].iloc[0]),
                "condition_count": group["condition"].nunique(),
                "body_site_count": group["body_site_group"].nunique(),
                "links": len(group),
                "n": int(group["focus_enriched_profile_count"].sum()),
                "s": joined_unique_count(group["focus_enriched_study_ids"]),
                "rank_score": float(group["rank_score"].max()),
            }
        )
    return (
        pd.DataFrame.from_records(records)
        .sort_values("trait_rank")
        .set_index("trait_id", drop=False)
    )


def aggregate_graph_edges(frame: pd.DataFrame) -> pd.DataFrame:
    """Aggregate displayed Sankey edges into condition-trait and trait-body edges."""
    records: list[dict[str, object]] = []
    condition_groups = frame.groupby(
        ["condition", "trait_id", "trait_label", "trait_color", "trait_rank"],
        sort=False,
    )
    for (condition, trait_id, trait_label, trait_color, trait_rank), group in condition_groups:
        records.append(
            {
                "edge_type": "condition_trait",
                "source_label": condition,
                "target_label": trait_label,
                "trait_id": trait_id,
                "trait_label": trait_label,
                "trait_color": trait_color,
                "trait_rank": int(trait_rank),
                "n": int(group["focus_enriched_profile_count"].sum()),
                "bg": int(group["background_no_enriched_profile_count"].sum()),
                "s": joined_unique_count(group["focus_enriched_study_ids"]),
                "links": len(group),
                "signal": weighted_signal(group),
                "rank_score": float(group["rank_score"].max()),
            }
        )
    body_site_groups = frame.groupby(
        ["trait_id", "trait_label", "trait_color", "trait_rank", "body_site_group"],
        sort=False,
    )
    for (trait_id, trait_label, trait_color, trait_rank, body_site), group in body_site_groups:
        records.append(
            {
                "edge_type": "trait_body_site",
                "source_label": trait_label,
                "target_label": body_site,
                "trait_id": trait_id,
                "trait_label": trait_label,
                "trait_color": trait_color,
                "trait_rank": int(trait_rank),
                "n": int(group["focus_enriched_profile_count"].sum()),
                "bg": int(group["background_no_enriched_profile_count"].sum()),
                "s": joined_unique_count(group["focus_enriched_study_ids"]),
                "links": len(group),
                "signal": weighted_signal(group),
                "rank_score": float(group["rank_score"].max()),
            }
        )
    return pd.DataFrame.from_records(records)


def trait_positions(traits: pd.DataFrame) -> dict[str, tuple[float, float]]:
    """Return stable positions for six central trait hub nodes."""
    coordinates = {
        0: (0.46, 0.72),
        1: (0.56, 0.72),
        2: (0.44, 0.52),
        3: (0.58, 0.52),
        4: (0.46, 0.32),
        5: (0.56, 0.32),
    }
    return {
        row.trait_label: coordinates[int(row.trait_rank)] for row in traits.itertuples(index=False)
    }


def text_color_for_fill(color: str) -> str:
    """Return readable text color for a filled trait node."""
    r, g, b = to_rgb(color)
    luminance = (0.2126 * r) + (0.7152 * g) + (0.0722 * b)
    return "white" if luminance < 0.47 else SLATE_900


def draw_curve(
    ax: plt.Axes,
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    color: str,
    linewidth: float,
    alpha: float,
    zorder: float,
    bend: float,
) -> None:
    """Draw one Cytoscape-like curved graph edge."""
    x0, y0 = start
    x1, y1 = end
    dx = x1 - x0
    control_offset = bend * abs(dx)
    path = MplPath(
        [
            (x0, y0),
            (x0 + dx * 0.45, y0 + control_offset),
            (x1 - dx * 0.45, y1 + control_offset),
            (x1, y1),
        ],
        [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4],
    )
    ax.add_patch(
        PathPatch(
            path,
            facecolor="none",
            edgecolor=color,
            linewidth=linewidth,
            alpha=alpha,
            capstyle="round",
            joinstyle="round",
            zorder=zorder,
        )
    )


def draw_circle_node(
    ax: plt.Axes,
    *,
    x: float,
    y: float,
    radius: float,
    facecolor: str,
    edgecolor: str,
    linewidth: float = 1.2,
    zorder: float = 5,
) -> None:
    """Draw a circular graph node."""
    ax.add_patch(
        Circle(
            (x, y),
            radius,
            facecolor=facecolor,
            edgecolor=edgecolor,
            linewidth=linewidth,
            zorder=zorder,
        )
    )


def render_static_graph(
    *,
    view_key: str,
    view_label: str,
    condition_mode: str,
    focus_direction: str,
    display_edges: pd.DataFrame,
    output_path: Path,
    edge_table_path: Path,
    node_table_path: Path,
    shared_signal_max: float,
) -> dict[str, object]:
    """Render one static graph PNG and write its audit node/edge tables."""
    frame = display_edges.copy()
    frame["color_value"] = frame.apply(value_for_color, axis=1)
    condition_frame = node_stats(frame, "condition")
    body_site_frame = node_stats(frame, "body_site_group")
    trait_frame = trait_stats(frame)
    graph_edges = aggregate_graph_edges(frame)
    graph_edges["color"] = [
        color_from_signal(row.trait_color, float(row.signal), shared_signal_max)
        for row in graph_edges.itertuples(index=False)
    ]
    graph_edges["width"] = scale_values(graph_edges["n"], output_min=0.45, output_max=5.8)

    condition_order = ordered_node_labels(condition_frame)
    body_site_order = ordered_node_labels(body_site_frame)
    condition_positions = evenly_spaced_positions(condition_order, x=0.15)
    body_site_positions = evenly_spaced_positions(body_site_order, x=0.86)
    center_positions = trait_positions(trait_frame)

    node_records: list[dict[str, object]] = []
    for label, row in condition_frame.iterrows():
        x, y = condition_positions[str(label)]
        node_records.append(
            {
                "node_type": "condition",
                "label": label,
                "x": x,
                "y": y,
                "n": int(row.n),
                "s": int(row.s),
                "links": int(row.links),
                "traits": int(row.traits),
            }
        )
    for row in trait_frame.itertuples(index=False):
        x, y = center_positions[row.trait_label]
        node_records.append(
            {
                "node_type": "trait",
                "label": row.trait_label,
                "trait_id": row.trait_id,
                "x": x,
                "y": y,
                "n": int(row.n),
                "s": int(row.s),
                "links": int(row.links),
                "condition_count": int(row.condition_count),
                "body_site_count": int(row.body_site_count),
            }
        )
    for label, row in body_site_frame.iterrows():
        x, y = body_site_positions[str(label)]
        node_records.append(
            {
                "node_type": "body_site",
                "label": label,
                "x": x,
                "y": y,
                "n": int(row.n),
                "s": int(row.s),
                "links": int(row.links),
                "traits": int(row.traits),
            }
        )
    pd.DataFrame.from_records(node_records).to_csv(node_table_path, sep="\t", index=False)
    graph_edges.to_csv(edge_table_path, sep="\t", index=False)

    fig_height = max(11.0, 4.0 + 0.36 * max(len(condition_order), len(body_site_order)))
    fig, ax = plt.subplots(figsize=(18.0, fig_height))
    fig.subplots_adjust(left=0.035, right=0.965, top=0.90, bottom=0.12)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.axis("off")

    for edge in graph_edges.sort_values("n", ascending=True).itertuples(index=False):
        if edge.edge_type == "condition_trait":
            start = condition_positions[str(edge.source_label)]
            end = center_positions[str(edge.target_label)]
            alpha = 0.36
            zorder = 1.0
        else:
            start = center_positions[str(edge.source_label)]
            end = body_site_positions[str(edge.target_label)]
            alpha = 0.58
            zorder = 1.4
        bend = (float(edge.trait_rank) - 2.5) * 0.010
        draw_curve(
            ax,
            start,
            end,
            color=str(edge.color),
            linewidth=float(edge.width),
            alpha=alpha,
            zorder=zorder,
            bend=bend,
        )

    condition_sizes = scale_values(condition_frame["n"], output_min=0.010, output_max=0.019)
    body_site_sizes = scale_values(body_site_frame["n"], output_min=0.011, output_max=0.021)
    trait_sizes = scale_values(trait_frame["n"], output_min=0.031, output_max=0.043)

    for label in condition_order:
        x, y = condition_positions[label]
        row = condition_frame.loc[label]
        draw_circle_node(
            ax,
            x=x,
            y=y,
            radius=float(condition_sizes.loc[label]),
            facecolor=NODE_FILL,
            edgecolor=SLATE_500,
            zorder=4.0,
        )
        ax.text(
            x - 0.022,
            y,
            (
                f"{wrap_node_label(label)}\n"
                f"links={int(row.links)}, traits={int(row.traits)}, n={int(row.n):,}"
            ),
            ha="right",
            va="center",
            fontsize=7.2,
            color=SLATE_900,
            linespacing=1.03,
        )

    for label in body_site_order:
        x, y = body_site_positions[label]
        row = body_site_frame.loc[label]
        draw_circle_node(
            ax,
            x=x,
            y=y,
            radius=float(body_site_sizes.loc[label]),
            facecolor=NODE_FILL,
            edgecolor=SLATE_500,
            zorder=4.0,
        )
        ax.text(
            x + 0.022,
            y,
            (
                f"{wrap_node_label(label)}\n"
                f"links={int(row.links)}, traits={int(row.traits)}, n={int(row.n):,}"
            ),
            ha="left",
            va="center",
            fontsize=7.2,
            color=SLATE_900,
            linespacing=1.03,
        )

    for row in trait_frame.itertuples(index=False):
        x, y = center_positions[row.trait_label]
        radius = float(trait_sizes.loc[row.trait_id])
        draw_circle_node(
            ax,
            x=x,
            y=y,
            radius=radius,
            facecolor=row.trait_color,
            edgecolor=SLATE_900,
            linewidth=1.4,
            zorder=6.0,
        )
        ax.text(
            x,
            y,
            wrap_node_label(row.trait_label, width=15),
            ha="center",
            va="center",
            fontsize=7.6,
            color=text_color_for_fill(row.trait_color),
            fontweight="bold",
            linespacing=1.02,
            zorder=7.0,
        )
        ax.text(
            x,
            y - radius - 0.018,
            (
                f"n={int(row.n):,}\n"
                f"conditions={int(row.condition_count)}, sites={int(row.body_site_count)}"
            ),
            ha="center",
            va="top",
            fontsize=6.8,
            color=SLATE_700,
            linespacing=1.05,
            zorder=7.0,
        )

    fig.text(
        0.035,
        0.965,
        f"Six shared traits Cytoscape-inspired graph: {focus_direction} profiles",
        ha="left",
        va="top",
        fontsize=17,
        fontweight="bold",
        color=SLATE_900,
    )
    fig.text(
        0.035,
        0.925,
        wrapped_plot_text(
            f"Graph view from displayed Sankey edges ({view_label}). Left nodes are "
            f"conditions using {condition_mode}; center nodes are the six q <= 0.001 "
            "trait hubs; right nodes are body-site groups. Edges are aggregated "
            "condition-to-trait and trait-to-body-site associations, not direct "
            "sankey flows. n = positive-z profile memberships; s = unique studies; "
            "bg = z <= 0 background profiles in the audit TSV.",
            width=156,
        ),
        ha="left",
        va="top",
        fontsize=9.5,
        color=SLATE_700,
        linespacing=1.15,
    )

    ax.text(0.15, 0.93, "Condition nodes", ha="center", va="center", fontsize=11.5, color=SLATE_700)
    ax.text(0.51, 0.93, "Six trait hubs", ha="center", va="center", fontsize=11.5, color=SLATE_700)
    ax.text(0.86, 0.93, "Body-site nodes", ha="center", va="center", fontsize=11.5, color=SLATE_700)

    width_examples = sorted(
        {
            int(graph_edges["n"].min()),
            int(graph_edges["n"].median()),
            int(graph_edges["n"].max()),
        }
    )
    width_lookup = scale_values(
        pd.Series(width_examples, index=width_examples),
        output_min=0.45,
        output_max=5.8,
    )
    legend_handles = [
        Line2D(
            [0],
            [0],
            color=SLATE_700,
            linewidth=float(width_lookup.loc[count]),
            alpha=0.72,
            solid_capstyle="round",
            label=f"edge n={count}",
        )
        for count in width_examples
    ]
    legend = ax.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.50, -0.085),
        ncol=len(legend_handles),
        frameon=False,
        fontsize=8,
        title="Edge widths; trait colors are central hubs",
        title_fontsize=9,
    )
    ax.add_artist(legend)

    fig.text(
        0.035,
        0.035,
        wrapped_plot_text(
            f"Static graph PNG generated from {view_key}. Edge TSV aggregates the "
            "displayed Sankey links into condition-trait and trait-body-site graph "
            "edges; pair-level condition/body-site detail remains in the source "
            "display-edge TSV.",
            width=155,
        ),
        ha="left",
        va="bottom",
        fontsize=8.0,
        color=SLATE_500,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=240, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    return {
        "view_key": view_key,
        "view_label": view_label,
        "focus_direction": focus_direction,
        "condition_mode": condition_mode,
        "condition_node_count": len(condition_order),
        "trait_node_count": len(trait_frame),
        "body_site_node_count": len(body_site_order),
        "graph_edge_count": len(graph_edges),
        "condition_trait_edge_count": int((graph_edges["edge_type"] == "condition_trait").sum()),
        "trait_body_site_edge_count": int((graph_edges["edge_type"] == "trait_body_site").sum()),
        "png_path": output_path.relative_to(REPO_ROOT).as_posix(),
        "node_table_path": node_table_path.relative_to(REPO_ROOT).as_posix(),
        "edge_table_path": edge_table_path.relative_to(REPO_ROOT).as_posix(),
    }


def append_readme_section(output_dir: Path, manifest: pd.DataFrame) -> None:
    """Append or replace the static graph PNG section in the README."""
    readme_path = output_dir / "README.md"
    if not readme_path.exists():
        return
    marker = "## Static Cytoscape-Inspired PNG Graphs"
    content = readme_path.read_text(encoding="utf-8").split(marker, maxsplit=1)[0].rstrip()
    lines = [
        "",
        "",
        marker,
        "",
        f"- Manifest: `{MANIFEST_FILENAME}`",
        "",
        "These PNGs are graph versions of the six shared trait Sankey views. "
        "They use six central trait hubs rather than condition-to-body-site flow "
        "bands, so they are meant to communicate the condition-trait and "
        "trait-body-site network structure in a Cytoscape-style static layout.",
        "",
        "| view | PNG | condition nodes | trait nodes | body-site nodes | graph edges |",
        "| --- | --- | ---: | ---: | ---: | ---: |",
    ]
    for row in manifest.itertuples(index=False):
        lines.append(
            f"| {row.view_label} | `{row.png_path}` | {row.condition_node_count} | "
            f"{row.trait_node_count} | {row.body_site_node_count} | {row.graph_edge_count} |"
        )
    lines.append("")
    readme_path.write_text(content + "\n".join(lines), encoding="utf-8")


def write_six_shared_trait_static_graphs(
    output_dir: Path = OUTPUT_DIR,
) -> pd.DataFrame:
    """Write Cytoscape-inspired static graph PNGs and audit tables."""
    frames = [read_display_edges(output_dir / view.display_edges_file) for view in GRAPH_VIEWS]
    signals = pd.concat(
        [frame.apply(value_for_color, axis=1) for frame in frames],
        ignore_index=True,
    )
    shared_signal_max = max(0.4, float(signals.max() if not signals.empty else 0.0))

    records: list[dict[str, object]] = []
    for view, frame in zip(GRAPH_VIEWS, frames, strict=True):
        stem = f"six_shared_traits_{view.key}_cytoscape_hub_graph"
        records.append(
            render_static_graph(
                view_key=view.key,
                view_label=view.label,
                condition_mode=view.condition_mode,
                focus_direction=view.focus_direction,
                display_edges=frame,
                output_path=output_dir / f"{stem}.png",
                edge_table_path=output_dir / f"{stem}_edges.tsv",
                node_table_path=output_dir / f"{stem}_nodes.tsv",
                shared_signal_max=shared_signal_max,
            )
        )
    manifest = pd.DataFrame.from_records(records)
    manifest.to_csv(output_dir / MANIFEST_FILENAME, sep="\t", index=False)
    append_readme_section(output_dir, manifest)
    return manifest


def main() -> None:
    """Write static Cytoscape-inspired graph outputs."""
    manifest = write_six_shared_trait_static_graphs()
    print(f"Wrote static graph PNGs: {len(manifest)}")
    print(f"Wrote manifest: {OUTPUT_DIR / MANIFEST_FILENAME}")


if __name__ == "__main__":
    main()
