"""Refresh PD trait contributor plots with specificity and co-carriage checks."""

from __future__ import annotations

import itertools
import math
from pathlib import Path
from textwrap import wrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from create_pd_cellulose_analysis import OUTPUT_DIR, load_profile_metadata, read_tsv, write_tsv
from create_pd_six_trait_taxon_heatmap import load_taxon_specificity

STEM = "pd_trait_contributors"
CELLS_PATH = OUTPUT_DIR / f"{STEM}_cells.tsv"
SELECTED_TAXA_PATH = OUTPUT_DIR / f"{STEM}_selected_taxa.tsv"
RANK_COVERAGE_PATH = OUTPUT_DIR / f"{STEM}_rank_coverage.tsv"

RANK_ORDER = {"species": 0, "family": 1}
DIRECTION_ORDER = {"inc": 0, "dec": 1}
DIRECTION_LABELS = {"inc": "increased", "dec": "decreased"}
DIRECTION_COLORS = {"inc": "#b91c1c", "dec": "#1d4ed8"}
GRID_COLOR = "#e5e7eb"
SLATE_500 = "#64748b"
SLATE_700 = "#334155"
SLATE_900 = "#111827"
SPECIFICITY_COLOR_LIMIT = 3.0


def base_trait_label(value: object) -> str:
    """Drop ontology suffixes from a display label."""
    return str(value).split(" (", maxsplit=1)[0]


def clean_existing_specificity(selected: pd.DataFrame) -> pd.DataFrame:
    """Remove previously generated specificity columns before a fresh merge."""
    generated = [
        "pd_feces_profiles",
        "other_feces_profiles",
        "pd_other_feces_counts",
        "pd_feces_studies",
        "other_feces_studies",
        "pd_specificity_fraction",
        "pd_vs_other_feces_log2_ratio",
    ]
    return selected.drop(columns=[column for column in generated if column in selected])


def augment_selected_taxa(selected: pd.DataFrame) -> pd.DataFrame:
    """Add exact reported taxon-label PD-vs-other feces specificity columns."""
    metadata = load_profile_metadata()
    specificity = load_taxon_specificity(pd.Index(selected["taxon"].astype(str)), metadata)
    specificity = specificity.rename(columns={"taxon_id": "taxon"})
    augmented = clean_existing_specificity(selected).merge(specificity, on="taxon", how="left")
    count_columns = [
        "pd_feces_profiles",
        "other_feces_profiles",
        "pd_feces_studies",
        "other_feces_studies",
    ]
    augmented[count_columns] = augmented[count_columns].fillna(0).astype(int)
    augmented["pd_other_feces_counts"] = (
        augmented["pd_feces_profiles"].astype(str)
        + "/"
        + augmented["other_feces_profiles"].astype(str)
    )
    return augmented


def trait_direction_columns(coverage: pd.DataFrame) -> list[tuple[str, str, str]]:
    """Return stable trait-direction columns as (trait_slug, trait_label, direction)."""
    traits = coverage[["trait_slug", "trait_label"]].drop_duplicates().reset_index(drop=True)
    columns: list[tuple[str, str, str]] = []
    for trait in traits.itertuples(index=False):
        for direction in DIRECTION_ORDER:
            columns.append((str(trait.trait_slug), str(trait.trait_label), direction))
    return columns


def plot_rank_panel(
    fig: plt.Figure,
    grid: plt.GridSpec,
    rank: str,
    cells: pd.DataFrame,
    selected: pd.DataFrame,
    coverage: pd.DataFrame,
    columns: list[tuple[str, str, str]],
    share_norm: Normalize,
) -> None:
    """Draw one species/family panel with unresolved bars and specificity strip."""
    taxa = selected[selected["rank"].eq(rank)].copy()
    taxa = taxa.sort_values(
        ["total_count", "filled_cells", "display_name"],
        ascending=[False, False, True],
    ).reset_index(drop=True)
    taxon_order = taxa["taxon"].tolist()
    taxon_y = {taxon: index for index, taxon in enumerate(taxon_order)}
    column_x = {
        (trait_slug, direction): index for index, (trait_slug, _, direction) in enumerate(columns)
    }

    nested = grid.subgridspec(
        2,
        2,
        height_ratios=[1.0, max(2.8, len(taxa) * 0.31)],
        width_ratios=[9.0, 1.35],
        hspace=0.04,
        wspace=0.04,
    )
    bar_ax = fig.add_subplot(nested[0, 0])
    spacer_ax = fig.add_subplot(nested[0, 1])
    ax = fig.add_subplot(nested[1, 0])
    specificity_ax = fig.add_subplot(nested[1, 1], sharey=ax)
    spacer_ax.axis("off")

    rank_coverage = coverage[coverage["rank"].eq(rank)].copy()
    coverage_lookup = {
        (row.trait_slug, row.direction): float(row.unresolved_or_ambiguous_share)
        for row in rank_coverage.itertuples(index=False)
    }
    unresolved = [
        coverage_lookup.get((trait_slug, direction), 0.0) for trait_slug, _, direction in columns
    ]
    x_positions = np.arange(len(columns))
    bar_ax.bar(x_positions, unresolved, color="#94a3b8", width=0.72)
    bar_ax.set_ylim(0, 1)
    bar_ax.set_xlim(-0.5, len(columns) - 0.5)
    bar_ax.set_ylabel("unresolved\nshare", fontsize=8, color=SLATE_700)
    bar_ax.tick_params(axis="x", labelbottom=False, length=0)
    bar_ax.tick_params(axis="y", labelsize=7)
    bar_ax.grid(axis="y", color=GRID_COLOR, linewidth=0.7)
    bar_ax.spines[["top", "right"]].set_visible(False)
    bar_ax.set_title(f"{rank.title()} contributors", loc="left", fontweight="bold")

    rank_cells = cells[cells["rank"].eq(rank)].copy()
    rank_cells["share"] = pd.to_numeric(rank_cells["share"], errors="coerce").fillna(0.0)
    rank_cells["count"] = pd.to_numeric(rank_cells["count"], errors="coerce").fillna(0.0)
    present = rank_cells[rank_cells["count"].gt(0) & rank_cells["taxon"].isin(taxon_order)]
    xs = [column_x[(row.trait_slug, row.direction)] for row in present.itertuples(index=False)]
    ys = [taxon_y[row.taxon] for row in present.itertuples(index=False)]
    sizes = [
        34.0 + 520.0 * math.sqrt(float(share_norm(row.share)))
        for row in present.itertuples(index=False)
    ]
    colors = [DIRECTION_COLORS[row.direction] for row in present.itertuples(index=False)]
    ax.scatter(xs, ys, s=sizes, c=colors, edgecolor="white", linewidth=0.8, alpha=0.88)

    ax.set_xlim(-0.5, len(columns) - 0.5)
    ax.set_ylim(len(taxa) - 0.5, -0.5)
    ax.set_xticks(
        x_positions,
        [
            f"{direction}\n" + "\n".join(wrap(base_trait_label(label), width=12))
            for trait_slug, label, direction in columns
        ],
        fontsize=7.5,
    )
    ax.set_yticks(np.arange(len(taxa)), taxa["display_name"], fontsize=7.5)
    ax.tick_params(length=0, pad=4)
    ax.set_xticks(np.arange(-0.5, len(columns), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(taxa), 1), minor=True)
    ax.grid(which="minor", color=GRID_COLOR, linewidth=0.7)
    ax.tick_params(which="minor", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    specificity_cmap = plt.get_cmap("YlOrBr").copy()
    specificity_cmap.set_bad("#f3f4f6")
    specificity_values = (
        taxa["pd_vs_other_feces_log2_ratio"]
        .replace([np.inf, -np.inf], np.nan)
        .clip(-SPECIFICITY_COLOR_LIMIT, SPECIFICITY_COLOR_LIMIT)
        .to_numpy(dtype=float, copy=True)
    )
    total_specificity_profiles = taxa["pd_feces_profiles"].astype(float) + taxa[
        "other_feces_profiles"
    ].astype(float)
    specificity_values[total_specificity_profiles.eq(0).to_numpy()] = np.nan
    specificity_ax.imshow(
        np.ma.masked_invalid(specificity_values[:, None]),
        cmap=specificity_cmap,
        vmin=-SPECIFICITY_COLOR_LIMIT,
        vmax=SPECIFICITY_COLOR_LIMIT,
        aspect="auto",
    )
    specificity_ax.set_xticks([0], ["PD specificity\nPD/other"], fontsize=8)
    specificity_ax.tick_params(left=False, labelleft=False, length=0, pad=4)
    specificity_ax.set_xticks([-0.5, 0.5], minor=True)
    specificity_ax.set_yticks(np.arange(-0.5, len(taxa), 1), minor=True)
    specificity_ax.grid(which="minor", color="white", linewidth=1)
    specificity_ax.tick_params(which="minor", length=0)
    for y, label in enumerate(taxa["pd_other_feces_counts"]):
        value = specificity_values[y]
        specificity_ax.text(
            0,
            y,
            label,
            ha="center",
            va="center",
            fontsize=7.5,
            color="white" if np.isfinite(value) and value >= 1.8 else SLATE_900,
        )
    for spine in specificity_ax.spines.values():
        spine.set_visible(False)


def plot_contributor_heatmap(
    cells: pd.DataFrame,
    selected: pd.DataFrame,
    coverage: pd.DataFrame,
    output_path: Path,
) -> None:
    """Plot the contributor dot heatmap with PD specificity row annotations."""
    columns = trait_direction_columns(coverage)
    ranks = [
        rank
        for rank in sorted(selected["rank"].unique(), key=lambda value: RANK_ORDER.get(value, 99))
        if rank in RANK_ORDER
    ]
    max_share = max(float(pd.to_numeric(cells["share"], errors="coerce").max()), 0.01)
    share_norm = Normalize(vmin=0.0, vmax=max_share)
    height = 3.2 + sum(max(4.1, len(selected[selected["rank"].eq(rank)]) * 0.37) for rank in ranks)
    fig = plt.figure(figsize=(18.5, height), facecolor="white")
    grid = fig.add_gridspec(
        len(ranks),
        1,
        left=0.22,
        right=0.96,
        top=0.91,
        bottom=0.7 / height,
        hspace=0.25,
        height_ratios=[max(4.1, len(selected[selected["rank"].eq(rank)]) * 0.37) for rank in ranks],
    )
    grids = [grid[index, 0] for index in range(len(ranks))]
    for rank, subgrid in zip(ranks, grids, strict=True):
        plot_rank_panel(fig, subgrid, rank, cells, selected, coverage, columns, share_norm)

    fig.text(
        0.03,
        0.986,
        "PD trait contributors with PD-vs-other feces specificity",
        fontsize=16,
        fontweight="bold",
        va="top",
    )
    subtitle = (
        "Dots mark selected species or family contributors for each trait-direction cell; dot area "
        "scales with the contributor's share of the resolved cell. The top bars show unresolved or "
        "ambiguous rank share for that cell. The right strip checks whether the same exact "
        "reported taxon label appears mostly in PD fecal profiles or in other fecal profiles."
    )
    fig.text(0.03, 0.958, "\n".join(wrap(subtitle, 135)), fontsize=9, color=SLATE_500)
    legend_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markerfacecolor=color,
            markeredgecolor="white",
            markersize=8,
            label=DIRECTION_LABELS[direction],
        )
        for direction, color in DIRECTION_COLORS.items()
    ]
    legend_handles.append(Patch(facecolor="#94a3b8", label="unresolved or ambiguous rank share"))
    fig.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.63, 0.936),
        ncol=3,
        frameon=False,
        fontsize=8.2,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with plt.rc_context({"svg.hashsalt": STEM}):
        fig.savefig(output_path, dpi=240, bbox_inches="tight", facecolor="white")
        fig.savefig(
            output_path.with_suffix(".svg"),
            bbox_inches="tight",
            facecolor="white",
            metadata={"Date": None},
        )
    plt.close(fig)


def build_trait_cocarriage_tables(
    cells: pd.DataFrame,
    selected: pd.DataFrame,
    coverage: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Collapse contributor cells into taxon trait presence and pairwise trait sharing."""
    trait_order = coverage[["trait_slug", "trait_label"]].drop_duplicates().reset_index(drop=True)
    summary = selected.copy()
    summary["taxon_key"] = summary["rank"].astype(str) + "|" + summary["taxon"].astype(str)

    present = cells.copy()
    present["count"] = pd.to_numeric(present["count"], errors="coerce").fillna(0.0)
    present["taxon_key"] = present["rank"].astype(str) + "|" + present["taxon"].astype(str)
    present = present[present["count"].gt(0)]

    for trait in trait_order.itertuples(index=False):
        slug = str(trait.trait_slug)
        keys = set(present.loc[present["trait_slug"].eq(slug), "taxon_key"])
        summary[f"{slug}_present"] = summary["taxon_key"].isin(keys)
    trait_flags = [f"{row.trait_slug}_present" for row in trait_order.itertuples(index=False)]
    summary["present_trait_count"] = summary[trait_flags].sum(axis=1)
    summary = summary.drop(columns=["taxon_key"])

    edges: list[dict[str, object]] = []
    trait_rows = list(trait_order.itertuples(index=False))
    for first, second in itertools.combinations(trait_rows, 2):
        first_slug = str(first.trait_slug)
        second_slug = str(second.trait_slug)
        shared = summary[
            summary[f"{first_slug}_present"] & summary[f"{second_slug}_present"]
        ].copy()
        shared = shared.sort_values(["present_trait_count", "total_count"], ascending=False)
        edges.append(
            {
                "trait_a_slug": first_slug,
                "trait_a_label": base_trait_label(first.trait_label),
                "trait_b_slug": second_slug,
                "trait_b_label": base_trait_label(second.trait_label),
                "shared_taxa": len(shared),
                "species_shared_taxa": int(shared["rank"].eq("species").sum()),
                "family_shared_taxa": int(shared["rank"].eq("family").sum()),
                "example_taxa": "; ".join(shared["display_name"].head(5).astype(str)),
            }
        )
    return summary, pd.DataFrame.from_records(edges)


def plot_trait_cocarriage_network(
    edges: pd.DataFrame,
    taxa_summary: pd.DataFrame,
    coverage: pd.DataFrame,
    output_path: Path,
) -> None:
    """Plot a compact trait-to-trait co-carriage network over selected taxa."""
    traits = coverage[["trait_slug", "trait_label"]].drop_duplicates().reset_index(drop=True)
    slugs = traits["trait_slug"].astype(str).tolist()
    labels = {row.trait_slug: base_trait_label(row.trait_label) for row in traits.itertuples()}
    counts = {
        slug: int(taxa_summary[f"{slug}_present"].sum())
        for slug in slugs
        if f"{slug}_present" in taxa_summary
    }
    angles = np.linspace(math.pi / 2, math.pi / 2 + 2 * math.pi, len(slugs), endpoint=False)
    positions = {
        slug: (math.cos(angle), math.sin(angle)) for slug, angle in zip(slugs, angles, strict=True)
    }

    fig, ax = plt.subplots(figsize=(9.2, 8.0), facecolor="white", constrained_layout=True)
    max_shared = max(int(edges["shared_taxa"].max()), 1)
    for row in edges.itertuples(index=False):
        if int(row.shared_taxa) <= 0:
            continue
        x1, y1 = positions[row.trait_a_slug]
        x2, y2 = positions[row.trait_b_slug]
        width = 0.7 + 5.5 * (int(row.shared_taxa) / max_shared)
        ax.plot(
            [x1, x2],
            [y1, y2],
            color="#94a3b8",
            linewidth=width,
            alpha=0.35 + 0.45 * (int(row.shared_taxa) / max_shared),
            solid_capstyle="round",
            zorder=1,
        )
        ax.text(
            (x1 + x2) / 2,
            (y1 + y2) / 2,
            str(int(row.shared_taxa)),
            ha="center",
            va="center",
            fontsize=8,
            color=SLATE_700,
            bbox={"boxstyle": "round,pad=0.2", "facecolor": "white", "edgecolor": GRID_COLOR},
            zorder=3,
        )

    max_count = max(counts.values(), default=1)
    for slug in slugs:
        x, y = positions[slug]
        size = 900 + 1800 * (counts.get(slug, 0) / max_count)
        ax.scatter(
            [x],
            [y],
            s=size,
            color="#f59e0b",
            edgecolor="white",
            linewidth=1.5,
            zorder=4,
        )
        label = "\n".join(wrap(labels[slug], width=15))
        ax.text(
            x * 1.25,
            y * 1.25,
            f"{label}\n{counts.get(slug, 0)} taxa",
            ha="center",
            va="center",
            fontsize=9,
            color=SLATE_900,
        )

    ax.set_title("Trait co-carriage among selected PD contributor taxa", fontweight="bold")
    ax.text(
        0,
        -1.38,
        "Edge labels count selected species/family labels carrying both trait annotations.",
        ha="center",
        va="center",
        fontsize=9,
        color=SLATE_500,
    )
    ax.set_xlim(-1.55, 1.55)
    ax.set_ylim(-1.48, 1.42)
    ax.set_aspect("equal")
    ax.axis("off")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with plt.rc_context({"svg.hashsalt": f"{STEM}_cocarriage"}):
        fig.savefig(output_path, dpi=240, bbox_inches="tight", facecolor="white")
        fig.savefig(
            output_path.with_suffix(".svg"),
            bbox_inches="tight",
            facecolor="white",
            metadata={"Date": None},
        )
    plt.close(fig)


def write_note(selected: pd.DataFrame, edges: pd.DataFrame, output_path: Path) -> None:
    """Write a concise Markdown note for the refreshed contributor outputs."""
    observed = selected[
        (
            selected["pd_feces_profiles"].astype(float)
            + selected["other_feces_profiles"].astype(float)
        )
        .gt(0)
        .to_numpy()
    ]
    specific = observed.sort_values("pd_vs_other_feces_log2_ratio", ascending=False).head(5)
    shared = edges.sort_values("shared_taxa", ascending=False).head(5)
    lines = [
        "# PD Trait Contributor Dot Heatmap",
        "",
        f"![Contributor heatmap]({STEM}_species_family_dot_heatmap.png)",
        "",
        "The refreshed heatmap keeps the species/family contributor view and adds a PD "
        "specificity strip beside each taxon. The strip reports exact reported taxon-label "
        "profiles in PD feces versus other feces, so it checks whether a contributor is "
        "specific to the PD fecal setting rather than simply common across the fecal corpus.",
        "",
        "Direction 5 is added as a trait co-carriage network:",
        "",
        f"![Trait co-carriage network]({STEM}_trait_cocarriage_network.png)",
        "",
        "Top PD-skewed reported taxon labels:",
        "",
        "| rank | taxon | PD/other feces profiles | log2 ratio |",
        "| --- | --- | ---: | ---: |",
    ]
    for row in specific.itertuples(index=False):
        lines.append(
            f"| {row.rank} | {row.display_name} | {row.pd_other_feces_counts} | "
            f"{row.pd_vs_other_feces_log2_ratio:.2f} |"
        )
    lines.extend(
        [
            "",
            "Largest trait co-carriage pairs among selected contributors:",
            "",
            "| trait pair | shared taxa | examples |",
            "| --- | ---: | --- |",
        ]
    )
    for row in shared.itertuples(index=False):
        lines.append(
            f"| {row.trait_a_label} + {row.trait_b_label} | {int(row.shared_taxa)} | "
            f"{row.example_taxa or 'none'} |"
        )
    lines.extend(
        [
            "",
            "Audit files:",
            "",
            f"- `{STEM}_cells.tsv`: long-form dot heatmap cells.",
            f"- `{STEM}_selected_taxa.tsv`: selected taxa with specificity and co-carriage counts.",
            f"- `{STEM}_taxon_specificity.tsv`: PD-vs-other feces specificity table.",
            f"- `{STEM}_trait_cocarriage_taxa.tsv`: taxon-by-trait presence table.",
            f"- `{STEM}_trait_cocarriage_edges.tsv`: pairwise trait sharing table.",
            f"- `{STEM}_species_family_dot_heatmap.svg`: vector contributor heatmap.",
            f"- `{STEM}_trait_cocarriage_network.svg`: vector co-carriage network.",
        ]
    )
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    """Refresh contributor artifacts from existing contributor tables."""
    cells = read_tsv(CELLS_PATH)
    selected = read_tsv(SELECTED_TAXA_PATH)
    coverage = read_tsv(RANK_COVERAGE_PATH)
    selected = augment_selected_taxa(selected)
    write_tsv(selected, SELECTED_TAXA_PATH)
    specificity_columns = [
        "scope",
        "rank",
        "taxon",
        "name",
        "display_name",
        "pd_feces_profiles",
        "other_feces_profiles",
        "pd_other_feces_counts",
        "pd_feces_studies",
        "other_feces_studies",
        "pd_specificity_fraction",
        "pd_vs_other_feces_log2_ratio",
    ]
    write_tsv(selected[specificity_columns], OUTPUT_DIR / f"{STEM}_taxon_specificity.tsv")
    plot_contributor_heatmap(
        cells,
        selected,
        coverage,
        OUTPUT_DIR / f"{STEM}_species_family_dot_heatmap.png",
    )
    cocarriage_taxa, cocarriage_edges = build_trait_cocarriage_tables(cells, selected, coverage)
    write_tsv(cocarriage_taxa, OUTPUT_DIR / f"{STEM}_trait_cocarriage_taxa.tsv")
    write_tsv(cocarriage_edges, OUTPUT_DIR / f"{STEM}_trait_cocarriage_edges.tsv")
    plot_trait_cocarriage_network(
        cocarriage_edges,
        cocarriage_taxa,
        coverage,
        OUTPUT_DIR / f"{STEM}_trait_cocarriage_network.png",
    )
    write_note(selected, cocarriage_edges, OUTPUT_DIR / f"{STEM}.md")
    print(f"Wrote {len(selected)} selected taxa and {len(cocarriage_edges)} trait-pair edges")
    print(OUTPUT_DIR / f"{STEM}_species_family_dot_heatmap.png")
    print(OUTPUT_DIR / f"{STEM}_trait_cocarriage_network.png")


if __name__ == "__main__":
    main()
