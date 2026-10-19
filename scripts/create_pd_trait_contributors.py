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

from create_pd_cellulose_analysis import (
    FECES,
    OUTPUT_DIR,
    PD_CONDITION,
    load_profile_metadata,
    read_tsv,
    write_tsv,
)
from create_pd_six_trait_taxon_heatmap import (
    load_taxon_specificity,
    require_current_trait_taxa_evidence,
)

STEM = "pd_trait_contributors"

RANK_ORDER = {"species": 0, "family": 1}
DIRECTION_ORDER = {"inc": 0, "dec": 1}
DIRECTION_LABELS = {"inc": "increased", "dec": "decreased"}
DIRECTION_COLORS = {"inc": "#b91c1c", "dec": "#1d4ed8"}
GRID_COLOR = "#e5e7eb"
SLATE_500 = "#64748b"
SLATE_700 = "#334155"
SLATE_900 = "#111827"
SPECIFICITY_PREVALENCE_DIFFERENCE_LIMIT = 0.05
MIN_LARGE_DIRECTION_SHARE_DIFFERENCE = 0.05
UNRESOLVED_COLOR = "#cbd5e1"
BACKGROUND_CIRCLE_AREA = 680.0
SHARE_TOLERANCE = 1e-12


def base_trait_label(value: object) -> str:
    """Drop ontology suffixes from a display label."""
    return str(value).split(" (", maxsplit=1)[0]


def strip_svg_trailing_whitespace(path: Path) -> None:
    """Normalize Matplotlib's multiline path output for clean generated diffs."""
    lines = path.read_text(encoding="utf-8").splitlines()
    path.write_text("\n".join(line.rstrip() for line in lines) + "\n", encoding="utf-8")


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
        "pd_feces_total_profiles",
        "other_feces_total_profiles",
        "pd_feces_prevalence",
        "other_feces_prevalence",
        "pd_vs_other_feces_log2_prevalence_ratio",
        "pd_minus_other_feces_prevalence_difference",
        "pd_other_feces_display",
    ]
    return selected.drop(columns=[column for column in generated if column in selected])


def augment_selected_taxa(
    selected: pd.DataFrame,
    metadata: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Add pool-normalized PD-vs-other feces occurrence context."""
    if metadata is None:
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
    feces = metadata[
        metadata["body_site_group"].eq(FECES) & metadata["inc_dec"].isin(DIRECTION_ORDER)
    ]
    pd_total = int(feces.loc[feces["condition"].eq(PD_CONDITION), "profile_id"].nunique())
    other_total = int(feces.loc[~feces["condition"].eq(PD_CONDITION), "profile_id"].nunique())
    if pd_total <= 0 or other_total <= 0:
        raise ValueError("PD and other-feces profile pools must both be non-empty")
    augmented["pd_feces_total_profiles"] = pd_total
    augmented["other_feces_total_profiles"] = other_total
    augmented["pd_feces_prevalence"] = augmented["pd_feces_profiles"] / pd_total
    augmented["other_feces_prevalence"] = augmented["other_feces_profiles"] / other_total
    augmented["pd_minus_other_feces_prevalence_difference"] = (
        augmented["pd_feces_prevalence"] - augmented["other_feces_prevalence"]
    )
    augmented["pd_other_feces_display"] = np.where(
        augmented["other_feces_profiles"].gt(0),
        augmented["pd_other_feces_counts"],
        "---",
    )
    pd_smoothed = (augmented["pd_feces_profiles"] + 0.5) / (pd_total + 1.0)
    other_smoothed = (augmented["other_feces_profiles"] + 0.5) / (other_total + 1.0)
    augmented["pd_vs_other_feces_log2_prevalence_ratio"] = np.log2(pd_smoothed / other_smoothed)
    return augmented


def trait_direction_columns(coverage: pd.DataFrame) -> list[tuple[str, str, str]]:
    """Return stable trait-direction columns as (trait_slug, trait_label, direction)."""
    traits = coverage[["trait_slug", "trait_label"]].drop_duplicates().reset_index(drop=True)
    columns: list[tuple[str, str, str]] = []
    for trait in traits.itertuples(index=False):
        for direction in DIRECTION_ORDER:
            columns.append((str(trait.trait_slug), str(trait.trait_label), direction))
    return columns


def build_direction_differences(
    cells: pd.DataFrame,
    coverage: pd.DataFrame,
) -> pd.DataFrame:
    """Keep taxon-trait pairs whose full-cell shares differ by direction."""
    trait_order = {
        str(row.trait_slug): index
        for index, row in enumerate(
            coverage[["trait_slug", "trait_label"]].drop_duplicates().itertuples(index=False)
        )
    }
    key_columns = [
        "scope",
        "rank",
        "taxon",
        "name",
        "display_name",
        "trait_key",
        "trait_slug",
        "trait_label",
    ]
    records: list[dict[str, object]] = []
    for keys, group in cells.groupby(key_columns, sort=False, dropna=False):
        if group["direction"].duplicated().any():
            raise ValueError(f"Duplicate direction rows for contributor pair {keys}")
        by_direction = group.set_index("direction")
        missing = [
            direction for direction in DIRECTION_ORDER if direction not in by_direction.index
        ]
        if missing:
            raise ValueError(f"Missing {missing} direction rows for contributor pair {keys}")
        increased = by_direction.loc["inc"]
        decreased = by_direction.loc["dec"]
        inc_share = float(increased["share"])
        dec_share = float(decreased["share"])
        share_difference = inc_share - dec_share
        if math.isclose(share_difference, 0.0, rel_tol=0.0, abs_tol=SHARE_TOLERANCE):
            continue
        inc_count = int(increased["count"])
        dec_count = int(decreased["count"])
        record = dict(zip(key_columns, keys, strict=True))
        record.update(
            {
                "increased_count": inc_count,
                "decreased_count": dec_count,
                "increased_share_of_full_cell": inc_share,
                "decreased_share_of_full_cell": dec_share,
                "increased_minus_decreased_share": share_difference,
                "increased_studies": int(increased["studies"]),
                "decreased_studies": int(decreased["studies"]),
                "direction_pattern": (
                    "larger increased share" if share_difference > 0 else "larger decreased share"
                ),
                "direction_exclusive": (inc_count > 0) != (dec_count > 0),
            }
        )
        records.append(record)

    differences = pd.DataFrame.from_records(records)
    if differences.empty:
        return differences
    differences["rank_order"] = differences["rank"].map(RANK_ORDER).fillna(99)
    differences["trait_order"] = differences["trait_slug"].map(trait_order).fillna(99)
    differences["absolute_share_difference"] = differences["increased_minus_decreased_share"].abs()
    differences = differences.sort_values(
        ["rank_order", "trait_order", "absolute_share_difference", "display_name"],
        ascending=[True, True, False, True],
    )
    return differences.drop(columns=["rank_order", "trait_order", "absolute_share_difference"])


def select_large_difference_taxa(
    selected: pd.DataFrame,
    differences: pd.DataFrame,
) -> pd.DataFrame:
    """Retain taxa with at least one direction-share difference of five points."""
    if differences.empty:
        raise ValueError("No direction differences are available for heatmap filtering")
    maxima = differences.assign(
        absolute_direction_share_difference=differences["increased_minus_decreased_share"].abs()
    )
    maxima = (
        maxima.groupby(["rank", "taxon"], as_index=False)["absolute_direction_share_difference"]
        .max()
        .rename(
            columns={
                "absolute_direction_share_difference": ("max_absolute_direction_share_difference")
            }
        )
    )
    heatmap_taxa = selected.merge(
        maxima,
        on=["rank", "taxon"],
        how="left",
        validate="one_to_one",
    )
    heatmap_taxa["max_absolute_direction_share_difference"] = heatmap_taxa[
        "max_absolute_direction_share_difference"
    ].fillna(0.0)
    heatmap_taxa = heatmap_taxa[
        heatmap_taxa["max_absolute_direction_share_difference"].ge(
            MIN_LARGE_DIRECTION_SHARE_DIFFERENCE
        )
    ].copy()
    if heatmap_taxa.empty:
        raise ValueError("No taxa pass the large direction-share difference threshold")
    return heatmap_taxa.sort_values(
        ["rank", "max_absolute_direction_share_difference", "display_name"],
        ascending=[True, False, True],
    ).reset_index(drop=True)


def direction_background_circles(
    coverage: pd.DataFrame,
    columns: list[tuple[str, str, str]],
) -> pd.DataFrame:
    """Size paired direction backgrounds within trait and retain unresolved mass in grey."""
    indexed = coverage.set_index(["trait_slug", "direction"])
    records: list[dict[str, object]] = []
    for trait_slug, _, direction in columns:
        row = indexed.loc[(trait_slug, direction)]
        pair = indexed.loc[trait_slug]
        pair_max = max(float(pair["total_count"].max()), 1.0)
        total_count = float(row["total_count"])
        resolved_share = float(row["resolved_share"])
        outer_area = BACKGROUND_CIRCLE_AREA * total_count / pair_max
        records.append(
            {
                "trait_slug": trait_slug,
                "direction": direction,
                "total_count": total_count,
                "resolved_share": resolved_share,
                "unresolved_or_ambiguous_share": float(row["unresolved_or_ambiguous_share"]),
                "outer_area": outer_area,
                "resolved_area": outer_area * resolved_share,
            }
        )
    return pd.DataFrame.from_records(records)


def plot_rank_panel(
    fig: plt.Figure,
    grid: plt.GridSpec,
    rank: str,
    cells: pd.DataFrame,
    differences: pd.DataFrame,
    selected: pd.DataFrame,
    coverage: pd.DataFrame,
    columns: list[tuple[str, str, str]],
    share_norm: Normalize,
) -> None:
    """Draw one rank panel with paired direction backgrounds and changed contributor pairs."""
    taxa = selected[selected["rank"].eq(rank)].copy()
    taxa = taxa.sort_values(
        ["total_count", "filled_cells", "display_name"],
        ascending=[False, False, True],
    ).reset_index(drop=True)
    taxon_order = taxa["taxon"].tolist()
    taxon_y = {taxon: index + 1 for index, taxon in enumerate(taxon_order)}
    column_x = {
        (trait_slug, direction): index for index, (trait_slug, _, direction) in enumerate(columns)
    }

    nested = grid.subgridspec(
        1,
        2,
        width_ratios=[9.0, 1.65],
        wspace=0.04,
    )
    ax = fig.add_subplot(nested[0, 0])
    specificity_ax = fig.add_subplot(nested[0, 1], sharey=ax)

    rank_coverage = coverage[coverage["rank"].eq(rank)].copy()
    x_positions = np.arange(len(columns))
    for trait_index in range(0, len(columns), 2):
        if (trait_index // 2) % 2:
            ax.axvspan(trait_index - 0.5, trait_index + 1.5, color="#f8fafc", zorder=0)

    backgrounds = direction_background_circles(rank_coverage, columns)
    ax.scatter(
        x_positions,
        np.zeros(len(columns)),
        s=backgrounds["outer_area"],
        color=UNRESOLVED_COLOR,
        edgecolor="white",
        linewidth=0.9,
        zorder=2,
    )
    ax.scatter(
        x_positions,
        np.zeros(len(columns)),
        s=backgrounds["resolved_area"],
        c=[DIRECTION_COLORS[direction] for _, _, direction in columns],
        edgecolor="white",
        linewidth=0.6,
        alpha=0.9,
        zorder=3,
    )
    for x, unresolved_share in zip(
        x_positions,
        backgrounds["unresolved_or_ambiguous_share"],
        strict=True,
    ):
        ax.text(
            x,
            0,
            f"{float(unresolved_share):.0%}",
            ha="center",
            va="center",
            fontsize=5.8,
            color="white" if float(unresolved_share) < 0.5 else SLATE_900,
            zorder=4,
        )

    rank_cells = cells[cells["rank"].eq(rank)].copy()
    rank_cells["share"] = pd.to_numeric(rank_cells["share"], errors="coerce").fillna(0.0)
    rank_cells["count"] = pd.to_numeric(rank_cells["count"], errors="coerce").fillna(0.0)
    changed_keys = set(
        zip(
            differences.loc[differences["rank"].eq(rank), "taxon"],
            differences.loc[differences["rank"].eq(rank), "trait_slug"],
            strict=True,
        )
    )
    rank_cells["changed_pair"] = [
        (taxon, trait_slug) in changed_keys
        for taxon, trait_slug in zip(rank_cells["taxon"], rank_cells["trait_slug"], strict=True)
    ]
    present = rank_cells[
        rank_cells["count"].gt(0)
        & rank_cells["taxon"].isin(taxon_order)
        & rank_cells["changed_pair"]
    ]
    xs = [column_x[(row.trait_slug, row.direction)] for row in present.itertuples(index=False)]
    ys = [taxon_y[row.taxon] for row in present.itertuples(index=False)]
    sizes = [28.0 + 560.0 * float(share_norm(row.share)) for row in present.itertuples(index=False)]
    colors = [DIRECTION_COLORS[row.direction] for row in present.itertuples(index=False)]
    ax.scatter(xs, ys, s=sizes, c=colors, edgecolor="white", linewidth=0.8, alpha=0.88)

    ax.set_xlim(-0.5, len(columns) - 0.5)
    ax.set_ylim(len(taxa) + 0.5, -0.65)
    ax.set_xticks(
        x_positions,
        [
            f"{direction}\n" + "\n".join(wrap(base_trait_label(label), width=12))
            for trait_slug, label, direction in columns
        ],
        fontsize=7.5,
    )
    ax.set_yticks(
        np.arange(len(taxa) + 1),
        ["direction background\n(grey = unresolved)", *taxa["display_name"].tolist()],
        fontsize=7.5,
    )
    ax.tick_params(length=0, pad=4)
    ax.set_xticks(np.arange(-0.5, len(columns), 1), minor=True)
    ax.set_yticks(np.arange(0.5, len(taxa) + 1, 1), minor=True)
    ax.grid(which="minor", color=GRID_COLOR, linewidth=0.7)
    ax.tick_params(which="minor", length=0)
    ax.set_title(f"{rank.title()} contributors", loc="left", fontweight="bold")
    for spine in ax.spines.values():
        spine.set_visible(False)

    specificity_cmap = plt.get_cmap("RdBu_r").copy()
    specificity_cmap.set_bad("#f3f4f6")
    taxon_specificity_values = (
        taxa["pd_minus_other_feces_prevalence_difference"]
        .replace([np.inf, -np.inf], np.nan)
        .clip(
            -SPECIFICITY_PREVALENCE_DIFFERENCE_LIMIT,
            SPECIFICITY_PREVALENCE_DIFFERENCE_LIMIT,
        )
        .to_numpy(dtype=float, copy=True)
    )
    total_specificity_profiles = taxa["pd_feces_profiles"].astype(float) + taxa[
        "other_feces_profiles"
    ].astype(float)
    taxon_specificity_values[total_specificity_profiles.eq(0).to_numpy()] = np.nan
    specificity_values = np.concatenate(([np.nan], taxon_specificity_values))
    specificity_ax.imshow(
        np.ma.masked_invalid(specificity_values[:, None]),
        cmap=specificity_cmap,
        vmin=-SPECIFICITY_PREVALENCE_DIFFERENCE_LIMIT,
        vmax=SPECIFICITY_PREVALENCE_DIFFERENCE_LIMIT,
        aspect="auto",
    )
    specificity_ax.set_xticks(
        [0],
        ["PD/other counts\nand prevalence difference"],
        fontsize=8,
    )
    specificity_ax.tick_params(left=False, labelleft=False, length=0, pad=4)
    specificity_ax.set_xticks([-0.5, 0.5], minor=True)
    specificity_ax.set_yticks(np.arange(-0.5, len(taxa) + 1, 1), minor=True)
    specificity_ax.grid(which="minor", color="white", linewidth=1)
    specificity_ax.tick_params(which="minor", length=0)
    for y, row in enumerate(taxa.itertuples(index=False), start=1):
        value = specificity_values[y]
        difference_pp = 100.0 * float(row.pd_minus_other_feces_prevalence_difference)
        difference_label = "0.0 pp" if abs(difference_pp) < 0.05 else f"{difference_pp:+.1f} pp"
        label = f"{row.pd_other_feces_display}\n{difference_label}"
        specificity_ax.text(
            0,
            y,
            label,
            ha="center",
            va="center",
            fontsize=6.8,
            color=(
                "white"
                if np.isfinite(value)
                and abs(float(value)) >= 0.72 * SPECIFICITY_PREVALENCE_DIFFERENCE_LIMIT
                else SLATE_900
            ),
        )
    for spine in specificity_ax.spines.values():
        spine.set_visible(False)


def plot_contributor_heatmap(
    cells: pd.DataFrame,
    differences: pd.DataFrame,
    selected: pd.DataFrame,
    coverage: pd.DataFrame,
    output_path: Path,
) -> None:
    """Plot changed contributor pairs with directional background and specificity context."""
    columns = trait_direction_columns(coverage)
    ranks = [
        rank
        for rank in sorted(selected["rank"].unique(), key=lambda value: RANK_ORDER.get(value, 99))
        if rank in RANK_ORDER
    ]
    max_share = max(float(pd.to_numeric(cells["share"], errors="coerce").max()), 0.01)
    share_norm = Normalize(vmin=0.0, vmax=max_share)
    height = 2.8 + sum(
        max(4.1, (len(selected[selected["rank"].eq(rank)]) + 1) * 0.34) for rank in ranks
    )
    fig = plt.figure(figsize=(18.5, height), facecolor="white")
    grid = fig.add_gridspec(
        len(ranks),
        1,
        left=0.22,
        right=0.96,
        top=0.83,
        bottom=0.7 / height,
        hspace=0.25,
        height_ratios=[
            max(4.1, (len(selected[selected["rank"].eq(rank)]) + 1) * 0.34) for rank in ranks
        ],
    )
    grids = [grid[index, 0] for index in range(len(ranks))]
    for rank, subgrid in zip(ranks, grids, strict=True):
        plot_rank_panel(
            fig,
            subgrid,
            rank,
            cells,
            differences,
            selected,
            coverage,
            columns,
            share_norm,
        )

    fig.text(
        0.03,
        0.992,
        "PD trait contributors with large direction differences",
        fontsize=16,
        fontweight="bold",
        va="top",
    )
    subtitle = (
        "Rows require at least one absolute increased-versus-decreased full-cell contributor-share "
        f"difference of {MIN_LARGE_DIRECTION_SHARE_DIFFERENCE:.0%}. Paired red/blue backgrounds "
        "show total propagated counts and unresolved mass; contributor-circle area scales with "
        "full-cell share. The right strip uses PD minus other-feces prevalence, not a raw ratio."
    )
    fig.text(
        0.03,
        0.955,
        "\n".join(wrap(subtitle, 180)),
        fontsize=9,
        color=SLATE_500,
        va="top",
    )
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
    legend_handles.append(
        Patch(facecolor=UNRESOLVED_COLOR, label="grey halo: unresolved/ambiguous background")
    )
    legend_handles.append(
        Patch(
            facecolor="white",
            edgecolor=SLATE_500,
            label=(
                "specificity: PD - other prevalence (clipped at +/-5 pp); "
                "raw x/0 counts display as ---"
            ),
        )
    )
    fig.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.63, 0.89),
        ncol=4,
        frameon=False,
        fontsize=8.2,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with plt.rc_context({"svg.hashsalt": STEM}):
        fig.savefig(output_path, dpi=240, bbox_inches="tight", facecolor="white")
        svg_path = output_path.with_suffix(".svg")
        fig.savefig(
            svg_path,
            bbox_inches="tight",
            facecolor="white",
            metadata={"Date": None},
        )
    strip_svg_trailing_whitespace(svg_path)
    plt.close(fig)


def build_trait_cocarriage_tables(
    cells: pd.DataFrame,
    selected: pd.DataFrame,
    coverage: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Build direction-specific trait sharing and retain changed pairs for display."""
    trait_order = coverage[["trait_slug", "trait_label"]].drop_duplicates().reset_index(drop=True)
    summary = selected.copy()
    summary["taxon_key"] = summary["rank"].astype(str) + "|" + summary["taxon"].astype(str)

    present = cells.copy()
    present["count"] = pd.to_numeric(present["count"], errors="coerce").fillna(0.0)
    present["taxon_key"] = present["rank"].astype(str) + "|" + present["taxon"].astype(str)
    present = present[present["count"].gt(0)]

    for trait in trait_order.itertuples(index=False):
        slug = str(trait.trait_slug)
        direction_columns = []
        for direction in DIRECTION_ORDER:
            keys = set(
                present.loc[
                    present["trait_slug"].eq(slug) & present["direction"].eq(direction),
                    "taxon_key",
                ]
            )
            column = f"{slug}_{direction}_present"
            summary[column] = summary["taxon_key"].isin(keys)
            direction_columns.append(column)
        summary[f"{slug}_present"] = summary[direction_columns].any(axis=1)
    trait_flags = [f"{row.trait_slug}_present" for row in trait_order.itertuples(index=False)]
    summary["present_trait_count"] = summary[trait_flags].sum(axis=1)

    edges: list[dict[str, object]] = []
    trait_rows = list(trait_order.itertuples(index=False))
    for first, second in itertools.combinations(trait_rows, 2):
        first_slug = str(first.trait_slug)
        second_slug = str(second.trait_slug)
        shared_by_direction: dict[str, pd.DataFrame] = {}
        for direction in DIRECTION_ORDER:
            shared_by_direction[direction] = summary[
                summary[f"{first_slug}_{direction}_present"]
                & summary[f"{second_slug}_{direction}_present"]
            ].copy()
        increased = shared_by_direction["inc"]
        decreased = shared_by_direction["dec"]
        increased_keys = set(increased["taxon_key"])
        decreased_keys = set(decreased["taxon_key"])
        increased_only_keys = increased_keys - decreased_keys
        decreased_only_keys = decreased_keys - increased_keys

        def display_names(keys: set[str]) -> str:
            changed = summary[summary["taxon_key"].isin(keys)].sort_values(
                ["present_trait_count", "total_count"], ascending=False
            )
            return "; ".join(changed["display_name"].astype(str))

        difference = len(increased) - len(decreased)
        edges.append(
            {
                "trait_a_slug": first_slug,
                "trait_a_label": base_trait_label(first.trait_label),
                "trait_b_slug": second_slug,
                "trait_b_label": base_trait_label(second.trait_label),
                "increased_shared_taxa": len(increased),
                "decreased_shared_taxa": len(decreased),
                "increased_minus_decreased_shared_taxa": difference,
                "absolute_direction_difference": abs(difference),
                "changes_across_directions": difference != 0,
                "increased_only_taxa": display_names(increased_only_keys),
                "decreased_only_taxa": display_names(decreased_only_keys),
                "species_increased_shared_taxa": int(increased["rank"].eq("species").sum()),
                "species_decreased_shared_taxa": int(decreased["rank"].eq("species").sum()),
                "family_increased_shared_taxa": int(increased["rank"].eq("family").sum()),
                "family_decreased_shared_taxa": int(decreased["rank"].eq("family").sum()),
            }
        )
    all_edges = pd.DataFrame.from_records(edges)
    changed_edges = all_edges[all_edges["changes_across_directions"]].copy()
    changed_edges = changed_edges.sort_values(
        ["absolute_direction_difference", "trait_a_label", "trait_b_label"],
        ascending=[False, True, True],
    ).reset_index(drop=True)
    return summary.drop(columns=["taxon_key"]), all_edges, changed_edges


def plot_trait_cocarriage_network(
    edges: pd.DataFrame,
    taxa_summary: pd.DataFrame,
    coverage: pd.DataFrame,
    output_path: Path,
) -> None:
    """Plot only trait pairs whose selected-taxon co-carriage changes by direction."""
    traits = coverage[["trait_slug", "trait_label"]].drop_duplicates().reset_index(drop=True)
    labels = {row.trait_slug: base_trait_label(row.trait_label) for row in traits.itertuples()}
    active_slugs = set(edges["trait_a_slug"]).union(edges["trait_b_slug"])
    slugs = [slug for slug in traits["trait_slug"].astype(str) if slug in active_slugs]
    angles = np.linspace(math.pi / 2, math.pi / 2 + 2 * math.pi, len(slugs), endpoint=False)
    positions = {
        slug: (math.cos(angle), math.sin(angle)) for slug, angle in zip(slugs, angles, strict=True)
    }

    fig, ax = plt.subplots(figsize=(9.2, 8.0), facecolor="white", constrained_layout=True)
    if edges.empty:
        ax.text(0, 0, "No trait-pair co-carriage counts change across directions.", ha="center")
    max_shared = (
        1
        if edges.empty
        else max(
            int(edges[["increased_shared_taxa", "decreased_shared_taxa"]].to_numpy().max()),
            1,
        )
    )
    for row in edges.itertuples(index=False):
        x1, y1 = positions[row.trait_a_slug]
        x2, y2 = positions[row.trait_b_slug]
        distance = max(math.hypot(x2 - x1, y2 - y1), 1e-9)
        normal_x = -(y2 - y1) / distance
        normal_y = (x2 - x1) / distance
        for direction, offset_sign, count in [
            ("inc", 1.0, int(row.increased_shared_taxa)),
            ("dec", -1.0, int(row.decreased_shared_taxa)),
        ]:
            offset = 0.027 * offset_sign
            width = 1.2 + 4.8 * count / max_shared
            ax.plot(
                [x1 + normal_x * offset, x2 + normal_x * offset],
                [y1 + normal_y * offset, y2 + normal_y * offset],
                color=DIRECTION_COLORS[direction],
                linewidth=width,
                alpha=0.78,
                solid_capstyle="round",
                zorder=1,
            )
            ax.text(
                (x1 + x2) / 2 + normal_x * 0.09 * offset_sign,
                (y1 + y2) / 2 + normal_y * 0.09 * offset_sign,
                str(count),
                ha="center",
                va="center",
                fontsize=8,
                color=DIRECTION_COLORS[direction],
                bbox={
                    "boxstyle": "round,pad=0.18",
                    "facecolor": "white",
                    "edgecolor": GRID_COLOR,
                },
                zorder=3,
            )

    for slug in slugs:
        x, y = positions[slug]
        inc_count = int(taxa_summary[f"{slug}_inc_present"].sum())
        dec_count = int(taxa_summary[f"{slug}_dec_present"].sum())
        ax.scatter(
            [x],
            [y],
            s=2200,
            color="white",
            edgecolor=SLATE_700,
            linewidth=1.5,
            zorder=4,
        )
        label = "\n".join(wrap(labels[slug], width=15))
        ax.text(
            x,
            y,
            f"{label}\n{inc_count} inc | {dec_count} dec",
            ha="center",
            va="center",
            fontsize=8.5,
            color=SLATE_900,
            zorder=5,
        )

    ax.set_title("Direction-changing trait co-carriage pairs", fontweight="bold")
    ax.text(
        0,
        -1.38,
        "Only pairs with different increased/decreased shared-taxon counts are shown.",
        ha="center",
        va="center",
        fontsize=9,
        color=SLATE_500,
    )
    ax.legend(
        handles=[
            Line2D([0], [0], color=DIRECTION_COLORS[direction], linewidth=4, label=label)
            for direction, label in DIRECTION_LABELS.items()
        ],
        loc="upper center",
        ncol=2,
        frameon=False,
    )
    ax.set_xlim(-1.55, 1.55)
    ax.set_ylim(-1.48, 1.42)
    ax.set_aspect("equal")
    ax.axis("off")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with plt.rc_context({"svg.hashsalt": f"{STEM}_cocarriage"}):
        fig.savefig(output_path, dpi=240, bbox_inches="tight", facecolor="white")
        svg_path = output_path.with_suffix(".svg")
        fig.savefig(
            svg_path,
            bbox_inches="tight",
            facecolor="white",
            metadata={"Date": None},
        )
    strip_svg_trailing_whitespace(svg_path)
    plt.close(fig)


def write_note(
    selected: pd.DataFrame,
    heatmap_taxa: pd.DataFrame,
    differences: pd.DataFrame,
    changed_edges: pd.DataFrame,
    all_edges: pd.DataFrame,
    coverage: pd.DataFrame,
    trait_intervals: pd.DataFrame,
    total_direction_pairs: int,
    output_path: Path,
) -> None:
    """Write interpretation, bias limits, and the change-only results."""
    observed = selected[
        (
            selected["pd_feces_profiles"].astype(float)
            + selected["other_feces_profiles"].astype(float)
        )
        .gt(0)
        .to_numpy()
    ]
    specific = observed.sort_values(
        "pd_minus_other_feces_prevalence_difference", ascending=False
    ).head(5)
    supported_patterns = trait_intervals[
        trait_intervals["overall_ci95_low"].gt(0) | trait_intervals["overall_ci95_high"].lt(0)
    ].copy()
    unresolved = coverage.sort_values("unresolved_or_ambiguous_share", ascending=False).head(6)
    lines = [
        "# PD Trait Contributor Dot Heatmap",
        "",
        f"![Contributor heatmap]({STEM}_species_family_dot_heatmap.png)",
        "",
        "The background row replaces the separate unresolved-share bars. Within each trait, "
        "paired red/blue circle area compares the increased/decreased propagated count; the grey "
        "outer area is the portion unresolved or ambiguous at the displayed rank. Contributor "
        "circle area uses the full cell denominator, so unresolved mass already reduces the "
        "colored contributor share rather than being discarded.",
        "",
        (
            f"The heatmap retains `{len(heatmap_taxa)}` of `{len(selected)}` candidate taxa. A "
            "row must have at least one absolute increased-versus-decreased full-cell contributor "
            f"share difference of `{MIN_LARGE_DIRECTION_SHARE_DIFFERENCE:.0%}`. This removes rows "
            "containing only small directional changes."
        ),
        "",
        f"The direction-difference table retains {len(differences)} of "
        f"{total_direction_pairs} taxon-trait pairs whose full-cell contributor share changes "
        "between directions. Equal pairs, including pairs absent in both directions, are omitted.",
        "",
        "The right strip colors taxa by the profile-pool-normalized prevalence difference "
        "(`PD prevalence - other-feces prevalence`) in percentage points. This finite difference "
        "is used instead of a count ratio, so a small or zero other-feces count cannot create an "
        "infinite ratio. Raw count displays with a zero denominator, such as `1/0`, are rendered "
        "as `---`. The metric remains descriptive and does not account for study clustering.",
        "",
        "## Direction-changing co-carriage pairs",
        "",
        f"![Trait co-carriage network]({STEM}_trait_cocarriage_network.png)",
        "",
        f"Only {len(changed_edges)} of {len(all_edges)} trait pairs have different selected-taxon "
        "co-carriage counts between directions and are displayed:",
        "",
        "| trait pair | increased | decreased | inc - dec | increased-only taxa | "
        "decreased-only taxa |",
        "| --- | ---: | ---: | ---: | --- | --- |",
    ]
    for row in changed_edges.itertuples(index=False):
        lines.append(
            f"| {row.trait_a_label} + {row.trait_b_label} | "
            f"{int(row.increased_shared_taxa)} | {int(row.decreased_shared_taxa)} | "
            f"{int(row.increased_minus_decreased_shared_taxa):+d} | "
            f"{row.increased_only_taxa or 'none'} | {row.decreased_only_taxa or 'none'} |"
        )
    lines.extend(
        [
            "",
            "## Aggregate direction patterns",
            "",
            "The contributor counts are not valid independent observations for a significance "
            "test. As external context, the existing paired, study-weighted aggregate analysis "
            "has 95% intervals excluding zero for these traits:",
            "",
            "| trait | mean increased - decreased z | 95% interval | pattern |",
            "| --- | ---: | ---: | --- |",
        ]
    )
    for row in supported_patterns.itertuples(index=False):
        lines.append(
            f"| {row.trait} | {float(row.overall_mean_z_contrast):.3f} | "
            f"[{float(row.overall_ci95_low):.3f}, {float(row.overall_ci95_high):.3f}] | "
            f"{row.overall_direction} |"
        )
    lines.extend(
        [
            "",
            "These intervals are aggregate directional context, not contributor-level tests, and "
            "the table does not add a new multiple-testing correction.",
            "",
            "## Largest unresolved or ambiguous shares",
            "",
            "| rank | trait | direction | unresolved share |",
            "| --- | --- | --- | ---: |",
        ]
    )
    for row in unresolved.itertuples(index=False):
        lines.append(
            f"| {row.rank} | {base_trait_label(row.trait_label)} | {row.direction} | "
            f"{float(row.unresolved_or_ambiguous_share):.1%} |"
        )
    lines.extend(
        [
            "",
            "High unresolved shares limit rank-specific attribution. They do not create extra "
            "colored contributor mass because every contributor share is divided by the full "
            "trait-direction count.",
            "",
            "## PD-versus-other-feces occurrence context",
            "",
            "Top profile-pool-normalized reported taxon labels (descriptive):",
            "",
            "| rank | taxon | PD/other profiles | PD/other prevalence | prevalence difference |",
            "| --- | --- | ---: | ---: | ---: |",
        ]
    )
    for row in specific.itertuples(index=False):
        lines.append(
            f"| {row.rank} | {row.display_name} | {row.pd_other_feces_display} | "
            f"{float(row.pd_feces_prevalence):.2%}/{float(row.other_feces_prevalence):.2%} | "
            f"{100.0 * float(row.pd_minus_other_feces_prevalence_difference):+.2f} pp |"
        )
    lines.extend(
        [
            "",
            "## Counting and selection bias",
            "",
            "- Counts are propagated taxon-trait assertions, not independent patients, samples, "
            "or original source assertions. Repeated experiments and descendant propagation can "
            "increase them.",
            "- The displayed taxa are the union of the top five resolved contributors per "
            "trait-direction cell. Broadly annotated or frequently reported taxa therefore have "
            "more opportunities to enter and remain in the selected set.",
            "- Species and family panels are separate resolution views. Parent and descendant "
            "labels must not be added as if they were disjoint contributors.",
            "- Co-carriage means that the selected label contributes to both annotation cells; it "
            "does not establish within-sample expression, mechanism, or causal PD specificity.",
            "",
            "Audit files:",
            "",
            f"- `{STEM}_cells.tsv`: long-form source dot cells.",
            f"- `{STEM}_direction_differences.tsv`: change-only taxon-trait direction pairs.",
            f"- `{STEM}_heatmap_taxa.tsv`: taxa passing the large-difference row filter.",
            f"- `{STEM}_selected_taxa.tsv`: selected taxa with normalized specificity context.",
            f"- `{STEM}_taxon_specificity.tsv`: PD-vs-other-feces counts and prevalence metrics.",
            f"- `{STEM}_trait_cocarriage_taxa.tsv`: direction-specific taxon-by-trait presence.",
            f"- `{STEM}_trait_cocarriage_all_edges.tsv`: all pairwise direction counts.",
            f"- `{STEM}_trait_cocarriage_edges.tsv`: only pairs that change by direction.",
            f"- `{STEM}_species_family_dot_heatmap.svg`: vector contributor heatmap.",
            f"- `{STEM}_trait_cocarriage_network.svg`: vector change-only network.",
        ]
    )
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    """Refresh contributor artifacts from existing contributor tables."""
    metadata = load_profile_metadata()
    require_current_trait_taxa_evidence(metadata, OUTPUT_DIR)
    # Every path is built from OUTPUT_DIR here, so redirecting it redirects every write.
    selected_taxa_path = OUTPUT_DIR / f"{STEM}_selected_taxa.tsv"
    cells = read_tsv(OUTPUT_DIR / f"{STEM}_cells.tsv")
    selected = read_tsv(selected_taxa_path)
    coverage = read_tsv(OUTPUT_DIR / f"{STEM}_rank_coverage.tsv")
    trait_intervals = read_tsv(OUTPUT_DIR / "pd_species_six_traits_influence_heatmap_traits.tsv")
    selected = augment_selected_taxa(selected, metadata)
    write_tsv(selected, selected_taxa_path)
    differences = build_direction_differences(cells, coverage)
    write_tsv(differences, OUTPUT_DIR / f"{STEM}_direction_differences.tsv")
    heatmap_taxa = select_large_difference_taxa(selected, differences)
    write_tsv(heatmap_taxa, OUTPUT_DIR / f"{STEM}_heatmap_taxa.tsv")
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
        "pd_feces_total_profiles",
        "other_feces_total_profiles",
        "pd_feces_prevalence",
        "other_feces_prevalence",
        "pd_vs_other_feces_log2_prevalence_ratio",
        "pd_minus_other_feces_prevalence_difference",
        "pd_other_feces_display",
    ]
    write_tsv(selected[specificity_columns], OUTPUT_DIR / f"{STEM}_taxon_specificity.tsv")
    plot_contributor_heatmap(
        cells,
        differences,
        heatmap_taxa,
        coverage,
        OUTPUT_DIR / f"{STEM}_species_family_dot_heatmap.png",
    )
    cocarriage_taxa, all_cocarriage_edges, changed_cocarriage_edges = build_trait_cocarriage_tables(
        cells, selected, coverage
    )
    write_tsv(cocarriage_taxa, OUTPUT_DIR / f"{STEM}_trait_cocarriage_taxa.tsv")
    write_tsv(all_cocarriage_edges, OUTPUT_DIR / f"{STEM}_trait_cocarriage_all_edges.tsv")
    write_tsv(
        changed_cocarriage_edges,
        OUTPUT_DIR / f"{STEM}_trait_cocarriage_edges.tsv",
    )
    plot_trait_cocarriage_network(
        changed_cocarriage_edges,
        cocarriage_taxa,
        coverage,
        OUTPUT_DIR / f"{STEM}_trait_cocarriage_network.png",
    )
    write_note(
        selected,
        heatmap_taxa,
        differences,
        changed_cocarriage_edges,
        all_cocarriage_edges,
        coverage,
        trait_intervals,
        cells.groupby(["rank", "taxon", "trait_slug"]).ngroups,
        OUTPUT_DIR / f"{STEM}.md",
    )
    print(
        f"Wrote {len(selected)} selected taxa, {len(heatmap_taxa)} large-difference heatmap taxa, "
        f"{len(differences)} direction differences, "
        f"and {len(changed_cocarriage_edges)}/{len(all_cocarriage_edges)} changed trait pairs"
    )
    print(OUTPUT_DIR / f"{STEM}_species_family_dot_heatmap.png")
    print(OUTPUT_DIR / f"{STEM}_trait_cocarriage_network.png")


if __name__ == "__main__":
    main()
