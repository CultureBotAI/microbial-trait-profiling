"""Create PD ecological-reporter specificity and condition-taxon-trait figures."""

from __future__ import annotations

import math
from pathlib import Path
from textwrap import wrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import FancyBboxPatch, PathPatch, Rectangle
from matplotlib.path import Path as MplPath
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[1]
VIS_DIR = REPO_ROOT / "trait-visualization"
INVESTIGATION_DIR = VIS_DIR / "trait-profile-investigation"
SIX_SHARED_DIR = INVESTIGATION_DIR / "six-shared-traits"
OUTPUT_DIR = INVESTIGATION_DIR / "PD-cellulose"

PROFILE_MATRIX_PATH = REPO_ROOT / "data/derived/study_trait_profiles_wide.tsv"
PROFILE_METADATA_PATH = VIS_DIR / "study_trait_profile_study_metadata.tsv"
TAXA_PATH = REPO_ROOT / "data/derived/taxa.tsv"
TAXA_TRAIT_MATRIX_PATH = REPO_ROOT / "data/derived/taxa_trait_matrix.tsv"

PD_CONDITION = "Parkinson's disease"
FECES = "Feces"
INCREASED = "increased"
DECREASED = "decreased"
DIRECTION_LABELS = {"inc": INCREASED, "dec": DECREASED}

CELLULOSE = "cellulose catabolic process"
POSITIVE_COMPOSITE = "three-trait increased-side composite"
REPORTER_COMPOSITE = "ecological-reporter composite"
POSITIVE_COMPOSITE_TRAITS = [
    "moderately halophilic",
    "facultatively aerobic",
    "arsenate reduction",
]
FIGURE_TRAIT_ORDER = [
    *POSITIVE_COMPOSITE_TRAITS,
    "coccus shaped",
    "sphere shaped",
    CELLULOSE,
]

COMBINED_CONDITIONS = [
    PD_CONDITION,
    "Chronic constipation",
    "Antimicrobial agent",
    "Diet",
    "Irritable bowel syndrome",
    "Colorectal cancer",
    "Type II diabetes mellitus",
    "Stroke",
    "COVID-19",
]
NETWORK_CONDITIONS = [
    PD_CONDITION,
    "Chronic constipation",
    "Antimicrobial agent",
    "Diet",
    "Colorectal cancer",
]
ALLOWED_TAXON_RANKS = {"genus", "species"}

TRAIT_COLORS = {
    "facultatively aerobic": "#0072b2",
    CELLULOSE: "#d55e00",
    "coccus shaped": "#009e73",
    "sphere shaped": "#cc79a7",
    "arsenate reduction": "#e69f00",
    "moderately halophilic": "#6a3d9a",
}
CONDITION_COLORS = {
    PD_CONDITION: "#c2410c",
    "Chronic constipation": "#7c3aed",
    "Antimicrobial agent": "#dc2626",
    "Diet": "#16a34a",
    "Irritable bowel syndrome": "#0891b2",
    "Colorectal cancer": "#9333ea",
    "Type II diabetes mellitus": "#ca8a04",
    "Stroke": "#475569",
    "COVID-19": "#0f766e",
}
INCREASED_COLOR = "#c2410c"
DECREASED_COLOR = "#2563eb"
SLATE_50 = "#f8fafc"
SLATE_100 = "#f1f5f9"
SLATE_200 = "#e2e8f0"
SLATE_300 = "#cbd5e1"
SLATE_500 = "#64748b"
SLATE_700 = "#334155"
SLATE_900 = "#0f172a"

COMBINED_STEM = "pd_ecological_reporter_combined"
NETWORK_STEM = "pd_condition_taxon_trait_network"


def read_tsv(path: Path, **kwargs: object) -> pd.DataFrame:
    """Read a TSV, accepting a gzip-compressed sibling when the plain file is absent."""
    source = path
    if not source.exists():
        compressed = path.with_suffix(f"{path.suffix}.gz")
        if not compressed.exists():
            raise FileNotFoundError(path)
        source = compressed
    return pd.read_csv(source, sep="\t", **kwargs)


def save_figure(fig: plt.Figure, output_path: Path) -> None:
    """Save a figure as high-resolution PNG and vector SVG."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(output_path.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_tsv(frame: pd.DataFrame, output_path: Path) -> None:
    """Write an auditable TSV with stable missing-value output."""
    frame.to_csv(output_path, sep="\t", index=False, na_rep="NA", float_format="%.8g")


def load_trait_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load metadata and recreate the globally normalized six-trait profile z-scores."""
    traits = read_tsv(SIX_SHARED_DIR / "six_shared_trait_summary.tsv")
    trait_ids = traits["trait_id"].tolist()
    metadata = read_tsv(PROFILE_METADATA_PATH, index_col=0)
    metadata = metadata.copy()
    metadata.index = metadata.index.astype(str)
    metadata["profile_id"] = metadata.index
    for column in [
        "study_id",
        "experiment_id",
        "inc_dec",
        "condition",
        "body_site_group",
        "condition_group",
    ]:
        metadata[column] = metadata[column].fillna("").astype(str).str.strip()
    metadata = metadata[metadata["inc_dec"].isin(DIRECTION_LABELS)].copy()

    profile_counts = read_tsv(
        PROFILE_MATRIX_PATH,
        usecols=["study_id", "experiment_id", "inc_dec", *trait_ids],
    )
    profile_counts["profile_id"] = (
        profile_counts["study_id"].astype(str)
        + "|"
        + profile_counts["experiment_id"].astype(str)
        + "|"
        + profile_counts["inc_dec"].astype(str)
    )
    values = profile_counts.groupby("profile_id", sort=False)[trait_ids].sum()
    missing_profiles = metadata.index.difference(values.index)
    if len(missing_profiles):
        raise ValueError(f"Profile matrix is missing {len(missing_profiles)} metadata profiles")
    values = values.loc[metadata.index]
    totals = pd.to_numeric(metadata["total_trait_signal"], errors="coerce").replace(0, pd.NA)
    composition = values.div(totals, axis=0).fillna(0.0)
    scaled = (composition * 1000).map(math.log1p)
    zscores = (scaled - scaled.mean()) / scaled.std(ddof=0).replace(0, pd.NA)
    return metadata, zscores.fillna(0.0), traits


def _confidence_interval(values: pd.Series) -> tuple[float, float]:
    """Return a descriptive two-sided 95% t interval for study-level contrasts."""
    clean = pd.to_numeric(values, errors="coerce").dropna()
    if len(clean) < 2:
        return math.nan, math.nan
    standard_error = float(clean.std(ddof=1) / math.sqrt(len(clean)))
    critical = float(stats.t.ppf(0.975, len(clean) - 1))
    mean = float(clean.mean())
    return mean - critical * standard_error, mean + critical * standard_error


def build_condition_trait_summaries(
    metadata: pd.DataFrame,
    zscores: pd.DataFrame,
    traits: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build paired within-study direction contrasts and condition-level summaries."""
    trait_ids = traits["trait_id"].tolist()
    trait_names = dict(zip(traits["trait_id"], traits["trait_label"], strict=True))
    trait_ids_by_name = {name: trait_id for trait_id, name in trait_names.items()}
    missing_composite = (set(POSITIVE_COMPOSITE_TRAITS) | {CELLULOSE}) - set(trait_ids_by_name)
    if missing_composite:
        raise ValueError(f"Missing composite traits: {sorted(missing_composite)}")

    profiles = zscores[trait_ids].join(
        metadata[
            [
                "study_id",
                "experiment_id",
                "inc_dec",
                "condition",
                "condition_group",
                "body_site_group",
            ]
        ]
    )
    profiles = profiles[profiles["body_site_group"].eq(FECES)].copy()
    records: list[dict[str, object]] = []
    for (condition, study_id), group in profiles.groupby(["condition", "study_id"], sort=False):
        direction_means = group.groupby("inc_dec")[trait_ids].mean()
        if not {"inc", "dec"}.issubset(direction_means.index):
            continue
        deltas: dict[str, float] = {}
        for trait_id in trait_ids:
            delta = float(
                direction_means.loc["inc", trait_id] - direction_means.loc["dec", trait_id]
            )
            deltas[trait_names[trait_id]] = delta
            records.append(
                {
                    "condition": condition,
                    "condition_group": group["condition_group"].iloc[0],
                    "study_id": study_id,
                    "trait_id": trait_id,
                    "trait": trait_names[trait_id],
                    "contrast": delta,
                    "contrast_definition": "mean z increased - mean z decreased",
                    "increased_profile_count": int(group["inc_dec"].eq("inc").sum()),
                    "decreased_profile_count": int(group["inc_dec"].eq("dec").sum()),
                }
            )
        positive_composite = float(np.mean([deltas[name] for name in POSITIVE_COMPOSITE_TRAITS]))
        reporter_composite = positive_composite - deltas[CELLULOSE]
        for trait, contrast, definition in [
            (
                POSITIVE_COMPOSITE,
                positive_composite,
                "mean increased-minus-decreased z for three increased-side traits",
            ),
            (
                REPORTER_COMPOSITE,
                reporter_composite,
                "three-trait composite minus cellulose increased-minus-decreased z",
            ),
        ]:
            records.append(
                {
                    "condition": condition,
                    "condition_group": group["condition_group"].iloc[0],
                    "study_id": study_id,
                    "trait_id": "composite",
                    "trait": trait,
                    "contrast": contrast,
                    "contrast_definition": definition,
                    "increased_profile_count": int(group["inc_dec"].eq("inc").sum()),
                    "decreased_profile_count": int(group["inc_dec"].eq("dec").sum()),
                }
            )
    study_contrasts = pd.DataFrame.from_records(records)
    if study_contrasts.empty:
        raise ValueError("No fecal studies contain both increased and decreased profiles")

    summary_records: list[dict[str, object]] = []
    for (condition, trait), group in study_contrasts.groupby(["condition", "trait"]):
        values = group["contrast"]
        ci_low, ci_high = _confidence_interval(values)
        positive_count = int(values.gt(0).sum())
        negative_count = int(values.lt(0).sum())
        summary_records.append(
            {
                "condition": condition,
                "condition_group": group["condition_group"].iloc[0],
                "trait": trait,
                "study_weighted_mean_contrast": float(values.mean()),
                "study_weighted_median_contrast": float(values.median()),
                "ci_low": ci_low,
                "ci_high": ci_high,
                "ci_excludes_zero": bool(ci_low > 0 or ci_high < 0),
                "paired_study_count": int(group["study_id"].nunique()),
                "positive_study_count": positive_count,
                "negative_study_count": negative_count,
                "zero_study_count": int(values.eq(0).sum()),
                "majority_direction_fraction": max(positive_count, negative_count) / len(group),
                "increased_profile_count": int(group["increased_profile_count"].sum()),
                "decreased_profile_count": int(group["decreased_profile_count"].sum()),
                "contrast_definition": group["contrast_definition"].iloc[0],
            }
        )
    summary = pd.DataFrame.from_records(summary_records)
    return summary, study_contrasts


def prepare_pd_taxon_panel(
    carriers: pd.DataFrame,
    trait_order: list[str],
    max_taxa: int = 11,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Select supported genus/species labels and create a complete taxon-trait grid."""
    selected = carriers[
        carriers["trait"].isin(trait_order) & carriers["taxonomic_rank"].isin(ALLOWED_TAXON_RANKS)
    ].copy()
    if selected.duplicated(["taxon_id", "trait"]).any():
        raise ValueError("Duplicate taxon/trait rows in PD carrier evidence")
    selected["total_studies"] = selected["increased_studies"] + selected["decreased_studies"]
    ranking = (
        selected.groupby(["taxon_id", "taxon_name", "taxonomic_rank"], as_index=False)
        .agg(
            annotated_trait_count=("trait", "nunique"),
            total_studies=("total_studies", "max"),
            increased_studies=("increased_studies", "max"),
            decreased_studies=("decreased_studies", "max"),
        )
        .assign(
            selection_score=lambda frame: (
                frame["total_studies"] * np.sqrt(frame["annotated_trait_count"])
            )
        )
        .sort_values(
            ["selection_score", "total_studies", "taxon_name"],
            ascending=[False, False, True],
        )
        .head(max_taxa)
        .reset_index(drop=True)
    )
    grid = pd.MultiIndex.from_product(
        [ranking["taxon_id"], trait_order], names=["taxon_id", "trait"]
    )
    cells = selected.set_index(["taxon_id", "trait"]).reindex(grid).reset_index()
    lookup = ranking.set_index("taxon_id")
    for column in ["taxon_name", "taxonomic_rank"]:
        cells[column] = cells["taxon_id"].map(lookup[column])
    cells["trait_annotated"] = cells["annotation_count"].notna()
    cells["direction_log2_ratio"] = np.nan
    annotated = cells["trait_annotated"]
    cells.loc[annotated, "direction_log2_ratio"] = np.log2(
        (cells.loc[annotated, "increased_studies"].astype(float) + 0.5)
        / (cells.loc[annotated, "decreased_studies"].astype(float) + 0.5)
    )
    cells["inc_dec_study_counts"] = ""
    cells.loc[annotated, "inc_dec_study_counts"] = (
        cells.loc[annotated, "increased_studies"].astype(int).astype(str)
        + "/"
        + cells.loc[annotated, "decreased_studies"].astype(int).astype(str)
    )
    cells["blank_reason"] = np.where(
        cells["trait_annotated"], "", "no annotation in local taxa-trait matrix"
    )
    return cells, ranking


def _selected_condition_order(summary: pd.DataFrame) -> list[str]:
    available = set(summary["condition"])
    missing = [condition for condition in COMBINED_CONDITIONS if condition not in available]
    if missing:
        raise ValueError(f"Combined-figure conditions missing from summary: {missing}")
    return COMBINED_CONDITIONS


def _plot_condition_heatmap(
    ax: plt.Axes,
    summary: pd.DataFrame,
    trait_order: list[str],
) -> None:
    """Draw the condition-by-trait study-weighted contrast heatmap."""
    conditions = _selected_condition_order(summary)
    matrix = summary.pivot(
        index="condition", columns="trait", values="study_weighted_mean_contrast"
    )
    study_counts = summary.pivot(index="condition", columns="trait", values="paired_study_count")
    significant = summary.pivot(index="condition", columns="trait", values="ci_excludes_zero")
    matrix = matrix.reindex(index=conditions, columns=trait_order)
    study_counts = study_counts.reindex(index=conditions, columns=trait_order)
    significant = significant.reindex(index=conditions, columns=trait_order)
    color_limit = max(0.75, float(np.nanquantile(np.abs(matrix.to_numpy(float)), 0.95)))
    norm = TwoSlopeNorm(vmin=-color_limit, vcenter=0, vmax=color_limit)
    image = ax.imshow(matrix.to_numpy(float), cmap="RdBu_r", norm=norm, aspect="auto")
    ax.set_xticks(
        range(len(trait_order)),
        [
            "\n".join(wrap(name.replace(" ecological-", "\necological-"), 17))
            for name in trait_order
        ],
        fontsize=8,
    )
    ax.set_yticks(range(len(conditions)), conditions, fontsize=9)
    ax.tick_params(length=0, pad=6)
    ax.set_xticks(np.arange(-0.5, len(trait_order), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(conditions), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.4)
    ax.tick_params(which="minor", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    for y, condition in enumerate(conditions):
        for x, trait in enumerate(trait_order):
            value = matrix.loc[condition, trait]
            if pd.isna(value):
                continue
            count = int(study_counts.loc[condition, trait])
            text_color = "white" if abs(value) >= color_limit * 0.52 else SLATE_900
            ax.text(
                x,
                y,
                f"{value:+.2f}\nn={count}",
                ha="center",
                va="center",
                fontsize=7.2,
                color=text_color,
                linespacing=1.05,
            )
            if bool(significant.loc[condition, trait]):
                ax.add_patch(
                    Rectangle(
                        (x - 0.47, y - 0.47),
                        0.94,
                        0.94,
                        fill=False,
                        edgecolor=SLATE_900,
                        linewidth=1.25,
                    )
                )
    pd_y = conditions.index(PD_CONDITION)
    ax.add_patch(
        Rectangle(
            (-0.54, pd_y - 0.54),
            len(trait_order) + 0.08,
            1.08,
            fill=False,
            edgecolor=INCREASED_COLOR,
            linewidth=2.2,
        )
    )
    bar = ax.figure.colorbar(image, ax=ax, pad=0.015, fraction=0.035)
    bar.set_label("Study-weighted mean delta z\n(increased - decreased)", fontsize=8)
    bar.ax.tick_params(labelsize=7)
    ax.set_title(
        "A  Disease fingerprints across the six shared traits",
        loc="left",
        fontweight="bold",
        fontsize=12,
        y=1.105,
        pad=0,
    )
    ax.text(
        0,
        1.045,
        "Cell text: effect and paired-study count; dark border: descriptive 95% t "
        "interval excludes 0",
        transform=ax.transAxes,
        fontsize=7.5,
        color=SLATE_500,
    )


def _plot_specificity_quadrant(ax: plt.Axes, summary: pd.DataFrame) -> None:
    """Plot cellulose depletion against the independent three-trait positive composite."""
    table = summary.pivot(index="condition", columns="trait", values="study_weighted_mean_contrast")
    counts = summary.pivot(index="condition", columns="trait", values="paired_study_count")
    table = table.dropna(subset=[CELLULOSE, POSITIVE_COMPOSITE])
    table = table[
        (counts.loc[table.index, CELLULOSE] >= 3)
        & (counts.loc[table.index, POSITIVE_COMPOSITE] >= 3)
    ].copy()
    table["cellulose_decreased_skew"] = -table[CELLULOSE]
    table["positive_trait_composite"] = table[POSITIVE_COMPOSITE]
    table["paired_studies"] = counts.loc[table.index, CELLULOSE]
    sizes = 18 + np.sqrt(table["paired_studies"].astype(float)) * 12
    ax.scatter(
        table["cellulose_decreased_skew"],
        table["positive_trait_composite"],
        s=sizes,
        color=SLATE_300,
        edgecolor="white",
        linewidth=0.45,
        alpha=0.62,
        zorder=2,
    )
    label_offsets = {
        PD_CONDITION: (-6, -14),
        "Chronic constipation": (6, -18),
        "Antimicrobial agent": (-6, 8),
        "Diet": (-6, -14),
        "Irritable bowel syndrome": (6, 6),
        "Colorectal cancer": (6, 9),
        "Type II diabetes mellitus": (-6, 11),
        "Stroke": (-6, -14),
        "COVID-19": (6, 6),
    }
    for condition in COMBINED_CONDITIONS:
        if condition not in table.index:
            continue
        row = table.loc[condition]
        color = CONDITION_COLORS[condition]
        marker_size = float(32 + math.sqrt(row["paired_studies"]) * 15)
        ax.scatter(
            row["cellulose_decreased_skew"],
            row["positive_trait_composite"],
            s=marker_size,
            color=color,
            edgecolor=SLATE_900 if condition == PD_CONDITION else "white",
            linewidth=1.5 if condition == PD_CONDITION else 0.7,
            zorder=4,
        )
        dx, dy = label_offsets[condition]
        ax.annotate(
            condition,
            (row["cellulose_decreased_skew"], row["positive_trait_composite"]),
            xytext=(dx, dy),
            textcoords="offset points",
            ha="left" if dx > 0 else "right",
            fontsize=7.2,
            color=color,
            fontweight="bold" if condition == PD_CONDITION else "normal",
            zorder=5,
        )
    ax.axhline(0, color=SLATE_700, linewidth=0.8)
    ax.axvline(0, color=SLATE_700, linewidth=0.8)
    ax.grid(color=SLATE_200, linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)
    ax.set_xlabel("Cellulose decreased-profile skew\n-(increased - decreased mean z)", fontsize=9)
    ax.set_ylabel("Three-trait increased-side composite\nmean(increased - decreased z)", fontsize=9)
    ax.set_title(
        "B  Specificity among fecal condition labels",
        loc="left",
        fontweight="bold",
        fontsize=12,
        y=1.105,
        pad=0,
    )
    ax.text(
        0,
        1.045,
        f"{len(table)} conditions with ≥3 paired studies; upper-right matches the proposed pattern",
        transform=ax.transAxes,
        fontsize=7.5,
        color=SLATE_500,
    )
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=8)


def _plot_pd_taxon_panel(
    ax: plt.Axes,
    cells: pd.DataFrame,
    ranking: pd.DataFrame,
    trait_order: list[str],
) -> None:
    """Draw KG-Microbe annotation glyphs colored by BugSigDB study-direction support."""
    lookup = cells.set_index(["taxon_id", "trait"])
    values = pd.to_numeric(cells["direction_log2_ratio"], errors="coerce")
    limit = max(2.0, float(np.nanquantile(np.abs(values), 0.95)))
    norm = TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit)
    cmap = plt.get_cmap("RdBu_r")
    for y, row in ranking.iterrows():
        for x, trait in enumerate(trait_order):
            cell = lookup.loc[(row["taxon_id"], trait)]
            if not bool(cell["trait_annotated"]):
                ax.scatter(x, y, s=26, marker="x", color=SLATE_300, linewidth=0.8)
                continue
            value = float(cell["direction_log2_ratio"])
            total = int(cell["increased_studies"] + cell["decreased_studies"])
            ax.scatter(
                x,
                y,
                s=30 + total * 7,
                color=cmap(norm(value)),
                edgecolor="white",
                linewidth=0.7,
                zorder=3,
            )
    rank_labels = {"genus": "G", "species": "S"}
    labels = [
        f"{row.taxon_name} [{rank_labels[row.taxonomic_rank]}]  "
        f"{int(row.increased_studies)}/{int(row.decreased_studies)} studies"
        for row in ranking.itertuples()
    ]
    ax.set_xticks(
        range(len(trait_order)),
        ["\n".join(wrap(name, 17)) for name in trait_order],
        fontsize=8,
    )
    ax.set_yticks(range(len(ranking)), labels, fontsize=8)
    ax.set_xlim(-0.6, len(trait_order) - 0.4)
    ax.set_ylim(len(ranking) - 0.5, -0.5)
    ax.tick_params(length=0, pad=5)
    ax.set_xticks(np.arange(-0.5, len(trait_order), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(ranking), 1), minor=True)
    ax.grid(which="minor", color=SLATE_200, linewidth=0.7)
    ax.tick_params(which="minor", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    scalar = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    bar = ax.figure.colorbar(scalar, ax=ax, pad=0.015, fraction=0.025)
    bar.set_label("log2((increased studies + 0.5) /\n(decreased studies + 0.5))", fontsize=8)
    bar.ax.tick_params(labelsize=7)
    ax.set_title(
        "C  PD reported taxa: KG-Microbe annotations x BugSigDB direction support",
        loc="left",
        fontweight="bold",
        fontsize=12,
        y=1.1,
        pad=0,
    )
    ax.text(
        0,
        1.04,
        "Circle = annotated trait; size = distinct-study support; x = no local "
        "annotation; genus/species only",
        transform=ax.transAxes,
        fontsize=7.5,
        color=SLATE_500,
    )


def plot_combined_figure(
    summary: pd.DataFrame,
    taxon_cells: pd.DataFrame,
    taxon_ranking: pd.DataFrame,
    trait_order: list[str],
    output_path: Path,
) -> None:
    """Create the three-panel ecological-reporter specificity figure."""
    display_traits = [*trait_order, REPORTER_COMPOSITE]
    fig = plt.figure(figsize=(20, 17), facecolor="white")
    grid = fig.add_gridspec(
        2,
        2,
        height_ratios=[1.04, 0.96],
        width_ratios=[1.25, 1.0],
        hspace=0.38,
        wspace=0.32,
    )
    heatmap_ax = fig.add_subplot(grid[0, 0])
    quadrant_ax = fig.add_subplot(grid[0, 1])
    taxon_ax = fig.add_subplot(grid[1, :])
    _plot_condition_heatmap(heatmap_ax, summary, display_traits)
    _plot_specificity_quadrant(quadrant_ax, summary)
    _plot_pd_taxon_panel(taxon_ax, taxon_cells, taxon_ranking, trait_order)
    fig.suptitle(
        "Parkinson's disease ecological-reporter hypothesis: specificity and provenance",
        fontsize=20,
        fontweight="bold",
        x=0.04,
        ha="left",
        y=0.995,
    )
    fig.text(
        0.04,
        0.972,
        "BugSigDB differential-abundance profiles provide direction and study support; "
        "KG-Microbe provides taxon-trait annotations. Values are annotation-derived "
        "ecological readouts, not measured activity.",
        fontsize=10,
        color=SLATE_700,
        ha="left",
    )
    fig.subplots_adjust(top=0.86, bottom=0.055)
    save_figure(fig, output_path)


def build_condition_taxon_trait_network(
    metadata: pd.DataFrame,
    taxa: pd.DataFrame,
    taxa_trait_matrix: pd.DataFrame,
    traits: pd.DataFrame,
    conditions: list[str],
    min_context_studies: int = 2,
    top_taxa_per_context: int = 2,
    max_taxa: int = 16,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Build a sparse evidence network without combining overlapping taxonomic ranks."""
    trait_ids = traits["trait_id"].tolist()
    trait_names = dict(zip(traits["trait_id"], traits["trait_label"], strict=True))
    profiles = metadata[
        metadata["condition"].isin(conditions) & metadata["body_site_group"].eq(FECES)
    ].copy()
    profiles["direction"] = profiles["inc_dec"].map(DIRECTION_LABELS)
    missing_conditions = set(conditions) - set(profiles["condition"])
    if missing_conditions:
        raise ValueError(
            f"Network conditions absent from fecal metadata: {sorted(missing_conditions)}"
        )
    occurrences = taxa.merge(
        profiles[["study_id", "experiment_id", "profile_id", "condition", "direction"]],
        on=["study_id", "experiment_id", "direction"],
        how="inner",
        validate="many_to_one",
    ).drop_duplicates(["condition", "direction", "profile_id", "taxon_id"])
    matrix = taxa_trait_matrix[["taxon_id", "taxon_name", "taxonomic_rank", *trait_ids]].copy()
    observed = occurrences.merge(
        matrix,
        on="taxon_id",
        how="inner",
        suffixes=("_reported", ""),
        validate="many_to_one",
    )
    observed = observed[observed["taxonomic_rank"].isin(ALLOWED_TAXON_RANKS)].copy()
    annotation_values = observed[trait_ids].apply(pd.to_numeric, errors="coerce").fillna(0)
    observed = observed[annotation_values.gt(0).any(axis=1)].copy()
    if observed.empty:
        raise ValueError("No annotated genus/species occurrences for network conditions")

    context_support = (
        observed.groupby(
            ["condition", "direction", "taxon_id", "taxon_name", "taxonomic_rank"],
            as_index=False,
        )
        .agg(study_count=("study_id", "nunique"), profile_count=("profile_id", "nunique"))
        .loc[lambda frame: frame["study_count"].ge(min_context_studies)]
        .sort_values(
            ["condition", "direction", "study_count", "profile_count", "taxon_name"],
            ascending=[True, True, False, False, True],
        )
        .groupby(["condition", "direction"], group_keys=False)
        .head(top_taxa_per_context)
    )
    if context_support.empty:
        raise ValueError("No network context meets the minimum study threshold")
    selected_ids = set(context_support["taxon_id"])
    selected_observed = observed[observed["taxon_id"].isin(selected_ids)].copy()
    selected_observed["context_study_key"] = (
        selected_observed["condition"].astype(str)
        + "|"
        + selected_observed["direction"].astype(str)
        + "|"
        + selected_observed["study_id"].astype(str)
    )
    taxon_ranking = (
        selected_observed.groupby(["taxon_id", "taxon_name", "taxonomic_rank"], as_index=False)
        .agg(
            context_study_support=("context_study_key", "nunique"),
            condition_count=("condition", "nunique"),
        )
        .sort_values(
            ["condition_count", "context_study_support", "taxon_name"],
            ascending=[False, False, True],
        )
        .head(max_taxa)
    )
    selected_ids = set(taxon_ranking["taxon_id"])
    context_support = context_support[context_support["taxon_id"].isin(selected_ids)].copy()
    selected_observed = selected_observed[selected_observed["taxon_id"].isin(selected_ids)].copy()

    nodes: list[dict[str, object]] = []
    edges: list[dict[str, object]] = []
    for condition in conditions:
        nodes.append(
            {
                "node_id": f"condition::{condition}",
                "label": condition,
                "layer": "condition",
                "node_type": "BugSigDB condition",
                "taxonomic_rank": "",
            }
        )
        for direction in [INCREASED, DECREASED]:
            context_id = f"context::{condition}::{direction}"
            context_rows = context_support[
                context_support["condition"].eq(condition)
                & context_support["direction"].eq(direction)
            ]
            if context_rows.empty:
                continue
            nodes.append(
                {
                    "node_id": context_id,
                    "label": f"{'↑' if direction == INCREASED else '↓'} {direction}",
                    "layer": "direction",
                    "node_type": "BugSigDB signature direction",
                    "taxonomic_rank": "",
                }
            )
            represented = selected_observed[
                selected_observed["condition"].eq(condition)
                & selected_observed["direction"].eq(direction)
                & selected_observed["taxon_id"].isin(context_rows["taxon_id"])
            ]
            edges.append(
                {
                    "source": f"condition::{condition}",
                    "target": context_id,
                    "source_layer": "condition",
                    "target_layer": "direction",
                    "study_support": int(represented["study_id"].nunique()),
                    "profile_support": int(represented["profile_id"].nunique()),
                    "annotation_count": pd.NA,
                    "evidence_type": "BugSigDB condition-direction support",
                }
            )
            for row in context_rows.itertuples():
                edges.append(
                    {
                        "source": context_id,
                        "target": f"taxon::{row.taxon_id}",
                        "source_layer": "direction",
                        "target_layer": "taxon",
                        "study_support": int(row.study_count),
                        "profile_support": int(row.profile_count),
                        "annotation_count": pd.NA,
                        "evidence_type": "BugSigDB exact-taxon report",
                    }
                )
    for row in taxon_ranking.itertuples():
        nodes.append(
            {
                "node_id": f"taxon::{row.taxon_id}",
                "label": row.taxon_name,
                "layer": "taxon",
                "node_type": "BugSigDB reported taxon",
                "taxonomic_rank": row.taxonomic_rank,
            }
        )
        matrix_row = matrix[matrix["taxon_id"].eq(row.taxon_id)].iloc[0]
        taxon_support = selected_observed[selected_observed["taxon_id"].eq(row.taxon_id)][
            "context_study_key"
        ].nunique()
        for trait_id in trait_ids:
            annotation_count = int(pd.to_numeric(matrix_row[trait_id], errors="coerce") or 0)
            if annotation_count <= 0:
                continue
            edges.append(
                {
                    "source": f"taxon::{row.taxon_id}",
                    "target": f"trait::{trait_id}",
                    "source_layer": "taxon",
                    "target_layer": "trait",
                    "study_support": int(taxon_support),
                    "profile_support": pd.NA,
                    "annotation_count": annotation_count,
                    "evidence_type": "KG-Microbe taxon-trait annotation",
                }
            )
    for trait_id in trait_ids:
        nodes.append(
            {
                "node_id": f"trait::{trait_id}",
                "label": trait_names[trait_id],
                "layer": "trait",
                "node_type": "KG-Microbe trait",
                "taxonomic_rank": "",
            }
        )
    nodes_frame = pd.DataFrame.from_records(nodes).drop_duplicates("node_id")
    edges_frame = pd.DataFrame.from_records(edges)
    return nodes_frame, edges_frame, context_support


def _even_positions(items: list[str], top: float = 0.92, bottom: float = 0.08) -> dict[str, float]:
    if len(items) == 1:
        return {items[0]: (top + bottom) / 2}
    values = np.linspace(top, bottom, len(items))
    return dict(zip(items, values, strict=True))


def _network_positions(
    nodes: pd.DataFrame,
    edges: pd.DataFrame,
    conditions: list[str],
    trait_order: list[str],
) -> dict[str, tuple[float, float]]:
    """Place four network layers and order taxa by their weighted source barycenter."""
    positions: dict[str, tuple[float, float]] = {}
    condition_ids = [f"condition::{condition}" for condition in conditions]
    for node_id, y in _even_positions(condition_ids).items():
        positions[node_id] = (0.0, y)
    direction_nodes = nodes[nodes["layer"].eq("direction")]["node_id"].tolist()
    direction_nodes.sort(
        key=lambda node_id: (
            conditions.index(node_id.split("::")[1]),
            0 if node_id.endswith(INCREASED) else 1,
        )
    )
    for node_id, y in _even_positions(direction_nodes, 0.96, 0.04).items():
        positions[node_id] = (1.0, y)
    taxon_nodes = nodes[nodes["layer"].eq("taxon")]["node_id"].tolist()
    taxon_scores: list[tuple[str, float]] = []
    for node_id in taxon_nodes:
        incoming = edges[edges["target"].eq(node_id) & edges["source_layer"].eq("direction")]
        weights = incoming["study_support"].astype(float).clip(lower=1)
        score = float(
            np.average([positions[source][1] for source in incoming["source"]], weights=weights)
        )
        taxon_scores.append((node_id, score))
    taxon_scores.sort(key=lambda item: item[1], reverse=True)
    for node_id, y in _even_positions([item[0] for item in taxon_scores], 0.96, 0.04).items():
        positions[node_id] = (2.0, y)
    trait_lookup = nodes[nodes["layer"].eq("trait")].set_index("label")["node_id"]
    trait_ids = [trait_lookup[trait] for trait in trait_order]
    for node_id, y in _even_positions(trait_ids, 0.9, 0.1).items():
        positions[node_id] = (3.0, y)
    return positions


def _network_edge_color(row: object, nodes: pd.DataFrame) -> str:
    if row.source_layer == "direction":
        return INCREASED_COLOR if INCREASED in row.source else DECREASED_COLOR
    if row.source_layer == "taxon":
        trait_label = nodes.set_index("node_id").loc[row.target, "label"]
        return TRAIT_COLORS[str(trait_label)]
    return SLATE_500


def plot_condition_taxon_trait_network(
    nodes: pd.DataFrame,
    edges: pd.DataFrame,
    conditions: list[str],
    trait_order: list[str],
    output_path: Path,
) -> None:
    """Render the sparse four-layer condition-direction-taxon-trait evidence network."""
    positions = _network_positions(nodes, edges, conditions, trait_order)
    fig, ax = plt.subplots(figsize=(20, 13), facecolor="white")
    ax.set_xlim(-0.38, 3.45)
    ax.set_ylim(-0.02, 1.02)
    ax.axis("off")
    for row in edges.itertuples(index=False):
        x0, y0 = positions[row.source]
        x1, y1 = positions[row.target]
        bend = 0.34 * (x1 - x0)
        path = MplPath(
            [(x0, y0), (x0 + bend, y0), (x1 - bend, y1), (x1, y1)],
            [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4],
        )
        support = max(1.0, float(row.study_support))
        ax.add_patch(
            PathPatch(
                path,
                facecolor="none",
                edgecolor=_network_edge_color(row, nodes),
                linewidth=0.35 + 0.8 * math.sqrt(support),
                alpha=0.23 if row.source_layer == "taxon" else 0.36,
                zorder=1,
            )
        )
    layer_widths = {"condition": 0.38, "direction": 0.2, "taxon": 0.27, "trait": 0.36}
    layer_heights = {"condition": 0.048, "direction": 0.035, "taxon": 0.028, "trait": 0.045}
    for row in nodes.itertuples(index=False):
        x, y = positions[row.node_id]
        width = layer_widths[row.layer]
        height = layer_heights[row.layer]
        if row.layer == "condition":
            color = CONDITION_COLORS.get(row.label, SLATE_500)
        elif row.layer == "direction":
            color = INCREASED_COLOR if INCREASED in row.label else DECREASED_COLOR
        elif row.layer == "trait":
            color = TRAIT_COLORS[row.label]
        else:
            color = "#475569" if row.taxonomic_rank == "genus" else "#94a3b8"
        ax.add_patch(
            FancyBboxPatch(
                (x - width / 2, y - height / 2),
                width,
                height,
                boxstyle="round,pad=0.008,rounding_size=0.009",
                facecolor=color,
                edgecolor="white",
                linewidth=0.8,
                zorder=3,
            )
        )
        label = row.label
        if row.layer == "taxon":
            label = f"{label} [{row.taxonomic_rank[0].upper()}]"
        wrap_width = 15 if row.layer == "trait" else 18
        ax.text(
            x,
            y,
            "\n".join(wrap(label, width=wrap_width)),
            ha="center",
            va="center",
            fontsize=7.2 if row.layer == "taxon" else 7.7,
            color="white",
            fontweight="bold" if row.layer in {"condition", "trait"} else "normal",
            zorder=4,
        )
    for x, title in enumerate(
        ["BugSigDB condition", "Signature direction", "Exact reported taxon", "KG-Microbe trait"]
    ):
        ax.text(x, 1.005, title, ha="center", va="bottom", fontsize=11, fontweight="bold")
    fig.suptitle(
        "Condition → direction → taxon → trait evidence network",
        x=0.04,
        y=0.99,
        ha="left",
        fontsize=20,
        fontweight="bold",
    )
    fig.text(
        0.04,
        0.958,
        "Fecal BugSigDB profiles; top two genus/species labels per condition-direction with "
        "support from ≥2 studies. Line width scales with distinct-study support. Taxon→trait "
        "links are KG-Microbe annotations and branch independently, so widths are not conserved.",
        fontsize=10,
        color=SLATE_700,
    )
    fig.text(
        0.04,
        0.025,
        "[G] genus; [S] species. Exact taxon IDs are used; ranks are not summed. Annotation links "
        "do not establish activity, strain-level co-occurrence, or causality.",
        fontsize=9,
        color=SLATE_500,
    )
    save_figure(fig, output_path)


def write_documentation(
    summary: pd.DataFrame,
    network_nodes: pd.DataFrame,
    network_edges: pd.DataFrame,
) -> None:
    """Write a compact interpretation and audit guide for both figures."""
    pd_summary = summary[summary["condition"].eq(PD_CONDITION)].set_index("trait")
    cellulose = pd_summary.loc[CELLULOSE]
    reporter = pd_summary.loc[REPORTER_COMPOSITE]
    condition_count = summary[summary["paired_study_count"].ge(3)]["condition"].nunique()
    lines = [
        "# PD Ecological-Reporter Specificity and Provenance",
        "",
        f"![Combined figure]({COMBINED_STEM}.png)",
        "",
        "The combined figure separates three evidence levels. Panel A uses paired, study-level "
        "contrasts of annotation-derived profile z-scores. Panel B places PD among all fecal "
        f"condition labels with at least three paired studies ({condition_count} in this dataset). "
        "Panel C intersects exact reported PD genus/species labels with the local KG-Microbe "
        "taxa-trait matrix and colors annotations by BugSigDB direction support.",
        "",
        f"For PD, the study-weighted cellulose increased-minus-decreased contrast is "
        f"`{cellulose.study_weighted_mean_contrast:+.3f}` across "
        f"`{int(cellulose.paired_study_count)}` paired studies. The exploratory reporter "
        f"composite is `{reporter.study_weighted_mean_contrast:+.3f}`. The composite equals "
        "the mean contrasts for moderately halophilic, facultatively aerobic, and arsenate "
        "reduction minus the cellulose contrast; it is not a validated diagnostic score.",
        "",
        f"[Condition/trait audit table]({COMBINED_STEM}_condition_trait_summary.tsv) | "
        f"[Study contrasts]({COMBINED_STEM}_study_contrasts.tsv) | "
        f"[Taxon-panel cells]({COMBINED_STEM}_taxon_panel.tsv) | "
        f"[Vector figure]({COMBINED_STEM}.svg)",
        "",
        "## Condition-Taxon-Trait Network",
        "",
        f"![Network]({NETWORK_STEM}.png)",
        "",
        f"The network contains `{len(network_nodes)}` nodes and `{len(network_edges)}` edges. "
        "BugSigDB provides the condition, signature direction, exact taxon report, and study "
        "support. KG-Microbe provides the taxon-to-trait annotations. Genus and species are "
        "shown together as separate exact IDs but are never summed. Taxon-to-trait links repeat "
        "the taxon's evidence support for each annotation and therefore are not conserved flows.",
        "",
        f"[Network nodes]({NETWORK_STEM}_nodes.tsv) | "
        f"[Network edges]({NETWORK_STEM}_edges.tsv) | "
        f"[Taxon selection]({NETWORK_STEM}_taxon_selection.tsv) | "
        f"[Vector figure]({NETWORK_STEM}.svg)",
        "",
        "## Interpretation Boundary",
        "",
        "These are ecological annotations of taxa found in curated differential-abundance "
        "signatures. They are not measurements of gene expression, metabolite flux, cellulose "
        "degradation, or longitudinal host physiology. The condition comparison is descriptive "
        "and is not an external diagnostic validation.",
        "",
    ]
    (OUTPUT_DIR / "pd_ecological_reporter_figures.md").write_text(
        "\n".join(lines), encoding="utf-8"
    )


def main() -> None:
    """Generate both requested figures and all audit tables."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    metadata, zscores, traits = load_trait_inputs()
    trait_order = FIGURE_TRAIT_ORDER
    missing_figure_traits = set(trait_order) - set(traits["trait_label"])
    if missing_figure_traits:
        raise ValueError(
            f"Figure traits missing from selected trait summary: {missing_figure_traits}"
        )
    summary, study_contrasts = build_condition_trait_summaries(metadata, zscores, traits)
    carriers = read_tsv(OUTPUT_DIR / "pd_six_trait_taxon_carriers.tsv")
    taxon_cells, taxon_ranking = prepare_pd_taxon_panel(carriers, trait_order)
    plot_combined_figure(
        summary,
        taxon_cells,
        taxon_ranking,
        trait_order,
        OUTPUT_DIR / f"{COMBINED_STEM}.png",
    )

    taxa = read_tsv(TAXA_PATH)
    taxa_trait_matrix = read_tsv(
        TAXA_TRAIT_MATRIX_PATH,
        usecols=["taxon_id", "taxon_name", "taxonomic_rank", *traits["trait_id"].tolist()],
    )
    network_nodes, network_edges, taxon_selection = build_condition_taxon_trait_network(
        metadata,
        taxa,
        taxa_trait_matrix,
        traits,
        NETWORK_CONDITIONS,
    )
    plot_condition_taxon_trait_network(
        network_nodes,
        network_edges,
        NETWORK_CONDITIONS,
        trait_order,
        OUTPUT_DIR / f"{NETWORK_STEM}.png",
    )
    write_tsv(summary, OUTPUT_DIR / f"{COMBINED_STEM}_condition_trait_summary.tsv")
    write_tsv(study_contrasts, OUTPUT_DIR / f"{COMBINED_STEM}_study_contrasts.tsv")
    write_tsv(taxon_cells, OUTPUT_DIR / f"{COMBINED_STEM}_taxon_panel.tsv")
    write_tsv(network_nodes, OUTPUT_DIR / f"{NETWORK_STEM}_nodes.tsv")
    write_tsv(network_edges, OUTPUT_DIR / f"{NETWORK_STEM}_edges.tsv")
    write_tsv(taxon_selection, OUTPUT_DIR / f"{NETWORK_STEM}_taxon_selection.tsv")
    write_documentation(summary, network_nodes, network_edges)
    print(f"Wrote combined figure and network to {OUTPUT_DIR}")
    print(
        f"Condition summaries: {len(summary)}; study contrasts: {len(study_contrasts)}; "
        f"network: {len(network_nodes)} nodes / {len(network_edges)} edges"
    )


if __name__ == "__main__":
    main()
