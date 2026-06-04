"""Investigate facultative-aerobe increased-profile signals."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Circle, PathPatch
from matplotlib.path import Path as MplPath

REPO_ROOT = Path(__file__).resolve().parents[1]
VIS_DIR = REPO_ROOT / "trait-visualization"
OUTPUT_DIR = VIS_DIR / "trait-profile-investigation" / "facultative-aerobe-trait-expansion"

GUILD_TRAITS = {
    "METPO:1000608": "facultatively aerobic",
    "METPO:1000602": "aerobic",
    "METPO:1000699": "gram negative",
    "METPO:1000702": "motile",
    "METPO:1000604": "microaerophilic",
    "METPO:2000103|pathways:nitrate_reduction": "nitrate_reduction",
}

TARGET_TAXA = {
    "NCBITaxon:33033": "Parvimonas micra",
    "NCBITaxon:573": "Klebsiella pneumoniae",
    "NCBITaxon:1328": "Streptococcus anginosus",
    "NCBITaxon:562": "Escherichia coli",
    "NCBITaxon:851": "Fusobacterium nucleatum",
    "NCBITaxon:208479": "Enterocloster bolteae",
    "NCBITaxon:84112": "Eggerthella lenta",
    "NCBITaxon:28131": "Prevotella intermedia",
    "NCBITaxon:837": "Porphyromonas gingivalis",
    "NCBITaxon:39778": "Veillonella dispar",
}
EXCLUDED_CONDITION_SUMMARY_LABELS = {
    "Body fat distribution",
    "Exposure",
    "Health trait",
    "Microbiome",
    "Population",
    "Physical activity",
}
MIN_PAIRED_SMD_PROFILE_PAIRS = 3


def read_tsv(path: Path, **kwargs: object) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", **kwargs)


def clean_text(value: object) -> str:
    text = str(value or "").strip()
    if not text or text.lower() == "nan" or text.upper() == "NA":
        return "unknown"
    return text


def safe_ratio(numerator: float, denominator: float) -> float:
    if denominator == 0:
        return math.nan
    return numerator / denominator


def paired_guild_enrichment_stats(group: pd.DataFrame) -> dict[str, float | int]:
    """Calculate paired inc-vs-dec guild enrichment within study/experiment."""
    paired = group.pivot_table(
        index=["study_id", "experiment_id"],
        columns="inc_dec",
        values="guild_z_mean",
        aggfunc="mean",
    )
    if not {"inc", "dec"}.issubset(paired.columns):
        return {
            "paired_profile_pair_count": 0,
            "paired_mean_guild_z_difference": math.nan,
            "paired_guild_z_difference_sd": math.nan,
            "paired_standardized_mean_difference": math.nan,
        }

    paired = paired.dropna(subset=["inc", "dec"])
    differences = paired["inc"] - paired["dec"]
    pair_count = len(differences)
    mean_difference = differences.mean() if pair_count else math.nan
    difference_sd = differences.std(ddof=1) if pair_count > 1 else math.nan
    smd = math.nan
    if (
        pair_count >= MIN_PAIRED_SMD_PROFILE_PAIRS
        and not pd.isna(difference_sd)
        and difference_sd > 0
    ):
        smd = mean_difference / difference_sd

    return {
        "paired_profile_pair_count": pair_count,
        "paired_mean_guild_z_difference": mean_difference,
        "paired_guild_z_difference_sd": difference_sd,
        "paired_standardized_mean_difference": smd,
    }


def profile_id_from_taxa_row(row: pd.Series) -> str:
    direction = "inc" if row["direction"] == "increased" else "dec"
    return f"{row['study_id']}|{row['experiment_id']}|{direction}"


def load_study_heatmap() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    matrix = read_tsv(VIS_DIR / "study_trait_profile_clustered_heatmap_matrix.tsv", index_col=0)
    metadata = read_tsv(VIS_DIR / "study_trait_profile_study_metadata.tsv", index_col=0)
    trait_metadata = read_tsv(VIS_DIR / "study_trait_profile_trait_metadata.tsv", index_col=0)

    matrix.index = metadata.index
    column_labels = {
        row["display_label"]: row["object_name"] for _, row in trait_metadata.iterrows()
    }
    matrix = matrix.rename(columns=column_labels)
    return matrix, metadata, trait_metadata


def add_guild_scores(matrix: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    available_traits = [trait for trait in GUILD_TRAITS.values() if trait in matrix.columns]
    scored = metadata.copy()
    scored["guild_z_mean"] = matrix[available_traits].mean(axis=1)
    scored["oxygen_tolerance_z_mean"] = matrix[
        ["facultatively aerobic", "aerobic", "microaerophilic", "nitrate_reduction"]
    ].mean(axis=1)
    scored["gram_negative_motility_z_mean"] = matrix[["gram negative", "motile"]].mean(axis=1)
    for trait_name in available_traits:
        scored[f"z_{trait_name}"] = matrix[trait_name]
    return scored


def summarize_direction(metadata: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for trait_name in GUILD_TRAITS.values():
        column = f"z_{trait_name}"
        inc_mean = metadata.loc[metadata["inc_dec"] == "inc", column].mean()
        dec_mean = metadata.loc[metadata["inc_dec"] == "dec", column].mean()
        rows.append(
            {
                "trait": trait_name,
                "increased_profile_mean_z": inc_mean,
                "decreased_profile_mean_z": dec_mean,
                "increased_minus_decreased_z": inc_mean - dec_mean,
            }
        )
    for score_name in [
        "guild_z_mean",
        "oxygen_tolerance_z_mean",
        "gram_negative_motility_z_mean",
    ]:
        inc_mean = metadata.loc[metadata["inc_dec"] == "inc", score_name].mean()
        dec_mean = metadata.loc[metadata["inc_dec"] == "dec", score_name].mean()
        rows.append(
            {
                "trait": score_name,
                "increased_profile_mean_z": inc_mean,
                "decreased_profile_mean_z": dec_mean,
                "increased_minus_decreased_z": inc_mean - dec_mean,
            }
        )
    return pd.DataFrame(rows)


def summarize_group(metadata: pd.DataFrame, group_column: str, min_profiles: int) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    for group_value, group in metadata.groupby(group_column, dropna=False):
        if len(group) < min_profiles:
            continue
        inc = group[group["inc_dec"] == "inc"]
        dec = group[group["inc_dec"] == "dec"]
        record = {
            group_column: clean_text(group_value),
            "profile_count": len(group),
            "study_count": group["study_id"].nunique(),
            "increased_profile_count": len(inc),
            "decreased_profile_count": len(dec),
            "mean_guild_z": group["guild_z_mean"].mean(),
            "increased_mean_guild_z": inc["guild_z_mean"].mean(),
            "decreased_mean_guild_z": dec["guild_z_mean"].mean(),
            "increased_minus_decreased_guild_z": inc["guild_z_mean"].mean()
            - dec["guild_z_mean"].mean(),
            "increased_mean_oxygen_tolerance_z": inc["oxygen_tolerance_z_mean"].mean(),
            "increased_mean_gram_negative_motility_z": inc["gram_negative_motility_z_mean"].mean(),
        }
        if group_column == "condition":
            record.update(paired_guild_enrichment_stats(group))
        records.append(record)
    return pd.DataFrame(records).sort_values(
        ["increased_mean_guild_z", "profile_count"],
        ascending=[False, False],
    )


def remove_excluded_condition_labels(condition_summary: pd.DataFrame) -> pd.DataFrame:
    """Remove broad study-design labels from condition summaries."""
    return condition_summary[
        ~condition_summary["condition"].map(is_excluded_condition_label)
    ].copy()


def is_excluded_condition_label(value: object) -> bool:
    """Return whether a condition label is too broad for disease-context summaries."""
    label = clean_text(value)
    broad_terms = ("exposure", "measurement", "population")
    lower_label = label.lower()
    return label in EXCLUDED_CONDITION_SUMMARY_LABELS or any(
        term in lower_label for term in broad_terms
    )


def select_top_conditions(condition_summary: pd.DataFrame, limit: int) -> pd.DataFrame:
    """Select top condition rows after removing broad excluded labels."""
    condition_summary = remove_excluded_condition_labels(condition_summary)
    top_conditions = (
        condition_summary[
            (condition_summary["increased_profile_count"] >= 3)
            & (condition_summary["study_count"] >= 2)
            & (condition_summary["increased_minus_decreased_guild_z"] > 0)
            & (condition_summary["increased_mean_guild_z"] > 0)
        ]
        .sort_values(
            ["increased_minus_decreased_guild_z", "increased_mean_guild_z"],
            ascending=[False, False],
        )
        .head(limit)
    )
    if top_conditions.empty:
        top_conditions = condition_summary[condition_summary["increased_profile_count"] >= 3].head(
            limit
        )
    return top_conditions


def paired_smd_label(row: Any) -> str:
    pair_count = int(row.paired_profile_pair_count)
    smd = row.paired_standardized_mean_difference
    if pd.isna(smd):
        return f"paired d=NA, n={pair_count}"
    return f"paired d={smd:.2f}, n={pair_count}"


def load_taxa_direction_summary() -> pd.DataFrame:
    taxa_metadata = read_tsv(VIS_DIR / "taxa_guild_discovery_taxa_metadata.tsv", index_col=0)
    summary = taxa_metadata.loc[taxa_metadata.index.intersection(TARGET_TAXA.keys())].copy()
    summary["target_taxon_name"] = summary.index.map(TARGET_TAXA)
    summary["increased_fraction"] = [
        safe_ratio(float(row["increased"]), float(row["bugsigdb_profile_count"]))
        for _, row in summary.iterrows()
    ]
    summary["increased_profile_tendency"] = [
        safe_ratio(
            float(row["increased"]) - float(row["decreased"]),
            float(row["bugsigdb_profile_count"]),
        )
        for _, row in summary.iterrows()
    ]
    summary["increased_to_decreased_ratio"] = [
        safe_ratio(float(row["increased"]), float(row["decreased"]))
        for _, row in summary.iterrows()
    ]
    return summary[
        [
            "target_taxon_name",
            "taxon_name",
            "taxonomic_rank",
            "bugsigdb_study_count",
            "bugsigdb_profile_count",
            "decreased",
            "increased",
            "dominant_direction",
            "increased_fraction",
            "increased_profile_tendency",
            "increased_to_decreased_ratio",
            "trait_count",
        ]
    ].sort_values("increased_profile_tendency", ascending=False)


def load_target_taxa_trait_counts() -> pd.DataFrame:
    usecols = ["taxon_id", "taxon_name", "taxonomic_rank", *GUILD_TRAITS.keys()]
    matrix = read_tsv(REPO_ROOT / "data/derived/taxa_trait_matrix.tsv", usecols=usecols)
    matrix = matrix[matrix["taxon_id"].isin(TARGET_TAXA)].copy()
    taxa_metadata = read_tsv(VIS_DIR / "taxa_guild_discovery_taxa_metadata.tsv", index_col=0)
    matrix["target_taxon_name"] = matrix["taxon_id"].map(TARGET_TAXA)
    matrix["heatmap_taxon_name"] = matrix["taxon_id"].map(taxa_metadata["taxon_name"])
    matrix["heatmap_taxonomic_rank"] = matrix["taxon_id"].map(taxa_metadata["taxonomic_rank"])

    for trait_id, trait_name in GUILD_TRAITS.items():
        matrix[f"{trait_name}_count"] = pd.to_numeric(matrix[trait_id], errors="coerce").fillna(0)
        matrix[f"{trait_name}_present"] = matrix[f"{trait_name}_count"] > 0

    present_columns = [f"{trait_name}_present" for trait_name in GUILD_TRAITS.values()]
    matrix["guild_trait_presence_count"] = matrix[present_columns].sum(axis=1)
    output_columns = [
        "taxon_id",
        "target_taxon_name",
        "heatmap_taxonomic_rank",
        "guild_trait_presence_count",
    ]
    for trait_name in GUILD_TRAITS.values():
        output_columns.extend([f"{trait_name}_present", f"{trait_name}_count"])
    output = matrix[output_columns].rename(columns={"heatmap_taxonomic_rank": "taxonomic_rank"})
    return output.sort_values(
        ["guild_trait_presence_count", "target_taxon_name"],
        ascending=[False, True],
    )


def load_target_taxa_occurrences(scored_profiles: pd.DataFrame) -> pd.DataFrame:
    taxa = read_tsv(
        REPO_ROOT / "data/derived/taxa.tsv",
        usecols=[
            "study_id",
            "experiment_id",
            "signature_id",
            "direction",
            "taxon_id",
            "taxon_name",
            "taxonomic_rank",
        ],
    )
    taxa = taxa[taxa["taxon_id"].isin(TARGET_TAXA)].copy()
    taxa["profile_id"] = taxa.apply(profile_id_from_taxa_row, axis=1)
    taxa = taxa.drop_duplicates(subset=["taxon_id", "profile_id"])
    taxa["target_taxon_name"] = taxa["taxon_id"].map(TARGET_TAXA)
    profile_columns = [
        "study_id",
        "experiment_id",
        "inc_dec",
        "direction_label",
        "title",
        "publication_year",
        "condition",
        "body_site",
        "body_site_group",
        "condition_group",
        "location",
        "host_species",
        "guild_z_mean",
        "oxygen_tolerance_z_mean",
        "gram_negative_motility_z_mean",
    ]
    profiles = scored_profiles[profile_columns].copy()
    return taxa.merge(
        profiles,
        left_on="profile_id",
        right_index=True,
        how="left",
        suffixes=("_taxa", ""),
    )


def summarize_taxa_conditions(occurrences: pd.DataFrame) -> pd.DataFrame:
    grouped = (
        occurrences.groupby(["target_taxon_name", "condition", "body_site_group", "direction"])[
            ["profile_id", "study_id"]
        ]
        .nunique()
        .rename(columns={"profile_id": "profile_count", "study_id": "study_count"})
        .reset_index()
    )
    if grouped.empty:
        return grouped

    profile_pivot = grouped.pivot_table(
        index=["target_taxon_name", "condition", "body_site_group"],
        columns="direction",
        values="profile_count",
        aggfunc="sum",
        fill_value=0,
    )
    study_pivot = grouped.pivot_table(
        index=["target_taxon_name", "condition", "body_site_group"],
        columns="direction",
        values="study_count",
        aggfunc="sum",
        fill_value=0,
    )
    for direction in ["decreased", "increased"]:
        if direction not in profile_pivot.columns:
            profile_pivot[direction] = 0
        if direction not in study_pivot.columns:
            study_pivot[direction] = 0

    output = profile_pivot.rename(
        columns={
            "decreased": "decreased_profile_count",
            "increased": "increased_profile_count",
        }
    ).join(
        study_pivot.rename(
            columns={
                "decreased": "decreased_study_count",
                "increased": "increased_study_count",
            }
        )
    )
    output = output.reset_index()
    output["total_profile_count"] = (
        output["decreased_profile_count"] + output["increased_profile_count"]
    )
    output["increased_profile_tendency"] = (
        output["increased_profile_count"] - output["decreased_profile_count"]
    ) / output["total_profile_count"].replace(0, pd.NA)
    return output.sort_values(
        ["increased_profile_tendency", "increased_study_count", "total_profile_count"],
        ascending=[False, False, False],
    )


def count_target_taxa_conditions(occurrences: pd.DataFrame, direction: str) -> pd.Series:
    """Count condition labels for target taxon-profile occurrences."""
    selected = occurrences[occurrences["direction"] == direction].copy()
    selected["condition"] = selected["condition"].map(clean_text)
    selected = selected[~selected["condition"].map(is_excluded_condition_label)].copy()
    return selected["condition"].value_counts()


def save_condition_breakdown_bar_chart(
    counts: pd.Series,
    *,
    output_path: Path,
    title_prefix: str,
    color: str,
) -> None:
    """Save a horizontal bar chart for condition counts."""
    if counts.empty:
        return

    bar_counts = counts.head(20)
    fig, ax = plt.subplots(figsize=(10, max(5, len(bar_counts) * 0.35)))
    ax.barh(bar_counts.index, bar_counts.to_numpy(), color=color)
    ax.invert_yaxis()
    ax.set_xlabel("Target taxon-profile occurrences")
    ax.set_title(f"{title_prefix}: top conditions")
    for index, value in enumerate(bar_counts.to_numpy()):
        ax.text(value, index, f" {int(value)}", va="center", fontsize=8)
    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def display_body_site(value: object) -> str:
    """Return concise display labels for body-site groups."""
    label = clean_text(value)
    if label.lower() == "unknown":
        return "unknown body site"
    return label


def is_other_body_site(value: object) -> bool:
    return clean_text(value).lower() == "other"


def anatomical_site_base(value: object) -> tuple[float, float]:
    """Map body-site group labels to approximate body-diagram coordinates."""
    label = clean_text(value).lower()
    if "nasopharynx" in label or "airway" in label:
        return 0.30, 0.86
    if "subgingival" in label or "dental" in label:
        return 0.275, 0.825
    if "saliva" in label:
        return 0.325, 0.825
    if "mouth" in label or "oral" in label:
        return 0.30, 0.825
    if "colon" in label:
        return 0.285, 0.47
    if "caecum" in label or "cecum" in label:
        return 0.33, 0.46
    if "feces" in label or "faeces" in label or "stool" in label:
        return 0.31, 0.37
    if "vagina" in label or "urogenital" in label:
        return 0.31, 0.25
    if "unknown" in label:
        return 0.16, 0.16
    return 0.16, 0.47


def add_anatomical_outline(ax: plt.Axes) -> None:
    """Draw an original front-facing body-template outline for site callouts."""
    face = "#f8fafc"
    edge = "#94a3b8"

    def add_path(
        start: tuple[float, float],
        segments: list[tuple[tuple[float, float], ...]],
    ) -> None:
        vertices = [start]
        codes = [MplPath.MOVETO]
        for segment in segments:
            if len(segment) == 1:
                vertices.append(segment[0])
                codes.append(MplPath.LINETO)
            else:
                vertices.extend(segment)
                codes.extend([MplPath.CURVE4] * 3)
        vertices.append(start)
        codes.append(MplPath.CLOSEPOLY)
        ax.add_patch(
            PathPatch(
                MplPath(vertices, codes),
                facecolor=face,
                edgecolor=edge,
                lw=1.4,
                joinstyle="round",
                zorder=1,
            )
        )

    def mirror_point(point: tuple[float, float]) -> tuple[float, float]:
        return 0.60 - point[0], point[1]

    def mirror_segments(
        segments: list[tuple[tuple[float, float], ...]],
    ) -> list[tuple[tuple[float, float], ...]]:
        return [tuple(mirror_point(point) for point in segment) for segment in segments]

    ax.add_patch(
        Circle(
            (0.30, 0.87),
            0.045,
            facecolor=face,
            edgecolor=edge,
            lw=1.4,
            zorder=1,
        )
    )
    torso_segments = [
        ((0.224, 0.800), (0.207, 0.775), (0.196, 0.735)),
        ((0.178, 0.660), (0.178, 0.565), (0.195, 0.480)),
        ((0.205, 0.420), (0.220, 0.355), (0.240, 0.320)),
        ((0.265, 0.295), (0.335, 0.295), (0.360, 0.320)),
        ((0.380, 0.355), (0.395, 0.420), (0.405, 0.480)),
        ((0.422, 0.565), (0.422, 0.660), (0.404, 0.735)),
        ((0.393, 0.775), (0.376, 0.800), (0.355, 0.805)),
        ((0.335, 0.812), (0.265, 0.812), (0.245, 0.805)),
    ]
    left_arm_segments = [
        ((0.175, 0.665), (0.160, 0.585), (0.142, 0.505)),
        ((0.127, 0.445), (0.112, 0.395), (0.105, 0.365)),
        ((0.113, 0.350), (0.135, 0.346), (0.150, 0.358)),
        ((0.162, 0.395), (0.178, 0.455), (0.198, 0.535)),
        ((0.207, 0.600), (0.213, 0.665), (0.224, 0.715)),
    ]
    left_leg_segments = [
        ((0.258, 0.245), (0.250, 0.170), (0.238, 0.090)),
        ((0.233, 0.060), (0.250, 0.045), (0.278, 0.050)),
        ((0.292, 0.100), (0.295, 0.225), (0.295, 0.300)),
    ]
    add_path((0.245, 0.805), torso_segments)
    add_path((0.205, 0.725), left_arm_segments)
    add_path(mirror_point((0.205, 0.725)), mirror_segments(left_arm_segments))
    add_path((0.265, 0.300), left_leg_segments)
    add_path(mirror_point((0.265, 0.300)), mirror_segments(left_leg_segments))
    ax.plot([0.280, 0.280], [0.805, 0.770], color=edge, lw=1.2, zorder=2)
    ax.plot([0.320, 0.320], [0.805, 0.770], color=edge, lw=1.2, zorder=2)


def save_anatomical_sites_diagram(
    body_site_summary: pd.DataFrame,
    *,
    output_path: Path,
    title: str,
    signal_column: str,
    signal_label: str,
    legend_note: str | None = None,
    max_sites: int = 12,
) -> None:
    """Save an anatomical-site callout diagram for body-site guild summaries."""
    if body_site_summary.empty:
        return

    plot_sites = (
        body_site_summary[~body_site_summary["body_site_group"].map(is_other_body_site)]
        .head(max_sites)
        .copy()
    )
    if plot_sites.empty:
        return
    plot_sites["profile_count"] = pd.to_numeric(
        plot_sites["profile_count"],
        errors="coerce",
    ).fillna(0)
    plot_sites[signal_column] = pd.to_numeric(
        plot_sites[signal_column],
        errors="coerce",
    ).fillna(0)
    points = [anatomical_site_base(site) for site in plot_sites["body_site_group"]]
    plot_sites["anchor_x"] = [point[0] for point in points]
    plot_sites["anchor_y"] = [point[1] for point in points]
    plot_sites = plot_sites.sort_values(
        ["anchor_y", "profile_count"],
        ascending=[False, False],
    )

    max_count = max(float(plot_sites["profile_count"].max()), 1.0)
    max_signal = max(float(plot_sites[signal_column].abs().max()), 0.1)
    norm = TwoSlopeNorm(vmin=-max_signal, vcenter=0, vmax=max_signal)
    cmap = plt.get_cmap("RdBu_r")

    fig_height = max(7.2, 1.0 + len(plot_sites) * 0.82)
    fig, ax = plt.subplots(figsize=(11, fig_height))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")
    add_anatomical_outline(ax)

    label_ys = [0.84]
    if len(plot_sites) > 1:
        step = 0.70 / (len(plot_sites) - 1)
        label_ys = [0.84 - index * step for index in range(len(plot_sites))]

    for label_y, row in zip(label_ys, plot_sites.itertuples(index=False), strict=False):
        count = int(row.profile_count)
        signal = float(getattr(row, signal_column))
        color = cmap(norm(signal))
        size = 95 + 900 * math.sqrt(count / max_count)
        ax.scatter(
            [row.anchor_x],
            [row.anchor_y],
            s=size,
            c=[color],
            edgecolors="#111827",
            linewidths=0.8,
            zorder=4,
        )
        ax.plot(
            [row.anchor_x, 0.52, 0.57],
            [row.anchor_y, label_y, label_y],
            color="#6b7280",
            lw=0.9,
            zorder=2,
        )
        bar_width = 0.30 * count / max_count
        ax.plot(
            [0.58, 0.58 + bar_width],
            [label_y, label_y],
            color=color,
            lw=7,
            solid_capstyle="round",
            zorder=3,
        )
        ax.text(
            0.58,
            label_y + 0.018,
            display_body_site(row.body_site_group),
            ha="left",
            va="bottom",
            fontsize=10,
            fontweight="bold",
        )
        ax.text(
            0.58,
            label_y - 0.018,
            f"{count:,} profiles; {signal_label} {signal:+.2f}",
            ha="left",
            va="top",
            fontsize=8.5,
            color="#374151",
        )

    ax.text(0.02, 0.965, title, ha="left", va="top", fontsize=17, fontweight="bold")
    footer = (
        "Bubble and bar length show profile count; color shows the direction-specific guild signal."
    )
    if legend_note:
        footer = f"{footer} {legend_note}"
    ax.text(
        0.02,
        0.035,
        footer,
        ha="left",
        va="bottom",
        fontsize=8.5,
        color="#374151",
    )
    scalar_mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    scalar_mappable.set_array([])
    colorbar = fig.colorbar(scalar_mappable, ax=ax, fraction=0.03, pad=0.02)
    colorbar.set_label(signal_label)

    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def write_figure_methods() -> None:
    """Write methods text for generated figure files."""
    lines = [
        "# Figure Methods",
        "",
        "These method notes describe how each PNG in this investigation folder was "
        "generated by `scripts/investigate_facultative_aerobe_expansion.py`.",
        "",
        "## Shared Inputs",
        "",
        "Study-profile z-scores were read from "
        "`trait-visualization/study_trait_profile_clustered_heatmap_matrix.tsv` and "
        "joined to study metadata from "
        "`trait-visualization/study_trait_profile_study_metadata.tsv`. The "
        "facultative-aerobe guild score was the row-wise mean of the selected guild "
        "traits: facultatively aerobic, aerobic, gram negative, motile, "
        "microaerophilic, and nitrate_reduction.",
        "Direction labels (`increased` and `decreased`) were taken from BugSigDB's "
        "curated `Abundance in Group 1` signature field. This script did not "
        "recalculate differential abundance; it grouped and filtered "
        "BugSigDB-reported study/experiment/direction profiles by condition and "
        "body-site metadata.",
        "",
        "The focal taxa were defined a priori in the script's `TARGET_TAXA` map "
        "using NCBI Taxonomy identifiers, which avoids ambiguity from synonym or "
        "strain-name differences. The selected taxa were Parvimonas micra "
        "(NCBITaxon:33033), Klebsiella pneumoniae (NCBITaxon:573), Streptococcus "
        "anginosus (NCBITaxon:1328), Escherichia coli (NCBITaxon:562), "
        "Fusobacterium nucleatum (NCBITaxon:851), Enterocloster bolteae "
        "(NCBITaxon:208479), Eggerthella lenta (NCBITaxon:84112), Prevotella "
        "intermedia (NCBITaxon:28131), Porphyromonas gingivalis (NCBITaxon:837), "
        "and Veillonella dispar (NCBITaxon:39778). They were chosen as focal "
        "examples for the facultative-aerobe increased-profile signal because "
        "they include taxa with oxygen-tolerance, gram-negative, motility, "
        "microaerophilic, or nitrate-reduction annotations in the derived "
        "taxa-trait matrix, and because they provide gut, oral, and airway "
        "contexts where facultative-aerobe or aerotolerant increased-profile "
        "patterns may be "
        "relevant. The script writes `candidate_taxa_trait_presence.tsv` from "
        "`data/derived/taxa_trait_matrix.tsv` as the audit trail for those trait "
        "annotations; some focal taxa carry broad multi-trait support, while "
        "others are retained as single-trait focal taxa for condition/body-site "
        "follow-up.",
        "",
        "## `candidate_taxa_increased_profile_tendency.png`",
        "",
        "This figure uses the fixed focal taxa listed above, then intersects their "
        "NCBI Taxonomy IDs with "
        "`trait-visualization/taxa_guild_discovery_taxa_metadata.tsv`. That "
        "metadata table supplies the total BugSigDB profile count and the number "
        "of profiles where each taxon was reported as increased or decreased. For "
        "each target taxon, increased-profile tendency was calculated as "
        "`(increased profile count - decreased profile count) / total BugSigDB "
        "profile count`. The horizontal bars show this signed value for each target "
        "taxon. No "
        "condition or body-site filter is applied to this figure; it summarizes "
        "the overall BugSigDB direction tendency of the preselected focal taxa.",
        "",
        "## `condition_increased_guild_signal.png`",
        "",
        "Study profiles from all body sites were grouped by BugSigDB condition "
        "label. Groups with fewer than five profiles were excluded. Broad "
        "non-condition or metadata labels (`Body fat distribution`, `Exposure`, "
        "`Health trait`, `Microbiome`, `Population`, `Physical activity`, and "
        "labels containing `exposure`, `measurement`, or `population`) were removed "
        "from the displayed condition "
        "summary. Those labels usually describe exposure, cohort, or measurement "
        "metadata rather than a disease or condition context. The figure keeps the "
        "top 12 eligible conditions with at least three increased profiles, at least "
        "two studies, positive increased-minus-decreased guild signal, and positive "
        "increased-profile mean guild z-score, using fallback filtering if needed. "
        "Conditions are selected and sorted by increased-minus-decreased guild signal. "
        "Bar length is this enrichment-oriented contrast, and the right-side labels "
        "show the paired standardized mean difference (`paired d`) calculated from "
        "matched increased/decreased profiles within the same study and experiment. "
        "`paired d=NA` means fewer than three matched profile pairs, or zero variance "
        "among the paired differences.",
        "",
        "## `guild_trait_direction_shift.png`",
        "",
        "For each selected guild trait, the script calculated the mean z-score in "
        "all increased study profiles and the mean z-score in all decreased study "
        "profiles. Bars show `increased mean z - decreased mean z`; positive values "
        "indicate traits that are higher in increased profiles.",
        "",
        "## `anatomical_sites_diagram.png`",
        "",
        "Study profiles were grouped by `body_site_group`, requiring at least 10 "
        "profiles per group. The display excludes the generic `Other` body-site "
        "group. Bubble and bar length encode the number of study profiles in the "
        "body-site group. Color encodes the mean facultative-aerobe guild z-score "
        "among increased profiles at that body site. Body-site labels were placed "
        "on a simplified, neutral front-view body outline using fixed anatomical "
        "coordinates for feces/stool, colon, caecum, oral, airway, urogenital, and "
        "unknown labels. Unknown body site means the BugSigDB profile did not have "
        "a specific body-site label; these profiles are placed off-body and should "
        "not be interpreted as a physical anatomical location.",
        "",
        "Oral-site increased-profile signals are biologically plausible because the "
        "oral cavity is oxygen-exposed and includes aerotolerant, facultatively "
        "aerobic, or microaerophilic taxa such as Streptococcus, Fusobacterium, "
        "Prevotella, Porphyromonas, and Veillonella. Those oral-site signals should "
        "be interpreted as body-site-specific microbial ecology, not as a gut "
        "facultative-aerobe increased-profile signal.",
        "",
        "## `target_taxa_condition_breakdown_bar.png`",
        "",
        "Target-taxon occurrences were loaded from `data/derived/taxa.tsv`, limited "
        "to the fixed focal NCBI Taxonomy IDs, deduplicated by taxon and study "
        "profile, and joined to the scored study profiles. This produces one row "
        "per focal taxon/profile occurrence with condition, body site, direction, "
        "and guild-score metadata. The bar chart counts increased focal-taxon/profile "
        "occurrences by condition after removing the broad non-condition and "
        "metadata labels described above. A condition with multiple focal taxa can "
        "therefore contribute multiple occurrences; the chart is an occurrence "
        "breakdown, not a unique-study breakdown. It displays the top 20 explicit "
        "condition labels as horizontal bars, with bar labels showing occurrence "
        "counts.",
        "",
    ]
    (OUTPUT_DIR / "figure_methods.md").write_text("\n".join(lines), encoding="utf-8")


def make_plots(
    taxa_summary: pd.DataFrame,
    condition_summary: pd.DataFrame,
    direction_summary: pd.DataFrame,
    body_site_summary: pd.DataFrame,
    occurrences: pd.DataFrame,
) -> None:
    plt.rcParams.update({"figure.dpi": 160, "savefig.dpi": 200})

    fig, ax = plt.subplots(figsize=(9, 5.4))
    plot_taxa = taxa_summary.sort_values("increased_profile_tendency")
    ax.barh(
        plot_taxa["target_taxon_name"],
        plot_taxa["increased_profile_tendency"],
        color="#b91c1c",
    )
    ax.axvline(0, color="#333333", linewidth=0.8)
    ax.set_xlabel("(increased profiles - decreased profiles) / total profiles")
    ax.set_title("Facultative-Aerobe Target Taxa: Increased-Profile Tendency")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "candidate_taxa_increased_profile_tendency.png")
    plt.close(fig)

    top_conditions = select_top_conditions(condition_summary, 12).sort_values(
        "increased_minus_decreased_guild_z"
    )
    fig, ax = plt.subplots(figsize=(10, 5.5))
    condition_values = top_conditions["increased_minus_decreased_guild_z"]
    ax.barh(
        top_conditions["condition"],
        condition_values,
        color="#dc2626",
    )
    if not top_conditions.empty:
        x_max = condition_values.max()
        x_pad = max(x_max * 0.025, 0.02)
        ax.set_xlim(0, x_max * 1.52)
        for y_pos, row in enumerate(top_conditions.itertuples(index=False)):
            ax.text(
                row.increased_minus_decreased_guild_z + x_pad,
                y_pos,
                paired_smd_label(row),
                va="center",
                fontsize=8,
                color="#334155",
            )
    ax.set_xlabel("Guild enrichment: increased mean z - decreased mean z")
    ax.set_title("Conditions With Increased-Profile Guild Enrichment")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "condition_increased_guild_signal.png")
    plt.close(fig)

    plot_direction = direction_summary[
        direction_summary["trait"].isin(GUILD_TRAITS.values())
    ].sort_values("increased_minus_decreased_z")
    fig, ax = plt.subplots(figsize=(8, 4.4))
    ax.barh(
        plot_direction["trait"],
        plot_direction["increased_minus_decreased_z"],
        color="#be123c",
    )
    ax.axvline(0, color="#333333", linewidth=0.8)
    ax.set_xlabel("increased profile mean z - decreased profile mean z")
    ax.set_title("Guild Traits Are Higher in Increased Profiles")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "guild_trait_direction_shift.png")
    plt.close(fig)

    condition_counts = count_target_taxa_conditions(occurrences, "increased")
    save_condition_breakdown_bar_chart(
        condition_counts,
        output_path=OUTPUT_DIR / "target_taxa_condition_breakdown_bar.png",
        title_prefix="Increased facultative-aerobe target taxa",
        color="#dc2626",
    )
    save_anatomical_sites_diagram(
        body_site_summary,
        output_path=OUTPUT_DIR / "anatomical_sites_diagram.png",
        title="Facultative-Aerobe Increased-Profile Signal: Anatomical Sites",
        signal_column="increased_mean_guild_z",
        signal_label="increased-profile guild z",
        legend_note=(
            "Unknown body site means no specific BugSigDB body-site label was available; "
            "it is plotted off-body."
        ),
    )


def write_summary_report(
    direction_summary: pd.DataFrame,
    body_site_summary: pd.DataFrame,
    condition_summary: pd.DataFrame,
    taxa_summary: pd.DataFrame,
    condition_taxa_summary: pd.DataFrame,
    top_increased_profiles: pd.DataFrame,
) -> None:
    direction_rows = direction_summary[
        direction_summary["trait"].isin(GUILD_TRAITS.values())
    ].sort_values("increased_minus_decreased_z", ascending=False)
    top_conditions = select_top_conditions(condition_summary, 8)
    top_taxa = taxa_summary.head(10)
    top_taxon_conditions = condition_taxa_summary[
        (condition_taxa_summary["increased_profile_count"] >= 3)
        & (condition_taxa_summary["increased_study_count"] >= 2)
    ].head(12)

    lines = [
        "# Facultative-Aerobe Increased-Profile Signal Investigation",
        "",
        "This investigation uses the existing study trait profile clustered heatmap matrix, "
        "study metadata, taxa guild heatmap metadata, and derived taxa/trait matrices.",
        "",
        "## Scope",
        "",
        "- Guild traits: "
        + ", ".join(f"{name} ({trait_id})" for trait_id, name in GUILD_TRAITS.items()),
        "- Target taxa: " + ", ".join(TARGET_TAXA.values()),
        "- Context: all body sites, with body-site and condition summaries to separate "
        "gut, oral, airway, and other patterns.",
        "",
        "## Main Signal",
        "",
        "The selected guild traits are higher in increased profiles than in decreased "
        "profiles, supporting a facultative-aerobe increased-profile signal.",
        "",
        "| trait | increased - decreased z |",
        "| --- | ---: |",
    ]
    for row in direction_rows.itertuples(index=False):
        lines.append(f"| {row.trait} | {row.increased_minus_decreased_z:.3f} |")

    lines.extend(
        [
            "",
            "## Target Taxa Most Often Reported In Increased Profiles",
            "",
            "| taxon | increased / decreased profiles | increased-profile tendency |",
            "| --- | ---: | ---: |",
        ]
    )
    for row in top_taxa.itertuples(index=False):
        lines.append(
            f"| {row.target_taxon_name} | {int(row.increased)} / {int(row.decreased)} | "
            f"{row.increased_profile_tendency:.3f} |"
        )

    lines.extend(
        [
            "",
            "## Recurrent Target-Taxon Increased Profiles By Condition",
            "",
            "| taxon | condition | body site | "
            "increased / decreased profiles | increased studies |",
            "| --- | --- | --- | ---: | ---: |",
        ]
    )
    for row in top_taxon_conditions.itertuples(index=False):
        lines.append(
            f"| {row.target_taxon_name} | {clean_text(row.condition)} | "
            f"{clean_text(row.body_site_group)} | "
            f"{int(row.increased_profile_count)} / {int(row.decreased_profile_count)} | "
            f"{int(row.increased_study_count)} |"
        )

    lines.extend(
        [
            "",
            "## Contexts With High Increased-Profile Guild Signal",
            "",
            "| condition | profiles | studies | increased - decreased z | paired d | "
            "matched pairs | increased mean guild z |",
            "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in top_conditions.itertuples(index=False):
        diff = row.increased_minus_decreased_guild_z
        diff_text = "" if pd.isna(diff) else f"{diff:.3f}"
        smd = row.paired_standardized_mean_difference
        smd_text = "" if pd.isna(smd) else f"{smd:.3f}"
        lines.append(
            f"| {row.condition} | {int(row.profile_count)} | {int(row.study_count)} | "
            f"{diff_text} | {smd_text} | {int(row.paired_profile_pair_count)} | "
            f"{row.increased_mean_guild_z:.3f} |"
        )

    lines.extend(
        [
            "",
            "Broad labels such as `Exposure`, `Population`, and labels containing "
            "`measurement` come from BugSigDB condition metadata but usually describe "
            "study design, cohort, exposure, or measurement context rather than a "
            "disease/condition context. They are excluded from displayed condition "
            "bar charts and condition summaries.",
            "",
            "Oral-site increased-profile signals are expected in this analysis because "
            "the oral cavity is oxygen-exposed and several focal taxa are common oral "
            "or airway-associated aerotolerant, facultatively aerobic, or "
            "microaerophilic lineages. These signals should be interpreted as "
            "body-site-specific oral microbial ecology, not as a gut "
            "increased-profile signal.",
        ]
    )

    lines.extend(
        [
            "",
            "## Body-Site Check",
            "",
            "| body_site_group | profiles | increased mean guild z | increased - decreased z |",
            "| --- | ---: | ---: | ---: |",
        ]
    )
    for row in body_site_summary.head(10).itertuples(index=False):
        diff = row.increased_minus_decreased_guild_z
        diff_text = "" if pd.isna(diff) else f"{diff:.3f}"
        lines.append(
            f"| {row.body_site_group} | {int(row.profile_count)} | "
            f"{row.increased_mean_guild_z:.3f} | {diff_text} |"
        )

    lines.extend(
        [
            "",
            "## Top Increased Profiles To Inspect",
            "",
            "| profile | condition | body site | guild z | title |",
            "| --- | --- | --- | ---: | --- |",
        ]
    )
    for row in top_increased_profiles.head(10).itertuples(index=False):
        title = clean_text(row.title).replace("|", "/")
        if len(title) > 90:
            title = f"{title[:87]}..."
        profile_id = str(row.profile_id).replace("|", "\\|")
        lines.append(
            f"| {profile_id} | {clean_text(row.condition)} | {clean_text(row.body_site)} | "
            f"{row.guild_z_mean:.3f} | {title} |"
        )

    lines.extend(
        [
            "",
            "## Suggested Follow-Up",
            "",
            "1. Split the signal by body site before interpreting it as a "
            "facultative-aerobe increased-profile signal; "
            "oral and airway profiles can dominate aerobic/microaerophilic traits.",
            "2. In gut-only analyses, test the target taxa against inflammation markers, "
            "nitrate respiration, oxygen tolerance, virulence genes, and resistome features.",
            "3. Compare increased-profile taxa against the anaerobe-depletion guild to see "
            "whether the same conditions show reciprocal decreased- and "
            "increased-profile patterns.",
            "",
            "## Files",
            "",
            "- `guild_trait_direction_shift.tsv`",
            "- `body_site_guild_summary.tsv`",
            "- `condition_guild_summary.tsv`",
            "- `candidate_taxa_direction_summary.tsv`",
            "- `candidate_taxa_trait_presence.tsv`",
            "- `candidate_taxa_condition_direction.tsv`",
            "- `candidate_taxa_profile_occurrences.tsv`",
            "- `top_increased_guild_profiles.tsv`",
            "- `figure_methods.md`",
            "- `candidate_taxa_increased_profile_tendency.png`",
            "- `condition_increased_guild_signal.png`",
            "- `guild_trait_direction_shift.png`",
            "- `anatomical_sites_diagram.png`",
            "- `target_taxa_condition_breakdown_bar.png`",
            "",
        ]
    )
    (OUTPUT_DIR / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    study_matrix, study_metadata, study_trait_metadata = load_study_heatmap()
    scored_profiles = add_guild_scores(study_matrix, study_metadata)

    direction_summary = summarize_direction(scored_profiles)
    body_site_summary = summarize_group(scored_profiles, "body_site_group", min_profiles=10)
    condition_summary = summarize_group(scored_profiles, "condition", min_profiles=5)
    condition_summary_output = remove_excluded_condition_labels(condition_summary)
    taxa_summary = load_taxa_direction_summary()
    trait_presence = load_target_taxa_trait_counts()
    occurrences = load_target_taxa_occurrences(scored_profiles)
    condition_taxa_summary = summarize_taxa_conditions(occurrences)

    trait_definitions = study_trait_metadata.loc[
        study_trait_metadata.index.intersection(GUILD_TRAITS.keys())
    ].copy()
    trait_definitions["investigation_trait_name"] = trait_definitions.index.map(GUILD_TRAITS)

    top_increased_profiles = scored_profiles[scored_profiles["inc_dec"] == "inc"].copy()
    top_increased_profiles = top_increased_profiles.sort_values(
        ["guild_z_mean", "total_trait_signal"],
        ascending=[False, False],
    ).head(100)
    top_increased_profiles = top_increased_profiles.reset_index(names="profile_id")

    direction_summary.to_csv(OUTPUT_DIR / "guild_trait_direction_shift.tsv", sep="\t", index=False)
    body_site_summary.to_csv(OUTPUT_DIR / "body_site_guild_summary.tsv", sep="\t", index=False)
    condition_summary_output.to_csv(
        OUTPUT_DIR / "condition_guild_summary.tsv",
        sep="\t",
        index=False,
    )
    taxa_summary.to_csv(OUTPUT_DIR / "candidate_taxa_direction_summary.tsv", sep="\t")
    trait_presence.to_csv(OUTPUT_DIR / "candidate_taxa_trait_presence.tsv", sep="\t", index=False)
    condition_taxa_summary.to_csv(
        OUTPUT_DIR / "candidate_taxa_condition_direction.tsv", sep="\t", index=False
    )
    occurrences.sort_values(["target_taxon_name", "profile_id"]).to_csv(
        OUTPUT_DIR / "candidate_taxa_profile_occurrences.tsv",
        sep="\t",
        index=False,
    )
    top_increased_profiles.to_csv(
        OUTPUT_DIR / "top_increased_guild_profiles.tsv",
        sep="\t",
        index=False,
    )
    trait_definitions.to_csv(OUTPUT_DIR / "guild_trait_definitions.tsv", sep="\t")

    make_plots(taxa_summary, condition_summary, direction_summary, body_site_summary, occurrences)
    write_figure_methods()
    write_summary_report(
        direction_summary,
        body_site_summary,
        condition_summary,
        taxa_summary,
        condition_taxa_summary,
        top_increased_profiles,
    )

    print(f"Wrote investigation outputs to: {OUTPUT_DIR}")
    print(f"Study profiles analyzed: {len(scored_profiles)}")
    print(f"Target taxa occurrence rows: {len(occurrences)}")
    print(f"Condition/taxon rows: {len(condition_taxa_summary)}")


if __name__ == "__main__":
    main()
