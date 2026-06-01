"""Investigate gut anaerobe and fermentation depletion signals."""

from __future__ import annotations

import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Circle, PathPatch
from matplotlib.path import Path as MplPath

REPO_ROOT = Path(__file__).resolve().parents[1]
VIS_DIR = REPO_ROOT / "trait-visualization"
OUTPUT_DIR = VIS_DIR / "trait-profile-investigation" / "gut-anaerobe-fermentation-depletion"

GUILD_TRAITS = {
    "METPO:1000603": "anaerobic",
    "METPO:1000607": "obligately anaerobic",
    "METPO:1002005": "Fermentation",
    "METPO:2000103|GO:0030245": "cellulose catabolic process",
    "METPO:2000103|GO:0045493": "xylan catabolic process",
}

TARGET_TAXA = {
    "NCBITaxon:853": "Faecalibacterium prausnitzii",
    "NCBITaxon:39488": "Anaerobutyricum hallii",
    "NCBITaxon:39485": "Lachnospira eligens",
    "NCBITaxon:39491": "Agathobacter rectalis",
    "NCBITaxon:166486": "Roseburia intestinalis",
    "NCBITaxon:40518": "Ruminococcus bromii",
    "NCBITaxon:360807": "Roseburia inulinivorans",
}

GUT_BODY_SITE_RE = re.compile(
    r"feces|faeces|stool|colon|caecum|cecum|gut|intestin|digestive",
    re.IGNORECASE,
)
EXCLUDED_CONDITION_SUMMARY_LABELS = {"Clinical treatment"}


def read_tsv(path: Path, **kwargs: object) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", **kwargs)


def clean_text(value: object) -> str:
    text = str(value or "").strip()
    if not text or text.lower() == "nan" or text.upper() == "NA":
        return "unknown"
    return text


def is_gut_body_site(value: object) -> bool:
    return bool(GUT_BODY_SITE_RE.search(clean_text(value)))


def safe_ratio(numerator: float, denominator: float) -> float:
    if denominator == 0:
        return math.nan
    return numerator / denominator


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
    scored["gut_context"] = scored["body_site"].map(is_gut_body_site) | scored[
        "body_site_group"
    ].map(is_gut_body_site)
    scored["guild_z_mean"] = matrix[available_traits].mean(axis=1)
    scored["anaerobe_fermentation_z_mean"] = matrix[
        ["anaerobic", "obligately anaerobic", "Fermentation"]
    ].mean(axis=1)
    scored["fiber_degradation_z_mean"] = matrix[
        ["cellulose catabolic process", "xylan catabolic process"]
    ].mean(axis=1)
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
                "decreased_minus_increased_z": dec_mean - inc_mean,
            }
        )
    for score_name in [
        "guild_z_mean",
        "anaerobe_fermentation_z_mean",
        "fiber_degradation_z_mean",
    ]:
        inc_mean = metadata.loc[metadata["inc_dec"] == "inc", score_name].mean()
        dec_mean = metadata.loc[metadata["inc_dec"] == "dec", score_name].mean()
        rows.append(
            {
                "trait": score_name,
                "increased_profile_mean_z": inc_mean,
                "decreased_profile_mean_z": dec_mean,
                "decreased_minus_increased_z": dec_mean - inc_mean,
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
        records.append(
            {
                group_column: clean_text(group_value),
                "profile_count": len(group),
                "study_count": group["study_id"].nunique(),
                "increased_profile_count": len(inc),
                "decreased_profile_count": len(dec),
                "mean_guild_z": group["guild_z_mean"].mean(),
                "increased_mean_guild_z": inc["guild_z_mean"].mean(),
                "decreased_mean_guild_z": dec["guild_z_mean"].mean(),
                "decreased_minus_increased_guild_z": dec["guild_z_mean"].mean()
                - inc["guild_z_mean"].mean(),
                "decreased_mean_fiber_z": dec["fiber_degradation_z_mean"].mean(),
                "decreased_mean_anaerobe_fermentation_z": dec[
                    "anaerobe_fermentation_z_mean"
                ].mean(),
            }
        )
    return pd.DataFrame(records).sort_values(
        ["decreased_mean_guild_z", "profile_count"],
        ascending=[False, False],
    )


def remove_excluded_condition_labels(condition_summary: pd.DataFrame) -> pd.DataFrame:
    """Remove broad study-design labels from disease/injury condition summaries."""
    return condition_summary[
        ~condition_summary["condition"].isin(EXCLUDED_CONDITION_SUMMARY_LABELS)
    ].copy()


def load_taxa_direction_summary() -> pd.DataFrame:
    taxa_metadata = read_tsv(VIS_DIR / "taxa_guild_discovery_taxa_metadata.tsv", index_col=0)
    summary = taxa_metadata.loc[taxa_metadata.index.intersection(TARGET_TAXA.keys())].copy()
    summary["target_taxon_name"] = summary.index.map(TARGET_TAXA)
    summary["decreased_fraction"] = [
        safe_ratio(float(row["decreased"]), float(row["bugsigdb_profile_count"]))
        for _, row in summary.iterrows()
    ]
    summary["depletion_bias"] = [
        safe_ratio(
            float(row["decreased"]) - float(row["increased"]),
            float(row["bugsigdb_profile_count"]),
        )
        for _, row in summary.iterrows()
    ]
    summary["decreased_to_increased_ratio"] = [
        safe_ratio(float(row["decreased"]), float(row["increased"]))
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
            "decreased_fraction",
            "depletion_bias",
            "decreased_to_increased_ratio",
            "trait_count",
        ]
    ].sort_values("depletion_bias", ascending=False)


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
        "gut_context",
        "guild_z_mean",
        "anaerobe_fermentation_z_mean",
        "fiber_degradation_z_mean",
    ]
    profiles = scored_profiles[profile_columns].copy()
    merged = taxa.merge(
        profiles,
        left_on="profile_id",
        right_index=True,
        how="left",
        suffixes=("_taxa", ""),
    )
    return merged


def summarize_taxa_conditions(occurrences: pd.DataFrame) -> pd.DataFrame:
    gut = occurrences[occurrences["gut_context"].fillna(False)].copy()
    grouped = (
        gut.groupby(["target_taxon_name", "condition", "body_site_group", "direction"])[
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
    output["depletion_bias"] = (
        output["decreased_profile_count"] - output["increased_profile_count"]
    ) / output["total_profile_count"].replace(0, pd.NA)
    return output.sort_values(
        ["depletion_bias", "decreased_study_count", "total_profile_count"],
        ascending=[False, False, False],
    )


def count_target_taxa_conditions(
    occurrences: pd.DataFrame,
    direction: str,
    *,
    gut_context_only: bool,
) -> pd.Series:
    """Count condition labels for target taxon-profile occurrences."""
    selected = occurrences[occurrences["direction"] == direction].copy()
    if gut_context_only:
        selected = selected[selected["gut_context"].fillna(False)].copy()
    selected["condition"] = selected["condition"].map(clean_text)
    selected = selected[
        ~selected["condition"].isin(EXCLUDED_CONDITION_SUMMARY_LABELS)
    ].copy()
    return selected["condition"].value_counts()


def counts_with_other_remainder(counts: pd.Series, max_direct_slices: int) -> pd.Series:
    """Keep top categories and collapse the remaining long tail."""
    top_counts = counts.head(max_direct_slices).copy()
    remainder = counts.iloc[max_direct_slices:]
    if remainder.empty:
        return top_counts
    other_count = pd.Series({"other conditions": int(remainder.sum())})
    return pd.concat([top_counts, other_count])


def percentage_label(percent: float) -> str:
    """Hide tiny pie labels to keep dense pies readable."""
    return f"{percent:.1f}%" if percent >= 2 else ""


def save_condition_breakdown_plots(
    counts: pd.Series,
    *,
    pie_output_path: Path,
    histogram_output_path: Path,
    title_prefix: str,
    color: str,
) -> None:
    """Save pie and horizontal histogram plots for condition counts."""
    if counts.empty:
        return

    pie_counts = counts_with_other_remainder(counts, max_direct_slices=12)
    fig, ax = plt.subplots(figsize=(10, 7))
    wedges, _texts, autotexts = ax.pie(
        pie_counts.to_numpy(),
        autopct=percentage_label,
        startangle=90,
        counterclock=False,
        pctdistance=0.78,
    )
    for text in autotexts:
        text.set_fontsize(8)
    ax.legend(
        wedges,
        pie_counts.index,
        title="Condition",
        loc="center left",
        bbox_to_anchor=(1.0, 0.5),
        frameon=False,
        fontsize=8,
        title_fontsize=9,
    )
    ax.set_title(f"{title_prefix}: condition share")
    fig.tight_layout()
    fig.savefig(pie_output_path, bbox_inches="tight")
    plt.close(fig)

    histogram_counts = counts.head(20)
    fig, ax = plt.subplots(figsize=(10, max(5, len(histogram_counts) * 0.35)))
    ax.barh(histogram_counts.index, histogram_counts.to_numpy(), color=color)
    ax.invert_yaxis()
    ax.set_xlabel("Target taxon-profile occurrences")
    ax.set_title(f"{title_prefix}: top conditions")
    for index, value in enumerate(histogram_counts.to_numpy()):
        ax.text(value, index, f" {int(value)}", va="center", fontsize=8)
    fig.tight_layout()
    fig.savefig(histogram_output_path, bbox_inches="tight")
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
    max_sites: int = 12,
) -> None:
    """Save an anatomical-site callout diagram for body-site guild summaries."""
    if body_site_summary.empty:
        return

    plot_sites = body_site_summary[
        ~body_site_summary["body_site_group"].map(is_other_body_site)
    ].head(max_sites).copy()
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
    ax.text(
        0.02,
        0.035,
        "Bubble and bar length show profile count; color shows the "
        "direction-specific guild signal.",
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
        "generated by `scripts/investigate_gut_anaerobe_depletion.py`.",
        "",
        "## Shared Inputs",
        "",
        "Study-profile z-scores were read from "
        "`trait-visualization/study_trait_profile_clustered_heatmap_matrix.tsv` and "
        "joined to study metadata from "
        "`trait-visualization/study_trait_profile_study_metadata.tsv`. The gut "
        "anaerobe/fermentation guild score was the row-wise mean of the selected "
        "guild traits: anaerobic, obligately anaerobic, Fermentation, cellulose "
        "catabolic process, and xylan catabolic process. Gut-context profiles were "
        "defined by body-site or grouped body-site labels matching "
        "`feces`, `faeces`, `stool`, `colon`, `caecum`, `cecum`, `gut`, "
        "`intestin`, or `digestive`.",
        "",
        "The focal taxa were defined a priori in the script's `TARGET_TAXA` map "
        "using NCBI Taxonomy identifiers, which avoids ambiguity from synonym or "
        "strain-name differences. The selected taxa were "
        "Faecalibacterium prausnitzii (NCBITaxon:853), Anaerobutyricum hallii "
        "(NCBITaxon:39488), Lachnospira eligens (NCBITaxon:39485), "
        "Agathobacter rectalis (NCBITaxon:39491), Roseburia intestinalis "
        "(NCBITaxon:166486), Ruminococcus bromii (NCBITaxon:40518), and "
        "Roseburia inulinivorans (NCBITaxon:360807). They were chosen as focal "
        "examples for the gut anaerobe/fermentation depletion hypothesis because "
        "they are species-level gut anaerobe or fiber-fermentation candidates in "
        "the selected guild, and the derived taxa-trait matrix gives each selected "
        "taxon at least one of the investigation guild-trait annotations. The "
        "script writes `candidate_taxa_trait_presence.tsv` from "
        "`data/derived/taxa_trait_matrix.tsv` as the audit trail for those trait "
        "annotations.",
        "",
        "## `candidate_taxa_depletion_bias.png`",
        "",
        "This figure uses the fixed focal taxa listed above, then intersects their "
        "NCBI Taxonomy IDs with "
        "`trait-visualization/taxa_guild_discovery_taxa_metadata.tsv`. That "
        "metadata table supplies the total BugSigDB profile count and the number "
        "of profiles where each taxon was reported as increased or decreased. For "
        "each target taxon, depletion bias was calculated as `(decreased profile "
        "count - increased profile count) / total BugSigDB profile count`. The "
        "horizontal bars show this depletion-bias value for each target taxon. No "
        "condition or body-site filter is applied to this figure; it summarizes "
        "the overall BugSigDB direction tendency of the preselected focal taxa.",
        "",
        "## `condition_decreased_guild_signal.png`",
        "",
        "Gut-context study profiles were grouped by BugSigDB condition label. Groups "
        "with fewer than five profiles were excluded, and the broad label "
        "`Clinical treatment` was removed from the displayed condition summary. The "
        "figure keeps the top 12 eligible conditions with at least three decreased "
        "profiles and at least two studies, using fallback filtering if needed. Bar "
        "height represents the mean gut anaerobe/fermentation guild z-score among "
        "decreased profiles for that condition.",
        "",
        "## `guild_trait_direction_shift.png`",
        "",
        "For each selected guild trait, the script calculated the mean z-score in "
        "all decreased study profiles and the mean z-score in all increased study "
        "profiles. Bars show `decreased mean z - increased mean z`; positive values "
        "indicate traits that are higher in decreased profiles.",
        "",
        "## `anatomical_sites_diagram.png`",
        "",
        "Gut-context profiles were grouped by `body_site_group`, requiring at least "
        "10 profiles per group. The display excludes the generic `Other` body-site "
        "group. Bubble and bar length encode the number of study profiles in the "
        "body-site group. Color encodes the mean gut anaerobe/fermentation guild "
        "z-score among decreased profiles at that body site. Body-site labels were "
        "placed on a simplified, neutral front-view body outline using fixed "
        "anatomical coordinates for feces/stool, colon, caecum, oral, airway, "
        "urogenital, and unknown labels.",
        "",
        "## `target_taxa_condition_breakdown_pie.png`",
        "",
        "Target-taxon occurrences were loaded from `data/derived/taxa.tsv`, limited "
        "to the fixed focal NCBI Taxonomy IDs, deduplicated by taxon and study "
        "profile, and joined to the scored study profiles. This produces one row "
        "per focal taxon/profile occurrence with condition, body site, direction, "
        "and guild-score metadata. The pie chart counts decreased "
        "target-taxon/profile occurrences in gut-context profiles by condition, "
        "after removing `Clinical treatment`. A condition with multiple focal "
        "taxa can therefore contribute multiple occurrences; the chart is an "
        "occurrence breakdown, not a unique-study breakdown. The top 12 condition "
        "labels are shown directly; all remaining labels are collapsed into "
        "`other conditions`. Percent labels below 2% are hidden for readability.",
        "",
        "## `target_taxa_condition_breakdown_histogram.png`",
        "",
        "The histogram uses the same decreased gut-context focal-taxon/profile "
        "condition counts as the pie chart, including the same target-taxon "
        "selection, gut-context filter, deduplication rule, and `Clinical "
        "treatment` exclusion. It displays the top 20 explicit condition labels "
        "as horizontal bars, with bar labels showing occurrence counts.",
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

    fig, ax = plt.subplots(figsize=(9, 4.8))
    plot_taxa = taxa_summary.sort_values("depletion_bias")
    ax.barh(plot_taxa["target_taxon_name"], plot_taxa["depletion_bias"], color="#8c2d04")
    ax.axvline(0, color="#333333", linewidth=0.8)
    ax.set_xlabel("(decreased profiles - increased profiles) / total profiles")
    ax.set_title("Gut Anaerobe Candidate Taxa Are Depletion-Biased")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "candidate_taxa_depletion_bias.png")
    plt.close(fig)

    top_conditions = condition_summary[
        (condition_summary["decreased_profile_count"] >= 3)
        & (condition_summary["study_count"] >= 2)
    ].head(12)
    if top_conditions.empty:
        top_conditions = condition_summary[condition_summary["decreased_profile_count"] >= 3].head(
            12
        )
    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.barh(
        top_conditions["condition"][::-1],
        top_conditions["decreased_mean_guild_z"][::-1],
        color="#2563eb",
    )
    ax.set_xlabel("Mean gut anaerobe/fermentation guild z-score in decreased profiles")
    ax.set_title("Conditions With High Decreased-Profile Guild Signal")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "condition_decreased_guild_signal.png")
    plt.close(fig)

    plot_direction = direction_summary[
        direction_summary["trait"].isin(GUILD_TRAITS.values())
    ].sort_values("decreased_minus_increased_z")
    fig, ax = plt.subplots(figsize=(8, 4.2))
    ax.barh(
        plot_direction["trait"],
        plot_direction["decreased_minus_increased_z"],
        color="#047857",
    )
    ax.axvline(0, color="#333333", linewidth=0.8)
    ax.set_xlabel("decreased profile mean z - increased profile mean z")
    ax.set_title("Guild Traits Are Higher in Decreased Profiles")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "guild_trait_direction_shift.png")
    plt.close(fig)

    condition_counts = count_target_taxa_conditions(
        occurrences,
        "decreased",
        gut_context_only=True,
    )
    save_condition_breakdown_plots(
        condition_counts,
        pie_output_path=OUTPUT_DIR / "target_taxa_condition_breakdown_pie.png",
        histogram_output_path=OUTPUT_DIR / "target_taxa_condition_breakdown_histogram.png",
        title_prefix="Decreased gut anaerobe target taxa",
        color="#2563eb",
    )
    save_anatomical_sites_diagram(
        body_site_summary,
        output_path=OUTPUT_DIR / "anatomical_sites_diagram.png",
        title="Gut Anaerobe/Fermentation Depletion: Anatomical Sites",
        signal_column="decreased_mean_guild_z",
        signal_label="decreased-profile guild z",
    )


def write_summary_report(
    direction_summary: pd.DataFrame,
    body_site_summary: pd.DataFrame,
    condition_summary: pd.DataFrame,
    taxa_summary: pd.DataFrame,
    condition_taxa_summary: pd.DataFrame,
    top_depletion_profiles: pd.DataFrame,
) -> None:
    direction_rows = direction_summary[
        direction_summary["trait"].isin(GUILD_TRAITS.values())
    ].sort_values("decreased_minus_increased_z", ascending=False)
    top_conditions = condition_summary[
        (condition_summary["decreased_profile_count"] >= 3)
        & (condition_summary["study_count"] >= 2)
    ].head(8)
    if top_conditions.empty:
        top_conditions = condition_summary[condition_summary["decreased_profile_count"] >= 3].head(
            8
        )
    top_taxa = taxa_summary.head(7)
    top_taxon_conditions = condition_taxa_summary[
        (condition_taxa_summary["decreased_profile_count"] >= 3)
        & (condition_taxa_summary["decreased_study_count"] >= 2)
    ].head(10)

    lines = [
        "# Gut Anaerobe / Fermentation Depletion Investigation",
        "",
        "This investigation uses the existing study trait profile clustered heatmap matrix, "
        "study metadata, taxa guild heatmap metadata, and derived taxa/trait matrices.",
        "",
        "## Scope",
        "",
        "- Guild traits: "
        + ", ".join(f"{name} ({trait_id})" for trait_id, name in GUILD_TRAITS.items()),
        "- Target taxa: " + ", ".join(TARGET_TAXA.values()),
        "- Gut context: body-site labels matching feces/faeces/stool/colon/caecum/cecum/gut/"
        "intestin/digestive, plus matching grouped body-site labels.",
        "",
        "## Main Signal",
        "",
        "The selected guild traits are consistently higher in decreased profiles than in "
        "increased profiles, which supports the gut anaerobe/fermentation depletion hypothesis.",
        "",
        "| trait | decreased - increased z |",
        "| --- | ---: |",
    ]
    for row in direction_rows.itertuples(index=False):
        lines.append(f"| {row.trait} | {row.decreased_minus_increased_z:.3f} |")

    lines.extend(
        [
            "",
            "## Strongest Depletion-Biased Taxa",
            "",
            "| taxon | decreased / increased profiles | depletion bias |",
            "| --- | ---: | ---: |",
        ]
    )
    for row in top_taxa.itertuples(index=False):
        lines.append(
            f"| {row.target_taxon_name} | {int(row.decreased)} / {int(row.increased)} | "
            f"{row.depletion_bias:.3f} |"
        )

    lines.extend(
        [
            "",
            "## Recurrent Target-Taxon Depletions By Condition",
            "",
            "| taxon | condition | body site | "
            "decreased / increased profiles | decreased studies |",
            "| --- | --- | --- | ---: | ---: |",
        ]
    )
    for row in top_taxon_conditions.itertuples(index=False):
        lines.append(
            f"| {row.target_taxon_name} | {clean_text(row.condition)} | "
            f"{clean_text(row.body_site_group)} | "
            f"{int(row.decreased_profile_count)} / {int(row.increased_profile_count)} | "
            f"{int(row.decreased_study_count)} |"
        )

    lines.extend(
        [
            "",
            "## Gut Contexts With High Decreased-Profile Guild Signal",
            "",
            "| condition | profiles | studies | decreased mean guild z | decreased - increased z |",
            "| --- | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in top_conditions.itertuples(index=False):
        diff = row.decreased_minus_increased_guild_z
        diff_text = "" if pd.isna(diff) else f"{diff:.3f}"
        lines.append(
            f"| {row.condition} | {int(row.profile_count)} | {int(row.study_count)} | "
            f"{row.decreased_mean_guild_z:.3f} | {diff_text} |"
        )

    lines.extend(
        [
            "",
            "## Body-Site Check",
            "",
            "| body_site_group | profiles | decreased mean guild z | decreased - increased z |",
            "| --- | ---: | ---: | ---: |",
        ]
    )
    for row in body_site_summary.head(10).itertuples(index=False):
        diff = row.decreased_minus_increased_guild_z
        diff_text = "" if pd.isna(diff) else f"{diff:.3f}"
        lines.append(
            f"| {row.body_site_group} | {int(row.profile_count)} | "
            f"{row.decreased_mean_guild_z:.3f} | {diff_text} |"
        )

    lines.extend(
        [
            "",
            "## Top Decreased Profiles To Inspect",
            "",
            "| profile | condition | body site | guild z | title |",
            "| --- | --- | --- | ---: | --- |",
        ]
    )
    for row in top_depletion_profiles.head(10).itertuples(index=False):
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
            "1. Treat the target taxa as a depletion guild and test enrichment by "
            "condition/body site "
            "rather than analyzing each taxon in isolation.",
            "2. Prioritize conditions where decreased profiles have high guild scores "
            "and target taxa "
            "are repeatedly decreased across multiple studies.",
            "3. Re-run the same summaries after adding explicit SCFA/butyrate "
            "production annotations, "
            "because the current selected heatmap traits capture fermentation and "
            "fiber degradation "
            "but do not fully encode butyrate production.",
            "",
            "## Files",
            "",
            "- `guild_trait_direction_shift.tsv`",
            "- `gut_body_site_guild_summary.tsv`",
            "- `gut_condition_guild_summary.tsv`",
            "- `candidate_taxa_direction_summary.tsv`",
            "- `candidate_taxa_trait_presence.tsv`",
            "- `candidate_taxa_condition_direction.tsv`",
            "- `candidate_taxa_profile_occurrences.tsv`",
            "- `top_decreased_guild_profiles.tsv`",
            "- `figure_methods.md`",
            "- `candidate_taxa_depletion_bias.png`",
            "- `condition_decreased_guild_signal.png`",
            "- `guild_trait_direction_shift.png`",
            "- `anatomical_sites_diagram.png`",
            "- `target_taxa_condition_breakdown_pie.png`",
            "- `target_taxa_condition_breakdown_histogram.png`",
            "",
        ]
    )
    (OUTPUT_DIR / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    study_matrix, study_metadata, study_trait_metadata = load_study_heatmap()
    scored_profiles = add_guild_scores(study_matrix, study_metadata)
    gut_profiles = scored_profiles[scored_profiles["gut_context"]].copy()

    direction_summary = summarize_direction(scored_profiles)
    body_site_summary = summarize_group(gut_profiles, "body_site_group", min_profiles=10)
    condition_summary = remove_excluded_condition_labels(
        summarize_group(gut_profiles, "condition", min_profiles=5)
    )
    taxa_summary = load_taxa_direction_summary()
    trait_presence = load_target_taxa_trait_counts()
    occurrences = load_target_taxa_occurrences(scored_profiles)
    condition_taxa_summary = summarize_taxa_conditions(occurrences)

    trait_definitions = study_trait_metadata.loc[
        study_trait_metadata.index.intersection(GUILD_TRAITS.keys())
    ].copy()
    trait_definitions["investigation_trait_name"] = trait_definitions.index.map(GUILD_TRAITS)

    top_depletion_profiles = gut_profiles[gut_profiles["inc_dec"] == "dec"].copy()
    top_depletion_profiles = top_depletion_profiles.sort_values(
        ["guild_z_mean", "total_trait_signal"],
        ascending=[False, False],
    ).head(100)
    top_depletion_profiles = top_depletion_profiles.reset_index(names="profile_id")

    direction_summary.to_csv(OUTPUT_DIR / "guild_trait_direction_shift.tsv", sep="\t", index=False)
    body_site_summary.to_csv(OUTPUT_DIR / "gut_body_site_guild_summary.tsv", sep="\t", index=False)
    condition_summary.to_csv(OUTPUT_DIR / "gut_condition_guild_summary.tsv", sep="\t", index=False)
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
    top_depletion_profiles.to_csv(
        OUTPUT_DIR / "top_decreased_guild_profiles.tsv",
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
        top_depletion_profiles,
    )

    print(f"Wrote investigation outputs to: {OUTPUT_DIR}")
    print(f"Study profiles analyzed: {len(scored_profiles)}")
    print(f"Gut-context profiles analyzed: {len(gut_profiles)}")
    print(f"Target taxa occurrence rows: {len(occurrences)}")
    print(f"Condition/taxon rows: {len(condition_taxa_summary)}")


if __name__ == "__main__":
    main()
