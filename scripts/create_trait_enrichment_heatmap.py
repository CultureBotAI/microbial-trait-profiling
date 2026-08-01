#!/usr/bin/env python3
"""Create a clustered heatmap from data/derived/trait_enrichment_report.txt."""

from __future__ import annotations

import argparse
import math
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import NamedTuple

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist

DEFAULT_REPORT_PATH = Path("data/derived/trait_enrichment_results.tsv")
DEFAULT_OUTPUT_PATH = Path("data/derived/trait_enrichment_clustered_heatmap.png")
DEFAULT_MATRIX_PATH = Path("data/derived/trait_enrichment_clustered_heatmap_matrix.tsv")
METPO_MAPPING_PATH = Path("metpo_phenotype_mapping.tsv")


@dataclass(frozen=True)
class ParsedHit:
    study_id: str
    experiment_id: str
    inc_dec: str
    trait: str
    pvalue: float
    fold: float


class ProfileParts(NamedTuple):
    study_id: str
    experiment_id: str
    inc_dec: str


def load_metpo_label_lookup() -> dict[str, str]:
    if not METPO_MAPPING_PATH.exists():
        return {}

    mapping = pd.read_csv(METPO_MAPPING_PATH, sep="\t", dtype=str).fillna("")
    if not {"METPO_ID", "Description"}.issubset(mapping.columns):
        return {}

    return {
        row["METPO_ID"].strip(): row["Description"].strip()
        for _, row in mapping.iterrows()
        if row["METPO_ID"].strip()
    }


def _clean_label(value: str) -> str:
    cleaned = str(value).strip()
    return cleaned


def _format_profile_label(row: pd.Series) -> str:
    disease = _clean_label(row.get("condition", ""))
    anatomy = _clean_label(row.get("body_site", ""))
    study_ref = _clean_label(row.get("study_ref", row.get("study_id", "")))
    experiment_id = _clean_label(row.get("experiment_id", ""))
    direction = _clean_label(row.get("inc_dec", ""))

    detail_parts: list[str] = []
    if disease:
        detail_parts.append(f"Disease: {disease}")
    if anatomy:
        detail_parts.append(f"Anatomy: {anatomy}")
    detail = "\n".join(detail_parts)
    support = f"Study: {study_ref} | Exp: {experiment_id} | {direction}"
    if detail:
        return f"{detail}\n{support}"
    return support


def trait_display_name(trait: str, lookup: dict[str, str]) -> str:
    return lookup.get(trait, trait)


def _euclidean_linkage(matrix: pd.DataFrame, axis: int) -> object | None:
    values = matrix.to_numpy() if axis == 0 else matrix.to_numpy().T
    if values.shape[0] < 2 or values.shape[1] == 0:
        return None
    distances = pdist(values, metric="euclidean")
    if len(distances) == 0:
        return None
    return linkage(distances, method="average")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a clustered heatmap from the trait enrichment report."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_REPORT_PATH,
        help="Input trait enrichment report path.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT_PATH,
        help="Output heatmap image path.",
    )
    parser.add_argument(
        "--matrix-output",
        type=Path,
        default=DEFAULT_MATRIX_PATH,
        help="Output TSV path for the plotted matrix.",
    )
    parser.add_argument(
        "--top-traits",
        type=int,
        default=0,
        help="Keep the most prevalent traits across profiles. Use 0 for all traits.",
    )
    parser.add_argument(
        "--min-trait-prevalence",
        type=int,
        default=14,
        help="Keep only traits enriched in at least this many profiles.",
    )
    parser.add_argument(
        "--max-profiles",
        type=int,
        default=0,
        help="Maximum number of profiles to plot. Use 0 for all profiles.",
    )
    parser.add_argument(
        "--min-profile-hits",
        type=int,
        default=12,
        help="Keep only profiles with at least this many enriched traits in the selected matrix.",
    )
    parser.add_argument(
        "--min-profile-signal",
        type=float,
        default=50.0,
        help="Keep only profiles with signal_sum >= this value.",
    )
    parser.add_argument(
        "--all-data",
        action="store_true",
        help="Shortcut for plotting all traits and all profiles.",
    )
    parser.add_argument(
        "--skip-heatmap",
        action="store_true",
        help="Only write matrix TSV and skip heatmap rendering.",
    )
    return parser.parse_args()


def _parse_profile_header(line: str) -> ProfileParts:
    parts = line.split("\t")
    if len(parts) != 3:
        msg = f"Unparseable profile header: {line!r}"
        raise ValueError(msg)
    return ProfileParts(parts[0], parts[1], parts[2])


def _parse_fold(value: str) -> float:
    cleaned = value.strip().lower()
    if cleaned in {"inf", "+inf"}:
        return math.inf
    if cleaned == "-inf":
        return -math.inf
    return float(value)


def parse_report(path: Path) -> list[ParsedHit]:
    """Parse either the new structured TSV or the old text report."""
    if path.suffix == ".tsv":
        df = pd.read_csv(path, sep="\t", comment="#")
        hits = []
        for _, row in df.iterrows():
            hits.append(
                ParsedHit(
                    study_id=str(row["study_id"]),
                    experiment_id=str(row["experiment_id"]),
                    inc_dec=str(row["inc_dec"]),
                    trait=str(row["trait"]),
                    pvalue=float(row.get("fdr_q", row["p_value"])),
                    fold=float(row["fold"]),
                )
            )
        return hits

    hits: list[ParsedHit] = []
    current_profile: ProfileParts | None = None

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line == "---" or line.startswith("p < ") or line.startswith("FDR q < ") or line.startswith("Min Fold"):
                continue

            if line.startswith("- "):
                profile = current_profile
                if profile is None:
                    msg = f"Found a trait row before any profile header: {line!r}"
                    raise ValueError(msg)

                payload = line[2:]
                parts = payload.split(", ")
                if len(parts) not in {3, 4}:
                    msg = f"Unparseable trait row: {line!r}"
                    raise ValueError(msg)

                trait, pvalue_text = parts[0], parts[1]
                if len(parts) == 4:
                    summary = parts[3]
                    if "fold=" not in summary:
                        msg = f"Unparseable summary in trait row: {line!r}"
                        raise ValueError(msg)
                    fold_text = summary.rsplit("fold=", 1)[-1]
                else:
                    summary = parts[2]
                    if not summary.startswith("fold="):
                        msg = f"Unparseable summary in trait row: {line!r}"
                        raise ValueError(msg)
                    fold_text = summary.split("=", 1)[-1]
                hits.append(
                    ParsedHit(
                        study_id=profile.study_id,
                        experiment_id=profile.experiment_id,
                        inc_dec=profile.inc_dec,
                        trait=trait,
                        pvalue=float(pvalue_text),
                        fold=_parse_fold(fold_text),
                    )
                )
                continue

            current_profile = _parse_profile_header(line)

    if not hits:
        msg = f"No trait rows were parsed from {path}"
        raise ValueError(msg)

    return hits


def build_matrix(
    hits: list[ParsedHit],
    top_traits: int,
    min_trait_prevalence: int,
    max_profiles: int,
    min_profile_hits: int,
    min_profile_signal: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    frame = pd.DataFrame([hit.__dict__ for hit in hits])
    frame["profile_id"] = (
        frame["study_id"] + "|" + frame["experiment_id"] + "|" + frame["inc_dec"]
    )
    frame["neglog10p"] = frame["pvalue"].map(lambda value: -math.log10(max(value, 1e-300)))

    trait_prevalence = frame.groupby("trait")["profile_id"].nunique().sort_values(ascending=False)
    selected_prevalence = trait_prevalence
    if min_trait_prevalence > 0:
        selected_prevalence = selected_prevalence[selected_prevalence >= min_trait_prevalence]
    if top_traits > 0:
        selected_traits = selected_prevalence.head(top_traits).index
    else:
        selected_traits = selected_prevalence.index

    filtered = frame[frame["trait"].isin(selected_traits)].copy()
    if filtered.empty:
        msg = "No traits remain after top-trait filtering"
        raise ValueError(msg)

    label_lookup = load_metpo_label_lookup()

    matrix = (
        filtered.pivot_table(
            index="profile_id",
            columns="trait",
            values="neglog10p",
            aggfunc="max",
            fill_value=0.0,
        )
        .astype(float)
        .fillna(0.0)
    )

    # Prune empty rows and columns to reduce white space
    # This ensures that we only plot traits and profiles that have at least one significant hit
    matrix = matrix.loc[(matrix > 0).any(axis=1)]
    matrix = matrix.loc[:, (matrix > 0).any(axis=0)]

    if matrix.empty:
        msg = "No data remaining after pruning empty rows/cols"
        raise ValueError(msg)

    ordered_traits = (
        selected_prevalence.loc[matrix.columns].sort_values(ascending=False).index
    )
    matrix = matrix.loc[:, ordered_traits]
    display_columns = [trait_display_name(trait, label_lookup) for trait in matrix.columns]
    if len(set(display_columns)) != len(display_columns):
        counts: dict[str, int] = {}
        for label in display_columns:
            counts[label] = counts.get(label, 0) + 1
        display_columns = [
            f"{label} ({trait})" if counts[label] > 1 else label
            for trait, label in zip(matrix.columns, display_columns, strict=True)
        ]
    matrix.columns = display_columns

    profile_meta = (
        filtered.groupby("profile_id")
        .agg(
            study_id=("study_id", "first"),
            experiment_id=("experiment_id", "first"),
            inc_dec=("inc_dec", "first"),
            hit_count=("trait", "nunique"),
            signal_sum=("neglog10p", "sum"),
        )
        .loc[matrix.index]
    )

    # Attempt to enrich profile metadata with experiment-level fields (condition, body_site)
    try:
        experiments_path = Path("data/raw/experiments.csv")
        if experiments_path.exists():
            exdf = pd.read_csv(experiments_path, dtype=str)
            # normalize column names
            exdf.columns = [c.strip().lower().replace(" ", "_") for c in exdf.columns]
            # build keys matching parsed profile ids (e.g. study_33947803 / exp1)
            exdf = exdf.assign(
                study_key=("study_" + exdf["study"].astype(str).str.strip()),
                experiment_key=(
                    "exp"
                    + exdf["experiment_page_name"]
                    .astype(str)
                    .str.extract(r"(\d+)", expand=False)
                    .astype(str)
                ),
            )

            exdf_lookup = exdf.set_index(["study_key", "experiment_key"])

            conds = []
            sites = []
            studies = []
            for _, row in profile_meta.iterrows():
                sk = row["study_id"].strip()
                ek = row["experiment_id"].strip()
                try:
                    entry = exdf_lookup.loc[(sk, ek)]
                    # if multiple rows match, take first
                    if isinstance(entry, pd.DataFrame):
                        entry = entry.iloc[0]
                    conds.append(entry.get("condition", ""))
                    sites.append(entry.get("body_site", ""))
                    studies.append(entry.get("study", ""))
                except Exception:
                    conds.append("")
                    sites.append("")
                    studies.append("")

            profile_meta = profile_meta.assign(condition=conds, body_site=sites, study_ref=studies)
    except Exception:
        # best-effort only; silently continue if anything goes wrong
        profile_meta = profile_meta.assign(condition="", body_site="", study_ref="")

    if min_profile_hits > 0:
        keep = profile_meta["hit_count"] >= min_profile_hits
        matrix = matrix.loc[keep]
        profile_meta = profile_meta.loc[keep]

    if min_profile_signal > 0:
        keep = profile_meta["signal_sum"] >= min_profile_signal
        matrix = matrix.loc[keep]
        profile_meta = profile_meta.loc[keep]

    if matrix.empty:
        msg = "No profiles remain after profile filtering"
        raise ValueError(msg)

    if max_profiles > 0 and len(matrix) > max_profiles:
        selected_profiles = profile_meta["signal_sum"].sort_values(ascending=False).head(
            max_profiles
        ).index
        matrix = matrix.loc[selected_profiles]
        profile_meta = profile_meta.loc[selected_profiles]

    matrix = matrix.sort_index()
    profile_meta = profile_meta.loc[matrix.index]
    # Map trait IDs (e.g. METPO:XXXX) to human phenotype names when a mapping
    # file exists. If multiple traits map to the same phenotype, disambiguate by
    # appending the original trait id in brackets so column names remain unique.
    try:
        mapping_path = Path("metpo_phenotype_mapping.tsv")
        if mapping_path.exists():
            mapdf = pd.read_csv(mapping_path, sep="\t", dtype=str)
            # expected columns: METPO_ID, Description
            if "METPO_ID" in mapdf.columns and "Description" in mapdf.columns:
                lookup = dict(
                    zip(mapdf["METPO_ID"], mapdf["Description"], strict=False)
                )
                orig_cols = list(matrix.columns)
                mapped = [lookup.get(c, c) for c in orig_cols]

                counts = Counter(mapped)
                if any(v > 1 for v in counts.values()):
                    new_cols = []
                    # disambiguate duplicates by appending original trait id
                    for orig, lab in zip(orig_cols, mapped, strict=False):
                        if counts[lab] == 1:
                            new_cols.append(lab)
                        else:
                            new_cols.append(f"{lab} [{orig}]")
                    matrix.columns = new_cols
                else:
                    matrix.columns = mapped
    except Exception:
        # best-effort only; do not fail plotting if mapping unavailable
        pass
    return matrix, profile_meta


def save_heatmap(matrix: pd.DataFrame, profile_meta: pd.DataFrame, output_path: Path) -> None:
    sns.set_theme(style="white")
    row_colors = profile_meta["inc_dec"].map({"inc": "#c0392b", "dec": "#2874a6"}).fillna(
        "#7f8c8d"
    )
    row_count = len(matrix)
    col_count = len(matrix.columns)
    show_row_labels = row_count <= 500
    show_col_labels = col_count <= 200

    # Clustering scales quadratically. For large full-dataset plots, disable it.
    row_cluster = row_count <= 2000
    col_cluster = col_count <= 1200

    nonzero_values = matrix.to_numpy().ravel()
    nonzero_values = nonzero_values[nonzero_values > 0]
    # Cap vmax to prevent extreme outliers from washing out the heatmap.
    # We'll use the 90th percentile for the top of the color scale.
    vmax = float(pd.Series(nonzero_values).quantile(0.90)) if len(nonzero_values) else 1.0
    vmax = min(max(vmax, 5.0), 40.0)

    # Use a mask for 0 values so they don't get colored.
    mask = matrix == 0

    row_linkage = _euclidean_linkage(matrix, axis=0) if row_cluster else None
    col_linkage = _euclidean_linkage(matrix, axis=1) if col_cluster else None

    fig_width = min(48, max(20, col_count * 0.08))
    fig_height = min(48, max(20, row_count * 0.02))

    # Set background color for masked cells to light gray so hits pop out
    plt.rcParams["axes.facecolor"] = "#f0f0f0"

    cluster = sns.clustermap(
        matrix,
        mask=mask,
        cmap="Blues",  # Traditional blue sequential palette
        method="average",
        metric="euclidean",
        row_colors=row_colors,
        linewidths=0,
        xticklabels=show_col_labels,
        yticklabels=show_row_labels,
        vmin=2.0,
        vmax=vmax,
        row_linkage=row_linkage,
        col_linkage=col_linkage,
        row_cluster=row_cluster,
        col_cluster=col_cluster,
        figsize=(fig_width, fig_height),
        cbar_pos=(0.02, 0.82, 0.02, 0.12),
        dendrogram_ratio=(0.14, 0.12),
        colors_ratio=(0.03, 0.03),
    )

    cluster.ax_heatmap.set_xlabel("Trait")
    cluster.ax_heatmap.set_ylabel("Study profile")
    cluster.ax_heatmap.set_title(
        "BugSigDB Trait Enrichment Clustered Heatmap\n"
        "Cell values are -log10(p) for significant traits in each "
        "study|experiment|direction profile",
        pad=18,
    )
    cluster.ax_heatmap.xaxis.label.set_size(11)
    cluster.ax_heatmap.yaxis.label.set_size(11)

    if show_col_labels:
        for label in cluster.ax_heatmap.get_xticklabels():
            label.set_rotation(60)
            label.set_ha("right")
            label.set_fontsize(7)

    if show_row_labels:
        # Build richer row labels including disease and anatomy when available.
        # We MUST use the reordered indices if clustering occurred to ensure labels match rows.
        reordered_indices = (
            cluster.dendrogram_row.reordered_ind if row_cluster else range(len(matrix))
        )
        reordered_meta = profile_meta.iloc[reordered_indices]
        new_labels = [_format_profile_label(r) for _, r in reordered_meta.iterrows()]

        cluster.ax_heatmap.set_yticklabels(new_labels)
        # Move y-axis tick labels to the right side of the heatmap
        cluster.ax_heatmap.yaxis.tick_right()
        cluster.ax_heatmap.yaxis.set_label_position("right")
        for label in cluster.ax_heatmap.get_yticklabels():
            label.set_rotation(0)
            # When labels are on the right, left-align them relative to the ticks
            label.set_ha("left")
            label.set_va("center")
            label.set_fontsize(8)
            label.set_clip_on(False)
        # Slight padding so labels don't touch the heatmap
        cluster.ax_heatmap.tick_params(axis="y", pad=6)

    from matplotlib.patches import Patch

    legend_handles = [
        Patch(facecolor="#c0392b", edgecolor="none", label="inc"),
        Patch(facecolor="#2874a6", edgecolor="none", label="dec"),
    ]
    cluster.fig.legend(
        handles=legend_handles,
        title="Row annotation",
        loc="upper right",
        bbox_to_anchor=(0.985, 0.985),
        frameon=False,
        fontsize=9,
        title_fontsize=10,
    )
    # Provide extra space on the right for y-axis labels moved to the right
    try:
        cluster.fig.subplots_adjust(top=0.94, right=0.88)
    except Exception:
        cluster.fig.subplots_adjust(top=0.94)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    # Lower DPI for very large matrices to keep runtime/memory manageable.
    cell_count = row_count * col_count
    dpi = 120 if cell_count > 2_000_000 else 300
    cluster.fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(cluster.fig)


def main() -> int:
    args = parse_args()
    if args.all_data:
        args.top_traits = 0
        args.max_profiles = 0

    hits = parse_report(args.input)
    matrix, profile_meta = build_matrix(
        hits,
        args.top_traits,
        args.min_trait_prevalence,
        args.max_profiles,
        args.min_profile_hits,
        args.min_profile_signal,
    )

    matrix_to_save = profile_meta.join(matrix)
    matrix_to_save.to_csv(args.matrix_output, sep="\t", index_label="profile_id")

    if not args.skip_heatmap:
        save_heatmap(matrix, profile_meta, args.output)

    print(f"Input rows: {len(hits)}")
    print(f"Profiles plotted: {len(matrix)}")
    print(f"Traits plotted: {len(matrix.columns)}")
    print(f"All-data mode: {args.all_data}")
    print(f"Skipped heatmap: {args.skip_heatmap}")
    print(f"Min trait prevalence: {args.min_trait_prevalence}")
    print(f"Min profile hits: {args.min_profile_hits}")
    print(f"Min profile signal: {args.min_profile_signal}")
    print(f"Wrote matrix: {args.matrix_output}")
    if not args.skip_heatmap:
        print(f"Wrote heatmap: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())