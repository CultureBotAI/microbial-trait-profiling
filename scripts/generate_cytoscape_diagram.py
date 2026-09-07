#!/usr/bin/env python3
"""Generate Cytoscape-inspired tripartite network diagrams for increased and decreased profiles."""

import argparse
from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.path import Path as MplPath

# Default file paths
DEFAULT_ENRICHMENT_PATH = Path("data/derived/trait_enrichment_results.tsv")
DEFAULT_EXPERIMENTS_PATH = Path("data/derived/experiments.tsv")
DEFAULT_RAW_EXPERIMENTS_PATH = Path("data/raw/full_dump.csv")
DEFAULT_TRAIT_META_PATH = Path("data/derived/taxa_trait_matrix_columns.tsv")
DEFAULT_PHENOTYPE_MAPPING_PATH = Path("metpo_phenotype_mapping.tsv")
DEFAULT_OUTPUT_DIR = Path("data/derived/Cytoscape_Diagrams")

# Six canonical traits matching the reference
DEFAULT_TRAITS = [
    "METPO:1000608",  # facultatively aerobic
    "METPO:2000103|GO:0030245",  # cellulose catabolic process
    "METPO:1000668",  # coccus shaped
    "METPO:1000683",  # sphere shaped
    "METPO:2000103|pathways:arsenate_reduction",  # arsenate reduction
    "METPO:1000623",  # moderately halophilic
]

TRAIT_COLORS = [
    "#0072b2",  # blue: facultatively aerobic
    "#d55e00",  # vermilion/orange: cellulose catabolic process
    "#009e73",  # bluish green: coccus shaped
    "#cc79a7",  # reddish purple/pink: sphere shaped
    "#e69f00",  # yellow-orange: arsenate reduction
    "#562b80",  # deep purple: moderately halophilic
]

BG_COLOR = "#d9d9d9"

# High-level anatomical display groups
BODY_SITE_MAP = {
    "feces": "Feces",
    "stool": "Feces",
    "saliva": "Oral cavity / dental",
    "oral cavity": "Oral cavity / dental",
    "tongue": "Oral cavity / dental",
    "dental plaque": "Oral cavity / dental",
    "supragingival plaque": "Oral cavity / dental",
    "subgingival": "Oral cavity / dental",
    "oropharynx": "Oral cavity / dental",
    "pharynx": "Oral cavity / dental",
    "gingival sulcus": "Oral cavity / dental",
    "vagina": "Vaginal / reproductive",
    "cervix": "Vaginal / reproductive",
    "cervicovaginal": "Vaginal / reproductive",
    "endometrium": "Vaginal / reproductive",
    "semen": "Vaginal / reproductive",
    "penis": "Vaginal / reproductive",
    "nasal cavity": "Upper airway / nasal",
    "nasopharynx": "Upper airway / nasal",
    "throat": "Upper airway / nasal",
    "nose": "Upper airway / nasal",
    "upper respiratory tract": "Upper airway / nasal",
    "colon": "Lower GI tract",
    "duodenum": "Lower GI tract",
    "ileum": "Lower GI tract",
    "rectum": "Lower GI tract",
    "intestine": "Lower GI tract",
    "cecum": "Lower GI tract",
    "small intestine": "Lower GI tract",
    "stomach": "Stomach / esophagus",
    "mucosa of body of stomach": "Stomach / esophagus",
    "gastric mucosa": "Stomach / esophagus",
    "esophagus": "Stomach / esophagus",
    "eye": "Ocular / eye",
    "conjunctiva": "Ocular / eye",
    "tear film": "Ocular / eye",
    "ocular surface": "Ocular / eye",
    "meconium": "Birth-associated samples",
    "placenta": "Birth-associated samples",
    "amniotic fluid": "Birth-associated samples",
    "breast milk": "Birth-associated samples",
    "cord blood": "Birth-associated samples",
    "lung": "Lower airway / lung",
    "bronchus": "Lower airway / lung",
    "lower respiratory tract": "Lower airway / lung",
    "sputum": "Lower airway / lung",
    "skin": "Skin",
    "scalp": "Skin",
    "forehead": "Skin",
    "colonic mucosa": "Intestinal / GI mucosa",
    "duodenal mucosa": "Intestinal / GI mucosa",
    "intestinal mucosa": "Intestinal / GI mucosa",
}


def map_body_site(raw_site: str | None) -> str:
    """Map a raw BugSigDB body site to high-level anatomical category."""
    if not raw_site or pd.isna(raw_site):
        return "Other"
    raw_lower = str(raw_site).strip().lower()
    for k in sorted(BODY_SITE_MAP.keys(), key=len, reverse=True):
        if k in raw_lower:
            return BODY_SITE_MAP[k]
    return raw_site.strip()


def load_trait_labels(trait_meta_path: Path, mapping_path: Path) -> dict[str, str]:
    """Load human-readable trait labels from column metadata or mapping files."""
    labels = {
        "METPO:1000608": "facultatively\naerobic",
        "METPO:2000103|GO:0030245": "cellulose\ncatabolic\nprocess",
        "METPO:1000668": "coccus shaped",
        "METPO:1000683": "sphere shaped",
        "METPO:2000103|pathways:arsenate_reduction": "arsenate_reduction",
        "METPO:1000623": "moderately\nhalophilic",
    }
    if trait_meta_path.exists():
        try:
            tm = pd.read_csv(trait_meta_path, sep="\t")
            for _, r in tm.iterrows():
                cid = str(r["column_id"])
                lbl = str(r.get("column_label", cid))
                lbl = (
                    lbl.split(" [METPO")[0]
                    .split(" METPO:")[0]
                    .split(" [GO:")[0]
                    .split(" GO:")[0]
                    .replace("capable of: ", "")
                )
                if "arsenate reduction" in lbl.lower() or "arsenate_reduction" in lbl.lower():
                    lbl = "arsenate_reduction"
                elif len(lbl.split()) > 2 and "\n" not in lbl:
                    words = lbl.split()
                    mid = len(words) // 2
                    lbl = " ".join(words[:mid]) + "\n" + " ".join(words[mid:])
                labels[cid] = lbl
        except Exception:
            pass
    if mapping_path.exists():
        try:
            mp = pd.read_csv(mapping_path, sep="\t", header=None)
            for _, r in mp.iterrows():
                cid = str(r[0])
                lbl = (
                    str(r[1])
                    .split(" [METPO")[0]
                    .split(" METPO:")[0]
                    .split(" [GO:")[0]
                    .split(" GO:")[0]
                )
                if "arsenate reduction" in lbl.lower() or "arsenate_reduction" in lbl.lower():
                    lbl = "arsenate_reduction"
                elif len(lbl.split()) > 2 and "\n" not in lbl:
                    words = lbl.split()
                    mid = len(words) // 2
                    lbl = " ".join(words[:mid]) + "\n" + " ".join(words[mid:])
                labels[cid] = lbl
        except Exception:
            pass
    # Force canonical labels
    labels["METPO:2000103|GO:0030245"] = "cellulose\ncatabolic\nprocess"
    labels["METPO:2000103|pathways:arsenate_reduction"] = "arsenate_reduction"
    return labels


def get_cubic_bezier_path(x0: float, y0: float, x1: float, y1: float, width: float) -> MplPath:
    """Create a horizontal cubic Bezier ribbon of specified vertical width."""
    hw = width / 2.0
    dx = (x1 - x0) * 0.5

    verts = [
        (x0, y0 - hw),
        (x0 + dx, y0 - hw),
        (x1 - dx, y1 - hw),
        (x1, y1 - hw),
        (x1, y1 + hw),
        (x1 - dx, y1 + hw),
        (x0 + dx, y0 + hw),
        (x0, y0 + hw),
        (x0, y0 - hw),
    ]
    codes = [
        MplPath.MOVETO,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.LINETO,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CLOSEPOLY,
    ]
    return MplPath(verts, codes)


def plot_cytoscape_diagram(
    enrich_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    selected_traits: list[str],
    trait_labels: dict[str, str],
    direction: str,
    output_png: Path,
    output_svg: Path | None = None,
    min_profiles: int = 3,
    min_studies: int = 2,
) -> None:
    """Generate a Cytoscape-inspired tripartite network diagram."""
    # 1. Filter enrichment results by selected traits and direction
    df = enrich_df[
        (enrich_df["trait"].isin(selected_traits)) & (enrich_df["inc_dec"] == direction)
    ].copy()

    if df.empty:
        print(f"Warning: No data found for direction={direction} with selected traits.")
        return

    # Merge metadata if needed
    if ("condition" not in df.columns or "body_site" not in df.columns) and not metadata_df.empty:
        df = df.merge(
            metadata_df[["study_id", "experiment_id", "condition", "body_site"]],
            on=["study_id", "experiment_id"],
            how="left",
        )

    df["display_site"] = df["body_site"].apply(map_body_site)
    df["display_condition"] = df["condition"].fillna("Unknown Condition")

    # 2. Trait color mapping
    trait_color_map = {
        t: TRAIT_COLORS[i % len(TRAIT_COLORS)] for i, t in enumerate(selected_traits)
    }

    # 3. Aggregate Condition-to-Trait relationships
    cond_trait_grp = (
        df.groupby(["display_condition", "trait"])
        .agg(
            n_profiles=("study_id", "count"),
            n_studies=("study_id", "nunique"),
        )
        .reset_index()
    )
    cond_trait_grp = cond_trait_grp[
        (cond_trait_grp["n_profiles"] >= min_profiles)
        & (cond_trait_grp["n_studies"] >= min_studies)
    ]

    # Aggregate Trait-to-BodySite relationships
    trait_site_grp = (
        df.groupby(["trait", "display_site"])
        .agg(
            n_profiles=("study_id", "count"),
            n_studies=("study_id", "nunique"),
        )
        .reset_index()
    )
    trait_site_grp = trait_site_grp[
        (trait_site_grp["n_profiles"] >= min_profiles)
        & (trait_site_grp["n_studies"] >= min_studies)
    ]

    # Qualifying conditions and sites
    qualifying_conditions = cond_trait_grp["display_condition"].unique()
    qualifying_sites = trait_site_grp["display_site"].unique()

    if len(qualifying_conditions) == 0 or len(qualifying_sites) == 0:
        # Relax constraints if needed
        cond_trait_grp = (
            df.groupby(["display_condition", "trait"])
            .agg(n_profiles=("study_id", "count"), n_studies=("study_id", "nunique"))
            .reset_index()
        )
        trait_site_grp = (
            df.groupby(["trait", "display_site"])
            .agg(n_profiles=("study_id", "count"), n_studies=("study_id", "nunique"))
            .reset_index()
        )
        qualifying_conditions = cond_trait_grp["display_condition"].unique()
        qualifying_sites = trait_site_grp["display_site"].unique()

    # Condition statistics
    cond_stats = (
        cond_trait_grp.groupby("display_condition")
        .agg(
            total_links=("trait", "count"),
            total_traits=("trait", "nunique"),
            total_n=("n_profiles", "sum"),
        )
        .reset_index()
    )
    cond_stats = cond_stats.sort_values(
        by=["total_n", "total_links", "total_traits"], ascending=False
    )
    sorted_conditions = cond_stats["display_condition"].tolist()

    # Body site statistics
    site_stats = (
        trait_site_grp.groupby("display_site")
        .agg(
            total_links=("trait", "count"),
            total_traits=("trait", "nunique"),
            total_n=("n_profiles", "sum"),
        )
        .reset_index()
    )
    site_stats = site_stats.sort_values(
        by=["total_n", "total_links", "total_traits"], ascending=False
    )
    sorted_sites = site_stats["display_site"].tolist()

    # Trait statistics
    trait_stats = {}
    for t in selected_traits:
        t_sub = df[df["trait"] == t]
        trait_stats[t] = {
            "n": len(t_sub),
            "conditions": t_sub["display_condition"].nunique(),
            "sites": t_sub["display_site"].nunique(),
        }

    # 4. Figure setup
    fig = plt.figure(figsize=(24, 15), facecolor="white")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")

    # Header
    dir_title = "decreased profiles" if direction == "dec" else "increased profiles"
    ax.text(
        3,
        97,
        f"Six shared traits Cytoscape-inspired graph: {dir_title}",
        fontsize=22,
        fontweight="bold",
        fontfamily="sans-serif",
        color="#111111",
        va="top",
        ha="left",
    )
    subtitle = (
        f"Graph view from displayed edges ({direction.capitalize()}, individual conditions). "
        "Left nodes are conditions using raw individual condition labels;\n"
        "center nodes are the six trait hubs; right nodes are body-site groups. "
        "Edges are aggregated condition-to-trait and trait-to-body-site associations.\n"
        "n = positive-z profile memberships."
    )
    ax.text(
        3,
        93.5,
        subtitle,
        fontsize=11,
        fontfamily="sans-serif",
        color="#444444",
        va="top",
        ha="left",
    )

    # Column titles
    ax.text(
        15,
        86,
        "Condition nodes",
        fontsize=15,
        fontweight="medium",
        fontfamily="sans-serif",
        color="#333333",
        ha="center",
    )
    ax.text(
        50,
        86,
        "Six trait hubs",
        fontsize=15,
        fontweight="medium",
        fontfamily="sans-serif",
        color="#333333",
        ha="center",
    )
    ax.text(
        85,
        86,
        "Body-site nodes",
        fontsize=15,
        fontweight="medium",
        fontfamily="sans-serif",
        color="#333333",
        ha="center",
    )

    # 5. Compute Node Layout Coordinates
    # Trait Hubs (2 columns x 3 rows centered at X=45.5 and X=56.5)
    trait_coords = {}
    trait_cols = [45.5, 56.5]
    trait_rows = [68.0, 50.0, 32.0]
    for i, t in enumerate(selected_traits[:6]):
        c_idx = i % 2
        r_idx = i // 2
        trait_coords[t] = (trait_cols[c_idx], trait_rows[r_idx])

    # Condition nodes (Left column at X=15)
    n_conds = len(sorted_conditions)
    y_top_c = 82.0
    y_bot_c = 15.0
    cond_y_spacing = (y_top_c - y_bot_c) / max(n_conds - 1, 1)
    cond_coords = {}
    for i, c in enumerate(sorted_conditions):
        y = y_top_c - i * cond_y_spacing
        cond_coords[c] = (15.0, y)

    # Body-site nodes (Right column at X=85)
    n_sites = len(sorted_sites)
    y_top_s = 80.0
    y_bot_s = 16.0
    site_y_spacing = (y_top_s - y_bot_s) / max(n_sites - 1, 1)
    site_coords = {}
    for i, s in enumerate(sorted_sites):
        y = y_top_s - i * site_y_spacing
        site_coords[s] = (85.0, y)

    # 6. Edge Width Scaling Function
    max_profile_n = max(
        cond_trait_grp["n_profiles"].max() if not cond_trait_grp.empty else 1,
        trait_site_grp["n_profiles"].max() if not trait_site_grp.empty else 1,
        10,
    )

    def calc_width(n_val: float) -> float:
        # Scale width between 0.2 and 1.8 coordinate units
        norm = np.sqrt(max(n_val, 1) / max_profile_n)
        return float(np.clip(0.25 + norm * 1.55, 0.25, 2.0))

    # 7. Draw Edges (Beziers)
    # Sort edges thin-to-thick so thick edges don't completely bury thin ones
    cond_trait_grp_sorted = cond_trait_grp.sort_values(by="n_profiles", ascending=True)
    for _, row in cond_trait_grp_sorted.iterrows():
        c = row["display_condition"]
        t = row["trait"]
        n_p = row["n_profiles"]
        if c not in cond_coords or t not in trait_coords:
            continue
        x0, y0 = cond_coords[c]
        x1, y1 = trait_coords[t]
        # Start from right edge of condition oval (x0+1.8) & end at left of trait ellipse (x1-3.8)
        w = calc_width(n_p)
        color = trait_color_map[t]
        bezier_p = get_cubic_bezier_path(x0 + 1.8, y0, x1 - 3.8, y1, w)
        patch = patches.PathPatch(
            bezier_p,
            facecolor=color,
            edgecolor=color,
            alpha=0.35,
            linewidth=0.2,
            zorder=2,
        )
        ax.add_patch(patch)

    trait_site_grp_sorted = trait_site_grp.sort_values(by="n_profiles", ascending=True)
    for _, row in trait_site_grp_sorted.iterrows():
        t = row["trait"]
        s = row["display_site"]
        n_p = row["n_profiles"]
        if t not in trait_coords or s not in site_coords:
            continue
        x0, y0 = trait_coords[t]
        x1, y1 = site_coords[s]
        # Start from right edge of trait ellipse (x0 + 3.8) and end at left of site oval (x1 - 1.8)
        w = calc_width(n_p)
        color = trait_color_map[t]
        bezier_p = get_cubic_bezier_path(x0 + 3.8, y0, x1 - 1.8, y1, w)
        patch = patches.PathPatch(
            bezier_p,
            facecolor=color,
            edgecolor=color,
            alpha=0.35,
            linewidth=0.2,
            zorder=2,
        )
        ax.add_patch(patch)

    # 8. Render Nodes
    # Left Condition Nodes
    for c, (x, y) in cond_coords.items():
        # Node ellipse
        c_patch = patches.Ellipse(
            (x, y),
            width=3.6,
            height=min(cond_y_spacing * 0.75, 1.8),
            facecolor="#fcfcfc",
            edgecolor="#555555",
            linewidth=1.2,
            zorder=4,
        )
        ax.add_patch(c_patch)

        # Label & stats on the left of node
        c_sub = cond_stats[cond_stats["display_condition"] == c].iloc[0]
        # Truncate long condition label if necessary
        disp_name = c if len(c) <= 28 else c[:25] + "..."
        l_cnt = c_sub["total_links"]
        t_cnt = c_sub["total_traits"]
        n_cnt = c_sub["total_n"]
        info_text = f"links={l_cnt}, traits={t_cnt}, n={n_cnt}"
        ax.text(
            x - 2.2,
            y + 0.3,
            disp_name,
            fontsize=8.5,
            fontweight="normal",
            fontfamily="sans-serif",
            color="#000000",
            ha="right",
            va="center",
            zorder=5,
        )
        ax.text(
            x - 2.2,
            y - 0.45,
            info_text,
            fontsize=7.5,
            fontfamily="sans-serif",
            color="#222222",
            ha="right",
            va="center",
            zorder=5,
        )

    # Right Body Site Nodes
    for s, (x, y) in site_coords.items():
        s_patch = patches.Ellipse(
            (x, y),
            width=3.6,
            height=min(site_y_spacing * 0.75, 1.8),
            facecolor="#fcfcfc",
            edgecolor="#555555",
            linewidth=1.2,
            zorder=4,
        )
        ax.add_patch(s_patch)

        s_sub = site_stats[site_stats["display_site"] == s].iloc[0]
        l_cnt = s_sub["total_links"]
        t_cnt = s_sub["total_traits"]
        n_cnt = s_sub["total_n"]
        info_text = f"links={l_cnt}, traits={t_cnt}, n={n_cnt}"
        ax.text(
            x + 2.2,
            y + 0.3,
            s,
            fontsize=8.5,
            fontweight="normal",
            fontfamily="sans-serif",
            color="#000000",
            ha="left",
            va="center",
            zorder=5,
        )
        ax.text(
            x + 2.2,
            y - 0.45,
            info_text,
            fontsize=7.5,
            fontfamily="sans-serif",
            color="#222222",
            ha="left",
            va="center",
            zorder=5,
        )

    # Center Six Trait Hubs
    for t, (x, y) in trait_coords.items():
        color = trait_color_map[t]
        t_patch = patches.Ellipse(
            (x, y),
            width=8.0,
            height=4.8,
            facecolor=color,
            edgecolor="#111111",
            linewidth=2.0,
            zorder=4,
        )
        ax.add_patch(t_patch)

        # Label inside
        lbl = trait_labels.get(t, t)
        # Determine text contrast color
        rgb = mcolors.to_rgb(color)
        luminance = 0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2]
        txt_color = "#ffffff" if luminance < 0.55 else "#111111"

        ax.text(
            x,
            y,
            lbl,
            fontsize=10.5,
            fontweight="bold",
            fontfamily="sans-serif",
            color=txt_color,
            ha="center",
            va="center",
            zorder=5,
        )

        # Stats below hub
        st = trait_stats.get(t, {"n": 0, "conditions": 0, "sites": 0})
        stats_text = f"n={st['n']}\nconditions={st['conditions']}, sites={st['sites']}"
        ax.text(
            x,
            y - 3.4,
            stats_text,
            fontsize=8.5,
            fontfamily="sans-serif",
            color="#333333",
            ha="center",
            va="top",
            zorder=5,
        )

    # 9. Bottom Legend
    ax.text(
        50,
        5.5,
        "Edge widths; trait colors are central hubs",
        fontsize=10,
        fontweight="bold",
        fontfamily="sans-serif",
        color="#222222",
        ha="center",
    )

    legend_n_vals = [
        int(max(1, max_profile_n * 0.03)),
        int(max(2, max_profile_n * 0.15)),
        int(max(3, max_profile_n * 0.85)),
    ]
    legend_xs = [42.0, 50.0, 59.0]
    for n_val, lx in zip(legend_n_vals, legend_xs, strict=False):
        lw = calc_width(n_val) * 3.0
        ax.plot(
            [lx - 1.8, lx - 0.4],
            [3.8, 3.8],
            color="#444444",
            linewidth=lw,
            solid_capstyle="round",
        )
        ax.text(
            lx - 0.1,
            3.8,
            f"edge n={n_val}",
            fontsize=9,
            fontfamily="sans-serif",
            color="#222222",
            va="center",
            ha="left",
        )

    # Footer note
    ax.text(
        0.5,
        1.0,
        f"Static graph PNG generated from {direction}_individual_conditions. "
        "Edge data aggregates the displayed links into condition-trait "
        "and trait-body-site graph edges.",
        fontsize=8.5,
        fontfamily="sans-serif",
        color="#777777",
        ha="left",
        va="bottom",
    )

    # Save outputs
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    print(f"Saved Cytoscape diagram to {output_png}")

    if output_svg:
        output_svg.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_svg, format="svg", bbox_inches="tight")
        print(f"Saved Cytoscape SVG diagram to {output_svg}")

    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate Cytoscape-inspired tripartite network diagrams."
    )
    parser.add_argument(
        "--enrichment",
        type=Path,
        default=DEFAULT_ENRICHMENT_PATH,
        help="Path to trait enrichment results TSV",
    )
    parser.add_argument(
        "--experiments",
        type=Path,
        default=DEFAULT_EXPERIMENTS_PATH,
        help="Path to experiments TSV",
    )
    parser.add_argument(
        "--raw-experiments",
        type=Path,
        default=DEFAULT_RAW_EXPERIMENTS_PATH,
        help="Path to raw full_dump.csv",
    )
    parser.add_argument(
        "--trait-meta",
        type=Path,
        default=DEFAULT_TRAIT_META_PATH,
        help="Path to taxa_trait_matrix_columns.tsv",
    )
    parser.add_argument(
        "--phenotype-mapping",
        type=Path,
        default=DEFAULT_PHENOTYPE_MAPPING_PATH,
        help="Path to metpo_phenotype_mapping.tsv",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory to save output diagrams",
    )
    parser.add_argument(
        "--direction",
        type=str,
        choices=["inc", "dec", "both"],
        default="both",
        help="Enrichment direction (inc, dec, or both)",
    )
    parser.add_argument(
        "--min-profiles",
        type=int,
        default=3,
        help="Minimum profile count per edge",
    )
    parser.add_argument(
        "--min-studies",
        type=int,
        default=2,
        help="Minimum study count per edge",
    )
    args = parser.parse_args()

    if not args.enrichment.exists():
        raise FileNotFoundError(f"Enrichment file not found: {args.enrichment}")

    enrich_df = pd.read_csv(args.enrichment, sep="\t", comment="#")

    # Load experiments metadata
    metadata_df = pd.DataFrame()
    if args.experiments.exists():
        metadata_df = pd.read_csv(args.experiments, sep="\t")
    elif args.raw_experiments.exists():
        raw_df = pd.read_csv(args.raw_experiments)
        meta_cols = {}
        if "Study" in raw_df.columns:
            meta_cols["Study"] = "study_id"
        if "Experiment" in raw_df.columns:
            meta_cols["Experiment"] = "experiment_id"
        if "Condition" in raw_df.columns:
            meta_cols["Condition"] = "condition"
        if "Body site" in raw_df.columns:
            meta_cols["Body site"] = "body_site"
        metadata_df = raw_df[list(meta_cols.keys())].rename(columns=meta_cols).drop_duplicates()

    trait_labels = load_trait_labels(args.trait_meta, args.phenotype_mapping)

    directions = ["dec", "inc"] if args.direction == "both" else [args.direction]

    for dir_code in directions:
        dir_name = "decreased" if dir_code == "dec" else "increased"
        # Auto-rank top 6 traits for this direction if needed
        sub_df = enrich_df[enrich_df["inc_dec"] == dir_code]
        top_traits_in_data = sub_df["trait"].value_counts().head(6).index.tolist()
        # Use default canonical traits if available in data, else use top 6
        selected_traits = [t for t in DEFAULT_TRAITS if t in sub_df["trait"].values]
        for t in top_traits_in_data:
            if len(selected_traits) >= 6:
                break
            if t not in selected_traits:
                selected_traits.append(t)

        out_png = args.output_dir / f"cytoscape_six_traits_{dir_name}.png"
        out_svg = args.output_dir / f"cytoscape_six_traits_{dir_name}.svg"

        plot_cytoscape_diagram(
            enrich_df=enrich_df,
            metadata_df=metadata_df,
            selected_traits=selected_traits[:6],
            trait_labels=trait_labels,
            direction=dir_code,
            output_png=out_png,
            output_svg=out_svg,
            min_profiles=args.min_profiles,
            min_studies=args.min_studies,
        )


if __name__ == "__main__":
    main()
