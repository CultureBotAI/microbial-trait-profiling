#!/usr/bin/env python3
"""Generate grouped six-trait Sankey diagrams for increased and decreased profiles."""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.path import Path as MplPath

# Default file paths
DEFAULT_ENRICHMENT_PATH = Path("data/derived/trait_enrichment_results.tsv")
DEFAULT_RANKED_PAIRS_PATH = Path("data/derived/trait_enrichment_trait_study_pairs_ranked.tsv")
DEFAULT_EXPERIMENTS_PATH = Path("data/derived/experiments.tsv")
DEFAULT_RAW_EXPERIMENTS_PATH = Path("data/raw/full_dump.csv")
DEFAULT_TRAIT_META_PATH = Path("data/derived/taxa_trait_matrix_columns.tsv")
DEFAULT_PHENOTYPE_MAPPING_PATH = Path("metpo_phenotype_mapping.tsv")
DEFAULT_OUTPUT_DIR = Path("data/derived/Sankey_Diagrams")

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
    "#1f77b4",  # blue: facultatively aerobic
    "#d95f02",  # orange: cellulose catabolic process
    "#2ca02c",  # green: coccus shaped
    "#e377c2",  # pink: sphere shaped
    "#e6ab02",  # yellow-orange: arsenate reduction
    "#9467bd",  # purple: moderately halophilic
]

BG_COLOR = "#d9d9d9"  # Light grey underlay for non-enriched / z <= 0 baseline

# Body site mapping dictionary to high-level anatomical display groups
BODY_SITE_MAP = {
    "feces": "Feces",
    "stool": "Feces",
    "saliva": "Oral cavity /\ndental",
    "oral cavity": "Oral cavity /\ndental",
    "tongue": "Oral cavity /\ndental",
    "dental plaque": "Oral cavity /\ndental",
    "supragingival plaque": "Oral cavity /\ndental",
    "subgingival": "Oral cavity /\ndental",
    "oropharynx": "Oral cavity /\ndental",
    "pharynx": "Oral cavity /\ndental",
    "gingival sulcus": "Oral cavity /\ndental",
    "vagina": "Vaginal /\nreproductive",
    "cervix": "Vaginal /\nreproductive",
    "cervicovaginal": "Vaginal /\nreproductive",
    "endometrium": "Vaginal /\nreproductive",
    "semen": "Vaginal /\nreproductive",
    "penis": "Vaginal /\nreproductive",
    "nasal cavity": "Upper airway /\nnasal",
    "nasopharynx": "Upper airway /\nnasal",
    "throat": "Upper airway /\nnasal",
    "nose": "Upper airway /\nnasal",
    "upper respiratory tract": "Upper airway /\nnasal",
    "colon": "Lower GI tract",
    "duodenum": "Lower GI tract",
    "ileum": "Lower GI tract",
    "rectum": "Lower GI tract",
    "intestine": "Lower GI tract",
    "cecum": "Lower GI tract",
    "small intestine": "Lower GI tract",
    "stomach": "Stomach /\nesophagus",
    "mucosa of body of stomach": "Stomach /\nesophagus",
    "gastric mucosa": "Stomach /\nesophagus",
    "esophagus": "Stomach /\nesophagus",
    "eye": "Ocular / eye",
    "conjunctiva": "Ocular / eye",
    "tear film": "Ocular / eye",
    "ocular surface": "Ocular / eye",
    "meconium": "Birth-associated\nsamples",
    "placenta": "Birth-associated\nsamples",
    "amniotic fluid": "Birth-associated\nsamples",
    "breast milk": "Birth-associated\nsamples",
    "cord blood": "Birth-associated\nsamples",
    "lung": "Lower airway /\nlung",
    "bronchus": "Lower airway /\nlung",
    "lower respiratory tract": "Lower airway /\nlung",
    "sputum": "Lower airway /\nlung",
    "skin": "Skin",
    "scalp": "Skin",
    "forehead": "Skin",
    "colonic mucosa": "Intestinal / GI\nmucosa",
    "duodenal mucosa": "Intestinal / GI\nmucosa",
    "intestinal mucosa": "Intestinal / GI\nmucosa",
}


def map_body_site(raw_site: str | None) -> str:
    if not raw_site or pd.isna(raw_site):
        return "Other"
    raw_lower = str(raw_site).strip().lower()
    # Sort keys by length descending to match more specific phrases first
    for k in sorted(BODY_SITE_MAP.keys(), key=len, reverse=True):
        if k in raw_lower:
            return BODY_SITE_MAP[k]
    return raw_site.strip()


def load_trait_labels(trait_meta_path: Path, mapping_path: Path) -> dict[str, str]:
    labels = {
        "METPO:1000608": "facultatively aerobic",
        "METPO:2000103|GO:0030245": "cellulose catabolic process",
        "METPO:1000668": "coccus shaped",
        "METPO:1000683": "sphere shaped",
        "METPO:2000103|pathways:arsenate_reduction": "arsenate reduction",
        "METPO:1000623": "moderately halophilic",
    }
    if trait_meta_path.exists():
        try:
            tm = pd.read_csv(trait_meta_path, sep="\t")
            for _, r in tm.iterrows():
                cid = str(r["column_id"])
                lbl = str(r.get("column_label", cid))
                lbl = lbl.split(" [METPO")[0].split(" METPO:")[0].replace("capable of: ", "")
                labels[cid] = lbl
        except Exception:
            pass
    if mapping_path.exists():
        try:
            mp = pd.read_csv(mapping_path, sep="\t", header=None)
            for _, r in mp.iterrows():
                cid = str(r[0])
                lbl = str(r[1]).split(" [METPO")[0].split(" METPO:")[0]
                labels[cid] = lbl
        except Exception:
            pass
    return labels


def compute_z_score_or_signal(df: pd.DataFrame) -> pd.Series:
    """Compute positive signal (z-score approximation from p-value or log2 fold)."""
    if "p_value" in df.columns:
        p_clamped = df["p_value"].clip(lower=1e-15, upper=1.0)
        # Standard normal inverse survival approximation: z ~ sqrt(2 * ln(1/p))
        z_approx = np.sqrt(np.maximum(0.0, 2 * np.log(1.0 / p_clamped)))
        return z_approx.clip(lower=0.0, upper=3.5)
    elif "fold" in df.columns:
        return np.log2(np.maximum(df["fold"], 1.0)).clip(lower=0.0, upper=3.5)
    return pd.Series([1.0] * len(df), index=df.index)


def get_cubic_bezier_ribbon(
    x0: float, y0_bot: float, y0_top: float, x1: float, y1_bot: float, y1_top: float
) -> MplPath:
    """Create a cubic Bezier ribbon with explicit top and bottom bounds at each end."""
    dx = (x1 - x0) * 0.5

    verts = [
        (x0, y0_bot),
        (x0 + dx, y0_bot),
        (x1 - dx, y1_bot),
        (x1, y1_bot),
        (x1, y1_top),
        (x1 - dx, y1_top),
        (x0 + dx, y0_top),
        (x0, y0_top),
        (x0, y0_bot),
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


def get_cubic_bezier_path(x0: float, y0: float, x1: float, y1: float, width: float) -> MplPath:
    """Create a cubic Bezier ribbon of specified width centered at y0 and y1."""
    hw = width / 2.0
    return get_cubic_bezier_ribbon(x0, y0 - hw, y0 + hw, x1, y1 - hw, y1 + hw)


def plot_grouped_six_trait_sankey(
    enrich_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    selected_traits: list[str],
    trait_labels: dict[str, str],
    direction: str,
    output_png: Path,
    output_svg: Path | None = None,
    min_profiles: int = 2,
    min_studies: int = 1,
) -> None:
    # 1. Prepare data
    df = enrich_df[
        (enrich_df["trait"].isin(selected_traits)) & (enrich_df["inc_dec"] == direction)
    ].copy()

    if df.empty:
        print(f"No enrichment data for direction={direction}")
        return

    # Map condition and body site
    if "condition" not in df.columns or "body_site" not in df.columns:
        lookup = metadata_df.drop_duplicates("experiment_id").set_index("experiment_id")
        df["condition"] = df["experiment_id"].map(lookup["condition"]).fillna("Unknown")
        df["raw_body_site"] = df["experiment_id"].map(lookup["body_site"]).fillna("Unknown")
    else:
        df["raw_body_site"] = df["body_site"]

    df["display_body_site"] = df["raw_body_site"].apply(map_body_site)
    df = df[(df["condition"] != "Unknown") & (df["display_body_site"] != "Other")].copy()

    df["z_score"] = compute_z_score_or_signal(df)

    # 2. Aggregate edges by condition, display_body_site, trait
    edge_agg = (
        df.groupby(["condition", "display_body_site", "trait"])
        .agg(
            n_profiles=("profile_id", "nunique"),
            n_studies=("study_id", "nunique"),
            mean_z=("z_score", "mean"),
        )
        .reset_index()
    )

    # Filter by minimum support
    edge_agg = edge_agg[
        (edge_agg["n_profiles"] >= min_profiles) & (edge_agg["n_studies"] >= min_studies)
    ].copy()

    if edge_agg.empty:
        msg = f"No edges passed thresholds (min_p={min_profiles}, min_s={min_studies})"
        print(msg)
        return

    # 3. Calculate Condition & Body Site totals and ordering
    cond_stats = (
        edge_agg.groupby("condition")
        .agg(
            total_n=("n_profiles", "sum"),
            n_links=("trait", "count"),
            n_traits=("trait", "nunique"),
        )
        .sort_values("total_n", ascending=False)
    )

    site_stats = (
        edge_agg.groupby("display_body_site")
        .agg(
            total_n=("n_profiles", "sum"),
            n_links=("trait", "count"),
            n_traits=("trait", "nunique"),
        )
        .sort_values("total_n", ascending=False)
    )

    conditions = cond_stats.index.tolist()
    body_sites = site_stats.index.tolist()

    # Layout coordinates (Canvas size 2400 x 2400)
    fig_w, fig_h = 24.0, 24.0
    _fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=200)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")

    # Header title and notes
    dir_word = "decreased" if direction == "dec" else "increased"
    title = f"Six shared traits, individual conditions: {dir_word} profile Sankey"
    ax.text(
        3.0,
        96.5,
        title,
        fontsize=22,
        fontweight="bold",
        fontfamily="sans-serif",
        color="#111111",
        va="top",
    )

    subtitle_lines = [
        f"Colored curves are {dir_word} profiles with trait z-score > 0; "
        "parallel trait colors distinguish the six traits.",
        "Width = positive-z profile count; shade = mean positive z-score "
        "on a shared numeric scale for all six traits.",
        "Within each parallel condition/body-site stack, curves are ordered "
        "thick-to-thin. Light grey underlays show z <= 0 background.",
        "Node key: links = plotted trait-specific links, traits = distinct traits, "
        "n = positive-z profile memberships summed across displayed trait links.",
        "Body-site trunks sum same-color links. Conditions use raw individual labels; "
        "body sites keep grouped labels.",
    ]
    subtitle = "\n".join(subtitle_lines)
    ax.text(
        3.0,
        94.0,
        subtitle,
        fontsize=10.5,
        color="#333333",
        fontfamily="sans-serif",
        va="top",
        linespacing=1.4,
    )

    # Column Titles
    ax.text(
        22.0,
        86.0,
        "Individual conditions",
        fontsize=14,
        fontweight="bold",
        color="#222222",
        ha="center",
    )
    ax.text(
        66.0,
        86.0,
        "Body-site display groups",
        fontsize=14,
        fontweight="bold",
        color="#222222",
        ha="center",
    )

    # Y-placement calculation for nodes
    top_y = 82.0
    bot_y = 10.0
    usable_h = top_y - bot_y

    n_cond = len(conditions)
    cond_spacing = usable_h / max(1, n_cond)
    cond_y_map = {c: top_y - (i + 0.5) * cond_spacing for i, c in enumerate(conditions)}

    n_sites = len(body_sites)
    site_spacing = usable_h / max(1, n_sites)
    site_y_map = {s: top_y - (i + 0.5) * site_spacing for i, s in enumerate(body_sites)}

    x_cond = 32.0
    x_site = 60.0

    # Function to scale ribbon width proportionally to profile count
    def calc_ribbon_width(n_profiles: int) -> float:
        return float(np.clip(0.04 + 0.035 * np.sqrt(n_profiles), 0.05, 0.55))

    # Compute stacked heights for condition nodes and body-site nodes
    # Sort edges exiting c: order by target site position (top to bottom), then by trait
    # Sort edges entering s: order by source condition position (top to bottom), then by trait
    cond_index_map = {c: i for i, c in enumerate(conditions)}
    site_index_map = {s: i for i, s in enumerate(body_sites)}

    edge_agg["cond_idx"] = edge_agg["condition"].map(cond_index_map)
    edge_agg["site_idx"] = edge_agg["display_body_site"].map(site_index_map)
    edge_agg["ribbon_w"] = edge_agg["n_profiles"].apply(calc_ribbon_width)

    # Compute total stacked height per condition and per body-site
    cond_total_w = edge_agg.groupby("condition")["ribbon_w"].sum().to_dict()
    site_total_w = edge_agg.groupby("display_body_site")["ribbon_w"].sum().to_dict()

    # Calculate exact vertical slices on condition side
    # Sort edges for condition: top-to-bottom target site, then trait index
    cond_edge_sort = edge_agg.sort_values(
        by=["cond_idx", "site_idx", "n_profiles"], ascending=[True, True, False]
    ).copy()

    # Track current offset per condition (starting from top edge down to bottom edge)
    cond_curr_y = {c: cond_y_map[c] + cond_total_w.get(c, 0.0) / 2.0 for c in conditions}
    edge_cond_bounds = {}
    for edge_id, row in cond_edge_sort.iterrows():
        c = row["condition"]
        w = row["ribbon_w"]
        top_y = cond_curr_y[c]
        bot_y = top_y - w
        cond_curr_y[c] = bot_y
        edge_cond_bounds[edge_id] = (bot_y, top_y)

    # Calculate exact vertical slices on body site side
    site_edge_sort = edge_agg.sort_values(
        by=["site_idx", "cond_idx", "n_profiles"], ascending=[True, True, False]
    ).copy()

    site_curr_y = {s: site_y_map[s] + site_total_w.get(s, 0.0) / 2.0 for s in body_sites}
    edge_site_bounds = {}
    for edge_id, row in site_edge_sort.iterrows():
        s = row["display_body_site"]
        w = row["ribbon_w"]
        top_y = site_curr_y[s]
        bot_y = top_y - w
        site_curr_y[s] = bot_y
        edge_site_bounds[edge_id] = (bot_y, top_y)

    # Colormaps for each trait
    trait_color_map = {
        t: TRAIT_COLORS[i % len(TRAIT_COLORS)] for i, t in enumerate(selected_traits)
    }
    trait_cmaps = {}
    for t, base_hex in trait_color_map.items():
        base_rgb = mcolors.to_rgb(base_hex)
        cmap = mcolors.LinearSegmentedColormap.from_list(
            f"cmap_{t}",
            [(0.95, 0.95, 0.95), base_rgb],
        )
        trait_cmaps[t] = cmap

    # Draw Nodes: node height matches total connection thickness (pill/capsule shape)
    node_w = 1.4
    for c in conditions:
        cy = cond_y_map[c]
        tot_h = max(cond_total_w.get(c, 0.6), 0.6)
        c_patch = patches.FancyBboxPatch(
            (x_cond - node_w / 2.0, cy - tot_h / 2.0),
            node_w,
            tot_h,
            boxstyle=patches.BoxStyle.Round(pad=0.1, rounding_size=0.3),
            facecolor="#f4f4f4",
            edgecolor="#555555",
            linewidth=1.2,
            zorder=10,
        )
        ax.add_patch(c_patch)

        # Label on left
        c_tot = cond_stats.loc[c, "total_n"]
        c_links = cond_stats.loc[c, "n_links"]
        c_traits = cond_stats.loc[c, "n_traits"]
        label_text = (
            f"{c}\n"
            rf"$\mathrm{{links}}={c_links}, \mathrm{{traits}}={c_traits}, n={c_tot}$"
        )
        ax.text(
            x_cond - 1.8,
            cy,
            label_text,
            ha="right",
            va="center",
            fontsize=9.5,
            color="#111111",
            fontfamily="sans-serif",
        )

    for s in body_sites:
        sy = site_y_map[s]
        tot_h = max(site_total_w.get(s, 0.6), 0.6)
        s_patch = patches.FancyBboxPatch(
            (x_site - node_w / 2.0, sy - tot_h / 2.0),
            node_w,
            tot_h,
            boxstyle=patches.BoxStyle.Round(pad=0.1, rounding_size=0.3),
            facecolor="#f4f4f4",
            edgecolor="#555555",
            linewidth=1.2,
            zorder=10,
        )
        ax.add_patch(s_patch)

        # Label on right
        s_tot = site_stats.loc[s, "total_n"]
        s_links = site_stats.loc[s, "n_links"]
        s_traits = site_stats.loc[s, "n_traits"]
        site_clean_label = s.replace("\n", " ")
        label_text = (
            f"{site_clean_label}\n"
            rf"$\mathrm{{links}}={s_links}, \mathrm{{traits}}={s_traits}, n={s_tot}$"
        )
        ax.text(
            x_site + 1.8,
            sy,
            label_text,
            ha="left",
            va="center",
            fontsize=9.5,
            color="#111111",
            fontfamily="sans-serif",
        )

    # 4. Draw Stacked Non-Overlapping Bezier Ribbons
    for edge_id, row in edge_agg.iterrows():
        t = row["trait"]
        mz = row["mean_z"]

        c_bot, c_top = edge_cond_bounds[edge_id]
        s_bot, s_top = edge_site_bounds[edge_id]

        # Shading by z-score (0 to 3.0 scale)
        norm_z = np.clip(mz / 2.9, 0.15, 1.0)
        color_rgba = trait_cmaps[t](norm_z)

        # Background underlay link (slight margin)
        bg_ribbon = get_cubic_bezier_ribbon(
            x_cond + node_w / 2.0,
            c_bot - 0.01,
            c_top + 0.01,
            x_site - node_w / 2.0,
            s_bot - 0.01,
            s_top + 0.01,
        )
        patch_bg = patches.PathPatch(
            bg_ribbon, facecolor=BG_COLOR, edgecolor="none", alpha=0.35, zorder=2
        )
        ax.add_patch(patch_bg)

        # Foreground colored ribbon
        fg_ribbon = get_cubic_bezier_ribbon(
            x_cond + node_w / 2.0,
            c_bot,
            c_top,
            x_site - node_w / 2.0,
            s_bot,
            s_top,
        )
        patch_fg = patches.PathPatch(
            fg_ribbon, facecolor=color_rgba, edgecolor="none", alpha=0.85, zorder=4
        )
        ax.add_patch(patch_fg)

    # 6. Draw Trait Color Scales on the Right
    x_cbar = 83.0
    cbar_w = 12.0
    cbar_h = 0.8
    y_cbar_top = 85.0
    cbar_spacing = 8.5

    ax.text(
        x_cbar,
        y_cbar_top + 6.0,
        "Trait color scales",
        fontsize=14,
        fontweight="bold",
        color="#222222",
    )
    ax.text(
        x_cbar,
        y_cbar_top + 3.2,
        "Shared numeric range for all traits.\nShade = mean positive z-score.",
        fontsize=9.0,
        color="#555555",
        linespacing=1.3,
    )

    for i, t in enumerate(selected_traits):
        yc = y_cbar_top - i * cbar_spacing
        lbl = trait_labels.get(t, t)

        ax.text(
            x_cbar,
            yc + 1.6,
            lbl,
            fontsize=10.0,
            fontweight="semibold",
            color="#222222",
        )

        # Gradient bar
        grad = np.linspace(0.1, 1.0, 256).reshape(1, -1)
        ax.imshow(
            grad,
            extent=(x_cbar, x_cbar + cbar_w, yc, yc + cbar_h),
            cmap=trait_cmaps[t],
            aspect="auto",
            zorder=10,
        )
        rect_border = patches.Rectangle(
            (x_cbar, yc),
            cbar_w,
            cbar_h,
            facecolor="none",
            edgecolor="#555555",
            linewidth=0.8,
            zorder=11,
        )
        ax.add_patch(rect_border)

        ax.text(x_cbar, yc - 1.2, "0", fontsize=8.5, color="#444444", ha="left")
        ax.text(x_cbar + cbar_w, yc - 1.2, "2.9", fontsize=8.5, color="#444444", ha="right")

    # 7. Draw Bottom Legend
    y_leg = 3.5
    ax.text(
        50.0,
        y_leg + 2.0,
        "Edge width: source n; body-site trunks sum n per trait",
        fontsize=11.0,
        fontweight="bold",
        ha="center",
        color="#222222",
    )

    leg_items = [
        ("n=3 profiles", 0.10, "#333333"),
        ("n=10 profiles", 0.22, "#333333"),
        ("n=67 profiles", 0.45, "#333333"),
        ("body-site trunk: summed n", 0.30, "#333333"),
        ("bg profiles: z <= 0", 0.15, BG_COLOR),
    ]

    x_leg_start = 28.0
    for label, width, color in leg_items:
        ax.plot(
            [x_leg_start, x_leg_start + 3.0],
            [y_leg, y_leg],
            color=color,
            linewidth=width * 15.0,
            solid_capstyle="round",
        )
        ax.text(
            x_leg_start + 3.5,
            y_leg,
            label,
            va="center",
            fontsize=9.0,
            color="#333333",
        )
        x_leg_start += 9.5

    # Footer notes
    footer_lines = [
        "Every line is trait-specific. Node counts are profile-trait memberships across "
        "displayed links, not deduplicated profile totals.",
        "Display filtering uses thresholds after aggregating by raw condition, "
        "grouped body-site label, and trait: at least 2 positive-z focus profiles.",
        "Body-site trunks sum all displayed condition links by trait color. "
        "Hybrid display mode keeps individual condition labels on the left "
        "and grouped body-site labels on the right.",
    ]
    footer_text = "\n".join(footer_lines)
    ax.text(
        3.0,
        0.8,
        footer_text,
        fontsize=8.5,
        color="#555555",
        fontfamily="sans-serif",
        va="bottom",
        linespacing=1.3,
    )

    output_png.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(output_png, dpi=250, bbox_inches="tight")
    print(f"Saved PNG to {output_png}")

    if output_svg:
        plt.savefig(output_svg, bbox_inches="tight")
        print(f"Saved SVG to {output_svg}")

    plt.close()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate grouped 6-trait Sankey diagrams for increased and decreased profiles."
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
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory to save output Sankey diagrams",
    )
    parser.add_argument(
        "--direction",
        choices=["inc", "dec", "both"],
        default="both",
        help="Direction to generate ('inc', 'dec', or 'both')",
    )
    parser.add_argument(
        "--min-profiles",
        type=int,
        default=2,
        help="Minimum profile count per edge",
    )
    parser.add_argument(
        "--min-studies",
        type=int,
        default=1,
        help="Minimum study count per edge",
    )

    args = parser.parse_args()

    # Load data
    enrich_path = args.enrichment
    if not enrich_path.exists():
        if DEFAULT_RANKED_PAIRS_PATH.exists():
            enrich_path = DEFAULT_RANKED_PAIRS_PATH
        else:
            raise FileNotFoundError(f"Cannot find enrichment TSV at {args.enrichment}")

    print(f"Loading enrichment data from {enrich_path}...")
    enrich_df = pd.read_csv(enrich_path, sep="\t", comment="#")

    meta_df = pd.DataFrame()
    if args.experiments.exists():
        print(f"Loading metadata from {args.experiments}...")
        meta_df = pd.read_csv(args.experiments, sep="\t")
    elif DEFAULT_RAW_EXPERIMENTS_PATH.exists():
        print(f"Loading metadata from {DEFAULT_RAW_EXPERIMENTS_PATH}...")
        raw_df = pd.read_csv(DEFAULT_RAW_EXPERIMENTS_PATH, comment="#", dtype=str)
        raw_df.columns = [c.strip().replace('"', "") for c in raw_df.columns]
        raw_df["study_key"] = "bsdb:" + raw_df["Study"].astype(str).str.strip()
        raw_df["exp_num"] = raw_df["Experiment"].str.extract(r"(\d+)", expand=False).fillna("1")
        raw_df["experiment_id"] = raw_df["study_key"] + "-" + raw_df["exp_num"]
        meta_df = raw_df.rename(columns={"Body site": "body_site", "Condition": "condition"})

    trait_labels = load_trait_labels(DEFAULT_TRAIT_META_PATH, DEFAULT_PHENOTYPE_MAPPING_PATH)

    directions = ["dec", "inc"] if args.direction == "both" else [args.direction]
    args.output_dir.mkdir(parents=True, exist_ok=True)

    for d in directions:
        dir_name = "decreased" if d == "dec" else "increased"
        out_png = args.output_dir / f"six_trait_sankey_{dir_name}.png"
        out_svg = args.output_dir / f"six_trait_sankey_{dir_name}.svg"

        print(f"Generating 6-trait Sankey diagram for {dir_name} profiles...")
        plot_grouped_six_trait_sankey(
            enrich_df=enrich_df,
            metadata_df=meta_df,
            selected_traits=DEFAULT_TRAITS,
            trait_labels=trait_labels,
            direction=d,
            output_png=out_png,
            output_svg=out_svg,
            min_profiles=args.min_profiles,
            min_studies=args.min_studies,
        )


if __name__ == "__main__":
    main()
