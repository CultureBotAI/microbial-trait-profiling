from pathlib import Path

import pandas as pd
import plotly.graph_objects as go

# Paths
ENRICHMENT_PATH = Path("data/derived/trait_enrichment_trait_study_pairs_ranked.tsv")
EXPERIMENTS_PATH = Path("data/raw/full_dump.csv")
OUTPUT_DIR = Path("data/derived/Sankey_Diagrams")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Configuration
GUILDS = {
    "Anaerobe": [
        "METPO:1000603",  # anaerobic
        "METPO:1000607",  # obligately anaerobic
    ],
    "Facultative-Aerobe": [
        "METPO:1000608",  # facultatively aerobic
        "METPO:1000602",  # aerobic
        "METPO:1000604",  # microaerophilic
    ],
}

MIN_PROFILES = 3
MIN_STUDIES = 2

# Colors
COLOR_INC = "rgba(214, 39, 40, 0.65)"  # Red for increased
COLOR_DEC = "rgba(31, 119, 180, 0.65)"  # Blue for decreased
COLOR_NONE = "rgba(180, 180, 180, 0.30)"  # Gray for no enrichment


def generate_single_guild_sankey_with_background(
    enrich_df, total_counts_df, guild_name, direction, title, max_conditions: int = 20
):
    guild_traits = GUILDS[guild_name]
    g_df = enrich_df[
        (enrich_df["trait"].isin(guild_traits)) & (enrich_df["inc_dec"] == direction)
    ].copy()

    if g_df.empty:
        print(f"No data for {guild_name} in {direction} direction.")
        return

    # Profile-level aggregation for enriched
    profile_level = (
        g_df.groupby(["condition", "body_site", "study_id", "profile_id"])
        .agg(max_fold=("fold", "max"))
        .reset_index()
    )

    enriched_edges = (
        profile_level.groupby(["condition", "body_site"])
        .agg(
            n_profiles=("profile_id", "nunique"),
            s_studies=("study_id", "nunique"),
            mean_fold=("max_fold", "mean"),
        )
        .reset_index()
    )

    # Filter by support thresholds for enriched edges
    qualifying_pairs = enriched_edges[
        (enriched_edges["n_profiles"] >= MIN_PROFILES)
        & (enriched_edges["s_studies"] >= MIN_STUDIES)
    ][["condition", "body_site"]].drop_duplicates()

    if qualifying_pairs.empty:
        print(f"No edges meet thresholds for {guild_name} ({direction}).")
        return

    # Merge with total counts to compute non-enriched parallel edges
    merged = pd.merge(qualifying_pairs, total_counts_df, on=["condition", "body_site"], how="left")
    merged = pd.merge(
        merged,
        enriched_edges[["condition", "body_site", "n_profiles", "s_studies", "mean_fold"]],
        on=["condition", "body_site"],
        how="left",
    )
    merged["n_profiles"] = merged["n_profiles"].fillna(0).astype(int)
    merged["n_no_enrich"] = (merged["total_profiles"] - merged["n_profiles"]).clip(lower=0)

    # Filter to top conditions by total profiles to ensure publication-quality spacing
    cond_order = (
        merged.groupby("condition")["total_profiles"]
        .sum()
        .sort_values(ascending=False)
        .head(max_conditions)
        .index.tolist()
    )
    merged = merged[merged["condition"].isin(cond_order)].copy()

    site_order = (
        merged.groupby("body_site")["total_profiles"]
        .sum()
        .sort_values(ascending=False)
        .index.tolist()
    )

    all_nodes = cond_order + site_order
    node_map = {node: i for i, node in enumerate(all_nodes)}

    sources = []
    targets = []
    values = []
    edge_colors = []

    for _, row in merged.iterrows():
        c_idx = node_map[row["condition"]]
        s_idx = node_map[row["body_site"]]

        # Enriched link
        if row["n_profiles"] > 0:
            sources.append(c_idx)
            targets.append(s_idx)
            values.append(int(row["n_profiles"]))
            edge_colors.append(COLOR_INC if direction == "inc" else COLOR_DEC)

        # Non-enriched parallel link
        if row["n_no_enrich"] > 0:
            sources.append(c_idx)
            targets.append(s_idx)
            values.append(int(row["n_no_enrich"]))
            edge_colors.append(COLOR_NONE)

    # Node Labels
    node_labels = []
    for node in all_nodes:
        if node in cond_order:
            n = merged[merged["condition"] == node]["total_profiles"].sum()
            s = merged[merged["condition"] == node]["total_studies"].max()
        else:
            n = merged[merged["body_site"] == node]["total_profiles"].sum()
            s = merged[merged["body_site"] == node]["total_studies"].max()
        node_labels.append(f"{node}<br>total n={n}, s={s}")

    fig = go.Figure(
        data=[
            go.Sankey(
                node={
                    "pad": 18,
                    "thickness": 20,
                    "line": {"color": "black", "width": 0.5},
                    "label": node_labels,
                    "color": "lightgray",
                    "x": [0.01] * len(cond_order) + [0.99] * len(site_order),
                },
                link={"source": sources, "target": targets, "value": values, "color": edge_colors},
                arrangement="snap",
            )
        ]
    )

    color_desc = "Red" if direction == "inc" else "Blue"
    subtitle = (
        f"<span style='font-size:10px'>Colored links: {direction}reased enrichment ({color_desc}). "
        "Gray links: no enrichment. Node n sums profile memberships.<br>"
        "Interpretation: Results measure annotation representation in curated taxon signatures "
        "(urn model, M = 5,184 background taxa), "
        "not direct patient-level organism abundance.</span>"
    )
    fig.update_layout(
        title_text=f"<b>{title}</b><br>{subtitle}",
        font_size=12,
        width=1200,
        height=1300,
        margin={"t": 140, "b": 40, "l": 40, "r": 40},
    )

    filename = OUTPUT_DIR / f"{guild_name.lower()}_sankey_{direction}_with_unenriched.png"
    fig.write_image(str(filename), scale=2)
    print(f"Saved {filename}")

    legacy_filename = OUTPUT_DIR / f"{guild_name.lower()}_sankey_{direction}.png"
    fig.write_image(str(legacy_filename), scale=2)
    print(f"Saved {legacy_filename}")


def generate_combined_sankey(enrich_df, total_counts_df, max_conditions: int = 20):
    """Generates a 3-column Sankey:

    Left: Aerobe (Conditions split by Inc, Dec, None)
    Center: Anatomy / Body Sites
    Right: Anaerobe (Conditions split by Inc, Dec, None)
    """
    aerobe_traits = GUILDS["Facultative-Aerobe"]
    anaerobe_traits = GUILDS["Anaerobe"]

    # Filter enrichments
    aero_inc_profiles = (
        enrich_df[(enrich_df["trait"].isin(aerobe_traits)) & (enrich_df["inc_dec"] == "inc")][
            ["condition", "body_site", "profile_id", "study_id"]
        ]
        .drop_duplicates()
        .groupby(["condition", "body_site"])
        .agg(n_inc=("profile_id", "nunique"), s_inc=("study_id", "nunique"))
        .reset_index()
    )

    aero_dec_profiles = (
        enrich_df[(enrich_df["trait"].isin(aerobe_traits)) & (enrich_df["inc_dec"] == "dec")][
            ["condition", "body_site", "profile_id", "study_id"]
        ]
        .drop_duplicates()
        .groupby(["condition", "body_site"])
        .agg(n_dec=("profile_id", "nunique"), s_dec=("study_id", "nunique"))
        .reset_index()
    )

    anaero_inc_profiles = (
        enrich_df[(enrich_df["trait"].isin(anaerobe_traits)) & (enrich_df["inc_dec"] == "inc")][
            ["condition", "body_site", "profile_id", "study_id"]
        ]
        .drop_duplicates()
        .groupby(["condition", "body_site"])
        .agg(n_inc=("profile_id", "nunique"), s_inc=("study_id", "nunique"))
        .reset_index()
    )

    anaero_dec_profiles = (
        enrich_df[(enrich_df["trait"].isin(anaerobe_traits)) & (enrich_df["inc_dec"] == "dec")][
            ["condition", "body_site", "profile_id", "study_id"]
        ]
        .drop_duplicates()
        .groupby(["condition", "body_site"])
        .agg(n_dec=("profile_id", "nunique"), s_dec=("study_id", "nunique"))
        .reset_index()
    )

    # Filter significant conditions based on support thresholds
    significant_pairs = set()
    for df in [aero_inc_profiles, aero_dec_profiles, anaero_inc_profiles, anaero_dec_profiles]:
        hits = df[(df.get("n_inc", 0) >= MIN_PROFILES) | (df.get("n_dec", 0) >= MIN_PROFILES)]
        for _, r in hits.iterrows():
            significant_pairs.add((r["condition"], r["body_site"]))

    if not significant_pairs:
        print("No significant pairs for combined sankey.")
        return

    sig_df = pd.DataFrame(list(significant_pairs), columns=["condition", "body_site"])
    base_df = pd.merge(sig_df, total_counts_df, on=["condition", "body_site"], how="left")

    # Filter to top conditions by total profiles to ensure legible publication layout
    left_conditions = (
        base_df.groupby("condition")["total_profiles"]
        .sum()
        .sort_values(ascending=False)
        .head(max_conditions)
        .index.tolist()
    )
    base_df = base_df[base_df["condition"].isin(left_conditions)].copy()
    center_sites = (
        base_df.groupby("body_site")["total_profiles"]
        .sum()
        .sort_values(ascending=False)
        .index.tolist()
    )
    right_conditions = list(left_conditions)

    # Node indexing
    # Left nodes: 0 .. len(left_conditions)-1
    # Center nodes: len(left_conditions) .. len(left_conditions)+len(center_sites)-1
    # Right nodes: len(left_conditions)+len(center_sites) ..
    left_node_map = {c: i for i, c in enumerate(left_conditions)}
    offset_center = len(left_conditions)
    center_node_map = {s: offset_center + i for i, s in enumerate(center_sites)}
    offset_right = offset_center + len(center_sites)
    right_node_map = {c: offset_right + i for i, c in enumerate(right_conditions)}

    node_labels = []
    for c in left_conditions:
        n = base_df[base_df["condition"] == c]["total_profiles"].sum()
        node_labels.append(f"Aerobes: {c}<br>(n={n})")
    for s in center_sites:
        n = base_df[base_df["body_site"] == s]["total_profiles"].sum()
        node_labels.append(f"<b>{s}</b><br>(n={n})")
    for c in right_conditions:
        n = base_df[base_df["condition"] == c]["total_profiles"].sum()
        node_labels.append(f"Anaerobes: {c}<br>(n={n})")

    sources = []
    targets = []
    values = []
    edge_colors = []

    # Map Aerobe edges (Left -> Center)
    for _, row in base_df.iterrows():
        c = row["condition"]
        s = row["body_site"]
        c_idx = left_node_map[c]
        s_idx = center_node_map[s]
        tot = int(row["total_profiles"])

        # Find aero inc / dec
        a_inc = aero_inc_profiles[
            (aero_inc_profiles["condition"] == c) & (aero_inc_profiles["body_site"] == s)
        ]
        a_dec = aero_dec_profiles[
            (aero_dec_profiles["condition"] == c) & (aero_dec_profiles["body_site"] == s)
        ]

        n_inc = int(a_inc["n_inc"].iloc[0]) if not a_inc.empty else 0
        n_dec = int(a_dec["n_dec"].iloc[0]) if not a_dec.empty else 0
        n_none = max(0, tot - n_inc - n_dec)

        if n_inc > 0:
            sources.append(c_idx)
            targets.append(s_idx)
            values.append(n_inc)
            edge_colors.append(COLOR_INC)

        if n_dec > 0:
            sources.append(c_idx)
            targets.append(s_idx)
            values.append(n_dec)
            edge_colors.append(COLOR_DEC)

        if n_none > 0:
            sources.append(c_idx)
            targets.append(s_idx)
            values.append(n_none)
            edge_colors.append(COLOR_NONE)

    # Map Anaerobe edges (Center -> Right)
    for _, row in base_df.iterrows():
        c = row["condition"]
        s = row["body_site"]
        s_idx = center_node_map[s]
        c_idx = right_node_map[c]
        tot = int(row["total_profiles"])

        # Find anaero inc / dec
        an_inc = anaero_inc_profiles[
            (anaero_inc_profiles["condition"] == c) & (anaero_inc_profiles["body_site"] == s)
        ]
        an_dec = anaero_dec_profiles[
            (anaero_dec_profiles["condition"] == c) & (anaero_dec_profiles["body_site"] == s)
        ]

        n_inc = int(an_inc["n_inc"].iloc[0]) if not an_inc.empty else 0
        n_dec = int(an_dec["n_dec"].iloc[0]) if not an_dec.empty else 0
        n_none = max(0, tot - n_inc - n_dec)

        if n_inc > 0:
            sources.append(s_idx)
            targets.append(c_idx)
            values.append(n_inc)
            edge_colors.append(COLOR_INC)

        if n_dec > 0:
            sources.append(s_idx)
            targets.append(c_idx)
            values.append(n_dec)
            edge_colors.append(COLOR_DEC)

        if n_none > 0:
            sources.append(s_idx)
            targets.append(c_idx)
            values.append(n_none)
            edge_colors.append(COLOR_NONE)

    fig = go.Figure(
        data=[
            go.Sankey(
                node={
                    "pad": 18,
                    "thickness": 20,
                    "line": {"color": "black", "width": 0.5},
                    "label": node_labels,
                    "color": "lightgray",
                    "x": (
                        [0.01] * len(left_conditions)
                        + [0.50] * len(center_sites)
                        + [0.99] * len(right_conditions)
                    ),
                },
                link={"source": sources, "target": targets, "value": values, "color": edge_colors},
                arrangement="snap",
            )
        ]
    )

    combined_subtitle = (
        "<span style='font-size:11px'>Left: Aerobes & Facultative Aerobes | "
        "Center: Anatomy / Body Sites | Right: Anaerobes<br>"
        "Red = Increased Enrichment | Blue = Decreased Enrichment | "
        "Gray = No Significant Enrichment<br>"
        "Interpretation: Results measure annotation representation in curated taxon signatures "
        "(urn model, M = 5,184 background taxa), "
        "not direct patient-level organism abundance.</span>"
    )
    title = "<b>Combined Guild Enrichment Across Anatomical Sites</b>"
    fig.update_layout(
        title_text=f"{title}<br>{combined_subtitle}",
        font_size=12,
        width=1600,
        height=1300,
        margin={"t": 140, "b": 40, "l": 40, "r": 40},
    )

    filename = OUTPUT_DIR / "combined_aerobe_anaerobe_sankey.png"
    fig.write_image(str(filename), scale=2)
    print(f"Saved {filename}")


def main():
    print("Loading enrichment data...")
    enrich_df = pd.read_csv(ENRICHMENT_PATH, sep="\t", comment="#")

    if "profile_id" not in enrich_df.columns:
        enrich_df["profile_id"] = (
            enrich_df["study_id"].astype(str)
            + "|"
            + enrich_df["experiment_id"].astype(str)
            + "|"
            + enrich_df["inc_dec"].astype(str)
        )

    print("Loading metadata...")
    derived_exp_path = Path("data/derived/experiments.tsv")
    if derived_exp_path.exists():
        print(f"Loading metadata from {derived_exp_path}...")
        meta_df = pd.read_csv(derived_exp_path, sep="\t", dtype=str)
        lookup = (
            meta_df.drop_duplicates(subset="experiment_id")
            .set_index("experiment_id")[["body_site", "condition"]]
            .to_dict("index")
        )
        enrich_df["body_site"] = enrich_df["experiment_id"].map(
            lambda x: lookup.get(x, {}).get("body_site", "Unknown")
        )
        enrich_df["condition"] = enrich_df["experiment_id"].map(
            lambda x: lookup.get(x, {}).get("condition", "Unknown")
        )
    elif EXPERIMENTS_PATH.exists():
        print(f"Loading metadata from {EXPERIMENTS_PATH}...")
        meta_df = pd.read_csv(EXPERIMENTS_PATH, comment="#", dtype=str)
        meta_df.columns = [c.strip().replace('"', "") for c in meta_df.columns]

        meta_df["study_key"] = "bsdb:" + meta_df["Study"].astype(str).str.strip()
        meta_df["exp_num"] = meta_df["Experiment"].str.extract(r"(\d+)", expand=False).fillna("1")
        meta_df["exp_key"] = meta_df["study_key"] + "-" + meta_df["exp_num"]

        lookup = (
            meta_df.drop_duplicates(subset="exp_key")
            .set_index("exp_key")[["Body site", "Condition"]]
            .to_dict("index")
        )

        print("Mapping metadata...")
        enrich_df["body_site"] = enrich_df["experiment_id"].map(
            lambda x: lookup.get(x, {}).get("Body site", "Unknown")
        )
        enrich_df["condition"] = enrich_df["experiment_id"].map(
            lambda x: lookup.get(x, {}).get("Condition", "Unknown")
        )

    # Clean and filter
    enrich_df = enrich_df[
        (enrich_df["body_site"] != "Unknown") & (enrich_df["condition"] != "Unknown")
    ].copy()

    # Calculate total baseline profile & study counts per (condition, body_site)
    total_counts_df = (
        enrich_df.groupby(["condition", "body_site"])
        .agg(
            total_profiles=("profile_id", "nunique"),
            total_studies=("study_id", "nunique"),
        )
        .reset_index()
    )

    print("Generating single guild Sankey diagrams with background non-enriched links...")
    # 1. Anaerobe Guild in Decreased Profiles
    generate_single_guild_sankey_with_background(
        enrich_df,
        total_counts_df,
        "Anaerobe",
        "dec",
        "Condition-to-Body-Site Map: Anaerobe Guild in Decreased Profiles",
    )

    # 2. Facultative-Aerobe Guild in Increased Profiles
    generate_single_guild_sankey_with_background(
        enrich_df,
        total_counts_df,
        "Facultative-Aerobe",
        "inc",
        "Condition-to-Body-Site Map: Facultative-Aerobe Guild in Increased Profiles",
    )

    print("Generating combined 3-column Sankey diagram...")
    generate_combined_sankey(enrich_df, total_counts_df)

    print("Done.")


if __name__ == "__main__":
    main()
