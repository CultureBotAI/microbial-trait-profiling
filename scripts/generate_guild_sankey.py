import pandas as pd
import plotly.graph_objects as go
import numpy as np
from pathlib import Path

# Paths
ENRICHMENT_PATH = Path("data/derived/trait_enrichment_trait_study_pairs_ranked_q0.001.tsv")
EXPERIMENTS_PATH = Path("bsdbng/data/raw/experiments.csv")
OUTPUT_DIR = Path("bsdbng/data/derived/Sankey_Diagrams")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Configuration
GUILDS = {
    "Anaerobe": [
        "METPO:1000603",  # anaerobic
        "METPO:1000607"   # obligately anaerobic
    ],
    "Facultative-Aerobe": [
        "METPO:1000608",  # facultatively aerobic
        "METPO:1000602",  # aerobic
        "METPO:1000604"   # microaerophilic
    ]
}

MIN_PROFILES = 3
MIN_STUDIES = 2

def generate_guild_sankey(df, guild_name, direction, title):
    # Filter for guild traits and direction
    guild_traits = GUILDS[guild_name]
    g_df = df[(df["trait"].isin(guild_traits)) & (df["inc_dec"] == direction)].copy()
    
    if g_df.empty:
        print(f"No data for {guild_name} in {direction} direction.")
        return

    # Aggregate by profile to handle multiple guild hits in one profile
    profile_level = g_df.groupby(["condition", "body_site", "study_id", "profile_id"]).agg(
        max_fold=("fold", "max")
    ).reset_index()

    # Aggregate for the Sankey edges
    edge_stats = profile_level.groupby(["condition", "body_site"]).agg(
        n_profiles=("profile_id", "nunique"),
        s_studies=("study_id", "nunique"),
        mean_fold=("max_fold", "mean")
    ).reset_index()
    
    # Filter by support thresholds
    edge_stats = edge_stats[
        (edge_stats["n_profiles"] >= MIN_PROFILES) & 
        (edge_stats["s_studies"] >= MIN_STUDIES)
    ]
    
    if edge_stats.empty:
        print(f"No edges meet thresholds for {guild_name} ({direction}).")
        return

    # Order Conditions: thickest (top) to thinnest (bottom)
    cond_order = edge_stats.groupby("condition")["n_profiles"].sum().sort_values(ascending=False).index.tolist()
    # Order Body Sites: thickest (top) to thinnest (bottom)
    site_order = edge_stats.groupby("body_site")["n_profiles"].sum().sort_values(ascending=False).index.tolist()
    
    all_nodes = cond_order + site_order
    node_map = {node: i for i, node in enumerate(all_nodes)}
    
    sources = [node_map[c] for c in edge_stats["condition"]]
    targets = [node_map[s] for s in edge_stats["body_site"]]
    values = edge_stats["n_profiles"].tolist()
    
    # Color scale (Red-Blue)
    max_f = edge_stats["mean_fold"].max()
    min_f = edge_stats["mean_fold"].min()
    
    def get_color(val):
        norm = (val - min_f) / (max_f - min_f) if max_f != min_f else 0.5
        alpha = 0.6
        intensity = 0.3 + (norm * 0.7)
        r = int(180 * intensity)
        g = int(50 * (1-intensity))
        b = int(50 * (1-intensity))
        return f"rgba({r}, {g}, {b}, {alpha})"

    edge_colors = [get_color(f) for f in edge_stats["mean_fold"]]
    
    # Node Labels
    node_labels = []
    for node in all_nodes:
        n = edge_stats[edge_stats["condition"] == node]["n_profiles"].sum() if node in cond_order else edge_stats[edge_stats["body_site"] == node]["n_profiles"].sum()
        s = edge_stats[edge_stats["condition"] == node]["s_studies"].sum() if node in cond_order else edge_stats[edge_stats["body_site"] == node]["s_studies"].sum()
        node_labels.append(f"{node}<br>n={n}, s={s}")

    fig = go.Figure(data=[go.Sankey(
        node = dict(
            pad = 15, 
            thickness = 20, 
            line = dict(color = "black", width = 0.5), 
            label = node_labels, 
            color = "lightgray"
        ),
        link = dict(source = sources, target = targets, value = values, color = edge_colors),
        arrangement = "fixed" # Crucial to maintain our calculated order
    )])

    # Add custom positions to force vertical order
    # x=0 for conditions, x=1 for sites
    # y ranges from 0 (top) to 1 (bottom)
    node_x = [0.01] * len(cond_order) + [0.99] * len(site_order)
    
    # Calculate spacing for y
    def get_y_pos(counts):
        total = sum(counts)
        y = []
        curr = 0
        for c in counts:
            # Center of the block
            y.append((curr + c/2) / total)
            curr += c
        return y

    cond_y = get_y_pos(edge_stats.groupby("condition")["n_profiles"].sum().loc[cond_order].tolist())
    site_y = get_y_pos(edge_stats.groupby("body_site")["n_profiles"].sum().loc[site_order].tolist())
    node_y = cond_y + site_y

    fig.data[0].node.x = node_x
    fig.data[0].node.y = node_y

    fig.update_layout(
        title_text=f"<b>{title}</b><br><span style='font-size:10px'>Edges connect BugSigDB conditions to anatomical body sites. Width = profile support; Color = mean guild fold-change; s = unique studies.</span>",
        font_size=12,
        width=1200, height=1200 # Increased height for more labels
    )
    
    filename = OUTPUT_DIR / f"{guild_name.lower()}_sankey_{direction}.png"
    fig.write_image(str(filename), scale=2)
    print(f"Saved {filename}")

def main():
    print("Loading enrichment data...")
    enrich_df = pd.read_csv(ENRICHMENT_PATH, sep="\t", comment="#")
    
    print("Loading metadata...")
    meta_df = pd.read_csv(EXPERIMENTS_PATH, dtype=str)
    meta_df.columns = [c.strip().replace('"', '') for c in meta_df.columns]
    
    meta_df['study_key'] = "bsdb:" + meta_df['Study'].astype(str).str.strip()
    meta_df['exp_num'] = meta_df['Experiment page name'].str.extract(r'(\d+)', expand=False).fillna("1")
    meta_df['exp_key'] = meta_df['study_key'] + "-" + meta_df['exp_num']
    
    lookup = meta_df.drop_duplicates(subset='exp_key').set_index('exp_key')[['Body site', 'Condition']].to_dict('index')
    
    print("Mapping metadata...")
    enrich_df['body_site'] = enrich_df['experiment_id'].map(lambda x: lookup.get(x, {}).get('Body site', 'Unknown'))
    enrich_df['condition'] = enrich_df['experiment_id'].map(lambda x: lookup.get(x, {}).get('Condition', 'Unknown'))
    
    # Clean and filter
    enrich_df = enrich_df[(enrich_df['body_site'] != 'Unknown') & (enrich_df['condition'] != 'Unknown')].copy()
    
    print("Generating Sankey diagrams...")
    # 1. Anaerobe Guild Enrichment in Decreased Profiles
    generate_guild_sankey(enrich_df, "Anaerobe", "dec", "Condition-to-Body-Site Map: Anaerobe Guild Enrichment in Decreased Profiles")
    
    # 2. Facultative-Aerobe Guild Enrichment in Increased Profiles
    generate_guild_sankey(enrich_df, "Facultative-Aerobe", "inc", "Condition-to-Body-Site Map: Facultative-Aerobe Guild Enrichment in Increased Profiles")
    
    print("Done.")

if __name__ == "__main__":
    main()
