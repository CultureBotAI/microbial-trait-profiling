"""Plot species-level PD taxa against the four positive-side lollipop traits."""

from __future__ import annotations

import math
from pathlib import Path
from textwrap import wrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUTPUT_DIR = (
    Path(__file__).resolve().parents[1]
    / "trait-visualization/trait-profile-investigation/PD-cellulose"
)
TRAITS = [
    "moderately halophilic",
    "facultatively aerobic",
    "arsenate reduction",
    "coccus shaped",
]
STEM = "pd_species_four_increased_traits_inc_dec_heatmap"
COLOR_LIMIT = 3.0


def prepare_cells(carriers: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Keep observed species; preserve unannotated, zero-denominator, and finite cells."""
    selected = carriers[carriers["trait"].isin(TRAITS)].copy()
    if selected.duplicated(["taxon_id", "trait"]).any():
        raise ValueError("Duplicate taxon/trait rows in carrier table")
    species = selected[
        selected["taxonomic_rank"].eq("species") & (selected["increased_profiles"] > 0)
    ].copy()
    if species.empty:
        raise ValueError("No species-level increased PD taxa carry the selected traits")
    row_order = (
        species.groupby(["taxon_id", "taxon_name"], as_index=False)
        .agg(
            annotated_trait_count=("trait", "nunique"),
            increased_profiles=("increased_profiles", "max"),
            decreased_profiles=("decreased_profiles", "max"),
        )
        .sort_values(
            ["annotated_trait_count", "increased_profiles", "taxon_name"],
            ascending=[False, False, True],
        )
    )
    grid = pd.MultiIndex.from_product([row_order["taxon_id"], TRAITS], names=["taxon_id", "trait"])
    cells = species.set_index(["taxon_id", "trait"]).reindex(grid).reset_index()
    cells["taxon_name"] = cells["taxon_id"].map(row_order.set_index("taxon_id")["taxon_name"])
    cells["taxonomic_rank"] = "species"
    cells["trait_annotated"] = cells["annotation_count"].notna()
    for column in ["increased_profiles", "decreased_profiles"]:
        cells[column] = cells[column].astype("Int64")
    cells["inc_dec_counts"] = ""
    cells["inc_dec_ratio"] = np.nan
    cells["log2_inc_dec_ratio"] = np.nan
    for index, row in cells[cells["trait_annotated"]].iterrows():
        inc, dec = int(row["increased_profiles"]), int(row["decreased_profiles"])
        ratio = inc / dec if dec else math.inf
        cells.loc[index, "inc_dec_counts"] = f"{inc}/{dec}"
        cells.loc[index, "inc_dec_ratio"] = ratio
        cells.loc[index, "log2_inc_dec_ratio"] = math.log2(ratio)
    cells["color_log2_ratio"] = cells["log2_inc_dec_ratio"].clip(-COLOR_LIMIT, COLOR_LIMIT)
    cells["blank_reason"] = np.where(
        cells["trait_annotated"], "", "no annotation in the local taxa-trait matrix"
    )
    return cells, row_order


def plot_heatmap(cells: pd.DataFrame, rows: pd.DataFrame, path: Path) -> None:
    """Color exact inc/dec ratios and label them with their unsmoothed count pairs."""
    values = cells.pivot(index="taxon_id", columns="trait", values="color_log2_ratio")
    labels = cells.pivot(index="taxon_id", columns="trait", values="inc_dec_counts")
    values, labels = [
        frame.reindex(index=rows["taxon_id"], columns=TRAITS) for frame in [values, labels]
    ]
    height = max(7.0, len(rows) * 0.32 + 2.8)
    fig, ax = plt.subplots(figsize=(12.2, height), facecolor="white")
    fig.subplots_adjust(left=0.37, right=0.89, top=1 - 1.3 / height, bottom=1.1 / height)
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#e5e7eb")
    image = ax.imshow(
        values.to_numpy(float), cmap=cmap, vmin=-COLOR_LIMIT, vmax=COLOR_LIMIT, aspect="auto"
    )
    ax.set_xticks(range(4), ["\n".join(wrap(trait, width=18)) for trait in TRAITS], fontsize=10)
    ax.set_yticks(range(len(rows)), rows["taxon_name"], fontsize=9)
    ax.tick_params(length=0, pad=7)
    ax.set_xticks(np.arange(-0.5, 4, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(rows), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1)
    ax.tick_params(which="minor", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    for y, taxon_id in enumerate(rows["taxon_id"]):
        for x, trait in enumerate(TRAITS):
            value = values.loc[taxon_id, trait]
            if pd.isna(value):
                continue
            ax.text(
                x,
                y,
                labels.loc[taxon_id, trait],
                ha="center",
                va="center",
                fontsize=9,
                color="white" if abs(value) >= 2 else "#111827",
            )
    bar = fig.colorbar(image, ax=ax, pad=0.025, fraction=0.04, shrink=0.5, extend="both")
    bar.set_ticks([-3, -2, -1, 0, 1, 2, 3], labels=["1/8", "1/4", "1/2", "1", "2", "4", "8"])
    bar.set_label("Inc/dec profile-count ratio (log2 color scale)")
    fig.text(
        0.03,
        1 - 0.12 / height,
        "PD species: four increased-direction traits",
        fontsize=17,
        fontweight="bold",
        va="top",
    )
    subtitle = (
        f"{len(rows)} species reported in increased PD fecal profiles. Cells = increased/decreased "
        "profile counts; red = inc > dec, blue = dec > inc. Counts include all PD profiles "
        "containing the annotated taxon, without a profile z-score filter. Gray = no trait "
        "annotation. Ratios above 8 (including inc/0) saturate red; no pseudocount is added."
    )
    fig.text(
        0.03,
        1 - 0.5 / height,
        "\n".join(wrap(subtitle, width=115)),
        fontsize=9,
        color="#4b5563",
        va="top",
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=240, bbox_inches="tight", facecolor="white")
    fig.savefig(path.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    """Generate the species heatmap and exact cell-value audit table."""
    carriers = pd.read_csv(OUTPUT_DIR / "pd_six_trait_taxon_carriers.tsv", sep="\t")
    cells, rows = prepare_cells(carriers)
    plot_heatmap(cells, rows, OUTPUT_DIR / f"{STEM}.png")
    cells.to_csv(OUTPUT_DIR / f"{STEM}_cells.tsv", sep="\t", index=False, na_rep="NA")
    note = (
        "# PD Species-by-Four-Trait Heatmap\n\n"
        f"![Heatmap]({STEM}.png)\n\n"
        f"Includes all {len(rows)} reported species occurring in at least one increased PD fecal "
        "profile and annotated for at least one of the four positive-side lollipop traits. "
        "Species is the lowest rank available in the carrier data. Genus and higher-rank "
        "observations are excluded, not reassigned to species. Labels such as 'sp.' retain "
        "the species rank and names supplied by the source.\n\n"
        "Cells display inc/dec counts of distinct profiles containing that annotated taxon, "
        "out of the source cohort of 130 increased and 123 decreased PD fecal profiles. "
        "These are raw occurrence counts, not abundance estimates or normalized prevalence "
        "ratios. No positive-z or p-value filter is applied. An inc/0 ratio is infinite; "
        "its color saturates at the upper scale limit without adding a pseudocount. "
        "Gray indicates an absent annotation, not a zero observation count "
        "or proven trait absence. "
        "Counts repeat across traits carried by the same species because the underlying "
        "taxon occurrences are the same.\n\n"
        "KG-Microbe annotations may originate at the reported taxon or its descendants. "
        "Multiple annotated traits need not be experimentally verified in the same strain. "
        "The coccus-shaped lollipop effect is nearly neutral (+0.017).\n\n"
        f"[Exact cell counts and ratios]({STEM}_cells.tsv) | "
        f"[Vector figure]({STEM}.svg)\n"
    )
    (OUTPUT_DIR / f"{STEM}.md").write_text(note, encoding="utf-8")
    print(f"Wrote {len(rows)} species x 4 traits; {cells['trait_annotated'].sum()} annotated cells")
    print(OUTPUT_DIR / f"{STEM}.png")


if __name__ == "__main__":
    main()
