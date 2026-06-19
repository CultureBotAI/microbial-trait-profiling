"""Create a summary diagram connecting profile enrichment to guild assignment."""

from __future__ import annotations

from pathlib import Path
from textwrap import fill

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

REPO_ROOT = Path(__file__).resolve().parents[1]
VIS_DIR = REPO_ROOT / "trait-visualization"
INVESTIGATION_DIR = VIS_DIR / "trait-profile-investigation"
OUTPUT_PNG = INVESTIGATION_DIR / "profile_enrichment_guild_assignment_summary.png"
OUTPUT_SVG = INVESTIGATION_DIR / "profile_enrichment_guild_assignment_summary.svg"

SLATE_50 = "#f8fafc"
SLATE_100 = "#f1f5f9"
SLATE_200 = "#e2e8f0"
SLATE_400 = "#94a3b8"
SLATE_500 = "#64748b"
SLATE_700 = "#334155"
SLATE_900 = "#0f172a"
RED_50 = "#fef2f2"
RED_600 = "#dc2626"
BLUE_50 = "#eff6ff"
BLUE_600 = "#2563eb"
GREEN_50 = "#ecfdf5"
GREEN_700 = "#047857"
AMBER_50 = "#fffbeb"
AMBER_600 = "#d97706"


def read_tsv(path: Path, **kwargs: object) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", **kwargs)


def count_rows(path: Path) -> int:
    with path.open(encoding="utf-8") as handle:
        return max(0, sum(1 for _ in handle) - 1)


def format_count(value: int) -> str:
    return f"{value:,}"


def get_summary_counts() -> dict[str, int]:
    metadata = read_tsv(
        VIS_DIR / "study_trait_profile_study_metadata.tsv",
        usecols=["study_id"],
    )
    return {
        "profile_count": len(metadata),
        "study_count": metadata["study_id"].nunique(),
        "facultative_universe_count": count_rows(
            INVESTIGATION_DIR / "facultative-aerobe-v2" / "kg_microbe_candidate_taxon_universe.tsv"
        ),
        "facultative_condition_count": count_rows(
            INVESTIGATION_DIR / "facultative-aerobe-v2" / "condition_guild_summary.tsv"
        ),
        "gut_condition_count": count_rows(
            INVESTIGATION_DIR / "gut-anaerobe-v2" / "gut_condition_guild_summary.tsv"
        ),
        "facultative_target_taxa_count": count_rows(
            INVESTIGATION_DIR / "facultative-aerobe-v2" / "candidate_taxa_trait_presence.tsv"
        ),
        "gut_target_taxa_count": count_rows(
            INVESTIGATION_DIR / "gut-anaerobe-v2" / "candidate_taxa_trait_presence.tsv"
        ),
    }


def wrapped(text: str, width: int = 31) -> str:
    return "\n".join(fill(line, width=width) for line in text.split("\n"))


def add_box(
    ax: plt.Axes,
    *,
    xy: tuple[float, float],
    width: float,
    height: float,
    title: str,
    body: str,
    facecolor: str = "white",
    edgecolor: str = SLATE_200,
    title_color: str = SLATE_900,
    body_color: str = SLATE_700,
    body_width: int = 31,
    title_size: float = 10.5,
    body_size: float = 8.7,
) -> None:
    x, y = xy
    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.55,rounding_size=1.25",
        linewidth=1.2,
        edgecolor=edgecolor,
        facecolor=facecolor,
        zorder=2,
    )
    ax.add_patch(patch)
    ax.text(
        x + 1.0,
        y + height - 1.25,
        title,
        ha="left",
        va="top",
        fontsize=title_size,
        fontweight="bold",
        color=title_color,
        zorder=3,
    )
    ax.text(
        x + 1.0,
        y + height - 4.15,
        wrapped(body, body_width),
        ha="left",
        va="top",
        fontsize=body_size,
        color=body_color,
        linespacing=1.24,
        zorder=3,
    )


def add_arrow(
    ax: plt.Axes,
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    color: str = SLATE_400,
    connectionstyle: str = "arc3,rad=0",
    label: str | None = None,
    label_xy: tuple[float, float] | None = None,
) -> None:
    arrow = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=14,
        linewidth=1.8,
        color=color,
        connectionstyle=connectionstyle,
        shrinkA=6,
        shrinkB=6,
        zorder=1,
    )
    ax.add_patch(arrow)
    if label and label_xy:
        ax.text(
            label_xy[0],
            label_xy[1],
            label,
            ha="center",
            va="center",
            fontsize=7.8,
            color=color,
            fontweight="bold",
            zorder=4,
        )


def add_badge(
    ax: plt.Axes,
    xy: tuple[float, float],
    text: str,
    *,
    facecolor: str,
    edgecolor: str,
    text_color: str,
) -> None:
    x, y = xy
    patch = FancyBboxPatch(
        (x, y),
        16.6,
        2.6,
        boxstyle="round,pad=0.25,rounding_size=1.3",
        linewidth=1.0,
        edgecolor=edgecolor,
        facecolor=facecolor,
        zorder=3,
    )
    ax.add_patch(patch)
    ax.text(
        x + 8.3,
        y + 1.3,
        text,
        ha="center",
        va="center",
        fontsize=8.0,
        color=text_color,
        fontweight="bold",
        zorder=4,
    )


def create_diagram() -> None:
    counts = get_summary_counts()
    INVESTIGATION_DIR.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(15.2, 8.6), dpi=180)
    fig.patch.set_facecolor(SLATE_50)
    ax.set_facecolor(SLATE_50)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 62)
    ax.axis("off")

    ax.text(
        3.5,
        59.4,
        "Profile enrichment connects to guild assignment through trait-level evidence",
        fontsize=21,
        fontweight="bold",
        color=SLATE_900,
        ha="left",
        va="top",
    )
    ax.text(
        3.5,
        56.4,
        (
            "BugSigDB curated profile directions are converted to trait z-scores, "
            "contrasted within metadata groups, then summarized as v2 guild signals."
        ),
        fontsize=10.2,
        color=SLATE_500,
        ha="left",
        va="top",
    )

    add_box(
        ax,
        xy=(3.0, 42.4),
        width=16.8,
        height=10.5,
        title="1. Profiles",
        body=(
            f"{format_count(counts['profile_count'])} profiles\n"
            f"{format_count(counts['study_count'])} studies\n"
            "curated BugSigDB inc/dec labels"
        ),
        facecolor="white",
        body_width=25,
    )
    add_box(
        ax,
        xy=(22.0, 42.4),
        width=16.8,
        height=10.5,
        title="2. Trait matrix",
        body=(
            "Taxon occurrences are aggregated with KG-Microbe traits and standardized "
            "as profile trait z-scores."
        ),
        facecolor="white",
        body_width=26,
    )
    add_box(
        ax,
        xy=(41.0, 42.4),
        width=16.8,
        height=10.5,
        title="3. Enrichment",
        body="Within each condition or body site, compare focus direction against its comparator.",
        facecolor="white",
        body_width=26,
    )
    add_box(
        ax,
        xy=(60.0, 42.4),
        width=16.8,
        height=10.5,
        title="4. Evidence",
        body="Keep traits with positive focus mean z and positive focus-minus-comparator contrast.",
        facecolor="white",
        body_width=27,
    )
    add_box(
        ax,
        xy=(79.0, 42.4),
        width=17.8,
        height=10.5,
        title="5. Guild call",
        body="Selected trait sets define a named guild score and condition-level signal.",
        facecolor="white",
        body_width=27,
    )

    for start_x, end_x in [(19.8, 22.0), (38.8, 41.0), (57.8, 60.0), (76.8, 79.0)]:
        add_arrow(ax, (start_x, 47.65), (end_x, 47.65))

    add_box(
        ax,
        xy=(8.0, 21.1),
        width=39.0,
        height=14.7,
        title="Facultative-aerobe v2 guild",
        body=(
            "Focus direction: increased profiles\n"
            "Guild score: mean of facultatively aerobic, aerobic, and microaerophilic z-scores\n"
            f"Condition summaries: {format_count(counts['facultative_condition_count'])}; "
            f"target taxa: {format_count(counts['facultative_target_taxa_count'])}; "
            f"candidate universe: {format_count(counts['facultative_universe_count'])} taxa"
        ),
        facecolor=RED_50,
        edgecolor="#fecaca",
        title_color="#991b1b",
        body_width=55,
    )
    add_box(
        ax,
        xy=(53.0, 21.1),
        width=39.0,
        height=14.7,
        title="Gut anaerobe v2 guild",
        body=(
            "Focus direction: decreased profiles in gut context\n"
            "Guild score: mean of anaerobic and obligately anaerobic z-scores\n"
            f"Condition summaries: {format_count(counts['gut_condition_count'])}; "
            f"target taxa: {format_count(counts['gut_target_taxa_count'])}; "
            "gut labels include feces, stool, colon, caecum, gut, intestine"
        ),
        facecolor=BLUE_50,
        edgecolor="#bfdbfe",
        title_color="#1d4ed8",
        body_width=55,
    )

    add_arrow(
        ax,
        (83.0, 42.3),
        (30.0, 35.9),
        color=RED_600,
        connectionstyle="arc3,rad=0.12",
        label="increased",
        label_xy=(67.8, 40.7),
    )
    add_arrow(
        ax,
        (88.0, 42.3),
        (72.0, 35.9),
        color=BLUE_600,
        connectionstyle="arc3,rad=-0.12",
        label="decreased gut",
        label_xy=(81.0, 37.6),
    )

    add_box(
        ax,
        xy=(5.0, 5.2),
        width=26.5,
        height=11.0,
        title="Guardrail",
        body=(
            "Target-taxon recurrence is an audit layer. It is not the same as a "
            "recomputed guild score."
        ),
        facecolor=AMBER_50,
        edgecolor="#fde68a",
        title_color="#92400e",
        body_width=34,
    )
    add_box(
        ax,
        xy=(35.0, 5.2),
        width=30.0,
        height=11.0,
        title="Taxon support assignment",
        body=(
            "Fixed focal taxa use NCBI IDs for direction and occurrence counts. "
            "Data-inferred support intersects enriched traits with observed taxa carrying "
            "KG-Microbe annotations."
        ),
        facecolor=GREEN_50,
        edgecolor="#a7f3d0",
        title_color=GREEN_700,
        body_width=42,
    )
    add_box(
        ax,
        xy=(68.5, 5.2),
        width=26.5,
        height=11.0,
        title="Generated outputs",
        body=(
            "Condition guild signals, guild-trait shifts, body-site diagrams, "
            "target-taxon condition counts, and inferred evidence TSVs."
        ),
        facecolor=SLATE_100,
        edgecolor=SLATE_200,
        title_color=SLATE_900,
        body_width=35,
    )

    add_arrow(
        ax,
        (29.0, 21.0),
        (47.0, 16.3),
        color=RED_600,
        connectionstyle="arc3,rad=-0.08",
    )
    add_arrow(
        ax,
        (72.0, 21.0),
        (54.0, 16.3),
        color=BLUE_600,
        connectionstyle="arc3,rad=0.08",
    )
    add_arrow(ax, (65.0, 10.7), (68.5, 10.7), color=SLATE_400)
    add_arrow(ax, (31.5, 10.7), (35.0, 10.7), color=AMBER_600)

    add_badge(
        ax,
        (8.5, 37.2),
        "guild-specific panels",
        facecolor=RED_50,
        edgecolor="#fecaca",
        text_color="#991b1b",
    )
    add_badge(
        ax,
        (27.0, 37.2),
        "taxon audit panels",
        facecolor=AMBER_50,
        edgecolor="#fde68a",
        text_color="#92400e",
    )
    add_badge(
        ax,
        (45.5, 37.2),
        "inferred evidence",
        facecolor=GREEN_50,
        edgecolor="#a7f3d0",
        text_color=GREEN_700,
    )

    ax.text(
        3.5,
        2.0,
        (
            "Summary diagram generated from current v2 investigation tables. "
            "Profile enrichment is the focus-minus-comparator trait z-score contrast; "
            "guild assignment is based on selected KG-Microbe trait sets."
        ),
        fontsize=8.4,
        color=SLATE_500,
        ha="left",
        va="bottom",
    )

    fig.savefig(OUTPUT_PNG, dpi=220, bbox_inches="tight", facecolor=fig.get_facecolor())
    fig.savefig(OUTPUT_SVG, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)


def main() -> None:
    create_diagram()
    print(f"Wrote: {OUTPUT_PNG}")
    print(f"Wrote: {OUTPUT_SVG}")


if __name__ == "__main__":
    main()
