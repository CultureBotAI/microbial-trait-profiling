"""Create a summary diagram connecting profile enrichment to guild assignment."""

from __future__ import annotations

from pathlib import Path
from textwrap import fill

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

REPO_ROOT = Path(__file__).resolve().parents[1]
VIS_DIR = REPO_ROOT / "trait-visualization"
INVESTIGATION_DIR = VIS_DIR / "trait-profile-investigation"
OUTPUT_PNG = INVESTIGATION_DIR / "profile_enrichment_guild_assignment_summary.png"
OUTPUT_SVG = INVESTIGATION_DIR / "profile_enrichment_guild_assignment_summary.svg"
ANAEROBE_GUILD_DEFINITION_TABLE_PNG = INVESTIGATION_DIR / "anaerobe_v2_guild_definition_table.png"
ANAEROBE_GUILD_DEFINITION_TABLE_SVG = INVESTIGATION_DIR / "anaerobe_v2_guild_definition_table.svg"
PRE_GRID_OUTPUT_PNG = (
    INVESTIGATION_DIR / "profile_enrichment_guild_assignment_pre_direction_grid_summary.png"
)
GRID_OUTPUT_TSV = INVESTIGATION_DIR / "profile_enrichment_guild_direction_grid.tsv"
ANAEROBE_GUILD_SELECTION_TSV = (
    INVESTIGATION_DIR / "anaerobe-depletion-v2" / "guild_trait_selection.tsv"
)

AEROBE_GUILD_TRAITS = {
    "METPO:1000608": "facultatively aerobic",
    "METPO:1000602": "aerobic",
    "METPO:1000604": "microaerophilic",
}
ANAEROBE_GUILD_TRAITS = {
    "METPO:1000603": "anaerobic",
    "METPO:1000607": "obligately anaerobic",
}
ANAEROBE_GUILD_TRAIT_DEFINITIONS = {
    "METPO:1000603": (
        "Taxa carrying anaerobic oxygen-tolerance annotations; broad counterpart "
        "to oxygen-positive traits."
    ),
    "METPO:1000607": (
        "Stricter anaerobe annotation for taxa expected to rely on oxygen-free growth contexts."
    ),
}
ANAEROBE_GUILD_TRAIT_ROLES = {
    "METPO:1000603": "Core trait for anaerobe enrichment in decreased profiles.",
    "METPO:1000607": (
        "Keeps the selected trait set specific to oxygen tolerance rather than "
        "fermentation metabolism."
    ),
}
EXCLUDED_CONDITION_LABELS = {
    "Body fat distribution",
    "Clinical treatment",
    "Exposure",
    "Health trait",
    "Microbiome",
    "Population",
    "Physical activity",
}
MIN_CONDITION_PROFILES = 5
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


def clean_text(value: object) -> str:
    text = str(value or "").strip()
    if not text or text.lower() == "nan" or text.upper() == "NA":
        return "unknown"
    return text


def is_excluded_condition_label(value: object) -> bool:
    label = clean_text(value)
    broad_terms = ("exposure", "measurement", "population")
    lower_label = label.lower()
    return label in EXCLUDED_CONDITION_LABELS or any(term in lower_label for term in broad_terms)


def count_rows(path: Path) -> int:
    if not path.exists():
        return 0
    with path.open(encoding="utf-8") as handle:
        return max(0, sum(1 for _ in handle) - 1)


def format_count(value: int) -> str:
    return f"{value:,}"


def count_unique_conditions(path: Path) -> int:
    if not path.exists():
        return 0
    return read_tsv(path, usecols=["condition"])["condition"].nunique()


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


def write_direction_grid_table() -> None:
    """Write a unified 2x2 trait-set-by-direction condition summary table."""
    matrix, metadata, _ = load_study_heatmap()
    metadata = metadata.copy()
    metadata["condition_clean"] = metadata["condition"].map(clean_text)

    guilds = [
        {
            "guild": "aerobe_v2",
            "traits": list(AEROBE_GUILD_TRAITS.values()),
            "parent_direction": "increased",
        },
        {
            "guild": "anaerobe_v2_depletion",
            "traits": list(ANAEROBE_GUILD_TRAITS.values()),
            "parent_direction": "decreased",
        },
    ]
    direction_to_inc_dec = {"increased": "inc", "decreased": "dec"}
    comparator_direction = {"increased": "decreased", "decreased": "increased"}
    condition_counts = metadata["condition_clean"].value_counts()
    conditions = sorted(
        condition
        for condition, profile_count in condition_counts.items()
        if profile_count >= MIN_CONDITION_PROFILES and not is_excluded_condition_label(condition)
    )

    records: list[dict[str, object]] = []
    for guild in guilds:
        profiles = metadata.copy()
        available_traits = [trait for trait in guild["traits"] if trait in matrix.columns]
        profiles["guild_z_mean"] = matrix.loc[profiles.index, available_traits].mean(axis=1)

        for condition in conditions:
            condition_profiles = profiles[profiles["condition_clean"] == condition].copy()
            for profile_direction in ["increased", "decreased"]:
                comparator = comparator_direction[profile_direction]
                direction_profiles = condition_profiles[
                    condition_profiles["inc_dec"] == direction_to_inc_dec[profile_direction]
                ]
                comparator_profiles = condition_profiles[
                    condition_profiles["inc_dec"] == direction_to_inc_dec[comparator]
                ]
                mean_guild_z = direction_profiles["guild_z_mean"].mean()
                comparator_mean_guild_z = comparator_profiles["guild_z_mean"].mean()
                direction_contrast = mean_guild_z - comparator_mean_guild_z
                records.append(
                    {
                        "guild": guild["guild"],
                        "profile_direction": profile_direction,
                        "condition": condition,
                        "profile_count": len(direction_profiles),
                        "study_count": direction_profiles["study_id"].nunique(),
                        "mean_guild_z": mean_guild_z,
                        "comparator_mean_guild_z": comparator_mean_guild_z,
                        "direction_contrast": direction_contrast,
                        "parent_focus": profile_direction == guild["parent_direction"],
                    }
                )

    pd.DataFrame(records).to_csv(GRID_OUTPUT_TSV, sep="\t", index=False)


def get_summary_counts() -> dict[str, int]:
    metadata = read_tsv(
        VIS_DIR / "study_trait_profile_study_metadata.tsv",
        usecols=["study_id"],
    )
    return {
        "profile_count": len(metadata),
        "study_count": metadata["study_id"].nunique(),
        "direction_grid_condition_count": count_unique_conditions(GRID_OUTPUT_TSV),
        "direction_grid_row_count": count_rows(GRID_OUTPUT_TSV),
        "facultative_universe_count": count_rows(
            INVESTIGATION_DIR / "facultative-aerobe-v2" / "kg_microbe_candidate_taxon_universe.tsv"
        ),
        "facultative_condition_count": count_rows(
            INVESTIGATION_DIR / "facultative-aerobe-v2" / "condition_guild_summary.tsv"
        ),
        "anaerobe_condition_count": count_rows(
            INVESTIGATION_DIR / "anaerobe-depletion-v2" / "condition_guild_summary.tsv"
        ),
        "anaerobe_universe_count": count_rows(
            INVESTIGATION_DIR / "anaerobe-depletion-v2" / "kg_microbe_candidate_taxon_universe.tsv"
        ),
        "facultative_target_taxa_count": count_rows(
            INVESTIGATION_DIR / "facultative-aerobe-v2" / "candidate_taxa_trait_presence.tsv"
        ),
        "anaerobe_target_taxa_count": count_rows(
            INVESTIGATION_DIR / "anaerobe-depletion-v2" / "candidate_taxa_trait_presence.tsv"
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


def add_compact_direction_grid(
    ax: plt.Axes, *, xy: tuple[float, float], width: float, height: float
) -> None:
    """Add the 2x2 direction grid with parent signals and guardrail cells."""
    x, y = xy
    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.45,rounding_size=1.25",
        linewidth=1.2,
        edgecolor=SLATE_200,
        facecolor="white",
        zorder=2,
    )
    ax.add_patch(patch)
    ax.text(
        x + 1.0,
        y + height - 1.15,
        "2x2 Direction Grid",
        ha="left",
        va="top",
        fontsize=10.2,
        fontweight="bold",
        color=SLATE_900,
        zorder=4,
    )

    grid_x = x + 1.0
    grid_y = y + 0.85
    grid_w = width - 2.0
    grid_h = height - 4.1
    header_h = 1.35
    row_label_w = 8.4
    cell_w = (grid_w - row_label_w) / 2
    cell_h = (grid_h - header_h) / 2
    ax.text(
        grid_x + row_label_w + cell_w / 2,
        grid_y + grid_h - 0.15,
        "inc",
        ha="center",
        va="top",
        fontsize=7.6,
        color=SLATE_500,
        fontweight="bold",
    )
    ax.text(
        grid_x + row_label_w + cell_w * 1.5,
        grid_y + grid_h - 0.15,
        "dec",
        ha="center",
        va="top",
        fontsize=7.6,
        color=SLATE_500,
        fontweight="bold",
    )
    row_labels = [("Aerobe", RED_600), ("Anaerobe", BLUE_600)]
    cells = [
        ("parent:\naerobe enriched", RED_50, "#fecaca", "#991b1b"),
        ("guardrail:\ngenerally high?", SLATE_100, SLATE_200, SLATE_700),
        ("guardrail:\nmixed signal?", SLATE_100, SLATE_200, SLATE_700),
        ("parent:\nanaerobe enriched", BLUE_50, "#bfdbfe", "#1d4ed8"),
    ]
    for row_index, (label, label_color) in enumerate(row_labels):
        cell_y = grid_y + (1 - row_index) * cell_h
        ax.text(
            grid_x + row_label_w - 0.5,
            cell_y + cell_h / 2,
            label,
            ha="right",
            va="center",
            fontsize=7.8,
            color=label_color,
            fontweight="bold",
            zorder=5,
        )
        for col_index in range(2):
            cell_index = row_index * 2 + col_index
            text, facecolor, edgecolor, text_color = cells[cell_index]
            cell_x = grid_x + row_label_w + col_index * cell_w
            rect = FancyBboxPatch(
                (cell_x + 0.25, cell_y + 0.22),
                cell_w - 0.5,
                cell_h - 0.44,
                boxstyle="round,pad=0.18,rounding_size=0.55",
                linewidth=0.8,
                edgecolor=edgecolor,
                facecolor=facecolor,
                zorder=3,
            )
            ax.add_patch(rect)
            ax.text(
                cell_x + cell_w / 2,
                cell_y + cell_h / 2,
                text,
                ha="center",
                va="center",
                fontsize=7.3,
                color=text_color,
                fontweight="bold",
                linespacing=1.1,
                zorder=4,
            )


def add_grid_cell(
    ax: plt.Axes,
    *,
    xy: tuple[float, float],
    width: float,
    height: float,
    title: str,
    body: str,
    facecolor: str,
    edgecolor: str,
    title_color: str,
    body_color: str = SLATE_700,
    linewidth: float = 1.1,
) -> None:
    x, y = xy
    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.25,rounding_size=0.75",
        linewidth=linewidth,
        edgecolor=edgecolor,
        facecolor=facecolor,
        zorder=2,
    )
    ax.add_patch(patch)
    ax.text(
        x + 0.7,
        y + height - 0.85,
        title,
        ha="left",
        va="top",
        fontsize=8.0,
        fontweight="bold",
        color=title_color,
        zorder=3,
    )
    ax.text(
        x + 0.7,
        y + height - 2.65,
        wrapped(body, 38),
        ha="left",
        va="top",
        fontsize=6.9,
        color=body_color,
        linespacing=1.12,
        zorder=3,
    )


def add_large_direction_grid(ax: plt.Axes, counts: dict[str, int]) -> None:
    grid_x = 8.0
    grid_y = 19.0
    grid_w = 84.0
    grid_h = 16.5
    header_h = 3.0
    label_w = 19.5
    row_h = (grid_h - header_h) / 2
    col_w = (grid_w - label_w) / 2
    top_row_y = grid_y + row_h
    header_y = grid_y + grid_h - header_h
    inc_x = grid_x + label_w
    dec_x = inc_x + col_w

    outer = FancyBboxPatch(
        (grid_x, grid_y),
        grid_w,
        grid_h,
        boxstyle="round,pad=0.35,rounding_size=1.0",
        linewidth=1.1,
        edgecolor=SLATE_200,
        facecolor="white",
        zorder=1,
    )
    ax.add_patch(outer)
    ax.text(
        grid_x + 1.0,
        header_y + 1.95,
        "2 x 2 direction grid",
        ha="left",
        va="center",
        fontsize=9.4,
        fontweight="bold",
        color=SLATE_900,
        zorder=3,
    )
    ax.text(
        inc_x + col_w / 2,
        header_y + 1.95,
        "increased profiles",
        ha="center",
        va="center",
        fontsize=8.4,
        fontweight="bold",
        color=SLATE_900,
        zorder=3,
    )
    ax.text(
        dec_x + col_w / 2,
        header_y + 1.95,
        "decreased profiles",
        ha="center",
        va="center",
        fontsize=8.4,
        fontweight="bold",
        color=SLATE_900,
        zorder=3,
    )

    ax.text(
        grid_x + 1.0,
        top_row_y + row_h / 2,
        "Aerobe\ntrait set",
        ha="left",
        va="center",
        fontsize=8.2,
        fontweight="bold",
        color="#991b1b",
        linespacing=1.15,
        zorder=3,
    )
    ax.text(
        grid_x + 1.0,
        grid_y + row_h / 2,
        "Anaerobe\ntrait set",
        ha="left",
        va="center",
        fontsize=8.2,
        fontweight="bold",
        color="#1d4ed8",
        linespacing=1.15,
        zorder=3,
    )

    add_grid_cell(
        ax,
        xy=(inc_x, top_row_y),
        width=col_w,
        height=row_h,
        title="Parent focus: aerobe enriched",
        body=(
            "Aerobe traits enriched in increased profiles.\n"
            f"{format_count(counts['direction_grid_condition_count'])} BugSigDB "
            "condition-label rows (not experiments)."
        ),
        facecolor=RED_50,
        edgecolor=RED_600,
        title_color="#991b1b",
        linewidth=1.8,
    )
    add_grid_cell(
        ax,
        xy=(dec_x, top_row_y),
        width=col_w,
        height=row_h,
        title="Off-diagonal check",
        body=(
            "Aerobe trait score in decreased profiles.\n"
            "Used as the comparator/reversal check."
        ),
        facecolor=SLATE_100,
        edgecolor=SLATE_200,
        title_color=SLATE_700,
    )
    add_grid_cell(
        ax,
        xy=(inc_x, grid_y),
        width=col_w,
        height=row_h,
        title="Off-diagonal check",
        body=(
            "Anaerobe trait score in increased profiles.\n"
            "Used as the comparator/reversal check."
        ),
        facecolor=SLATE_100,
        edgecolor=SLATE_200,
        title_color=SLATE_700,
    )
    add_grid_cell(
        ax,
        xy=(dec_x, grid_y),
        width=col_w,
        height=row_h,
        title="Parent focus: anaerobe enriched",
        body=(
            "Anaerobe traits enriched in decreased profiles.\n"
            f"{format_count(counts['direction_grid_condition_count'])} BugSigDB "
            "condition-label rows (not experiments)."
        ),
        facecolor=BLUE_50,
        edgecolor=BLUE_600,
        title_color="#1d4ed8",
        linewidth=1.8,
    )

    ax.text(
        grid_x + 1.0,
        grid_y - 1.0,
        f"Unified table: {GRID_OUTPUT_TSV.name}; rows are guild x direction x condition.",
        ha="left",
        va="top",
        fontsize=7.2,
        color=SLATE_500,
        zorder=3,
    )


def format_optional_count(value: object) -> str:
    try:
        return format_count(int(float(str(value))))
    except ValueError:
        return "n/a"


def load_anaerobe_guild_definition_rows() -> list[dict[str, str]]:
    selection = pd.DataFrame()
    if ANAEROBE_GUILD_SELECTION_TSV.exists():
        selection = read_tsv(ANAEROBE_GUILD_SELECTION_TSV).set_index("trait_id")

    rows: list[dict[str, str]] = []
    for trait_id, trait_name in ANAEROBE_GUILD_TRAITS.items():
        source = selection.loc[trait_id] if trait_id in selection.index else pd.Series(dtype=str)
        taxon_count = source.get("kg_microbe_taxa_with_trait_count", "")
        rows.append(
            {
                "trait": f"{trait_name}\n{trait_id}",
                "definition": ANAEROBE_GUILD_TRAIT_DEFINITIONS[trait_id],
                "role": ANAEROBE_GUILD_TRAIT_ROLES[trait_id],
                "taxa": format_optional_count(taxon_count),
            }
        )
    return rows


def add_definition_table_cell(
    ax: plt.Axes,
    *,
    x: float,
    y: float,
    width: float,
    height: float,
    text: str,
    facecolor: str,
    edgecolor: str,
    color: str = SLATE_700,
    fontsize: float = 7.2,
    fontweight: str = "normal",
    ha: str = "left",
    wrap_width: int = 34,
) -> None:
    rect = Rectangle(
        (x, y),
        width,
        height,
        linewidth=0.65,
        edgecolor=edgecolor,
        facecolor=facecolor,
        zorder=2,
    )
    ax.add_patch(rect)
    text_x = x + width / 2 if ha == "center" else x + 0.75
    ax.text(
        text_x,
        y + height / 2,
        wrapped(text, wrap_width),
        ha=ha,
        va="center",
        fontsize=fontsize,
        fontweight=fontweight,
        color=color,
        linespacing=1.1,
        zorder=3,
    )


def create_anaerobe_guild_definition_table(
    *,
    output_png: Path = ANAEROBE_GUILD_DEFINITION_TABLE_PNG,
    output_svg: Path | None = ANAEROBE_GUILD_DEFINITION_TABLE_SVG,
) -> None:
    """Create a compact diagram-ready table for the anaerobe v2 selected traits."""
    rows = load_anaerobe_guild_definition_rows()

    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(15.2, 4.2), dpi=180)
    fig.patch.set_facecolor(SLATE_50)
    ax.set_facecolor(SLATE_50)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 30)
    ax.axis("off")

    ax.text(
        3.5,
        27.4,
        "Anaerobe v2 trait-enrichment guild definitions",
        fontsize=18,
        fontweight="bold",
        color=SLATE_900,
        ha="left",
        va="top",
    )
    ax.text(
        3.5,
        24.9,
        (
            "Parent focus: decreased profiles. Trait score is the row-wise mean of "
            "the selected oxygen-tolerance trait z-scores."
        ),
        fontsize=9.4,
        color=SLATE_500,
        ha="left",
        va="top",
    )
    add_badge(
        ax,
        (78.5, 24.3),
        "anaerobe enriched in dec",
        facecolor=BLUE_50,
        edgecolor="#bfdbfe",
        text_color="#1d4ed8",
    )

    table_x = 3.5
    table_y = 5.8
    table_w = 93.0
    header_h = 3.1
    row_h = 6.35
    table_h = header_h + row_h * len(rows)
    columns = [
        ("Trait", 18.6, 18),
        ("Diagram definition", 34.0, 43),
        ("Guild role", 30.8, 39),
        ("KG taxa", 9.6, 9),
    ]

    outer = FancyBboxPatch(
        (table_x, table_y),
        table_w,
        table_h,
        boxstyle="round,pad=0.35,rounding_size=0.9",
        linewidth=1.1,
        edgecolor=SLATE_200,
        facecolor="white",
        zorder=1,
    )
    ax.add_patch(outer)

    x = table_x
    for header, width, wrap_width in columns:
        add_definition_table_cell(
            ax,
            x=x,
            y=table_y + row_h * len(rows),
            width=width,
            height=header_h,
            text=header,
            facecolor=BLUE_50,
            edgecolor=SLATE_200,
            color="#1d4ed8",
            fontsize=7.4,
            fontweight="bold",
            ha="center" if header == "KG taxa" else "left",
            wrap_width=wrap_width,
        )
        x += width

    for row_index, row in enumerate(rows):
        y = table_y + row_h * (len(rows) - 1 - row_index)
        facecolor = "white" if row_index == 0 else SLATE_50
        x = table_x
        for key, width, wrap_width in [
            ("trait", columns[0][1], columns[0][2]),
            ("definition", columns[1][1], columns[1][2]),
            ("role", columns[2][1], columns[2][2]),
            ("taxa", columns[3][1], columns[3][2]),
        ]:
            add_definition_table_cell(
                ax,
                x=x,
                y=y,
                width=width,
                height=row_h,
                text=row[key],
                facecolor=facecolor,
                edgecolor=SLATE_200,
                color="#1d4ed8" if key == "trait" else SLATE_700,
                fontsize=7.0 if key != "trait" else 7.2,
                fontweight="bold" if key == "trait" else "normal",
                ha="center" if key == "taxa" else "left",
                wrap_width=wrap_width,
            )
            x += width

    ax.text(
        3.5,
        2.6,
        (
            "This table defines the v2 trait-enrichment score. Focal target taxa are "
            "recurrence audit examples and do not define the formal trait set."
        ),
        fontsize=8.2,
        color=SLATE_500,
        ha="left",
        va="bottom",
    )

    fig.savefig(output_png, dpi=220, bbox_inches="tight", facecolor=fig.get_facecolor())
    if output_svg is not None:
        fig.savefig(output_svg, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)


def create_diagram(
    *,
    use_direction_grid: bool = True,
    output_png: Path = OUTPUT_PNG,
    output_svg: Path | None = OUTPUT_SVG,
) -> None:
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
            "contrasted within metadata groups, then summarized as v2 trait enrichment."
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
            "curated BugSigDB increased/decreased labels"
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
        body=(
            "Keep traits whose z-score is above zero in the focus direction and enriched "
            "relative to the comparator."
        ),
        facecolor="white",
        body_width=27,
    )
    add_box(
        ax,
        xy=(79.0, 40.8),
        width=17.8,
        height=13.8,
        title="5. Guild call",
        body=(
            "Trait sets are chosen a priori from KG-Microbe oxygen-tolerance "
            "annotations, not from focal taxa or top conditions. Each trait-enrichment score "
            "is the row-wise mean of its selected trait z-scores; selection TSVs "
            "list IDs, rationale, and taxon counts."
        ),
        facecolor="white",
        body_width=29,
        body_size=7.3,
    )

    for start_x, end_x in [(19.8, 22.0), (38.8, 41.0), (57.8, 60.0), (76.8, 79.0)]:
        add_arrow(ax, (start_x, 47.65), (end_x, 47.65))

    if use_direction_grid:
        add_large_direction_grid(ax, counts)
        add_arrow(
            ax,
            (83.0, 42.3),
            (43.6, 31.6),
            color=RED_600,
            connectionstyle="arc3,rad=0.12",
        )
        add_arrow(
            ax,
            (88.0, 42.3),
            (76.0, 24.2),
            color=BLUE_600,
            connectionstyle="arc3,rad=-0.12",
        )
    else:
        add_box(
            ax,
            xy=(8.0, 21.1),
            width=39.0,
            height=14.7,
            title="Facultative-aerobe v2 guild",
            body=(
                "Focus direction: increased profiles\n"
                "Trait score: mean of facultatively aerobic, aerobic, and "
                "microaerophilic z-scores\n"
                "Condition summary table: "
                f"{format_count(counts['facultative_condition_count'])} BugSigDB "
                "condition-label rows (not experiments)\n"
                f"target taxa audit: {format_count(counts['facultative_target_taxa_count'])}; "
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
            title="Anaerobe v2 decreased-profile trait set",
            body=(
                "Focus direction: decreased profiles across all body sites\n"
                "Trait score: mean of anaerobic and obligately anaerobic z-scores\n"
                "Condition summary table: "
                f"{format_count(counts['anaerobe_condition_count'])} BugSigDB "
                "condition-label rows (not experiments)\n"
                f"target taxa audit: {format_count(counts['anaerobe_target_taxa_count'])}; "
                f"candidate universe: {format_count(counts['anaerobe_universe_count'])} taxa"
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
        )
        add_arrow(
            ax,
            (88.0, 42.3),
            (72.0, 35.9),
            color=BLUE_600,
            connectionstyle="arc3,rad=-0.12",
        )

    add_box(
        ax,
        xy=(5.0, 5.2),
        width=26.5,
        height=11.0,
        title="Recurrence audit",
        body=(
            "Target-taxon recurrence counts how often fixed taxa reappear in the "
            "focus direction across profiles/studies; it does not define the "
            "trait-enrichment score."
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
        title="Inferred taxon support",
        body=(
            "For each enriched trait, list observed candidate-universe taxa that carry "
            "the matching KG-Microbe annotation."
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
            "Unified direction-grid TSV, condition trait enrichment, body-site diagrams, "
            "target-taxon counts, and inferred evidence TSVs."
        ),
        facecolor=SLATE_100,
        edgecolor=SLATE_200,
        title_color=SLATE_900,
        body_width=35,
    )

    if not use_direction_grid:
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

    fig.savefig(output_png, dpi=220, bbox_inches="tight", facecolor=fig.get_facecolor())
    if output_svg is not None:
        fig.savefig(output_svg, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)


def main() -> None:
    write_direction_grid_table()
    create_diagram()
    create_anaerobe_guild_definition_table()
    create_diagram(
        use_direction_grid=False,
        output_png=PRE_GRID_OUTPUT_PNG,
        output_svg=None,
    )
    print(f"Wrote: {GRID_OUTPUT_TSV}")
    print(f"Wrote: {OUTPUT_PNG}")
    print(f"Wrote: {OUTPUT_SVG}")
    print(f"Wrote: {ANAEROBE_GUILD_DEFINITION_TABLE_PNG}")
    print(f"Wrote: {ANAEROBE_GUILD_DEFINITION_TABLE_SVG}")
    print(f"Wrote: {PRE_GRID_OUTPUT_PNG}")


if __name__ == "__main__":
    main()
