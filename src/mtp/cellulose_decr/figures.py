"""Draw the ``cellulose_decr`` figures from the analysis tables (#552).

Every number, count and example on a figure is read from the TSV and JSON files
that :mod:`mtp.cellulose_decr.analysis` writes to ``data/analysis/cellulose_decr/``,
and axis limits follow the data. Writes the statistical reference plates to
``figures/`` and the presentation plates to ``figures_graphical/`` as SVG, PDF
and PNG, and a combined PDF of each set. These are regenerated, not committed
(CONTRIBUTING.md; #560); reruns on the same tables give identical bytes.

Trait names are the ``object_name`` labels joined from
``data/derived/taxa_trait_matrix_columns.tsv`` by the analysis. Captions keep
to methods/biological_interpretation.md: these are annotation counts of
curated signature taxa, not measurements of cellulose degradation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import textwrap
from collections.abc import Callable, Mapping, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Literal, cast

import matplotlib
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import to_rgba
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyBboxPatch
from scipy import stats

from mtp.cellulose_decr.analysis import (
    ALL_CONTRASTS,
    ALPHA,
    ANALYSIS_DIR,
    NON_FECAL_GI,
    ORAL,
    ORIENTED,
    OTHER,
    POOLED_STRATIFIED,
)

# Condition hues: categorical slots validated all-pairs for colour-vision
# deficiency (orange, aqua, violet). Direction: a diverging blue/red pair.
CONDITION_COLOR = {"PD": "#eb6834", "HIV": "#1baf7a", "CRC": "#6250d6"}
DOWN = "#256abf"
UP = "#e34948"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8a8984"
GRID = "#e4e3df"
FAINT = "#f4f3f0"
SHORT = ("PD", "HIV", "CRC")
FULL = {"PD": "Parkinson's disease", "HIV": "HIV infection", "CRC": "Colorectal cancer"}
PAIRS = (("PD", "HIV"), ("PD", "CRC"), ("HIV", "CRC"))
COMPARTMENT_LABEL = {
    "Feces": "Stool",
    NON_FECAL_GI: "GI mucosa / lumen",
    ORAL: "Oral",
    OTHER: "Other / multi-site",
}
PNG_DPI = 150
# Neutral highlight for p < 0.05 cells, distinct from the direction colours.
SIGNIFICANT_FILL = "#ebe4cf"

# SVG text is drawn as paths, so it renders without the font installed, and a
# fixed hash salt and no date stamps make SVG and PDF bytes repeatable.
matplotlib.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.edgecolor": MUTED,
        "axes.labelcolor": INK2,
        "axes.titlecolor": INK,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.color": INK2,
        "ytick.color": INK2,
        "svg.fonttype": "path",
        "svg.hashsalt": "cellulose_decr",
        "pdf.fonttype": 42,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    }
)


@dataclass(frozen=True)
class Results:
    summary: dict[str, Any]
    pairwise: pd.DataFrame
    sensitivity: pd.DataFrame
    compartments: pd.DataFrame
    interaction: pd.DataFrame
    catalog: pd.DataFrame
    benchmark: pd.DataFrame
    scenarios: pd.DataFrame
    scenario_interaction: pd.DataFrame


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    # Whole option names only, so the recipe's guards see every option (#617).
    parser = argparse.ArgumentParser(
        description=__doc__.splitlines()[0] if __doc__ else None, allow_abbrev=False
    )
    parser.add_argument("--analysis-dir", type=Path, default=ANALYSIS_DIR)
    parser.add_argument(
        "--code-commit",
        default="unrecorded",
        help="Revision that draws the figures, e.g. $(git rev-parse HEAD) with +dirty when the "
        "code or tables are modified. Decks link the figures' code and tables at it (#556).",
    )
    return parser.parse_args(argv)


# Every table the figures read, recorded in the figure provenance.
FIGURE_INPUTS = (
    "cellulose_decr_summary.json",
    "pairwise_jaccard_significance.tsv",
    "sensitivity_definitions_jaccard.tsv",
    "compartment_jaccard_and_direction_stats.tsv",
    "compartment_direction_interaction.tsv",
    "cellulose_decreased_taxa_catalog.tsv",
    "trait_convergence_benchmark.tsv",
    "condition_definition_sensitivity.tsv",
    "condition_definition_interaction.tsv",
)


def load(directory: Path) -> Results:
    def tsv(name: str) -> pd.DataFrame:
        if name not in FIGURE_INPUTS:
            raise ValueError(f"{name} is not declared in FIGURE_INPUTS")
        frame: pd.DataFrame = pd.read_csv(directory / name, sep="\t")
        return frame

    summary = json.loads((directory / FIGURE_INPUTS[0]).read_text(encoding="utf-8"))
    return Results(
        summary=summary,
        pairwise=tsv("pairwise_jaccard_significance.tsv"),
        sensitivity=tsv("sensitivity_definitions_jaccard.tsv"),
        compartments=tsv("compartment_jaccard_and_direction_stats.tsv"),
        interaction=tsv("compartment_direction_interaction.tsv"),
        catalog=tsv("cellulose_decreased_taxa_catalog.tsv"),
        benchmark=tsv("trait_convergence_benchmark.tsv"),
        scenarios=tsv("condition_definition_sensitivity.tsv"),
        scenario_interaction=tsv("condition_definition_interaction.tsv"),
    )


def scenario_row(
    results: Results, scenario: str, population: str, short: str, unit: str
) -> dict[str, Any]:
    table = results.scenarios
    return only(
        table[
            table["scenario"].eq(scenario)
            & table["contrast_population"].eq(population)
            & table["short_condition"].eq(short)
            & table["unit"].eq(unit)
        ]
    )


def interaction_row(
    results: Results, population: str, short: str, unit: str, compartment: str = NON_FECAL_GI
) -> dict[str, Any]:
    table = results.interaction
    return only(
        table[
            table["contrast_population"].eq(population)
            & table["short_condition"].eq(short)
            & table["unit"].eq(unit)
            & table["compartment_b"].eq(compartment)
        ]
    )


def not_estimable(row: Mapping[str, Any]) -> str:
    """Why an odds ratio of the down share (stool vs other site) cannot be estimated."""
    if int(row["b_down"]) + int(row["b_up"]) == 0:
        return "not estimable (no untied pairs outside stool)"
    if int(row["a_down"]) + int(row["a_up"]) == 0:
        return "not estimable (no untied pairs in stool)"
    if int(row["a_up"]) + int(row["b_up"]) == 0:
        return "not estimable (no ↑ in either site)"
    return "not estimable (no ↓ in either site)"


def zero_cell(row: Mapping[str, Any]) -> str:
    """Why an odds ratio of the down share (stool vs other site) is infinite or zero."""
    if int(row["a_up"]) == 0 or int(row["b_down"]) == 0:
        reason = "no ↑ in stool" if int(row["a_up"]) == 0 else "no ↓ outside stool"
        return f"OR unbounded ({reason})"
    reason = "no ↓ in stool" if int(row["a_down"]) == 0 else "no ↑ outside stool"
    return f"OR 0 ({reason})"


def records(frame: pd.DataFrame) -> list[dict[str, Any]]:
    """Rows as plain dicts, so values are read without pandas scalar types."""
    return [{str(key): value for key, value in row.items()} for row in frame.to_dict("records")]


def only(frame: pd.DataFrame) -> dict[str, Any]:
    rows = records(frame)
    if len(rows) != 1:
        raise ValueError(f"expected exactly one row, found {len(rows)}")
    return rows[0]


# --------------------------------------------------------------------------- text helpers

SUPERSCRIPT = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")
TIMES = "\N{MULTIPLICATION SIGN}"
DASH = "\N{EN DASH}"


def sci(value: float) -> str:
    """Readable p-value: two significant digits from 0.001, a x 10^k form below."""
    if value is None or math.isnan(value):
        return "n/a"
    if value >= 0.1:
        return f"{value:.2f}"
    if value >= 0.001:
        return f"{value:#.2g}"
    exponent = math.floor(math.log10(value))
    mantissa = round(value / 10**exponent, 1)
    if mantissa >= 10:
        mantissa, exponent = mantissa / 10, exponent + 1
    return f"{mantissa:.1f}{TIMES}10{str(exponent).translate(SUPERSCRIPT)}"


def mc_p(value: float, iterations: int) -> str:
    """Monte Carlo p; at the floor (no draw as extreme) say so instead of a number."""
    floor = 1 / (iterations + 1)
    if value <= floor * 1.0001:
        return f"p ≤ {sci(floor)}"
    return f"p = {sci(value)}"


def wrap(text: str, width: int) -> str:
    return "\n".join(textwrap.wrap(text, width=width))


def header(fig: Figure, title: str, subtitle: str) -> None:
    """Title wrapped to the figure width, with the subtitle placed below its last line."""
    lines = wrap(title, round(fig.get_figwidth() * 8.4))
    fig.text(0.015, 0.975, lines, ha="left", va="top", fontsize=15, fontweight="bold", color=INK)
    drop = (lines.count("\n") + 1) * 15 * 1.35 / (fig.get_figheight() * 72)
    fig.text(
        0.015,
        0.975 - drop - 0.012,
        wrap(subtitle, 175),
        ha="left",
        va="top",
        fontsize=10.5,
        color=INK2,
    )


def footnote(fig: Figure, text: str, width: int = 210) -> None:
    fig.text(0.015, 0.012, wrap(text, width), ha="left", va="bottom", fontsize=7.5, color=MUTED)


def panel_title(ax: Axes, text: str) -> None:
    ax.set_title(text, loc="left", fontsize=11, fontweight="bold", color=INK, pad=10)


def recessive_grid(ax: Axes, axis: Literal["both", "x", "y"] = "x") -> None:
    ax.grid(axis=axis, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def save(fig: Figure, directory: Path, stem: str, book: PdfPages) -> list[Path]:
    directory.mkdir(parents=True, exist_ok=True)
    written = [directory / f"{stem}{suffix}" for suffix in (".svg", ".pdf", ".png")]
    fig.savefig(written[0], metadata={"Date": None})
    fig.savefig(written[1], metadata={"CreationDate": None})
    fig.savefig(written[2], dpi=PNG_DPI)
    add_page = cast(Callable[[Figure], None], book.savefig)
    add_page(fig)
    return written


def nice_limit(value: float, step: float) -> float:
    """The smallest multiple of ``step`` at or above ``value``."""
    return step * math.ceil(value / step)


def vintage_note(results: Results) -> str:
    return (
        "PD, Parkinson's disease; HIV, HIV infection; CRC, colorectal cancer. "
        f"BugSigDB {results.summary['bugsigdb_vintage']}; human profiles; condition and body "
        "site from each experiment's own record."
    )


# --------------------------------------------------------------------------- shared pieces


def regions(catalog: pd.DataFrame) -> list[tuple[str, tuple[bool, bool, bool], list[str]]]:
    """The seven Venn regions as (label, membership in PD/HIV/CRC, taxon names)."""
    order = (
        ("PD, HIV and CRC", (True, True, True)),
        ("PD and CRC only", (True, False, True)),
        ("PD and HIV only", (True, True, False)),
        ("HIV and CRC only", (False, True, True)),
        ("PD only", (True, False, False)),
        ("CRC only", (False, False, True)),
        ("HIV only", (False, True, False)),
    )
    found = []
    for label, member in order:
        mask = (
            catalog["in_PD_feces_dec"].eq(member[0])
            & catalog["in_HIV_feces_dec"].eq(member[1])
            & catalog["in_CRC_feces_dec"].eq(member[2])
        )
        found.append((label, member, list(catalog.loc[mask, "taxon_name"])))
    return found


NULL_ROWS = (
    ("U1", "U1  uniform draw from all 189 cellulose-annotated matrix taxa"),
    ("U2", "U2  uniform draw from the taxa reported in stool, either direction"),
    ("U3", "U3  uniform draw from the taxa ever reported decreased in stool"),
    ("U4", "U4  random trait placed on all stool-decreased taxa of the three conditions"),
    ("U5", "U5  U3 taxa weighted by how often other stool studies report them decreased"),
    ("U6", "U6  whole study signatures permuted (co-reported taxa stay together)"),
    ("U6_size_matched", "U6  as above, keeping draws whose set sizes match the observed"),
)


def null_label(results: Results, key: str, text: str) -> str:
    universes = results.summary["universes"]
    text = text.replace("189", str(universes["U1_all_cellulose_taxa_in_matrix"]))
    if key == "U6_size_matched":
        iterations = results.summary["null_models"][key]["iterations"]
        text = f"{text} ({iterations:,} draws)"
    return text


def stool_rows(results: Results, population: str) -> dict[str, dict[str, Any]]:
    table = results.compartments
    chosen = table[table["contrast_population"].eq(population) & table["compartment"].eq("Feces")]
    return {str(row["short_condition"]): row for row in records(chosen)}


def compartment_row(
    results: Results, population: str, compartment: str, short: str
) -> dict[str, Any] | None:
    table = results.compartments
    chosen = table[
        table["contrast_population"].eq(population)
        & table["compartment"].eq(compartment)
        & table["short_condition"].eq(short)
    ]
    return only(chosen) if len(chosen) else None


# --------------------------------------------------------------------------- graphical plates


def graphical_convergence(results: Results) -> Figure:
    observed = results.summary["observed"]
    mean_d = observed["mean_pairwise_jaccard_distance"]
    u6 = results.summary["null_models"]["U6"]
    fig = Figure(figsize=(15, 7.2))
    header(
        fig,
        "Cellulose-annotated taxa reported decreased in stool overlap across PD, HIV and CRC",
        (
            f"Mean pairwise Jaccard distance d_J = {mean_d:.3f}; {observed['three_way_core']} of "
            f"{observed['three_way_union']} taxa appear in all three. Randomly chosen stool "
            "studies "
            "overlap about as much when their co-reported taxa are kept together "
            f"(U6, {mc_p(u6['mean_pairwise_jaccard_distance']['p'], u6['iterations'])})."
        ),
    )
    pairs = results.summary["four_highlighted_pairs"]
    contributing = {
        short: int(
            pairs[f"{index}_{short}_Feces_Cellulose_dec"]["n_supporting_studies_in_target_dir"]
        )
        for index, short in enumerate(SHORT, start=1)
    }

    # A. Network
    ax = fig.add_axes((0.015, 0.12, 0.30, 0.70))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.axis("off")
    panel_title(ax, "A. Shared taxa between each pair")
    centers = {"PD": (0.2, 0.74), "CRC": (0.8, 0.74), "HIV": (0.5, 0.2)}
    radius = 0.15
    for a, b in PAIRS:
        name = f"{a} vs. {b}"
        similarity = 1 - observed["jaccard_distance"][name]
        (x1, y1), (x2, y2) = centers[a], centers[b]
        length = math.hypot(x2 - x1, y2 - y1)
        ux, uy = (x2 - x1) / length, (y2 - y1) / length
        ax.plot(
            [x1 + ux * radius, x2 - ux * radius],
            [y1 + uy * radius, y2 - uy * radius],
            color="#b9b8b2",
            linewidth=2 + 16 * similarity,
            solid_capstyle="butt",
            zorder=1,
        )
        mid_x, mid_y = (x1 + x2) / 2, (y1 + y2) / 2
        offset = {"PD vs. CRC": (0, 0.0), "PD vs. HIV": (-0.13, 0.0), "HIV vs. CRC": (0.13, 0.0)}
        dx, dy = offset[name]
        ax.text(
            mid_x + dx,
            mid_y + dy,
            f"{observed['intersections'][name]} of {observed['unions'][name]} shared\n"
            f"d_J {observed['jaccard_distance'][name]:.2f}",
            ha="center",
            va="center",
            fontsize=9,
            color=INK,
            bbox={"boxstyle": "round,pad=0.3", "fc": "white", "ec": GRID},
            zorder=3,
        )
    for short, (x, y) in centers.items():
        ax.add_patch(
            Circle(
                (x, y),
                radius,
                facecolor=to_rgba(CONDITION_COLOR[short], 0.16),
                edgecolor=CONDITION_COLOR[short],
                linewidth=2,
                zorder=2,
            )
        )
        ax.text(x, y + 0.045, short, ha="center", va="center", fontsize=13, fontweight="bold")
        ax.text(
            x,
            y - 0.035,
            f"{observed['sizes'][short]} taxa\nfrom {contributing[short]} studies",
            ha="center",
            va="center",
            fontsize=8.5,
            color=INK2,
        )
    ax.text(
        0.5,
        -0.02,
        wrap(
            "Line width: Jaccard similarity of the pair. Studies: stool studies reporting at "
            "least one of the condition's taxa decreased.",
            62,
        ),
        ha="center",
        va="top",
        fontsize=7.5,
        color=MUTED,
    )

    # B. Venn regions
    ax = fig.add_axes((0.33, 0.12, 0.30, 0.70))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.axis("off")
    panel_title(ax, f"B. Where the {observed['three_way_union']} taxa fall")
    circles = {"PD": (0.37, 0.60), "CRC": (0.63, 0.60), "HIV": (0.50, 0.38)}
    for short, (x, y) in circles.items():
        ax.add_patch(
            Circle(
                (x, y),
                0.25,
                facecolor=CONDITION_COLOR[short] + "1f",
                edgecolor=CONDITION_COLOR[short],
                linewidth=1.8,
            )
        )
    label_at = {
        "PD": (0.13, 0.86, "left"),
        "CRC": (0.87, 0.86, "right"),
        "HIV": (0.50, 0.075, "center"),
    }
    for short, (x, y, align) in label_at.items():
        ax.text(
            x,
            y,
            f"{short} · {observed['sizes'][short]}",
            ha=align,
            va="center",
            fontsize=10.5,
            fontweight="bold",
            color=INK,
        )
    spot = {
        "PD, HIV and CRC": (0.50, 0.54),
        "PD and CRC only": (0.50, 0.74),
        "PD and HIV only": (0.36, 0.43),
        "HIV and CRC only": (0.64, 0.43),
        "PD only": (0.24, 0.66),
        "CRC only": (0.76, 0.66),
        "HIV only": (0.50, 0.22),
    }
    for label, _, names in regions(results.catalog):
        x, y = spot[label]
        core = label == "PD, HIV and CRC"
        ax.text(
            x,
            y,
            str(len(names)),
            ha="center",
            va="center",
            fontsize=17 if core else 12,
            fontweight="bold",
            color=INK,
        )
    ax.text(
        0.5,
        -0.02,
        "Circle areas are not to scale.",
        ha="center",
        va="top",
        fontsize=7.5,
        color=MUTED,
    )

    # C. Null models
    ax = fig.add_axes((0.80, 0.16, 0.185, 0.62))
    panel_title(ax, "C. Is the overlap unusual?")
    nulls = results.summary["null_models"]
    ys = np.arange(len(NULL_ROWS))[::-1]
    lowest = min(
        mean_d, *(nulls[key]["mean_pairwise_jaccard_distance"]["p2_5"] for key, _ in NULL_ROWS)
    )
    for y, (key, text) in zip(ys, NULL_ROWS, strict=True):
        record = nulls[key]["mean_pairwise_jaccard_distance"]
        ax.plot([record["p2_5"], record["p97_5"]], [y, y], color=MUTED, linewidth=3, zorder=1)
        ax.plot(record["expected"], y, "o", color=INK2, markersize=7, zorder=2)
        label = wrap(null_label(results, key, text), 34)
        ax.text(
            -0.03,
            y,
            f"{label}\n{mc_p(record['p'], nulls[key]['iterations'])}",
            ha="right",
            va="center",
            fontsize=7.4,
            color=INK2,
            transform=ax.get_yaxis_transform(),
        )
    ax.axvline(mean_d, color=INK, linewidth=1.6, zorder=3)
    ax.text(
        mean_d + 0.01,
        len(NULL_ROWS) - 0.35,
        f"observed {mean_d:.3f}",
        ha="left",
        va="center",
        fontsize=8.5,
        fontweight="bold",
    )
    ax.set_xlim(math.floor((lowest - 0.03) * 20) / 20, 1.0)
    ax.set_ylim(-0.6, len(NULL_ROWS) - 0.1)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("mean pairwise d_J under the null\n(dot: mean; bar: 95% of draws)", fontsize=8)
    recessive_grid(ax)

    footnote(
        fig,
        "Taxa: KG-Microbe matrix taxa annotated METPO:2000103|GO:0030245 (cellulose catabolic "
        "process; any rank, annotations rolled up from descendants) in at least one decreased "
        "stool signature of the condition. Overlap counts reported taxonomic labels, not "
        "consistent loss: a taxon counts once it is reported decreased in any study. Monte Carlo "
        "p from "
        f"{u6['iterations']:,} draws. {vintage_note(results)}",
    )
    return fig


def direct_core_text(results: Results) -> str:
    direct = int(results.summary["core_taxa_with_direct_trait_edge"])
    core = int(results.summary["observed"]["three_way_core"])
    if direct == 0:
        return f"none of the {core} core taxa has its own annotation edge"
    if direct == core:
        return f"all {core} core taxa have their own annotation edge"
    verb = "has its" if direct == 1 else "have their"
    return f"{direct} of the {core} core taxa {verb} own annotation edge"


# Ranks at or below genus form the first group of the core figures; the rest the second.
LOW_RANKS = ("genus", "species", "subspecies", "strain")
FAMILY_TO_PHYLUM = {"family", "order", "class", "phylum"}


def core_taxa(catalog: pd.DataFrame) -> pd.DataFrame:
    core = catalog[catalog["is_3way_core"].astype(bool)].copy()
    core["total_dec"] = core[[f"{s}_dec_studies" for s in SHORT]].sum(axis=1)
    core["group"] = np.where(core["taxonomic_rank"].isin(LOW_RANKS), 0, 1)
    return core.sort_values(["group", "total_dec", "taxon_name"], ascending=[True, False, True])


def rank_groups(core: pd.DataFrame) -> tuple[str, str]:
    """The two rank groups of the core, counted and named by the ranks they hold."""
    low = set(core.loc[core["group"].eq(0), "taxonomic_rank"])
    high = set(core.loc[core["group"].eq(1), "taxonomic_rank"])
    n_low, n_high = int(core["group"].eq(0).sum()), int(core["group"].eq(1).sum())
    if low <= {"genus"}:
        low_text = f"{n_low} genera"
    elif low <= {"genus", "species"}:
        low_text = f"{n_low} genus- or species-rank taxa"
    else:
        low_text = f"{n_low} genus-or-lower-rank taxa"
    nodes = "family-to-phylum" if high <= FAMILY_TO_PHYLUM else "higher-rank"
    return low_text, f"{n_high} {nodes} nodes"


SENSITIVITY_LABELS = {
    "1_feces_dec_all_ranks": "Primary: all ranks",
    "3_feces_dec_genus_species": "Genus and species ranks only",
    "2_feces_dec_z_gt_0": "Profiles with cellulose z > 0 only",
    "7_feces_dec_disease_vs_healthy_oriented": "Disease-vs-healthy contrasts only",
    "8_feces_dec_disease_vs_healthy_strong": "Disease-vs-healthy, strong evidence",
    "10_feces_dec_synonyms_merged": "Synonym condition labels counted",
    "9_feces_dec_without_site_miscurations": "Without the rectal-mucosa study",
    "5_all_body_sites_dec": "All body sites, not only stool",
    "4_feces_net_decreased": "Net decreased (more dec than inc profiles)",
    "4b_feces_net_decreased_by_studies": "Net decreased (more dec than inc studies)",
    "6_baseline_all_feces_dec_taxa": "Reference, not a null: all stool-decreased taxa",
}


def reversed_examples(core: pd.DataFrame, limit: int = 2) -> str:
    """Core taxa reported increased in more studies than decreased, largest margins first."""
    examples = []
    for row in records(core):
        for short in SHORT:
            margin = int(row[f"{short}_inc_studies"]) - int(row[f"{short}_dec_studies"])
            if margin > 0:
                examples.append((-margin, str(row["taxon_name"]), short))
    if not examples:
        return ""
    shown = ", ".join(f"{name} in {short}" for _, name, short in sorted(examples)[:limit])
    return f" (for example {shown})" if len(examples) > limit else f" ({shown})"


def graphical_core(results: Results) -> Figure:
    catalog = results.catalog
    core = core_taxa(catalog)
    n_direct = int(core["has_direct_trait_edge"].astype(bool).sum())
    net = core[core["is_strictly_net_decreased_all3"].astype(bool)]
    sens = results.sensitivity
    mean_rows = {
        str(row["definition_slug"]): row
        for row in records(sens[sens["comparison"].str.startswith("Mean")])
    }
    gs_mean = float(mean_rows["3_feces_dec_genus_species"]["jaccard_distance"])
    fig = Figure(figsize=(15, 7.6))
    if n_direct == 0:
        direct_text = (
            "none has its own annotation edge; each inherits it from annotated descendants"
        )
    else:
        direct_text = (
            f"{n_direct} {'has its' if n_direct == 1 else 'have their'} own annotation edge"
        )
    low_text, high_text = rank_groups(core)
    header(
        fig,
        (
            f"The {len(core)} taxa reported decreased in stool in all three conditions: "
            f"{low_text} and {high_text}"
        ),
        (
            f"Of the {len(core)}, {direct_text}. {len(net)} had more decreased than increased "
            "stool studies in every condition. "
            f"Keeping only genus and species ranks gives mean d_J = {gs_mean:.3f}, close to the "
            f"primary {results.summary['observed']['mean_pairwise_jaccard_distance']:.3f}."
        ),
    )

    # A. Sensitivity of the mean distance to the set definition
    ax = fig.add_axes((0.215, 0.14, 0.26, 0.66))
    panel_title(ax, "A. Mean pairwise d_J under each set definition")
    order = list(SENSITIVITY_LABELS)
    values = [float(mean_rows[slug]["jaccard_distance"]) for slug in order]
    ys = np.arange(len(order))[::-1]
    for y, slug, value in zip(ys, order, values, strict=True):
        primary = slug.startswith("1_")
        baseline = slug.startswith("6_")
        if baseline:
            ax.barh(y, value, height=0.62, color="white", edgecolor=MUTED, hatch="///")
        else:
            ax.barh(y, value, height=0.62, color=INK if primary else "#9a9993")
        core_n = int(mean_rows[slug]["intersection"])
        union_n = int(mean_rows[slug]["union"])
        ax.text(
            value + 0.01,
            y,
            f"{value:.3f}   {core_n}/{union_n} in all three",
            va="center",
            fontsize=8,
            color=INK,
        )
    ax.set_yticks(ys)
    ax.set_yticklabels([SENSITIVITY_LABELS[slug] for slug in order], fontsize=8.5)
    ax.set_xlim(0, 1.15)
    ax.set_xlabel("mean pairwise Jaccard distance (lower = more overlap)")
    recessive_grid(ax)

    # B. Dot table of the core taxa
    ax = fig.add_axes((0.62, 0.17, 0.36, 0.65))
    panel_title(ax, "B. Stool studies reporting each core taxon ↓ or ↑")
    ax.set_xlim(-0.5, 4.2)
    rows = list(core.itertuples(index=False))
    n_rows = len(rows) + 1
    ax.set_ylim(-0.2, n_rows)
    ax.axis("off")
    columns = {"PD": 0.4, "HIV": 1.3, "CRC": 2.2}
    for short, x in columns.items():
        ax.text(x, n_rows - 0.3, short, ha="center", fontsize=10, fontweight="bold")
    ax.text(3.25, n_rows - 0.3, "net ↓ in all 3", ha="center", fontsize=8.5, fontweight="bold")
    y = n_rows - 1.2
    previous_group = None
    for row in rows:
        if previous_group is not None and row.group != previous_group:
            ax.axhline(y + 0.55, xmin=0.0, xmax=1.0, color=GRID, linewidth=1)
        previous_group = row.group
        rank = str(row.taxonomic_rank)
        ax.text(
            -0.5,
            y,
            f"{row.taxon_name}",
            ha="right",
            va="center",
            fontsize=9,
            style="italic" if rank in ("genus", "species") else "normal",
        )
        ax.text(-0.45, y, rank, ha="left", va="center", fontsize=7, color=MUTED)
        for short, x in columns.items():
            dec = int(getattr(row, f"{short}_dec_studies"))
            inc = int(getattr(row, f"{short}_inc_studies"))
            color = DOWN if dec > inc else (UP if inc > dec else MUTED)
            ax.plot(x - 0.3, y, "o", color=color, markersize=6)
            ax.text(x - 0.18, y, f"{dec}↓ {inc}↑", va="center", fontsize=8.5, color=INK)
        if row.is_strictly_net_decreased_all3:
            ax.text(3.25, y, "✓", ha="center", va="center", fontsize=11, color=DOWN)
        y -= 1
    fig.text(
        0.62,
        0.155,
        wrap(
            "Dot: blue when ↓ studies outnumber ↑, red when ↑ outnumber ↓, grey when equal. "
            "Upper block: genus rank; lower block: family, order, class and phylum. A taxon "
            "bears the annotation when it or a descendant has the KG edge. Net ↓ = more "
            "decreased than increased studies in every condition.",
            92,
        ),
        ha="left",
        va="top",
        fontsize=7.5,
        color=MUTED,
    )
    footnote(
        fig,
        "A taxon is in a condition's set when at least one decreased stool signature reports it; "
        "the counts in B show that some core taxa are also reported increased, sometimes more "
        f"often{reversed_examples(core)}. Hatched bar in A: every stool-decreased taxon whatever "
        "its traits, for scale only; it is not a null model (U4 in the null-model figure tests "
        f"the cellulose set against it). {vintage_note(results)}",
    )
    return fig


def untied(row: Mapping[str, Any], unit: str) -> tuple[int, int]:
    if unit == "experiment":
        return int(row["paired_exp_down_dec_gt_inc"]), int(row["paired_exp_up_inc_gt_dec"])
    return int(row["study_down"]), int(row["study_up"])


def graphical_compartments(results: Results) -> Figure:
    fig = Figure(figsize=(15, 8.8))
    pooled = compartment_row(results, "all_contrasts", "Feces", "Pooled")
    if pooled is None:
        raise ValueError("no pooled stool row in the compartment table")
    crc_gi = compartment_row(results, "all_contrasts", NON_FECAL_GI, "CRC")
    other_sites = [
        (row["short_condition"], row["compartment"], sum(untied(row, "experiment")))
        for row in records(
            results.compartments[
                results.compartments["contrast_population"].eq("all_contrasts")
                & ~results.compartments["compartment"].eq("Feces")
                & results.compartments["short_condition"].isin(SHORT)
            ]
        )
    ]
    largest = max(other_sites, key=lambda item: item[2])
    lean = ""
    if crc_gi is not None and float(crc_gi["share_down_untied"]) < 0.5:
        gi_down, gi_up = untied(crc_gi, "experiment")
        lean = f" CRC GI mucosa / lumen leans the other way ({gi_down}↓ {gi_up}↑)."
    exp_down, exp_up = untied(pooled, "experiment")
    study_down, study_up = untied(pooled, "study")
    header(
        fig,
        "The cellulose annotation share is larger in the decreased signature, most clearly in "
        "stool",
        (
            f"Across PD, HIV and CRC stool, {exp_down} of {exp_down + exp_up} untied paired "
            f"experiments (p = {sci(pooled['binom_p_two_sided'])}) and {study_down} of "
            f"{study_down + study_up} untied studies "
            f"(p = {sci(pooled['study_binom_p_two_sided'])}) "
            f"had the larger cellulose annotation share in the decreased signature.{lean} Outside "
            f"stool the largest single-condition stratum is {largest[0]} "
            f"{COMPARTMENT_LABEL[str(largest[1])]} with {largest[2]} untied pairs."
        ),
    )

    # A. Butterfly of experiment counts
    ax = fig.add_axes((0.165, 0.21, 0.30, 0.56))
    panel_title(ax, "A. Paired experiments with the larger share in ↓ or ↑ signature")
    strata: list[tuple[str, str, dict[str, Any]]] = []
    for compartment in ("Feces", NON_FECAL_GI, ORAL, OTHER):
        names = (
            ("PD", "HIV", "CRC", "Pooled")
            if compartment in ("Feces", NON_FECAL_GI)
            else ("Pooled",)
        )
        for short in names:
            found = compartment_row(results, "all_contrasts", compartment, short)
            if found is not None:
                strata.append((compartment, short, found))
    ys = np.arange(len(strata))[::-1].astype(float)
    gap = 0.0
    positions = []
    previous = None
    for y, (compartment, _, _) in zip(ys, strata, strict=True):
        if previous is not None and compartment != previous:
            gap += 0.7
        positions.append(y - gap)
        previous = compartment
    widest = max(max(untied(row, "experiment")) for _, _, row in strata)
    limit = nice_limit(widest * 1.12, 10)
    for y, (compartment, short, row) in zip(positions, strata, strict=True):
        down, up = untied(row, "experiment")
        ax.barh(y, -down, height=0.6, color=DOWN)
        ax.barh(y, up, height=0.6, color=UP)
        ax.text(-down - limit / 60, y, str(down), ha="right", va="center", fontsize=8, color=INK)
        ax.text(up + limit / 60, y, str(up), ha="left", va="center", fontsize=8, color=INK)
        label = "All three" if short == "Pooled" else short
        ax.text(
            -limit * 1.03,
            y,
            f"{COMPARTMENT_LABEL[compartment]} · {label}",
            ha="right",
            va="center",
            fontsize=8.5,
            fontweight="bold" if short == "Pooled" else "normal",
        )
        s_down, s_up = untied(row, "study")
        lines = []
        if float(row["min_attainable_binom_p_two_sided"]) > 0.05:
            lines.append(f"studies {s_down}↓ {s_up}↑ · too few pairs for p < 0.05")
        else:
            lines.append(f"experiments p = {sci(row['binom_p_two_sided'])}")
            study_p = (
                "too few studies for p < 0.05"
                if s_down + s_up == 0 or 2 * 0.5 ** (s_down + s_up) > 0.05
                else f"p = {sci(row['study_binom_p_two_sided'])}"
            )
            lines.append(f"studies {s_down}↓ {s_up}↑ · {study_p}")
        if compartment == "Feces":
            healthy = compartment_row(results, "disease_vs_healthy_oriented", "Feces", short)
            if healthy is not None:
                h_down, h_up = untied(healthy, "study")
                lines.append(
                    f"disease vs healthy studies {h_down}↓ {h_up}↑, "
                    f"p = {sci(healthy['study_binom_p_two_sided'])}"
                )
        ax.text(
            limit * 1.03,
            y,
            "\n".join(lines),
            ha="left",
            va="center",
            fontsize=7.2,
            color=INK2,
            linespacing=1.15,
        )
    ax.axvline(0, color=INK, linewidth=0.8)
    ax.set_xlim(-limit, limit)
    ax.set_ylim(min(positions) - 0.7, max(positions) + 0.7)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    step = 50 if limit >= 100 else (20 if limit >= 40 else 10)
    reach = step * math.floor(limit / step)
    ticks = [float(tick) for tick in np.arange(-reach, reach + step, step)]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{abs(tick):.0f}" for tick in ticks])
    ax.set_xlabel("untied paired experiments (← larger in decreased · larger in increased →)")
    recessive_grid(ax)

    # B. Stool versus GI mucosa/lumen
    ax = fig.add_axes((0.845, 0.50, 0.14, 0.30))
    panel_title(ax, "B. Stool vs GI\nmucosa / lumen")
    chosen = [
        ("CRC", "experiment", "all_contrasts", "CRC · all · experiments"),
        ("CRC", "study", "all_contrasts", "CRC · all · studies"),
        ("CRC", "experiment", "disease_vs_healthy_oriented", "CRC · disease vs healthy · exp."),
        ("CRC", "study", "disease_vs_healthy_oriented", "CRC · disease vs healthy · studies"),
        (POOLED_STRATIFIED, "experiment", "all_contrasts", "All three, stratified · exp."),
        (POOLED_STRATIFIED, "study", "all_contrasts", "All three, stratified · studies"),
    ]
    picked = [
        interaction_row(results, population, short, unit) for short, unit, population, _ in chosen
    ]
    finite_high = [float(row["conditional_or_ci95_high"]) for row in picked]
    finite_high = [value for value in finite_high if math.isfinite(value)]
    finite_low = [float(row["conditional_or_ci95_low"]) for row in picked]
    x_max = 10 ** math.ceil(math.log10(max([*finite_high, 10.0]) * 1.2))
    x_min = min(0.5, 10 ** math.floor(math.log10(max(min([*finite_low, 1.0]), 1e-3))))
    ys = np.arange(len(chosen))[::-1]
    for y, (_, _, _, label), row in zip(ys, chosen, picked, strict=True):
        low = float(row["conditional_or_ci95_low"])
        high = float(row["conditional_or_ci95_high"])
        estimate = float(row["conditional_or"])
        note = f"p = {sci(row['fisher_p_two_sided'])}"
        if math.isnan(estimate):
            note = f"{note} · {not_estimable(row)}"
        else:
            ax.plot([max(low, x_min), min(high, x_max)], [y, y], color=MUTED, linewidth=2.5)
            if high > x_max:
                ax.plot(x_max, y, ">", color=MUTED, markersize=5, clip_on=False)
            if math.isinf(estimate) or estimate == 0:
                edge = x_max if math.isinf(estimate) else x_min
                marker = ">" if math.isinf(estimate) else "<"
                ax.plot(edge, y, marker, color=INK, markersize=6, clip_on=False)
                note = f"{note} · {zero_cell(row)}"
            else:
                ax.plot(estimate, y, "o", color=INK, markersize=6)
        ax.text(
            -0.05,
            y,
            f"{label}\n{note}",
            ha="right",
            va="center",
            fontsize=7.4,
            transform=ax.get_yaxis_transform(),
        )
    ax.axvline(1, color=INK2, linewidth=0.8, linestyle="--")
    ax.set_xscale("log")
    ax.set_xlim(x_min, x_max)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel(
        "odds ratio of the ↓ share, stool vs GI\n(log scale; exact 95% interval)", fontsize=8
    )
    recessive_grid(ax)

    # C. Sampling effort
    ax = fig.add_axes((0.75, 0.07, 0.22, 0.22))
    panel_title(ax, "C. Decreased-signature studies by site")
    ax.axis("off")
    compartments = ("Feces", NON_FECAL_GI, ORAL, OTHER)
    for j, short in enumerate(SHORT):
        ax.text(0.45 + j * 0.2, 0.95, short, ha="center", fontsize=9, fontweight="bold")
    for i, compartment in enumerate(compartments):
        y = 0.75 - i * 0.2
        ax.text(0.33, y, COMPARTMENT_LABEL[compartment], ha="right", va="center", fontsize=8.5)
        for j, short in enumerate(SHORT):
            found = compartment_row(results, "all_contrasts", compartment, short)
            value = int(found["n_studies_dec"]) if found is not None else 0
            ax.text(
                0.45 + j * 0.2,
                y,
                str(value),
                ha="center",
                va="center",
                fontsize=9.5,
                color=INK if value >= 10 else MUTED,
            )
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    exact_hiv = scenario_row(results, "exact_labels", "disease_vs_healthy_oriented", "HIV", "study")
    merged_hiv = scenario_row(
        results, "synonyms_merged", "disease_vs_healthy_oriented", "HIV", "study"
    )
    art_hiv = scenario_row(
        results, "synonyms_hiv_minus_vs_plus_and_art", "disease_vs_healthy_oriented", "HIV", "study"
    )
    gi_all = compartment_row(results, "all_contrasts", NON_FECAL_GI, "CRC")
    gi_healthy = compartment_row(results, "disease_vs_healthy_oriented", NON_FECAL_GI, "CRC")
    gi_note = ""
    if gi_all is not None and gi_healthy is not None:
        table = results.scenario_interaction
        without = only(
            table[
                table["scenario"].eq("crc_within_patient_excluded")
                & table["contrast_population"].eq("disease_vs_healthy_oriented")
                & table["short_condition"].eq("CRC")
                & table["unit"].eq("study")
            ]
        )
        gi_note = (
            f"only {int(gi_healthy['paired_experiments'])} of {int(gi_all['paired_experiments'])} "
            "paired CRC mucosa/lumen experiments are classified disease vs healthy, and one of "
            "those compares tumour with adjacent tissue (without it the disease-vs-healthy "
            f"study row of B gives p = {sci(float(without['fisher_p_two_sided']))}), "
        )
    footnote(
        fig,
        "Paired experiment: one BugSigDB experiment with both an increased and a decreased "
        "signature; the share compared is cellulose annotations over all retained non-MediaDive "
        "annotations. Ties excluded; p-values are two-sided sign tests, nominal. Study rows take "
        "the sign of each study's mean paired difference. HIV depends on which experiments count: "
        f"disease vs healthy studies {int(exact_hiv['down'])}↓ {int(exact_hiv['up'])}↑ "
        f"(p = {sci(exact_hiv['binom_p_two_sided'])}) with exact labels, "
        f"{int(merged_hiv['down'])}↓ {int(merged_hiv['up'])}↑ "
        f"(p = {sci(merged_hiv['binom_p_two_sided'])}) with synonym labels, "
        f"{int(art_hiv['down'])}↓ {int(art_hiv['up'])}↑ (p = {sci(art_hiv['binom_p_two_sided'])}) "
        "when HIV-negative-vs-positive and treated contrasts are added "
        "(condition_definition_sensitivity.tsv). The B tests treat units as independent, "
        f"{gi_note}and the sites differ in sampling effort (C), so a difference between "
        f"compartments is not evidence of anatomical specificity. {vintage_note(results)}",
        width=178,
    )
    return fig


ORDINAL = {1: "lowest", 2: "second-lowest", 3: "third-lowest", 4: "fourth-lowest"}


def trait_display(row: Mapping[str, Any]) -> str:
    """Canonical label, with its predicate when that is not plain capability."""
    trait_id = str(row["trait_id"])
    label = str(row["trait_label"])
    predicate = row["trait_predicate"]
    if isinstance(predicate, str) and predicate and predicate != "capable of":
        label = f"{predicate}: {label}"
    short_id = trait_id.split("|")[-1]
    if short_id.startswith("pathways:"):
        short_id = "pathway"
    return f"{label} ({short_id})"


def graphical_benchmark(results: Results) -> Figure:
    bench = results.benchmark
    biological = bench[~bench["is_growth_medium"].astype(bool)].reset_index(drop=True)
    rows = records(biological)
    n_traits = len(rows)
    cell = next(row for row in rows if row["trait_id"] == results.summary["trait_id"])
    first = rows[0]
    rank = int(cell["rank_excluding_growth_media"])
    fig = Figure(figsize=(15, 6.8))
    header(
        fig,
        (
            f"Among {n_traits} non-MediaDive traits, {cell['trait_label']} has the "
            f"{ORDINAL.get(rank, f'rank-{rank}')} cross-condition Jaccard distance"
        ),
        (
            f"Traits with at least 10 stool-decreased taxa in each condition. Rank {rank} of "
            f"{n_traits}, after {first['trait_label']} ({first['trait_id']}). Ranks describe this "
            "corpus; they are not a test, and traits share taxa."
        ),
    )
    ax = fig.add_axes((0.05, 0.16, 0.33, 0.62))
    panel_title(ax, f"A. Mean pairwise d_J of all {n_traits} traits")
    values = [float(row["mean_pairwise_jaccard_distance"]) for row in rows]
    start = math.floor(min(values) * 50) / 50
    stop = math.ceil(max(values) * 50) / 50
    bins = [round(start + 0.02 * step, 2) for step in range(round((stop - start) / 0.02) + 1)]
    ax.hist(values, bins=bins, color="#c9c8c2", edgecolor="white")
    # Headroom so the labels at the top clear the tallest bar.
    ax.set_ylim(0, ax.get_ylim()[1] * 1.15)
    median = float(np.median(values))
    ax.axvline(median, color=INK2, linestyle="--", linewidth=1)
    top_y = ax.get_ylim()[1]
    ax.text(
        0.98,
        0.97,
        f"dashed: median {median:.3f}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8.5,
        color=INK2,
    )
    for row, height in ((first, 0.55), (cell, 0.97)):
        value = float(row["mean_pairwise_jaccard_distance"])
        ax.axvline(value, color=INK, linewidth=1.6)
        ax.text(
            value + 0.006,
            top_y * height,
            f"rank {int(row['rank_excluding_growth_media'])}\n"
            f"{wrap(str(row['trait_label']), 20)}\n{value:.3f}",
            ha="left",
            va="top",
            fontsize=8.5,
        )
    ax.set_xlabel("mean pairwise Jaccard distance across PD, HIV and CRC")
    ax.set_ylabel("traits")
    recessive_grid(ax, "y")

    ax = fig.add_axes((0.62, 0.16, 0.36, 0.62))
    panel_title(ax, "B. The ten lowest distances")
    top = rows[:10]
    ys = np.arange(len(top))[::-1]
    for y, row in zip(ys, top, strict=True):
        value = float(row["mean_pairwise_jaccard_distance"])
        highlight = row["trait_id"] == results.summary["trait_id"]
        color = INK if highlight else "#9a9993"
        ax.plot([0, value], [y, y], color=color, linewidth=2.5)
        ax.plot(value, y, "o", color=color, markersize=7)
        ax.text(
            value + 0.012,
            y,
            f"{value:.3f} · {int(row['three_way_intersection'])}/{int(row['three_way_union'])} "
            "in all three",
            va="center",
            fontsize=8,
        )
    ax.set_yticks(ys)
    ax.set_yticklabels([trait_display(row) for row in top], fontsize=8.5)
    for label, row in zip(ax.get_yticklabels(), top, strict=True):
        if row["trait_id"] == results.summary["trait_id"]:
            label.set_fontweight("bold")
    ax.set_xlim(
        0, nice_limit(max(float(row["mean_pairwise_jaccard_distance"]) for row in top) + 0.2, 0.1)
    )
    ax.set_xlabel("mean pairwise Jaccard distance (lower = more overlap)")
    recessive_grid(ax)
    footnote(
        fig,
        "Trait names are the object_name labels in data/derived/taxa_trait_matrix_columns.tsv; "
        "a predicate other than 'capable of' is shown before the name. Excluded: MediaDive "
        "growth-medium columns (METPO:2000517|mediadive.medium:*). "
        f"{vintage_note(results)}",
    )
    return fig


# --------------------------------------------------------------------------- detailed plates


TABLE_NULLS = (
    ("u1", "U1", "U1\nmatrix\nuniform"),
    ("u2", "U2", "U2\nstool\nuniform"),
    ("u3", "U3", "U3\nstool dec\nuniform"),
    ("u4_shared_trait", "U4", "U4\nrandom\ntrait"),
    ("u5_emp_weighted", "U5", "U5\nprevalence\nweighted"),
    ("u6_signature_perm", "U6", "U6\nsignature\npermutation"),
    ("u6_size_matched", "U6_size_matched", "U6\nsize-\nmatched"),
)


def detailed_nulls(results: Results) -> Figure:
    pairwise = results.pairwise
    nulls = results.summary["null_models"]
    observed = results.summary["observed"]
    fig = Figure(figsize=(15, 9.8))

    def mean_p(key: str) -> str:
        return mc_p(nulls[key]["mean_pairwise_jaccard_distance"]["p"], nulls[key]["iterations"])

    def core_p(key: str) -> str:
        return mc_p(nulls[key]["three_way_core"]["p"], nulls[key]["iterations"])

    independent = max(
        float(nulls[key]["mean_pairwise_jaccard_distance"]["p"]) for key in ("U1", "U2", "U3", "U4")
    )
    header(
        fig,
        "Overlap of stool-decreased cellulose-annotated taxa against six null models",
        (
            f"Against nulls that draw taxa independently (U1{DASH}U4) the mean overlap is "
            f"extreme (every p ≤ {sci(independent)}). Weighting taxa by how often other stool "
            f"studies report them (U5) weakens it (mean {mean_p('U5')}; core {core_p('U5')}); "
            "permuting whole study signatures (U6) removes it "
            f"(mean {mean_p('U6')}; core {core_p('U6')})."
        ),
    )
    # A. Table graphic: expected d_J and p for each comparison and null
    ax = fig.add_axes((0.125, 0.50, 0.52, 0.34))
    panel_title(ax, "A. Expected Jaccard distance and p under each null")
    ax.axis("off")
    comparisons = [f"{a} vs. {b}" for a, b in PAIRS] + ["Mean of 3 pairs"]
    columns = ["observed", *[key for key, _, _ in TABLE_NULLS]]
    n_cols = len(columns)
    for j, (_, _, label) in enumerate([("observed", "", "observed\nd_J"), *TABLE_NULLS]):
        ax.text((j + 0.5) / n_cols, 0.97, label, ha="center", va="top", fontsize=8.5)
    for i, name in enumerate(comparisons):
        y = 0.62 - i * 0.18
        ax.text(-0.01, y, name, ha="right", va="center", fontsize=9, fontweight="bold")
        if name.startswith("Mean"):
            row = only(pairwise[pairwise["endpoint"].eq("mean_pairwise_jaccard_distance")])
        else:
            row = only(pairwise[pairwise["comparison"].eq(name)])
        null_keys = {column: key for column, key, _ in TABLE_NULLS}
        for j, column in enumerate(columns):
            x = (j + 0.5) / n_cols
            if column == "observed":
                ax.text(
                    x,
                    y,
                    f"{row['jaccard_distance']:.3f}",
                    ha="center",
                    va="center",
                    fontsize=10,
                    fontweight="bold",
                )
                continue
            key = null_keys[column]
            expected = float(row.get(f"{column}_expected_jaccard_distance", float("nan")))
            if column in ("u1", "u2", "u3"):
                exact = not name.startswith("Mean")
                p = float(row[f"{column}_hypergeom_p" if exact else f"{column}_mc_p"])
                p_text = f"exact {sci(p)}" if exact else mc_p(p, nulls[key]["iterations"])
            else:
                p = float(row[f"{column}_p"])
                p_text = mc_p(p, nulls[key]["iterations"])
            strong = not math.isnan(p) and p < 0.05
            ax.add_patch(
                FancyBboxPatch(
                    (x - 0.45 / n_cols, y - 0.075),
                    0.9 / n_cols,
                    0.15,
                    boxstyle="round,pad=0.002",
                    facecolor=SIGNIFICANT_FILL if strong else FAINT,
                    edgecolor="none",
                    transform=ax.transAxes,
                )
            )
            ax.text(x, y + 0.025, f"{expected:.3f}", ha="center", va="center", fontsize=9)
            ax.text(
                x,
                y - 0.035,
                p_text.replace("p = ", "p "),
                ha="center",
                va="center",
                fontsize=7.4,
                color=INK2,
            )
    fig.text(
        0.125,
        0.475,
        wrap(
            "Observed values are in the first column. Shaded sand: p < 0.05 (nominal). "
            f"U1{DASH}U3 pair p-values are exact hypergeometric tails; all other p are Monte "
            "Carlo. U4 places a random trait on the taxa reported decreased in stool in the "
            "three conditions; U5 draws from the U3 taxa weighted by "
            "their reporting frequency in other stool studies; U6 permutes whole study signatures.",
            150,
        ),
        ha="left",
        va="top",
        fontsize=7.5,
        color=MUTED,
    )

    # B. 3-way core against the nulls
    ax = fig.add_axes((0.80, 0.50, 0.18, 0.32))
    panel_title(ax, "B. Taxa in all three sets")
    keys = [
        ("U1", "U1"),
        ("U2", "U2"),
        ("U3", "U3"),
        ("U4", "U4"),
        ("U5", "U5"),
        ("U6", "U6"),
        ("U6_size_matched", "U6 size-matched"),
    ]
    ys = np.arange(len(keys))[::-1]
    for y, (key, label) in zip(ys, keys, strict=True):
        record = nulls[key]["three_way_core"]
        ax.plot([record["p2_5"], record["p97_5"]], [y, y], color=MUTED, linewidth=3)
        ax.plot(record["expected"], y, "o", color=INK2, markersize=6)
        ax.text(
            -0.04,
            y,
            f"{label}\n{mc_p(record['p'], nulls[key]['iterations'])}",
            ha="right",
            va="center",
            fontsize=7.8,
            transform=ax.get_yaxis_transform(),
        )
    core = observed["three_way_core"]
    ax.axvline(core, color=INK, linewidth=1.6)
    ax.text(
        core + 0.5,
        len(keys) - 0.15,
        f"observed {core}",
        ha="left",
        va="center",
        fontsize=8.5,
        fontweight="bold",
    )
    widest = max(core, *(nulls[key]["three_way_core"]["p97_5"] for key, _ in keys))
    ax.set_xlim(0, nice_limit(widest * 1.25, 5))
    ax.set_ylim(-0.7, len(keys) + 0.2)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("taxa in all three (dot: null mean; bar: 95% of draws)", fontsize=8)
    recessive_grid(ax)

    # C. Partition with taxon names
    ax = fig.add_axes((0.05, 0.07, 0.93, 0.31))
    panel_title(
        ax,
        f"C. The {observed['three_way_union']} taxa by the conditions that report them "
        "decreased in stool",
    )
    ax.axis("off")
    found = regions(results.catalog)
    width = 1 / len(found)
    for i, (label, _, names) in enumerate(found):
        x = i * width
        ax.text(
            x + width / 2,
            0.98,
            f"{label}\n{len(names)}",
            ha="center",
            va="top",
            fontsize=9,
            fontweight="bold",
        )
        ax.text(
            x + 0.008,
            0.80,
            "\n".join(names),
            ha="left",
            va="top",
            fontsize=6.9,
            color=INK2,
            linespacing=1.25,
        )
    footnote(
        fig,
        "Names are KG-Microbe matrix taxon names. A taxon bears the annotation when it or a "
        f"descendant has the KG edge; {direct_core_text(results)}. {vintage_note(results)}",
    )
    return fig


def share_interval(down: int, up: int) -> tuple[float, float, float]:
    total = down + up
    if total == 0:
        return float("nan"), float("nan"), float("nan")
    interval = stats.binomtest(down, total).proportion_ci(method="exact")
    return down / total, float(interval.low), float(interval.high)


def detailed_compartments(results: Results) -> Figure:
    fig = Figure(figsize=(15, 10.2))
    header(
        fig,
        "Direction of the cellulose annotation shift by body site, inference unit and contrast "
        "type",
        (
            "Share of untied units whose decreased signature had the larger cellulose annotation "
            "share, with exact 95% intervals; 0.5 means no preference."
        ),
    )
    ax = fig.add_axes((0.22, 0.42, 0.40, 0.45))
    panel_title(ax, "A. Share of untied units whose decreased signature has the larger share")
    rows = []
    for compartment in ("Feces", NON_FECAL_GI, ORAL, OTHER):
        for short in ("PD", "HIV", "CRC", "Pooled"):
            rows.append((compartment, short))
    ys = np.arange(len(rows))[::-1].astype(float)
    # Series differ by marker shape and fill, not by the direction colours.
    markers = (
        ("all_contrasts", "experiment", "o", INK, INK, "all contrasts · experiments"),
        ("all_contrasts", "study", "s", "#9a9993", "#9a9993", "all contrasts · studies"),
        (
            "disease_vs_healthy_oriented",
            "experiment",
            "D",
            INK,
            "white",
            "disease vs healthy · experiments",
        ),
        (
            "disease_vs_healthy_oriented",
            "study",
            "^",
            "#9a9993",
            "white",
            "disease vs healthy · studies",
        ),
    )
    for y, (compartment, short) in zip(ys, rows, strict=True):
        label = "All three" if short == "Pooled" else short
        ax.text(
            -0.04,
            y,
            f"{COMPARTMENT_LABEL[compartment]} · {label}",
            ha="right",
            va="center",
            fontsize=8.3,
            transform=ax.get_yaxis_transform(),
        )
        for offset, (population, unit, marker, color, fill, _) in zip(
            (0.27, 0.09, -0.09, -0.27), markers, strict=True
        ):
            row = compartment_row(results, population, compartment, short)
            if row is None:
                continue
            down, up = untied(row, unit)
            share, low, high = share_interval(down, up)
            if math.isnan(share):
                continue
            ax.plot(
                [low, high],
                [y + offset, y + offset],
                color=color,
                linewidth=1.2,
                alpha=0.8,
                clip_on=False,
            )
            ax.plot(
                share,
                y + offset,
                marker,
                color=color,
                markerfacecolor=fill,
                markersize=4.5,
                clip_on=False,
            )
            if population == "all_contrasts" and unit == "experiment":
                ax.text(
                    1.03,
                    y,
                    f"{down}↓ {up}↑",
                    ha="left",
                    va="center",
                    fontsize=7.6,
                    color=INK2,
                    transform=ax.get_yaxis_transform(),
                )
    ax.axvline(0.5, color=INK2, linestyle="--", linewidth=0.8)
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.7, len(rows) - 0.3)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel(
        "share of untied units in which the decreased signature has the larger "
        "cellulose annotation share"
    )
    recessive_grid(ax)
    handles = [
        Line2D(
            [], [], marker=marker, color=color, markerfacecolor=fill, linestyle="none", label=text
        )
        for _, _, marker, color, fill, text in markers
    ]
    ax.legend(
        handles=handles,
        loc="lower left",
        fontsize=7.6,
        frameon=False,
        bbox_to_anchor=(0.0, -0.22),
        ncol=2,
    )

    # B. Sensitivity table
    sens = results.sensitivity
    ax = fig.add_axes((0.69, 0.42, 0.30, 0.45))
    panel_title(ax, "B. Set definitions: pairwise and mean d_J")
    ax.axis("off")
    order = list(SENSITIVITY_LABELS)
    headers = (f"PD{DASH}\nHIV", f"PD{DASH}\nCRC", f"HIV{DASH}\nCRC", "mean", "core/\nunion")
    xs = (0.50, 0.61, 0.72, 0.83, 0.95)
    for x, text in zip(xs, headers, strict=True):
        ax.text(x, 0.99, text, ha="center", va="top", fontsize=8, fontweight="bold")
    step = 0.86 / len(order)
    for i, slug in enumerate(order):
        y = 0.88 - i * step
        part = sens[sens["definition_slug"].eq(slug)]
        ax.text(0.40, y, wrap(SENSITIVITY_LABELS[slug], 34), ha="right", va="center", fontsize=7.6)
        for x, comparison in zip(xs[:3], ("PD vs. HIV", "PD vs. CRC", "HIV vs. CRC"), strict=True):
            value = part.loc[part["comparison"].eq(comparison), "jaccard_distance"].iloc[0]
            ax.text(x, y, f"{value:.3f}", ha="center", va="center", fontsize=8)
        mean_row = part[part["comparison"].str.startswith("Mean")].iloc[0]
        ax.text(
            xs[3],
            y,
            f"{mean_row['jaccard_distance']:.3f}",
            ha="center",
            va="center",
            fontsize=8,
            fontweight="bold",
        )
        ax.text(
            xs[4],
            y,
            f"{int(mean_row['intersection'])}/{int(mean_row['union'])}",
            ha="center",
            va="center",
            fontsize=8,
        )
    compared = ["n_a", "n_b", "intersection", "union"]
    same_sets = (
        sens.loc[sens["definition_slug"].eq("7_feces_dec_disease_vs_healthy_oriented"), compared]
        .to_numpy()
        .tolist()
        == sens.loc[sens["definition_slug"].eq("8_feces_dec_disease_vs_healthy_strong"), compared]
        .to_numpy()
        .tolist()
    )
    if same_sets:
        ax.text(
            0.0,
            0.0,
            "The two disease-vs-healthy definitions give the same sets.",
            ha="left",
            va="top",
            fontsize=7.4,
            color=MUTED,
        )
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    # C. Interaction tests
    inter = results.interaction
    ax = fig.add_axes((0.22, 0.115, 0.40, 0.18))
    panel_title(ax, "C. Stool against each other compartment, all contrasts")
    picked = inter[
        inter["contrast_population"].eq("all_contrasts")
        & inter["short_condition"].isin(["CRC", POOLED_STRATIFIED])
    ]
    entries = records(picked)
    ys = np.arange(len(entries))[::-1]
    x_min, x_max = 0.01, 1000.0
    for y, row in zip(ys, entries, strict=True):
        name = "All three, stratified" if row["short_condition"] == POOLED_STRATIFIED else "CRC"
        ax.text(
            -0.04,
            y,
            f"{name} · {COMPARTMENT_LABEL[str(row['compartment_b'])]} · {row['unit']}",
            ha="right",
            va="center",
            fontsize=7.6,
            transform=ax.get_yaxis_transform(),
        )
        counts = (
            f"{int(row['a_down'])}↓{int(row['a_up'])}↑ vs "
            f"{int(row['b_down'])}↓{int(row['b_up'])}↑ · "
            f"p = {sci(float(row['fisher_p_two_sided']))}"
        )
        estimate = float(row["conditional_or"])
        if math.isnan(estimate):
            counts = f"{counts} · {not_estimable(row)}"
        else:
            low = float(row["conditional_or_ci95_low"])
            high = float(row["conditional_or_ci95_high"])
            ax.plot([max(low, x_min), min(high, x_max)], [y, y], color=MUTED, linewidth=2)
            if high > x_max:
                ax.plot(x_max, y, ">", color=MUTED, markersize=5, clip_on=False)
            if low < x_min:
                ax.plot(x_min, y, "<", color=MUTED, markersize=5, clip_on=False)
            if math.isinf(estimate) or estimate == 0:
                edge = x_max if math.isinf(estimate) else x_min
                ax.plot(
                    edge,
                    y,
                    ">" if math.isinf(estimate) else "<",
                    color=INK,
                    markersize=5,
                    clip_on=False,
                )
                counts += f" · {zero_cell(row)}"
            else:
                ax.plot(estimate, y, "o", color=INK, markersize=4.5)
        ax.text(
            1.02,
            y,
            counts,
            ha="left",
            va="center",
            fontsize=7.4,
            color=INK2,
            transform=ax.get_yaxis_transform(),
        )
    ax.axvline(1, color=INK2, linestyle="--", linewidth=0.8)
    ax.set_xscale("log")
    ax.set_xlim(x_min, x_max)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel(
        "odds ratio of the ↓ share, stool vs the other site (log scale; exact 95% interval; "
        "arrows: beyond the axis)",
        fontsize=8,
    )
    recessive_grid(ax)
    footnote(
        fig,
        "Experiments within a study are not independent; the study rows collapse them to one "
        "vote per study. Disease-vs-healthy rows keep experiments the contrast manifest classifies "
        "disease_vs_healthy, relabelled so decreased means lower in disease. Intervals are exact "
        f"(Clopper{DASH}Pearson for shares; conditional maximum likelihood for odds ratios). "
        "Stratified rows test a common odds ratio across PD, HIV and CRC, so the conditions' "
        "different sampling of the sites does not enter the estimate. "
        f"{vintage_note(results)}",
    )
    return fig


def detailed_core_and_benchmark(results: Results) -> Figure:
    catalog = results.catalog
    shared = catalog[catalog["n_conditions_feces_dec"] >= 2].copy()
    shared["total_dec"] = shared[[f"{s}_dec_studies" for s in SHORT]].sum(axis=1)
    shared["block"] = np.where(shared["is_3way_core"].astype(bool), 0, 1)
    shared = shared.sort_values(["block", "total_dec", "taxon_name"], ascending=[True, False, True])
    fig = Figure(figsize=(15, 11))
    header(
        fig,
        "Taxa shared by at least two conditions, and the trait-wide convergence benchmark",
        (
            f"{len(shared)} cellulose-annotated taxa are reported decreased in stool in two or "
            f"more conditions; cells give decreased/increased stool studies."
        ),
    )
    ax = fig.add_axes((0.20, 0.09, 0.30, 0.79))
    panel_title(ax, "A. Decreased / increased stool studies")
    ax.axis("off")
    rows = records(shared)
    n = len(rows)
    ax.set_xlim(0, 4.6)
    ax.set_ylim(-1, n + 1)
    for j, short in enumerate(SHORT):
        ax.text(0.5 + j, n + 0.3, short, ha="center", fontsize=10, fontweight="bold")
    ax.text(3.6, n + 0.3, "flags", ha="center", fontsize=9, fontweight="bold")
    y = n - 0.5
    for row in rows:
        ax.text(
            -0.05,
            y,
            f"{row['taxon_name']} ({row['taxonomic_rank']})",
            ha="right",
            va="center",
            fontsize=8,
        )
        for j, short in enumerate(SHORT):
            dec = int(row[f"{short}_dec_studies"])
            inc = int(row[f"{short}_inc_studies"])
            balance = dec - inc
            shade = min(abs(balance), 10) / 10
            color = DOWN if balance > 0 else UP
            face = to_rgba(color, 0.12 + 0.6 * shade) if balance else FAINT
            ax.add_patch(
                FancyBboxPatch(
                    (j + 0.08, y - 0.42),
                    0.84,
                    0.84,
                    boxstyle="round,pad=0.0",
                    facecolor=face,
                    edgecolor="white",
                )
            )
            ax.text(j + 0.5, y, f"{dec}/{inc}", ha="center", va="center", fontsize=7.8, color=INK)
        flags = []
        if row["is_3way_core"]:
            flags.append("core")
        if row["has_direct_trait_edge"]:
            flags.append("direct")
        if row["taxonomic_rank"] in ("genus", "species"):
            flags.append("g/s")
        if row["is_strictly_net_decreased_all3"]:
            flags.append("net↓")
        ax.text(3.15, y, " · ".join(flags), ha="left", va="center", fontsize=7.4, color=INK2)
        y -= 1
    core_rows = int(shared["block"].eq(0).sum())
    ax.axhline(n - core_rows, color=INK2, linewidth=0.8, xmin=0, xmax=0.66)
    fig.text(
        0.58,
        0.20,
        wrap(
            "Panel A fill: blue where decreased studies outnumber increased, red where the "
            "reverse, darker with a larger difference. core = in all three sets; direct = the "
            "taxon has its own annotation edge (otherwise it inherits it from descendants); "
            "g/s = genus or species rank; net↓ = more decreased than increased studies in all "
            "three conditions.",
            80,
        ),
        ha="left",
        va="top",
        fontsize=7.4,
        color=MUTED,
    )
    bench = results.benchmark
    biological = bench[~bench["is_growth_medium"].astype(bool)].reset_index(drop=True)
    ax = fig.add_axes((0.70, 0.30, 0.28, 0.56))
    panel_title(
        ax,
        f"B. Fifteen lowest of {len(biological)} non-MediaDive traits\n"
        "(descriptive, not a test; traits share taxa)",
    )
    top = records(biological.head(15))
    ys = np.arange(len(top))[::-1]
    for y, row in zip(ys, top, strict=True):
        value = float(row["mean_pairwise_jaccard_distance"])
        highlight = row["trait_id"] == results.summary["trait_id"]
        color = INK if highlight else "#9a9993"
        ax.plot([0, value], [y, y], color=color, linewidth=2)
        ax.plot(value, y, "o", color=color, markersize=5)
        ax.text(value + 0.012, y, f"{value:.3f}", va="center", fontsize=7.6)
    ax.set_yticks(ys)
    ax.set_yticklabels([trait_display(row) for row in top], fontsize=7.8)
    ax.set_xlim(
        0, nice_limit(max(float(row["mean_pairwise_jaccard_distance"]) for row in top) + 0.1, 0.1)
    )
    ax.set_xlabel("mean pairwise Jaccard distance")
    recessive_grid(ax)
    footnote(
        fig,
        "Trait names come from taxa_trait_matrix_columns.tsv (object_name). A taxon bears an "
        "annotation when it or a descendant has the KG edge. "
        f"{vintage_note(results)}",
    )
    return fig


# --------------------------------------------------------------------------- slide captions


@dataclass(frozen=True)
class Caption:
    """Slide text for one figure, worded from the same tables as the figure (#556).

    Decks quote these instead of copying numbers by hand, so a rerun of the analysis
    updates the slides with the figures. Each verdict follows the values, including
    direction, unit and population, so the text stays true after a re-pin (#589).
    """

    title: str
    takeaway: str
    caveat: str
    evidence: str


# How captions name each condition row of the tables.
CAPTION_NAMES = {**{short: short for short in SHORT}, "Pooled": "all three pooled"}


def listing(names: Sequence[str]) -> str:
    """'A', 'A and B' or 'A, B and C'."""
    if len(names) < 2:
        return "".join(names)
    return f"{', '.join(names[:-1])} and {names[-1]}"


def sentence(text: str) -> str:
    """``text`` with a capital first letter, to start a sentence."""
    return text[:1].upper() + text[1:]


@dataclass(frozen=True)
class Verdict:
    """One condition's disease-vs-healthy stool study result across the condition definitions.

    ``kind`` is "every" or "some" when the decreased signature more often has the larger
    share at p < 0.05 under every definition or only some, "opposite" or "opposite some"
    for the increased signature, "mixed" when the direction changes with the definition,
    "none" when no definition reaches p < 0.05, and "no data" without rows.
    """

    kind: str
    low: float = math.nan
    high: float = math.nan
    readings: int = 0
    definitions: int = 0


# Verdicts that change with the condition definition.
DEFINITION_DEPENDENT = ("some", "opposite some", "mixed")


#: How a p range names its direction where it is not the decreased signature's (#589, #608).
SPAN_DIRECTION = {
    "mixed": " in opposite directions",
    "opposite some": " toward the increased signature",
}


def p_span(verdict: Verdict, before: str, between: str) -> str:
    """A verdict's p range, naming its direction unless that is the decreased signature's."""
    direction = SPAN_DIRECTION.get(verdict.kind, "")
    return f"{before}{sci(verdict.low)}{between}{sci(verdict.high)}{direction}"


def direction_verdicts(results: Results) -> dict[str, Verdict]:
    """The verdict for PD, HIV, CRC and the pooled three (condition_definition_sensitivity)."""
    table = results.scenarios
    verdicts = {}
    for short in CAPTION_NAMES:
        rows = records(
            table[
                table["contrast_population"].eq(ORIENTED)
                & table["short_condition"].eq(short)
                & table["unit"].eq("study")
            ]
        )
        if not rows:
            verdicts[short] = Verdict("no data")
            continue
        p_values = [float(row["binom_p_two_sided"]) for row in rows]
        down = [
            p < ALPHA and int(r["down"]) > int(r["up"]) for p, r in zip(p_values, rows, strict=True)
        ]
        up = [
            p < ALPHA and int(r["up"]) > int(r["down"]) for p, r in zip(p_values, rows, strict=True)
        ]
        if all(down):
            kind = "every"
        elif all(up):
            kind = "opposite"
        elif any(down) and any(up):
            kind = "mixed"
        elif any(down):
            kind = "some"
        elif any(up):
            kind = "opposite some"
        else:
            kind = "none"
        readings = {(int(row["down"]), int(row["up"]), int(row["tied"])) for row in rows}
        verdicts[short] = Verdict(kind, min(p_values), max(p_values), len(readings), len(rows))
    return verdicts


def fragile_statement(results: Results) -> str:
    """The conditions whose result depends on the condition definition, or ''."""
    shown = {
        short: verdict
        for short, verdict in direction_verdicts(results).items()
        if verdict.kind in DEFINITION_DEPENDENT
    }
    if not shown:
        return ""
    names = {short: "pooled" if short == "Pooled" else short for short in shown}
    spans = {short: p_span(verdict, "", DASH) for short, verdict in shown.items()}
    if len(shown) == 1:
        ranges = "".join(spans.values())
        noun = "result depends"
    else:
        ranges = "; ".join(f"{names[short]} {span}" for short, span in spans.items())
        noun = "results depend"
    return (
        f"the {listing(list(names.values()))} {noun} on the condition definition "
        f"(study-level p {ranges})"
    )


def overlap_statement(results: Results) -> str:
    """The three-way overlap against the whole-signature null (U6), from both endpoints.

    U6's mean-distance p is a lower tail (more overlap) and its core p an upper tail (more
    shared taxa), so p > 1 - alpha on both means less overlap than random.
    """
    observed = results.summary["observed"]
    u6 = results.summary["null_models"]["U6"]
    mean_p = float(u6["mean_pairwise_jaccard_distance"]["p"])
    core_p = float(u6["three_way_core"]["p"])
    shared = (
        f"{observed['three_way_core']} of {observed['three_way_union']} decreased "
        "cellulose-annotated stool taxa are shared by PD, HIV and CRC"
    )
    if mean_p < ALPHA and core_p < ALPHA:
        return f"{shared}; they overlap more than random sets of whole stool signatures"
    if mean_p > 1 - ALPHA and core_p > 1 - ALPHA:
        return f"{shared}; they overlap less than random sets of whole stool signatures"
    if ALPHA <= mean_p <= 1 - ALPHA and ALPHA <= core_p <= 1 - ALPHA:
        return f"{shared}; random sets of whole stool signatures overlap about as much"
    return (
        f"{shared}; against random sets of whole stool signatures, mean-distance "
        f"p = {sci(mean_p)} and core p = {sci(core_p)}"
    )


def crc_compartment_statement(results: Results) -> str:
    """CRC stool against GI mucosa / lumen, by unit and population.

    The tumour-vs-adjacent-tissue experiment is named only while excluding it changes the
    disease-vs-healthy study result, as it did on the v2 contrast manifest (#569, #580).
    """
    crc = interaction_row(results, ALL_CONTRASTS, "CRC", "experiment")
    experiments = interaction_row(results, ORIENTED, "CRC", "experiment")
    studies = interaction_row(results, ORIENTED, "CRC", "study")
    table = results.scenario_interaction
    chosen = table[
        table["scenario"].eq("crc_within_patient_excluded")
        & table["contrast_population"].eq(ORIENTED)
        & table["short_condition"].eq("CRC")
        & table["unit"].eq("study")
        & table["compartment_b"].eq(NON_FECAL_GI)
    ]
    study_p = float(studies["fisher_p_two_sided"])
    without = ""
    if len(chosen):
        without_p = float(only(chosen)["fisher_p_two_sided"])
        if not math.isclose(without_p, study_p):
            without = f" ({sci(without_p)} without the tumour-vs-adjacent-tissue experiment)"
    return (
        f"CRC stool vs GI mucosa / lumen, all contrasts: {int(crc['a_down'])}↓ "
        f"{int(crc['a_up'])}↑ vs {int(crc['b_down'])}↓ {int(crc['b_up'])}↑ experiments "
        f"(odds ratio {float(crc['conditional_or']):.1f}, "
        f"p = {sci(float(crc['fisher_p_two_sided']))}); disease vs healthy: "
        f"p = {sci(float(experiments['fisher_p_two_sided']))} by experiment, {sci(study_p)} "
        f"by study{without}. The sites differ in contrast mix and sampling, so this does not "
        "show anatomical specificity"
    )


def null_endpoints(results: Results, key: str) -> tuple[str, str]:
    """A null's p for the mean distance and for the three-way core."""
    null = results.summary["null_models"][key]
    return (
        sci(float(null["mean_pairwise_jaccard_distance"]["p"])),
        sci(float(null["three_way_core"]["p"])),
    )


def independent_nulls(results: Results) -> str:
    """'p ≤ x' for the weakest mean-distance p of the independent-draw nulls (U1-U4)."""
    nulls = results.summary["null_models"]
    worst = max(
        float(nulls[key]["mean_pairwise_jaccard_distance"]["p"]) for key in ("U1", "U2", "U3", "U4")
    )
    return f"p ≤ {sci(worst)}"


def caption_overlap(results: Results) -> Caption:
    observed = results.summary["observed"]
    sizes = observed["sizes"]
    mean, core = null_endpoints(results, "U6")
    return Caption(
        title=overlap_statement(results),
        takeaway=(
            f"PD {sizes['PD']}, HIV {sizes['HIV']} and CRC {sizes['CRC']} taxa; mean pairwise "
            f"Jaccard distance {float(observed['mean_pairwise_jaccard_distance']):.3f}. "
            f"Independent-draw nulls (U1{DASH}U4): {independent_nulls(results)}; whole-signature "
            f"permutation (U6): p = {mean} (mean), {core} (core)."
        ),
        caveat=(
            "A taxon enters a set once any stool study reports it decreased. Many taxa are "
            "reported together in the same signatures whatever the condition, which "
            "independent-draw nulls ignore."
        ),
        evidence="pairwise_jaccard_significance.tsv",
    )


def caption_shared_taxa(results: Results) -> Caption:
    core = core_taxa(results.catalog)
    n_core = len(core)
    n_direct = int(core["has_direct_trait_edge"].astype(bool).sum())
    net = int(core["is_strictly_net_decreased_all3"].astype(bool).sum())
    low_text, high_text = rank_groups(core)
    sens = results.sensitivity
    means = {
        str(row["definition_slug"]): float(row["jaccard_distance"])
        for row in records(sens[sens["comparison"].str.startswith("Mean")])
    }
    if n_direct == 0:
        direct = "none annotated directly"
        inherit = "Each inherits the cellulose annotation from annotated descendants"
    elif n_direct == n_core:
        direct = "all annotated directly"
        inherit = "All carry their own cellulose annotation edge"
    else:
        direct = f"{n_direct} annotated directly"
        inherit = (
            f"{n_direct} {'carries its' if n_direct == 1 else 'carry their'} own cellulose "
            "annotation edge and the others inherit it from annotated descendants"
        )
    examples = reversed_examples(core)
    reported_up = (
        f"Some are reported increased more often in one condition{examples}. " if examples else ""
    )
    return Caption(
        title=f"The {n_core} taxa in all three sets are {low_text} and {high_text}, {direct}",
        takeaway=(
            f"{inherit}; {net} of the {n_core} are reported decreased more often than increased "
            "in all three conditions."
        ),
        caveat=(
            f"{reported_up}Genus and species ranks alone give mean d_J "
            f"{means['3_feces_dec_genus_species']:.3f}, against "
            f"{means['1_feces_dec_all_ranks']:.3f} for all ranks."
        ),
        evidence="cellulose_decreased_taxa_catalog.tsv",
    )


def caption_direction(results: Results) -> Caption:
    verdicts = direction_verdicts(results)
    groups: dict[str, list[str]] = {}
    for short, verdict in verdicts.items():
        groups.setdefault(verdict.kind, []).append(CAPTION_NAMES[short])
    lead = (
        "Disease-vs-healthy stool studies: the decreased signature more often carries the "
        "larger cellulose annotation share"
    )
    clauses = []
    if groups.get("every"):
        clauses.append(f"in {listing(groups['every'])} under every condition definition")
    if groups.get("some"):
        clauses.append(f"in {listing(groups['some'])} only under some")
    title = f"{lead} {'; '.join(clauses)}" if clauses else f"{lead} in no condition at p < 0.05"
    others = []
    reversed_clauses = []
    if groups.get("opposite"):
        reversed_clauses.append(f"in {listing(groups['opposite'])} under every definition")
    if groups.get("opposite some"):
        reversed_clauses.append(f"in {listing(groups['opposite some'])} only under some")
    if reversed_clauses:
        others.append(f"the increased signature does {' and '.join(reversed_clauses)}")
    if groups.get("mixed"):
        changes = "changes" if len(groups["mixed"]) == 1 else "change"
        others.append(f"{listing(groups['mixed'])} {changes} direction with the definition")
    if groups.get("none") and clauses:
        reaches = "reaches" if len(groups["none"]) == 1 else "reach"
        others.append(f"{listing(groups['none'])} {reaches} p < 0.05 under no definition")
    if groups.get("no data"):
        others.append(f"no data for {listing(groups['no data'])}")
    if others:
        title = f"{title}; {'; '.join(others)}"

    healthy = stool_rows(results, ORIENTED)
    counts = ", ".join(
        f"{'all three' if short == 'Pooled' else short} {int(healthy[short]['study_down'])} : "
        f"{int(healthy[short]['study_up'])} : {int(healthy[short]['study_tied'])} "
        f"(p = {sci(float(healthy[short]['study_binom_p_two_sided']))})"
        for short in CAPTION_NAMES
        if short in healthy
    )
    ranges = "".join(
        f"{sentence(CAPTION_NAMES[short])}: {p_span(v, 'p from ', ' to ')} across "
        f"{v.readings} distinct readings of the {v.definitions} sensitivity scenarios. "
        for short, v in verdicts.items()
        if v.kind in DEFINITION_DEPENDENT
    )
    return Caption(
        title=title,
        takeaway=f"Exact condition labels, down : up : tied studies: {counts}.",
        caveat=f"{ranges}{crc_compartment_statement(results)}.",
        evidence="condition_definition_sensitivity.tsv",
    )


def caption_ranking(results: Results) -> Caption:
    bench = results.benchmark
    rows = records(bench[~bench["is_growth_medium"].astype(bool)].reset_index(drop=True))
    cell = next(row for row in rows if row["trait_id"] == results.summary["trait_id"])
    first = rows[0]
    median = float(np.median([float(row["mean_pairwise_jaccard_distance"]) for row in rows]))
    leader = (
        ""
        if first["trait_id"] == cell["trait_id"]
        else f"{first['trait_label']} ({first['trait_id']}) "
        f"{float(first['mean_pairwise_jaccard_distance']):.3f}, "
    )
    label = str(cell["trait_label"])
    return Caption(
        title=(
            f"{label[:1].upper()}{label[1:]} ranks {int(cell['rank_excluding_growth_media'])} of "
            f"{len(rows)} non-growth-medium traits by overlap of stool-decreased taxa across "
            "PD, HIV and CRC; the ranking describes this corpus and is not a test"
        ),
        takeaway=(
            f"Mean pairwise Jaccard distance: {leader}{label} "
            f"{float(cell['mean_pairwise_jaccard_distance']):.3f}; median of the {len(rows)} "
            f"traits {median:.3f}."
        ),
        caveat=(
            "Traits share taxa, so their distances are not independent, and the trait was picked "
            "from this corpus as an exploratory lead."
        ),
        evidence="trait_convergence_benchmark.tsv",
    )


def caption_nulls(results: Results) -> Caption:
    iterations = int(results.summary["null_models"]["U1"]["iterations"])
    matched = int(results.summary["null_models"]["U6_size_matched"]["iterations"])
    u5 = null_endpoints(results, "U5")
    u6 = null_endpoints(results, "U6")
    size_matched = null_endpoints(results, "U6_size_matched")
    return Caption(
        title="Overlap of stool-decreased cellulose-annotated taxa against six null models",
        takeaway=(
            f"Independent-draw nulls (U1{DASH}U4): every mean {independent_nulls(results)}. "
            f"U5 (prevalence-weighted): {u5[0]} (mean), {u5[1]} (core). U6 (whole signatures): "
            f"{u6[0]}, {u6[1]}; size-matched: {size_matched[0]}, {size_matched[1]}."
        ),
        caveat=(
            f"U1{DASH}U3 pair p-values are exact hypergeometric tails, the others Monte Carlo "
            f"with {iterations:,} draws ({matched:,} kept for the size-matched U6); only U6 "
            "keeps co-reported taxa together."
        ),
        evidence="pairwise_jaccard_significance.tsv",
    )


def caption_body_sites(results: Results) -> Caption:
    def odds(population: str, unit: str) -> str:
        row = interaction_row(results, population, POOLED_STRATIFIED, unit)
        return f"{float(row['conditional_or']):.2f} (p = {sci(float(row['fisher_p_two_sided']))})"

    table = results.compartments
    elsewhere = [
        (str(row["short_condition"]), str(row["compartment"]), sum(untied(row, "experiment")))
        for row in records(
            table[
                table["contrast_population"].eq(ALL_CONTRASTS)
                & ~table["compartment"].eq("Feces")
                & table["short_condition"].isin(SHORT)
            ]
        )
    ]
    crc_gi = sum(n for short, site, n in elsewhere if (short, site) == ("CRC", NON_FECAL_GI))
    others = [n for short, site, n in elsewhere if (short, site) != ("CRC", NON_FECAL_GI)]
    return Caption(
        title=(
            "Direction of the cellulose annotation shift by body site, inference unit and "
            "contrast type"
        ),
        takeaway=(
            "Stool vs GI mucosa / lumen, common odds ratio stratified by condition: all contrasts "
            f"{odds(ALL_CONTRASTS, 'experiment')} by experiment and "
            f"{odds(ALL_CONTRASTS, 'study')} by study; disease vs healthy "
            f"{odds(ORIENTED, 'study')} by study."
        ),
        caveat=(
            "Experiments within a study are not independent. Outside stool, single conditions "
            f"have at most {max(others, default=0)} untied pairs, except CRC GI mucosa / lumen "
            f"({crc_gi})."
        ),
        evidence="compartment_direction_interaction.tsv",
    )


def caption_shared_by_two(results: Results) -> Caption:
    observed = results.summary["observed"]
    bench = results.benchmark
    n_traits = int((~bench["is_growth_medium"].astype(bool)).sum())
    direct = direct_core_text(results)
    return Caption(
        title=(
            f"{observed['two_or_more_way']} cellulose-annotated taxa are reported decreased in "
            f"stool in two or more conditions, {observed['three_way_core']} in all three"
        ),
        takeaway=(
            "Cells give decreased / increased stool studies per condition; flags mark the core, "
            "direct annotation edges, genus or species rank and net-decreased taxa."
        ),
        caveat=(
            f"{direct[:1].upper()}{direct[1:]}. The right panel ranks {n_traits} "
            "non-growth-medium traits descriptively; it is not a test."
        ),
        evidence="cellulose_decreased_taxa_catalog.tsv",
    )


# --------------------------------------------------------------------------- main


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    results = load(args.analysis_dir)
    # The stems describe what each figure shows, not a conclusion (#570).
    plates: Mapping[
        str, list[tuple[str, Callable[[Results], Figure], Callable[[Results], Caption]]]
    ] = {
        "figures_graphical": [
            ("fig1_overlap_and_null_models", graphical_convergence, caption_overlap),
            ("fig2_shared_taxa_and_set_definitions", graphical_core, caption_shared_taxa),
            ("fig3_direction_by_body_site", graphical_compartments, caption_direction),
            ("fig4_trait_overlap_ranking", graphical_benchmark, caption_ranking),
        ],
        "figures": [
            ("fig1_pairwise_and_3way_jaccard_nulls", detailed_nulls, caption_nulls),
            (
                "fig2_direction_by_body_site_and_definition",
                detailed_compartments,
                caption_body_sites,
            ),
            (
                "fig3_core_taxa_and_trait_benchmark",
                detailed_core_and_benchmark,
                caption_shared_by_two,
            ),
        ],
    }
    combined = {
        "figures_graphical": "cellulose_decr_graphical_figures.pdf",
        "figures": "cellulose_decr_figures.pdf",
    }
    outputs: list[str] = []
    captions: dict[str, dict[str, str]] = {}
    for folder, entries in plates.items():
        with PdfPages(
            args.analysis_dir / combined[folder], metadata={"CreationDate": None}
        ) as book:
            for stem, draw, caption in entries:
                fig = draw(results)
                written = save(fig, args.analysis_dir / folder, stem, book)
                outputs.extend(f"{folder}/{path.name}" for path in written)
                captions[f"{folder}/{stem}.svg"] = asdict(caption(results))
                print(f"wrote {', '.join(path.name for path in written)}")
        outputs.append(combined[folder])
    # Slide text for decks, from the same tables (scripts/build_slide_deck.py; #556).
    captions_path = args.analysis_dir / "figures" / "figure_captions.json"
    statements = {
        "overlap": overlap_statement(results),
        "crc_compartments": crc_compartment_statement(results),
        "fragile": fragile_statement(results),
    }
    captions_path.write_text(
        json.dumps({"figures": captions, "statements": statements}, indent=1, ensure_ascii=False)
        + "\n",
        encoding="utf-8",
    )
    outputs.append("figures/figure_captions.json")
    provenance = {
        "script": "src/mtp/cellulose_decr/figures.py",
        "script_sha256": sha256(Path(__file__).resolve()),
        "code_commit": args.code_commit,
        "matplotlib": matplotlib.__version__,
        "inputs": {name: sha256(args.analysis_dir / name) for name in FIGURE_INPUTS},
        "outputs": sorted(outputs),
    }
    (args.analysis_dir / "figures" / "figures_provenance.json").write_text(
        json.dumps(provenance, indent=1) + "\n", encoding="utf-8"
    )
    return 0


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


if __name__ == "__main__":
    sys.exit(main())
