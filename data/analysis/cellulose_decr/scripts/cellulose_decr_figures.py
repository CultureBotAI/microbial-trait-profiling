#!/usr/bin/env python3
"""Draw the ``cellulose_decr`` figures from the analysis tables (#552).

Every number on a figure is read from the TSV and JSON files that
``cellulose_decr_analysis.py`` writes to ``data/analysis/cellulose_decr/``;
none is typed in. Writes the statistical reference plates to ``figures/`` and
the presentation plates to ``figures_graphical/`` as SVG, PDF and PNG, and a
combined PDF of each set.

Trait names are the ``object_name`` labels joined from
``data/derived/taxa_trait_matrix_columns.tsv`` by the analysis. Captions keep
to methods/biological_interpretation.md: these are annotation counts of
curated signature taxa, not measurements of cellulose degradation.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import textwrap
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
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

ANALYSIS_DIR = Path(__file__).resolve().parents[1]

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
NON_FECAL_GI = "Non-fecal GI (mucosa/lumen)"
ORAL = "Oral / Oropharyngeal"
OTHER = "Other non-GI"
COMPARTMENT_LABEL = {
    "Feces": "Stool",
    NON_FECAL_GI: "GI mucosa / lumen",
    ORAL: "Oral",
    OTHER: "Other non-GI",
}
PNG_DPI = 150

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
        "svg.fonttype": "none",
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


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0] if __doc__ else None)
    parser.add_argument("--analysis-dir", type=Path, default=ANALYSIS_DIR)
    return parser.parse_args(argv)


def load(directory: Path) -> Results:
    def tsv(name: str) -> pd.DataFrame:
        frame: pd.DataFrame = pd.read_csv(directory / name, sep="\t")
        return frame

    summary = json.loads((directory / "cellulose_decr_summary.json").read_text(encoding="utf-8"))
    return Results(
        summary=summary,
        pairwise=tsv("pairwise_jaccard_significance.tsv"),
        sensitivity=tsv("sensitivity_definitions_jaccard.tsv"),
        compartments=tsv("compartment_jaccard_and_direction_stats.tsv"),
        interaction=tsv("compartment_direction_interaction.tsv"),
        catalog=tsv("cellulose_decreased_taxa_catalog.tsv"),
        benchmark=tsv("trait_convergence_benchmark.tsv"),
    )


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
    """Readable p-value: three decimals above 0.001, a x 10^k form below."""
    if value is None or math.isnan(value):
        return "n/a"
    if value >= 0.001:
        return f"{value:.3f}" if value < 0.1 else f"{value:.2f}"
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
    fig.text(0.015, 0.975, title, ha="left", va="top", fontsize=15, fontweight="bold", color=INK)
    fig.text(0.015, 0.928, wrap(subtitle, 175), ha="left", va="top", fontsize=10.5, color=INK2)


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
    fig.savefig(written[0])
    fig.savefig(written[1])
    fig.savefig(written[2], dpi=PNG_DPI)
    add_page = cast(Callable[[Figure], None], book.savefig)
    add_page(fig)
    return written


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
    ("U3", "U3  uniform draw from the taxa ever reported decreased in stool"),
    ("U4", "U4  random trait placed on all stool-decreased taxa of the three conditions"),
    ("U5", "U5  taxa weighted by how often other stool studies report them decreased"),
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
    studies = stool_rows(results, "all_contrasts")

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
        n_studies = int(studies[short]["n_studies_dec"])
        ax.text(x, y + 0.045, short, ha="center", va="center", fontsize=13, fontweight="bold")
        ax.text(
            x,
            y - 0.035,
            f"{observed['sizes'][short]} taxa\n{n_studies} studies",
            ha="center",
            va="center",
            fontsize=8.5,
            color=INK2,
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
    for y, (key, text) in zip(ys, NULL_ROWS, strict=True):
        record = nulls[key]["mean_pairwise_jaccard_distance"]
        ax.plot([record["p2_5"], record["p97_5"]], [y, y], color=MUTED, linewidth=3, zorder=1)
        ax.plot(record["expected"], y, "o", color=INK2, markersize=7, zorder=2)
        ax.text(
            1.0,
            y - 0.33,
            mc_p(record["p"], nulls[key]["iterations"]),
            ha="right",
            va="top",
            fontsize=8,
            color=INK,
        )
        ax.text(
            -0.03,
            y,
            wrap(null_label(results, key, text), 34),
            ha="right",
            va="center",
            fontsize=7.6,
            color=INK2,
            transform=ax.get_yaxis_transform(),
        )
    ax.axvline(mean_d, color=INK, linewidth=1.6, zorder=3)
    ax.text(
        mean_d + 0.01,
        len(NULL_ROWS) - 0.15,
        f"observed {mean_d:.3f}",
        ha="left",
        va="center",
        fontsize=8.5,
        fontweight="bold",
    )
    ax.set_xlim(0.38, 1.0)
    ax.set_ylim(-0.8, len(NULL_ROWS) + 0.2)
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


def core_taxa(catalog: pd.DataFrame) -> pd.DataFrame:
    core = catalog[catalog["is_3way_core"].astype(bool)].copy()
    core["total_dec"] = core[[f"{s}_dec_studies" for s in SHORT]].sum(axis=1)
    core["group"] = np.where(core["is_direct_or_genus_species"].astype(bool), 0, 1)
    return core.sort_values(["group", "total_dec", "taxon_name"], ascending=[True, False, True])


SENSITIVITY_LABELS = {
    "1_feces_dec_all_ranks": "Primary: all ranks",
    "3_feces_dec_genus_species": "Genus and species ranks only",
    "2_feces_dec_z_gt_0": "Profiles with cellulose z > 0 only",
    "7_feces_dec_disease_vs_healthy_oriented": "Disease-vs-healthy contrasts only",
    "8_feces_dec_disease_vs_healthy_strong": "Disease-vs-healthy, strong evidence",
    "5_all_body_sites_dec": "All body sites, not only stool",
    "4_feces_net_decreased": "Net decreased (more dec than inc profiles)",
    "4b_feces_net_decreased_by_studies": "Net decreased (more dec than inc studies)",
    "6_baseline_all_feces_dec_taxa": "Baseline: all stool-decreased taxa, any trait",
}


def graphical_core(results: Results) -> Figure:
    catalog = results.catalog
    core = core_taxa(catalog)
    n_genus = int(core["group"].eq(0).sum())
    n_higher = int(core["group"].eq(1).sum())
    net = core[core["is_strictly_net_decreased_all3"].astype(bool)]
    sens = results.sensitivity
    mean_rows = {
        str(row["definition_slug"]): row
        for row in records(sens[sens["comparison"].str.startswith("Mean")])
    }
    gs_mean = float(mean_rows["3_feces_dec_genus_species"]["jaccard_distance"])
    fig = Figure(figsize=(15, 7.6))
    header(
        fig,
        (
            f"The {len(core)} taxa reported decreased in stool in all three conditions: "
            f"{n_genus} genera and {n_higher} higher-rank nodes"
        ),
        (
            f"{len(net)} of them had more decreased than increased stool studies in every "
            "condition. "
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
        color = INK if primary else (MUTED if baseline else "#9a9993")
        ax.barh(y, value, height=0.62, color=color)
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
            "Upper block: genera; lower block: family, order, class and phylum nodes, whose "
            "annotations are rolled up from descendant taxa. Net ↓ = more decreased than "
            "increased studies in every condition.",
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
        "often (for example Oscillospiraceae in PD and Bacteroidota in CRC). "
        f"{vintage_note(results)}",
    )
    return fig


def graphical_compartments(results: Results) -> Figure:
    fig = Figure(figsize=(15, 7.4))
    pooled = compartment_row(results, "all_contrasts", "Feces", "Pooled")
    if pooled is None:
        raise ValueError("no pooled stool row in the compartment table")
    header(
        fig,
        "The decreased-over-increased shift is clearest in stool; other sites have too few "
        "comparisons to judge",
        (
            f"Across PD, HIV and CRC stool, {int(pooled['paired_exp_down_dec_gt_inc'])} of "
            f"{int(pooled['paired_exp_down_dec_gt_inc'] + pooled['paired_exp_up_inc_gt_dec'])} "
            "untied paired experiments had the larger cellulose annotation share in the decreased "
            f"signature (p = {sci(pooled['binom_p_two_sided'])}); outside stool most strata have "
            "fewer than six untied pairs."
        ),
    )

    # A. Butterfly of experiment counts
    ax = fig.add_axes((0.15, 0.17, 0.38, 0.60))
    panel_title(ax, "A. Paired experiments: decreased signature larger (↓) or smaller (↑)")
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
    limit = 120
    for y, (compartment, short, row) in zip(positions, strata, strict=True):
        down = int(row["paired_exp_down_dec_gt_inc"])
        up = int(row["paired_exp_up_inc_gt_dec"])
        ax.barh(y, -down, height=0.6, color=DOWN)
        ax.barh(y, up, height=0.6, color=UP)
        ax.text(-down - 2, y, str(down), ha="right", va="center", fontsize=8, color=INK)
        ax.text(up + 2, y, str(up), ha="left", va="center", fontsize=8, color=INK)
        label = "All three" if short == "Pooled" else short
        ax.text(
            -limit - 4,
            y,
            f"{COMPARTMENT_LABEL[compartment]} · {label}",
            ha="right",
            va="center",
            fontsize=8.5,
            fontweight="bold" if short == "Pooled" else "normal",
        )
        study = f"studies {int(row['study_down'])}↓ {int(row['study_up'])}↑"
        detail = f"p = {sci(row['binom_p_two_sided'])} · {study}"
        if row["min_attainable_binom_p_two_sided"] > 0.05:
            detail = f"{study} · too few pairs for p < 0.05"
        ax.text(limit + 4, y, detail, ha="left", va="center", fontsize=7.8, color=INK2)
    ax.axvline(0, color=INK, linewidth=0.8)
    ax.set_xlim(-limit, limit)
    ax.set_ylim(min(positions) - 0.7, max(positions) + 0.7)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xticks([-100, -50, 0, 50, 100])
    ax.set_xticklabels(["100", "50", "0", "50", "100"])
    ax.set_xlabel("untied paired experiments (← decreased larger · increased larger →)")
    recessive_grid(ax)

    # B. Stool versus GI mucosa/lumen
    inter = results.interaction
    ax = fig.add_axes((0.845, 0.51, 0.14, 0.24))
    panel_title(ax, "B. Stool vs GI\nmucosa / lumen")
    chosen = [
        ("CRC", "experiment", "all_contrasts", "CRC · experiments"),
        ("CRC", "study", "all_contrasts", "CRC · studies"),
        ("Pooled", "experiment", "all_contrasts", "All three · experiments"),
        ("Pooled", "study", "all_contrasts", "All three · studies"),
        ("CRC", "experiment", "disease_vs_healthy_oriented", "CRC · disease vs healthy"),
    ]
    ys = np.arange(len(chosen))[::-1]
    for y, (short, unit, population, label) in zip(ys, chosen, strict=True):
        row = only(
            inter[
                inter["short_condition"].eq(short)
                & inter["unit"].eq(unit)
                & inter["contrast_population"].eq(population)
                & inter["compartment_b"].eq(NON_FECAL_GI)
            ]
        )
        low = float(row["conditional_or_ci95_low"])
        high = float(row["conditional_or_ci95_high"])
        ax.plot([low, min(high, 400)], [y, y], color=MUTED, linewidth=2.5)
        ax.plot(row["conditional_or"], y, "o", color=INK, markersize=6)
        ax.text(
            -0.05,
            y,
            f"{label}\np = {sci(row['fisher_p_two_sided'])}",
            ha="right",
            va="center",
            fontsize=7.8,
            transform=ax.get_yaxis_transform(),
        )
    ax.axvline(1, color=INK2, linewidth=0.8, linestyle="--")
    ax.set_xscale("log")
    ax.set_xlim(0.5, 400)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("odds ratio, ↓ share stool vs GI\n(exact 95% interval)", fontsize=8)
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

    footnote(
        fig,
        "Paired experiment: one BugSigDB experiment with both an increased and a decreased "
        "signature; the share compared is cellulose annotations over all retained non-MediaDive "
        "annotations. Ties excluded. Study rows take the sign of each study's mean paired "
        "difference. The Fisher tests treat units as independent; a few studies sampled both "
        "compartments, and the sites differ in sampling effort (C), so a difference between "
        "compartments is not evidence of anatomical specificity. "
        f"{vintage_note(results)}",
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
    bins = [round(0.30 + 0.02 * step, 2) for step in range(31)]
    ax.hist(values, bins=bins, color="#c9c8c2", edgecolor="white")
    median = float(np.median(values))
    ax.axvline(median, color=INK2, linestyle="--", linewidth=1)
    top_y = ax.get_ylim()[1]
    ax.text(0.735, top_y * 0.95, f"dashed: median {median:.3f}", fontsize=8.5, color=INK2)
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
    ax.set_xlim(0, 0.85)
    ax.set_xlabel("mean pairwise Jaccard distance (lower = more overlap)")
    recessive_grid(ax)
    footnote(
        fig,
        "Trait names are the object_name labels in data/derived/taxa_trait_matrix_columns.tsv. "
        "Excluded: MediaDive growth-medium columns (METPO:2000517|mediadive.medium:*). Several "
        "listed traits are morphology or temperature classes, not metabolic processes. "
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
    header(
        fig,
        "Overlap of stool-decreased cellulose-annotated taxa against six null models",
        (
            f"Nulls that draw taxa independently (U1{DASH}U4) make the overlap look extreme; nulls "
            "that respect how often taxa are reported (U5) and which taxa are reported together "
            "(U6) do not."
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
                    facecolor="#dde9f8" if strong else FAINT,
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
            "Observed values are in the first column. Shaded blue: p < 0.05. "
            f"U1{DASH}U3 pair p-values are exact hypergeometric tails; all other p are Monte "
            "Carlo. U4 places a random trait on the taxa reported decreased in stool in the "
            "three conditions; U5 weights taxa by "
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
    ax.set_xlim(0, 30)
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
        "Names are KG-Microbe matrix taxon names; family-and-above nodes carry annotations rolled "
        f"up from descendants. {vintage_note(results)}",
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
    panel_title(ax, "A. Share of untied units with the decreased signature larger")
    rows = []
    for compartment in ("Feces", NON_FECAL_GI, ORAL, OTHER):
        for short in ("PD", "HIV", "CRC", "Pooled"):
            rows.append((compartment, short))
    ys = np.arange(len(rows))[::-1].astype(float)
    markers = (
        ("all_contrasts", "experiment", "o", INK, "all contrasts · experiments"),
        ("all_contrasts", "study", "s", "#9a9993", "all contrasts · studies"),
        (
            "disease_vs_healthy_oriented",
            "experiment",
            "D",
            DOWN,
            "disease vs healthy · experiments",
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
        for offset, (population, unit, marker, color, _) in zip(
            (0.22, 0.0, -0.22), markers, strict=True
        ):
            row = compartment_row(results, population, compartment, short)
            if row is None:
                continue
            if unit == "experiment":
                down, up = (
                    int(row["paired_exp_down_dec_gt_inc"]),
                    int(row["paired_exp_up_inc_gt_dec"]),
                )
            else:
                down, up = int(row["study_down"]), int(row["study_up"])
            share, low, high = share_interval(down, up)
            if math.isnan(share):
                continue
            ax.plot([low, high], [y + offset, y + offset], color=color, linewidth=1.2, alpha=0.8)
            ax.plot(share, y + offset, marker, color=color, markersize=4.5)
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
    ax.set_xlabel("share of untied units with decreased signature larger")
    recessive_grid(ax)
    handles = [
        Line2D([], [], marker=marker, color=color, linestyle="none", label=text)
        for _, _, marker, color, text in markers
    ]
    ax.legend(
        handles=handles,
        loc="lower left",
        fontsize=7.6,
        frameon=False,
        bbox_to_anchor=(0.0, -0.2),
        ncol=3,
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
    for i, slug in enumerate(order):
        y = 0.88 - i * 0.095
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
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    # C. Interaction tests
    inter = results.interaction
    ax = fig.add_axes((0.22, 0.115, 0.40, 0.18))
    panel_title(
        ax,
        "C. Stool against each other compartment (odds ratio of the ↓ share, exact 95% interval)",
    )
    picked = inter[
        inter["contrast_population"].eq("all_contrasts")
        & inter["short_condition"].isin(["CRC", "Pooled"])
    ]
    entries = records(picked)
    ys = np.arange(len(entries))[::-1]
    for y, row in zip(ys, entries, strict=True):
        name = "All three" if row["short_condition"] == "Pooled" else "CRC"
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
            ax.text(
                0.6, y, "not estimable (no untied pairs)", va="center", fontsize=7.4, color=MUTED
            )
            continue
        low = max(float(row["conditional_or_ci95_low"]), 0.011)
        high = min(float(row["conditional_or_ci95_high"]), 990.0)
        ax.plot([low, high], [y, y], color=MUTED, linewidth=2)
        if math.isinf(estimate):
            ax.plot(990, y, ">", color=INK, markersize=5)
            counts += " · OR unbounded (no ↑ outside stool)"
        elif estimate == 0:
            ax.plot(0.011, y, "<", color=INK, markersize=5)
            counts += " · OR 0"
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
    ax.set_xlim(0.01, 1000)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    recessive_grid(ax)
    footnote(
        fig,
        "Experiments within a study are not independent; the study rows collapse them to one "
        "vote per study. Disease-vs-healthy rows keep experiments the contrast manifest classifies "
        "disease_vs_healthy, relabelled so decreased means lower in disease. Intervals are exact "
        f"(Clopper{DASH}Pearson for shares; conditional maximum likelihood for odds ratios). "
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
        if row["is_direct_or_genus_species"]:
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
            "reverse, darker with a larger difference. core = in all three sets; g/s = genus or "
            "species rank; net↓ = more decreased than increased studies in all three conditions.",
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
    panel_title(ax, f"B. Fifteen lowest of {len(biological)} non-MediaDive traits")
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
    ax.set_xlim(0, 0.8)
    ax.set_xlabel("mean pairwise Jaccard distance")
    recessive_grid(ax)
    footnote(
        fig,
        "Trait names come from taxa_trait_matrix_columns.tsv (object_name). Family-and-above "
        f"nodes carry rolled-up annotations. {vintage_note(results)}",
    )
    return fig


# --------------------------------------------------------------------------- main


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    results = load(args.analysis_dir)
    plates: Mapping[str, list[tuple[str, Any]]] = {
        "figures_graphical": [
            ("fig1_surprising_convergence_network_venn", graphical_convergence),
            ("fig2_genus_species_core_guild", graphical_core),
            ("fig3_feces_compartment_specificity", graphical_compartments),
            ("fig4_top_functional_pathway_podium", graphical_benchmark),
        ],
        "figures": [
            ("fig1_pairwise_and_3way_jaccard_nulls", detailed_nulls),
            ("fig2_compartment_specificity_and_sensitivity", detailed_compartments),
            ("fig3_core_taxa_and_trait_benchmark", detailed_core_and_benchmark),
        ],
    }
    combined = {
        "figures_graphical": "cellulose_decr_graphical_figures.pdf",
        "figures": "cellulose_decr_figures.pdf",
    }
    for folder, entries in plates.items():
        with PdfPages(args.analysis_dir / combined[folder]) as book:
            for stem, draw in entries:
                fig = draw(results)
                written = save(fig, args.analysis_dir / folder, stem, book)
                print(f"wrote {', '.join(path.name for path in written)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
