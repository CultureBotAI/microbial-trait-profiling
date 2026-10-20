#!/usr/bin/env python3
"""Draw the deck's two methods plates and seven Tufte plates from the three-arm rerun data.

The figures are computed from the TSV outputs of the rerun on main (``RERUN_COMMIT``)
under ``trait-visualization/`` and ``data/derived/``, many of them gitignored, and are
written as PNG and SVG to ``data/derived/tufte_plates/``. Plates IV and V also quote this
checkout's cellulose_decr tables and the ``figure_captions.json`` that ``just
cellulose-decr-figures`` writes beside them. ``plates_provenance.json`` records the
drawing commit and each plate's inputs with their hashes. ``just figure-review-plates``
runs this script; ``scripts/build_slide_deck.py`` places the plates in the deck (#573).
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import re
import textwrap
from collections.abc import Mapping, Sequence
from pathlib import Path
from statistics import median
from typing import Any

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

from mtp.paths import file_sha256
from mtp.quoted_results import (
    CALIBRATION_CRITERIA,
    GUILD_CORRELATION,
    ORA_DEPLETION,
    ORA_ENRICHMENT,
    McNemar,
    OraCalls,
    screened_rejection_rate,
)
from mtp.table1_directions import (
    ALL_STUDIES,
    CELLULOSE,
    COMPARISON,
    CONTRAST_PROVENANCE,
    ORIENTED,
    Population,
    populations,
)
from mtp.trait_guilds import OXYGEN_GUILDS

REPO_ROOT = Path(__file__).resolve().parents[1]
PLATES_DIR = Path("data/derived/tufte_plates")
CELLULOSE_DIR = Path("data/analysis/cellulose_decr")
TRAIT_COLUMNS = Path("data/derived/taxa_trait_matrix_columns.tsv")
# The rerun the plates' tables come from. The deck stamps and links the plate slides at the
# same revision and refuses inputs that differ from it (#637); the two move together (#579).
RERUN_COMMIT = "622e5c71d8c8e4190981a98c4220823b5e6060df"
RERUN_SOURCE = f"Source: Rerun on main ({RERUN_COMMIT[:7]})"

# Tufte palette: muted, high-contrast, colorblind-safe ink on warm paper/clean white
BG_COLOR = "#fcfbf8"
PANEL_BG = "#ffffff"
INK = "#162433"
MUTED = "#526273"
LIGHT_RULE = "#dce3ea"
MID_RULE = "#aab7c4"
INC_COLOR = "#b8332a"  # Vermilion / brick red for increased-side shift
INC_LIGHT = "#f5d7d4"
DEC_COLOR = "#1b6ca8"  # Deep Prussian blue for decreased-side shift
DEC_LIGHT = "#d4e6f4"
NEUTRAL_INK = "#6c7a89"
GOLD_ACCENT = "#b07d12"
TEAL_ACCENT = "#157a6e"
PURPLE_ACCENT = "#6a4c93"


# The files the plate being drawn reads; main() records them with their hashes (#637).
INPUTS_READ: set[Path] = set()
# Counts the plates draw that the deck's slide text quotes (#640).
PLATE_COUNTS: dict[str, int] = {}


def read_tsv(path: Path) -> list[dict[str, str]]:
    INPUTS_READ.add(path)
    with path.open(newline="", encoding="utf-8") as f:
        lines = [line for line in f if not line.startswith("#")]
        return list(csv.DictReader(lines, delimiter="\t"))


def read_json(path: Path) -> Any:
    INPUTS_READ.add(path)
    return json.loads(path.read_text(encoding="utf-8"))


def fmt_p(value: float) -> str:
    """A p-value to two significant digits, or more where fewer would round across 0.05."""
    for digits in range(2, 8):
        text = f"{value:#.{digits}g}" if value >= 0.001 else f"{value:.{digits - 1}e}"
        if (float(text) < 0.05) == (value < 0.05):
            return text
    return repr(value)


SUPERSCRIPT = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def sci(value: float) -> str:
    """A p-value as fmt_p gives it, with any power of ten written out as 10⁻⁶."""
    text = fmt_p(value)
    if "e" not in text:
        return text
    mantissa, exponent = text.split("e")
    return f"{mantissa} \N{MULTIPLICATION SIGN} 10{str(int(exponent)).translate(SUPERSCRIPT)}"


def ora_calls(path: Path) -> OraCalls:
    """An ORA table's tests and calls, from the lines heading it (#579)."""
    INPUTS_READ.add(path)
    with path.open(encoding="utf-8") as f:
        return OraCalls.from_header(f)


def screening_bias(repo: Path) -> str:
    """The screened sign test's rejection rate under the complete null, "0.091" (#579)."""
    return f"{screened_rejection_rate(read_tsv(repo / CALIBRATION_CRITERIA)):.3f}"


def stratum_size(rows: Mapping[str, Mapping[str, str]]) -> tuple[int, int]:
    """A paired-contrast stratum's studies and pairs: the most any trait is scored on."""
    studies = max(int(float(r["study_count"])) for r in rows.values())
    pairs = max(int(float(r["pair_count"])) for r in rows.values())
    return studies, pairs


def stratum_text(rows: Mapping[str, Mapping[str, str]]) -> str:
    """'959 studies · 2,335 pairs', from the stratum's table (#579)."""
    studies, pairs = stratum_size(rows)
    return f"{studies:,} studies · {pairs:,} pairs"


def pd_split(rows: Sequence[Mapping[str, str]], per_side: int = 2) -> str:
    """Plate V's note on Parkinson's disease: the genera in the panel with the most PD
    profiles on each side, from the switch table's top-condition columns (#680)."""

    def pd_profiles(conditions: str) -> int:
        for part in conditions.split(";"):
            found = re.match(r"\s*Parkinson's disease \((\d+) profiles", part)
            if found:
                return int(found.group(1))
        return 0

    sides = []
    for column, arrow in (("top_increased_conditions", "↑"), ("top_decreased_conditions", "↓")):
        counts = sorted(
            ((pd_profiles(r[column]), r["taxon_name"]) for r in rows), key=lambda c: -c[0]
        )
        named = [f"{name} ({n}{arrow})" for n, name in counts[:per_side] if n]
        sides.append(" & ".join(named) or "none")
    return " vs ".join(sides)


def table1_signs(repo: Path) -> dict[str, Population]:
    """The Table 1 sign tests in all studies and the disease-vs-healthy oriented subset,
    which Methods Plate II and Plates III and VII quote (#579)."""
    return populations(read_tsv(repo / COMPARISON))


def increased_side(population: Population) -> list[str]:
    """The Table 1 traits whose share is larger in the increased signature."""
    return [trait for trait, test in population.fractions.items() if test.share_up > 0.5]


def survivor_codes(population: Population) -> str:
    """'FA/CEL, CEL/SPH, CEL/ARS': the log-ratios at BY q <= 0.05, or 'none'."""
    return ", ".join(test.label for test in population.ratio_survivors()) or "none"


def passing_pairs(population: Population, pairs: Sequence[tuple[str, str]]) -> int:
    """How many of the given log-ratios pass BY q <= 0.05."""
    return sum(population.ratio(*pair).passes() for pair in pairs)


def increased_pairs_text(
    every: Population, oriented: Population, pairs: Sequence[tuple[str, str]]
) -> str:
    """Plate III's label for a group of pairs: those passing BY in either cohort, or the
    lowest BY q among them (#579, #686)."""
    if not pairs:
        return "no pairs"
    passing = sorted(
        {
            population.ratio(*pair).label
            for population in (every, oriented)
            for pair in pairs
            if population.ratio(*pair).passes()
        }
    )
    if passing:
        verb = "passes" if len(passing) == 1 else "pass"
        return f"{', '.join(passing)} {verb} BY in all or oriented"
    lowest = min(
        population.ratio(*pair).q_value for population in (every, oriented) for pair in pairs
    )
    return f"all BY q ≥ {math.floor(lowest * 100) / 100:.2f} in both all & oriented"


def oriented_cellulose_pairs_text(
    every: Population, oriented: Population, pairs: Sequence[tuple[str, str]]
) -> str:
    """Plate III's line on what the oriented subset changes among the cellulose pairs (#579)."""
    before = {every.ratio(*pair).label for pair in pairs if every.ratio(*pair).passes()}
    after = [oriented.ratio(*pair) for pair in pairs]
    passing = [test for test in after if test.passes()]
    added = [test.label for test in passing if test.label not in before]
    lost = sorted(before - {test.label for test in passing})
    verb = "passes" if len(passing) == 1 else "pass"
    counted = f"{len(passing)} of {len(pairs)} Cellulose pairs {verb} BY"
    if added or lost:
        changes = [f"adds {', '.join(added)}"] if added else []
        changes += [f"loses {', '.join(lost)}"] if lost else []
        return f"Oriented subset {' and '.join(changes)} ({counted})"
    failing = [test for test in after if not test.passes()]
    closest = min(failing, key=lambda test: test.q_value) if failing else None
    nearest = f"; next {closest.label} at BY q = {fmt_p(closest.q_value)}" if closest else ""
    # Short enough for the panel's title line under any classifier (#704)
    return f"Oriented subset: the same {counted}{nearest}"


def study_bh_count(rows: dict[str, dict[str, str]]) -> int:
    """Traits passing study-level BH q ≤ 0.05 in one stratum's paired-contrast table."""
    return sum(
        1 for r in rows.values() if r.get("study_q_value") and float(r["study_q_value"]) <= 0.05
    )


def check_plate_two_rows(trait_ids: Sequence[str], pooled: dict[str, dict[str, str]]) -> None:
    """Every Plate II row must be in the pooled table, whose counts the plate and slide quote.

    A row the stratified screen did not test elsewhere is drawn as "n.t." there (#616).
    """
    missing = [trait for trait in trait_ids if trait not in pooled]
    if missing:
        raise ValueError(f"Plate II rows missing from the pooled table: {missing}")


def check_plate_two_directions(
    rows: Sequence[tuple[str, str, str]], pooled: Mapping[str, Mapping[str, str]]
) -> None:
    """Plate II groups its rows by direction and tags the decreased groups "↓", so a row's
    pooled study contrast must have its group's sign. A trait that changed sides in a rerun
    stops the plate rather than sitting under the wrong tag (#579)."""
    wrong = []
    for trait, _, group in rows:
        value = pooled[trait].get("study_mean_increased_minus_decreased")
        if value and (float(value) < 0) != group.endswith("↓"):
            wrong.append(f"{trait} in {group} at {float(value):+.3f}")
    if wrong:
        raise ValueError(f"Plate II rows whose pooled contrast has the other sign: {wrong}")


def check_classifiers(repo_root: Path, input_root: Path) -> None:
    """Refuse to draw plates whose cellulose_decr notes and Table 1 text define disease vs
    healthy with different contrast classifiers (#679)."""
    provenance = repo_root / CELLULOSE_DIR / "cellulose_decr_provenance.json"
    cellulose = json.loads(provenance.read_text(encoding="utf-8")).get(
        "contrast_classifier_version"
    )
    contrasts = input_root / CONTRAST_PROVENANCE
    table1 = json.loads(contrasts.read_text(encoding="utf-8"))["classifier_version"]
    if cellulose != table1:
        raise SystemExit(
            f"the cellulose_decr tables use {cellulose} and the Table 1 tables {table1}; re-pin "
            "cellulose_decr (just cellulose-decr) or rerun just contrast-orientation (#679)"
        )


def input_record(repo: Path) -> dict[str, str]:
    """The files the last plate read, relative to the root they came from, with hashes.

    The cellulose_decr files come from this checkout and the rest from the input root,
    so the deck can check each against the revision it cites (#637).
    """
    record = {}
    for path in sorted(INPUTS_READ):
        root = REPO_ROOT if path.is_relative_to(REPO_ROOT / CELLULOSE_DIR) else repo
        record[path.relative_to(root).as_posix()] = file_sha256(path)
    return record


def write_plates_provenance(
    out_dir: Path,
    inputs: Mapping[str, Mapping[str, str]],
    code_commit: str,
    outputs: Mapping[str, str] | None = None,
    counts: Mapping[str, int] | None = None,
) -> Path:
    """Record the commit that drew the plates, each plate's inputs and outputs with their
    hashes, and the counts the slide text quotes.

    The deck links the plate code at that commit (#613), checks the inputs against the
    revisions its slides cite, and checks that the plates it places are the ones drawn
    (#637).
    """
    path = out_dir / "plates_provenance.json"
    record = {
        "script": "scripts/generate_tufte_plates.py",
        "code_commit": code_commit,
        "rerun_commit": RERUN_COMMIT,
        "plates": sorted(inputs),
        "inputs": {plate: dict(sorted(files.items())) for plate, files in sorted(inputs.items())},
        "outputs": dict(sorted((outputs or {}).items())),
        "counts": dict(sorted((counts or {}).items())),
    }
    path.write_text(json.dumps(record, indent=1) + "\n", encoding="utf-8")
    return path


def trait_labels(repo: Path) -> dict[str, str]:
    """Each trait column's label as the plates show it (#556, #591, #607).

    The label is the column's object_name. A compound the taxon uses as a carbon source
    adds "(carbon source)", since its object_name alone ("raffinose") names no trait.
    """
    labels = {}
    for row in read_tsv(repo / TRAIT_COLUMNS):
        name = row["object_name"].replace("_", " ")
        if row.get("predicate_name") == "uses as carbon source":
            name = f"{name} (carbon source)"
        labels[row["column_id"]] = name
    return labels


def canonical_display(display: str, label: str | None) -> str:
    """A row label with its trait name replaced by the canonical one; tags and marks kept."""
    if label is None:
        return display
    match = re.match(r"^(\[[^\]]*\]\s+)?(.*?)((?:\s+\[[^\]]*\])?(?:\s+★)?)$", display)
    if match is None:
        return display
    prefix, _, suffix = match.groups()
    return f"{prefix or ''}{label[:1].upper()}{label[1:]}{suffix}"


def cellulose_decr_notes(repo: Path) -> dict[str, str]:
    """Plate annotations for the cellulose_decr results, read from its committed tables.

    The plates used to carry numbers copied from a prompt, which #556 found to be wrong;
    reading them keeps the plates in step with the tables. The statement of which results
    depend on the condition definition comes from ``just cellulose-decr-figures``.
    """
    captions = repo / CELLULOSE_DIR / "figures" / "figure_captions.json"
    if not captions.exists():
        raise FileNotFoundError(f"missing {captions}; run `just cellulose-decr-figures` first")
    fragile = read_json(captions)["statements"]["fragile"]
    rows = read_tsv(repo / CELLULOSE_DIR / "compartment_direction_interaction.tsv")

    def crc(population: str, unit: str) -> dict[str, str]:
        return next(
            row
            for row in rows
            if row["contrast_population"] == population
            and row["short_condition"] == "CRC"
            and row["unit"] == unit
            and row["compartment_b"] == "Non-fecal GI (mucosa/lumen)"
        )

    stool = crc("all_contrasts", "experiment")
    healthy = crc("disease_vs_healthy_oriented", "study")
    summary = read_json(repo / CELLULOSE_DIR / "cellulose_decr_summary.json")
    observed = summary["observed"]
    u6 = summary["null_models"]["U6"]
    counts = f"stool {stool['a_down']}↓:{stool['a_up']}↑ vs GI {stool['b_down']}↓:{stool['b_up']}↑"
    odds = f"OR {float(stool['conditional_or']):.1f}, p={fmt_p(float(stool['fisher_p_two_sided']))}"
    healthy_p = f"disease-vs-healthy studies p={fmt_p(float(healthy['fisher_p_two_sided']))}"
    bullets = [
        f"CRC cellulose, paired experiments: {counts} · {odds} · {healthy_p}",
        f"Stool overlap: {observed['three_way_core']} of {observed['three_way_union']} taxa in "
        "all three · whole-signature null "
        f"p={fmt_p(float(u6['mean_pairwise_jaccard_distance']['p']))} (mean d_J), "
        f"{fmt_p(float(u6['three_way_core']['p']))} (core)",
    ]
    if fragile:
        # Plate IV stars HIV cellulose from the profile-level Fisher test (#590).
        bullets.append(f"Cellulose deep-dive: {fragile}")
    return {
        # Wrapped to stay inside Plate IV's panel B.
        "plate4_callout": "\n".join(
            textwrap.fill(f"• {bullet}", width=88, subsequent_indent="  ") for bullet in bullets
        ),
        "crc_source": f"CRC cellulose stool vs GI: {odds}; {healthy_p} (cellulose_decr).",
    }


def setup_tufte_rc() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
            "axes.edgecolor": MID_RULE,
            "axes.linewidth": 0.8,
            "xtick.color": INK,
            "ytick.color": INK,
            "text.color": INK,
            "figure.facecolor": BG_COLOR,
            "axes.facecolor": PANEL_BG,
            "savefig.facecolor": BG_COLOR,
            "savefig.dpi": 220,
        }
    )


def despine(
    ax: plt.Axes, top: bool = True, right: bool = True, left: bool = False, bottom: bool = False
) -> None:
    if top:
        ax.spines["top"].set_visible(False)
    if right:
        ax.spines["right"].set_visible(False)
    if left:
        ax.spines["left"].set_visible(False)
    if bottom:
        ax.spines["bottom"].set_visible(False)


# ==============================================================================
# METHODS PLATE 1: High-Level Data & Taxonomic Integration Architecture
# ==============================================================================
def guild_text(labels: dict[str, str]) -> str:
    """The aerobe and anaerobe guild terms of mtp.trait_guilds, by canonical label."""
    lines = []
    for guild, terms in OXYGEN_GUILDS.items():
        first, *rest = terms
        others = ", ".join(f"{labels[term]} ({term.removeprefix('METPO:')})" for term in rest)
        lines.append(
            f"• {guild.capitalize()} Guild ({len(terms)} terms): {labels[first]} ({first}),\n"
            f"  {others}"
        )
    return "\n".join(lines)


def build_methods_plate1(repo: Path, out_dir: Path) -> Path:
    fig = plt.figure(figsize=(16.0, 9.0), facecolor=BG_COLOR)
    ax = fig.add_axes([0.03, 0.08, 0.94, 0.80])
    ax.axis("off")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)

    fig.suptitle(
        "METHODS OVERVIEW 1 · BugSigDB \N{MULTIPLICATION SIGN} KG-Microbe Data Architecture & "
        "Taxonomic Resolution",
        fontsize=15.6,
        fontweight="bold",
        x=0.03,
        ha="left",
        y=0.955,
        color=INK,
    )
    fig.text(
        0.03,
        0.912,
        "How curated human microbiome abundance signatures map onto organism-level microbiological "
        "trait annotations via NCBI Taxonomy subtree rollup",
        fontsize=9.8,
        color=MUTED,
    )

    def draw_card(
        x: float,
        y: float,
        w: float,
        h: float,
        title: str,
        header_bg: str,
        edge: str = MID_RULE,
        lw: float = 1.3,
    ) -> None:
        body = mpatches.FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle="round,pad=0.6,rounding_size=1.4",
            facecolor=PANEL_BG,
            edgecolor=edge,
            linewidth=lw,
            zorder=1,
        )
        ax.add_patch(body)
        hdr = mpatches.FancyBboxPatch(
            (x, y + h - 7.5),
            w,
            7.5,
            boxstyle="round,pad=0.6,rounding_size=1.4",
            facecolor=header_bg,
            edgecolor=edge,
            linewidth=lw,
            zorder=2,
        )
        ax.add_patch(hdr)
        ax.text(
            x + w / 2,
            y + h - 3.7,
            title,
            fontsize=10.6,
            fontweight="bold",
            color="#ffffff",
            ha="center",
            va="center",
            zorder=3,
        )

    # --- Column 1: BugSigDB Corpus Funnel ---
    draw_card(
        1.5,
        6.0,
        28.5,
        88.0,
        "1. BugSigDB Literature Corpus (v1.3.1)",
        DEC_COLOR,
        edge=DEC_COLOR,
        lw=1.5,
    )

    paired_dir = repo / "trait-visualization/trait-profile-investigation"
    pooled = read_tsv(paired_dir / "trait_profile_paired_contrast_trait_stats_pooled.tsv")
    feces = read_tsv(paired_dir / "trait_profile_paired_contrast_trait_stats_feces.tsv")
    feces_excluded = read_tsv(
        paired_dir / "trait_profile_paired_contrast_trait_stats_feces_excluded.tsv"
    )
    pairs = max(int(r["pair_count"]) for r in pooled)
    studies = max(int(float(r["study_count"])) for r in pooled)
    feces_studies = max(int(float(r["study_count"])) for r in feces)
    non_feces_studies = max(int(float(r["study_count"])) for r in feces_excluded)
    ind_cells = [
        r
        for r in read_tsv(
            paired_dir
            / "six-shared-traits"
            / "six_shared_trait_condition_body_site_direction_contrasts.tsv"
        )
        if r.get("view_key") == "individual_conditions_direction_contrast"
    ]
    cond_count = len({r["condition"] for r in ind_cells})
    site_count = len({r["body_site_group"] for r in ind_cells})

    funnel_steps = [
        (
            74.0,
            "1,039 Curated Studies",
            "Published human microbiome case\N{EN DASH}control studies",
            "#eef5fb",
            DEC_COLOR,
        ),
        (
            56.0,
            f"{pairs:,} Paired Experiments",
            f"{studies:,} unique studies with both INC & DEC taxon sets\n"
            f"{feces_studies:,} Feces studies  ·  {non_feces_studies:,} Non-Feces studies",
            "#e1effa",
            DEC_COLOR,
        ),
        (
            37.0,
            f"{2 * pairs:,} Directional Signatures",
            f"{pairs:,} Increased (↑) + {pairs:,} Decreased (↓) profiles\n"
            f"Spanning {cond_count:,} conditions \N{MULTIPLICATION SIGN} {site_count:,} anatomical "
            "sites",
            "#d4e6f4",
            DEC_COLOR,
        ),
        (
            18.0,
            "48,370 Taxon Occurrences",
            "6,356 unique NCBI TaxIDs across all ranks\nMedian signature size: 7 taxa (IQR "
            "4\N{EN DASH}13)",
            "#f5d7d4",
            INC_COLOR,
        ),
    ]
    for idx_f, (fy, ftitle, fsub, fbg, fcol) in enumerate(funnel_steps):
        box = mpatches.FancyBboxPatch(
            (3.5, fy),
            24.5,
            11.8,
            boxstyle="round,pad=0.4,rounding_size=1.0",
            facecolor=fbg,
            edgecolor=fcol,
            linewidth=1.2,
            zorder=3,
        )
        ax.add_patch(box)
        ax.text(
            15.75,
            fy + 7.8,
            ftitle,
            fontsize=10.2,
            fontweight="bold",
            color=INK,
            ha="center",
            va="center",
            zorder=4,
        )
        ax.text(
            15.75, fy + 3.2, fsub, fontsize=8.0, color=MUTED, ha="center", va="center", zorder=4
        )
        if idx_f < len(funnel_steps) - 1:
            ax.annotate(
                "",
                xy=(15.75, funnel_steps[idx_f + 1][0] + 12.4),
                xytext=(15.75, fy - 0.4),
                arrowprops={
                    "arrowstyle": "-|>,head_width=0.6,head_length=0.6",
                    "color": DEC_COLOR,
                    "lw": 1.8,
                },
                zorder=5,
            )

    ax.text(
        15.75,
        10.5,
        "Reported abundance shifts (not cell density)",
        fontsize=8.0,
        style="italic",
        color=MUTED,
        ha="center",
    )

    # Arrow Column 1 -> Column 2
    ax.annotate(
        "",
        xy=(34.2, 50.0),
        xytext=(30.6, 50.0),
        arrowprops={"arrowstyle": "-|>,head_width=0.8,head_length=0.8", "color": INK, "lw": 2.2},
        zorder=6,
    )

    # --- Column 2: KG-Microbe Trait Graph & Subtree Rollup ---
    draw_card(
        34.8,
        6.0,
        30.4,
        88.0,
        "2. KG-Microbe Trait Graph (20260120)",
        TEAL_ACCENT,
        edge=TEAL_ACCENT,
        lw=1.5,
    )

    kg_top = mpatches.FancyBboxPatch(
        (36.8, 71.5),
        26.4,
        13.5,
        boxstyle="round,pad=0.4,rounding_size=1.0",
        facecolor="#e8f5f3",
        edgecolor=TEAL_ACCENT,
        linewidth=1.2,
        zorder=3,
    )
    ax.add_patch(kg_top)
    ax.text(
        50.0,
        80.8,
        "5,829 Trait-Bearing Taxa \N{MULTIPLICATION SIGN} 317 Traits",
        fontsize=10.1,
        fontweight="bold",
        color=INK,
        ha="center",
        va="center",
        zorder=4,
    )
    ax.text(
        50.0,
        75.2,
        "Filtered from 2,650 raw features by excluding\n2,333 MediaDive growth-media recipes",
        fontsize=8.1,
        color=MUTED,
        ha="center",
        va="center",
        zorder=4,
    )

    # 4 Ontology pills in 2x2 grid
    ont_specs = [
        (
            36.8,
            56.5,
            "METPO Phenotypes",
            "O₂ tolerance, cell shape,\nGram stain, halophily",
            "#fef9ec",
            GOLD_ACCENT,
        ),
        (
            50.4,
            56.5,
            "GO Biological Process",
            "Cellulose, xylan, lignin &\nurea catabolism",
            "#eef5fb",
            DEC_COLOR,
        ),
        (
            36.8,
            42.5,
            "ChEBI Substrates",
            "Carbon, sugar, organic acid\n& nucleoside utilization",
            "#f3eef8",
            PURPLE_ACCENT,
        ),
        (
            50.4,
            42.5,
            "Metabolic Pathways",
            "Arsenate, nitrate, sulfate &\nthiosulfate redox",
            "#f5d7d4",
            INC_COLOR,
        ),
    ]
    for ox, oy, otitle, osub, obg, ocol in ont_specs:
        op = mpatches.FancyBboxPatch(
            (ox, oy),
            12.8,
            11.2,
            boxstyle="round,pad=0.3,rounding_size=0.8",
            facecolor=obg,
            edgecolor=ocol,
            linewidth=1.1,
            zorder=3,
        )
        ax.add_patch(op)
        ax.text(
            ox + 6.4,
            oy + 7.8,
            otitle,
            fontsize=8.3,
            fontweight="bold",
            color=INK,
            ha="center",
            va="center",
            zorder=4,
        )
        ax.text(
            ox + 6.4, oy + 3.2, osub, fontsize=7.3, color=MUTED, ha="center", va="center", zorder=4
        )

    # Subtree Rollup Box
    rollup_box = mpatches.FancyBboxPatch(
        (36.8, 12.0),
        26.4,
        26.5,
        boxstyle="round,pad=0.4,rounding_size=1.0",
        facecolor="#f7f9fa",
        edgecolor=INK,
        linewidth=1.2,
        zorder=3,
    )
    ax.add_patch(rollup_box)
    ax.text(
        50.0,
        34.8,
        "Taxonomic Subtree Trait Integration",
        fontsize=9.6,
        fontweight="bold",
        color=INK,
        ha="center",
        va="center",
        zorder=4,
    )
    ax.text(
        50.0,
        28.2,
        "Species/Strain traits propagate up NCBI Taxonomy\nto Genus, Family, Order, Class & Phylum "
        "nodes",
        fontsize=8.0,
        color=MUTED,
        ha="center",
        va="center",
        zorder=4,
    )
    ax.text(
        50.0,
        19.8,
        "• Arm 1 & Arm 3: Continuous / Sign profile matrix\n  weighted by descendant trait "
        "assertion fractions\n• Arm 2 (Taxon ORA): Binary presence/absence\n  over M = 5,829 "
        "distinct organisms",
        fontsize=7.8,
        color=INK,
        ha="center",
        va="center",
        zorder=4,
    )

    # Arrow Column 2 -> Column 3
    ax.annotate(
        "",
        xy=(69.4, 50.0),
        xytext=(65.8, 50.0),
        arrowprops={"arrowstyle": "-|>,head_width=0.8,head_length=0.8", "color": INK, "lw": 2.2},
        zorder=6,
    )

    # --- Column 3: Coverage by Rank & Guild Definitions ---
    draw_card(
        70.0,
        6.0,
        28.5,
        88.0,
        "3. Coverage & Functional Guilds",
        GOLD_ACCENT,
        edge=GOLD_ACCENT,
        lw=1.5,
    )

    ax.text(
        84.25,
        82.5,
        "Signature Occurrence Coverage by Rank",
        fontsize=9.5,
        fontweight="bold",
        color=INK,
        ha="center",
    )
    # Read from the table the plate cites, which Plate VI also draws (#640).
    coverage = {
        r["taxonomic_rank"]: r
        for r in read_tsv(
            repo / "trait-visualization/trait-profile-investigation/taxa-coverage/"
            "taxa_coverage_by_rank.tsv"
        )
    }
    rank_bars = [
        (
            rank.capitalize(),
            float(coverage[rank]["trait_bearing_occurrences_pct"]),
            f"{int(coverage[rank]['trait_bearing_occurrences']):,} / "
            f"{int(coverage[rank]['total_occurrences']):,} occ.",
            INC_COLOR if rank == "species" else TEAL_ACCENT,
        )
        for rank in ("phylum", "order", "genus", "family", "species")
    ]
    species_gap = float(coverage["species"]["zero_trait_occurrences_pct"])
    for idx_b, (rname, rpct, rlbl, rcol) in enumerate(rank_bars):
        by = 76.5 - idx_b * 5.2
        ax.text(71.8, by, rname, fontsize=8.2, fontweight="bold", color=INK, ha="left", va="center")
        # bar track
        track = mpatches.Rectangle(
            (77.8, by - 1.5), 12.0, 3.0, facecolor="#eaeef2", edgecolor="none", zorder=2
        )
        ax.add_patch(track)
        fill = mpatches.Rectangle(
            (77.8, by - 1.5),
            12.0 * (rpct / 100.0),
            3.0,
            facecolor=rcol,
            alpha=0.88,
            edgecolor="none",
            zorder=3,
        )
        ax.add_patch(fill)
        ax.text(
            90.3,
            by,
            f"{rpct:.1f}%",
            fontsize=8.0,
            fontweight="bold",
            color=rcol,
            ha="left",
            va="center",
        )
        ax.text(
            97.5,
            by,
            f"({rlbl.split('/')[0].strip()})",
            fontsize=7.1,
            color=MUTED,
            ha="right",
            va="center",
        )

    ax.text(
        84.25,
        50.8,
        f"Species gap ({species_gap:.0f}% of occurrences): uncultivated MAGs, slash-bins &\n"
        "2 NCBI TaxID desyncs (Oscillospiraceae, Rothia)",
        fontsize=7.6,
        color=INC_COLOR,
        ha="center",
        va="center",
    )

    # Guild Definition Box (replaces the standalone guild table slides!)
    guild_box = mpatches.FancyBboxPatch(
        (71.8, 28.5),
        24.9,
        18.5,
        boxstyle="round,pad=0.4,rounding_size=0.9",
        facecolor="#fef9ec",
        edgecolor=GOLD_ACCENT,
        linewidth=1.1,
        zorder=3,
    )
    ax.add_patch(guild_box)
    ax.text(
        84.25,
        44.2,
        "OR-Rule Functional Guilds (1 if ≥1 member term)",
        fontsize=8.6,
        fontweight="bold",
        color=INK,
        ha="center",
    )
    ax.text(
        72.6,
        36.5,
        guild_text(trait_labels(repo))
        + "\n• Gram Envelope: Gram+ (1000698) vs Gram\N{MINUS SIGN} (1000699)",
        fontsize=7.4,
        color=INK,
        ha="left",
        va="center",
    )

    # Canonical Table 1 Box
    t1_box = mpatches.FancyBboxPatch(
        (71.8, 9.5),
        24.9,
        16.5,
        boxstyle="round,pad=0.4,rounding_size=0.9",
        facecolor="#f4f7f9",
        edgecolor=INK,
        linewidth=1.1,
        zorder=3,
    )
    ax.add_patch(t1_box)
    ax.text(
        84.25,
        23.3,
        "6 Canonical Table 1 Focal Traits",
        fontsize=8.8,
        fontweight="bold",
        color=INK,
        ha="center",
    )
    ax.text(
        72.6,
        15.8,
        "• 5 Increased (↑): Facultatively aerobic, Moderately\n  halophilic, Arsenate reduction, "
        "Coccus, Sphere shaped\n• 1 Decreased (↓): Cellulose catabolic process\n  (plus Ecological "
        "Composite = meanInc(5) \N{MINUS SIGN} Cellulose)",
        fontsize=7.5,
        color=INK,
        ha="left",
        va="center",
    )

    fig.text(
        0.03,
        0.025,
        f"{RERUN_SOURCE} · BugSigDB v1.3.1 (Zenodo 10627578) "
        "\N{MULTIPLICATION SIGN} KG-Microbe 20260120 · taxa_coverage_by_rank.tsv, "
        "condition_trait_module_guild_definition_table.tsv.",
        fontsize=8.0,
        color=MUTED,
    )

    out_path = out_dir / "methods_plate1_data_architecture_funnel.png"
    fig.savefig(out_path, dpi=220)
    fig.savefig(out_path.with_suffix(".svg"))
    plt.close(fig)
    return out_path


# ==============================================================================
# METHODS PLATE 2: High-Level Three-Arm Analytical Architecture & Convergence
# ==============================================================================
def build_methods_plate2(repo: Path, out_dir: Path) -> Path:
    # Arm 1's survivors, read from the tables Plates I and II draw (#598, #623)
    stats = repo / "trait-visualization/trait-profile-investigation"
    paired = {
        stratum: {
            r["trait"]: r
            for r in read_tsv(stats / f"trait_profile_paired_contrast_trait_stats_{stratum}.tsv")
        }
        for stratum in ("pooled", "feces")
    }
    contrasts = read_tsv(
        stats / "six-shared-traits/six_shared_trait_condition_body_site_direction_contrasts.tsv"
    )
    fisher = sum(
        r.get("selected_after_bh") == "True"
        for r in contrasts
        if r.get("view_key") == "individual_conditions_direction_contrast"
    )
    survivors = (
        f"{study_bh_count(paired['pooled'])} Pooled & {study_bh_count(paired['feces'])} Feces "
        f"traits\npass Study BH q ≤ 0.05 · {fisher} Cond"
        "\N{MULTIPLICATION SIGN}Site Fisher BH survivors"
    )
    # The test families of Arm 1's stratified screens and Arm 3's 317-trait screen (#624, #640)
    views = [r.get("view_key") for r in contrasts]
    individual = views.count("individual_conditions_direction_contrast")
    grouped = views.count("grouped_direction_contrast")
    wilcoxon = len(
        read_tsv(stats / "PD-cellulose/pd_condition_taxon_trait_network_condition_trait_tests.tsv")
    )
    screen = read_tsv(repo / "data/derived/trait_screen/trait_direction_screen.tsv")
    # Arm 2's tests and calls and the oxygen-guild McNemar test, read from their tables (#579)
    enriched, depleted = ora_calls(repo / ORA_ENRICHMENT), ora_calls(repo / ORA_DEPLETION)
    enrichment_rows = read_tsv(repo / ORA_ENRICHMENT)
    by_side = {side: sum(r["inc_dec"] == side for r in enrichment_rows) for side in ("inc", "dec")}
    enriched_traits = {r["trait"] for r in enrichment_rows}
    mcnemar = McNemar.from_rows(read_tsv(repo / GUILD_CORRELATION))
    concordant, strict = three_arm_counts(
        three_arm_traits(
            paired["pooled"],
            {
                r["trait"]: r
                for r in read_tsv(repo / "data/derived/trait_enrichment_traits_overall_ranked.tsv")
            },
            {r["trait"]: r for r in screen},
            {tid for tid, *_ in TABLE1_TRAITS},
        )
    )
    # The slide's caveat quotes these, from the plates' provenance (#640)
    PLATE_COUNTS.update(
        {
            "traits_observed": len(paired["pooled"]),
            "three_arm_concordant": concordant,
            "three_arm_strict_by": strict,
        }
    )

    def passing(column: str) -> int:
        return sum(has_value(r.get(column)) and float(r[column]) <= 0.05 for r in screen)

    family = str(len(screen)).translate(str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉"))
    signs = table1_signs(repo)
    every, oriented = signs[ALL_STUDIES], signs[ORIENTED]
    arm3 = (
        f"Arm 3 Survivors: {passing('bh_q_nominal')} traits pass BH{family} ≤ 0.05 "
        f"({passing('by_q_nominal')} BY)"
    )
    fig = plt.figure(figsize=(16.0, 9.0), facecolor=BG_COLOR)
    ax = fig.add_axes([0.03, 0.08, 0.94, 0.80])
    ax.axis("off")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)

    fig.suptitle(
        "METHODS OVERVIEW 2 · The Three Independent Analytical Arms & Cross-Arm Convergence Funnel",
        fontsize=15.6,
        fontweight="bold",
        x=0.03,
        ha="left",
        y=0.955,
        color=INK,
    )
    fig.text(
        0.03,
        0.912,
        "Triangulating trait shifts across (1) Continuous Paired Contrasts, (2) Taxon-Calibrated "
        "Hypergeometric ORA, and (3) Non-Parametric Sign & Log-Ratio Tests",
        fontsize=9.8,
        color=MUTED,
    )

    def draw_arm_card(
        x: float, y: float, w: float, h: float, title: str, subtitle: str, hdr_col: str
    ) -> None:
        card = mpatches.FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle="round,pad=0.5,rounding_size=1.4",
            facecolor=PANEL_BG,
            edgecolor=hdr_col,
            linewidth=1.6,
            zorder=1,
        )
        ax.add_patch(card)
        hdr = mpatches.FancyBboxPatch(
            (x, y + h - 10.5),
            w,
            10.5,
            boxstyle="round,pad=0.5,rounding_size=1.4",
            facecolor=hdr_col,
            edgecolor=hdr_col,
            linewidth=1.6,
            zorder=2,
        )
        ax.add_patch(hdr)
        ax.text(
            x + w / 2,
            y + h - 3.6,
            title,
            fontsize=10.4,
            fontweight="bold",
            color="#ffffff",
            ha="center",
            va="center",
            zorder=3,
        )
        ax.text(
            x + w / 2,
            y + h - 7.8,
            subtitle,
            fontsize=8.0,
            color="#f0f4f8",
            ha="center",
            va="center",
            zorder=3,
        )

    # --- ARM 1 (Left Column) ---
    draw_arm_card(
        1.5,
        25.0,
        30.5,
        69.0,
        "ARM 1 · Continuous Paired Contrasts",
        "Row/Col Standardized z-Scores & Stratified Contrasts",
        DEC_COLOR,
    )

    a1_blocks = [
        (
            68.5,
            "1. Dual Normalization (4,670 \N{MULTIPLICATION SIGN} 317)",
            "log1p(assertion count) → profile row z-score\n→ trait column z-score across corpus",
            "#eef5fb",
        ),
        (
            53.5,
            "2. Within-Study Paired Shift (Δz)",
            "Δz = z_increased \N{MINUS SIGN} z_decreased per experiment\nAggregated across N = 959 "
            "unique studies",
            "#e1effa",
        ),
        (
            38.5,
            "3. Condition \N{MULTIPLICATION SIGN} Site Stratification",
            f"• Positive-z share Fisher ({individual:,} ind. / {grouped:,} grouped cells)\n"
            f"• Fecal Study Wilcoxon signed-rank ({wilcoxon:,} tests)",
            "#d4e6f4",
        ),
    ]
    for by, btitle, bsub, bbg in a1_blocks:
        b = mpatches.FancyBboxPatch(
            (3.3, by),
            26.9,
            11.8,
            boxstyle="round,pad=0.35,rounding_size=0.9",
            facecolor=bbg,
            edgecolor=DEC_COLOR,
            lw=1.0,
            zorder=3,
        )
        ax.add_patch(b)
        ax.text(
            16.75,
            by + 8.0,
            btitle,
            fontsize=9.0,
            fontweight="bold",
            color=INK,
            ha="center",
            va="center",
            zorder=4,
        )
        ax.text(
            16.75, by + 3.2, bsub, fontsize=7.7, color=MUTED, ha="center", va="center", zorder=4
        )

    ax.text(
        16.75,
        30.5,
        f"Arm 1 Survivors: {survivors}",
        fontsize=8.2,
        fontweight="bold",
        color=DEC_COLOR,
        ha="center",
        va="center",
    )

    # --- ARM 2 (Center Column) ---
    draw_arm_card(
        34.75,
        25.0,
        30.5,
        69.0,
        "ARM 2 · Taxon-Calibrated ORA",
        "2\N{MULTIPLICATION SIGN}2 Hypergeometric Over-Representation (M = 5,829 Taxa)",
        TEAL_ACCENT,
    )

    a2_top = mpatches.FancyBboxPatch(
        (36.55, 67.5),
        26.9,
        12.8,
        boxstyle="round,pad=0.35,rounding_size=0.9",
        facecolor="#e8f5f3",
        edgecolor=TEAL_ACCENT,
        lw=1.0,
        zorder=3,
    )
    ax.add_patch(a2_top)
    ax.text(
        50.0,
        76.2,
        "Distinct-Taxon Hypergeometric Urn",
        fontsize=9.1,
        fontweight="bold",
        color=INK,
        ha="center",
        va="center",
        zorder=4,
    )
    ax.text(
        50.0,
        70.8,
        f"{enriched.tests:,} profile \N{MULTIPLICATION SIGN} trait tests (BH q ≤ "
        f"{enriched.p_threshold})\nEach taxon counts once (prevents genus hub inflation)",
        fontsize=7.7,
        color=MUTED,
        ha="center",
        va="center",
        zorder=4,
    )

    # Visual comparison of Legacy vs Taxon-Calibrated ORA (absorbs old figure2!)
    ora_box = mpatches.FancyBboxPatch(
        (36.55, 38.5),
        26.9,
        25.8,
        boxstyle="round,pad=0.35,rounding_size=0.9",
        facecolor="#f7f9fa",
        edgecolor=INK,
        lw=1.1,
        zorder=3,
    )
    ax.add_patch(ora_box)
    ax.text(
        50.0,
        61.2,
        "Calibration Eliminates False Depletions",
        fontsize=8.8,
        fontweight="bold",
        color=INK,
        ha="center",
        zorder=5,
    )

    # Mini bars inside ora_box (zorder=5 so they sit on top of ora_box!)
    ax.text(
        37.6,
        55.6,
        "Legacy Assertion Rollup:",
        fontsize=7.6,
        fontweight="bold",
        color=MUTED,
        ha="left",
        zorder=5,
    )
    ax.add_patch(
        mpatches.Rectangle((37.6, 51.4), 8.2, 2.7, facecolor=TEAL_ACCENT, alpha=0.70, zorder=5)
    )
    ax.text(46.3, 52.75, "63,292 Enrich.", fontsize=7.3, color=TEAL_ACCENT, va="center", zorder=6)
    ax.add_patch(
        mpatches.Rectangle((37.6, 47.6), 12.6, 2.7, facecolor=INC_COLOR, alpha=0.88, zorder=5)
    )
    ax.text(
        50.7,
        48.95,
        "100,999 False Depl.",
        fontsize=7.3,
        fontweight="bold",
        color=INC_COLOR,
        va="center",
        zorder=6,
    )

    ax.text(
        37.6,
        44.0,
        "Taxon-Calibrated ORA:",
        fontsize=7.6,
        fontweight="bold",
        color=INK,
        ha="left",
        zorder=5,
    )
    ax.add_patch(
        mpatches.Rectangle((37.6, 40.0), 8.8, 2.7, facecolor=TEAL_ACCENT, alpha=0.95, zorder=5)
    )
    ax.text(
        46.9,
        41.35,
        f"{enriched.calls:,} Enrich. vs {depleted.calls:,} Depl.!",
        fontsize=7.3,
        fontweight="bold",
        color=TEAL_ACCENT,
        va="center",
        zorder=6,
    )

    ax.text(
        50.0,
        30.5,
        f"Arm 2 Survivors: {by_side['inc']:,} Inc + {by_side['dec']:,} Dec\nenrichments across "
        f"{len(enriched_traits)} traits (fold ≥ {enriched.fold})",
        fontsize=8.2,
        fontweight="bold",
        color=TEAL_ACCENT,
        ha="center",
        va="center",
    )

    # --- ARM 3 (Right Column) ---
    draw_arm_card(
        68.0,
        25.0,
        30.5,
        69.0,
        "ARM 3 · Non-Parametric Sign & Log-Ratio",
        "Exact Binomial, McNemar & Compositional Closure",
        PURPLE_ACCENT,
    )

    a3_blocks = [
        (
            68.5,
            "1. Study Exact Binomial Sign Test",
            f"Majority vote per study ({every.tested_studies:,} all / "
            f"{oriented.tested_studies:,} oriented)\n"
            f"{len(screen)}-trait BH/BY screen (fixes {screening_bias(repo)} screening bias)",
            "#f3eef8",
        ),
        (
            53.5,
            "2. Coupled Oxygen-Guild McNemar Test",
            f"{mcnemar.discordant:,} discordant studies: {mcnemar.aerobe_up_anaerobe_down} "
            f"Aerobe↑/Anaerobe↓\nvs {mcnemar.anaerobe_up_aerobe_down} opposite (OR = "
            f"{mcnemar.odds_ratio:.2f}, p = {sci(mcnemar.p_value)})",
            "#fef9ec",
        ),
        (
            38.5,
            "3. Reference-Explicit Log-Ratios",
            "Δlog₂(Trait A / Trait B) across 15 Table 1 pairs\n"
            f"{len(every.ratio_survivors())} BY hits (all {every.studies:,}) → "
            f"{len(oriented.ratio_survivors())} BY in oriented ({oriented.studies:,} st.)",
            "#f5d7d4",
        ),
    ]
    for by, btitle, bsub, bbg in a3_blocks:
        b = mpatches.FancyBboxPatch(
            (69.8, by),
            26.9,
            11.8,
            boxstyle="round,pad=0.35,rounding_size=0.9",
            facecolor=bbg,
            edgecolor=PURPLE_ACCENT,
            lw=1.0,
            zorder=3,
        )
        ax.add_patch(b)
        ax.text(
            83.25,
            by + 8.0,
            btitle,
            fontsize=9.0,
            fontweight="bold",
            color=INK,
            ha="center",
            va="center",
            zorder=4,
        )
        ax.text(
            83.25, by + 3.2, bsub, fontsize=7.7, color=MUTED, ha="center", va="center", zorder=4
        )

    ax.text(
        83.25,
        30.5,
        f"{arm3}\nOriented ({oriented.tested_studies:,} st.): Cellulose ↓ "
        f"p={fmt_p(oriented.fractions[CELLULOSE].p_value)} · "
        f"{len(oriented.ratio_survivors())} LR BY hits",
        fontsize=8.0,
        fontweight="bold",
        color=PURPLE_ACCENT,
        ha="center",
        va="center",
    )

    # Downward arrows from all 3 arms into the Bottom Convergence Ribbon
    for ax_x in (16.75, 50.0, 83.25):
        ax.annotate(
            "",
            xy=(ax_x, 18.8),
            xytext=(ax_x, 24.6),
            arrowprops={
                "arrowstyle": "-|>,head_width=0.7,head_length=0.7",
                "color": INK,
                "lw": 2.0,
            },
            zorder=5,
        )

    # --- Bottom Convergence Ribbon ---
    ribbon = mpatches.FancyBboxPatch(
        (1.5, 2.5),
        97.0,
        15.5,
        boxstyle="round,pad=0.4,rounding_size=1.2",
        facecolor="#fef9ec",
        edgecolor=GOLD_ACCENT,
        linewidth=1.8,
        zorder=2,
    )
    ax.add_patch(ribbon)
    ax.text(
        50.0,
        15.2,
        "CROSS-ARM CONVERGENCE FUNNEL (ALL THREE ARMS CONCORDANT IN DIRECTION)",
        fontsize=9.6,
        fontweight="bold",
        color=INK,
        ha="center",
    )

    conv_steps = [
        (
            3.5,
            "317 Biological Traits",
            f"Full KG-Microbe phenome\n({len(paired['pooled'])} observed in corpus)",
            "#ffffff",
            MID_RULE,
        ),
        (
            27.5,
            f"{concordant} Three-Arm Traits",
            "Arm 1 Study BH q ≤ 0.05\n+ Arm 3 BH₃₁₇ ≤ 0.05 + Arm 2 ORA",
            "#e8f5f3",
            TEAL_ACCENT,
        ),
        (
            51.5,
            f"{strict} Strict BY Survivors",
            "Pass conservative 317-trait\nBenjamini\N{EN DASH}Yekutieli q ≤ 0.05",
            "#eef5fb",
            DEC_COLOR,
        ),
        (
            75.5,
            "6 Table 1 + Oriented Audit",
            f"{len(increased_side(every))} Inc (↑) vs Cellulose (↓)\nOriented "
            f"{oriented.tested_studies:,} st.: Cell. BH₆ "
            f"q={fmt_p(oriented.fractions[CELLULOSE].q_value)}",
            "#f5d7d4",
            INC_COLOR,
        ),
    ]
    for idx_c, (cx, ctitle, csub, cbg, ccol) in enumerate(conv_steps):
        cb = mpatches.FancyBboxPatch(
            (cx, 4.2),
            21.0,
            8.8,
            boxstyle="round,pad=0.3,rounding_size=0.8",
            facecolor=cbg,
            edgecolor=ccol,
            lw=1.3,
            zorder=3,
        )
        ax.add_patch(cb)
        ax.text(
            cx + 10.5,
            10.2,
            ctitle,
            fontsize=8.8,
            fontweight="bold",
            color=INK,
            ha="center",
            va="center",
            zorder=4,
        )
        ax.text(cx + 10.5, 6.4, csub, fontsize=7.4, color=MUTED, ha="center", va="center", zorder=4)
        if idx_c < len(conv_steps) - 1:
            ax.annotate(
                "",
                xy=(conv_steps[idx_c + 1][0] - 0.4, 8.6),
                xytext=(cx + 21.4, 8.6),
                arrowprops={
                    "arrowstyle": "-|>,head_width=0.6,head_length=0.6",
                    "color": GOLD_ACCENT,
                    "lw": 1.8,
                },
                zorder=5,
            )

    fig.text(
        0.03,
        0.025,
        f"{RERUN_SOURCE} · Synthesizes Arm 1 (just trait-profile-investigation), Arm 2 "
        "(just analyze), Arm 3 (just table1-trait-maps, just trait-direction-screen, just "
        "table1-compositional-stats), just contrast-orientation and calibration.",
        fontsize=7.8,
        color=MUTED,
    )

    out_path = out_dir / "methods_plate2_three_arm_architecture.png"
    fig.savefig(out_path, dpi=220)
    fig.savefig(out_path.with_suffix(".svg"))
    plt.close(fig)
    return out_path


# ==============================================================================
# PLATE 1: Multi-Scale Inference & Attenuation Architecture (Visual Funnel + Slopegraph)
# ==============================================================================
# The six Table 1 traits Plate I follows: column id, label, display name and the direction
# the Table 1 study test found.
TABLE1_TRAITS = (
    ("METPO:1000608", "facultatively aerobic", "Facultatively aerobic", "inc"),
    ("METPO:2000103|GO:0030245", "cellulose catabolic process", "Cellulose catabolism", "dec"),
    ("METPO:1000668", "coccus shaped", "Coccus shaped", "inc"),
    (
        "METPO:2000103|pathways:arsenate_reduction",
        "arsenate reduction",
        "Arsenate reduction",
        "inc",
    ),
    ("METPO:1000683", "sphere shaped", "Sphere shaped", "inc"),
    ("METPO:1000623", "moderately halophilic", "Moderately halophilic", "inc"),
)
TRAIT_SHORT = {
    "facultatively aerobic": "facultative aerobe",
    "cellulose catabolic process": "cellulose",
    "coccus shaped": "coccus",
    "arsenate reduction": "arsenate",
    "sphere shaped": "sphere",
    "moderately halophilic": "halophily",
}
CONDITION_SHORT = {
    "Parkinson's disease": "PD",
    "HIV infection": "HIV",
    "Colorectal cancer": "CRC",
}


def has_value(value: str | None) -> bool:
    return value not in (None, "", "NA", "nan")


def plate1_summary(
    pooled: Mapping[str, dict[str, str]],
    feces: Mapping[str, dict[str, str]],
    sign: Mapping[str, dict[str, str]],
    contrasts: Sequence[dict[str, str]],
    wilcoxon: Sequence[dict[str, str]],
) -> dict[str, Any]:
    """Plate I's sample sizes, test families, survivors and points, read from its tables.

    The plate typed them and claimed survivors its own tables do not have (#636). A trait
    survives Scales 1-4 at BH q <= 0.05, and Scales 5 and 6 when the screen selects one of
    its cells after BH. Scale 5 plots each trait's lowest q among the cells that pass the
    screen's support rule, the cells BH can select.
    """
    cells = [
        r for r in contrasts if r.get("view_key") == "individual_conditions_direction_contrast"
    ]
    supported = [
        r for r in cells if r.get("passes_support") == "True" and has_value(r.get("fisher_q_value"))
    ]
    survivors = [0] * 6
    points: dict[str, list[float]] = {}
    selected: dict[int, list[str]] = {4: [], 5: []}
    for tid, label, _display, _direction in TABLE1_TRAITS:
        tests = [r for r in wilcoxon if r["trait"] == label and has_value(r.get("q_value"))]
        fisher = [float(r["fisher_q_value"]) for r in supported if r["trait_label"] == label]
        points[label] = [
            float(pooled[tid]["q_value"]),
            float(pooled[tid]["study_q_value"]),
            float(feces[tid]["study_q_value"]),
            float(sign[label]["bh_q"]),
            min(fisher, default=1.0),
            min((float(r["q_value"]) for r in tests), default=1.0),
        ]
        chosen = {
            4: [
                r
                for r in cells
                if r["trait_label"] == label and r.get("selected_after_bh") == "True"
            ],
            5: [r for r in tests if r.get("selected_after_bh") == "True"],
        }
        for index in range(4):
            survivors[index] += points[label][index] <= 0.05
        for index, rows in chosen.items():
            if rows:
                survivors[index] += 1
                places = sorted(
                    CONDITION_SHORT.get(r["condition"], r["condition"])
                    + (
                        ""
                        if r.get("body_site_group", "Feces") == "Feces"
                        else f" {r['body_site_group']}"
                    )
                    for r in rows
                )
                selected[index].append(f"{TRAIT_SHORT[label]}: {', '.join(places)}")
    studies = [int(r["paired_study_count"]) for r in wilcoxon]
    lowest = min(
        (r for r in wilcoxon if has_value(r.get("q_value"))), key=lambda r: float(r["q_value"])
    )
    sign_studies = max(
        int(r["studies_up"]) + int(r["studies_down"]) + int(r["studies_tied"])
        for r in sign.values()
    )
    return {
        "pairs": max(int(r["pair_count"]) for r in pooled.values()),
        "studies": max(int(float(r["study_count"])) for r in pooled.values()),
        "feces_studies": max(int(float(r["study_count"])) for r in feces.values()),
        "sign_studies": sign_studies,
        "m": [
            sum(has_value(r.get("q_value")) for r in pooled.values()),
            sum(has_value(r.get("study_q_value")) for r in pooled.values()),
            sum(has_value(r.get("study_q_value")) for r in feces.values()),
            int(next(iter(sign.values()))["family_size"]),
            len(cells),
            len(wilcoxon),
        ],
        "cell_profiles": round(
            median(
                int(r["increased_profile_count"]) + int(r["decreased_profile_count"])
                for r in supported
            )
        ),
        "wilcoxon_studies": (min(studies), max(studies)),
        "survivors": survivors,
        "selected": selected,
        "lowest_wilcoxon": (
            float(lowest["q_value"]),
            f"{TRAIT_SHORT.get(lowest['trait'], lowest['trait'])} in "
            f"{CONDITION_SHORT.get(lowest['condition'], lowest['condition'])}",
        ),
        "points": points,
    }


def plate1_subtitle(summary: Mapping[str, Any]) -> str:
    """Panel B's verdict on each scale, e.g. 'All 6 pass Scales 1-4; 2 pass Scale 5 (...)'."""
    survivors = summary["survivors"]
    if all(count == 6 for count in survivors[:4]):
        parts = ["All 6 pass Scales 1\N{EN DASH}4"]
    else:
        parts = [
            f"{count} {'passes' if count == 1 else 'pass'} Scale {index + 1}"
            for index, count in enumerate(survivors[:4])
        ]
    for index in (4, 5):
        if survivors[index]:
            verb = "passes" if survivors[index] == 1 else "pass"
            parts.append(
                f"{survivors[index]} {verb} Scale {index + 1} "
                f"({'; '.join(summary['selected'][index])})"
            )
        elif index == 5:
            q, where = summary["lowest_wilcoxon"]
            parts.append(f"none passes Scale 6 (lowest q = {fmt_p(q)}, {where})")
        else:
            parts.append("none passes Scale 5")
    return "; ".join(parts)


def build_plate1(repo: Path, out_dir: Path) -> Path:
    stats = repo / "trait-visualization/trait-profile-investigation"
    pooled_rows = {
        r["trait"]: r
        for r in read_tsv(stats / "trait_profile_paired_contrast_trait_stats_pooled.tsv")
    }
    feces_rows = {
        r["trait"]: r
        for r in read_tsv(stats / "trait_profile_paired_contrast_trait_stats_feces.tsv")
    }
    table1_sign = {
        r["trait"].replace("_", " "): r
        for r in read_tsv(repo / "data/derived/table1_traits/table1_traits_direction_shift.tsv")
    }
    contrasts = read_tsv(
        stats / "six-shared-traits/six_shared_trait_condition_body_site_direction_contrasts.tsv"
    )
    wilcoxon_rows = read_tsv(
        stats / "PD-cellulose/pd_condition_taxon_trait_network_condition_trait_tests.tsv"
    )
    summary = plate1_summary(pooled_rows, feces_rows, table1_sign, contrasts, wilcoxon_rows)
    m1, m2, m3, m4, m5, m6 = summary["m"]
    low, high = summary["wilcoxon_studies"]
    survivors = summary["survivors"]

    fig = plt.figure(figsize=(16.0, 9.0), facecolor=BG_COLOR)
    gs = fig.add_gridspec(
        1, 2, width_ratios=[0.78, 1.42], wspace=0.18, left=0.145, right=0.985, top=0.84, bottom=0.13
    )

    ax_funnel = fig.add_subplot(gs[0, 0])
    ax_slope = fig.add_subplot(gs[0, 1])

    fig.suptitle(
        "TUFTE PLATE 1 · Multi-Scale Statistical Attenuation Across Granularity & Multiplicity "
        "Regimes",
        fontsize=15.2,
        fontweight="bold",
        x=0.03,
        ha="left",
        y=0.96,
        color=INK,
    )
    fig.text(
        0.03,
        0.920,
        f"Why pooled signals (N = {summary['pairs']:,} pairs; N = {summary['studies']:,} studies) "
        "are overwhelmingly significant while Condition \N{MULTIPLICATION SIGN} Site study "
        "screens attenuate under sample partitioning",
        fontsize=9.2,
        color=MUTED,
    )

    # --- Panel A: Visual Sample-Size vs Multiplicity Burden Bars Across the 6 Scales ---
    def passed(index: int, noun: str = "Table 1 pass") -> str:
        return f"{survivors[index]} / 6 {noun}"

    scale_rows = [
        (
            "Scale 1 · Experiment Paired t",
            summary["pairs"],
            passed(0),
            TEAL_ACCENT,
            f"{summary['pairs']:,} pairs (m = {m1:,})",
        ),
        (
            "Scale 2 · Pooled Study Paired t",
            summary["studies"],
            passed(1),
            TEAL_ACCENT,
            f"{summary['studies']:,} studies (m = {m2:,})",
        ),
        (
            "Scale 3 · Feces Study Paired t",
            summary["feces_studies"],
            passed(2),
            TEAL_ACCENT,
            f"{summary['feces_studies']:,} studies (m = {m3:,})",
        ),
        (
            "Scale 4 · Pooled Study Sign Test",
            summary["sign_studies"],
            passed(3),
            TEAL_ACCENT,
            f"≤ {summary['sign_studies']:,} studies (m = {m4:,})",
        ),
        (
            "Scale 5 · Cond\N{MULTIPLICATION SIGN}Site Profile Fisher",
            summary["cell_profiles"],
            passed(4, "with a BH-selected cell"),
            GOLD_ACCENT,
            f"median {summary['cell_profiles']} profiles/supported cell (m = {m5:,})",
        ),
        (
            "Scale 6 · Cond\N{MULTIPLICATION SIGN}Site Study Wilcoxon",
            high,
            passed(5, "with a BH-selected cell"),
            INC_COLOR,
            f"{low}\N{EN DASH}{high} studies/cell (m = {m6:,})",
        ),
    ]
    y_f = np.arange(len(scale_rows))[::-1]
    ax_funnel.set_ylim(-0.6, len(scale_rows) - 0.2)
    ax_funnel.set_xlim(0, 4100)

    for idx_s, (_sname, n_unit, surv_txt, scol, n_lbl) in enumerate(scale_rows):
        yp = y_f[idx_s]
        ax_funnel.barh(yp, n_unit, height=0.56, color=scol, alpha=0.86)
        ax_funnel.text(
            n_unit + 55,
            yp + 0.11,
            f"{n_lbl}",
            fontsize=7.9,
            fontweight="bold",
            color=INK,
            va="center",
        )
        ax_funnel.text(
            n_unit + 55,
            yp - 0.16,
            f"Survivors: {surv_txt}",
            fontsize=7.4,
            fontweight="bold",
            color=scol,
            va="center",
        )

    ax_funnel.set_yticks(y_f)
    ax_funnel.set_yticklabels([r[0] for r in scale_rows], fontsize=8.0, fontweight="bold")
    ax_funnel.set_xlabel(
        "Effective Sample Size per Test Cell (N pairs, profiles, or studies)",
        fontsize=8.2,
        fontweight="bold",
    )
    ax_funnel.set_title(
        "A. Sample Partitioning & Multiplicity Funnel\nPartitioning "
        f"{summary['studies']:,} studies into Condition \N{MULTIPLICATION SIGN} Site cells "
        f"drops N to {low}\N{EN DASH}{high}",
        fontsize=9.5,
        fontweight="bold",
        loc="left",
        pad=10,
    )
    despine(ax_funnel)

    # --- Panel B: Spacious Multi-Scale -log10(q) Slopegraph ---
    scales = [
        f"Scale 1\nExp. Paired t\n({summary['pairs']:,} pairs\nm = {m1:,})",
        f"Scale 2\nPooled Study t\n({summary['studies']:,} studies\nm = {m2:,})",
        f"Scale 3\nFeces Study t\n({summary['feces_studies']:,} studies\nm = {m3:,})",
        f"Scale 4\nStudy Sign\n(≤ {summary['sign_studies']:,} studies\nm = {m4:,})",
        f"Scale 5\nBest supported\nCond\N{MULTIPLICATION SIGN}Site\nFisher (m={m5:,})",
        f"Scale 6\nBest Cond\N{MULTIPLICATION SIGN}Site\nWilcoxon (m={m6:,})",
    ]
    x_pos = np.arange(len(scales))

    ax_slope.set_ylim(-0.5, 9.8)

    ax_slope.axhspan(0, -math.log10(0.05), color="#f3ece6", alpha=0.65, zorder=0)
    ax_slope.axhline(-math.log10(0.05), color=INC_COLOR, linestyle="--", linewidth=1.1, zorder=1)
    ax_slope.axhline(-math.log10(0.001), color=TEAL_ACCENT, linestyle=":", linewidth=1.0, zorder=1)
    ax_slope.text(
        2.3,
        -math.log10(0.05) + 0.16,
        "BH FDR q = 0.05 threshold (-log10 q = 1.30)",
        fontsize=8.0,
        color=INC_COLOR,
        ha="center",
        fontweight="bold",
        bbox={
            "boxstyle": "round,pad=0.15",
            "facecolor": PANEL_BG,
            "edgecolor": "none",
            "alpha": 0.85,
        },
    )
    ax_slope.text(
        2.3,
        -math.log10(0.001) + 0.16,
        "BH q = 0.001 (-log10 q = 3.00)",
        fontsize=7.8,
        color=TEAL_ACCENT,
        ha="center",
    )

    for xi in x_pos:
        ax_slope.axvline(xi, color=LIGHT_RULE, linewidth=0.9, zorder=0)

    trait_colors = {
        "facultatively aerobic": "#b8332a",
        "cellulose catabolic process": "#1b6ca8",
        "coccus shaped": "#d96b27",
        "arsenate reduction": "#7b4f9d",
        "sphere shaped": "#2e8b57",
        "moderately halophilic": "#b07d12",
    }
    right_y_Explicit = {
        "moderately halophilic": 2.20,
        "cellulose catabolic process": 1.62,
        "arsenate reduction": 1.05,
        "facultatively aerobic": 0.52,
        "coccus shaped": 0.02,
        "sphere shaped": -0.38,
    }

    for _tid, tname, tdisp, tdir in TABLE1_TRAITS:
        q_vals = summary["points"][tname]
        wilc_q = q_vals[-1]
        y_vals = [min(9.2, -math.log10(max(q, 1e-12))) for q in q_vals]
        c = trait_colors[tname]
        marker = "o" if tdir == "inc" else "s"

        ax_slope.plot(
            x_pos,
            y_vals,
            color=c,
            linestyle="-",
            linewidth=2.3,
            marker=marker,
            markersize=6.2,
            alpha=0.92,
            zorder=3,
        )

        # Clean right-hand legend callouts with leader lines from Scale 6
        ry = right_y_Explicit[tname]
        ax_slope.annotate(
            f"{tdisp} (q={fmt_p(wilc_q)})",
            xy=(len(scales) - 1, y_vals[-1]),
            xytext=(len(scales) - 0.85, ry),
            fontsize=7.7,
            fontweight="bold" if wilc_q <= 0.10 else "normal",
            color=c,
            va="center",
            ha="left",
            arrowprops={"arrowstyle": "-", "color": c, "lw": 0.8, "alpha": 0.7},
        )

    ax_slope.set_xlim(-0.20, len(scales) + 0.92)
    ax_slope.set_xticks(x_pos)
    ax_slope.set_xticklabels(scales, fontsize=7.9, fontweight="bold")
    ax_slope.set_ylabel(
        "Multiplicity-Adjusted Significance  $-\\log_{10}(\\mathrm{BH}\\ q)$  [capped at 9.2]",
        fontsize=8.8,
        fontweight="bold",
    )
    ax_slope.set_title(
        "B. Significance Trajectory of the 6 Table 1 Traits Across Scales 1\N{EN DASH}6\n"
        + textwrap.fill(plate1_subtitle(summary), width=96),
        fontsize=9.6,
        fontweight="bold",
        loc="left",
        pad=10,
    )
    despine(ax_slope)

    fig.text(
        0.05,
        0.025,
        f"{RERUN_SOURCE} · "
        "trait_profile_paired_contrast_trait_stats_{pooled,feces}.tsv, "
        "table1_traits_direction_shift.tsv, "
        "six_shared_trait_condition_body_site_direction_contrasts.tsv, "
        "pd_condition_taxon_trait_network_condition_trait_tests.tsv.",
        fontsize=8.0,
        color=MUTED,
    )

    out_path = out_dir / "tufte_plate1_scale_attenuation_slopegraph.png"
    fig.savefig(out_path, dpi=220)
    fig.savefig(out_path.with_suffix(".svg"))
    plt.close(fig)
    return out_path


# ==============================================================================
# PLATE 2: Beyond Table 1 — The 28-Trait Functional Guild Dichotomy & Anatomy
# ==============================================================================
def build_plate2(repo: Path, out_dir: Path) -> Path:
    pooled = {
        r["trait"]: r
        for r in read_tsv(
            repo / "trait-visualization/trait-profile-investigation/"
            "trait_profile_paired_contrast_trait_stats_pooled.tsv"
        )
    }
    feces = {
        r["trait"]: r
        for r in read_tsv(
            repo / "trait-visualization/trait-profile-investigation/"
            "trait_profile_paired_contrast_trait_stats_feces.tsv"
        )
    }
    nonfeces = {
        r["trait"]: r
        for r in read_tsv(
            repo / "trait-visualization/trait-profile-investigation/"
            "trait_profile_paired_contrast_trait_stats_feces_excluded.tsv"
        )
    }
    saliva = {
        r["trait"]: r
        for r in read_tsv(
            repo / "trait-visualization/trait-profile-investigation/"
            "trait_profile_paired_contrast_trait_stats_saliva.tsv"
        )
    }

    labels = trait_labels(repo)
    trait_rows_spec = [
        # Increased Guilds
        ("METPO:1000608", "[O₂ Guild]  Facultatively aerobic ★", "O₂ / Aerotolerance"),
        ("METPO:1000606", "Microaerophilic", "O₂ / Aerotolerance"),
        ("METPO:1000602", "Aerobic", "O₂ / Aerotolerance"),
        ("METPO:1000623", "[Halophily]  Moderately halophilic ★", "Halotolerance / Osmotic"),
        ("METPO:1000467", "NaCl optimum mid2", "Halotolerance / Osmotic"),
        ("METPO:1000627", "Euryhaline", "Halotolerance / Osmotic"),
        (
            "METPO:2000103|pathways:arsenate_reduction",
            "[Redox]  Arsenate reduction ★",
            "Redox",
        ),
        (
            "METPO:2000103|pathways:trithionate_oxidation",
            "Trithionate oxidation",
            "Redox",
        ),
        (
            "METPO:2000103|pathways:thiosulfate_oxidation",
            "Thiosulfate oxidation",
            "Redox",
        ),
        (
            "METPO:2000103|pathways:thiosulfate_oxidation_dark",
            "Thiosulfate oxidation (dark)",
            "Redox",
        ),
        ("METPO:2000103|pathways:nitrate_reduction", "Nitrate reduction", "Redox"),
        (
            "METPO:2000103|pathways:methane_oxidation",
            "Methane oxidation",
            "Redox",
        ),
        ("METPO:1000668", "[Morphology]  Coccus shaped ★", "Morphology"),
        ("METPO:1000683", "Sphere shaped ★", "Morphology"),
        # Carbon sources and lignin catabolism; the tag sits on the shortest label so the
        # row fits the left margin.
        ("METPO:2000006|CHEBI:25115", "[Carbon use]  Malate", "Carbon use"),
        ("METPO:2000006|CHEBI:17596", "Inosine", "Carbon use"),
        ("METPO:2000006|CHEBI:16808", "2-Dehydro-D-gluconate", "Carbon use"),
        ("METPO:2000006|CHEBI:320055", "Methyl β-D-glucopyranoside", "Carbon use"),
        ("METPO:2000103|GO:0046274", "Lignin catabolic process", "Carbon use"),
        # Decreased Guilds
        (
            "METPO:2000103|GO:0030245",
            "[Plant glycan ↓]  Cellulose catabolic process ★",
            "Plant glycan ↓",
        ),
        ("METPO:2000103|GO:0045493", "Xylan catabolic process", "Plant glycan ↓"),
        ("METPO:2000006|CHEBI:16634", "Raffinose", "Plant glycan ↓"),
        ("METPO:1000429", "[Anaerobic / low GC ↓]  GC low", "Anaerobic / low GC ↓"),
        ("METPO:1000603", "Anaerobic", "Anaerobic / low GC ↓"),
        ("METPO:1000607", "Obligately anaerobic", "Anaerobic / low GC ↓"),
        ("METPO:1000845", "Acetogenesis", "Anaerobic / low GC ↓"),
    ]
    check_plate_two_rows([trait for trait, _, _ in trait_rows_spec], pooled)
    check_plate_two_directions(trait_rows_spec, pooled)
    # Every row takes the label its column has in taxa_trait_matrix_columns.tsv (#607).
    trait_rows_spec = [
        (trait_id, canonical_display(display, labels.get(trait_id)), group)
        for trait_id, display, group in trait_rows_spec
    ]

    fig = plt.figure(figsize=(16.0, 9.2), facecolor=BG_COLOR)
    gs = fig.add_gridspec(
        1,
        4,
        width_ratios=[1.0, 1.0, 1.0, 1.0],
        wspace=0.12,
        left=0.23,
        right=0.98,
        top=0.84,
        bottom=0.11,
    )

    fig.suptitle(
        f"TUFTE PLATE 2 · Beyond Table 1: {len(trait_rows_spec)} Traits by Direction of "
        "Annotation Shift & Body Site",
        fontsize=15.0,
        fontweight="bold",
        x=0.03,
        ha="left",
        y=0.96,
        color=INK,
    )
    fig.text(
        0.03,
        0.922,
        "Study-level paired contrasts (mean Δz ± 95% CI) · ** study BH q ≤ 0.005 · * q ≤ 0.05 · "
        f"† nominal p ≤ 0.05 · Feces ({stratum_size(feces)[0]:,} studies) carries most of the "
        "pooled signal; "
        f"{study_bh_count(nonfeces)} traits pass BH without Feces and {study_bh_count(saliva)} "
        "in Saliva",
        fontsize=9.0,
        color=MUTED,
    )

    strata = [
        (
            f"Pooled All Sites\n({stratum_text(pooled)}\n"
            f"{study_bh_count(pooled)} traits study BH q ≤ 0.05)",
            pooled,
            (-0.38, 0.72),
        ),
        (
            f"Feces Only\n({stratum_text(feces)}\n"
            f"{study_bh_count(feces)} traits study BH q ≤ 0.05)",
            feces,
            (-0.44, 0.78),
        ),
        (
            f"Feces Excluded\n({stratum_text(nonfeces)}\n"
            f"{study_bh_count(nonfeces)} traits study BH q ≤ 0.05)",
            nonfeces,
            (-0.40, 0.56),
        ),
        (
            f"Saliva Only\n({stratum_text(saliva)}\n"
            f"{study_bh_count(saliva)} traits study BH q ≤ 0.05)",
            saliva,
            (-0.78, 0.92),
        ),
    ]

    n_rows = len(trait_rows_spec)
    y_positions = np.arange(n_rows)[::-1]
    # The rule between the increased groups and the first decreased (↓) group
    first_decreased = next(i for i, row in enumerate(trait_rows_spec) if row[2].endswith("↓"))
    div_y = y_positions[first_decreased] + 0.5

    for col_idx, (col_title, data_dict, xlim) in enumerate(strata):
        ax = fig.add_subplot(gs[0, col_idx])
        ax.set_ylim(-0.7, n_rows - 0.3)
        ax.set_xlim(*xlim)
        ax.axvline(0.0, color=INK, linestyle="-", linewidth=0.9, zorder=1)
        ax.axhline(div_y, color=INK, linestyle="--", linewidth=1.0, zorder=2)

        for idx_r, yp in enumerate(y_positions):
            if idx_r % 2 == 1:
                ax.axhspan(yp - 0.45, yp + 0.45, color="#f4f1ea", alpha=0.65, zorder=0)

        for idx_r, (tid, tlabel, _guild_lbl) in enumerate(trait_rows_spec):
            yp = y_positions[idx_r]
            row = data_dict.get(tid)
            if not row or not row.get("study_mean_increased_minus_decreased"):
                # Not tested in this stratum (see the footer)
                ax.text(
                    0.02 * (xlim[1] - xlim[0]),
                    yp,
                    "n.t.",
                    fontsize=7.0,
                    color=MUTED,
                    ha="left",
                    va="center",
                )
                continue
            mean_v = float(row["study_mean_increased_minus_decreased"])
            ci_lo = float(row["study_ci95_lower"])
            ci_hi = float(row["study_ci95_upper"])
            st_t = float(row["study_t_statistic"])
            st_p = float(row["study_p_value"])
            st_q = float(row["study_q_value"])

            is_sig_q = st_q <= 0.05
            is_nom_p = st_p <= 0.05
            c = INC_COLOR if mean_v >= 0 else DEC_COLOR
            alpha = 0.95 if (is_sig_q or is_nom_p) else 0.42

            ax.plot(
                [ci_lo, ci_hi],
                [yp, yp],
                color=c,
                linewidth=1.9 if is_sig_q else 1.1,
                alpha=alpha,
                zorder=3,
            )
            ax.scatter(
                [mean_v],
                [yp],
                color=c if is_sig_q else (PANEL_BG if not is_nom_p else c),
                edgecolor=c,
                linewidth=1.3,
                s=36 if is_sig_q else 22,
                zorder=4,
            )

            sig_tag = (
                "**" if st_q <= 0.005 else ("*" if st_q <= 0.05 else ("†" if st_p <= 0.05 else ""))
            )
            if sig_tag:
                ax.text(
                    xlim[1] - 0.02 * (xlim[1] - xlim[0]),
                    yp,
                    f"t={st_t:+.1f}{sig_tag}",
                    fontsize=7.3,
                    fontweight="bold" if is_sig_q else "normal",
                    color=c,
                    ha="right",
                    va="center",
                )

            if col_idx == 0:
                ax.text(
                    xlim[0] - 0.03 * (xlim[1] - xlim[0]),
                    yp,
                    tlabel,
                    fontsize=8.2,
                    fontweight="bold"
                    if ("★" in tlabel or "Saliva" in tlabel or "[" in tlabel)
                    else "normal",
                    color=INK,
                    ha="right",
                    va="center",
                )

        ax.set_yticks([])
        ax.set_title(col_title, fontsize=9.5, fontweight="bold", pad=8)
        ax.set_xlabel("Study Mean Δz (Inc \N{MINUS SIGN} Dec)", fontsize=8.4, fontweight="bold")
        ax.tick_params(axis="x", labelsize=7.8)
        despine(ax, left=True)

    fig.text(
        0.03,
        0.012,
        f"{RERUN_SOURCE} · "
        "trait_profile_paired_contrast_trait_stats_{pooled,feces,feces_excluded,saliva}.tsv · "
        "★ Table 1 trait · rows show the object_name of taxa_trait_matrix_columns.tsv, with "
        '"(carbon source)" for compounds used as carbon sources.\nn.t. = not tested in that '
        "stratum: the stratified screen tests a trait carried by at least 6 and at most 95% of "
        "the stratum's profiles, with nonzero variance.",
        fontsize=7.9,
        color=MUTED,
        va="bottom",
        linespacing=1.5,
    )

    out_path = out_dir / "tufte_plate2_317_trait_guild_dichotomy.png"
    fig.savefig(out_path, dpi=220)
    fig.savefig(out_path.with_suffix(".svg"))
    plt.close(fig)
    return out_path


# ==============================================================================
# PLATE 3: Compositional Closure & Reference-Explicit Paired Log-Ratio Matrix
# ==============================================================================
def build_plate3(repo: Path, out_dir: Path) -> Path:
    lr_rows = read_tsv(repo / "data/derived/table1_composition/trait_log_ratio_stats.tsv")
    orient_rows = read_tsv(repo / COMPARISON)
    signs = table1_signs(repo)
    every, oriented = signs[ALL_STUDIES], signs[ORIENTED]
    version = read_json(repo / CONTRAST_PROVENANCE)["classifier_version"]
    orient_lr = {}
    for r in orient_rows:
        if (
            r.get("comparison") == "log_ratio"
            and r.get("population") == "disease_vs_healthy_oriented"
        ):
            orient_lr[(r["numerator"], r["reference"])] = (r, 1.0)
            orient_lr[(r["reference"], r["numerator"])] = (r, -1.0)

    short_names = {
        "facultatively_aerobic": "Fac. aerobic (FA)",
        "cellulose_catabolic_process": "Cellulose cat. (CEL)",
        "coccus_shaped": "Coccus shaped (COC)",
        "sphere_shaped": "Sphere shaped (SPH)",
        "arsenate_reduction": "Arsenate red. (ARS)",
        "moderately_halophilic": "Mod. halophilic (HAL)",
    }
    codes = {
        "facultatively_aerobic": "FA",
        "cellulose_catabolic_process": "CEL",
        "coccus_shaped": "COC",
        "sphere_shaped": "SPH",
        "arsenate_reduction": "ARS",
        "moderately_halophilic": "HAL",
    }
    trait_order = [
        "facultatively_aerobic",
        "coccus_shaped",
        "sphere_shaped",
        "moderately_halophilic",
        "arsenate_reduction",
        "cellulose_catabolic_process",
    ]

    fig = plt.figure(figsize=(16.0, 9.0), facecolor=BG_COLOR)
    gs = fig.add_gridspec(
        1, 2, width_ratios=[0.90, 1.40], wspace=0.19, left=0.12, right=0.985, top=0.84, bottom=0.11
    )

    fig.suptitle(
        "TUFTE PLATE 3 · Compositional Closure: Reference-Explicit Paired Log-Ratios Across the 15 "
        "Table 1 Trait Pairs",
        fontsize=15.0,
        fontweight="bold",
        x=0.04,
        ha="left",
        y=0.96,
        color=INK,
    )
    fig.text(
        0.04,
        0.922,
        "Within-experiment Δlog₂(Trait A / Trait B) where signature size T cancels: "
        f"{len(every.ratio_survivors())} BY survivors in all {every.studies:,} studies "
        f"({survivor_codes(every)}) → {len(oriented.ratio_survivors())} in Disease vs. Healthy "
        f"Oriented ({oriented.studies:,} studies; {survivor_codes(oriented)})",
        fontsize=8.9,
        color=MUTED,
    )

    others = [trait for trait in every.fractions if trait != CELLULOSE]
    cellulose_pairs = [(CELLULOSE, trait) for trait in others]
    increased_pairs = list(itertools.combinations(others, 2))

    # --- Panel A: 6x6 Matrix of 15 Pairs ---
    ax_mat = fig.add_subplot(gs[0, 0])
    n_t = len(trait_order)
    ax_mat.set_xlim(-0.5, n_t - 0.5)
    ax_mat.set_ylim(-0.5, n_t - 0.5)
    ax_mat.invert_yaxis()

    pair_lookup = {}
    for r in lr_rows:
        a, b = r["numerator"], r["reference"]
        pair_lookup[(a, b)] = (r, 1.0)
        pair_lookup[(b, a)] = (r, -1.0)

    for i, tr_row in enumerate(trait_order):
        for j, tr_col in enumerate(trait_order):
            if i == j:
                rect = mpatches.Rectangle(
                    (j - 0.48, i - 0.48),
                    0.96,
                    0.96,
                    facecolor="#e8ecef",
                    edgecolor=MID_RULE,
                    lw=0.8,
                )
                ax_mat.add_patch(rect)
                ax_mat.text(
                    j,
                    i,
                    codes[tr_row],
                    ha="center",
                    va="center",
                    fontsize=9.8,
                    fontweight="bold",
                    color=INK,
                )
                continue

            r, sign = pair_lookup[(tr_row, tr_col)]
            or_r, _ = orient_lr.get((tr_row, tr_col), ({}, 1.0))
            or_by_q = (
                float(or_r["q_value"])
                if or_r and or_r.get("q_value") not in (None, "", "NA")
                else 1.0
            )

            mean_lr = float(r["mean_log2_ratio_change"]) * sign
            share_up = (
                float(r["share_ratio_up"]) if sign > 0 else (1.0 - float(r["share_ratio_up"]))
            )
            by_q = float(r["by_q"])
            up_n = int(r["studies_up"]) if sign > 0 else int(r["studies_down"])
            dn_n = int(r["studies_down"]) if sign > 0 else int(r["studies_up"])

            if i < j:
                bg = INC_LIGHT if mean_lr > 0.15 else (DEC_LIGHT if mean_lr < -0.15 else "#f7f9fa")
                both_sig = by_q <= 0.05
                or_only_sig = (by_q > 0.05) and (or_by_q <= 0.05)
                edge_c = (
                    GOLD_ACCENT
                    if or_only_sig
                    else (
                        INC_COLOR
                        if (both_sig and mean_lr > 0)
                        else (DEC_COLOR if (both_sig and mean_lr < 0) else MID_RULE)
                    )
                )
                lw = 2.4 if both_sig else (2.2 if or_only_sig else 0.7)
                ls = "--" if or_only_sig else "-"
                rect = mpatches.Rectangle(
                    (j - 0.47, i - 0.47),
                    0.94,
                    0.94,
                    facecolor=bg,
                    edgecolor=edge_c,
                    lw=lw,
                    linestyle=ls,
                )
                ax_mat.add_patch(rect)
                if both_sig:
                    q_str = f"q={fmt_p(by_q)}*"
                elif or_only_sig:
                    q_str = f"Or.q={fmt_p(or_by_q)}★"
                else:
                    q_str = f"q={fmt_p(by_q)}"
                ax_mat.text(
                    j,
                    i - 0.12,
                    f"{mean_lr:+.2f}",
                    ha="center",
                    va="center",
                    fontsize=8.8,
                    fontweight="bold",
                    color=INK,
                )
                ax_mat.text(
                    j,
                    i + 0.18,
                    q_str,
                    ha="center",
                    va="center",
                    fontsize=6.9 if or_only_sig else 7.4,
                    color="#92400e" if or_only_sig else (INC_COLOR if both_sig else MUTED),
                    fontweight="bold" if (both_sig or or_only_sig) else "normal",
                )
            else:
                bg = "#f4f6f8" if (by_q > 0.05 and or_by_q > 0.05) else "#fff6e5"
                rect = mpatches.Rectangle(
                    (j - 0.47, i - 0.47), 0.94, 0.94, facecolor=bg, edgecolor=MID_RULE, lw=0.7
                )
                ax_mat.add_patch(rect)
                ax_mat.text(
                    j,
                    i - 0.12,
                    f"{share_up * 100:.1f}%",
                    ha="center",
                    va="center",
                    fontsize=8.5,
                    fontweight="bold",
                    color=INK,
                )
                ax_mat.text(
                    j,
                    i + 0.18,
                    f"{up_n}:{dn_n}",
                    ha="center",
                    va="center",
                    fontsize=7.4,
                    color=MUTED,
                )

    ax_mat.set_xticks(range(n_t))
    ax_mat.set_xticklabels([codes[t] for t in trait_order], fontsize=9.0, fontweight="bold")
    ax_mat.set_yticks(range(n_t))
    ax_mat.set_yticklabels([short_names[t] for t in trait_order], fontsize=8.6, fontweight="bold")
    ax_mat.set_title(
        "A. 6\N{MULTIPLICATION SIGN}6 Pairwise Log-Ratio Matrix (Row / Col)\nUpper: Mean Δlog₂ & "
        "BY q (dashed ★ = Oriented BY≤0.05) · Lower: % Up",
        fontsize=9.6,
        fontweight="bold",
        loc="left",
        pad=10,
    )
    despine(ax_mat, top=True, right=True, left=True, bottom=True)

    # --- Panel B: Spacious Forest Plot of Share of Studies with Ratio > 0 (All 959 vs Oriented 569)
    # ---
    ax_for = fig.add_subplot(gs[0, 1])
    cel_pairs = [
        r for r in lr_rows if "cellulose_catabolic_process" in (r["numerator"], r["reference"])
    ]
    other_pairs = [
        r for r in lr_rows if "cellulose_catabolic_process" not in (r["numerator"], r["reference"])
    ]
    cel_pairs.sort(key=lambda r: float(r["by_q"]))
    other_pairs.sort(key=lambda r: float(r["sign_p"]))
    sorted_lr = cel_pairs + other_pairs

    y_pos = np.arange(len(sorted_lr))[::-1]
    ax_for.set_ylim(-0.6, len(sorted_lr) + 0.4)
    ax_for.set_xlim(22.0, 118.0)
    ax_for.set_xticks([30, 40, 50, 60, 70])
    ax_for.axvline(50.0, color=INK, linestyle="--", linewidth=1.1, zorder=1)
    ax_for.axhline(y_pos[4] - 0.5, color=INK, linestyle=":", linewidth=1.1, zorder=1)

    ax_for.text(
        50.8,
        len(sorted_lr) - 0.05,
        "Null = 50% (co-moving traits cancel)",
        fontsize=7.9,
        color=MUTED,
        style="italic",
    )
    # The labels name the sides only while the table puts every other trait on the
    # increased side and cellulose on the decreased one (#704)
    others_up = all(every.fractions[trait].share_up > 0.5 for trait in others)
    cellulose_down = every.fractions[CELLULOSE].share_up < 0.5
    cellulose_group = "Increased Traits" if others_up and cellulose_down else "the Other Traits"
    passing_all = passing_pairs(every, cellulose_pairs)
    passing_oriented = passing_pairs(oriented, cellulose_pairs)
    ax_for.text(
        23.5,
        y_pos[0] + 0.42,
        f"Cellulose vs {cellulose_group} ({len(cellulose_pairs)} pairs: "
        f"{passing_all} {'passes' if passing_all == 1 else 'pass'} BY in all "
        f"{every.studies:,} [●] → {passing_oriented} "
        f"{'passes' if passing_oriented == 1 else 'pass'} in {oriented.studies:,} oriented [◆])",
        fontsize=8.0,
        fontweight="bold",
        color=DEC_COLOR,
    )
    ax_for.text(
        23.5,
        y_pos[5] + 0.48,
        f"{'Increased-Side vs Increased-Side Traits' if others_up else 'Pairs Without Cellulose'} "
        f"({len(increased_pairs)} pairs: "
        f"{increased_pairs_text(every, oriented, increased_pairs)})",
        fontsize=8.0,
        fontweight="bold",
        color=MUTED,
        bbox={"facecolor": BG_COLOR, "edgecolor": "none", "pad": 1.5},
        zorder=5,
    )

    for idx_r, r in enumerate(sorted_lr):
        yp = y_pos[idx_r]
        num_c = codes[r["numerator"]]
        ref_c = codes[r["reference"]]
        lbl = f"{num_c} / {ref_c}"
        sh = float(r["share_ratio_up"]) * 100.0
        lo = float(r["share_ci_low"]) * 100.0
        hi = float(r["share_ci_high"]) * 100.0
        by_q = float(r["by_q"])
        mean_lr = float(r["mean_log2_ratio_change"])
        up_n = int(r["studies_up"])
        dn_n = int(r["studies_down"])
        elig = int(r["studies_eligible"])

        or_r, _ = orient_lr.get((r["numerator"], r["reference"]), ({}, 1.0))
        or_sh = (
            float(or_r["share_up"]) * 100.0
            if or_r and or_r.get("share_up") not in (None, "", "NA")
            else None
        )
        or_by_q = (
            float(or_r["q_value"]) if or_r and or_r.get("q_value") not in (None, "", "NA") else 1.0
        )
        or_up = (
            int(or_r["studies_up"])
            if or_r and or_r.get("studies_up") not in (None, "", "NA")
            else 0
        )
        or_dn = (
            int(or_r["studies_down"])
            if or_r and or_r.get("studies_down") not in (None, "", "NA")
            else 0
        )

        is_sig = by_q <= 0.05
        or_sig = or_by_q <= 0.05
        col = (INC_COLOR if sh > 50 else DEC_COLOR) if (is_sig or or_sig) else NEUTRAL_INK
        ax_for.plot(
            [lo, hi], [yp, yp], color=col, lw=2.3 if is_sig else (1.8 if or_sig else 1.2), zorder=3
        )
        ax_for.scatter(
            [sh],
            [yp],
            color=col if is_sig else PANEL_BG,
            edgecolor=col,
            s=42 if is_sig else 26,
            lw=1.4,
            zorder=4,
        )

        if or_sh is not None:
            ax_for.scatter(
                [or_sh],
                [yp],
                marker="D",
                s=46 if or_sig else 24,
                facecolor=GOLD_ACCENT if or_sig else "#fef9ec",
                edgecolor=INK if or_sig else col,
                lw=1.1,
                zorder=5,
            )

        ax_for.text(
            21.0,
            yp,
            lbl,
            fontsize=8.5,
            fontweight="bold" if (is_sig or or_sig) else "normal",
            color=INK,
            ha="right",
            va="center",
        )
        if is_sig and or_sig:
            q_badge = f"BY q={fmt_p(by_q)}* (Or:{fmt_p(or_by_q)}★)"
        elif or_sig:
            q_badge = (
                f"All q={fmt_p(by_q)} → Or:{or_up}:{or_dn} ({or_sh:.1f}%) BY q={fmt_p(or_by_q)}★"
            )
        else:
            q_badge = f"BY q={fmt_p(by_q)}  [N={elig}]"

        if is_sig:
            row_txt = f"{up_n}:{dn_n} ({sh:.1f}%) · Δlog₂={mean_lr:+.2f} · {q_badge}"
        elif or_sig:
            row_txt = f"All {up_n}:{dn_n} ({sh:.1f}%) · {q_badge}"
        else:
            row_txt = f"{up_n}:{dn_n} ({sh:.1f}%) · Δlog₂={mean_lr:+.2f} · {q_badge}"

        ax_for.text(
            117.0,
            yp,
            row_txt,
            fontsize=7.15 if (is_sig or or_sig) else 7.5,
            fontweight="bold" if (is_sig or or_sig) else "normal",
            color=col if (is_sig or or_sig) else MUTED,
            ha="right",
            va="center",
        )

    ax_for.set_yticks([])
    ax_for.set_xlabel(
        f"% of Eligible Studies with Δlog₂(Numerator / Reference) > 0  (● All {every.studies:,} "
        f"st. ±95% CI  |  ◆ Oriented Disease-vs-Healthy {oriented.studies:,} st.)",
        fontsize=8.3,
        fontweight="bold",
    )
    ax_for.set_title(
        f"B. Directional Consistency Across All 15 Pairs (Unoriented All {every.studies:,} vs. "
        f"Oriented {oriented.studies:,} Studies)\n"
        f"{oriented_cellulose_pairs_text(every, oriented, cellulose_pairs)}",
        fontsize=9.6,
        fontweight="bold",
        loc="left",
        pad=10,
    )
    despine(ax_for, left=True)

    fig.text(
        0.04,
        0.025,
        f"{RERUN_SOURCE} · trait_log_ratio_stats.tsv & "
        f"disease_control_comparison.tsv (just contrast-orientation, {version}). *All-studies BY "
        f"q ≤ 0.05; ★Oriented disease-vs-healthy ({oriented.studies:,} st.) BY q ≤ 0.05.",
        fontsize=7.9,
        color=MUTED,
    )

    out_path = out_dir / "tufte_plate3_compositional_log_ratio_matrix.png"
    fig.savefig(out_path, dpi=220)
    fig.savefig(out_path.with_suffix(".svg"))
    plt.close(fig)
    return out_path


# ==============================================================================
# PLATE 4: Clean Visual Condition x Trait Matrix + Lead Fecal Cohort Forest Plot
# ==============================================================================
def build_plate4(repo: Path, out_dir: Path) -> Path:
    cellulose = cellulose_decr_notes(REPO_ROOT)
    eco_summary = read_tsv(
        repo / "trait-visualization/trait-profile-investigation/PD-cellulose/"
        "pd_ecological_reporter_combined_condition_trait_summary.tsv"
    )
    wilc_tests = read_tsv(
        repo / "trait-visualization/trait-profile-investigation/PD-cellulose/"
        "pd_condition_taxon_trait_network_condition_trait_tests.tsv"
    )
    fisher_rows = read_tsv(
        repo / "trait-visualization/trait-profile-investigation/six-shared-traits/"
        "six_shared_trait_condition_body_site_direction_contrasts.tsv"
    )

    by_eco = {(r["condition"], r["trait"]): r for r in eco_summary}
    by_wilc = {(r["condition"], r["trait"]): r for r in wilc_tests}
    by_fish = {
        (r["condition"], r["trait_label"]): r
        for r in fisher_rows
        if r.get("view_key") == "individual_conditions_direction_contrast"
        and r.get("body_site_group") == "Feces"
    }

    # Each row names its paired studies, read from the reporter table (#579)
    studies = {
        r["condition"]: int(r["paired_study_count"])
        for r in eco_summary
        if r["trait"] == "ecological-reporter composite"
    }
    conditions = [
        (condition, f"{shown} ({studies[condition]} st.)")
        for condition, shown in (
            ("Parkinson's disease", "Parkinson's disease"),
            ("Colorectal cancer", "Colorectal cancer"),
            ("HIV infection", "HIV infection"),
            ("Irritable bowel syndrome", "Irritable bowel syndrome"),
            ("COVID-19", "COVID-19"),
            ("Gestational diabetes", "Gestational diabetes"),
            ("Age", "Age / aging"),
            ("Acute lymphoblastic leukemia", "Acute lymphoblastic leukemia"),
            ("Atopic eczema", "Atopic eczema"),
            ("Autism", "Autism"),
            ("Stroke", "Stroke"),
            ("Ulcerative colitis", "Ulcerative colitis"),
            ("Clinical treatment", "Clinical treatment"),
            ("Socioeconomic status", "Socioeconomic status"),
        )
    ]

    traits_cols = [
        ("cellulose catabolic process", "Cellulose\ncatabolism"),
        ("facultatively aerobic", "Facultatively\naerobic"),
        ("moderately halophilic", "Moderately\nhalophilic"),
        ("arsenate reduction", "Arsenate\nreduction"),
        ("coccus shaped", "Coccus\nshaped"),
        ("sphere shaped", "Sphere\nshaped"),
        ("ecological-reporter composite", "Ecological\nComposite"),
    ]

    fig = plt.figure(figsize=(16.0, 9.1), facecolor=BG_COLOR)
    gs = fig.add_gridspec(
        1, 2, width_ratios=[1.12, 0.98], wspace=0.48, left=0.17, right=0.98, top=0.80, bottom=0.11
    )
    ax = fig.add_subplot(gs[0, 0])
    ax_lead = fig.add_subplot(gs[0, 1])

    fig.suptitle(
        "TUFTE PLATE 4 · Fecal Condition \N{MULTIPLICATION SIGN} Trait Paired-Shift Matrix & Lead "
        "Cohort Triangulation",
        fontsize=15.0,
        fontweight="bold",
        x=0.03,
        ha="left",
        y=0.96,
        color=INK,
    )
    fig.text(
        0.03,
        0.920,
        f"Left: Study-weighted mean paired shift Δz across {len(conditions)} Fecal cohorts "
        "\N{MULTIPLICATION SIGN} "
        "6 Table 1 traits + Ecological Composite · Right: Study win/loss consistency & 95% CIs for "
        "the lead disease cohorts",
        fontsize=8.9,
        color=MUTED,
    )

    # --- Panel A: Clean Visual Heatmap (1 concise value per cell!) ---
    n_r = len(conditions)
    n_c = len(traits_cols)
    ax.set_xlim(-0.5, n_c - 0.5)
    ax.set_ylim(-0.5, n_r - 0.5)
    ax.invert_yaxis()

    for r_idx, (cond_raw, _cond_disp) in enumerate(conditions):
        for c_idx, (tr_raw, _tr_disp) in enumerate(traits_cols):
            er = by_eco.get((cond_raw, tr_raw), {})
            wr = by_wilc.get((cond_raw, tr_raw), {})
            fr = by_fish.get((cond_raw, tr_raw), {})

            if not er or er.get("study_weighted_mean_contrast") in (None, "", "NA"):
                continue

            mean_z = float(er["study_weighted_mean_contrast"])
            ci_excl = er.get("ci_excludes_zero") == "True"
            w_p = float(wr["p_value"]) if wr and wr.get("p_value") else 1.0
            w_q = float(wr["q_value"]) if wr and wr.get("q_value") else 1.0
            f_q = (
                float(fr["fisher_q_value"])
                if fr and fr.get("fisher_q_value") not in (None, "", "NA")
                else 1.0
            )
            f_sig = f_q <= 0.05

            norm_v = max(-1.6, min(1.6, mean_z)) / 1.6
            if norm_v > 0.06:
                bg_col = (1.0, 1.0 - 0.34 * norm_v, 1.0 - 0.38 * norm_v)
            elif norm_v < -0.06:
                bg_col = (1.0 - 0.38 * abs(norm_v), 1.0 - 0.18 * abs(norm_v), 1.0)
            else:
                bg_col = (0.975, 0.975, 0.975)

            edge_col = GOLD_ACCENT if f_sig else (INK if (w_p <= 0.05 or ci_excl) else MID_RULE)
            edge_lw = 2.6 if f_sig else (1.5 if (w_p <= 0.05 or ci_excl) else 0.6)

            rect = mpatches.Rectangle(
                (c_idx - 0.46, r_idx - 0.44),
                0.92,
                0.88,
                facecolor=bg_col,
                edgecolor=edge_col,
                lw=edge_lw,
            )
            ax.add_patch(rect)

            sig_mark = ""
            if f_sig:
                sig_mark += "★"
            if w_q <= 0.05:
                sig_mark += "**"
            elif w_p <= 0.05 or ci_excl:
                sig_mark += "*"

            txt_col = INC_COLOR if mean_z >= 0.20 else (DEC_COLOR if mean_z <= -0.20 else INK)
            ax.text(
                c_idx,
                r_idx,
                f"{mean_z:+.2f}{sig_mark}",
                fontsize=8.3,
                fontweight="bold" if sig_mark else "normal",
                color=txt_col,
                ha="center",
                va="center",
            )

    ax.set_xticks(range(n_c))
    ax.set_xticklabels([t[1] for t in traits_cols], fontsize=8.2, fontweight="bold")
    ax.xaxis.tick_top()
    ax.set_yticks(range(n_r))
    ax.set_yticklabels([c[1] for c in conditions], fontsize=8.4, fontweight="bold")
    ax.set_title(
        "A. Fecal Cohort \N{MULTIPLICATION SIGN} Trait Mean Paired Shift Δz\n★ Gold border = "
        "Fisher BH q≤0.05 · **Wilcoxon BH q≤0.05 · *p≤0.05 / CI≠0",
        fontsize=9.3,
        fontweight="bold",
        loc="left",
        pad=12,
    )
    despine(ax, top=True, right=True, left=True, bottom=True)

    # --- Panel B: Lead Fecal Cohort Forest Plot (Mean Δz ± 95% CI + Study Win/Loss + Fisher q) ---
    # Each label's arrow is the sign of the study-weighted mean it plots (#579)
    lead_specs = [
        (condition, trait, f"{short} · {shown} {lead_arrow(by_eco.get((condition, trait), {}))}")
        for condition, short, trait, shown in (
            ("Parkinson's disease", "PD", "ecological-reporter composite", "Ecological Composite"),
            ("Parkinson's disease", "PD", "moderately halophilic", "Moderately halophilic"),
            ("Parkinson's disease", "PD", "arsenate reduction", "Arsenate reduction"),
            ("Parkinson's disease", "PD", "facultatively aerobic", "Facultatively aerobic"),
            ("Colorectal cancer", "CRC", "ecological-reporter composite", "Ecological Composite"),
            ("Colorectal cancer", "CRC", "moderately halophilic", "Moderately halophilic"),
            ("HIV infection", "HIV", "ecological-reporter composite", "Ecological Composite"),
            ("COVID-19", "COVID-19", "cellulose catabolic process", "Cellulose catabolism"),
            ("Parkinson's disease", "PD", "cellulose catabolic process", "Cellulose catabolism"),
            ("Colorectal cancer", "CRC", "cellulose catabolic process", "Cellulose catabolism"),
            ("HIV infection", "HIV", "cellulose catabolic process", "Cellulose catabolism"),
        )
    ]
    # The sizes its slide quotes (#704)
    PLATE_COUNTS.update(
        {
            "plate4_conditions": len(conditions),
            "plate4_traits": len(traits_cols),
            "plate4_leads": len(lead_specs),
        }
    )
    y_l = np.arange(len(lead_specs))[::-1]
    # Room below the rows for the cellulose callout
    ax_lead.set_ylim(-3.4, len(lead_specs) - 0.15)
    ax_lead.set_xlim(-1.55, 4.45)
    ax_lead.set_xticks([-1.5, -0.5, 0.0, 0.5, 1.5])
    ax_lead.axvline(0.0, color=INK, linestyle="--", linewidth=1.0, zorder=1)
    ax_lead.axhline(y_l[6] - 0.5, color=MID_RULE, linestyle=":", linewidth=1.0, zorder=1)

    for idx_l, (cond_raw, tr_raw, _disp_lbl) in enumerate(lead_specs):
        yp = y_l[idx_l]
        er = by_eco.get((cond_raw, tr_raw), {})
        wr = by_wilc.get((cond_raw, tr_raw), {})
        fr = by_fish.get((cond_raw, tr_raw), {})
        if not er:
            continue
        mean_z = float(er["study_weighted_mean_contrast"])
        ci_lo = float(er["ci_low"]) if er.get("ci_low") not in (None, "", "NA") else mean_z
        ci_hi = float(er["ci_high"]) if er.get("ci_high") not in (None, "", "NA") else mean_z
        pos_st = int(er.get("positive_study_count", 0))
        neg_st = int(er.get("negative_study_count", 0))
        w_p = float(wr["p_value"]) if wr and wr.get("p_value") else 1.0
        w_q = float(wr["q_value"]) if wr and wr.get("q_value") else 1.0
        f_q = (
            float(fr["fisher_q_value"])
            if fr and fr.get("fisher_q_value") not in (None, "", "NA")
            else None
        )

        col = INC_COLOR if mean_z >= 0 else DEC_COLOR
        ax_lead.plot([ci_lo, ci_hi], [yp, yp], color=col, lw=2.3, zorder=3)
        ax_lead.scatter(
            [mean_z],
            [yp],
            color=col,
            edgecolor=GOLD_ACCENT if (f_q and f_q <= 0.05) else col,
            s=48,
            lw=1.6,
            zorder=4,
        )

        stat_parts = [f"{pos_st}↑/{neg_st}↓ st."]
        if f_q is not None and f_q <= 0.05:
            stat_parts.append(f"Fisher q={f_q:.3f}★")
        if w_q <= 0.10:
            stat_parts.append(f"Wilc q={w_q:.3f}")
        elif w_p <= 0.05:
            stat_parts.append(f"Wilc p={w_p:.3f}")

        ax_lead.text(
            4.40,
            yp,
            " · ".join(stat_parts),
            fontsize=7.6,
            fontweight="bold" if (w_q <= 0.10 or (f_q and f_q <= 0.05)) else "normal",
            color=col,
            ha="right",
            va="center",
        )

    # Bottom verification callout badge inside Panel B
    ax_lead.text(
        0.5,
        0.012,
        cellulose["plate4_callout"],
        transform=ax_lead.transAxes,
        fontsize=6.05,
        fontweight="bold",
        color=INK,
        ha="center",
        va="bottom",
        multialignment="left",
        bbox={
            "boxstyle": "round,pad=0.26",
            "facecolor": "#fef9ec",
            "edgecolor": GOLD_ACCENT,
            "lw": 1.0,
        },
        zorder=5,
    )

    ax_lead.set_yticks(y_l)
    ax_lead.set_yticklabels([s[2] for s in lead_specs], fontsize=8.0, fontweight="bold")
    ax_lead.set_xlabel(
        "Study-Weighted Mean Paired Shift Δz (±95% CI across unique studies)",
        fontsize=8.5,
        fontweight="bold",
    )
    ax_lead.set_title(
        "B. Lead Fecal Disease \N{MULTIPLICATION SIGN} Trait Contrasts (PD, CRC, HIV, "
        "COVID-19)\nStudy win/loss counts, Profile Fisher FDR (★), Study Wilcoxon p/q",
        fontsize=8.35,
        fontweight="bold",
        loc="left",
        pad=12,
    )
    despine(ax_lead, left=False)

    fig.text(
        0.03,
        0.025,
        f"{RERUN_SOURCE} · "
        "pd_ecological_reporter_combined_condition_trait_summary.tsv, "
        "pd_condition_taxon_trait_network_condition_trait_tests.tsv, "
        "six_shared_trait_condition_body_site_direction_contrasts.tsv.",
        fontsize=7.8,
        color=MUTED,
    )

    out_path = out_dir / "tufte_plate4_diverging_condition_trait_triangulation.png"
    fig.savefig(out_path, dpi=220)
    fig.savefig(out_path.with_suffix(".svg"))
    plt.close(fig)
    return out_path


# ==============================================================================
# PLATE 5: Anatomical Oxygen-Switch Polarity Inversion & Cross-Disease Switchers
# ==============================================================================
def lead_arrow(row: Mapping[str, str]) -> str:
    """'↑' or '↓' for the sign of a reporter row's study-weighted mean, '' without one."""
    value = row.get("study_weighted_mean_contrast")
    if not has_value(value):
        return ""
    mean = float(str(value))
    return "↑" if mean > 0 else "↓" if mean < 0 else "="


def plate5_strata(
    rows: Sequence[dict[str, str]],
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Plate V's strata: the 10 highest-ranked scoring above 0.20 (rank_score weighs in
    study support) and the 10 lowest-scoring below -0.15. The deck's oxygen_switch_strata
    must select the same; a test ties the two (#639)."""

    def score(row: dict[str, str]) -> float:
        return float(row["aerobe_up_anaerobe_down_score"])

    ranked = sorted(rows, key=lambda r: -float(r["rank_score"]))
    up = [r for r in ranked if score(r) > 0.20][:10]
    down = sorted((r for r in rows if score(r) < -0.15), key=score)[:10]
    return up, down


def crc_oxygen_split(rows: Sequence[dict[str, str]]) -> str:
    """'Feces +0.68 (aerobe↑) vs Lower GI tract -0.79 (anaerobe↑)', or '' without both rows."""
    scores = {
        r["body_site_display"]: float(r["aerobe_up_anaerobe_down_score"])
        for r in rows
        if r["condition"] == "Colorectal cancer"
    }
    if "Feces" not in scores or "Lower GI tract" not in scores:
        return ""

    def part(site: str) -> str:
        value = scores[site]
        return f"{site} {value:+.2f} ({'aerobe' if value > 0 else 'anaerobe'}↑)"

    return f"{part('Feces')} vs {part('Lower GI tract')}"


def build_plate5(repo: Path, out_dir: Path) -> Path:
    cellulose = cellulose_decr_notes(REPO_ROOT)
    mcnemar = McNemar.from_rows(read_tsv(repo / GUILD_CORRELATION))
    oxy_rows = read_tsv(
        repo / "trait-visualization/trait-profile-investigation/taxa-oxygen-switches/"
        "condition_oxygen_switch_scores.tsv"
    )
    sw_rows = read_tsv(
        repo / "trait-visualization/trait-profile-investigation/taxa-oxygen-switches/"
        "taxa_cross_disease_switch_summary.tsv"
    )

    fig = plt.figure(figsize=(16.0, 9.1), facecolor=BG_COLOR)
    gs = fig.add_gridspec(
        1, 2, width_ratios=[1.05, 1.15], wspace=0.16, left=0.05, right=0.98, top=0.84, bottom=0.10
    )

    fig.suptitle(
        "TUFTE PLATE 5 · Oxygen-Switch Scores by Body Site & Cross-Disease Genus Switchers",
        fontsize=15.0,
        fontweight="bold",
        x=0.05,
        ha="left",
        y=0.96,
        color=INK,
    )
    # --- Panel A: Anatomical Oxygen-Switch Polarity Inversion ---
    ax_oxy = fig.add_subplot(gs[0, 0])

    pos_strata, neg_strata = plate5_strata(oxy_rows)
    combined_oxy = pos_strata + neg_strata[::-1]
    combined_oxy.sort(key=lambda r: float(r["aerobe_up_anaerobe_down_score"]))

    fecal = [r for r in combined_oxy if r["body_site_display"] == "Feces"]
    other = [r for r in combined_oxy if r["body_site_display"] != "Feces"]
    fecal_up = sum(float(r["aerobe_up_anaerobe_down_score"]) > 0 for r in fecal)
    other_down = sum(float(r["aerobe_up_anaerobe_down_score"]) < 0 for r in other)
    fig.text(
        0.05,
        0.922,
        f"Of the {len(combined_oxy)} strata shown, {fecal_up} of {len(fecal)} fecal strata score "
        f"aerobe-up and {other_down} of {len(other)} others anaerobe-up; major genera switch "
        "direction across diseases. st. = the larger of the aerobe-only and anaerobe-only study "
        "counts.",
        fontsize=8.8,
        color=MUTED,
    )

    ax_oxy.set_ylim(-0.7, len(combined_oxy) - 0.2)
    ax_oxy.set_xlim(-1.45, 1.50)
    ax_oxy.axvline(0.0, color=INK, lw=1.0, zorder=2)

    for idx_r, r in enumerate(combined_oxy):
        score = float(r["aerobe_up_anaerobe_down_score"])
        cond = r["condition"]
        cond_disp = "COPD" if cond == "Chronic obstructive pulmonary disease" else cond
        site = r["body_site_display"]
        st_a = int(float(r["aerobic_only_study_count"]))
        st_an = int(float(r["anaerobic_only_study_count"]))
        is_crc = cond == "Colorectal cancer"

        col = INC_COLOR if score > 0 else DEC_COLOR
        ax_oxy.barh(
            idx_r,
            score,
            height=0.58,
            color=col,
            alpha=0.85 if not is_crc else 1.0,
            edgecolor=GOLD_ACCENT if is_crc else "none",
            lw=1.8 if is_crc else 0,
        )

        lbl = f"{cond_disp} · {site}"
        if is_crc:
            lbl = f"★ {lbl}"

        if score >= 0:
            ax_oxy.text(
                -0.03,
                idx_r,
                lbl,
                fontsize=7.9,
                fontweight="bold" if is_crc else "normal",
                color=INK,
                ha="right",
                va="center",
            )
            ax_oxy.text(
                score + 0.03,
                idx_r,
                f"{score:+.2f} ({max(st_a, st_an)} st.)",
                fontsize=7.5,
                fontweight="bold",
                color=col,
                ha="left",
                va="center",
            )
        else:
            ax_oxy.text(
                0.03,
                idx_r,
                lbl,
                fontsize=7.9,
                fontweight="bold" if is_crc else "normal",
                color=INK,
                ha="left",
                va="center",
            )
            ax_oxy.text(
                score - 0.03,
                idx_r,
                f"{score:+.2f} ({max(st_a, st_an)} st.)",
                fontsize=7.5,
                fontweight="bold",
                color=col,
                ha="right",
                va="center",
            )

    ax_oxy.set_yticks([])
    ax_oxy.set_xlabel(
        "Oxygen-Switch Score:  ← Anaerobe-Up / Aerobe-Down   |   Aerobe-Up / Anaerobe-Down →",
        fontsize=8.4,
        fontweight="bold",
    )
    crc = crc_oxygen_split(oxy_rows)
    ax_oxy.set_title(
        f"A. Oxygen-Switch Score: {len(pos_strata)} Top-Ranked Aerobe-Up & {len(neg_strata)} "
        "Lowest-Scoring Strata" + (f"\n★ CRC: {crc}" if crc else ""),
        fontsize=9.1,
        fontweight="bold",
        loc="left",
        pad=10,
    )
    despine(ax_oxy, left=True)

    # --- Panel B: Cross-Disease Genus Switchers (Top 12 with generous spacing) ---
    ax_sw = fig.add_subplot(gs[0, 1])
    top_sw = [r for r in sw_rows if r["taxonomic_rank"] == "genus"][:12]
    y_sw = np.arange(len(top_sw))[::-1]

    ax_sw.set_ylim(-0.6, len(top_sw) - 0.1)
    ax_sw.set_xlim(-185, 235)
    ax_sw.axvline(0, color=INK, lw=1.0, zorder=2)

    def short_cond_list(raw: str, max_items: int = 2) -> str:
        parts = [p.strip() for p in raw.split(";") if p.strip()][:max_items]
        cleaned = []
        for p in parts:
            name = p.split("(")[0].strip()
            cnt = p.split("(")[1].split(" ")[0] if "(" in p else ""
            short_map = {
                "Parkinson's disease": "PD",
                "Colorectal cancer": "CRC",
                "Irritable bowel syndrome": "IBS",
                "HIV infection": "HIV",
                "Human immunodeficiency virus": "HIV",
                "Chronic constipation": "Constip.",
                "Gestational diabetes": "GDM",
                "Type II diabetes mellitus": "T2DM",
                "Cervical glandular intraepithelial neoplasia": "CGIN",
                "Human papilloma virus infection": "HPV",
                "Rheumatoid arthritis": "RA",
                "Gastric cancer": "Gastric Ca",
                "Crohn's disease": "Crohn's",
            }
            name = short_map.get(name, name)
            cleaned.append(f"{name} ({cnt})")
        return ", ".join(cleaned)

    for idx_r, r in enumerate(top_sw):
        yp = y_sw[idx_r]
        gname = r["taxon_name"]
        inc_p = int(r["increased_side_profiles"])
        dec_p = int(r["decreased_side_profiles"])
        top_inc = short_cond_list(r["top_increased_conditions"], 2)
        top_dec = short_cond_list(r["top_decreased_conditions"], 2)

        ax_sw.barh(yp - 0.14, -dec_p, height=0.36, color=DEC_COLOR, alpha=0.84, zorder=3)
        ax_sw.barh(yp - 0.14, inc_p, height=0.36, color=INC_COLOR, alpha=0.84, zorder=3)

        ax_sw.text(
            0,
            yp + 0.10,
            f"{gname}  [switch={float(r['switch_score']):.0f} · {dec_p}↓ vs {inc_p}↑]",
            fontsize=7.5,
            fontweight="bold",
            color=INK,
            ha="center",
            va="bottom",
            bbox={
                "boxstyle": "round,pad=0.08",
                "facecolor": BG_COLOR,
                "edgecolor": "none",
                "alpha": 0.92,
            },
            zorder=5,
        )
        ax_sw.text(
            -dec_p - 4,
            yp - 0.14,
            f"↓ {top_dec}",
            fontsize=7.2,
            color=DEC_COLOR,
            ha="right",
            va="center",
        )
        ax_sw.text(
            inc_p + 4,
            yp - 0.14,
            f"↑ {top_inc}",
            fontsize=7.2,
            color=INC_COLOR,
            ha="left",
            va="center",
        )

    ax_sw.set_yticks([])
    ax_sw.set_xlabel(
        "← Decreased-Signature Profiles (top conditions)   |   Increased-Signature Profiles (top "
        "conditions) →",
        fontsize=8.4,
        fontweight="bold",
    )
    ax_sw.set_title(
        f"B. Top {len(top_sw)} Cross-Disease Switcher Genera (Context-Dependent Polarity)\n"
        f"Note PD split: {pd_split(top_sw)}",
        fontsize=9.5,
        fontweight="bold",
        loc="left",
        pad=10,
    )
    despine(ax_sw, left=True)

    fig.text(
        0.05,
        0.025,
        f"{RERUN_SOURCE} · condition_oxygen_switch_scores.tsv, "
        "taxa_cross_disease_switch_summary.tsv. Pooled McNemar OR = "
        f"{mcnemar.odds_ratio:.3f} ({mcnemar.aerobe_up_anaerobe_down} vs "
        f"{mcnemar.anaerobe_up_aerobe_down} st., p = {fmt_p(mcnemar.p_value)}); "
        f"{cellulose['crc_source']}",
        fontsize=7.8,
        color=MUTED,
    )

    out_path = out_dir / "tufte_plate5_anatomical_oxygen_polarity_and_switchers.png"
    fig.savefig(out_path, dpi=220)
    fig.savefig(out_path.with_suffix(".svg"))
    plt.close(fig)
    return out_path


# ==============================================================================
# PLATE 6: Taxonomic Resolution vs Trait Coverage & Root-Cause Diagnostic
# ==============================================================================
def build_plate6(repo: Path, out_dir: Path) -> Path:
    rank_rows = read_tsv(
        repo
        / "trait-visualization/trait-profile-investigation/taxa-coverage/taxa_coverage_by_rank.tsv"
    )
    zero_rows = read_tsv(
        repo / "trait-visualization/trait-profile-investigation/taxa-coverage/"
        "top_zero_trait_bugsigdb_taxa.tsv"
    )

    fig = plt.figure(figsize=(16.0, 9.0), facecolor=BG_COLOR)
    gs = fig.add_gridspec(
        2,
        2,
        width_ratios=[0.90, 1.30],
        height_ratios=[1.05, 0.95],
        wspace=0.48,
        hspace=0.35,
        left=0.13,
        right=0.98,
        top=0.84,
        bottom=0.11,
    )

    ax_rank = fig.add_subplot(gs[0, 0])
    ax_cat = fig.add_subplot(gs[1, 0])
    ax_top = fig.add_subplot(gs[:, 1])

    fig.suptitle(
        "TUFTE PLATE 6 · Taxonomic Resolution vs KG-Microbe Trait Coverage & Zero-Trait Root-Cause "
        "Audit",
        fontsize=15.0,
        fontweight="bold",
        x=0.04,
        ha="left",
        y=0.96,
        color=INK,
    )
    fig.text(
        0.04,
        0.922,
        "Why Genus/Family signatures achieve 96.5% trait coverage while 45.9% of distinct Species "
        "(22.0% of species occurrences) have zero traits in KG-Microbe — and how TaxID remapping "
        "and LCA propagation recover them.",
        fontsize=8.7,
        color=MUTED,
    )

    # --- Panel A: Rank-Stratified Coverage Cliff ---
    ranks_keep = [r for r in rank_rows if r["taxonomic_rank"] != "superkingdom"]
    y_rk = np.arange(len(ranks_keep))[::-1]
    ax_rank.set_ylim(-0.5, len(ranks_keep) - 0.3)
    ax_rank.set_xlim(38, 132)
    ax_rank.set_xticks([40, 60, 80, 100])

    for idx_r, r in enumerate(ranks_keep):
        yp = y_rk[idx_r]
        u_cov = float(r["trait_bearing_pct"])
        o_cov = float(r["trait_bearing_occurrences_pct"])

        ax_rank.plot([u_cov, o_cov], [yp, yp], color=INK, lw=2.0, zorder=2)
        ax_rank.scatter(
            [u_cov],
            [yp],
            color=DEC_COLOR,
            s=46,
            zorder=3,
            label="Unique Taxa Covered (%)" if idx_r == 0 else "",
        )
        ax_rank.scatter(
            [o_cov],
            [yp],
            color=TEAL_ACCENT,
            s=52,
            marker="D",
            zorder=3,
            label="Signature Occurrences Covered (%)" if idx_r == 0 else "",
        )

        ax_rank.text(
            u_cov - 1.5,
            yp,
            f"{u_cov:.1f}%",
            fontsize=7.5,
            fontweight="bold",
            color=DEC_COLOR,
            ha="right",
            va="center",
        )
        ax_rank.text(
            o_cov + 1.5,
            yp,
            f"{o_cov:.1f}% "
            f"({int(r['trait_bearing_occurrences']):,}/{int(r['total_occurrences']):,})",
            fontsize=7.4,
            fontweight="bold",
            color=TEAL_ACCENT,
            ha="left",
            va="center",
        )

    ax_rank.set_yticks(y_rk)
    ax_rank.set_yticklabels(
        [f"{r['display_rank']}\n({int(r['total_taxa']):,} taxa)" for r in ranks_keep],
        fontsize=8.0,
        fontweight="bold",
    )
    ax_rank.set_xlabel("% with ≥1 KG-Microbe Trait Annotation", fontsize=8.2, fontweight="bold")
    ax_rank.set_title(
        "A. The Species-Level Coverage Cliff\nSpecies: 54.1% unique / 78.0% occurrences vs Genus: "
        "76.0% / 96.5%",
        fontsize=9.4,
        fontweight="bold",
        loc="left",
        pad=8,
    )
    ax_rank.legend(
        loc="lower left", frameon=True, facecolor=PANEL_BG, edgecolor=MID_RULE, fontsize=7.2
    )
    despine(ax_rank)

    # --- Panel B: Root-Cause Pareto Breakdown of Top 100 Zero-Trait BugSigDB Taxa ---
    cat_totals: dict[str, int] = {}
    cat_taxa_cnt: dict[str, int] = {}
    for r in zero_rows:
        cat = r["diagnosis_category"]
        occ = int(r["signature_occurrences"])
        cat_totals[cat] = cat_totals.get(cat, 0) + occ
        cat_taxa_cnt[cat] = cat_taxa_cnt.get(cat, 0) + 1

    sorted_cats = sorted(cat_totals.items(), key=lambda x: x[1], reverse=True)
    y_c = np.arange(len(sorted_cats))[::-1]
    cat_colors = {
        "Uncultivated / Metagenomic MAG": "#3b7ea1",
        "Rare / Novel Isolate (Lacks KG Phenotype)": "#6c7a89",
        "Ambiguous Taxon / Bin Slash": "#d96b27",
        "Non-Bacterial / Fungal Eukaryote": "#7b4f9d",
        "NCBI TaxID Desynchronization": "#b8332a",
        "Candidatus Lineage": "#157a6e",
    }
    cat_short = {
        "Uncultivated / Metagenomic MAG": "Uncultivated / MAG",
        "Rare / Novel Isolate (Lacks KG Phenotype)": "Rare / Novel Isolate",
        "Ambiguous Taxon / Bin Slash": "Ambiguous Slash Bin",
        "Non-Bacterial / Fungal Eukaryote": "Fungal Eukaryote",
        "NCBI TaxID Desynchronization": "NCBI TaxID Desync",
        "Candidatus Lineage": "Candidatus Lineage",
    }

    ax_cat.set_ylim(-0.5, len(sorted_cats) - 0.3)
    max_cat_occ = max(v for _, v in sorted_cats)
    ax_cat.set_xlim(0, max_cat_occ * 1.45)

    for idx_c, (cat_name, tot_occ) in enumerate(sorted_cats):
        yp = y_c[idx_c]
        col = cat_colors.get(cat_name, TEAL_ACCENT)
        ax_cat.barh(yp, tot_occ, height=0.58, color=col, alpha=0.88)
        ax_cat.text(
            tot_occ + 10,
            yp,
            f"{tot_occ:,} occ. ({cat_taxa_cnt[cat_name]} taxa)",
            fontsize=7.7,
            fontweight="bold",
            color=col,
            va="center",
        )

    ax_cat.set_yticks(y_c)
    ax_cat.set_yticklabels(
        [cat_short.get(c[0], c[0]) for c in sorted_cats], fontsize=8.0, fontweight="bold"
    )
    ax_cat.set_xlabel(
        "Zero-Trait Signature Occurrences (Top 100 Zero-Trait Taxa)",
        fontsize=8.2,
        fontweight="bold",
    )
    ax_cat.set_title(
        "B. Root-Cause Taxonomy of Top 100 Zero-Trait BugSigDB Taxa\nActionable fixes: TaxID "
        "remapping, slash-bin LCA & family rollup",
        fontsize=9.3,
        fontweight="bold",
        loc="left",
        pad=8,
    )
    despine(ax_cat)

    # --- Panel C: Top 18 Zero-Trait Taxa Diagnostic Bars (wider xlim so zero clipping) ---
    top18 = zero_rows[:18]
    y_t = np.arange(len(top18))[::-1]
    ax_top.set_ylim(-0.6, len(top18) - 0.2)
    ax_top.set_xlim(0, 385)
    ax_top.set_xticks([0, 50, 100, 150, 200])

    for idx_t, r in enumerate(top18):
        yp = y_t[idx_t]
        occ = int(r["signature_occurrences"])
        st = int(r["study_count"])
        cat = r["diagnosis_category"]
        col = cat_colors.get(cat, INK)
        rem = r["suggested_remediation"]
        if "NCBITaxon:" in rem:
            rem_short = "Remap → " + rem.split("traits:")[-1].strip()
        elif "lowest common ancestor" in rem:
            rem_short = "Multi-taxon LCA resolution"
        elif "Propagate higher-rank" in rem:
            rem_short = "Propagate family/order traits"
        elif "FungiTraits" in rem:
            rem_short = "Fungal eukaryote (non-bacterial)"
        elif "Candidatus" in cat:
            rem_short = "MAG metabolic model inference"
        else:
            rem_short = "Curate in BacDive / KG-Microbe"

        ax_top.barh(yp, occ, height=0.58, color=col, alpha=0.88)
        ax_top.text(
            occ + 3,
            yp,
            f"{occ} occ ({st} st.) · {rem_short}",
            fontsize=7.4,
            fontweight="bold" if cat == "NCBI TaxID Desynchronization" or occ >= 60 else "normal",
            color=col if cat == "NCBI TaxID Desynchronization" else INK,
            va="center",
        )

    ax_top.set_yticks(y_t)
    ax_top.set_yticklabels(
        [f"{r['taxon_name'][:26]} [{r['taxon_id'].split(':')[-1]}]" for r in top18],
        fontsize=7.6,
        fontweight="bold",
    )
    ax_top.set_xlabel(
        "BugSigDB Signature Occurrences with Zero KG-Microbe Traits",
        fontsize=8.2,
        fontweight="bold",
    )
    ax_top.set_title(
        "C. Top 18 Zero-Trait BugSigDB Taxa & Remediation Path\nHighlighting Oscillospiraceae (60 "
        "st.) & Rothia (23 st.) TaxID desyncs",
        fontsize=9.3,
        fontweight="bold",
        loc="left",
        pad=10,
    )
    despine(ax_top)

    fig.text(
        0.04,
        0.025,
        f"{RERUN_SOURCE} · taxa_coverage_by_rank.tsv, "
        "top_zero_trait_bugsigdb_taxa.tsv (create_taxa_coverage_histogram.py). Remapping "
        "NCBITaxon:541000 → 216572 recovers 162 traits across 60 studies; NCBITaxon:508215 → 32207 "
        "recovers 87 traits across 23 studies.",
        fontsize=7.7,
        color=MUTED,
    )

    out_path = out_dir / "tufte_plate6_taxonomic_coverage_blindspots.png"
    fig.savefig(out_path, dpi=220)
    fig.savefig(out_path.with_suffix(".svg"))
    plt.close(fig)
    return out_path


# ==============================================================================
# PLATE 7: Three-Arm 317-Trait Concordance & Expanded Table 1 Forest Plot
# ==============================================================================
def three_arm_traits(
    arm1_rows: Mapping[str, dict[str, str]],
    arm2_rows: Mapping[str, dict[str, str]],
    arm3_screen_rows: Mapping[str, dict[str, str]],
    t1_ids: set[str],
) -> list[dict[str, Any]]:
    """The traits Plate VII draws, and Methods Plate II counts (#640): the Table 1 traits
    plus every trait passing Arm 1 study BH and Arm 3 BH, with each arm's direction."""
    # The traits passing both Arm 1 Study FDR <= 0.05 and Arm 3 Study BH FDR <= 0.05, which
    # include traits whose arms disagree (paired-down / ORA-inverted), plus the Table 1 traits
    selected_traits = []
    for tid, s_row in arm3_screen_rows.items():
        a1 = arm1_rows.get(tid)
        if not a1:
            continue
        a1_st_q = (
            float(a1["study_q_value"]) if a1.get("study_q_value") not in (None, "", "NA") else 1.0
        )
        a1_diff = (
            float(a1["mean_increased_minus_decreased"])
            if a1.get("mean_increased_minus_decreased") not in (None, "", "NA")
            else 0.0
        )
        a3_bh_q = (
            float(s_row["bh_q_nominal"])
            if s_row.get("bh_q_nominal") not in (None, "", "NA")
            else 1.0
        )
        a3_by_q = (
            float(s_row["by_q_nominal"])
            if s_row.get("by_q_nominal") not in (None, "", "NA")
            else 1.0
        )
        if (a1_st_q <= 0.05 and a3_bh_q <= 0.05) or (tid in t1_ids):
            raw_lbl = a1["trait_label"]
            lbl = re.sub(
                r"\s+(?:METPO:\S+|CHEBI:\S+|GO:\S+|EC:\S+|[a-z0-9_]+_[a-z0-9_]+)$", "", raw_lbl
            )
            lbl = (
                lbl.replace("metabolite utilization: ", "Metabolite: ")
                .replace("enzyme activity: ", "Enzyme: ")
                .replace("process: ", "")
            )
            a2 = arm2_rows.get(tid) or arm2_rows.get(s_row.get("term_id", ""), {})
            inc_h = int(a2.get("inc_hits", 0)) if a2 else 0
            dec_h = int(a2.get("dec_hits", 0)) if a2 else 0
            tot_h = int(a2.get("hit_rows", 0)) if a2 else 0
            up_st = int(s_row["studies_up"])
            dn_st = int(s_row["studies_down"])
            n_st = up_st + dn_st
            sh_up = (
                float(s_row["share_up"]) if s_row.get("share_up") not in (None, "", "NA") else 0.5
            )

            paired_dir = "inc" if a1_diff > 0 else "dec"
            arm3_dir = "inc" if sh_up > 0.5 else "dec"
            arm2_dir = "inc" if inc_h > dec_h else ("dec" if dec_h > inc_h else "tie")
            all3_concordant = (paired_dir == arm3_dir == arm2_dir) and tot_h > 0
            # Arms 1 and 3 agree and ORA points the other way, or the arms differ otherwise:
            # no ORA call, an ORA tie, or Arms 1 and 3 disagreeing (#685)
            if tid in t1_ids:
                group = "table1"
            elif all3_concordant:
                group = "concordant"
            elif paired_dir == arm3_dir and arm2_dir not in (paired_dir, "tie"):
                group = "ora_inverted"
            else:
                group = "not_concordant"

            selected_traits.append(
                {
                    "tid": tid,
                    "label": lbl,
                    "in_t1": tid in t1_ids,
                    "a1_st_q": a1_st_q,
                    "a1_diff": a1_diff,
                    "inc_h": inc_h,
                    "dec_h": dec_h,
                    "up_st": up_st,
                    "dn_st": dn_st,
                    "n_st": n_st,
                    "sh_up": sh_up,
                    "a3_bh_q": a3_bh_q,
                    "a3_by_q": a3_by_q,
                    "paired_dir": paired_dir,
                    "arm2_dir": arm2_dir,
                    "all3_concordant": all3_concordant,
                    "group": group,
                }
            )

    return selected_traits


# Plate VII's markers, label badges and notes for its trait groups (#685)
MARKERS = {"table1": "D", "concordant": "o", "ora_inverted": "^", "not_concordant": "s"}
PREFIXES = {"table1": "★ ", "concordant": "", "ora_inverted": "▲ ", "not_concordant": "■ "}
NOTES = {"table1": "", "concordant": "", "ora_inverted": " (ORA inv.)", "not_concordant": ""}
CELLULOSE_ID = "METPO:2000103|GO:0030245"


def three_arm_counts(selected: Sequence[dict[str, Any]]) -> tuple[int, int]:
    """Traits passing Arm 1 study BH and Arm 3 BH with all three arms concordant, and how
    many of those also pass Arm 3 BY."""
    concordant = [
        t
        for t in selected
        if t["all3_concordant"] and t["a1_st_q"] <= 0.05 and t["a3_bh_q"] <= 0.05
    ]
    return len(concordant), sum(t["a3_by_q"] <= 0.05 for t in concordant)


def plate_seven_groups(selected: Sequence[dict[str, Any]]) -> dict[str, int]:
    """The Table 1, additional concordant, ORA-inverted and otherwise discordant traits
    Plate VII draws (#579, #685)."""
    concordant = [
        t
        for t in selected
        if t["all3_concordant"] and t["a1_st_q"] <= 0.05 and t["a3_bh_q"] <= 0.05
    ]
    return {
        "table1_traits": sum(t["in_t1"] for t in selected),
        "table1_concordant": sum(t["in_t1"] for t in concordant),
        "three_arm_additional": sum(not t["in_t1"] for t in concordant),
        "ora_inverted": sum(t["group"] == "ora_inverted" for t in selected),
        "not_concordant": sum(t["group"] == "not_concordant" for t in selected),
        "forest_traits": len(selected),
    }


def build_plate7(repo: Path, out_dir: Path) -> Path:
    arm1_rows = {
        r["trait"]: r
        for r in read_tsv(
            repo / "trait-visualization/trait-profile-investigation/"
            "trait_profile_paired_contrast_trait_stats_pooled.tsv"
        )
    }
    arm2_rows = {
        r["trait"]: r
        for r in read_tsv(repo / "data/derived/trait_enrichment_traits_overall_ranked.tsv")
    }
    arm3_screen_rows = {
        r["trait"]: r
        for r in read_tsv(repo / "data/derived/trait_screen/trait_direction_screen.tsv")
    }
    orient_rows = read_tsv(repo / COMPARISON)
    signs = table1_signs(repo)
    every, oriented = signs[ALL_STUDIES], signs[ORIENTED]
    cellulose, cellulose_all = oriented.fractions[CELLULOSE], every.fractions[CELLULOSE]
    orient_t1 = {
        r["numerator"]: r
        for r in orient_rows
        if r.get("comparison") == "fraction"
        and r.get("population") == "disease_vs_healthy_oriented"
    }
    t1_metric_map = {
        "METPO:1000608": "facultatively_aerobic",
        "METPO:2000103|GO:0030245": "cellulose_catabolic_process",
        "METPO:1000668": "coccus_shaped",
        "METPO:1000683": "sphere_shaped",
        "METPO:2000103|pathways:arsenate_reduction": "arsenate_reduction",
        "METPO:1000623": "moderately_halophilic",
    }

    t1_ids = set(t1_metric_map.keys())
    selected_traits = three_arm_traits(arm1_rows, arm2_rows, arm3_screen_rows, t1_ids)

    # Sort: 1) Original 6 Table 1 traits (by share_up desc), 2) Additional 3-Arm Concordant
    # traits (by share_up desc), 3) ORA-inverted traits
    order = ("table1", "concordant", "ora_inverted", "not_concordant")
    selected_traits.sort(key=lambda x: (order.index(x["group"]), -x["sh_up"]))
    group_of = {t["tid"]: t["group"] for t in selected_traits}
    arm1_studies = max(
        int(float(r["study_count"])) for r in arm1_rows.values() if r.get("study_count")
    )
    # The groups the plate and its slide count (#579, #640)
    concordant, _ = three_arm_counts(selected_traits)
    groups = plate_seven_groups(selected_traits)
    PLATE_COUNTS.update({**groups, "screened_traits": len(arm3_screen_rows)})

    fig = plt.figure(figsize=(16.0, 9.0), facecolor=BG_COLOR)
    gs = fig.add_gridspec(
        1, 2, width_ratios=[0.88, 1.32], wspace=0.44, left=0.05, right=0.985, top=0.84, bottom=0.11
    )
    ax_scat = fig.add_subplot(gs[0, 0])
    ax_for = fig.add_subplot(gs[0, 1])

    fig.suptitle(
        "TUFTE PLATE 7 · Three-Arm Concordance & Expanded Table 1 Forest Plot Across All "
        f"{len(arm3_screen_rows)} Traits",
        fontsize=15.0,
        fontweight="bold",
        x=0.05,
        y=0.955,
        ha="left",
        color=INK,
    )
    fig.text(
        0.05,
        0.916,
        f"Full {len(arm3_screen_rows)}-trait BH/BY screen (correcting {screening_bias(repo)} "
        f"post-selection bias) yields {concordant} three-arm concordant traits; Cellulose is "
        f"{cellulose.share_down:.1%}↓ in the {oriented.tested_studies:,} oriented studies (p = "
        f"{fmt_p(cellulose.p_value)}, BH₆ q = {fmt_p(cellulose.q_value)}) and "
        f"{cellulose_all.share_down:.1%}↓ in all {every.tested_studies:,}.",
        fontsize=8.9,
        color=MUTED,
        ha="left",
    )

    # --- Panel A (Left): Cross-Arm Concordance Scatter Plot Across Every Tested Trait ---
    drawn = 0
    for tid, s_row in arm3_screen_rows.items():
        a1 = arm1_rows.get(tid)
        if (
            not a1
            or a1.get("study_t_statistic") in (None, "", "NA")
            or s_row.get("share_up") in (None, "", "NA")
        ):
            continue
        t_st = float(a1["study_t_statistic"])
        sh_u = (float(s_row["share_up"]) - 0.5) * 100.0
        a1_sq = (
            float(a1["study_q_value"]) if a1.get("study_q_value") not in (None, "", "NA") else 1.0
        )
        a3_bq = (
            float(s_row["bh_q_nominal"])
            if s_row.get("bh_q_nominal") not in (None, "", "NA")
            else 1.0
        )

        if tid in t1_ids:
            drawn += 1
            col = INC_COLOR if t_st > 0 else DEC_COLOR
            ax_scat.scatter(
                t_st, sh_u, s=90, color=col, edgecolor=INK, linewidth=1.3, zorder=5, marker="D"
            )
            m_id = t1_metric_map[tid]
            or_r = orient_t1.get(m_id, {})
            if or_r and or_r.get("share_up") not in (None, "", "NA"):
                or_sh_u = (float(or_r["share_up"]) - 0.5) * 100.0
                ax_scat.scatter(
                    t_st,
                    or_sh_u,
                    s=42,
                    facecolor="#ffffff",
                    edgecolor=INK,
                    linewidth=1.4,
                    zorder=6,
                    marker="o",
                )
        elif a1_sq <= 0.05 and a3_bq <= 0.05:
            drawn += 1
            is_inv = group_of.get(tid) in ("ora_inverted", "not_concordant")
            col = GOLD_ACCENT if is_inv else (TEAL_ACCENT if t_st > 0 else PURPLE_ACCENT)
            mk = MARKERS[group_of.get(tid, "concordant")]
            ax_scat.scatter(
                t_st,
                sh_u,
                s=62 if is_inv else 54,
                color=col,
                edgecolor=INK if is_inv else "white",
                linewidth=0.8,
                zorder=4,
                marker=mk,
            )
        else:
            drawn += 1
            ax_scat.scatter(t_st, sh_u, s=18, color="#cbd5e1", alpha=0.65, zorder=2)

    ax_scat.set_xlim(-5.8, 7.2)
    ax_scat.set_ylim(-21.5, 22.0)
    ax_scat.axhline(0, color=MID_RULE, lw=0.9, ls="--")
    ax_scat.axvline(0, color=MID_RULE, lw=0.9, ls="--")

    # Where each label sits; its badge, arrow and note come from the tables (#685)
    annot_specs = [
        ("METPO:2000103|GO:0046274", "Lignin catabolism", (-5.5, 18.2)),
        ("METPO:1000608", "Facultatively aerobic", (-5.5, 14.2)),
        ("METPO:1000623", "Moderately halophilic", (-5.5, 10.6)),
        ("METPO:2000103|pathways:arsenate_reduction", "Arsenate reduction", (-5.5, 7.0)),
        ("METPO:1000683", "Sphere shaped", (-5.5, 3.6)),
        ("METPO:1000668", "Coccus shaped", (1.6, 1.4)),
        ("METPO:1000429", "GC low", (-5.6, -4.5)),
        ("METPO:2000103|GO:0030245", "Cellulose", (-5.6, -13.3)),
        ("METPO:1000603", "Anaerobic", (-5.6, -16.0)),
    ]
    traits_by_id = {t["tid"]: t for t in selected_traits}
    for tid, name, (tx_pos, ty_pos) in annot_specs:
        if tid in arm1_rows and tid in arm3_screen_rows:
            paired = float(arm1_rows[tid]["study_t_statistic"]) > 0
            shared = float(arm3_screen_rows[tid]["share_up"]) > 0.5
            text_lbl = f"{name} {'↑' if paired else '↓'}"
            if shared != paired:
                text_lbl += f" (Arm 3 {'↑' if shared else '↓'})"
            trait = traits_by_id.get(tid)
            if trait is not None:
                text_lbl = f"{PREFIXES[trait['group']]}{text_lbl}{NOTES[trait['group']]}"
                if trait["group"] == "concordant" and trait["a3_by_q"] <= 0.05:
                    text_lbl += " (3-arm BY ✓)"
            if tid == CELLULOSE_ID:
                text_lbl += (
                    f" (Or.{oriented.tested_studies}: {(cellulose.share_up - 0.5) * 100:+.1f}%, "
                    f"BH₆ q={fmt_p(cellulose.q_value)})"
                )
            tx = float(arm1_rows[tid]["study_t_statistic"])
            ty = (float(arm3_screen_rows[tid]["share_up"]) - 0.5) * 100.0
            ax_scat.annotate(
                text_lbl,
                xy=(tx, ty),
                xytext=(tx_pos, ty_pos),
                fontsize=7.2,
                fontweight="bold" if "★" in text_lbl or "▲" in text_lbl else "normal",
                color=INK,
                bbox={
                    "boxstyle": "round,pad=0.14",
                    "facecolor": BG_COLOR,
                    "edgecolor": "none",
                    "alpha": 0.88,
                },
                arrowprops={"arrowstyle": "-", "color": MID_RULE, "lw": 0.7},
                zorder=6,
            )

    # Compact visual legend badge in bottom-right of scatter plot
    ax_scat.text(
        0.97,
        0.04,
        f"◆ ★ Original 6 Table 1 traits ({groups['table1_concordant']} with all 3 arms ✓)\n"
        f"○ Oriented {oriented.tested_studies}-st. Table 1 (Cell. ↓ p={fmt_p(cellulose.p_value)})\n"
        f"● {groups['three_arm_additional']} Additional 3-Arm Study-FDR traits (✓)\n"
        f"▲ {groups['ora_inverted']} with Arms 1 & 3 agreeing, ORA opposite"
        + (
            f"\n■ {groups['not_concordant']} with the arms differing otherwise"
            if groups["not_concordant"]
            else ""
        ),
        transform=ax_scat.transAxes,
        fontsize=7.0,
        fontweight="bold",
        va="bottom",
        ha="right",
        color=INK,
        bbox={
            "boxstyle": "round,pad=0.35",
            "facecolor": "#fef9ec",
            "edgecolor": GOLD_ACCENT,
            "lw": 1.0,
        },
        zorder=7,
    )

    ax_scat.set_xlabel(
        f"Arm 1 Study-Level Paired Contrast t-Statistic (N = {arm1_studies:,} studies)",
        fontsize=8.6,
        fontweight="bold",
    )
    ax_scat.set_ylabel(
        "Arm 3 Study Sign-Test Shift (% Studies Up \N{MINUS SIGN} 50%)",
        fontsize=8.6,
        fontweight="bold",
    )
    ax_scat.set_title(
        f"A. Cross-Arm Alignment Across All {drawn} Tested Traits\nStrong concordance between "
        "continuous Δz (Arm 1) and study sign-test share (Arm 3)",
        fontsize=9.3,
        fontweight="bold",
        loc="left",
        pad=10,
    )
    despine(ax_scat)

    # --- Panel B (Right): Forest Plot of All 19 Study-FDR Multi-Arm Traits ---
    n_tr = len(selected_traits)
    y_pos = np.arange(n_tr)[::-1]
    ax_for.set_ylim(-0.7, n_tr - 0.1)
    ax_for.set_xlim(30.0, 117.0)
    ax_for.axvline(50.0, color=INK, lw=1.1, ls="--", zorder=2)

    # A rule above the first trait of each group after the Table 1 traits
    boundaries = [
        idx
        for idx in range(1, n_tr)
        if selected_traits[idx]["group"] != selected_traits[idx - 1]["group"]
    ]
    y_labels_for = []
    for idx, tr in enumerate(selected_traits):
        yp = y_pos[idx]
        pct = tr["sh_up"] * 100.0
        n_s = max(tr["n_st"], 1)
        se_pct = math.sqrt(tr["sh_up"] * (1.0 - tr["sh_up"]) / n_s) * 100.0
        ci_lo = max(31.0, pct - 1.96 * se_pct)
        ci_hi = min(73.0, pct + 1.96 * se_pct)

        if idx in boundaries:
            first = idx == boundaries[0]
            ax_for.axhline(yp + 0.5, color=GOLD_ACCENT if first else MID_RULE, lw=1.0, ls=":")

        if tr["in_t1"]:
            ax_for.axhspan(yp - 0.44, yp + 0.44, color="#fef9ec", alpha=0.85, zorder=0)

        col = INC_COLOR if tr["paired_dir"] == "inc" else DEC_COLOR
        if not tr["all3_concordant"]:
            col = GOLD_ACCENT

        ax_for.plot([ci_lo, ci_hi], [yp, yp], color=col, lw=2.2, zorder=3)
        mk = MARKERS[tr["group"]]
        ax_for.scatter(
            [pct],
            [yp],
            color=col,
            s=55 if tr["in_t1"] else 44,
            marker=mk,
            edgecolor=INK if tr["in_t1"] else "white",
            lw=0.9,
            zorder=4,
        )

        or_badge = ""
        if tr["in_t1"]:
            m_id = t1_metric_map[tr["tid"]]
            or_r = orient_t1.get(m_id, {})
            if or_r and or_r.get("share_up") not in (None, "", "NA"):
                or_pct = float(or_r["share_up"]) * 100.0
                or_p = float(or_r["p_value"])
                ax_for.scatter(
                    [or_pct],
                    [yp],
                    marker="o",
                    s=38,
                    facecolor="#ffffff",
                    edgecolor=INK,
                    lw=1.4,
                    zorder=5,
                )
                or_badge = f" · Or.{oriented.tested_studies}:{or_pct:.1f}% (p={or_p:.1e})"

        prefix = PREFIXES[tr["group"]] or "• "
        arrow = "↑" if tr["paired_dir"] == "inc" else "↓"
        y_labels_for.append(f"{prefix}{arrow} {tr['label'][:28]}")

        by_badge = "BY✓" if tr["a3_by_q"] <= 0.05 else "BH✓"
        ora_badge = f"ORA:{tr['inc_h']}↑/{tr['dec_h']}↓"
        if tr["in_t1"]:
            stat_str = f"{pct:.1f}% ({tr['up_st']}↑/{tr['dn_st']}↓){or_badge} · {by_badge}"
        else:
            stat_str = (
                f"{pct:.1f}% ({tr['up_st']}↑/{tr['dn_st']}↓) · A1 q={tr['a1_st_q']:.1e} · "
                f"{ora_badge} · {by_badge}"
            )
        ax_for.text(
            73.8,
            yp,
            stat_str,
            fontsize=6.85 if tr["in_t1"] else 7.0,
            fontweight="bold" if tr["in_t1"] else "normal",
            color=col if tr["all3_concordant"] else "#92400e",
            ha="left",
            va="center",
        )

    ax_for.set_yticks(y_pos)
    ax_for.set_yticklabels(y_labels_for, fontsize=7.5, fontweight="bold")
    ax_for.set_xticks([35, 45, 50, 55, 65])
    ax_for.set_xticklabels(["35%", "45%", "50% (Null)", "55%", "65%"], fontsize=8.0)
    ax_for.set_xlabel(
        f"Arm 3 Study Sign-Test % Increased (◆/● All {every.tested_studies:,} st. ±95% CI  |  "
        f"○ Oriented {oriented.tested_studies:,} st. for Table 1)",
        fontsize=8.4,
        fontweight="bold",
    )
    ax_for.set_title(
        f"B. The {len(selected_traits)} Study-FDR Multi-Arm Traits ({groups['table1_traits']} "
        f"Table 1 ★ + {groups['three_arm_additional']} Expanded 3-Arm ✓ + "
        f"{groups['ora_inverted']} ORA-Inv. ▲"
        + (f" + {groups['not_concordant']} Discordant ■" if groups["not_concordant"] else "")
        + f")\nOriented {oriented.tested_studies:,}-study filter (○): "
        f"Cellulose {cellulose.share_up:.1%}↑ / {cellulose.share_down:.1%}↓ (p = "
        f"{fmt_p(cellulose.p_value)}, BH₆ q = {fmt_p(cellulose.q_value)})",
        fontsize=8.35,
        fontweight="bold",
        loc="left",
        pad=10,
    )
    despine(ax_for, left=False)

    fig.text(
        0.05,
        0.025,
        f"{RERUN_SOURCE} · trait_profile_paired_contrast_trait_stats_pooled.tsv, "
        "trait_enrichment_traits_overall_ranked.tsv, trait_direction_screen.tsv, "
        "disease_control_comparison.tsv.",
        fontsize=7.7,
        color=MUTED,
    )

    out_path = out_dir / "tufte_plate7_three_arm_317_trait_concordance.png"
    fig.savefig(out_path, dpi=220)
    fig.savefig(out_path.with_suffix(".svg"))
    plt.close(fig)
    return out_path


def main(argv: Sequence[str] | None = None) -> int:
    # Whole option names only, so the recipe's guards see every option (#617).
    parser = argparse.ArgumentParser(
        description=__doc__.splitlines()[0] if __doc__ else None, allow_abbrev=False
    )
    parser.add_argument(
        "--input-root",
        type=Path,
        default=REPO_ROOT,
        help="Checkout holding the rerun tables (default: this one).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / PLATES_DIR,
        help=f"Where the plates go (default: this checkout's {PLATES_DIR}, whatever "
        "--input-root is).",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing plates. They are gitignored outputs with no other copy.",
    )
    parser.add_argument(
        "--code-commit",
        default="unrecorded",
        help="Revision that draws the plates, e.g. $(git rev-parse HEAD) with +dirty when the "
        "code is modified; the deck links the plate code at it (#613).",
    )
    args = parser.parse_args(argv)
    repo = args.input_root.resolve()
    out_dir = args.output_dir
    if any(out_dir.glob("*.png")) and not args.overwrite:
        raise SystemExit(f"{out_dir} has plates; back them up and pass --overwrite to replace them")
    check_classifiers(REPO_ROOT, repo)
    out_dir.mkdir(parents=True, exist_ok=True)
    # A run that fails partway must not leave new plates beside the previous record (#621).
    (out_dir / "plates_provenance.json").unlink(missing_ok=True)

    setup_tufte_rc()
    builders = (
        build_methods_plate1,
        build_methods_plate2,
        build_plate1,
        build_plate2,
        build_plate3,
        build_plate4,
        build_plate5,
        build_plate6,
        build_plate7,
    )
    inputs: dict[str, dict[str, str]] = {}
    outputs: dict[str, str] = {}
    PLATE_COUNTS.clear()
    for build in builders:
        INPUTS_READ.clear()
        path = build(repo, out_dir)
        inputs[path.name] = input_record(repo)
        for drawn in (path, path.with_suffix(".svg")):
            outputs[drawn.name] = file_sha256(drawn)
        print(f"Generated: {path} ({path.stat().st_size:,} bytes)")
    recorded = write_plates_provenance(out_dir, inputs, args.code_commit, outputs, PLATE_COUNTS)
    print(f"Recorded: {recorded}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
