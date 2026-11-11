"""Build the 30 NCASM Tufte plates: plates I-XX and the ten data-story plots (21-30).

Ported from ``build_tufte_plates_20.py`` in the slide-wizard deck folder
``slides/NCASM-disease-microbiome-traits/`` (sha256 ``95c07e5a5e13``, 2026-10-22; #2264).
The plate builders are called in the original order. The differences:

- the output folder is an argument (default ``data/derived/ncasm_tufte_plates/``), where
  the original wrote into the slide-wizard deck folder;
- every plate is drawn from matplotlib's default style, restored afterwards, so one
  plate's settings do not reach the next and the caller's settings are not changed;
- the data-story plots are written under ``data_story_plots/`` in the output folder and
  copied to their plate names, as before;
- one combined PDF, ``NCASM_Tufte_Plates_30.pdf``, replaces three identical 30-page
  copies (two of them named ``..._20.pdf``), and is written without a date;
- ``manifest.json`` gives paths relative to the output folder, the code revision, and
  for every plate where its numbers and data marks come from (``data_status``);
- the data-story plots draw their points from a random generator (#775), so the build
  refuses to draw them unless ``--audit-random-marks`` is given, and then each carries
  a banner saying so.
"""

from __future__ import annotations

import argparse
import json
import shutil
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import TypedDict

import matplotlib as mpl
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from mtp.ncasm_tufte_plates.data_status import (
    DataStatus,
    random_marks_status,
    typed_literals_status,
)
from mtp.ncasm_tufte_plates.plates_common import (
    PLATE_DPI,
    PLATE_HEIGHT_IN,
    PLATE_WIDTH_IN,
    ROMAN,
    TUFTE_PAPER,
    PlateInfo,
)
from mtp.ncasm_tufte_plates.plates_part1 import (
    build_plate_01,
    build_plate_02,
    build_plate_03,
    build_plate_04,
    build_plate_05,
    build_plate_06,
    build_plate_07,
    build_plate_08,
    build_plate_09,
    build_plate_10,
)
from mtp.ncasm_tufte_plates.plates_part2 import (
    build_plate_11,
    build_plate_12,
    build_plate_13,
    build_plate_14,
    build_plate_15,
    build_plate_16,
    build_plate_17,
    build_plate_18,
    build_plate_19,
    build_plate_20,
)
from mtp.ncasm_tufte_plates.story_plots import (
    PLOT_BUILDERS,
    PLOT_SPEAKER_NOTES,
    RANDOM_MARK_PLOTS,
)

__all__ = [
    "COMBINED_PDF",
    "DEFAULT_OUTPUT_DIR",
    "MANIFEST",
    "PLATE_BUILDERS",
    "SPEAKER_NOTES",
    "STORY_DIR",
    "Manifest",
    "ManifestEntry",
    "RandomMarksError",
    "build_plates",
    "main",
    "parse_args",
]

DEFAULT_OUTPUT_DIR = Path("data/derived/ncasm_tufte_plates")
COMBINED_PDF = "NCASM_Tufte_Plates_30.pdf"
MANIFEST = "manifest.json"
MANIFEST_SCHEMA_VERSION = 1
STORY_DIR = "data_story_plots"

SPEAKER_NOTES: dict[int, str] = {
    1: (
        "Plate I (Module 1 · Background — Hook & Paradigm): Why do microbiome meta-analyses "
        "across diseases look "
        "like a phonebook trap? At the taxon-name level, PD, HIV, and CRC share only 21.8% of "
        "depleted gut taxa (21.2\u201322.6%). "
        "When we map taxa to KG-Microbe ecophysiological guilds (METPO + GO), sharing jumps 2.3x "
        "to 50.6% (p=3.4e-6). "
        "Epistemological guardrail: we measure inherited genomic/physiological capacity of "
        "differentially abundant taxa, "
        "not single-sample metatranscriptomics."
    ),
    2: (
        "Plate II (Module 1 · Background — Rollup & Two-Part Curation Cascade): Reference traits "
        "in BacDive, Madin, and GO "
        "live mostly at strain/species level, whereas 16S studies report genera and families. "
        "Transitive NCBITaxon "
        "rdfs:subClassOf+ closure bridges this gap. Panel B shows the two-part curation cascade: "
        "Part A (KG-Microbe Trait "
        "Curation, Steps 1\u20133) filters 2,650 predicates to 317 non-MediaDive features (135 "
        "METPO + 102 CHEBI + 68 pathways + 12 GO) "
        "and 199 non-redundant traits; Part B (BugSigDB Study Curation, Steps 1\u20133) narrows "
        "5,142 signatures to 1,870 directional "
        "signatures across 935 Tier 1+2 subject-level case-vs-control contrasts in 548 studies."
    ),
    3: (
        "Plate III (Module 1 · Results — Rank Coverage & Pooled Shifts): Across 5,829 unique "
        "BugSigDB taxa (2,216 trait-bearing "
        "in Tier 1+2 contrasts), 45.9% of species (1,764/3,846) lack direct trait assertions, "
        "whereas Genus (76.0% mapped) and Family "
        "(73.5% mapped) provide dense trait coverage. Across 935 Tier 1+2 contrasts in 548 "
        "studies (Slide 4), facultatively "
        "aerobic (+0.16, 56.3% Up, study sign BY q=0.0408), halophilic (+0.13, 57.6% Up), coccus "
        "(+0.14, 55.4% Up), and arsenate "
        "reduction (+0.10, 54.4% Up) lead disease gains, while cellulose catabolism (-0.17, 60.7% "
        "Down, study sign BY q=0.0408) leads losses."
    ),
    4: (
        "Plate IV (Module 2 · Background — Three Statistical Arms): One shared taxon-trait matrix "
        "supports three "
        "distinct statistical questions: Arm 1 tests continuous within-contrast composition "
        "shifts (n=935 Tier 1+2 contrasts), "
        "Arm 2 tests finite-urn unique-taxon hypergeometric enrichment (N=2,216 unique Tier 1+2 "
        "taxa), and Arm 3 tests "
        "study-level majority votes and paired log2 ratios (n=548 studies). Because sample sizes "
        "span orders of "
        "magnitude over the same data, cross-arm p-values are never pooled."
    ),
    5: (
        "Plate V (Module 2 · Background — Direction & Evidence Contract): 'Direction' has four "
        "distinct mathematical "
        "definitions: marginal prevalence, continuous paired shift Delta, unweighted study "
        "majority vote %, and scale-invariant "
        "log-ratios. On the Condition x Body-Site grid (1,872 cells), explicit display-state "
        "partitioning separates the "
        "316 drawn cells (16.9%) from 209 low-study masked cells (11.2%), 7 "
        "sampled-without-signal cells (0.4%), and "
        "1,340 unsampled cells (71.6%)."
    ),
    6: (
        "Plate VI (Module 3 · Background — Four Threats & Near-Twin Redundancy): Knowledge-graph "
        "trait profiling faces "
        "four statistical threats: (1) curation-density bias toward model pathogens, (2) "
        "correlated-trait FDR runaway from "
        "near-twin terms (up to 95.9% Jaccard overlap), (3) unit-sum compositional closure, and "
        "(4) multi-contrast "
        "pseudoreplication. Each has a dedicated structural control."
    ),
    7: (
        "Plate VII (Module 3 · Results — Density Neutralization & Taxon-Urn ORA): Within-contrast "
        "Inc-Dec pairing collapses "
        "the annotation-density regression slope 20-fold (from -0.237 to -0.011, R^2=0.0016), and "
        "PC1 density residualization "
        "preserves effect sizes (r=0.994). In Arm 2, deduplicating the hypergeometric urn from "
        "26,412 assertions to "
        "1,246 unique background taxa (2,216 unique Tier 1+2 taxa) restores a uniform null "
        "(median p=0.543)."
    ),
    8: (
        "Plate VIII (Module 3 · Results — Permutation Nulls & Log-Ratios): Flipping "
        "Increased/Decreased labels within "
        "studies 100 times across 246 Condition x Feces cells (24,600 tests) produces zero false "
        "discoveries at BH q <= 0.10 "
        "(0.00% FPR vs 3 observed hits at q<=0.10 plus CRC at p=0.015). Scale-invariant paired "
        "log2 trait ratios across 548 "
        "studies show 2 of 15 pairs pass BY q<=0.05 (Facultative/Cellulose +1.05, BY q=0.0029; "
        "Arsenate/Cellulose +0.91, BY q=0.0029)."
    ),
    9: (
        "Plate IX (Module 4 · Background — Redox Hypothesis & 317-Trait Screen): Under the Winter "
        "& Byndloss epithelial "
        "oxygenation hypothesis, healthy colonocytes consume O2 via butyrate beta-oxidation to "
        "maintain luminal anoxia, "
        "whereas inflamed mucosa leaks O2 and nitrate into the lumen. We test this unbiasedly "
        "across all 317 non-MediaDive "
        "traits in 935 Tier 1+2 contrasts (548 studies), yielding 22 study t-test BH hits "
        "(q<=0.05), 6 study t-test BY hits "
        "(q<=0.05), and 2 study sign BY hits at m=317 (q=0.0408; 3 at m_eff=75)."
    ),
    10: (
        "Plate X (Module 4 · Results — Global 317-Trait Screen & Study Votes): In the unbiased "
        "317-trait Wilcoxon screen across "
        "935 Tier 1+2 contrasts, cellulose catabolism is the #1 depleted trait (-0.17, 60.7% "
        "Down, q=3.9e-7; study sign 181:265, "
        "BY q=0.0408), while facultative anaerobe (+0.16, 56.3% Up, study sign 295:219, BY "
        "q=0.0408) leads study sign gains. "
        "Arm 3 one-vote-per-study sign tests across 548 studies confirm the 6 Table 1 focal traits."
    ),
    11: (
        "Plate XI (Module 4 · Results — Arm 1 & Arm 2 Disease Atlases across Seven Major "
        "Conditions): Across 246 Condition x Feces "
        "cells, three disease-trait pairs pass multi-cell BH FDR (q<=0.10) and CRC x Cellulose is "
        "nominally significant: "
        "PD x Cellulose (n=71 stool, 75.4% Dec, p=8.4e-5), HIV x Cellulose (n=20 stool, 70.6% "
        "Dec, p=0.072), PD x Moderately halophilic "
        "(+24.4 pp, q=0.088), and CRC x Cellulose (n=44 stool, 69.7% Dec, p=0.035). All four are "
        "independently validated by Arm 2 ORA."
    ),
    12: (
        "Plate XII (Module 4 · Results — Paired McNemar Oxygen-Guild Rewiring): Testing whether "
        "aerotolerant gains and "
        "anaerobe losses co-occur within the exact same contrast, the paired McNemar test across "
        "935 Tier 1+2 contrasts (897 active) "
        "finds 343 canonical redox-swap pairs vs 231 reverse pairs (OR=1.48, p=3.4e-6), and "
        "strict Stool (n=547) shows 261 vs 134 (OR=1.95, p=1.9e-10)."
    ),
    13: (
        "Plate XIII (Module 5 · Background — Venn Partitions & Six Nulls): Why does an exclusive "
        "2-way Venn slice like "
        "'HIV & CRC only = 1' look small when full HIV & CRC overlap is 15/35 (42.9%)? Because 14 "
        "of the 15 shared taxa "
        "are ALSO depleted in Parkinson's disease (the 3-way core of 14/41 = 34.1%). Benchmarking "
        "against U1-U6 null "
        "models shows observed convergence (50.6%) is 2.3x the strict fecal depleted-taxon urn U4 "
        "(21.8% across 402 taxa, p=3.4e-6; U6 whole-study p=0.742)."
    ),
    14: (
        "Plate XIV (Module 5 · Results — 142-Trait Pathway Benchmark & Four Pairs): Across 142 "
        "qualifying non-MediaDive "
        "traits in PD, CRC, and HIV stool contrasts (n=135), cellulose catabolism ranks #1 among "
        "GO biological processes "
        "(#2 overall at 50.6% taxon / 52.1% guild after acetogenesis 66.8%). Panel B shows the "
        "four headline Condition x Feces shifts."
    ),
    15: (
        "Plate XV (Module 5 · Results — Three-Disease Slopegraph & 41-Taxon Partition): Pairwise "
        "Jaccard overlap across Tier 1+2 "
        "PD (n=71), CRC (n=44), and HIV (n=20) stool contrasts rises from 21.8% across all 402 "
        "depleted stool taxa (21.2-22.6%) to "
        "50.6% (42.9-57.5%, p=3.4e-6) among the 41 fiber-digesting taxa. Over two-thirds (28/41 = "
        "68.3%) are shared in >=2 diseases, "
        "and 14/41 (34.1%) in all 3."
    ),
    16: (
        "Plate XVI (Module 5 · Results — Nine Core Genera & Rank Sensitivity): Nine butyrate- and "
        "fiber-degrading genera "
        "(Roseburia [34 contrasts], Clostridium [19], Bacteroides [17], Blautia [16], Lachnospira "
        "[12], Ruminococcus [12], "
        "Anaerostipes [6], Eubacterium [6], Butyrivibrio [3]) plus five named cellulolytic "
        "species anchor the guild. Restricting "
        "to Genus + Species only (21 taxa) preserves 42.9% overlap (2.1x U5 null of 20.0%, "
        "p=7.0e-4)."
    ),
    17: (
        "Plate XVII (Module 6 · Background — Contributor Math & LOTO Architecture): Unconditional "
        "full-cell contributor "
        "decomposition divides each taxon's carrier frequency by all N contrasts in the cell so "
        "shares sum additively. "
        "Leave-One-Taxon-Out (LOTO, Issue #190) across 2,216 unique Tier 1+2 taxa produces 0 "
        "direction flips, and clade-level "
        "ablation separates multi-family guilds from single-clade hitchhikers."
    ),
    18: (
        "Plate XVIII (Module 6 · Results — Slide 21 PD Species & Family Contributors): "
        "Decomposing Tier 1+2 Parkinson's disease "
        "shifts across the 11 taxa with |Delta| >= 5 pp (Slide 21: 3 species and 8 families; 80 "
        "PD vs 855 other contrasts) shows "
        "Lachnospiraceae (36 PD / 160 other, +12.4 pp PD specificity) driving cellulose loss "
        "(-28.6 pp net), while Enterobacteriaceae "
        "(+1.9 pp, +14.3 pp FA net) and Streptococcaceae (+1.4 pp, +4.8 pp coccus net) drive gains."
    ),
    19: (
        "Plate XIX (Module 6 · Results — LOTO Sensitivity & Assay/Literature Audit): Single-taxon "
        "LOTO across 2,216 Tier 1+2 taxa "
        "yields 0 direction flips; under clade ablation, cellulose loss (p=0.0256), facultative "
        "anaerobiosis (p=0.0089), and "
        "halophily (p=0.0180) survive as multi-family guilds, whereas coccus and arsenate "
        "hitchhike. Both 16S rRNA (n=805, 62.0% down) "
        "and Shotgun WMS (n=128, 69.0% down) confirm cellulose loss."
    ),
    20: (
        "Plate XX (Module 7 · Results & Synthesis — Anatomical Specificity & Luminal Dysbiosis "
        "Model): Shared fiber-guild "
        "loss is strictly a luminal stool phenomenon (73.0% decreased in PD+CRC+HIV stool, "
        "81/111, p=1.4e-6; 63.7% across all 496 "
        "Tier 1+2 stool contrasts, p=1.1e-9; 50.6% Jaccard overlap), opposite oral cavity (34.6% "
        "decreased, 18/52; Stool vs Oral "
        "OR=5.10, p=5.3e-6) and absent in CRC colon mucosa (37.5%, 3/8, p=0.73)."
    ),
}

PLATE_BUILDERS: tuple[Callable[[Path], PlateInfo], ...] = (
    build_plate_01,
    build_plate_02,
    build_plate_03,
    build_plate_04,
    build_plate_05,
    build_plate_06,
    build_plate_07,
    build_plate_08,
    build_plate_09,
    build_plate_10,
    build_plate_11,
    build_plate_12,
    build_plate_13,
    build_plate_14,
    build_plate_15,
    build_plate_16,
    build_plate_17,
    build_plate_18,
    build_plate_19,
    build_plate_20,
)


class ManifestEntry(TypedDict):
    """One plate in ``manifest.json``; paths are relative to the output folder."""

    plate_num: int
    roman: str
    stem: str
    png_path: str
    svg_path: str
    pdf_path: str | None
    out_of_bounds: list[str]
    overlaps: list[str]
    speaker_notes: str
    data_status: DataStatus


class Manifest(TypedDict):
    """``manifest.json``: the build's record of what it drew and from what."""

    schema_version: int
    generator: str
    code_commit: str | None
    combined_pdf: str
    plates: list[ManifestEntry]


class RandomMarksError(ValueError):
    """Plates that draw random points were requested without ``audit_random_marks``."""


def _report(entry: ManifestEntry) -> None:
    oob, ovl = entry["out_of_bounds"], entry["overlaps"]
    status = "OK" if (not oob and not ovl) else f"WARN(oob={len(oob)}, ovl={len(ovl)})"
    print(f"Plate {entry['plate_num']:02d} ({entry['stem']}): {status}")
    if oob:
        print(f"   Out-of-bounds: {oob}")
    if ovl:
        print(f"   Overlaps: {ovl}")


def _random_mark_plates() -> list[int]:
    first_story_plate = len(PLATE_BUILDERS) + 1
    return [
        number
        for number, plot_builder in enumerate(PLOT_BUILDERS, start=first_story_plate)
        if plot_builder in RANDOM_MARK_PLOTS
    ]


def build_plates(
    out_dir: Path, *, audit_random_marks: bool = False, code_commit: str | None = None
) -> Manifest:
    """Draw all 30 plates into ``out_dir`` and write the manifest and the combined PDF.

    Raises ``RandomMarksError`` before drawing anything when a plate's points come from a
    random generator and ``audit_random_marks`` is false (#775)."""
    random_plates = _random_mark_plates()
    if random_plates and not audit_random_marks:
        numbers = ", ".join(str(number) for number in random_plates)
        raise RandomMarksError(
            f"plates {numbers} draw their points from a random generator, not from data "
            "(#775); pass --audit-random-marks to draw them with a banner saying so"
        )
    out_dir.mkdir(parents=True, exist_ok=True)
    plates: list[ManifestEntry] = []

    for builder in PLATE_BUILDERS:
        with mpl.style.context("default"):
            info = builder(out_dir)
        pnum = info["plate_num"]
        entry: ManifestEntry = {
            "plate_num": pnum,
            "roman": ROMAN[pnum - 1],
            "stem": info["stem"],
            "png_path": Path(info["png_path"]).name,
            "svg_path": Path(info["svg_path"]).name,
            "pdf_path": Path(info["pdf_path"]).name,
            "out_of_bounds": info["out_of_bounds"],
            "overlaps": info["overlaps"],
            "speaker_notes": SPEAKER_NOTES[pnum],
            "data_status": typed_literals_status(),
        }
        plates.append(entry)
        _report(entry)

    story_dir = out_dir / STORY_DIR
    for idx, plot_builder in enumerate(PLOT_BUILDERS, start=len(PLATE_BUILDERS) + 1):
        if plot_builder not in RANDOM_MARK_PLOTS:
            raise NotImplementedError(
                f"plate {idx} is not listed as drawing random points, and no builder that "
                "reads run tables records its inputs yet (#2264)"
            )
        with mpl.style.context("default"):
            story = plot_builder(story_dir)
        plate_stem = f"plate_{idx:02d}_{story['stem']}"
        shutil.copy2(story["png_path"], out_dir / f"{plate_stem}.png")
        shutil.copy2(story["svg_path"], out_dir / f"{plate_stem}.svg")
        entry = {
            "plate_num": idx,
            "roman": ROMAN[idx - 1],
            "stem": plate_stem,
            "png_path": f"{plate_stem}.png",
            "svg_path": f"{plate_stem}.svg",
            "pdf_path": None,
            "out_of_bounds": story["out_of_bounds"],
            "overlaps": [f"{left!r} <-> {right!r}" for left, right in story["overlaps"]],
            "speaker_notes": PLOT_SPEAKER_NOTES[idx - len(PLATE_BUILDERS)],
            "data_status": random_marks_status(),
        }
        plates.append(entry)
        _report(entry)

    manifest: Manifest = {
        "schema_version": MANIFEST_SCHEMA_VERSION,
        "generator": "python -m mtp.ncasm_tufte_plates",
        "code_commit": code_commit,
        "combined_pdf": COMBINED_PDF,
        "plates": plates,
    }
    (out_dir / MANIFEST).write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    # One page per plate, each the plate's PNG on a full-bleed canvas.
    combined_pdf = out_dir / COMBINED_PDF
    with (
        mpl.style.context("default"),
        PdfPages(combined_pdf, metadata={"CreationDate": None}) as pdf,
    ):
        for item in plates:
            img = mpimg.imread(out_dir / item["png_path"])
            fig = plt.figure(
                figsize=(PLATE_WIDTH_IN, PLATE_HEIGHT_IN), dpi=PLATE_DPI, facecolor=TUFTE_PAPER
            )
            ax = fig.add_axes((0.0, 0.0, 1.0, 1.0))
            ax.imshow(img)
            ax.axis("off")
            pdf.savefig(fig, facecolor=TUFTE_PAPER)
            plt.close(fig)

    total_oob = sum(len(item["out_of_bounds"]) for item in plates)
    total_overlaps = sum(len(item["overlaps"]) for item in plates)
    print(f"\nCompleted {len(plates)} Tufte plates -> {out_dir}")
    print(f"Combined PDF -> {combined_pdf}")
    print(f"Geometry summary: total_oob={total_oob}, total_overlaps={total_overlaps}")
    if random_plates:
        numbers = ", ".join(str(number) for number in random_plates)
        print(f"Plates drawn from random points, with the banner (#775): {numbers}")
    return manifest


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Command-line arguments of ``python -m mtp.ncasm_tufte_plates``."""
    parser = argparse.ArgumentParser(
        prog="python -m mtp.ncasm_tufte_plates",
        description="Draw the 30 NCASM Tufte plates, their manifest and the combined PDF.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"folder to write into (default: {DEFAULT_OUTPUT_DIR})",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="write into an output folder that already holds files",
    )
    parser.add_argument(
        "--audit-random-marks",
        action="store_true",
        help=(
            "also draw the plates whose points come from a random generator, each with a "
            "banner saying so (#775); without it the build stops before drawing anything"
        ),
    )
    parser.add_argument(
        "--code-commit",
        help="revision of the code drawing the plates, recorded in the manifest",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    """Build the plates; refuse a non-empty output folder unless ``--overwrite``."""
    args = parse_args(argv)
    out_dir: Path = args.output_dir
    if out_dir.is_dir() and any(out_dir.iterdir()) and not args.overwrite:
        print(f"{out_dir} is not empty; pass --overwrite to replace its plates")
        return 1
    try:
        build_plates(
            out_dir, audit_random_marks=args.audit_random_marks, code_commit=args.code_commit
        )
    except RandomMarksError as error:
        print(error)
        return 1
    return 0
