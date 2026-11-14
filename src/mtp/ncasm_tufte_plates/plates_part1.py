"""Part 1 of the 20-slide NCASM Tufte Plate generator (Plates I through X)."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from mtp.ncasm_tufte_plates.plates_common import (
    TUFTE_AMBER,
    TUFTE_CRIMSON,
    TUFTE_FOREST,
    TUFTE_INK,
    TUFTE_LIGHT,
    TUFTE_MUTED,
    TUFTE_PAPER,
    TUFTE_SAND,
    TUFTE_SLATE,
    PlateInfo,
    SlopeSeries,
    apply_range_frame,
    draw_labeled_hline,
    draw_marginal_note,
    draw_slopegraph_series,
    draw_tufte_table,
    new_plate,
    save_and_verify_plate,
    setup_canvas_panel,
    style_tufte_axis,
    title_panel,
    two_panels,
)


def build_plate_01(out_dir: Path) -> PlateInfo:
    """Plate I: Hook & Paradigm Shift — From Taxonomic Phonebooks to Ecophysiological Guilds."""
    fig = new_plate(
        plate_num=1,
        module_tag="Module 1 — Background",
        role_tag="Hook & Paradigm Shift",
        title="Five Names, One Missing Job: Why Taxonomic Lists Miss Shared Dysbiosis",
        subtitle=(
            "Cross-disease microbiome meta-analyses appear discordant at the taxon-list level "
            "(20.2\u201325.4% overlap) "
            "yet show 2.3\u00d7 overlap (50.6%, U\u2083 p = 1.1\u00d710\u207b\u00b9\u2070; "
            "U\u2086 p = 0.144) among GO:0030245-annotated taxa."
        ),
        footer=(
            "Source: BugSigDB v1.3.1 (548 Tier 1+2 studies; 935 subject-control contrasts) "
            "integrated with KG-Microbe (1.4M nodes, 5.2M edges; 2,216 mapped Tier 1+2 taxa, 317 "
            "traits)."
        ),
    )
    ax_a, ax_b = two_panels(fig)

    # Panel A: 2-column slopegraph of Taxon-List Overlap vs Functional-Guild Overlap
    style_tufte_axis(ax_a, keep_bottom=True, keep_left=True)
    x_cols = [0.0, 1.0]
    series: list[SlopeSeries] = [
        {
            "name": "PD \u2229 CRC",
            "values": [25.4, 57.5],
            "color": TUFTE_CRIMSON,
            "lw": 2.4,
            "marker_size": 6.8,
            "zorder": 5,
            "left_label": "PD \u2229 CRC: 25.4%",
            "right_label": "57.5% (23/40 cellulose taxa)",
            "label_weight": "bold",
        },
        {
            "name": "PD \u2229 HIV",
            "values": [20.2, 51.4],
            "color": TUFTE_SLATE,
            "lw": 2.1,
            "marker_size": 6.2,
            "zorder": 4,
            "left_label": "PD \u2229 HIV: 20.2%",
            "right_label": "51.4% (18/35 cellulose taxa)",
            "label_weight": "bold",
        },
        {
            "name": "Mean Pairwise",
            "values": [21.8, 50.6],
            "color": TUFTE_INK,
            "lw": 1.8,
            "ls": "--",
            "marker_size": 5.5,
            "zorder": 4,
            "left_label": "Mean across pairs: 21.8%",
            "right_label": "50.6% (2.3\u00d7 U\u2084 · U\u2086 p = 0.144)",
            "label_weight": "bold",
        },
        {
            "name": "HIV \u2229 CRC",
            "values": [20.3, 42.9],
            "color": TUFTE_FOREST,
            "lw": 2.0,
            "marker_size": 6.0,
            "zorder": 4,
            "left_label": "HIV \u2229 CRC: 20.3%",
            "right_label": "42.9% (15/35 cellulose taxa)",
            "label_weight": "normal",
        },
    ]
    draw_slopegraph_series(
        ax_a,
        x_cols,
        series,
        min_dy=3.8,
        left_pad=0.05,
        annotation_size=8.4,
    )
    draw_labeled_hline(
        ax_a,
        20.0,
        "Genus+Species baseline (U5 = 20.0%)",
        x_text=0.48,
        dy=-3.6,
        color=TUFTE_MUTED,
        size=8.6,
    )
    ax_a.set_xlim(-0.26, 2.05)
    ax_a.set_ylim(10.0, 62.0)
    ax_a.set_xticks(x_cols)
    ax_a.set_xticklabels(
        [
            "All Depleted Stool Taxa\n(U4 Baseline, N = 402)",
            "Cellulose-Annotated Taxa\n(N = 41 Stool Taxa)",
        ],
        fontsize=10.0,
        fontweight="bold",
    )
    ax_a.set_yticks([15, 25, 35, 45, 57.5])
    ax_a.set_ylabel("Pairwise Cross-Disease Jaccard Overlap (%)", fontsize=10.2, color=TUFTE_INK)
    apply_range_frame(ax_a, x_bounds=(0.0, 1.0), y_bounds=(15.0, 57.5))
    title_panel(
        ax_a,
        "A. Stool Overlap: All Taxa (21.8%) vs. Cellulose Guild (50.6%)",
        "Pairwise Jaccard overlap across Parkinson's (PD), HIV & Colorectal Cancer (CRC)",
    )

    # Panel B: 5-row Foundation Table + Epistemological Guardrail
    setup_canvas_panel(ax_b)
    title_panel(
        ax_b,
        "B. Knowledge-Graph Architecture (548 Studies · 935 Contrasts)",
        "Five building blocks linking 2,216 mapped Tier 1+2 taxa to 317 curated traits",
    )
    headers = [
        (0.01, "Component", "left"),
        (0.30, "Scale / Standard", "left"),
        (0.61, "Role in Trait Profiling", "left"),
    ]
    rows = [
        (
            "1. BugSigDB Corpus",
            "548 Tier 1+2 studies\n935 paired contrasts",
            "Subject-level case-vs-control differential-\nabundance contrasts (897 active pairs)",
        ),
        (
            "2. KG-Microbe",
            "2,216 mapped taxa\n(3,400 Tier 1+2 TaxIDs)",
            "Integrates BacDive, MediaDive, Madin\net al. (2020), GOLD, and UniProt GO",
        ),
        (
            "3. METPO + CHEBI + GO",
            "317 curated features\n(135M + 102C + 68P + 12G)",
            "Standardized oxygen, morphology, salt,\nand carbohydrate catabolism terms",
        ),
        (
            "4. Ancestral Closure",
            "NCBITaxon\nrdfs:subClassOf+",
            "Propagates strain/species traits to\n16S genus (94.8% occ.) & family lists",
        ),
        (
            "5. Redox Hypothesis",
            "Winter & Byndloss\n(2013, 2017)",
            "Inflamed epithelium leaks O\u2082/NO\u2083\u207b,\nreplacing anaerobes with "
            "aerotolerants",
        ),
    ]
    y_bot = draw_tufte_table(
        ax_b,
        headers,
        rows,
        y_top=0.91,
        row_step=0.112,
        font_size=8.9,
        header_size=9.0,
        row_accents={0: TUFTE_SLATE, 4: TUFTE_CRIMSON},
    )
    draw_marginal_note(
        ax_b,
        x=0.01,
        y_top=y_bot - 0.025,
        kicker="Epistemological Guardrail: Reference Capacity \u2260 In Situ Transcription",
        body_lines=[
            "\u2022 Trait profiles measure inherited genomic/physiological capacity of reported "
            "taxa.",
            "\u2022 They test ecological guild selection across hosts, not single-sample "
            "metatranscriptomics.",
        ],
        accent_color=TUFTE_CRIMSON,
        line_step=0.045,
        font_size=8.9,
    )

    return save_and_verify_plate(fig, 1, "hook_paradigm_shift", out_dir)


def build_plate_02(out_dir: Path) -> PlateInfo:
    """Plate II: Ancestral NCBITaxon Closure & Two-Part Curation Cascade."""
    fig = new_plate(
        plate_num=2,
        module_tag="Module 1 — Background",
        role_tag="Rollup & Filter Cascade",
        title="Filtering 2,650 Raw Terms to 317 Traits Across 935 Tier 1+2 Contrasts",
        subtitle=(
            "Two independent curation cascades: (A) filtering 2,650 KG-Microbe terms to 317 "
            "non-MediaDive features, and "
            "(B) curating 1,237 BugSigDB studies (959 human) to 548 Tier 1+2 case-control studies "
            "(935 contrasts)."
        ),
        footer=(
            "Source: KG-Microbe predicate audit (2,333 MediaDive growth-medium terms excluded -> "
            "317 features) and "
            "BugSigDB Tier 1+2 curation (1,237 -> 959 -> 548 studies; 2,571 -> 1,825 -> 935 "
            "contrasts)."
        ),
    )
    ax_a, ax_b = two_panels(fig)

    # Panel A: Taxonomic Rollup & Trait Column Filtering Table
    setup_canvas_panel(ax_a)
    title_panel(
        ax_a,
        "A. Bridging the Taxonomic Rank Gap via Ancestral Closure",
        "42.2% of unique species lack direct records; genus rollup covers 76.3% unique (94.8% "
        "occ.)",
    )
    headers = [
        (0.01, "Pipeline Stage", "left"),
        (0.32, "Input Representation", "left"),
        (0.64, "Resolved Output", "left"),
    ]
    rows = [
        (
            "1. Source Assertions",
            "BacDive / Madin / GO\nat Strain & Species",
            "Direct records on cultured strains;\n852 / 2,018 species (42.2%) unmapped",
        ),
        (
            "2. Ancestral Rollup",
            "NCBITaxon hierarchy\n(rdfs:subClassOf+)",
            "Genus: 76.3% unique (94.8% occ.,\n8,139/8,582) · Family: 95.3% occ.",
        ),
        (
            "3. Trait Curation\n(Part A: Steps 1\u20133)",
            "2,650 raw KG-Microbe\ntrait columns",
            "Drop 2,333 growth media (-88.0%)\n\u2192 317 features (135M+102C+68P+12G)",
        ),
        (
            "4. Cohort Curation\n(Part B: Steps 1\u20133)",
            "1,237 export \u2192 959\nhuman-host studies",
            "2,571 \u2192 1,825 \u2192 935 Tier 1+2 pairs\nacross 548 case-control studies",
        ),
    ]
    y_bot = draw_tufte_table(
        ax_a,
        headers,
        rows,
        y_top=0.90,
        row_step=0.145,
        font_size=9.0,
        header_size=9.1,
        row_accents={1: TUFTE_SLATE, 2: TUFTE_CRIMSON, 3: TUFTE_FOREST},
    )
    draw_marginal_note(
        ax_a,
        x=0.01,
        y_top=y_bot - 0.035,
        kicker="Why Separate Trait Curation (A) from Study Curation (B)?",
        body_lines=[
            "\u2022 Part A removes 2,333 DSMZ growth-medium tautologies to retain 317 biological "
            "features.",
            "\u2022 Part B restricts BugSigDB to 548 Tier 1+2 subject-level case-vs-control "
            "studies (935 pairs).",
        ],
        accent_color=TUFTE_FOREST,
        line_step=0.048,
        font_size=9.0,
    )

    # Panel B: Two-Part Curation Cascade with Separate 1-3 Numbering & Demarcation at y=2.75
    style_tufte_axis(ax_b, keep_bottom=True, keep_left=False)
    ax_b.set_xlim(-2050, 5600)
    ax_b.set_ylim(-1.35, 6.85)
    ax_b.set_yticks([])

    sec_x = -2000.0
    y_label_x = -90.0

    # Section A header: Trait Curation
    ax_b.text(
        sec_x,
        6.35,
        "A. TRAIT CURATION (KG-MICROBE \u00b7 STEPS 1\u20133)",
        fontsize=8.6,
        fontweight="bold",
        color=TUFTE_CRIMSON,
        ha="left",
        va="center",
    )

    # Dashed horizontal demarcation rule at y = 2.75 between Part A and Part B
    ax_b.plot(
        [sec_x, 5300],
        [2.75, 2.75],
        color=TUFTE_INK,
        lw=1.1,
        ls="--",
        zorder=1,
    )

    # Section B header: Signature & Study Curation
    ax_b.text(
        sec_x,
        2.15,
        "B. SIGNATURE & STUDY CURATION (BUGSIGDB \u00b7 STEPS 1\u20133)",
        fontsize=8.6,
        fontweight="bold",
        color=TUFTE_SLATE,
        ha="left",
        va="center",
    )

    items = [
        (5.5, "1. Raw KG-Microbe Terms", 2650, "2,650 total terms", TUFTE_LIGHT),
        (
            4.5,
            "2. Exclude Growth Media (-2,333)",
            2333,
            "2,333 excluded (88.0% media)",
            TUFTE_AMBER,
        ),
        (
            3.5,
            "3. Retained Biological Features",
            317,
            "317 features (135M+102C+68P+12G)",
            TUFTE_CRIMSON,
        ),
        (1.3, "1. Directional Signatures", 1870, "1,870 sigs (2,571->1,825->935)", TUFTE_SLATE),
        (0.3, "2. Tier 1+2 Paired Contrasts", 935, "935 contrasts (897 active pairs)", TUFTE_SLATE),
        (-0.7, "3. Curated Human Studies", 548, "548 studies (1,237->959->548)", TUFTE_FOREST),
    ]

    for y, lab, val, vtxt, col in items:
        ax_b.barh(y, val, height=0.52, color=col, alpha=0.90, edgecolor="none", zorder=2)
        ax_b.text(
            y_label_x,
            y,
            lab,
            fontsize=8.2,
            fontweight="bold" if col in (TUFTE_CRIMSON, TUFTE_SLATE, TUFTE_FOREST) else "normal",
            color=TUFTE_INK,
            ha="right",
            va="center",
        )
        ax_b.text(
            val + 90.0,
            y,
            vtxt,
            fontsize=8.2,
            fontweight="bold",
            color=col if col != TUFTE_LIGHT else TUFTE_INK,
            ha="left",
            va="center",
        )

    ax_b.set_xticks([0, 1000, 2000, 3000])
    ax_b.set_xticklabels(["0", "1,000", "2,000", "3,000"])
    ax_b.set_xlabel(
        "Count of Traits (Part A) or Signatures / Contrasts / Studies (Part B)",
        fontsize=9.4,
        color=TUFTE_INK,
    )
    apply_range_frame(ax_b, x_bounds=(0, 3000), y_bounds=None, keep_labels=True)
    title_panel(
        ax_b,
        "B. Two-Part Curation: (A) Traits & (B) Case-Control Contrasts",
        "Separate 1\u20133 numbering for KG-Microbe trait filtering and BugSigDB Tier 1+2 cohort "
        "curation",
    )

    return save_and_verify_plate(fig, 2, "rollup_and_filter_cascade", out_dir)


def build_plate_03(out_dir: Path) -> PlateInfo:
    """Plate III: Rank-Stratified Trait Coverage & Six Focal Pooled Shifts."""
    fig = new_plate(
        plate_num=3,
        module_tag="Module 1 — Results",
        role_tag="Coverage & Pooled Shifts",
        title="Genus & Family Annotations Cover ~95% of 16S Occurrences & Reveal Six Shifts",
        subtitle=(
            "Across 3,400 unique Tier 1+2 TaxIDs (19,388 occurrences), only 57.8% of unique "
            "species have direct traits "
            "(42.2% unmapped), vs. 76.3% of unique genera (94.8% occ.) and 75.3% of families "
            "(95.3% occ.)."
        ),
        footer=(
            "Source: Taxon-trait coverage across 3,400 Tier 1+2 BugSigDB TaxIDs (2,216 "
            "trait-bearing; 19,388 occurrences) "
            "and pooled within-experiment trait contrasts across 935 Tier 1+2 contrasts in 548 "
            "studies."
        ),
    )
    ax_a, ax_b = two_panels(fig)

    # Panel A: Stacked Horizontal Coverage Bar by Rank (3,400 Unique TaxIDs · 19,388 Occurrences)
    style_tufte_axis(ax_a, keep_bottom=True, keep_left=False)
    ranks = [
        "Species (n = 2,018)",
        "Family (n = 247)",
        "Genus (n = 973)",
        "Order (n = 112)",
        "Phylum (n = 49)",
        "Class (n = 53)",
    ]
    mapped_uniq = np.array([57.8, 75.3, 76.3, 81.2, 93.9, 94.3])
    unmapped = 100.0 - mapped_uniq
    y_pos = np.arange(len(ranks))

    ax_a.barh(y_pos, mapped_uniq, height=0.55, color=TUFTE_SLATE, label="Mapped unique TaxIDs (%)")
    ax_a.barh(
        y_pos,
        unmapped,
        left=mapped_uniq,
        height=0.55,
        color=TUFTE_SAND,
        edgecolor=TUFTE_LIGHT,
        lw=0.6,
        label="Unmapped unique (%)",
    )

    right_notes = [
        "57.8% uniq · 75.8% occ. (852/2,018 = 42.2% gap)",
        "75.3% uniq · 95.3% occ. coverage",
        "76.3% uniq · 94.8% occ. (8,139 / 8,582)",
        "81.2% uniq · 97.6% occ. coverage",
        "93.9% uniq · 99.6% occ. coverage",
        "94.3% uniq · 99.5% occ. coverage",
    ]
    for i, (m_u, r_note) in enumerate(zip(mapped_uniq, right_notes, strict=False)):
        ax_a.text(
            m_u / 2.0,
            i,
            f"{m_u:.1f}%",
            va="center",
            ha="center",
            fontsize=8.6,
            fontweight="bold",
            color=TUFTE_PAPER,
        )
        ax_a.text(
            101.5,
            i,
            r_note,
            va="center",
            ha="left",
            fontsize=8.4,
            fontweight="bold" if i in (0, 2) else "normal",
            color=TUFTE_CRIMSON if i == 0 else TUFTE_INK,
        )

    ax_a.set_yticks(y_pos)
    ax_a.set_yticklabels(ranks, fontsize=9.4, color=TUFTE_INK)
    ax_a.set_xlim(0, 168)
    ax_a.set_ylim(-0.6, 5.7)
    ax_a.set_xticks([0, 25, 50, 75, 100])
    ax_a.set_xlabel("Percentage of 3,400 Tier 1+2 Unique TaxIDs (%)", fontsize=9.8, color=TUFTE_INK)
    ax_a.legend(
        loc="upper left",
        bbox_to_anchor=(0.0, 1.01),
        ncol=2,
        frameon=False,
        fontsize=8.6,
        handlelength=1.2,
    )
    apply_range_frame(ax_a, x_bounds=(0, 100), y_bounds=None)
    title_panel(
        ax_a,
        "A. Rank Coverage: 3,400 Unique TaxIDs · 19,388 Occurrences",
        "852 / 2,018 species (42.2%) lack direct traits; Genus covers 94.8% of occurrences",
    )

    # Panel B: Diverging Lollipop + 95% CI of Six Focal Pooled Contrasts
    style_tufte_axis(ax_b, keep_bottom=True, keep_left=False)
    traits = [
        "Cellulose catabolism\n(n = 678 · q = 1.9e-7 · Depleted)",
        "Moderately halophilic\n(n = 824 · q = 6.9e-4)",
        "Sphere-shaped\n(n = 768 · q = 0.0040)",
        "Coccus-shaped\n(n = 780 · q = 0.025)",
        "Facultatively aerobic\n(n = 853 · q = 0.011)",
        "Arsenate reduction\n(n = 773 · q = 0.044)",
    ]
    means = np.array([-0.25, 0.13, 0.14, 0.17, 0.21, 0.24])
    ci_lo = np.array([-0.33, 0.05, 0.06, 0.08, 0.12, 0.15])
    ci_hi = np.array([-0.16, 0.22, 0.23, 0.26, 0.30, 0.33])
    y_b = np.arange(len(traits))

    ax_b.axvline(0.0, color=TUFTE_INK, lw=0.9, ls="--", zorder=1)
    for i, (m, lo, hi) in enumerate(zip(means, ci_lo, ci_hi, strict=False)):
        col = TUFTE_SLATE if m < 0 else TUFTE_CRIMSON
        ax_b.plot([lo, hi], [i, i], color=col, lw=2.6, zorder=3)
        ax_b.scatter([m], [i], color=col, s=55, zorder=4)
        if m < 0:
            ax_b.text(
                0.012,
                i,
                f"{m:+.2f} [{lo:+.2f}, {hi:+.2f}]",
                va="center",
                ha="left",
                fontsize=9.1,
                fontweight="bold",
                color=col,
            )
        else:
            ax_b.text(
                hi + 0.012,
                i,
                f"{m:+.2f} [{lo:+.2f}, {hi:+.2f}]",
                va="center",
                ha="left",
                fontsize=9.1,
                fontweight="bold" if i in (4, 5) else "normal",
                color=col,
            )

    ax_b.set_yticks(y_b)
    ax_b.set_yticklabels(traits, fontsize=9.1, color=TUFTE_INK)
    ax_b.set_xlim(-0.36, 0.56)
    ax_b.set_ylim(-0.6, 5.6)
    ax_b.set_xticks([-0.30, -0.15, 0.0, 0.15, 0.30])
    ax_b.set_xlabel(
        "Pooled Within-Experiment Trait Shift \u0394z (Increased \u2212 Decreased)",
        fontsize=9.8,
        color=TUFTE_INK,
    )
    apply_range_frame(ax_b, x_bounds=(-0.30, 0.30), y_bounds=None)
    title_panel(
        ax_b,
        "B. Six Focal Trait Shifts Across 935 Contrasts (548 Studies)",
        "Aerotolerance, salt tolerance & coccus morphology rise; cellulose catabolism falls",
    )

    return save_and_verify_plate(fig, 3, "rank_coverage_and_pooled_contrasts", out_dir)


def build_plate_04(out_dir: Path) -> PlateInfo:
    """Plate IV: Three Complementary Statistical Arms, Sample Sizes & Rank Loss Cascade."""
    fig = new_plate(
        plate_num=4,
        module_tag="Module 2 — Background",
        role_tag="Three Arms & Rank Cascade",
        title="Testing the Same Taxon\u2013Trait Data Across Three Biological Scales (548 Studies)",
        subtitle=(
            "We test shifts within experiments (Arm 1: 935 contrasts), across unique taxa (Arm 2: "
            "5,829 full / 2,216 Tier 1+2 taxa; 1,870 sigs), "
            "and across studies (Arm 3: 548 studies; 55 FDR \u2264 0.05 signals vs. 12 at "
            "genus-only)."
        ),
        footer=(
            "Source: MTP Three-Arm Design on Tier 1+2 (Arm 1: 935 contrasts; Arm 2: 5,829 matrix "
            "taxa, 592,790 tests; Arm 3: 548 studies; "
            "all-rank coverage 89.2% unique / 96.4% occurrence-weighted across 317 traits)."
        ),
    )
    ax_a, ax_b = two_panels(fig)

    # Panel A: Clean 3-Arm Specification Table
    setup_canvas_panel(ax_a)
    title_panel(
        ax_a,
        "A. Three Statistical Arms Over 548 Studies & 935 Contrasts",
        "Each arm guards against a specific failure mode of microbiome meta-analysis",
    )
    headers = [
        (0.01, "Statistical Arm", "left"),
        (0.28, "Sampling Unit & Test", "left"),
        (0.62, "Primary Question Answered", "left"),
    ]
    rows = [
        (
            "Arm 1: Paired\nComposition Shift",
            "n = 935 contrasts (897 active)\nWilcoxon + (1|study_id)",
            "Does the continuous trait score\nshift between Inc and Dec lists?",
        ),
        (
            "Arm 2: Taxon-Urn\nHypergeometric ORA",
            "M = 5,829 matrix taxa\n1,870 sigs · 592,790 tests",
            "Is each signature enriched above\nthe finite taxon-background urn?",
        ),
        (
            "Arm 3a: Study Vote\nSign Test",
            "n = 548 studies (1 vote)\nExact binomial + CI",
            "Do a majority of independent\npapers agree on shift direction?",
        ),
        (
            "Arm 3b: Paired\nLog-Ratio Test",
            "n = 548 studies / 935 pairs\nPaired log\u2082(Trait A / B)",
            "Does the shift survive scale-\ninvariant compositional closure?",
        ),
    ]
    y_bot = draw_tufte_table(
        ax_a,
        headers,
        rows,
        y_top=0.90,
        row_step=0.145,
        font_size=9.1,
        header_size=9.1,
        row_accents={0: TUFTE_CRIMSON, 1: TUFTE_SLATE, 2: TUFTE_FOREST, 3: TUFTE_AMBER},
    )
    draw_marginal_note(
        ax_a,
        x=0.01,
        y_top=y_bot - 0.035,
        kicker="Rank Loss Cascade Across 317 Traits (All Ranks vs. Genus-Only)",
        body_lines=[
            "\u2022 All ranks (Class/Phylum/Order/Genus/Family/Species): 89.2% unique / 96.4% "
            "occ. coverage.",
            "\u2022 All-rank rollup detects 55 FDR \u2264 0.05 signals vs. only 12 signals when "
            "restricted to Genus-only.",
        ],
        accent_color=TUFTE_CRIMSON,
        line_step=0.048,
        font_size=8.9,
    )

    # Panel B: Log10 Effective Sample Size Scale
    style_tufte_axis(ax_b, keep_bottom=True, keep_left=False)
    units = [
        "Arm 3 Curated Studies\n(One collapsed vote per study, N = 548)",
        "Arm 1 Paired Contrasts\n(935 Tier 1+2 contrasts · 897 active)",
        "Arm 1 Directional Signatures\n(1,870 Increased + Decreased lists)",
        "Arm 2 Unique Taxon Urn\n(Deduplicated Tier 1+2 taxa, N = 2,216)",
        "Arm 2 Total ORA Grid Tests\n(1,870 signatures \u00d7 317 traits = 592,790)",
    ]
    n_vals = np.array([548, 935, 1870, 2216, 592790], dtype=float)
    log_vals = np.log10(n_vals)
    cols = [TUFTE_FOREST, TUFTE_CRIMSON, TUFTE_CRIMSON, TUFTE_SLATE, TUFTE_AMBER]
    y_pos = np.arange(len(units))

    ax_b.hlines(y_pos, 2.0, log_vals, color=cols, lw=2.6, zorder=3)
    ax_b.scatter(log_vals, y_pos, color=cols, s=65, zorder=4)
    for y, lv, nv, col in zip(y_pos, log_vals, n_vals, cols, strict=False):
        ax_b.text(
            lv + 0.12,
            float(y),
            f"N = {int(nv):,}  (10^{lv:.2f})",
            va="center",
            ha="left",
            fontsize=9.3,
            fontweight="bold",
            color=col,
        )

    ax_b.set_yticks(y_pos)
    ax_b.set_yticklabels(units, fontsize=9.1, color=TUFTE_INK)
    ax_b.set_xlim(1.8, 7.6)
    ax_b.set_ylim(-0.6, 4.6)
    ax_b.set_xticks([2, 3, 4, 5, 6])
    ax_b.set_xticklabels(
        ["10\u00b2 (100)", "10\u00b3 (1k)", "10\u00b4 (10k)", "10\u00b5 (100k)", "10\u2076 (1M)"]
    )
    ax_b.set_xlabel(
        "Effective Sample Size / Multiplicity Scale (log\u2081\u2080 N)",
        fontsize=9.8,
        color=TUFTE_INK,
    )
    apply_range_frame(ax_b, x_bounds=(2.0, 6.0), y_bounds=None)
    title_panel(
        ax_b,
        "B. Sample Size by Test Across 548 Studies (10\u00b2 to 10\u2076)",
        "From 548 studies and 935 paired contrasts to 592,790 hypergeometric ORA tests",
    )

    return save_and_verify_plate(fig, 4, "three_statistical_arms", out_dir)


def build_plate_05(out_dir: Path) -> PlateInfo:
    """Plate V: Four Meanings of 'Direction' and the 1,872-Cell Evidence Contract."""
    fig = new_plate(
        plate_num=5,
        module_tag="Module 2 — Background",
        role_tag="Evidence Contract & Direction",
        title="Measuring Trait Direction & Separating Missing Data from Nulls (935 Contrasts)",
        subtitle=(
            "Across 935 Tier 1+2 contrasts in 548 studies (897 non-zero active pairs; 94.8% genus "
            "occ. coverage), "
            "four complementary direction estimators and a 1,872-cell display contract separate "
            "sparsity from nulls."
        ),
        footer=(
            "Source: Condition \u00d7 Body-Site stratification grid (36 Conditions \u00d7 52 Body "
            "Sites = 1,872 total cells) "
            "across 935 Tier 1+2 contrasts in 548 studies (897 non-zero active pairs)."
        ),
    )
    ax_a, ax_b = two_panels(fig)

    # Panel A: 4 Meanings of Direction Table
    setup_canvas_panel(ax_a)
    title_panel(
        ax_a,
        "A. Four Ways to Quantify 'Direction' (897 Active Pairs)",
        "Each metric compares against its own null baseline across 548 Tier 1+2 studies",
    )
    headers = [
        (0.01, "Directional Metric", "left"),
        (0.31, "Mathematical Form", "left"),
        (0.62, "Null Baseline & Scope", "left"),
    ]
    rows = [
        (
            "1. Margin Prevalence",
            "P(trait > 0 in Inc)\nvs P(trait > 0 in Dec)",
            "Trait density baseline across\n1,870 directional signatures",
        ),
        (
            "2. Paired Shift \u0394",
            "Mean(f_Inc \u2212 f_Dec)\nwithin contrast i",
            "\u0394 = 0.00 null across 935 pairs\n(897 non-zero active pairs)",
        ),
        (
            "3. Study Vote %",
            "Count(\u0394_study > 0) /\nNon-tied studies",
            "50.0% coin flip across 548\nTier 1+2 case-control studies",
        ),
        (
            "4. Paired Log\u2082 Ratio",
            "Mean log\u2082((f_A+\u03b5) /\n(f_B+\u03b5)) per pair",
            "0.0 log-ratio; closure-free\n(e.g., FA / Cellulose = +1.05)",
        ),
    ]
    y_bot = draw_tufte_table(
        ax_a,
        headers,
        rows,
        y_top=0.90,
        row_step=0.145,
        font_size=9.1,
        header_size=9.1,
        row_accents={0: TUFTE_MUTED, 1: TUFTE_CRIMSON, 2: TUFTE_SLATE, 3: TUFTE_FOREST},
    )
    draw_marginal_note(
        ax_a,
        x=0.01,
        y_top=y_bot - 0.035,
        kicker="Active Pair Coverage Across 935 Tier 1+2 Contrasts",
        body_lines=[
            "\u2022 897 of 935 Tier 1+2 contrasts (95.9%) form non-zero active pairs with mapped "
            "traits.",
            "\u2022 Anchored by 94.8% genus occurrence coverage (8,139 / 8,582 genus occurrences).",
        ],
        accent_color=TUFTE_SLATE,
        line_step=0.048,
        font_size=9.0,
    )

    # Panel B: 1,872-Cell Display-State Partition Bar Chart
    style_tufte_axis(ax_b, keep_bottom=True, keep_left=False)
    states = [
        "Tested Null (\u0394 = 0)\n(Sampled without signal, |\u0394| \u2264 0.01)",
        "Masked (<2 studies)\n(Low-study cells suppressed from plot)",
        "Drawn (\u22652 studies)\n(Active Condition \u00d7 Body-Site signal)",
        "Unsampled (0 sigs)\n(Zero BugSigDB experiments in cell)",
    ]
    cell_counts = [7, 209, 316, 1340]
    cell_pcts = [0.4, 11.2, 16.9, 71.6]
    cols = [TUFTE_AMBER, TUFTE_MUTED, TUFTE_CRIMSON, TUFTE_SAND]
    y_pos = np.arange(len(states))

    ax_b.barh(y_pos, cell_counts, height=0.52, color=cols, edgecolor=TUFTE_LIGHT, lw=0.6)
    for y, c, p, col in zip(y_pos, cell_counts, cell_pcts, cols, strict=False):
        txt_col = TUFTE_INK if col == TUFTE_SAND else col
        ax_b.text(
            c + 30,
            float(y),
            f"{c:,} cells  ({p:.1f}%)",
            va="center",
            ha="left",
            fontsize=9.5,
            fontweight="bold" if y == 2 else "normal",
            color=txt_col,
        )

    ax_b.set_yticks(y_pos)
    ax_b.set_yticklabels(states, fontsize=9.3, color=TUFTE_INK)
    ax_b.set_xlim(0, 1780)
    ax_b.set_ylim(-0.6, 3.6)
    ax_b.set_xticks([0, 316, 600, 1000, 1340])
    ax_b.set_xlabel(
        "Number of Condition \u00d7 Body-Site Grid Cells (36 \u00d7 52 = 1,872)",
        fontsize=9.8,
        color=TUFTE_INK,
    )
    apply_range_frame(ax_b, x_bounds=(0, 1340), y_bounds=None)
    title_panel(
        ax_b,
        "B. Coverage of 1,872 Grid Cells (36 Conditions \u00d7 52 Sites)",
        "Prevents conflating unsampled anatomical pairs (71.6%) with tested nulls (0.4%)",
    )

    return save_and_verify_plate(fig, 5, "direction_and_evidence_contract", out_dir)


def build_plate_06(out_dir: Path) -> PlateInfo:
    """Plate VI: Four Statistical Threats & Near-Twin Semantic Redundancy."""
    fig = new_plate(
        plate_num=6,
        module_tag="Module 3 — Background",
        role_tag="Threats & Controls",
        title="Four Statistical Threats in Knowledge-Graph Trait Analysis & Their Controls",
        subtitle=(
            "Why naive enrichment over literature-curated knowledge graphs fails\u2014and how "
            "pairing, "
            "Benjamini\u2013Yekutieli FDR, log-ratios, and mixed models neutralize each threat "
            "across 935 contrasts."
        ),
        footer=(
            "Source: microbial-trait-profiling methodological audit (Fig S3 annotation density; "
            "trait-pair Jaccard "
            "co-occurrence matrix across 317 non-MediaDive features in 548 studies / 935 Tier 1+2 "
            "contrasts)."
        ),
    )
    ax_a, ax_b = two_panels(fig)

    # Panel A: 4 Threats -> Controls Table
    setup_canvas_panel(ax_a)
    title_panel(
        ax_a,
        "A. Four Methodological Threats & Built-In Remedies",
        "Every statistical claim across 935 Tier 1+2 contrasts (548 studies) survives all four "
        "controls",
    )
    headers = [
        (0.01, "Statistical Threat", "left"),
        (0.30, "Failure Mechanism", "left"),
        (0.64, "Implemented Control", "left"),
    ]
    rows = [
        (
            "1. Curation-Density\nBias (Fig S3)",
            "Pathogens have 10\u00d7 more\nBacDive edges than commensals",
            "Within-study Inc\u2013Dec pairing\n(57\u00d7 R\u00b2 drop: 9.06%\u21920.16%)",
        ),
        (
            "2. Correlated-Trait\nFDR Runaway",
            "Ontology parent-child & near-\ntwins violate BH independence",
            "Benjamini\u2013Yekutieli (BY) FDR\n+ 100\u00d7 within-study label flips",
        ),
        (
            "3. Compositional\nClosure",
            "Unit-sum normalization couples\ngains to apparent losses",
            "Scale-invariant paired log\u2082\nratios across trait pairs",
        ),
        (
            "4. Multi-Contrast\nPseudoreplication",
            "548 studies contribute 935\ncontrasts (up to 18 per paper)",
            "Mixed-effects (1|study_id) +\nArm 3 548-study sign test",
        ),
    ]
    y_bot = draw_tufte_table(
        ax_a,
        headers,
        rows,
        y_top=0.90,
        row_step=0.145,
        font_size=8.8,
        header_size=9.0,
        row_accents={0: TUFTE_CRIMSON, 1: TUFTE_SLATE, 2: TUFTE_FOREST, 3: TUFTE_AMBER},
    )
    draw_marginal_note(
        ax_a,
        x=0.01,
        y_top=y_bot - 0.035,
        kicker="Why Standard Benjamini\u2013Hochberg (BH) Is Not Enough",
        body_lines=[
            "\u2022 Across 317 traits and 2,216 Tier 1+2 taxa, ontology hierarchy and co-curation "
            "create overlap.",
            "\u2022 Benjamini\u2013Yekutieli (BY) guarantees FDR control under arbitrary "
            "dependency.",
        ],
        accent_color=TUFTE_SLATE,
        line_step=0.048,
        font_size=9.2,
    )

    # Panel B: Near-Twin Semantic Redundancy Bar Chart
    style_tufte_axis(ax_b, keep_bottom=True, keep_left=False)
    twin_pairs = [
        "Spore-forming \u2194 Gram-positive\n(38.5% shared carrier taxa)",
        "Obligate anaerobe \u2194 Catalase-negative\n(44.1% co-regulated oxygen sensitivity)",
        "Cellulose \u2194 Xylan catabolism\n(48.2% plant-polysaccharide co-carriage)",
        "Aerobe \u2194 Facultative aerobe\n(62.6% overlapping respiratory guilds)",
        "Coccus shaped \u2194 Sphere shaped\n(88.4% METPO / BacDive synonymy)",
        "pH range mid1 \u2194 pH range mid2\n(95.9% adjacent neutrophilic bins)",
    ]
    jaccards = [38.5, 44.1, 48.2, 62.6, 88.4, 95.9]
    cols = [TUFTE_FOREST, TUFTE_SLATE, TUFTE_SLATE, TUFTE_AMBER, TUFTE_CRIMSON, TUFTE_CRIMSON]
    y_pos = np.arange(len(twin_pairs))

    ax_b.barh(y_pos, jaccards, height=0.52, color=cols, alpha=0.88)
    for y, j, col in zip(y_pos, jaccards, cols, strict=False):
        ax_b.text(
            j + 1.8,
            float(y),
            f"{j:.1f}% shared taxa",
            va="center",
            ha="left",
            fontsize=9.5,
            fontweight="bold" if j > 85 else "normal",
            color=col,
        )

    ax_b.set_yticks(y_pos)
    ax_b.set_yticklabels(twin_pairs, fontsize=9.2, color=TUFTE_INK)
    ax_b.set_xlim(0, 122)
    ax_b.set_ylim(-0.6, 5.6)
    ax_b.set_xticks([0, 25, 50, 75, 100])
    ax_b.set_xlabel("Pairwise Taxon-Carrier Jaccard Overlap (%)", fontsize=10.0, color=TUFTE_INK)
    apply_range_frame(ax_b, x_bounds=(0, 100), y_bounds=None)
    title_panel(
        ax_b,
        "B. Measured Near-Twin Redundancy in KG-Microbe (317 Traits)",
        "Up to 95.9% carrier overlap between adjacent bins motivates BY FDR & guild grouping",
    )

    return save_and_verify_plate(fig, 6, "statistical_threats_and_redundancy", out_dir)


def build_plate_07(out_dir: Path) -> PlateInfo:
    """Plate VII: Curation-Density Neutralization & Unique-Taxon ORA Calibration."""
    fig = new_plate(
        plate_num=7,
        module_tag="Module 3 — Results",
        role_tag="Density & ORA Calibration",
        title=(
            "Within-Study Pairing Collapses Curation Bias 57\u00d7; Taxon Urn Restores Uniform Null"
        ),
        subtitle=(
            "Pairing Increased and Decreased signatures within each contrast collapses "
            "annotation-density "
            "variance 57\u00d7 (from 9.06% [slope \u22120.237] to 0.16%), while unique-taxon ORA "
            "(M = 5,829) restores a uniform null."
        ),
        footer=(
            "Source: Slide 30 (v4_slide_030_fig_s3_annotation_density; 935 Tier 1+2 contrasts in "
            "548 studies) "
            "and Arm 2 ORA calibration (592,790 tests across M = 5,829 full-matrix / 2,216 Tier "
            "1+2 unique taxa)."
        ),
    )
    ax_a, ax_b = two_panels(fig)

    # Panel A: Dumbbell comparison of Raw vs Standardized/Residualized shifts + 57x R^2 collapse
    style_tufte_axis(ax_a, keep_bottom=True, keep_left=False)
    traits = [
        "Cellulose catabolism (\u2193 61.1%, 414:264, q = 1.9e-7)",
        "Obligate anaerobe (\u2193 52.8%, 444:397, q = 0.24)",
        "Arsenate reduction (\u2191 54.6%, 422:351, q = 0.044)",
        "Moderately halophilic (\u2191 57.9%, 477:347, q = 6.9e-4)",
        "Coccus shaped (\u2191 55.5%, 433:347, q = 0.025)",
        "Facultatively aerobic (\u2191 56.3%, 480:373, q = 0.011)",
    ]
    raw_delta = np.array([-0.18, -0.15, 0.10, 0.13, 0.14, 0.16])
    res_delta = np.array([-0.25, -0.10, 0.24, 0.13, 0.17, 0.21])
    y_pos = np.arange(len(traits))

    ax_a.axvline(0.0, color=TUFTE_INK, lw=0.8, ls="--")
    for y, r_val, p_val in zip(y_pos, raw_delta, res_delta, strict=False):
        ax_a.plot([r_val, p_val], [y, y], color=TUFTE_INK, lw=2.0, zorder=2)
    ax_a.scatter(raw_delta, y_pos, color=TUFTE_CRIMSON, s=60, zorder=4, label="Paired \u0394 (Raw)")
    ax_a.scatter(
        res_delta,
        y_pos,
        facecolors=TUFTE_PAPER,
        edgecolors=TUFTE_SLATE,
        lw=2.0,
        s=60,
        zorder=5,
        label="Standardized / Residualized \u0394z",
    )

    for y, r_val, p_val in zip(y_pos, raw_delta, res_delta, strict=False):
        hi = max(r_val, p_val)
        if r_val < 0:
            ax_a.text(
                0.012,
                float(y),
                f"\u0394 {r_val:+.2f} / \u0394z {p_val:+.2f}",
                va="center",
                ha="left",
                fontsize=9.0,
                fontweight="bold" if y == 0 else "normal",
                color=TUFTE_SLATE,
            )
        else:
            ax_a.text(
                hi + 0.014,
                float(y),
                f"\u0394 {r_val:+.2f} / \u0394z {p_val:+.2f}",
                va="center",
                ha="left",
                fontsize=9.0,
                color=TUFTE_INK,
            )

    ax_a.text(
        0.05,
        5.65,
        "Unpaired R\u00b2 = 9.06% (slope \u22120.237) \u2192 Paired \u0394z R\u00b2 = 0.16% "
        "(57\u00d7 drop; unpaired z = 2.25%)",
        fontsize=8.3,
        fontweight="bold",
        color=TUFTE_CRIMSON,
        va="bottom",
        ha="center",
    )

    ax_a.set_yticks(y_pos)
    ax_a.set_yticklabels(traits, fontsize=8.9, color=TUFTE_INK)
    ax_a.set_xlim(-0.30, 0.46)
    ax_a.set_ylim(-0.6, 6.15)
    ax_a.set_xticks([-0.25, -0.10, 0.0, 0.10, 0.25])
    ax_a.set_xlabel(
        "Paired Trait Shift (\u0394 and \u0394z) Across 935 Contrasts in 548 Studies",
        fontsize=9.6,
        color=TUFTE_INK,
    )
    ax_a.legend(loc="lower right", frameon=False, fontsize=8.8)
    apply_range_frame(ax_a, x_bounds=(-0.25, 0.25), y_bounds=None)
    title_panel(
        ax_a,
        "A. 57\u00d7 Curation-Density Variance Collapse via Pairing",
        "Pairing neutralizes curation bias (9.06% \u2192 0.16% R\u00b2) across all 935 Tier 1+2 "
        "contrasts",
    )

    # Panel B: ORA Urn Calibration Table
    setup_canvas_panel(ax_b)
    title_panel(
        ax_b,
        "B. Hypergeometric ORA Background Urn Calibration (Arm 2)",
        "Deduplicating the urn from 26,412 assertions to M = 5,829 unique matrix taxa restores "
        "calibration",
    )
    headers = [
        (0.01, "Calibration Metric", "left"),
        (0.38, "Legacy Assertion Urn", "left"),
        (0.69, "Unique-Taxon Urn (Fixed)", "left"),
    ]
    rows = [
        (
            "1. Urn Denominator (M)",
            "26,412 trait assertions\n(multi-counted taxa)",
            "5,829 full-matrix taxa\n(2,216 Tier 1+2 unique taxa)",
        ),
        (
            "2. Null FPR & Median p",
            "4.80\u00d7 FPR inflation\n(p_median = 0.082)",
            "1.00\u00d7 calibrated null\n(p_median = 0.543)",
        ),
        (
            "3. Evaluated ORA Space",
            "Unbounded assertion\nbackground rows",
            "592,790 tests (1,870\nsigs \u00d7 317 traits)",
        ),
        (
            "4. Disease\u00d7Feces Hits",
            "Spurious enrichment in\nwell-studied pathogens",
            "4 concordant cells pass\nstrict FDR (q = 0.0077\u20130.034)",
        ),
    ]
    y_bot = draw_tufte_table(
        ax_b,
        headers,
        rows,
        y_top=0.90,
        row_step=0.145,
        font_size=9.2,
        header_size=9.2,
        row_accents={1: TUFTE_CRIMSON, 2: TUFTE_SLATE, 3: TUFTE_FOREST},
        col_colors={1: TUFTE_MUTED, 2: TUFTE_INK},
        bold_cols={0, 2},
    )
    draw_marginal_note(
        ax_b,
        x=0.01,
        y_top=y_bot - 0.035,
        kicker="Why Assertion-Count Urns Fail in Knowledge Graphs",
        body_lines=[
            "\u2022 Counting edges (26,412) instead of unique organisms (M = 5,829) inflates FPR "
            "by 4.80\u00d7.",
            "\u2022 Deduplicating to M = 5,829 unique taxa across 1,870 signatures restores "
            "1.00\u00d7 null.",
        ],
        accent_color=TUFTE_FOREST,
        line_step=0.048,
        font_size=9.2,
    )

    return save_and_verify_plate(fig, 7, "density_neutralization_and_ora_urn", out_dir)


def build_plate_08(out_dir: Path) -> PlateInfo:
    """Plate VIII: 100x Within-Study Label-Flip Nulls & Scale-Invariant Log-Ratios."""
    fig = new_plate(
        plate_num=8,
        module_tag="Module 3 — Results",
        role_tag="Permutation & Closure",
        title="Zero Permutation False Positives (0 / 43,200) & Scale-Invariant Log-Ratio Shifts",
        subtitle=(
            "Randomly flipping Increased and Decreased labels within studies 100 times produces "
            "zero false discoveries "
            "at BH q \u2264 0.10, while paired log\u2082 trait ratios across 935 contrasts (548 "
            "studies) rule out unit-sum closure."
        ),
        footer=(
            "Source: microbial-trait-profiling permutation audit (100 iterations \u00d7 432 "
            "Condition \u00d7 Feces cells "
            "= 43,200 null tests) and Arm 3b paired log\u2082(Trait A / Trait B) closure tests "
            "across 548 studies (935 contrasts)."
        ),
    )
    ax_a, ax_b = two_panels(fig)

    # Panel A: Permutation Null vs Observed Discovery Comparison
    setup_canvas_panel(ax_a)
    title_panel(
        ax_a,
        "A. 100\u00d7 Within-Study Inc/Dec Label-Flip Benchmark",
        "Preserves exact trait covariance and study clustering across 548 studies while "
        "destroying direction",
    )
    headers = [
        (0.01, "Evaluation Regime", "left"),
        (0.36, "Cells / Traits Tested", "left"),
        (0.65, "Discoveries", "left"),
    ]
    rows = [
        (
            "1. Permuted Null\n(100\u00d7 sign flips)",
            "43,200 cell tests\n(100 \u00d7 432 Feces cells)",
            "0 / 43,200 at BH q \u2264 0.10\n(0.00% FPR; min q = 0.278)",
        ),
        (
            "2. Observed Feces\n(BH FDR \u2264 0.10)",
            "432 Condition \u00d7 Feces\ncells evaluated",
            "4 significant cells\n(PD, HIV, CRC + PD Halophily)",
        ),
        (
            "3. Observed Feces\n(Stringent BY \u2264 0.10)",
            "432 Condition \u00d7 Feces\ncells evaluated",
            "3 core cells survive BY\narbitrary-dependence FDR",
        ),
        (
            "4. Global 317-Trait\nTier 1+2 Screen",
            "317 non-MediaDive\ntraits (n = 935 contrasts)",
            "6 study-level t-test BY hits\n(2 binomial sign BY hits at m=317)",
        ),
    ]
    y_bot = draw_tufte_table(
        ax_a,
        headers,
        rows,
        y_top=0.90,
        row_step=0.145,
        font_size=9.2,
        header_size=9.2,
        row_accents={0: TUFTE_SLATE, 1: TUFTE_CRIMSON, 2: TUFTE_FOREST},
    )
    draw_marginal_note(
        ax_a,
        x=0.01,
        y_top=y_bot - 0.035,
        kicker="Why Within-Study Sign Flipping Is the Gold-Standard Null",
        body_lines=[
            "\u2022 Keeps the exact 317\u00d7317 trait correlation matrix intact inside every "
            "paper.",
            "\u2022 Proves that near-twin trait correlations cannot manufacture false positives "
            "under pairing.",
        ],
        accent_color=TUFTE_SLATE,
        line_step=0.048,
        font_size=9.2,
    )

    # Panel B: Paired Log2 Ratio Forest Plot
    # (Compositional Closure Control across 548 studies / 935 contrasts)
    style_tufte_axis(ax_b, keep_bottom=True, keep_left=False)
    ratio_labels = [
        "Anaerobe / Facultative\n(60.5% Down [300:196] · p=3.3e-6, q=2.0e-5)",
        "Bacillus / Coccus shaped\n(57.4% Down [272:202] · p=0.0015, q=0.0059)",
        "Sphere / Cellulose catabolism\n(59.9% Up [230:154] · BY q = 0.070)",
        "Halophilic / Cellulose catabolism\n(60.5% Up [228:149] · BY q = 0.070)",
        "Coccus / Cellulose catabolism\n(56.0% Up [219:172] · BY q = 0.594)",
        "Arsenate / Cellulose catabolism\n(62.5% Up [235:141] · BY q = 0.0029)",
        "Facultative / Cellulose catabolism\n(62.8% Up [240:142] · BY q = 0.0029)",
    ]
    means = np.array([-0.63, -0.45, 0.65, 0.72, 0.75, 0.91, 1.05])
    ci_lo = np.array([-0.89, -0.72, 0.25, 0.30, 0.24, 0.49, 0.62])
    ci_hi = np.array([-0.37, -0.18, 1.05, 1.14, 1.26, 1.33, 1.48])
    y_pos = np.arange(len(ratio_labels))

    ax_b.axvline(0.0, color=TUFTE_INK, lw=0.9, ls="--")
    for i, (m, lo, hi) in enumerate(zip(means, ci_lo, ci_hi, strict=False)):
        col = TUFTE_SLATE if m < 0 else TUFTE_CRIMSON
        ax_b.plot([lo, hi], [i, i], color=col, lw=2.4, zorder=3)
        ax_b.scatter([m], [i], color=col, s=50, zorder=4)
        x_lbl = 0.05 if m < 0 else hi + 0.05
        ax_b.text(
            x_lbl,
            i,
            f"{m:+.2f} [{lo:+.2f}, {hi:+.2f}]",
            va="center",
            ha="left",
            fontsize=8.7,
            fontweight="bold" if i in (0, 5, 6) else "normal",
            color=col,
        )

    ax_b.set_yticks(y_pos)
    ax_b.set_yticklabels(ratio_labels, fontsize=8.6, color=TUFTE_INK)
    ax_b.set_xlim(-1.15, 2.35)
    ax_b.set_ylim(-0.6, 6.6)
    ax_b.set_xticks([-1.0, -0.5, 0.0, 0.5, 1.0, 1.5])
    ax_b.set_xlabel(
        "Paired Shift in log\u2082(Trait A / Trait B) Across 548 Studies (95% CI)",
        fontsize=9.6,
        color=TUFTE_INK,
    )
    apply_range_frame(ax_b, x_bounds=(-1.0, 1.5), y_bounds=None)
    title_panel(
        ax_b,
        "B. Scale-Invariant Paired Log\u2082 Trait Ratios (548 Studies)",
        "2 of 15 all-pairs pass BY q\u22640.05 (Facultative/Cellulose +1.05, Arsenate/Cellulose "
        "+0.91, q=0.0029)",
    )

    return save_and_verify_plate(fig, 8, "study_flip_null_and_log_ratios", out_dir)


def build_plate_09(out_dir: Path) -> PlateInfo:
    """Plate IX: Unbiased 317-Trait Screening & the Epithelial Oxygenation Model."""
    fig = new_plate(
        plate_num=9,
        module_tag="Module 4 — Background",
        role_tag="Redox Hypothesis & Screen",
        title="The Epithelial Oxygenation Model & Multi-Arm 317-Trait Convergence Funnel",
        subtitle=(
            "Testing the Winter & Byndloss luminal oxygenation hypothesis across all 317 "
            "non-MediaDive features "
            "in 935 Tier 1+2 contrasts (548 studies) without pre-selecting candidate pathways."
        ),
        footer=(
            "Source: Winter & Byndloss (Science 2013; Nat Rev Microbiol 2017) redox framework "
            "tested across "
            "317 KG-Microbe features in 935 Tier 1+2 subject-level case-vs-control contrasts (548 "
            "studies)."
        ),
    )
    ax_a, ax_b = two_panels(fig)

    # Panel A: Healthy vs Inflamed Lumen Ecological Comparison Table
    setup_canvas_panel(ax_a)
    title_panel(
        ax_a,
        "A. Ecological Mechanism: Luminal Oxygen & Nitrate Leak",
        "How host colonocyte metabolism governs the balance between fiber fermenters and "
        "pathobionts",
    )
    headers = [
        (0.01, "Ecological Dimension", "left"),
        (0.30, "Healthy Homeostatic Colon", "left"),
        (0.64, "Inflamed / Dysbiotic Colon", "left"),
    ]
    rows = [
        (
            "1. Colonocyte Energy",
            "Mitochondrial \u03b2-oxidation\nof microbiota butyrate",
            "Shift to anaerobic glycolysis;\nepithelial hypoxia lost",
        ),
        (
            "2. Luminal Electron\nAcceptors",
            "Strict luminal anoxia\n(pO\u2082 < 1 mmHg; low NO\u2083\u207b)",
            "O\u2082 diffuses from mucosa;\niNOS generates NO\u2083\u207b / S-oxides",
        ),
        (
            "3. Winning Guild",
            "Obligate anaerobes &\nplant-fiber degraders",
            "Facultative anaerobes,\nhalophiles & Enterobacteriaceae",
        ),
        (
            "4. Predicted Trait\nSignature",
            "High cellulose/xylan catabolism\n& catalase-negative anaerobes",
            "Rise in facultative respiration,\narsenate/nitrate redox & cocci",
        ),
    ]
    y_bot = draw_tufte_table(
        ax_a,
        headers,
        rows,
        y_top=0.90,
        row_step=0.145,
        font_size=8.7,
        header_size=8.9,
        row_accents={1: TUFTE_SLATE, 2: TUFTE_FOREST, 3: TUFTE_CRIMSON},
        col_colors={1: TUFTE_FOREST, 2: TUFTE_CRIMSON},
        bold_cols={0},
    )
    draw_marginal_note(
        ax_a,
        x=0.01,
        y_top=y_bot - 0.035,
        kicker="Coupled Oxygen McNemar Test Across 935 Tier 1+2 Contrasts",
        body_lines=[
            "\u2022 Across 935 contrasts (897 active): 343 Aerobe-Up & Anaerobe-Down vs 231 "
            "reverse (OR = 1.48, p = 3.4e-6).",
            "\u2022 Feces subset (n = 547): 261 vs 134 (OR = 1.95, p = 1.9e-10, polarity +0.71).",
        ],
        accent_color=TUFTE_CRIMSON,
        line_step=0.048,
        font_size=9.0,
    )

    # Panel B: 5-Stage Cross-Arm Trait Convergence Funnel
    style_tufte_axis(ax_b, keep_bottom=True, keep_left=False)
    tiers = [
        "Stage 5: Study Sign BY \u2264 0.05 (m = 317)\n(2 core: Cellulose & Facultative, q=0.041)",
        "Stage 4: Study Sign BY \u2264 0.05 (m_eff = 75)\n(3 core traits survive eigenvalue FDR)",
        "Stage 3: Study Paired t-Test BY \u2264 0.05\n(6 traits pass study continuous BY FDR)",
        "Stage 2: Study Paired t-Test BH \u2264 0.05\n(22 traits pass study BH FDR, 548 studies)",
        "Stage 1: Non-MediaDive Feature Space\n(317 features: 135M + 102C + 68P + 12G)",
    ]
    counts = np.array([2, 3, 6, 22, 317], dtype=float)
    cols = [TUFTE_CRIMSON, TUFTE_CRIMSON, TUFTE_SLATE, TUFTE_FOREST, TUFTE_MUTED]
    y_pos = np.arange(len(tiers))

    ax_b.barh(y_pos, counts, height=0.52, color=cols, alpha=0.88)
    for y, c, col in zip(y_pos, counts, cols, strict=False):
        pct_317 = (c / 317.0) * 100.0
        ax_b.text(
            c + 6.0,
            float(y),
            f"{int(c):,} traits  ({pct_317:.1f}% of 317)",
            va="center",
            ha="left",
            fontsize=9.3,
            fontweight="bold" if y <= 2 else "normal",
            color=col,
        )

    ax_b.set_yticks(y_pos)
    ax_b.set_yticklabels(tiers, fontsize=9.1, color=TUFTE_INK)
    ax_b.set_xlim(0.0, 435.0)
    ax_b.set_ylim(-0.6, 4.6)
    ax_b.set_xticks([0, 22, 75, 160, 317])
    ax_b.set_xlabel(
        "Number of Features Retained Across Stringency Tiers (out of 317)",
        fontsize=9.8,
        color=TUFTE_INK,
    )
    apply_range_frame(ax_b, x_bounds=(0.0, 317.0), y_bounds=None)
    title_panel(
        ax_b,
        "B. Study-Level Convergence Funnel (548 Studies / 935 Contrasts)",
        "317 tested \u2192 22 Study BH \u2192 6 Study t-test BY \u2192 3 Sign BY (m_eff=75) "
        "\u2192 2 Sign BY (m=317)",
    )

    return save_and_verify_plate(fig, 9, "global_screen_and_redox_model", out_dir)


def build_plate_10(out_dir: Path) -> PlateInfo:
    """Plate X: Global 317-Trait Screen Results & Arm 3 Study Majority Votes."""
    fig = new_plate(
        plate_num=10,
        module_tag="Module 4 — Results",
        role_tag="317-Trait Screen & Study Votes",
        title="Cellulose Breakdown Ranks #1 in Global Losses; Study Votes Confirm Focal Shifts",
        subtitle=(
            "Across all 317 non-MediaDive features in 935 Tier 1+2 contrasts (548 studies), "
            "cellulose catabolism is the #1 "
            "depleted trait (61.1% down, 414:264, p = 9.7\u00d710\u207b\u2079, q = "
            "1.9\u00d710\u207b\u2077), confirmed by 548-study majority votes."
        ),
        footer=(
            "Source: Global 317-trait paired screen across 935 Tier 1+2 contrasts (897 active "
            "pairs in 548 studies) and "
            "Arm 3 study-level majority-vote sign tests across 548 independent studies."
        ),
    )
    ax_a, ax_b = two_panels(fig)

    # Panel A: Six Focal Signals in the Unbiased 317-Trait Screen (-log10 p)
    style_tufte_axis(ax_a, keep_bottom=True, keep_left=False)
    top_traits = [
        "Obligate anaerobe\n(\u2193 52.8% [444:397] · \u0394=\u22120.15, q=0.24)",
        "Arsenate reduction\n(\u2191 54.6% [422:351] · \u0394=+0.10, q=0.044)",
        "Coccus shaped\n(\u2191 55.5% [433:347] · \u0394=+0.14, q=0.025)",
        "Facultatively aerobic\n(\u2191 56.3% [480:373] · \u0394=+0.16, q=0.011)",
        "Moderately halophilic\n(\u2191 57.9% [477:347] · \u0394=+0.13, q=6.9e-4)",
        "Cellulose catabolism\n(\u2193 #1 61.1% [414:264] · \u0394=\u22120.18, q=1.9e-7)",
    ]
    neg_log_p = np.array([0.96, 1.96, 2.64, 3.55, 5.18, 8.01])
    p_strs = [
        "p = 0.11 (444:397 \u2193)",
        "p = 0.011 (422:351 \u2191)",
        "p = 0.0023 (433:347 \u2191)",
        "p = 2.8e-4 (480:373 \u2191)",
        "p = 6.6e-6 (477:347 \u2191)",
        "p = 9.7e-9 (414:264 \u2193 #1)",
    ]
    cols = [
        TUFTE_SLATE,
        TUFTE_CRIMSON,
        TUFTE_CRIMSON,
        TUFTE_CRIMSON,
        TUFTE_CRIMSON,
        TUFTE_SLATE,
    ]
    y_pos = np.arange(len(top_traits))

    ax_a.axvline(2.0, color=TUFTE_MUTED, lw=0.9, ls=":", zorder=1)
    ax_a.text(
        2.05, 0.1, "BH q \u2264 0.05", fontsize=8.4, color=TUFTE_MUTED, ha="left", va="bottom"
    )

    ax_a.hlines(y_pos, 0.0, neg_log_p, color=cols, lw=2.6, zorder=3)
    ax_a.scatter(neg_log_p, y_pos, color=cols, s=58, zorder=4)
    for y, v, ps, col in zip(y_pos, neg_log_p, p_strs, cols, strict=False):
        ax_a.text(
            v + 0.16,
            float(y),
            ps,
            va="center",
            ha="left",
            fontsize=8.9,
            fontweight="bold" if y in (4, 5) else "normal",
            color=col,
        )

    ax_a.set_yticks(y_pos)
    ax_a.set_yticklabels(top_traits, fontsize=8.9, color=TUFTE_INK)
    ax_a.set_xlim(0.0, 11.4)
    ax_a.set_ylim(-0.6, 5.6)
    ax_a.set_xticks([0, 2, 4, 6, 8])
    ax_a.set_xlabel(
        "Paired Sign / Wilcoxon Strength  \u2212log\u2081\u2080(p)  Across 935 Contrasts",
        fontsize=9.6,
        color=TUFTE_INK,
    )
    apply_range_frame(ax_a, x_bounds=(0.0, 8.0), y_bounds=None)
    title_panel(
        ax_a,
        "A. Six Focal Guilds in 317-Trait Screen (n = 935 Contrasts)",
        "Cellulose catabolism (#1 global signal, q = 1.9e-7) and halophilic/facultative gains lead",
    )

    # Panel B: Arm 3 Study-Level Majority Vote Forest Plot across 548 Studies
    style_tufte_axis(ax_b, keep_bottom=True, keep_left=False)
    arm3_traits = [
        "Cellulose catabolism\n(256:153 [62.6% Down, q=1.5e-4] / 181:265 [BY q=0.0408])",
        "Obligate anaerobe\n(264:235 [52.9% Down] · paired ratio 60.5% Down, q=2.0e-5)",
        "Coccus shaped\n(240:196 [55.0%, p=0.039] / 280:215 [56.6% Up])",
        "Facultatively aerobic\n(285:210 [57.6% Up, q=0.0022] / 295:219 [BY q=0.0408])",
        "Moderately halophilic\n(272:194 [58.4% Up, q=0.0014])",
        "Arsenate reduction\n(243:197 [55.2%, p=0.032] / 206:141 [59.4% Up])",
    ]
    pct_up = np.array([37.4, 47.1, 56.6, 57.6, 58.4, 59.4])
    ci_lo = np.array([32.7, 42.6, 52.1, 53.1, 53.8, 54.0])
    ci_hi = np.array([42.3, 51.6, 61.0, 62.0, 62.9, 64.6])
    y_b = np.arange(len(arm3_traits))

    ax_b.axvline(50.0, color=TUFTE_INK, lw=0.9, ls="--", zorder=1)
    for i, (m, lo, hi) in enumerate(zip(pct_up, ci_lo, ci_hi, strict=False)):
        col = TUFTE_SLATE if m < 50.0 else TUFTE_CRIMSON
        ax_b.plot([lo, hi], [i, i], color=col, lw=2.6, zorder=3)
        ax_b.scatter([m], [i], color=col, s=58, zorder=4)
        ax_b.text(
            hi + 0.8,
            i,
            f"{m:.1f}% [{lo:.1f}%, {hi:.1f}%]",
            va="center",
            ha="left",
            fontsize=8.9,
            fontweight="bold" if i in (0, 3, 4, 5) else "normal",
            color=col,
        )

    ax_b.set_yticks(y_b)
    ax_b.set_yticklabels(arm3_traits, fontsize=8.6, color=TUFTE_INK)
    ax_b.set_xlim(30.0, 77.0)
    ax_b.set_ylim(-0.6, 5.6)
    ax_b.set_xticks([35, 40, 45, 50, 55, 60, 65])
    ax_b.set_xlabel(
        "Percentage of Independent Studies Shifting Increased (%) [95% CI]",
        fontsize=9.6,
        color=TUFTE_INK,
    )
    apply_range_frame(ax_b, x_bounds=(35.0, 65.0), y_bounds=None)
    title_panel(
        ax_b,
        "B. Arm 3 Study Majority-Vote Sign Test (548 Studies)",
        "One vote per study eliminates multi-contrast pseudoreplication; Bacillus/Coccus 57.4% "
        "down (q=0.0059)",
    )

    return save_and_verify_plate(fig, 10, "global_317_results_and_arm3_votes", out_dir)
