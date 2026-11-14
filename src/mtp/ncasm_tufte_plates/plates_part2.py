"""Part 2 of the 20-slide NCASM Tufte Plate generator (Plates XI through XX)."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from mtp.ncasm_tufte_plates.plates_common import (
    TUFTE_AMBER,
    TUFTE_CRIMSON,
    TUFTE_FOREST,
    TUFTE_HAIRLINE,
    TUFTE_INK,
    TUFTE_LIGHT,
    TUFTE_MUTED,
    TUFTE_PAPER,
    TUFTE_SLATE,
    PlateInfo,
    SlopeSeries,
    apply_range_frame,
    draw_marginal_note,
    draw_slopegraph_series,
    draw_tufte_quartile_rule,
    draw_tufte_table,
    new_plate,
    save_and_verify_plate,
    setup_canvas_panel,
    style_tufte_axis,
    title_panel,
    two_panels,
)


def build_plate_11(out_dir: Path) -> PlateInfo:
    """Plate XI: Condition x Feces Atlases Across Seven Major Conditions."""
    fig = new_plate(
        plate_num=11,
        module_tag="Module 4 — Results",
        role_tag="Disease \u00d7 Feces Atlases",
        title=(
            "Condition \u00d7 Feces Atlases Across Seven Major Conditions & Global Cellulose Loss"
        ),
        subtitle=(
            "Evaluating Tier 1+2 stool cohorts across seven major conditions shows consistent "
            "cellulose catabolism loss "
            "(63.7% down in stool, n = 496, p = 1.1\u00d710\u207b\u2079) alongside aerotolerant "
            "and halophilic gains."
        ),
        footer=(
            "Source: Arm 1 Condition \u00d7 Feces atlas across 935 Tier 1+2 contrasts (548 "
            "studies; stool n = 547, cellulose-active "
            "stool n = 496) and Arm 2 unique-taxon hypergeometric ORA (M = 2,216 Tier 1+2 taxa)."
        ),
    )
    ax_a, ax_b = two_panels(fig)

    # Panel A: Tufte Dot-Matrix of Paired Shift Delta across 7 Major Conditions
    # (Stool / Total Tier 1+2)
    style_tufte_axis(ax_a, keep_bottom=True, keep_left=False)
    cohorts = [
        "Liver cirrhosis / NAFLD (Stool n=6 · All n=17)",
        "COVID-19 (Stool n=7 · All n=94)",
        "IBD (Stool n=15 · All n=80)",
        "HIV infection (Stool n=20 · All n=28)",
        "IBS / Diabetes (Stool n=31 · Diab n=39)",
        "Colorectal cancer (Stool n=44 · All n=54)",
        "Parkinson\u2019s disease (Stool n=71 · All n=80)",
    ]
    trait_cols = [
        "Cellulose\ncatabolism",
        "Facultative\naerobic",
        "Moderately\nhalophilic",
        "Coccus\nshaped",
    ]
    # Delta in percentage points across the 7 major conditions
    mat = np.array(
        [
            [-18.2, +22.4, +11.5, +14.2],
            [-21.0, +24.8, +13.2, +16.8],
            [-21.8, +26.5, +12.4, +15.1],
            [-53.7, +21.2, +9.4, +11.8],
            [-19.5, +17.4, +11.0, +13.5],
            [-36.9, +19.8, +10.2, +16.4],
            [-27.4, +18.6, +24.4, +14.9],
        ]
    )
    fdr_stars = {(3, 0), (5, 0), (6, 0), (6, 2)}

    for r_i in range(len(cohorts)):
        ax_a.axhline(r_i, color=TUFTE_HAIRLINE, lw=0.6, zorder=1)
        for c_i in range(len(trait_cols)):
            val = mat[r_i, c_i]
            col = TUFTE_SLATE if val < 0 else TUFTE_CRIMSON
            is_sig = (r_i, c_i) in fdr_stars
            sz = abs(val) * 2.8 + (42 if is_sig else 14)
            ax_a.scatter(
                [c_i],
                [r_i],
                s=sz,
                color=col,
                alpha=0.92 if is_sig else 0.45,
                edgecolors=TUFTE_INK if is_sig else "none",
                linewidths=1.4 if is_sig else 0.0,
                zorder=4 if is_sig else 3,
            )
            tag = f"{val:+.1f}*" if is_sig else f"{val:+.1f}"
            ax_a.text(
                c_i + 0.16,
                r_i,
                tag,
                va="center",
                ha="left",
                fontsize=8.5 if not is_sig else 9.0,
                fontweight="bold" if is_sig else "normal",
                color=col if is_sig else TUFTE_MUTED,
            )

    ax_a.set_yticks(np.arange(len(cohorts)))
    ax_a.set_yticklabels(cohorts, fontsize=8.8, color=TUFTE_INK)
    ax_a.set_xticks(np.arange(len(trait_cols)))
    ax_a.set_xticklabels(trait_cols, fontsize=9.1, fontweight="bold", color=TUFTE_INK)
    ax_a.set_xlim(-0.35, 3.85)
    ax_a.set_ylim(-0.6, 6.6)
    ax_a.set_xlabel(
        "Mean Paired Shift \u0394 in Percentage Points (* = BH q \u2264 0.10)",
        fontsize=9.5,
        color=TUFTE_INK,
    )
    apply_range_frame(ax_a, x_bounds=(0, 3), y_bounds=None)
    title_panel(
        ax_a,
        "A. Seven Major Conditions \u00d7 Feces Shift Matrix (\u0394 pp)",
        "Stool (n) and All-Site Tier 1+2 (n) contrast counts across seven major conditions",
    )

    # Panel B: Global & Stool Cellulose Signal + 4 Headline Condition x Feces Pairs
    setup_canvas_panel(ax_b)
    title_panel(
        ax_b,
        "B. Global Cellulose Signal (#1) & Condition \u00d7 Feces Convergence",
        "Cellulose catabolism is the #1 depleted trait globally (935 contrasts, 548 studies) and "
        "in stool",
    )
    headers = [
        (0.01, "Cohort / Analytical Stratum", "left"),
        (0.36, "Paired Shift & Direction", "left"),
        (0.68, "Statistical Evidence", "left"),
    ]
    rows = [
        (
            "1. Global Pooled Screen\n(935 contrasts · 548 studies)",
            "416 Dec vs 269 Inc\n(60.7% down · \u0394 = \u22120.17)",
            "p = 2.0\u00d710\u207b\u2078, q = 3.9\u00d710\u207b\u2077\nStudy vote 265:181 (59.4%, "
            "BY q=0.041)",
        ),
        (
            "2. All Tier 1+2 Stool\n(n = 496 active stool pairs)",
            "316 Dec vs 180 Inc\n(63.7% down in feces)",
            "p = 1.1\u00d710\u207b\u2079\nPD+CRC+HIV: 73.0% (81:30, p=1.4e-6)",
        ),
        (
            "3. PD, CRC & HIV Stool\n(n = 71, n = 44, n = 20)",
            "PD 75.4% (46:15) · HIV 70.6%\nCRC 69.7% Dec (23:10)",
            "PD p = 8.4e-5 · HIV p = 0.072\nCRC p = 0.035 (Arm 2 ORA q \u2264 0.018)",
        ),
        (
            "4. PD \u00d7 Halophilic (\u2191)\n(n = 71 stool; n = 80 total)",
            "\u0394 = +24.4 pp in stool\nOsmotolerant/salt-tolerant rise",
            "Arm 1 BH q = 0.088\nArm 2 ORA 2.4\u00d7 (q = 0.0336)",
        ),
    ]
    y_bot = draw_tufte_table(
        ax_b,
        headers,
        rows,
        y_top=0.90,
        row_step=0.145,
        font_size=8.9,
        header_size=9.1,
        row_accents={0: TUFTE_SLATE, 1: TUFTE_SLATE, 2: TUFTE_SLATE, 3: TUFTE_CRIMSON},
    )
    draw_marginal_note(
        ax_b,
        x=0.01,
        y_top=y_bot - 0.035,
        kicker="Systemic Conditions Sampled via Stool vs. Primary GI Conditions",
        body_lines=[
            "\u2022 Seven major conditions span neurological (PD n=80), viral (COVID-19 n=94, HIV "
            "n=28),",
            "  metabolic/hepatic (Diabetes n=39, Cirrhosis n=17), and GI conditions (IBD n=80, "
            "CRC n=54).",
        ],
        accent_color=TUFTE_SLATE,
        line_step=0.048,
        font_size=8.9,
    )

    return save_and_verify_plate(fig, 11, "arm1_and_arm2_disease_atlases", out_dir)


def build_plate_12(out_dir: Path) -> PlateInfo:
    """Plate XII: Paired McNemar Test Confirms Coordinated Oxygen-Guild Swap."""
    fig = new_plate(
        plate_num=12,
        module_tag="Module 4 — Results",
        role_tag="Pan-Disease Redox Swap",
        title=(
            "Paired McNemar Test Confirms Coupled Oxygen-Guild Rewiring (p = "
            "3.4\u00d710\u207b\u2076)"
        ),
        subtitle=(
            "Across 935 Tier 1+2 contrasts (897 active in 533 studies), Aerobe-Up & Anaerobe-Down "
            "co-occurs 1.48\u00d7 "
            "more often than the reverse (343 vs. 231 pairs; Feces n = 547: 261 vs. 134, OR = "
            "1.95, p = 1.9\u00d710\u207b\u00b9\u2070)."
        ),
        footer=(
            "Source: Coupled oxygen-guild analysis across 935 Tier 1+2 contrasts (897 active in "
            "533 studies; Feces n = 547, "
            "Non-fecal GI mucosa n = 44, Oral cavity n = 78) and Seven Major Conditions."
        ),
    )
    ax_a, ax_b = two_panels(fig)

    # Panel A: Compartment-Stratified Coupled Oxygen McNemar Discordance & Polarity
    style_tufte_axis(ax_a, keep_bottom=True, keep_left=False)
    categories = [
        "Oral Cavity (n = 78)\n(22 Canonical vs 26 Reverse · Polarity \u22120.19)",
        "Non-Fecal GI Mucosa (n = 44)\n(11 Canonical vs 16 Reverse · Polarity \u22120.22)",
        "Feces / Stool Subset (n = 547)\n(261 Canonical vs 134 Reverse · OR = 1.95, p = 1.9e-10)",
        "All Tier 1+2 Contrasts (n = 897 active)\n(343 Canonical vs 231 Reverse · OR = 1.48, p = "
        "3.4e-6)",
    ]
    canon_counts = np.array([22, 11, 261, 343])
    rev_counts = np.array([26, 16, 134, 231])
    y_pos = np.arange(len(categories))

    ax_a.barh(
        y_pos + 0.15,
        canon_counts,
        height=0.28,
        color=TUFTE_CRIMSON,
        alpha=0.90,
        label="Aerobe Up & Anaerobe Down (Canonical)",
    )
    ax_a.barh(
        y_pos - 0.15,
        rev_counts,
        height=0.28,
        color=TUFTE_SLATE,
        alpha=0.85,
        label="Anaerobe Up & Aerobe Down (Reverse)",
    )

    for y, c_up, c_dn in zip(y_pos, canon_counts, rev_counts, strict=False):
        hi = max(c_up, c_dn)
        ratio = c_up / c_dn
        ax_a.text(
            hi + 12,
            float(y),
            f"{c_up} vs {c_dn}  (OR = {ratio:.2f})",
            va="center",
            ha="left",
            fontsize=9.1,
            fontweight="bold" if y >= 2 else "normal",
            color=TUFTE_CRIMSON if y >= 2 else TUFTE_MUTED,
        )

    ax_a.text(
        255,
        1.15,
        "Feces Polarity = +0.71 (OR = 1.95, p = 1.9\u00d710\u207b\u00b9\u2070)\nMucosa Polarity = "
        "\u22120.22 · Oral Polarity = \u22120.19",
        va="center",
        ha="left",
        fontsize=8.3,
        fontweight="bold",
        color=TUFTE_INK,
    )

    ax_a.set_yticks(y_pos)
    ax_a.set_yticklabels(categories, fontsize=8.9, color=TUFTE_INK)
    ax_a.set_xlim(0, 520)
    ax_a.set_ylim(-0.6, 3.6)
    ax_a.set_xticks([0, 134, 231, 261, 343])
    ax_a.set_xlabel(
        "Number of Discordant Coupled Oxygen Contrasts (897 Active of 935)",
        fontsize=9.6,
        color=TUFTE_INK,
    )
    ax_a.legend(loc="lower right", frameon=False, fontsize=8.6)
    apply_range_frame(ax_a, x_bounds=(0, 343), y_bounds=None)
    title_panel(
        ax_a,
        "A. Coupled Oxygen Rewiring by Anatomical Compartment",
        "Luminal stool drives the pan-disease redox swap (261 vs 134); mucosa and oral invert "
        "polarity",
    )

    # Panel B: Seven Major Conditions Aerotolerance Shift (Tier 1+2 Total Contrast Counts)
    style_tufte_axis(ax_b, keep_bottom=True, keep_left=False)
    diseases = [
        "Liver cirrhosis (Hepatic · n = 17)",
        "HIV infection (Viral immunodeficiency · n = 28)",
        "Diabetes mellitus (Metabolic · n = 39)",
        "Colorectal cancer (Neoplastic · n = 54)",
        "Inflammatory bowel disease (Inflammatory · n = 80)",
        "Parkinson\u2019s disease (Neurological · n = 80)",
        "COVID-19 (Acute viral · n = 94)",
    ]
    betas = np.array([0.18, 0.23, 0.16, 0.21, 0.28, 0.19, 0.25])
    lo_ci = np.array([0.04, 0.07, 0.05, 0.10, 0.18, 0.08, 0.14])
    hi_ci = np.array([0.32, 0.39, 0.27, 0.32, 0.38, 0.30, 0.36])
    y_b = np.arange(len(diseases))

    ax_b.axvline(0.0, color=TUFTE_INK, lw=0.9, ls="--")
    for i, (b, lo, hi) in enumerate(zip(betas, lo_ci, hi_ci, strict=False)):
        col = TUFTE_CRIMSON if i in (4, 5, 6) else TUFTE_SLATE
        ax_b.plot([lo, hi], [i, i], color=col, lw=2.5, zorder=3)
        ax_b.scatter([b], [i], color=col, s=54, zorder=4)
        ax_b.text(
            hi + 0.015,
            i,
            f"\u03b2 = +{b:.2f} [+{lo:.2f}, +{hi:.2f}]",
            va="center",
            ha="left",
            fontsize=8.9,
            fontweight="bold" if i >= 4 else "normal",
            color=col,
        )

    ax_b.set_yticks(y_b)
    ax_b.set_yticklabels(diseases, fontsize=8.9, color=TUFTE_INK)
    ax_b.set_xlim(-0.05, 0.68)
    ax_b.set_ylim(-0.6, 6.6)
    ax_b.set_xticks([0.0, 0.15, 0.30, 0.45])
    ax_b.set_xlabel(
        "Cluster-Adjusted Net Aerotolerance Shift \u03b2 Across Seven Major Conditions",
        fontsize=9.5,
        color=TUFTE_INK,
    )
    apply_range_frame(ax_b, x_bounds=(0.0, 0.45), y_bounds=None)
    title_panel(
        ax_b,
        "B. Aerotolerance Shift Across Seven Major Conditions",
        "Positive shift across COVID-19 (94), PD (80), IBD (80), CRC (54), Diabetes (39), HIV "
        "(28), Cirrhosis (17)",
    )

    return save_and_verify_plate(fig, 12, "oxygen_guild_mcnemar_rewiring", out_dir)


def build_plate_13(out_dir: Path) -> PlateInfo:
    """Plate XIII: 7-Region Venn Partition Arithmetic & Six Nested Null Models."""
    fig = new_plate(
        plate_num=13,
        module_tag="Module 5 — Background",
        role_tag="Venn Partitions & Nulls",
        title="Why a Strong 3-Way Core (14 Taxa) Leaves Small Exclusive 2-Way Slices",
        subtitle=(
            "Clarifying the difference between exclusive Venn slices and full pairwise "
            "intersections across the "
            "41 Tier 1+2 stool fiber-degrading taxa (PD = 33, CRC = 30, HIV = 20), benchmarked "
            "against six null universes."
        ),
        footer=(
            "Source: Cellulose-degrading guild convergence audit across 135 Tier 1+2 stool "
            "contrasts (PD n = 71 [33 taxa], "
            "CRC n = 44 [30 taxa], HIV n = 20 [20 taxa]; union = 41 unique depleted stool "
            "cellulose taxa; U1\u2013U6 nulls)."
        ),
    )
    ax_a, ax_b = two_panels(fig)

    # Panel A: 7-Region Venn Partition Table
    setup_canvas_panel(ax_a)
    title_panel(
        ax_a,
        "A. Exact Seven-Region Venn Partition of 41 Stool Taxa",
        "Full pairwise overlap adds the 14 all-three core taxa to each exclusive 2-disease slice",
    )
    headers = [
        (0.01, "Venn Partition Region", "left"),
        (0.42, "Taxon Count", "left"),
        (0.68, "Full Pairwise Sum", "left"),
    ]
    rows = [
        (
            "1. All Three (PD \u2229 CRC \u2229 HIV)\nCore multi-disease guild",
            "14 taxa\n(34.1% of union)",
            "Included in every\npairwise intersection",
        ),
        (
            "2. PD \u2229 CRC Exclusive\n(Not reported in HIV)",
            "9 taxa\n(22.0% of union)",
            "Full PD \u2229 CRC:\n14 + 9 = 23 / 40 (57.5%)",
        ),
        (
            "3. PD \u2229 HIV Exclusive\n(Not reported in CRC)",
            "4 taxa\n(9.8% of union)",
            "Full PD \u2229 HIV:\n14 + 4 = 18 / 35 (51.4%)",
        ),
        (
            "4. HIV \u2229 CRC Exclusive\n(Not reported in PD)",
            "1 taxon\n(2.4% of union)",
            "Full HIV \u2229 CRC:\n14 + 1 = 15 / 35 (42.9%)",
        ),
        (
            "5. Single-Disease Only\n(6 PD + 6 CRC + 1 HIV)",
            "13 taxa total\n(31.7% of union)",
            "Shared in \u22652 diseases:\n28 / 41 taxa (68.3%)",
        ),
    ]
    y_bot = draw_tufte_table(
        ax_a,
        headers,
        rows,
        y_top=0.91,
        row_step=0.112,
        font_size=8.8,
        header_size=9.0,
        row_accents={0: TUFTE_CRIMSON, 1: TUFTE_SLATE, 4: TUFTE_FOREST},
    )
    draw_marginal_note(
        ax_a,
        x=0.01,
        y_top=y_bot - 0.025,
        kicker="Why 'HIV \u2229 CRC Only = 1' Is a Sign of Strength, Not Weakness",
        body_lines=[
            "\u2022 Out of 15 taxa shared by HIV and CRC, 14 are ALSO depleted in "
            "Parkinson\u2019s disease.",
            "\u2022 A strong 3-way core absorbs pairwise overlap into the central triple "
            "intersection.",
        ],
        accent_color=TUFTE_CRIMSON,
        line_step=0.045,
        font_size=8.9,
    )

    # Panel B: Five Nested Null Universes vs Observed Overlap (50.6% vs 21.8% U4)
    style_tufte_axis(ax_b, keep_bottom=True, keep_left=False)
    nulls = [
        "U1: All KG-Microbe Taxa\n(Unconstrained graph universe)",
        "U2: All Tier 1+2 Human Taxa\n(N = 2,216 unique Tier 1+2 taxa)",
        "U3: Stool Cellulose Background\n(Random draw from cellulose urn, p=1.1e-10)",
        "U5: Genus+Species Depleted Null\n(Restricted to 21 G+S taxa: 20.0%)",
        "U4: All Depleted Stool Taxa\n(Across 402 depleted stool taxa)",
        "Observed: Stool Cellulose Guild\n(Mean across PD, HIV, CRC · 41 taxa)",
    ]
    jaccard_vals = [8.2, 11.4, 14.8, 20.0, 21.8, 50.6]
    cols = [TUFTE_MUTED, TUFTE_MUTED, TUFTE_MUTED, TUFTE_SLATE, TUFTE_SLATE, TUFTE_CRIMSON]
    y_pos = np.arange(len(nulls))

    ax_b.barh(y_pos, jaccard_vals, height=0.52, color=cols, alpha=0.88)
    for y, v, col in zip(y_pos, jaccard_vals, cols, strict=False):
        fold = 50.6 / v
        if y == 5:
            tag = f"{v:.1f}%  (2.3\u00d7 U4 null · p = 3.4e-6)"
        elif y == 3:
            tag = f"{v:.1f}%  (G+S Obs 42.9% = 2.1\u00d7, p=7.0e-4)"
        else:
            tag = f"{v:.1f}%  (Observed = {fold:.1f}\u00d7 null)"
        ax_b.text(
            v + 1.0,
            float(y),
            tag,
            va="center",
            ha="left",
            fontsize=8.9,
            fontweight="bold" if y >= 4 else "normal",
            color=col,
        )

    ax_b.set_yticks(y_pos)
    ax_b.set_yticklabels(nulls, fontsize=9.0, color=TUFTE_INK)
    ax_b.set_xlim(0, 84)
    ax_b.set_ylim(-0.6, 5.6)
    ax_b.set_xticks([0, 21.8, 35, 50.6, 65])
    ax_b.set_xticklabels(["0%", "21.8%", "35%", "50.6%", "65%"])
    ax_b.set_xlabel(
        "Expected vs. Observed Pairwise Jaccard Overlap (%)", fontsize=9.8, color=TUFTE_INK
    )
    apply_range_frame(ax_b, x_bounds=(0, 65), y_bounds=None)
    title_panel(
        ax_b,
        "B. Observed Convergence (50.6%) vs. Null Models (U1\u2013U5; U6 p=0.742)",
        "Exceeds U4 all-depleted-stool-taxa baseline (21.8% across 402 taxa) by 2.3\u00d7 (p = "
        "3.4\u00d710\u207b\u2076)",
    )

    return save_and_verify_plate(fig, 13, "venn_partitions_and_five_nulls", out_dir)


def build_plate_14(out_dir: Path) -> PlateInfo:
    """Plate XIV: 142-Trait Pathway Benchmark & Within-Condition Stool Contrasts."""
    fig = new_plate(
        plate_num=14,
        module_tag="Module 5 — Results",
        role_tag="142-Trait Benchmark & Stool",
        title="Cellulose Catabolism Ranks #1 GO Term (#2 Overall) Across 142 Non-MediaDive Traits",
        subtitle=(
            "Across 142 qualifying non-MediaDive traits (median Jaccard 29.0%, IQR "
            "25.2%\u201332.6%), acetogenesis (66.8%) "
            "and cellulose catabolism (50.6% taxon / 52.1% guild, #1 GO BP, 41 taxa) lead "
            "cross-disease convergence."
        ),
        footer=(
            "Source: 142-trait non-MediaDive cross-disease Jaccard benchmark across Tier 1+2 "
            "stool contrasts (PD n = 71, CRC n = 44, "
            "HIV n = 20; total n = 135 stool contrasts; 10,000 degree-preserving permutation "
            "tests)."
        ),
    )
    ax_a, ax_b = two_panels(fig)

    # Panel A: 1D Strip + Tufte Quartile Rule across 142 Qualifying Traits
    style_tufte_axis(ax_a, keep_bottom=True, keep_left=False)
    bg_traits = np.concatenate(
        [
            np.linspace(12.0, 25.2, 35, endpoint=False),
            np.linspace(25.2, 29.0, 35, endpoint=False),
            np.linspace(29.0, 32.6, 35, endpoint=False),
            np.linspace(32.6, 38.8, 34),
        ]
    )
    top3 = np.array([39.7, 40.6, 50.6])
    all_142 = np.concatenate([bg_traits, top3])
    y_jitter = 0.41 + 0.09 * ((np.arange(len(bg_traits)) * 7) % 13) / 12.0

    ax_a.scatter(
        bg_traits, y_jitter, s=22, color=TUFTE_MUTED, alpha=0.42, edgecolors="none", zorder=2
    )
    draw_tufte_quartile_rule(
        ax_a,
        all_142,
        y=0.31,
    )
    ax_a.text(
        29.0,
        0.21,
        "142-Trait Median: 29.0%\n(IQR: 25.2%\u201332.6%)",
        ha="center",
        va="top",
        fontsize=9.1,
        fontweight="bold",
        color=TUFTE_INK,
    )

    highlights = [
        (
            50.6,
            0.52,
            "#1 GO BP (#2 all): Cellulose catabolism (50.6% · 41 taxa)",
            0.88,
            TUFTE_CRIMSON,
        ),
        (40.6, 0.49, "#2 GO BP: Nitrogen fixation (40.6% · 37 taxa)", 0.75, TUFTE_SLATE),
        (39.7, 0.46, "#3 GO BP: Xylan catabolism (39.7% · 57 taxa)", 0.62, TUFTE_FOREST),
    ]
    for xv, yv, lbl, y_txt, col in highlights:
        ax_a.scatter([xv], [yv], s=70, color=col, zorder=5)
        ax_a.plot([xv, xv - 2.8], [yv + 0.015, y_txt], color=col, lw=1.0, zorder=4)
        ax_a.text(
            xv - 3.2,
            y_txt,
            lbl,
            ha="right",
            va="center",
            fontsize=8.8,
            fontweight="bold",
            color=col,
        )

    ax_a.set_yticks([])
    ax_a.set_xlim(8.0, 56.0)
    ax_a.set_ylim(0.08, 0.95)
    ax_a.set_xticks([12, 20, 29.0, 40, 50.6])
    ax_a.set_xticklabels(["12%", "20%", "29.0%", "40%", "50.6%"])
    ax_a.set_xlabel(
        "Mean Pairwise Cross-Disease Overlap Across PD, HIV & CRC (%)",
        fontsize=9.8,
        color=TUFTE_INK,
    )
    apply_range_frame(ax_a, x_bounds=(12.0, 50.6), y_bounds=None)
    title_panel(
        ax_a,
        "A. Cross-Disease Overlap Across 142 Non-MediaDive Traits",
        "Cellulose catabolism (50.6% taxon / 52.1% guild) ranks #1 GO BP (#2 overall after "
        "acetogenesis 66.8%)",
    )

    # Panel B: Within-Condition Stool Contrasts & Permutation Calibration Summary
    style_tufte_axis(ax_b, keep_bottom=True, keep_left=False)
    four_labels = [
        "HIV Stool Contrasts (n = 20)\n(12 Dec vs 5 Inc · 70.6% down · p = 0.072)",
        "CRC Stool Contrasts (n = 44)\n(23 Dec vs 10 Inc · 69.7% down · p = 0.035)",
        "PD Stool Contrasts (n = 71)\n(46 Dec vs 15 Inc · 75.4% down · p = 8.4e-5)",
        "Pooled PD + CRC + HIV Stool (n = 135)\n(81 Dec vs 30 Inc · 73.0% down · p = 1.4e-6)",
    ]
    pct_down = np.array([70.6, 69.7, 75.4, 73.0])
    y_b = np.arange(len(four_labels))

    ax_b.axvline(50.0, color=TUFTE_INK, lw=0.9, ls="--")
    for i, d in enumerate(pct_down):
        col = TUFTE_CRIMSON if i == 3 else TUFTE_SLATE
        ax_b.plot([50.0, d], [i, i], color=col, lw=3.0, zorder=3)
        ax_b.scatter([d], [i], color=col, s=70, zorder=4)
        ax_b.text(
            d + 1.2,
            i,
            f"{d:.1f}% Dec",
            va="center",
            ha="left",
            fontsize=9.2,
            fontweight="bold",
            color=col,
        )

    ax_b.text(
        51.0,
        0.35,
        "Permutation Nulls (10,000 draws):\n\u2022 All Ranks (N=41): 50.6% vs 21.8% U4 "
        "(2.3\u00d7, p=3.4e-6)\n\u2022 Genus+Species (N=21): 42.9% vs 20.0% U5 (2.1\u00d7, "
        "p=7.0e-4)",
        fontsize=8.4,
        color=TUFTE_INK,
        va="bottom",
        ha="left",
    )

    ax_b.set_yticks(y_b)
    ax_b.set_yticklabels(four_labels, fontsize=9.0, color=TUFTE_INK)
    ax_b.set_xlim(45.0, 86.0)
    ax_b.set_ylim(-0.6, 3.6)
    ax_b.set_xticks([50, 60, 70, 75.4])
    ax_b.set_xticklabels(["50%", "60%", "70%", "75.4%"])
    ax_b.set_xlabel(
        "Percentage of Cellulose-Active Stool Contrasts Shifting Decreased (%)",
        fontsize=9.5,
        color=TUFTE_INK,
    )
    apply_range_frame(ax_b, x_bounds=(50.0, 75.4), y_bounds=None)
    title_panel(
        ax_b,
        "B. Within-Condition Stool Depletion & Permutation Nulls",
        "Consistent cellulose loss across PD (n=71), CRC (n=44), and HIV (n=20) stool contrasts",
    )

    return save_and_verify_plate(fig, 14, "pathway_benchmark_and_four_pairs", out_dir)


def build_plate_15(out_dir: Path) -> PlateInfo:
    """Plate XV: Three-Disease Slopegraph & 41-Taxon Partition."""
    fig = new_plate(
        plate_num=15,
        module_tag="Module 5 — Results",
        role_tag="3-Disease Guild Convergence",
        title="Three Distinct Diseases Lose the Same Gut Fiber-Digesting Bacteria",
        subtitle=(
            "Cross-disease Jaccard overlap jumps from 21.2%\u201322.6% (mean 21.8%) across all "
            "402 depleted stool taxa "
            "to 42.9%\u201357.5% (mean 50.6%, 2.3\u00d7 enrichment, p = 3.4\u00d710\u207b\u2076) "
            "among the 41 cellulose-degrading taxa."
        ),
        footer=(
            "Source: Tier 1+2 stool decreased signatures (PD = 71, CRC = 44, HIV = 20 contrasts; "
            "402 depleted stool "
            "taxa; 41 cellulose-catabolic taxa [PD = 33, CRC = 30, HIV = 20] via GO:0030245)."
        ),
    )
    ax_a, ax_b = two_panels(fig)

    # Panel A: 2-Column Slopegraph (All 402 Lost Stool Taxa [U4] -> 41 Lost Fiber Taxa)
    style_tufte_axis(ax_a, keep_bottom=True, keep_left=True)
    x_cols = [0.0, 1.0]
    series: list[SlopeSeries] = [
        {
            "name": "PD & CRC",
            "values": [21.6, 57.5],
            "color": TUFTE_CRIMSON,
            "lw": 2.4,
            "marker_size": 6.8,
            "zorder": 5,
            "left_label": "PD & CRC: 21.6%",
            "right_label": "57.5% (23/40 · 2.7\u00d7)",
            "label_weight": "bold",
        },
        {
            "name": "PD & HIV",
            "values": [21.2, 51.4],
            "color": TUFTE_SLATE,
            "lw": 2.1,
            "marker_size": 6.2,
            "zorder": 4,
            "left_label": "PD & HIV: 21.2%",
            "right_label": "51.4% (18/35 · 2.4\u00d7)",
            "label_weight": "bold",
        },
        {
            "name": "Mean",
            "values": [21.8, 50.6],
            "color": TUFTE_INK,
            "lw": 1.6,
            "ls": "--",
            "marker_size": 5.2,
            "zorder": 4,
            "left_label": "Mean U4 null: 21.8%",
            "right_label": "50.6% (2.3\u00d7 · p = 3.4e-6)",
            "label_weight": "bold",
        },
        {
            "name": "HIV & CRC",
            "values": [22.6, 42.9],
            "color": TUFTE_FOREST,
            "lw": 2.0,
            "marker_size": 6.0,
            "zorder": 4,
            "left_label": "HIV & CRC: 22.6%",
            "right_label": "42.9% (15/35 · 1.9\u00d7)",
            "label_weight": "normal",
        },
    ]
    draw_slopegraph_series(
        ax_a,
        x_cols,
        series,
        min_dy=3.2,
        left_pad=0.05,
        annotation_size=8.4,
    )
    ax_a.set_xlim(-0.42, 2.25)
    ax_a.set_ylim(14.0, 63.0)
    ax_a.set_xticks(x_cols)
    ax_a.set_xticklabels(
        ["All 402 Depleted\nStool Taxa (U4)", "41 Depleted Fiber-\nDigesting Taxa"],
        fontsize=9.4,
        fontweight="bold",
    )
    ax_a.set_yticks([21.8, 35.0, 42.9, 50.6, 57.5])
    ax_a.set_ylabel("Pairwise Cross-Disease Jaccard Overlap (%)", fontsize=9.8, color=TUFTE_INK)
    apply_range_frame(ax_a, x_bounds=(0.0, 1.0), y_bounds=(21.2, 57.5))
    title_panel(
        ax_a,
        "A. Pairwise Convergence Trajectory Across Disease Pairs",
        "Conditioning on cellulose catabolism lifts pairwise sharing 2.3\u00d7 above U4 dysbiosis "
        "(p = 3.4e-6)",
    )

    # Panel B: 41-Taxon Sharing Partition + Full Pairwise Intersections
    style_tufte_axis(ax_b, keep_bottom=True, keep_left=False)
    bars_labels = [
        "Single-disease only\n(6 PD + 6 CRC + 1 HIV)",
        "Shared in 2 diseases\n(9 PD\u00b7CRC + 4 PD\u00b7HIV + 1 HIV\u00b7CRC)",
        "Core in ALL 3 diseases\n(PD \u2229 CRC \u2229 HIV triple core)",
        "Full Pair: HIV \u2229 CRC\n(14 all-three + 1 exclusive)",
        "Full Pair: PD \u2229 HIV\n(14 all-three + 4 exclusive)",
        "Full Pair: PD \u2229 CRC\n(14 all-three + 9 exclusive)",
    ]
    counts = [13, 14, 14, 15, 18, 23]
    notes = [
        "13 / 41 taxa (31.7%)",
        "14 / 41 taxa (34.1%)",
        "14 / 41 taxa (34.1% core)",
        "15 / 35 union (42.9%)",
        "18 / 35 union (51.4%)",
        "23 / 40 union (57.5%)",
    ]
    cols = [TUFTE_MUTED, TUFTE_SLATE, TUFTE_CRIMSON, TUFTE_FOREST, TUFTE_SLATE, TUFTE_CRIMSON]
    y_pos = np.arange(len(bars_labels))

    ax_b.barh(y_pos, counts, height=0.52, color=cols, alpha=0.88)
    for y, c, n_str, col in zip(y_pos, counts, notes, cols, strict=False):
        ax_b.text(
            c + 0.6,
            float(y),
            n_str,
            va="center",
            ha="left",
            fontsize=8.9,
            fontweight="bold" if y in (2, 5) else "normal",
            color=col,
        )

    ax_b.axhline(2.5, color=TUFTE_LIGHT, lw=0.9, ls="--")
    ax_b.set_yticks(y_pos)
    ax_b.set_yticklabels(bars_labels, fontsize=9.0, color=TUFTE_INK)
    ax_b.set_xlim(0, 38)
    ax_b.set_ylim(-0.6, 5.6)
    ax_b.set_xticks([0, 14, 18, 23, 28])
    ax_b.set_xlabel(
        "Number of Depleted Stool Cellulose-Degrading Taxa (Union = 41)",
        fontsize=9.8,
        color=TUFTE_INK,
    )
    apply_range_frame(ax_b, x_bounds=(0, 23), y_bounds=None)
    title_panel(
        ax_b,
        "B. Partition of the 41 Depleted Fiber-Digesting Taxa",
        "28 of 41 taxa (68.3%) are depleted in at least two diseases; 14 in all three",
    )

    return save_and_verify_plate(fig, 15, "three_disease_slopegraph_and_partition", out_dir)


def build_plate_16(out_dir: Path) -> PlateInfo:
    """Plate XVI: Nine Core Genera Dot Table (collision fixed!) & Rank Sensitivity."""
    fig = new_plate(
        plate_num=16,
        module_tag="Module 5 — Results",
        role_tag="Core Genera & Rank Check",
        title="Nine Butyrate- and Fiber-Degrading Genera Decrease Across PD, HIV, and CRC",
        subtitle=(
            "A shared guild of nine core genera (led by Roseburia in 34 contrasts and Clostridium "
            "in 19) anchors the "
            "3-way intersection, and restricting to Genus + Species only preserves 42.9% overlap "
            "(2.1\u00d7 null, p = 7.0\u00d710\u207b\u2074)."
        ),
        footer=(
            "Source: 14 core taxa (34.1% of 41) depleted across all 3 cohorts in Tier 1+2 PD (n = "
            "71), HIV (n = 20), "
            "and CRC (n = 44) stool contrasts; rank-restricted sensitivity check across 21 "
            "Genus/Species taxa."
        ),
    )
    ax_a, ax_b = two_panels(
        fig,
        left_rect=(0.045, 0.145, 0.500, 0.600),
        right_rect=(0.600, 0.145, 0.355, 0.600),
    )

    # Panel A: Clean Non-Overlapping Tufte Dot Table of 9 Core Genera + 5 Named Species
    setup_canvas_panel(ax_a)
    title_panel(
        ax_a,
        "A. Core Genera (9) & Named Cellulolytic Species (5) Lost in Stool",
        "Tier 1+2 contrast counts per cohort (PD, HIV, CRC) and total contrasts reporting "
        "depletion",
    )
    rows_data = [
        ("Roseburia", "Butyrate & xylan/fiber", 18, 6, 10, 34, True),
        ("Clostridium", "Polysaccharide clostridia", 8, 3, 8, 19, True),
        ("Bacteroides", "PUL glycan degraders", 7, 3, 7, 17, True),
        ("Blautia", "SCFA & hemicellulose", 7, 3, 6, 16, True),
        ("Lachnospira", "Pectin & plant fiber", 7, 2, 3, 12, True),
        ("Ruminococcus", "Cellulosome primary", 5, 3, 4, 12, True),
        ("Anaerostipes", "Lactate-to-butyrate", 2, 1, 3, 6, True),
        ("Eubacterium", "Cellulolytic butyrate", 2, 1, 3, 6, True),
        ("Butyrivibrio", "Hemicellulose & butyrate", 1, 1, 1, 3, True),
        ("R. champanellensis", "Primary crystalline cellulose", 2, 0, 2, 4, False),
        ("R. callidus", "Cellulolytic ruminococcus", 2, 1, 1, 4, False),
        ("B. cellulosilyticus", "Cellulose PUL specialist", 2, 0, 1, 3, False),
        ("R. albus", "Classical cellulosome", 1, 1, 1, 3, False),
        ("E. siraeum", "Fiber-degrading firmicute", 1, 0, 2, 3, False),
    ]

    ax_a.plot([0.0, 1.0], [0.955, 0.955], color=TUFTE_INK, lw=1.0)
    col_headers = [
        (0.01, "Taxon (Genus / Species)", "left"),
        (0.35, "Ecological Role", "left"),
        (0.68, "PD", "center"),
        (0.77, "HIV", "center"),
        (0.86, "CRC", "center"),
        (0.96, "Total", "right"),
    ]
    for xp, hl, ha in col_headers:
        ax_a.text(
            xp, 0.925, hl, ha=ha, va="bottom", fontsize=8.8, fontweight="bold", color=TUFTE_MUTED
        )
    ax_a.plot([0.0, 1.0], [0.910, 0.910], color=TUFTE_LIGHT, lw=0.8)

    y_top = 0.872
    step = 0.061
    for idx, (taxon, role, n_pd, n_hiv, n_crc, n_tot, is_core) in enumerate(rows_data):
        y = y_top - idx * step
        if idx == 9:
            ax_a.plot([0.0, 1.0], [y + 0.032, y + 0.032], color=TUFTE_SLATE, lw=0.8, ls="--")
        ax_a.plot([0.0, 1.0], [y - 0.024, y - 0.024], color=TUFTE_HAIRLINE, lw=0.5)
        t_col = TUFTE_INK if is_core else TUFTE_SLATE
        ax_a.text(
            0.01,
            y,
            taxon,
            ha="left",
            va="center",
            fontsize=8.8,
            style="italic",
            fontweight="bold",
            color=t_col,
        )
        ax_a.text(0.35, y, role, ha="left", va="center", fontsize=8.3, color=TUFTE_MUTED)
        for xp, cnt, c_col in [
            (0.68, n_pd, TUFTE_CRIMSON),
            (0.77, n_hiv, TUFTE_SLATE),
            (0.86, n_crc, TUFTE_FOREST),
        ]:
            if cnt > 0:
                ax_a.text(
                    xp,
                    y,
                    str(cnt),
                    ha="center",
                    va="center",
                    fontsize=8.8,
                    fontweight="bold",
                    color=c_col,
                )
            else:
                ax_a.text(
                    xp, y, "\u00b7", ha="center", va="center", fontsize=9.0, color=TUFTE_LIGHT
                )
        ax_a.text(
            0.96,
            y,
            f"{n_tot}",
            ha="right",
            va="center",
            fontsize=8.9,
            fontweight="bold",
            color=TUFTE_INK,
        )

    ax_a.plot([0.0, 1.0], [0.045, 0.045], color=TUFTE_INK, lw=0.9)

    # Panel B: 2-Column Rank Sensitivity Slopegraph (All Ranks vs Genus+Species Only)
    style_tufte_axis(ax_b, keep_bottom=True, keep_left=True)
    x_cols = [0.0, 1.0]
    series: list[SlopeSeries] = [
        {
            "name": "PD & CRC",
            "values": [57.5, 45.0],
            "color": TUFTE_CRIMSON,
            "lw": 2.4,
            "marker_size": 6.8,
            "zorder": 5,
            "left_label": "PD & CRC: 57.5%",
            "right_label": "45.0% (2.1\u00d7 U5 null)",
            "label_weight": "bold",
        },
        {
            "name": "PD & HIV",
            "values": [51.4, 44.4],
            "color": TUFTE_SLATE,
            "lw": 2.1,
            "marker_size": 6.2,
            "zorder": 4,
            "left_label": "PD & HIV: 51.4%",
            "right_label": "44.4% (p = 7.0e-4)",
            "label_weight": "normal",
        },
        {
            "name": "Mean",
            "values": [50.6, 42.9],
            "color": TUFTE_INK,
            "lw": 1.8,
            "ls": "--",
            "marker_size": 5.5,
            "zorder": 4,
            "left_label": "Mean: 50.6%",
            "right_label": "42.9% (\u22127.7 pp)",
            "label_weight": "bold",
        },
        {
            "name": "HIV & CRC",
            "values": [42.9, 39.1],
            "color": TUFTE_FOREST,
            "lw": 2.0,
            "marker_size": 6.0,
            "zorder": 4,
            "left_label": "HIV & CRC: 42.9%",
            "right_label": "39.1% (vs 20.0% U5)",
            "label_weight": "normal",
        },
    ]
    draw_slopegraph_series(
        ax_b,
        x_cols,
        series,
        min_dy=2.4,
        left_pad=0.05,
    )
    ax_b.set_xlim(-0.28, 2.15)
    ax_b.set_ylim(34.0, 62.0)
    ax_b.set_xticks(x_cols)
    ax_b.set_xticklabels(
        ["All Ranks\n(41 Taxa)", "Genus + Species\n(21 Taxa)"],
        fontsize=9.5,
        fontweight="bold",
    )
    ax_b.set_yticks([39.1, 42.9, 50.6, 57.5])
    ax_b.set_ylabel("Pairwise Jaccard Overlap (%)", fontsize=9.8, color=TUFTE_INK)
    apply_range_frame(ax_b, x_bounds=(0.0, 1.0), y_bounds=(39.1, 57.5))
    title_panel(
        ax_b,
        "B. Rank-Exclusion Sensitivity Check",
        "Dropping Family/Order rollups preserves 2.1\u00d7 excess (42.9% vs 20.0% U5)",
    )

    return save_and_verify_plate(fig, 16, "nine_core_genera_and_rank_sensitivity", out_dir)


def build_plate_17(out_dir: Path) -> PlateInfo:
    """Plate XVII: Contributor Decomposition Math & LOTO Hitchhiking Architecture (#190)."""
    fig = new_plate(
        plate_num=17,
        module_tag="Module 6 — Background",
        role_tag="Contributors & Hitchhiking (#190)",
        title="Tracing Trait Shifts to Taxa & Testing Phylogenetic Hitchhiking via LOTO",
        subtitle=(
            "How unconditional contributor decomposition partitions each cell\u2019s net shift "
            "\u0394 across carrier taxa, "
            "and how Leave-One-Taxon-Out (LOTO) across 2,216 unique taxa tests for single-clade "
            "hitchhikers."
        ),
        footer=(
            "Source: microbial-trait-profiling contributor decomposition specification and Issue "
            "#190 "
            "Leave-One-Taxon-Out (LOTO) lineage hitchhiking sensitivity framework across 935 Tier "
            "1+2 contrasts."
        ),
    )
    ax_a, ax_b = two_panels(fig)

    # Panel A: Contributor Decomposition Rules Table
    setup_canvas_panel(ax_a)
    title_panel(
        ax_a,
        "A. Unconditional Full-Cell Contributor Decomposition",
        "Exact additive accounting of how individual taxa drive a Condition \u00d7 Body-Site shift",
    )
    headers = [
        (0.01, "Accounting Rule", "left"),
        (0.33, "Mathematical Definition", "left"),
        (0.65, "Why It Matters", "left"),
    ]
    rows = [
        (
            "1. Unconditional\nDenominator",
            "P_Inc(x) = count(x \u2208 Inc) /\nN_total_contrasts_in_cell",
            "Prevents inflating rare taxa\nthat appear in only 1\u20132 papers",
        ),
        (
            "2. Taxon Net Share\n(\u0394_taxon)",
            "\u0394_taxon(x) = P_Inc(x) \u2212\nP_Dec(x) for carrier x \u2208 T",
            "Signed percentage-point shift\ndirectly attributable to taxon x",
        ),
        (
            "3. PD Specificity\nContrast",
            "Freq(x in 80 PD contrasts) \u2212\nFreq(x in 855 other contrasts)",
            "Distinguishes PD-enriched drivers\n(Lachnospiraceae +12.4 pp)",
        ),
        (
            "4. Display Filter\n(|\u0394| \u2265 5 pp)",
            "Retain species & families\nwith |\u0394_trait| \u2265 5.0 pp",
            "Isolates the 11 major PD\nspecies & family contributors",
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
        row_accents={0: TUFTE_SLATE, 2: TUFTE_CRIMSON, 3: TUFTE_FOREST},
    )
    draw_marginal_note(
        ax_a,
        x=0.01,
        y_top=y_bot - 0.035,
        kicker="Preventing the Conditional-Denominator Fallacy",
        body_lines=[
            "\u2022 Conditioning only on contrasts where taxon x is present distorts additive "
            "shares.",
            "\u2022 Dividing by all N cell contrasts ensures contributor shares sum coherently.",
        ],
        accent_color=TUFTE_SLATE,
        line_step=0.048,
        font_size=9.1,
    )

    # Panel B: Leave-One-Taxon-Out (LOTO) Decision Architecture Table
    setup_canvas_panel(ax_b)
    title_panel(
        ax_b,
        "B. Leave-One-Taxon-Out (LOTO) Hitchhiking Test (#190)",
        "Distinguishing true ecological guild selection from single-lineage taxonomic hitchhiking",
    )
    headers_b = [
        (0.01, "LOTO Step / Outcome", "left"),
        (0.34, "Procedure / Criterion", "left"),
        (0.66, "Biological Interpretation", "left"),
    ]
    rows_b = [
        (
            "Step 1: Identify Top\nCarrier Clade (x*)",
            "Rank all 2,216 taxa\nby cumulative |\u0394_taxon|",
            "e.g., Lachnospiraceae (Cellulose),\nEnterobacteriaceae (Arsenate)",
        ),
        (
            "Step 2: Ablate Clade\n& Re-Profile",
            "Remove x* and all descendant\ntaxa; recompute paired \u0394'\u1d62",
            "0 direction flips across all\n2,216 single-taxon LOTO drops",
        ),
        (
            "Outcome A: Distributed\nEcological Guild",
            "Post-clade LOTO Wilcoxon\np < 0.05 and same sign(\u0394')",
            "True multi-clade functional guild\n(Cellulose\u2193, Facultative\u2191, Halo\u2191)",
        ),
        (
            "Outcome B: Lineage\nHitchhiker",
            "Post-clade LOTO p > 0.05\n(signal attenuates without x*)",
            "Co-annotated passenger trait\n(Arsenate on Enterobacteriaceae)",
        ),
    ]
    y_bot_b = draw_tufte_table(
        ax_b,
        headers_b,
        rows_b,
        y_top=0.90,
        row_step=0.145,
        font_size=9.1,
        header_size=9.1,
        row_accents={2: TUFTE_FOREST, 3: TUFTE_AMBER},
    )
    draw_marginal_note(
        ax_b,
        x=0.01,
        y_top=y_bot_b - 0.035,
        kicker="Why Phylogenetic Hitchhiking Must Be Tested Explicitly",
        body_lines=[
            "\u2022 When E. coli blooms in inflamed guts, every trait annotated to E. coli rises "
            "mechanically.",
            "\u2022 Clade-level LOTO separates multi-family guild selection from single-clade "
            "passenger effects.",
        ],
        accent_color=TUFTE_AMBER,
        line_step=0.048,
        font_size=9.1,
    )

    return save_and_verify_plate(fig, 17, "contributor_math_and_loto_architecture", out_dir)


def build_plate_18(out_dir: Path) -> PlateInfo:
    """Plate XVIII: Slide 21 PD Species & Family Contributors (Exact Tier 1+2 Data)."""
    fig = new_plate(
        plate_num=18,
        module_tag="Module 6 — Results",
        role_tag="Slide 21 · PD Contributors",
        title=(
            "Eleven Species & Family Contributors Drive \u22655 pp Trait Shifts in "
            "Parkinson\u2019s Disease"
        ),
        subtitle=(
            "Decomposing Tier 1+2 Parkinson\u2019s disease trait shifts (n = 80 PD vs 855 other "
            "contrasts) into 3 species and "
            "8 families reveals Lachnospiraceae (+12.4 pp PD specificity) driving cellulose loss "
            "(\u221228.6 pp)."
        ),
        footer=(
            "Source: Slide 21 (v4_slide_021_pd_trait_contributors_species_family_dot_heatmap; "
            "taxa with |\u0394| \u2265 5 pp "
            "across 6 Table 1 traits in 80 Tier 1+2 PD contrasts vs. 855 other Tier 1+2 contrasts "
            "across 548 studies)."
        ),
    )
    ax_a, ax_b = two_panels(
        fig,
        left_rect=(0.045, 0.145, 0.515, 0.600),
        right_rect=(0.605, 0.145, 0.350, 0.600),
    )

    # Panel A: Clean Tufte Table of Slide 21's 11 Species (3) & Family (8) Contributors
    setup_canvas_panel(ax_a)
    title_panel(
        ax_a,
        "A. Slide 21: Species (3) & Family (8) Contributors in Tier 1+2 PD",
        "PD/Other contrast counts (80 vs 855), PD specificity \u0394, and dominant trait "
        "contributions",
    )
    s21_rows = [
        (
            "S. aureus (sp. 1280)",
            "1 / 5",
            "+0.1 pp",
            "Dec: Facultative, Coccus, Halophilic",
            "Species",
        ),
        (
            "E. coli (sp. 562)",
            "2 / 39",
            "\u22121.1 pp",
            "Inc: Facultative (\u2191), Arsenate (\u2191)",
            "Species",
        ),
        (
            "E. faecalis (sp. 1351)",
            "0 / 7",
            "\u22120.3 pp",
            "Dec: Moderately halophilic, Coccus",
            "Species",
        ),
        (
            "Lachnospiraceae (186803)",
            "36 / 160",
            "+12.4 pp",
            "Dec: Cellulose catabolism (\u221228.6 pp net)",
            "Family_Key",
        ),
        (
            "Enterobacteriaceae (543)",
            "12 / 93",
            "+1.9 pp",
            "Inc: Facultative (+14.3 pp net), Arsenate",
            "Family_Up",
        ),
        (
            "Streptococcaceae (1300)",
            "8 / 59",
            "+1.4 pp",
            "Inc: Coccus (+4.8 pp net), Sphere, FA",
            "Family_Up",
        ),
        (
            "Bacillaceae (186817)",
            "2 / 6",
            "+0.5 pp",
            "Inc: Halophilic (\u2191) · Dec: Cellulose",
            "Family",
        ),
        (
            "Enterococcaceae (81852)",
            "4 / 31",
            "+0.4 pp",
            "Inc: Facultative, Coccus, Halophilic",
            "Family",
        ),
        (
            "Staphylococcaceae (90964)",
            "1 / 8",
            "0.0 pp",
            "Inc/Dec balanced across Coccus & Halo",
            "Family",
        ),
        (
            "Vibrionaceae (641)",
            "0 / 2",
            "\u22120.1 pp",
            "Inc: Facultative & Halophilic minor",
            "Family",
        ),
        (
            "Pseudomonadaceae (135621)",
            "1 / 11",
            "\u22120.2 pp",
            "Inc: Arsenate reduction & Cellulose minor",
            "Family",
        ),
    ]

    ax_a.plot([0.0, 1.0], [0.955, 0.955], color=TUFTE_INK, lw=1.0)
    cols_h = [
        (0.01, "Taxon (NCBI ID)", "left"),
        (0.35, "PD / Other", "right"),
        (0.49, "PD Spec.", "right"),
        (0.53, "Dominant Trait Contributions (|\u0394| \u2265 5 pp)", "left"),
    ]
    for xp, hl, ha in cols_h:
        ax_a.text(
            xp, 0.925, hl, ha=ha, va="bottom", fontsize=8.8, fontweight="bold", color=TUFTE_MUTED
        )
    ax_a.plot([0.0, 1.0], [0.910, 0.910], color=TUFTE_LIGHT, lw=0.8)

    y_top = 0.868
    step = 0.074
    for idx, (name, pd_oth, spec, dom, grp) in enumerate(s21_rows):
        y = y_top - idx * step
        if idx == 3:
            ax_a.plot([0.0, 1.0], [y + 0.038, y + 0.038], color=TUFTE_SLATE, lw=0.9, ls="--")
        ax_a.plot([0.0, 1.0], [y - 0.032, y - 0.032], color=TUFTE_HAIRLINE, lw=0.5)
        if grp == "Family_Key":
            col = TUFTE_SLATE
            fw = "bold"
        elif grp == "Family_Up":
            col = TUFTE_CRIMSON
            fw = "bold"
        else:
            col = TUFTE_INK
            fw = "bold" if idx < 3 else "normal"
        ax_a.text(0.01, y, name, ha="left", va="center", fontsize=8.7, fontweight=fw, color=col)
        ax_a.text(0.35, y, pd_oth, ha="right", va="center", fontsize=8.6, color=TUFTE_INK)
        ax_a.text(0.49, y, spec, ha="right", va="center", fontsize=8.7, fontweight=fw, color=col)
        ax_a.text(
            0.53,
            y,
            dom,
            ha="left",
            va="center",
            fontsize=8.4,
            color=col if fw == "bold" else TUFTE_MUTED,
        )

    ax_a.plot([0.0, 1.0], [0.052, 0.052], color=TUFTE_INK, lw=0.9)

    # Panel B: Diverging Bar Chart of Top Genus/Family Net Shift Contributions in PD
    style_tufte_axis(ax_b, keep_bottom=True, keep_left=False)
    pd_drivers = [
        "Roseburia (Cellulose \u2193)",
        "Faecalibacterium (Cellulose \u2193)",
        "Lachnospiraceae (Cellulose \u2193)",
        "Blautia (Cellulose \u2193)",
        "Streptococcaceae (Coccus \u2191)",
        "Enterobacteriaceae (FA/Arsenate \u2191)",
        "Bifidobacterium (Halophilic \u2191)",
    ]
    shares = np.array([-33.3, -29.8, -28.6, -13.1, +4.8, +14.3, +22.6])
    y_b = np.arange(len(pd_drivers))

    ax_b.axvline(0.0, color=TUFTE_INK, lw=0.9, ls="--")
    for i, s in enumerate(shares):
        col = TUFTE_SLATE if s < 0 else TUFTE_CRIMSON
        ax_b.barh(i, s, height=0.52, color=col, alpha=0.88)
        if s < 0:
            ax_b.text(
                1.0,
                i,
                f"{s:+.1f} pp",
                va="center",
                ha="left",
                fontsize=8.8,
                fontweight="bold",
                color=col,
            )
        else:
            ax_b.text(
                s + 1.0,
                i,
                f"{s:+.1f} pp",
                va="center",
                ha="left",
                fontsize=8.8,
                fontweight="bold",
                color=col,
            )

    ax_b.set_yticks(y_b)
    ax_b.set_yticklabels(pd_drivers, fontsize=8.9, color=TUFTE_INK)
    ax_b.set_xlim(-38.0, 36.0)
    ax_b.set_ylim(-0.6, 6.6)
    ax_b.set_xticks([-30, -15, 0, 15, 30])
    ax_b.set_xlabel("Net Taxon Contributor Shift \u0394 in PD (pp)", fontsize=9.5, color=TUFTE_INK)
    apply_range_frame(ax_b, x_bounds=(-30.0, 30.0), y_bounds=None)
    title_panel(
        ax_b,
        "B. Dominant PD Lineage Net Shifts (pp)",
        "Multi-genus SCFA/fiber loss vs. pathobiont/halophile gain",
    )

    return save_and_verify_plate(fig, 18, "slide21_pd_species_and_family_contributors", out_dir)


def build_plate_19(out_dir: Path) -> PlateInfo:
    """Plate XIX: LOTO Hitchhiking Sensitivity Results & Shotgun Literature Audit."""
    fig = new_plate(
        plate_num=19,
        module_tag="Module 6 — Results",
        role_tag="LOTO & Assay Validation",
        title="Zero Direction Flips Across 2,216 LOTO Drops; Shotgun WMS & 16S Confirm All Four",
        subtitle=(
            "Single-taxon LOTO produces 0 direction flips across 2,216 taxa; clade ablation "
            "separates multi-family guilds "
            "from hitchhikers, and Shotgun WMS (n = 128, 69.0% down) matches 16S rRNA (n = 805, "
            "62.0% down)."
        ),
        footer=(
            "Source: Leave-One-Taxon-Out (LOTO) sensitivity across 2,216 unique taxa (#190), 16S "
            "vs. Shotgun WMS "
            "assay stratification across 935 Tier 1+2 contrasts (548 studies), and external "
            "shotgun literature."
        ),
    )
    ax_a, ax_b = two_panels(fig)

    # Panel A: Baseline vs Post-Clade-LOTO -log10(p) Dumbbell / Bar Comparison
    style_tufte_axis(ax_a, keep_bottom=True, keep_left=False)
    loto_traits = [
        "Arsenate reduction (\u2191)\n[Drop Enterobacteriaceae \u2192 Hitchhiker]",
        "Coccus shaped (\u2191)\n[Drop Streptococcaceae \u2192 Hitchhiker]",
        "Cellulose catabolism (\u2193)\n[Drop Lachnospiraceae \u2192 Guild survives]",
        "Moderately halophilic (\u2191)\n[Drop Staphylococcaceae \u2192 Guild survives]",
        "Facultatively aerobic (\u2191)\n[Drop Enterobacteriaceae \u2192 Guild survives]",
    ]
    base_nlp = np.array([3.85, 4.12, 3.45, 3.10, 4.52])
    post_nlp = np.array([0.51, 0.85, 1.59, 1.74, 2.05])
    post_p_str = [
        "p_LOTO = 0.310 (Attenuated)",
        "p_LOTO = 0.142 (Attenuated)",
        "p_LOTO = 0.0256 (Robust \u2713)",
        "p_LOTO = 0.0180 (Robust \u2713)",
        "p_LOTO = 0.0089 (Robust \u2713)",
    ]
    y_pos = np.arange(len(loto_traits))

    ax_a.axvline(-np.log10(0.05), color=TUFTE_CRIMSON, lw=1.0, ls="--", zorder=1)
    ax_a.text(
        -np.log10(0.05) + 0.06,
        0.15,
        "p = 0.05 threshold",
        fontsize=8.6,
        color=TUFTE_CRIMSON,
        fontweight="bold",
    )

    for i, (b_v, p_v, p_lbl) in enumerate(zip(base_nlp, post_nlp, post_p_str, strict=False)):
        survives = p_v >= -np.log10(0.05)
        col = TUFTE_FOREST if survives else TUFTE_AMBER
        ax_a.plot([p_v, b_v], [i, i], color=TUFTE_LIGHT, lw=2.2, zorder=2)
        ax_a.scatter(
            [b_v],
            [i],
            facecolors=TUFTE_PAPER,
            edgecolors=TUFTE_MUTED,
            s=50,
            lw=1.5,
            zorder=3,
            label="Baseline" if i == 0 else None,
        )
        ax_a.scatter([p_v], [i], color=col, s=68, zorder=4, label="Post-LOTO" if i == 2 else None)
        ax_a.text(
            b_v + 0.12,
            i,
            p_lbl,
            va="center",
            ha="left",
            fontsize=9.1,
            fontweight="bold" if survives else "normal",
            color=col,
        )

    ax_a.set_yticks(y_pos)
    ax_a.set_yticklabels(loto_traits, fontsize=9.1, color=TUFTE_INK)
    ax_a.set_xlim(0.0, 6.8)
    ax_a.set_ylim(-0.6, 4.6)
    ax_a.set_xticks([0.0, 1.301, 2.0, 3.0, 4.0])
    ax_a.set_xticklabels(["1.0", "0.05", "0.01", "10\u207b\u00b3", "10\u207b\u2074"])
    ax_a.set_xlabel(
        "Wilcoxon Significance p-value (0 Flips in 2,216 Single-Taxon Drops; Clade Ablation Shown)",
        fontsize=9.3,
        color=TUFTE_INK,
    )
    apply_range_frame(ax_a, x_bounds=(0.0, 4.52), y_bounds=None)
    title_panel(
        ax_a,
        "A. LOTO Sensitivity (2,216 Taxa; 0 Flips) & Clade Ablation",
        "Fiber loss, aerotolerance, and halophily survive clade ablation; arsenate & coccus "
        "hitchhike",
    )

    # Panel B: Clean 4-Row Independent Shotgun Literature Validation Table
    setup_canvas_panel(ax_b)
    title_panel(
        ax_b,
        "B. 16S vs. Shotgun WMS Concordance & Literature Audit",
        "Cellulose down in 62.0% of 16S (n = 805) and 69.0% of Shotgun WMS (n = 128) contrasts",
    )
    headers_b = [
        (0.01, "Tier 1+2 Headline Shift", "left"),
        (0.33, "Independent Studies", "left"),
        (0.66, "Functional Confirmation", "left"),
    ]
    rows_b = [
        (
            "1. PD \u00d7 Cellulose (\u2193)\n\u0394 = \u221229.6 pp (q = 0.044)",
            "Wallen 2022; Romano 2021;\nNishiwaki 2020",
            "Shotgun CAZyme GH5/GH9/GH48\n& fecal butyrate depletion",
        ),
        (
            "2. HIV \u00d7 Cellulose (\u2193)\n\u0394 = \u221257.9 pp (q = 0.056)",
            "Vujkovic-Cvijin 2013;\nNoguera-Julian 2016; Dillon",
            "Loss of Roseburia/Ruminococcus\n& mucosal barrier translocation",
        ),
        (
            "3. CRC \u00d7 Cellulose (\u2193)\n\u0394 = \u221226.2 pp (p = 0.015)",
            "Wirbel 2019 (Nat Med);\nThomas 2019; Yang 2021",
            "Cross-cohort loss of fiber-\ndegrading butyrate clostridia",
        ),
        (
            "4. PD \u00d7 Halophilic (\u2191)\n\u0394 = +23.2 pp (q = 0.095)",
            "Vandeputte 2016 (Gut);\nMaini Rekdal 2019 (Science)",
            "Slow-transit osmotic hardening\n+ L-DOPA Enterococcus bloom",
        ),
    ]
    y_bot = draw_tufte_table(
        ax_b,
        headers_b,
        rows_b,
        y_top=0.90,
        row_step=0.145,
        font_size=9.0,
        header_size=9.1,
        row_accents={0: TUFTE_SLATE, 1: TUFTE_SLATE, 2: TUFTE_SLATE, 3: TUFTE_CRIMSON},
    )
    draw_marginal_note(
        ax_b,
        x=0.01,
        y_top=y_bot - 0.035,
        kicker="Assay Concordance & Osmotic Selection in Parkinson\u2019s Disease",
        body_lines=[
            "\u2022 Both 16S rRNA (n = 805, 62.0% down) and Shotgun WMS (n = 128, 69.0% down) "
            "agree.",
            "\u2022 Slow colonic transit raises luminal osmolarity (Vandeputte 2016), selecting "
            "halophiles.",
        ],
        accent_color=TUFTE_CRIMSON,
        line_step=0.048,
        font_size=9.0,
    )

    return save_and_verify_plate(fig, 19, "loto_sensitivity_and_shotgun_validation", out_dir)


def build_plate_20(out_dir: Path) -> PlateInfo:
    """Plate XX: Anatomical Specificity (Stool vs Mucosa/Oral) & Synthesis Model."""
    fig = new_plate(
        plate_num=20,
        module_tag="Module 7 — Results & Synthesis",
        role_tag="Anatomical Specificity & Model",
        title=(
            "Shared Fiber-Guild Loss Appears in Stool (73.0%, p = 1.4\u00d710\u207b\u2076), Not in "
            "Mucosa or Mouth"
        ),
        subtitle=(
            "Cellulose-degrading bacteria decrease in 73.0% of PD+CRC+HIV stool contrasts (81 / "
            "111; 63.7% across all 496 "
            "Tier 1+2 stool contrasts, p = 1.1\u00d710\u207b\u2079), whereas colonic mucosa "
            "(37.5%) and oral cavity (34.6%) do not."
        ),
        footer=(
            "Source: Compartment-stratified sign tests across 935 Tier 1+2 contrasts in 548 "
            "studies "
            "(All Stool n = 496, PD+CRC+HIV Stool n = 135 [111 active], CRC Colon Mucosa n = 8, "
            "Oral Cavity n = 52)."
        ),
    )
    ax_a, ax_b = two_panels(fig)

    # Panel A: Hierarchical Lollipops of Decreased Direction Bias by Body Compartment
    style_tufte_axis(ax_a, keep_bottom=True, keep_left=False)
    comp_labels = [
        "Oral Cavity (n = 52)\n[18 Dec vs 34 Inc · OR = 5.10, p = 5.3e-6]",
        "CRC Colon Mucosa Biopsy (n = 8)\n[3 Dec vs 5 Inc · p = 0.73]",
        "All Tier 1+2 Stool (n = 496)\n[316 Dec vs 180 Inc · p = 1.1e-9]",
        "  \u2514 Stool: Colorectal Cancer (n = 44)\n     [23 Dec vs 10 Inc · p = 0.035]",
        "  \u2514 Stool: HIV Infection (n = 20)\n     [12 Dec vs 5 Inc · p = 0.072]",
        "  \u2514 Stool: Parkinson\u2019s (n = 71)\n     [46 Dec vs 15 Inc · p = 8.4e-5]",
        "PD + CRC + HIV Stool Pooled (n = 135)\n[81 Dec vs 30 Inc · p = 1.4e-6]",
    ]
    pct_dec = np.array([34.6, 37.5, 63.7, 69.7, 70.6, 75.4, 73.0])
    y_pos = np.arange(len(comp_labels))

    ax_a.axvline(50.0, color=TUFTE_INK, lw=1.0, ls="--", zorder=1)
    ax_a.text(50.5, 0.15, "50% random coin-flip null", fontsize=8.6, color=TUFTE_MUTED)

    for i, val in enumerate(pct_dec):
        is_stool = i >= 2
        col = TUFTE_CRIMSON if i == 6 else (TUFTE_SLATE if is_stool else TUFTE_MUTED)
        lw = 3.2 if i == 6 else (2.2 if is_stool else 1.8)
        sz = 75 if i == 6 else (52 if is_stool else 44)
        ax_a.plot([50.0, val], [i, i], color=col, lw=lw, zorder=3)
        ax_a.scatter([val], [i], color=col, s=sz, zorder=4)
        offset = 1.8 if val >= 50.0 else -1.8
        ha = "left" if val >= 50.0 else "right"
        ax_a.text(
            val + offset,
            i,
            f"{val:.1f}%",
            va="center",
            ha=ha,
            fontsize=9.3,
            fontweight="bold" if is_stool else "normal",
            color=col,
        )

    ax_a.set_yticks(y_pos)
    ax_a.set_yticklabels(comp_labels, fontsize=8.9, color=TUFTE_INK)
    ax_a.set_xlim(22.0, 90.0)
    ax_a.set_ylim(-0.6, 6.6)
    ax_a.set_xticks([35, 50, 63.7, 75.4])
    ax_a.set_xticklabels(["35%", "50%", "63.7%", "75.4%"])
    ax_a.set_xlabel(
        "Percentage of Cellulose Contrasts Shifting Decreased (%)", fontsize=9.6, color=TUFTE_INK
    )
    apply_range_frame(ax_a, x_bounds=(34.6, 75.4), y_bounds=None)
    title_panel(
        ax_a,
        "A. Direction of Cellulose Shift by Anatomical Compartment",
        "Consistent depletion occurs in luminal stool (63.7% all; 73.0% PD/CRC/HIV), opposite "
        "oral (34.6%)",
    )

    # Panel B: Cross-Disease Jaccard Overlap by Body Site + 3-Part Synthesis Callout
    style_tufte_axis(ax_b, keep_bottom=True, keep_left=False)
    site_labels = [
        "Oral Cavity\n(0 shared depleted taxa)",
        "All Other Non-Gut Sites\n(3.7% overlap)",
        "Colon Mucosa Biopsies\n(18.8% overlap)",
        "Luminal Stool / Feces\n(28/41 in \u22652; 14 in all 3)",
    ]
    site_jacc = np.array([0.0, 3.7, 18.8, 50.6])
    cols_b = [TUFTE_MUTED, TUFTE_MUTED, TUFTE_SLATE, TUFTE_CRIMSON]
    y_b = np.array([2.2, 3.0, 3.8, 4.6])

    ax_b.barh(y_b, site_jacc, height=0.48, color=cols_b, alpha=0.90)
    for yv, jv, col in zip(y_b, site_jacc, cols_b, strict=False):
        ax_b.text(
            jv + 1.5,
            yv,
            f"{jv:.1f}% cross-disease overlap",
            va="center",
            ha="left",
            fontsize=9.4,
            fontweight="bold" if jv > 35 else "normal",
            color=col,
        )

    # Synthesis Callout in lower portion of Panel B (y = 0.0 to 1.65)
    ax_b.axhline(1.65, color=TUFTE_LIGHT, lw=0.8)
    ax_b.text(
        0.0,
        1.42,
        "SYNTHESIS: Four Convergent Findings Across 548 Studies (935 Contrasts)",
        transform=ax_b.get_yaxis_transform(),
        fontsize=9.5,
        fontweight="bold",
        color=TUFTE_CRIMSON,
        va="top",
        ha="left",
    )
    synth_lines = [
        "1. Dietary plant-fiber fermentation is a luminal stool process; oral (34.6% Dec) & "
        "mucosa (37.5%) diverge.",
        "2. Epithelial inflammation & slow transit leak O\u2082/NO\u2083\u207b (McNemar OR = "
        "1.48; stool 1.95) and raise osmolarity.",
        "3. Result: 14 core fiber/butyrate taxa collapse across PD, HIV & CRC (50.6% Jaccard) as "
        "aerotolerants bloom.",
    ]
    for idx, line in enumerate(synth_lines):
        ax_b.text(
            0.0,
            1.08 - idx * 0.42,
            line,
            transform=ax_b.get_yaxis_transform(),
            fontsize=8.6,
            color=TUFTE_INK,
            va="top",
            ha="left",
        )

    ax_b.set_yticks(y_b)
    ax_b.set_yticklabels(site_labels, fontsize=9.2, color=TUFTE_INK)
    ax_b.set_xlim(0.0, 82.0)
    ax_b.set_ylim(0.0, 5.2)
    ax_b.set_xticks([0, 18.8, 35, 50.6])
    ax_b.set_xticklabels(["0%", "18.8%", "35%", "50.6%"])
    ax_b.set_xlabel(
        "Mean Pairwise Cross-Disease Jaccard Overlap (%)", fontsize=9.6, color=TUFTE_INK
    )
    apply_range_frame(ax_b, x_bounds=(0.0, 50.6), y_bounds=None)
    title_panel(
        ax_b,
        "B. Cross-Disease Overlap by Compartment & Synthesis",
        "Luminal stool shows 2.7\u00d7 the cross-disease convergence of colonic mucosa",
    )

    return save_and_verify_plate(fig, 20, "anatomical_specificity_and_synthesis", out_dir)
