"""Audit trait correlation, density adjustment, near-twin collinearity,
and ontology redundancy (#264, #352).

Addresses Issue #264 and Issue #352:
- Computes baseline pairwise correlation distribution and PC1 variance across non-MediaDive traits.
- Evaluates PC1 regression adjustment and profile composition normalization (Task 1).
- Audits METPO ontology redundancy, states branch collapsing rules, and bounds
  distinct trait space (Task 2).
- Tests collapse of facultatively aerobic x moderately halophilic vs.
  survival of convergent six.
- Tests conditional contrast and VIF for coccus shaped vs. sphere shaped (Task 3).
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path
from typing import Any

from mtp.density_adjustment import (
    DEFAULT_COLUMNS_PATH,
    DEFAULT_OUTPUT_DIR,
    DEFAULT_PROFILES_PATH,
    SIGNIFICANT_POST_ADJUSTMENT_STATUSES,
    CorrelationMetrics,
    add_trait_list_argument,
    format_correlation_metric,
    format_metric_diff,
    format_p_and_q,
    format_paired_regression_slope,
    format_paired_regression_test,
    format_profile_regression_slope,
    format_q_value,
    format_report_number,
    nominal_test_count_statement,
    q_value_statement,
    reported_traits_heading,
    run_density_audit,
)
from mtp.ontology_audit import run_ontology_redundancy_audit
from mtp.trait_lists import DeclaredTraitList, load_trait_list, output_dir_argument


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Audit trait correlation, annotation density adjustment, conditional near-twins, "
            "and METPO ontology redundancy."
        )
    )
    parser.add_argument(
        "--profiles",
        type=Path,
        default=DEFAULT_PROFILES_PATH,
        help=f"Path to study_trait_profiles_wide.tsv (default: {DEFAULT_PROFILES_PATH})",
    )
    parser.add_argument(
        "--columns",
        type=Path,
        default=DEFAULT_COLUMNS_PATH,
        help=f"Path to taxa_trait_matrix_columns.tsv (default: {DEFAULT_COLUMNS_PATH})",
    )
    parser.add_argument(
        "--output-dir",
        type=output_dir_argument,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Directory to write summary TSVs and reports (default: {DEFAULT_OUTPUT_DIR})",
    )
    parser.add_argument(
        "--code-commit",
        type=str,
        default="",
        help="Code commit hash for provenance tracking",
    )
    add_trait_list_argument(parser)
    return parser.parse_args()


def _q_text(q: float) -> str:
    """A BH q-value for the conclusions; NA when the trait is untestable or not reported.

    Three significant digits, or more when a q next to 0.05 would otherwise round onto
    the other side of the threshold (#1325).
    """
    return f"q={format_q_value(q)}"


def _reported_family_size(t1_effects: list[Any]) -> int:
    """The q-value family, which must be the reported rows: never empty or mismatched (#1324)."""
    if not t1_effects:
        raise ValueError("The density results report no traits, so there is no q-value family")
    sizes = {eff.bh_family_size for eff in t1_effects}
    if sizes != {len(t1_effects)}:
        raise ValueError(
            f"The q-value family sizes {sorted(sizes)} differ from the "
            f"{len(t1_effects)} reported traits"
        )
    return len(t1_effects)


def _conclusion_number(
    value: float, format_spec: str, *, allow_positive_infinity: bool = False
) -> str:
    """Render section 6 summaries without treating undefined values as zero (#1328)."""
    if math.isfinite(value):
        return format(value, format_spec)
    if allow_positive_infinity and value == math.inf:
        return "inf"
    return "NA"


def _conclusion_metrics(metrics: CorrelationMetrics) -> tuple[str, str]:
    """Return median r and PC1 percent; the producer's no-variation zeros are unavailable."""
    return (
        format_correlation_metric(metrics, "median_pairwise_r", "+.4f"),
        format_correlation_metric(metrics, "pc1_variance_pct", ".2f"),
    )


def _plural(count: int, singular: str, plural: str) -> str:
    """Return singular when count == 1, otherwise plural (#1407)."""
    return singular if count == 1 else plural


def write_consolidated_report(
    output_dir: Path,
    density_res: dict[str, Any],
    ontology_res: dict[str, Any],
    trait_list: DeclaredTraitList | None = None,
) -> Path:
    """Generate a comprehensive consolidated report for Issue #264 covering all tasks.

    With a declared ``trait_list`` the conclusions of section 6 are not written: they
    name Table 1 traits (#1054). They quote BH q-values and statuses from the reported
    effects, and read NA for a trait that is not reported (#1321).
    Section 6 describes the computed correlations and VIFs without inferring their
    cause, a severity category, or joint-model coefficient stability (#1328).
    """
    t1_effects = density_res["table1_effects"]
    family_size = _reported_family_size(t1_effects)
    raw = density_res["raw_pearson"]
    raw_sp = density_res["raw_spearman"]
    adj_pc1 = density_res["adj_pc1_pearson"]
    adj_norm = density_res.get("adj_norm_pearson", density_res.get("adj_dw_pearson"))
    twin = density_res["twin_summary"]

    n_pairs = ontology_res["n_redundant_pairs"]
    n_clusters = ontology_res["n_clusters"]
    n_collapsed = ontology_res["n_collapsed_features"]
    bounded_space = ontology_res["bounded_trait_space"]
    nat_m = ontology_res["native_metrics"]
    max_m = ontology_res["max_metrics"]
    adj_col_m = ontology_res["adj_collapsed_metrics"]

    t1_by_id = {eff.trait_id: eff for eff in t1_effects}
    aerobic_eff = t1_by_id.get("METPO:1000608")
    aerobic_label = aerobic_eff.trait_label if aerobic_eff is not None else "facultatively aerobic"
    nan = float("nan")
    aerobic_status = aerobic_eff.status if aerobic_eff is not None else "not reported"
    aerobic_q_raw_str = _q_text(aerobic_eff.raw_wilcoxon_bh_q if aerobic_eff is not None else nan)
    aerobic_q_adj_str = _q_text(aerobic_eff.adj_wilcoxon_bh_q if aerobic_eff is not None else nan)
    survivors = ", ".join(
        eff.trait_label for eff in t1_effects if eff.status in SIGNIFICANT_POST_ADJUSTMENT_STATUSES
    )
    report_out = output_dir / "density_adjustment_report.txt"
    lines = [
        "=" * 80,
        "TRAIT CORRELATION, ANNOTATION DENSITY ADJUSTMENT, AND REDUNDANCY AUDIT",
        "=" * 80,
        *(
            []
            if trait_list is None
            else [
                f"Trait list: {trait_list.source_name} (sha256 {trait_list.sha256}, "
                f"{len(trait_list.trait_ids)} traits)"
            ]
        ),
        "",
        "1. BASELINE DIMENSIONALITY AND CORRELATION STRUCTURE",
        "-" * 80,
        f"Input profiles evaluated: {raw.n_samples:,}",
        f"Non-MediaDive trait columns: {raw.n_traits}",
        "Pearson pairwise correlation: "
        f"median {format_correlation_metric(raw, 'median_pairwise_r', '+.4f')}, "
        f"mean {format_correlation_metric(raw, 'mean_pairwise_r', '+.4f')}",
        "Pairs exceeding +0.5: "
        f"{format_correlation_metric(raw, 'pct_pairs_gt_05', '.1f', suffix='%')} | "
        "Pairs exceeding +0.9: "
        f"{format_correlation_metric(raw, 'pct_pairs_gt_09', '.1f', suffix='%')}",
        "PC1 variance explained: "
        f"{format_correlation_metric(raw, 'pc1_variance_pct', '.2f', suffix='%')} | "
        "PC1-PC5 variance explained: "
        f"{format_correlation_metric(raw, 'pc1_5_variance_pct', '.2f', suffix='%')}",
        "Effective number of independent tests: "
        f"Li & Ji (2005) = {raw.m_eff_li_ji:.1f}, "
        f"Galwey (2009) = {raw.m_eff_galwey:.1f}, "
        f"Cheverud/Nyholt = {raw.m_eff_nyholt:.1f}",
        "",
        "Spearman baseline verification:",
        "  Median pairwise rho: "
        f"{format_correlation_metric(raw_sp, 'median_pairwise_r', '+.4f')}, "
        "Mean pairwise rho: "
        f"{format_correlation_metric(raw_sp, 'mean_pairwise_r', '+.4f')}",
        "  Pairs > 0.5: "
        f"{format_correlation_metric(raw_sp, 'pct_pairs_gt_05', '.1f', suffix='%')} | "
        "Pairs > 0.9: "
        f"{format_correlation_metric(raw_sp, 'pct_pairs_gt_09', '.1f', suffix='%')}",
        "  PC1 variance: "
        f"{format_correlation_metric(raw_sp, 'pc1_variance_pct', '.2f', suffix='%')}",
        "",
        "2. DENSITY AXIS ADJUSTMENT AND CORRELATION CENTERING",
        "-" * 80,
        f"Test Case Pair: {aerobic_label} x moderately halophilic",
        f"  Baseline Pearson r: {format_correlation_metric(raw, 'test_pair_r', '+.4f')}",
        "  Post-PC1 regression r: "
        f"{format_correlation_metric(adj_pc1, 'test_pair_r', '+.4f')} "
        f"(change: {format_metric_diff(adj_pc1, raw)})",
    ]
    if adj_norm is not None:
        lines.append(
            f"  Post-normalization r: {format_correlation_metric(adj_norm, 'test_pair_r', '+.4f')}"
        )

    lines.extend(
        [
            "",
            "Global Correlation Matrix Centering (Method A - PC1 Regression):",
            "  Median pairwise r: "
            f"{format_correlation_metric(adj_pc1, 'median_pairwise_r', '+.4f')}",
            f"  Mean pairwise r: {format_correlation_metric(adj_pc1, 'mean_pairwise_r', '+.4f')}",
            "  Pairs > 0.5: "
            f"{format_correlation_metric(adj_pc1, 'pct_pairs_gt_05', '.1f', suffix='%')}",
            "",
            (
                f"3. {reported_traits_heading(trait_list)} "
                f"(PAIRED EXPERIMENTS N={twin.n_paired_experiments:,})"
            ),
            "-" * 80,
            *q_value_statement(family_size),
            nominal_test_count_statement(t1_effects),
        ]
    )

    for eff in t1_effects:
        lines.append(
            f"  {eff.trait_label:<30s} | "
            f"Raw mean: {eff.raw_mean_diff:+.3f} "
            f"({format_p_and_q(eff.raw_wilcoxon_p, eff.raw_wilcoxon_bh_q)}) | "
            f"Adj mean: {eff.adj_mean_diff:+.3f} "
            f"({format_p_and_q(eff.adj_wilcoxon_p, eff.adj_wilcoxon_bh_q)}) | "
            f"{eff.status}"
        )

    lines.extend(
        [
            "",
            "4. NEAR-TWIN CONDITIONAL TEST AND COLLINEARITY (COCCUS VS. SPHERE)",
            "-" * 80,
            f"Candidate pair: {twin.trait_a_label} ({twin.trait_a}) vs. "
            f"{twin.trait_b_label} ({twin.trait_b})",
            "Profile-level correlation: r = "
            f"{format_report_number(twin.profile_correlation, '+.4f')} | "
            "Variance Inflation Factor (VIF): "
            f"{format_report_number(twin.vif_profile, '.2f', allow_positive_infinity=True)}",
            "Profile OLS regression slope (coccus ~ sphere): "
            f"{format_profile_regression_slope(twin)}",
            f"Paired contrasts (N={twin.n_paired_experiments}): "
            f"r = {format_report_number(twin.paired_correlation, '+.4f')} | Paired VIF: "
            f"{format_report_number(twin.vif_paired, '.2f', allow_positive_infinity=True)}",
            f"Paired regression slope: "
            f"{format_paired_regression_slope(twin.paired_slope_a_on_b)} "
            f"({format_paired_regression_test(twin.paired_slope_pvalue)})",
            "",
            "5. METPO ONTOLOGY REDUNDANCY AUDIT AND BRANCH COLLAPSING (TASK 2)",
            "-" * 80,
            f"Audited non-MediaDive traits: {raw.n_traits}",
            f"Redundant pairs identified: {n_pairs} {_plural(n_pairs, 'pair', 'pairs')} "
            f"across {n_clusters} concept {_plural(n_clusters, 'cluster', 'clusters')}",
            f"Redundant features collapsed: {n_collapsed}",
            f"Bounded distinct trait space: {bounded_space} "
            f"{_plural(bounded_space, 'feature', 'features')}",
            "",
            "Stated Collapsing Rule:",
            "  1. Survivor Selection (Max-Support Evidence Rule): The term with the greatest",
            "     total annotation evidence / profile counts survives (e.g. coccus over sphere;",
            "     motile over flagellated/gliding/axially filamented).",
            "  2. Count Semantics:",
            "     - Max Member Count (lower bound on the union): C_i = max_{m} C_im per profile.",
            "       Avoids artificial inflation from double-counting co-annotated taxa.",
            "     - Native Surviving Representative Count: C_i = C_{i, survivor}.",
            "",
            "Dimensionality Comparison (Raw vs Collapsed vs PC1 Adjusted):",
            f"  - Raw Uncollapsed (P={raw.n_traits}):       "
            f"PC1 var = {format_correlation_metric(raw, 'pc1_variance_pct', '.2f', suffix='%')} | "
            f"Med r = {format_correlation_metric(raw, 'median_pairwise_r', '+.4f')} | "
            f"M_eff (Galwey) = {raw.m_eff_galwey:.1f}",
            f"  - Collapsed Native (P={bounded_space}):      "
            "PC1 var = "
            f"{format_correlation_metric(nat_m, 'pc1_variance_pct', '.2f', suffix='%')} | "
            f"Med r = {format_correlation_metric(nat_m, 'median_pairwise_r', '+.4f')} | "
            f"M_eff (Galwey) = {nat_m.m_eff_galwey:.1f}",
            f"  - Collapsed Deduplicated (P={bounded_space}): "
            "PC1 var = "
            f"{format_correlation_metric(max_m, 'pc1_variance_pct', '.2f', suffix='%')} | "
            f"Med r = {format_correlation_metric(max_m, 'median_pairwise_r', '+.4f')} | "
            f"M_eff (Galwey) = {max_m.m_eff_galwey:.1f}",
            f"  - Collapsed + PC1 Adj (P={bounded_space}):   "
            "PC1 var = "
            f"{format_correlation_metric(adj_col_m, 'pc1_variance_pct', '.2f', suffix='%')} | "
            f"Med r = {format_correlation_metric(adj_col_m, 'median_pairwise_r', '+.4f')} | "
            f"M_eff (Galwey) = {adj_col_m.m_eff_galwey:.1f}",
            "",
            "Redundancy vs. PC1 Summary:",
            f"  Branch collapsing bounds distinct trait space to {bounded_space} "
            f"{_plural(bounded_space, 'feature', 'features')}",
            "  (raw PC1 variance explained = "
            f"{format_correlation_metric(raw, 'pc1_variance_pct', '.2f', suffix='%')}; "
            "collapsed deduplicated PC1 variance explained = "
            f"{format_correlation_metric(max_m, 'pc1_variance_pct', '.2f', suffix='%')};",
            "  collapsed + PC1 regression PC1 variance explained = "
            f"{format_correlation_metric(adj_col_m, 'pc1_variance_pct', '.2f', suffix='%')}).",
            "  These descriptive comparisons summarize the annotation matrices without",
            "  attributing the correlation axis to curation effort or ontology redundancy.",
            "",
            "6. CONCLUSION & METHODOLOGICAL IMPLICATIONS",
            "-" * 80,
        ]
    )
    if trait_list is not None:
        lines.extend(
            [
                "Not written for a declared trait list; section 3 reports each listed",
                "trait's computed effects.",
                "=" * 80,
            ]
        )
        report_out.write_text("\n".join(lines), encoding="utf-8")
        return report_out

    raw_r, raw_pc1 = _conclusion_metrics(raw)
    adjusted_r, _ = _conclusion_metrics(adj_pc1)
    native_r, native_pc1 = _conclusion_metrics(nat_m)
    max_member_r, max_member_pc1 = _conclusion_metrics(max_m)
    collapsed_adjusted_r, collapsed_adjusted_pc1 = _conclusion_metrics(adj_col_m)
    profile_r = _conclusion_number(twin.profile_correlation, "+.4f")
    paired_r = _conclusion_number(twin.paired_correlation, "+.4f")
    profile_vif = _conclusion_number(twin.vif_profile, ".2f", allow_positive_infinity=True)
    paired_vif = _conclusion_number(twin.vif_paired, ".2f", allow_positive_infinity=True)
    n_orphan = int(ontology_res.get("n_orphan_features", 0))
    n_collapsed_in_clusters = n_collapsed - n_orphan
    cluster_summary = (
        f"{n_collapsed_in_clusters} "
        f"{_plural(n_collapsed_in_clusters, 'feature', 'features')} "
        f"collapsed into {n_clusters} {_plural(n_clusters, 'cluster', 'clusters')}"
    )
    orphan_summary = (
        f"{n_orphan} {_plural(n_orphan, 'feature', 'features')} "
        f"removed because "
        f"{_plural(n_orphan, 'its cluster survivor is', 'their cluster survivors are')} absent"
    )
    if n_orphan > 0 and n_collapsed_in_clusters > 0:
        collapsed_summary = f"{cluster_summary}; {orphan_summary}"
    elif n_orphan > 0:
        collapsed_summary = orphan_summary
    else:
        collapsed_summary = cluster_summary
    lines.extend(
        [
            "1. Median pairwise Pearson r across profile trait annotations: "
            f"raw = {raw_r}; post-PC1 regression = {adjusted_r}.",
            f"   Raw PC1 variance explained (%) = {raw_pc1}.",
            f"2. After PC1 regression the {aerobic_label} paired contrast is labelled "
            f"'{aerobic_status}'",
            f"   (raw BH {aerobic_q_raw_str}, adjusted BH {aerobic_q_adj_str}).",
            "3. Significant before and after PC1 regression in the same direction "
            f"(BH q <= 0.05; survived or enhanced): {survivors or 'none'}.",
            f"4. {twin.trait_a_label} and {twin.trait_b_label}: "
            f"profile r = {profile_r}, VIF = {profile_vif}.",
            "   Paired increased-minus-decreased experiment contrasts "
            f"(N={twin.n_paired_experiments:,}): r = {paired_r}, VIF = {paired_vif}.",
            "   Joint-model coefficient stability was not evaluated.",
            f"5. Ontology redundancy accounts for {n_pairs} "
            f"redundant {_plural(n_pairs, 'pair', 'pairs')} ({collapsed_summary}),",
            f"   bounding distinct trait space to {bounded_space} "
            f"{_plural(bounded_space, 'feature', 'features')}.",
            f"   PC1 variance explained (%): raw = {raw_pc1}; collapsed native = {native_pc1};",
            f"   collapsed max-member = {max_member_pc1}; "
            f"collapsed + PC1 regression = {collapsed_adjusted_pc1}.",
            f"   Median pairwise Pearson r: collapsed native = {native_r}; "
            f"collapsed max-member = {max_member_r};",
            f"   collapsed + PC1 regression = {collapsed_adjusted_r}.",
            "These are descriptive comparisons of annotation matrices; "
            "they do not identify the source of the observed associations.",
            "NA denotes an unavailable summary; inf denotes a reported positive-infinite VIF.",
            "=" * 80,
        ]
    )

    report_text = "\n".join(lines)
    report_out.write_text(report_text, encoding="utf-8")
    return report_out


def main() -> int:
    args = parse_args()
    try:
        trait_list = None if args.trait_list is None else load_trait_list(args.trait_list)
        density_res = run_density_audit(
            profiles_path=args.profiles,
            columns_path=args.columns,
            output_dir=args.output_dir,
            code_commit=args.code_commit,
            trait_list=trait_list,
        )
        ontology_res = run_ontology_redundancy_audit(
            profiles_path=args.profiles,
            columns_path=args.columns,
            output_dir=args.output_dir,
        )
        report_path = write_consolidated_report(
            args.output_dir, density_res, ontology_res, trait_list
        )
    except Exception as exc:
        sys.stderr.write(f"Error during density & redundancy audit: {exc}\n")
        return 1

    print("\n" + "=" * 80)
    print("DENSITY ADJUSTMENT & ONTOLOGY REDUNDANCY AUDIT COMPLETED")
    print("=" * 80)
    print(f"Summary metrics written to:          {density_res['summary_tsv']}")
    print(f"Table 1 adjusted effects to:        {density_res['table1_tsv']}")
    print(f"Ontology redundant pairs to:        {ontology_res['pairs_tsv']}")
    print(f"Ontology redundancy summary to:     {ontology_res['summary_tsv']}")
    print(f"Ontology collapsing comparison to:  {ontology_res['comparison_tsv']}")
    print(f"Consolidated comprehensive report:  {report_path}")
    print("\nKey Findings:")
    raw = density_res["raw_pearson"]
    adj = density_res["adj_pc1_pearson"]
    twin = density_res["twin_summary"]
    t1_by_id = {eff.trait_id: eff for eff in density_res["table1_effects"]}
    aerobic_label = (
        t1_by_id["METPO:1000608"].trait_label
        if "METPO:1000608" in t1_by_id
        else "facultatively aerobic"
    )
    print(
        f"  - Baseline PC1 explains "
        f"{format_correlation_metric(raw, 'pc1_variance_pct', '.2f', suffix='%')} of variance "
        f"across {raw.n_traits} traits."
    )
    print(
        f"  - Post-adjustment median pairwise r: "
        f"{format_correlation_metric(adj, 'median_pairwise_r', '+.4f')} "
        f"(raw: {format_correlation_metric(raw, 'median_pairwise_r', '+.4f')})."
    )
    print(
        f"  - {aerobic_label} x moderately halophilic r: "
        f"{format_correlation_metric(raw, 'test_pair_r', '+.4f')} "
        f"-> {format_correlation_metric(adj, 'test_pair_r', '+.4f')} "
        f"(change: {format_metric_diff(adj, raw)})."
    )
    print(
        f"  - Near-twins (coccus x sphere): r = "
        f"{format_report_number(twin.profile_correlation, '+.4f')}, "
        f"VIF = {format_report_number(twin.vif_profile, '.2f', allow_positive_infinity=True)}."
    )
    n_clusters_label = _plural(ontology_res["n_clusters"], "cluster", "clusters")
    print(
        f"  - Ontology Redundancy: {ontology_res['n_redundant_pairs']} redundant "
        f"{_plural(ontology_res['n_redundant_pairs'], 'pair', 'pairs')} across "
        f"{ontology_res['n_clusters']} {n_clusters_label}."
    )
    print(
        f"  - Distinct Trait Space Bounded: {ontology_res['bounded_trait_space']} "
        f"{_plural(ontology_res['bounded_trait_space'], 'feature', 'features')} "
        f"(collapsed {ontology_res['n_collapsed_features']} "
        f"{_plural(ontology_res['n_collapsed_features'], 'feature', 'features')} from "
        f"{raw.n_traits})."
    )
    max_pc1 = format_correlation_metric(
        ontology_res["max_metrics"], "pc1_variance_pct", ".2f", suffix="%"
    )
    raw_pc1 = format_correlation_metric(raw, "pc1_variance_pct", ".2f", suffix="%")
    adj_col_pc1 = format_correlation_metric(
        ontology_res["adj_collapsed_metrics"], "pc1_variance_pct", ".2f", suffix="%"
    )
    print(
        f"  - Collapsed Matrix PC1 Variance: {max_pc1} "
        f"(raw: {raw_pc1}; collapsed + PC1 regression: {adj_col_pc1})."
    )
    print("=" * 80 + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
