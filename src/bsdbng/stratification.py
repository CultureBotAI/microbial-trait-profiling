"""Background stratification and annotation-density sensitivity evaluation.

Addresses Issue #118:
1. Partition study trait profiles into stratified universes (pooled, feces, feces_excluded,
   or specific user-selected anatomical site strata).
2. Run hypergeometric ORA enrichment per stratum, computing background universes (M, K)
   within each stratum.
3. Export ranked tables under stable filenames:
   - trait_enrichment_traits_overall_ranked_<stratum>.tsv
   - trait_enrichment_trait_study_pairs_ranked_<stratum>.tsv
   carrying n_profiles_tested, trait_total, and body_site_group.
4. Generate stratification sensitivity summary (TSV and text report) evaluating:
   - Annotation-density diagnostic deciles.
   - Cross-stratum concordance, retained counts, effect heterogeneity, and support states.
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from bsdbng.paths import DERIVED_DATA_DIR, REPO_ROOT

DEFAULT_PROFILES_PATH = (
    DERIVED_DATA_DIR / "study_trait_profiles_ora_wide.tsv"
    if (DERIVED_DATA_DIR / "study_trait_profiles_ora_wide.tsv").exists()
    else DERIVED_DATA_DIR / "study_trait_profiles_wide.tsv"
)
DEFAULT_TAXA_MATRIX_PATH = DERIVED_DATA_DIR / "taxa_trait_matrix.tsv"
DEFAULT_EXPERIMENTS_PATH = DERIVED_DATA_DIR / "experiments.tsv"
DEFAULT_OUTPUT_DIR = DERIVED_DATA_DIR

PRIMARY_STRATA = ("pooled", "feces", "feces_excluded")

PRECAUTIONARY_NOTE = (
    "Avoid equating a nonsignificant non-fecal result with evidence that the association "
    "is absent or does not generalize; effect heterogeneity and uncertainty need direct "
    "assessment. Significant-set comparisons alone cannot separate background changes from "
    "reduced sample size."
)

ANNOTATION_DENSITY_NOTE = (
    "The annotation-density diagnostic is descriptive. Background counts alone do not establish "
    "the direction or magnitude of calibration bias, and a proposed density-matched permutation "
    "needs a justified null and exchangeability assumptions."
)


def is_feces_site(raw_site: str | None) -> bool:
    """Return True if the anatomical site corresponds to feces or stool."""
    if not raw_site or pd.isna(raw_site):
        return False
    site_str = str(raw_site).strip().lower()
    tokens = [t.strip() for t in site_str.split(",")]
    return "feces" in tokens or "stool" in tokens or site_str in ("feces", "stool")


def classify_site_stratum(raw_site: str | None) -> str:
    """Classify a raw body site string into 'feces' or 'feces_excluded'."""
    return "feces" if is_feces_site(raw_site) else "feces_excluded"


def partition_profiles(
    profiles_path: Path,
    experiments_path: Path,
    output_dir: Path,
    strata: list[str] | tuple[str, ...] = PRIMARY_STRATA,
) -> dict[str, Path]:
    """Partition wide study trait profiles by anatomical site stratum.

    Writes:
      output_dir / f"study_trait_profiles_{stratum}.tsv"
    for each stratum in strata. Returns a dictionary mapping stratum name to Path.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    profiles_df = pd.read_csv(profiles_path, sep="\t")
    experiments_df = pd.read_csv(experiments_path, sep="\t")

    meta_cols = ["study_id", "experiment_id", "body_site"]
    exp_subset = experiments_df[meta_cols].drop_duplicates()
    merged = profiles_df.merge(exp_subset, on=["study_id", "experiment_id"], how="left")

    partition_paths: dict[str, Path] = {}
    profile_cols = [c for c in profiles_df.columns if c != "body_site"]

    for stratum in strata:
        if stratum == "pooled":
            stratum_df = profiles_df[profile_cols]
        elif stratum == "feces":
            mask = merged["body_site"].apply(is_feces_site)
            stratum_df = merged[mask][profile_cols]
        elif stratum == "feces_excluded":
            mask = ~merged["body_site"].apply(is_feces_site)
            stratum_df = merged[mask][profile_cols]
        else:
            # Custom stratum: match exact body_site or comma-separated token
            target_stratum = stratum.strip().lower()

            def matches_custom(site: Any, target: str = target_stratum) -> bool:
                if not site or pd.isna(site):
                    return False
                tokens = [t.strip().lower() for t in str(site).split(",")]
                return target in tokens or str(site).strip().lower() == target

            mask = merged["body_site"].apply(matches_custom)
            stratum_df = merged[mask][profile_cols]

        if stratum_df.empty:
            raise ValueError(
                f"Stratum '{stratum}' has 0 matching profiles in {profiles_path} "
                f"using experiments from {experiments_path}."
            )

        out_path = output_dir / f"study_trait_profiles_{stratum}.tsv"
        stratum_df.to_csv(out_path, sep="\t", index=False)
        partition_paths[stratum] = out_path

    return partition_paths


@dataclass(frozen=True)
class StratumOutputs:
    """Output files for a single stratum."""

    stratum: str
    profiles_tsv: Path
    results_tsv: Path
    report_txt: Path
    overall_ranked_tsv: Path
    pairs_ranked_tsv: Path
    all_results_tsv: Path | None = None


def _run_script(script_path: Path, args: list[str]) -> None:
    mod_name = f"_strat_{script_path.stem}"
    if mod_name in sys.modules:
        module = sys.modules[mod_name]
    else:
        spec = importlib.util.spec_from_file_location(mod_name, script_path)
        if spec is None or spec.loader is None:
            raise ImportError(f"Cannot load script: {script_path}")
        module = importlib.util.module_from_spec(spec)
        sys.modules[mod_name] = module
        spec.loader.exec_module(module)
    ret = module.main(args)
    if ret is not None and ret != 0:
        raise RuntimeError(f"{script_path.name} exited with status {ret}")


def run_stratum_enrichment(
    stratum: str,
    profiles_tsv: Path,
    output_dir: Path,
    p_threshold: float = 0.001,
    min_fold: float = 1.5,
    fdr: bool = True,
    fdr_scope: str = "global",
    taxa_matrix: Path | None = None,
) -> StratumOutputs:
    """Execute ORA enrichment and rankings export for one stratum."""
    output_dir.mkdir(parents=True, exist_ok=True)
    results_tsv = output_dir / f"trait_enrichment_results_{stratum}.tsv"
    report_txt = output_dir / f"trait_enrichment_report_{stratum}.txt"
    all_results_tsv = output_dir / f"trait_enrichment_all_results_{stratum}.tsv"
    overall_ranked_tsv = output_dir / f"trait_enrichment_traits_overall_ranked_{stratum}.tsv"
    pairs_ranked_tsv = output_dir / f"trait_enrichment_trait_study_pairs_ranked_{stratum}.tsv"

    resolved_matrix = taxa_matrix or DEFAULT_TAXA_MATRIX_PATH
    if not resolved_matrix.exists():
        raise FileNotFoundError(
            f"Taxa-trait matrix not found at '{resolved_matrix}'. "
            "Pass --taxa-matrix or build it via `just taxa-trait-matrix`."
        )

    analyzer_args: list[str] = [
        "--input",
        str(profiles_tsv),
        "--output",
        str(report_txt),
        "--tsv-output",
        str(results_tsv),
        "--all-results-output",
        str(all_results_tsv),
        "--taxa-matrix",
        str(resolved_matrix),
        "--p-threshold",
        str(p_threshold),
        "--min-fold",
        str(min_fold),
        "--fdr-scope",
        fdr_scope,
    ]
    if not fdr:
        analyzer_args.append("--no-fdr")

    _run_script(REPO_ROOT / "scripts" / "trait_enrichment.py", analyzer_args)

    exporter_args: list[str] = [
        "--input",
        str(report_txt),
        "--results",
        str(results_tsv),
        "--profiles",
        str(profiles_tsv),
        "--overall-output",
        str(overall_ranked_tsv),
        "--pair-output",
        str(pairs_ranked_tsv),
        "--body-site-group",
        stratum,
    ]
    _run_script(REPO_ROOT / "scripts" / "export_trait_enrichment_rankings.py", exporter_args)

    return StratumOutputs(
        stratum=stratum,
        profiles_tsv=profiles_tsv,
        results_tsv=results_tsv,
        report_txt=report_txt,
        overall_ranked_tsv=overall_ranked_tsv,
        pairs_ranked_tsv=pairs_ranked_tsv,
        all_results_tsv=all_results_tsv,
    )


def compute_annotation_density_deciles(
    pooled_overall_df: pd.DataFrame, n_deciles: int = 10
) -> pd.DataFrame:
    """Bin traits into background-annotation-density deciles based on trait_total."""
    df = pooled_overall_df.copy()
    empty_deciles = pd.DataFrame(
        columns=[
            "decile",
            "min_trait_total",
            "max_trait_total",
            "total_traits",
            "significant_traits",
            "pct_significant",
        ]
    )
    if df.empty or "trait_total" not in df.columns:
        return empty_deciles

    actual_bins = min(n_deciles, len(df))
    if actual_bins <= 0:
        return empty_deciles

    # Compute rank-based quantile bins to handle duplicate counts gracefully
    ranks = df["trait_total"].rank(method="first")
    df["density_decile"] = pd.qcut(ranks, q=actual_bins, labels=list(range(1, actual_bins + 1)))

    rows: list[dict[str, Any]] = []
    for decile, group in df.groupby("density_decile", observed=False):
        t_total = group["trait_total"]
        n_traits = len(group)
        stat_col = "best_qvalue" if "best_qvalue" in group.columns else "best_pvalue"
        sig_count = (group[stat_col] <= 0.05).sum() if stat_col in group.columns else n_traits
        rows.append(
            {
                "decile": int(decile),
                "min_trait_total": int(t_total.min()) if n_traits > 0 else 0,
                "max_trait_total": int(t_total.max()) if n_traits > 0 else 0,
                "total_traits": n_traits,
                "significant_traits": int(sig_count),
                "pct_significant": (float(sig_count) / n_traits * 100.0) if n_traits > 0 else 0.0,
            }
        )

    return pd.DataFrame(rows)


def build_sensitivity_summary_table(
    stratum_outputs: dict[str, StratumOutputs],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build the cross-stratum per-trait sensitivity table and stratum metadata table.

    Returns:
        (per_trait_df, strata_summary_df)
    """
    overall_dfs: dict[str, pd.DataFrame] = {}
    meta_records: list[dict[str, Any]] = []

    for stratum, outputs in stratum_outputs.items():
        if outputs.overall_ranked_tsv.exists():
            df = pd.read_csv(outputs.overall_ranked_tsv, sep="\t", comment="#")
            overall_dfs[stratum] = df
        else:
            overall_dfs[stratum] = pd.DataFrame()

        # Read comment header metadata
        n_testable = 0
        n_ranked = 0
        if outputs.overall_ranked_tsv.exists():
            with outputs.overall_ranked_tsv.open("r", encoding="utf-8") as f:
                for line in f:
                    if line.startswith("# Total-Testable-Profiles:"):
                        n_testable = int(line.split(":", 1)[1].strip())
                    elif line.startswith("# Ranked-Traits:"):
                        n_ranked = int(line.split(":", 1)[1].strip())

        meta_records.append(
            {
                "stratum": stratum,
                "total_testable_profiles": n_testable,
                "ranked_traits": n_ranked,
                "correction_family": "Benjamini-Hochberg (global)",
            }
        )

    strata_summary_df = pd.DataFrame(meta_records)

    pooled_df = overall_dfs.get("pooled", pd.DataFrame())
    feces_df = overall_dfs.get("feces", pd.DataFrame())
    feces_ex_df = overall_dfs.get("feces_excluded", pd.DataFrame())

    all_traits = sorted(
        set().union(*(set(df["trait"]) for df in overall_dfs.values() if "trait" in df.columns))
    )

    def _get_row(df: pd.DataFrame, trait_name: str) -> pd.Series | None:
        if df.empty or "trait" not in df.columns:
            return None
        match = df[df["trait"] == trait_name]
        if match.empty:
            return None
        return match.iloc[0]

    trait_rows: list[dict[str, Any]] = []
    for trait in all_traits:
        p_row = _get_row(pooled_df, trait)
        f_row = _get_row(feces_df, trait)
        fe_row = _get_row(feces_ex_df, trait)

        p_trait_tot = (
            int(p_row["trait_total"]) if p_row is not None and "trait_total" in p_row else 0
        )
        p_tested = (
            int(p_row["n_profiles_tested"])
            if p_row is not None and "n_profiles_tested" in p_row
            else 0
        )
        f_tested = (
            int(f_row["n_profiles_tested"])
            if f_row is not None and "n_profiles_tested" in f_row
            else 0
        )
        fe_tested = (
            int(fe_row["n_profiles_tested"])
            if fe_row is not None and "n_profiles_tested" in fe_row
            else 0
        )

        p_prev = float(p_row["prevalence"]) if p_row is not None and "prevalence" in p_row else 0.0
        f_prev = float(f_row["prevalence"]) if f_row is not None and "prevalence" in f_row else 0.0
        fe_prev = (
            float(fe_row["prevalence"]) if fe_row is not None and "prevalence" in fe_row else 0.0
        )

        p_fold = float(p_row["mean_fold"]) if p_row is not None and "mean_fold" in p_row else np.nan
        f_fold = float(f_row["mean_fold"]) if f_row is not None and "mean_fold" in f_row else np.nan
        fe_fold = (
            float(fe_row["mean_fold"]) if fe_row is not None and "mean_fold" in fe_row else np.nan
        )

        def _extract_q_or_p(row: Any) -> float:
            if row is None:
                return np.nan
            if "best_qvalue" in row and pd.notna(row["best_qvalue"]):
                return float(row["best_qvalue"])
            if "best_pvalue" in row and pd.notna(row["best_pvalue"]):
                return float(row["best_pvalue"])
            return np.nan

        p_q = _extract_q_or_p(p_row)
        f_q = _extract_q_or_p(f_row)
        fe_q = _extract_q_or_p(fe_row)

        p_sig = bool(p_row is not None and (pd.isna(p_q) or p_q <= 0.05))
        f_sig = bool(f_row is not None and (pd.isna(f_q) or f_q <= 0.05))
        fe_sig = bool(fe_row is not None and (pd.isna(fe_q) or fe_q <= 0.05))

        # Determine support status
        if f_tested < 2 or fe_tested < 2:
            support_status = "insufficient_support"
        elif f_sig and fe_sig:
            support_status = "both_strata_significant"
        elif f_sig and not fe_sig:
            support_status = "feces_only_significant"
        elif fe_sig and not f_sig:
            support_status = "feces_excluded_only_significant"
        elif p_sig:
            support_status = "pooled_only_significant"
        else:
            support_status = "neither_significant"

        delta_fold = (f_fold - fe_fold) if (pd.notna(f_fold) and pd.notna(fe_fold)) else np.nan
        delta_prev = (f_prev - fe_prev) if (pd.notna(f_prev) and pd.notna(fe_prev)) else np.nan

        trait_rows.append(
            {
                "trait": trait,
                "trait_total_pooled": p_trait_tot,
                "pooled_n_tested": p_tested,
                "pooled_prevalence": p_prev,
                "pooled_best_stat": p_q,
                "pooled_mean_fold": p_fold,
                "pooled_significant": p_sig,
                "feces_n_tested": f_tested,
                "feces_prevalence": f_prev,
                "feces_best_stat": f_q,
                "feces_mean_fold": f_fold,
                "feces_significant": f_sig,
                "feces_excluded_n_tested": fe_tested,
                "feces_excluded_prevalence": fe_prev,
                "feces_excluded_best_stat": fe_q,
                "feces_excluded_mean_fold": fe_fold,
                "feces_excluded_significant": fe_sig,
                "support_status": support_status,
                "delta_mean_fold": delta_fold,
                "delta_prevalence": delta_prev,
            }
        )

    per_trait_df = pd.DataFrame(trait_rows)
    return per_trait_df, strata_summary_df


def generate_sensitivity_reports(
    stratum_outputs: dict[str, StratumOutputs],
    output_dir: Path,
) -> tuple[Path, Path]:
    """Generate stratification sensitivity summary TSV and text reports."""
    output_dir.mkdir(parents=True, exist_ok=True)
    tsv_path = output_dir / "stratification_sensitivity_summary.tsv"
    txt_path = output_dir / "stratification_sensitivity_summary.txt"

    per_trait_df, strata_summary_df = build_sensitivity_summary_table(stratum_outputs)

    # Compute annotation density deciles from pooled overall table
    pooled_df = pd.DataFrame()
    if "pooled" in stratum_outputs and stratum_outputs["pooled"].overall_ranked_tsv.exists():
        pooled_df = pd.read_csv(stratum_outputs["pooled"].overall_ranked_tsv, sep="\t", comment="#")
    decile_df = compute_annotation_density_deciles(pooled_df)

    # Write TSV with metadata comments
    with tsv_path.open("w", encoding="utf-8", newline="") as f:
        f.write(f"# Note: {PRECAUTIONARY_NOTE}\n")
        f.write(f"# Diagnostic: {ANNOTATION_DENSITY_NOTE}\n")
        f.write(f"# Evaluated-Strata: {','.join(stratum_outputs.keys())}\n")
        per_trait_df.to_csv(f, sep="\t", index=False)

    # Write comprehensive text report
    support_counts = (
        per_trait_df["support_status"].value_counts().to_dict()
        if "support_status" in per_trait_df.columns
        else {}
    )

    txt_lines: list[str] = [
        "=" * 80,
        "ORA BACKGROUND STRATIFICATION AND ANNOTATION-DENSITY SENSITIVITY REPORT",
        "=" * 80,
        "",
        "1. METHODOLOGICAL FRAMEWORK & PRECAUTIONARY INTERPRETATION",
        "-" * 80,
        f"IMPORTANT: {PRECAUTIONARY_NOTE}",
        "",
        f"DIAGNOSTIC NOTE: {ANNOTATION_DENSITY_NOTE}",
        "",
        "2. STRATUM PROFILE AND REPAIR COUNTS",
        "-" * 80,
    ]

    for _, row in strata_summary_df.iterrows():
        txt_lines.append(
            f"Stratum: {row['stratum']:<15} | "
            f"Testable Profiles: {row['total_testable_profiles']:<6} | "
            f"Ranked Traits: {row['ranked_traits']:<5} | "
            f"Correction: {row['correction_family']}"
        )

    txt_lines.extend(
        [
            "",
            "3. ANNOTATION-DENSITY SENSITIVITY DIAGNOSTIC (POOLED BACKGROUND K)",
            "-" * 80,
            "Traits binned by background annotation count (trait_total) deciles:",
            (
                f"{'Decile':<8} {'Min K':<8} {'Max K':<8} {'Total Traits':<14} "
                f"{'Sig Traits':<12} {'% Sig':<8}"
            ),
        ]
    )

    for _, d in decile_df.iterrows():
        txt_lines.append(
            f"{int(d['decile']):<8} {int(d['min_trait_total']):<8} {int(d['max_trait_total']):<8} "
            f"{int(d['total_traits']):<14} {int(d['significant_traits']):<12} "
            f"{d['pct_significant']:<8.1f}%"
        )

    txt_lines.extend(
        [
            "",
            "4. CROSS-STRATUM CONCORDANCE & SUPPORT STATES",
            "-" * 80,
            f"Total evaluated traits across all strata: {len(per_trait_df)}",
        ]
    )
    for status, count in sorted(support_counts.items(), key=lambda x: -x[1]):
        pct = count / len(per_trait_df) * 100.0 if len(per_trait_df) > 0 else 0.0
        txt_lines.append(f"  - {status:<32}: {count:>4} ({pct:>5.1f}%)")

    txt_lines.extend(
        [
            "",
            "5. OUTPUT ARTIFACTS",
            "-" * 80,
            f"  - Stratification sensitivity TSV : {tsv_path}",
            f"  - Stratification sensitivity text: {txt_path}",
        ]
    )
    for stratum, outs in stratum_outputs.items():
        txt_lines.append(f"  - Stratum '{stratum}' overall ranking  : {outs.overall_ranked_tsv}")
        txt_lines.append(f"  - Stratum '{stratum}' pair rankings     : {outs.pairs_ranked_tsv}")

    txt_lines.append("=" * 80 + "\n")

    txt_path.write_text("\n".join(txt_lines), encoding="utf-8")
    return tsv_path, txt_path


def run_stratification_pipeline(
    profiles_path: Path = DEFAULT_PROFILES_PATH,
    experiments_path: Path = DEFAULT_EXPERIMENTS_PATH,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    p_threshold: float = 0.001,
    min_fold: float = 1.5,
    fdr: bool = True,
    fdr_scope: str = "global",
    strata: list[str] | tuple[str, ...] = PRIMARY_STRATA,
    taxa_matrix: Path | None = None,
) -> tuple[dict[str, StratumOutputs], Path, Path]:
    """Execute the full background stratification pipeline."""
    print("=" * 80)
    print("Running ORA Background Stratification Pipeline (#118)")
    print(f"Strata: {', '.join(strata)}")
    print(f"Input profiles: {profiles_path}")
    print(f"Experiments: {experiments_path}")
    print(f"Output directory: {output_dir}")
    print("=" * 80)

    # 1. Partition profiles
    print("\n[Step 1/3] Partitioning profiles by anatomical site stratum...")
    partitions = partition_profiles(
        profiles_path=profiles_path,
        experiments_path=experiments_path,
        output_dir=output_dir,
        strata=strata,
    )
    for stratum, path in partitions.items():
        print(f"  - Stratum '{stratum}': {path}")

    # 2. Run enrichment and ranking export per stratum
    print("\n[Step 2/3] Executing hypergeometric enrichment per stratum...")
    stratum_outputs: dict[str, StratumOutputs] = {}
    for stratum, prof_tsv in partitions.items():
        print(f"\n--- Running stratum: {stratum} ---")
        outputs = run_stratum_enrichment(
            stratum=stratum,
            profiles_tsv=prof_tsv,
            output_dir=output_dir,
            p_threshold=p_threshold,
            min_fold=min_fold,
            fdr=fdr,
            fdr_scope=fdr_scope,
            taxa_matrix=taxa_matrix,
        )
        stratum_outputs[stratum] = outputs

    # 3. Generate sensitivity comparison and diagnostics
    print("\n[Step 3/3] Generating sensitivity comparison summary & density diagnostics...")
    tsv_path, txt_path = generate_sensitivity_reports(
        stratum_outputs=stratum_outputs,
        output_dir=output_dir,
    )
    print(f"Wrote sensitivity summary TSV : {tsv_path}")
    print(f"Wrote sensitivity summary text: {txt_path}")
    print("Stratification pipeline complete.\n")

    return stratum_outputs, tsv_path, txt_path


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Evaluate ORA background stratification and annotation-density sensitivity (#118)"
        )
    )
    parser.add_argument(
        "--profiles",
        type=Path,
        default=DEFAULT_PROFILES_PATH,
        help="Path to study trait profiles wide TSV.",
    )
    parser.add_argument(
        "--experiments",
        type=Path,
        default=DEFAULT_EXPERIMENTS_PATH,
        help="Path to experiments TSV containing body_site metadata.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory to write stratified outputs and sensitivity summaries.",
    )
    parser.add_argument(
        "--p-threshold",
        type=float,
        default=0.001,
        help="Significance threshold (applied to FDR if used, else p-value).",
    )
    parser.add_argument(
        "--min-fold",
        type=float,
        default=1.5,
        help="Minimum fold enrichment threshold.",
    )
    parser.add_argument(
        "--fdr",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Apply Benjamini-Hochberg FDR correction (default: True).",
    )
    parser.add_argument(
        "--fdr-scope",
        choices=["global", "profile"],
        default="global",
        help="Scope for FDR correction (default: global).",
    )
    parser.add_argument(
        "--strata",
        nargs="+",
        default=list(PRIMARY_STRATA),
        help="List of strata to evaluate (default: pooled feces feces_excluded).",
    )
    parser.add_argument(
        "--site",
        type=str,
        default=None,
        help="Evaluate a single anatomical site stratum (overrides --strata).",
    )
    parser.add_argument(
        "--taxa-matrix",
        type=Path,
        default=DEFAULT_TAXA_MATRIX_PATH,
        help="Path to taxa_trait_matrix.tsv (default: data/derived/taxa_trait_matrix.tsv).",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    strata = [args.site] if args.site is not None else args.strata
    run_stratification_pipeline(
        profiles_path=args.profiles,
        experiments_path=args.experiments,
        output_dir=args.output_dir,
        p_threshold=args.p_threshold,
        min_fold=args.min_fold,
        fdr=args.fdr,
        fdr_scope=args.fdr_scope,
        strata=strata,
        taxa_matrix=args.taxa_matrix,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
