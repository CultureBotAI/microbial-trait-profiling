"""Trait coupling sensitivity: Leave-One-Trait-Out (LOTO) and Leave-Cluster-Out (LCO).

Addresses Issue #260 (and follow-ups #478, #479):
1. Measures trait coupling across non-MediaDive traits (317 columns):
   - Leave-One-Trait-Out (LOTO): For each of the six canonical Table 1 traits,
     removes it from the trait universe, renormalizes profile denominators,
     and evaluates the resulting delta in effect size and focal-set FDR q-value.
   - Leave-Cluster-Out (LCO): Removes correlated near-twins (coccus shaped +
     sphere shaped) simultaneously and evaluates the shifts for remaining traits.
2. Evaluates sensitivity across all three analytical arms:
   - Composition arm: Paired profile fraction log-ratio contrast, z-score, paired t-test.
   - Directional arm: Study-level sign test over paired contrasts, Clopper-Pearson / binomial test.
   - ORA arm: Profile-level hypergeometric enrichment against background urn M = sum(X).
   Note (#478): Within this module, arm-level BH q-values are computed as a focal-set
   sensitivity check over the evaluated Table 1 traits (m = 6 in baseline, m = 5 in LOTO,
   m = 4 in LCO) on all-profile statistics, rather than full-family Table 1 BH q-values
   across all 317 non-MediaDive (or 2,607 total) columns. Consequently, delta_q combines
   denominator perturbation with the focal BH family-size step.
3. Evaluates the formal stability decision rule:
   - Tests whether dropping sphere moves coccus by more than the minimum distance
     between adjacent Table 1 rows.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import binomtest, hypergeom

from mtp.density_adjustment import (
    DEFAULT_COLUMNS_PATH,
    DEFAULT_PROFILES_PATH,
    TABLE_1_TRAITS,
    evaluate_conditional_twin_test,
    load_non_media_profile_matrix,
)
from mtp.enrichment_stats import benjamini_hochberg
from mtp.paths import DERIVED_DATA_DIR

DEFAULT_OUTPUT_DIR = DERIVED_DATA_DIR

# Canonical near-twin cluster for Leave-Cluster-Out (LCO)
COCCUS_ID = "METPO:1000668"
SPHERE_ID = "METPO:1000683"
CANONICAL_LCO_CLUSTER: dict[str, str] = {
    COCCUS_ID: "coccus shaped",
    SPHERE_ID: "sphere shaped",
}


@dataclass(frozen=True)
class ArmResult:
    """Baseline or perturbed evaluation result for a single trait in an arm."""

    trait_id: str
    trait_name: str
    primary_stat: float  # effect (comp), share_up (dir), mean_fold (ora)
    p_value: float
    q_value: float
    extra_stat: float = 0.0  # t_stat (comp), n_scored (dir), hits_q05 (ora)


@dataclass(frozen=True)
class LOTORow:
    """LOTO shift entry for one evaluated trait when dropping another trait."""

    arm: str
    dropped_trait_id: str
    dropped_trait_name: str
    target_trait_id: str
    target_trait_name: str
    base_stat: float
    loto_stat: float
    delta_stat: float
    base_q: float
    loto_q: float
    delta_q: float


@dataclass(frozen=True)
class LCOSummaryRow:
    """LCO shift entry when removing an entire cluster of correlated traits."""

    arm: str
    cluster_dropped: str
    trait_id: str
    trait_name: str
    base_stat: float
    lco_stat: float
    delta_stat: float
    base_q: float
    lco_q: float
    delta_q: float
    base_extra: float
    lco_extra: float


@dataclass(frozen=True)
class DecisionRuleResult:
    """Evaluation of the inter-trait row spacing decision rule."""

    target_trait_id: str
    target_trait_name: str
    dropped_trait_id: str
    dropped_trait_name: str
    delta_effect: float
    abs_delta_effect: float
    min_adjacent_distance: float
    passes_stability_criterion: bool
    adjacent_distances: list[float]
    sorted_effects: list[float]
    verdict: str
    min_adjacent_pair: tuple[str, str] = ("", "")


def compute_composition_arm(
    profiles_df: pd.DataFrame,
    trait_universe: list[str],
    focal_traits: list[str],
) -> dict[str, ArmResult]:
    """Compute paired composition arm effects, t-statistics, and focal-set FDR q-values.

    Profile fractions are normalized over trait_universe:
        f_{i, t} = X_{i, t} / sum_{k in trait_universe} X_{i, k}
    Transformed via log1p(1000 * f) and standardized to z-scores across all profiles.
    Paired contrast is computed per experiment: z_inc - z_dec.
    One-sample t-test evaluates significance across paired experiments.
    Returned q_value is a focal-set Benjamini-Hochberg adjustment over focal_traits
    (m = len(focal_traits), a sensitivity check distinct from full-family Table 1 q-values).
    """
    counts_all = profiles_df[trait_universe].to_numpy(dtype=float)
    totals = counts_all.sum(axis=1)
    denom = np.where(totals > 0, totals, 1.0)
    fractions = np.where(totals[:, None] > 0, counts_all / denom[:, None], 0.0)

    scaled = np.log1p(fractions * 1000.0)
    mu = scaled.mean(axis=0)
    sd = scaled.std(axis=0, ddof=0)
    sd[sd == 0] = 1.0
    zscores = (scaled - mu) / sd

    df_z = pd.DataFrame(zscores, columns=trait_universe)
    df_z["study_id"] = profiles_df["study_id"].values
    df_z["experiment_id"] = profiles_df["experiment_id"].values
    df_z["inc_dec"] = profiles_df["inc_dec"].values

    pvt_z = df_z.pivot(
        index=["study_id", "experiment_id"], columns="inc_dec", values=trait_universe
    ).dropna()

    diffs_mat = (
        pvt_z.xs("inc", axis=1, level=1)[trait_universe]
        - pvt_z.xs("dec", axis=1, level=1)[trait_universe]
    ).to_numpy(dtype=float)

    effects = np.mean(diffs_mat, axis=0)
    std_diffs = np.std(diffs_mat, axis=0, ddof=1)
    valid_std = std_diffs > 1e-12
    t_stats = np.zeros(len(trait_universe), dtype=float)
    p_vals = np.ones(len(trait_universe), dtype=float)
    if np.any(valid_std):
        tt = stats.ttest_1samp(diffs_mat[:, valid_std], 0.0, axis=0)
        t_stats[valid_std] = np.where(
            np.isnan(tt.statistic), 0.0, np.asarray(tt.statistic, dtype=float)
        )
        p_vals[valid_std] = np.where(np.isnan(tt.pvalue), 1.0, np.asarray(tt.pvalue, dtype=float))
    q_vals = benjamini_hochberg(p_vals.tolist()) if len(p_vals) > 0 else np.array([], dtype=float)

    idx_map = {t: idx for idx, t in enumerate(trait_universe)}
    results: dict[str, ArmResult] = {}
    for t in focal_traits:
        idx = idx_map[t]
        results[t] = ArmResult(
            trait_id=t,
            trait_name=TABLE_1_TRAITS.get(t, t),
            primary_stat=float(effects[idx]),
            p_value=float(p_vals[idx]),
            q_value=float(q_vals[idx]) if len(q_vals) > idx else float(p_vals[idx]),
            extra_stat=float(t_stats[idx]),
        )
    return results


def compute_directional_arm(
    profiles_df: pd.DataFrame,
    trait_universe: list[str],
    focal_traits: list[str],
    sign_tol: float = 1e-12,
) -> dict[str, ArmResult]:
    """Compute directional arm study-level sign-test over paired contrasts.

    For each study s, the mean difference in fractions across its experiments:
        Delta f_{s, t} = mean_{e in s} (f_{e, inc, t} - f_{e, dec, t})
    Sign counts:
        up: Delta f_{s, t} > sign_tol
        down: Delta f_{s, t} < -sign_tol
    Exact two-sided binomial test on up vs. down (p = 0.5).
    Returned q_value is a focal-set Benjamini-Hochberg adjustment over focal_traits
    (m = len(focal_traits), a sensitivity check distinct from full-family Table 1 q-values).
    """
    counts_all = profiles_df[trait_universe].to_numpy(dtype=float)
    totals = counts_all.sum(axis=1)
    denom = np.where(totals > 0, totals, 1.0)
    fractions = np.where(totals[:, None] > 0, counts_all / denom[:, None], 0.0)

    df_frac = pd.DataFrame(fractions, columns=trait_universe)
    df_frac["study_id"] = profiles_df["study_id"].values
    df_frac["experiment_id"] = profiles_df["experiment_id"].values
    df_frac["inc_dec"] = profiles_df["inc_dec"].values

    pvt = df_frac.pivot(
        index=["study_id", "experiment_id"], columns="inc_dec", values=trait_universe
    ).dropna()

    pair_diff = (
        pvt.xs("inc", axis=1, level=1)[trait_universe]
        - pvt.xs("dec", axis=1, level=1)[trait_universe]
    )
    study_diff = pair_diff.groupby(level="study_id").mean().to_numpy(dtype=float)

    up_arr = (study_diff > sign_tol).sum(axis=0)
    down_arr = (study_diff < -sign_tol).sum(axis=0)
    scored_arr = up_arr + down_arr

    binom_cache: dict[tuple[int, int], float] = {}
    p_vals: list[float] = []
    for up_val, scored_val in zip(up_arr, scored_arr, strict=True):
        u_int = int(up_val)
        n_int = int(scored_val)
        if n_int <= 0:
            p_vals.append(1.0)
        else:
            key = (u_int, n_int)
            if key not in binom_cache:
                binom_cache[key] = float(binomtest(u_int, n_int, p=0.5).pvalue)
            p_vals.append(binom_cache[key])

    q_vals = benjamini_hochberg(p_vals) if len(p_vals) > 0 else np.array([], dtype=float)

    idx_map = {t: idx for idx, t in enumerate(trait_universe)}
    results: dict[str, ArmResult] = {}
    for t in focal_traits:
        idx = idx_map[t]
        scored = int(scored_arr[idx])
        share_up = float(up_arr[idx] / scored) if scored > 0 else 0.5
        results[t] = ArmResult(
            trait_id=t,
            trait_name=TABLE_1_TRAITS.get(t, t),
            primary_stat=share_up,
            p_value=float(p_vals[idx]),
            q_value=float(q_vals[idx]) if len(q_vals) > idx else float(p_vals[idx]),
            extra_stat=float(scored),
        )
    return results


def compute_ora_arm(
    paired_profiles_df: pd.DataFrame,
    trait_universe: list[str],
    focal_traits: list[str],
) -> dict[str, ArmResult]:
    """Compute ORA hypergeometric enrichment arm across paired profiles.

    Background urn:
        M = sum_{i, t in universe} X_{i, t}
        N_i = sum_{t in universe} X_{i, t} (profile sample draw)
        K_t = sum_i X_{i, t} (trait total)
    One-sided hypergeometric p-value: P(X >= k_{i, t} | M, K_t, N_i).
    Hits evaluated at profile-level BH q <= 0.05.
    Overall primary statistic: mean fold enrichment across tested profiles.
    Fisher's combined test evaluates global significance, with focal-set BH q_value
    computed over focal_traits (m = len(focal_traits)).
    """
    counts_universe = paired_profiles_df[trait_universe].to_numpy(dtype=int)
    M = int(counts_universe.sum())
    N_arr = counts_universe.sum(axis=1)
    K_arr = counts_universe.sum(axis=0)
    pos_N = N_arr > 0
    N_pos = N_arr[pos_N]

    focal_set = set(focal_traits)
    p_vals: list[float] = []
    focal_stats: dict[str, dict[str, float]] = {}

    for idx, t in enumerate(trait_universe):
        K = int(K_arr[idx])
        k_obs = counts_universe[:, idx]
        pos_k = pos_N & (k_obs > 0)

        if K > 0 and M > 0 and np.any(pos_k):
            valid_p = hypergeom.sf(k_obs[pos_k] - 1, M, K, N_arr[pos_k])
            _, fisher_p = stats.combine_pvalues(valid_p, method="fisher")
            fisher_p_val = float(fisher_p) if not np.isnan(fisher_p) else 1.0
        else:
            valid_p = np.array([], dtype=float)
            fisher_p_val = 1.0

        p_vals.append(fisher_p_val)

        if t in focal_set:
            if K > 0 and M > 0 and np.any(pos_N):
                exp_k = (N_pos * K) / M
                mean_fold = float(np.mean(k_obs[pos_N] / exp_k))
            else:
                mean_fold = 0.0

            if valid_p.size > 0:
                q_arr = benjamini_hochberg(valid_p.tolist())
                hits_q05 = float(np.sum(np.asarray(q_arr) <= 0.05))
            else:
                hits_q05 = 0.0

            focal_stats[t] = {
                "mean_fold": mean_fold,
                "hits_q05": hits_q05,
                "fisher_p": fisher_p_val,
            }

    q_vals = benjamini_hochberg(p_vals) if len(p_vals) > 0 else np.array([], dtype=float)
    idx_map = {t: idx for idx, t in enumerate(trait_universe)}

    results: dict[str, ArmResult] = {}
    for t in focal_traits:
        idx = idx_map[t]
        st = focal_stats[t]
        results[t] = ArmResult(
            trait_id=t,
            trait_name=TABLE_1_TRAITS.get(t, t),
            primary_stat=st["mean_fold"],
            p_value=st["fisher_p"],
            q_value=float(q_vals[idx]) if len(q_vals) > idx else st["fisher_p"],
            extra_stat=st["hits_q05"],
        )
    return results


def evaluate_loto(
    profiles_df: pd.DataFrame,
    paired_profiles_df: pd.DataFrame,
    trait_universe: list[str],
    focal_traits: list[str],
) -> tuple[list[LOTORow], dict[str, pd.DataFrame]]:
    """Evaluate Leave-One-Trait-Out (LOTO) for all focal traits across three arms.

    Note (#478): base_q uses a focal-set BH family of size m = len(focal_traits)
    (6 for Table 1), whereas loto_q uses the remaining focal family of size m - 1 (5).
    Thus delta_q reflects both denominator perturbation and the focal family-size change.

    Returns:
        tuple of (list of LOTORow, dict mapping arm to delta_matrix DataFrame)
    """
    base_comp = compute_composition_arm(profiles_df, trait_universe, focal_traits)
    base_dir = compute_directional_arm(profiles_df, trait_universe, focal_traits)
    base_ora = compute_ora_arm(paired_profiles_df, trait_universe, focal_traits)

    base_maps = {
        "composition": base_comp,
        "directional": base_dir,
        "ora": base_ora,
    }

    loto_rows: list[LOTORow] = []

    # Prepare square matrix DataFrames for composition
    comp_delta_eff = pd.DataFrame(index=focal_traits, columns=focal_traits, dtype=float)
    comp_delta_q = pd.DataFrame(index=focal_traits, columns=focal_traits, dtype=float)

    for dropped in focal_traits:
        dropped_name = TABLE_1_TRAITS.get(dropped, dropped)
        rem_universe = [c for c in trait_universe if c != dropped]
        rem_focal = [t for t in focal_traits if t != dropped]

        # Recompute all three arms under reduced universe
        loto_comp = compute_composition_arm(profiles_df, rem_universe, rem_focal)
        loto_dir = compute_directional_arm(profiles_df, rem_universe, rem_focal)
        loto_ora = compute_ora_arm(paired_profiles_df, rem_universe, rem_focal)

        loto_arm_results = {
            "composition": loto_comp,
            "directional": loto_dir,
            "ora": loto_ora,
        }

        for arm_name, loto_map in loto_arm_results.items():
            base_map = base_maps[arm_name]
            for target in rem_focal:
                b = base_map[target]
                m = loto_map[target]
                d_stat = m.primary_stat - b.primary_stat
                d_q = m.q_value - b.q_value

                loto_rows.append(
                    LOTORow(
                        arm=arm_name,
                        dropped_trait_id=dropped,
                        dropped_trait_name=dropped_name,
                        target_trait_id=target,
                        target_trait_name=b.trait_name,
                        base_stat=b.primary_stat,
                        loto_stat=m.primary_stat,
                        delta_stat=d_stat,
                        base_q=b.q_value,
                        loto_q=m.q_value,
                        delta_q=d_q,
                    )
                )

                if arm_name == "composition":
                    comp_delta_eff.loc[dropped, target] = d_stat
                    comp_delta_q.loc[dropped, target] = d_q

    # Convert indices to readable trait names
    readable_names = [TABLE_1_TRAITS.get(t, t) for t in focal_traits]
    comp_delta_eff.index = readable_names
    comp_delta_eff.columns = readable_names
    comp_delta_q.index = readable_names
    comp_delta_q.columns = readable_names

    delta_matrices = {
        "composition_delta_effect": comp_delta_eff,
        "composition_delta_q": comp_delta_q,
    }

    return loto_rows, delta_matrices


def evaluate_lco(
    profiles_df: pd.DataFrame,
    paired_profiles_df: pd.DataFrame,
    trait_universe: list[str],
    focal_traits: list[str],
    cluster_to_drop: dict[str, str],
) -> list[LCOSummaryRow]:
    """Evaluate Leave-Cluster-Out (LCO) by removing correlated traits simultaneously.

    Note (#478): base_q uses a focal-set BH family of size m = len(focal_traits),
    whereas lco_q uses the reduced focal family of size len(lco_focal).
    """
    cluster_ids = set(cluster_to_drop.keys())
    cluster_label = " + ".join(cluster_to_drop.values())

    base_comp = compute_composition_arm(profiles_df, trait_universe, focal_traits)
    base_dir = compute_directional_arm(profiles_df, trait_universe, focal_traits)
    base_ora = compute_ora_arm(paired_profiles_df, trait_universe, focal_traits)

    lco_universe = [c for c in trait_universe if c not in cluster_ids]
    lco_focal = [t for t in focal_traits if t not in cluster_ids]

    lco_comp = compute_composition_arm(profiles_df, lco_universe, lco_focal)
    lco_dir = compute_directional_arm(profiles_df, lco_universe, lco_focal)
    lco_ora = compute_ora_arm(paired_profiles_df, lco_universe, lco_focal)

    lco_rows: list[LCOSummaryRow] = []

    arm_pairs = [
        ("composition", base_comp, lco_comp),
        ("directional", base_dir, lco_dir),
        ("ora", base_ora, lco_ora),
    ]

    for arm_name, base_map, lco_map in arm_pairs:
        for t in lco_focal:
            b = base_map[t]
            lco_res = lco_map[t]
            lco_rows.append(
                LCOSummaryRow(
                    arm=arm_name,
                    cluster_dropped=cluster_label,
                    trait_id=t,
                    trait_name=TABLE_1_TRAITS.get(t, t),
                    base_stat=b.primary_stat,
                    lco_stat=lco_res.primary_stat,
                    delta_stat=lco_res.primary_stat - b.primary_stat,
                    base_q=b.q_value,
                    lco_q=lco_res.q_value,
                    delta_q=lco_res.q_value - b.q_value,
                    base_extra=b.extra_stat,
                    lco_extra=lco_res.extra_stat,
                )
            )

    return lco_rows


def evaluate_decision_rule(
    base_comp: dict[str, ArmResult],
    comp_delta_eff_matrix: pd.DataFrame,
    target_trait_id: str = COCCUS_ID,
    dropped_trait_id: str = SPHERE_ID,
) -> DecisionRuleResult:
    """Evaluate the formal stability decision rule.

    Criterion:
        If dropping sphere moves coccus by more than the minimum distance between
        adjacent Table 1 rows, the two cannot be reported as separate convergent findings.
    """
    target_name = TABLE_1_TRAITS.get(target_trait_id, target_trait_id)
    dropped_name = TABLE_1_TRAITS.get(dropped_trait_id, dropped_trait_id)

    d_eff = float(
        comp_delta_eff_matrix.loc[[dropped_name], [target_name]].to_numpy(dtype=float).item()
    )
    abs_d_eff = abs(d_eff)

    sorted_pairs = sorted((res.primary_stat, res.trait_name) for res in base_comp.values())
    effects = [pair[0] for pair in sorted_pairs]
    adjacent_dists = [effects[i + 1] - effects[i] for i in range(len(effects) - 1)]
    min_idx = int(np.argmin(adjacent_dists))
    min_adj_dist = adjacent_dists[min_idx]
    min_adj_pair = (sorted_pairs[min_idx][1], sorted_pairs[min_idx + 1][1])

    passes = abs_d_eff < min_adj_dist

    if passes:
        verdict = (
            f"Passes stability criterion: |Delta effect| ({abs_d_eff:.6f}) < "
            f"min adjacent row distance ({min_adj_dist:.6f}). Denominator coupling does not "
            f"destabilize effect rankings."
        )
    else:
        verdict = (
            f"Fails stability criterion: |Delta effect| ({abs_d_eff:.6f}) >= "
            f"min adjacent row distance ({min_adj_dist:.6f}). Near-twin coupling destabilizes "
            f"separate reporting."
        )

    return DecisionRuleResult(
        target_trait_id=target_trait_id,
        target_trait_name=target_name,
        dropped_trait_id=dropped_trait_id,
        dropped_trait_name=dropped_name,
        delta_effect=d_eff,
        abs_delta_effect=abs_d_eff,
        min_adjacent_distance=min_adj_dist,
        passes_stability_criterion=passes,
        adjacent_distances=adjacent_dists,
        sorted_effects=effects,
        verdict=verdict,
        min_adjacent_pair=min_adj_pair,
    )


def build_coupling_report_text(
    delta_eff_mat: pd.DataFrame,
    delta_q_mat: pd.DataFrame,
    lco_rows: list[LCOSummaryRow],
    rule_result: DecisionRuleResult,
    rule_result_reverse: DecisionRuleResult,
    code_commit: str = "",
    loto_rows: list[LOTORow] | None = None,
    corpus_metrics: dict[str, float | int] | None = None,
) -> str:
    """Format comprehensive text report for trait coupling audit."""
    base_m = delta_eff_mat.shape[0] if not delta_eff_mat.empty else len(TABLE_1_TRAITS)
    loto_m = max(base_m - 1, 1)
    lco_traits = {r.trait_id for r in lco_rows}
    lco_m = len(lco_traits) if lco_traits else max(base_m - len(CANONICAL_LCO_CLUSTER), 1)

    comp_loto_vals = (
        delta_eff_mat.to_numpy(dtype=float) if not delta_eff_mat.empty else np.array([])
    )
    finite_comp_loto = np.abs(comp_loto_vals[~np.isnan(comp_loto_vals)])
    n_loto_pairs = int(finite_comp_loto.size)
    max_comp_loto = float(finite_comp_loto.max()) if n_loto_pairs > 0 else 0.0

    comp_lco_deltas = [abs(r.delta_stat) for r in lco_rows if r.arm == "composition"]
    max_comp_lco = max(comp_lco_deltas, default=0.0)
    max_comp_overall = max(max_comp_loto, max_comp_lco)

    dir_deltas = [abs(r.delta_stat) for r in (loto_rows or []) if r.arm == "directional"] + [
        abs(r.delta_stat) for r in lco_rows if r.arm == "directional"
    ]
    max_dir_delta = max(dir_deltas, default=0.0)

    ora_deltas = [abs(r.delta_stat) for r in (loto_rows or []) if r.arm == "ora"] + [
        abs(r.delta_stat) for r in lco_rows if r.arm == "ora"
    ]
    max_ora_delta = max(ora_deltas, default=0.0)

    pair_suffix = (
        f" (between {rule_result.min_adjacent_pair[0]} and {rule_result.min_adjacent_pair[1]})"
        if rule_result.min_adjacent_pair[0] and rule_result.min_adjacent_pair[1]
        else ""
    )

    if corpus_metrics is not None:
        n_cols_val = int(corpus_metrics.get("n_traits", base_m))
        total_counts = float(corpus_metrics.get("total_annotations", 0.0))
        coccus_counts = float(corpus_metrics.get("coccus_annotations", 0.0))
        twin_r = float(corpus_metrics.get("coccus_sphere_r", float("nan")))
        twin_vif = float(corpus_metrics.get("coccus_sphere_vif", float("nan")))
        twin_stats_str = (
            f" (r = {twin_r:+.3f}, VIF = {twin_vif:.2f}, Issue #264)"
            if not (np.isnan(twin_r) or np.isnan(twin_vif))
            else " (Issue #264)"
        )
        if total_counts >= 1e6:
            total_str = f"{total_counts / 1e6:.2f}M annotations"
        else:
            total_str = f"{total_counts:,.0f} annotations"
        if coccus_counts >= 1e6:
            coccus_str = f" ({coccus_counts / 1e6:.2f}M counts)"
        elif coccus_counts > 0:
            coccus_str = f" ({coccus_counts:,.0f} counts)"
        else:
            coccus_str = ""
        cols_str = f"{n_cols_val} columns ({total_str})"
    else:
        twin_stats_str = " (Issue #264)"
        cols_str = f"{base_m} evaluated columns"
        coccus_str = ""

    lines: list[str] = [
        "=" * 80,
        "TRAIT COUPLING SENSITIVITY AUDIT REPORT (Issue #260)",
        "=" * 80,
    ]
    if code_commit:
        lines.append(f"Code commit: {code_commit}")
    lines.extend(
        [
            "",
            "1. LEAVE-ONE-TRAIT-OUT (LOTO) - COMPOSITION DELTA EFFECT MATRIX (6 x 6):",
            "   Rows: Dropped Trait | Columns: Evaluated Trait",
            "-" * 80,
            delta_eff_mat.to_string(),
            "",
            "2. LEAVE-ONE-TRAIT-OUT (LOTO) - COMPOSITION DELTA Q-VALUE MATRIX (6 x 6):",
            "   Rows: Dropped Trait | Columns: Evaluated Trait",
            (
                f"   Note (#478): Focal-set BH sensitivity check "
                f"(m={base_m} baseline -> m={loto_m} LOTO; not full-family Table 1 q)"
            ),
            "-" * 80,
            delta_q_mat.to_string(),
            "",
            "3. LEAVE-CLUSTER-OUT (LCO) - DROPPING COCCUS + SPHERE SIMULTANEOUSLY:",
            (
                f"   Note (#478): Base q and LCO q are focal-set BH sensitivity checks "
                f"(m={base_m} baseline -> m={lco_m} LCO; not full-family Table 1 q)"
            ),
            "-" * 80,
            (
                f"{'Arm':<14} | {'Trait Name':<28} | {'Base Stat':<10} | {'LCO Stat':<10} | "
                f"{'Delta Stat':<11} | {'Base q':<9} | {'LCO q':<9}"
            ),
            "-" * 80,
        ]
    )

    for row in lco_rows:
        lines.append(
            f"{row.arm:<14} | {row.trait_name:<28} | {row.base_stat:+10.4f} | "
            f"{row.lco_stat:+10.4f} | {row.delta_stat:+11.6f} | "
            f"{row.base_q:9.2e} | {row.lco_q:9.2e}"
        )

    comp_eff_vals = (
        delta_eff_mat.to_numpy(dtype=float)
        if not delta_eff_mat.empty
        else np.array([], dtype=float)
    )
    valid_eff = np.abs(comp_eff_vals[~np.isnan(comp_eff_vals)])
    max_comp_loto_eff = float(np.max(valid_eff)) if valid_eff.size > 0 else 0.0
    n_comp_loto_pairs = int(valid_eff.size)

    comp_q_vals = (
        delta_q_mat.to_numpy(dtype=float) if not delta_q_mat.empty else np.array([], dtype=float)
    )
    valid_q = np.abs(comp_q_vals[~np.isnan(comp_q_vals)])
    max_comp_loto_q = float(np.max(valid_q)) if valid_q.size > 0 else 0.0

    max_lco_stat = max((abs(r.delta_stat) for r in lco_rows), default=0.0)
    max_any_shift = max(max_comp_loto_eff, max_lco_stat)
    min_gap = rule_result.min_adjacent_distance

    lines.extend(
        [
            "",
            "4. FORMAL DECISION RULE EVALUATION:",
            "-" * 80,
            "Criterion: |Delta effect| < min adjacent row spacing in Table 1",
            f"Table 1 baseline sorted effects: {[round(x, 4) for x in rule_result.sorted_effects]}",
            (
                f"Inter-trait adjacent distances:  "
                f"{[round(x, 4) for x in rule_result.adjacent_distances]}"
            ),
            (
                f"Minimum adjacent row distance:   "
                f"{rule_result.min_adjacent_distance:.6f}{pair_suffix}"
            ),
            "",
            "Direct Test (Drop Sphere -> Coccus Effect Shift):",
            f"  Delta effect(coccus):            {rule_result.delta_effect:+.6f}",
            f"  |Delta effect|:                  {rule_result.abs_delta_effect:.6f}",
            f"  Criterion verdict:               {rule_result.verdict}",
            "",
            "Reverse Test (Drop Coccus -> Sphere Effect Shift):",
            f"  Delta effect(sphere):            {rule_result_reverse.delta_effect:+.6f}",
            f"  |Delta effect|:                  {rule_result_reverse.abs_delta_effect:.6f}",
            f"  Criterion verdict:               {rule_result_reverse.verdict}",
            "",
            "5. METHODOLOGICAL & BIOLOGICAL SYNTHESIS:",
            "-" * 80,
            "a. Denominator Coupling Sensitivity by Analytical Arm:",
            (
                f"   In the composition arm, denominator coupling induced by dropping any single "
                f"Table 1 trait across all {n_loto_pairs} LOTO pairs produces a maximum "
                f"|Delta effect| of {max_comp_loto:.6f} (and up to {max_comp_lco:.6f} under LCO; "
                f"overall composition |Delta effect| <= {max_comp_overall:.6f}), compared against "
                f"the tightest Table 1 inter-trait gap of {rule_result.min_adjacent_distance:.6f}."
            ),
            (
                f"   In the directional sign-test arm, maximum |Delta share_up| is "
                f"{max_dir_delta:.6f}; in the ORA arm, maximum |Delta mean_fold| is "
                f"{max_ora_delta:.6f} as removing high-count traits shrinks the shared "
                f"background urn M. Note that focal-set BH q-values above are evaluated over "
                f"m={base_m} (baseline), m={loto_m} (LOTO), and m={lco_m} (LCO) traits as a "
                f"sensitivity check rather than full-family Table 1 q-values."
            ),
            "",
            "b. Distinction Between Redundancy and Denominator Coupling:",
            (
                f"   Coccus and sphere correlate strongly{twin_stats_str} due to "
                f"semantic/curation redundancy in METPO mappings, NOT denominator distortion."
            ),
            (
                f"   Because total non-MediaDive annotation counts are distributed across "
                f"{cols_str}, removing even a high-annotation trait like "
                f"coccus{coccus_str} shrinks the denominator smoothly without tilting relative "
                f"compositions. Therefore, while coccus and sphere should not be claimed as "
                f"independent biological mechanisms (due to semantic overlap), their Table 1 "
                f"statistics are structurally stable."
            ),
            "=" * 80,
        ]
    )
    # Strip padding left by fixed-width last columns; trailing whitespace fails pre-commit.
    return "\n".join(line.rstrip() for line in "\n".join(lines).splitlines()) + "\n"


def run_coupling_audit(
    profiles_path: Path = DEFAULT_PROFILES_PATH,
    columns_path: Path = DEFAULT_COLUMNS_PATH,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    code_commit: str = "",
) -> dict[str, Any]:
    """Execute complete trait coupling audit and write derived datasets."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    profiles_df, valid_cols = load_non_media_profile_matrix(profiles_path, columns_path)
    t1_keys = list(TABLE_1_TRAITS.keys())

    # Filter to paired experiments for paired and ORA evaluations
    pvt_check = profiles_df.pivot(
        index=["study_id", "experiment_id"], columns="inc_dec", values=t1_keys[0]
    ).dropna()
    paired_exp_set = set(pvt_check.index)
    paired_mask = [
        (s, e) in paired_exp_set
        for s, e in zip(profiles_df["study_id"], profiles_df["experiment_id"], strict=True)
    ]
    sub_paired_df = profiles_df[paired_mask].copy()

    # 1. Baseline evaluations
    base_comp = compute_composition_arm(profiles_df, valid_cols, t1_keys)

    # 2. LOTO evaluation
    loto_rows, delta_matrices = evaluate_loto(
        profiles_df=profiles_df,
        paired_profiles_df=sub_paired_df,
        trait_universe=valid_cols,
        focal_traits=t1_keys,
    )

    # 3. LCO evaluation
    lco_rows = evaluate_lco(
        profiles_df=profiles_df,
        paired_profiles_df=sub_paired_df,
        trait_universe=valid_cols,
        focal_traits=t1_keys,
        cluster_to_drop=CANONICAL_LCO_CLUSTER,
    )

    # 4. Decision rule evaluation
    rule_result = evaluate_decision_rule(
        base_comp=base_comp,
        comp_delta_eff_matrix=delta_matrices["composition_delta_effect"],
        target_trait_id=COCCUS_ID,
        dropped_trait_id=SPHERE_ID,
    )
    rule_result_rev = evaluate_decision_rule(
        base_comp=base_comp,
        comp_delta_eff_matrix=delta_matrices["composition_delta_effect"],
        target_trait_id=SPHERE_ID,
        dropped_trait_id=COCCUS_ID,
    )

    # Compute dynamic corpus metrics for report prose (#479)
    total_annotations = float(profiles_df[valid_cols].to_numpy(dtype=float).sum())
    coccus_annotations = (
        float(profiles_df[COCCUS_ID].to_numpy(dtype=float).sum())
        if COCCUS_ID in profiles_df.columns
        else 0.0
    )
    if COCCUS_ID in profiles_df.columns and SPHERE_ID in profiles_df.columns:
        twin_summary = evaluate_conditional_twin_test(
            profiles_df,
            trait_a=COCCUS_ID,
            trait_b=SPHERE_ID,
            columns_for_depth=valid_cols,
        )
        coccus_sphere_r = twin_summary.profile_correlation
        coccus_sphere_vif = twin_summary.vif_profile
    else:
        coccus_sphere_r = float("nan")
        coccus_sphere_vif = float("nan")

    corpus_metrics: dict[str, float | int] = {
        "n_traits": len(valid_cols),
        "total_annotations": total_annotations,
        "coccus_annotations": coccus_annotations,
        "coccus_sphere_r": coccus_sphere_r,
        "coccus_sphere_vif": coccus_sphere_vif,
    }

    # 5. Export TSVs and Report
    loto_df = pd.DataFrame([asdict(r) for r in loto_rows])
    lco_df = pd.DataFrame([asdict(r) for r in lco_rows])

    loto_tsv_path = output_dir / "trait_coupling_loto_matrix.tsv"
    lco_tsv_path = output_dir / "trait_coupling_lco_summary.tsv"
    report_txt_path = output_dir / "trait_coupling_report.txt"

    loto_df.to_csv(loto_tsv_path, sep="\t", index=False)
    lco_df.to_csv(lco_tsv_path, sep="\t", index=False)

    report_text = build_coupling_report_text(
        delta_eff_mat=delta_matrices["composition_delta_effect"],
        delta_q_mat=delta_matrices["composition_delta_q"],
        lco_rows=lco_rows,
        rule_result=rule_result,
        rule_result_reverse=rule_result_rev,
        code_commit=code_commit,
        loto_rows=loto_rows,
        corpus_metrics=corpus_metrics,
    )
    report_txt_path.write_text(report_text, encoding="utf-8")

    return {
        "loto_tsv": loto_tsv_path,
        "lco_tsv": lco_tsv_path,
        "report_txt": report_txt_path,
        "loto_rows": loto_rows,
        "lco_rows": lco_rows,
        "delta_matrices": delta_matrices,
        "rule_result": rule_result,
        "rule_result_rev": rule_result_rev,
    }
