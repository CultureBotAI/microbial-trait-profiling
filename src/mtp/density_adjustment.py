"""Trait correlation, annotation density adjustment, and conditional near-twin evaluation.

Addresses Issue #264, Issue #352, and follow-ups #375-#383:
1. Quantifies the baseline trait correlation structure across non-MediaDive traits:
   - Pairwise correlation distribution (median, mean, proportion > 0.5, proportion > 0.9).
   - Dominant annotation density axis: PC1 and top principal components variance explained.
   - Effective number of independent tests: Li & Ji (2005), Galwey (2009), Cheverud/Nyholt (2004).
     (Note: Cheverud/Nyholt is reported for completeness but excluded as unreliable per #325, #357).
2. Adjusts for the curation density axis:
   - Method A: Profile-level PC1 regression / residualization.
   - Method B: Profile composition normalization (Total Sum Scaling / profile
     annotation depth normalization).
   - Re-evaluates pairwise correlations and tests the diagnostic test pair
     (facultatively aerobic x moderately halophilic) across profile representations.
3. Evaluates conditional test and collinearity for near-twins (coccus shaped vs. sphere shaped):
   - Computes bivariate correlation, paired contrast correlation, and
     Variance Inflation Factor (VIF).
   - Evaluates partial correlation and paired contrast regression.
"""

from __future__ import annotations

import argparse
import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pandas as pd
from scipy import stats

from mtp.composition_denominators import MEDIADIVE_MARKER, trait_fractions
from mtp.paths import DERIVED_DATA_DIR
from mtp.trait_key_policy import display_trait_label, read_recorded_member_ids
from mtp.trait_lists import (
    DeclaredTraitList,
    file_sha256,
    load_trait_list,
    output_dir_argument,
    trait_labels,
)
from mtp.trait_screen import effective_family_size

DEFAULT_PROFILES_PATH = DERIVED_DATA_DIR / "study_trait_profiles_wide.tsv"
DEFAULT_COLUMNS_PATH = DERIVED_DATA_DIR / "taxa_trait_matrix_columns.tsv"
DEFAULT_OUTPUT_DIR = DERIVED_DATA_DIR / "density_adjustment"
TRAIT_LIST_PROVENANCE_FILENAME = "trait_list_provenance.json"

# The six canonical Table 1 traits
TABLE_1_TRAITS: dict[str, str] = {
    "METPO:1000608": "facultatively aerobic",
    "METPO:2000103|GO:0030245": "cellulose catabolic process",
    "METPO:1000668": "coccus shaped",
    "METPO:1000683": "sphere shaped",
    "METPO:2000103|pathways:arsenate_reduction": "arsenate reduction",
    "METPO:1000623": "moderately halophilic",
}


def table1_traits_for_columns(columns_path: Path | None) -> dict[str, str]:
    """Return Table 1 trait labels resolved against recorded column member metadata."""
    members = read_recorded_member_ids(columns_path)
    return {
        trait_id: display_trait_label(trait_id, label, members.get(trait_id))
        for trait_id, label in TABLE_1_TRAITS.items()
    }


# The specific test case pair identified in Issue #264
TEST_PAIR_A = "METPO:1000608"  # facultatively aerobic
TEST_PAIR_B = "METPO:1000623"  # moderately halophilic

# The near-twin candidate pair
TWIN_A = "METPO:1000668"  # coccus shaped
TWIN_B = "METPO:1000683"  # sphere shaped


@dataclass(frozen=True)
class CorrelationMetrics:
    """Summary of correlation and dimensionality metrics across a trait matrix."""

    n_samples: int
    n_traits: int
    constant_traits_excluded: int
    median_pairwise_r: float
    mean_pairwise_r: float
    pct_pairs_gt_05: float
    pct_pairs_gt_09: float
    pc1_variance_pct: float
    pc1_5_variance_pct: float
    m_eff_li_ji: float
    m_eff_galwey: float
    m_eff_nyholt: float
    test_pair_r: float
    near_twin_r: float


@dataclass(frozen=True)
class ConditionalTwinSummary:
    """Collinearity and conditional contrast summary for near-twin traits."""

    trait_a: str
    trait_b: str
    trait_a_label: str
    trait_b_label: str
    profile_correlation: float
    profile_partial_r_depth: float
    vif_profile: float
    profile_slope_a_on_b: float
    profile_intercept_a_on_b: float
    n_paired_experiments: int
    paired_correlation: float
    vif_paired: float
    paired_slope_a_on_b: float
    paired_slope_pvalue: float


@dataclass(frozen=True)
class Table1TraitEffect:
    """Paired contrast effect for one Table 1 trait before and after density adjustment.

    ``relative_change_pct`` is the signed change (adj - raw) / |raw| x 100, kept for
    continuity. It reads as attenuation only for a positive raw mean: for a negative
    one, a negative value means the effect grew. ``magnitude_change_pct`` is
    (|adj| - |raw|) / |raw| x 100, negative exactly when the effect shrank, and
    ``sign_flipped`` marks an adjusted mean of the opposite sign (#1084). A mean
    within 1e-9 of zero has no sign, and a raw one leaves the magnitude change NaN.
    """

    trait_id: str
    trait_label: str
    raw_mean_diff: float
    raw_wilcoxon_p: float
    adj_mean_diff: float
    adj_wilcoxon_p: float
    relative_change_pct: float
    status: str
    magnitude_change_pct: float
    sign_flipped: bool


# Statuses whose raw and adjusted paired contrasts are both significant (Wilcoxon p <= 0.05)
# in the same direction (#1084, #1097, #1098); deck figure 4 draws only these as surviving.
SIGNIFICANT_POST_ADJUSTMENT_STATUSES = frozenset(
    {"enhanced (significant, larger magnitude)", "survived (significant post-adjustment)"}
)

# A paired mean difference within this of zero has no sign (#1097).
_ZERO_MEAN = 1e-9


def _sign_flipped(raw_mean: float, adj_mean: float) -> bool:
    """Whether two means clearly away from zero have opposite signs."""
    return min(abs(raw_mean), abs(adj_mean)) > _ZERO_MEAN and raw_mean * adj_mean < 0


def adjustment_status(raw_mean: float, raw_p: float, adj_mean: float, adj_p: float) -> str:
    """Classify a trait by raw and adjusted significance and by the change in effect magnitude.

    Significance comes from the Wilcoxon p-values; "enhanced" and "attenuated" compare
    |adj| with |raw| (#1084). A reversed sign and a raw mean of zero, where the
    magnitude change is undefined, each get their own status (#1097), and a significant
    adjusted contrast whose raw contrast was nonsignificant (or had <= 10 non-zero pairs)
    is classified as "significant only post-adjustment" rather than "survived" or
    "enhanced" (#1098).
    """
    if math.isnan(adj_p):
        return "insufficient non-zero pairs"
    if abs(raw_mean) <= _ZERO_MEAN:
        if adj_p <= 0.05:
            return "undefined magnitude change (raw mean zero, significant post-adjustment)"
        return "undefined magnitude change (raw mean zero, nonsignificant post-adjustment)"
    if _sign_flipped(raw_mean, adj_mean):
        if adj_p <= 0.05:
            return "reversed (significant post-adjustment)"
        return "reversed (nonsignificant post-adjustment)"
    raw_significant = not math.isnan(raw_p) and raw_p <= 0.05
    larger = abs(adj_mean) > abs(raw_mean)
    if adj_p <= 0.05:
        if not raw_significant:
            return "significant only post-adjustment"
        if larger:
            return "enhanced (significant, larger magnitude)"
        return "survived (significant post-adjustment)"
    if raw_significant:
        if larger:
            return "collapsed (nonsignificant post-adjustment, larger magnitude)"
        return "collapsed (attenuated to nonsignificant)"
    if larger:
        return "enhanced magnitude (nonsignificant post-adjustment)"
    return "attenuated (nonsignificant post-adjustment)"


def read_non_media_trait_columns(columns_path: Path) -> pd.DataFrame:
    """Return the columns-table rows of the non-MediaDive traits, in table order.

    Non-MediaDive traits are defined as columns where object_category != 'METPO:1004005'.
    """
    if not columns_path.exists():
        gz_candidate = columns_path.with_name(f"{columns_path.name}.gz")
        if gz_candidate.exists():
            columns_path = gz_candidate
        else:
            raise FileNotFoundError(f"Columns file not found: {columns_path} (or {gz_candidate})")

    df = pd.read_csv(columns_path, sep="\t")
    if "object_category" not in df.columns or "column_id" not in df.columns:
        raise ValueError(f"Columns file missing required columns: {columns_path}")
    return df[
        (df["object_category"] != "METPO:1004005")
        & (~df["column_id"].astype(str).str.lower().str.contains(MEDIADIVE_MARKER, regex=False))
    ]


def get_non_media_trait_columns(columns_path: Path) -> list[str]:
    """Return ordered list of non-MediaDive trait column identifiers."""
    return [str(c) for c in read_non_media_trait_columns(columns_path)["column_id"].tolist()]


def load_non_media_profile_matrix(
    profiles_path: Path,
    columns_path: Path,
) -> tuple[pd.DataFrame, list[str]]:
    """Load wide study trait profiles restricted to the non-MediaDive traits."""
    if not profiles_path.exists():
        gz_candidate = profiles_path.with_name(f"{profiles_path.name}.gz")
        if gz_candidate.exists():
            profiles_path = gz_candidate
        else:
            raise FileNotFoundError(f"Profiles file not found: {profiles_path} (or {gz_candidate})")

    if not columns_path.exists():
        gz_candidate = columns_path.with_name(f"{columns_path.name}.gz")
        if gz_candidate.exists():
            columns_path = gz_candidate
        else:
            raise FileNotFoundError(f"Columns file not found: {columns_path} (or {gz_candidate})")

    non_media_cols = get_non_media_trait_columns(columns_path)
    profiles_df = pd.read_csv(profiles_path, sep="\t")

    valid_cols = [c for c in non_media_cols if c in profiles_df.columns]
    if len(valid_cols) != len(non_media_cols):
        missing = set(non_media_cols) - set(valid_cols)
        raise ValueError(f"Profiles table missing {len(missing)} expected trait columns: {missing}")

    return profiles_df, valid_cols


def compute_effective_tests(
    eigenvalues: Sequence[float] | np.ndarray,
) -> tuple[float, float, float]:
    """Compute effective number of tests: Li & Ji (2005), Galwey (2009), Cheverud/Nyholt.

    Delegates to canonical trait_screen.effective_family_size (#378).
    Handles edge cases:
      - p == 0: returns (0.0, 0.0, 0.0)
      - p == 1: returns (min(abs(val), 1.0), min(abs(val), 1.0), 1.0)
      - p >= 2: calls effective_family_size and caps at p (#377).

    Args:
        eigenvalues: 1D array or sequence of eigenvalues of the correlation matrix.

    Returns:
        tuple of (M_eff_li_ji, M_eff_galwey, M_eff_nyholt)
    """
    p = len(eigenvalues)
    if p == 0:
        return 0.0, 0.0, 0.0
    if p < 2:
        val = float(abs(eigenvalues[0]))
        return min(val, 1.0), min(val, 1.0), 1.0

    eval_list = [float(x) for x in eigenvalues]
    if not any(x > 0.0 for x in eval_list):
        return 0.0, 0.0, float(p)

    li_ji, galwey, nyholt = effective_family_size(eval_list)
    return min(li_ji, float(p)), min(galwey, float(p)), min(max(nyholt, 1.0), float(p))


def compute_correlation_metrics(
    trait_matrix: np.ndarray,
    column_ids: Sequence[str],
    *,
    method: Literal["pearson", "spearman"] = "pearson",
) -> CorrelationMetrics:
    """Compute pairwise correlation distribution, PC1 variance, and effective test numbers.

    Constant columns (zero variance) are excluded before correlation estimation (#378),
    preventing spurious eigenvalues of 1.0 from artificially inflating effective test counts.
    """
    if trait_matrix.ndim != 2:
        raise ValueError(f"trait_matrix must be 2D, got shape {trait_matrix.shape}")
    n_samples, n_traits = trait_matrix.shape
    if n_traits < 2:
        raise ValueError("Need at least 2 traits to compute correlation metrics.")
    if len(column_ids) != n_traits:
        raise ValueError(
            f"column_ids length ({len(column_ids)}) does not match "
            f"trait_matrix columns ({n_traits})"
        )

    # Drop constant columns before correlating (#378)
    stds = np.std(trait_matrix, axis=0)
    varying_mask = stds > 0.0
    constant_traits_excluded = int(np.sum(~varying_mask))
    varying_indices = np.where(varying_mask)[0]

    if len(varying_indices) < 2:
        return CorrelationMetrics(
            n_samples=n_samples,
            n_traits=n_traits,
            constant_traits_excluded=constant_traits_excluded,
            median_pairwise_r=0.0,
            mean_pairwise_r=0.0,
            pct_pairs_gt_05=0.0,
            pct_pairs_gt_09=0.0,
            pc1_variance_pct=0.0,
            pc1_5_variance_pct=0.0,
            m_eff_li_ji=float(len(varying_indices)),
            m_eff_galwey=float(len(varying_indices)),
            m_eff_nyholt=float(len(varying_indices)),
            test_pair_r=0.0,
            near_twin_r=float("nan"),
        )

    reduced_matrix = trait_matrix[:, varying_indices]
    reduced_cols = [column_ids[i] for i in varying_indices]
    n_varying = len(varying_indices)

    if method == "pearson":
        corr = np.corrcoef(reduced_matrix.T)
    elif method == "spearman":
        res = stats.spearmanr(reduced_matrix)
        stat = getattr(res, "statistic", None)
        if stat is None:
            stat = res[0]
        if n_varying == 2 or np.ndim(stat) == 0:
            s_val = float(stat)
            corr = np.array([[1.0, s_val], [s_val, 1.0]], dtype=float)
        else:
            corr = np.asarray(stat, dtype=float)
    else:
        raise ValueError(f"Unsupported correlation method: {method}")

    # Ensure symmetric and non-nan
    corr = np.nan_to_num(corr, nan=0.0)
    np.fill_diagonal(corr, 1.0)

    iu = np.triu_indices(n_varying, k=1)
    upper_vals = corr[iu]

    median_r = float(np.median(upper_vals))
    mean_r = float(np.mean(upper_vals))
    pct_05 = float(np.mean(upper_vals > 0.5) * 100.0)
    pct_09 = float(np.mean(upper_vals > 0.9) * 100.0)

    # Eigenvalues of the correlation matrix
    evals = np.sort(np.linalg.eigvalsh(corr))[::-1]
    evals_sum = float(np.sum(evals))
    pc1_pct = float(evals[0] / evals_sum * 100.0) if evals_sum > 0 else 0.0
    pc1_5_pct = float(np.sum(evals[:5]) / evals_sum * 100.0) if evals_sum > 0 else 0.0

    li_ji, galwey, nyholt = compute_effective_tests(evals)

    # Specific pair correlations
    test_pair_r = 0.0
    if TEST_PAIR_A in reduced_cols and TEST_PAIR_B in reduced_cols:
        idx_a = reduced_cols.index(TEST_PAIR_A)
        idx_b = reduced_cols.index(TEST_PAIR_B)
        test_pair_r = float(corr[idx_a, idx_b])

    near_twin_r = float("nan")
    if TWIN_A in reduced_cols and TWIN_B in reduced_cols:
        idx_twin_a = reduced_cols.index(TWIN_A)
        idx_twin_b = reduced_cols.index(TWIN_B)
        near_twin_r = float(corr[idx_twin_a, idx_twin_b])

    return CorrelationMetrics(
        n_samples=n_samples,
        n_traits=n_traits,
        constant_traits_excluded=constant_traits_excluded,
        median_pairwise_r=median_r,
        mean_pairwise_r=mean_r,
        pct_pairs_gt_05=pct_05,
        pct_pairs_gt_09=pct_09,
        pc1_variance_pct=pc1_pct,
        pc1_5_variance_pct=pc1_5_pct,
        m_eff_li_ji=li_ji,
        m_eff_galwey=galwey,
        m_eff_nyholt=nyholt,
        test_pair_r=test_pair_r,
        near_twin_r=near_twin_r,
    )


def regress_out_pc1(trait_matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Regress out the first principal component from each trait column.

    Args:
        trait_matrix: 2D array of shape (N, P) containing trait profile counts.

    Returns:
        tuple of (adjusted_matrix, pc1_scores, pc1_variance_pct)
    """
    if trait_matrix.ndim != 2 or trait_matrix.shape[0] < 2 or trait_matrix.shape[1] < 1:
        return (
            trait_matrix.astype(float),
            np.zeros(trait_matrix.shape[0], dtype=float),
            0.0,
        )

    x = np.nan_to_num(trait_matrix.astype(float), nan=0.0)
    mu = x.mean(axis=0)
    std = x.std(axis=0, ddof=1)
    std = np.nan_to_num(std, nan=1.0)
    std[std == 0.0] = 1.0
    z = (x - mu) / std
    z = np.nan_to_num(z, nan=0.0)

    # SVD on standardized matrix
    u, s, _ = np.linalg.svd(z, full_matrices=False)
    if len(s) == 0:
        return x.copy(), np.zeros(x.shape[0], dtype=float), 0.0

    pc1_scores = u[:, 0] * s[0]
    total_var = float(np.sum(s**2))
    pc1_var_pct = float((s[0] ** 2) / total_var * 100.0) if total_var > 0 else 0.0

    # Center PC1 scores for OLS regression
    pc1_c = pc1_scores - pc1_scores.mean()
    pc1_denom = float(np.dot(pc1_c, pc1_c))

    if pc1_denom > 0:
        betas = np.dot(x.T, pc1_c) / pc1_denom
        residuals = x - np.outer(pc1_c, betas)
    else:
        residuals = x.copy()

    return residuals, pc1_scores, pc1_var_pct


def normalize_traits_by_profile_composition(
    df: pd.DataFrame,
    cols: Sequence[str] | pd.Index[str],
) -> np.ndarray:
    """Normalize trait counts per profile by total annotation depth (profile composition / TSS).

    Method B adjusts for overall variation in profile annotation richness by dividing
    each profile's trait counts by its row sum across the candidate trait columns (#352, #379):
        x_norm = x_i / sum_j(x_ij)
    This accounts for profiles with higher total numbers of annotated traits without
    conflating taxon count with sequencing depth.
    """
    cols_list = list(cols)
    if not cols_list:
        raise ValueError("cols must not be empty")

    missing = [c for c in cols_list if c not in df.columns]
    if missing:
        raise ValueError(f"DataFrame missing {len(missing)} requested columns: {missing[:5]}")

    fractions = trait_fractions(df[cols_list])
    return fractions.to_numpy(dtype=float)


def evaluate_table1_paired_effects(
    profiles_df: pd.DataFrame,
    adj_matrix: np.ndarray,
    non_media_cols: list[str],
    traits: Mapping[str, str] = TABLE_1_TRAITS,
) -> list[Table1TraitEffect]:
    """Evaluate paired experiment contrasts (inc - dec) before and after PC1 adjustment.

    ``traits`` maps the reported trait identifiers to labels, in report order. PC1 is
    fit on every non-MediaDive trait whatever is reported, so a trait's row does not
    depend on which other traits are listed (#1054).
    """
    t1_keys = [k for k in traits if k in non_media_cols]
    has_paired_directions = (
        not profiles_df.empty
        and "study_id" in profiles_df.columns
        and "experiment_id" in profiles_df.columns
        and "inc_dec" in profiles_df.columns
        and "inc" in profiles_df["inc_dec"].values
        and "dec" in profiles_df["inc_dec"].values
        and len(t1_keys) > 0
        and adj_matrix.ndim == 2
        and adj_matrix.shape[0] == len(profiles_df)
        and adj_matrix.shape[1] == len(non_media_cols)
    )

    pvt_raw = pd.DataFrame()
    pvt_adj = pd.DataFrame()

    if has_paired_directions:
        adj_df = pd.DataFrame(adj_matrix, columns=non_media_cols)
        adj_df["study_id"] = profiles_df["study_id"].values
        adj_df["experiment_id"] = profiles_df["experiment_id"].values
        adj_df["inc_dec"] = profiles_df["inc_dec"].values

        pvt_raw = profiles_df.pivot_table(
            index=["study_id", "experiment_id"],
            columns="inc_dec",
            values=t1_keys,
            aggfunc="mean",
        ).dropna()

        pvt_adj = adj_df.pivot_table(
            index=["study_id", "experiment_id"],
            columns="inc_dec",
            values=t1_keys,
            aggfunc="mean",
        ).dropna()

    results: list[Table1TraitEffect] = []
    for trait_id in t1_keys:
        label = traits[trait_id]

        can_eval = (
            has_paired_directions
            and not pvt_raw.empty
            and not pvt_adj.empty
            and (trait_id, "inc") in pvt_raw.columns
            and (trait_id, "dec") in pvt_raw.columns
            and (trait_id, "inc") in pvt_adj.columns
            and (trait_id, "dec") in pvt_adj.columns
        )

        if can_eval:
            raw_diff = (pvt_raw[(trait_id, "inc")] - pvt_raw[(trait_id, "dec")]).to_numpy(
                dtype=float
            )
            adj_diff = (pvt_adj[(trait_id, "inc")] - pvt_adj[(trait_id, "dec")]).to_numpy(
                dtype=float
            )

            raw_mean = float(np.mean(raw_diff)) if len(raw_diff) > 0 else 0.0
            adj_mean = float(np.mean(adj_diff)) if len(adj_diff) > 0 else 0.0

            # Wilcoxon signed-rank test on finite non-zero differences
            nz_raw = raw_diff[np.isfinite(raw_diff) & (raw_diff != 0.0)]
            nz_adj = adj_diff[np.isfinite(adj_diff) & (adj_diff != 0.0)]

            # Use NaN when <= 10 non-zero pairs instead of fabricating 1.0 (#383)
            raw_p = float(stats.wilcoxon(nz_raw).pvalue) if len(nz_raw) > 10 else float("nan")
            adj_p = float(stats.wilcoxon(nz_adj).pvalue) if len(nz_adj) > 10 else float("nan")

            rel_change = (
                float((adj_mean - raw_mean) / abs(raw_mean) * 100.0)
                if abs(raw_mean) > 1e-9
                else 0.0
            )
        else:
            raw_mean = 0.0
            adj_mean = 0.0
            raw_p = float("nan")
            adj_p = float("nan")
            rel_change = 0.0

        magnitude_change = (
            float((abs(adj_mean) - abs(raw_mean)) / abs(raw_mean) * 100.0)
            if abs(raw_mean) > _ZERO_MEAN
            else float("nan")
        )

        results.append(
            Table1TraitEffect(
                trait_id=trait_id,
                trait_label=label,
                raw_mean_diff=raw_mean,
                raw_wilcoxon_p=raw_p,
                adj_mean_diff=adj_mean,
                adj_wilcoxon_p=adj_p,
                relative_change_pct=rel_change,
                # General statistical criteria rather than hardcoded trait IDs (#383)
                status=adjustment_status(raw_mean, raw_p, adj_mean, adj_p),
                magnitude_change_pct=magnitude_change,
                sign_flipped=_sign_flipped(raw_mean, adj_mean),
            )
        )

    return results


def evaluate_conditional_twin_test(
    profiles_df: pd.DataFrame,
    trait_a: str = TWIN_A,
    trait_b: str = TWIN_B,
    columns_for_depth: Sequence[str] | None = None,
) -> ConditionalTwinSummary:
    """Evaluate collinearity, regression, and partial correlation for near-twins (#380).

    Computes:
      - Bivariate profile correlation r, R^2, and Variance Inflation Factor (VIF).
      - Partial correlation controlling for profile total annotation depth.
      - Paired experiment contrast correlation, paired VIF, and paired OLS regression.
    """
    if trait_a not in profiles_df.columns or trait_b not in profiles_df.columns:
        raise ValueError(f"Missing twin columns: {trait_a} or {trait_b}")

    vec_a = profiles_df[trait_a].to_numpy(dtype=float)
    vec_b = profiles_df[trait_b].to_numpy(dtype=float)

    std_a = float(np.std(vec_a)) if len(vec_a) > 0 else 0.0
    std_b = float(np.std(vec_b)) if len(vec_b) > 0 else 0.0

    # Profile-level bivariate correlation
    if len(vec_a) >= 2 and std_a > 0.0 and std_b > 0.0:
        corr_prof = float(np.corrcoef(vec_a, vec_b)[0, 1])
        r2_prof = corr_prof**2
        vif_prof = float(1.0 / (1.0 - r2_prof)) if (1.0 - r2_prof) > 1e-12 else float("inf")
    else:
        corr_prof = float("nan")
        r2_prof = float("nan")
        vif_prof = float("nan")

    # Partial correlation controlling for profile annotation depth
    if columns_for_depth:
        valid_depth_cols = [c for c in columns_for_depth if c in profiles_df.columns]
        if valid_depth_cols:
            depth = profiles_df[valid_depth_cols].sum(axis=1).to_numpy(dtype=float)
        else:
            depth = vec_a + vec_b
    else:
        depth = vec_a + vec_b

    std_depth = float(np.std(depth)) if len(depth) > 0 else 0.0
    if (
        len(vec_a) >= 2
        and std_a > 0.0
        and std_b > 0.0
        and std_depth > 0.0
        and not math.isnan(corr_prof)
    ):
        r_az = float(np.corrcoef(vec_a, depth)[0, 1])
        r_bz = float(np.corrcoef(vec_b, depth)[0, 1])
        denom = math.sqrt(max(0.0, (1.0 - r_az**2) * (1.0 - r_bz**2)))
        partial_r = (
            float(np.clip((corr_prof - r_az * r_bz) / denom, -1.0, 1.0)) if denom > 1e-12 else 0.0
        )
    else:
        partial_r = float("nan")

    # Profile OLS regression: A ~ B
    if len(vec_b) >= 2 and std_b > 0.0:
        res_prof = stats.linregress(vec_b, vec_a)
        slope_prof = float(res_prof.slope)
        intercept_prof = float(res_prof.intercept)
    else:
        slope_prof = 0.0
        intercept_prof = float(np.mean(vec_a)) if len(vec_a) > 0 else 0.0

    # Paired experiment contrasts
    has_paired_directions = (
        not profiles_df.empty
        and "study_id" in profiles_df.columns
        and "experiment_id" in profiles_df.columns
        and "inc_dec" in profiles_df.columns
        and "inc" in profiles_df["inc_dec"].values
        and "dec" in profiles_df["inc_dec"].values
    )
    if has_paired_directions:
        pvt = profiles_df.pivot_table(
            index=["study_id", "experiment_id"],
            columns="inc_dec",
            values=[trait_a, trait_b],
            aggfunc="mean",
        ).dropna()
    else:
        pvt = pd.DataFrame()

    can_eval_paired = (
        not pvt.empty
        and (trait_a, "inc") in pvt.columns
        and (trait_a, "dec") in pvt.columns
        and (trait_b, "inc") in pvt.columns
        and (trait_b, "dec") in pvt.columns
    )

    if can_eval_paired:
        delta_a = (pvt[(trait_a, "inc")] - pvt[(trait_a, "dec")]).to_numpy(dtype=float)
        delta_b = (pvt[(trait_b, "inc")] - pvt[(trait_b, "dec")]).to_numpy(dtype=float)
    else:
        delta_a = np.array([], dtype=float)
        delta_b = np.array([], dtype=float)

    n_paired = len(delta_a)
    std_delta_a = float(np.std(delta_a)) if n_paired > 0 else 0.0
    std_delta_b = float(np.std(delta_b)) if n_paired > 0 else 0.0

    if n_paired >= 2 and std_delta_a > 0.0 and std_delta_b > 0.0:
        corr_paired = float(np.corrcoef(delta_a, delta_b)[0, 1])
        r2_paired = corr_paired**2
        vif_paired = float(1.0 / (1.0 - r2_paired)) if (1.0 - r2_paired) > 1e-12 else float("inf")
    else:
        corr_paired = float("nan")
        r2_paired = float("nan")
        vif_paired = float("nan")

    # Paired regression: delta_A ~ delta_B
    if n_paired >= 2 and std_delta_b > 0.0:
        res_paired = stats.linregress(delta_b, delta_a)
        slope_paired = float(res_paired.slope)
        pval_paired_slope = float(res_paired.pvalue)
    else:
        slope_paired = 0.0
        pval_paired_slope = float("nan")

    return ConditionalTwinSummary(
        trait_a=trait_a,
        trait_b=trait_b,
        trait_a_label=TABLE_1_TRAITS.get(trait_a, trait_a),
        trait_b_label=TABLE_1_TRAITS.get(trait_b, trait_b),
        profile_correlation=corr_prof,
        profile_partial_r_depth=partial_r,
        vif_profile=vif_prof,
        profile_slope_a_on_b=slope_prof,
        profile_intercept_a_on_b=intercept_prof,
        n_paired_experiments=n_paired,
        paired_correlation=corr_paired,
        vif_paired=vif_paired,
        paired_slope_a_on_b=slope_paired,
        paired_slope_pvalue=pval_paired_slope,
    )


def declared_density_traits(
    columns_path: Path,
    trait_list: DeclaredTraitList,
) -> dict[str, str]:
    """Label a declared list; every trait must be a non-MediaDive column (#1054)."""
    non_media = read_non_media_trait_columns(columns_path)
    non_media_ids = set(non_media["column_id"].astype(str))
    excluded = [trait_id for trait_id in trait_list.trait_ids if trait_id not in non_media_ids]
    if excluded:
        raise ValueError(f"Declared traits are not non-MediaDive trait columns: {excluded}")
    return trait_labels(non_media, trait_list.trait_ids)


def write_trait_list_provenance(
    output_dir: Path,
    trait_list: DeclaredTraitList,
    *,
    profiles_path: Path,
    columns_path: Path,
    code_commit: str,
) -> Path:
    """Record the declared list, its hash and the inputs it was evaluated on."""
    # The loaders above fall back to a .gz sibling; hash the file they actually read.
    resolved_profiles, resolved_columns = (
        path if path.exists() else path.with_name(f"{path.name}.gz")
        for path in (profiles_path, columns_path)
    )
    provenance = {
        **trait_list.provenance(),
        "reported_table": "table1_density_adjusted_effects.tsv",
        "pc1_fit_traits": "every non-MediaDive trait column, independent of the list",
        "profiles_file": resolved_profiles.name,
        "profiles_sha256": file_sha256(resolved_profiles),
        "columns_file": resolved_columns.name,
        "columns_sha256": file_sha256(resolved_columns),
        "code_commit": code_commit,
    }
    path = output_dir / TRAIT_LIST_PROVENANCE_FILENAME
    path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    return path


def reported_traits_heading(trait_list: DeclaredTraitList | None) -> str:
    """Name the reported trait set in report headings; the default text is unchanged."""
    if trait_list is None:
        return "SURVIVAL OF CONVERGENT SIX TRAIT SIGNALS"
    return f"DECLARED TRAIT LIST SIGNALS ({len(trait_list.trait_ids)} TRAITS)"


def _has_density_directory_component(path: Path) -> bool:
    """Inspect every symlink target, including reserved names hidden between aliases."""
    pending = list(path.absolute().parts)
    current = Path(pending.pop(0))
    seen: set[tuple[Path, tuple[str, ...]]] = set()
    while pending:
        component = pending.pop(0)
        if component.casefold() == "density_adjustment":
            return True
        if component == "..":
            current = current.parent
            continue
        current /= component
        if current.is_symlink():
            state = (current, tuple(pending))
            if state in seen:
                raise ValueError("Cyclic alias in declared-list --output-dir")
            seen.add(state)
            target = current.readlink()
            expanded = target if target.is_absolute() else current.parent / target
            current = Path(expanded.anchor)
            pending = [*expanded.parts[1:], *pending]
    return False


def _validate_declared_density_output(output_dir: Path) -> None:
    """Protect six-trait tables before declared-list input reads or writes (#1093).

    Reserve arm directory names in both requested and resolved paths, including
    nested destinations and aliases. An existing canonical table name is reserved
    regardless of its contents: partial files and dangling links must not be
    replaced. Declared-list reruns therefore need a fresh output directory.
    """
    resolved = output_dir.resolve()
    if (
        resolved.is_relative_to(DEFAULT_OUTPUT_DIR.resolve())
        or any(component.casefold() == "density_adjustment" for component in resolved.parts)
        or _has_density_directory_component(output_dir)
    ):
        raise ValueError(
            "A trait-list run needs its own --output-dir outside the six-trait "
            "density_adjustment/ directories, including their nested paths and aliases"
        )
    table = output_dir / "table1_density_adjusted_effects.tsv"
    if table.exists() or table.is_symlink():
        raise ValueError(
            "A trait-list run needs a fresh --output-dir: the reserved "
            "table1_density_adjusted_effects.tsv already exists. "
            "Use a new directory for each declared-list run."
        )


def run_density_audit(
    profiles_path: Path = DEFAULT_PROFILES_PATH,
    columns_path: Path = DEFAULT_COLUMNS_PATH,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    code_commit: str = "",
    trait_list: DeclaredTraitList | None = None,
) -> dict[str, Any]:
    """Execute the complete trait density adjustment and near-twin audit.

    Writes:
      output_dir / 'density_adjustment_summary.tsv'
      output_dir / 'table1_density_adjusted_effects.tsv'
      output_dir / 'density_adjustment_report.txt'

    With a declared ``trait_list`` the effects table reports those traits instead of
    the six Table 1 traits, and ``trait_list_provenance.json`` records the list (#1054).
    A list run never writes into default or arm density directories, or replaces
    an existing canonical density table. Its destination is checked before any
    input reads or directory creation (#1083, #1093).
    """
    if trait_list is not None:
        _validate_declared_density_output(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    profiles_df, non_media_cols = load_non_media_profile_matrix(profiles_path, columns_path)
    table1_traits = table1_traits_for_columns(columns_path)
    reported_traits = (
        table1_traits if trait_list is None else declared_density_traits(columns_path, trait_list)
    )
    test_pair_a_label = reported_traits.get(TEST_PAIR_A, table1_traits[TEST_PAIR_A])
    raw_matrix = np.nan_to_num(profiles_df[non_media_cols].to_numpy(dtype=float), nan=0.0)

    # 1. Baseline metrics (Pearson & Spearman)
    raw_pearson = compute_correlation_metrics(raw_matrix, non_media_cols, method="pearson")
    raw_spearman = compute_correlation_metrics(raw_matrix, non_media_cols, method="spearman")

    # 2. Method A: PC1 Regression Residualization
    adj_matrix_pc1, _pc1_scores, _pc1_var_pct = regress_out_pc1(raw_matrix)
    adj_pc1_pearson = compute_correlation_metrics(adj_matrix_pc1, non_media_cols, method="pearson")

    # 3. Method B: Profile Composition / Annotation Normalization (#352, #379)
    comp_matrix = normalize_traits_by_profile_composition(profiles_df, non_media_cols)
    adj_comp_pearson = compute_correlation_metrics(comp_matrix, non_media_cols, method="pearson")

    # 4. Table 1 Trait Effects under PC1 adjustment
    t1_effects = evaluate_table1_paired_effects(
        profiles_df, adj_matrix_pc1, non_media_cols, reported_traits
    )

    # 5. Near-twin conditional test
    twin_summary = evaluate_conditional_twin_test(
        profiles_df, TWIN_A, TWIN_B, columns_for_depth=non_media_cols
    )

    # 6. Save Table 1 Effects TSV
    t1_df = pd.DataFrame([asdict(e) for e in t1_effects])
    t1_out = output_dir / "table1_density_adjusted_effects.tsv"
    t1_df.to_csv(t1_out, sep="\t", index=False)

    # 7. Save Summary Metrics TSV
    summary_rows = [
        {"metric_set": "raw_pearson", **asdict(raw_pearson)},
        {"metric_set": "raw_spearman", **asdict(raw_spearman)},
        {"metric_set": "adjusted_pc1_pearson", **asdict(adj_pc1_pearson)},
        {"metric_set": "adjusted_profile_composition_pearson", **asdict(adj_comp_pearson)},
    ]
    summary_df = pd.DataFrame(summary_rows)
    summary_out = output_dir / "density_adjustment_summary.tsv"
    summary_df.to_csv(summary_out, sep="\t", index=False)

    # 8. Save Comprehensive Text Report (fully dynamic, no hard-coded values)
    report_out = output_dir / "density_adjustment_report.txt"
    report_lines = [
        "=" * 80,
        "TRAIT CORRELATION, ANNOTATION DENSITY ADJUSTMENT, AND NEAR-TWIN AUDIT",
        "=" * 80,
    ]
    if code_commit:
        report_lines.append(f"Code commit: {code_commit}")
    report_lines.extend(
        [
            f"Profiles source: {profiles_path}",
            f"Columns source: {columns_path}",
        ]
    )
    if trait_list is not None:
        report_lines.append(
            f"Trait list: {trait_list.source_name} (sha256 {trait_list.sha256}, "
            f"{len(trait_list.trait_ids)} traits)"
        )
    report_lines.extend(
        [
            "",
            "1. BASELINE DIMENSIONALITY AND CORRELATION STRUCTURE",
            "-" * 80,
            f"Input profiles evaluated: {raw_pearson.n_samples:,}",
            f"Non-MediaDive trait columns: {raw_pearson.n_traits} "
            f"({raw_pearson.constant_traits_excluded} constant excluded)",
            f"Pearson pairwise correlation: median {raw_pearson.median_pairwise_r:+.4f}, "
            f"mean {raw_pearson.mean_pairwise_r:+.4f}",
            f"Pairs exceeding +0.5: {raw_pearson.pct_pairs_gt_05:.1f}% | "
            f"Pairs exceeding +0.9: {raw_pearson.pct_pairs_gt_09:.1f}%",
            f"PC1 variance explained: {raw_pearson.pc1_variance_pct:.2f}% | "
            f"PC1-PC5 variance explained: {raw_pearson.pc1_5_variance_pct:.2f}%",
            f"Effective number of independent tests: "
            f"Li & Ji (2005) = {raw_pearson.m_eff_li_ji:.1f}, "
            f"Galwey (2009) = {raw_pearson.m_eff_galwey:.1f}, "
            f"Cheverud/Nyholt = {raw_pearson.m_eff_nyholt:.1f} (unreliable, #325, #357)",
            "",
            "Spearman baseline verification:",
            f"  Median pairwise rho: {raw_spearman.median_pairwise_r:+.4f}, "
            f"Mean pairwise rho: {raw_spearman.mean_pairwise_r:+.4f}",
            f"  Pairs > 0.5: {raw_spearman.pct_pairs_gt_05:.1f}% | "
            f"Pairs > 0.9: {raw_spearman.pct_pairs_gt_09:.1f}%",
            f"  PC1 variance: {raw_spearman.pc1_variance_pct:.2f}%",
            "",
            "2. DENSITY AXIS ADJUSTMENT: SEPARATING BIOLOGY FROM CURATION EFFORT",
            "-" * 80,
            f"Test Case Pair: {test_pair_a_label} x moderately halophilic",
            f"  Baseline Pearson r: {raw_pearson.test_pair_r:+.4f}",
            f"  Post-PC1 regression r: {adj_pc1_pearson.test_pair_r:+.4f} "
            f"(attenuation: {adj_pc1_pearson.test_pair_r - raw_pearson.test_pair_r:+.4f})",
            f"  Post-profile-composition r: {adj_comp_pearson.test_pair_r:+.4f}",
            "",
            "Global Correlation Matrix Centering (Method A - PC1 Regression):",
            f"  Median pairwise r: {adj_pc1_pearson.median_pairwise_r:+.4f} (centered around 0)",
            f"  Mean pairwise r: {adj_pc1_pearson.mean_pairwise_r:+.4f}",
            f"  Pairs > 0.5: {adj_pc1_pearson.pct_pairs_gt_05:.1f}%",
            "",
            f"3. {reported_traits_heading(trait_list)} "
            f"(PAIRED EXPERIMENTS N={twin_summary.n_paired_experiments:,})",
            "-" * 80,
        ]
    )

    for eff in t1_effects:
        p_raw_str = f"p={eff.raw_wilcoxon_p:.2e}" if not math.isnan(eff.raw_wilcoxon_p) else "p=NA"
        p_adj_str = f"p={eff.adj_wilcoxon_p:.2e}" if not math.isnan(eff.adj_wilcoxon_p) else "p=NA"
        report_lines.append(
            f"  {eff.trait_label:<30s} | "
            f"Raw mean: {eff.raw_mean_diff:+.3f} ({p_raw_str}) | "
            f"Adj mean: {eff.adj_mean_diff:+.3f} ({p_adj_str}) | "
            f"{eff.status}"
        )

    report_lines.extend(
        [
            "",
            "4. NEAR-TWIN CONDITIONAL TEST AND COLLINEARITY (COCCUS VS. SPHERE)",
            "-" * 80,
            f"Candidate pair: {twin_summary.trait_a_label} ({twin_summary.trait_a}) vs. "
            f"{twin_summary.trait_b_label} ({twin_summary.trait_b})",
            f"Profile-level correlation: r = {twin_summary.profile_correlation:+.4f} | "
            f"Partial r (controlling for depth): r = {twin_summary.profile_partial_r_depth:+.4f}",
            f"Variance Inflation Factor (VIF): {twin_summary.vif_profile:.2f}",
            f"Profile OLS regression slope (coccus ~ sphere): "
            f"{twin_summary.profile_slope_a_on_b:.4f}",
            f"Paired contrasts (N={twin_summary.n_paired_experiments}): "
            f"r = {twin_summary.paired_correlation:+.4f} | "
            f"Paired VIF: {twin_summary.vif_paired:.2f}",
            f"Paired regression slope: {twin_summary.paired_slope_a_on_b:.4f} "
            f"(p={twin_summary.paired_slope_pvalue:.2e})",
            "",
            "5. METHODOLOGICAL FINDINGS SUMMARY",
            "-" * 80,
            f"1. A substantial positive correlation component among traits is driven by PC1 "
            f"({raw_pearson.pc1_variance_pct:.2f}% variance explained).",
            f"2. Method A (PC1 regression) centers median pairwise correlation from "
            f"{raw_pearson.median_pairwise_r:+.4f} to {adj_pc1_pearson.median_pairwise_r:+.4f}.",
            f"3. Test case pair ({test_pair_a_label} x moderately halophilic) "
            f"correlation attenuates from r={raw_pearson.test_pair_r:+.4f} "
            f"to r={adj_pc1_pearson.test_pair_r:+.4f}.",
            f"4. Near-twins ({twin_summary.trait_a_label} and {twin_summary.trait_b_label}) show "
            f"profile VIF = {twin_summary.vif_profile:.2f} "
            f"(paired VIF = {twin_summary.vif_paired:.2f}), "
            "indicating strong collinearity.",
            "=" * 80,
        ]
    )

    report_text = "\n".join(report_lines) + "\n"
    report_out.write_text(report_text, encoding="utf-8")

    provenance_out = None
    if trait_list is None:
        # A default rerun into a reused directory must not keep a list it did not use.
        (output_dir / TRAIT_LIST_PROVENANCE_FILENAME).unlink(missing_ok=True)
    else:
        provenance_out = write_trait_list_provenance(
            output_dir,
            trait_list,
            profiles_path=profiles_path,
            columns_path=columns_path,
            code_commit=code_commit,
        )

    return {
        "raw_pearson": raw_pearson,
        "raw_spearman": raw_spearman,
        "adj_pc1_pearson": adj_pc1_pearson,
        "adj_comp_pearson": adj_comp_pearson,
        "table1_effects": t1_effects,
        "twin_summary": twin_summary,
        "summary_tsv": summary_out,
        "table1_tsv": t1_out,
        "report_txt": report_out,
        "trait_list_provenance": provenance_out,
    }


def add_trait_list_argument(parser: argparse.ArgumentParser) -> None:
    """Declare the optional trait list on each density command line (#1054)."""
    parser.add_argument(
        "--trait-list",
        type=Path,
        default=None,
        help=(
            "File of trait identifiers (one per line) to report instead of the six "
            "Table 1 traits; its SHA-256 is recorded in trait_list_provenance.json"
        ),
    )


def main() -> None:
    """CLI entrypoint for trait density audit."""
    parser = argparse.ArgumentParser(
        description="Audit trait correlation, density adjustment, and conditional near-twins."
    )
    parser.add_argument(
        "--profiles",
        type=Path,
        default=DEFAULT_PROFILES_PATH,
        help="Path to study_trait_profiles_wide.tsv",
    )
    parser.add_argument(
        "--columns",
        type=Path,
        default=DEFAULT_COLUMNS_PATH,
        help="Path to taxa_trait_matrix_columns.tsv",
    )
    parser.add_argument(
        "--output-dir",
        type=output_dir_argument,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory to write output TSVs and reports",
    )
    parser.add_argument(
        "--code-commit",
        type=str,
        default="",
        help="Code commit hash for provenance tracking",
    )
    add_trait_list_argument(parser)

    args = parser.parse_args()
    if args.trait_list is not None:
        _validate_declared_density_output(args.output_dir)
    results = run_density_audit(
        profiles_path=args.profiles,
        columns_path=args.columns,
        output_dir=args.output_dir,
        code_commit=args.code_commit,
        trait_list=None if args.trait_list is None else load_trait_list(args.trait_list),
    )
    print(f"Density adjustment audit completed successfully. Outputs written to {args.output_dir}")
    print(f"Summary TSV: {results['summary_tsv']}")
    print(f"Table 1 Effects TSV: {results['table1_tsv']}")
    print(f"Report: {results['report_txt']}")


if __name__ == "__main__":
    main()
