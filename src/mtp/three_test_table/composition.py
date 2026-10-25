"""Full-family observed composition and secondary study signed-rank test (#1044).

Reuse the calibration adapter's observed inference: all retained trait fractions,
log1p(1000*f), population-SD column z scores, within-experiment inc-minus-dec,
then equal-weight study means. Zero-total profiles remain zero fractions. Adapter
SD thresholds of 1e-12 remain unchanged; its untestable p=1 placeholders enter
full-family correction, while this export leaves their inferential fields blank.

The secondary test is two-sided Wilcoxon signed rank on those study means,
zero_method='wilcox', method='approx', correction=False. Exact zero study means
are excluded from its ranks. All-zero or empty vectors are untestable, held at
p=1 for correction over the full catalog. This test never decides a headline.
No ranks choose an estimand: both signed means and their absolute magnitudes have
explicit competition ranks (minimum rank for ties). Inputs must already carry
the verified arm's human corpus scope; this function does not infer host species.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd
from numpy.typing import NDArray
from scipy.stats import wilcoxon

from mtp.composition_calibration.adapters.jed_composition import (
    PROFILE_NONTRAIT_COLUMNS,
    JedCompositionAdapter,
    compute_composition_contrasts,
)
from mtp.multiple_testing import adjusted_probabilities, load_policy
from mtp.study_trait_profiles import direction_to_inc_dec
from mtp.trait_guilds import MEDIADIVE_MARKER

COMPOSITION_CONTRACT = "full-catalog-study-delta-z-v1"
WILCOXON_CONTRACT = "study-mean-delta-z-wilcoxon-v1"
PROFILE_KEYS = ["study_id", "experiment_id", "inc_dec"]
MAX_EXACT_FLOAT_COUNT = 2**53 - 1


@dataclass(frozen=True)
class CompositionResult:
    """One row per declared trait and the exact experiment/study contrast matrices."""

    table: pd.DataFrame
    pair_contrasts: pd.DataFrame
    study_contrasts: pd.DataFrame


def _validated_profiles(profiles: pd.DataFrame, traits: Sequence[str]) -> pd.DataFrame:
    """Reject lossy inputs before the historical adapter's permissive coercion."""
    if not traits or any(not isinstance(trait, str) or not trait.strip() for trait in traits):
        raise ValueError("The retained trait catalog must contain nonblank trait identifiers")
    if len(set(traits)) != len(traits):
        raise ValueError("Duplicate trait identifiers in retained catalog")
    if not profiles.columns.is_unique:
        raise ValueError("Duplicate profile/count columns")
    if any(not isinstance(column, str) for column in profiles.columns):
        raise ValueError("Profile/count column identifiers must be strings")
    if not set(PROFILE_KEYS).issubset(profiles.columns):
        raise ValueError("Missing profile identity columns")
    count_columns = [str(column) for column in profiles if column not in PROFILE_NONTRAIT_COLUMNS]
    retained = [column for column in count_columns if MEDIADIVE_MARKER not in column.lower()]
    if set(retained) != set(traits):
        raise ValueError("Declared trait catalog differs from the complete retained count columns")
    data = profiles[[*PROFILE_KEYS, *count_columns]].copy()
    for column in PROFILE_KEYS:
        if any(not isinstance(value, str) or not value.strip() for value in data[column]):
            raise ValueError(f"Missing or blank profile identity: {column}")
    data["inc_dec"] = data.inc_dec.map(direction_to_inc_dec)
    if data.inc_dec.isna().any():
        raise ValueError("Unrecognized profile direction")
    if data.duplicated(PROFILE_KEYS).any():
        raise ValueError("Duplicate profile identity after direction normalization")
    try:
        counts = data[count_columns].apply(pd.to_numeric, errors="raise")
        values = counts.to_numpy(dtype=float)
    except (ValueError, TypeError) as error:
        raise ValueError("Trait counts must be finite nonnegative integers") from error
    if not np.isfinite(values).all() or np.any(values < 0) or np.any(values != np.floor(values)):
        raise ValueError("Trait counts must be finite nonnegative integers")
    # Bound sums before the adapter's integer aggregation and float conversion.
    # Beyond this limit a one-count distinction is not guaranteed in float64.
    retained_positions = [count_columns.index(column) for column in retained]
    if np.any(values > MAX_EXACT_FLOAT_COUNT) or np.any(
        values[:, retained_positions].sum(axis=1) > MAX_EXACT_FLOAT_COUNT
    ):
        raise ValueError("Trait counts and retained count totals exceed exact float64 range")
    data[count_columns] = counts
    return data


def _t_reason(values: NDArray[np.float64]) -> str:
    """Expose the adapter's study testability without duplicating its t statistic."""
    if len(values) < 2:
        return "insufficient_support"
    if float(np.std(values, ddof=1)) <= 1e-12:
        return "all_tied" if np.all(values == 0.0) else "insufficient_variance"
    return ""


def _add_ranks(table: pd.DataFrame) -> None:
    """Rank defined effect estimates, including defined zeros; never rank missing effects."""
    for effect in ("experiment_mean_delta_z", "study_mean_delta_z"):
        values = pd.to_numeric(table[effect])
        for name, ranked, ascending in (
            ("ascending", values, True),
            ("descending", values, False),
            ("absolute_descending", values.abs(), False),
        ):
            table[f"{effect}_rank_{name}"] = ranked.rank(method="min", ascending=ascending).astype(
                "Int64"
            )


def _native_number(table: pd.DataFrame, trait: str, column: str) -> float:
    """Require a finite scalar from the adapter for a declared testable trait."""
    value = table.at[trait, column]
    if not isinstance(value, (float, int, np.floating, np.integer)) or not np.isfinite(value):
        raise ValueError(f"Nonfinite or nonnumeric native statistic: {trait}, {column}")
    return float(value)


def composition_table(profiles: pd.DataFrame, traits: Sequence[str]) -> CompositionResult:
    """Build the declared all-trait table without calibration simulations or native edits.

    Returned inferential missing values are NaN for TSV export. Descriptive means
    remain defined for constant/all-zero study vectors. Unpaired input profiles
    are explicitly counted and excluded by the shared adapter's pairing logic.
    """
    traits = tuple(traits)
    data = _validated_profiles(profiles, traits)
    contrasts = compute_composition_contrasts(data, traits)
    pairs, studies = contrasts.pair_contrasts, contrasts.study_contrasts
    if (
        not np.isfinite(pairs.to_numpy(dtype=float)).all()
        or not np.isfinite(studies.to_numpy(dtype=float)).all()
    ):
        raise ValueError("Nonfinite transformed composition contrast")
    native: dict[str, pd.DataFrame] = {}
    adapter = JedCompositionAdapter()
    for correction in ("bh", "by"):
        native[correction] = adapter.run_inference(
            data,
            pd.DataFrame(),
            {"traits": traits, "correction": correction, "study_aggregation": True, "alpha": 0.05},
        ).results_df.set_index("trait")
    paired_keys = set(pairs.index)
    in_pair = [
        key in paired_keys
        for key in data[["study_id", "experiment_id"]].itertuples(index=False, name=None)
    ]
    paired_profiles = data.loc[in_pair]
    zero_profiles = int(paired_profiles[list(traits)].sum(axis=1).eq(0).sum())
    support = paired_profiles[list(traits)].gt(0).sum(axis=0)
    rows: list[dict[str, object]] = []
    wilcoxon_probabilities: list[float] = []
    for trait in traits:
        values = studies[trait].to_numpy(dtype=float)
        reason = _t_reason(values)
        nonzero = int(np.count_nonzero(values))
        wilcoxon_reason = (
            "insufficient_support" if not len(values) else "all_tied" if not nonzero else ""
        )
        if not reason and (trait not in native["bh"].index or trait not in native["by"].index):
            raise ValueError("Native adapter omitted a testable full-catalog trait")
        row: dict[str, object] = {
            "trait": trait,
            "composition_contract": COMPOSITION_CONTRACT,
            "composition_family_size": len(traits),
            "composition_input_profiles": len(data),
            "composition_profiles": len(paired_profiles),
            "composition_profile_support": int(support[trait]),
            "composition_profile_support_fraction": (
                int(support[trait]) / len(paired_profiles) if len(paired_profiles) else np.nan
            ),
            "composition_unpaired_profiles_excluded": len(data) - len(paired_profiles),
            "composition_zero_total_profiles": zero_profiles,
            "composition_pairs": len(pairs),
            "composition_studies": len(studies),
            "composition_effect_units": "standardized_log1p_1000_retained_fraction_difference",
            "composition_test_unit": "study",
            "composition_zero_total_policy": "retained_as_zero_fraction",
            "experiment_mean_delta_z": float(pairs[trait].mean()) if len(pairs) else np.nan,
            "study_mean_delta_z": float(values.mean()) if len(values) else np.nan,
            "study_delta_z_sd": float(np.std(values, ddof=1)) if len(values) >= 2 else np.nan,
            "study_t_status": "untestable" if reason else "tested",
            "study_t_reason": reason,
            "study_t_primary_correction": "BY",
            "study_wilcoxon_contract": WILCOXON_CONTRACT,
            "study_wilcoxon_status": "untestable" if wilcoxon_reason else "tested",
            "study_wilcoxon_reason": wilcoxon_reason,
            "study_wilcoxon_nonzero_studies": nonzero,
            "study_wilcoxon_zero_studies": len(values) - nonzero,
        }
        for output, column in (
            ("study_t_statistic", "t_statistic"),
            ("study_t_p", "p_value"),
            ("study_t_q_bh", "q_value"),
            ("study_mean_delta_z_ci95_low", "ci95_lower"),
            ("study_mean_delta_z_ci95_high", "ci95_upper"),
        ):
            row[output] = _native_number(native["bh"], trait, column) if not reason else np.nan
        row["study_t_q_by"] = (
            _native_number(native["by"], trait, "q_value") if not reason else np.nan
        )
        if wilcoxon_reason:
            row["study_wilcoxon_statistic"] = row["study_wilcoxon_p"] = np.nan
            wilcoxon_probabilities.append(1.0)
        else:
            signed_rank = wilcoxon(
                values,
                alternative="two-sided",
                zero_method="wilcox",
                method="approx",
                correction=False,
            )
            row["study_wilcoxon_statistic"] = float(signed_rank.statistic)
            row["study_wilcoxon_p"] = float(signed_rank.pvalue)
            wilcoxon_probabilities.append(float(signed_rank.pvalue))
        rows.append(row)
    by, bh = adjusted_probabilities(wilcoxon_probabilities, policy=load_policy())
    for row, q_by, q_bh in zip(rows, by, bh, strict=True):
        tested = row["study_wilcoxon_status"] == "tested"
        row["study_wilcoxon_q_by"] = float(q_by) if tested else np.nan
        row["study_wilcoxon_q_bh"] = float(q_bh) if tested else np.nan
    table = pd.DataFrame.from_records(rows)
    _add_ranks(table)
    return CompositionResult(table, pairs, studies)


__all__ = [
    "COMPOSITION_CONTRACT",
    "MAX_EXACT_FLOAT_COUNT",
    "PROFILE_KEYS",
    "WILCOXON_CONTRACT",
    "CompositionResult",
    "composition_table",
]
