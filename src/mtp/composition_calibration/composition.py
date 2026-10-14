"""Paired composition arm calibration against the study-flip complete null (#263).

Evaluates empirical Type I error control, family-wise error rates, and Table 1
canonical trait calibration for the study-aggregated paired compositional inference
arm (JedCompositionAdapter) under a label-swapping null benchmark.

For a random half of the 959 paired studies, the increased and decreased direction
labels are exchanged. In paired composition inference, exchanging direction labels
inverts the study-level mean contrast (Delta'_sj = -Delta_sj) while holding all
profile trait counts, normalizations, and column z-score scalings invariant.

Produces:
1. composition_calibration_summary.json: Full metadata, Type I errors, and family metrics.
2. composition_calibration_summary.tsv: Tabular comparison of observed vs null calls.
3. composition_calibration_replicates.tsv: Per-replicate call counts by threshold.
4. composition_null_pvalues.tsv.gz: All null p-values across all simulation replicates.
5. table1_composition_calibration.tsv: Table 1 traits observed vs null distribution.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import platform
import sys
import time
from collections.abc import Sequence
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from importlib import metadata
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

from mtp.composition_calibration.adapters.jed_composition import (
    DEFAULT_PSEUDOCOUNT_SCALE,
    JedCompositionAdapter,
    compute_composition_contrasts,
)
from mtp.composition_calibration.scenarios import replicate_rng
from mtp.density_adjustment import TABLE_1_TRAITS
from mtp.enrichment_stats import benjamini_hochberg
from mtp.paths import DERIVED_DATA_DIR

logger = logging.getLogger(__name__)

STREAM = "complete_null_study_flip"
UNIT = "study"
SEED = 182
DEFAULT_REPLICATES = 100
DEFAULT_PROFILES_PATH = DERIVED_DATA_DIR / "study_trait_profiles_wide.tsv"
DEFAULT_COLUMNS_PATH = DERIVED_DATA_DIR / "taxa_trait_matrix_columns.tsv"
DEFAULT_OUTPUT_DIR = DERIVED_DATA_DIR / "calibration" / "composition"

THRESHOLDS = (0.05, 0.01, 0.001)
CORRECTIONS = ("BH", "BY")

__all__ = [
    "CORRECTIONS",
    "DEFAULT_COLUMNS_PATH",
    "DEFAULT_OUTPUT_DIR",
    "DEFAULT_PROFILES_PATH",
    "DEFAULT_REPLICATES",
    "SEED",
    "STREAM",
    "THRESHOLDS",
    "UNIT",
    "ReplicateCallCount",
    "SummaryRow",
    "Table1CalibrationRow",
    "build_parser",
    "load_composition_trait_columns",
    "main",
    "run_composition_calibration",
    "save_composition_calibration_artifacts",
]


@dataclass(frozen=True)
class ReplicateCallCount:
    """One replicate's discovery count for one correction at one threshold."""

    seed: int
    stream: str
    unit: str
    replicate: int
    correction: str
    threshold: float
    n_tests: int
    significant: int


@dataclass(frozen=True)
class SummaryRow:
    """Consolidated observed vs null comparison for one correction and threshold."""

    correction: str
    threshold: float
    n_tests: int
    replicates: int
    mean_null_calls: float
    max_null_calls: int
    replicates_with_any_call: int
    family_error: float
    family_error_ci_low: float
    family_error_ci_high: float
    observed_calls: int
    replicates_reaching_observed: int


@dataclass(frozen=True)
class Table1CalibrationRow:
    """Table 1 trait calibration metrics under the complete null."""

    trait_id: str
    trait_name: str
    observed_effect: float
    observed_t: float
    observed_p: float
    observed_bh_q: float
    observed_by_q: float
    null_mean_t: float
    null_min_t: float
    null_max_t: float
    null_type_i_alpha_05: float
    null_reps_reaching_observed: int
    survives_null_calibration: bool


def _file_sha256(path: Path) -> str:
    """Compute SHA-256 hash of a file."""
    if not path.is_file():
        return ""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _software_versions() -> dict[str, str]:
    """Capture runtime software versions for reproducibility provenance."""
    versions = {"python": platform.python_version()}
    for pkg in ("mtp", "numpy", "pandas", "scipy"):
        try:
            versions[pkg] = metadata.version(pkg)
        except metadata.PackageNotFoundError:
            versions[pkg] = "unknown"
    return versions


def load_composition_trait_columns(
    columns_path: Path = DEFAULT_COLUMNS_PATH,
    profiles_df: pd.DataFrame | None = None,
) -> list[str]:
    """Load non-MediaDive trait columns retained for composition calibration."""
    if columns_path.is_file():
        cols_df = pd.read_csv(columns_path, sep="\t")
        if "object_category" in cols_df.columns and "column_id" in cols_df.columns:
            non_media = cols_df[cols_df["object_category"] != "METPO:1004005"]["column_id"]
            traits = [str(c) for c in non_media.tolist()]
            if profiles_df is not None:
                return [t for t in traits if t in profiles_df.columns]
            return traits

    if profiles_df is not None:
        meta_cols = {
            "_norm_dir",
            "body_site",
            "condition",
            "experiment_id",
            "inc_dec",
            "signature_size",
            "study_id",
            "taxa_count",
        }
        return [
            c for c in profiles_df.columns if c not in meta_cols and "mediadive" not in c.lower()
        ]

    raise ValueError(f"Could not resolve trait columns: {columns_path} not found")


def run_composition_calibration(
    profiles_df: pd.DataFrame,
    trait_columns: Sequence[str],
    replicates: int = DEFAULT_REPLICATES,
    seed: int = SEED,
    stream: str = STREAM,
    pseudocount_scale: float = DEFAULT_PSEUDOCOUNT_SCALE,
) -> tuple[
    dict[str, Any],
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
]:
    """Run full calibration battery for the paired composition arm."""
    traits = list(trait_columns)
    m = len(traits)
    if m == 0:
        raise ValueError("Trait list cannot be empty")

    contrast_data = compute_composition_contrasts(
        profiles_df, traits, pseudocount_scale=pseudocount_scale
    )
    study_contrasts = contrast_data.study_contrasts.copy()
    if study_contrasts.empty:
        raise ValueError("No paired studies available for composition inference")

    study_keys = sorted(study_contrasts.index.astype(str), key=repr)
    study_contrasts = study_contrasts.loc[study_keys]
    n_studies = len(study_keys)
    harmonic_sum = float(np.sum(1.0 / np.arange(1, m + 1)))

    logger.info(
        "Paired composition calibration: %d studies, %d pairs, %d traits, %d replicates",
        n_studies,
        contrast_data.n_pairs,
        m,
        replicates,
    )

    # 1. Observed Inference
    adapter = JedCompositionAdapter()
    obs_res = adapter.run_inference(
        profiles_df,
        pd.DataFrame(),
        config={"traits": traits, "pseudocount_scale": pseudocount_scale},
    )
    obs_df = obs_res.results_df.set_index("trait")

    obs_p = obs_df.loc[traits, "p_value"].to_numpy(dtype=float)
    obs_bh = obs_df.loc[traits, "q_value"].to_numpy(dtype=float)
    obs_by = np.clip(obs_bh * harmonic_sum, 0.0, 1.0)
    obs_t = obs_df.loc[traits, "t_statistic"].to_numpy(dtype=float)
    obs_eff = obs_df.loc[traits, "effect_size"].to_numpy(dtype=float)

    observed_calls: dict[tuple[str, float], int] = {}
    for thresh in THRESHOLDS:
        observed_calls[("BH", thresh)] = int(np.sum(obs_bh <= thresh))
        observed_calls[("BY", thresh)] = int(np.sum(obs_by <= thresh))

    # 2. Null Replicates via Shared Study-Flip Stream
    replicate_records: list[ReplicateCallCount] = []
    null_p_records: list[dict[str, Any]] = []
    null_t_mat = np.zeros((replicates, m), dtype=float)
    null_p_mat = np.zeros((replicates, m), dtype=float)

    study_mat_values = study_contrasts[traits].to_numpy(dtype=float)

    for rep in range(replicates):
        rng = replicate_rng(seed, stream, rep)
        draws = np.array([rng.random() for _ in study_keys])
        signs = np.where(draws < 0.5, -1.0, 1.0)[:, np.newaxis]
        sim_mat = study_mat_values * signs

        mean_r = np.mean(sim_mat, axis=0)
        std_r = np.std(sim_mat, axis=0, ddof=1)
        # Avoid division by zero if std is negligible
        valid_std = std_r > 1e-12
        t_r = np.zeros(m, dtype=float)
        p_r = np.ones(m, dtype=float)

        t_r[valid_std] = mean_r[valid_std] / (std_r[valid_std] / np.sqrt(n_studies))
        p_r[valid_std] = 2.0 * stats.t.sf(np.abs(t_r[valid_std]), df=n_studies - 1)

        null_t_mat[rep, :] = t_r
        null_p_mat[rep, :] = p_r

        bh_r = benjamini_hochberg(p_r)
        by_r = np.clip(bh_r * harmonic_sum, 0.0, 1.0)

        for thresh in THRESHOLDS:
            replicate_records.append(
                ReplicateCallCount(
                    seed=seed,
                    stream=stream,
                    unit=UNIT,
                    replicate=rep,
                    correction="BH",
                    threshold=thresh,
                    n_tests=m,
                    significant=int(np.sum(bh_r <= thresh)),
                )
            )
            replicate_records.append(
                ReplicateCallCount(
                    seed=seed,
                    stream=stream,
                    unit=UNIT,
                    replicate=rep,
                    correction="BY",
                    threshold=thresh,
                    n_tests=m,
                    significant=int(np.sum(by_r <= thresh)),
                )
            )

        for t_idx, trait_name in enumerate(traits):
            null_p_records.append(
                {
                    "null_p_value": p_r[t_idx],
                    "null_t_statistic": t_r[t_idx],
                    "replicate": rep,
                    "trait": trait_name,
                }
            )

    # 3. Summarize Replicates
    summary_rows: list[SummaryRow] = []
    for corr in CORRECTIONS:
        for thresh in THRESHOLDS:
            here = [
                r
                for r in replicate_records
                if r.correction == corr and abs(r.threshold - thresh) < 1e-6
            ]
            vals = [r.significant for r in here]
            any_call = sum(1 for v in vals if v > 0)
            ci = stats.binomtest(any_call, len(vals)).proportion_ci(
                confidence_level=0.95, method="exact"
            )
            seen = observed_calls.get((corr, thresh), 0)
            reached = sum(1 for v in vals if v >= seen)

            summary_rows.append(
                SummaryRow(
                    correction=corr,
                    threshold=thresh,
                    n_tests=m,
                    replicates=len(vals),
                    mean_null_calls=float(np.mean(vals)),
                    max_null_calls=int(np.max(vals)),
                    replicates_with_any_call=any_call,
                    family_error=float(any_call / len(vals)),
                    family_error_ci_low=float(ci.low),
                    family_error_ci_high=float(ci.high),
                    observed_calls=seen,
                    replicates_reaching_observed=reached,
                )
            )

    summary_df = pd.DataFrame([asdict(row) for row in summary_rows])
    replicates_df = pd.DataFrame([asdict(rec) for rec in replicate_records])
    null_pvals_df = pd.DataFrame(null_p_records)

    # 4. Overall Type I Error Metrics
    flat_null_p = null_p_mat.flatten()
    type_i_metrics = {
        f"alpha_{alpha}": {
            "empirical_type_i": float(np.mean(flat_null_p <= alpha)),
            "nominal_alpha": alpha,
            "total_null_tests": len(flat_null_p),
        }
        for alpha in (0.05, 0.01, 0.001)
    }

    # 5. Table 1 Trait Evaluation
    table1_rows: list[Table1CalibrationRow] = []
    trait_idx_map = {t: i for i, t in enumerate(traits)}

    for t_id, t_name in TABLE_1_TRAITS.items():
        if t_id not in trait_idx_map:
            continue
        idx = trait_idx_map[t_id]
        t_obs_val = float(obs_t[idx])
        p_obs_val = float(obs_p[idx])
        bh_obs_val = float(obs_bh[idx])
        by_obs_val = float(obs_by[idx])
        eff_obs_val = float(obs_eff[idx])

        null_ts = null_t_mat[:, idx]
        null_ps = null_p_mat[:, idx]

        reps_reached = int(np.sum(np.abs(null_ts) >= np.abs(t_obs_val)))
        null_type_i_05 = float(np.mean(null_ps <= 0.05))
        survives = bool(bh_obs_val <= 0.05 and reps_reached == 0)

        table1_rows.append(
            Table1CalibrationRow(
                trait_id=t_id,
                trait_name=t_name,
                observed_effect=eff_obs_val,
                observed_t=t_obs_val,
                observed_p=p_obs_val,
                observed_bh_q=bh_obs_val,
                observed_by_q=by_obs_val,
                null_mean_t=float(np.mean(null_ts)),
                null_min_t=float(np.min(null_ts)),
                null_max_t=float(np.max(null_ts)),
                null_type_i_alpha_05=null_type_i_05,
                null_reps_reaching_observed=reps_reached,
                survives_null_calibration=survives,
            )
        )

    table1_df = pd.DataFrame([asdict(r) for r in table1_rows])

    # Metadata & summary object
    summary_dict: dict[str, Any] = {
        "method_id": adapter.method_id,
        "owner": adapter.owner,
        "issue": "https://github.com/CultureBotAI/microbial-trait-profiling/issues/263",
        "seed": seed,
        "stream": stream,
        "unit": UNIT,
        "replicates_requested": replicates,
        "replicates_completed": replicates,
        "n_studies": n_studies,
        "n_pairs": contrast_data.n_pairs,
        "n_traits": m,
        "pseudocount_scale": pseudocount_scale,
        "type_i_error_metrics": type_i_metrics,
        "observed_calls": {f"{k[0]}_q_{k[1]}": v for k, v in observed_calls.items()},
        "summary": [asdict(row) for row in summary_rows],
        "table1_traits": [asdict(r) for r in table1_rows],
    }

    return summary_dict, summary_df, replicates_df, null_pvals_df, table1_df


def save_composition_calibration_artifacts(
    summary_dict: dict[str, Any],
    summary_df: pd.DataFrame,
    replicates_df: pd.DataFrame,
    null_pvals_df: pd.DataFrame,
    table1_df: pd.DataFrame,
    output_dir: Path,
) -> dict[str, Path]:
    """Persist calibration outputs to destination directory."""
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "null_pvalues_tsv": output_dir / "composition_null_pvalues.tsv.gz",
        "replicates_tsv": output_dir / "composition_calibration_replicates.tsv",
        "summary_json": output_dir / "composition_calibration_summary.json",
        "summary_tsv": output_dir / "composition_calibration_summary.tsv",
        "table1_tsv": output_dir / "table1_composition_calibration.tsv",
    }

    summary_df.to_csv(paths["summary_tsv"], sep="\t", index=False)
    replicates_df.to_csv(paths["replicates_tsv"], sep="\t", index=False)
    null_pvals_df.to_csv(paths["null_pvalues_tsv"], sep="\t", index=False, compression="gzip")
    table1_df.to_csv(paths["table1_tsv"], sep="\t", index=False)

    paths["summary_json"].write_text(
        json.dumps(summary_dict, indent=2) + "\n",
        encoding="utf-8",
    )

    return paths


def build_parser() -> argparse.ArgumentParser:
    """CLI argument parser for composition arm calibration."""
    parser = argparse.ArgumentParser(
        prog="python -m mtp.composition_calibration.composition",
        description="Run paired composition arm calibration against the study-flip null (#263).",
    )
    parser.add_argument(
        "--profiles",
        type=Path,
        default=DEFAULT_PROFILES_PATH,
        help="Wide study trait profiles TSV.",
    )
    parser.add_argument(
        "--columns",
        type=Path,
        default=DEFAULT_COLUMNS_PATH,
        help="Taxa-trait matrix columns TSV (for non-MediaDive trait filtering).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Destination directory for output artifacts.",
    )
    parser.add_argument(
        "--replicates",
        type=int,
        default=DEFAULT_REPLICATES,
        help="Number of study-flip null replicates (default: 100).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=SEED,
        help="Master seed for reproducible stream generation (default: 182).",
    )
    parser.add_argument(
        "--pseudocount-scale",
        type=float,
        default=DEFAULT_PSEUDOCOUNT_SCALE,
        help="Scale factor for log1p transformation (default: 1000.0).",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """CLI entrypoint for composition arm calibration."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = build_parser().parse_args(argv)
    t0 = time.monotonic()
    started = datetime.now(UTC)

    if not args.profiles.is_file():
        logger.error("Profiles file not found: %s", args.profiles)
        return 1

    profiles_df = pd.read_csv(args.profiles, sep="\t")
    trait_columns = load_composition_trait_columns(args.columns, profiles_df)

    summary_dict, summary_df, replicates_df, null_pvals_df, table1_df = run_composition_calibration(
        profiles_df=profiles_df,
        trait_columns=trait_columns,
        replicates=args.replicates,
        seed=args.seed,
        pseudocount_scale=args.pseudocount_scale,
    )

    provenance = {
        "elapsed_seconds": round(time.monotonic() - t0, 3),
        "finished_utc": datetime.now(UTC).isoformat(),
        "inputs": {
            "columns": {
                "path": str(args.columns),
                "sha256": _file_sha256(args.columns),
            },
            "profiles": {
                "path": str(args.profiles),
                "rows": len(profiles_df),
                "sha256": _file_sha256(args.profiles),
            },
        },
        "software": _software_versions(),
        "started_utc": started.isoformat(),
    }
    summary_dict["provenance"] = provenance

    saved = save_composition_calibration_artifacts(
        summary_dict, summary_df, replicates_df, null_pvals_df, table1_df, args.output_dir
    )

    print("\n--- Paired Composition Arm Calibration Summary (#263) ---")
    for row in summary_dict["summary"]:
        print(
            f"  {row['correction']:4} q<={row['threshold']:<6}: "
            f"observed {row['observed_calls']:>3}  |  "
            f"null mean {row['mean_null_calls']:.2f}  max {row['max_null_calls']:>2}  |  "
            f"any-call share {row['family_error']:.3f} "
            f"(95% CI {row['family_error_ci_low']:.3f}-{row['family_error_ci_high']:.3f}) | "
            f"reaching observed: {row['replicates_reaching_observed']}/{row['replicates']}"
        )

    print("\n--- Table 1 Canonical Traits Calibration ---")
    for t in summary_dict["table1_traits"]:
        print(
            f"  {t['trait_name']:<28}: "
            f"obs_t={t['observed_t']:+6.3f} (q={t['observed_bh_q']:.2e}) | "
            f"null_t range=[{t['null_min_t']:+6.3f}, {t['null_max_t']:+6.3f}] | "
            f"reaching |t_obs|={t['null_reps_reaching_observed']}/100 | "
            f"survives={t['survives_null_calibration']}"
        )

    print("\nSaved artifacts to:")
    for path in saved.values():
        print(f"  {path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
