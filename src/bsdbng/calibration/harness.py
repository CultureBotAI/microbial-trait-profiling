"""Monte Carlo statistical calibration harness for BugSigDB ORA and Depletion (#348, #351, #343).

Evaluates empirical Type I error control, discrete ECDF coverage, and False Discovery
Proportion (FDR) under:
1. Uniform hypergeometric urn draws without replacement (synthetic harness benchmark check)
2. Shared study-level direction flip null (#263 / #343)
3. Mixed signal profiles with known ground-truth spikes

Writes machine-readable calibration artifacts:
- ``calibration_summary.json``: parameters, provenance, metrics, and calibration verdict
- ``calibration_replicates.tsv``: per-replicate execution status, tests, and discoveries
- ``calibration_null_pvalues.tsv``: pooled null p-values across replicates
- ``calibration_ecdf.tsv``: empirical vs theoretical uniform discrete cumulative distribution
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
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from bsdbng.calibration.adapters.ilyas_ora import IlyasORAAdapter
from bsdbng.calibration.metrics import (
    FDRMetric,
    TypeIErrorMetric,
    compute_discrete_ecdf,
    compute_empirical_fdr,
    compute_type_i_error,
)
from bsdbng.calibration.protocol import MethodAdapter
from bsdbng.calibration.scenarios import (
    generate_complete_null_ora_profiles,
    generate_mixed_signal_ora_profiles,
    generate_study_flip_null_ora_profiles,
    resolve_signature_sizes,
    shared_trait_columns,
    trait_population_counts,
)
from bsdbng.calibration.study_flip import replicate_rng
from bsdbng.paths import DERIVED_DATA_DIR

logger = logging.getLogger(__name__)

STATUS_SUCCESS = "SUCCESS"
STATUS_PARTIAL = "PARTIAL"
STATUS_FAILED = "FAILED"

EXIT_SUCCESS = 0
EXIT_FAILURE = 1
EXIT_PARTIAL = 2

SCENARIOS = ("all", "null", "mixed")

SUMMARY_FILENAME = "calibration_summary.json"
REPLICATES_FILENAME = "calibration_replicates.tsv"
NULL_PVALUES_FILENAME = "calibration_null_pvalues.tsv"
ECDF_FILENAME = "calibration_ecdf.tsv"


@dataclass
class ReplicateRecord:
    """Per-replicate execution record across evaluated scenarios."""

    replicate_id: int
    scenario: str
    status: str
    n_tests: int = 0
    n_discoveries: int = 0
    type_i_rate: float | None = None
    fdp: float | None = None
    error_message: str | None = None


@dataclass
class CalibrationProvenance:
    """Full execution provenance including software environment and parameters."""

    started_utc: str
    finished_utc: str
    elapsed_seconds: float
    command_args: dict[str, Any]
    inputs: dict[str, dict[str, Any]]
    output_dir: str
    software: dict[str, str]
    requested_replicates: int
    completed_replicates: int
    failed_replicates: int
    invalid_replicates: int
    status: str
    error_message: str | None = None


@dataclass
class CalibrationBenchmarkReport:
    """Aggregated results, metrics, and provenance of a calibration benchmark run."""

    status: str
    requested_replicates: int
    completed_replicates: int
    failed_replicates: int
    invalid_replicates: int
    replicate_records: list[ReplicateRecord] = field(default_factory=list)
    null_pvalues: pd.DataFrame = field(default_factory=lambda: pd.DataFrame(columns=["p_value"]))
    complete_null_ecdf: pd.DataFrame = field(default_factory=pd.DataFrame)
    complete_null_type_i: TypeIErrorMetric | None = None
    study_flip_null_type_i: TypeIErrorMetric | None = None
    mixed_signal_fdr: FDRMetric | None = None
    provenance: CalibrationProvenance | None = None

    @property
    def harness_check_passed(self) -> bool:
        """Verdict: True if null benchmark was conducted and is conservative (#389)."""
        if self.complete_null_type_i is not None and not self.complete_null_type_i.is_conservative:
            return False
        if (
            self.study_flip_null_type_i is not None
            and not self.study_flip_null_type_i.is_conservative
        ):
            return False
        return (self.complete_null_type_i is not None) or (self.study_flip_null_type_i is not None)

    @property
    def is_calibrated(self) -> bool:
        """Backward-compatible alias for harness_check_passed (#389)."""
        return self.harness_check_passed

    def to_dict(self) -> dict[str, Any]:
        """Convert report to JSON-serializable dictionary."""
        return {
            "status": self.status,
            "harness_check_passed": self.harness_check_passed,
            "is_calibrated": self.is_calibrated,
            "requested_replicates": self.requested_replicates,
            "completed_replicates": self.completed_replicates,
            "failed_replicates": self.failed_replicates,
            "invalid_replicates": self.invalid_replicates,
            "n_null_pvalues": len(self.null_pvalues),
            "complete_null_type_i": (
                asdict(self.complete_null_type_i) if self.complete_null_type_i else None
            ),
            "study_flip_null_type_i": (
                asdict(self.study_flip_null_type_i) if self.study_flip_null_type_i else None
            ),
            "mixed_signal_fdr": (asdict(self.mixed_signal_fdr) if self.mixed_signal_fdr else None),
            "provenance": asdict(self.provenance) if self.provenance else None,
        }


def _clear_previous_outputs(output_dir: Path) -> None:
    for filename in (SUMMARY_FILENAME, REPLICATES_FILENAME, NULL_PVALUES_FILENAME, ECDF_FILENAME):
        p = output_dir / filename
        if p.exists():
            p.unlink()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _software_versions() -> dict[str, str]:
    return {
        "python": platform.python_version(),
        "numpy": str(np.__version__),
        "pandas": str(pd.__version__),
        "platform": platform.platform(),
    }


def select_spike_traits(
    profiles_df: pd.DataFrame, taxa_matrix_df: pd.DataFrame, n_spikes: int = 3
) -> list[str]:
    """Select the most prevalent shared traits in the taxa matrix as spike targets."""
    traits = shared_trait_columns(profiles_df, taxa_matrix_df)
    if not traits:
        raise ValueError("No shared trait columns between profiles and taxa matrix")
    k_counts = trait_population_counts(taxa_matrix_df, traits)
    ranked_indices = np.argsort(-k_counts)
    selected = [traits[i] for i in ranked_indices[: min(n_spikes, len(traits))]]
    return selected


def run_ora_calibration_benchmark(
    adapter: MethodAdapter,
    empirical_profiles_df: pd.DataFrame,
    taxa_matrix_df: pd.DataFrame,
    spike_traits: list[str],
    replicates: int = 50,
    seed: int = 42,
    alpha: float = 0.05,
    scenario: str = "all",
    config: dict[str, Any] | None = None,
    fdr_q: float | None = None,
    signal_strength: float = 0.6,
) -> CalibrationBenchmarkReport:
    """Run a Monte Carlo calibration benchmark on an ORA adapter (#348, #351, #343)."""
    if replicates <= 0:
        raise ValueError(f"Replicates count must be positive integer, got {replicates}")
    if not 0.0 < alpha < 1.0:
        raise ValueError(f"Alpha threshold must be in (0, 1), got {alpha}")
    fdr_q = alpha if fdr_q is None else fdr_q
    if not 0.0 < fdr_q < 1.0:
        raise ValueError(f"FDR target q must be in (0, 1), got {fdr_q}")
    if scenario not in (*SCENARIOS, "study_flip"):
        raise ValueError(f"Unknown scenario '{scenario}'; expected one of {SCENARIOS}")
    if empirical_profiles_df.empty:
        raise ValueError("Empirical profiles DataFrame cannot be empty")
    if taxa_matrix_df.empty:
        raise ValueError("Taxa matrix DataFrame cannot be empty")

    run_uniform_null = scenario in ("all", "null")
    # dropped from active ORA calibration scenarios (#388)
    run_study_flip = scenario == "study_flip"
    run_mixed = scenario in ("all", "mixed")
    test_type = config.get("test_type", "enrichment") if config else "enrichment"

    replicate_records: list[ReplicateRecord] = []
    null_pvals_all: list[float] = []
    study_flip_pvals_all: list[float] = []
    null_rep_rates: list[float] = []
    study_flip_rep_rates: list[float] = []
    mixed_rep_discoveries: list[set[tuple[str, str]]] = []
    true_spikes_pooled: set[tuple[str, str]] = set()

    n_completed = 0
    n_failed = 0
    n_invalid = 0

    for rep_idx in range(replicates):
        rep_success = True

        # 1. Uniform draw null scenario (synthetic harness benchmark)
        if run_uniform_null:
            null_rng = replicate_rng(seed, "ora_uniform_null", rep_idx)
            try:
                null_df = generate_complete_null_ora_profiles(
                    empirical_profiles_df, taxa_matrix_df, null_rng
                )
                null_res = adapter.run_inference(null_df, taxa_matrix_df, config=config)
                pvals = null_res.results_df["p_value"].dropna().to_numpy(dtype=float)
                if len(pvals) == 0:
                    n_invalid += 1
                    replicate_records.append(
                        ReplicateRecord(
                            replicate_id=rep_idx,
                            scenario="null",
                            status="INVALID",
                            error_message="Zero tests evaluated",
                        )
                    )
                    rep_success = False
                else:
                    sig_cnt = int(np.sum(pvals <= alpha))
                    rate = sig_cnt / len(pvals)
                    null_rep_rates.append(rate)
                    null_pvals_all.extend(pvals.tolist())
                    replicate_records.append(
                        ReplicateRecord(
                            replicate_id=rep_idx,
                            scenario="null",
                            status="SUCCESS",
                            n_tests=len(pvals),
                            n_discoveries=sig_cnt,
                            type_i_rate=rate,
                        )
                    )
            except Exception as exc:
                n_failed += 1
                replicate_records.append(
                    ReplicateRecord(
                        replicate_id=rep_idx,
                        scenario="null",
                        status="FAILED",
                        error_message=str(exc),
                    )
                )
                rep_success = False

        # 2. Shared study-level flip null scenario (#263, #343, #388)
        if run_study_flip:
            flip_rng = replicate_rng(seed, "ora_study_flip", rep_idx)
            try:
                flip_df = generate_study_flip_null_ora_profiles(empirical_profiles_df, flip_rng)
                flip_res = adapter.run_inference(flip_df, taxa_matrix_df, config=config)
                pvals = flip_res.results_df["p_value"].dropna().to_numpy(dtype=float)
                if len(pvals) == 0:
                    n_invalid += 1
                    replicate_records.append(
                        ReplicateRecord(
                            replicate_id=rep_idx,
                            scenario="study_flip",
                            status="INVALID",
                            error_message="Zero tests evaluated",
                        )
                    )
                    rep_success = False
                else:
                    sig_cnt = int(np.sum(pvals <= alpha))
                    rate = sig_cnt / len(pvals)
                    study_flip_rep_rates.append(rate)
                    study_flip_pvals_all.extend(pvals.tolist())
                    replicate_records.append(
                        ReplicateRecord(
                            replicate_id=rep_idx,
                            scenario="study_flip",
                            status="SUCCESS",
                            n_tests=len(pvals),
                            n_discoveries=sig_cnt,
                            type_i_rate=rate,
                        )
                    )
            except Exception as exc:
                n_failed += 1
                replicate_records.append(
                    ReplicateRecord(
                        replicate_id=rep_idx,
                        scenario="study_flip",
                        status="FAILED",
                        error_message=str(exc),
                    )
                )
                rep_success = False

        # 3. Mixed signal scenario (empirical FDR)
        if run_mixed:
            mixed_rng = replicate_rng(seed, "ora_mixed_signal", rep_idx)
            try:
                mixed_df, true_spikes = generate_mixed_signal_ora_profiles(
                    empirical_profiles_df,
                    taxa_matrix_df,
                    spike_traits=spike_traits,
                    signal_strength=signal_strength,
                    rng=mixed_rng,
                    test_type=test_type,
                )
                true_spikes_pooled.update(true_spikes)
                mixed_res = adapter.run_inference(mixed_df, taxa_matrix_df, config=config)
                sig_hits = mixed_res.results_df[mixed_res.results_df["is_significant"]]
                discoveries = set(zip(sig_hits["profile_id"], sig_hits["trait"], strict=True))
                mixed_rep_discoveries.append(discoveries)
                replicate_records.append(
                    ReplicateRecord(
                        replicate_id=rep_idx,
                        scenario="mixed",
                        status="SUCCESS",
                        n_tests=len(mixed_res.results_df),
                        n_discoveries=len(discoveries),
                    )
                )
            except Exception as exc:
                n_failed += 1
                replicate_records.append(
                    ReplicateRecord(
                        replicate_id=rep_idx,
                        scenario="mixed",
                        status="FAILED",
                        error_message=str(exc),
                    )
                )
                rep_success = False

        if rep_success:
            n_completed += 1

    overall_status = (
        STATUS_SUCCESS
        if n_completed == replicates
        else (STATUS_PARTIAL if n_completed > 0 else STATUS_FAILED)
    )

    t1_null = (
        compute_type_i_error(null_pvals_all, alpha=alpha, replicate_rates=null_rep_rates)
        if run_uniform_null and null_rep_rates
        else None
    )
    t1_flip = (
        compute_type_i_error(
            study_flip_pvals_all, alpha=alpha, replicate_rates=study_flip_rep_rates
        )
        if run_study_flip and study_flip_rep_rates
        else None
    )
    fdr_metric = (
        compute_empirical_fdr(mixed_rep_discoveries, true_spikes_pooled, target_q=fdr_q)
        if run_mixed and mixed_rep_discoveries
        else None
    )

    # Use null p-values for discrete ECDF check
    pvals_for_ecdf = null_pvals_all if null_pvals_all else study_flip_pvals_all
    ecdf_df = compute_discrete_ecdf(pvals_for_ecdf)
    null_pvals_df = pd.DataFrame({"p_value": pvals_for_ecdf})

    return CalibrationBenchmarkReport(
        status=overall_status,
        requested_replicates=replicates,
        completed_replicates=n_completed,
        failed_replicates=n_failed,
        invalid_replicates=n_invalid,
        replicate_records=replicate_records,
        null_pvalues=null_pvals_df,
        complete_null_ecdf=ecdf_df,
        complete_null_type_i=t1_null,
        study_flip_null_type_i=t1_flip,
        mixed_signal_fdr=fdr_metric,
    )


def save_calibration_results(
    report: CalibrationBenchmarkReport,
    output_dir: Path,
) -> dict[str, Path]:
    """Save machine-readable summary, per-replicate status, ECDF, and null p-values."""
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "summary_json": output_dir / SUMMARY_FILENAME,
        "replicates_tsv": output_dir / REPLICATES_FILENAME,
        "ecdf_tsv": output_dir / ECDF_FILENAME,
        "null_pvalues_tsv": output_dir / NULL_PVALUES_FILENAME,
    }
    pd.DataFrame([asdict(rec) for rec in report.replicate_records]).to_csv(
        paths["replicates_tsv"], sep="\t", index=False
    )
    report.complete_null_ecdf.to_csv(paths["ecdf_tsv"], sep="\t", index=False)
    report.null_pvalues.to_csv(paths["null_pvalues_tsv"], sep="\t", index=False)
    paths["summary_json"].write_text(
        json.dumps(report.to_dict(), indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    return paths


def _positive_int(val: str) -> int:
    iv = int(val)
    if iv <= 0:
        raise argparse.ArgumentTypeError(f"Must be positive integer: {val}")
    return iv


def _probability(val: str) -> float:
    fv = float(val)
    if not 0.0 < fv < 1.0:
        raise argparse.ArgumentTypeError(f"Must be in (0, 1): {val}")
    return fv


def build_parser() -> argparse.ArgumentParser:
    """Build command-line parser for the calibration benchmark harness."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--profiles",
        type=Path,
        default=DERIVED_DATA_DIR / "study_trait_profiles_ora_wide.tsv",
        help="Input profiles TSV.",
    )
    parser.add_argument(
        "--taxa-matrix",
        type=Path,
        default=DERIVED_DATA_DIR / "taxa_trait_matrix.tsv",
        help="Taxa x trait matrix TSV (one row per taxon in the background universe).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DERIVED_DATA_DIR / "calibration",
        help="Destination directory for calibration outputs.",
    )
    parser.add_argument("--replicates", type=_positive_int, default=10)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--alpha",
        type=_probability,
        default=0.05,
        help="Nominal level for Type I error check.",
    )
    parser.add_argument(
        "--fdr-q",
        type=_probability,
        default=None,
        help="BH target q for discoveries in the mixed scenario (default: --alpha).",
    )
    parser.add_argument("--scenario", choices=SCENARIOS, default="all")
    parser.add_argument("--test-type", choices=["enrichment", "depletion"], default="enrichment")
    parser.add_argument(
        "--min-fold",
        type=float,
        default=1.5,
        help="Minimum fold enrichment for significance (default: 1.5).",
    )
    parser.add_argument(
        "--max-fold",
        type=float,
        default=0.67,
        help="Maximum fold depletion for significance (default: 0.67).",
    )
    parser.add_argument(
        "--signal-strength",
        type=float,
        default=0.6,
        help="Fraction of each signature replaced by spike-trait taxa in the mixed scenario.",
    )
    spikes = parser.add_mutually_exclusive_group()
    spikes.add_argument("--spike-traits", nargs="+", default=None, metavar="TRAIT")
    spikes.add_argument(
        "--n-spike-traits",
        type=_positive_int,
        default=3,
        help="If --spike-traits is not given, spike the N most prevalent traits.",
    )
    parser.add_argument(
        "--max-profiles",
        type=_positive_int,
        default=None,
        help="Use only the first N profile rows (for small reproducible runs).",
    )
    return parser


def parse_args(args: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse CLI arguments for the calibration harness."""
    return build_parser().parse_args(args)


def _load_inputs(
    parsed: argparse.Namespace,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    for label, path in (("Profiles", parsed.profiles), ("Taxa matrix", parsed.taxa_matrix)):
        if not path.is_file():
            raise ValueError(f"{label} input file not found: {path}")
    profiles_df = pd.read_csv(parsed.profiles, sep="\t")
    taxa_matrix_df = pd.read_csv(parsed.taxa_matrix, sep="\t", index_col=0)
    if parsed.max_profiles is not None:
        profiles_df = profiles_df.head(parsed.max_profiles)
    if profiles_df.empty:
        raise ValueError(f"Profiles TSV has no rows: {parsed.profiles}")
    if taxa_matrix_df.empty:
        raise ValueError(f"Taxa matrix TSV has no rows: {parsed.taxa_matrix}")
    if taxa_matrix_df.index.has_duplicates:
        raise ValueError(f"Taxa matrix has duplicate taxon IDs: {parsed.taxa_matrix}")
    shared = shared_trait_columns(profiles_df, taxa_matrix_df)
    if not shared:
        raise ValueError("No trait columns shared between profiles and taxa matrix")
    sizes = resolve_signature_sizes(profiles_df)
    if not np.any(sizes > 0):
        raise ValueError("All profiles have signature size 0; nothing to test")

    # Binarize taxa matrix trait columns so cells never exceed presence (> 0) (#348)
    for col in shared:
        taxa_matrix_df[col] = (
            pd.to_numeric(taxa_matrix_df[col], errors="coerce").fillna(0) > 0
        ).astype(np.int64)

    return profiles_df, taxa_matrix_df


def main(args: Sequence[str] | None = None) -> int:
    """CLI entry point for the calibration benchmark. Returns a process exit code."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    parsed = parse_args(args)
    started = datetime.now(UTC)
    t0 = time.monotonic()
    output_dir: Path = parsed.output_dir

    command_args = {k: (str(v) if isinstance(v, Path) else v) for k, v in vars(parsed).items()}
    inputs: dict[str, dict[str, Any]] = {}
    report: CalibrationBenchmarkReport | None = None
    error_message: str | None = None

    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        _clear_previous_outputs(output_dir)
        profiles_df, taxa_matrix_df = _load_inputs(parsed)
        inputs = {
            "profiles": {
                "path": str(parsed.profiles),
                "sha256": _sha256(parsed.profiles),
                "rows_used": len(profiles_df),
            },
            "taxa_matrix": {
                "path": str(parsed.taxa_matrix),
                "sha256": _sha256(parsed.taxa_matrix),
                "rows": len(taxa_matrix_df),
            },
        }
        spike_traits: list[str] = []
        if parsed.scenario in ("all", "mixed"):
            spike_traits = parsed.spike_traits or select_spike_traits(
                profiles_df, taxa_matrix_df, parsed.n_spike_traits
            )
        command_args["spike_traits"] = spike_traits
        logger.info(
            "Calibration: %d replicates, scenario=%s, test_type=%s, profiles=%d, taxa=%d",
            parsed.replicates,
            parsed.scenario,
            parsed.test_type,
            len(profiles_df),
            len(taxa_matrix_df),
        )
        report = run_ora_calibration_benchmark(
            adapter=IlyasORAAdapter(),
            empirical_profiles_df=profiles_df,
            taxa_matrix_df=taxa_matrix_df,
            spike_traits=spike_traits,
            replicates=parsed.replicates,
            seed=parsed.seed,
            alpha=parsed.alpha,
            fdr_q=parsed.fdr_q,
            scenario=parsed.scenario,
            signal_strength=parsed.signal_strength,
            config={
                "test_type": parsed.test_type,
                "min_fold": parsed.min_fold,
                "max_fold": parsed.max_fold,
            },
        )
    except Exception as err:
        error_message = f"{type(err).__name__}: {err}"
        logger.error("Calibration run failed: %s", error_message)

    finished = datetime.now(UTC)
    status = report.status if report is not None else STATUS_FAILED
    provenance = CalibrationProvenance(
        started_utc=started.isoformat(),
        finished_utc=finished.isoformat(),
        elapsed_seconds=round(time.monotonic() - t0, 3),
        command_args=command_args,
        inputs=inputs,
        output_dir=str(output_dir),
        software=_software_versions(),
        requested_replicates=parsed.replicates,
        completed_replicates=report.completed_replicates if report else 0,
        failed_replicates=report.failed_replicates if report else 0,
        invalid_replicates=report.invalid_replicates if report else 0,
        status=status,
        error_message=error_message,
    )

    try:
        if report is not None:
            report.provenance = provenance
            paths = save_calibration_results(report, output_dir)
        else:
            output_dir.mkdir(parents=True, exist_ok=True)
            paths = {"summary_json": output_dir / SUMMARY_FILENAME}
            paths["summary_json"].write_text(
                json.dumps({"status": STATUS_FAILED, "provenance": asdict(provenance)}, indent=2)
                + "\n",
                encoding="utf-8",
            )
    except OSError as err:
        logger.error("Could not write calibration outputs to %s: %s", output_dir, err)
        return EXIT_FAILURE

    if report is None:
        return EXIT_FAILURE

    logger.info(
        "Status %s: %d/%d replicates completed (%d failed, %d invalid). Summary: %s",
        status,
        report.completed_replicates,
        report.requested_replicates,
        report.failed_replicates,
        report.invalid_replicates,
        paths["summary_json"],
    )
    if report.complete_null_type_i is not None:
        t1 = report.complete_null_type_i
        logger.info(
            "Uniform null Type I at alpha=%.3f: %.4f (MC SE %.4f, %d p-values)",
            t1.alpha,
            t1.empirical_rate,
            t1.monte_carlo_se,
            t1.n_tests,
        )
    if report.study_flip_null_type_i is not None:
        t1_flip = report.study_flip_null_type_i
        logger.info(
            "Study-flip null Type I at alpha=%.3f: %.4f (MC SE %.4f, %d p-values)",
            t1_flip.alpha,
            t1_flip.empirical_rate,
            t1_flip.monte_carlo_se,
            t1_flip.n_tests,
        )
    if report.mixed_signal_fdr is not None:
        fdr = report.mixed_signal_fdr
        logger.info(
            "Mixed-signal empirical FDR at q=%.3f: %.4f (FDP sd %.4f, mean discoveries %.1f)",
            fdr.target_q,
            fdr.empirical_fdr,
            fdr.fdp_std,
            fdr.mean_discoveries,
        )
    return {STATUS_SUCCESS: EXIT_SUCCESS, STATUS_PARTIAL: EXIT_PARTIAL}.get(status, EXIT_FAILURE)


if __name__ == "__main__":
    sys.exit(main())
