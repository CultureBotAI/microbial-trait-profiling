"""Monte Carlo calibration harness for evaluating statistical pipelines.

Run with ``uv run python -m bsdbng.calibration.harness --help`` or ``just calibrate``.

Exit codes:
    0  every requested replicate completed and produced valid output
    1  invalid inputs/arguments, or the run failed / produced no valid replicate
    2  partial run: results were saved, but some replicates failed or were invalid

A zero exit code means the benchmark ran and its outputs were written. It does
not mean the method is calibrated; that is decided from the reported metrics
(see #263 and #182).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import platform
import sys
import time
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from importlib import metadata
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
    resolve_signature_sizes,
    shared_trait_columns,
    trait_population_counts,
)
from bsdbng.paths import DERIVED_DATA_DIR

logger = logging.getLogger(__name__)

EXIT_SUCCESS = 0
EXIT_FAILURE = 1
EXIT_PARTIAL = 2

SCENARIOS = ("all", "null", "mixed")
STATUS_SUCCESS = "SUCCESS"
STATUS_PARTIAL = "PARTIAL"
STATUS_FAILED = "FAILED"

REPLICATE_COMPLETED = "completed"
REPLICATE_FAILED = "failed"
REPLICATE_INVALID = "invalid"

SUMMARY_FILENAME = "calibration_summary.json"
REPLICATES_FILENAME = "calibration_replicates.tsv"
ECDF_FILENAME = "calibration_ecdf.tsv"
NULL_PVALUES_FILENAME = "calibration_null_pvalues.tsv.gz"
OUTPUT_FILENAMES = (SUMMARY_FILENAME, REPLICATES_FILENAME, ECDF_FILENAME, NULL_PVALUES_FILENAME)


@dataclass
class ReplicateRecord:
    """Outcome of one independent Monte Carlo replicate."""

    replicate_index: int
    seed_spawn_key: int
    status: str
    n_null_tests: int = 0
    n_mixed_tests: int = 0
    n_mixed_discoveries: int = 0
    n_true_signal_pairs: int = 0
    error_message: str | None = None


@dataclass
class CalibrationProvenance:
    """Provenance metadata for calibration execution."""

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
    """Consolidated calibration evaluation results."""

    method_id: str
    scenario: str
    requested_replicates: int
    replicate_records: list[ReplicateRecord]
    complete_null_type_i: TypeIErrorMetric | None
    complete_null_ecdf: pd.DataFrame
    mixed_signal_fdr: FDRMetric | None
    null_pvalues: pd.DataFrame
    spike_traits: list[str]
    adapter_metadata: dict[str, Any] = field(default_factory=dict)
    provenance: CalibrationProvenance | None = None

    def _count(self, status: str) -> int:
        return sum(1 for rec in self.replicate_records if rec.status == status)

    @property
    def completed_replicates(self) -> int:
        return self._count(REPLICATE_COMPLETED)

    @property
    def failed_replicates(self) -> int:
        return self._count(REPLICATE_FAILED)

    @property
    def invalid_replicates(self) -> int:
        return self._count(REPLICATE_INVALID)

    @property
    def status(self) -> str:
        if self.completed_replicates == 0:
            return STATUS_FAILED
        if self.completed_replicates < self.requested_replicates:
            return STATUS_PARTIAL
        return STATUS_SUCCESS

    def to_dict(self) -> dict[str, Any]:
        """Convert report to a JSON-serializable dictionary."""
        return {
            "status": self.status,
            "method_id": self.method_id,
            "scenario": self.scenario,
            "requested_replicates": self.requested_replicates,
            "completed_replicates": self.completed_replicates,
            "failed_replicates": self.failed_replicates,
            "invalid_replicates": self.invalid_replicates,
            "spike_traits": self.spike_traits,
            "n_null_pvalues": len(self.null_pvalues),
            "complete_null_type_i": (
                asdict(self.complete_null_type_i) if self.complete_null_type_i else None
            ),
            "mixed_signal_fdr": asdict(self.mixed_signal_fdr) if self.mixed_signal_fdr else None,
            "adapter_metadata": self.adapter_metadata,
            "replicates": [asdict(rec) for rec in self.replicate_records],
            "provenance": asdict(self.provenance) if self.provenance else None,
        }


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
    """Run a Monte Carlo calibration benchmark on an ORA adapter.

    Each replicate gets an independent generator spawned from ``seed`` via
    ``numpy.random.SeedSequence``, so replicate ``i`` is reproducible on its own.
    A replicate is ``invalid`` if it ran without error but produced no tests for
    a requested scenario, and ``failed`` if it raised.
    """
    if replicates <= 0:
        raise ValueError(f"Replicates count must be positive integer, got {replicates}")
    if not 0.0 < alpha < 1.0:
        raise ValueError(f"Alpha threshold must be in (0, 1), got {alpha}")
    fdr_q = alpha if fdr_q is None else fdr_q
    if not 0.0 < fdr_q < 1.0:
        raise ValueError(f"FDR target q must be in (0, 1), got {fdr_q}")
    if scenario not in SCENARIOS:
        raise ValueError(f"Unknown scenario '{scenario}'; expected one of {SCENARIOS}")
    if empirical_profiles_df.empty:
        raise ValueError("Empirical profiles DataFrame cannot be empty")
    if taxa_matrix_df.empty:
        raise ValueError("Taxa trait matrix DataFrame cannot be empty")

    run_null = scenario in ("all", "null")
    run_mixed = scenario in ("all", "mixed")
    if run_mixed and not spike_traits:
        raise ValueError("Mixed-signal scenario requires at least one spike trait")

    min_fold = config.get("min_fold", 1.5) if config and "min_fold" in config else 1.5
    adapter_config = {
        "alpha": fdr_q,
        "min_fold": min_fold,
        **(config or {}),
    }
    child_seeds = np.random.SeedSequence(seed).spawn(replicates)

    records: list[ReplicateRecord] = []
    null_frames: list[pd.DataFrame] = []
    replicate_null_rates: list[float] = []
    mixed_discoveries: list[set[tuple[str, str]]] = []
    all_true_signal_pairs: set[tuple[str, str]] = set()
    adapter_metadata: dict[str, Any] = {}

    for rep_idx, child in enumerate(child_seeds):
        rng = np.random.default_rng(child)
        record = ReplicateRecord(
            replicate_index=rep_idx,
            seed_spawn_key=int(child.spawn_key[-1]),
            status=REPLICATE_COMPLETED,
        )
        rep_null: pd.DataFrame | None = None
        rep_null_rate: float = 0.0
        rep_disc: set[tuple[str, str]] | None = None
        rep_truth: set[tuple[str, str]] = set()
        try:
            logger.info("Replicate %d/%d", rep_idx + 1, replicates)
            if run_null:
                null_sim_df = generate_complete_null_ora_profiles(
                    empirical_profiles_df, taxa_matrix_df, rng
                )
                null_res = adapter.run_inference(null_sim_df, taxa_matrix_df, adapter_config)
                adapter_metadata = null_res.metadata
                record.n_null_tests = len(null_res.results_df)
                rep_null_pvals = null_res.results_df["p_value"].to_numpy(dtype=float)
                rep_null_rate = (
                    float(np.mean(rep_null_pvals <= alpha)) if len(rep_null_pvals) > 0 else 0.0
                )
                rep_null = pd.DataFrame(
                    {
                        "replicate_index": rep_idx,
                        "null_p_value": rep_null_pvals,
                    }
                )
            if run_mixed:
                mixed_sim_df, rep_truth = generate_mixed_signal_ora_profiles(
                    empirical_profiles_df,
                    taxa_matrix_df,
                    spike_traits=spike_traits,
                    signal_strength=signal_strength,
                    rng=rng,
                )
                mixed_res = adapter.run_inference(mixed_sim_df, taxa_matrix_df, adapter_config)
                adapter_metadata = mixed_res.metadata
                record.n_mixed_tests = len(mixed_res.results_df)
                record.n_true_signal_pairs = len(rep_truth)
                sig = mixed_res.results_df[mixed_res.results_df["is_significant"].astype(bool)]
                rep_disc = set(
                    zip(sig["profile_id"].astype(str), sig["trait"].astype(str), strict=True)
                )
                record.n_mixed_discoveries = len(rep_disc)
        except Exception as err:
            logger.warning("Replicate %d failed: %s", rep_idx, err)
            record.status = REPLICATE_FAILED
            record.error_message = f"{type(err).__name__}: {err}"
            records.append(record)
            continue

        problems = []
        if run_null and record.n_null_tests == 0:
            problems.append("null scenario produced zero tests")
        if run_mixed and record.n_mixed_tests == 0:
            problems.append("mixed scenario produced zero tests")
        if run_mixed and record.n_true_signal_pairs == 0:
            problems.append("mixed scenario produced zero true signal pairs")
        if problems:
            record.status = REPLICATE_INVALID
            record.error_message = "; ".join(problems)
            logger.warning("Replicate %d invalid: %s", rep_idx, record.error_message)
        else:
            if rep_null is not None:
                null_frames.append(rep_null)
                replicate_null_rates.append(rep_null_rate)
            if rep_disc is not None:
                # Key pairs by replicate so each replicate's FDP is computed
                # against its own ground truth, not the union across replicates.
                tag = f"{rep_idx}|"
                mixed_discoveries.append({(tag + p, t) for p, t in rep_disc})
                all_true_signal_pairs.update((tag + p, t) for p, t in rep_truth)
        records.append(record)

    null_pvalues = (
        pd.concat(null_frames, ignore_index=True)
        if null_frames
        else pd.DataFrame(columns=["replicate_index", "null_p_value"])
    )
    pvals = null_pvalues["null_p_value"].to_numpy(dtype=float)

    type_i = (
        compute_type_i_error(pvals, alpha=alpha, replicate_rates=replicate_null_rates)
        if run_null and len(pvals)
        else None
    )
    ecdf_df = compute_discrete_ecdf(pvals)
    fdr_metric = (
        compute_empirical_fdr(mixed_discoveries, all_true_signal_pairs, target_q=fdr_q)
        if run_mixed and mixed_discoveries
        else None
    )

    return CalibrationBenchmarkReport(
        method_id=adapter.method_id,
        scenario=scenario,
        requested_replicates=replicates,
        replicate_records=records,
        complete_null_type_i=type_i,
        complete_null_ecdf=ecdf_df,
        mixed_signal_fdr=fdr_metric,
        null_pvalues=null_pvalues,
        spike_traits=list(spike_traits),
        adapter_metadata=adapter_metadata,
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
    # Summary last, so its presence implies the other files were written.
    paths["summary_json"].write_text(
        json.dumps(report.to_dict(), indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    return paths


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _software_versions() -> dict[str, str]:
    versions = {"python": platform.python_version()}
    for pkg in ("bsdbng", "numpy", "pandas", "scipy"):
        try:
            versions[pkg] = metadata.version(pkg)
        except metadata.PackageNotFoundError:
            versions[pkg] = "unknown"
    return versions


def _clear_previous_outputs(output_dir: Path) -> None:
    """Remove outputs of earlier runs so stale results cannot pass as current ones."""
    for name in OUTPUT_FILENAMES:
        (output_dir / name).unlink(missing_ok=True)


def select_spike_traits(
    profiles_df: pd.DataFrame, taxa_matrix_df: pd.DataFrame, n_traits: int
) -> list[str]:
    """Pick the ``n_traits`` most prevalent shared traits (ties broken by name)."""
    traits = shared_trait_columns(profiles_df, taxa_matrix_df)
    counts = trait_population_counts(taxa_matrix_df, traits).tolist()
    ranked = sorted(
        ((int(k), trait) for trait, k in zip(traits, counts, strict=True) if k > 0),
        key=lambda item: (-item[0], item[1]),
    )
    return [trait for _, trait in ranked[:n_traits]]


def _probability(value: str) -> float:
    prob = float(value)
    if not 0.0 < prob < 1.0:
        raise argparse.ArgumentTypeError(f"must be in (0, 1), got {value}")
    return prob


def _positive_int(value: str) -> int:
    num = int(value)
    if num <= 0:
        raise argparse.ArgumentTypeError(f"must be a positive integer, got {value}")
    return num


def build_parser() -> argparse.ArgumentParser:
    """Build the CLI parser for the calibration harness."""
    parser = argparse.ArgumentParser(
        prog="python -m bsdbng.calibration.harness",
        description="Run the Monte Carlo calibration benchmark for the ORA adapter.",
    )
    parser.add_argument(
        "--profiles",
        type=Path,
        default=DERIVED_DATA_DIR / "study_trait_profiles_wide.tsv",
        help="Empirical study trait profiles TSV (must include taxa_count).",
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
        help="Nominal level for the complete-null Type I error check.",
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


def parse_args(args: list[str] | None = None) -> argparse.Namespace:
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
    if not shared_trait_columns(profiles_df, taxa_matrix_df):
        raise ValueError("No trait columns shared between profiles and taxa matrix")
    sizes = resolve_signature_sizes(profiles_df)
    if not np.any(sizes > 0):
        raise ValueError("All profiles have signature size 0; nothing to test")
    return profiles_df, taxa_matrix_df


def main(args: list[str] | None = None) -> int:
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
            config={"test_type": parsed.test_type, "min_fold": parsed.min_fold},
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
            "Complete-null Type I at alpha=%.3f: %.4f (MC SE %.4f, %d p-values)",
            t1.alpha,
            t1.empirical_rate,
            t1.monte_carlo_se,
            t1.replicates,
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
