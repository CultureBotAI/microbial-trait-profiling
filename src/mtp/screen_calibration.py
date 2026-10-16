"""How many traits does the widened directional screen call when nothing is there? (#263)

The trait screen of #255 tests all 317 retained traits and corrects across the
whole family, twice: Benjamini-Yekutieli and Benjamini-Hochberg at the nominal
family size, and Benjamini-Yekutieli again at the effective family size estimated
from the study contrasts (#265). On the pinned corpus it calls 18 traits at BY
q <= 0.05 with m = 317 and 56 with m = 112. This module asks what those numbers
would be if no trait differed between the directions.

The null is the directional arm's complete null of #182, reused rather than
re-implemented: for a random half of the studies, the increased and decreased
labels of every profile in the study are exchanged, drawn from the same
``complete_null_study_flip`` stream with the same seed. Replicate ``i`` here is
therefore replicate ``i`` of the shared flip-assignment artefact of #263, and the
counts are measured on the identical null the other arms are offered.

Each replicate runs the production screen unchanged, including its re-estimate of
the effective family size, so the effective-size column is calibrated as the
procedure is actually run. Every call under this null is a false one, so the
per-replicate count is the false-discovery count, and the share of replicates
with any call is the family error of the whole screen. Neither is a verdict on
the observed 18 or 56; they are what the same procedure produces with no signal.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from statistics import fmean

from scipy.stats import binomtest

from mtp.calibration import SEED, SHARED_FLIP_REPLICATES, flip_directions, replicate_rng
from mtp.paths import DERIVED_DATA_DIR, file_sha256
from mtp.trait_guilds import (
    DEFAULT_EXPERIMENTS_PATH,
    DEFAULT_FULL_DUMP_PATH,
    DEFAULT_PROFILE_PATH,
    Profile,
    load_profiles,
    read_bugsigdb_vintage,
)
from mtp.trait_screen import (
    ScreenResult,
    candidate_slots,
    run_screen,
    screen_candidates,
    significant,
)

VERSION = "screen-calibration-v1"
DEFAULT_OUTPUT_DIR = DERIVED_DATA_DIR / "calibration" / "screen"
DEFAULT_SHARED_FLIPS_PATH = DERIVED_DATA_DIR / "calibration" / "calibration_flip_assignments.tsv"
#: The stream whose replicates these are; the shared artefact of #263 names it too.
STREAM = "complete_null_study_flip"
UNIT = "study"
#: Multi-threshold calibration grid across nominal and effective q columns (#341).
THRESHOLDS = (0.10, 0.05, 0.01, 0.001)
#: The q columns of the screen, each with the correction it applies and how the
#: family size behind it is read off the result.
CORRECTIONS = (
    ("by_q_nominal", "Benjamini-Yekutieli", "nominal"),
    ("bh_q_nominal", "Benjamini-Hochberg", "nominal"),
    ("by_q_effective", "Benjamini-Yekutieli", "effective"),
)


@dataclass(frozen=True)
class ReplicateCount:
    """One replicate's calls under one correction at one threshold."""

    replicate: int
    q_column: str
    correction: str
    family_basis: str
    family_size: int | None
    threshold: float
    significant: int


def _family_size(result: ScreenResult, basis: str) -> int | None:
    return result.family_size if basis == "nominal" else result.effective_family_size_used


def count_calls(result: ScreenResult, replicate: int) -> list[ReplicateCount]:
    """The number of traits each q column calls at each threshold in one screen."""
    counts: list[ReplicateCount] = []
    for column, correction, basis in CORRECTIONS:
        size = _family_size(result, basis)
        for threshold in THRESHOLDS:
            counts.append(
                ReplicateCount(
                    replicate=replicate,
                    q_column=column,
                    correction=correction,
                    family_basis=basis,
                    family_size=size,
                    threshold=threshold,
                    # An undefined effective size leaves the column blank, which
                    # is no call rather than a call of everything.
                    significant=len(significant(result.rows, column, threshold))
                    if size is not None
                    else 0,
                )
            )
    return counts


def replicate_flipped_profiles(
    profiles: Sequence[Profile],
    replicate: int,
    seed: int = SEED,
) -> list[Profile]:
    """The i-th replicate's flipped corpus under the production null."""
    return flip_directions(profiles, replicate_rng(seed, STREAM, replicate), unit=UNIT)


def run_null_screens(
    profiles: Sequence[Profile],
    columns: Sequence[str],
    replicates: int,
    seed: int = SEED,
) -> list[ReplicateCount]:
    """The production screen on ``replicates`` study-flipped copies of the corpus."""
    counts: list[ReplicateCount] = []
    for index in range(replicates):
        flipped = replicate_flipped_profiles(profiles, index, seed)
        counts.extend(count_calls(run_screen(flipped, columns), index))
    return counts


def summarize(
    counts: Sequence[ReplicateCount], observed: Mapping[tuple[str, float], int]
) -> list[dict[str, object]]:
    """Per correction and threshold: calls under the null beside the observed calls.

    ``family_error`` is one Bernoulli per replicate - any call at all - so its
    exact interval is over independent replicates. The mean count is the expected
    number of false calls; the share of replicates calling at least the observed
    number says how often the null alone reaches what the corpus shows.
    """
    rows: list[dict[str, object]] = []
    for column, correction, basis in CORRECTIONS:
        for threshold in THRESHOLDS:
            here = [c for c in counts if c.q_column == column and c.threshold == threshold]
            if not here:
                continue
            values = [c.significant for c in here]
            any_call = sum(1 for v in values if v > 0)
            interval = binomtest(any_call, len(values)).proportion_ci(
                confidence_level=0.95, method="exact"
            )
            seen = observed.get((column, threshold))
            reached = sum(1 for v in values if seen is not None and v >= seen)
            sizes = sorted({c.family_size for c in here if c.family_size is not None})
            rows.append(
                {
                    "q_column": column,
                    "correction": correction,
                    "family_basis": basis,
                    "family_size_min": sizes[0] if sizes else None,
                    "family_size_max": sizes[-1] if sizes else None,
                    "threshold": threshold,
                    "replicates": len(values),
                    "mean_null_calls": fmean(values),
                    "max_null_calls": max(values),
                    "replicates_with_any_call": any_call,
                    "family_error": any_call / len(values),
                    "family_error_ci_low": float(interval.low),
                    "family_error_ci_high": float(interval.high),
                    "observed_calls": seen,
                    "replicates_reaching_observed": (
                        reached if seen is not None and seen > 0 else None
                    ),
                }
            )
    return rows


def _format(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def _write(path: Path, rows: Sequence[Mapping[str, object]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        if rows:
            header = list(rows[0])
            writer.writerow(header)
            for row in rows:
                writer.writerow([_format(row[name]) for name in header])
    return path


def shared_flip_artefact_provenance(seed: int, replicates: int) -> dict[str, int | bool]:
    """Verify and record shared flip artefact seed and replicate count (#342)."""
    return {
        "shared_flip_artefact_seed": SEED,
        "shared_flip_artefact_replicates": SHARED_FLIP_REPLICATES,
        "matches_shared_flip_artefact": (seed == SEED and 0 < replicates <= SHARED_FLIP_REPLICATES),
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile-input", type=Path, default=DEFAULT_PROFILE_PATH)
    parser.add_argument("--experiments-input", type=Path, default=DEFAULT_EXPERIMENTS_PATH)
    parser.add_argument("--full-dump-input", type=Path, default=DEFAULT_FULL_DUMP_PATH)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--replicates",
        type=int,
        default=SHARED_FLIP_REPLICATES,
        help="Study-flip replicates; the default matches the shared artefact of #263",
    )
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument(
        "--shared-flips",
        type=Path,
        default=DEFAULT_SHARED_FLIPS_PATH,
        help="Shared flip-assignment TSV from #263 for provenance recording (#342)",
    )
    parser.add_argument("--code-commit", help="Source revision, supplied by the Justfile recipe")
    args = parser.parse_args(argv)

    columns = screen_candidates(args.profile_input)
    profiles, stats = load_profiles(
        args.profile_input,
        args.experiments_input,
        candidate_slots(columns),
        reject_missing_counts=True,
    )
    observed_result = run_screen(profiles, columns)
    observed = {(c.q_column, c.threshold): c.significant for c in count_calls(observed_result, -1)}
    counts = run_null_screens(profiles, columns, args.replicates, args.seed)

    output = args.output_dir
    per_replicate = [
        {
            "seed": args.seed,
            "stream": STREAM,
            "unit": UNIT,
            "replicate": c.replicate,
            "q_column": c.q_column,
            "correction": c.correction,
            "family_basis": c.family_basis,
            "family_size": c.family_size,
            "threshold": c.threshold,
            "significant": c.significant,
        }
        for c in counts
    ]
    summary = summarize(counts, observed)
    written = [
        _write(output / "screen_null_replicates.tsv", per_replicate),
        _write(output / "screen_null_summary.tsv", summary),
    ]
    replicates_completed = len({c.replicate for c in counts})
    provenance = {
        "version": VERSION,
        "issue": "https://github.com/CultureBotAI/microbial-trait-profiling/issues/263",
        "seed": args.seed,
        "stream": STREAM,
        "unit": UNIT,
        "replicates_completed": replicates_completed,
        **shared_flip_artefact_provenance(args.seed, replicates_completed),
        "rng": "numpy.random.default_rng([seed, sha256_64(stream), replicate])",
        "candidates": len(columns),
        "studies": observed_result.studies,
        "corpus": asdict(stats),
        "bugsigdb_vintage": read_bugsigdb_vintage(args.full_dump_input),
        "inputs": {
            "profile_sha256": file_sha256(args.profile_input),
            "experiments_sha256": file_sha256(args.experiments_input),
            "full_dump_sha256": file_sha256(args.full_dump_input),
            "shared_flip_assignments_sha256": (
                file_sha256(args.shared_flips) if args.shared_flips.is_file() else None
            ),
            "analysis_source_sha256": file_sha256(Path(__file__)),
        },
        "code_commit": args.code_commit,
    }
    provenance_path = output / "screen_null_provenance.json"
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    written.append(provenance_path)
    for row in summary:
        mean = float(str(row["mean_null_calls"]))
        share = float(str(row["family_error"]))
        print(
            f"  {row['q_column']:16} q<={row['threshold']:<6} "
            f"observed {row['observed_calls']!s:>4}  null mean {mean:.2f}  "
            f"max {row['max_null_calls']}  any-call share {share:.3f}"
        )
    for path in written:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
