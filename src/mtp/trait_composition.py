"""Paired trait log-ratios with one directional observation per study.

The reference is part of every hypothesis. Ratios cancel the profile total;
they do not identify absolute changes or remove annotation/selection biases.
Zero counts are excluded as complete experiment pairs, never imputed.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from collections.abc import Sequence
from dataclasses import asdict, dataclass, replace
from itertools import combinations
from pathlib import Path
from statistics import fmean, median

from scipy.stats import binomtest, false_discovery_control

from mtp.corpus_identity import CorpusIdentity, corpus_identity, summary_lines
from mtp.paths import DERIVED_DATA_DIR
from mtp.paths import file_sha256 as file_sha256
from mtp.trait_guilds import (
    DEFAULT_EXPERIMENTS_PATH,
    DEFAULT_FULL_DUMP_PATH,
    DEFAULT_PROFILE_PATH,
    DIRECTIONS,
    TABLE1_TRAITS,
    CorpusStats,
    Profile,
    load_profiles,
    read_bugsigdb_vintage,
    write_stats_report,
)

ANALYSIS_VERSION = "paired-log-ratio-sign-v1"
LOG_RATIO_TIE_TOLERANCE = 1e-12
CONFIDENCE_LEVEL = 0.95


@dataclass(frozen=True)
class StudyLogRatio:
    """Coverage and mean eligible paired log-ratio change for one study."""

    numerator: str
    reference: str
    study_id: str
    paired_experiments: int
    eligible_experiments: int
    excluded_zero_experiments: int
    mean_log2_ratio_change: float | None


@dataclass(frozen=True)
class LogRatioSummary:
    """A relative directional hypothesis, including its conditional population."""

    numerator: str
    reference: str
    studies_total: int
    studies_eligible: int
    studies_excluded: int
    studies_partially_observed: int
    studies_up: int
    studies_down: int
    studies_tied: int
    studies_scored: int
    paired_experiments: int
    eligible_experiments: int
    excluded_zero_experiments: int
    mean_log2_ratio_change: float | None
    median_log2_ratio_change: float | None
    share_ratio_up: float | None
    share_ci_low: float | None
    share_ci_high: float | None
    sign_p: float | None
    by_q: float | None
    status: str


@dataclass(frozen=True)
class LogRatioResult:
    summaries: list[LogRatioSummary]
    studies: list[StudyLogRatio]
    family_size: int


def _matched_profiles(
    profiles: Sequence[Profile], traits: Sequence[str]
) -> list[tuple[Profile, Profile]]:
    """Reject ambiguous/missing keys before forming experiment-level contrasts."""
    by_experiment: defaultdict[tuple[str, str], dict[str, Profile]] = defaultdict(dict)
    for profile in profiles:
        key = (profile.study_id, profile.experiment_id)
        if profile.inc_dec not in DIRECTIONS:
            raise ValueError(f"Unexpected direction for {key}: {profile.inc_dec}")
        if profile.inc_dec in by_experiment[key]:
            raise ValueError(f"Duplicate profile for {key}, {profile.inc_dec}")
        for trait in traits:
            if trait not in profile.guild_counts:
                raise ValueError(f"Missing trait {trait} for {key}; missing is not zero")
            if profile.guild_counts[trait] < 0:
                raise ValueError(f"Negative count for {trait} in {key}")
        by_experiment[key][profile.inc_dec] = profile
    pairs: list[tuple[Profile, Profile]] = []
    for key, directions in sorted(by_experiment.items()):
        if set(directions) != set(DIRECTIONS):
            raise ValueError(f"Unpaired experiment {key}; filter complete pairs first")
        inc, dec = directions["inc"], directions["dec"]
        if (inc.condition, inc.body_site) != (dec.condition, dec.body_site):
            raise ValueError(f"Mismatched experiment metadata for {key}")
        pairs.append((inc, dec))
    return pairs


def _summarize(studies: list[StudyLogRatio]) -> LogRatioSummary:
    values = [s.mean_log2_ratio_change for s in studies if s.mean_log2_ratio_change is not None]
    up = sum(v > LOG_RATIO_TIE_TOLERANCE for v in values)
    down = sum(v < -LOG_RATIO_TIE_TOLERANCE for v in values)
    scored = up + down
    share: float | None = None
    low: float | None = None
    high: float | None = None
    p_value: float | None = None
    if scored:
        test = binomtest(up, scored, p=0.5, alternative="two-sided")
        interval = test.proportion_ci(confidence_level=CONFIDENCE_LEVEL, method="exact")
        share, low, high, p_value = (
            up / scored,
            float(interval.low),
            float(interval.high),
            float(test.pvalue),
        )
    status = "tested" if scored else "all_tied" if values else "no_positive_pairs"
    first = studies[0]
    return LogRatioSummary(
        numerator=first.numerator,
        reference=first.reference,
        studies_total=len(studies),
        studies_eligible=len(values),
        studies_excluded=len(studies) - len(values),
        studies_partially_observed=sum(
            0 < s.eligible_experiments < s.paired_experiments for s in studies
        ),
        studies_up=up,
        studies_down=down,
        studies_tied=len(values) - scored,
        studies_scored=scored,
        paired_experiments=sum(s.paired_experiments for s in studies),
        eligible_experiments=sum(s.eligible_experiments for s in studies),
        excluded_zero_experiments=sum(s.excluded_zero_experiments for s in studies),
        mean_log2_ratio_change=fmean(values) if values else None,
        median_log2_ratio_change=median(values) if values else None,
        share_ratio_up=share,
        share_ci_low=low,
        share_ci_high=high,
        sign_p=p_value,
        by_q=None,
        status=status,
    )


def paired_log2_ratio_change(
    numerator_inc: int,
    reference_inc: int,
    numerator_dec: int,
    reference_dec: int,
) -> float | None:
    """Return the increased-minus-decreased log2 ratio for four positive counts.

    A nonpositive count excludes the whole pair, without a pseudocount. Taking
    logs before subtraction avoids overflow or underflow in intermediate ratios.
    """
    counts = (numerator_inc, reference_inc, numerator_dec, reference_dec)
    if min(counts) <= 0:
        return None
    a_inc, b_inc, a_dec, b_dec = (math.log2(count) for count in counts)
    return math.fsum((a_inc, -b_inc, -a_dec, b_dec))


def compute_log_ratio_shifts(
    profiles: Sequence[Profile], traits: Sequence[str] = tuple(TABLE1_TRAITS)
) -> LogRatioResult:
    """Test all unordered trait pairs, using positive counts on both directions.

    For A/B, calculate log2(A_inc/B_inc) - log2(A_dec/B_dec) per
    experiment, then average within study. No profile total enters the ratio.
    BY adjusts the complete declared pair family, including untestable pairs
    internally represented by p=1. Undefined results stay missing in outputs.
    """
    if len(traits) < 2 or len(set(traits)) != len(traits):
        raise ValueError("Provide at least two distinct traits")
    pairs = _matched_profiles(profiles, traits)
    if not pairs:
        raise ValueError("No paired profiles available for log-ratio analysis")
    summaries: list[LogRatioSummary] = []
    all_studies: list[StudyLogRatio] = []
    for numerator, reference in combinations(traits, 2):
        contrasts: defaultdict[str, list[float]] = defaultdict(list)
        totals: defaultdict[str, int] = defaultdict(int)
        for inc, dec in pairs:
            study = inc.study_id
            totals[study] += 1
            contrast = paired_log2_ratio_change(
                inc.guild_counts[numerator],
                inc.guild_counts[reference],
                dec.guild_counts[numerator],
                dec.guild_counts[reference],
            )
            if contrast is not None:
                contrasts[study].append(contrast)
        studies = [
            StudyLogRatio(
                numerator=numerator,
                reference=reference,
                study_id=study,
                paired_experiments=total,
                eligible_experiments=len(contrasts[study]),
                excluded_zero_experiments=total - len(contrasts[study]),
                mean_log2_ratio_change=fmean(contrasts[study]) if contrasts[study] else None,
            )
            for study, total in sorted(totals.items())
        ]
        all_studies.extend(studies)
        summaries.append(_summarize(studies))
    adjusted = false_discovery_control(
        [s.sign_p if s.sign_p is not None else 1.0 for s in summaries], method="by"
    )
    summaries = [
        replace(summary, by_q=float(q) if summary.sign_p is not None else None)
        for summary, q in zip(summaries, adjusted, strict=True)
    ]
    return LogRatioResult(summaries, all_studies, len(summaries))


def _display(value: float | None) -> str:
    return "NA" if value is None else f"{value:.4g}"


def write_log_ratio_outputs(
    result: LogRatioResult,
    output_dir: Path,
    *,
    profile_path: Path,
    experiments_path: Path,
    full_dump_path: Path,
    body_site: str | None = None,
    corpus_stats: CorpusStats | None = None,
    code_commit: str | None = None,
    identity: CorpusIdentity | None = None,
) -> list[Path]:
    """Write statistics, study-level exclusions, source fingerprints and a report."""
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = output_dir / "trait_log_ratio_stats.tsv"
    study_path = output_dir / "trait_log_ratio_studies.tsv"
    provenance_path = output_dir / "trait_log_ratio_provenance.json"
    report_path = output_dir / "trait_log_ratio_report.md"
    scope = {
        "analysis_version": ANALYSIS_VERSION,
        "inference_level": "body_site" if body_site is not None else "pooled",
        "body_site": body_site,
        "correction_method": "Benjamini-Yekutieli",
        "family_size": result.family_size,
    }
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        rows = [{**asdict(summary), **scope} for summary in result.summaries]
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    with study_path.open("w", encoding="utf-8", newline="") as handle:
        study_rows = [asdict(study) for study in result.studies]
        writer = csv.DictWriter(
            handle, fieldnames=list(study_rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(study_rows)
    resolved_vintage = (
        identity.bugsigdb_vintage if identity is not None else read_bugsigdb_vintage(full_dump_path)
    )
    resolved_identity = (
        identity
        if identity is not None
        else corpus_identity(profile_path, bugsigdb_vintage=resolved_vintage)
    )
    provenance = {
        **scope,
        "profile_path": str(profile_path),
        "profile_sha256": resolved_identity.matrix.sha256,
        "experiments_path": str(experiments_path),
        "experiments_sha256": file_sha256(experiments_path),
        "analysis_source_sha256": file_sha256(Path(__file__)),
        "shared_loader_source_sha256": file_sha256(Path(__file__).with_name("trait_guilds.py")),
        "code_commit": code_commit,
        "bugsigdb_vintage": resolved_vintage,
        "corpus_identity": resolved_identity.as_record(),
        "input_filter_counts_before_site_selection": asdict(corpus_stats) if corpus_stats else None,
        "trait_ids": TABLE1_TRAITS,
        "host": "Homo sapiens",
        "pairing": "study_id, experiment_id; both inc and dec required",
        "zero_policy": "exclude whole experiment pair if either trait is zero on either side",
        "study_weighting": "one mean of eligible experiment log-ratio contrasts per study",
        "null": "P(study mean log2 ratio change > 0 | non-tied, eligible) = 0.5",
        "tie_tolerance_log2": LOG_RATIO_TIE_TOLERANCE,
        "confidence_level": CONFIDENCE_LEVEL,
        "interval": "exact binomial interval for share_ratio_up, not for mean effect",
        "untestable_family_policy": "p=1 internally for BY; undefined p/q left blank in TSV",
    }
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    lines = [
        "# Paired trait log-ratio analysis",
        "",
        f"Scope: {body_site if body_site is not None else 'all body sites pooled'}. "
        f"Version: {ANALYSIS_VERSION}; {result.family_size} pair hypotheses, "
        "Benjamini-Yekutieli adjustment.",
        "",
        "Each study contributes one mean paired log2(A/B) change. Positive means "
        "A/B is higher in increased signatures. This is relative annotation "
        "representation, not an absolute trait change or independent replication.",
        "",
        "Only experiments with both traits positive on both sides are eligible; "
        "zero counts are not imputed. Studies with partial coverage use their "
        "eligible experiments only. Coverage can differ between comparisons.",
        "",
        "The confidence interval describes the share of non-tied studies shifting "
        "up; it is not an interval for the mean log-ratio effect. NA means the "
        "directional test was undefined, not a negative or nonsignificant result.",
        "",
        "| Numerator / reference | Eligible / total studies | Eligible / total pairs | "
        "Up / down / tied studies | Mean log2 change | Share up (95% CI) | BY q | Status |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for s in result.summaries:
        lines.append(
            f"| {s.numerator} / {s.reference} | {s.studies_eligible}/{s.studies_total} | "
            f"{s.eligible_experiments}/{s.paired_experiments} | "
            f"{s.studies_up}/{s.studies_down}/{s.studies_tied} | "
            f"{_display(s.mean_log2_ratio_change)} | {_display(s.share_ratio_up)} "
            f"({_display(s.share_ci_low)}, {_display(s.share_ci_high)}) | "
            f"{_display(s.by_q)} | {s.status} |"
        )
    lines.extend(
        [
            "",
            "The six Table 1 traits were selected using overlapping data. These "
            "are exploratory reference-dependent comparisons; BY does not remove "
            "that selection bias, shared annotations, or dependence between studies.",
            "",
            "Study-level coverage and effects: [TSV](trait_log_ratio_studies.tsv). "
            "Complete statistics: [TSV](trait_log_ratio_stats.tsv). "
            "Source hashes and analysis rules: [JSON](trait_log_ratio_provenance.json).",
        ]
    )
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return [summary_path, study_path, provenance_path, report_path]


def main(argv: Sequence[str] | None = None) -> int:
    from mtp.workflow.consumers import (
        bind_count_arguments,
        core_corpus_identity,
        preflight_consumer_outputs,
        write_consumer_record,
    )

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--core-run", type=Path, help="Verified core workflow run manifest")
    parser.add_argument("--profile-input", type=Path)
    parser.add_argument("--experiments-input", type=Path)
    parser.add_argument("--full-dump-input", type=Path)
    parser.add_argument("--output-dir", type=Path, default=DERIVED_DATA_DIR / "table1_composition")
    parser.add_argument("--body-site", help="Exact source label; no display cuts or Other pooling")
    parser.add_argument("--code-commit", help="Source revision, supplied by the Justfile recipe")
    args = parser.parse_args(argv)
    bundle = bind_count_arguments(args)
    if bundle is not None:
        preflight_consumer_outputs(
            bundle,
            args.output_dir / "core_consumption.json",
            output_directory=args.output_dir,
        )
    args.profile_input = args.profile_input or DEFAULT_PROFILE_PATH
    args.experiments_input = args.experiments_input or DEFAULT_EXPERIMENTS_PATH
    args.full_dump_input = args.full_dump_input or DEFAULT_FULL_DUMP_PATH
    vintage = read_bugsigdb_vintage(args.full_dump_input)
    identity = (
        core_corpus_identity(bundle, bugsigdb_vintage=vintage)
        if bundle is not None
        else corpus_identity(args.profile_input, bugsigdb_vintage=vintage)
    )
    print("\n".join(summary_lines(identity)))
    profiles, stats = load_profiles(
        args.profile_input, args.experiments_input, TABLE1_TRAITS, reject_missing_counts=True
    )
    write_stats_report(stats, vintage)
    if args.body_site is not None:
        profiles = [p for p in profiles if p.body_site == args.body_site]
    result = compute_log_ratio_shifts(profiles)
    paths = write_log_ratio_outputs(
        result,
        args.output_dir,
        profile_path=args.profile_input,
        experiments_path=args.experiments_input,
        full_dump_path=args.full_dump_input,
        body_site=args.body_site,
        corpus_stats=stats,
        code_commit=args.code_commit,
        identity=identity,
    )
    if bundle is not None:
        paths.append(
            write_consumer_record(
                bundle,
                args.output_dir / "core_consumption.json",
                consumer="mtp.trait_composition",
                outputs=paths,
            )
        )
    for path in paths:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
