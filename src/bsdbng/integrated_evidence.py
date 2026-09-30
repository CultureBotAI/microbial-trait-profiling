"""One table holding every method's result without merging them (#178).

Three frameworks report on the same matrix and mean different things by the word
"significant". A blank cell can mean untested, all tied, no annotation, masked from a
figure, or tested and not significant, and those are not the same claim. This module
builds a long table in which every result keeps its own row, its own family and its
own units, together with the states that a join would otherwise flatten.

What it deliberately does not do:

* It does not pool p-values, average effects across methods, or claim a false
  discovery rate for a set defined as "significant under any method". Separate
  corrections do not combine into one guarantee.
* It does not turn an absent source row into a negative result. A method whose
  output was not supplied is recorded as ``source_not_supplied`` for every feature it
  would have covered.
* It does not mark anything confirmed. The six Table 1 traits were selected on this
  same corpus, which every row carries in ``selection_note``.

Rows are produced for the analyses this repository runs, and read from files for the
other arms and for the assessments tracked in #182, #183, #184, #188 and #190. Each
of those is a status, not a result: ``pending`` until its own deliverable exists.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import asdict, dataclass, fields, replace
from pathlib import Path

from bsdbng.guild_maps import build_matrix, compute_correlation
from bsdbng.paths import DERIVED_DATA_DIR
from bsdbng.table1_trait_maps import compute_direction_shifts
from bsdbng.trait_composition import ANALYSIS_VERSION as LOG_RATIO_VERSION
from bsdbng.trait_composition import compute_log_ratio_shifts
from bsdbng.trait_guilds import (
    DEFAULT_EXPERIMENTS_PATH,
    DEFAULT_FULL_DUMP_PATH,
    DEFAULT_PROFILE_PATH,
    OXYGEN_GUILDS,
    TABLE1_TRAITS,
    Profile,
    drop_zero_total_pairs,
    load_profiles,
    read_bugsigdb_vintage,
    write_stats_report,
)

VERSION = "integrated-evidence-v1"
DEFAULT_OUTPUT_DIR = DERIVED_DATA_DIR / "integrated_evidence"
#: Where each assessment writes its own status. Absent means pending, never clear.
DEFAULT_CALIBRATION_STATUS = DERIVED_DATA_DIR / "calibration" / "calibration_status.json"
DEFAULT_SENSITIVITY_STATUS = (
    DERIVED_DATA_DIR / "table1_composition_sensitivity" / "sensitivity_status.json"
)
DEFAULT_ORIENTATION_STATUS = (
    DERIVED_DATA_DIR / "contrast_orientation" / "contrast_orientation_provenance.json"
)

#: Test states. They are recorded separately from the decision and from display.
TESTED = "tested"
UNTESTED = "untested"
INSUFFICIENT = "insufficient_support"
ALL_TIED = "all_tied"
ZERO_ANNOTATION = "zero_annotation_signal"
NOT_SAMPLED = "not_sampled"
SOURCE_MISSING = "source_not_supplied"

#: Decisions, which exist only where a threshold was declared and a test was run.
SIGNIFICANT = "significant"
NOT_SIGNIFICANT = "not_significant"
UNDECIDED = "undecided"

#: Display states, which never change a decision.
DRAWN = "drawn"
MASKED = "masked_below_min_studies"
NOT_DRAWN = "not_drawn"
NOT_APPLICABLE = "not_applicable"

SELECTION_NOTE = (
    "the six Table 1 traits were selected on this corpus by the composition arm; "
    "correcting within them is not independent confirmation"
)
#: No single cutoff is native to all three families collected here, so each row
#: is decided at both and each threshold names where it comes from (#294). The
#: composition arm publishes q with no cutoff at all, which is why neither of
#: these is described as the arm's own rule.
DEFAULT_OVERLAY_Q = 0.05
DEFAULT_MANUSCRIPT_Q = 0.001
OVERLAY_Q_SOURCE = (
    "chosen by this integration for display, matching the q <= 0.05 overlay of PR #173"
)
MANUSCRIPT_Q_SOURCE = (
    "the manuscript's q <= 0.001 rule, as the enrichment arm's committed run applies it"
)
#: Assessments that belong to other issues. Absent input means pending, not clear.
ASSESSMENTS = {
    "calibration_status": "#182",
    "sensitivity_status": "#183",
    "coverage_status": "#184",
    "orientation_status": "#188",
    "taxon_influence_status": "#190",
}


@dataclass(frozen=True)
class Record:
    """One evidence record. Every field means one thing and is never reused."""

    method: str
    method_version: str
    inference_level: str
    feature_type: str
    trait: str
    trait_qualifier: str
    reference: str
    stratum_kind: str
    stratum: str
    population: str
    denominator: str
    background: str
    eligibility_rule: str
    unique_studies: int | None
    paired_experiments: int | None
    profiles: int | None
    testable_units: int | None
    ties: int | None
    effect_value: float | None
    effect_units: str
    interval_low: float | None
    interval_high: float | None
    interval_type: str
    p_value: float | None
    q_value: float | None
    correction_method: str
    family_id: str
    family_size: int | None
    q_threshold: float | None
    effect_threshold: float | None
    test_status: str
    decision: str
    display_status: str
    zero_annotation_state: str
    #: The same row decided at the stricter threshold, so the manuscript's
    #: q <= 0.001 rule and the q <= 0.05 overlay can be compared per row rather
    #: than by rerunning the whole table (#294). Each threshold names its source,
    #: because no single cutoff is native to all three families here.
    strict_q_threshold: float | None = None
    strict_decision: str = UNDECIDED
    q_threshold_source: str = ""
    strict_q_threshold_source: str = ""
    selection_note: str = SELECTION_NOTE
    calibration_status: str = "pending #182"
    sensitivity_status: str = "pending #183"
    coverage_status: str = "pending #184"
    orientation_status: str = "pending #188"
    taxon_influence_status: str = "pending #190"
    source: str = ""
    notes: str = ""


def untested(record: Record) -> Record:
    """Blank every estimate on a row that was not tested.

    A source can carry a placeholder for an untested feature - a share of 0 and a
    q-value derived from p = 1 are what the fraction analysis reports for an all-tied
    trait before #165 - and a placeholder read as an effect is exactly the confusion
    this table exists to prevent.
    """
    if record.test_status == TESTED:
        return record
    return replace(
        record,
        effect_value=None,
        interval_low=None,
        interval_high=None,
        interval_type="unavailable",
        p_value=None,
        q_value=None,
        decision=UNDECIDED,
    )


def _decision(q_value: float | None, threshold: float | None, status: str) -> str:
    if status != TESTED or q_value is None or threshold is None:
        return UNDECIDED
    return SIGNIFICANT if q_value <= threshold else NOT_SIGNIFICANT


def fraction_sign_records(
    profiles: Sequence[Profile],
    *,
    source: str,
    q_threshold: float | None,
    strict_q_threshold: float | None,
) -> list[Record]:
    """The six Table 1 traits' study-level directional consistency test."""
    result = compute_direction_shifts(list(profiles))
    studies = len({p.study_id for p in profiles})
    experiments = len({(p.study_id, p.experiment_id) for p in profiles})
    records = []
    for shift in result.shifts:
        if shift.studies_scored:
            status = TESTED
        elif shift.studies_tied:
            status = ALL_TIED
        else:
            status = UNTESTED
        records.append(
            Record(
                method="directional_consistency_fraction_sign_test",
                method_version="table1-direction-shift",
                inference_level="corpus",
                feature_type="trait",
                trait=shift.trait,
                trait_qualifier=";".join(TABLE1_TRAITS[shift.trait]),
                reference="",
                stratum_kind="pooled",
                stratum="all body sites and conditions",
                population="human, paired experiments",
                denominator="summed non-MediaDive trait counts of the profile",
                background="",
                eligibility_rule="every study with a defined contrast; ties excluded from the test",
                unique_studies=studies,
                paired_experiments=experiments,
                profiles=len(profiles),
                testable_units=shift.studies_scored,
                ties=shift.studies_tied,
                effect_value=shift.share_up,
                effect_units="share of non-tied studies shifting toward the increased side",
                interval_low=shift.ci_low,
                interval_high=shift.ci_high,
                interval_type=(
                    "exact binomial 95% interval for the directional proportion"
                    if shift.share_up is not None
                    else "unavailable"
                ),
                p_value=shift.p_value,
                q_value=shift.q_value,
                correction_method="Benjamini-Hochberg",
                family_id="table1_six_fractions",
                family_size=len(result.shifts),
                q_threshold=q_threshold,
                strict_q_threshold=strict_q_threshold,
                q_threshold_source=OVERLAY_Q_SOURCE,
                strict_q_threshold_source=MANUSCRIPT_Q_SOURCE,
                effect_threshold=None,
                test_status=status,
                decision=_decision(shift.q_value, q_threshold, status),
                strict_decision=_decision(shift.q_value, strict_q_threshold, status),
                display_status=NOT_APPLICABLE,
                zero_annotation_state="profiles with no retained annotation score share 0",
                source=source,
            )
        )
    return [untested(record) for record in records]


def log_ratio_records(
    profiles: Sequence[Profile],
    *,
    source: str,
    q_threshold: float | None,
    strict_q_threshold: float | None,
) -> list[Record]:
    """The 15 paired log-ratio direction tests, each on its own eligible population."""
    result = compute_log_ratio_shifts(list(profiles))
    records = []
    for summary in result.summaries:
        status = {
            "tested": TESTED,
            "all_tied": ALL_TIED,
            "no_positive_pairs": ZERO_ANNOTATION,
        }.get(summary.status, UNTESTED)
        records.append(
            Record(
                method="directional_consistency_log_ratio_sign_test",
                method_version=LOG_RATIO_VERSION,
                inference_level="corpus",
                feature_type="trait_pair",
                trait=summary.numerator,
                trait_qualifier=";".join(TABLE1_TRAITS[summary.numerator]),
                reference=summary.reference,
                stratum_kind="pooled",
                stratum="all body sites and conditions",
                population="human, paired experiments with all four counts positive",
                denominator="none; the profile total cancels in the ratio",
                background="",
                eligibility_rule=(
                    "experiment excluded when either trait is zero on either side; no pseudocount"
                ),
                unique_studies=summary.studies_eligible,
                paired_experiments=summary.eligible_experiments,
                profiles=None,
                testable_units=summary.studies_scored,
                ties=summary.studies_tied,
                effect_value=summary.share_ratio_up,
                effect_units="share of non-tied eligible studies with a higher log2 ratio",
                interval_low=summary.share_ci_low,
                interval_high=summary.share_ci_high,
                interval_type=(
                    "exact binomial 95% interval for the directional proportion, not for the "
                    "mean log2 change"
                    if summary.share_ratio_up is not None
                    else "unavailable"
                ),
                p_value=summary.sign_p,
                q_value=summary.by_q,
                correction_method="Benjamini-Yekutieli",
                family_id="table1_fifteen_ratios",
                family_size=result.family_size,
                q_threshold=q_threshold,
                strict_q_threshold=strict_q_threshold,
                q_threshold_source=OVERLAY_Q_SOURCE,
                strict_q_threshold_source=MANUSCRIPT_Q_SOURCE,
                effect_threshold=None,
                test_status=status,
                decision=_decision(summary.by_q, q_threshold, status),
                strict_decision=_decision(summary.by_q, strict_q_threshold, status),
                display_status=NOT_APPLICABLE,
                zero_annotation_state=(
                    f"{summary.excluded_zero_experiments} experiments excluded for a zero count"
                ),
                source=source,
                notes=(
                    f"mean log2 change {summary.mean_log2_ratio_change:.4g}"
                    if summary.mean_log2_ratio_change is not None
                    else "mean log2 change undefined"
                ),
            )
        )
    return [untested(record) for record in records]


def guild_records(
    profiles: Sequence[Profile],
    *,
    source: str,
    q_threshold: float | None,
    strict_q_threshold: float | None,
) -> list[Record]:
    """The paired aerobe-against-anaerobe McNemar on study signs."""
    matrix = build_matrix(list(profiles), max_conditions=60, max_body_sites=8, min_studies=3)
    correlation = compute_correlation(list(profiles), matrix)
    up_down, down_up = correlation.discordant
    status = TESTED if up_down + down_up else INSUFFICIENT
    records = [
        Record(
            method="oxygen_guild_mcnemar",
            method_version="guild-maps-mcnemar",
            inference_level="corpus",
            feature_type="guild_pair",
            trait="aerobe",
            trait_qualifier=";".join(OXYGEN_GUILDS["aerobe"]),
            reference="anaerobe",
            stratum_kind="pooled",
            stratum="all body sites and conditions",
            population="human, paired experiments",
            denominator="summed non-MediaDive trait counts of the profile",
            background="",
            eligibility_rule="studies with a non-tied sign for both guilds",
            unique_studies=correlation.studies,
            paired_experiments=len({(p.study_id, p.experiment_id) for p in profiles}),
            profiles=len(profiles),
            testable_units=up_down + down_up,
            ties=None,
            effect_value=correlation.odds_ratio,
            effect_units="odds ratio of discordant study signs",
            interval_low=None,
            interval_high=None,
            interval_type="unavailable",
            p_value=correlation.mcnemar_p,
            q_value=correlation.mcnemar_p,
            correction_method="none; a single test",
            family_id="guild_mcnemar",
            family_size=1,
            q_threshold=q_threshold,
            strict_q_threshold=strict_q_threshold,
            q_threshold_source=OVERLAY_Q_SOURCE,
            strict_q_threshold_source=MANUSCRIPT_Q_SOURCE,
            effect_threshold=None,
            test_status=status,
            decision=_decision(correlation.mcnemar_p, q_threshold, status),
            strict_decision=_decision(correlation.mcnemar_p, strict_q_threshold, status),
            display_status=NOT_APPLICABLE,
            zero_annotation_state="profiles with no retained annotation score share 0",
            selection_note="the oxygen guilds were fixed a priori, not selected on this corpus",
            source=source,
        )
    ]
    return [untested(record) for record in records]


def placeholder_records(
    method: str, level: str, reason: str, features: Iterable[str]
) -> list[Record]:
    """Rows for a method whose output was not supplied.

    Without these the table would say nothing about the method, and a reader would
    have to guess whether it found nothing or was never run.
    """
    return [
        Record(
            method=method,
            method_version="",
            inference_level=level,
            feature_type="trait",
            trait=feature,
            trait_qualifier=";".join(TABLE1_TRAITS.get(feature, ())),
            reference="",
            stratum_kind="",
            stratum="",
            population="",
            denominator="",
            background="",
            eligibility_rule="",
            unique_studies=None,
            paired_experiments=None,
            profiles=None,
            testable_units=None,
            ties=None,
            effect_value=None,
            effect_units="",
            interval_low=None,
            interval_high=None,
            interval_type="unavailable",
            p_value=None,
            q_value=None,
            correction_method="",
            family_id="",
            family_size=None,
            q_threshold=None,
            effect_threshold=None,
            test_status=SOURCE_MISSING,
            decision=UNDECIDED,
            display_status=NOT_DRAWN,
            zero_annotation_state="unknown",
            notes=reason,
        )
        for feature in features
    ]


def read_status(path: Path | None, key: str, issue: str) -> str:
    """Read one assessment's status from its own output, or say it is pending."""
    if path is None or not path.exists():
        return f"pending {issue}"
    try:
        content = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return f"unreadable {issue}"
    if key == "calibration_status":
        verdicts = content.get("calibration_status", {})
        return (
            "; ".join(f"{name}: {verdict}" for name, verdict in sorted(verdicts.items()))
            or f"pending {issue}"
        )
    if key == "sensitivity_status":
        scenarios = content.get("scenario_status", {})
        return f"{len(scenarios)} scenarios reported" if scenarios else f"pending {issue}"
    if key == "orientation_status":
        counts = content.get("class_counts", [])
        kept = sum(
            int(row.get("experiments", 0))
            for row in counts
            if row.get("contrast_class") == "disease_vs_healthy"
        )
        return f"{kept} disease-vs-healthy experiments classified" if counts else f"pending {issue}"
    return content.get(key) or f"pending {issue}"


def apply_statuses(records: Sequence[Record], statuses: Mapping[str, str]) -> list[Record]:
    """Attach the assessment statuses to every row, so no row implies more than is known."""
    missing = set(ASSESSMENTS).difference(statuses)
    if missing:
        raise ValueError(f"no status supplied for {sorted(missing)}")
    return [
        replace(
            record,
            calibration_status=statuses["calibration_status"],
            sensitivity_status=statuses["sensitivity_status"],
            coverage_status=statuses["coverage_status"],
            orientation_status=statuses["orientation_status"],
            taxon_influence_status=statuses["taxon_influence_status"],
        )
        for record in records
    ]


def summarize(records: Sequence[Record]) -> list[dict[str, object]]:
    """Counts per method and state, so a reader can see what each method contributed."""
    rows: list[dict[str, object]] = []
    for method in sorted({r.method for r in records}):
        subset = [r for r in records if r.method == method]
        rows.append(
            {
                "method": method,
                "records": len(subset),
                "families": ";".join(sorted({r.family_id for r in subset if r.family_id})),
                "tested": sum(r.test_status == TESTED for r in subset),
                "untested_or_unsupplied": sum(r.test_status != TESTED for r in subset),
                "significant": sum(r.decision == SIGNIFICANT for r in subset),
                "not_significant": sum(r.decision == NOT_SIGNIFICANT for r in subset),
                "undecided": sum(r.decision == UNDECIDED for r in subset),
                "effect_units": ";".join(
                    sorted({r.effect_units for r in subset if r.effect_units})
                ),
                "combined_guarantee": "none; each family is corrected on its own",
            }
        )
    return rows


def _format(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def write_table(records: Sequence[Record], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    names = [field.name for field in fields(Record)]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(names)
        for record in records:
            row = asdict(record)
            writer.writerow([_format(row[name]) for name in names])
    return path


DICTIONARY: tuple[tuple[str, str], ...] = (
    ("method", "which procedure produced the row; never merged with another"),
    ("method_version", "the procedure's own version string"),
    ("inference_level", "what the result is about: corpus, body site, condition and site, profile"),
    ("feature_type", "trait, trait_pair or guild_pair"),
    ("trait / reference", "the feature, and the reference where a comparison names one"),
    ("trait_qualifier", "the raw identifiers behind the trait name"),
    ("population / eligibility_rule", "who is included, and the rule that decides"),
    ("denominator / background", "what the effect is relative to; empty when none applies"),
    (
        "unique_studies / paired_experiments / profiles / testable_units / ties",
        "distinct counts, never interchangeable",
    ),
    ("effect_value / effect_units", "the estimate and what it measures; units differ by method"),
    (
        "interval_low / interval_high / interval_type",
        "the interval and what it is an interval for; 'unavailable' where none was computed",
    ),
    (
        "p_value / q_value / correction_method / family_id / family_size",
        "the test and the family it was corrected in; families are not comparable",
    ),
    ("q_threshold / effect_threshold", "the thresholds this row's decision used, if any"),
    (
        "test_status",
        f"{TESTED}, {UNTESTED}, {INSUFFICIENT}, {ALL_TIED}, {ZERO_ANNOTATION}, "
        f"{NOT_SAMPLED} or {SOURCE_MISSING}",
    ),
    ("decision", f"{SIGNIFICANT}, {NOT_SIGNIFICANT} or {UNDECIDED}; only ever within one family"),
    ("display_status", f"{DRAWN}, {MASKED}, {NOT_DRAWN} or {NOT_APPLICABLE}; never a result"),
    ("zero_annotation_state", "how zero or absent annotation was treated"),
    ("selection_note", "whether the feature was selected on this corpus"),
    (
        "calibration_status / sensitivity_status / coverage_status / orientation_status / "
        "taxon_influence_status",
        "the state of #182, #183, #184, #188 and #190 for this row; 'pending' means not assessed",
    ),
    ("source", "the input the row was computed from"),
)


def write_view(
    records: Sequence[Record], summary: Sequence[Mapping[str, object]], path: Path
) -> Path:
    """A short Markdown view: the data dictionary, the per-method summary, the rows."""
    lines = [
        "# Integrated evidence table (generated)",
        "",
        f"Version {VERSION}. Every method keeps its own row, family and units. Nothing here",
        "combines them: there is no pooled p-value, no averaged effect and no false discovery",
        "rate for a set defined as significant under any method.",
        "",
        "## Data dictionary",
        "",
        "| Field | Meaning |",
        "| --- | --- |",
    ]
    lines += [f"| `{name}` | {meaning} |" for name, meaning in DICTIONARY]
    lines += [
        "",
        "## What each method contributed",
        "",
        "| Method | Records | Families | Tested | Not tested | Significant | Not significant "
        "| Undecided |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in summary:
        lines.append(
            f"| {row['method']} | {row['records']} | {row['families']} | {row['tested']} | "
            f"{row['untested_or_unsupplied']} | {row['significant']} | {row['not_significant']} | "
            f"{row['undecided']} |"
        )
    lines += [
        "",
        "## Records",
        "",
        "| Method | Feature | Level | Effect | Units | q (family) | Decision | Test state |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for record in records:
        feature = f"{record.trait}/{record.reference}" if record.reference else record.trait
        effect = "NA" if record.effect_value is None else f"{record.effect_value:.4g}"
        q_value = "NA" if record.q_value is None else f"{record.q_value:.3g}"
        lines.append(
            f"| {record.method} | {feature} | {record.inference_level} | {effect} | "
            f"{record.effect_units or 'NA'} | {q_value} ({record.family_id or 'NA'}) | "
            f"{record.decision} | {record.test_status} |"
        )
    lines += [
        "",
        "Reading rules: a blank q is not a negative result; `source_not_supplied` means the",
        "method's output was not given to this run; a decision holds only inside its own family;",
        "and the six Table 1 traits were selected on this corpus, so no row confirms them.",
        "",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _sha256(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile-input", type=Path, default=DEFAULT_PROFILE_PATH)
    parser.add_argument("--experiments-input", type=Path, default=DEFAULT_EXPERIMENTS_PATH)
    parser.add_argument("--full-dump-input", type=Path, default=DEFAULT_FULL_DUMP_PATH)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--q-threshold",
        type=float,
        default=DEFAULT_OVERLAY_Q,
        help="Display threshold for the primary decision column.",
    )
    parser.add_argument(
        "--strict-q-threshold",
        type=float,
        default=DEFAULT_MANUSCRIPT_Q,
        help="Second threshold decided per row, for the manuscript comparison (#294).",
    )
    parser.add_argument("--calibration-status", type=Path, default=DEFAULT_CALIBRATION_STATUS)
    parser.add_argument("--sensitivity-status", type=Path, default=DEFAULT_SENSITIVITY_STATUS)
    parser.add_argument("--orientation-provenance", type=Path, default=DEFAULT_ORIENTATION_STATUS)
    parser.add_argument("--coverage-status", type=Path)
    parser.add_argument("--influence-status", type=Path)
    parser.add_argument("--code-commit", help="Source revision, supplied by the Justfile recipe")
    args = parser.parse_args(argv)

    guilds = {**TABLE1_TRAITS, **OXYGEN_GUILDS}
    profiles, stats = load_profiles(
        args.profile_input, args.experiments_input, guilds, reject_missing_counts=True
    )
    vintage = read_bugsigdb_vintage(args.full_dump_input)
    write_stats_report(stats, vintage)
    source = f"{args.profile_input.name}@{_sha256(args.profile_input)[:12]}"

    # The shared comparison basis excludes a pair whose retained-annotation total
    # is zero on either side; scoring it as share 0 makes its contrast a tie and
    # merges zero-annotation studies into the all-tied count (#294, #178).
    exclusion = drop_zero_total_pairs(profiles)
    profiles = exclusion.profiles
    print(
        f"Excluded {exclusion.zero_total_experiments} experiments with a zero "
        f"retained-annotation total on a side; studies {exclusion.studies_before} "
        f"-> {exclusion.studies_after}"
    )

    thresholds = {"q_threshold": args.q_threshold, "strict_q_threshold": args.strict_q_threshold}
    records: list[Record] = []
    records += fraction_sign_records(profiles, source=source, **thresholds)
    records += log_ratio_records(profiles, source=source, **thresholds)
    records += guild_records(profiles, source=source, **thresholds)
    records += placeholder_records(
        "composition_arm_native_transformed_contrast",
        "body_site",
        "the composition arm's native output was not supplied to this run (#174, #197)",
        TABLE1_TRAITS,
    )
    records += placeholder_records(
        "enrichment_arm_profile_hypergeometric",
        "profile",
        "the enrichment arm's native output was not supplied to this run (#131)",
        TABLE1_TRAITS,
    )
    records += placeholder_records(
        "aggregated_profile_support",
        "condition_body_site",
        "the cell support rule lives on the #173 branch and was not supplied (#177)",
        TABLE1_TRAITS,
    )

    statuses = {
        "calibration_status": read_status(args.calibration_status, "calibration_status", "#182"),
        "sensitivity_status": read_status(args.sensitivity_status, "sensitivity_status", "#183"),
        "coverage_status": read_status(args.coverage_status, "coverage_status", "#184"),
        "orientation_status": read_status(
            args.orientation_provenance, "orientation_status", "#188"
        ),
        "taxon_influence_status": read_status(
            args.influence_status, "taxon_influence_status", "#190"
        ),
    }
    records = apply_statuses(records, statuses)
    summary = summarize(records)

    output = args.output_dir
    written = [
        write_table(records, output / "integrated_evidence.tsv"),
        write_view(records, summary, output / "integrated_evidence.md"),
    ]
    provenance = {
        "version": VERSION,
        "bugsigdb_vintage": vintage,
        "code_commit": args.code_commit,
        "profile_sha256": _sha256(args.profile_input),
        "experiments_sha256": _sha256(args.experiments_input),
        "analysis_source_sha256": _sha256(Path(__file__)),
        "corpus": asdict(stats),
        "q_threshold": args.q_threshold,
        "assessment_status": statuses,
        "methods": [dict(row) for row in summary],
        "no_combined_guarantee": (
            "each family is corrected on its own; a union of decisions has no false discovery rate"
        ),
        "selection_note": SELECTION_NOTE,
    }
    provenance_path = output / "integrated_evidence_provenance.json"
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    written.append(provenance_path)
    for row in summary:
        print(row)
    for path in written:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
