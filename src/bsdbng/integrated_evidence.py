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

from bsdbng.guild_maps import (
    CELL_BELOW_MIN_STUDIES,
    CELL_NO_SIGNAL,
    CELL_NOT_SAMPLED,
    CELL_VALUE,
    GuildMatrix,
    build_matrix,
    cell_status_rows,
    compute_correlation,
)
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
    format_pair_exclusion,
    load_profiles,
    read_bugsigdb_vintage,
    resolve_trait_columns,
    write_stats_report,
)

#: Declared term ids and the raw matrix columns they resolved to, per feature.
#: The raw column is the only key that joins a result back to the matrix: two of
#: the six Table 1 traits are carried as METPO-qualified composites (#294).
FeatureColumns = Mapping[str, Sequence[tuple[str, str]]]

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
#: The method's eligibility rule left nothing to test. Not the same claim as
#: ``ZERO_ANNOTATION``: two traits can both be well annotated across the corpus
#: and still never both be positive on both sides of one experiment (#294).
NO_ELIGIBLE_PAIRS = "no_eligible_pairs"
#: Sampled, but no taxon in the profile carried any term of the feature.
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
#: The guild rows carried "fixed a priori", which nothing in this repository
#: establishes and which is not true of the aerobe guild in any case:
#: METPO:1000608 is both a member of it and one of the six traits the
#: composition arm selected on this corpus (#294).
GUILD_SELECTION_NOTE = (
    "the two oxygen guilds are the groupings both arms use, with no record in this "
    "repository of where or when they were declared; the aerobe guild contains "
    "METPO:1000608, which is also one of the six traits selected on this corpus"
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
#: The oxygen guild figure's aggregation parameters. The McNemar row and the cell
#: support inventory share one matrix, so both describe the same figure and the
#: masking threshold a cell row reports is the one that masked it.
GUILD_MAX_CONDITIONS = 60
GUILD_MAX_BODY_SITES = 8
GUILD_MIN_STUDIES = 3

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
    #: The declared term ids behind ``trait``, and the raw matrix columns they
    #: resolved to. The declared id is not a join key: ``GO:0030245`` is carried
    #: in the matrix as ``METPO:2000103|GO:0030245`` (#294).
    trait_qualifier: str
    trait_columns: str
    reference: str
    #: The same two identifiers for the reference of a comparison. Dropping them
    #: left a ratio row naming only half of its own hypothesis.
    reference_qualifier: str
    reference_columns: str
    stratum_kind: str
    stratum: str
    #: The stratum as separate identifiers rather than one packed label, so a
    #: cell row can be joined on any of them. Blank where the row is pooled.
    condition: str
    body_site: str
    direction: str
    population: str
    denominator: str
    background: str
    eligibility_rule: str
    unique_studies: int | None
    #: The ratio population, which is conditional and differs per hypothesis:
    #: how many studies the method could have used, how many its eligibility rule
    #: removed, and how many contributed some but not all of their experiments.
    studies_total: int | None
    studies_excluded: int | None
    studies_partially_observed: int | None
    paired_experiments: int | None
    experiments_eligible: int | None
    experiments_excluded_zero_count: int | None
    profiles: int | None
    testable_units: int | None
    ties: int | None
    effect_value: float | None
    effect_units: str
    #: The distribution of the per-study values the test was run over, which is a
    #: description of the input rather than the tested estimate. It is kept on an
    #: untested row, where ``effect_value`` is blanked, because "every study tied"
    #: and "no study was eligible" are different facts about that distribution.
    study_effect_mean: float | None
    study_effect_median: float | None
    study_effect_units: str
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
    #: The minimum number of studies a figure cell needs before it is drawn
    #: rather than masked. It is a display rule, not an inference threshold.
    support_threshold: int | None
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


def is_cell_row(record: Record) -> bool:
    """Whether the row describes a figure cell rather than a method result.

    ``not_drawn`` is shared with a method whose output was not supplied, so the
    display state cannot identify the inventory; the support threshold can,
    because only a cell is subject to one.
    """
    return record.support_threshold is not None


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


def declared_terms(columns: FeatureColumns, feature: str) -> str:
    """The declared term ids of a feature, as the analysis asked for them."""
    return ";".join(term for term, _ in columns.get(feature, ()))


def matrix_columns(columns: FeatureColumns, feature: str) -> str:
    """The raw matrix columns a feature resolved to, in header order.

    This is the join key. ``trait_qualifier`` is what was asked for; this is what
    was found, and the two differ for every METPO-qualified composite (#294).
    """
    return ";".join(column for _, column in columns.get(feature, ()))


def fraction_sign_records(
    profiles: Sequence[Profile],
    *,
    source: str,
    columns: FeatureColumns,
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
                trait_qualifier=declared_terms(columns, shift.trait),
                trait_columns=matrix_columns(columns, shift.trait),
                reference="",
                reference_qualifier="",
                reference_columns="",
                stratum_kind="pooled",
                stratum="all body sites and conditions",
                condition="",
                body_site="",
                direction="",
                population="human, paired experiments",
                denominator="summed non-MediaDive trait counts of the profile",
                background="",
                eligibility_rule="every study with a defined contrast; ties excluded from the test",
                unique_studies=studies,
                studies_total=studies,
                studies_excluded=0,
                studies_partially_observed=None,
                paired_experiments=experiments,
                experiments_eligible=experiments,
                experiments_excluded_zero_count=0,
                profiles=len(profiles),
                testable_units=shift.studies_scored,
                ties=shift.studies_tied,
                effect_value=shift.share_up,
                effect_units="share of non-tied studies shifting toward the increased side",
                # The source run publishes the per-study contrast tallies, not a
                # summary of the contrast distribution, so none is invented here.
                study_effect_mean=None,
                study_effect_median=None,
                study_effect_units="",
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
                support_threshold=None,
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
    columns: FeatureColumns,
    q_threshold: float | None,
    strict_q_threshold: float | None,
) -> list[Record]:
    """The 15 paired log-ratio direction tests, each on its own eligible population."""
    result = compute_log_ratio_shifts(list(profiles))
    records = []
    for summary in result.summaries:
        # ``no_positive_pairs`` means the pair's eligibility rule left nothing to
        # test, which is not the same claim as "the annotation is absent": the
        # traits can be well annotated across the corpus and still never both be
        # positive on both sides of one experiment (#294).
        status = {
            "tested": TESTED,
            "all_tied": ALL_TIED,
            "no_positive_pairs": NO_ELIGIBLE_PAIRS,
        }.get(summary.status, UNTESTED)
        records.append(
            Record(
                method="directional_consistency_log_ratio_sign_test",
                method_version=LOG_RATIO_VERSION,
                inference_level="corpus",
                feature_type="trait_pair",
                trait=summary.numerator,
                trait_qualifier=declared_terms(columns, summary.numerator),
                trait_columns=matrix_columns(columns, summary.numerator),
                reference=summary.reference,
                reference_qualifier=declared_terms(columns, summary.reference),
                reference_columns=matrix_columns(columns, summary.reference),
                stratum_kind="pooled",
                stratum="all body sites and conditions",
                condition="",
                body_site="",
                direction="",
                population="human, paired experiments with all four counts positive",
                denominator="none; the profile total cancels in the ratio",
                background="",
                eligibility_rule=(
                    "experiment excluded when either trait is zero on either side; no pseudocount"
                ),
                unique_studies=summary.studies_eligible,
                studies_total=summary.studies_total,
                studies_excluded=summary.studies_excluded,
                studies_partially_observed=summary.studies_partially_observed,
                paired_experiments=summary.paired_experiments,
                experiments_eligible=summary.eligible_experiments,
                experiments_excluded_zero_count=summary.excluded_zero_experiments,
                profiles=None,
                testable_units=summary.studies_scored,
                ties=summary.studies_tied,
                effect_value=summary.share_ratio_up,
                effect_units="share of non-tied eligible studies with a higher log2 ratio",
                study_effect_mean=summary.mean_log2_ratio_change,
                study_effect_median=summary.median_log2_ratio_change,
                study_effect_units="per-study mean log2 change in the numerator/reference ratio",
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
                support_threshold=None,
                test_status=status,
                decision=_decision(summary.by_q, q_threshold, status),
                strict_decision=_decision(summary.by_q, strict_q_threshold, status),
                display_status=NOT_APPLICABLE,
                # The count lives in experiments_excluded_zero_count; this names
                # the rule, so no column packs a label and a number together.
                zero_annotation_state="whole experiment excluded on a zero count; no pseudocount",
                source=source,
            )
        )
    return [untested(record) for record in records]


def build_guild_matrix(profiles: Sequence[Profile]) -> GuildMatrix:
    """The oxygen guild figure's aggregation, at the figure's own parameters."""
    return build_matrix(
        list(profiles),
        max_conditions=GUILD_MAX_CONDITIONS,
        max_body_sites=GUILD_MAX_BODY_SITES,
        min_studies=GUILD_MIN_STUDIES,
    )


def guild_records(
    profiles: Sequence[Profile],
    *,
    matrix: GuildMatrix,
    source: str,
    columns: FeatureColumns,
    q_threshold: float | None,
    strict_q_threshold: float | None,
) -> list[Record]:
    """The paired aerobe-against-anaerobe McNemar on study signs.

    The test is uncorrected, so the row carries a p and leaves q blank. A
    threshold decision in this table is a decision on a corrected value, and
    there is no corrected value here; writing the uncorrected p into ``q_value``
    made a single uncorrected test look like a one-member FDR family (#294).
    """
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
            trait_qualifier=declared_terms(columns, "aerobe"),
            trait_columns=matrix_columns(columns, "aerobe"),
            reference="anaerobe",
            reference_qualifier=declared_terms(columns, "anaerobe"),
            reference_columns=matrix_columns(columns, "anaerobe"),
            stratum_kind="pooled",
            stratum="all body sites and conditions",
            condition="",
            body_site="",
            direction="",
            population="human, paired experiments",
            denominator="summed non-MediaDive trait counts of the profile",
            background="",
            eligibility_rule="studies with a non-tied sign for both guilds",
            unique_studies=correlation.studies,
            studies_total=correlation.studies,
            studies_excluded=None,
            studies_partially_observed=None,
            paired_experiments=len({(p.study_id, p.experiment_id) for p in profiles}),
            experiments_eligible=None,
            experiments_excluded_zero_count=0,
            profiles=len(profiles),
            testable_units=up_down + down_up,
            ties=None,
            effect_value=correlation.odds_ratio,
            effect_units="odds ratio of discordant study signs",
            study_effect_mean=None,
            study_effect_median=None,
            study_effect_units="",
            interval_low=None,
            interval_high=None,
            interval_type="unavailable",
            p_value=correlation.mcnemar_p,
            q_value=None,
            correction_method="none; a single uncorrected test, so no q is computed",
            family_id="guild_mcnemar",
            family_size=1,
            q_threshold=q_threshold,
            strict_q_threshold=strict_q_threshold,
            q_threshold_source=OVERLAY_Q_SOURCE,
            strict_q_threshold_source=MANUSCRIPT_Q_SOURCE,
            effect_threshold=None,
            support_threshold=None,
            test_status=status,
            decision=_decision(None, q_threshold, status),
            strict_decision=_decision(None, strict_q_threshold, status),
            display_status=NOT_APPLICABLE,
            zero_annotation_state="profiles with no retained annotation score share 0",
            selection_note=GUILD_SELECTION_NOTE,
            source=source,
        )
    ]
    return [untested(record) for record in records]


#: How a matrix cell's render state maps onto the table's own two state columns.
#: A cell is never tested here - the tests in this table are at corpus level - so
#: its test_status reports the support it has and its display_status reports what
#: the figure does with it. They are separate columns because a cell held back
#: for too little support is not a negative result (#294).
CELL_STATES: dict[str, tuple[str, str]] = {
    CELL_VALUE: (UNTESTED, DRAWN),
    CELL_NO_SIGNAL: (ZERO_ANNOTATION, DRAWN),
    CELL_BELOW_MIN_STUDIES: (INSUFFICIENT, MASKED),
    CELL_NOT_SAMPLED: (NOT_SAMPLED, NOT_DRAWN),
}


def cell_support_records(
    matrix: GuildMatrix,
    *,
    source: str,
    columns: FeatureColumns,
    min_studies: int,
    method: str = "oxygen_guild_cell_support",
) -> list[Record]:
    """One row per figure cell, saying why the cell is drawn or is not.

    ``drawn``, ``masked_below_min_studies`` and ``not_sampled`` only exist at
    cell level: a corpus-level test is not drawn anywhere, so its display state
    is ``not_applicable``, and before these rows existed no code path could set
    the other three at all - the enum named states the table never emitted.

    The classification is ``guild_maps.cell_status``, the same function the
    figure writer uses, so a cell cannot be reported as drawn by one artifact
    and masked by the other. The inventory is complete, every cell of the
    matrix including the unsampled ones, because a partial inventory would put
    a reader back to guessing what an absent row meant.
    """
    statuses = cell_status_rows(matrix)
    records: list[Record] = []
    for row_index, condition in enumerate(matrix.conditions):
        for column_index, column in enumerate(matrix.columns):
            site, suffix = column.split("|", 1)
            guild, direction = suffix.rsplit("_", 1)
            test_status, display_status = CELL_STATES[statuses[row_index][column_index]]
            records.append(
                Record(
                    method=method,
                    method_version="guild-maps-build-matrix",
                    inference_level="condition_body_site",
                    feature_type="guild",
                    trait=guild,
                    trait_qualifier=declared_terms(columns, guild),
                    trait_columns=matrix_columns(columns, guild),
                    reference="",
                    reference_qualifier="",
                    reference_columns="",
                    stratum_kind="condition x body site x direction",
                    stratum=f"{condition}, {site}, {direction}",
                    condition=condition,
                    body_site=site,
                    direction=direction,
                    population="human, paired experiments",
                    denominator="summed non-MediaDive trait counts of the profile",
                    background="",
                    eligibility_rule=(
                        "no test is run on a cell; the cell is drawn when enough studies "
                        "back it and any of them carries the annotation"
                    ),
                    unique_studies=matrix.study_counts[row_index][column_index],
                    studies_total=None,
                    studies_excluded=None,
                    studies_partially_observed=None,
                    paired_experiments=None,
                    experiments_eligible=None,
                    experiments_excluded_zero_count=None,
                    profiles=None,
                    testable_units=None,
                    ties=None,
                    effect_value=None,
                    effect_units="",
                    study_effect_mean=None,
                    study_effect_median=None,
                    study_effect_units="",
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
                    support_threshold=min_studies,
                    test_status=test_status,
                    decision=UNDECIDED,
                    display_status=display_status,
                    zero_annotation_state=(
                        "sampled, and no taxon carried any term of the guild"
                        if test_status == ZERO_ANNOTATION
                        else "profiles with no retained annotation score share 0"
                    ),
                    selection_note=GUILD_SELECTION_NOTE,
                    source=source,
                )
            )
    return records


def placeholder_records(
    method: str,
    level: str,
    reason: str,
    features: Iterable[str],
    *,
    columns: FeatureColumns | None = None,
) -> list[Record]:
    """Rows for a method whose output was not supplied.

    Without these the table would say nothing about the method, and a reader would
    have to guess whether it found nothing or was never run. The rows cover the
    six Table 1 traits only: they are not an inventory of what the method would
    report, which for the enrichment arm is every retained trait in the matrix.
    """
    resolved: FeatureColumns = {} if columns is None else columns
    return [
        Record(
            method=method,
            method_version="",
            inference_level=level,
            feature_type="trait",
            trait=feature,
            trait_qualifier=declared_terms(resolved, feature),
            trait_columns=matrix_columns(resolved, feature),
            reference="",
            reference_qualifier="",
            reference_columns="",
            stratum_kind="",
            stratum="",
            condition="",
            body_site="",
            direction="",
            population="",
            denominator="",
            background="",
            eligibility_rule="",
            unique_studies=None,
            studies_total=None,
            studies_excluded=None,
            studies_partially_observed=None,
            paired_experiments=None,
            experiments_eligible=None,
            experiments_excluded_zero_count=None,
            profiles=None,
            testable_units=None,
            ties=None,
            effect_value=None,
            effect_units="",
            study_effect_mean=None,
            study_effect_median=None,
            study_effect_units="",
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
            support_threshold=None,
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


def family_inventory(records: Sequence[Record]) -> list[dict[str, object]]:
    """Every correction family in the table, with its declared size (#236, #289).

    A q-value means nothing without the family it was corrected in and how large
    that family was, so each family is exported with both, and the declared size
    is checked against the rows actually present. A mismatch means the exported
    table is not the family the correction was computed over, which would make
    every q in it unreadable - so it is an error, not a note.
    """
    rows: list[dict[str, object]] = []
    for family in sorted({r.family_id for r in records if r.family_id}):
        subset = [r for r in records if r.family_id == family]
        declared = {r.family_size for r in subset}
        if len(declared) != 1:
            msg = f"family {family!r} reports more than one size: {sorted(map(str, declared))}"
            raise ValueError(msg)
        size = subset[0].family_size
        if size != len(subset):
            msg = (
                f"family {family!r} declares {size} members but the table exports "
                f"{len(subset)} rows for it; a q corrected over a different set of "
                "tests than the one exported cannot be read"
            )
            raise ValueError(msg)
        rows.append(
            {
                "family_id": family,
                "family_size": size,
                "correction_method": ";".join(sorted({r.correction_method for r in subset})),
                "methods": ";".join(sorted({r.method for r in subset})),
                "rows_exported": len(subset),
                "tested": sum(r.test_status == TESTED for r in subset),
                "with_q_value": sum(r.q_value is not None for r in subset),
                "significant_at_q_threshold": sum(r.decision == SIGNIFICANT for r in subset),
                "significant_at_strict_q_threshold": sum(
                    r.strict_decision == SIGNIFICANT for r in subset
                ),
            }
        )
    unfamilied = [r for r in records if r.q_value is not None and not r.family_id]
    if unfamilied:
        detail = "; ".join(f"{r.method}/{r.trait}" for r in unfamilied)
        msg = f"rows carry a q-value with no correction family: {detail}"
        raise ValueError(msg)
    return rows


def feature_column_rows(records: Sequence[Record], columns: FeatureColumns) -> list[dict[str, str]]:
    """One tidy row per (feature, declared term, raw matrix column).

    The main table lists a feature's columns in one cell, which is not atomic
    when a guild pools several terms. This is the joinable form: one value per
    cell, keyed on the raw column name the matrix actually carries (#294).
    """
    used = {r.trait for r in records if r.trait} | {r.reference for r in records if r.reference}
    return [
        {"feature": feature, "term_id": term, "matrix_column": column}
        for feature in sorted(used)
        for term, column in columns.get(feature, ())
    ]


def unmatched_terms(columns: FeatureColumns, declared: Mapping[str, Sequence[str]]) -> list[str]:
    """Declared terms that no matrix column carries.

    A guild with no column at all is rejected by the loader, but one declared
    term of several can go missing without any error, and the guild sum would
    then quietly cover fewer terms than it names.
    """
    return sorted(
        f"{feature}/{term}"
        for feature, terms in declared.items()
        for term in terms
        if term not in {found for found, _ in columns.get(feature, ())}
    )


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


def write_tidy_tsv(
    rows: Sequence[Mapping[str, object]], fieldnames: Sequence[str], path: Path
) -> Path:
    """Write identifier columns plus one value per cell, header always present."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(list(fieldnames))
        for row in rows:
            writer.writerow([_format(row.get(name)) for name in fieldnames])
    return path


DICTIONARY: tuple[tuple[str, str], ...] = (
    ("method", "which procedure produced the row; never merged with another"),
    ("method_version", "the procedure's own version string"),
    ("inference_level", "what the result is about: corpus, body site, condition and site, profile"),
    ("feature_type", "trait, trait_pair or guild_pair"),
    ("trait / reference", "the feature, and the reference where a comparison names one"),
    (
        "trait_qualifier / reference_qualifier",
        "the term ids the analysis declared for each side of the hypothesis",
    ),
    (
        "trait_columns / reference_columns",
        "the raw matrix columns those terms resolved to; this is the join key, and it "
        "differs from the declared term for every METPO-qualified composite",
    ),
    ("population / eligibility_rule", "who is included, and the rule that decides"),
    ("denominator / background", "what the effect is relative to; empty when none applies"),
    (
        "unique_studies / paired_experiments / profiles / testable_units / ties",
        "distinct counts, never interchangeable",
    ),
    (
        "studies_total / studies_excluded / studies_partially_observed",
        "the conditional population: studies the method could have used, studies its "
        "eligibility rule removed, and studies contributing only some experiments. "
        "Blank where the method defines no such rule",
    ),
    (
        "experiments_eligible / experiments_excluded_zero_count",
        "experiments the method's own rule kept and removed for a zero count; the "
        "corpus-wide zero-total pair exclusion is in the provenance record, not here",
    ),
    ("effect_value / effect_units", "the estimate and what it measures; units differ by method"),
    (
        "study_effect_mean / study_effect_median / study_effect_units",
        "the distribution of per-study values the test ran over, which describes the "
        "input rather than the tested estimate, so it survives on an untested row",
    ),
    (
        "interval_low / interval_high / interval_type",
        "the interval and what it is an interval for; 'unavailable' where none was computed",
    ),
    (
        "p_value / q_value / correction_method / family_id / family_size",
        "the test and the family it was corrected in; families are not comparable",
    ),
    (
        "q_threshold / effect_threshold / support_threshold",
        "the thresholds this row used: two for the decision, and the minimum studies "
        "a figure cell needs to be drawn. The support threshold is a display rule and "
        "is set only on a cell row, which is what identifies the cell inventory",
    ),
    (
        "stratum_kind / stratum / condition / body_site / direction",
        "what the row is about within the corpus, as one label and as separate "
        "identifiers; blank on a pooled row",
    ),
    (
        "test_status",
        f"{TESTED}, {UNTESTED}, {INSUFFICIENT}, {ALL_TIED}, {NO_ELIGIBLE_PAIRS}, "
        f"{ZERO_ANNOTATION}, {NOT_SAMPLED} or {SOURCE_MISSING}. The last three are "
        "different claims: nothing was eligible to test, the feature was sampled "
        "but carried no annotation, and nothing was sampled at all",
    ),
    ("decision", f"{SIGNIFICANT}, {NOT_SIGNIFICANT} or {UNDECIDED}; only ever within one family"),
    (
        "display_status",
        f"{DRAWN}, {MASKED}, {NOT_DRAWN} or {NOT_APPLICABLE}; never a result. The first "
        "two belong to a figure cell, which is what the cell rows carry; a corpus-level "
        f"test is drawn nowhere, so it is {NOT_APPLICABLE}",
    ),
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
    cells = [record for record in records if is_cell_row(record)]
    if cells:
        lines += [
            "",
            "## Figure cells, by what the figure does with them",
            "",
            "Cell rows report support and display, never a test. They are listed here as",
            "counts rather than one row each; the whole inventory is in the TSV.",
            "",
            "| Display | Test state | Cells |",
            "| --- | --- | --- |",
        ]
        states = sorted({(record.display_status, record.test_status) for record in cells})
        for display, test_status in states:
            count = sum(
                record.display_status == display and record.test_status == test_status
                for record in cells
            )
            lines.append(f"| {display} | {test_status} | {count} |")

    lines += [
        "",
        "## Records",
        "",
        "| Method | Feature | Level | Effect | Units | p | q (family) | Decision | Test state |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for record in records:
        if is_cell_row(record):
            continue
        feature = f"{record.trait}/{record.reference}" if record.reference else record.trait
        effect = "NA" if record.effect_value is None else f"{record.effect_value:.4g}"
        p_value = "NA" if record.p_value is None else f"{record.p_value:.3g}"
        q_value = "NA" if record.q_value is None else f"{record.q_value:.3g}"
        lines.append(
            f"| {record.method} | {feature} | {record.inference_level} | {effect} | "
            f"{record.effect_units or 'NA'} | {p_value} | {q_value} "
            f"({record.family_id or 'NA'}) | {record.decision} | {record.test_status} |"
        )
    lines += [
        "",
        "Reading rules: a blank q is not a negative result; `source_not_supplied` means the",
        "method's output was not given to this run; a decision holds only inside its own family;",
        "and the six Table 1 traits were selected on this corpus, so no row confirms them.",
        "A row with a p and no q was not corrected for multiplicity, so it carries no",
        "threshold decision either.",
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

    # The declared term is not a join key: GO:0030245 is carried in the matrix as
    # METPO:2000103|GO:0030245, so every row exports the raw column too (#294).
    columns = resolve_trait_columns(args.profile_input, guilds)
    unmatched = unmatched_terms(columns, guilds)
    if unmatched:
        print(f"Declared terms with no matrix column: {', '.join(unmatched)}")

    # The shared comparison basis excludes a pair whose retained-annotation total
    # is zero on either side; scoring it as share 0 makes its contrast a tie and
    # merges zero-annotation studies into the all-tied count (#294, #178).
    exclusion = drop_zero_total_pairs(profiles)
    profiles = exclusion.profiles
    print(format_pair_exclusion(exclusion))

    thresholds = {"q_threshold": args.q_threshold, "strict_q_threshold": args.strict_q_threshold}
    records: list[Record] = []
    records += fraction_sign_records(profiles, source=source, columns=columns, **thresholds)
    records += log_ratio_records(profiles, source=source, columns=columns, **thresholds)
    matrix = build_guild_matrix(profiles)
    records += guild_records(profiles, matrix=matrix, source=source, columns=columns, **thresholds)
    records += cell_support_records(
        matrix, source=source, columns=columns, min_studies=GUILD_MIN_STUDIES
    )
    records += placeholder_records(
        "composition_arm_native_transformed_contrast",
        "body_site",
        "the composition arm's native output was not supplied to this run (#174, #197)",
        TABLE1_TRAITS,
        columns=columns,
    )
    records += placeholder_records(
        "enrichment_arm_profile_hypergeometric",
        "profile",
        "the enrichment arm's native output was not supplied to this run (#131)",
        TABLE1_TRAITS,
        columns=columns,
    )
    records += placeholder_records(
        "aggregated_profile_support",
        "condition_body_site",
        "the cell support rule lives on the #173 branch and was not supplied (#177)",
        TABLE1_TRAITS,
        columns=columns,
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
    families = family_inventory(records)
    feature_rows = feature_column_rows(records, columns)

    output = args.output_dir
    written = [
        write_table(records, output / "integrated_evidence.tsv"),
        write_tidy_tsv(
            feature_rows,
            ("feature", "term_id", "matrix_column"),
            output / "integrated_evidence_features.tsv",
        ),
        write_tidy_tsv(
            families,
            (
                "family_id",
                "family_size",
                "correction_method",
                "methods",
                "rows_exported",
                "tested",
                "with_q_value",
                "significant_at_q_threshold",
                "significant_at_strict_q_threshold",
            ),
            output / "integrated_evidence_families.tsv",
        ),
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
        # The exclusion that defines the shared comparison basis, recorded here
        # as well as printed, so a reader of the files alone can see it (#294).
        "pair_exclusion": {
            "zero_total_experiments": exclusion.zero_total_experiments,
            "studies_before": exclusion.studies_before,
            "studies_after": exclusion.studies_after,
            "studies_lost": exclusion.studies_lost,
            "experiments_kept": len({(p.study_id, p.experiment_id) for p in exclusion.profiles}),
            "rule": (
                "a whole experiment pair is dropped when its retained-annotation total "
                "is zero on either side; the share is undefined, not zero, and scoring "
                "it as zero would count the pair as a tie"
            ),
        },
        "thresholds": {
            "q_threshold": args.q_threshold,
            "q_threshold_source": OVERLAY_Q_SOURCE,
            "strict_q_threshold": args.strict_q_threshold,
            "strict_q_threshold_source": MANUSCRIPT_Q_SOURCE,
        },
        "assessment_status": statuses,
        "methods": [dict(row) for row in summary],
        "families": families,
        "feature_columns": {
            feature: [list(pair) for pair in pairs] for feature, pairs in sorted(columns.items())
        },
        "declared_terms_without_a_matrix_column": unmatched,
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
