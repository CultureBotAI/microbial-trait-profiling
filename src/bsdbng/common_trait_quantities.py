"""Common study-weighted fraction and log-ratio quantities for issue #256.

This module implements the descriptive contract in ``three-framework-comparison-v1``.
It deliberately does not add a hypothesis test or correction family.  Full-population
fractions, fractions restricted to each ratio's complete-case population, and the
ratio on that identical population are separate rows so that a denominator change,
a population change, and a change of estimand cannot be silently conflated.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from collections.abc import Callable, Mapping, Sequence
from dataclasses import asdict, dataclass, fields
from decimal import Decimal, InvalidOperation
from itertools import combinations
from pathlib import Path
from statistics import fmean, median

from scipy.stats import binomtest

from bsdbng.paths import DERIVED_DATA_DIR
from bsdbng.trait_composition import _matched_profiles
from bsdbng.trait_guilds import (
    DEFAULT_EXPERIMENTS_PATH,
    DEFAULT_FULL_DUMP_PATH,
    DEFAULT_PROFILE_PATH,
    MEDIADIVE_MARKER,
    PROFILE_ID_COLUMNS,
    TABLE1_TRAIT_LABELS,
    TABLE1_TRAITS,
    Profile,
    column_term_id,
    load_profiles,
    read_bugsigdb_vintage,
    write_stats_report,
)

VERSION = "common-study-weighted-trait-quantities-v1"
COMPARISON_SPEC = "three-framework-comparison-v1"
DENOMINATOR_DEFINITION = "sum of all retained non-MediaDive trait columns before selection"
FRACTION_TIE_TOLERANCE = 1e-12
LOG_RATIO_TIE_TOLERANCE = 1e-12
CONFIDENCE_LEVEL = 0.95
DEFAULT_OUTPUT_DIR = DERIVED_DATA_DIR / "table1_common_quantities"

FRACTION_FULL = "fraction_contrast_full_population"
FRACTION_MATCHED = "fraction_contrast_ratio_matched_population"
LOG_RATIO = "paired_log2_ratio_contrast"


@dataclass(frozen=True)
class Pair:
    """One validated increased/decreased experiment pair."""

    study_id: str
    experiment_id: str
    condition: str
    body_site: str
    inc: Profile
    dec: Profile

    @property
    def pair_id(self) -> str:
        return f"{self.study_id}|{self.experiment_id}"

    def positive(self, traits: Sequence[str]) -> bool:
        return all(self.inc.guild_counts[t] > 0 and self.dec.guild_counts[t] > 0 for t in traits)


@dataclass(frozen=True)
class Population:
    """An exact eligible-pair population used by one or more quantities."""

    population_id: str
    population_type: str
    numerator: str
    reference: str
    eligibility_rule: str
    pairs: tuple[Pair, ...]


@dataclass(frozen=True)
class QuantitySummary:
    """One common descriptive estimate and its population/uncertainty contract."""

    row_id: str
    comparison_spec: str
    quantity_type: str
    numerator: str
    reference: str
    numerator_source_column: str
    reference_source_column: str
    population_id: str
    population_numerator: str
    population_reference: str
    population_definition: str
    denominator_definition: str
    eligibility_rule: str
    independent_unit: str
    study_weighting: str
    base_unique_studies: int
    eligible_unique_studies: int
    entirely_excluded_studies: int
    partially_observed_studies: int
    base_experiment_pairs: int
    eligible_experiment_pairs: int
    excluded_zero_total_pairs: int
    excluded_zero_named_count_pairs: int
    excluded_missing_or_malformed_pairs: int
    positive_studies: int
    negative_studies: int
    tied_studies: int
    non_tied_studies: int
    mean_study_effect: float | None
    median_study_effect: float | None
    effect_units: str
    proportion_positive_non_tied: float | None
    proportion_ci_low: float | None
    proportion_ci_high: float | None
    uncertainty_parameter: str
    mean_effect_uncertainty: str
    status: str
    calibration_status: str


@dataclass(frozen=True)
class StudyEffect:
    """Coverage and one equally weighted effect for a study/quantity row."""

    row_id: str
    population_id: str
    quantity_type: str
    numerator: str
    reference: str
    study_id: str
    base_experiment_pairs: int
    eligible_experiment_pairs: int
    excluded_experiment_pairs: int
    mean_study_effect: float | None
    effect_units: str
    direction_status: str


@dataclass(frozen=True)
class TraitIdentity:
    """A declared Table 1 trait resolved to one exact input column."""

    order: int
    trait: str
    display_label: str
    declared_term_ids: str
    source_column: str
    source_term_id: str


@dataclass(frozen=True)
class DenominatorAudit:
    """Comparison of the shared loader with PR #197's committed parity artifact."""

    reference_revision: str
    summary: tuple[dict[str, object], ...]
    mismatches: tuple[dict[str, object], ...]

    @property
    def ok(self) -> bool:
        return not any(bool(row["blocking"]) for row in self.mismatches)


@dataclass(frozen=True)
class CommonResult:
    """All summaries plus exact study and pair populations."""

    summaries: tuple[QuantitySummary, ...]
    study_effects: tuple[StudyEffect, ...]
    base_pairs: tuple[Pair, ...]
    base_population_id: str
    populations: tuple[Population, ...]


def complete_pairs(profiles: Sequence[Profile], traits: Sequence[str]) -> list[Pair]:
    """Build complete pairs after the shared duplicate/key/metadata checks."""
    pairs = [
        Pair(inc.study_id, inc.experiment_id, inc.condition, inc.body_site, inc, dec)
        for inc, dec in _matched_profiles(profiles, traits)
    ]
    for pair in pairs:
        if pair.inc.total < 0 or pair.dec.total < 0:
            raise ValueError(f"Negative retained total for {pair.pair_id}")
    return pairs


def _pair_digest(pairs: Sequence[Pair]) -> str:
    payload = "".join(f"{p.study_id}\t{p.experiment_id}\n" for p in sorted(pairs, key=_pair_key))
    return hashlib.sha256(payload.encode()).hexdigest()


def _pair_key(pair: Pair) -> tuple[str, str]:
    return pair.study_id, pair.experiment_id


def _population(
    population_type: str,
    numerator: str,
    reference: str,
    eligibility_rule: str,
    pairs: Sequence[Pair],
) -> Population:
    prefix = population_type
    if numerator:
        prefix += f":{numerator}"
    if reference:
        prefix += f":{reference}"
    population_id = f"{prefix}:sha256:{_pair_digest(pairs)}"
    return Population(
        population_id,
        population_type,
        numerator,
        reference,
        eligibility_rule,
        tuple(sorted(pairs, key=_pair_key)),
    )


def fraction_contrast(pair: Pair, trait: str) -> float:
    """Increased-minus-decreased fraction in raw fraction units."""
    if pair.inc.total <= 0 or pair.dec.total <= 0:
        raise ValueError(f"Fraction requested for an ineligible zero-total pair: {pair.pair_id}")
    return (
        pair.inc.guild_counts[trait] / pair.inc.total
        - pair.dec.guild_counts[trait] / pair.dec.total
    )


def log_ratio_contrast(pair: Pair, numerator: str, reference: str) -> float:
    """Paired log2(A/B) contrast; all four named counts must be positive."""
    counts = (
        pair.inc.guild_counts[numerator],
        pair.inc.guild_counts[reference],
        pair.dec.guild_counts[numerator],
        pair.dec.guild_counts[reference],
    )
    if min(counts) <= 0:
        raise ValueError(f"Log-ratio requested for an ineligible zero count: {pair.pair_id}")
    a_inc, b_inc, a_dec, b_dec = (math.log2(count) for count in counts)
    return math.fsum((a_inc, -b_inc, -a_dec, b_dec))


def _fraction_contrast_for(trait: str) -> Callable[[Pair], float]:
    def contrast(pair: Pair) -> float:
        return fraction_contrast(pair, trait)

    return contrast


def _log_ratio_contrast_for(numerator: str, reference: str) -> Callable[[Pair], float]:
    def contrast(pair: Pair) -> float:
        return log_ratio_contrast(pair, numerator, reference)

    return contrast


def _summarize_quantity(
    *,
    base_pairs: Sequence[Pair],
    population: Population,
    quantity_type: str,
    numerator: str,
    reference: str,
    source_columns: Mapping[str, str],
    contrast: Callable[[Pair], float],
    scale: float,
    effect_units: str,
    tie_tolerance: float,
    denominator_definition: str,
    population_definition: str,
    excluded_zero_total_pairs: int,
    excluded_zero_named_count_pairs: int,
) -> tuple[QuantitySummary, list[StudyEffect]]:
    by_study_base: defaultdict[str, int] = defaultdict(int)
    by_study_values: defaultdict[str, list[float]] = defaultdict(list)
    for pair in base_pairs:
        by_study_base[pair.study_id] += 1
    for pair in population.pairs:
        by_study_values[pair.study_id].append(contrast(pair))

    raw_study_effects: dict[str, float] = {
        study: fmean(values) for study, values in by_study_values.items()
    }
    positive = sum(value > tie_tolerance for value in raw_study_effects.values())
    negative = sum(value < -tie_tolerance for value in raw_study_effects.values())
    tied = len(raw_study_effects) - positive - negative
    non_tied = positive + negative
    share: float | None = None
    ci_low: float | None = None
    ci_high: float | None = None
    if non_tied:
        interval = binomtest(positive, non_tied).proportion_ci(
            confidence_level=CONFIDENCE_LEVEL, method="exact"
        )
        share = positive / non_tied
        ci_low, ci_high = float(interval.low), float(interval.high)

    values = list(raw_study_effects.values())
    status = "estimated" if non_tied else "all_tied" if values else "no_eligible_pairs"
    partial = sum(0 < len(by_study_values[study]) < total for study, total in by_study_base.items())
    if quantity_type == FRACTION_MATCHED:
        row_id = (
            f"{quantity_type}:{numerator}:population:{population.numerator}/{population.reference}"
        )
    else:
        row_id = (
            f"{quantity_type}:{numerator}:{reference}"
            if reference
            else f"{quantity_type}:{numerator}"
        )
    summary = QuantitySummary(
        row_id=row_id,
        comparison_spec=COMPARISON_SPEC,
        quantity_type=quantity_type,
        numerator=numerator,
        reference=reference,
        numerator_source_column=source_columns[numerator],
        reference_source_column=source_columns.get(reference, ""),
        population_id=population.population_id,
        population_numerator=population.numerator,
        population_reference=population.reference,
        population_definition=population_definition,
        denominator_definition=denominator_definition,
        eligibility_rule=population.eligibility_rule,
        independent_unit="study",
        study_weighting=(
            "mean eligible experiment contrasts within study; equal weight per eligible study"
        ),
        base_unique_studies=len(by_study_base),
        eligible_unique_studies=len(values),
        entirely_excluded_studies=len(by_study_base) - len(values),
        partially_observed_studies=partial,
        base_experiment_pairs=len(base_pairs),
        eligible_experiment_pairs=len(population.pairs),
        excluded_zero_total_pairs=excluded_zero_total_pairs,
        excluded_zero_named_count_pairs=excluded_zero_named_count_pairs,
        excluded_missing_or_malformed_pairs=0,
        positive_studies=positive,
        negative_studies=negative,
        tied_studies=tied,
        non_tied_studies=non_tied,
        mean_study_effect=fmean(values) * scale if values else None,
        median_study_effect=median(values) * scale if values else None,
        effect_units=effect_units,
        proportion_positive_non_tied=share,
        proportion_ci_low=ci_low,
        proportion_ci_high=ci_high,
        uncertainty_parameter=(
            "exact 95% binomial interval for the proportion of non-tied eligible studies "
            "with positive contrast"
        ),
        mean_effect_uncertainty="unavailable; repeated-study effect uncertainty is tracked in #133",
        status=status,
        calibration_status="descriptive estimate; not evaluated by full-pipeline calibration #182",
    )
    study_rows: list[StudyEffect] = []
    for study, total in sorted(by_study_base.items()):
        eligible = len(by_study_values[study])
        raw_effect = raw_study_effects.get(study)
        if raw_effect is None:
            direction = "ineligible"
        elif raw_effect > tie_tolerance:
            direction = "positive"
        elif raw_effect < -tie_tolerance:
            direction = "negative"
        else:
            direction = "tied"
        study_rows.append(
            StudyEffect(
                row_id=row_id,
                population_id=population.population_id,
                quantity_type=quantity_type,
                numerator=numerator,
                reference=reference,
                study_id=study,
                base_experiment_pairs=total,
                eligible_experiment_pairs=eligible,
                excluded_experiment_pairs=total - eligible,
                mean_study_effect=raw_effect * scale if raw_effect is not None else None,
                effect_units=effect_units,
                direction_status=direction,
            )
        )
    return summary, study_rows


def build_common_quantities(
    profiles: Sequence[Profile],
    traits: Sequence[str],
    source_columns: Mapping[str, str] | None = None,
) -> CommonResult:
    """Calculate 6 full fractions, 30 matched fractions, and 15 ratios."""
    if len(traits) < 2 or len(set(traits)) != len(traits):
        raise ValueError("Provide at least two distinct traits")
    sources = dict(source_columns or {trait: trait for trait in traits})
    missing_sources = set(traits).difference(sources)
    if missing_sources:
        raise ValueError(f"Missing source columns for {sorted(missing_sources)}")
    pairs = complete_pairs(profiles, traits)
    if not pairs:
        raise ValueError("No paired profiles available for common quantities")

    base_id = f"base:sha256:{_pair_digest(pairs)}"
    positive_totals = [pair for pair in pairs if pair.inc.total > 0 and pair.dec.total > 0]
    full_population = _population(
        "positive-retained-totals",
        "",
        "",
        "both retained non-MediaDive profile totals are positive; focal zero is valid",
        positive_totals,
    )
    populations = [full_population]
    summaries: list[QuantitySummary] = []
    study_effects: list[StudyEffect] = []
    zero_total_pairs = len(pairs) - len(positive_totals)

    for trait in traits:
        summary, studies = _summarize_quantity(
            base_pairs=pairs,
            population=full_population,
            quantity_type=FRACTION_FULL,
            numerator=trait,
            reference="",
            source_columns=sources,
            contrast=_fraction_contrast_for(trait),
            scale=100.0,
            effect_units="percentage-point increased-minus-decreased fraction contrast",
            tie_tolerance=FRACTION_TIE_TOLERANCE,
            denominator_definition=DENOMINATOR_DEFINITION,
            population_definition="full fraction-eligible base population",
            excluded_zero_total_pairs=zero_total_pairs,
            excluded_zero_named_count_pairs=0,
        )
        summaries.append(summary)
        study_effects.extend(studies)

    for numerator, reference in combinations(traits, 2):
        ratio_pairs = [pair for pair in pairs if pair.positive((numerator, reference))]
        population = _population(
            "ratio-complete-case",
            numerator,
            reference,
            "all four named counts are positive; no pseudocount; whole pair excluded otherwise",
            ratio_pairs,
        )
        populations.append(population)
        zero_named = len(pairs) - len(ratio_pairs)
        for trait in (numerator, reference):
            summary, studies = _summarize_quantity(
                base_pairs=pairs,
                population=population,
                quantity_type=FRACTION_MATCHED,
                numerator=trait,
                reference="",
                source_columns=sources,
                contrast=_fraction_contrast_for(trait),
                scale=100.0,
                effect_units="percentage-point increased-minus-decreased fraction contrast",
                tie_tolerance=FRACTION_TIE_TOLERANCE,
                denominator_definition=DENOMINATOR_DEFINITION,
                population_definition=(
                    f"fraction for {trait} on the exact {numerator}/{reference} ratio population"
                ),
                excluded_zero_total_pairs=0,
                excluded_zero_named_count_pairs=zero_named,
            )
            summaries.append(summary)
            study_effects.extend(studies)
        summary, studies = _summarize_quantity(
            base_pairs=pairs,
            population=population,
            quantity_type=LOG_RATIO,
            numerator=numerator,
            reference=reference,
            source_columns=sources,
            contrast=_log_ratio_contrast_for(numerator, reference),
            scale=1.0,
            effect_units="log2 increased-minus-decreased change in numerator/reference",
            tie_tolerance=LOG_RATIO_TIE_TOLERANCE,
            denominator_definition=(
                "no profile total in estimand; the shared total cancels, but the named "
                "reference remains"
            ),
            population_definition=f"complete-case population for {numerator}/{reference}",
            excluded_zero_total_pairs=0,
            excluded_zero_named_count_pairs=zero_named,
        )
        summaries.append(summary)
        study_effects.extend(studies)

    return CommonResult(
        tuple(summaries),
        tuple(study_effects),
        tuple(sorted(pairs, key=_pair_key)),
        base_id,
        tuple(populations),
    )


def inspect_trait_columns(
    profile_path: Path, trait_map: Mapping[str, Sequence[str]]
) -> tuple[list[str], list[TraitIdentity]]:
    """Return the exact retained universe and uniquely resolved Table 1 columns."""
    with profile_path.open(encoding="utf-8", newline="") as handle:
        header = next(csv.reader(handle, delimiter="\t"), None)
    if header is None:
        raise ValueError(f"{profile_path} is empty")
    if header[:3] != list(PROFILE_ID_COLUMNS):
        raise ValueError(f"{profile_path} must start with {PROFILE_ID_COLUMNS}")
    if len(header) != len(set(header)):
        raise ValueError("Duplicate profile-matrix column names are ambiguous")
    retained = [
        column
        for column in header
        if column not in PROFILE_ID_COLUMNS and MEDIADIVE_MARKER not in column.lower()
    ]
    identities: list[TraitIdentity] = []
    for order, (trait, ids) in enumerate(trait_map.items(), start=1):
        wanted = set(ids)
        matches = [
            column for column in retained if column in wanted or column_term_id(column) in wanted
        ]
        if len(matches) != 1:
            raise ValueError(f"Trait {trait!r} resolved to {len(matches)} columns: {matches}")
        source = matches[0]
        identities.append(
            TraitIdentity(
                order=order,
                trait=trait,
                display_label=TABLE1_TRAIT_LABELS.get(trait, trait),
                declared_term_ids=";".join(ids),
                source_column=source,
                source_term_id=column_term_id(source),
            )
        )
    return retained, identities


def _whole_number(raw: str, *, field: str, profile_id: str) -> int:
    try:
        value = Decimal(raw)
    except InvalidOperation as error:
        raise ValueError(f"Malformed {field} for {profile_id}: {raw!r}") from error
    if not value.is_finite() or value != value.to_integral_value():
        raise ValueError(f"Non-integral {field} for {profile_id}: {raw!r}")
    return int(value)


def _optional_whole_number(raw: str, *, field: str, profile_id: str) -> int | None:
    if raw.strip().upper() in {"", "NA"}:
        return None
    return _whole_number(raw, field=field, profile_id=profile_id)


def _required_reference_total(row: Mapping[str, int | None], field: str, profile_id: str) -> int:
    value = row[field]
    if value is None:
        raise ValueError(f"Required {field} is missing for {profile_id}")
    return value


def audit_denominators(
    profiles: Sequence[Profile], reference_path: Path, reference_revision: str
) -> DenominatorAudit:
    """Compare current matched keys/totals with the pinned PR #197 audit.

    The PR artifact contains every raw profile plus a ``paired`` flag.  Denominator
    parity is checked on the current shared base population against those raw rows;
    its paired flag is compared separately so a population mismatch cannot be
    mistaken for a denominator mismatch.
    """
    required = {
        "profile_id",
        "legacy_total",
        "harmonized_total",
        "mediadive_total",
        "ilyas_total",
        "noel_total",
        "paired",
    }
    reference: dict[str, dict[str, int | None]] = {}
    reference_paired: set[str] = set()
    with reference_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Reference denominator audit is missing {sorted(missing)}")
        for row in reader:
            paired = row["paired"].strip().casefold()
            if paired not in {"true", "false"}:
                raise ValueError(f"Unexpected paired value: {row['paired']!r}")
            profile_id = row["profile_id"]
            if profile_id in reference:
                raise ValueError(f"Duplicate reference profile ID: {profile_id}")
            if paired == "true":
                reference_paired.add(profile_id)
            reference[profile_id] = {
                "legacy_total": _whole_number(
                    row["legacy_total"], field="legacy_total", profile_id=profile_id
                ),
                "harmonized_total": _whole_number(
                    row["harmonized_total"], field="harmonized_total", profile_id=profile_id
                ),
                "mediadive_total": _whole_number(
                    row["mediadive_total"], field="mediadive_total", profile_id=profile_id
                ),
                "ilyas_total": _whole_number(
                    row["ilyas_total"], field="ilyas_total", profile_id=profile_id
                ),
                "noel_total": _optional_whole_number(
                    row["noel_total"], field="noel_total", profile_id=profile_id
                ),
            }

    current: dict[str, int] = {}
    for profile in profiles:
        profile_id = f"{profile.study_id}|{profile.experiment_id}|{profile.inc_dec}"
        if profile_id in current:
            raise ValueError(f"Duplicate current profile ID: {profile_id}")
        current[profile_id] = profile.total

    mismatches: list[dict[str, object]] = []
    missing_from_artifact = sorted(set(current).difference(reference))
    for profile_id in missing_from_artifact:
        mismatches.append(
            {
                "profile_id": profile_id,
                "check": "current_base_key_missing_from_pr197_artifact",
                "current": current[profile_id],
                "reference": "missing",
                "blocking": True,
            }
        )
    current_not_pr_paired = sorted(set(current).difference(reference_paired))
    pr_paired_not_current = sorted(reference_paired.difference(current))
    for profile_id in current_not_pr_paired:
        mismatches.append(
            {
                "profile_id": profile_id,
                "check": "current_base_profile_not_pr197_paired",
                "current": "paired",
                "reference": "not paired",
                "blocking": False,
            }
        )
    for profile_id in pr_paired_not_current:
        mismatches.append(
            {
                "profile_id": profile_id,
                "check": "pr197_paired_profile_not_current_base",
                "current": "not paired",
                "reference": reference[profile_id]["harmonized_total"],
                "blocking": False,
            }
        )

    checks = {
        "current_loader_vs_pr197_harmonized": 0,
        "pr197_ilyas_vs_harmonized": 0,
        "pr197_noel_vs_harmonized": 0,
        "pr197_legacy_minus_harmonized_vs_mediadive": 0,
    }
    legacy_changed = 0
    common = sorted(set(current).intersection(reference))
    for profile_id in common:
        reference_row = reference[profile_id]
        legacy_total = _required_reference_total(reference_row, "legacy_total", profile_id)
        harmonized_total = _required_reference_total(reference_row, "harmonized_total", profile_id)
        mediadive_total = _required_reference_total(reference_row, "mediadive_total", profile_id)
        ilyas_total = _required_reference_total(reference_row, "ilyas_total", profile_id)
        noel_total = reference_row["noel_total"]
        comparisons = {
            "current_loader_vs_pr197_harmonized": (
                current[profile_id],
                harmonized_total,
            ),
            "pr197_ilyas_vs_harmonized": (ilyas_total, harmonized_total),
            "pr197_noel_vs_harmonized": (noel_total, harmonized_total),
            "pr197_legacy_minus_harmonized_vs_mediadive": (
                legacy_total - harmonized_total,
                mediadive_total,
            ),
        }
        legacy_changed += legacy_total != harmonized_total
        for check, (left, right) in comparisons.items():
            if left != right:
                checks[check] += 1
                mismatches.append(
                    {
                        "profile_id": profile_id,
                        "check": check,
                        "current": left,
                        "reference": right,
                        "blocking": True,
                    }
                )

    paired_noel_missing = sorted(
        profile_id for profile_id in reference_paired if reference[profile_id]["noel_total"] is None
    )
    for profile_id in paired_noel_missing:
        mismatches.append(
            {
                "profile_id": profile_id,
                "check": "pr197_paired_noel_total_missing",
                "current": current.get(profile_id, "not in current base"),
                "reference": "NA",
                "blocking": False,
            }
        )

    summary: list[dict[str, object]] = [
        {
            "check": "current_base_profile_keys_in_pr197_artifact",
            "profiles_compared": len(common),
            "mismatches": len(missing_from_artifact),
            "status": "pass" if not missing_from_artifact else "fail",
            "reference_revision": reference_revision,
        },
        {
            "check": "pr197_paired_membership_vs_current_shared_base",
            "profiles_compared": len(current),
            "mismatches": len(current_not_pr_paired) + len(pr_paired_not_current),
            "status": (
                "pass"
                if not current_not_pr_paired and not pr_paired_not_current
                else "population_discrepancy"
            ),
            "reference_revision": reference_revision,
        },
        {
            "check": "pr197_paired_rows_missing_noel_total",
            "profiles_compared": len(reference_paired),
            "mismatches": len(paired_noel_missing),
            "status": ("pass" if not paired_noel_missing else "reference_field_incomplete"),
            "reference_revision": reference_revision,
        },
    ]
    summary.extend(
        {
            "check": check,
            "profiles_compared": len(common),
            "mismatches": count,
            "status": ("pass" if count == 0 else "fail"),
            "reference_revision": reference_revision,
        }
        for check, count in checks.items()
    )
    summary.append(
        {
            "check": "profiles_changed_by_legacy_MediaDive_inclusion",
            "profiles_compared": len(common),
            "mismatches": legacy_changed,
            "status": "diagnostic_expected_difference",
            "reference_revision": reference_revision,
        }
    )
    return DenominatorAudit(reference_revision, tuple(summary), tuple(mismatches))


def _sha256(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def _format(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return format(value, ".17g")
    return str(value)


def _write_dicts(path: Path, rows: Sequence[Mapping[str, object]], names: Sequence[str]) -> Path:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(names), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({name: _format(row.get(name)) for name in names})
    return path


def _display(value: float | None) -> str:
    return "NA" if value is None else f"{value:.4g}"


def _write_report(result: CommonResult, audit: DenominatorAudit, path: Path) -> Path:
    full = [row for row in result.summaries if row.quantity_type == FRACTION_FULL]
    ratios = [row for row in result.summaries if row.quantity_type == LOG_RATIO]
    matched = [row for row in result.summaries if row.quantity_type == FRACTION_MATCHED]
    lines = [
        "# Common study-weighted fraction and ratio quantities (generated)",
        "",
        f"Version `{VERSION}` under `{COMPARISON_SPEC}`. These are descriptive estimates, not a",
        "new test or correction family. Every study contributes one mean of its eligible",
        "experiment contrasts. The exact 95% interval is for the positive-direction",
        "proportion among non-tied studies, not for the mean effect.",
        "",
        "## Full-population fraction contrasts",
        "",
        "| Trait | Studies | Pairs | Mean pp | Median pp | Positive share (95% CI) | Ties |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for summary in full:
        lines.append(
            f"| {summary.numerator} | "
            f"{summary.eligible_unique_studies}/{summary.base_unique_studies} | "
            f"{summary.eligible_experiment_pairs}/{summary.base_experiment_pairs} | "
            f"{_display(summary.mean_study_effect)} | "
            f"{_display(summary.median_study_effect)} | "
            f"{_display(summary.proportion_positive_non_tied)} "
            f"({_display(summary.proportion_ci_low)}, "
            f"{_display(summary.proportion_ci_high)}) | {summary.tied_studies} |"
        )
    lines += [
        "",
        "## Ratios and the two fractions on each exact ratio population",
        "",
        "| Ratio | Studies | Pairs | Numerator fraction mean pp | Reference fraction mean pp | "
        "Ratio mean log2 | Ratio median log2 | Positive share (95% CI) | Ties |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for ratio in ratios:
        same = [row for row in matched if row.population_id == ratio.population_id]
        numerator_fraction = next(row for row in same if row.numerator == ratio.numerator)
        reference_fraction = next(row for row in same if row.numerator == ratio.reference)
        lines.append(
            f"| {ratio.numerator}/{ratio.reference} | "
            f"{ratio.eligible_unique_studies}/{ratio.base_unique_studies} | "
            f"{ratio.eligible_experiment_pairs}/{ratio.base_experiment_pairs} | "
            f"{_display(numerator_fraction.mean_study_effect)} | "
            f"{_display(reference_fraction.mean_study_effect)} | "
            f"{_display(ratio.mean_study_effect)} | {_display(ratio.median_study_effect)} | "
            f"{_display(ratio.proportion_positive_non_tied)} "
            f"({_display(ratio.proportion_ci_low)}, {_display(ratio.proportion_ci_high)}) | "
            f"{ratio.tied_studies} |"
        )
    lines += [
        "",
        "## Denominator audit",
        "",
        "| Check | Profiles | Mismatches | Status |",
        "| --- | --- | --- | --- |",
    ]
    for audit_row in audit.summary:
        lines.append(
            f"| {audit_row['check']} | {audit_row['profiles_compared']} | "
            f"{audit_row['mismatches']} | {audit_row['status']} |"
        )
    lines += [
        "",
        "The full TSV retains population IDs, partial exclusions, ties, validity state, exact",
        "source columns, uncertainty parameter, and calibration status. Pair and study-effect",
        "manifests make equality of the matched fraction/ratio populations auditable.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_outputs(
    result: CommonResult,
    output_dir: Path,
    *,
    retained_columns: Sequence[str],
    identities: Sequence[TraitIdentity],
    audit: DenominatorAudit,
    provenance: Mapping[str, object],
) -> list[Path]:
    """Write common rows, exact manifests, parity audit, provenance, and a view."""
    output_dir.mkdir(parents=True, exist_ok=True)
    written = [
        _write_dicts(
            output_dir / "common_quantities.tsv",
            [asdict(summary) for summary in result.summaries],
            [field.name for field in fields(QuantitySummary)],
        ),
        _write_dicts(
            output_dir / "common_quantity_study_effects.tsv",
            [asdict(study) for study in result.study_effects],
            [field.name for field in fields(StudyEffect)],
        ),
        _write_dicts(
            output_dir / "trait_manifest.tsv",
            [asdict(identity) for identity in identities],
            [field.name for field in fields(TraitIdentity)],
        ),
    ]
    selected = {identity.source_column: identity.trait for identity in identities}
    retained_rows = [
        {
            "order": order,
            "source_column": column,
            "source_term_id": column_term_id(column),
            "table1_trait": selected.get(column, ""),
        }
        for order, column in enumerate(retained_columns, start=1)
    ]
    written.append(
        _write_dicts(
            output_dir / "retained_columns.tsv",
            retained_rows,
            ("order", "source_column", "source_term_id", "table1_trait"),
        )
    )
    base_rows = [
        {
            "base_population_id": result.base_population_id,
            "pair_id": pair.pair_id,
            "study_id": pair.study_id,
            "experiment_id": pair.experiment_id,
            "condition": pair.condition,
            "body_site": pair.body_site,
            "increased_profile_id": f"{pair.pair_id}|inc",
            "decreased_profile_id": f"{pair.pair_id}|dec",
            "increased_retained_total": pair.inc.total,
            "decreased_retained_total": pair.dec.total,
        }
        for pair in result.base_pairs
    ]
    written.append(
        _write_dicts(
            output_dir / "base_experiment_pairs.tsv",
            base_rows,
            (
                "base_population_id",
                "pair_id",
                "study_id",
                "experiment_id",
                "condition",
                "body_site",
                "increased_profile_id",
                "decreased_profile_id",
                "increased_retained_total",
                "decreased_retained_total",
            ),
        )
    )
    population_rows = [
        {
            "population_id": population.population_id,
            "population_type": population.population_type,
            "numerator": population.numerator,
            "reference": population.reference,
            "eligibility_rule": population.eligibility_rule,
            "pair_id": pair.pair_id,
            "study_id": pair.study_id,
            "experiment_id": pair.experiment_id,
        }
        for population in result.populations
        for pair in population.pairs
    ]
    written.append(
        _write_dicts(
            output_dir / "eligible_experiment_pairs.tsv",
            population_rows,
            (
                "population_id",
                "population_type",
                "numerator",
                "reference",
                "eligibility_rule",
                "pair_id",
                "study_id",
                "experiment_id",
            ),
        )
    )
    written.append(
        _write_dicts(
            output_dir / "denominator_audit.tsv",
            audit.summary,
            ("check", "profiles_compared", "mismatches", "status", "reference_revision"),
        )
    )
    written.append(
        _write_dicts(
            output_dir / "denominator_mismatches.tsv",
            audit.mismatches,
            ("profile_id", "check", "current", "reference", "blocking"),
        )
    )
    written.append(_write_report(result, audit, output_dir / "common_quantity_report.md"))

    fingerprints = {path.name: _sha256(path) for path in written}
    provenance_path = output_dir / "common_quantity_provenance.json"
    provenance_path.write_text(
        json.dumps(
            {
                **provenance,
                "version": VERSION,
                "comparison_spec": COMPARISON_SPEC,
                "denominator_definition": DENOMINATOR_DEFINITION,
                "retained_column_count": len(retained_columns),
                "base_population_id": result.base_population_id,
                "base_experiment_pairs": len(result.base_pairs),
                "base_unique_studies": len({pair.study_id for pair in result.base_pairs}),
                "quantity_rows": len(result.summaries),
                "population_manifests": {
                    population.population_id: len(population.pairs)
                    for population in result.populations
                },
                "fraction_zero_policy": (
                    "exclude whole pair if either retained total is zero; a focal zero with "
                    "positive totals is valid"
                ),
                "ratio_zero_policy": (
                    "require all four named counts positive; exclude whole pair; no pseudocount"
                ),
                "count_validity_policy": (
                    "blank, malformed, and negative retained counts fail before estimation; "
                    "they are not observed zeros or eligibility exclusions"
                ),
                "study_weighting": (
                    "mean eligible experiment contrasts within study; equal weight per study"
                ),
                "uncertainty": (
                    "exact 95% binomial interval for the non-tied positive-study proportion; "
                    "no interval for the mean effect"
                ),
                "testing": (
                    "none; common quantities are descriptive and create no p-value, q-value, "
                    "correction family, calibrated decision, or independent replication"
                ),
                "calibration_status": "not evaluated; route through #182",
                "sensitivity_boundary": (
                    "only required full and ratio-matched populations are emitted; broader "
                    "reference/zero/common-support sensitivity belongs to #183"
                ),
                "integration_target": "#178",
                "output_sha256": fingerprints,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    written.append(provenance_path)
    return written


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile-input", type=Path, default=DEFAULT_PROFILE_PATH)
    parser.add_argument("--experiments-input", type=Path, default=DEFAULT_EXPERIMENTS_PATH)
    parser.add_argument("--full-dump-input", type=Path, default=DEFAULT_FULL_DUMP_PATH)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--composition-audit-input", type=Path, required=True)
    parser.add_argument("--composition-reference-revision", required=True)
    parser.add_argument("--profile-git-blob", default="")
    parser.add_argument("--expected-profile-sha256")
    parser.add_argument("--code-commit", help="Source revision, supplied by the Justfile recipe")
    args = parser.parse_args(argv)

    profile_sha256 = _sha256(args.profile_input)
    if args.expected_profile_sha256 and profile_sha256 != args.expected_profile_sha256:
        raise ValueError(
            f"Profile matrix SHA256 changed: {profile_sha256}; "
            f"expected {args.expected_profile_sha256}"
        )
    retained, identities = inspect_trait_columns(args.profile_input, TABLE1_TRAITS)
    source_columns = {identity.trait: identity.source_column for identity in identities}
    profiles, stats = load_profiles(
        args.profile_input,
        args.experiments_input,
        TABLE1_TRAITS,
        reject_missing_counts=True,
    )
    vintage = read_bugsigdb_vintage(args.full_dump_input)
    write_stats_report(stats, vintage)
    result = build_common_quantities(profiles, tuple(TABLE1_TRAITS), source_columns)
    audit = audit_denominators(
        profiles, args.composition_audit_input, args.composition_reference_revision
    )
    if not audit.ok:
        blocking = [row for row in audit.mismatches if bool(row["blocking"])]
        examples = "; ".join(f"{row['profile_id']}:{row['check']}" for row in blocking[:3])
        raise ValueError(f"Denominator parity audit failed ({len(blocking)}): {examples}")

    provenance: dict[str, object] = {
        "code_commit": args.code_commit,
        "profile_path": str(args.profile_input),
        "profile_sha256": profile_sha256,
        "profile_git_blob": args.profile_git_blob,
        "experiments_path": str(args.experiments_input),
        "experiments_sha256": _sha256(args.experiments_input),
        "full_dump_path": str(args.full_dump_input),
        "full_dump_sha256": _sha256(args.full_dump_input),
        "bugsigdb_vintage": vintage,
        "analysis_source_sha256": _sha256(Path(__file__)),
        "shared_loader_source_sha256": _sha256(Path(__file__).with_name("trait_guilds.py")),
        "composition_reference_revision": args.composition_reference_revision,
        "composition_audit_path": str(args.composition_audit_input),
        "composition_audit_sha256": _sha256(args.composition_audit_input),
        "composition_audit_status": (
            "denominator parity is required on every current shared-base profile; PR #197 "
            "paired-membership and missing-field discrepancies remain separately recorded"
        ),
        "input_filter_counts": asdict(stats),
        "trait_identities": [asdict(identity) for identity in identities],
    }
    paths = write_outputs(
        result,
        args.output_dir,
        retained_columns=retained,
        identities=identities,
        audit=audit,
        provenance=provenance,
    )
    for path in paths:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
