"""How the paired log-ratio results depend on who is eligible (#183).

The log-ratio analysis (:mod:`bsdbng.trait_composition`) excludes an experiment
from A/B whenever either trait is zero on either side, so each ratio is measured on
its own complete-case population. This module runs the scenarios declared in
``methods/table1_log_ratio_sensitivity_protocol.json`` before any result was seen.
It separates three things a single comparison of fraction and ratio results would
mix up: a change of population (the same fraction question on fewer pairs), a
change of estimand (fractions against ratios on the same pairs), and dependence on
particular strata or studies.

Every scenario keeps its own correction family. Nothing here selects a preferred
scenario, and stability across these scenarios is robustness to these choices on
shared annotation data, not independent replication or absolute abundance.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from collections.abc import Callable, Iterable, Sequence
from dataclasses import asdict, dataclass, fields, replace
from itertools import combinations
from pathlib import Path
from statistics import fmean, median

from scipy.stats import binomtest, false_discovery_control

from bsdbng.paths import DERIVED_DATA_DIR, REPO_ROOT
from bsdbng.trait_composition import _matched_profiles
from bsdbng.trait_guilds import (
    DEFAULT_EXPERIMENTS_PATH,
    DEFAULT_FULL_DUMP_PATH,
    DEFAULT_PROFILE_PATH,
    TABLE1_TRAITS,
    Profile,
    load_profiles,
    read_bugsigdb_vintage,
    write_stats_report,
)

PROTOCOL_PATH = REPO_ROOT / "methods" / "table1_log_ratio_sensitivity_protocol.json"
DEFAULT_TAXA_PATH = DERIVED_DATA_DIR / "bugsigdb_taxa.tsv"
DEFAULT_OUTPUT_DIR = DERIVED_DATA_DIR / "table1_composition_sensitivity"
VERSION = "log-ratio-sensitivity-v1"

TIE_TOLERANCE = 1e-12
CONFIDENCE_LEVEL = 0.95
MIN_STUDIES = 30
REPORT_Q = 0.05

TESTED = "tested"
ALL_TIED = "all_tied"
INSUFFICIENT = "insufficient_support"
NO_PAIRS = "no_eligible_pairs"

S0 = "S0_historical_fraction"
S1 = "S1_fraction_baseline"
S2 = "S2_ratio_complete_case"
S3 = "S3_ratio_matched_fractions"
S4 = "S4_six_trait_common_support"
S5 = "S5_triple_common_support"
S6 = "S6_feces_only"
S7 = "S7_excluding_feces"

FECES = "Feces"
Contrast = Callable[["Pair"], float | None]


@dataclass(frozen=True)
class Pair:
    """One complete experiment pair with its metadata."""

    study_id: str
    experiment_id: str
    condition: str
    body_site: str
    inc: Profile
    dec: Profile

    def positive(self, traits: Iterable[str]) -> bool:
        return all(self.inc.guild_counts[t] > 0 and self.dec.guild_counts[t] > 0 for t in traits)


def complete_pairs(profiles: Sequence[Profile], traits: Sequence[str]) -> list[Pair]:
    """Validated complete pairs, using the log-ratio module's own key checks."""
    return [
        Pair(inc.study_id, inc.experiment_id, inc.condition, inc.body_site, inc, dec)
        for inc, dec in _matched_profiles(profiles, traits)
    ]


def fraction_contrast(trait: str, *, zero_total_as_zero: bool) -> Contrast:
    """Increased minus decreased share of ``trait``.

    With ``zero_total_as_zero`` a profile with no retained annotation scores 0, as
    the historical fraction analysis does; otherwise its pair is ineligible.
    """

    def contrast(pair: Pair) -> float | None:
        if not zero_total_as_zero and (pair.inc.total <= 0 or pair.dec.total <= 0):
            return None
        return pair.inc.fraction(trait) - pair.dec.fraction(trait)

    return contrast


def restricted(contrast: Contrast, eligible: Callable[[Pair], bool]) -> Contrast:
    """The same contrast, defined only on pairs another rule admits."""

    def inner(pair: Pair) -> float | None:
        return contrast(pair) if eligible(pair) else None

    return inner


def log_ratio_contrast(numerator: str, reference: str) -> Contrast:
    """Change in log2(numerator/reference); undefined unless all four counts are positive."""

    def contrast(pair: Pair) -> float | None:
        if not pair.positive((numerator, reference)):
            return None
        return math.fsum(
            (
                math.log2(pair.inc.guild_counts[numerator]),
                -math.log2(pair.inc.guild_counts[reference]),
                -math.log2(pair.dec.guild_counts[numerator]),
                math.log2(pair.dec.guild_counts[reference]),
            )
        )

    return contrast


@dataclass(frozen=True)
class Row:
    """One scenario result. Undefined quantities are None, never zero."""

    scenario: str
    comparison: str
    numerator: str
    reference: str
    population: str
    effect_units: str
    studies_total: int
    studies_eligible: int
    studies_partially_observed: int
    pairs_total: int
    pairs_eligible: int
    studies_up: int
    studies_down: int
    studies_tied: int
    mean_effect: float | None
    median_effect: float | None
    share_up: float | None
    ci_low: float | None
    ci_high: float | None
    p_value: float | None
    status: str
    family: str = ""
    family_size: int = 0
    correction: str = ""
    q_value: float | None = None


def study_means(pairs: Sequence[Pair], contrast: Contrast) -> dict[str, float]:
    """One mean of eligible experiment contrasts per study."""
    by_study: defaultdict[str, list[float]] = defaultdict(list)
    for pair in pairs:
        value = contrast(pair)
        if value is not None:
            by_study[pair.study_id].append(value)
    return {study: fmean(values) for study, values in by_study.items()}


def sign_test(values: Iterable[float]) -> tuple[int, int, int, float | None]:
    """Up, down and tied counts, and the two-sided exact binomial p of the non-tied."""
    listed = list(values)
    up = sum(value > TIE_TOLERANCE for value in listed)
    down = sum(value < -TIE_TOLERANCE for value in listed)
    scored = up + down
    p_value = float(binomtest(up, scored, p=0.5).pvalue) if scored else None
    return up, down, len(listed) - scored, p_value


def summarize(
    *,
    scenario: str,
    comparison: str,
    numerator: str,
    reference: str,
    population: str,
    pairs: Sequence[Pair],
    contrast: Contrast,
    min_studies: int = MIN_STUDIES,
) -> Row:
    """Average eligible experiment contrasts within study, then test study signs."""
    per_study_pairs = Counter(pair.study_id for pair in pairs)
    eligible = Counter(pair.study_id for pair in pairs if contrast(pair) is not None)
    means = study_means(pairs, contrast)
    up, down, tied, p_value = sign_test(means.values())
    status = TESTED
    share = low = high = None
    if not means:
        status, p_value = NO_PAIRS, None
    elif len(means) < min_studies:
        status, p_value = INSUFFICIENT, None
    elif p_value is None:
        status = ALL_TIED
    else:
        interval = binomtest(up, up + down, p=0.5).proportion_ci(
            confidence_level=CONFIDENCE_LEVEL, method="exact"
        )
        share, low, high = up / (up + down), float(interval.low), float(interval.high)
    values = list(means.values())
    return Row(
        scenario=scenario,
        comparison=comparison,
        numerator=numerator,
        reference=reference,
        population=population,
        effect_units=(
            "fraction contrast, increased minus decreased"
            if comparison == "fraction"
            else "change in log2(numerator/reference)"
        ),
        studies_total=len(per_study_pairs),
        studies_eligible=len(means),
        studies_partially_observed=sum(
            1 for study, n in eligible.items() if n < per_study_pairs[study]
        ),
        pairs_total=len(pairs),
        pairs_eligible=sum(eligible.values()),
        studies_up=up,
        studies_down=down,
        studies_tied=tied,
        mean_effect=fmean(values) if values else None,
        median_effect=median(values) if values else None,
        share_up=share,
        ci_low=low,
        ci_high=high,
        p_value=p_value,
        status=status,
    )


def matched_reference(numerator: str, reference: str) -> str:
    """The label an S3 row carries, naming the ratio whose pairs define its population."""
    return f"on {numerator}/{reference} pairs"


def correct(rows: Sequence[Row], family: str, method: str) -> list[Row]:
    """Adjust one declared family. Untested rows hold their place at p = 1."""
    if not rows:
        return []
    adjusted = false_discovery_control(
        [1.0 if row.p_value is None else row.p_value for row in rows], method=method
    )
    return [
        replace(
            row,
            family=family,
            family_size=len(rows),
            correction="Benjamini-Hochberg" if method == "bh" else "Benjamini-Yekutieli",
            q_value=None if row.p_value is None else float(q),
        )
        for row, q in zip(rows, adjusted, strict=True)
    ]


def _fractions(
    scenario: str,
    pairs: Sequence[Pair],
    traits: Sequence[str],
    population: str,
    *,
    min_studies: int,
    zero_total_as_zero: bool = False,
) -> list[Row]:
    return [
        summarize(
            scenario=scenario,
            comparison="fraction",
            numerator=trait,
            reference="",
            population=population,
            pairs=pairs,
            contrast=fraction_contrast(trait, zero_total_as_zero=zero_total_as_zero),
            min_studies=min_studies,
        )
        for trait in traits
    ]


def _ratios(
    scenario: str,
    pairs: Sequence[Pair],
    ratio_pairs: Sequence[tuple[str, str]],
    population: str,
    *,
    min_studies: int,
) -> list[Row]:
    return [
        summarize(
            scenario=scenario,
            comparison="log_ratio",
            numerator=numerator,
            reference=reference,
            population=population,
            pairs=pairs,
            contrast=log_ratio_contrast(numerator, reference),
            min_studies=min_studies,
        )
        for numerator, reference in ratio_pairs
    ]


def run_scenarios(
    pairs: Sequence[Pair], traits: Sequence[str], *, min_studies: int = MIN_STUDIES
) -> list[Row]:
    """Every declared scenario except leave-one-study-out, each in its own family."""
    ratio_pairs = list(combinations(traits, 2))
    rows: list[Row] = []
    rows += correct(
        _fractions(
            S0,
            pairs,
            traits,
            "all complete pairs; zero-total profiles score fraction 0",
            zero_total_as_zero=True,
            min_studies=min_studies,
        ),
        "S0:six_fractions",
        "bh",
    )
    rows += correct(
        _fractions(
            S1,
            pairs,
            traits,
            "complete pairs with positive totals both sides",
            min_studies=min_studies,
        ),
        "S1:six_fractions",
        "bh",
    )
    rows += correct(
        _ratios(
            S2,
            pairs,
            ratio_pairs,
            "pairs with both traits positive",
            min_studies=min_studies,
        ),
        "S2:fifteen_ratios",
        "by",
    )
    matched: list[Row] = []
    for numerator, reference in ratio_pairs:
        admitted = log_ratio_contrast(numerator, reference)

        def in_ratio(pair: Pair, admitted: Contrast = admitted) -> bool:
            return admitted(pair) is not None

        for trait in (numerator, reference):
            matched.append(
                summarize(
                    scenario=S3,
                    comparison="fraction",
                    numerator=trait,
                    reference=matched_reference(numerator, reference),
                    population=f"pairs eligible for {numerator}/{reference}",
                    pairs=pairs,
                    contrast=restricted(
                        fraction_contrast(trait, zero_total_as_zero=False), in_ratio
                    ),
                    min_studies=min_studies,
                )
            )
    rows += correct(matched, "S3:thirty_matched_fractions", "by")

    common = [pair for pair in pairs if pair.positive(traits)]
    population = "pairs with all six traits positive on both sides"
    rows += correct(
        _fractions(S4, common, traits, population, min_studies=min_studies),
        "S4:six_fractions",
        "bh",
    )
    rows += correct(
        _ratios(
            S4,
            common,
            ratio_pairs,
            population,
            min_studies=min_studies,
        ),
        "S4:fifteen_ratios",
        "by",
    )

    triples: list[Row] = []
    for triple in combinations(traits, 3):
        shared = [pair for pair in pairs if pair.positive(triple)]
        triples += _ratios(
            S5,
            shared,
            list(combinations(triple, 2)),
            f"pairs with {', '.join(triple)} positive on both sides",
            min_studies=min_studies,
        )
    rows += correct(triples, "S5:sixty_triple_ratios", "by")

    feces = [pair for pair in pairs if pair.body_site == FECES]
    others = [pair for pair in pairs if pair.body_site != FECES]
    rows += correct(
        _ratios(
            S6,
            feces,
            ratio_pairs,
            "S2 pairs in exact body site Feces",
            min_studies=min_studies,
        ),
        "S6:fifteen_ratios",
        "by",
    )
    rows += correct(
        _ratios(
            S7,
            others,
            ratio_pairs,
            "S2 pairs outside exact body site Feces",
            min_studies=min_studies,
        ),
        "S7:fifteen_ratios",
        "by",
    )
    return rows


@dataclass(frozen=True)
class PopulationComparison:
    """One trait's fraction result at S1 against the same result on one ratio's pairs.

    S1 and S3 differ in three ways at once: the eligible studies, the family size and
    the correction method. This row keeps all three visible, compares the share and
    its interval rather than the q status, and splits a changed estimate from a
    detection lost to the larger family.
    """

    trait: str
    ratio_numerator: str
    ratio_reference: str
    baseline_scenario: str
    baseline_family: str
    baseline_family_size: int
    baseline_correction: str
    baseline_studies_eligible: int
    baseline_share_up: float | None
    baseline_ci_low: float | None
    baseline_ci_high: float | None
    baseline_p_value: float | None
    baseline_q_value: float | None
    restricted_scenario: str
    restricted_family: str
    restricted_family_size: int
    restricted_correction: str
    restricted_studies_eligible: int
    restricted_share_up: float | None
    restricted_ci_low: float | None
    restricted_ci_high: float | None
    restricted_p_value: float | None
    restricted_q_value: float | None
    studies_retained: float | None
    share_shift: float | None
    restricted_share_in_baseline_ci: bool | None
    intervals_overlap: bool | None
    estimate_agreement: str
    detection_change: str


def _intervals_overlap(baseline: Row, restricted: Row) -> bool | None:
    low_b, high_b = baseline.ci_low, baseline.ci_high
    low_r, high_r = restricted.ci_low, restricted.ci_high
    if low_b is None or high_b is None or low_r is None or high_r is None:
        return None
    return low_b <= high_r and low_r <= high_b


def _estimate_agreement(baseline: Row, restricted: Row, overlap: bool | None) -> str:
    """Whether the restricted share is compatible with the baseline share.

    Interval overlap is reported first, so a row whose interval has genuinely moved is
    never recorded only as a change of side.
    """
    if baseline.share_up is None or restricted.share_up is None or overlap is None:
        return "not comparable"
    reversed_side = (baseline.share_up - 0.5) * (restricted.share_up - 0.5) < 0
    if overlap:
        return "side reversed, intervals overlap" if reversed_side else "consistent"
    return "shifted and side reversed" if reversed_side else "shifted"


def _detection_change(baseline: Row, restricted: Row) -> str:
    """Split a lost detection into one lost at the correction and one lost before it."""
    detected = baseline.q_value is not None and baseline.q_value <= REPORT_Q
    restricted_detected = restricted.q_value is not None and restricted.q_value <= REPORT_Q
    restricted_raw = restricted.p_value is not None and restricted.p_value <= REPORT_Q
    if detected and restricted_detected:
        return "retained"
    if detected and restricted_raw:
        return "lost at correction"
    if detected:
        return "lost before correction"
    return "gained" if restricted_detected else "absent in both"


def population_comparison(rows: Sequence[Row], traits: Sequence[str]) -> list[PopulationComparison]:
    """Pair every S3 matched-fraction row with its S1 full-population counterpart."""
    baseline = {row.numerator: row for row in rows if row.scenario == S1}
    restricted = {(row.numerator, row.reference): row for row in rows if row.scenario == S3}
    comparisons: list[PopulationComparison] = []
    for numerator, reference in combinations(traits, 2):
        label = matched_reference(numerator, reference)
        for trait in (numerator, reference):
            base, rest = baseline[trait], restricted[(trait, label)]
            overlap = _intervals_overlap(base, rest)
            in_ci = (
                None
                if rest.share_up is None or base.ci_low is None or base.ci_high is None
                else base.ci_low <= rest.share_up <= base.ci_high
            )
            comparisons.append(
                PopulationComparison(
                    trait=trait,
                    ratio_numerator=numerator,
                    ratio_reference=reference,
                    baseline_scenario=base.scenario,
                    baseline_family=base.family,
                    baseline_family_size=base.family_size,
                    baseline_correction=base.correction,
                    baseline_studies_eligible=base.studies_eligible,
                    baseline_share_up=base.share_up,
                    baseline_ci_low=base.ci_low,
                    baseline_ci_high=base.ci_high,
                    baseline_p_value=base.p_value,
                    baseline_q_value=base.q_value,
                    restricted_scenario=rest.scenario,
                    restricted_family=rest.family,
                    restricted_family_size=rest.family_size,
                    restricted_correction=rest.correction,
                    restricted_studies_eligible=rest.studies_eligible,
                    restricted_share_up=rest.share_up,
                    restricted_ci_low=rest.ci_low,
                    restricted_ci_high=rest.ci_high,
                    restricted_p_value=rest.p_value,
                    restricted_q_value=rest.q_value,
                    studies_retained=(
                        rest.studies_eligible / base.studies_eligible
                        if base.studies_eligible
                        else None
                    ),
                    share_shift=(
                        None
                        if base.share_up is None or rest.share_up is None
                        else rest.share_up - base.share_up
                    ),
                    restricted_share_in_baseline_ci=in_ci,
                    intervals_overlap=overlap,
                    estimate_agreement=_estimate_agreement(base, rest, overlap),
                    detection_change=_detection_change(base, rest),
                )
            )
    return comparisons


@dataclass(frozen=True)
class Influence:
    """Leave-one-study-out summary for one S2 ratio (S8)."""

    numerator: str
    reference: str
    studies: int
    full_share_up: float | None
    full_q: float | None
    min_share_up: float | None
    max_share_up: float | None
    min_p: float | None
    max_p: float | None
    removals_changing_q_status: int
    most_influential_studies: str


def leave_one_study_out(pairs: Sequence[Pair], traits: Sequence[str]) -> list[Influence]:
    """Remove each eligible study in turn and re-adjust the 15-ratio BY family.

    The other 14 ratios keep their full-data p-values, so a change in status is
    attributable to the removed study's effect on this ratio alone.
    """
    ratio_pairs = list(combinations(traits, 2))
    means = {pair: study_means(pairs, log_ratio_contrast(*pair)) for pair in ratio_pairs}
    full_p = [sign_test(values.values())[3] for values in means.values()]
    base = [1.0 if p is None else p for p in full_p]
    full_q = [float(q) for q in false_discovery_control(base, method="by")]

    results: list[Influence] = []
    for index, ratio in enumerate(ratio_pairs):
        values = means[ratio]
        up, down, _tied, _p = sign_test(values.values())
        full_share = up / (up + down) if up + down else None
        shares: list[float] = []
        p_values: list[float] = []
        flips = 0
        shifts: list[tuple[float, str]] = []
        for study in sorted(values):
            rest = [value for other, value in values.items() if other != study]
            r_up, r_down, _t, r_p = sign_test(rest)
            if r_p is None:
                continue
            share = r_up / (r_up + r_down)
            shares.append(share)
            p_values.append(r_p)
            family = list(base)
            family[index] = r_p
            q = float(false_discovery_control(family, method="by")[index])
            flips += (q <= REPORT_Q) != (full_q[index] <= REPORT_Q)
            if full_share is not None:
                shifts.append((abs(share - full_share), study))
        top = ", ".join(study for _delta, study in sorted(shifts, reverse=True)[:3])
        results.append(
            Influence(
                numerator=ratio[0],
                reference=ratio[1],
                studies=len(values),
                full_share_up=full_share,
                full_q=full_q[index] if full_p[index] is not None else None,
                min_share_up=min(shares, default=None),
                max_share_up=max(shares, default=None),
                min_p=min(p_values, default=None),
                max_p=max(p_values, default=None),
                removals_changing_q_status=flips,
                most_influential_studies=top,
            )
        )
    return results


@dataclass(frozen=True)
class ReportedTaxa:
    """Reported taxa for one experiment, both directions together."""

    taxa: int
    species: int


def read_reported_taxa(path: Path) -> dict[tuple[str, str], ReportedTaxa]:
    """Count reported taxa and species-rank taxa per experiment from the taxa TSV."""
    taxa: Counter[tuple[str, str]] = Counter()
    species: Counter[tuple[str, str]] = Counter()
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["study_id"], row["experiment_id"])
            taxa[key] += 1
            species[key] += row["taxonomic_rank"] == "species"
    return {key: ReportedTaxa(taxa[key], species[key]) for key in taxa}


def _median(values: Sequence[float]) -> float | None:
    return median(values) if values else None


def eligibility_tables(
    pairs: Sequence[Pair],
    traits: Sequence[str],
    reported: dict[tuple[str, str], ReportedTaxa],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Who each ratio keeps: per ratio overall, and per ratio by body site and condition.

    Annotation coverage here is what the profile matrix records - the retained
    annotation total and the numbers of reported taxa - since the per-profile
    coverage sidecar of #184 is not yet available.
    """
    overall: list[dict[str, object]] = []
    strata: list[dict[str, object]] = []
    for numerator, reference in combinations(traits, 2):
        eligible = [p for p in pairs if p.positive((numerator, reference))]
        excluded = [p for p in pairs if not p.positive((numerator, reference))]
        kept_studies = {p.study_id for p in eligible}
        all_studies = {p.study_id for p in pairs}
        partial = {p.study_id for p in excluded} & kept_studies

        def describe(group: Sequence[Pair]) -> dict[str, float | None]:
            known = [
                reported[(p.study_id, p.experiment_id)]
                for p in group
                if (p.study_id, p.experiment_id) in reported
            ]
            return {
                "median_total_annotations": _median(
                    [float(p.inc.total + p.dec.total) for p in group]
                ),
                "median_reported_taxa": _median([float(r.taxa) for r in known]),
                "median_species_share": _median([r.species / r.taxa for r in known if r.taxa]),
                "pairs_without_reported_taxa_record": float(len(group) - len(known)),
            }

        kept, dropped = describe(eligible), describe(excluded)
        overall.append(
            {
                "numerator": numerator,
                "reference": reference,
                "studies_total": len(all_studies),
                "studies_eligible": len(kept_studies),
                "studies_partially_observed": len(partial),
                "studies_entirely_excluded": len(all_studies - kept_studies),
                "pairs_total": len(pairs),
                "pairs_eligible": len(eligible),
                "pairs_zero_numerator": sum(
                    1
                    for p in excluded
                    if min(p.inc.guild_counts[numerator], p.dec.guild_counts[numerator]) == 0
                ),
                "pairs_zero_reference": sum(
                    1
                    for p in excluded
                    if min(p.inc.guild_counts[reference], p.dec.guild_counts[reference]) == 0
                ),
                **{f"eligible_{k}": v for k, v in kept.items()},
                **{f"excluded_{k}": v for k, v in dropped.items()},
                "coverage_sidecar_status": "pending #184",
            }
        )
        for kind in ("body_site", "condition"):
            grouped: defaultdict[str, list[Pair]] = defaultdict(list)
            for pair in pairs:
                grouped[pair.body_site if kind == "body_site" else pair.condition].append(pair)
            for label, members in sorted(grouped.items()):
                kept_members = [p for p in members if p.positive((numerator, reference))]
                strata.append(
                    {
                        "numerator": numerator,
                        "reference": reference,
                        "stratum_kind": kind,
                        "stratum": label,
                        "studies_total": len({p.study_id for p in members}),
                        "studies_eligible": len({p.study_id for p in kept_members}),
                        "pairs_total": len(members),
                        "pairs_eligible": len(kept_members),
                    }
                )
    return overall, strata


def _format(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def _write_rows(path: Path, rows: Sequence[dict[str, object]]) -> Path:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        if rows:
            writer.writerow(list(rows[0]))
            writer.writerows([_format(value) for value in row.values()] for row in rows)
    return path


def _sha256(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def _display(value: float | None, spec: str = ".3g") -> str:
    return "NA" if value is None else format(value, spec)


def write_outputs(
    rows: Sequence[Row],
    influence: Sequence[Influence],
    overall: Sequence[dict[str, object]],
    strata: Sequence[dict[str, object]],
    output_dir: Path,
    provenance: dict[str, object],
    *,
    population: Sequence[PopulationComparison],
) -> list[Path]:
    """Write every scenario row, the influence and eligibility tables, status and a report."""
    output_dir.mkdir(parents=True, exist_ok=True)
    written = [
        _write_rows(output_dir / "sensitivity_results.tsv", [asdict(row) for row in rows]),
        _write_rows(output_dir / "leave_one_study_out.tsv", [asdict(item) for item in influence]),
        _write_rows(output_dir / "eligibility_by_ratio.tsv", list(overall)),
        _write_rows(output_dir / "eligibility_by_stratum.tsv", list(strata)),
        _write_rows(
            output_dir / "population_sensitivity.tsv", [asdict(item) for item in population]
        ),
    ]
    by_scenario: defaultdict[str, Counter[str]] = defaultdict(Counter)
    for row in rows:
        by_scenario[row.scenario][row.status] += 1
        if row.q_value is not None and row.q_value <= REPORT_Q:
            by_scenario[row.scenario]["q_at_or_below_0.05"] += 1
    status = {
        **provenance,
        "version": VERSION,
        "protocol": str(PROTOCOL_PATH.relative_to(REPO_ROOT)),
        "protocol_sha256": _sha256(PROTOCOL_PATH),
        "row_fields": [field.name for field in fields(Row)],
        "scenario_status": {name: dict(counts) for name, counts in sorted(by_scenario.items())},
        "population_comparison_summary": dict(
            Counter(item.detection_change for item in population)
        ),
        "population_estimate_agreement": dict(
            Counter(item.estimate_agreement for item in population)
        ),
        "calibration_status": "not assessed (#182)",
        "orientation_status": "not verified as disease versus control (#188)",
        "coverage_status": "profile totals and reported taxa only; #184 sidecar pending",
    }
    status_path = output_dir / "sensitivity_status.json"
    status_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    written.append(status_path)

    lines = [
        "# Log-ratio sensitivity (generated)",
        "",
        f"Version {VERSION}; scenarios as declared in `{status['protocol']}`. Each scenario "
        "keeps its own correction family. Effects are study-level means; shares are of "
        "non-tied studies with an exact 95% interval for that share only. NA is undefined, "
        "not zero.",
        "",
        "| Scenario | Comparison | Numerator | Reference | Studies eligible/total | Pairs eligible "
        "| Up/down/tied | Mean effect | Share up (95% CI) | q (family) | Status |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in rows:
        lines.append(
            f"| {row.scenario} | {row.comparison} | {row.numerator} | {row.reference} | "
            f"{row.studies_eligible}/{row.studies_total} | {row.pairs_eligible} | "
            f"{row.studies_up}/{row.studies_down}/{row.studies_tied} | "
            f"{_display(row.mean_effect)} | {_display(row.share_up)} "
            f"({_display(row.ci_low)}, {_display(row.ci_high)}) | "
            f"{_display(row.q_value)} ({row.family}) | {row.status} |"
        )
    report = output_dir / "sensitivity_report.md"
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    written.append(report)
    return written


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile-input", type=Path, default=DEFAULT_PROFILE_PATH)
    parser.add_argument("--experiments-input", type=Path, default=DEFAULT_EXPERIMENTS_PATH)
    parser.add_argument("--full-dump-input", type=Path, default=DEFAULT_FULL_DUMP_PATH)
    parser.add_argument("--taxa-input", type=Path, default=DEFAULT_TAXA_PATH)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--code-commit", help="Source revision, supplied by the Justfile recipe")
    args = parser.parse_args(argv)

    traits = tuple(TABLE1_TRAITS)
    profiles, stats = load_profiles(
        args.profile_input, args.experiments_input, TABLE1_TRAITS, reject_missing_counts=True
    )
    vintage = read_bugsigdb_vintage(args.full_dump_input)
    write_stats_report(stats, vintage)
    pairs = complete_pairs(profiles, traits)
    rows = run_scenarios(pairs, traits)
    influence = leave_one_study_out(pairs, traits)
    population = population_comparison(rows, traits)
    overall, strata = eligibility_tables(pairs, traits, read_reported_taxa(args.taxa_input))
    provenance: dict[str, object] = {
        "bugsigdb_vintage": vintage,
        "code_commit": args.code_commit,
        "profile_sha256": _sha256(args.profile_input),
        "experiments_sha256": _sha256(args.experiments_input),
        "taxa_sha256": _sha256(args.taxa_input),
        "analysis_source_sha256": _sha256(Path(__file__)),
        "corpus": asdict(stats),
        "complete_pairs": len(pairs),
        "studies": len({pair.study_id for pair in pairs}),
    }
    for path in write_outputs(
        rows, influence, overall, strata, args.output_dir, provenance, population=population
    ):
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
