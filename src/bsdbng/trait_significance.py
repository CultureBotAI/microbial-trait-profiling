"""Evidence overlays for the Table 1 condition-by-body-site heatmap (#172).

The heatmap shows how large a trait's annotation share is. This module attaches
three separately defined kinds of evidence to it and keeps them apart in every
table, caption and legend (#177):

* **Aggregated profile-level support.** The enrichment arm's per-profile
  hypergeometric test (Analysis B, ``scripts/trait_enrichment.py`` at
  ``c75744a``) is re-run here and reproduces that arm's per-trait counts. A cell's
  profile calls are then summarised by a descriptive study-support rule
  (:data:`SUPPORT_RULE`). That rule is not a statistical test, and no cell-level
  p- or q-value is produced for it.
* **Body-site composition contrast.** A paired increased-minus-decreased test per
  exact source body-site label, modelled on the composition arm (Analysis A,
  ``scripts/create_trait_profile_infographic.py`` at ``3afeecd``). It is a
  reimplementation, not a transcription: :data:`COMPOSITION_DEPARTURES` lists
  every difference, and the figure and methods text repeat them.
* **Condition-by-body-site composition contrast.** The same test inside one exact
  condition and body-site label.

Inference does not depend on the figure (#176). Every test runs on exact source
labels over the corpus :func:`bsdbng.trait_guilds.load_profiles` returns, and
correction families are fixed before anything is drawn. The heatmap's top-N cuts,
masking and pooled Other groups only decide where a result is shown afterwards.
A study contributes one value to a composition test however many experiments or
conditions it has.

The six Table 1 traits were chosen using this same corpus, so none of these
results confirms them independently.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import defaultdict
from collections.abc import Callable, Collection, Iterable, Mapping, Sequence
from dataclasses import dataclass, field, replace
from pathlib import Path
from statistics import fmean

from scipy.stats import false_discovery_control, hypergeom, ttest_1samp

from bsdbng.trait_guilds import (
    DIRECTIONS,
    MEDIADIVE_MARKER,
    PROFILE_ID_COLUMNS,
    Profile,
    column_term_id,
)

#: Thresholds copied from ``scripts/trait_enrichment.py``'s defaults at ``c75744a``.
#: The arm has no depletion call, so none is made here.
ENRICHMENT_Q = 0.001
ENRICHMENT_MIN_FOLD = 1.5

#: The composition arm publishes q-values without a cutoff. This one is chosen by
#: this integration for display, and the figure says so.
COMPOSITION_Q = 0.05
#: The composition arm refuses a column with fewer than three values; here the
#: values are unique studies.
MIN_STUDIES = 3
#: The composition arm's stratum floor (``MIN_STRATIFIED_PAIRS``), in paired
#: experiments with a positive total on both sides.
MIN_SITE_PAIRS = 50
#: The composition arm's trait selection within a stratum
#: (``select_stratum_traits``): present in at least this many profiles, in at most
#: this share of them, and with nonzero variance of its share.
MIN_PROFILE_SUPPORT = 6
MAX_PROFILE_SUPPORT_FRACTION = 0.95

#: The descriptive support rule. A study supports a cell when at least one of its
#: profiles in the cell is enriched; a cell has majority support when more than
#: half of at least two eligible studies do.
SUPPORT_SHARE = 0.5
MIN_SUPPORT_STUDIES = 2
SUPPORT_RULE = (
    "strict majority of eligible studies (at least 2) have at least one profile in the "
    f"cell enriched under the profile-level test (q <= {ENRICHMENT_Q}, "
    f"fold >= {ENRICHMENT_MIN_FOLD})"
)
#: The sensitivity alternative from #177: a study supports a cell only when more
#: than half of its own eligible profiles there are enriched, so a study with many
#: experiments gets no more chances than a study with one.
ALTERNATIVE_SUPPORT_RULE = (
    "strict majority of eligible studies (at least 2) have a strict majority of their "
    "own eligible profiles in the cell enriched"
)

#: Evidence types. They are never merged into one another.
SUPPORT = "aggregated_profile_support"
SITE_COMPOSITION = "site_composition"
CELL_COMPOSITION = "condition_site_composition"

#: Support states for a cell.
MAJORITY_SUPPORT = "majority_support"
NO_MAJORITY = "no_majority"
INSUFFICIENT_STUDIES = "insufficient_studies"

#: Composition test states.
TESTED = "tested"
NOT_SELECTED = "trait_not_selected_in_site"
NO_VARIATION = "no_variation"
SITE_BELOW_MIN_PAIRS = "site_below_min_pairs"

#: Every way this composition test differs from Analysis A at ``3afeecd``, as
#: (aspect, the arm's committed run, this reimplementation).
COMPOSITION_DEPARTURES: tuple[tuple[str, str, str], ...] = (
    (
        "sampling unit",
        "each matched experiment pair is one observation",
        "experiment contrasts are averaged to one value per unique study, across all of "
        "its conditions in the stratum (#133, #176)",
    ),
    (
        "strata",
        "body_site_group labels with at least 50 pairs, plus Feces and Feces-excluded",
        f"exact source body-site labels with at least {MIN_SITE_PAIRS} eligible pairs; no "
        "pooled Other stratum, since that pool is defined by the figure's top-N cut",
    ),
    (
        "denominator",
        "sum of every column, MediaDive growth media included (#174)",
        "sum of the non-MediaDive trait columns, as the shared comparison basis specifies",
    ),
    (
        "zero totals",
        "a zero-total profile is scored as share 0",
        "a pair with a zero total on either side is excluded; its shares are undefined",
    ),
    (
        "correction family",
        "Benjamini-Hochberg within each stratum over every trait the stratum selects",
        "the same selection and family per body site; the six-trait family is reported "
        "only as a labelled sensitivity column",
    ),
    (
        "cutoff",
        "q-values published without a cutoff",
        f"q <= {COMPOSITION_Q} chosen here for display",
    ),
)

#: Keyed by (condition, body site, trait, direction), all exact source labels.
CellKey = tuple[str, str, str, str]
ProfileKey = tuple[str, str, str]


@dataclass(frozen=True)
class CountTable:
    """Validated counts for the retained trait columns of the profiles in use."""

    columns: tuple[str, ...]
    counts: dict[ProfileKey, tuple[int, ...]]
    totals: dict[ProfileKey, int] = field(init=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "totals", {key: sum(v) for key, v in self.counts.items()})

    def total(self, key: ProfileKey) -> int:
        return self.totals[key]


def read_count_table(profile_path: Path, keys: Collection[ProfileKey]) -> CountTable:
    """Read the non-MediaDive trait counts of ``keys``, rejecting malformed input.

    A blank, non-integer or negative cell raises rather than becoming zero, and so
    does a row of the wrong width or a repeated key: a silent zero would change
    every denominator and background downstream (#184).
    """
    wanted = set(keys)
    with profile_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if header is None:
            raise ValueError(f"{profile_path} is empty")
        if header[: len(PROFILE_ID_COLUMNS)] != list(PROFILE_ID_COLUMNS):
            raise ValueError(f"{profile_path} must start with {PROFILE_ID_COLUMNS}")
        positions = [
            index
            for index, column in enumerate(header)
            if column not in PROFILE_ID_COLUMNS and MEDIADIVE_MARKER not in column.lower()
        ]
        counts: dict[ProfileKey, tuple[int, ...]] = {}
        for line, row in enumerate(reader, start=2):
            if not row:
                continue
            if len(row) != len(header):
                raise ValueError(
                    f"{profile_path}:{line} has {len(row)} fields; the header has {len(header)}"
                )
            key = (row[0], row[1], row[2])
            if key not in wanted:
                continue
            if key in counts:
                raise ValueError(f"{profile_path}:{line} repeats profile {key}")
            values: list[int] = []
            for position in positions:
                raw = row[position].strip()
                if not (raw.isascii() and raw.isdigit()):
                    raise ValueError(
                        f"{profile_path}:{line} column {header[position]!r} holds {raw!r}; "
                        "counts must be non-negative integers and missing is not zero"
                    )
                values.append(int(raw))
            counts[key] = tuple(values)
    missing = wanted.difference(counts)
    if missing:
        raise ValueError(f"{len(missing)} profiles are absent from {profile_path}")
    return CountTable(columns=tuple(header[p] for p in positions), counts=counts)


def resolve_trait_slots(
    columns: Sequence[str], traits: Mapping[str, Sequence[str]]
) -> dict[str, int]:
    """Map each named trait to its column slot.

    Two Table 1 traits are carried behind a METPO qualifier, so the segment after
    the final pipe is matched as well as the whole column name.
    """
    slots: dict[str, int] = {}
    for slot, column in enumerate(columns):
        term = column_term_id(column)
        for name, ids in traits.items():
            if column in ids or term in ids:
                if name in slots:
                    raise ValueError(f"trait {name!r} matches more than one column")
                slots[name] = slot
    return slots


@dataclass(frozen=True)
class EnrichmentCall:
    """One profile's result for one trait under the profile-level enrichment test."""

    study_id: str
    experiment_id: str
    inc_dec: str
    trait: str
    count: int
    expected: float
    fold: float
    p_value: float
    q_value: float

    @property
    def enriched(self) -> bool:
        return self.q_value <= ENRICHMENT_Q and self.fold >= ENRICHMENT_MIN_FOLD


@dataclass(frozen=True)
class EnrichmentRun:
    calls: list[EnrichmentCall]
    untested_traits: list[str]
    profiles_tested: int
    tests_performed: int


def enrichment_calls(table: CountTable, traits: Mapping[str, Sequence[str]]) -> EnrichmentRun:
    """Run the enrichment arm's hypergeometric test over every profile in ``table``.

    The background is every profile in the table with a positive total, fixed
    before any display choice, as the shared comparison basis requires. Every
    observed positive profile-trait count is tested and Benjamini-Hochberg runs
    across all of them at once, so the family is the arm's full family rather than
    the six drawn traits. A profile with no retained annotation at all enters
    neither the background nor the tests.
    """
    slots = resolve_trait_slots(table.columns, traits)
    drawn = {slot: name for name, slot in slots.items()}
    width = len(table.columns)
    background = [0] * width
    rows: list[tuple[ProfileKey, int, list[tuple[int, int]]]] = []
    for key, values in table.counts.items():
        observed = [(slot, count) for slot, count in enumerate(values) if count > 0]
        if not observed:
            continue
        for slot, count in observed:
            background[slot] += count
        rows.append((key, sum(values), observed))
    background_total = sum(background)
    if not rows:
        return EnrichmentRun([], sorted(traits), 0, 0)

    records: list[tuple[ProfileKey, int, int, float, float]] = []
    p_values: list[float] = []
    for key, row_total, observed in rows:
        for slot, count in observed:
            successes = background[slot]
            # scipy's sf(k - 1, M, K, N) is P[X >= k], the arm's one-tailed test.
            p_value = float(hypergeom.sf(count - 1, background_total, successes, row_total))
            p_values.append(p_value)
            records.append((key, slot, count, p_value, row_total * successes / background_total))

    q_values = false_discovery_control(p_values, method="bh")
    calls = [
        EnrichmentCall(
            study_id=key[0],
            experiment_id=key[1],
            inc_dec=key[2],
            trait=drawn[slot],
            count=count,
            expected=expected,
            fold=count / expected,
            p_value=p_value,
            q_value=float(q_value),
        )
        for (key, slot, count, p_value, expected), q_value in zip(records, q_values, strict=True)
        if slot in drawn
    ]
    seen = {call.trait for call in calls}
    return EnrichmentRun(
        calls=calls,
        untested_traits=sorted(name for name in traits if name not in seen),
        profiles_tested=len(rows),
        tests_performed=len(p_values),
    )


@dataclass(frozen=True)
class SupportRecord:
    """Aggregated profile-level support for one cell. It carries no p- or q-value."""

    condition: str
    body_site: str
    trait: str
    direction: str
    eligible_studies: int
    supporting_studies: int
    alternative_supporting_studies: int
    eligible_profiles: int
    enriched_profiles: int
    zero_total_profiles: int
    studies_with_repeated_profiles: int
    max_profiles_per_study: int

    @staticmethod
    def _status(supporting: int, eligible: int) -> str:
        if eligible < MIN_SUPPORT_STUDIES:
            return INSUFFICIENT_STUDIES
        return MAJORITY_SUPPORT if supporting > SUPPORT_SHARE * eligible else NO_MAJORITY

    @property
    def status(self) -> str:
        return self._status(self.supporting_studies, self.eligible_studies)

    @property
    def alternative_status(self) -> str:
        return self._status(self.alternative_supporting_studies, self.eligible_studies)

    @property
    def support_share(self) -> float | None:
        return self.supporting_studies / self.eligible_studies if self.eligible_studies else None


def aggregate_support(
    profiles: Sequence[Profile], table: CountTable, run: EnrichmentRun, traits: Sequence[str]
) -> dict[CellKey, SupportRecord]:
    """Summarise profile calls per exact (condition, body site, trait, direction).

    A profile is eligible when its total is positive. One whose trait count is
    zero is eligible but cannot be enriched, which is what the arm's test implies.
    """
    enriched = {
        (call.study_id, call.experiment_id, call.inc_dec, call.trait)
        for call in run.calls
        if call.enriched
    }
    # cell -> study -> (eligible profiles, enriched profiles)
    per_study: defaultdict[CellKey, defaultdict[str, list[int]]] = defaultdict(
        lambda: defaultdict(lambda: [0, 0])
    )
    zero_totals: defaultdict[CellKey, int] = defaultdict(int)
    for profile in profiles:
        key = (profile.study_id, profile.experiment_id, profile.inc_dec)
        positive = table.total(key) > 0
        for trait in traits:
            cell = (profile.condition, profile.body_site, trait, profile.inc_dec)
            if not positive:
                zero_totals[cell] += 1
                continue
            tally = per_study[cell][profile.study_id]
            tally[0] += 1
            tally[1] += (*key, trait) in enriched

    records: dict[CellKey, SupportRecord] = {}
    for cell in sorted(set(per_study) | set(zero_totals)):
        tallies: list[list[int]] = list(per_study[cell].values()) if cell in per_study else []
        records[cell] = SupportRecord(
            condition=cell[0],
            body_site=cell[1],
            trait=cell[2],
            direction=cell[3],
            eligible_studies=len(tallies),
            supporting_studies=sum(1 for _, hits in tallies if hits > 0),
            alternative_supporting_studies=sum(1 for n, hits in tallies if hits > n / 2),
            eligible_profiles=sum(n for n, _ in tallies),
            enriched_profiles=sum(hits for _, hits in tallies),
            zero_total_profiles=zero_totals.get(cell, 0),
            studies_with_repeated_profiles=sum(1 for n, _ in tallies if n > 1),
            max_profiles_per_study=max((n for n, _ in tallies), default=0),
        )
    return records


@dataclass(frozen=True)
class CompositionResult:
    """One paired composition contrast. Site and cell results stay separate types."""

    evidence_type: str
    body_site: str
    condition: str | None
    trait: str
    column: str
    status: str
    studies: int
    paired_experiments: int
    conditions: int
    zero_total_pairs_excluded: int
    mean_contrast_z: float | None
    t_statistic: float | None
    p_value: float | None
    family: str
    family_size: int
    q_value: float | None
    table1_family_size: int
    table1_q_value: float | None

    @property
    def significant(self) -> bool:
        return self.q_value is not None and self.q_value <= COMPOSITION_Q

    @property
    def direction(self) -> str | None:
        if self.mean_contrast_z is None:
            return None
        return "increased" if self.mean_contrast_z >= 0 else "decreased"


def _transformed(count: int, total: int) -> float:
    """The composition arm's scale for one share: per 1,000, then log1p."""
    return math.log1p(1000.0 * count / total)


def _selected_slots(
    stratum: Sequence[tuple[ProfileKey, ProfileKey]],
    table: CountTable,
    min_profile_support: int = MIN_PROFILE_SUPPORT,
) -> dict[int, float]:
    """Traits the composition arm would select in a stratum, with their z-score scale.

    Selection mirrors ``select_stratum_traits``; the returned value is the
    population standard deviation of the transformed share over the stratum's
    profiles, which is the z-score divisor.
    """
    keys = [key for pair in stratum for key in pair]
    totals = [table.total(key) for key in keys]
    ceiling = max(2, int(len(keys) * MAX_PROFILE_SUPPORT_FRACTION))
    scales: dict[int, float] = {}
    for slot in range(len(table.columns)):
        counts = [table.counts[key][slot] for key in keys]
        support = sum(1 for count in counts if count > 0)
        if support < min_profile_support or support > ceiling:
            continue
        values = [_transformed(count, total) for count, total in zip(counts, totals, strict=True)]
        mean = fmean(values)
        deviation = math.sqrt(fmean([(value - mean) ** 2 for value in values]))
        if deviation > 0.0:
            scales[slot] = deviation
    return scales


def _study_values(
    pairs: Sequence[tuple[str, ProfileKey, ProfileKey]],
    table: CountTable,
    slot: int,
    scale: float,
) -> list[float]:
    """One mean z-scale contrast per unique study, pooling all of its pairs."""
    by_study: defaultdict[str, list[float]] = defaultdict(list)
    for study, inc, dec in pairs:
        up = _transformed(table.counts[inc][slot], table.total(inc))
        down = _transformed(table.counts[dec][slot], table.total(dec))
        by_study[study].append((up - down) / scale)
    return [fmean(values) for values in by_study.values()]


def _t_test(values: Sequence[float]) -> tuple[str, float | None, float | None, float | None]:
    if len(values) < MIN_STUDIES:
        return INSUFFICIENT_STUDIES, None, None, None
    mean = fmean(values)
    if all(value == values[0] for value in values):
        return NO_VARIATION, mean, None, None
    result = ttest_1samp(values, 0.0)
    return TESTED, mean, float(result.statistic), float(result.pvalue)


def _bh(
    results: list[CompositionResult], family_of: Callable[[CompositionResult], str]
) -> list[CompositionResult]:
    """Benjamini-Hochberg within each family over tested results only."""
    families: defaultdict[str, list[int]] = defaultdict(list)
    for index, result in enumerate(results):
        if result.p_value is not None:
            families[family_of(result)].append(index)
    adjusted = list(results)
    for family, indices in families.items():
        q_values = false_discovery_control([results[i].p_value for i in indices], method="bh")
        for index, q_value in zip(indices, q_values, strict=True):
            adjusted[index] = replace(
                adjusted[index], family=family, family_size=len(indices), q_value=float(q_value)
            )
    return adjusted


StudyPair = tuple[str, ProfileKey, ProfileKey]


@dataclass(frozen=True)
class _Stratum:
    """The eligible pairs of one exact site, or of one condition within it."""

    body_site: str
    condition: str | None
    pairs: list[StudyPair]
    conditions: int
    zero_total_pairs_excluded: int

    @property
    def studies(self) -> int:
        return len({study for study, _, _ in self.pairs})

    def result(
        self, status: str, slot: int, table: CountTable, names: Mapping[int, str]
    ) -> CompositionResult:
        return CompositionResult(
            evidence_type=SITE_COMPOSITION if self.condition is None else CELL_COMPOSITION,
            body_site=self.body_site,
            condition=self.condition,
            trait=names.get(slot, table.columns[slot]),
            column=table.columns[slot],
            status=status,
            studies=self.studies,
            paired_experiments=len(self.pairs),
            conditions=self.conditions,
            zero_total_pairs_excluded=self.zero_total_pairs_excluded,
            mean_contrast_z=None,
            t_statistic=None,
            p_value=None,
            family="",
            family_size=0,
            q_value=None,
            table1_family_size=0,
            table1_q_value=None,
        )

    def tested(
        self, slot: int, scale: float, table: CountTable, names: Mapping[int, str]
    ) -> CompositionResult:
        status, mean, statistic, p_value = _t_test(_study_values(self.pairs, table, slot, scale))
        return replace(
            self.result(status, slot, table, names),
            mean_contrast_z=mean,
            t_statistic=statistic,
            p_value=p_value,
        )


def composition_contrasts(
    profiles: Sequence[Profile],
    table: CountTable,
    traits: Mapping[str, Sequence[str]],
    *,
    min_site_pairs: int = MIN_SITE_PAIRS,
    min_profile_support: int = MIN_PROFILE_SUPPORT,
) -> list[CompositionResult]:
    """Run the paired composition contrast per exact site and per exact cell.

    Pairs are formed from the loader's paired corpus and grouped by exact source
    labels only. For a site, every eligible pair in that site enters its study's
    mean, whatever the pair's condition, so a study with several conditions is
    still one observation. Cells use the site's trait selection and z-score scale,
    and one Benjamini-Hochberg family spans every tested cell.
    """
    table1_slots = resolve_trait_slots(table.columns, traits)
    table1_names = {slot: name for name, slot in table1_slots.items()}

    by_experiment: defaultdict[tuple[str, str], dict[str, Profile]] = defaultdict(dict)
    for profile in profiles:
        by_experiment[(profile.study_id, profile.experiment_id)][profile.inc_dec] = profile

    # site -> list of (condition, study, inc key, dec key); excluded zero-total pairs.
    site_pairs: defaultdict[str, list[tuple[str, str, ProfileKey, ProfileKey]]] = defaultdict(list)
    zero_pairs: defaultdict[tuple[str, str | None], int] = defaultdict(int)
    for (study, experiment), sides in sorted(by_experiment.items()):
        if set(sides) != set(DIRECTIONS):
            raise ValueError(f"unpaired experiment {(study, experiment)}; load a paired corpus")
        inc, dec = sides["inc"], sides["dec"]
        if (inc.condition, inc.body_site) != (dec.condition, dec.body_site):
            raise ValueError(f"experiment {(study, experiment)} has mismatched metadata")
        inc_key = (study, experiment, "inc")
        dec_key = (study, experiment, "dec")
        if table.total(inc_key) == 0 or table.total(dec_key) == 0:
            zero_pairs[(inc.body_site, None)] += 1
            zero_pairs[(inc.body_site, inc.condition)] += 1
            continue
        site_pairs[inc.body_site].append((inc.condition, study, inc_key, dec_key))

    results: list[CompositionResult] = []
    for site, pairs in sorted(site_pairs.items()):
        by_condition: defaultdict[str, list[StudyPair]] = defaultdict(list)
        for condition, study, inc_key, dec_key in pairs:
            by_condition[condition].append((study, inc_key, dec_key))
        flat = [(study, inc, dec) for _, study, inc, dec in pairs]
        site_row = _Stratum(site, None, flat, len(by_condition), zero_pairs.get((site, None), 0))

        if len(pairs) < min_site_pairs:
            # Recorded for the six traits only, so the state is visible without
            # listing every untested trait of every small site.
            results.extend(
                site_row.result(SITE_BELOW_MIN_PAIRS, slot, table, table1_names)
                for slot in sorted(table1_names)
            )
            continue

        scales = _selected_slots(
            [(inc, dec) for _, _, inc, dec in pairs], table, min_profile_support
        )
        cells = [
            _Stratum(site, condition, cell_pairs, 1, zero_pairs.get((site, condition), 0))
            for condition, cell_pairs in sorted(by_condition.items())
        ]
        for slot in sorted(set(scales) | set(table1_names)):
            if slot not in scales:
                results.append(site_row.result(NOT_SELECTED, slot, table, table1_names))
                continue
            results.append(site_row.tested(slot, scales[slot], table, table1_names))
            for cell in cells:
                if cell.studies < MIN_STUDIES and slot not in table1_names:
                    continue
                results.append(cell.tested(slot, scales[slot], table, table1_names))

    results = _bh(
        results,
        lambda r: f"site:{r.body_site}" if r.evidence_type == SITE_COMPOSITION else "all_cells",
    )
    # The six-trait families, kept only as labelled sensitivity columns.
    table1_columns = {table.columns[slot] for slot in table1_names}
    six = [i for i, r in enumerate(results) if r.p_value is not None and r.column in table1_columns]
    for evidence_type in (SITE_COMPOSITION, CELL_COMPOSITION):
        indices = [i for i in six if results[i].evidence_type == evidence_type]
        if not indices:
            continue
        q_values = false_discovery_control([results[i].p_value for i in indices], method="bh")
        for index, q_value in zip(indices, q_values, strict=True):
            results[index] = replace(
                results[index], table1_family_size=len(indices), table1_q_value=float(q_value)
            )
    return results


@dataclass(frozen=True)
class EvidenceOverlay:
    """All three evidence types, keyed by exact source labels, before any display."""

    support: dict[CellKey, SupportRecord]
    composition: list[CompositionResult]
    enrichment: EnrichmentRun
    traits: tuple[str, ...]

    def table1_composition(self) -> list[CompositionResult]:
        return [r for r in self.composition if r.trait in self.traits]

    def cell_marks(self) -> dict[CellKey, frozenset[str]]:
        """Evidence types meeting their rule per cell, for cells anywhere in the corpus.

        Support marks the direction it was aggregated over. A condition-level
        composition result is one contrast between the two directions, so it marks
        only the side it favours. Site results never mark a cell.
        """
        marks: defaultdict[CellKey, set[str]] = defaultdict(set)
        for cell, record in self.support.items():
            if record.status == MAJORITY_SUPPORT:
                marks[cell].add(SUPPORT)
        for result in self.table1_composition():
            if (
                result.evidence_type == CELL_COMPOSITION
                and result.significant
                and result.condition is not None
                and result.direction is not None
            ):
                side = "inc" if result.direction == "increased" else "dec"
                marks[(result.condition, result.body_site, result.trait, side)].add(
                    CELL_COMPOSITION
                )
        return {cell: frozenset(kinds) for cell, kinds in marks.items()}

    def site_marks(self) -> dict[tuple[str, str], str]:
        """(exact body site, trait) -> the direction a significant site contrast favours."""
        return {
            (result.body_site, result.trait): result.direction
            for result in self.table1_composition()
            if result.evidence_type == SITE_COMPOSITION
            and result.significant
            and result.direction is not None
        }


def build_overlay(
    profiles: Sequence[Profile],
    traits: Mapping[str, Sequence[str]],
    profile_path: Path,
) -> EvidenceOverlay:
    """Compute every evidence type over ``profiles`` with no reference to a figure."""
    keys = [(p.study_id, p.experiment_id, p.inc_dec) for p in profiles]
    table = read_count_table(profile_path, keys)
    run = enrichment_calls(table, traits)
    return EvidenceOverlay(
        support=aggregate_support(profiles, table, run, list(traits)),
        composition=composition_contrasts(profiles, table, traits),
        enrichment=run,
        traits=tuple(traits),
    )


def _cell(value: float | int | str | None) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def write_evidence_tables(
    overlay: EvidenceOverlay,
    output_dir: Path,
    prefix: str,
    *,
    display_status: Callable[[CellKey], str],
    site_drawn: Callable[[str], bool],
    provenance: Mapping[str, object],
) -> list[Path]:
    """Write support, the six traits' contrasts, the full screen, and provenance.

    ``display_status`` and ``site_drawn`` describe the combined figure; they are
    written alongside the evidence and never feed back into it.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    support_path = output_dir / f"{prefix}_cell_support.tsv"
    contrasts_path = output_dir / f"{prefix}_composition_contrasts.tsv"
    screen_path = output_dir / f"{prefix}_composition_screen.tsv"
    provenance_path = output_dir / f"{prefix}_evidence_provenance.json"

    with support_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "condition",
                "body_site",
                "trait",
                "direction",
                "evidence_type",
                "support_status",
                "eligible_studies",
                "supporting_studies",
                "support_share",
                "alternative_support_status",
                "alternative_supporting_studies",
                "eligible_profiles",
                "enriched_profiles",
                "zero_total_profiles",
                "studies_with_repeated_profiles",
                "max_profiles_per_study",
                "display_status",
            ]
        )
        for cell, record in overlay.support.items():
            writer.writerow(
                [
                    record.condition,
                    record.body_site,
                    record.trait,
                    record.direction,
                    SUPPORT,
                    record.status,
                    record.eligible_studies,
                    record.supporting_studies,
                    _cell(record.support_share),
                    record.alternative_status,
                    record.alternative_supporting_studies,
                    record.eligible_profiles,
                    record.enriched_profiles,
                    record.zero_total_profiles,
                    record.studies_with_repeated_profiles,
                    record.max_profiles_per_study,
                    display_status(cell),
                ]
            )

    composition_header = [
        "evidence_type",
        "body_site",
        "condition",
        "trait",
        "column",
        "status",
        "studies",
        "paired_experiments",
        "conditions",
        "zero_total_pairs_excluded",
        "mean_contrast_z",
        "t_statistic",
        "p_value",
        "family",
        "family_size",
        "q_value",
        "significant",
        "table1_family_size",
        "table1_q_value",
        "display_status",
    ]

    def composition_row(result: CompositionResult) -> list[str]:
        if result.condition is None:
            shown = "site_drawn" if site_drawn(result.body_site) else "site_not_drawn"
        else:
            shown = display_status((result.condition, result.body_site, result.trait, "inc"))
        return [
            result.evidence_type,
            result.body_site,
            result.condition or "",
            result.trait,
            result.column,
            result.status,
            _cell(result.studies),
            _cell(result.paired_experiments),
            _cell(result.conditions),
            _cell(result.zero_total_pairs_excluded),
            _cell(result.mean_contrast_z),
            _cell(result.t_statistic),
            _cell(result.p_value),
            result.family,
            _cell(result.family_size),
            _cell(result.q_value),
            "true" if result.significant else "false",
            _cell(result.table1_family_size),
            _cell(result.table1_q_value),
            shown,
        ]

    ordered = sorted(
        overlay.composition,
        key=lambda r: (r.evidence_type, r.body_site, r.condition or "", r.trait),
    )
    for path, rows in (
        (contrasts_path, [r for r in ordered if r.trait in overlay.traits]),
        (screen_path, ordered),
    ):
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(composition_header)
            writer.writerows(composition_row(result) for result in rows)

    provenance_record: dict[str, object] = {
        **provenance,
        "evidence_types": {
            SUPPORT: {
                "profile_test": "one-tailed hypergeometric per observed positive profile-trait "
                "count; Benjamini-Hochberg across every test",
                "profile_thresholds": {"q": ENRICHMENT_Q, "min_fold": ENRICHMENT_MIN_FOLD},
                "background": "every profile in the analysed corpus with a positive "
                "non-MediaDive total",
                "profiles_tested": overlay.enrichment.profiles_tested,
                "tests_performed": overlay.enrichment.tests_performed,
                "rule": SUPPORT_RULE,
                "alternative_rule": ALTERNATIVE_SUPPORT_RULE,
                "cell_p_or_q": None,
                "source": "scripts/trait_enrichment.py at c75744a",
            },
            SITE_COMPOSITION: {
                "test": "one-sample t-test against 0 of one mean z-scale contrast per study",
                "family": "Benjamini-Hochberg per exact body site over the traits selected there",
                "min_site_pairs": MIN_SITE_PAIRS,
                "min_studies": MIN_STUDIES,
                "display_q": COMPOSITION_Q,
                "source": "scripts/create_trait_profile_infographic.py at 3afeecd",
            },
            CELL_COMPOSITION: {
                "test": "as site_composition, within one exact condition and body site",
                "family": "one Benjamini-Hochberg family over every tested cell",
                "min_studies": MIN_STUDIES,
                "display_q": COMPOSITION_Q,
            },
        },
        "composition_departures": [
            {"aspect": aspect, "composition_arm": theirs, "here": ours}
            for aspect, theirs, ours in COMPOSITION_DEPARTURES
        ],
        "selection_note": "the six Table 1 traits were selected on this corpus; no result "
        "here is independent confirmation of them",
    }
    provenance_path.write_text(json.dumps(provenance_record, indent=2) + "\n", encoding="utf-8")
    return [support_path, contrasts_path, screen_path, provenance_path]


def sha256(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def support_summary(records: Iterable[SupportRecord]) -> dict[str, int]:
    """Counts of cells per support state under both rules, for the run report."""
    summary: defaultdict[str, int] = defaultdict(int)
    for record in records:
        summary[record.status] += 1
        summary[f"alternative_{record.alternative_status}"] += 1
        if record.status != record.alternative_status:
            summary["changed_by_alternative_rule"] += 1
    return dict(summary)
