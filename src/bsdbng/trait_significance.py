"""Significance calls for the condition-by-body-site heatmap cells (#172).

The heatmap says how much of a profile's trait signal a trait accounts for. It
does not say whether that share is more than chance, and the project already has
two arms that answer exactly that question in two different ways. This module
runs both of them over the same corpus the figure is drawn from, so a cell can be
outlined when either arm calls it significant and the figure carries all three
approaches at once.

Neither test is invented here. Both are transcribed from the arm that owns it:

* **Enrichment (Analysis B, Ilyas Abdolcader,**
  ``scripts/trait_enrichment.py`` **on** ``origin/enrichment-analysis-heatmap-cluster``\ **).**
  One profile at a time, a one-tailed hypergeometric test asks whether that
  profile drew more assertions of a trait than a draw of its size from the corpus
  pool would give. Benjamini-Hochberg is applied once across every test, and a
  call needs ``q <= 0.001`` together with a fold of at least 1.5 (enriched) or at
  most 0.5 (depleted). Because it is per profile, it aggregates to a heatmap cell
  directly.

* **Composition (Analysis A, Jed Kim-Ozaeta,**
  ``scripts/create_trait_profile_infographic.py`` **on**
  ``origin/clustered-heatmap-generation``\ **).** Trait shares are scaled per
  1,000, ``log1p``-transformed, and z-scored within the stratum; the increased
  side of a matched pair minus its decreased side is tested against zero with a
  one-sample t-test. It compares the two directions against each other rather
  than against a corpus background, so it produces one verdict for a trait's
  increased and decreased columns together, not one per column.

Two deliberate departures from the committed runs, both stated on the figure:

* Pairs are averaged to one value per study before the composition t-test, which
  is the correction issue #133 asks for. 56.2% of studies in this corpus carry
  more than one experiment, so treating experiments as independent overstates the
  evidence.
* The composition test is run per condition inside a body site as well as over
  the whole site. The per-condition form is reported for completeness and is
  almost always underpowered, which is the finding, not a defect: see
  :func:`composition_calls`.

The enrichment reimplementation is checked against the committed table it was
transcribed from - see ``tests/test_trait_significance.py`` - and reproduces its
per-trait significant-profile counts exactly.
"""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from scipy.stats import false_discovery_control, hypergeom, ttest_1samp

from bsdbng.trait_guilds import (
    DEFAULT_PROFILE_PATH,
    DIRECTIONS,
    MEDIADIVE_MARKER,
    PROFILE_ID_COLUMNS,
    Profile,
    column_term_id,
)

#: Thresholds copied from ``scripts/trait_enrichment.py``'s defaults.
ENRICHMENT_Q = 0.001
ENRICHMENT_MIN_FOLD = 1.5
DEPLETION_MAX_FOLD = 0.5

#: The composition arm publishes q-values without naming a cut-off, so the
#: conventional one is used and stated on the figure.
COMPOSITION_Q = 0.05
#: ``compute_trait_stats`` refuses a column with fewer than three values.
MIN_PAIRS = 3

#: A cell is outlined when more than this share of the studies behind it were
#: called significant on their own. A simple majority keeps one study out of
#: several from speaking for the cell, and a cell backed by a single study is
#: never outlined at all.
CELL_SUPPORT_SHARE = 0.5
MIN_CELL_STUDIES = 2

#: What a border means. ``both`` is drawn where a cell is called by the
#: enrichment arm and its body site is also called by the composition arm.
ENRICHMENT = "enrichment"
COMPOSITION = "composition"
BOTH = "both"

#: Keyed by (condition, body site, trait, direction).
CellKey = tuple[str, str, str, str]


@dataclass(frozen=True)
class EnrichmentCall:
    """One profile's verdict on one trait, under the enrichment arm."""

    study_id: str
    experiment_id: str
    inc_dec: str
    trait: str
    p_value: float
    q_value: float
    fold: float
    count: int
    expected: float

    @property
    def verdict(self) -> str | None:
        """``enriched``, ``depleted``, or nothing at this profile's evidence."""
        if self.q_value > ENRICHMENT_Q:
            return None
        if self.fold >= ENRICHMENT_MIN_FOLD:
            return "enriched"
        if self.fold <= DEPLETION_MAX_FOLD:
            return "depleted"
        return None


@dataclass(frozen=True)
class CompositionCall:
    """One trait's paired increased-minus-decreased contrast within a stratum."""

    stratum_kind: str
    body_site: str
    condition: str | None
    trait: str
    studies: int
    mean_contrast: float
    t_statistic: float
    p_value: float
    q_value: float

    @property
    def significant(self) -> bool:
        return self.q_value <= COMPOSITION_Q

    @property
    def direction(self) -> str:
        return "increased" if self.mean_contrast >= 0 else "decreased"


@dataclass(frozen=True)
class SignificanceOverlay:
    """Everything the renderer and the emitted table need.

    The two arms land on different things, so they are kept apart. ``marks`` is
    per cell and is what gets a border. ``column_marks`` is per body site and
    trait, which is as fine as the composition arm resolves; painting it onto
    every condition row of a block would outline 587 of the 1,284 drawn cells and
    say nothing about any one of them, so it is drawn on the column label instead.
    """

    #: Cell -> ``enrichment`` / ``composition`` / ``both``.
    marks: dict[CellKey, str]
    #: (body site, trait) -> ``increased`` / ``decreased``, pooled over conditions.
    column_marks: dict[tuple[str, str], str]
    #: Cell -> ``enriched`` / ``depleted``, from the enrichment arm alone.
    enrichment_verdicts: dict[CellKey, str]
    #: Cell -> studies behind it, and how many were called on their own.
    support: dict[CellKey, tuple[int, int]]
    composition: list[CompositionCall]
    #: Traits with no enrichment call anywhere, so nothing could be outlined.
    untested_traits: list[str]
    profiles_tested: int
    tests_performed: int

    @property
    def both_directions(self) -> int:
        """Cells whose increased and decreased columns were both called.

        Marcin's read of the first cuts was that this happens but is rare, and
        the figure reports the count so the claim is checkable.
        """
        seen: defaultdict[tuple[str, str, str], int] = defaultdict(int)
        for condition, site, trait, _direction in self.enrichment_verdicts:
            seen[(condition, site, trait)] += 1
        return sum(1 for count in seen.values() if count == len(DIRECTIONS))


def _resolve_columns(header: Sequence[str], traits: Mapping[str, Sequence[str]]) -> dict[int, str]:
    """Map profile-matrix column positions to the trait group that claims them.

    Two of the six Table 1 traits are carried as METPO-qualified columns, so the
    segment after the final pipe is matched as well as the whole column name.
    """
    wanted = {name: set(ids) for name, ids in traits.items()}
    resolved: dict[int, str] = {}
    for index, column in enumerate(header):
        term = column_term_id(column)
        for name, ids in wanted.items():
            if column in ids or term in ids:
                resolved[index] = name
    return resolved


def enrichment_calls(
    keys: Iterable[tuple[str, str, str]],
    traits: Mapping[str, Sequence[str]],
    profile_path: Path = DEFAULT_PROFILE_PATH,
) -> tuple[list[EnrichmentCall], list[str], int, int]:
    """Run the enrichment arm's hypergeometric test over ``keys``.

    ``keys`` are the ``(study, experiment, direction)`` triples the figure kept,
    so the corpus filters are whatever :func:`bsdbng.trait_guilds.load_profiles`
    applied. The background pool is built from those profiles alone: a test that
    scored a profile against a pool it is not part of would be asking a different
    question from the one the figure's cells answer.

    Every non-MediaDive trait column enters the test, not only the six being
    drawn, because Benjamini-Hochberg is applied across the whole family the arm
    tests. Restricting the family to six traits would make each of them look far
    more significant than the arm reports.
    """
    wanted = set(keys)
    with profile_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if header is None:
            msg = f"{profile_path} is empty"
            raise ValueError(msg)

        positions = [
            index
            for index, column in enumerate(header)
            if column not in PROFILE_ID_COLUMNS and MEDIADIVE_MARKER not in column.lower()
        ]
        slot_of_position = {position: slot for slot, position in enumerate(positions)}
        drawn = {
            slot_of_position[position]: name
            for position, name in _resolve_columns(header, traits).items()
            if position in slot_of_position
        }

        rows: list[tuple[tuple[str, str, str], list[int], list[int]]] = []
        background = [0] * len(positions)
        background_total = 0
        for row in reader:
            if not row:
                continue
            key = (row[0], row[1], row[2])
            if key not in wanted:
                continue
            slots: list[int] = []
            counts: list[int] = []
            for slot, position in enumerate(positions):
                raw = row[position].strip()
                if not raw or raw == "0":
                    continue
                count = int(raw)
                if count <= 0:
                    continue
                slots.append(slot)
                counts.append(count)
                background[slot] += count
                background_total += count
            if counts:
                rows.append((key, slots, counts))

    if not rows:
        return [], sorted(traits), 0, 0

    records: list[tuple[tuple[str, str, str], int, int, float, float]] = []
    p_values: list[float] = []
    for key, slots, counts in rows:
        row_total = sum(counts)
        for slot, count in zip(slots, counts, strict=True):
            successes = background[slot]
            # scipy's sf(k - 1, M, K, N) is P[X >= k], the arm's one-tailed test.
            p_value = float(hypergeom.sf(count - 1, background_total, successes, row_total))
            expected = row_total * successes / background_total
            p_values.append(p_value)
            records.append((key, slot, count, p_value, expected))

    q_values = false_discovery_control(p_values, method="bh")
    calls = [
        EnrichmentCall(
            study_id=key[0],
            experiment_id=key[1],
            inc_dec=key[2],
            trait=drawn[slot],
            p_value=p_value,
            q_value=float(q_value),
            fold=count / expected if expected > 0 else math.inf,
            count=count,
            expected=expected,
        )
        for (key, slot, count, p_value, expected), q_value in zip(records, q_values, strict=True)
        if slot in drawn
    ]
    seen = {call.trait for call in calls}
    untested = sorted(name for name in traits if name not in seen)
    return calls, untested, len(rows), len(p_values)


def _z_scored(values: Sequence[float]) -> list[float]:
    """The composition arm's transform: share per 1,000, log1p, then z-score."""
    scaled = [math.log1p(value * 1000.0) for value in values]
    mean = sum(scaled) / len(scaled)
    variance = sum((value - mean) ** 2 for value in scaled) / len(scaled)
    deviation = math.sqrt(variance)
    if deviation == 0.0:
        return [0.0] * len(scaled)
    return [(value - mean) / deviation for value in scaled]


def _contrast(values: Sequence[float]) -> tuple[float, float, float] | None:
    """One trait column's t-test against zero, or nothing if it cannot be run."""
    if len(values) < MIN_PAIRS:
        return None
    mean = sum(values) / len(values)
    spread = sum((value - mean) ** 2 for value in values)
    if spread == 0.0:
        return None
    result = ttest_1samp(values, 0.0)
    return mean, float(result.statistic), float(result.pvalue)


def composition_calls(
    profiles: Sequence[Profile],
    traits: Sequence[str],
    *,
    bin_condition: Mapping[str, str],
    bin_body_site: Mapping[str, str],
) -> list[CompositionCall]:
    """Run the composition arm's paired contrast, per body site and per cell.

    The z-scores are taken within a body site, as the arm takes them within a
    stratum, and the pairs are then partitioned by condition. Doing it the other
    way - z-scoring inside a single condition - would standardise three or four
    values against themselves and mean nothing.

    Both grainings are returned. The per-body-site call is the one the arm itself
    publishes and the one with the power to say anything; the per-condition call
    is almost always untestable, because a cell rarely has the three matched
    study pairs the t-test needs.
    """
    by_site: defaultdict[str, list[Profile]] = defaultdict(list)
    for profile in profiles:
        site = bin_body_site.get(profile.body_site)
        if site is not None:
            by_site[site].append(profile)

    calls: list[CompositionCall] = []
    site_records: list[tuple[str, str, int, float, float, float]] = []
    cell_records: list[tuple[str, str, str, int, float, float, float]] = []

    for site, site_profiles in sorted(by_site.items()):
        z_scores = {
            trait: _z_scored([profile.fraction(trait) for profile in site_profiles])
            for trait in traits
        }
        # (condition, study) -> direction -> per-trait z-scores of its profiles.
        grouped: defaultdict[tuple[str, str], defaultdict[str, list[list[float]]]] = defaultdict(
            lambda: defaultdict(list)
        )
        for index, profile in enumerate(site_profiles):
            condition = bin_condition.get(profile.condition)
            if condition is None:
                continue
            grouped[(condition, profile.study_id)][profile.inc_dec].append(
                [z_scores[trait][index] for trait in traits]
            )

        per_condition: defaultdict[str, list[list[float]]] = defaultdict(list)
        pooled: list[list[float]] = []
        for (condition, _study), directions in grouped.items():
            if not all(direction in directions for direction in DIRECTIONS):
                continue
            # One value per study before the test, per issue #133.
            increased = [
                sum(values) / len(values) for values in zip(*directions["inc"], strict=True)
            ]
            decreased = [
                sum(values) / len(values) for values in zip(*directions["dec"], strict=True)
            ]
            contrast = [up - down for up, down in zip(increased, decreased, strict=True)]
            per_condition[condition].append(contrast)
            pooled.append(contrast)

        for position, trait in enumerate(traits):
            column = [row[position] for row in pooled]
            result = _contrast(column)
            if result is not None:
                mean, statistic, p_value = result
                site_records.append((site, trait, len(column), mean, statistic, p_value))
            for condition, rows in sorted(per_condition.items()):
                cell_column = [row[position] for row in rows]
                cell_result = _contrast(cell_column)
                if cell_result is not None:
                    mean, statistic, p_value = cell_result
                    cell_records.append(
                        (site, condition, trait, len(cell_column), mean, statistic, p_value)
                    )

    # Benjamini-Hochberg within each family, as the arm corrects within one table.
    if site_records:
        q_values = false_discovery_control([row[5] for row in site_records], method="bh")
        calls.extend(
            CompositionCall(
                stratum_kind="body_site",
                body_site=site,
                condition=None,
                trait=trait,
                studies=studies,
                mean_contrast=mean,
                t_statistic=statistic,
                p_value=p_value,
                q_value=float(q_value),
            )
            for (site, trait, studies, mean, statistic, p_value), q_value in zip(
                site_records, q_values, strict=True
            )
        )
    if cell_records:
        q_values = false_discovery_control([row[6] for row in cell_records], method="bh")
        calls.extend(
            CompositionCall(
                stratum_kind="condition_body_site",
                body_site=site,
                condition=condition,
                trait=trait,
                studies=studies,
                mean_contrast=mean,
                t_statistic=statistic,
                p_value=p_value,
                q_value=float(q_value),
            )
            for (site, condition, trait, studies, mean, statistic, p_value), q_value in zip(
                cell_records, q_values, strict=True
            )
        )
    return calls


def build_overlay(
    profiles: Sequence[Profile],
    traits: Mapping[str, Sequence[str]],
    *,
    bin_condition: Mapping[str, str],
    bin_body_site: Mapping[str, str],
    profile_path: Path = DEFAULT_PROFILE_PATH,
) -> SignificanceOverlay:
    """Run both arms and fold their verdicts into one per-cell overlay."""
    trait_order = list(traits)
    keys = [(p.study_id, p.experiment_id, p.inc_dec) for p in profiles]
    calls, untested, profiles_tested, tests = enrichment_calls(keys, traits, profile_path)

    where: dict[tuple[str, str, str], tuple[str, str]] = {}
    for profile in profiles:
        condition = bin_condition.get(profile.condition)
        site = bin_body_site.get(profile.body_site)
        if condition is not None and site is not None:
            where[(profile.study_id, profile.experiment_id, profile.inc_dec)] = (condition, site)

    # Cell -> study -> the verdicts that study's profiles carried.
    per_study: defaultdict[CellKey, defaultdict[str, set[str]]] = defaultdict(
        lambda: defaultdict(set)
    )
    for call in calls:
        placed = where.get((call.study_id, call.experiment_id, call.inc_dec))
        if placed is None:
            continue
        condition, site = placed
        verdict = call.verdict
        per_study[(condition, site, call.trait, call.inc_dec)][call.study_id].add(verdict or "")

    # Every cell the figure can draw, so a cell whose studies were all tested and
    # none called is distinguishable from one that was never tested.
    totals: defaultdict[CellKey, set[str]] = defaultdict(set)
    for profile in profiles:
        condition = bin_condition.get(profile.condition)
        site = bin_body_site.get(profile.body_site)
        if condition is None or site is None:
            continue
        for trait in trait_order:
            totals[(condition, site, trait, profile.inc_dec)].add(profile.study_id)

    verdicts: dict[CellKey, str] = {}
    support: dict[CellKey, tuple[int, int]] = {}
    for cell, studies in totals.items():
        called: Mapping[str, set[str]] = per_study.get(cell, {})
        supporting = {
            study: sorted(v for v in seen if v)
            for study, seen in called.items()
            if any(v for v in seen)
        }
        support[cell] = (len(studies), len(supporting))
        if len(studies) < MIN_CELL_STUDIES:
            continue
        if len(supporting) <= CELL_SUPPORT_SHARE * len(studies):
            continue
        enriched = sum(1 for seen in supporting.values() if "enriched" in seen)
        depleted = sum(1 for seen in supporting.values() if "depleted" in seen)
        verdicts[cell] = "enriched" if enriched >= depleted else "depleted"

    composition = composition_calls(
        profiles, trait_order, bin_condition=bin_condition, bin_body_site=bin_body_site
    )
    cell_significant = {
        (call.condition, call.body_site, call.trait)
        for call in composition
        if call.stratum_kind == "condition_body_site" and call.significant and call.condition
    }

    marks: dict[CellKey, str] = {}
    for cell in totals:
        condition, site, trait, _direction = cell
        by_enrichment = cell in verdicts
        # Only the condition-resolved contrast belongs on a cell. The pooled
        # body-site contrast is a statement about the whole block and goes on the
        # column label.
        by_composition = (condition, site, trait) in cell_significant
        if by_enrichment and by_composition:
            marks[cell] = BOTH
        elif by_enrichment:
            marks[cell] = ENRICHMENT
        elif by_composition:
            marks[cell] = COMPOSITION

    return SignificanceOverlay(
        marks=marks,
        column_marks={
            (call.body_site, call.trait): call.direction
            for call in composition
            if call.stratum_kind == "body_site" and call.significant
        },
        enrichment_verdicts=verdicts,
        support=support,
        composition=composition,
        untested_traits=untested,
        profiles_tested=profiles_tested,
        tests_performed=tests,
    )


SIGNIFICANCE_HEADER = (
    "condition",
    "body_site",
    "trait",
    "direction",
    "mark",
    "enrichment_verdict",
    "studies",
    "studies_called",
)

COMPOSITION_HEADER = (
    "stratum_kind",
    "body_site",
    "condition",
    "trait",
    "studies",
    "mean_increased_minus_decreased",
    "t_statistic",
    "p_value",
    "q_value",
    "significant",
)


def write_significance_tables(
    overlay: SignificanceOverlay, cell_path: Path, composition_path: Path
) -> tuple[Path, Path]:
    """Write the per-cell overlay and the composition contrasts as wide TSVs."""
    cell_path.parent.mkdir(parents=True, exist_ok=True)
    with cell_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(SIGNIFICANCE_HEADER)
        for cell in sorted(overlay.support):
            condition, site, trait, direction = cell
            studies, called = overlay.support[cell]
            writer.writerow(
                [
                    condition,
                    site,
                    trait,
                    direction,
                    overlay.marks.get(cell, ""),
                    overlay.enrichment_verdicts.get(cell, ""),
                    studies,
                    called,
                ]
            )

    with composition_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(COMPOSITION_HEADER)
        for call in sorted(
            overlay.composition,
            key=lambda c: (c.stratum_kind, c.body_site, c.condition or "", c.trait),
        ):
            writer.writerow(
                [
                    call.stratum_kind,
                    call.body_site,
                    call.condition or "",
                    call.trait,
                    call.studies,
                    f"{call.mean_contrast:.6g}",
                    f"{call.t_statistic:.6g}",
                    f"{call.p_value:.6g}",
                    f"{call.q_value:.6g}",
                    "true" if call.significant else "false",
                ]
            )
    return cell_path, composition_path
