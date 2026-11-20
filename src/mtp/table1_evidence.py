"""Evidence marks for the Table 1 condition-by-body-site heatmaps (#172).

The heatmap shows how large a trait's annotation share is. This module attaches
three separately defined kinds of evidence to it and keeps them apart in every
table, caption and legend (#177):

* **Aggregated profile-level support.** The enrichment arm's own per-profile calls
  (``is_significant`` in the unfiltered results that ``scripts/trait_enrichment.py``
  writes with ``--all-results-output``) are read, not recomputed, and summarised
  per cell by a descriptive study-support rule (:data:`SUPPORT_RULE`). That rule
  is not a statistical test, and no cell-level p- or q-value is produced for it.
* **Body-site composition contrast.** The composition arm's study-level paired
  test, run per exact source body-site label through the arm's own selection,
  transformation and study collapse (:mod:`mtp.native_composition`).
  :data:`COMPOSITION_DEPARTURES` lists every way it differs from the arm, and the
  figure and the methods text render that constant.
* **Condition-by-body-site composition contrast.** The same test inside one exact
  condition of a tested body site, with the site's trait selection and z scale.

Inference does not depend on the figure (#176). Every test runs on exact source
labels over the corpus the maps load, and correction families are fixed before
anything is drawn. The heatmap's top-N cuts, masking and pooled Other groups only
decide where a result is shown afterwards. A study contributes one value to a
composition test however many experiments or conditions it has.

The six Table 1 traits were chosen using this same corpus, so none of these
results confirms them independently.
"""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, ttest_1samp

from mtp.count_input import profile_count_reader
from mtp.enrichment_results import load_results_tsv
from mtp.guild_maps import OTHER_LABEL, EvidenceMarks, GuildMatrix, MarkStyle
from mtp.native_composition import (
    build_pair_matrix_from_profiles,
    build_zscore_matrix,
    select_stratum_traits,
    study_contrasts,
    study_moments,
)
from mtp.paths import file_sha256
from mtp.trait_guilds import DIRECTIONS, Profile, column_term_id, retained_trait_columns

#: The composition arm publishes q-values without a cutoff. This one is chosen by
#: this integration for display, and the figure says so.
COMPOSITION_Q = 0.05
#: The composition arm's stratum floor, ``MIN_STRATIFIED_PAIRS`` in
#: ``scripts/create_trait_profile_infographic.py``, in experiment pairs.
MIN_SITE_PAIRS = 50
#: A condition inside a tested site enters the condition-level family only with at
#: least this many unique studies; ``study_moments`` cannot test fewer.
MIN_CELL_STUDIES = 2

#: The descriptive support rule. A study supports a cell when at least one of its
#: tested profiles there is called significant by the enrichment arm; a cell has
#: majority support when more than half of at least two eligible studies do.
SUPPORT_SHARE = 0.5
MIN_SUPPORT_STUDIES = 2
SUPPORT_RULE = (
    "a strict majority of the cell's eligible studies (at least 2) have at least one "
    "profile there that the enrichment arm calls significant"
)
#: The sensitivity alternative from #177: a study supports a cell only when more
#: than half of its own tested profiles there are called, so a study with many
#: experiments gets no more chances than a study with one.
ALTERNATIVE_SUPPORT_RULE = (
    "a strict majority of eligible studies (at least 2) have a strict majority of their "
    "own tested profiles in the cell called significant"
)

#: Evidence types. They are never merged into one another.
SUPPORT = "aggregated_profile_support"
SITE_COMPOSITION = "site_composition"
CELL_COMPOSITION = "condition_site_composition"

#: How each cell-level evidence type is outlined. Site results never mark a cell.
MARK_STYLES: dict[str, MarkStyle] = {
    SUPPORT: MarkStyle(dash="solid", inset=0.0),
    CELL_COMPOSITION: MarkStyle(dash="dot", inset=0.18),
}

#: Support states for a cell.
MAJORITY_SUPPORT = "majority_support"
NO_MAJORITY = "no_majority"
INSUFFICIENT_STUDIES = "insufficient_studies"

#: Composition test states. ``tested`` and the two ``study_moments`` refusals
#: come from the arm's own helper; the rest are this module's.
TESTED = "tested"
NOT_SELECTED = "trait_not_selected_in_site"
SITE_BELOW_MIN_PAIRS = "site_below_min_pairs"
CELL_BELOW_MIN_STUDIES = "condition_below_min_studies"

ARM_SCRIPT = "scripts/create_trait_profile_infographic.py"

#: Every way the composition marks differ from the composition arm as it stands on
#: main, as (aspect, the arm, these marks). The figure legend and
#: methods/table1_trait_condition_body_site_maps.md both render these rows.
COMPOSITION_DEPARTURES: tuple[tuple[str, str, str], ...] = (
    (
        "strata",
        "body_site_group from the heatmap pair metadata: the top body sites plus a pooled "
        f"Other, each kept at {MIN_SITE_PAIRS} or more pairs, plus Feces and Feces-excluded",
        f"exact source body-site labels with {MIN_SITE_PAIRS} or more eligible pairs; no pooled "
        "Other or Feces-excluded stratum, because those pools are defined by a display cut",
    ),
    (
        "zero totals",
        "zero-total profiles stay in, with every share 0, through selection, z-scoring and "
        "the test",
        "a pair with a zero retained-annotation total on either side is left out first, as "
        "the maps leave it out, because its shares are undefined",
    ),
    (
        "condition contrasts",
        "none; conditions are only summarised descriptively",
        "the site test repeated inside each exact condition with at least "
        f"{MIN_CELL_STUDIES} studies, on the site's selection and z scale, in one "
        "Benjamini-Hochberg family over every such condition, site and trait",
    ),
    (
        "cutoff",
        "q-values published without a cutoff",
        f"q <= {COMPOSITION_Q}, chosen here for display",
    ),
)

#: What the marks take from the arm unchanged, for the methods text and provenance.
COMPOSITION_SHARED = (
    "denominator: the sum of every retained non-MediaDive trait column",
    "trait selection, log1p(per-mille share) transform and z-scoring: mtp.native_composition",
    "sampling unit: each pair's increased-minus-decreased z contrast, averaged to one value "
    "per unique study in the stratum, whatever its conditions",
    "test: two-sided one-sample t-test of the study values against 0",
    "site family: Benjamini-Hochberg over every trait the stratum selects, with untestable "
    "traits entering at p = 1",
)

#: Keyed by (condition, body site, trait, direction), all exact source labels.
CellKey = tuple[str, str, str, str]
ProfileKey = tuple[str, str, str]


def _profile_id(key: ProfileKey) -> str:
    return "|".join(key)


def read_count_matrix(profile_path: Path, keys: Iterable[ProfileKey]) -> pd.DataFrame:
    """Retained (non-MediaDive) trait counts of ``keys``, validated by the shared reader.

    A blank, fractional or negative count raises rather than becoming zero (#629).
    """
    wanted = set(keys)
    reader = profile_count_reader(profile_path)
    positions, columns, _ = retained_trait_columns(reader.header)
    offset = len(reader.leading_columns)
    index: list[str] = []
    rows: list[list[int]] = []
    for row in reader:
        key = (row.key[0], row.key[1], row.key[2])
        if key not in wanted:
            continue
        index.append(_profile_id(key))
        rows.append([row.counts[position - offset] for position in positions])
    missing = len(wanted) - len(index)
    if missing:
        raise ValueError(f"{missing} corpus profiles are absent from {profile_path}")
    return pd.DataFrame(np.asarray(rows, dtype=float), index=index, columns=columns)


def resolve_trait_columns(
    columns: Iterable[str], traits: Mapping[str, Sequence[str]]
) -> dict[str, str]:
    """Map each named trait to the raw matrix column carrying it.

    Two Table 1 traits are carried behind a METPO qualifier, so the segment after
    the final pipe is matched as well as the whole column name.
    """
    resolved: dict[str, str] = {}
    for column in columns:
        term = column_term_id(column)
        for name, ids in traits.items():
            if column in ids or term in ids:
                if name in resolved:
                    raise ValueError(f"trait {name!r} matches more than one column")
                resolved[name] = column
    return resolved


def _trait_labels(columns: Sequence[str]) -> pd.DataFrame:
    """The label frame ``select_stratum_traits`` joins; selection never reads the labels."""
    return pd.DataFrame(
        {
            "trait_label": list(columns),
            "display_label": list(columns),
            "column_label": "",
            "object_name": "",
            "object_category": "",
            "display_category": "unknown",
        },
        index=list(columns),
    )


def pair_metadata(profiles: Sequence[Profile]) -> pd.DataFrame:
    """One row per experiment pair, in the arm's pair-metadata shape."""
    by_experiment: defaultdict[tuple[str, str], dict[str, Profile]] = defaultdict(dict)
    for profile in profiles:
        by_experiment[(profile.study_id, profile.experiment_id)][profile.inc_dec] = profile
    records: list[dict[str, str]] = []
    for (study, experiment), sides in sorted(by_experiment.items()):
        if set(sides) != set(DIRECTIONS):
            raise ValueError(f"unpaired experiment {(study, experiment)}; load a paired corpus")
        inc, dec = sides["inc"], sides["dec"]
        if (inc.condition, inc.body_site) != (dec.condition, dec.body_site):
            raise ValueError(f"experiment {(study, experiment)} has mismatched metadata")
        records.append(
            {
                "pair_id": f"{study}|{experiment}|paired_inc_minus_dec",
                "study_id": study,
                "experiment_id": experiment,
                "condition": inc.condition,
                "body_site": inc.body_site,
                "increased_profile_ids": _profile_id((study, experiment, "inc")),
                "decreased_profile_ids": _profile_id((study, experiment, "dec")),
            }
        )
    columns = ["pair_id", "study_id", "experiment_id", "condition", "body_site"]
    columns += ["increased_profile_ids", "decreased_profile_ids"]
    return pd.DataFrame.from_records(records, columns=columns).set_index("pair_id")


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
    mean_contrast_z: float | None = None
    t_statistic: float | None = None
    p_value: float | None = None
    family: str = ""
    family_size: int = 0
    q_value: float | None = None
    table1_family_size: int = 0
    table1_q_value: float | None = None

    @property
    def significant(self) -> bool:
        return self.q_value is not None and self.q_value <= COMPOSITION_Q

    @property
    def direction(self) -> str | None:
        if self.mean_contrast_z is None:
            return None
        return "increased" if self.mean_contrast_z >= 0 else "decreased"

    @property
    def in_family(self) -> bool:
        return self.status not in {NOT_SELECTED, SITE_BELOW_MIN_PAIRS, CELL_BELOW_MIN_STUDIES}


def _study_test(template: CompositionResult, values: np.ndarray) -> CompositionResult:
    """The arm's study-level test of one trait over one value per unique study."""
    moments = study_moments(values)
    if moments.status != TESTED:
        return replace(template, status=moments.status, mean_contrast_z=moments.mean)
    result = ttest_1samp(values, 0.0)
    return replace(
        template,
        status=TESTED,
        mean_contrast_z=moments.mean,
        t_statistic=float(result.statistic),
        p_value=float(result.pvalue),
    )


def _p_or_one(p_value: float | None) -> float:
    return 1.0 if p_value is None else p_value


def _correct(
    results: list[CompositionResult],
    family_of: Callable[[CompositionResult], str | None],
    *,
    sensitivity: bool = False,
) -> list[CompositionResult]:
    """Benjamini-Hochberg per family, untestable members at p = 1 as the arm does."""
    families: defaultdict[str, list[int]] = defaultdict(list)
    for index, result in enumerate(results):
        family = family_of(result)
        if family is not None and result.in_family:
            families[family].append(index)
    adjusted = list(results)
    for family, indices in families.items():
        p_values = [_p_or_one(results[i].p_value) for i in indices]
        for index, q_value in zip(
            indices, false_discovery_control(p_values, method="bh"), strict=True
        ):
            q = float(q_value) if results[index].p_value is not None else None
            adjusted[index] = (
                replace(adjusted[index], table1_family_size=len(indices), table1_q_value=q)
                if sensitivity
                else replace(adjusted[index], family=family, family_size=len(indices), q_value=q)
            )
    return adjusted


def _template(
    site: str,
    column: str,
    names: Mapping[str, str],
    status: str,
    rows: pd.DataFrame,
    *,
    condition: str | None = None,
) -> CompositionResult:
    """An untested result for ``column`` over ``rows``, carrying its stratum's counts."""
    return CompositionResult(
        evidence_type=SITE_COMPOSITION if condition is None else CELL_COMPOSITION,
        body_site=site,
        condition=condition,
        trait=names.get(column, column),
        column=column,
        status=status,
        studies=int(rows["study_id"].nunique()),
        paired_experiments=len(rows),
        conditions=int(rows["condition"].nunique()),
    )


def composition_contrasts(
    pairs: pd.DataFrame,
    counts: pd.DataFrame,
    table1_columns: Mapping[str, str],
    *,
    min_site_pairs: int = MIN_SITE_PAIRS,
) -> list[CompositionResult]:
    """Run the composition arm's study-level test per exact site and per exact cell.

    ``pairs`` is the analysed corpus in pair-metadata form and ``counts`` its
    retained trait counts. Nothing here reads a figure, a top-N cut or a pooled
    label: strata are the exact ``body_site`` values of ``pairs``.
    """
    names = {column: name for name, column in table1_columns.items()}
    labels = _trait_labels(list(counts.columns))
    results: list[CompositionResult] = []
    for site, group in pairs.groupby("body_site", sort=True):
        site_name = str(site)
        if len(group) < min_site_pairs:
            results.extend(
                _template(site_name, column, names, SITE_BELOW_MIN_PAIRS, group)
                for column in sorted(names)
            )
            continue
        profile_ids = [*group["increased_profile_ids"], *group["decreased_profile_ids"]]
        stratum = counts.loc[profile_ids]
        selected = select_stratum_traits(stratum, labels).index
        pair_matrix = build_pair_matrix_from_profiles(build_zscore_matrix(stratum, selected), group)
        studies = group.loc[pair_matrix.index, "study_id"]
        results.extend(
            _template(site_name, column, names, NOT_SELECTED, group)
            for column in sorted(names)
            if column not in selected
        )
        site_values = study_contrasts(pair_matrix, studies)
        by_condition = [
            (
                _template(site_name, "", names, TESTED, rows, condition=str(condition)),
                study_contrasts(pair_matrix.loc[rows.index], studies),
            )
            for condition, rows in group.groupby("condition", sort=True)
        ]
        site_base = _template(site_name, "", names, TESTED, group)
        for column in selected:
            trait = names.get(column, column)
            site_row = replace(site_base, trait=trait, column=column)
            results.append(_study_test(site_row, site_values[column].to_numpy(dtype=float)))
            for base, values in by_condition:
                cell = replace(base, trait=trait, column=column)
                if cell.studies < MIN_CELL_STUDIES:
                    if column in names:
                        results.append(replace(cell, status=CELL_BELOW_MIN_STUDIES))
                    continue
                results.append(_study_test(cell, values[column].to_numpy(dtype=float)))

    def family(result: CompositionResult) -> str:
        if result.evidence_type == SITE_COMPOSITION:
            return f"site:{result.body_site}"
        return "all_condition_site_cells"

    results = _correct(results, family)
    # The six-trait families, kept only as labelled sensitivity columns (#215).
    return _correct(
        results,
        lambda r: family(r) if r.column in names else None,
        sensitivity=True,
    )


@dataclass(frozen=True)
class EnrichmentCalls:
    """The enrichment arm's own per-profile calls for the six traits."""

    path: Path
    sha256: str
    metadata: dict[str, str]
    tests_performed: int
    # (study, experiment, direction, trait name) -> the arm's is_significant.
    called: dict[tuple[str, str, str, str], bool]


def read_enrichment_calls(path: Path, table1_columns: Mapping[str, str]) -> EnrichmentCalls:
    """Read the arm's unfiltered results and keep the six traits' calls.

    The file is validated by :func:`mtp.enrichment_results.load_results_tsv`, and
    only an upper-tail, FDR-corrected run is accepted: support is a summary of
    enrichment calls, not of nominal or depletion ones.
    """
    table = load_results_tsv(path, allow_unfiltered=True)
    if table.metadata.get("Tail") != "upper" or table.metadata.get("FDR") != "True":
        raise ValueError(f"{path} is not an upper-tail FDR-corrected enrichment run")
    if "is_significant" not in table.frame.columns:
        raise ValueError(f"{path} has no is_significant column; pass the --all-results-output file")
    flags = table.frame["is_significant"]
    unknown = sorted(set(flags) - {"True", "False"})
    if unknown:
        raise ValueError(f"{path} has is_significant values {unknown}")
    names = {column: name for name, column in table1_columns.items()}
    called: dict[tuple[str, str, str, str], bool] = {}
    frame = table.frame[table.frame["trait"].isin(names)]
    for study, experiment, direction, trait, flag in frame[
        ["study_id", "experiment_id", "inc_dec", "trait", "is_significant"]
    ].itertuples(index=False, name=None):
        called[(study, experiment, direction, names[trait])] = flag == "True"
    return EnrichmentCalls(
        path=path,
        sha256=file_sha256(path),
        metadata=table.metadata,
        tests_performed=table.tests_performed or len(table.frame),
        called=called,
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
    tested_profiles: int
    called_profiles: int
    untested_profiles: int
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
    profiles: Sequence[Profile], calls: EnrichmentCalls, traits: Sequence[str]
) -> dict[CellKey, SupportRecord]:
    """Summarise the arm's profile calls per exact (condition, body site, trait, direction).

    A profile is eligible for a trait when the arm tested it for that trait. A
    corpus profile with no result row is counted as untested, never as a miss.
    """
    # cell -> study -> [tested profiles, called profiles]
    per_study: defaultdict[CellKey, defaultdict[str, list[int]]] = defaultdict(
        lambda: defaultdict(lambda: [0, 0])
    )
    untested: defaultdict[CellKey, int] = defaultdict(int)
    for profile in profiles:
        for trait in traits:
            cell = (profile.condition, profile.body_site, trait, profile.inc_dec)
            flag = calls.called.get(
                (profile.study_id, profile.experiment_id, profile.inc_dec, trait)
            )
            if flag is None:
                untested[cell] += 1
                continue
            tally = per_study[cell][profile.study_id]
            tally[0] += 1
            tally[1] += flag
    records: dict[CellKey, SupportRecord] = {}
    for cell in sorted(set(per_study) | set(untested)):
        tallies = list(per_study[cell].values()) if cell in per_study else []
        records[cell] = SupportRecord(
            condition=cell[0],
            body_site=cell[1],
            trait=cell[2],
            direction=cell[3],
            eligible_studies=len(tallies),
            supporting_studies=sum(1 for _, hits in tallies if hits > 0),
            alternative_supporting_studies=sum(1 for n, hits in tallies if hits > n / 2),
            tested_profiles=sum(n for n, _ in tallies),
            called_profiles=sum(hits for _, hits in tallies),
            untested_profiles=untested.get(cell, 0),
            studies_with_repeated_profiles=sum(1 for n, _ in tallies if n > 1),
            max_profiles_per_study=max((n for n, _ in tallies), default=0),
        )
    return records


@dataclass(frozen=True)
class Table1Evidence:
    """All three evidence types, keyed by exact source labels, before any display."""

    traits: tuple[str, ...]
    composition: list[CompositionResult]
    support: dict[CellKey, SupportRecord] | None
    enrichment: EnrichmentCalls | None

    def table1_composition(self) -> list[CompositionResult]:
        return [r for r in self.composition if r.trait in self.traits]

    def cell_marks(self) -> dict[CellKey, frozenset[str]]:
        """Evidence types meeting their rule per cell, anywhere in the corpus.

        Support marks the direction it was aggregated over. A condition-level
        composition result is one contrast between the two directions, so it marks
        only the side it favours. Site results never mark a cell.
        """
        marks: defaultdict[CellKey, set[str]] = defaultdict(set)
        for cell, record in (self.support or {}).items():
            if record.status == MAJORITY_SUPPORT:
                marks[cell].add(SUPPORT)
        for result in self.table1_composition():
            if (
                result.evidence_type == CELL_COMPOSITION
                and result.significant
                and result.condition is not None
            ):
                side = "inc" if result.direction == "increased" else "dec"
                marks[(result.condition, result.body_site, result.trait, side)].add(
                    CELL_COMPOSITION
                )
        return {cell: frozenset(kinds) for cell, kinds in marks.items()}

    def site_marks(self) -> dict[tuple[str, str], str]:
        """(exact body site, trait) -> the direction a significant site contrast favours."""
        return {
            (result.body_site, result.trait): str(result.direction)
            for result in self.table1_composition()
            if result.evidence_type == SITE_COMPOSITION and result.significant
        }


def build_evidence(
    profiles: Sequence[Profile],
    profile_path: Path,
    traits: Mapping[str, Sequence[str]],
    *,
    enrichment_results: Path | None = None,
    min_site_pairs: int = MIN_SITE_PAIRS,
) -> Table1Evidence:
    """Compute every evidence type over ``profiles`` with no reference to a figure.

    ``profiles`` must already be the analysed population: human, paired, and
    without zero-total pairs, as the maps draw it.
    """
    keys = [(p.study_id, p.experiment_id, p.inc_dec) for p in profiles]
    counts = read_count_matrix(profile_path, keys)
    table1_columns = resolve_trait_columns(counts.columns, traits)
    absent = sorted(set(traits) - set(table1_columns))
    if absent:
        raise ValueError(f"traits with no retained column: {', '.join(absent)}")
    calls = (
        None
        if enrichment_results is None
        else read_enrichment_calls(enrichment_results, table1_columns)
    )
    if calls is not None and not calls.called:
        raise ValueError(f"{enrichment_results} holds no result for the six traits")
    return Table1Evidence(
        traits=tuple(traits),
        composition=composition_contrasts(
            pair_metadata(profiles), counts, table1_columns, min_site_pairs=min_site_pairs
        ),
        support=None if calls is None else aggregate_support(profiles, calls, list(traits)),
        enrichment=calls,
    )


def evidence_marks(evidence: Table1Evidence) -> EvidenceMarks:
    """The figure's view of the evidence, with every legend line in reader's terms.

    Each line says which evidence type it describes, the rule that produced it and
    what it does not establish, so no mark reads as a calibrated cell test.
    """
    cells = evidence.cell_marks()
    totals = dict.fromkeys(MARK_STYLES, 0)
    for kinds in cells.values():
        for kind in kinds:
            totals[kind] += 1
    legend = {
        CELL_COMPOSITION: (
            "Dotted inner outline - condition-by-body-site composition contrast: the "
            "increased-minus-decreased annotation-share contrast within that exact condition and "
            f"body site, q <= {COMPOSITION_Q} in one family over every tested condition, site "
            "and trait, drawn on the side it favours."
        ),
    }
    if evidence.enrichment is not None:
        meta = evidence.enrichment.metadata
        legend[SUPPORT] = (
            "Solid outline - aggregated profile-level support: "
            f"{SUPPORT_RULE} (one-tailed hypergeometric, FDR q <= {meta.get('P-Threshold')}, "
            f"fold >= {meta.get('Min-Fold')}, over the arm's "
            f"{evidence.enrichment.tests_performed:,} profile-trait tests). A descriptive "
            "support rule with no cell-level p- or q-value."
        )
    notes = [
        "Composition marks reuse the composition arm's study-level test with these departures: "
        + "; ".join(f"{aspect}: {ours}" for aspect, _, ours in COMPOSITION_DEPARTURES)
        + ". The six traits were selected on this corpus, so no mark confirms them.",
    ]
    if evidence.enrichment is None:
        notes.append(
            "Aggregated profile-level support is not drawn: no enrichment results were supplied."
        )
    return EvidenceMarks(
        cells=cells,
        styles={kind: MARK_STYLES[kind] for kind in legend},
        legend=legend,
        corpus_totals={kind: totals[kind] for kind in legend},
        columns=evidence.site_marks(),
        column_legend=(
            "Triangle above a trait's two columns - body-site composition contrast pooled over "
            f"every condition in that exact body site, q <= {COMPOSITION_Q} within the site's "
            "family, pointing to the side with the higher share. It tests no single condition."
        ),
        notes=tuple(notes),
    )


def _cell(value: float | int | str | None) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


SUPPORT_COLUMNS = (
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
    "tested_profiles",
    "called_profiles",
    "untested_profiles",
    "studies_with_repeated_profiles",
    "max_profiles_per_study",
    "display_status",
)
COMPOSITION_COLUMNS = (
    "evidence_type",
    "body_site",
    "condition",
    "trait",
    "column",
    "status",
    "studies",
    "paired_experiments",
    "conditions",
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
)


def write_evidence_tables(
    evidence: Table1Evidence,
    output_dir: Path,
    prefix: str,
    *,
    display_status: Callable[[CellKey], str],
    site_drawn: Callable[[str], bool],
    provenance: Mapping[str, object],
) -> list[Path]:
    """Write support, the six traits' contrasts, the full screen and provenance.

    ``display_status`` and ``site_drawn`` describe the combined figure. They are
    written beside the evidence and never feed back into it.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    if evidence.support is not None:
        support_path = output_dir / f"{prefix}_cell_support.tsv"
        with support_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(SUPPORT_COLUMNS)
            for cell, record in evidence.support.items():
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
                        record.tested_profiles,
                        record.called_profiles,
                        record.untested_profiles,
                        record.studies_with_repeated_profiles,
                        record.max_profiles_per_study,
                        display_status(cell),
                    ]
                )
        written.append(support_path)

    def row(result: CompositionResult) -> list[str]:
        if result.condition is None:
            shown = "site_drawn" if site_drawn(result.body_site) else "site_not_drawn"
        else:
            side = "inc" if result.direction != "decreased" else "dec"
            shown = display_status((result.condition, result.body_site, result.trait, side))
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
        evidence.composition,
        key=lambda r: (r.evidence_type, r.body_site, r.condition or "", r.trait),
    )
    for name, rows in (
        ("composition_contrasts", [r for r in ordered if r.trait in evidence.traits]),
        ("composition_screen", ordered),
    ):
        path = output_dir / f"{prefix}_{name}.tsv"
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(COMPOSITION_COLUMNS)
            writer.writerows(row(result) for result in rows)
        written.append(path)

    support_record: dict[str, object] | None = None
    if evidence.enrichment is not None:
        support_record = {
            "source": "scripts/trait_enrichment.py --all-results-output",
            "results_file": evidence.enrichment.path.name,
            "results_sha256": evidence.enrichment.sha256,
            "results_metadata": evidence.enrichment.metadata,
            "tests_performed": evidence.enrichment.tests_performed,
            "rule": SUPPORT_RULE,
            "alternative_rule": ALTERNATIVE_SUPPORT_RULE,
            "cell_p_or_q": None,
        }
    provenance_record = {
        **provenance,
        "evidence_types": {
            SUPPORT: support_record,
            SITE_COMPOSITION: {
                "source": ARM_SCRIPT,
                "shared_with_arm": list(COMPOSITION_SHARED),
                "family": "Benjamini-Hochberg per exact body site over the traits selected there",
                "min_site_pairs": MIN_SITE_PAIRS,
                "display_q": COMPOSITION_Q,
            },
            CELL_COMPOSITION: {
                "test": "as site_composition, within one exact condition of a tested site",
                "family": "one Benjamini-Hochberg family over every tested condition-site-trait",
                "min_studies": MIN_CELL_STUDIES,
                "display_q": COMPOSITION_Q,
            },
        },
        "composition_departures": [
            {"aspect": aspect, "composition_arm": theirs, "here": ours}
            for aspect, theirs, ours in COMPOSITION_DEPARTURES
        ],
        "excluded_from_marks": f"the pooled {OTHER_LABEL} row and column, which are not one label",
        "selection_note": "the six Table 1 traits were selected on this corpus; no result "
        "here is independent confirmation of them",
    }
    path = output_dir / f"{prefix}_evidence_provenance.json"
    path.write_text(
        json.dumps(provenance_record, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    written.append(path)
    return written


def support_summary(records: Iterable[SupportRecord]) -> dict[str, int]:
    """Counts of cells per support state under both rules, for the run report."""
    summary: defaultdict[str, int] = defaultdict(int)
    for record in records:
        summary[record.status] += 1
        summary[f"alternative_{record.alternative_status}"] += 1
        if record.status != record.alternative_status:
            summary["changed_by_alternative_rule"] += 1
    return dict(summary)


def cell_display_status(matrix: GuildMatrix, cell: CellKey) -> str:
    """Where an evidence cell, in exact source labels, lands in a drawn matrix.

    A pooled Other row or column is never a match, because it is not one label.
    """
    condition, site, trait, direction = cell
    column = f"{site}|{trait}_{direction}"
    if site == OTHER_LABEL or column not in matrix.columns:
        return "site_not_drawn"
    if condition == OTHER_LABEL or condition not in matrix.conditions:
        return "condition_not_drawn"
    row = matrix.conditions.index(condition)
    index = matrix.columns.index(column)
    if matrix.values[row][index] is not None:
        return "drawn"
    if matrix.sampled_no_signal[row][index]:
        return "no_signal"
    return "below_min_studies" if matrix.study_counts[row][index] else "not_sampled"


def report_lines(evidence: Table1Evidence, matrix: GuildMatrix) -> list[str]:
    """Each evidence type on its own lines, with drawn and corpus-wide counts."""
    lines: list[str] = []
    marks = evidence.cell_marks()
    drawn = {cell for cell in marks if cell_display_status(matrix, cell) == "drawn"}
    if evidence.support is None or evidence.enrichment is None:
        lines.append("Aggregated profile support: not computed (no enrichment results).")
    else:
        summary = support_summary(evidence.support.values())
        lines.append(
            f"Aggregated profile support over {len(evidence.support)} exact cells, from "
            f"{evidence.enrichment.tests_performed} profile-trait tests: {summary}"
        )
        supported = {cell for cell, kinds in marks.items() if SUPPORT in kinds}
        lines.append(
            f"  cells meeting the support rule: {len(supported)}; drawn in the combined "
            f"figure: {len(supported & drawn)}"
        )
        both: defaultdict[tuple[str, str, str], int] = defaultdict(int)
        for condition, site, trait, _ in supported & drawn:
            both[(condition, site, trait)] += 1
        lines.append(
            "  drawn condition-site-trait combinations supported in both directions: "
            f"{sum(1 for n in both.values() if n == 2)} of {len(both)}"
        )
    for result in sorted(evidence.table1_composition(), key=lambda r: _p_or_one(r.q_value)):
        if result.evidence_type == SITE_COMPOSITION and result.significant:
            lines.append(
                f"  site composition {result.body_site} / {result.trait}: mean "
                f"{result.mean_contrast_z:+.3f} z over {result.studies} studies "
                f"({result.paired_experiments} pairs, {result.conditions} conditions), "
                f"q {result.q_value:.3g} in {result.family} (n={result.family_size}), "
                f"six-trait q {_p_or_one(result.table1_q_value):.3g}"
            )
    for evidence_type in (SITE_COMPOSITION, CELL_COMPOSITION):
        rows = [r for r in evidence.composition if r.evidence_type == evidence_type]
        table1 = [r for r in rows if r.trait in evidence.traits]
        lines.append(
            f"{evidence_type}: {sum(r.p_value is not None for r in rows)} tested, "
            f"{sum(r.significant for r in rows)} at q <= {COMPOSITION_Q}; six traits: "
            f"{sum(r.p_value is not None for r in table1)} tested, "
            f"{sum(r.significant for r in table1)} significant"
        )
    cell_marks = {cell for cell, kinds in marks.items() if CELL_COMPOSITION in kinds}
    lines.append(
        f"  condition-site composition marks: {len(cell_marks)}; drawn in the combined "
        f"figure: {len(cell_marks & drawn)}"
    )
    return lines
