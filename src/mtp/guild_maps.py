"""Disease-by-body-site oxygen guild heatmap and aerobe/anaerobe correlation.

Rows are conditions, columns are body sites crossed with the four
guild-by-direction combinations, so each body site contributes the 2x2
`{aerobe, anaerobe} x {inc, dec}` table unrolled into four adjacent columns.

Cell values are mean guild trait fractions. Because 55.4% of the figure's 922
studies contribute more than one experiment (issue #133; measured on the
committed corpus), profiles are averaged within a study before being averaged
across studies, and the companion n-matrix counts studies rather than experiments.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass, field, replace
from pathlib import Path

from scipy.cluster.hierarchy import fcluster, leaves_list, linkage
from scipy.spatial.distance import pdist
from scipy.stats import binomtest, spearmanr

from mtp.corpus_identity import corpus_identity, summary_lines, write_identity
from mtp.paths import DERIVED_DATA_DIR
from mtp.trait_guilds import (
    DEFAULT_EXPERIMENTS_PATH,
    DEFAULT_FULL_DUMP_PATH,
    DEFAULT_PROFILE_PATH,
    DIRECTIONS,
    OXYGEN_GUILDS,
    PairExclusion,
    Profile,
    drop_zero_total_pairs,
    format_pair_exclusion,
    load_profiles,
    read_bugsigdb_vintage,
    write_stats_report,
)

OTHER_LABEL = "Other"
GUILD_ORDER = ("aerobe", "anaerobe")

#: The oxygen guild figure's aggregation parameters, which are its command-line
#: defaults. The integrated evidence table inventories the figure's cells through
#: ``build_figure_matrix`` at these values. It used to keep its own copy, which
#: stayed at 8 body sites after the figure moved to 20 (#770).
FIGURE_MAX_CONDITIONS = 60
FIGURE_MAX_BODY_SITES = 20
FIGURE_MIN_STUDIES = 3
FIGURE_ROW_CLUSTERS = 6


@dataclass(frozen=True)
class PooledLabels:
    """How much of the corpus a top-N cut left out of its own row or column."""

    kind: str
    kept: int
    pooled_labels: int
    pooled_profiles: int
    total_profiles: int
    # Labels actually kept, which is below ``kept`` when the corpus has fewer
    # distinct labels than the cap allows.
    kept_labels: int = 0
    # Whether the remainder went into an Other bucket or was dropped outright.
    bucketed: bool = True
    # Exact labels outside the top-N cut, even when the Other row is later masked.
    pooled_names: tuple[str, ...] = ()

    @property
    def pooled_share(self) -> float:
        """Percentage of profiles that fell outside the kept labels."""
        return 100.0 * self.pooled_profiles / self.total_profiles if self.total_profiles else 0.0

    def describe(self) -> str:
        verb = "Other pools" if self.bucketed else "dropped"
        return (
            f"{self.kind}: kept top {self.kept}; {verb} {self.pooled_labels} labels "
            f"and {self.pooled_profiles} profiles ({self.pooled_share:.1f}%)"
        )


@dataclass(frozen=True)
class GuildMatrix:
    """Mean guild fractions and supporting study counts per condition/body site."""

    conditions: list[str]
    body_sites: list[str]
    columns: list[str]
    values: list[list[float | None]]
    study_counts: list[list[int]]
    # Cluster id per condition row, from Euclidean clustering of the row vectors.
    row_clusters: list[int]
    pooled: list[PooledLabels]
    # Conditions with no visible cells under the study/no-signal display policy,
    # excluding those made empty only by the body-site cuts.
    dropped_empty_conditions: list[str]
    # Body sites whose column block was backed by too few conditions to read as a
    # site-level pattern, removed with all of their columns.
    dropped_body_sites: list[str]
    # True where the body site was sampled for that condition but no taxon in the
    # profile carried any trait of that guild. Without this, an all-zero cell is
    # drawn in the same colour as a genuinely low fraction.
    sampled_no_signal: list[list[bool]]
    # Keep the build's visibility policy when cutting a focus panel. The masks
    # remain populated even when their shading is disabled.
    mark_no_signal: bool = True
    # Conditions with drawable cells before body-site exclusion, but none after.
    dropped_site_conditions: list[str] = field(default_factory=list)


def _top_labels(counts: dict[str, int], limit: int) -> list[str]:
    """Return the most frequent labels, ranked deterministically."""
    ranked = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    return [label for label, _ in ranked[:limit]]


def _binned(label: str, keep: set[str], *, bucket_rest: bool = True) -> str | None:
    """Map a label to itself, to the Other bucket, or out of the figure entirely."""
    if not label:
        return None
    if label in keep:
        return label
    return OTHER_LABEL if bucket_rest else None


def _pooling(
    counts: dict[str, int], keep: set[str], limit: int, kind: str, *, bucketed: bool = True
) -> PooledLabels:
    pooled = {k: v for k, v in counts.items() if k not in keep}
    return PooledLabels(
        kind=kind,
        kept=limit,
        pooled_labels=len(pooled),
        pooled_profiles=sum(pooled.values()),
        pooled_names=tuple(sorted(pooled)),
        total_profiles=sum(counts.values()),
        kept_labels=len(keep),
        bucketed=bucketed,
    )


def _cluster_rows(values: list[list[float | None]], n_clusters: int) -> tuple[list[int], list[int]]:
    """Order condition rows by Euclidean average-linkage clustering.

    Masked cells are treated as zero for the distance calculation only. That means
    the clustering responds to which body sites a condition was sampled in as well
    as to the guild fractions themselves, which is the intended behaviour here:
    conditions studied in the same anatomical sites should group together.
    """
    if len(values) < 3:
        return list(range(len(values))), [1] * len(values)

    matrix = [[0.0 if v is None else v for v in row] for row in values]
    distances = pdist(matrix, metric="euclidean")
    if not distances.any():
        return list(range(len(values))), [1] * len(values)

    # Ward on the Euclidean distances. Average linkage chains here, leaving one
    # 52-condition blob and five near-singletons; Ward gives usable groups.
    linkage_matrix = linkage(distances, method="ward")
    order = [int(i) for i in leaves_list(linkage_matrix)]
    n = max(1, min(n_clusters, len(values)))
    labels = [int(x) for x in fcluster(linkage_matrix, t=n, criterion="maxclust")]
    return order, labels


def build_matrix(
    profiles: list[Profile],
    *,
    max_conditions: int,
    max_body_sites: int,
    min_studies: int,
    n_clusters: int = 6,
    guild_order: Sequence[str] = GUILD_ORDER,
    include_other_conditions: bool = True,
    include_other_body_sites: bool = True,
    mark_no_signal: bool = True,
    min_site_conditions: int = 1,
) -> GuildMatrix:
    """Aggregate profiles into the condition x (group x direction) grid per body site.

    ``guild_order`` names the trait groups to unroll into columns. It is the two
    oxygen guilds by default and the six Table 1-named traits for the
    matched-corpus analog figure; the aggregation, masking, and clustering are
    identical either way.

    The two ``include_other_*`` flags decide, per axis, what happens to the labels
    the top-N cuts leave out: pooled into an Other row or column block, or dropped
    from the figure. Pooling keeps the corpus total in view but the pooled row or
    block mixes unrelated labels, so no cell in it reads as one condition or one
    anatomical site; dropping makes every drawn cell a single label at the cost of
    showing only part of the corpus. The axes are separate because the trade lands
    differently on each. ``mark_no_signal`` controls whether a
    condition kept solely because some site was sampled without any trait of the
    group is retained; those cells are shaded dark grey by the figure, so with the
    shading off such a row would render blank and is dropped instead.
    ``min_site_conditions`` drops a body site's whole column block when fewer than
    that many conditions are drawn in it: a block carrying one condition is a
    single row, not a site-level pattern, and it still costs one column per group
    and direction.
    """
    # Columns are encoded as "site|group_direction", so a body site containing the
    # separator would be split in the wrong place and silently mislabel its cells.
    unsafe = sorted({p.body_site for p in profiles if "|" in p.body_site})
    if unsafe:
        raise ValueError(f"body-site labels contain '|': {', '.join(unsafe)}")

    condition_counts: defaultdict[str, int] = defaultdict(int)
    body_site_counts: defaultdict[str, int] = defaultdict(int)
    for profile in profiles:
        if profile.condition:
            condition_counts[profile.condition] += 1
        if profile.body_site:
            body_site_counts[profile.body_site] += 1

    keep_conditions = set(_top_labels(dict(condition_counts), max_conditions))
    keep_body_sites = set(_top_labels(dict(body_site_counts), max_body_sites))
    pooled = [
        _pooling(
            dict(condition_counts),
            keep_conditions,
            max_conditions,
            "condition",
            bucketed=include_other_conditions,
        ),
        _pooling(
            dict(body_site_counts),
            keep_body_sites,
            max_body_sites,
            "body_site",
            bucketed=include_other_body_sites,
        ),
    ]

    # study -> cell -> list of per-profile fractions, so a study is averaged once.
    per_study: defaultdict[tuple[str, str, str], defaultdict[str, list[float]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for profile in profiles:
        condition = _binned(
            profile.condition, keep_conditions, bucket_rest=include_other_conditions
        )
        body_site = _binned(
            profile.body_site, keep_body_sites, bucket_rest=include_other_body_sites
        )
        if condition is None or not profile.body_site:
            continue
        # Preserve evidence at sites excluded by the top-N cut, without pooling
        # different excluded sites. It determines why a condition has no row.
        if body_site is None:
            body_site = profile.body_site
        for guild in guild_order:
            key = (condition, body_site, f"{guild}_{profile.inc_dec}")
            per_study[key][profile.study_id].append(profile.fraction(guild))

    other_condition = [OTHER_LABEL] if include_other_conditions and condition_counts else []
    other_site = [OTHER_LABEL] if include_other_body_sites and body_site_counts else []
    conditions = sorted(keep_conditions) + other_condition
    body_sites = sorted(keep_body_sites) + other_site
    columns = [
        f"{site}|{guild}_{direction}"
        for site in body_sites
        for guild in guild_order
        for direction in DIRECTIONS
    ]

    def summarize_cell(by_study: dict[str, list[float]]) -> tuple[float | None, int, bool]:
        """Apply the same study support and signal rules at kept and excluded sites."""
        n_studies = len(by_study)
        study_means = [sum(v) / len(v) for v in by_study.values()]
        zero_signal = n_studies > 0 and not any(m > 0 for m in study_means)
        value = (
            sum(study_means) / n_studies
            if n_studies > 0 and n_studies >= min_studies and not zero_signal
            else None
        )
        return value, n_studies, zero_signal

    visible_before_site_cuts: set[str] = set()
    for (condition, site, _), excluded_studies in per_study.items():
        if site not in body_sites:
            value, _, zero_signal = summarize_cell(excluded_studies)
            if value is not None or (mark_no_signal and zero_signal):
                visible_before_site_cuts.add(condition)

    values: list[list[float | None]] = []
    study_counts: list[list[int]] = []
    sampled_no_signal: list[list[bool]] = []
    for condition in conditions:
        value_row: list[float | None] = []
        count_row: list[int] = []
        blank_row: list[bool] = []
        for site in body_sites:
            for guild in guild_order:
                for direction in DIRECTIONS:
                    by_study: dict[str, list[float]] = per_study.get(
                        (condition, site, f"{guild}_{direction}"), {}
                    )
                    value, n_studies, zero_signal = summarize_cell(by_study)
                    count_row.append(n_studies)
                    blank_row.append(zero_signal)
                    value_row.append(value)
        values.append(value_row)
        study_counts.append(count_row)
        sampled_no_signal.append(blank_row)

    def drawn(row_index: int, column_indices: Sequence[int]) -> bool:
        """Whether a row shows anything at all within the given columns."""
        return any(values[row_index][i] is not None for i in column_indices) or (
            mark_no_signal and any(sampled_no_signal[row_index][i] for i in column_indices)
        )

    # Drop body-site blocks that too few conditions reach. Done before the rows so
    # a condition drawn only inside a dropped block goes with it.
    visible_before_site_cuts.update(
        condition for i, condition in enumerate(conditions) if drawn(i, range(len(columns)))
    )
    site_columns = {
        site: [i for i, column in enumerate(columns) if _site_of_column(column) == site]
        for site in body_sites
    }
    dropped_sites = [
        site
        for site, indices in site_columns.items()
        if sum(drawn(row, indices) for row in range(len(conditions))) < min_site_conditions
    ]
    if dropped_sites:
        keep_columns = [
            i
            for i, column in enumerate(columns)
            if _site_of_column(column) not in set(dropped_sites)
        ]
        body_sites = [site for site in body_sites if site not in set(dropped_sites)]
        columns = [columns[i] for i in keep_columns]
        values = [[row[i] for i in keep_columns] for row in values]
        study_counts = [[row[i] for i in keep_columns] for row in study_counts]
        sampled_no_signal = [[row[i] for i in keep_columns] for row in sampled_no_signal]

    # Drop rows that would render entirely blank. While the figure shades them, a
    # row keeping only dark-grey cells is retained: "sampled, no guild trait" is a
    # result, not an absence. With the shading off it is an invisible row.
    all_columns = list(range(len(columns)))
    keep_rows = [i for i in range(len(conditions)) if drawn(i, all_columns)]
    dropped = [conditions[i] for i in range(len(conditions)) if i not in set(keep_rows)]
    dropped_with_sites = [name for name in dropped if name in visible_before_site_cuts]
    dropped_blank = [name for name in dropped if name not in visible_before_site_cuts]
    conditions = [conditions[i] for i in keep_rows]
    values = [values[i] for i in keep_rows]
    study_counts = [study_counts[i] for i in keep_rows]
    sampled_no_signal = [sampled_no_signal[i] for i in keep_rows]

    order, labels = _cluster_rows(values, n_clusters)
    return GuildMatrix(
        conditions=[conditions[i] for i in order],
        body_sites=body_sites,
        columns=columns,
        values=[values[i] for i in order],
        study_counts=[study_counts[i] for i in order],
        row_clusters=[labels[i] for i in order],
        pooled=pooled,
        dropped_empty_conditions=dropped_blank,
        dropped_body_sites=dropped_sites,
        sampled_no_signal=[sampled_no_signal[i] for i in order],
        mark_no_signal=mark_no_signal,
        dropped_site_conditions=dropped_with_sites,
    )


def build_figure_matrix(
    profiles: Sequence[Profile],
    *,
    max_conditions: int = FIGURE_MAX_CONDITIONS,
    max_body_sites: int = FIGURE_MAX_BODY_SITES,
    min_studies: int = FIGURE_MIN_STUDIES,
    n_clusters: int = FIGURE_ROW_CLUSTERS,
    guild_order: Sequence[str] = GUILD_ORDER,
    include_other_conditions: bool = True,
    include_other_body_sites: bool = True,
    mark_no_signal: bool = True,
    min_site_conditions: int = 1,
) -> GuildMatrix:
    """The matrix the oxygen guild figure draws, at the figure's parameters by default.

    ``main`` draws the figure and writes its cell-status table from this call, and
    the integrated evidence table builds its cell inventory from it too, so the two
    cannot describe different grids (#770). The figure's population is the shared
    comparison basis: an experiment pair whose retained-annotation total is zero
    on a side is excluded here, whatever the caller passes. Its share is
    undefined, not zero, and would otherwise be averaged into a cell as a measured
    0 or shade one as sampled without signal (#783).

    The Table 1-trait maps build their matrix here too, passing their own traits
    and display options, so both condition-by-body-site figures share one
    population (#811).
    """
    return build_matrix(
        drop_zero_total_pairs(profiles).profiles,
        max_conditions=max_conditions,
        max_body_sites=max_body_sites,
        min_studies=min_studies,
        n_clusters=n_clusters,
        guild_order=guild_order,
        include_other_conditions=include_other_conditions,
        include_other_body_sites=include_other_body_sites,
        mark_no_signal=mark_no_signal,
        min_site_conditions=min_site_conditions,
    )


#: Contrasts smaller than this in magnitude are floating-point residue of
#: summing fractions in different orders, not a real direction.
SIGN_TOLERANCE = 1e-12


@dataclass(frozen=True)
class CorrelationResult:
    """Guild-by-direction association and cross-cell rank correlations."""

    table: dict[str, int]
    #: ``None`` when the discordant ratio is undefined, never ``inf``: with no
    #: study shifting the guilds the other way, 0/0 has no value to report.
    odds_ratio: float | None
    #: ``None`` when no study is discordant, so no McNemar test exists. Writing
    #: 1.0 there would read as "tested, no difference" rather than "not tested".
    mcnemar_p: float | None
    discordant: tuple[int, int]
    spearman: dict[str, tuple[float, float, int]]
    studies: int
    #: ``tested``, ``all_tied`` (studies scored, none discordant) or
    #: ``no_studies``. The same vocabulary the Table 1 arm uses.
    status: str
    #: Studies held out of the rank correlations for carrying no guild
    #: annotation on either side. The McNemar already drops them as ties.
    excluded_no_annotation: int = 0


def require_plotly() -> None:
    """Fail with the install command when the figure extra is absent (#288).

    plotly and kaleido are a ``viz`` extra rather than runtime dependencies, so
    an ingest-only install does not pull a headless-Chrome driver. Every figure
    function imports plotly lazily and calls this first, so the missing extra
    reads as a setup step rather than as a stack trace.
    """
    if importlib.util.find_spec("plotly") is not None:
        return
    msg = (
        "Rendering a figure needs the 'viz' extra (plotly and kaleido). Install it "
        "with `uv sync --all-extras`, or `uv sync --extra viz`. PNG export also "
        "needs a headless Chrome: `just get-chrome`."
    )
    raise ModuleNotFoundError(msg)


#: Per-cell render states of an aggregated matrix. Once the value is blank, a
#: cell shaded for carrying no annotation, a cell held back for too little
#: support, and a cell nothing was sampled in all look identical, and they are
#: three different claims. The figure writer and the integrated evidence table
#: both classify cells here, so no cell can be drawn in one artifact and masked
#: in the other (#294).
CELL_VALUE = "value"
CELL_NO_SIGNAL = "no_guild_signal"
CELL_BELOW_MIN_STUDIES = "below_min_studies"
CELL_NOT_SAMPLED = "not_sampled"


def cell_status(value: float | None, sampled_no_signal: bool, study_count: int) -> str:
    """Why a cell carries a value, or which of three reasons it does not."""
    if value is not None:
        return CELL_VALUE
    if sampled_no_signal:
        return CELL_NO_SIGNAL
    return CELL_BELOW_MIN_STUDIES if study_count else CELL_NOT_SAMPLED


def cell_status_rows(matrix: GuildMatrix) -> list[list[str]]:
    """``cell_status`` over the whole matrix, in its own row and column order."""
    return [
        [
            cell_status(value, blank, count)
            for value, blank, count in zip(v_row, b_row, c_row, strict=True)
        ]
        for v_row, b_row, c_row in zip(
            matrix.values, matrix.sampled_no_signal, matrix.study_counts, strict=True
        )
    ]


def _odds_ratio_label(odds_ratio: float | None) -> str:
    """Render the discordant ratio, saying so when it is undefined."""
    return "undefined (no reverse-discordant study)" if odds_ratio is None else f"{odds_ratio:.3f}"


def study_direction_means(
    profiles: list[Profile],
    guild_order: Sequence[str] = GUILD_ORDER,
    *,
    require_paired: bool = True,
) -> dict[str, dict[str, float]]:
    """Mean guild fraction per study and direction, collapsing repeated experiments.

    Pooling all increased and all decreased profiles of a study equals the mean of
    per-experiment paired differences only when every experiment contributes both
    sides, so by default the per-study inc and dec profile counts must match;
    callers on unpaired corpora must pass require_paired=False and must not
    describe the result as paired.
    """
    sums: defaultdict[str, defaultdict[str, float]] = defaultdict(lambda: defaultdict(float))
    counts: defaultdict[str, defaultdict[str, int]] = defaultdict(lambda: defaultdict(int))
    for profile in profiles:
        for guild in guild_order:
            key = f"{guild}_{profile.inc_dec}"
            sums[profile.study_id][key] += profile.fraction(guild)
            counts[profile.study_id][key] += 1
    if require_paired:
        for study, per_key in counts.items():
            for guild in guild_order:
                inc, dec = per_key.get(f"{guild}_inc", 0), per_key.get(f"{guild}_dec", 0)
                if inc != dec:
                    raise ValueError(
                        f"study {study} has {inc} increased vs {dec} decreased profiles; "
                        "study means are only paired contrasts on a paired corpus "
                        "(load_profiles(paired_only=True), the default)"
                    )
    return {
        study: {key: total / counts[study][key] for key, total in per_key.items()}
        for study, per_key in sums.items()
    }


def _sign(value: float) -> int:
    if abs(value) < SIGN_TOLERANCE:
        return 0
    return 1 if value > 0 else -1


def compute_correlation(
    profiles: list[Profile], _matrix: GuildMatrix | None = None
) -> CorrelationResult:
    """Test whether aerobes shift toward the increased side more often than anaerobes.

    Scoring counts studies, not experiments or trait assertions. A study is
    scored once per guild by the sign of (mean increased fraction - mean decreased
    fraction). Because the SAME studies provide both the aerobe and the anaerobe
    sign, the two guilds form matched pairs, not independent samples: the test is
    an exact McNemar (a binomial test on the discordant studies — aerobe up while
    anaerobe down, versus the reverse), and the odds ratio is the discordant
    ratio. An earlier version ran Fisher's exact on the marginal 2x2, which
    assumes independent rows and was anticonservative here by roughly three
    orders of magnitude. The marginal table is still reported for display.
    """
    study_means = study_direction_means(profiles)

    table = {f"{guild}_{d}": 0 for guild in GUILD_ORDER for d in DIRECTIONS}
    signs: dict[str, dict[str, int]] = defaultdict(dict)
    for study, per_key in study_means.items():
        for guild in GUILD_ORDER:
            inc = per_key.get(f"{guild}_inc")
            dec = per_key.get(f"{guild}_dec")
            if inc is None or dec is None:
                continue
            sign = _sign(inc - dec)
            if sign == 0:
                continue
            signs[study][guild] = sign
            table[f"{guild}_inc" if sign > 0 else f"{guild}_dec"] += 1

    aerobe_up_anaerobe_down = anaerobe_up_aerobe_down = 0
    for per_guild in signs.values():
        aerobe, anaerobe = per_guild.get("aerobe"), per_guild.get("anaerobe")
        if aerobe is None or anaerobe is None or aerobe == anaerobe:
            continue
        if aerobe > 0:
            aerobe_up_anaerobe_down += 1
        else:
            anaerobe_up_aerobe_down += 1

    # An undefined result is reported as undefined, not as a number that reads
    # like one (#289). This is the fault the Table 1 arm already fixed; the guild
    # arm carried it until now. Neither branch fires on the pinned corpus.
    n_discordant = aerobe_up_anaerobe_down + anaerobe_up_aerobe_down
    mcnemar_p = (
        float(binomtest(aerobe_up_anaerobe_down, n_discordant, p=0.5).pvalue)
        if n_discordant
        else None
    )
    odds_ratio = (
        aerobe_up_anaerobe_down / anaerobe_up_aerobe_down if anaerobe_up_aerobe_down else None
    )
    if not signs:
        status = "no_studies"
    elif not n_discordant:
        status = "all_tied"
    else:
        status = "tested"

    # Rank correlations across studies, one observation per study (#165). An
    # earlier version correlated condition x body-site cells, which counted a
    # study once per cell it reached and so treated dependent cells as independent.
    # Each pair uses the studies carrying both series; a constant series has no
    # defined rank correlation and its pair is omitted rather than written as nan.
    #
    # A study whose every guild fraction is zero carries no oxygen annotation on
    # either side. That zero is missingness, not a measured share, so it is
    # dropped rather than ranked (#289). The McNemar above already drops these as
    # ties; excluding them here makes both statistics describe one population.
    # Keeping them as identical (0, 0) points attenuates the within-side
    # correlations. On the committed corpus all 5 of 932 were zero-total pairs,
    # which the figure, the integrated table and calibration now exclude before
    # calling this (#783, #730). The rule still covers a study annotated without
    # any oxygen term.
    direction_keys = [f"{guild}_{d}" for guild in GUILD_ORDER for d in DIRECTIONS]
    per_study: dict[str, dict[str, float]] = {}
    excluded_no_annotation = 0
    for study, per_key in study_means.items():
        if all(per_key.get(key, 0.0) == 0.0 for key in direction_keys):
            excluded_no_annotation += 1
            continue
        row = dict(per_key)
        for guild in GUILD_ORDER:
            inc, dec = per_key.get(f"{guild}_inc"), per_key.get(f"{guild}_dec")
            if inc is not None and dec is not None:
                row[f"{guild}_contrast"] = inc - dec
        per_study[study] = row

    directions = [f"{guild}_{d}" for guild in GUILD_ORDER for d in DIRECTIONS]
    pairs = [(left, right) for i, left in enumerate(directions) for right in directions[i + 1 :]]
    # The Table 1 template: the two guilds' per-study contrasts against each other.
    pairs.append(("aerobe_contrast", "anaerobe_contrast"))
    spearman: dict[str, tuple[float, float, int]] = {}
    for left, right in pairs:
        paired = [
            (row[left], row[right]) for row in per_study.values() if left in row and right in row
        ]
        if len(paired) < 3:
            continue
        left_values = [x for x, _ in paired]
        right_values = [y for _, y in paired]
        if len(set(left_values)) < 2 or len(set(right_values)) < 2:
            continue
        rho, p_value = spearmanr(left_values, right_values)
        spearman[f"{left}~{right}"] = (float(rho), float(p_value), len(paired))

    return CorrelationResult(
        table=table,
        odds_ratio=odds_ratio,
        mcnemar_p=mcnemar_p,
        discordant=(aerobe_up_anaerobe_down, anaerobe_up_aerobe_down),
        spearman=spearman,
        studies=len(study_means),
        status=status,
        excluded_no_annotation=excluded_no_annotation,
    )


def _write_matrix_tsv(path: Path, matrix: GuildMatrix, rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["condition", *matrix.columns])
        for condition, row in zip(matrix.conditions, rows, strict=True):
            writer.writerow([condition, *row])


def write_matrix_outputs(
    matrix: GuildMatrix, output_dir: Path, prefix: str = "guild"
) -> list[Path]:
    """Write the value matrix, the study-count matrix, cell status, and row clusters.

    Shared by the oxygen guild figure and the matched-corpus Table 1-trait
    analog figure so the two read from one implementation; only the file-name
    prefix differs.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    value_rows = [["" if v is None else f"{v:.6f}" for v in row] for row in matrix.values]
    count_rows = [[str(c) for c in row] for row in matrix.study_counts]

    status_rows = cell_status_rows(matrix)

    values_path = output_dir / f"{prefix}_condition_body_site_matrix.tsv"
    status_path = output_dir / f"{prefix}_condition_body_site_cell_status.tsv"
    clusters_path = output_dir / f"{prefix}_condition_row_clusters.tsv"
    counts_path = output_dir / f"{prefix}_condition_body_site_study_counts.tsv"

    _write_matrix_tsv(values_path, matrix, value_rows)
    _write_matrix_tsv(counts_path, matrix, count_rows)
    _write_matrix_tsv(status_path, matrix, status_rows)

    with clusters_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["row_order", "condition", "cluster"])
        for order, (condition, cluster) in enumerate(
            zip(matrix.conditions, matrix.row_clusters, strict=True), start=1
        ):
            writer.writerow([order, condition, cluster])

    return [values_path, counts_path, status_path, clusters_path]


def write_outputs(
    matrix: GuildMatrix,
    correlation: CorrelationResult,
    output_dir: Path,
    vintage: str | None,
    exclusion: PairExclusion | None = None,
) -> list[Path]:
    """Write the value matrix, the study-count matrix, and the statistics table.

    ``exclusion`` is the zero-total pair exclusion the statistics ran after. It is
    written as the table's first row, so the table states its own population
    (#783).
    """
    written = write_matrix_outputs(matrix, output_dir)
    stats_path = output_dir / "guild_aerobe_anaerobe_correlation.tsv"

    with stats_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["statistic", "comparison", "value", "p_value", "n", "bugsigdb_vintage"])
        if exclusion is not None:
            writer.writerow(
                [
                    "excluded_experiment_pairs",
                    "zero_retained_annotation_total_on_a_side",
                    exclusion.zero_total_experiments,
                    "",
                    exclusion.studies_after,
                    vintage or "",
                ]
            )
        for key, count in sorted(correlation.table.items()):
            writer.writerow(["studies_shifted", key, count, "", correlation.studies, vintage or ""])
        up_down, down_up = correlation.discordant
        writer.writerow(
            [
                "discordant_studies",
                "aerobe_up_anaerobe_down",
                up_down,
                "",
                up_down + down_up,
                vintage or "",
            ]
        )
        writer.writerow(
            [
                "discordant_studies",
                "anaerobe_up_aerobe_down",
                down_up,
                "",
                up_down + down_up,
                vintage or "",
            ]
        )
        writer.writerow(
            [
                "mcnemar_odds_ratio",
                f"aerobe_vs_anaerobe_paired_signs ({correlation.status})",
                "" if correlation.odds_ratio is None else f"{correlation.odds_ratio:.6f}",
                "" if correlation.mcnemar_p is None else f"{correlation.mcnemar_p:.6g}",
                up_down + down_up,
                vintage or "",
            ]
        )
        for pair, (rho, p_value, n_studies) in sorted(correlation.spearman.items()):
            writer.writerow(
                [
                    "spearman_rho_across_studies",
                    pair,
                    f"{rho:.6f}",
                    f"{p_value:.6g}",
                    n_studies,
                    vintage or "",
                ]
            )

    return [*written, stats_path]


GUILD_COLORSCALES = {"aerobe": "OrRd", "anaerobe": "YlGn"}
# Colour-bar ticks as percentages with trailing zeros trimmed, so 0.0004 reads
# "0.04%" rather than plotly's default SI-prefixed "400u".
PERCENT_TICKS = ".3~%"
# Room for a panel figure's rotated column labels below it. The room above it is
# sized from its subtitle, as the combined figure's is.
PANEL_BOTTOM_MARGIN = 260
# A panel's colour bars lie flat under the figure, three across and two deep, in
# the band this reserves below the column labels. Stacking them vertically ties
# the legend's height to the number of conditions - a four-row cut would either
# pad itself out with empty rows or collapse its bars - and measuring them against
# the figure rather than the plotting area breaks that tie.
# Bars are laid out this far apart, and the band is sized from however many rows
# of them a figure's width leaves room for.
LEGEND_ROW_PITCH = 140
# A plotting area only a few rows tall makes the renderer fail outright
# ("axis scaling"), so a cut shorter than this keeps a little slack below it.
MIN_PANEL_ROWS = 8
# Every heatmap cell is this tall, in every figure, so a cell in a four-row cut
# is the same size as a cell in the combined figure.
ROW_HEIGHT = 26
NO_SIGNAL_COLOR = "#4d4d4d"
CLUSTER_GAP_ROWS = 1
# The combined figure's title line with the gaps around it, and one subtitle line,
# in pixels at plotly's default fonts. The top margin is sized from them.
TITLE_BLOCK_HEIGHT = 80
SUBTITLE_LINE_HEIGHT = 22
# The body-site names above the plotting area: a combined figure's block headers,
# at most two lines of 18-pixel bold text (``_wrapped``) lifted 12 pixels, or a
# panel figure's one-line subplot titles.
SITE_HEADER_HEIGHT = 60


def _guild_of_column(column: str) -> str:
    """Return the guild encoded in a `site|guild_direction` column label."""
    return column.split("|", 1)[1].rsplit("_", 1)[0]


def _site_of_column(column: str) -> str:
    """Return the body site encoded in a `site|guild_direction` column label."""
    return column.split("|", 1)[0]


def _wrapped(label: str, width: int = 16) -> str:
    """Break a long body-site name at the space nearest its middle.

    A header wider than its column block runs into the next one; the names that
    overflow are all two or three words, so one break is enough.
    """
    if len(label) <= width or " " not in label:
        return label
    middle = len(label) / 2
    split = min((i for i, char in enumerate(label) if char == " "), key=lambda i: abs(i - middle))
    return f"{label[:split]}<br>{label[split + 1 :]}"


def _site_blocks(columns: Sequence[str]) -> list[tuple[str, int, int]]:
    """Contiguous `(site, first index, last index)` runs over the column labels."""
    blocks: list[tuple[str, int, int]] = []
    for index, column in enumerate(columns):
        site = _site_of_column(column)
        if blocks and blocks[-1][0] == site:
            blocks[-1] = (site, blocks[-1][1], index)
        else:
            blocks.append((site, index, index))
    return blocks


def _upper_range(values: list[float], quantile: float) -> float:
    """Top of a colour range: the quantile of observed values, or their maximum.

    With one heavily sampled body site and a long tail of lightly sampled ones, a
    single extreme cell sets the maximum and flattens every other cell to nearly
    the same pale tint. Cutting the ramp at a high quantile restores contrast in
    the body of the distribution; cells above the cut saturate at the top colour
    and the cut is stated on the figure.
    """
    if not values:
        return 1.0
    if quantile >= 1.0:
        return max(values)
    ordered = sorted(values)
    index = min(len(ordered) - 1, max(0, round(quantile * (len(ordered) - 1))))
    return ordered[index] or max(ordered)


def _insert_cluster_gaps(
    matrix: GuildMatrix,
) -> tuple[list[str], list[list[float | None]], list[list[bool]]]:
    """Add blank spacer rows between clusters so the cut tree is visible.

    The spacers exist only in the rendered figure; the emitted tables keep one row
    per condition and carry the cluster id in a separate column.
    """
    width = len(matrix.columns)
    labels: list[str] = []
    values: list[list[float | None]] = []
    blanks: list[list[bool]] = []
    previous: int | None = None
    spacer = 0
    for index, condition in enumerate(matrix.conditions):
        cluster = matrix.row_clusters[index]
        if previous is not None and cluster != previous:
            for _ in range(CLUSTER_GAP_ROWS):
                spacer += 1
                # Figure-space padding keeps every category label unique and invisible.
                labels.append("\u2007" * spacer)
                values.append([None] * width)
                blanks.append([False] * width)
        labels.append(condition)
        values.append(matrix.values[index])
        blanks.append(matrix.sampled_no_signal[index])
        previous = cluster
    return labels, values, blanks


def _coverage_note(matrix: GuildMatrix) -> list[str]:
    """State how much of the corpus the top-N cuts left off the figure.

    Without this the reader has no way to tell from the figure that the rows are
    a heavily sampled minority of BugSigDB's conditions rather than all of them.
    """
    pooling = {p.kind: p for p in matrix.pooled}
    condition = pooling.get("condition")
    site = pooling.get("body_site")
    if condition is None or site is None:
        return []

    total_conditions = condition.kept_labels + condition.pooled_labels
    total_sites = site.kept_labels + site.pooled_labels
    drawn = [c for c in matrix.conditions if c != OTHER_LABEL]
    drawn_sites = [b for b in matrix.body_sites if b != OTHER_LABEL]
    condition_part = f"{len(drawn)} of {total_conditions} conditions"
    if condition.bucketed:
        condition_part += " plus a pooled Other row"
    site_part = f"{len(drawn_sites)} of {total_sites} body sites"
    if site.bucketed:
        site_part += " plus a pooled Other column block"

    def fate(pool: PooledLabels) -> str:
        verb = "pooled into Other" if pool.bucketed else "not drawn"
        return f"{verb} ({pool.pooled_share:.0f}% of profiles)"

    lines = [
        f"Rows and columns are a cut of the corpus: {condition_part} and {site_part}.",
        f"The remaining {condition.pooled_labels} conditions are {fate(condition)} and the "
        f"remaining {site.pooled_labels} body sites are {fate(site)}.",
    ]
    if condition.bucketed or site.bucketed:
        lines.append(
            "A pooled Other mixes unrelated labels, so no cell in it is one condition or one "
            "anatomical site."
        )
    return lines


def pair_exclusion_note(exclusion: PairExclusion) -> str:
    """The subtitle line that states a figure's population (#783).

    The oxygen guild figure and the Table 1-trait maps both draw from
    ``build_figure_matrix``, so both state the exclusion in these words (#811).
    """
    return (
        f"{exclusion.zero_total_experiments} experiment pairs with a zero "
        "retained-annotation total on a side are excluded (studies "
        f"{exclusion.studies_before} to {exclusion.studies_after}): their share is "
        "undefined, not zero."
    )


def write_figure(
    matrix: GuildMatrix,
    output_path: Path,
    vintage: str | None,
    *,
    guild_order: Sequence[str] = GUILD_ORDER,
    colorscales: dict[str, object] | None = None,
    title: str = "Oxygen guilds by condition and body site",
    group_noun: str = "guild",
    value_phrase: str = "guild trait fraction",
    shared_range: bool = True,
    colorbar_labels: dict[str, str] | None = None,
    range_quantile: float = 1.0,
    separate_body_sites: bool = False,
    mark_no_signal: bool = True,
    note: str | None = None,
) -> Path:
    """Render the heatmap with one colour scale per trait group.

    ``note`` is one more subtitle line, which the oxygen guild figure uses to state
    its population (#783).

    Each group is drawn as its own overlaid trace so it keeps its own hue. With
    ``shared_range`` the traces share the range 0 to the largest observed value,
    which makes a cell in one group directly comparable in magnitude to a cell in
    another; that only holds where the groups are commensurable, as the two oxygen
    guilds are. The Table 1-named traits differ in corpus-wide prevalence by more
    than an order of magnitude, so they are drawn on per-trait ranges instead and each
    colour bar carries its own maximum.

    Cells where the site was sampled but carried no trait of the group are filled
    dark grey, so an absent trait is distinguishable from a small fraction. With
    ``mark_no_signal`` off they are left blank instead: the shading dominates the
    figure where the traits are rare, and the state is still recorded in the
    emitted cell-status table.
    """
    require_plotly()
    import plotly.graph_objects as go

    scales: dict[str, object] = dict(GUILD_COLORSCALES) if colorscales is None else colorscales
    labels = colorbar_labels or {}
    row_labels, row_values, row_blanks = _insert_cluster_gaps(matrix)
    observed = [v for row in row_values for v in row if v is not None]
    z_max = _upper_range(observed, range_quantile)

    per_group: dict[str, list[list[float | None]]] = {}
    for guild in guild_order:
        per_group[guild] = [
            [
                value if _guild_of_column(column) == guild else None
                for column, value in zip(matrix.columns, row, strict=True)
            ]
            for row in row_values
        ]

    figure = go.Figure()
    count = len(guild_order)
    # Even slices of the plotting area, one per group. Sizing the bars and then
    # spacing them apart needs more than the whole height once there are six of
    # them, which silently clips the top and bottom bars.
    bar_slice = 1.0 / count
    bar_len = min(0.42, bar_slice * 0.8)
    for position, guild in enumerate(guild_order):
        z_guild = per_group[guild]
        group_values = [v for row in z_guild for v in row if v is not None]
        scale_max = z_max if shared_range else _upper_range(group_values, range_quantile)
        # Evenly space the colour bars down the right margin, top group first.
        centre = 1.0 - bar_slice * (position + 0.5)
        figure.add_trace(
            go.Heatmap(
                z=z_guild,
                x=matrix.columns,
                y=row_labels,
                colorscale=scales[guild],
                zmin=0.0,
                zmax=scale_max,
                hoverongaps=False,
                name=guild,
                colorbar={
                    "title": f"{labels.get(guild, guild)}<br>% of trait assertions",
                    "len": bar_len,
                    "y": centre,
                    "yanchor": "middle",
                    "thickness": 14,
                    # Percentages, not plotly's default SI prefixes: a tick reading
                    # "400u" for 0.0004 looks like a unit rather than a share.
                    "tickformat": PERCENT_TICKS,
                    # Six bars share the height; without a cap their ticks collide.
                    "nticks": 3,
                },
            )
        )

    shapes = [
        {
            "type": "rect",
            "xref": "x",
            "yref": "y",
            "x0": col - 0.5,
            "x1": col + 0.5,
            "y0": row - 0.5,
            "y1": row + 0.5,
            "fillcolor": NO_SIGNAL_COLOR,
            "line": {"color": NO_SIGNAL_COLOR, "width": 0},
            "layer": "above",
        }
        for row, flags in enumerate(row_blanks)
        for col, sampled in enumerate(flags)
        if sampled and mark_no_signal
    ]

    annotations: list[dict[str, object]] = []
    if separate_body_sites:
        # One rule between body-site blocks. With six traits times two directions
        # per site the blocks otherwise run together and the x labels are the only
        # thing telling the reader where one body site ends.
        boundaries = [
            index
            for index in range(1, len(matrix.columns))
            if _site_of_column(matrix.columns[index]) != _site_of_column(matrix.columns[index - 1])
        ]
        shapes.extend(
            {
                "type": "line",
                "xref": "x",
                "yref": "paper",
                "x0": index - 0.5,
                "x1": index - 0.5,
                "y0": 0,
                "y1": 1,
                "line": {"color": "#666666", "width": 1},
                "layer": "above",
            }
            for index in boundaries
        )
        # A header per block. The rotated per-column labels are unreadable at a
        # glance and repeat the site name twelve times; the anatomical site is the
        # first thing a reader needs, so it is named once, large, above the block.
        annotations.extend(
            {
                "x": (start + end) / 2.0,
                "y": 1.0,
                "xref": "x",
                "yref": "paper",
                "yshift": 12,
                "text": f"<b>{_wrapped(site)}</b>",
                "showarrow": False,
                "font": {"size": 18},
                "xanchor": "center",
                "yanchor": "bottom",
            }
            for site, start, end in _site_blocks(matrix.columns)
        )

    n_clusters = len(set(matrix.row_clusters))
    range_note = (
        f"All scales span 0 to {z_max:.3f}."
        if shared_range
        else f"Each {group_noun} is on its own scale; read magnitudes within a {group_noun}, "
        "not across them."
    )
    if range_quantile < 1.0:
        range_note += (
            f" Ramps are cut at the {range_quantile:.0%} quantile of observed cells; "
            "cells above it saturate."
        )
    cluster_note = (
        f"Rows clustered by Euclidean distance into {n_clusters} groups, separated by blank rows."
    )
    if mark_no_signal:
        cluster_note += (
            f" Dark grey = site sampled for that condition but no taxon carried the {group_noun}."
        )
    # Broken into explicit lines: as one run the subtitle overflows the plotting
    # area and is clipped by the right-hand colour bars.
    lines = [
        f"Mean {value_phrase}, averaged within study then across studies, shown as a "
        "percentage of the profile's summed trait-assertion counts (taxonomic "
        "rollups can count one reported taxon many times), not of its taxa.",
        range_note,
        cluster_note,
        *_coverage_note(matrix),
        *([note] if note else []),
        f"BugSigDB {vintage or 'vintage unknown'}.",
    ]
    title_layout: dict[str, object] = {"text": f"{title}<br><sub>{'<br>'.join(lines)}</sub>"}
    # The title block hangs from the top of the figure, and the plotting area starts
    # below its last line and below the site headers when they sit above it. A fixed
    # 220-pixel margin with the title centred in it let the last subtitle lines run
    # over the top rows (#783). A fixed 420 pixels with site headers let the
    # Table 1-trait maps' population line run into a two-line header (#811).
    top = TITLE_BLOCK_HEIGHT + SUBTITLE_LINE_HEIGHT * len(lines)
    if separate_body_sites:
        top += SITE_HEADER_HEIGHT
    height = ROW_HEIGHT * len(row_labels) + 300 + top
    title_layout.update({"yref": "container", "yanchor": "top", "y": 1 - 20 / height})
    figure.update_layout(
        title=title_layout,
        xaxis={"title": f"body site | {group_noun} _ direction", "tickangle": -60},
        yaxis={"title": "condition", "autorange": "reversed"},
        shapes=shapes,
        annotations=annotations,
        width=max(900, 26 * len(matrix.columns)),
        height=height,
        margin={"l": 300, "b": 300, "t": top},
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.write_image(str(output_path), scale=2)
    return output_path


def subset_matrix(
    matrix: GuildMatrix,
    *,
    conditions: Sequence[str] | None = None,
    sites: Sequence[str] | None = None,
) -> GuildMatrix:
    """Cut the matrix down to a set of conditions and body sites.

    One large heatmap answers "what was measured where"; it does not read as a
    result, because the block that carries most of the corpus sets the eye's
    scale and the sparse blocks beside it are mostly empty. Decomposing it into
    a few small heatmaps, each holding conditions that were sampled in the same
    places, makes each one legible on its own.

    Either axis may be left as ``None`` to be derived from the other: with only
    ``conditions`` given the body sites become the ones those conditions were
    actually drawn in, and with only ``sites`` given the conditions become the
    ones drawn in those sites. Rows and column blocks that end up carrying no
    visible cells are dropped, using the build's ``mark_no_signal`` policy so
    shaded no-signal cells survive the cut. A named condition absent from the
    matrix is skipped and recorded in
    ``dropped_empty_conditions``; the row order and cluster ids are the ones the
    full matrix was clustered into, so a subset stays comparable to it.
    """
    if conditions is None and sites is None:
        return matrix

    requested_sites = list(matrix.body_sites) if sites is None else list(sites)
    unknown = [site for site in requested_sites if site not in matrix.body_sites]
    if unknown:
        raise ValueError(f"body sites absent from the matrix: {', '.join(unknown)}")

    if conditions is None:
        row_indices = list(range(len(matrix.conditions)))
        missing: list[str] = []
    else:
        wanted = set(conditions)
        row_indices = [i for i, name in enumerate(matrix.conditions) if name in wanted]
        missing = [name for name in conditions if name not in matrix.conditions]

    column_indices = [
        i
        for site in requested_sites
        for i, column in enumerate(matrix.columns)
        if _site_of_column(column) == site
    ]

    def is_drawn(row: int, cols: Sequence[int]) -> bool:
        return any(matrix.values[row][i] is not None for i in cols) or (
            matrix.mark_no_signal and any(matrix.sampled_no_signal[row][i] for i in cols)
        )

    row_indices = [i for i in row_indices if is_drawn(i, column_indices)]

    if sites is None:
        # Derive the blocks from the rows that survived: a site nothing in this
        # subset was drawn in would be an empty band of columns.
        def block(site: str) -> list[int]:
            return [i for i in column_indices if _site_of_column(matrix.columns[i]) == site]

        requested_sites = [
            site
            for site in requested_sites
            if any(is_drawn(row, block(site)) for row in row_indices)
        ]
        column_indices = [
            i
            for site in requested_sites
            for i, column in enumerate(matrix.columns)
            if _site_of_column(column) == site
        ]

    dropped_sites = [site for site in matrix.body_sites if site not in set(requested_sites)]
    kept = {matrix.conditions[i] for i in row_indices}
    dropped_conditions = missing + [name for name in matrix.conditions if name not in kept]

    # A cut that leaves the Other row or column behind must stop claiming to pool
    # into it, or its subtitle describes the figure it was cut from.
    kept_labels = {
        "condition": [matrix.conditions[i] for i in row_indices],
        "body_site": requested_sites,
    }
    pooled = [
        replace(pool, bucketed=pool.bucketed and OTHER_LABEL in kept_labels[pool.kind])
        if pool.kind in kept_labels
        else pool
        for pool in matrix.pooled
    ]

    return replace(
        matrix,
        pooled=pooled,
        conditions=[matrix.conditions[i] for i in row_indices],
        body_sites=requested_sites,
        columns=[matrix.columns[i] for i in column_indices],
        values=[[matrix.values[r][c] for c in column_indices] for r in row_indices],
        study_counts=[[matrix.study_counts[r][c] for c in column_indices] for r in row_indices],
        row_clusters=[matrix.row_clusters[i] for i in row_indices],
        sampled_no_signal=[
            [matrix.sampled_no_signal[r][c] for c in column_indices] for r in row_indices
        ],
        dropped_empty_conditions=dropped_conditions,
        dropped_body_sites=dropped_sites,
    )


def write_body_site_panels(
    matrix: GuildMatrix,
    output_path: Path,
    vintage: str | None,
    *,
    sites: Sequence[str],
    guild_order: Sequence[str] = GUILD_ORDER,
    colorscales: dict[str, object] | None = None,
    title: str = "Body site detail",
    group_noun: str = "guild",
    value_phrase: str = "guild trait fraction",
    colorbar_labels: dict[str, str] | None = None,
    range_quantile: float = 1.0,
    mark_no_signal: bool = True,
    scale_reference: GuildMatrix | None = None,
    note: str | None = None,
) -> Path:
    """Draw one panel per named body site, at a width where the labels are readable.

    The combined figure fits every site by making each column narrow enough that
    the rotated `site | trait _ direction` labels run together. This pulls a few
    blocks out into their own panels: same rows in the same clustered order, same
    per-trait colour ranges as the combined figure so a cell here is the same
    colour there, but with the site named once above the panel and the columns
    labelled by trait and direction alone.
    """
    require_plotly()
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    scales: dict[str, object] = dict(GUILD_COLORSCALES) if colorscales is None else colorscales
    labels = colorbar_labels or {}
    row_labels, row_values, row_blanks = _insert_cluster_gaps(matrix)
    # Slack under a very short cut is added as blank rows, not by stretching the
    # cells: every cell is ROW_HEIGHT tall in every figure. A distinct invisible
    # character keeps each label unique without colliding with the cluster spacers.
    blank_width = len(matrix.columns)
    for index in range(MIN_PANEL_ROWS - len(row_labels)):
        row_labels.append("\u2008" * (index + 1))
        row_values.append([None] * blank_width)
        row_blanks.append([False] * blank_width)

    missing = [site for site in sites if site not in matrix.body_sites]
    if missing:
        raise ValueError(f"body sites absent from the matrix: {', '.join(missing)}")

    # Colour ranges come from ``scale_reference`` when one is given - the full
    # matrix - rather than from the panels, so a cell keeps the same colour in
    # every figure cut from it and the small heatmaps stay comparable.
    reference = scale_reference or matrix
    scale_max = {
        guild: _upper_range(
            [
                value
                for row in reference.values
                for column, value in zip(reference.columns, row, strict=True)
                if value is not None and _guild_of_column(column) == guild
            ],
            range_quantile,
        )
        for guild in guild_order
    }

    figure = make_subplots(
        rows=1,
        cols=len(sites),
        shared_yaxes=True,
        horizontal_spacing=0.06,
        subplot_titles=[f"<b>{site}</b>" for site in sites],
    )

    range_note = f"Each {group_noun} is on its own scale, the same one the combined figure uses"
    if range_quantile < 1.0:
        range_note += f", cut at the {range_quantile:.0%} quantile of observed cells"
    grey_note = (
        f" Dark grey = site sampled for that condition but no taxon carried the {group_noun}."
        if mark_no_signal
        else ""
    )
    lines = [
        f"Mean {value_phrase}, averaged within study then across studies, shown as a "
        "percentage of the profile's summed trait-assertion counts (taxonomic "
        "rollups can count one reported taxon many times), not of its taxa.",
        f"{range_note}.",
        "Rows are the conditions of the combined figure in the same clustered order, so a "
        f"blank row here was drawn in another body site.{grey_note}",
        *_coverage_note(matrix),
        *([note] if note else []),
        f"BugSigDB {vintage or 'vintage unknown'}.",
    ]
    # Sized from the subtitle and the site titles under it, as in ``write_figure``.
    # A fixed 400 pixels with the title centred in it let a panel's last subtitle
    # lines run into its site titles and top row (#811).
    top = TITLE_BLOCK_HEIGHT + SUBTITLE_LINE_HEIGHT * len(lines) + SITE_HEADER_HEIGHT

    count = len(guild_order)
    # All six bars on one line. Plotly places a horizontal bar's title a fixed
    # distance below it that does not track the row spacing set here, so a second
    # row of bars lands under the first row's titles; one row cannot collide with
    # itself. The width floor is what a six-bar legend and the subtitle need.
    width = max(1600, 300 + 480 * len(sites))
    legend_columns = count
    legend_depth = -(-count // legend_columns)
    legend_height = 50 + legend_depth * LEGEND_ROW_PITCH
    bottom = PANEL_BOTTOM_MARGIN + legend_height
    height = top + bottom + ROW_HEIGHT * len(row_labels)
    # Bar rows measured up from the foot of the figure, top row first.
    legend_rows = [
        (legend_height - 80.0 - row * LEGEND_ROW_PITCH) / height for row in range(legend_depth)
    ]
    shapes: list[dict[str, object]] = []
    for panel, site in enumerate(sites, start=1):
        indices = [i for i, column in enumerate(matrix.columns) if _site_of_column(column) == site]
        # The site is named once above the panel, so a column only has to carry its
        # trait and direction.
        tick_text = []
        for i in indices:
            guild = _guild_of_column(matrix.columns[i])
            direction = matrix.columns[i].rsplit("_", 1)[1]
            tick_text.append(f"{labels.get(guild, guild)} {direction}")
        for position, guild in enumerate(guild_order):
            z_guild = [
                [row[i] if _guild_of_column(matrix.columns[i]) == guild else None for i in indices]
                for row in row_values
            ]
            column, legend_row = position % legend_columns, position // legend_columns
            figure.add_trace(
                go.Heatmap(
                    z=z_guild,
                    x=tick_text,
                    y=row_labels,
                    colorscale=scales[guild],
                    zmin=0.0,
                    zmax=scale_max[guild],
                    hoverongaps=False,
                    name=guild,
                    # One colour bar per trait for the whole figure, not per panel.
                    showscale=panel == len(sites),
                    colorbar={
                        # The unit is in the subtitle: three bars across a narrow
                        # figure leave no room to repeat it under each one.
                        "title": {"text": labels.get(guild, guild), "side": "bottom"},
                        "orientation": "h",
                        "x": 0.07 + column * (0.88 / legend_columns),
                        "xanchor": "left",
                        "xref": "container",
                        "y": legend_rows[legend_row],
                        "yanchor": "bottom",
                        "yref": "container",
                        "len": 0.78 / legend_columns,
                        "lenmode": "fraction",
                        "thickness": 12,
                        "tickformat": PERCENT_TICKS,
                        # Three ticks per bar; more collide at this width.
                        "nticks": 3,
                    },
                ),
                row=1,
                col=panel,
            )
        if mark_no_signal:
            axis = "" if panel == 1 else str(panel)
            shapes.extend(
                {
                    "type": "rect",
                    "xref": f"x{axis}",
                    "yref": f"y{axis}",
                    "x0": col - 0.5,
                    "x1": col + 0.5,
                    "y0": row - 0.5,
                    "y1": row + 0.5,
                    "fillcolor": NO_SIGNAL_COLOR,
                    "line": {"color": NO_SIGNAL_COLOR, "width": 0},
                    "layer": "above",
                }
                for row, flags in enumerate(row_blanks)
                for col, i in enumerate(indices)
                if flags[i]
            )

    figure.update_layout(
        title={
            "text": f"{title}<br><sub>{'<br>'.join(lines)}</sub>",
            "yref": "container",
            "yanchor": "top",
            "y": 1 - 20 / height,
        },
        shapes=shapes,
        width=width,
        # The margins are text - the subtitle and site titles above, rotated column
        # labels below - so the plotting area has to be sized on top of them, not
        # carved out of a fixed total. A four-row cut would otherwise render as a
        # sliver with the colour bars stacked on each other.
        height=height,
        margin={"l": 300, "b": bottom, "t": top},
        plot_bgcolor="white",
    )
    figure.update_xaxes(tickangle=-60)
    figure.update_yaxes(autorange="reversed")
    figure.update_yaxes(title_text="condition", row=1, col=1)
    for annotation in figure.layout.annotations[: len(sites)]:
        annotation.font = {"size": 18}

    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.write_image(str(output_path), scale=2)
    return output_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Generate a condition x body-site oxygen guild heatmap with the 2x2 "
            "guild/direction table unrolled into four columns per body site."
        )
    )
    parser.add_argument("--core-run", type=Path, help="Verified core workflow run manifest")
    parser.add_argument("--profile-input", type=Path)
    parser.add_argument("--experiments-input", type=Path)
    parser.add_argument("--full-dump-input", type=Path)
    parser.add_argument("--output-dir", type=Path, default=DERIVED_DATA_DIR / "guild_maps")
    parser.add_argument("--max-conditions", type=int, default=FIGURE_MAX_CONDITIONS)
    parser.add_argument("--max-body-sites", type=int, default=FIGURE_MAX_BODY_SITES)
    parser.add_argument(
        "--row-clusters",
        type=int,
        default=FIGURE_ROW_CLUSTERS,
        help="Number of condition clusters to cut the Euclidean tree into.",
    )
    parser.add_argument(
        "--min-studies",
        type=int,
        default=FIGURE_MIN_STUDIES,
        help="Cells backed by fewer studies are masked rather than drawn.",
    )
    parser.add_argument(
        "--include-non-human",
        action="store_true",
        help="Keep non-human host species (issue #129 excludes them by default).",
    )
    parser.add_argument("--code-commit", help="Source revision, supplied by the Justfile recipe")
    return parser


def main() -> int:
    from mtp.workflow.consumers import (
        bind_count_arguments,
        core_corpus_identity,
        preflight_consumer_outputs,
        write_consumer_record,
    )

    args = build_parser().parse_args()
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
    profiles, stats = load_profiles(
        args.profile_input,
        args.experiments_input,
        OXYGEN_GUILDS,
        human_only=not args.include_non_human,
    )
    vintage = read_bugsigdb_vintage(args.full_dump_input)
    write_stats_report(stats, vintage)
    identity = (
        core_corpus_identity(bundle, bugsigdb_vintage=vintage)
        if bundle is not None
        else corpus_identity(args.profile_input, bugsigdb_vintage=vintage)
    )
    print("\n".join(summary_lines(identity)))

    # One population for the figure, its cell table and its McNemar: the shared
    # comparison basis, which the integrated evidence table's tests use too (#783).
    exclusion = drop_zero_total_pairs(profiles)
    print(format_pair_exclusion(exclusion))
    profiles = exclusion.profiles

    matrix = build_figure_matrix(
        profiles,
        max_conditions=args.max_conditions,
        max_body_sites=args.max_body_sites,
        min_studies=args.min_studies,
        n_clusters=args.row_clusters,
    )
    for pooled in matrix.pooled:
        print(pooled.describe())
    if matrix.dropped_body_sites:
        print(f"body sites dropped as entirely blank: {', '.join(matrix.dropped_body_sites)}")
    if matrix.dropped_empty_conditions:
        print(
            f"conditions dropped as entirely blank: {len(matrix.dropped_empty_conditions)} "
            f"({', '.join(sorted(matrix.dropped_empty_conditions))})"
        )
    if matrix.dropped_site_conditions:
        print(
            "conditions dropped with excluded body-site blocks: "
            f"{len(matrix.dropped_site_conditions)} "
            f"({', '.join(sorted(matrix.dropped_site_conditions))})"
        )
    correlation = compute_correlation(profiles)
    written = write_outputs(matrix, correlation, args.output_dir, vintage, exclusion=exclusion)
    written.append(
        write_figure(
            matrix,
            args.output_dir / "guild_condition_body_site_heatmap.png",
            vintage,
            note=pair_exclusion_note(exclusion),
        )
    )

    up_down, down_up = correlation.discordant
    ratio = "undefined" if correlation.odds_ratio is None else f"{correlation.odds_ratio:.3f}"
    p_text = "not tested" if correlation.mcnemar_p is None else f"p={correlation.mcnemar_p:.3g}"
    print(
        f"\nMcNemar (exact, paired signs): discordant {up_down}/{down_up}, "
        f"odds ratio {ratio} ({p_text}) over {correlation.studies} studies "
        f"[{correlation.status}]"
    )
    if correlation.excluded_no_annotation:
        print(
            f"rank correlations exclude {correlation.excluded_no_annotation} studies "
            "carrying no guild annotation on either side"
        )
    written.append(write_identity(identity, args.output_dir, code_commit=args.code_commit))
    if bundle is not None:
        written.append(
            write_consumer_record(
                bundle,
                args.output_dir / "core_consumption.json",
                consumer="mtp.guild_maps",
                outputs=written,
            )
        )
    for path in written:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
