"""Condition-by-body-site oxygen annotation guild heatmap and correlation.

Rows are conditions, columns are body sites crossed with the four
guild-by-direction combinations, so each body site contributes the 2x2
`{aerobe, anaerobe} x {inc, dec}` table unrolled into four adjacent columns.

Cell values are mean guild trait fractions. Because 56.2% of studies contribute
more than one experiment (issue #133; measured on the pinned v1.3.1 corpus),
profiles are averaged within a study before being averaged across studies, and
the companion n-matrix counts studies rather than experiments.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from scipy.cluster.hierarchy import fcluster, leaves_list, linkage
from scipy.spatial.distance import pdist
from scipy.stats import binomtest, spearmanr

from bsdbng.paths import DERIVED_DATA_DIR
from bsdbng.trait_guilds import (
    DEFAULT_EXPERIMENTS_PATH,
    DEFAULT_FULL_DUMP_PATH,
    DEFAULT_PROFILE_PATH,
    DIRECTIONS,
    OXYGEN_GUILDS,
    Profile,
    load_profiles,
    read_bugsigdb_vintage,
    write_stats_report,
)

OTHER_LABEL = "Other"
GUILD_ORDER = ("aerobe", "anaerobe")


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
    # Conditions removed because no cell survived the study threshold and none was
    # sampled-without-signal, so the row would render entirely blank.
    dropped_empty_conditions: list[str]
    # Body sites whose column block was backed by too few conditions to read as a
    # site-level pattern, removed with all of their columns.
    dropped_body_sites: list[str]
    # True where the body site was sampled but profiles have no recorded guild
    # annotation signal. Without this, an all-zero cell is
    # drawn in the same colour as a genuinely low fraction.
    sampled_no_signal: list[list[bool]]


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
    oxygen guilds by default and the six Table 1 traits for the Table 1 figure;
    the aggregation, masking, and clustering are identical either way.

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
        if condition is None or body_site is None:
            continue
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
                    n_studies = len(by_study)
                    count_row.append(n_studies)
                    study_means = [sum(v) / len(v) for v in by_study.values()]
                    zero_signal = n_studies > 0 and not any(m > 0 for m in study_means)
                    blank_row.append(zero_signal)
                    if n_studies < min_studies or zero_signal:
                        value_row.append(None)
                        continue
                    value_row.append(sum(study_means) / len(study_means))
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
        dropped_empty_conditions=dropped,
        dropped_body_sites=dropped_sites,
        sampled_no_signal=[sampled_no_signal[i] for i in order],
    )


#: Contrasts smaller than this in magnitude are floating-point residue of
#: summing fractions in different orders, not a real direction.
SIGN_TOLERANCE = 1e-12


@dataclass(frozen=True)
class CorrelationResult:
    """Guild-by-direction association and cross-cell rank correlations."""

    table: dict[str, int]
    odds_ratio: float
    mcnemar_p: float
    discordant: tuple[int, int]
    spearman: dict[str, tuple[float, float, int]]
    studies: int


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


def compute_correlation(profiles: list[Profile], matrix: GuildMatrix) -> CorrelationResult:
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

    n_discordant = aerobe_up_anaerobe_down + anaerobe_up_aerobe_down
    mcnemar_p = (
        float(binomtest(aerobe_up_anaerobe_down, n_discordant, p=0.5).pvalue)
        if n_discordant
        else 1.0
    )
    odds_ratio = (
        aerobe_up_anaerobe_down / anaerobe_up_aerobe_down
        if anaerobe_up_aerobe_down
        else float("inf")
    )

    # Rank correlation across condition x body-site cells. Only cells where both
    # members of a pair survived masking are used, so the series stay aligned.
    per_cell: defaultdict[tuple[str, str], dict[str, float]] = defaultdict(dict)
    for row_index, row in enumerate(matrix.values):
        condition = matrix.conditions[row_index]
        for index, column in enumerate(matrix.columns):
            site, suffix = column.split("|", 1)
            value = row[index]
            if value is not None:
                per_cell[(condition, site)][suffix] = value

    suffixes = [f"{guild}_{d}" for guild in GUILD_ORDER for d in DIRECTIONS]
    spearman: dict[str, tuple[float, float, int]] = {}
    for i, left in enumerate(suffixes):
        for right in suffixes[i + 1 :]:
            paired = [
                (cell[left], cell[right])
                for cell in per_cell.values()
                if left in cell and right in cell
            ]
            if len(paired) < 3:
                continue
            rho, p_value = spearmanr([x for x, _ in paired], [y for _, y in paired])
            spearman[f"{left}~{right}"] = (float(rho), float(p_value), len(paired))

    return CorrelationResult(
        table=table,
        odds_ratio=float(odds_ratio),
        mcnemar_p=mcnemar_p,
        discordant=(aerobe_up_anaerobe_down, anaerobe_up_aerobe_down),
        spearman=spearman,
        studies=len(study_means),
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

    Shared by the oxygen guild figure and the Table 1 trait figure so the two
    read from one implementation; only the file-name prefix differs.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    value_rows = [["" if v is None else f"{v:.6f}" for v in row] for row in matrix.values]
    count_rows = [[str(c) for c in row] for row in matrix.study_counts]

    status_rows = [
        [
            "value"
            if value is not None
            else "no_guild_signal"
            if blank
            else ("below_min_studies" if count else "not_sampled")
            for value, blank, count in zip(v_row, b_row, c_row, strict=True)
        ]
        for v_row, b_row, c_row in zip(
            matrix.values, matrix.sampled_no_signal, matrix.study_counts, strict=True
        )
    ]

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
) -> list[Path]:
    """Write the value matrix, the study-count matrix, and the statistics table."""
    written = write_matrix_outputs(matrix, output_dir)
    stats_path = output_dir / "guild_aerobe_anaerobe_correlation.tsv"

    with stats_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["statistic", "comparison", "value", "p_value", "n", "bugsigdb_vintage"])
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
                "aerobe_vs_anaerobe_paired_signs",
                f"{correlation.odds_ratio:.6f}",
                f"{correlation.mcnemar_p:.6g}",
                up_down + down_up,
                vintage or "",
            ]
        )
        for pair, (rho, p_value, n_cells) in sorted(correlation.spearman.items()):
            writer.writerow(
                [
                    "spearman_rho",
                    pair,
                    f"{rho:.6f}",
                    f"{p_value:.6g}",
                    n_cells,
                    vintage or "",
                ]
            )

    return [*written, stats_path]


GUILD_COLORSCALES = {"aerobe": "OrRd", "anaerobe": "YlGn"}
# Colour-bar ticks as percentages with trailing zeros trimmed, so 0.0004 reads
# "0.04%" rather than plotly's default SI-prefixed "400u".
PERCENT_TICKS = ".3~%"
NO_SIGNAL_COLOR = "#4d4d4d"
CLUSTER_GAP_ROWS = 1


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

    caveat = (
        "A pooled Other mixes unrelated labels, so no cell in it is one condition or one "
        "anatomical site."
        if condition.bucketed or site.bucketed
        else ""
    )
    return [
        f"Rows and columns are a cut of the corpus: {condition_part} and {site_part}.",
        f"The remaining {condition.pooled_labels} conditions are "
        f"{fate(condition)} and the remaining {site.pooled_labels} body sites are "
        f"{fate(site)}.",
        *([caveat] if caveat else []),
    ]


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
) -> Path:
    """Render the heatmap with one colour scale per trait group.

    Each group is drawn as its own overlaid trace so it keeps its own hue. With
    ``shared_range`` the traces share the range 0 to the largest observed value,
    which makes a cell in one group directly comparable in magnitude to a cell in
    another; that only holds where the groups are commensurable, as the two oxygen
    guilds are. The Table 1 traits differ in corpus-wide prevalence by more than
    an order of magnitude, so they are drawn on per-trait ranges instead and each
    colour bar carries its own maximum.

    Cells where the site was sampled but no group annotations were recorded are filled
    dark grey, so zero mapped signal is distinguishable from a small fraction. With
    ``mark_no_signal`` off they are left blank instead: the shading dominates the
    figure where the traits are rare, and the state is still recorded in the
    emitted cell-status table.
    """
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
    bar_len = min(0.42, 0.92 / count)
    for position, guild in enumerate(guild_order):
        z_guild = per_group[guild]
        group_values = [v for row in z_guild for v in row if v is not None]
        scale_max = z_max if shared_range else _upper_range(group_values, range_quantile)
        # Evenly space the colour bars down the right margin, top group first.
        centre = 0.5 + (count - 1) / 2.0 * (bar_len + 0.04) - position * (bar_len + 0.04)
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
        cluster_note += " Dark grey = sampled, with no recorded annotation signal."
    # Broken into explicit lines: as one run the subtitle overflows the plotting
    # area and is clipped by the right-hand colour bars.
    subtitle = "<br>".join(
        [
            f"Mean {value_phrase}: % of retained annotation counts; averaged within studies, "
            "then across studies.",
            "Taxonomic rollup can count one reported taxon many times; values are not "
            "percentages of taxa or patients.",
            range_note,
            cluster_note,
            *_coverage_note(matrix),
            "Directions are curated signature labels; annotation shares do not measure "
            "patient traits or exposure. See methods/biological_interpretation.md.",
            f"BugSigDB {vintage or 'vintage unknown'}.",
        ]
    )
    figure.update_layout(
        title={"text": f"{title}<br><sub>{subtitle}</sub>"},
        xaxis={"title": f"body site | {group_noun} _ direction", "tickangle": -60},
        yaxis={"title": "condition", "autorange": "reversed"},
        shapes=shapes,
        annotations=annotations,
        width=max(900, 26 * len(matrix.columns)),
        height=max(600, 26 * len(row_labels) + 380),
        # The site headers sit above the plotting area, so they need their own room.
        # Caption lines plus a two-line site header need the room.
        margin={"l": 300, "b": 300, "t": 445},
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.write_image(str(output_path), scale=2)
    return output_path


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
) -> Path:
    """Draw one panel per named body site, at a width where the labels are readable.

    The combined figure fits every site by making each column narrow enough that
    the rotated `site | trait _ direction` labels run together. This pulls a few
    blocks out into their own panels: same rows in the same clustered order, same
    per-trait colour ranges as the combined figure so a cell here is the same
    colour there, but with the site named once above the panel and the columns
    labelled by trait and direction alone.
    """
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    scales: dict[str, object] = dict(GUILD_COLORSCALES) if colorscales is None else colorscales
    labels = colorbar_labels or {}
    row_labels, row_values, row_blanks = _insert_cluster_gaps(matrix)

    missing = [site for site in sites if site not in matrix.body_sites]
    if missing:
        raise ValueError(f"body sites absent from the matrix: {', '.join(missing)}")

    # Colour ranges come from the whole matrix, not from the panels, so the
    # panels are readable at the same colours as the combined figure.
    scale_max = {
        guild: _upper_range(
            [
                value
                for row in row_values
                for column, value in zip(matrix.columns, row, strict=True)
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

    count = len(guild_order)
    bar_len = min(0.42, 0.92 / count)
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
            centre = 0.5 + (count - 1) / 2.0 * (bar_len + 0.04) - position * (bar_len + 0.04)
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
                        "title": f"{labels.get(guild, guild)}<br>% of trait assertions",
                        "len": bar_len,
                        "y": centre,
                        "yanchor": "middle",
                        "thickness": 14,
                        "tickformat": PERCENT_TICKS,
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

    range_note = f"Each {group_noun} is on its own scale, the same one the combined figure uses"
    if range_quantile < 1.0:
        range_note += f", cut at the {range_quantile:.0%} quantile of observed cells"
    grey_note = (
        " Dark grey = sampled, with no recorded annotation signal." if mark_no_signal else ""
    )
    subtitle = "<br>".join(
        [
            f"Mean {value_phrase}: % of retained annotation counts; averaged within studies, "
            "then across studies.",
            "Taxonomic rollup can count one reported taxon many times; values are not "
            "percentages of taxa or patients.",
            f"{range_note}.",
            "Rows are the conditions of the combined figure in the same clustered order, so a "
            f"blank row here was drawn in another body site.{grey_note}",
            *_coverage_note(matrix),
            "Directions are curated signature labels; annotation shares do not measure "
            "patient traits or exposure. See methods/biological_interpretation.md.",
            f"BugSigDB {vintage or 'vintage unknown'}.",
        ]
    )

    figure.update_layout(
        title={"text": f"{title}<br><sub>{subtitle}</sub>"},
        shapes=shapes,
        width=max(900, 150 + 620 * len(sites)),
        height=max(600, 26 * len(row_labels) + 400),
        # Caption lines above the panel titles.
        margin={"l": 300, "b": 260, "t": 425},
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
    parser.add_argument("--profile-input", type=Path, default=DEFAULT_PROFILE_PATH)
    parser.add_argument("--experiments-input", type=Path, default=DEFAULT_EXPERIMENTS_PATH)
    parser.add_argument("--full-dump-input", type=Path, default=DEFAULT_FULL_DUMP_PATH)
    parser.add_argument("--output-dir", type=Path, default=DERIVED_DATA_DIR / "guild_maps")
    parser.add_argument("--max-conditions", type=int, default=60)
    parser.add_argument("--max-body-sites", type=int, default=20)
    parser.add_argument(
        "--row-clusters",
        type=int,
        default=6,
        help="Number of condition clusters to cut the Euclidean tree into.",
    )
    parser.add_argument(
        "--min-studies",
        type=int,
        default=3,
        help="Cells backed by fewer studies are masked rather than drawn.",
    )
    parser.add_argument(
        "--include-non-human",
        action="store_true",
        help="Keep non-human host species (issue #129 excludes them by default).",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    profiles, stats = load_profiles(
        args.profile_input,
        args.experiments_input,
        OXYGEN_GUILDS,
        human_only=not args.include_non_human,
    )
    vintage = read_bugsigdb_vintage(args.full_dump_input)
    write_stats_report(stats, vintage)

    matrix = build_matrix(
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
    correlation = compute_correlation(profiles, matrix)
    written = write_outputs(matrix, correlation, args.output_dir, vintage)
    written.append(
        write_figure(matrix, args.output_dir / "guild_condition_body_site_heatmap.png", vintage)
    )

    up_down, down_up = correlation.discordant
    print(
        f"\nMcNemar (exact, paired signs): discordant {up_down}/{down_up}, "
        f"odds ratio {correlation.odds_ratio:.3f} "
        f"(p={correlation.mcnemar_p:.3g}) over {correlation.studies} studies"
    )
    for path in written:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
