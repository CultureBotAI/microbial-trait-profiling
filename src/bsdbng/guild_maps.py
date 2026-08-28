"""Disease-by-body-site oxygen guild heatmap and aerobe/anaerobe correlation.

Rows are conditions, columns are body sites crossed with the four
guild-by-direction combinations, so each body site contributes the 2x2
`{aerobe, anaerobe} x {inc, dec}` table unrolled into four adjacent columns.

Cell values are mean guild trait fractions. Because 60.2% of studies contribute
more than one experiment (issue #133), profiles are averaged within a study
before being averaged across studies, and the companion n-matrix counts studies
rather than experiments.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from scipy.cluster.hierarchy import fcluster, leaves_list, linkage
from scipy.spatial.distance import pdist
from scipy.stats import fisher_exact, spearmanr

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
    """How much of the corpus a top-N cut pushed into the Other bucket."""

    kind: str
    kept: int
    pooled_labels: int
    pooled_profiles: int
    total_profiles: int

    def describe(self) -> str:
        share = 100.0 * self.pooled_profiles / self.total_profiles if self.total_profiles else 0.0
        return (
            f"{self.kind}: kept top {self.kept}; Other pools {self.pooled_labels} labels "
            f"and {self.pooled_profiles} profiles ({share:.1f}%)"
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
    # True where the body site was sampled for that condition but no taxon in the
    # profile carried any trait of that guild. Without this, an all-zero cell is
    # drawn in the same colour as a genuinely low fraction.
    sampled_no_signal: list[list[bool]]


def _top_labels(counts: dict[str, int], limit: int) -> list[str]:
    """Return the most frequent labels, ranked deterministically."""
    ranked = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    return [label for label, _ in ranked[:limit]]


def _binned(label: str, keep: set[str]) -> str | None:
    if not label:
        return None
    return label if label in keep else OTHER_LABEL


def _pooling(counts: dict[str, int], keep: set[str], limit: int, kind: str) -> PooledLabels:
    pooled = {k: v for k, v in counts.items() if k not in keep}
    return PooledLabels(
        kind=kind,
        kept=limit,
        pooled_labels=len(pooled),
        pooled_profiles=sum(pooled.values()),
        total_profiles=sum(counts.values()),
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
) -> GuildMatrix:
    """Aggregate profiles into the condition x (body site x guild x direction) grid."""
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
        _pooling(dict(condition_counts), keep_conditions, max_conditions, "condition"),
        _pooling(dict(body_site_counts), keep_body_sites, max_body_sites, "body_site"),
    ]

    # study -> cell -> list of per-profile fractions, so a study is averaged once.
    per_study: defaultdict[tuple[str, str, str], defaultdict[str, list[float]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for profile in profiles:
        condition = _binned(profile.condition, keep_conditions)
        body_site = _binned(profile.body_site, keep_body_sites)
        if condition is None or body_site is None:
            continue
        for guild in GUILD_ORDER:
            key = (condition, body_site, f"{guild}_{profile.inc_dec}")
            per_study[key][profile.study_id].append(profile.fraction(guild))

    conditions = sorted(keep_conditions) + ([OTHER_LABEL] if condition_counts else [])
    body_sites = sorted(keep_body_sites) + ([OTHER_LABEL] if body_site_counts else [])
    columns = [
        f"{site}|{guild}_{direction}"
        for site in body_sites
        for guild in GUILD_ORDER
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
            for guild in GUILD_ORDER:
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

    # Drop rows that would render entirely blank. A row keeping only dark-grey
    # cells is retained: "sampled, no guild trait" is a result, not an absence.
    keep_rows = [
        i
        for i in range(len(conditions))
        if any(v is not None for v in values[i]) or any(sampled_no_signal[i])
    ]
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
        sampled_no_signal=[sampled_no_signal[i] for i in order],
    )


@dataclass(frozen=True)
class CorrelationResult:
    """Guild-by-direction association and cross-cell rank correlations."""

    table: dict[str, int]
    odds_ratio: float
    fisher_p: float
    spearman: dict[str, tuple[float, float, int]]
    studies: int


def _study_direction_means(profiles: list[Profile]) -> dict[str, dict[str, float]]:
    """Mean guild fraction per study and direction, collapsing repeated experiments."""
    sums: defaultdict[str, defaultdict[str, float]] = defaultdict(lambda: defaultdict(float))
    counts: defaultdict[str, defaultdict[str, int]] = defaultdict(lambda: defaultdict(int))
    for profile in profiles:
        for guild in GUILD_ORDER:
            key = f"{guild}_{profile.inc_dec}"
            sums[profile.study_id][key] += profile.fraction(guild)
            counts[profile.study_id][key] += 1
    return {
        study: {key: total / counts[study][key] for key, total in per_key.items()}
        for study, per_key in sums.items()
    }


def compute_correlation(profiles: list[Profile], matrix: GuildMatrix) -> CorrelationResult:
    """Test whether aerobes shift toward the increased side more often than anaerobes.

    The 2x2 table counts studies, not experiments or trait assertions. A study is
    scored once per guild by the sign of (mean increased fraction - mean decreased
    fraction), so the cells are `{aerobe, anaerobe} x {shifted up, shifted down}`.
    Counting non-zero signal instead would put nearly every study in every cell and
    drive the odds ratio to 1 regardless of the data.
    """
    study_means = _study_direction_means(profiles)

    table = {f"{guild}_{d}": 0 for guild in GUILD_ORDER for d in DIRECTIONS}
    for per_key in study_means.values():
        for guild in GUILD_ORDER:
            inc = per_key.get(f"{guild}_inc")
            dec = per_key.get(f"{guild}_dec")
            if inc is None or dec is None or inc == dec:
                continue
            table[f"{guild}_inc" if inc > dec else f"{guild}_dec"] += 1

    a, b = table["aerobe_inc"], table["aerobe_dec"]
    c, d = table["anaerobe_inc"], table["anaerobe_dec"]
    odds_ratio, fisher_p = fisher_exact([[a, b], [c, d]], alternative="two-sided")

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
        fisher_p=float(fisher_p),
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


def write_outputs(
    matrix: GuildMatrix,
    correlation: CorrelationResult,
    output_dir: Path,
    vintage: str | None,
) -> list[Path]:
    """Write the value matrix, the study-count matrix, and the statistics table."""
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

    values_path = output_dir / "guild_condition_body_site_matrix.tsv"
    status_path = output_dir / "guild_condition_body_site_cell_status.tsv"
    clusters_path = output_dir / "guild_condition_row_clusters.tsv"
    counts_path = output_dir / "guild_condition_body_site_study_counts.tsv"
    stats_path = output_dir / "guild_aerobe_anaerobe_correlation.tsv"

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

    with stats_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["statistic", "comparison", "value", "p_value", "bugsigdb_vintage"])
        for key, count in sorted(correlation.table.items()):
            writer.writerow(["studies_shifted", key, count, "", vintage or ""])
        writer.writerow(
            [
                "fisher_odds_ratio",
                "aerobe_vs_anaerobe_x_inc_vs_dec",
                f"{correlation.odds_ratio:.6f}",
                f"{correlation.fisher_p:.6g}",
                vintage or "",
                correlation.studies,
            ]
        )
        for pair, (rho, p_value, n_cells) in sorted(correlation.spearman.items()):
            writer.writerow(
                [
                    "spearman_rho",
                    pair,
                    f"{rho:.6f}",
                    f"{p_value:.6g}",
                    vintage or "",
                    n_cells,
                ]
            )

    return [values_path, counts_path, status_path, clusters_path, stats_path]


GUILD_COLORSCALES = {"aerobe": "OrRd", "anaerobe": "YlGn"}
NO_SIGNAL_COLOR = "#4d4d4d"
CLUSTER_GAP_ROWS = 1


def _guild_of_column(column: str) -> str:
    """Return the guild encoded in a `site|guild_direction` column label."""
    return column.split("|", 1)[1].rsplit("_", 1)[0]


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


def write_figure(matrix: GuildMatrix, output_path: Path, vintage: str | None) -> Path:
    """Render the heatmap with one colour scale per guild over a shared range.

    Aerobe and anaerobe are drawn as two overlaid traces so each guild keeps its
    own hue, but both use the range 0 to the larger of the two maxima, so a cell
    in one guild is directly comparable in magnitude to a cell in the other.
    Cells where the site was sampled but carried no guild trait are filled dark
    grey, so an absent guild is distinguishable from a small fraction.
    """
    import plotly.graph_objects as go

    row_labels, row_values, row_blanks = _insert_cluster_gaps(matrix)
    observed = [v for row in row_values for v in row if v is not None]
    z_max = max(observed) if observed else 1.0

    figure = go.Figure()
    for position, guild in enumerate(GUILD_ORDER):
        z_guild: list[list[float | None]] = [
            [
                value if _guild_of_column(column) == guild else None
                for column, value in zip(matrix.columns, row, strict=True)
            ]
            for row in row_values
        ]
        figure.add_trace(
            go.Heatmap(
                z=z_guild,
                x=matrix.columns,
                y=row_labels,
                colorscale=GUILD_COLORSCALES[guild],
                zmin=0.0,
                zmax=z_max,
                hoverongaps=False,
                name=guild,
                colorbar={
                    "title": f"{guild}<br>fraction",
                    "len": 0.42,
                    "y": 0.79 - 0.46 * position,
                    "yanchor": "middle",
                    "thickness": 14,
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
        if sampled
    ]

    n_clusters = len(set(matrix.row_clusters))
    subtitle = (
        "Mean guild trait fraction, averaged within study then across studies. "
        f"Both scales span 0 to {z_max:.3f}. Rows clustered by Euclidean distance into "
        f"{n_clusters} groups, separated by blank rows. "
        "Dark grey = site sampled for that condition but no taxon carried the guild trait. "
        f"BugSigDB {vintage or 'vintage unknown'}."
    )
    figure.update_layout(
        title={"text": f"Oxygen guilds by condition and body site<br><sub>{subtitle}</sub>"},
        xaxis={"title": "body site | guild _ direction", "tickangle": -60},
        yaxis={"title": "condition", "autorange": "reversed"},
        shapes=shapes,
        width=max(900, 26 * len(matrix.columns)),
        height=max(600, 26 * len(row_labels) + 300),
        margin={"l": 300, "b": 300},
    )
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

    print(
        f"\nFisher odds ratio {correlation.odds_ratio:.3f} "
        f"(p={correlation.fisher_p:.3g}) over {correlation.studies} studies"
    )
    for path in written:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
