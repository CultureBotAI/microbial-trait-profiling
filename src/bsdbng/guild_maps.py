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
class GuildMatrix:
    """Mean guild fractions and supporting study counts per condition/body site."""

    conditions: list[str]
    body_sites: list[str]
    columns: list[str]
    values: list[list[float | None]]
    study_counts: list[list[int]]


def _top_labels(counts: dict[str, int], limit: int) -> list[str]:
    """Return the most frequent labels, ranked deterministically."""
    ranked = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    return [label for label, _ in ranked[:limit]]


def _binned(label: str, keep: set[str]) -> str | None:
    if not label:
        return None
    return label if label in keep else OTHER_LABEL


def build_matrix(
    profiles: list[Profile],
    *,
    max_conditions: int,
    max_body_sites: int,
    min_studies: int,
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
    for condition in conditions:
        value_row: list[float | None] = []
        count_row: list[int] = []
        for site in body_sites:
            for guild in GUILD_ORDER:
                for direction in DIRECTIONS:
                    by_study: dict[str, list[float]] = per_study.get(
                        (condition, site, f"{guild}_{direction}"), {}
                    )
                    n_studies = len(by_study)
                    count_row.append(n_studies)
                    if n_studies < min_studies:
                        value_row.append(None)
                        continue
                    study_means = [sum(v) / len(v) for v in by_study.values()]
                    value_row.append(sum(study_means) / len(study_means))
        values.append(value_row)
        study_counts.append(count_row)

    return GuildMatrix(
        conditions=conditions,
        body_sites=body_sites,
        columns=columns,
        values=values,
        study_counts=study_counts,
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

    values_path = output_dir / "guild_condition_body_site_matrix.tsv"
    counts_path = output_dir / "guild_condition_body_site_study_counts.tsv"
    stats_path = output_dir / "guild_aerobe_anaerobe_correlation.tsv"

    _write_matrix_tsv(values_path, matrix, value_rows)
    _write_matrix_tsv(counts_path, matrix, count_rows)

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

    return [values_path, counts_path, stats_path]


def write_figure(matrix: GuildMatrix, output_path: Path, vintage: str | None) -> Path:
    """Render the heatmap, leaving masked cells blank rather than drawing a zero."""
    import plotly.graph_objects as go

    subtitle = (
        "Mean guild trait fraction, averaged within study then across studies. "
        f"BugSigDB {vintage or 'vintage unknown'}."
    )
    figure = go.Figure(
        data=go.Heatmap(
            z=matrix.values,
            x=matrix.columns,
            y=matrix.conditions,
            colorscale="Viridis",
            hoverongaps=False,
            colorbar={"title": "mean fraction"},
        )
    )
    figure.update_layout(
        title={"text": f"Oxygen guilds by condition and body site<br><sub>{subtitle}</sub>"},
        xaxis={"title": "body site | guild _ direction", "tickangle": -60},
        yaxis={"title": "condition", "autorange": "reversed"},
        width=max(900, 26 * len(matrix.columns)),
        height=max(600, 26 * len(matrix.conditions) + 260),
        margin={"l": 260, "b": 260},
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
    parser.add_argument("--max-conditions", type=int, default=25)
    parser.add_argument("--max-body-sites", type=int, default=8)
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
