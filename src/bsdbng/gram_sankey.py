"""Gram-positive / gram-negative directional Sankey over conditions and body sites.

Gram positive (METPO:1000698) and gram negative (METPO:1000699) are near
exhaustive and near equal across KG-Microbe, so raw per-direction counts mirror
each other by construction -- the compositional artefact recorded in issue #132.
Link weight is therefore the paired contrast, mean(increased fraction) minus
mean(decreased fraction) within a study, and only the studies whose contrast runs
in the flow's direction contribute to it.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from bsdbng.paths import DERIVED_DATA_DIR
from bsdbng.trait_guilds import (
    DEFAULT_EXPERIMENTS_PATH,
    DEFAULT_FULL_DUMP_PATH,
    DEFAULT_PROFILE_PATH,
    GRAM_GUILDS,
    Profile,
    load_profiles,
    read_bugsigdb_vintage,
    write_stats_report,
)

GRAM_CLASSES = ("gram positive", "gram negative")
DIRECTION_LABELS = {"inc": "increased", "dec": "decreased"}


@dataclass(frozen=True)
class SankeyEdge:
    """One condition -> body site flow for a gram class and contrast direction."""

    gram_class: str
    direction: str
    condition: str
    body_site: str
    studies: int
    mean_contrast: float


def _study_contrasts(profiles: list[Profile]) -> dict[tuple[str, str, str], dict[str, float]]:
    """Mean inc-minus-dec fraction per (study, condition, body site) and gram class."""
    sums: defaultdict[tuple[str, str, str], defaultdict[str, float]] = defaultdict(
        lambda: defaultdict(float)
    )
    experiments: defaultdict[tuple[str, str, str], set[str]] = defaultdict(set)
    for profile in profiles:
        if not profile.condition or not profile.body_site:
            continue
        key = (profile.study_id, profile.condition, profile.body_site)
        experiments[key].add(profile.experiment_id)
        sign = 1.0 if profile.inc_dec == "inc" else -1.0
        for gram_class in GRAM_CLASSES:
            sums[key][gram_class] += sign * profile.fraction(gram_class)

    # Profiles arrive as matched inc/dec pairs, so dividing the signed sum by the
    # number of experiments gives the mean per-experiment contrast for that study.
    return {
        key: {gram_class: total / len(experiments[key]) for gram_class, total in per_class.items()}
        for key, per_class in sums.items()
    }


def build_edges(profiles: list[Profile], *, min_studies: int) -> list[SankeyEdge]:
    """Aggregate per-study contrasts into condition -> body site edges."""
    contrasts = _study_contrasts(profiles)
    grouped: defaultdict[tuple[str, str, str, str], list[float]] = defaultdict(list)
    for (_study, condition, body_site), per_class in contrasts.items():
        for gram_class, contrast in per_class.items():
            if contrast == 0.0:
                continue
            direction = "inc" if contrast > 0 else "dec"
            grouped[(gram_class, direction, condition, body_site)].append(abs(contrast))

    edges = [
        SankeyEdge(
            gram_class=gram_class,
            direction=direction,
            condition=condition,
            body_site=body_site,
            studies=len(values),
            mean_contrast=sum(values) / len(values),
        )
        for (gram_class, direction, condition, body_site), values in grouped.items()
        if len(values) >= min_studies
    ]
    edges.sort(key=lambda e: (e.gram_class, e.direction, -e.studies, e.condition))
    return edges


def write_edge_table(edges: list[SankeyEdge], path: Path, vintage: str | None) -> Path:
    """Write every rendered flow so the figure can be audited without reading pixels."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "gram_class",
                "direction",
                "condition",
                "body_site",
                "studies",
                "mean_abs_contrast",
                "bugsigdb_vintage",
            ]
        )
        for edge in edges:
            writer.writerow(
                [
                    edge.gram_class,
                    DIRECTION_LABELS[edge.direction],
                    edge.condition,
                    edge.body_site,
                    edge.studies,
                    f"{edge.mean_contrast:.6f}",
                    vintage or "",
                ]
            )
    return path


def write_sankey(
    edges: list[SankeyEdge],
    gram_class: str,
    direction: str,
    output_path: Path,
    vintage: str | None,
) -> Path | None:
    """Render one gram class and direction as a condition -> body site Sankey."""
    import plotly.graph_objects as go

    selected = [e for e in edges if e.gram_class == gram_class and e.direction == direction]
    if not selected:
        return None

    conditions = sorted({e.condition for e in selected})
    body_sites = sorted({e.body_site for e in selected})
    nodes = conditions + body_sites
    index = {name: i for i, name in enumerate(nodes)}

    label = DIRECTION_LABELS[direction]
    figure = go.Figure(
        data=go.Sankey(
            node={
                "label": nodes,
                "pad": 14,
                "thickness": 14,
                "color": ["#4C78A8"] * len(conditions) + ["#72B7B2"] * len(body_sites),
            },
            link={
                "source": [index[e.condition] for e in selected],
                "target": [index[e.body_site] for e in selected],
                "value": [e.studies for e in selected],
                "label": [f"{e.studies} studies" for e in selected],
            },
        )
    )
    figure.update_layout(
        title={
            "text": (
                f"{gram_class.title()} {label} — condition to body site"
                f"<br><sub>Link width = studies whose paired inc-minus-dec contrast runs this "
                f"way. BugSigDB {vintage or 'vintage unknown'}.</sub>"
            )
        },
        width=1100,
        height=max(600, 20 * len(nodes) + 240),
        font={"size": 11},
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.write_image(str(output_path), scale=2)
    return output_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Generate gram-positive and gram-negative Sankey diagrams for increased "
            "and decreased paired contrasts across conditions and body sites."
        )
    )
    parser.add_argument("--profile-input", type=Path, default=DEFAULT_PROFILE_PATH)
    parser.add_argument("--experiments-input", type=Path, default=DEFAULT_EXPERIMENTS_PATH)
    parser.add_argument("--full-dump-input", type=Path, default=DEFAULT_FULL_DUMP_PATH)
    parser.add_argument("--output-dir", type=Path, default=DERIVED_DATA_DIR / "gram_sankey")
    parser.add_argument(
        "--min-studies",
        type=int,
        default=3,
        help="Flows supported by fewer studies are dropped rather than drawn.",
    )
    parser.add_argument("--include-non-human", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    profiles, stats = load_profiles(
        args.profile_input,
        args.experiments_input,
        GRAM_GUILDS,
        human_only=not args.include_non_human,
    )
    vintage = read_bugsigdb_vintage(args.full_dump_input)
    write_stats_report(stats, vintage)

    edges = build_edges(profiles, min_studies=args.min_studies)
    written = [
        write_edge_table(edges, args.output_dir / "gram_sankey_edges.tsv", vintage),
    ]
    for gram_class in GRAM_CLASSES:
        for direction in ("inc", "dec"):
            slug = f"{gram_class.replace(' ', '_')}_{direction}"
            path = write_sankey(
                edges, gram_class, direction, args.output_dir / f"{slug}_sankey.png", vintage
            )
            if path is None:
                print(f"no flows met the study threshold for {gram_class} {direction}")
                continue
            written.append(path)

    print(f"\nEdges rendered: {len(edges)}")
    for path in written:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
