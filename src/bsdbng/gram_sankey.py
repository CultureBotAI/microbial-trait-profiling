"""Gram-positive / gram-negative directional Sankey over conditions and body sites.

Gram positive (METPO:1000698) and gram negative (METPO:1000699) are near
exhaustive and near equal across KG-Microbe, so raw per-direction counts mirror
each other by construction - the compositional artefact recorded in issue #132.
Each (study, condition, body site) unit therefore contributes one paired contrast,
its mean increased fraction minus its mean decreased fraction.

That contrast does not escape the artefact one level up (#223): the two classes
are near-complementary shares of one denominator, so a unit's gram-positive and
gram-negative contrasts usually have opposite signs, and the gram-negative
diagrams are close to mirrors of the gram-positive ones. Directionality is
therefore tested rather than asserted. Every flow carries an exact sign test of
its up-versus-down unit counts, and every condition and body site carries an
exact McNemar test of the paired gram-positive and gram-negative signs, the same
test the guild figure uses for aerobes against anaerobes. A flow is drawn once,
in the direction most of its units take, with width equal to the net number of
units, and only flows passing their sign test are drawn at full strength.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from scipy.stats import binomtest, false_discovery_control

from bsdbng.guild_maps import SIGN_TOLERANCE, grouped_direction_means
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
#: A flow is drawn at full strength when its sign test passes this, within its
#: gram class's Benjamini-Hochberg family; otherwise it is drawn faded.
SIGN_TEST_Q = 0.05

#: (study, condition, body site): the unit one paired contrast is reported for.
Unit = tuple[str, str, str]


@dataclass(frozen=True)
class PairFilter:
    """Experiments kept for the contrast, and why the rest were not."""

    profiles: list[Profile]
    zero_total_experiments: int
    missing_metadata_experiments: int


def drop_undefined_pairs(profiles: list[Profile]) -> PairFilter:
    """Remove whole experiments whose shares are undefined or unplaceable.

    A profile with no retained annotation has no fraction; scoring it as zero
    would pull its study's contrast toward zero and can manufacture the exact-zero
    contrasts that carry no direction. Its partner goes with it, so the corpus
    stays paired.
    """
    by_experiment: defaultdict[tuple[str, str], list[Profile]] = defaultdict(list)
    for profile in profiles:
        by_experiment[(profile.study_id, profile.experiment_id)].append(profile)
    kept: list[Profile] = []
    zero_total = missing_metadata = 0
    for sides in by_experiment.values():
        if any(not p.condition or not p.body_site for p in sides):
            missing_metadata += 1
        elif any(p.total <= 0 for p in sides):
            zero_total += 1
        else:
            kept.extend(sides)
    return PairFilter(kept, zero_total, missing_metadata)


def unit_contrasts(profiles: list[Profile]) -> dict[Unit, dict[str, float]]:
    """Mean inc-minus-dec fraction per (study, condition, body site) and gram class.

    Built on the shared paired helper, which raises if a unit's increased and
    decreased profile counts differ, so the mean is a paired contrast by
    construction rather than by assumption.
    """
    means = grouped_direction_means(
        profiles,
        GRAM_CLASSES,
        lambda profile: (profile.study_id, profile.condition, profile.body_site),
    )
    return {
        unit: {
            gram_class: per_key[f"{gram_class}_inc"] - per_key[f"{gram_class}_dec"]
            for gram_class in GRAM_CLASSES
        }
        for unit, per_key in means.items()
    }


def _sign(value: float) -> int:
    if abs(value) < SIGN_TOLERANCE:
        return 0
    return 1 if value > 0 else -1


@dataclass(frozen=True)
class SankeyEdge:
    """One gram class's condition -> body site flow, with the test behind it."""

    gram_class: str
    condition: str
    body_site: str
    units_up: int
    units_down: int
    units_tied: int
    mean_contrast: float
    sign_test_p: float
    sign_test_q: float
    mcnemar_pos_up_neg_down: int
    mcnemar_neg_up_pos_down: int
    mcnemar_p: float
    mcnemar_q: float

    @property
    def units_scored(self) -> int:
        return self.units_up + self.units_down

    @property
    def net(self) -> int:
        return self.units_up - self.units_down

    @property
    def direction(self) -> str | None:
        """The side most scored units take, or None when they split evenly."""
        if self.net == 0:
            return None
        return "inc" if self.net > 0 else "dec"

    @property
    def significant(self) -> bool:
        return self.sign_test_q <= SIGN_TEST_Q


@dataclass(frozen=True)
class Complementarity:
    """How often a unit's two gram contrasts point opposite ways (#223)."""

    units: int
    opposite: int
    same: int
    either_tied: int


def complementarity(contrasts: dict[Unit, dict[str, float]]) -> Complementarity:
    opposite = same = tied = 0
    for per_class in contrasts.values():
        pos, neg = _sign(per_class["gram positive"]), _sign(per_class["gram negative"])
        if pos == 0 or neg == 0:
            tied += 1
        elif pos != neg:
            opposite += 1
        else:
            same += 1
    return Complementarity(len(contrasts), opposite, same, tied)


def _bh(p_values: list[float]) -> list[float]:
    if not p_values:
        return []
    return [float(q) for q in false_discovery_control(p_values, method="bh")]


def build_edges(profiles: list[Profile], *, min_studies: int) -> list[SankeyEdge]:
    """Test every (gram class, condition, body site) flow backed by enough units.

    ``min_studies`` counts scored units - studies within that condition and body
    site whose contrast is not a tie. The sign-test family is every such flow of
    one gram class; the McNemar family is every such condition and body site.
    """
    contrasts = unit_contrasts(profiles)
    signs: defaultdict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    values: defaultdict[tuple[str, str, str], list[float]] = defaultdict(list)
    for (_study, condition, body_site), per_class in contrasts.items():
        signs[(condition, body_site)].append(
            (_sign(per_class["gram positive"]), _sign(per_class["gram negative"]))
        )
        for gram_class in GRAM_CLASSES:
            values[(gram_class, condition, body_site)].append(per_class[gram_class])

    tallies: dict[tuple[str, str, str], tuple[int, int, int]] = {}
    for (gram_class, condition, body_site), contrast_values in values.items():
        scored = [_sign(v) for v in contrast_values]
        up, down = scored.count(1), scored.count(-1)
        if up + down >= min_studies:
            tallies[(gram_class, condition, body_site)] = (up, down, scored.count(0))

    mcnemar: dict[tuple[str, str], tuple[int, int, float]] = {}
    for site_key in sorted({(c, b) for _g, c, b in tallies}):
        pos_up = sum(1 for pos, neg in signs[site_key] if pos > 0 > neg)
        neg_up = sum(1 for pos, neg in signs[site_key] if neg > 0 > pos)
        discordant = pos_up + neg_up
        p_value = float(binomtest(pos_up, discordant, p=0.5).pvalue) if discordant else 1.0
        mcnemar[site_key] = (pos_up, neg_up, p_value)
    mcnemar_keys = sorted(mcnemar)
    mcnemar_q = dict(zip(mcnemar_keys, _bh([mcnemar[k][2] for k in mcnemar_keys]), strict=True))

    sign_p = {
        key: float(binomtest(up, up + down, p=0.5).pvalue)
        for key, (up, down, _t) in tallies.items()
    }
    sign_q: dict[tuple[str, str, str], float] = {}
    for gram_class in GRAM_CLASSES:
        keys = sorted(k for k in tallies if k[0] == gram_class)
        sign_q.update(zip(keys, _bh([sign_p[k] for k in keys]), strict=True))

    edges = []
    for key, (up, down, tied) in tallies.items():
        gram_class, condition, body_site = key
        pos_up, neg_up, m_p = mcnemar[(condition, body_site)]
        edges.append(
            SankeyEdge(
                gram_class=gram_class,
                condition=condition,
                body_site=body_site,
                units_up=up,
                units_down=down,
                units_tied=tied,
                mean_contrast=sum(values[key]) / len(values[key]),
                sign_test_p=sign_p[key],
                sign_test_q=sign_q[key],
                mcnemar_pos_up_neg_down=pos_up,
                mcnemar_neg_up_pos_down=neg_up,
                mcnemar_p=m_p,
                mcnemar_q=mcnemar_q[(condition, body_site)],
            )
        )
    edges.sort(key=lambda e: (e.gram_class, -abs(e.net), e.condition, e.body_site))
    return edges


EDGE_HEADER = (
    "gram_class",
    "condition",
    "body_site",
    "drawn_direction",
    "units_up",
    "units_down",
    "units_tied",
    "net_units",
    "mean_contrast",
    "sign_test_p",
    "sign_test_q",
    "sign_test_family",
    "mcnemar_pos_up_neg_down",
    "mcnemar_neg_up_pos_down",
    "mcnemar_p",
    "mcnemar_q",
    "bugsigdb_vintage",
)


def write_edge_table(edges: list[SankeyEdge], path: Path, vintage: str | None) -> Path:
    """Write every tested flow and its statistics, drawn or not.

    A flow whose units split evenly has no drawn direction and is written with an
    empty ``drawn_direction`` rather than left out.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(EDGE_HEADER)
        for edge in edges:
            writer.writerow(
                [
                    edge.gram_class,
                    edge.condition,
                    edge.body_site,
                    DIRECTION_LABELS[edge.direction] if edge.direction else "",
                    edge.units_up,
                    edge.units_down,
                    edge.units_tied,
                    edge.net,
                    f"{edge.mean_contrast:.6g}",
                    f"{edge.sign_test_p:.6g}",
                    f"{edge.sign_test_q:.6g}",
                    f"bh:{edge.gram_class}",
                    edge.mcnemar_pos_up_neg_down,
                    edge.mcnemar_neg_up_pos_down,
                    f"{edge.mcnemar_p:.6g}",
                    f"{edge.mcnemar_q:.6g}",
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
    balance: Complementarity,
) -> Path | None:
    """Render one gram class and net direction as a condition -> body site Sankey."""
    import plotly.graph_objects as go

    selected = [e for e in edges if e.gram_class == gram_class and e.direction == direction]
    if not selected:
        return None

    conditions = sorted({e.condition for e in selected})
    body_sites = sorted({e.body_site for e in selected})
    nodes = conditions + body_sites
    index = {name: i for i, name in enumerate(nodes)}
    passing = sum(e.significant for e in selected)

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
                "value": [abs(e.net) for e in selected],
                "color": [
                    "rgba(60,60,60,0.55)" if e.significant else "rgba(160,160,160,0.18)"
                    for e in selected
                ],
                "label": [
                    f"{e.units_up} up / {e.units_down} down, sign-test q={e.sign_test_q:.2g}"
                    for e in selected
                ],
            },
        )
    )
    subtitle = (
        f"Flows where more units shift {label} than the other way; width = net units "
        f"(up minus down). A unit is one study within one condition and body site, so a "
        f"study spanning several sites is counted once per site. Dark: sign-test "
        f"q <= {SIGN_TEST_Q} ({passing} of {len(selected)}); faded: not."
        f"<br>The two gram classes are near-complementary shares: {balance.opposite} of "
        f"{balance.opposite + balance.same} scored units shift them in opposite directions, "
        "so each gram-negative diagram largely mirrors a gram-positive one."
        f"<br>Annotation shares, not abundance. BugSigDB {vintage or 'vintage unknown'}."
    )
    figure.update_layout(
        title={
            "text": (
                f"{gram_class.title()}, net {label} - condition to body site"
                f"<br><sub>{subtitle}</sub>"
            )
        },
        width=1100,
        height=max(600, 20 * len(nodes) + 280),
        margin={"t": 150},
        font={"size": 11},
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.write_image(str(output_path), scale=2)
    return output_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Generate gram-positive and gram-negative Sankey diagrams of tested net "
            "increased and decreased paired contrasts across conditions and body sites."
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
        help="Flows with fewer non-tied studies are not tested or drawn.",
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
        reject_missing_counts=True,
    )
    vintage = read_bugsigdb_vintage(args.full_dump_input)
    write_stats_report(stats, vintage)
    filtered = drop_undefined_pairs(profiles)
    print(f"Dropped, zero total on a side: {filtered.zero_total_experiments} experiments")
    print(f"Dropped, no condition/site:    {filtered.missing_metadata_experiments} experiments")

    contrasts = unit_contrasts(filtered.profiles)
    balance = complementarity(contrasts)
    print(
        f"Study-condition-site units: {balance.units}; gram classes shift in opposite "
        f"directions in {balance.opposite}, the same direction in {balance.same}, "
        f"and at least one is tied in {balance.either_tied}"
    )

    edges = build_edges(filtered.profiles, min_studies=args.min_studies)
    written = [write_edge_table(edges, args.output_dir / "gram_sankey_edges.tsv", vintage)]
    for gram_class in GRAM_CLASSES:
        for direction in ("inc", "dec"):
            slug = f"{gram_class.replace(' ', '_')}_{direction}"
            path = write_sankey(
                edges,
                gram_class,
                direction,
                args.output_dir / f"{slug}_sankey.png",
                vintage,
                balance,
            )
            if path is None:
                print(f"no net flows for {gram_class} {direction}")
                continue
            written.append(path)

    drawn = [e for e in edges if e.direction is not None]
    print(
        f"\nFlows tested: {len(edges)}; drawn with a net direction: {len(drawn)}; "
        f"sign-test q <= {SIGN_TEST_Q}: {sum(e.significant for e in edges)}; "
        f"McNemar q <= {SIGN_TEST_Q}: "
        f"{len({(e.condition, e.body_site) for e in edges if e.mcnemar_q <= SIGN_TEST_Q})} "
        "condition-site pairs"
    )
    for path in written:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
