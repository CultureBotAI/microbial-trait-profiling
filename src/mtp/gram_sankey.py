"""Gram-positive and gram-negative directional Sankey over conditions and body sites.

Gram positive (METPO:1000698) and gram negative (METPO:1000699) are near
exhaustive and near equal across KG-Microbe, so raw per-direction hit counts
mirror each other by construction (#132). Each study within one condition and
body site therefore contributes one paired contrast per gram class: its mean
increased fraction minus its mean decreased fraction, built by the same
study-level averaging as the oxygen guild figure.

The contrast does not escape the mirror one level up (#223). The two classes are
near-complementary shares of one denominator, so a study's gram-positive and
gram-negative contrasts usually have opposite signs and the gram-negative
diagrams largely mirror the gram-positive ones. Directionality is tested rather
than asserted:

- every condition -> body site flow of a gram class carries an exact two-sided
  binomial sign test of its up-versus-down study counts, corrected within that
  gram class;
- every condition and body site carries an exact McNemar test on the studies
  whose two gram contrasts disagree, the paired test the guild figure uses for
  aerobes against anaerobes.

A test needs at least ``min_studies`` informative studies: non-tied studies for
a sign test, discordant studies for a McNemar. That support rule is fixed before
any test and does not look at direction. Each family is the selected members,
and the edge table records the unfiltered inventory beside it (methods/
multiple_testing_policy.md R6). A flow is drawn once, in the direction most of
its studies take, with width equal to the net number of studies, and only flows
passing their sign test are drawn at full strength.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from scipy.stats import binomtest

from mtp.corpus_identity import corpus_identity, summary_lines, write_identity
from mtp.enrichment_stats import benjamini_hochberg
from mtp.guild_maps import contrast_sign, exact_mcnemar_p, require_plotly
from mtp.paths import DERIVED_DATA_DIR
from mtp.table1_trait_maps import study_trait_contrasts, tally_direction_signs
from mtp.trait_guilds import (
    DEFAULT_EXPERIMENTS_PATH,
    DEFAULT_FULL_DUMP_PATH,
    DEFAULT_PROFILE_PATH,
    GRAM_GUILDS,
    PairExclusion,
    Profile,
    drop_zero_total_pairs,
    format_pair_exclusion,
    load_profiles,
    read_bugsigdb_vintage,
    write_stats_report,
)

GRAM_POSITIVE, GRAM_NEGATIVE = tuple(GRAM_GUILDS)
GRAM_CLASSES = (GRAM_POSITIVE, GRAM_NEGATIVE)
DIRECTION_LABELS = {"inc": "increased", "dec": "decreased"}

#: Tests need at least this many informative studies; the rest are not tested or drawn.
DEFAULT_MIN_STUDIES = 3
#: A flow is drawn at full strength when its BH q is at or below this.
Q_THRESHOLD = 0.05
PRIMARY_CORRECTION = "BH"
SIGN_FAMILY_PREFIX = "gram_sankey_flow_sign"
MCNEMAR_FAMILY_ID = "gram_sankey_site_mcnemar"
OUTPUT_DIR = DERIVED_DATA_DIR / "gram_sankey"
EDGES_FILENAME = "gram_sankey_edges.tsv"
SUMMARY_FILENAME = "gram_sankey_summary.json"

STATUS_TESTED = "tested"
STATUS_BELOW_MIN = "below_min_studies"
STATUS_BELOW_MIN_DISCORDANT = "below_min_discordant_studies"

#: A flow: one condition and one body site.
Site = tuple[str, str]


def sign_family_id(gram_class: str) -> str:
    """Name the sign-test family of one gram class."""
    return f"{SIGN_FAMILY_PREFIX}:{gram_class.replace(' ', '_')}"


@dataclass(frozen=True)
class SiteFilter:
    """Profiles placed on a condition and body site, and the experiments that were not."""

    profiles: list[Profile]
    missing_site_experiments: int


def drop_unplaced_experiments(profiles: Sequence[Profile]) -> SiteFilter:
    """Remove whole experiments with no recorded condition or body site.

    Such an experiment has no flow to join. Both of its profiles go, so the
    corpus stays paired.
    """
    unplaced = {
        (p.study_id, p.experiment_id) for p in profiles if not p.condition or not p.body_site
    }
    kept = [p for p in profiles if (p.study_id, p.experiment_id) not in unplaced]
    return SiteFilter(kept, len(unplaced))


def site_contrasts(profiles: Sequence[Profile]) -> dict[Site, dict[str, dict[str, float]]]:
    """Per condition and body site, each study's inc-minus-dec mean fraction per gram class.

    Within one flow a study contributes one contrast, so a sign test's n counts
    distinct studies. Averaging goes through the guild figure's paired helper,
    which raises when a study's increased and decreased profile counts differ.
    Zero-total pairs must already be excluded (``drop_zero_total_pairs``).
    """
    by_site: defaultdict[Site, list[Profile]] = defaultdict(list)
    for profile in profiles:
        by_site[(profile.condition, profile.body_site)].append(profile)
    return {
        site: study_trait_contrasts(group, GRAM_CLASSES, exclude_zero_totals=False)
        for site, group in sorted(by_site.items())
    }


@dataclass(frozen=True)
class Complementarity:
    """How often one study's two gram contrasts point opposite ways (#223)."""

    units: int
    opposite: int
    same: int
    either_tied: int


def complementarity(contrasts: Mapping[Site, Mapping[str, Mapping[str, float]]]) -> Complementarity:
    """Count study-within-site units by agreement of their gram-positive and gram-negative signs."""
    opposite = same = tied = units = 0
    for per_study in contrasts.values():
        for per_class in per_study.values():
            units += 1
            pos = contrast_sign(per_class[GRAM_POSITIVE])
            neg = contrast_sign(per_class[GRAM_NEGATIVE])
            if pos == 0 or neg == 0:
                tied += 1
            elif pos != neg:
                opposite += 1
            else:
                same += 1
    return Complementarity(units, opposite, same, tied)


@dataclass(frozen=True)
class SiteMcNemar:
    """Exact McNemar on one condition and body site's discordant studies."""

    positive_up_negative_down: int
    negative_up_positive_down: int
    p_nominal: float | None
    bh_q: float | None
    by_q: float | None
    family_size: int
    unfiltered_size: int

    @property
    def status(self) -> str:
        return STATUS_BELOW_MIN_DISCORDANT if self.p_nominal is None else STATUS_TESTED


@dataclass(frozen=True)
class SankeyFlow:
    """One gram class's condition -> body site flow, with the tests behind it."""

    gram_class: str
    condition: str
    body_site: str
    studies_up: int
    studies_down: int
    studies_tied: int
    mean_contrast: float
    sign_p_nominal: float | None
    sign_bh_q: float | None
    sign_by_q: float | None
    sign_family_size: int
    sign_unfiltered_size: int
    mcnemar: SiteMcNemar

    @property
    def studies_scored(self) -> int:
        return self.studies_up + self.studies_down

    @property
    def net(self) -> int:
        return self.studies_up - self.studies_down

    @property
    def status(self) -> str:
        return STATUS_BELOW_MIN if self.sign_p_nominal is None else STATUS_TESTED

    @property
    def direction(self) -> str | None:
        """The side most scored studies take; None when untested or split evenly."""
        if self.status != STATUS_TESTED or self.net == 0:
            return None
        return "inc" if self.net > 0 else "dec"

    @property
    def significant(self) -> bool:
        return self.sign_bh_q is not None and self.sign_bh_q <= Q_THRESHOLD


def _corrected(
    nominal: Sequence[float | None],
) -> tuple[list[float | None], list[float | None], int]:
    """BH and BY q-values over the tested members only; untested members stay None."""
    tested = [p for p in nominal if p is not None]
    bh_iter = iter(benjamini_hochberg(tested, method="bh"))
    by_iter = iter(benjamini_hochberg(tested, method="by"))
    bh = [None if p is None else next(bh_iter) for p in nominal]
    by = [None if p is None else next(by_iter) for p in nominal]
    return bh, by, len(tested)


def site_mcnemar_tests(
    contrasts: Mapping[Site, Mapping[str, Mapping[str, float]]],
    *,
    min_studies: int = DEFAULT_MIN_STUDIES,
) -> dict[Site, SiteMcNemar]:
    """McNemar-test each condition and body site with at least ``min_studies`` discordant studies.

    The family is those sites; the rest are reported untested, with the
    unfiltered site count beside the family size.
    """
    sites = sorted(contrasts)
    counts: list[tuple[int, int]] = []
    for site in sites:
        signs = [
            (contrast_sign(row[GRAM_POSITIVE]), contrast_sign(row[GRAM_NEGATIVE]))
            for row in contrasts[site].values()
        ]
        counts.append(
            (
                sum(1 for pos, neg in signs if pos > 0 > neg),
                sum(1 for pos, neg in signs if neg > 0 > pos),
            )
        )
    nominal = [
        exact_mcnemar_p(pos_up, neg_up) if pos_up + neg_up >= min_studies else None
        for pos_up, neg_up in counts
    ]
    bh_q, by_q, family_size = _corrected(nominal)
    return {
        site: SiteMcNemar(pos_up, neg_up, p, bh, by, family_size, len(sites))
        for site, (pos_up, neg_up), p, bh, by in zip(
            sites, counts, nominal, bh_q, by_q, strict=True
        )
    }


def build_flows(
    contrasts: Mapping[Site, Mapping[str, Mapping[str, float]]],
    *,
    min_studies: int = DEFAULT_MIN_STUDIES,
) -> list[SankeyFlow]:
    """Sign-test every flow of each gram class, corrected within that gram class.

    The family of a gram class is its flows with at least ``min_studies`` non-tied
    studies. Every other flow is returned untested, and each flow carries the
    unfiltered flow count beside the family size.
    """
    if min_studies < 1:
        raise ValueError("min_studies must be at least 1")
    mcnemar = site_mcnemar_tests(contrasts, min_studies=min_studies)
    sites = sorted(contrasts)
    flows: list[SankeyFlow] = []
    for gram_class in GRAM_CLASSES:
        tallies: list[tuple[int, int, int]] = []
        nominal: list[float | None] = []
        for site in sites:
            up, down, tied = tally_direction_signs(contrasts[site], (gram_class,))[gram_class]
            tallies.append((up, down, tied))
            scored = up + down
            nominal.append(
                float(binomtest(up, scored, p=0.5).pvalue) if scored >= min_studies else None
            )
        bh_q, by_q, family_size = _corrected(nominal)
        for index, site in enumerate(sites):
            up, down, tied = tallies[index]
            values = [row[gram_class] for row in contrasts[site].values()]
            flows.append(
                SankeyFlow(
                    gram_class=gram_class,
                    condition=site[0],
                    body_site=site[1],
                    studies_up=up,
                    studies_down=down,
                    studies_tied=tied,
                    mean_contrast=sum(values) / len(values),
                    sign_p_nominal=nominal[index],
                    sign_bh_q=bh_q[index],
                    sign_by_q=by_q[index],
                    sign_family_size=family_size,
                    sign_unfiltered_size=len(sites),
                    mcnemar=mcnemar[site],
                )
            )
    return flows


EDGE_HEADER = (
    "gram_class",
    "condition",
    "body_site",
    "status",
    "drawn_direction",
    "studies_up",
    "studies_down",
    "studies_tied",
    "net_studies",
    "mean_contrast",
    "sign_test_p_nominal",
    "sign_test_bh_q",
    "sign_test_by_q",
    "sign_test_family_id",
    "sign_test_family_size",
    "sign_test_unfiltered_flows",
    "mcnemar_status",
    "mcnemar_positive_up_negative_down",
    "mcnemar_negative_up_positive_down",
    "mcnemar_p_nominal",
    "mcnemar_bh_q",
    "mcnemar_by_q",
    "mcnemar_family_id",
    "mcnemar_family_size",
    "mcnemar_unfiltered_sites",
    "primary_correction",
    "q_threshold",
    "min_studies",
    "excluded_experiment_pairs",
    "bugsigdb_vintage",
)


def _number(value: float | None) -> str:
    return "" if value is None else f"{value:.6g}"


def write_edge_table(
    flows: Sequence[SankeyFlow],
    path: Path,
    *,
    vintage: str | None,
    min_studies: int,
    exclusion: PairExclusion,
) -> Path:
    """Write every flow and its statistics, drawn or not.

    An untested flow or site has blank p and q cells and is outside its family.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(EDGE_HEADER)
        for flow in flows:
            writer.writerow(
                [
                    flow.gram_class,
                    flow.condition,
                    flow.body_site,
                    flow.status,
                    DIRECTION_LABELS[flow.direction] if flow.direction else "",
                    flow.studies_up,
                    flow.studies_down,
                    flow.studies_tied,
                    flow.net,
                    f"{flow.mean_contrast:.6g}",
                    _number(flow.sign_p_nominal),
                    _number(flow.sign_bh_q),
                    _number(flow.sign_by_q),
                    sign_family_id(flow.gram_class),
                    flow.sign_family_size,
                    flow.sign_unfiltered_size,
                    flow.mcnemar.status,
                    flow.mcnemar.positive_up_negative_down,
                    flow.mcnemar.negative_up_positive_down,
                    _number(flow.mcnemar.p_nominal),
                    _number(flow.mcnemar.bh_q),
                    _number(flow.mcnemar.by_q),
                    MCNEMAR_FAMILY_ID,
                    flow.mcnemar.family_size,
                    flow.mcnemar.unfiltered_size,
                    PRIMARY_CORRECTION,
                    Q_THRESHOLD,
                    min_studies,
                    exclusion.zero_total_experiments,
                    vintage or "",
                ]
            )
    return path


def _pass_counts(
    nominal: Sequence[float | None], bh: Sequence[float | None], by: Sequence[float | None]
) -> dict[str, object]:
    tested_bh = [q for q in bh if q is not None]
    return {
        "nominal_p_below_0_05": sum(1 for p in nominal if p is not None and p < 0.05),
        "bh_q_passing": sum(1 for q in tested_bh if q <= Q_THRESHOLD),
        "by_q_passing": sum(1 for q in by if q is not None and q <= Q_THRESHOLD),
        "smallest_bh_q": min(tested_bh, default=None),
    }


def summarize(
    flows: Sequence[SankeyFlow], balance: Complementarity, *, min_studies: int
) -> dict[str, object]:
    """Family sizes, pass counts, and the two-class complementarity, for the summary JSON."""
    per_class: dict[str, dict[str, object]] = {}
    for gram_class in GRAM_CLASSES:
        members = [f for f in flows if f.gram_class == gram_class]
        per_class[gram_class] = {
            "family_id": sign_family_id(gram_class),
            "family_size": sum(f.status == STATUS_TESTED for f in members),
            "unfiltered_flows": len(members),
            "drawn": sum(f.direction is not None for f in members),
            **_pass_counts(
                [f.sign_p_nominal for f in members],
                [f.sign_bh_q for f in members],
                [f.sign_by_q for f in members],
            ),
        }
    sites = list({(f.condition, f.body_site): f.mcnemar for f in flows}.values())
    return {
        "min_studies": min_studies,
        "primary_correction": PRIMARY_CORRECTION,
        "q_threshold": Q_THRESHOLD,
        "sign_test": per_class,
        "mcnemar": {
            "family_id": MCNEMAR_FAMILY_ID,
            "family_size": sum(m.status == STATUS_TESTED for m in sites),
            "unfiltered_sites": len(sites),
            **_pass_counts(
                [m.p_nominal for m in sites], [m.bh_q for m in sites], [m.by_q for m in sites]
            ),
        },
        "complementarity": {
            "units": balance.units,
            "opposite": balance.opposite,
            "same": balance.same,
            "either_tied": balance.either_tied,
        },
    }


def figure_subtitle(
    selected: Sequence[SankeyFlow],
    gram_class: str,
    direction: str,
    balance: Complementarity,
    vintage: str | None,
) -> str:
    """State what the widths mean, what passed, and why the classes mirror each other."""
    passing = sum(f.significant for f in selected)
    family = sign_family_id(gram_class)
    if passing:
        verdict = (
            f"Dark: sign-test BH q <= {Q_THRESHOLD} in {family} "
            f"({passing} of {len(selected)}); faded: not."
        )
    else:
        verdict = (
            f"No flow passes its sign test (BH q <= {Q_THRESHOLD} in {family}); every link "
            "is faded. The diagram shows where study counts lean, not tested shifts."
        )
    return (
        f"Flows where more studies report {gram_class} annotations "
        f"{DIRECTION_LABELS[direction]} than the other way; width = net studies "
        "(up minus down) within one condition and body site."
        f"<br>{verdict}<br>The two gram classes are near-complementary shares: "
        f"{balance.opposite} of {balance.opposite + balance.same} scored units shift them in "
        "opposite directions, so each gram-negative diagram largely mirrors a gram-positive one."
        "<br>Trait-annotation shares of reported taxa, not abundance. Human 16S and shotgun "
        f"experiments. BugSigDB {vintage or 'vintage unknown'}."
    )


def write_sankey(
    flows: Sequence[SankeyFlow],
    gram_class: str,
    direction: str,
    output_path: Path,
    *,
    vintage: str | None,
    balance: Complementarity,
) -> Path | None:
    """Render one gram class and net direction as a condition -> body site Sankey."""
    selected = [f for f in flows if f.gram_class == gram_class and f.direction == direction]
    if not selected:
        return None
    require_plotly()
    import plotly.graph_objects as go

    conditions = sorted({f.condition for f in selected})
    body_sites = sorted({f.body_site for f in selected})
    nodes = conditions + body_sites
    condition_index = {name: i for i, name in enumerate(conditions)}
    site_index = {name: len(conditions) + i for i, name in enumerate(body_sites)}
    figure = go.Figure(
        data=go.Sankey(
            node={
                "label": nodes,
                "pad": 14,
                "thickness": 14,
                "color": ["#4C78A8"] * len(conditions) + ["#72B7B2"] * len(body_sites),
            },
            link={
                "source": [condition_index[f.condition] for f in selected],
                "target": [site_index[f.body_site] for f in selected],
                "value": [abs(f.net) for f in selected],
                "color": [
                    "rgba(60,60,60,0.55)" if f.significant else "rgba(160,160,160,0.25)"
                    for f in selected
                ],
                "label": [
                    f"{f.studies_up} up / {f.studies_down} down, sign-test BH q="
                    f"{_number(f.sign_bh_q)}"
                    for f in selected
                ],
            },
        )
    )
    label = DIRECTION_LABELS[direction]
    subtitle = figure_subtitle(selected, gram_class, direction, balance, vintage)
    figure.update_layout(
        title={
            "text": (
                f"{gram_class.capitalize()}, net {label}: condition to body site"
                f"<br><sub>{subtitle}</sub>"
            ),
            "x": 0.01,
            "y": 1,
            "yref": "container",
            "yanchor": "top",
            "pad": {"t": 20},
        },
        width=1400,
        height=max(700, 22 * max(len(conditions), len(body_sites)) + 300),
        margin={"t": 190, "l": 20, "r": 20, "b": 20},
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
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument(
        "--min-studies",
        type=int,
        default=DEFAULT_MIN_STUDIES,
        help="Tests need at least this many non-tied (sign) or discordant (McNemar) studies.",
    )
    parser.add_argument(
        "--no-figures",
        action="store_true",
        help="Write the edge table and summary only, without rendering PNGs.",
    )
    parser.add_argument("--code-commit", help="Source revision, supplied by the Justfile recipe")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    profiles, stats = load_profiles(args.profile_input, args.experiments_input, GRAM_GUILDS)
    vintage = read_bugsigdb_vintage(args.full_dump_input)
    write_stats_report(stats, vintage)
    identity = corpus_identity(args.profile_input, bugsigdb_vintage=vintage)
    print("\n".join(summary_lines(identity)))

    exclusion = drop_zero_total_pairs(profiles)
    print(format_pair_exclusion(exclusion))
    placed = drop_unplaced_experiments(exclusion.profiles)
    print(f"Excluded {placed.missing_site_experiments} experiments with no condition or body site")

    contrasts = site_contrasts(placed.profiles)
    balance = complementarity(contrasts)
    flows = build_flows(contrasts, min_studies=args.min_studies)
    summary = summarize(flows, balance, min_studies=args.min_studies)
    summary["excluded_zero_total_experiment_pairs"] = exclusion.zero_total_experiments
    summary["excluded_unplaced_experiments"] = placed.missing_site_experiments
    summary["bugsigdb_vintage"] = vintage

    output_dir: Path = args.output_dir
    written = [
        write_edge_table(
            flows,
            output_dir / EDGES_FILENAME,
            vintage=vintage,
            min_studies=args.min_studies,
            exclusion=exclusion,
        )
    ]
    summary_path = output_dir / SUMMARY_FILENAME
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    written.append(summary_path)
    if not args.no_figures:
        for gram_class in GRAM_CLASSES:
            for direction in DIRECTION_LABELS:
                slug = f"{gram_class.replace(' ', '_')}_{direction}"
                path = write_sankey(
                    flows,
                    gram_class,
                    direction,
                    output_dir / f"{slug}_sankey.png",
                    vintage=vintage,
                    balance=balance,
                )
                if path is None:
                    print(f"no tested net flows for {gram_class} {direction}")
                else:
                    written.append(path)
    written.append(write_identity(identity, output_dir, code_commit=args.code_commit))

    print(json.dumps(summary, indent=2))
    for path in written:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
