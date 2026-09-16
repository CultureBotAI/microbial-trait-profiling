"""Manuscript Table 1 traits: condition-by-body-site heatmap and directional shift.

This is the oxygen guild pair of figures (`bsdbng.guild_maps`) re-run over the six
individual traits of manuscript Table 1 rather than over the two pooled oxygen
guilds, with all six combined into a single figure each:

1. facultatively aerobic (METPO:1000608)
2. cellulose catabolic process (GO:0030245)
3. coccus shaped (METPO:1000668)
4. sphere shaped (METPO:1000683)
5. arsenate reduction (pathways:arsenate_reduction)
6. moderately halophilic (METPO:1000623)

The corpus definition, the study-level averaging, the masking rules, and the row
clustering are the same code path as the guild figure, so the two are directly
comparable. Two things necessarily differ. The six traits are single terms rather
than pooled guilds, and their corpus-wide prevalence spans more than an order of
magnitude, so the heatmap gives each trait its own colour range instead of one
shared range. And there is no aerobe-versus-anaerobe contrast to test, so the
companion figure reports each trait's own paired directional shift.
"""

from __future__ import annotations

import argparse
import csv
from collections.abc import Collection, Sequence
from dataclasses import dataclass
from pathlib import Path

from scipy.stats import binomtest, false_discovery_control, spearmanr

from bsdbng.guild_maps import (
    OTHER_LABEL,
    SIGN_TOLERANCE,
    build_matrix,
    study_direction_means,
    subset_matrix,
    write_body_site_panels,
    write_figure,
    write_matrix_outputs,
)
from bsdbng.paths import DERIVED_DATA_DIR
from bsdbng.trait_guilds import (
    DEFAULT_EXPERIMENTS_PATH,
    DEFAULT_FULL_DUMP_PATH,
    DEFAULT_PROFILE_PATH,
    TABLE1_TRAIT_COLORS,
    TABLE1_TRAIT_LABELS,
    TABLE1_TRAITS,
    Profile,
    load_profiles,
    read_bugsigdb_vintage,
    write_stats_report,
)

TRAIT_ORDER: tuple[str, ...] = tuple(TABLE1_TRAITS)
OUTPUT_PREFIX = "table1_traits"
CONFIDENCE_LEVEL = 0.95


@dataclass(frozen=True)
class FocusPanel:
    """One small heatmap cut from the combined figure.

    Exactly one axis is fixed; the other is derived from the data by
    `subset_matrix`, so a panel never asserts that a condition was studied in a
    body site it was not.
    """

    slug: str
    subtitle: str
    sites: tuple[str, ...] | None = None
    conditions: tuple[str, ...] | None = None


# The combined figure shows what was measured where, but the Feces block sets
# the eye's scale and everything beside it reads as empty. These four cuts group
# conditions that were sampled in the same places, so each is legible on its own:
# the gut block alone, the gut against the pooled remainder, the two cervical
# conditions, and the oral cluster the fourth cut was asked for. On v1.3.1 that
# last panel draws rheumatoid arthritis and smoking behavior only: gingivitis
# and periodontitis are in the corpus but fall outside the top-60 condition cut,
# so the figure names them rather than redrawing the group under the same title.
FOCUS_PANELS: tuple[FocusPanel, ...] = (
    FocusPanel(slug="feces", subtitle="Feces", sites=("Feces",)),
    FocusPanel(
        slug="feces_other",
        subtitle=f"Feces and {OTHER_LABEL}",
        sites=("Feces", OTHER_LABEL),
    ),
    FocusPanel(
        slug="cervical",
        subtitle="cervical neoplasia and papilloma virus infection",
        conditions=(
            "Cervical glandular intraepithelial neoplasia",
            "Human papilloma virus infection",
        ),
    ),
    FocusPanel(
        slug="oral",
        subtitle="the oral cluster, with rheumatoid arthritis",
        conditions=(
            "Gingivitis",
            "Rheumatoid arthritis",
            "Periodontitis",
            "Smoking behavior",
        ),
    ),
)


@dataclass(frozen=True)
class MissingConditions:
    """Requested panel conditions that are not drawn, split by cause."""

    absent_from_corpus: tuple[str, ...]
    cut_by_display: tuple[str, ...]


def missing_conditions(
    requested: Sequence[str],
    *,
    corpus_conditions: Collection[str],
    drawn_conditions: Sequence[str],
) -> MissingConditions:
    """Split undrawn panel conditions into the two causes that look alike.

    A requested condition can be missing because the export does not hold it, or
    because it is held and the top-N cut dropped it. Testing membership of the
    drawn conditions alone cannot tell those apart, since that set is already the
    cut one — and reporting the second as the first states something false about
    the export.
    """
    drawn = set(drawn_conditions)
    return MissingConditions(
        absent_from_corpus=tuple(name for name in requested if name not in corpus_conditions),
        cut_by_display=tuple(
            name for name in requested if name in corpus_conditions and name not in drawn
        ),
    )


def missing_condition_note(missing: MissingConditions, *, max_conditions: int) -> str | None:
    """Phrase each cause separately, or return None when everything is drawn."""
    parts = []
    if missing.absent_from_corpus:
        parts.append(
            "Not drawn because this BugSigDB export does not contain them: "
            f"{', '.join(missing.absent_from_corpus)}."
        )
    if missing.cut_by_display:
        parts.append(
            f"Sampled in this export but outside the top-{max_conditions} condition "
            f"cut, so not drawn: {', '.join(missing.cut_by_display)}."
        )
    return " ".join(parts) if parts else None


def trait_colorscale(trait: str) -> list[list[float | str]]:
    """Build a white-to-trait-colour ramp so each trait keeps a distinct hue.

    Plotly's named sequential scales are not distinguishable six at a time, and
    the trait colours are already fixed by the shared selection table, so the
    ramp is generated from that colour rather than chosen independently here.
    """
    return [[0.0, "#ffffff"], [1.0, TABLE1_TRAIT_COLORS[trait]]]


@dataclass(frozen=True)
class TraitShift:
    """One trait's paired directional shift, counted over studies."""

    trait: str
    studies_up: int
    studies_down: int
    studies_tied: int
    share_up: float
    ci_low: float
    ci_high: float
    p_value: float
    q_value: float

    @property
    def studies_scored(self) -> int:
        return self.studies_up + self.studies_down


@dataclass(frozen=True)
class ShiftResult:
    """Per-trait directional shifts plus the trait-by-trait contrast correlations."""

    shifts: list[TraitShift]
    # trait pair -> (rho, p, n studies contributing to both contrasts)
    spearman: dict[str, tuple[float, float, int]]
    studies: int


def _study_contrasts(
    profiles: list[Profile], trait_order: Sequence[str]
) -> dict[str, dict[str, float]]:
    """Per study, the increased minus decreased mean fraction for each trait.

    Collapsing to one value per study first is the fix issue #133 asks for: 56.2%
    of studies contribute more than one experiment on the pinned v1.3.1 corpus,
    so treating experiments as independent overstates the sample size.
    """
    means = study_direction_means(profiles, trait_order)
    contrasts: dict[str, dict[str, float]] = {}
    for study, per_key in means.items():
        row: dict[str, float] = {}
        for trait in trait_order:
            inc = per_key.get(f"{trait}_inc")
            dec = per_key.get(f"{trait}_dec")
            if inc is None or dec is None:
                continue
            row[trait] = inc - dec
        if row:
            contrasts[study] = row
    return contrasts


def compute_direction_shifts(
    profiles: list[Profile], trait_order: Sequence[str] = TRAIT_ORDER
) -> ShiftResult:
    """Score every study once per trait by the sign of its paired contrast.

    A study counts as shifted up for a trait when its mean increased fraction
    exceeds its mean decreased fraction, and down when the reverse holds; ties
    (including sub-tolerance float residue, see ``SIGN_TOLERANCE``) are recorded
    separately and excluded from the test. Each trait then gets a two-sided exact
    binomial test against a half-and-half null, and the six p-values are
    Benjamini-Hochberg corrected together because the six traits are one family
    of tests.
    """
    contrasts = _study_contrasts(profiles, trait_order)

    tallies: dict[str, tuple[int, int, int]] = {}
    for trait in trait_order:
        up = down = tied = 0
        for row in contrasts.values():
            value = row.get(trait)
            if value is None:
                continue
            if value > SIGN_TOLERANCE:
                up += 1
            elif value < -SIGN_TOLERANCE:
                down += 1
            else:
                tied += 1
        tallies[trait] = (up, down, tied)

    p_values: list[float] = []
    for trait in trait_order:
        up, down, _ = tallies[trait]
        scored = up + down
        p_values.append(float(binomtest(up, scored, p=0.5).pvalue) if scored else 1.0)
    q_values = [float(q) for q in false_discovery_control(p_values, method="bh")]

    shifts: list[TraitShift] = []
    for index, trait in enumerate(trait_order):
        up, down, tied = tallies[trait]
        scored = up + down
        if scored:
            test = binomtest(up, scored, p=0.5)
            interval = test.proportion_ci(confidence_level=CONFIDENCE_LEVEL)
            share, low, high = up / scored, float(interval.low), float(interval.high)
        else:
            share = low = high = 0.0
        shifts.append(
            TraitShift(
                trait=trait,
                studies_up=up,
                studies_down=down,
                studies_tied=tied,
                share_up=share,
                ci_low=low,
                ci_high=high,
                p_value=p_values[index],
                q_value=q_values[index],
            )
        )

    # Rank correlation between traits over the per-study contrasts. Unlike the
    # within-profile guild correlations in the oxygen figure, this compares a
    # trait against another trait across studies rather than two shares of one
    # profile total, so it is not forced negative by compositionality.
    spearman: dict[str, tuple[float, float, int]] = {}
    for i, left in enumerate(trait_order):
        for right in trait_order[i + 1 :]:
            paired = [
                (row[left], row[right])
                for row in contrasts.values()
                if left in row and right in row
            ]
            if len(paired) < 3:
                continue
            left_values = [x for x, _ in paired]
            right_values = [y for _, y in paired]
            # Spearman is undefined against a constant series and would return
            # nan; omit the pair rather than write nan into the table.
            if len(set(left_values)) < 2 or len(set(right_values)) < 2:
                continue
            rho, p_value = spearmanr(left_values, right_values)
            spearman[f"{left}~{right}"] = (float(rho), float(p_value), len(paired))

    return ShiftResult(shifts=shifts, spearman=spearman, studies=len(contrasts))


def write_shift_table(result: ShiftResult, output_dir: Path, vintage: str | None) -> list[Path]:
    """Write the per-trait shift table and the trait-pair correlation table.

    Both are long (tidy): fixed identifier columns, then one statistic per column
    with a single atomic value per cell.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    shifts_path = output_dir / f"{OUTPUT_PREFIX}_direction_shift.tsv"
    spearman_path = output_dir / f"{OUTPUT_PREFIX}_contrast_correlation.tsv"

    with shifts_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "trait",
                "trait_label",
                "studies_up",
                "studies_down",
                "studies_tied",
                "studies_scored",
                "share_up",
                "ci_low",
                "ci_high",
                "binomial_p",
                "bh_q",
                "bugsigdb_vintage",
            ]
        )
        for shift in result.shifts:
            writer.writerow(
                [
                    shift.trait,
                    TABLE1_TRAIT_LABELS[shift.trait],
                    shift.studies_up,
                    shift.studies_down,
                    shift.studies_tied,
                    shift.studies_scored,
                    f"{shift.share_up:.6f}",
                    f"{shift.ci_low:.6f}",
                    f"{shift.ci_high:.6f}",
                    f"{shift.p_value:.6g}",
                    f"{shift.q_value:.6g}",
                    vintage or "",
                ]
            )

    with spearman_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["trait_a", "trait_b", "spearman_rho", "p_value", "studies", "bugsigdb_vintage"]
        )
        for pair, (rho, p_value, n_studies) in sorted(result.spearman.items()):
            left, right = pair.split("~", 1)
            writer.writerow([left, right, f"{rho:.6f}", f"{p_value:.6g}", n_studies, vintage or ""])

    return [shifts_path, spearman_path]


def write_shift_figure(
    result: ShiftResult, output_path: Path, vintage: str | None, trait_order: Sequence[str]
) -> Path:
    """Render all six traits' directional shifts and their correlations in one figure.

    Left panel: the share of studies whose paired contrast is positive, with an
    exact binomial 95% interval, against the half-and-half null. Right panel: the
    rank correlation between traits across per-study contrasts, which is what
    shows that coccus shaped and sphere shaped are close to the same measurement.
    """
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    by_trait = {shift.trait: shift for shift in result.shifts}
    ordered = [by_trait[t] for t in trait_order]
    # Plot top-to-bottom in Table 1 order; the categorical axis is reversed below.
    labels = [TABLE1_TRAIT_LABELS[s.trait] for s in ordered]

    figure = make_subplots(
        rows=1,
        cols=2,
        column_widths=[0.56, 0.44],
        horizontal_spacing=0.18,
        subplot_titles=(
            "Share of studies shifted toward the increased side",
            "Rank correlation between per-study trait contrasts",
        ),
    )

    figure.add_trace(
        go.Scatter(
            x=[s.share_up for s in ordered],
            y=labels,
            mode="markers",
            marker={
                "size": 13,
                "color": [TABLE1_TRAIT_COLORS[s.trait] for s in ordered],
                "line": {"color": "#333333", "width": 1},
            },
            error_x={
                "type": "data",
                "symmetric": False,
                "array": [s.ci_high - s.share_up for s in ordered],
                "arrayminus": [s.share_up - s.ci_low for s in ordered],
                "color": "#666666",
                "thickness": 1.5,
                "width": 6,
            },
            hovertemplate="%{y}<br>share=%{x:.3f}<extra></extra>",
            showlegend=False,
        ),
        row=1,
        col=1,
    )
    figure.add_vline(x=0.5, line={"color": "#999999", "dash": "dash"}, row=1, col=1)
    # One annotation column to the right of the widest interval, so the counts
    # line up instead of tracking each trait's own error bar.
    annotation_x = max(s.ci_high for s in ordered) + 0.012
    for shift, label in zip(ordered, labels, strict=True):
        figure.add_annotation(
            x=annotation_x,
            y=label,
            xref="x",
            yref="y",
            text=f"{shift.studies_up}/{shift.studies_scored} studies, q={shift.q_value:.2g}",
            showarrow=False,
            xanchor="left",
            font={"size": 12},
        )

    lookup = {tuple(sorted(k.split("~", 1))): v[0] for k, v in result.spearman.items()}
    heat: list[list[float | None]] = []
    text: list[list[str]] = []
    for left in trait_order:
        row: list[float | None] = []
        row_text: list[str] = []
        for right in trait_order:
            # The diagonal is 1 by construction and would be the only saturated
            # cell on the panel, so it is left empty rather than drawn.
            value = None if left == right else lookup.get(tuple(sorted((left, right))))
            row.append(value)
            row_text.append("" if value is None else f"{value:.2f}")
        heat.append(row)
        text.append(row_text)

    short = [t.replace("_", " ") for t in trait_order]
    figure.add_trace(
        go.Heatmap(
            z=heat,
            x=short,
            y=short,
            text=text,
            texttemplate="%{text}",
            textfont={"size": 12},
            colorscale="RdBu",
            reversescale=True,
            zmid=0.0,
            zmin=-1.0,
            zmax=1.0,
            hoverongaps=False,
            colorbar={"title": "rho", "len": 0.6, "thickness": 14, "x": 1.02},
        ),
        row=1,
        col=2,
    )

    subtitle = (
        "Each study contributes one value per trait: the mean increased fraction minus the "
        f"mean decreased fraction over its paired experiments. {result.studies} studies; ties "
        "are excluded from a trait's test, so the denominators differ between traits."
        "<br>Error bars are exact binomial 95% intervals; q is Benjamini-Hochberg across the "
        f"six traits. Correlations are between traits over those per-study contrasts. "
        f"BugSigDB {vintage or 'vintage unknown'}."
    )
    figure.update_layout(
        title={"text": f"Table 1 traits: paired directional shift<br><sub>{subtitle}</sub>"},
        width=1650,
        height=640,
        margin={"l": 290, "r": 120, "t": 170, "b": 90},
        plot_bgcolor="white",
    )
    figure.update_xaxes(title_text="share of studies with increased > decreased", row=1, col=1)
    figure.update_yaxes(autorange="reversed", row=1, col=1)
    figure.update_xaxes(tickangle=-40, row=1, col=2)
    figure.update_yaxes(autorange="reversed", row=1, col=2)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.write_image(str(output_path), scale=2)
    return output_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Generate the manuscript Table 1 trait pair: a condition x body-site "
            "heatmap with all six traits combined, and their paired directional shifts."
        )
    )
    parser.add_argument("--profile-input", type=Path, default=DEFAULT_PROFILE_PATH)
    parser.add_argument("--experiments-input", type=Path, default=DEFAULT_EXPERIMENTS_PATH)
    parser.add_argument("--full-dump-input", type=Path, default=DEFAULT_FULL_DUMP_PATH)
    parser.add_argument("--output-dir", type=Path, default=DERIVED_DATA_DIR / "table1_traits")
    parser.add_argument("--max-conditions", type=int, default=60)
    parser.add_argument(
        "--max-body-sites",
        type=int,
        default=8,
        help=(
            "Fewer body sites than the guild figure keeps: six traits times two "
            "directions is twelve columns per site, so 20 sites would be unreadable."
        ),
    )
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
    parser.add_argument(
        "--min-site-conditions",
        type=int,
        default=2,
        help=(
            "Drop a body site's whole column block when fewer than this many "
            "conditions are drawn in it. A block carrying one condition costs "
            "twelve columns to show a single row."
        ),
    )
    parser.add_argument(
        "--include-other-conditions",
        action="store_true",
        help=(
            "Pool the conditions outside the top-N cut into an Other row. Off by "
            "default: that row mixes hundreds of unrelated conditions, so no cell in "
            "it reads as one condition."
        ),
    )
    parser.add_argument(
        "--drop-other-body-sites",
        action="store_true",
        help=(
            "Drop the Other body-site column block instead of keeping it. It is kept "
            "by default because it carries most of the sampling outside the top few "
            "sites, but it mixes anatomical sites and is not one of them."
        ),
    )
    parser.add_argument(
        "--mark-no-signal",
        action="store_true",
        help=(
            "Shade cells sampled without any Table 1 trait dark grey. Off by default: "
            "these traits are rare enough that the shading covers most of the figure. "
            "The state stays in the emitted cell-status table either way."
        ),
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    profiles, stats = load_profiles(
        args.profile_input,
        args.experiments_input,
        TABLE1_TRAITS,
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
        guild_order=TRAIT_ORDER,
        include_other_conditions=args.include_other_conditions,
        include_other_body_sites=not args.drop_other_body_sites,
        mark_no_signal=args.mark_no_signal,
        min_site_conditions=args.min_site_conditions,
    )
    for pooled in matrix.pooled:
        print(pooled.describe())
    if matrix.dropped_body_sites:
        print(
            f"body sites dropped below {args.min_site_conditions} conditions: "
            f"{', '.join(matrix.dropped_body_sites)}"
        )
    if matrix.dropped_empty_conditions:
        print(
            f"conditions dropped as entirely blank: {len(matrix.dropped_empty_conditions)} "
            f"({', '.join(sorted(matrix.dropped_empty_conditions))})"
        )

    corpus_conditions = {profile.condition for profile in profiles if profile.condition}

    result = compute_direction_shifts(profiles, TRAIT_ORDER)

    written = write_matrix_outputs(matrix, args.output_dir, OUTPUT_PREFIX)
    written.extend(write_shift_table(result, args.output_dir, vintage))
    written.append(
        write_figure(
            matrix,
            args.output_dir / f"{OUTPUT_PREFIX}_condition_body_site_heatmap.png",
            vintage,
            guild_order=TRAIT_ORDER,
            colorscales={trait: trait_colorscale(trait) for trait in TRAIT_ORDER},
            title="Manuscript Table 1 traits by condition and body site",
            group_noun="trait",
            value_phrase="trait fraction",
            shared_range=False,
            colorbar_labels={t: t.replace("_", " ") for t in TRAIT_ORDER},
            range_quantile=0.98,
            separate_body_sites=True,
            mark_no_signal=args.mark_no_signal,
        )
    )
    # The combined figure squeezes every body site in, which leaves each column
    # too narrow to label and the sparse blocks unreadable next to the gut. The
    # same matrix is cut into the four smaller heatmaps of FOCUS_PANELS, all on
    # the combined figure's colour ranges so a cell is the same colour in each.
    for panel in FOCUS_PANELS:
        cut = subset_matrix(matrix, conditions=panel.conditions, sites=panel.sites)
        if not cut.conditions or not cut.body_sites:
            print(f"focus panel {panel.slug}: nothing drawn, skipped")
            continue
        missing = missing_conditions(
            panel.conditions or (),
            corpus_conditions=corpus_conditions,
            drawn_conditions=matrix.conditions,
        )
        note = missing_condition_note(missing, max_conditions=args.max_conditions)
        if missing.absent_from_corpus:
            print(
                f"focus panel {panel.slug}: absent from the corpus: "
                f"{', '.join(missing.absent_from_corpus)}"
            )
        if missing.cut_by_display:
            print(
                f"focus panel {panel.slug}: in the corpus but outside the top-"
                f"{args.max_conditions} cut: {', '.join(missing.cut_by_display)}"
            )
        print(
            f"focus panel {panel.slug}: {len(cut.conditions)} conditions x "
            f"{len(cut.body_sites)} body sites ({', '.join(cut.body_sites)})"
        )
        written.append(
            write_body_site_panels(
                cut,
                args.output_dir / f"{OUTPUT_PREFIX}_focus_{panel.slug}.png",
                vintage,
                sites=cut.body_sites,
                guild_order=TRAIT_ORDER,
                colorscales={trait: trait_colorscale(trait) for trait in TRAIT_ORDER},
                title=f"Manuscript Table 1 traits: {panel.subtitle}",
                group_noun="trait",
                value_phrase="trait fraction",
                colorbar_labels={t: t.replace("_", " ") for t in TRAIT_ORDER},
                range_quantile=0.98,
                mark_no_signal=args.mark_no_signal,
                scale_reference=matrix,
                note=note,
            )
        )
    written.append(
        write_shift_figure(
            result,
            args.output_dir / f"{OUTPUT_PREFIX}_direction_shift.png",
            vintage,
            TRAIT_ORDER,
        )
    )

    print(f"\nPaired directional shift over {result.studies} studies:")
    for shift in result.shifts:
        print(
            f"  {shift.trait:<28} up {shift.studies_up:>5} / down {shift.studies_down:>5}"
            f"  share {shift.share_up:.3f}  q {shift.q_value:.3g}"
        )
    for path in written:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
