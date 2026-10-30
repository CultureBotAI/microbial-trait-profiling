"""Matched-corpus Table 1-trait analog: heatmap and directional shift.

This is the oxygen guild pair of figures (`mtp.guild_maps`) re-run over the six
individual traits named in manuscript Table 1 rather than over the two pooled
oxygen guilds, with all six combined into a single figure each:

1. facultatively aerobic (METPO:1000608)
2. cellulose catabolic process (GO:0030245)
3. coccus shaped (METPO:1000668)
4. sphere shaped (METPO:1000683)
5. arsenate reduction (pathways:arsenate_reduction)
6. moderately halophilic (METPO:1000623)

The corpus definition, the study-level averaging, the masking rules, and the row
clustering are the same code path as the guild figure, so the two are directly
comparable as a matched-corpus analog. That includes the population: both
matrices come from ``build_figure_matrix``, which leaves out every experiment pair
whose retained-annotation total is zero on a side (#783, #811). This does not
claim to reproduce the published Table 1 baseline. Two things necessarily differ. The six traits are
single terms rather than pooled guilds, and their corpus-wide prevalence spans
more than an order of magnitude, so the heatmap gives each trait its own colour
range instead of one shared range. And there is no aerobe-versus-anaerobe
contrast to test, so the companion figure reports each trait's own paired
directional shift.
"""

from __future__ import annotations

import argparse
import csv
from collections.abc import Collection, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

from scipy.stats import binomtest, false_discovery_control, spearmanr

from mtp.corpus_identity import corpus_identity, summary_lines, write_identity
from mtp.guild_maps import (
    OTHER_LABEL,
    SIGN_TOLERANCE,
    build_figure_matrix,
    pair_exclusion_note,
    require_plotly,
    study_direction_means,
    subset_matrix,
    write_body_site_panels,
    write_figure,
    write_matrix_outputs,
)
from mtp.paths import DERIVED_DATA_DIR
from mtp.trait_guilds import (
    DEFAULT_EXPERIMENTS_PATH,
    DEFAULT_FULL_DUMP_PATH,
    DEFAULT_PROFILE_PATH,
    TABLE1_TRAIT_COLORS,
    TABLE1_TRAITS,
    PairExclusion,
    Profile,
    drop_zero_total_pairs,
    format_pair_exclusion,
    load_profiles,
    read_bugsigdb_vintage,
    table1_trait_label,
    write_stats_report,
)

if TYPE_CHECKING:
    # plotly is the optional "viz" extra, so it is imported for annotations only;
    # the figure builders import it at call time behind require_plotly().
    from plotly.graph_objects import Figure

TRAIT_ORDER: tuple[str, ...] = tuple(TABLE1_TRAITS)
OUTPUT_PREFIX = "table1_traits"
CONFIDENCE_LEVEL = 0.95
ANALYSIS_LABEL = "Matched-corpus Table 1-trait analog"
ANALYSIS_NOTE = (
    "Uses the pinned human paired KG-Microbe profile corpus; it is not a "
    "reproduction of the published Table 1 baseline."
)


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
# last panel draws rheumatoid arthritis and smoking behavior only. Gingivitis and
# periodontitis are both in the corpus but fail different display rules:
# periodontitis falls just outside the top-60 cut, while gingivitis has at most
# one study per cell and is dropped as blank at any cut width. The figure names
# each with its own reason rather than redrawing the group under the same title.
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
    """Requested conditions without their own drawn row, split by display state."""

    absent_from_corpus: tuple[str, ...]
    blank_after_masking: tuple[str, ...]
    cut_by_display: tuple[str, ...]
    pooled_in_other: tuple[str, ...] = ()
    pooled_not_drawn: tuple[str, ...] = ()
    left_out_of_panel: tuple[str, ...] = ()
    # The panel names the pooled Other row itself, and it is not drawn.
    other_not_drawn: bool = False
    dropped_with_body_sites: tuple[str, ...] = ()
    # In the export, but only in experiment pairs the zero-total exclusion drops.
    excluded_zero_total: tuple[str, ...] = ()


def missing_conditions(
    requested: Sequence[str],
    *,
    corpus_conditions: Collection[str],
    drawn_conditions: Sequence[str],
    blank_conditions: Collection[str] = (),
    pooled_conditions: Collection[str] = (),
    combined_conditions: Collection[str] = (),
    site_dropped_conditions: Collection[str] = (),
    zero_total_conditions: Collection[str] = (),
) -> MissingConditions:
    """Split undrawn panel conditions by why they have no row of their own.

    A requested condition can be undrawn because the export does not hold it,
    because study-count masking and no-signal display settings leave its row
    blank, because its visible cells belonged to excluded body-site blocks,
    because the top-N cut never selected it, because it
    was assigned to the pooled Other row, or because the combined figure draws it
    but the panel's own cut does not. Testing
    membership of the drawn conditions alone cannot tell these apart, since that
    set is already the cut one — and naming the wrong cause states something
    false about the data. `blank_conditions` is the matrix's own record of rows it
    dropped as entirely blank, while ``site_dropped_conditions`` records rows
    lost with their body-site blocks. These causes are read from the build rather
    than guessed at. ``pooled_conditions`` lists actual Other membership from the
    full matrix build; the presence of an Other row alone cannot establish it.
    ``drawn_conditions`` belongs to the figure being annotated, which may omit
    the combined figure's Other row. ``combined_conditions`` are the rows the
    combined figure draws: one of those missing from a panel was removed by the
    panel's own cut, not by the top-N width cut, so it is reported apart from
    ``cut_by_display``. Other membership is assignment, not contribution: site
    binning or masking can leave a pooled condition with no profile in any drawn
    Other cell, so the note says "assigned to" and claims no more. A panel may
    name ``OTHER_LABEL`` itself; it is never a corpus condition, so it is
    reported only as drawn or not. Repeated requests are counted once.

    ``corpus_conditions`` are the conditions of the matrix's population, the
    shared comparison basis. ``zero_total_conditions`` are those the export holds
    only in experiment pairs with a zero retained-annotation total on a side.
    The basis leaves those pairs out, so such a condition is neither absent from
    the export nor cut by the display, and is reported apart from both (#811).
    """
    drawn = set(drawn_conditions)
    blank = set(blank_conditions)
    pooled = set(pooled_conditions)
    combined = set(combined_conditions)
    site_dropped = set(site_dropped_conditions)
    zero_total = set(zero_total_conditions) - set(corpus_conditions)
    requested = tuple(dict.fromkeys(requested))
    held = [name for name in requested if name in corpus_conditions and name not in drawn]
    left_out = tuple(name for name in held if name in combined)
    held = [name for name in held if name not in combined]
    pooled_missing = tuple(name for name in held if name in pooled)
    return MissingConditions(
        absent_from_corpus=tuple(
            name
            for name in requested
            if name not in corpus_conditions and name not in zero_total and name != OTHER_LABEL
        ),
        blank_after_masking=tuple(
            name
            for name in held
            if name in blank and name not in pooled and name not in site_dropped
        ),
        cut_by_display=tuple(
            name
            for name in held
            if name not in blank and name not in pooled and name not in site_dropped
        ),
        pooled_in_other=pooled_missing if OTHER_LABEL in drawn else (),
        pooled_not_drawn=pooled_missing if OTHER_LABEL not in drawn else (),
        left_out_of_panel=left_out,
        other_not_drawn=OTHER_LABEL in requested and OTHER_LABEL not in drawn,
        dropped_with_body_sites=tuple(
            name for name in held if name in site_dropped and name not in pooled
        ),
        excluded_zero_total=tuple(name for name in requested if name in zero_total),
    )


def missing_condition_note(
    missing: MissingConditions, *, max_conditions: int, min_studies: int
) -> str | None:
    """Phrase each cause separately, or return None when everything is drawn."""
    parts = []
    if missing.absent_from_corpus:
        parts.append(
            "Not drawn because this BugSigDB export does not contain them: "
            f"{', '.join(missing.absent_from_corpus)}."
        )
    if missing.excluded_zero_total:
        parts.append(
            "In this export, but only in experiment pairs with a zero retained-annotation "
            "total on a side, which are excluded, so not drawn: "
            f"{', '.join(missing.excluded_zero_total)}."
        )
    if missing.blank_after_masking:
        parts.append(
            "Sampled in this export, but blank after study-count masking "
            f"(minimum {min_studies} studies) and no-signal display settings, so not drawn: "
            f"{', '.join(missing.blank_after_masking)}."
        )
    if missing.dropped_with_body_sites:
        parts.append(
            "Sampled in this export, but its visible cells were in body-site blocks "
            "excluded from this figure, so not drawn: "
            f"{', '.join(missing.dropped_with_body_sites)}."
        )
    if missing.cut_by_display:
        parts.append(
            f"Sampled in this export but outside the top-{max_conditions} condition "
            f"cut, so not drawn: {', '.join(missing.cut_by_display)}."
        )
    if missing.pooled_in_other:
        parts.append(
            "Assigned to the pooled Other row rather than drawn as a row of its own: "
            f"{', '.join(missing.pooled_in_other)}."
        )
    if missing.pooled_not_drawn:
        parts.append(
            "Assigned to the pooled Other row, which is not drawn in this panel: "
            f"{', '.join(missing.pooled_not_drawn)}."
        )
    if missing.other_not_drawn:
        parts.append("The pooled Other row is not drawn in this panel.")
    if missing.left_out_of_panel:
        parts.append(
            "Drawn in the combined figure but left out of this panel's cut: "
            f"{', '.join(missing.left_out_of_panel)}."
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
    # None when no study is scored: the share and its interval are undefined, and
    # the trait keeps its place in the six-trait family at p = 1 without a p or q
    # being reported for it.
    share_up: float | None
    ci_low: float | None
    ci_high: float | None
    p_value: float | None
    q_value: float | None

    @property
    def studies_scored(self) -> int:
        return self.studies_up + self.studies_down

    @property
    def status(self) -> str:
        if self.studies_scored:
            return "tested"
        return "all_tied" if self.studies_tied else "no_studies"


@dataclass(frozen=True)
class ShiftResult:
    """Per-trait directional shifts plus the trait-by-trait contrast correlations."""

    shifts: list[TraitShift]
    # trait pair -> (rho, p, n studies contributing to both contrasts)
    spearman: dict[str, tuple[float, float, int]]
    studies: int


def study_trait_contrasts(
    profiles: list[Profile],
    trait_order: Sequence[str],
    *,
    exclude_zero_totals: bool = True,
) -> dict[str, dict[str, float]]:
    """Per study, the increased minus decreased mean fraction for each trait.

    Collapsing to one value per study first is the fix issue #133 asks for: 55.4%
    of the 922 studies left on the committed corpus after the zero-total
    exclusion contribute more than one experiment, so treating experiments as
    independent overstates the sample size.

    By default, experiments whose retained-annotation total is zero on either
    side are excluded via ``drop_zero_total_pairs`` to satisfy the shared
    comparison basis of methods/three_arm_statistics.md (#228). Passing
    ``exclude_zero_totals=False`` reproduces the historical S0 behavior where
    zero-total profiles were evaluated as 0.0 fraction.
    """
    if exclude_zero_totals:
        profiles = drop_zero_total_pairs(profiles).profiles
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


def tally_direction_signs(
    contrasts: Mapping[str, Mapping[str, float]], trait_order: Sequence[str]
) -> dict[str, tuple[int, int, int]]:
    """Per trait, the number of studies shifted up, shifted down, and tied.

    A study counts as shifted up for a trait when its mean increased fraction
    exceeds its mean decreased fraction, and down when the reverse holds; ties
    (including sub-tolerance float residue, see ``SIGN_TOLERANCE``) are counted
    separately and are excluded from any test built on this tally.

    Shared with ``mtp.trait_screen`` so that the widened trait family of #255
    votes by exactly this rule rather than by a second copy of it.
    """
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
    return tallies


def compute_direction_shifts(
    profiles: list[Profile],
    trait_order: Sequence[str] = TRAIT_ORDER,
    *,
    exclude_zero_totals: bool = True,
) -> ShiftResult:
    """Score every study once per trait by the sign of its paired contrast.

    A study counts as shifted up for a trait when its mean increased fraction
    exceeds its mean decreased fraction, and down when the reverse holds; ties
    (including sub-tolerance float residue, see ``SIGN_TOLERANCE``) are recorded
    separately and excluded from the test. Each trait then gets a two-sided exact
    binomial test against a half-and-half null, and the p-values of
    ``trait_order`` are Benjamini-Hochberg corrected together as one family.
    The widened family of #255 is corrected separately in ``mtp.trait_screen``.

    By default, experiments whose retained-annotation total is zero on either
    side are excluded via ``drop_zero_total_pairs`` to satisfy the shared
    comparison basis of methods/three_arm_statistics.md (#228). Passing
    ``exclude_zero_totals=False`` reproduces the historical S0 behavior.
    """
    contrasts = study_trait_contrasts(
        profiles, trait_order, exclude_zero_totals=exclude_zero_totals
    )
    tallies = tally_direction_signs(contrasts, trait_order)

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
        share: float | None = None
        low: float | None = None
        high: float | None = None
        if scored:
            test = binomtest(up, scored, p=0.5)
            interval = test.proportion_ci(confidence_level=CONFIDENCE_LEVEL)
            share, low, high = up / scored, float(interval.low), float(interval.high)
        shifts.append(
            TraitShift(
                trait=trait,
                studies_up=up,
                studies_down=down,
                studies_tied=tied,
                share_up=share,
                ci_low=low,
                ci_high=high,
                p_value=p_values[index] if scored else None,
                q_value=q_values[index] if scored else None,
            )
        )

    # These descriptive correlations still depend on profile denominators.
    # Taking differences and averaging by study does not remove that dependence;
    # closure permits positive as well as negative component correlations.
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


#: The correction family every q in the direction-shift table belongs to. Named so a
#: row read on its own says which family it was corrected in and how large that
#: family was (#289, carried to #178).
SHIFT_FAMILY_ID = "table1_six_fractions"
SHIFT_CORRECTION = "Benjamini-Hochberg"


def write_shift_table(
    result: ShiftResult,
    output_dir: Path,
    vintage: str | None,
    exclusion: PairExclusion | None = None,
    *,
    column_metadata_path: Path | None = None,
) -> list[Path]:
    """Write the per-trait shift table and the trait-pair correlation table.

    Both are long (tidy): fixed identifier columns, then one statistic per column
    with a single atomic value per cell.

    ``exclusion`` is the zero-total pair exclusion the shifts and the maps ran
    after. Its count is written on every shift row as ``excluded_experiment_pairs``,
    so a row read on its own states its population, as the guild figure's
    statistics table does (#783, #811). It is blank when no exclusion is given.
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
                "family_id",
                "family_size",
                "correction_method",
                "status",
                "excluded_experiment_pairs",
                "bugsigdb_vintage",
            ]
        )

        def cell(value: float | None, spec: str) -> str:
            return "" if value is None else format(value, spec)

        excluded = "" if exclusion is None else exclusion.zero_total_experiments
        for shift in result.shifts:
            writer.writerow(
                [
                    shift.trait,
                    table1_trait_label(shift.trait, column_metadata_path=column_metadata_path),
                    shift.studies_up,
                    shift.studies_down,
                    shift.studies_tied,
                    shift.studies_scored,
                    cell(shift.share_up, ".6f"),
                    cell(shift.ci_low, ".6f"),
                    cell(shift.ci_high, ".6f"),
                    cell(shift.p_value, ".6g"),
                    cell(shift.q_value, ".6g"),
                    SHIFT_FAMILY_ID,
                    len(result.shifts),
                    SHIFT_CORRECTION,
                    shift.status,
                    excluded,
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


def build_shift_figure(
    result: ShiftResult,
    vintage: str | None,
    trait_order: Sequence[str],
    *,
    column_metadata_path: Path | None = None,
) -> Figure:
    """Build all six traits' directional shifts and their correlations as one figure.

    Left panel: the share of studies whose paired contrast is positive, with an
    exact binomial 95% interval, against the half-and-half null. Right panel: the
    rank correlation between traits across per-study contrasts, which is what
    shows that coccus shaped and sphere shaped are close to the same measurement.

    Separate from ``write_shift_figure`` so the undefined-share branches below -
    the absent marker, the absent error bar, the fallback annotation column and
    the "share undefined" text - are reachable in a test without a kaleido image
    export (#289). Nothing between the two functions depends on the output path.
    """
    require_plotly()
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    by_trait = {shift.trait: shift for shift in result.shifts}
    ordered = [by_trait[t] for t in trait_order]
    # Plot top-to-bottom in Table 1 order; the categorical axis is reversed below.
    labels = [
        table1_trait_label(s.trait, column_metadata_path=column_metadata_path) for s in ordered
    ]

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

    # A trait with no scored study has no share. It keeps its row but no marker: a
    # marker at 0 would read as every study shifting down.
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
                "array": [
                    None if s.ci_high is None or s.share_up is None else s.ci_high - s.share_up
                    for s in ordered
                ],
                "arrayminus": [
                    None if s.ci_low is None or s.share_up is None else s.share_up - s.ci_low
                    for s in ordered
                ],
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
    annotation_x = max((s.ci_high for s in ordered if s.ci_high is not None), default=0.5) + 0.012
    for shift, label in zip(ordered, labels, strict=True):
        figure.add_annotation(
            x=annotation_x,
            y=label,
            xref="x",
            yref="y",
            text=(
                f"{shift.studies_up}/{shift.studies_scored} studies, q={shift.q_value:.2g}"
                if shift.q_value is not None
                else f"no scored study ({shift.studies_tied} tied); share undefined"
            ),
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
        "Both fractions and their contrasts remain denominator-dependent. "
        f"BugSigDB {vintage or 'vintage unknown'}."
    )
    figure.update_layout(
        title={
            "text": (
                f"{ANALYSIS_LABEL}: paired directional shift<br>"
                f"<sub>{ANALYSIS_NOTE} {subtitle}</sub>"
            )
        },
        width=1650,
        height=640,
        margin={"l": 290, "r": 120, "t": 170, "b": 90},
        plot_bgcolor="white",
    )
    figure.update_xaxes(title_text="share of studies with increased > decreased", row=1, col=1)
    figure.update_yaxes(autorange="reversed", row=1, col=1)
    figure.update_xaxes(tickangle=-40, row=1, col=2)
    figure.update_yaxes(autorange="reversed", row=1, col=2)
    return figure


def write_shift_figure(
    result: ShiftResult,
    output_path: Path,
    vintage: str | None,
    trait_order: Sequence[str],
    *,
    column_metadata_path: Path | None = None,
) -> Path:
    """Render the directional-shift figure to an image file."""
    figure = build_shift_figure(
        result, vintage, trait_order, column_metadata_path=column_metadata_path
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.write_image(str(output_path), scale=2)
    return output_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Generate the matched-corpus Table 1-trait analog: a condition x "
            "body-site heatmap with all six traits combined, and their paired "
            "directional shifts. This is not a reproduction of the published "
            "Table 1 baseline."
        )
    )
    parser.add_argument("--core-run", type=Path, help="Verified core workflow run manifest")
    parser.add_argument("--profile-input", type=Path)
    parser.add_argument("--columns", type=Path)
    parser.add_argument("--experiments-input", type=Path)
    parser.add_argument("--full-dump-input", type=Path)
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
    parser.add_argument("--code-commit", help="Source revision, supplied by the Justfile recipe")
    return parser


def _absent_site_reason(
    site: str,
    *,
    dropped: bool,
    pooled: bool,
    other_drawn: bool,
    max_body_sites: int,
    min_site_conditions: int,
) -> str:
    """Why a body site a focus panel names is not in the figure (#324, #413)."""
    if dropped:
        return f"{site} (dropped: fewer than {min_site_conditions} conditions)"
    if pooled and other_drawn:
        return f"{site} (pooled into {OTHER_LABEL} by the top-{max_body_sites} body-site cut)"
    if pooled:
        return f"{site} (outside the top-{max_body_sites} body-site cut)"
    return f"{site} (not in this figure's body sites)"


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
    args.columns = args.columns or args.profile_input.with_name("taxa_trait_matrix_columns.tsv")
    args.experiments_input = args.experiments_input or DEFAULT_EXPERIMENTS_PATH
    args.full_dump_input = args.full_dump_input or DEFAULT_FULL_DUMP_PATH
    profiles, stats = load_profiles(
        args.profile_input,
        args.experiments_input,
        TABLE1_TRAITS,
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

    # One population for the maps and the shift figure: the shared comparison
    # basis, as the oxygen guild figure has it (#783). A zero-total pair's share is
    # undefined, so it must not enter a map cell as a measured 0 (#811). The shift
    # drops it on its own too (#228).
    loaded_conditions = {profile.condition for profile in profiles if profile.condition}
    exclusion = drop_zero_total_pairs(profiles)
    print(format_pair_exclusion(exclusion))
    profiles = exclusion.profiles
    population_note = pair_exclusion_note(exclusion)

    matrix = build_figure_matrix(
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
    if matrix.dropped_site_conditions:
        print(
            "conditions dropped with excluded body-site blocks: "
            f"{len(matrix.dropped_site_conditions)} "
            f"({', '.join(sorted(matrix.dropped_site_conditions))})"
        )

    corpus_conditions = {profile.condition for profile in profiles if profile.condition}
    pooled_conditions = {
        name
        for pool in matrix.pooled
        if pool.kind == "condition" and pool.bucketed
        for name in pool.pooled_names
    }

    result = compute_direction_shifts(profiles, TRAIT_ORDER)

    written = write_matrix_outputs(matrix, args.output_dir, OUTPUT_PREFIX)
    written.extend(
        write_shift_table(
            result,
            args.output_dir,
            vintage,
            exclusion,
            column_metadata_path=args.columns,
        )
    )
    written.append(
        write_figure(
            matrix,
            args.output_dir / f"{OUTPUT_PREFIX}_condition_body_site_heatmap.png",
            vintage,
            guild_order=TRAIT_ORDER,
            colorscales={trait: trait_colorscale(trait) for trait in TRAIT_ORDER},
            title=f"{ANALYSIS_LABEL}: condition and body site",
            group_noun="trait",
            value_phrase="trait fraction",
            shared_range=False,
            colorbar_labels={t: t.replace("_", " ") for t in TRAIT_ORDER},
            range_quantile=0.98,
            separate_body_sites=True,
            mark_no_signal=args.mark_no_signal,
            note=population_note,
        )
    )
    # The combined figure squeezes every body site in, which leaves each column
    # too narrow to label and the sparse blocks unreadable next to the gut. The
    # same matrix is cut into the four smaller heatmaps of FOCUS_PANELS, all on
    # the combined figure's colour ranges so a cell is the same colour in each.
    for panel in FOCUS_PANELS:
        # A panel naming a body site the figure does not draw (e.g. Other under
        # --drop-other-body-sites, or a site too few conditions reach) cannot be
        # cut; drawing the remaining sites would mislabel it, so skip it (#324).
        absent_sites = [site for site in panel.sites or () if site not in matrix.body_sites]
        if absent_sites:
            pooled_sites = {
                name
                for pool in matrix.pooled
                if pool.kind == "body_site"
                for name in pool.pooled_names
            }
            reasons = [
                _absent_site_reason(
                    site,
                    dropped=site in matrix.dropped_body_sites,
                    pooled=site in pooled_sites,
                    other_drawn=OTHER_LABEL in matrix.body_sites,
                    max_body_sites=args.max_body_sites,
                    min_site_conditions=args.min_site_conditions,
                )
                for site in absent_sites
            ]
            print(
                f"focus panel {panel.slug}: names body sites this figure does not draw, "
                f"skipped: {', '.join(reasons)}"
            )
            continue
        cut = subset_matrix(matrix, conditions=panel.conditions, sites=panel.sites)
        missing = missing_conditions(
            panel.conditions or (),
            corpus_conditions=corpus_conditions,
            drawn_conditions=cut.conditions,
            blank_conditions=matrix.dropped_empty_conditions,
            pooled_conditions=pooled_conditions,
            combined_conditions=matrix.conditions,
            site_dropped_conditions=matrix.dropped_site_conditions,
            zero_total_conditions=loaded_conditions - corpus_conditions,
        )
        note = missing_condition_note(
            missing, max_conditions=args.max_conditions, min_studies=args.min_studies
        )
        if missing.absent_from_corpus:
            print(
                f"focus panel {panel.slug}: absent from the corpus: "
                f"{', '.join(missing.absent_from_corpus)}"
            )
        if missing.excluded_zero_total:
            print(
                f"focus panel {panel.slug}: only in excluded zero-total pairs: "
                f"{', '.join(missing.excluded_zero_total)}"
            )
        if missing.blank_after_masking:
            print(
                f"focus panel {panel.slug}: blank after study-count masking "
                f"(minimum {args.min_studies} studies) and no-signal display settings: "
                f"{', '.join(missing.blank_after_masking)}"
            )
        if missing.dropped_with_body_sites:
            print(
                f"focus panel {panel.slug}: body-site blocks excluded, so not drawn: "
                f"{', '.join(missing.dropped_with_body_sites)}"
            )
        if missing.cut_by_display:
            print(
                f"focus panel {panel.slug}: in the corpus but outside the top-"
                f"{args.max_conditions} cut: {', '.join(missing.cut_by_display)}"
            )
        if missing.pooled_in_other:
            print(
                f"focus panel {panel.slug}: assigned to Other: {', '.join(missing.pooled_in_other)}"
            )
        if missing.pooled_not_drawn:
            print(
                f"focus panel {panel.slug}: assigned to Other, not drawn in this panel: "
                f"{', '.join(missing.pooled_not_drawn)}"
            )
        if missing.other_not_drawn:
            print(f"focus panel {panel.slug}: the pooled Other row is not drawn")
        if missing.left_out_of_panel:
            print(
                f"focus panel {panel.slug}: drawn in the combined figure, not in this panel: "
                f"{', '.join(missing.left_out_of_panel)}"
            )
        if not cut.conditions or not cut.body_sites:
            print(f"focus panel {panel.slug}: nothing drawn, skipped")
            continue
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
                title=f"{ANALYSIS_LABEL}: {panel.subtitle}",
                group_noun="trait",
                value_phrase="trait fraction",
                colorbar_labels={t: t.replace("_", " ") for t in TRAIT_ORDER},
                range_quantile=0.98,
                mark_no_signal=args.mark_no_signal,
                scale_reference=matrix,
                # A cut is shown on its own, so it states the population on a line
                # of its own too.
                note="<br>".join(line for line in (note, population_note) if line),
            )
        )
    written.append(
        write_shift_figure(
            result,
            args.output_dir / f"{OUTPUT_PREFIX}_direction_shift.png",
            vintage,
            TRAIT_ORDER,
            column_metadata_path=args.columns,
        )
    )

    print(f"\nPaired directional shift over {result.studies} studies:")
    for shift in result.shifts:
        print(
            f"  {shift.trait:<28} up {shift.studies_up:>5} / down {shift.studies_down:>5}"
            + (
                f"  share {shift.share_up:.3f}  q {shift.q_value:.3g}"
                if shift.share_up is not None and shift.q_value is not None
                else f"  {shift.status}"
            )
        )
    written.append(write_identity(identity, args.output_dir, code_commit=args.code_commit))
    if bundle is not None:
        written.append(
            write_consumer_record(
                bundle,
                args.output_dir / "core_consumption.json",
                consumer="mtp.table1_trait_maps",
                outputs=written,
            )
        )
    for path in written:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
