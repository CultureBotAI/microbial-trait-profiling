"""Evidence marks on the Table 1 heatmaps (#172) and the review issues on #173.

#176: a site test counts each study once and ignores the display cuts.
#177: profile support is descriptive and carries no cell-level p or q.
#215: every departure from the composition arm is one constant, rendered on the
figure and in the methods text.
#216: a legend counts only the marks its figure draws.
#217: site glyphs are drawn once per trait, and every colour bar survives a cut.
"""

from __future__ import annotations

import csv
import re
from dataclasses import replace
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import pytest

from mtp import table1_trait_maps
from mtp.guild_maps import (
    COLUMN_MARK_GLYPH,
    MARK_LINE_COLOR,
    OTHER_LABEL,
    EvidenceMarks,
    GuildMatrix,
    MarkStyle,
    build_body_site_panels,
    build_figure,
)
from mtp.paths import REPO_ROOT
from mtp.table1_evidence import (
    CELL_COMPOSITION,
    COMPOSITION_DEPARTURES,
    INSUFFICIENT_STUDIES,
    MAJORITY_SUPPORT,
    NO_MAJORITY,
    SITE_COMPOSITION,
    SUPPORT,
    SUPPORT_COLUMNS,
    TESTED,
    CompositionResult,
    EnrichmentCalls,
    Table1Evidence,
    _correct,
    aggregate_support,
    cell_display_status,
    composition_contrasts,
    evidence_marks,
    pair_metadata,
    read_enrichment_calls,
)
from mtp.table1_trait_maps import TRAIT_ORDER, trait_colorscale
from mtp.trait_guilds import Profile

if TYPE_CHECKING:
    from plotly.graph_objects import Figure

TABLE1_COLUMNS = {
    "facultatively_aerobic": "METPO:1000608",
    "cellulose_catabolic_process": "METPO:2000103|GO:0030245",
    "coccus_shaped": "METPO:1000668",
    "sphere_shaped": "METPO:1000683",
    "arsenate_reduction": "METPO:2000103|pathways:arsenate_reduction",
    "moderately_halophilic": "METPO:1000623",
}
EXTRA_COLUMNS = [f"METPO:90000{i:02d}" for i in range(6)]
COLUMNS = [*TABLE1_COLUMNS.values(), *EXTRA_COLUMNS]
SCALES: dict[str, object] = {trait: trait_colorscale(trait) for trait in TRAIT_ORDER}
METHODS_DOC = REPO_ROOT / "methods" / "table1_trait_condition_body_site_maps.md"

Pair = tuple[str, str, str, str]  # study, experiment, condition, body site


def _profiles(pairs: list[Pair]) -> list[Profile]:
    return [
        Profile(
            study_id=study,
            experiment_id=experiment,
            inc_dec=side,
            condition=condition,
            body_site=site,
            guild_counts={},
            total=1,
        )
        for study, experiment, condition, site in pairs
        for side in ("inc", "dec")
    ]


def _counts(pairs: list[Pair], seed: int, *, shift: dict[str, int] | None = None) -> pd.DataFrame:
    """Random counts per profile; ``shift`` adds to the increased side of a study."""
    rng = np.random.default_rng(seed)
    rows: dict[str, list[float]] = {}
    for study, experiment, _, _ in pairs:
        for side in ("inc", "dec"):
            values = rng.integers(1, 9, size=len(COLUMNS)).astype(float)
            values[rng.random(len(COLUMNS)) < 0.15] = 0.0
            if side == "inc" and shift is not None:
                values[0] += shift.get(study, 0)
            rows[f"{study}|{experiment}|{side}"] = list(values)
    return pd.DataFrame.from_dict(rows, orient="index", columns=COLUMNS)


def _site_corpus() -> list[Pair]:
    """Ten Feces studies; study s0 spans three conditions with unequal experiments."""
    pairs: list[Pair] = [
        ("s0", "s0-1", "A", "Feces"),
        ("s0", "s0-2", "A", "Feces"),
        ("s0", "s0-3", "A", "Feces"),
        ("s0", "s0-4", "B", "Feces"),
        ("s0", "s0-5", "C", "Feces"),
    ]
    for index in range(1, 10):
        pairs.append((f"s{index}", f"s{index}-1", "AB"[index % 2], "Feces"))
    return pairs


def _site_rows(results: list[CompositionResult]) -> dict[str, CompositionResult]:
    return {r.trait: r for r in results if r.evidence_type == SITE_COMPOSITION}


def test_a_study_with_several_conditions_is_one_site_observation() -> None:
    """#176: grouping by (condition, study) entered s0 three times; by study, once."""
    pairs = _site_corpus()
    # s0 shifts hard one way; opposing studies the other way.
    counts = _counts(pairs, 1, shift={"s0": 40})
    results = composition_contrasts(
        pair_metadata(_profiles(pairs)), counts, TABLE1_COLUMNS, min_site_pairs=5
    )
    site = _site_rows(results)["facultatively_aerobic"]
    assert site.status == TESTED
    assert site.studies == 10
    assert site.paired_experiments == len(pairs)
    assert site.conditions == 3

    # Relabelling s0's conditions as one changes nothing at site level: the site
    # value is one mean per study, whatever conditions the study spans.
    merged = [(s, e, "A" if s == "s0" else c, b) for s, e, c, b in pairs]
    again = composition_contrasts(
        pair_metadata(_profiles(merged)), counts, TABLE1_COLUMNS, min_site_pairs=5
    )
    for trait, result in _site_rows(results).items():
        other = _site_rows(again)[trait]
        assert (result.mean_contrast_z, result.p_value, result.q_value, result.family_size) == (
            other.mean_contrast_z,
            other.p_value,
            other.q_value,
            other.family_size,
        )


def test_conditions_keep_their_own_contrasts() -> None:
    pairs = _site_corpus()
    results = composition_contrasts(
        pair_metadata(_profiles(pairs)), _counts(pairs, 2), TABLE1_COLUMNS, min_site_pairs=5
    )
    cells = [
        r for r in results if r.evidence_type == CELL_COMPOSITION and r.trait == "coccus_shaped"
    ]
    by_condition = {r.condition: r for r in cells}
    # A: s0 plus s2, s4, s6, s8. B: s0 plus s1, s3, s5, s7, s9. C: s0 alone.
    assert by_condition["A"].studies == 5
    assert by_condition["B"].studies == 6
    assert by_condition["C"].status == "condition_below_min_studies"
    assert by_condition["C"].q_value is None
    assert {r.family for r in cells if r.q_value is not None} == {"all_condition_site_cells"}


def test_untestable_members_stay_in_the_family_at_p_one() -> None:
    """The arm fills an untestable trait's p with 1 before correcting; so does this."""
    base = CompositionResult(SITE_COMPOSITION, "Feces", None, "t", "t", TESTED, 3, 3, 1)
    results = [
        replace(base, trait="a", column="a", p_value=0.01),
        replace(base, trait="b", column="b", p_value=0.04),
        replace(base, trait="c", column="c", status="all_tied"),
    ]
    corrected = _correct(results, lambda _result: "site:Feces")
    assert [r.family_size for r in corrected] == [3, 3, 3]
    assert corrected[0].q_value == pytest.approx(0.03)
    assert corrected[2].q_value is None


def _enrichment_file(path: Path, rows: list[tuple[str, str, str, str, bool]]) -> Path:
    header = [
        "study_id",
        "experiment_id",
        "inc_dec",
        "trait",
        "p_value",
        "fdr_q",
        "fold",
        "count",
        "expected",
        "testable_profiles",
        "background_m",
        "background_k",
        "is_significant",
        "correction_family",
        "eligibility_status",
    ]
    hits = sum(flag for *_, flag in rows)
    lines = [
        "# Tail: upper",
        "# FDR: True",
        "# FDR-Scope: global",
        "# P-Threshold: 0.001",
        "# Min-Fold: 1.5",
        f"# Total-Tests-Performed: {len(rows)}",
        f"# Significant-Hits: {hits}",
        f"# Result-Status: {'discoveries' if hits else 'zero_discoveries'}",
        "\t".join(header),
    ]
    for study, experiment, side, column, flag in rows:
        q = "0.0001" if flag else "0.5"
        stats = [q, q, "2.0", "3", "1.5", "10", "100", "20"]
        extra = [str(flag), "global", "testable"]
        lines.append("\t".join([study, experiment, side, column, *stats, *extra]))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


FAC = TABLE1_COLUMNS["facultatively_aerobic"]


def _calls(tmp_path: Path, rows: list[tuple[str, str, str, str, bool]]) -> EnrichmentCalls:
    return read_enrichment_calls(_enrichment_file(tmp_path / "all.tsv", rows), TABLE1_COLUMNS)


def test_one_study_never_speaks_for_a_cell(tmp_path: Path) -> None:
    pairs: list[Pair] = [("s1", "e1", "X", "Feces")]
    calls = _calls(tmp_path, [("s1", "e1", "inc", FAC, True)])
    record = aggregate_support(_profiles(pairs), calls, ["facultatively_aerobic"])[
        ("X", "Feces", "facultatively_aerobic", "inc")
    ]
    assert (record.eligible_studies, record.supporting_studies) == (1, 1)
    assert record.status == INSUFFICIENT_STUDIES


def test_support_counts_unique_studies_and_the_alternative_removes_extra_chances(
    tmp_path: Path,
) -> None:
    """#177: a study with three experiments and one call supports under the any-profile
    rule but not under the per-study majority rule."""
    pairs: list[Pair] = [
        ("s1", "e1", "X", "Feces"),
        ("s1", "e2", "X", "Feces"),
        ("s1", "e3", "X", "Feces"),
        ("s2", "e1", "X", "Feces"),
        ("s3", "e1", "X", "Feces"),
    ]
    calls = _calls(
        tmp_path,
        [
            ("s1", "e1", "inc", FAC, True),
            ("s1", "e2", "inc", FAC, False),
            ("s1", "e3", "inc", FAC, False),
            ("s2", "e1", "inc", FAC, True),
            ("s3", "e1", "inc", FAC, False),
        ],
    )
    record = aggregate_support(_profiles(pairs), calls, ["facultatively_aerobic"])[
        ("X", "Feces", "facultatively_aerobic", "inc")
    ]
    assert record.eligible_studies == 3
    assert record.supporting_studies == 2
    assert record.alternative_supporting_studies == 1
    assert record.status == MAJORITY_SUPPORT
    assert record.alternative_status == NO_MAJORITY
    assert (record.tested_profiles, record.called_profiles) == (5, 2)
    assert record.studies_with_repeated_profiles == 1
    # The dec side has no result rows at all: untested, never a miss.
    dec = aggregate_support(_profiles(pairs), calls, ["facultatively_aerobic"])[
        ("X", "Feces", "facultatively_aerobic", "dec")
    ]
    assert (dec.eligible_studies, dec.untested_profiles) == (0, 5)


def test_only_an_fdr_corrected_enrichment_run_is_accepted(tmp_path: Path) -> None:
    path = _enrichment_file(tmp_path / "all.tsv", [("s1", "e1", "inc", FAC, True)])
    path.write_text(
        path.read_text(encoding="utf-8").replace("# Tail: upper", "# Tail: lower"),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="upper-tail"):
        read_enrichment_calls(path, TABLE1_COLUMNS)


def _matrix() -> GuildMatrix:
    columns = [
        f"{site}|{trait}_{side}"
        for site in ("Feces", OTHER_LABEL)
        for trait in TRAIT_ORDER
        for side in ("inc", "dec")
    ]
    width = len(columns)
    values: list[list[float | None]] = [
        [0.1] * width,
        [None] * 2 + [0.1] * (width - 2),
    ]
    return GuildMatrix(
        conditions=["Colorectal cancer", "Asthma"],
        body_sites=["Feces", OTHER_LABEL],
        columns=columns,
        values=values,
        study_counts=[[5] * width, [1] * 2 + [5] * (width - 2)],
        row_clusters=[1, 1],
        pooled=[],
        dropped_empty_conditions=[],
        dropped_body_sites=[],
        sampled_no_signal=[[False] * width, [False] * width],
    )


def _marks(
    cells: dict[tuple[str, str, str, str], frozenset[str]],
    columns: dict[tuple[str, str], str] | None = None,
    totals: dict[str, int] | None = None,
) -> EvidenceMarks:
    styles = {SUPPORT: MarkStyle("solid", 0.0), CELL_COMPOSITION: MarkStyle("dot", 0.18)}
    return EvidenceMarks(
        cells=cells,
        styles=styles,
        legend={kind: f"<{kind}>" for kind in styles},
        corpus_totals=totals or {kind: sum(kind in k for k in cells.values()) for kind in styles},
        columns=columns or {},
        column_legend="<glyph>",
    )


def _outlines(figure: Figure) -> list[str]:
    """The dash style of every evidence outline the figure draws."""
    return [s.line.dash for s in figure.layout.shapes if s.line.color == MARK_LINE_COLOR]


def test_legend_counts_only_the_marks_a_figure_draws() -> None:
    """#216: a mark on a masked cell is neither drawn nor counted as drawn."""
    matrix = _matrix()
    drawn = ("Colorectal cancer", "Feces", "facultatively_aerobic", "inc")
    masked = ("Asthma", "Feces", "facultatively_aerobic", "inc")
    assert cell_display_status(matrix, drawn) == "drawn"
    assert cell_display_status(matrix, masked) == "below_min_studies"
    marks = _marks({drawn: frozenset({SUPPORT}), masked: frozenset({SUPPORT})})
    figure = build_figure(matrix, None, guild_order=TRAIT_ORDER, colorscales=SCALES, evidence=marks)
    assert len(_outlines(figure)) == 1
    subtitle = str(figure.layout.title.text).replace("<br>", " ")
    assert f"<{SUPPORT}> Drawn here: 1 cells. Across the corpus 2 cells meet the rule" in subtitle


def test_every_mark_type_is_explained_even_when_none_is_drawn() -> None:
    figure = build_figure(
        _matrix(), None, guild_order=TRAIT_ORDER, colorscales=SCALES, evidence=_marks({})
    )
    subtitle = str(figure.layout.title.text).replace("<br>", " ")
    assert f"<{SUPPORT}> Drawn here: 0 cells" in subtitle
    assert f"<{CELL_COMPOSITION}> Drawn here: 0 cells" in subtitle


def test_an_unknown_mark_type_is_an_error_not_a_default_style() -> None:
    cell = ("Colorectal cancer", "Feces", "coccus_shaped", "inc")
    with pytest.raises(KeyError):
        build_figure(
            _matrix(),
            None,
            guild_order=TRAIT_ORDER,
            colorscales=SCALES,
            evidence=_marks({cell: frozenset({"both"})}),
        )


def test_two_evidence_types_on_one_cell_are_both_visible() -> None:
    cell = ("Colorectal cancer", "Feces", "coccus_shaped", "inc")
    figure = build_figure(
        _matrix(),
        None,
        guild_order=TRAIT_ORDER,
        colorscales=SCALES,
        evidence=_marks({cell: frozenset({SUPPORT, CELL_COMPOSITION})}),
    )
    assert sorted(_outlines(figure)) == ["dot", "solid"]


def test_the_pooled_other_block_is_never_marked() -> None:
    cell = ("Colorectal cancer", OTHER_LABEL, "coccus_shaped", "inc")
    assert cell_display_status(_matrix(), cell) == "site_not_drawn"
    figure = build_figure(
        _matrix(),
        None,
        guild_order=TRAIT_ORDER,
        colorscales=SCALES,
        evidence=_marks(
            {cell: frozenset({SUPPORT})}, columns={(OTHER_LABEL, "coccus_shaped"): "increased"}
        ),
    )
    assert not _outlines(figure)
    assert not [a for a in figure.layout.annotations if a.text in COLUMN_MARK_GLYPH.values()]


def test_a_site_glyph_is_drawn_once_per_trait_over_both_columns() -> None:
    """#217: the glyph is one statement about the trait, centred on its two columns."""
    marks = _marks({}, columns={("Feces", "coccus_shaped"): "decreased"})
    for figure in (
        build_figure(_matrix(), None, guild_order=TRAIT_ORDER, colorscales=SCALES, evidence=marks),
        build_body_site_panels(
            _matrix(),
            None,
            sites=["Feces"],
            guild_order=TRAIT_ORDER,
            colorscales=SCALES,
            evidence=marks,
        ),
    ):
        glyphs = [a for a in figure.layout.annotations if a.text in COLUMN_MARK_GLYPH.values()]
        assert len(glyphs) == 1
        assert glyphs[0].text == COLUMN_MARK_GLYPH["decreased"]
        # coccus shaped is the third trait: columns 4 and 5 of the Feces block.
        assert glyphs[0].x == 4.5
        assert glyphs[0].textangle is None


def test_every_trait_keeps_its_colour_bar_when_the_last_panel_lacks_it() -> None:
    """#217: a trait with no cell in the last panel lost its bar."""
    matrix = _matrix()
    missing = {i for i, c in enumerate(matrix.columns) if c.startswith(f"{OTHER_LABEL}|arsenate")}
    matrix = replace(
        matrix,
        values=[[None if i in missing else v for i, v in enumerate(row)] for row in matrix.values],
    )
    figure = build_body_site_panels(
        matrix,
        None,
        sites=["Feces", OTHER_LABEL],
        guild_order=TRAIT_ORDER,
        colorscales=SCALES,
    )
    # Plotly drops the colour bar of a trace with no values, so each bar must sit on
    # a trace that draws at least one cell.
    shown = [
        trace.name
        for trace in figure.data
        if trace.showscale and any(v is not None for row in trace.z for v in row)
    ]
    assert sorted(shown) == sorted(TRAIT_ORDER)


def test_site_results_are_never_promoted_to_cells() -> None:
    site = CompositionResult(
        SITE_COMPOSITION,
        "Feces",
        None,
        "coccus_shaped",
        FAC,
        TESTED,
        9,
        9,
        3,
        mean_contrast_z=0.3,
        p_value=0.001,
        q_value=0.001,
    )
    evidence = Table1Evidence(traits=TRAIT_ORDER, composition=[site], support=None, enrichment=None)
    assert evidence.cell_marks() == {}
    assert evidence.site_marks() == {("Feces", "coccus_shaped"): "increased"}
    marks = evidence_marks(evidence)
    assert SUPPORT not in marks.legend
    assert any("profile-level support is not drawn" in note for note in marks.notes)


def test_a_condition_contrast_marks_only_the_side_it_favours() -> None:
    cell = CompositionResult(
        CELL_COMPOSITION,
        "Feces",
        "X",
        "coccus_shaped",
        FAC,
        TESTED,
        4,
        4,
        1,
        mean_contrast_z=-0.4,
        p_value=0.001,
        q_value=0.01,
    )
    evidence = Table1Evidence(traits=TRAIT_ORDER, composition=[cell], support=None, enrichment=None)
    assert evidence.cell_marks() == {
        ("X", "Feces", "coccus_shaped", "dec"): frozenset({CELL_COMPOSITION})
    }


def test_support_is_labelled_descriptive_with_no_cell_p_or_q(tmp_path: Path) -> None:
    """#177: the legend says what the rule is, and the table has no p or q column."""
    pairs: list[Pair] = [("s1", "e1", "X", "Feces"), ("s2", "e1", "X", "Feces")]
    calls = _calls(tmp_path, [("s1", "e1", "inc", FAC, True), ("s2", "e1", "inc", FAC, True)])
    evidence = Table1Evidence(
        traits=TRAIT_ORDER,
        composition=[],
        support=aggregate_support(_profiles(pairs), calls, list(TRAIT_ORDER)),
        enrichment=calls,
    )
    marks = evidence_marks(evidence)
    assert "A descriptive support rule with no cell-level p- or q-value." in marks.legend[SUPPORT]
    assert marks.corpus_totals[SUPPORT] == 1
    assert not [c for c in SUPPORT_COLUMNS if re.search(r"(^|_)(p|q)(_|$)", c)]


def _methods_rows() -> set[tuple[str, str, str]]:
    rows = set()
    for line in METHODS_DOC.read_text(encoding="utf-8").splitlines():
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        if len(cells) == 3 and line.startswith("|"):
            rows.add((cells[0], cells[1], cells[2]))
    return rows


def test_departures_are_one_constant_rendered_in_figure_and_methods() -> None:
    """#215: the figure note and the methods table both come from the same rows."""
    evidence = Table1Evidence(traits=TRAIT_ORDER, composition=[], support=None, enrichment=None)
    note = evidence_marks(evidence).notes[0]
    rows = _methods_rows()
    for aspect, arm, ours in COMPOSITION_DEPARTURES:
        assert f"{aspect}: {ours}" in note
        assert (aspect, arm.replace("|", "\\|"), ours.replace("|", "\\|")) in rows


def _run_main(
    corpus: tuple[Path, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch, *extra: str
) -> tuple[Path, list[EvidenceMarks]]:
    marks: list[EvidenceMarks] = []

    def record(*_args: object, **kwargs: object) -> Path:
        evidence = kwargs["evidence"]
        assert isinstance(evidence, EvidenceMarks)
        marks.append(evidence)
        return tmp_path / "figure.png"

    monkeypatch.setattr(table1_trait_maps, "write_figure", record)
    monkeypatch.setattr(table1_trait_maps, "write_body_site_panels", record)
    monkeypatch.setattr(table1_trait_maps, "write_shift_figure", lambda *_a, **_k: tmp_path / "s")
    monkeypatch.setattr(table1_trait_maps, "read_bugsigdb_vintage", lambda _path: "fixture")
    output = tmp_path / f"out{len(list(tmp_path.glob('out*')))}"
    monkeypatch.setattr(
        "sys.argv",
        [
            "table1-trait-maps",
            "--profile-input",
            str(corpus[0]),
            "--experiments-input",
            str(corpus[1]),
            "--output-dir",
            str(output),
            "--min-studies",
            "1",
            "--min-site-conditions",
            "1",
            *extra,
        ],
    )
    assert table1_trait_maps.main() == 0
    return output, marks


@pytest.fixture
def corpus(tmp_path: Path) -> tuple[Path, Path]:
    """Sixty Feces pairs over twelve studies and four conditions, plus some Saliva."""
    rng = np.random.default_rng(7)
    profiles = tmp_path / "profiles.tsv"
    experiments = tmp_path / "experiments.tsv"
    header = ["study_id", "experiment_id", "inc_dec", *COLUMNS]
    lines = ["\t".join(header)]
    meta = ["\t".join(["study_id", "experiment_id", "condition", "body_site", "host_species"])]
    for index in range(70):
        study = f"bsdb:{index % 12}"
        experiment = f"{study}-{index}"
        site = "Feces" if index < 60 else "Saliva"
        meta.append("\t".join([study, experiment, f"C{index % 4}", site, "Homo sapiens"]))
        for side in ("inc", "dec"):
            counts = rng.integers(1, 9, size=len(COLUMNS))
            counts[rng.random(len(COLUMNS)) < 0.15] = 0
            if side == "inc":
                counts[0] += 5
            lines.append("\t".join([study, experiment, side, *map(str, counts)]))
    profiles.write_text("\n".join(lines) + "\n", encoding="utf-8")
    experiments.write_text("\n".join(meta) + "\n", encoding="utf-8")
    return profiles, experiments


def _evidence_rows(path: Path) -> list[list[str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    drop = rows[0].index("display_status")
    return [row[:drop] + row[drop + 1 :] for row in rows]


def test_evidence_does_not_depend_on_the_display_cut(
    corpus: tuple[Path, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """#176: narrowing the drawn rows and columns moves marks, never results."""
    wide, wide_marks = _run_main(corpus, tmp_path, monkeypatch)
    narrow, narrow_marks = _run_main(
        corpus, tmp_path, monkeypatch, "--max-conditions", "1", "--max-body-sites", "1"
    )
    for name in ("composition_contrasts", "composition_screen"):
        left = _evidence_rows(wide / f"table1_traits_{name}.tsv")
        right = _evidence_rows(narrow / f"table1_traits_{name}.tsv")
        assert left == right
    feces = [row for row in left if row[1] == "Feces" and row[0] == SITE_COMPOSITION]
    assert feces
    assert all(row[6] == "12" for row in feces)
    assert wide_marks[0].cells == narrow_marks[0].cells
    assert wide_marks[0].corpus_totals == narrow_marks[0].corpus_totals
    assert wide_marks[0].columns == narrow_marks[0].columns
