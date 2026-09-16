"""Table 1 trait selection, the paired directional shift, and its FDR control."""

from pathlib import Path

import pytest
from scipy.stats import binomtest, false_discovery_control

from bsdbng.guild_maps import (
    COLORBAR_TITLE_CHARS,
    COLORBAR_TOP_TICK_SHARE,
    COLUMN_MARK_GLYPH,
    MARK_STYLES,
    OTHER_LABEL,
    SIGN_TOLERANCE,
    EvidenceMarks,
    GuildMatrix,
    body_site_panels_figure,
    build_matrix,
    condition_body_site_figure,
    subset_matrix,
)
from bsdbng.table1_trait_maps import (
    TRAIT_ORDER,
    FocusPanel,
    MissingConditions,
    cell_display_status,
    compute_direction_shifts,
    evidence_marks,
    missing_condition_note,
    missing_conditions,
    trait_colorscale,
    uncovered_conditions,
)
from bsdbng.trait_guilds import (
    TABLE1_TRAIT_COLORS,
    TABLE1_TRAIT_LABELS,
    TABLE1_TRAITS,
    Profile,
    column_term_id,
    load_profiles,
)
from bsdbng.trait_significance import CELL_COMPOSITION, SUPPORT, build_overlay

# Two of the six traits reach the profile matrix behind a METPO qualifier, so the
# fixture uses the qualified spelling to prove the suffix lookup is what resolves
# them. The other four are bare CURIEs, as in the real header.
QUALIFIED_CELLULOSE = "METPO:2000103|GO:0030245"
QUALIFIED_ARSENATE = "METPO:2000103|pathways:arsenate_reduction"
TRAIT_COLUMNS = [
    "METPO:1000608",
    QUALIFIED_CELLULOSE,
    "METPO:1000668",
    "METPO:1000683",
    QUALIFIED_ARSENATE,
    "METPO:1000623",
]
OTHER_COLUMN = "METPO:1000615"
MEDIA_COLUMN = "METPO:2000517|mediadive.medium:J14"


def _write_profiles(path: Path, rows: list[tuple[str, str, str, list[int]]]) -> None:
    header = ["study_id", "experiment_id", "inc_dec", *TRAIT_COLUMNS, OTHER_COLUMN, MEDIA_COLUMN]
    lines = ["\t".join(header)]
    for study_id, experiment_id, inc_dec, values in rows:
        lines.append("\t".join([study_id, experiment_id, inc_dec, *[str(v) for v in values]]))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_experiments(path: Path, rows: list[tuple[str, str, str, str, str]]) -> None:
    header = ["study_id", "experiment_id", "condition", "body_site", "host_species"]
    lines = ["\t".join(header)]
    for row in rows:
        lines.append("\t".join(row))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


@pytest.fixture
def corpus(tmp_path: Path) -> tuple[Path, Path]:
    """Three paired human studies with a deliberate per-trait direction pattern.

    Cells are counts and the figure works in fractions, so the balancing
    non-guild column is set to give each study the same total on both sides.
    Otherwise a trait with identical counts either side would still register a
    shift, purely because the two denominators differ. That gives, per trait:
    facultatively aerobic up in all three studies, cellulose catabolic process
    down in all three, coccus shaped two up and one down, sphere shaped tied in
    every study, arsenate reduction one up and two ties, moderately halophilic
    one up, one down, and one tie.
    """
    profiles = tmp_path / "profiles.tsv"
    experiments = tmp_path / "experiments.tsv"
    _write_profiles(
        profiles,
        [
            # fac_aer, cellulose, coccus, sphere, arsenate, halophilic, other, media
            ("bsdb:1", "bsdb:1-1", "inc", [8, 1, 6, 2, 0, 3, 5, 100]),
            ("bsdb:1", "bsdb:1-1", "dec", [2, 7, 2, 2, 0, 3, 9, 100]),
            ("bsdb:2", "bsdb:2-1", "inc", [9, 0, 7, 1, 4, 2, 4, 100]),
            ("bsdb:2", "bsdb:2-1", "dec", [3, 6, 3, 1, 1, 5, 8, 100]),
            ("bsdb:3", "bsdb:3-1", "inc", [7, 2, 1, 3, 2, 4, 6, 100]),
            ("bsdb:3", "bsdb:3-1", "dec", [1, 9, 5, 3, 2, 1, 4, 100]),
            # Unpaired and non-human rows, dropped by the shared corpus filters.
            ("bsdb:4", "bsdb:4-1", "inc", [5, 5, 5, 5, 5, 5, 5, 100]),
            ("bsdb:5", "bsdb:5-1", "inc", [9, 0, 9, 0, 9, 0, 1, 100]),
            ("bsdb:5", "bsdb:5-1", "dec", [0, 9, 0, 9, 0, 9, 1, 100]),
        ],
    )
    _write_experiments(
        experiments,
        [
            ("bsdb:1", "bsdb:1-1", "Colorectal cancer", "Feces", "Homo sapiens"),
            ("bsdb:2", "bsdb:2-1", "Colorectal cancer", "Feces", "Homo sapiens"),
            ("bsdb:3", "bsdb:3-1", "Asthma", "Feces", "Homo sapiens"),
            ("bsdb:4", "bsdb:4-1", "Asthma", "Feces", "Homo sapiens"),
            ("bsdb:5", "bsdb:5-1", "Colitis", "Caecum", "Mus musculus"),
        ],
    )
    return profiles, experiments


def test_table1_definitions_cover_the_six_manuscript_traits() -> None:
    assert TRAIT_ORDER == (
        "facultatively_aerobic",
        "cellulose_catabolic_process",
        "coccus_shaped",
        "sphere_shaped",
        "arsenate_reduction",
        "moderately_halophilic",
    )
    assert set(TABLE1_TRAIT_LABELS) == set(TRAIT_ORDER)
    assert set(TABLE1_TRAIT_COLORS) == set(TRAIT_ORDER)
    # Every trait is a single term: these are traits, not pooled guilds.
    assert all(len(ids) == 1 for ids in TABLE1_TRAITS.values())
    assert trait_colorscale("coccus_shaped") == [[0.0, "#ffffff"], [1.0, "#009e73"]]


def test_qualified_columns_resolve_by_term_suffix(corpus: tuple[Path, Path]) -> None:
    """GO:0030245 reaches the matrix as METPO:2000103|GO:0030245 and must still match."""
    assert column_term_id(QUALIFIED_CELLULOSE) == "GO:0030245"
    assert column_term_id("METPO:1000608") == "METPO:1000608"

    profiles, stats = load_profiles(*corpus, TABLE1_TRAITS)
    inc = next(p for p in profiles if p.study_id == "bsdb:2" and p.inc_dec == "inc")

    assert inc.guild_counts["cellulose_catabolic_process"] == 0
    assert inc.guild_counts["arsenate_reduction"] == 4
    assert inc.guild_counts["facultatively_aerobic"] == 9
    # The MediaDive column is excluded from the denominator, the other trait is not.
    assert inc.total == 9 + 0 + 7 + 1 + 4 + 2 + 4
    assert inc.fraction("arsenate_reduction") == pytest.approx(4 / 27)
    assert stats.mediadive_columns_excluded == 1
    assert stats.dropped_non_human == 2
    assert stats.dropped_unpaired == 1


def test_ambiguous_term_lookup_raises(tmp_path: Path) -> None:
    """A term matching two columns would be summed twice; that must fail loudly."""
    profiles = tmp_path / "profiles.tsv"
    experiments = tmp_path / "experiments.tsv"
    header = ["study_id", "experiment_id", "inc_dec", "METPO:1000608", "METPO:9999|METPO:1000608"]
    profiles.write_text(
        "\t".join(header) + "\nbsdb:1\tbsdb:1-1\tinc\t3\t3\n",
        encoding="utf-8",
    )
    _write_experiments(experiments, [("bsdb:1", "bsdb:1-1", "Asthma", "Feces", "Homo sapiens")])

    with pytest.raises(ValueError, match="more than one profile column"):
        load_profiles(profiles, experiments, {"facultatively_aerobic": ("METPO:1000608",)})


def test_matrix_unrolls_six_traits_per_body_site(corpus: tuple[Path, Path]) -> None:
    profiles, _ = load_profiles(*corpus, TABLE1_TRAITS)
    matrix = build_matrix(
        profiles,
        max_conditions=10,
        max_body_sites=10,
        min_studies=1,
        guild_order=TRAIT_ORDER,
    )

    # Two body sites survive (Feces and the Other bucket), each contributing six
    # traits times two directions.
    assert len(matrix.columns) == len(matrix.body_sites) * len(TRAIT_ORDER) * 2
    assert "Feces|facultatively_aerobic_inc" in matrix.columns
    assert "Feces|cellulose_catabolic_process_dec" in matrix.columns


def test_direction_shift_counts_studies_and_matches_scipy(corpus: tuple[Path, Path]) -> None:
    profiles, _ = load_profiles(*corpus, TABLE1_TRAITS)
    result = compute_direction_shifts(profiles, TRAIT_ORDER)
    by_trait = {shift.trait: shift for shift in result.shifts}

    assert result.studies == 3

    aerobic = by_trait["facultatively_aerobic"]
    assert (aerobic.studies_up, aerobic.studies_down, aerobic.studies_tied) == (3, 0, 0)
    assert aerobic.share_up == pytest.approx(1.0)

    cellulose = by_trait["cellulose_catabolic_process"]
    assert (cellulose.studies_up, cellulose.studies_down) == (0, 3)
    assert cellulose.share_up == pytest.approx(0.0)

    coccus = by_trait["coccus_shaped"]
    assert (coccus.studies_up, coccus.studies_down) == (2, 1)

    # Sphere shaped holds the same fraction on both sides of every experiment, so
    # every study ties and none is scored; the trait must not silently vanish.
    sphere = by_trait["sphere_shaped"]
    assert (sphere.studies_up, sphere.studies_down, sphere.studies_tied) == (0, 0, 3)
    assert sphere.studies_scored == 0
    assert sphere.p_value == 1.0

    # Ties are excluded from the test rather than split between the two sides.
    arsenate = by_trait["arsenate_reduction"]
    assert (arsenate.studies_up, arsenate.studies_down, arsenate.studies_tied) == (1, 0, 2)
    assert arsenate.studies_scored == 1
    halophilic = by_trait["moderately_halophilic"]
    assert (halophilic.studies_up, halophilic.studies_down, halophilic.studies_tied) == (1, 1, 1)

    # Golden values: the p-values must equal scipy's exact binomial test, and the
    # q-values scipy's Benjamini-Hochberg over the same six p-values in order.
    expected_p = [
        float(binomtest(by_trait[t].studies_up, by_trait[t].studies_scored, p=0.5).pvalue)
        if by_trait[t].studies_scored
        else 1.0
        for t in TRAIT_ORDER
    ]
    expected_q = [float(q) for q in false_discovery_control(expected_p, method="bh")]
    assert [by_trait[t].p_value for t in TRAIT_ORDER] == pytest.approx(expected_p)
    assert [by_trait[t].q_value for t in TRAIT_ORDER] == pytest.approx(expected_q)
    assert aerobic.p_value == pytest.approx(0.25)


def test_contrast_correlation_is_computed_across_studies(corpus: tuple[Path, Path]) -> None:
    """Trait pairs correlate over per-study contrasts, not within a profile total."""
    profiles, _ = load_profiles(*corpus, TABLE1_TRAITS)
    result = compute_direction_shifts(profiles, TRAIT_ORDER)

    key = "facultatively_aerobic~cellulose_catabolic_process"
    assert key in result.spearman
    rho, _, n_studies = result.spearman[key]
    assert n_studies == 3
    assert -1.0 <= rho <= 1.0

    # Sphere shaped has the same contrast in every study, so Spearman is
    # undefined against it; the pair is omitted rather than written as nan.
    assert not [k for k in result.spearman if "sphere_shaped" in k]


def test_panel_figure_rejects_a_body_site_it_cannot_draw(corpus: tuple[Path, Path]) -> None:
    """Naming an absent site would otherwise render a silently empty panel."""
    profiles, _ = load_profiles(*corpus, TABLE1_TRAITS)
    matrix = build_matrix(
        profiles,
        max_conditions=10,
        max_body_sites=10,
        min_studies=1,
        guild_order=TRAIT_ORDER,
    )

    with pytest.raises(ValueError, match="Nowhere"):
        body_site_panels_figure(matrix, None, sites=["Nowhere"], guild_order=TRAIT_ORDER)


def _paired_profile(
    study: str, experiment: str, inc_dec: str, count: int, total: int = 10
) -> Profile:
    counts = dict.fromkeys(TRAIT_ORDER, 0)
    counts["facultatively_aerobic"] = count
    return Profile(
        study_id=study,
        experiment_id=experiment,
        inc_dec=inc_dec,
        condition="X",
        body_site="Feces",
        guild_counts=counts,
        total=total,
    )


def test_studies_are_scored_once_not_once_per_experiment() -> None:
    """The #133 collapse itself: a multi-experiment study contributes ONE sign.

    The two experiments disagree (+0.8 and -0.2), so per-experiment scoring
    would tally 1 up + 1 down; the study-level collapse averages first
    (mean inc 0.55, mean dec 0.25) and tallies a single up.
    """
    study = [
        _paired_profile("bsdb:1", "bsdb:1-1", "inc", 9),
        _paired_profile("bsdb:1", "bsdb:1-1", "dec", 1),
        _paired_profile("bsdb:1", "bsdb:1-2", "inc", 2),
        _paired_profile("bsdb:1", "bsdb:1-2", "dec", 4),
    ]
    result = compute_direction_shifts(study, ("facultatively_aerobic",))
    (shift,) = result.shifts
    assert (shift.studies_up, shift.studies_down, shift.studies_tied) == (1, 0, 0)
    assert result.studies == 1


def test_float_residue_contrasts_are_ties_not_signs() -> None:
    """A contrast below SIGN_TOLERANCE is summation noise, not a direction."""
    residue = SIGN_TOLERANCE / 10
    study = [
        _paired_profile("bsdb:1", "bsdb:1-1", "inc", 1, total=int(1 / residue)),
        _paired_profile("bsdb:1", "bsdb:1-1", "dec", 0, total=10),
    ]
    result = compute_direction_shifts(study, ("facultatively_aerobic",))
    (shift,) = result.shifts
    assert (shift.studies_up, shift.studies_down, shift.studies_tied) == (0, 0, 1)
    assert shift.p_value == 1.0


def test_undrawn_panel_conditions_separate_all_three_causes() -> None:
    missing = missing_conditions(
        ("Gingivitis", "Periodontitis", "Rheumatoid arthritis", "Never studied"),
        corpus_conditions={"Gingivitis", "Periodontitis", "Rheumatoid arthritis"},
        drawn_conditions=["Rheumatoid arthritis"],
        blank_conditions=["Gingivitis"],
    )
    assert missing.absent_from_corpus == ("Never studied",)
    assert missing.blank_after_masking == ("Gingivitis",)
    assert missing.cut_by_display == ("Periodontitis",)


def test_a_cut_condition_is_never_reported_as_missing_from_the_export() -> None:
    # v1.3.1 holds periodontitis; only the top-N cut drops it. Reporting that as
    # an absent condition states something false about the export, which is the
    # defect this pair of helpers exists to prevent.
    note = missing_condition_note(
        missing_conditions(
            ("Periodontitis",),
            corpus_conditions={"Periodontitis"},
            drawn_conditions=[],
        ),
        max_conditions=60,
        min_studies=3,
    )
    assert note is not None
    assert "does not contain" not in note
    assert "outside the top-60 condition cut, so not drawn: Periodontitis" in note


def test_a_blank_row_is_not_blamed_on_the_width_of_the_cut() -> None:
    # Gingivitis is inside the cut at this width and still undrawn, because every
    # cell is below min_studies. Calling that a top-N exclusion would be false --
    # at max_conditions=451 nothing is outside the cut at all.
    note = missing_condition_note(
        missing_conditions(
            ("Gingivitis",),
            corpus_conditions={"Gingivitis"},
            drawn_conditions=[],
            blank_conditions=["Gingivitis"],
        ),
        max_conditions=451,
        min_studies=3,
    )
    assert note is not None
    assert "outside the top-451" not in note
    assert "fewer than 3 studies" in note
    assert "Gingivitis" in note


def test_the_note_states_every_cause_in_its_own_words() -> None:
    note = missing_condition_note(
        MissingConditions(
            absent_from_corpus=("Never studied",),
            blank_after_masking=("Gingivitis",),
            cut_by_display=("Periodontitis",),
        ),
        max_conditions=60,
        min_studies=3,
    )
    assert note is not None
    assert "does not contain them: Never studied" in note
    assert "fewer than 3 studies, so the row is blank and not drawn: Gingivitis" in note
    assert "outside the top-60 condition cut, so not drawn: Periodontitis" in note


def test_no_note_when_every_requested_condition_is_drawn() -> None:
    assert (
        missing_condition_note(
            missing_conditions(
                ("Rheumatoid arthritis",),
                corpus_conditions={"Rheumatoid arthritis"},
                drawn_conditions=["Rheumatoid arthritis"],
            ),
            max_conditions=60,
            min_studies=3,
        )
        is None
    )


CellMarks = dict[tuple[str, str, str, str], frozenset[str]]
SCALES: dict[str, object] = {trait: trait_colorscale(trait) for trait in TRAIT_ORDER}


def _matrix(
    corpus: tuple[Path, Path], *, max_conditions: int = 10, min_studies: int = 1
) -> tuple[list[Profile], GuildMatrix]:
    profiles, _ = load_profiles(*corpus, TABLE1_TRAITS)
    matrix = build_matrix(
        profiles,
        max_conditions=max_conditions,
        max_body_sites=10,
        min_studies=min_studies,
        guild_order=TRAIT_ORDER,
        include_other_conditions=False,
        mark_no_signal=False,
    )
    return profiles, matrix


def _evidence(cells: CellMarks, columns: dict[tuple[str, str], str] | None = None) -> EvidenceMarks:
    totals: dict[str, int] = {}
    for kinds in cells.values():
        for kind in kinds:
            totals[kind] = totals.get(kind, 0) + 1
    return EvidenceMarks(
        cells=cells,
        totals=totals,
        legend={kind: f"<{kind}>" for kind in MARK_STYLES},
        columns=columns or {},
        column_legend="<glyph>",
    )


def test_legend_counts_only_the_marks_a_figure_draws(corpus: tuple[Path, Path]) -> None:
    """#216: a mark on a masked cell is counted as not drawn, never as drawn."""
    _, matrix = _matrix(corpus, min_studies=2)
    drawn_cell = ("Colorectal cancer", "Feces", "facultatively_aerobic", "inc")
    masked_cell = ("Asthma", "Feces", "facultatively_aerobic", "inc")
    assert cell_display_status(matrix, drawn_cell) == "drawn"
    assert cell_display_status(matrix, masked_cell) == "condition_not_drawn"
    evidence = _evidence({drawn_cell: frozenset({SUPPORT}), masked_cell: frozenset({SUPPORT})})

    figure = condition_body_site_figure(
        matrix, None, guild_order=TRAIT_ORDER, colorscales=SCALES, evidence=evidence
    )
    outlines = [shape for shape in figure.layout.shapes if shape.line.color == "#111111"]
    assert len(outlines) == 1
    subtitle = str(figure.layout.title.text)
    assert "Drawn here: 1 cells. Across the corpus 2 cells meet the rule" in subtitle


def test_every_mark_type_present_is_explained_even_when_none_is_drawn(
    corpus: tuple[Path, Path],
) -> None:
    _, matrix = _matrix(corpus, min_studies=2)
    undrawn = ("Asthma", "Feces", "coccus_shaped", "dec")
    evidence = _evidence({undrawn: frozenset({CELL_COMPOSITION})})
    subtitle = str(
        condition_body_site_figure(
            matrix, None, guild_order=TRAIT_ORDER, colorscales=SCALES, evidence=evidence
        ).layout.title.text
    )
    assert f"<{CELL_COMPOSITION}> Drawn here: 0 cells" in subtitle
    assert f"<{SUPPORT}>" not in subtitle


def test_two_evidence_types_on_one_cell_are_both_visible(corpus: tuple[Path, Path]) -> None:
    _, matrix = _matrix(corpus)
    cell = ("Colorectal cancer", "Feces", "coccus_shaped", "inc")
    figure = condition_body_site_figure(
        matrix,
        None,
        guild_order=TRAIT_ORDER,
        colorscales=SCALES,
        evidence=_evidence({cell: frozenset({SUPPORT, CELL_COMPOSITION})}),
    )
    outlines = [shape for shape in figure.layout.shapes if shape.line.color == "#111111"]
    assert sorted(shape.line.dash for shape in outlines) == ["dot", "solid"]
    widths = sorted(shape.x1 - shape.x0 for shape in outlines)
    assert widths[0] < widths[1]


def test_an_unknown_mark_type_is_an_error_not_a_default_style(corpus: tuple[Path, Path]) -> None:
    _, matrix = _matrix(corpus)
    cell = ("Colorectal cancer", "Feces", "coccus_shaped", "inc")
    with pytest.raises(KeyError):
        condition_body_site_figure(
            matrix,
            None,
            guild_order=TRAIT_ORDER,
            colorscales=SCALES,
            evidence=_evidence({cell: frozenset({"both"})}),
        )


def test_a_site_glyph_is_drawn_once_per_trait_over_both_columns(corpus: tuple[Path, Path]) -> None:
    """#217: the glyph is one statement about the trait, centred on its two columns."""
    _, matrix = _matrix(corpus)
    evidence = _evidence({}, columns={("Feces", "coccus_shaped"): "decreased"})
    for figure in (
        condition_body_site_figure(
            matrix, None, guild_order=TRAIT_ORDER, colorscales=SCALES, evidence=evidence
        ),
        body_site_panels_figure(
            matrix,
            None,
            sites=["Feces"],
            guild_order=TRAIT_ORDER,
            colorscales=SCALES,
            evidence=evidence,
        ),
    ):
        glyphs = [a for a in figure.layout.annotations if a.text in COLUMN_MARK_GLYPH.values()]
        assert len(glyphs) == 1
        assert glyphs[0].text == COLUMN_MARK_GLYPH["decreased"]
        # coccus shaped is the third trait: columns 4 and 5 within the block.
        assert glyphs[0].x == 4.5
        # Not rotated with the column labels.
        assert glyphs[0].textangle is None


def test_every_trait_keeps_its_colour_bar_when_a_panel_lacks_it(corpus: tuple[Path, Path]) -> None:
    """#217: the arsenate reduction bar vanished when its trace had no cell in the last panel."""
    _, matrix = _matrix(corpus)
    figure = body_site_panels_figure(
        matrix, None, sites=["Feces"], guild_order=TRAIT_ORDER, colorscales=SCALES
    )
    bars = [trace for trace in figure.data if trace.showscale]
    assert len(bars) == len(TRAIT_ORDER)
    positions = [bar.colorbar.x for bar in bars]
    assert len(set(positions)) == len(TRAIT_ORDER)
    for bar in bars:
        assert all(
            len(line) <= COLORBAR_TITLE_CHARS for line in bar.colorbar.title.text.split("<br>")
        )
        assert bar.colorbar.tickvals[-1] <= COLORBAR_TOP_TICK_SHARE * bar.zmax


def test_a_short_cut_is_as_tall_as_its_rows(corpus: tuple[Path, Path]) -> None:
    _, matrix = _matrix(corpus)
    cut = subset_matrix(matrix, conditions=["Asthma"])
    figure = body_site_panels_figure(
        cut, None, sites=cut.body_sites, guild_order=TRAIT_ORDER, colorscales=SCALES
    )
    rows = {label for trace in figure.data if not trace.showscale for label in trace.y}
    assert rows == {"Asthma"}


def test_a_cut_records_only_requested_conditions_as_dropped(corpus: tuple[Path, Path]) -> None:
    """#214: a condition that was never asked for is not reported as dropped."""
    _, matrix = _matrix(corpus)
    by_name = subset_matrix(matrix, conditions=["Asthma", "Nowhere"])
    assert by_name.dropped_empty_conditions == ["Nowhere"]
    by_site = subset_matrix(matrix, sites=["Feces"])
    assert by_site.dropped_empty_conditions == []


def test_uncovered_conditions_are_reported(corpus: tuple[Path, Path]) -> None:
    _, matrix = _matrix(corpus)
    panels = (FocusPanel(slug="one", subtitle="one", conditions=("Asthma",)),)
    assert uncovered_conditions(matrix, panels) == ["Colorectal cancer"]
    everything = (FocusPanel(slug="all", subtitle="all", sites=("Feces",)),)
    assert uncovered_conditions(matrix, everything) == []


def test_evidence_does_not_depend_on_the_display_cut(corpus: tuple[Path, Path]) -> None:
    """#176: narrowing the drawn conditions changes where marks land, not what is tested."""
    profiles, wide = _matrix(corpus, max_conditions=10)
    _, narrow = _matrix(corpus, max_conditions=1)
    assert set(narrow.conditions) < set(wide.conditions)
    overlay = build_overlay(profiles, TABLE1_TRAITS, corpus[0])
    evidence = evidence_marks(overlay, overlay.enrichment.tests_performed)
    # The overlay takes no display input at all; the marks and their totals are
    # the same object whichever matrix they are drawn on.
    subtitles = [
        str(
            condition_body_site_figure(
                matrix, None, guild_order=TRAIT_ORDER, colorscales=SCALES, evidence=evidence
            ).layout.title.text
        )
        for matrix in (wide, narrow)
    ]
    for kind in evidence.totals:
        total = evidence.totals[kind]
        assert all(f"Across the corpus {total} cells" in text for text in subtitles)
    dropped = next(c for c in wide.conditions if c not in narrow.conditions)
    cell = (dropped, "Feces", "facultatively_aerobic", "inc")
    assert cell_display_status(narrow, cell) == "condition_not_drawn"


def test_the_pooled_other_block_is_never_a_match() -> None:
    matrix = GuildMatrix(
        conditions=["A"],
        body_sites=[OTHER_LABEL],
        columns=[f"{OTHER_LABEL}|coccus_shaped_inc"],
        values=[[0.1]],
        study_counts=[[5]],
        row_clusters=[1],
        pooled=[],
        dropped_empty_conditions=[],
        dropped_body_sites=[],
        sampled_no_signal=[[False]],
    )
    assert (
        cell_display_status(matrix, ("A", OTHER_LABEL, "coccus_shaped", "inc")) == "site_not_drawn"
    )


@pytest.mark.integration
def test_a_one_row_cut_renders(tmp_path: Path, corpus: tuple[Path, Path]) -> None:
    """The padding rows existed because short cuts were believed not to render."""
    from bsdbng.guild_maps import save_figure

    _, matrix = _matrix(corpus)
    cut = subset_matrix(matrix, conditions=["Asthma"])
    path = save_figure(
        body_site_panels_figure(
            cut, None, sites=cut.body_sites, guild_order=TRAIT_ORDER, colorscales=SCALES
        ),
        tmp_path / "one_row.png",
    )
    assert path.stat().st_size > 0
