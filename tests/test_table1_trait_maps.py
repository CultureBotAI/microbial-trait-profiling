"""Table 1 trait selection, the paired directional shift, and its FDR control."""

import csv
from pathlib import Path

import pytest
from scipy.stats import binomtest, false_discovery_control

from bsdbng.guild_maps import SIGN_TOLERANCE, build_matrix, write_body_site_panels
from bsdbng.table1_trait_maps import (
    TRAIT_ORDER,
    MissingConditions,
    compute_direction_shifts,
    missing_condition_note,
    missing_conditions,
    trait_colorscale,
    write_shift_table,
)
from bsdbng.trait_guilds import (
    TABLE1_TRAIT_COLORS,
    TABLE1_TRAIT_LABELS,
    TABLE1_TRAITS,
    Profile,
    column_term_id,
    load_profiles,
)

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


def test_every_shift_row_names_its_correction_family(
    corpus: tuple[Path, Path], tmp_path: Path
) -> None:
    """Each exported q says which family it was corrected in and its size (#289)."""
    profiles, _ = load_profiles(*corpus, TABLE1_TRAITS)
    result = compute_direction_shifts(profiles, TRAIT_ORDER)
    shifts_path = write_shift_table(result, tmp_path, "v")[0]
    with shifts_path.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == len(TRAIT_ORDER)
    assert {(r["family_id"], r["family_size"], r["correction_method"]) for r in rows} == {
        ("table1_six_fractions", str(len(TRAIT_ORDER)), "Benjamini-Hochberg")
    }


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


def test_panel_figure_rejects_a_body_site_it_cannot_draw(
    tmp_path: Path, corpus: tuple[Path, Path]
) -> None:
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
        write_body_site_panels(
            matrix,
            tmp_path / "panels.png",
            None,
            sites=["Nowhere"],
            guild_order=TRAIT_ORDER,
        )


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
