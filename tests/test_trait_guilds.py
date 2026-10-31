"""Corpus filters, guild aggregation, and the guild-map statistics."""

import gzip
import re
from collections.abc import Callable
from math import isnan
from pathlib import Path
from statistics import fmean

import pytest
from scipy.stats import binomtest, spearmanr

from bsdbng.guild_maps import (
    CELL_BELOW_MIN_STUDIES,
    CELL_NO_SIGNAL,
    CELL_NOT_SAMPLED,
    CELL_VALUE,
    OTHER_LABEL,
    build_matrix,
    cell_status,
    compute_correlation,
    require_plotly,
    study_direction_means,
    subset_matrix,
    write_figure,
)
from bsdbng.trait_guilds import (
    GRAM_GUILDS,
    OXYGEN_GUILDS,
    Profile,
    load_profiles,
    read_bugsigdb_vintage,
    resolve_trait_columns,
)

AEROBE_A, AEROBE_B, AEROBE_C = OXYGEN_GUILDS["aerobe"]
ANAEROBE_A, ANAEROBE_B = OXYGEN_GUILDS["anaerobe"]
MEDIA_COLUMN = "METPO:2000517|mediadive.medium:J14"


def _write_profiles(path: Path, rows: list[tuple[str, str, str, list[int | str]]]) -> None:
    header = [
        "study_id",
        "experiment_id",
        "inc_dec",
        AEROBE_A,
        AEROBE_B,
        AEROBE_C,
        ANAEROBE_A,
        ANAEROBE_B,
        "METPO:1000615",
        MEDIA_COLUMN,
    ]
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
    profiles = tmp_path / "profiles.tsv"
    experiments = tmp_path / "experiments.tsv"
    _write_profiles(
        profiles,
        [
            # Paired human experiment: aerobes up, anaerobes down.
            ("bsdb:1", "bsdb:1-1", "inc", [6, 0, 0, 1, 0, 3, 100]),
            ("bsdb:1", "bsdb:1-1", "dec", [1, 0, 0, 6, 0, 3, 100]),
            # Paired human experiment in a second condition.
            ("bsdb:2", "bsdb:2-1", "inc", [5, 1, 0, 2, 0, 2, 100]),
            ("bsdb:2", "bsdb:2-1", "dec", [1, 0, 0, 5, 1, 4, 100]),
            # Unpaired: increased only, must be dropped by the pair filter.
            ("bsdb:3", "bsdb:3-1", "inc", [4, 0, 0, 1, 0, 5, 100]),
            # Mouse study, must be dropped by the host filter.
            ("bsdb:4", "bsdb:4-1", "inc", [9, 0, 0, 1, 0, 1, 100]),
            ("bsdb:4", "bsdb:4-1", "dec", [1, 0, 0, 9, 0, 1, 100]),
        ],
    )
    _write_experiments(
        experiments,
        [
            ("bsdb:1", "bsdb:1-1", "Colorectal cancer", "Feces", "Homo sapiens"),
            ("bsdb:2", "bsdb:2-1", "Asthma", "Feces", "Homo sapiens"),
            ("bsdb:3", "bsdb:3-1", "Asthma", "Feces", "Homo sapiens"),
            ("bsdb:4", "bsdb:4-1", "Colitis", "Caecum", "Mus musculus"),
        ],
    )
    return profiles, experiments


def test_filters_drop_non_human_unpaired_and_mediadive(corpus: tuple[Path, Path]) -> None:
    profiles, stats = load_profiles(*corpus, OXYGEN_GUILDS)

    assert stats.dropped_non_human == 2
    assert stats.dropped_unpaired == 1
    assert stats.mediadive_columns_excluded == 1
    assert stats.kept_experiments == 2
    assert stats.kept_profiles == 4
    assert {p.study_id for p in profiles} == {"bsdb:1", "bsdb:2"}


def test_guild_counts_sum_member_terms_and_exclude_media(corpus: tuple[Path, Path]) -> None:
    profiles, _ = load_profiles(*corpus, OXYGEN_GUILDS)
    inc = next(p for p in profiles if p.study_id == "bsdb:2" and p.inc_dec == "inc")

    # aerobe = 5 + 1 + 0; anaerobe = 2 + 0; total excludes the MediaDive column.
    assert inc.guild_counts["aerobe"] == 6
    assert inc.guild_counts["anaerobe"] == 2
    assert inc.total == 5 + 1 + 0 + 2 + 0 + 2
    assert inc.fraction("aerobe") == pytest.approx(6 / 10)


def test_unmatched_experiments_are_counted_not_silently_dropped(
    tmp_path: Path, corpus: tuple[Path, Path]
) -> None:
    profiles_path, _ = corpus
    experiments = tmp_path / "sparse_experiments.tsv"
    _write_experiments(
        experiments, [("bsdb:1", "bsdb:1-1", "Colorectal cancer", "Feces", "Homo sapiens")]
    )
    _, stats = load_profiles(profiles_path, experiments, OXYGEN_GUILDS)

    assert stats.unmatched_experiments == 5
    assert stats.kept_experiments == 1


def test_missing_guild_column_raises(corpus: tuple[Path, Path]) -> None:
    profiles_path, experiments_path = corpus
    with pytest.raises(ValueError, match="gram positive"):
        load_profiles(profiles_path, experiments_path, GRAM_GUILDS)


def test_rows_with_no_surviving_cells_are_dropped(corpus: tuple[Path, Path]) -> None:
    """A condition whose every cell is masked would render as a blank row."""
    profiles, _ = load_profiles(*corpus, OXYGEN_GUILDS)

    # Each condition here is backed by a single study, so at a threshold of two
    # nothing survives and every row is dropped rather than drawn empty.
    masked = build_matrix(profiles, max_conditions=10, max_body_sites=10, min_studies=2)
    assert masked.conditions == []
    assert set(masked.dropped_empty_conditions) == {"Colorectal cancer", "Asthma", OTHER_LABEL}

    # At a threshold of one the same conditions carry values and are kept. Other
    # is still dropped: the cap is above the number of conditions, so it pooled
    # nothing and its row would be blank.
    kept = build_matrix(profiles, max_conditions=10, max_body_sites=10, min_studies=1)
    assert set(kept.conditions) == {"Colorectal cancer", "Asthma"}
    assert kept.dropped_empty_conditions == [OTHER_LABEL]
    assert any(value is not None for row in kept.values for value in row)


def test_correlation_counts_studies_by_shift_direction(corpus: tuple[Path, Path]) -> None:
    profiles, _ = load_profiles(*corpus, OXYGEN_GUILDS)
    result = compute_correlation(profiles)

    # Both surviving studies put aerobes higher on the increased side and
    # anaerobes higher on the decreased side.
    assert result.table == {
        "aerobe_inc": 2,
        "aerobe_dec": 0,
        "anaerobe_inc": 0,
        "anaerobe_dec": 2,
    }
    assert result.studies == 2

    # The same studies provide both guild signs, so the test is a paired one:
    # both studies are discordant in the same direction (aerobe up while
    # anaerobe down) and the exact McNemar is a two-sided binomial on the
    # discordant studies. Fisher's exact on the marginal table would treat the
    # two rows as independent samples and be anticonservative.
    assert result.discordant == (2, 0)
    assert result.mcnemar_p == pytest.approx(float(binomtest(2, 2, p=0.5).pvalue))
    assert result.mcnemar_p == pytest.approx(0.5)
    assert result.status == "tested"
    # The paired odds ratio is the discordant ratio. With no reverse-discordant
    # study the denominator is zero, so the ratio has no value - reported as
    # undefined rather than as infinity, which would print as a number (#289).
    assert result.odds_ratio is None


def test_a_guild_result_with_no_discordant_study_is_undefined_not_unity() -> None:
    """No discordant study means no McNemar test, not a test that found nothing.

    Writing ``p = 1.0`` here would read as "tested, no difference"; the truth is
    that the test does not exist. This is the fault the Table 1 arm already
    fixed, checked here for the guild arm (#289). It does not fire on the pinned
    corpus, which is why it needs a fixture.
    """
    profiles: list[Profile] = []
    for index in range(4):
        study = f"bsdb:{index}"
        # Both guilds shift the same way in every study, so no study is
        # discordant and the discordant table is empty on both sides.
        profiles.append(
            Profile(
                study_id=study,
                experiment_id=f"{study}-1",
                inc_dec="inc",
                guild_counts={"aerobe": 8, "anaerobe": 8},
                total=16,
                condition="c",
                body_site="s",
            )
        )
        profiles.append(
            Profile(
                study_id=study,
                experiment_id=f"{study}-1",
                inc_dec="dec",
                guild_counts={"aerobe": 4, "anaerobe": 4},
                total=16,
                condition="c",
                body_site="s",
            )
        )

    result = compute_correlation(profiles)

    assert result.discordant == (0, 0)
    assert result.mcnemar_p is None
    assert result.odds_ratio is None
    assert result.status == "all_tied"


def test_a_study_with_no_guild_annotation_is_not_ranked_as_a_zero() -> None:
    """An unannotated study is missing data, not a measured share of zero (#289).

    Ranking it as (0, 0) would put every such study at one identical point and
    pull the rank correlation toward that block. The McNemar already drops them
    as ties, so the two statistics would otherwise describe different populations.
    """

    def study(name: str, aerobe_inc: int, anaerobe_inc: int) -> list[Profile]:
        return [
            Profile(
                study_id=name,
                experiment_id=f"{name}-1",
                inc_dec="inc",
                guild_counts={"aerobe": aerobe_inc, "anaerobe": anaerobe_inc},
                total=20,
                condition="c",
                body_site="s",
            ),
            Profile(
                study_id=name,
                experiment_id=f"{name}-1",
                inc_dec="dec",
                guild_counts={"aerobe": anaerobe_inc, "anaerobe": aerobe_inc},
                total=20,
                condition="c",
                body_site="s",
            ),
        ]

    annotated: list[Profile] = []
    for index, (aerobe, anaerobe) in enumerate([(9, 1), (7, 3), (3, 7), (1, 9)]):
        annotated.extend(study(f"bsdb:a{index}", aerobe, anaerobe))
    unannotated = [p for index in range(3) for p in study(f"bsdb:z{index}", 0, 0)]

    with_blanks = compute_correlation(annotated + unannotated)
    without_blanks = compute_correlation(annotated)

    assert with_blanks.excluded_no_annotation == 3
    assert without_blanks.excluded_no_annotation == 0
    # The blank studies are held out, so both runs rank the same four studies and
    # agree on every correlation - including the count behind each one.
    assert with_blanks.spearman == without_blanks.spearman
    for _rho, _p, n_studies in with_blanks.spearman.values():
        assert n_studies == 4


def test_a_guild_result_with_no_scored_study_reports_no_studies() -> None:
    """A corpus where every guild contrast ties yields no scored study at all."""
    profiles = [
        Profile(
            study_id="bsdb:1",
            experiment_id="bsdb:1-1",
            inc_dec=direction,
            guild_counts={"aerobe": 5, "anaerobe": 5},
            total=20,
            condition="c",
            body_site="s",
        )
        for direction in ("inc", "dec")
    ]

    result = compute_correlation(profiles)

    assert result.status == "no_studies"
    assert result.mcnemar_p is None
    assert result.odds_ratio is None


def test_study_means_reject_an_unpaired_corpus() -> None:
    """Pooled means only equal paired contrasts when inc and dec counts match."""
    lonely = Profile(
        study_id="bsdb:9",
        experiment_id="bsdb:9-1",
        inc_dec="inc",
        condition="X",
        body_site="Y",
        guild_counts={"aerobe": 3, "anaerobe": 1},
        total=10,
    )
    with pytest.raises(ValueError, match="bsdb:9"):
        study_direction_means([lonely])
    # An explicitly unpaired caller can still opt in.
    assert "bsdb:9" in study_direction_means([lonely], require_paired=False)


def test_vintage_is_read_from_the_dump_header(tmp_path: Path) -> None:
    dump = tmp_path / "full_dump.csv"
    dump.write_text(
        "# BugSigDB 2026-04-24_00:41_UTC, License: CC BY 4.0\nBSDB ID,Study\n",
        encoding="utf-8",
    )
    assert read_bugsigdb_vintage(dump) == "2026-04-24_00:41_UTC"
    assert read_bugsigdb_vintage(tmp_path / "absent.csv") is None


def test_fraction_is_zero_when_a_profile_has_no_signal() -> None:
    empty = Profile(
        study_id="bsdb:9",
        experiment_id="bsdb:9-1",
        inc_dec="inc",
        condition="X",
        body_site="Y",
        guild_counts={"aerobe": 0},
        total=0,
    )
    assert empty.fraction("aerobe") == 0.0


def test_sampled_cells_without_guild_signal_are_flagged_not_drawn_as_zero(
    tmp_path: Path,
) -> None:
    """A site sampled for a condition but carrying no guild trait is its own state."""
    profiles_path = tmp_path / "profiles.tsv"
    experiments_path = tmp_path / "experiments.tsv"
    # Both guilds are absent; only the non-guild trait carries any count.
    _write_profiles(
        profiles_path,
        [
            ("bsdb:1", "bsdb:1-1", "inc", [0, 0, 0, 0, 0, 7, 100]),
            ("bsdb:1", "bsdb:1-1", "dec", [0, 0, 0, 0, 0, 7, 100]),
        ],
    )
    _write_experiments(
        experiments_path,
        [("bsdb:1", "bsdb:1-1", "Periodontitis", "Saliva", "Homo sapiens")],
    )
    profiles, _ = load_profiles(profiles_path, experiments_path, OXYGEN_GUILDS)
    matrix = build_matrix(profiles, max_conditions=10, max_body_sites=10, min_studies=1)

    saliva = [i for i, c in enumerate(matrix.columns) if c.startswith("Saliva|")]
    assert saliva, "expected Saliva columns"
    # The row has no values at all, but is kept: "sampled, no guild trait" is a
    # result, so only entirely blank rows are dropped.
    assert "Periodontitis" in matrix.conditions
    assert "Periodontitis" not in matrix.dropped_empty_conditions
    row = matrix.conditions.index("Periodontitis")

    for i in saliva:
        # Sampled, so a study is counted, but no value is drawn and the cell is
        # flagged so the figure can shade it distinctly from a small fraction.
        assert matrix.study_counts[row][i] == 1
        assert matrix.values[row][i] is None
        assert matrix.sampled_no_signal[row][i] is True


def test_rows_are_clustered_and_pooling_is_reported(corpus: tuple[Path, Path]) -> None:
    """Clustering reorders conditions, and the Other bucket reports what it hides."""
    profiles, _ = load_profiles(*corpus, OXYGEN_GUILDS)
    matrix = build_matrix(
        profiles,
        max_conditions=1,
        max_body_sites=10,
        min_studies=1,
        n_clusters=2,
    )

    assert len(matrix.row_clusters) == len(matrix.conditions)
    assert OTHER_LABEL in matrix.conditions

    condition_pooling = next(p for p in matrix.pooled if p.kind == "condition")
    # Two conditions survive the filters; keeping one pools the other.
    assert condition_pooling.kept == 1
    assert condition_pooling.pooled_labels == 1
    assert condition_pooling.pooled_profiles == 2
    assert "Other pools 1 labels" in condition_pooling.describe()


def test_other_bucket_can_be_excluded_instead_of_pooled(corpus: tuple[Path, Path]) -> None:
    """Dropping Other removes the mixed row and column instead of drawing them."""
    profiles, _ = load_profiles(*corpus, OXYGEN_GUILDS)
    matrix = build_matrix(
        profiles,
        max_conditions=1,
        max_body_sites=10,
        min_studies=1,
        n_clusters=2,
        include_other_conditions=False,
        include_other_body_sites=False,
    )

    assert OTHER_LABEL not in matrix.conditions
    assert OTHER_LABEL not in matrix.body_sites
    assert all(not column.startswith(f"{OTHER_LABEL}|") for column in matrix.columns)

    # The pooling report still names what the cut left out, now as dropped.
    condition_pooling = next(p for p in matrix.pooled if p.kind == "condition")
    assert condition_pooling.pooled_labels == 1
    assert condition_pooling.bucketed is False
    assert "dropped 1 labels" in condition_pooling.describe()


def test_rows_kept_only_by_no_signal_cells_go_with_the_shading(tmp_path: Path) -> None:
    """Without the dark-grey shading such a row would be an invisible blank line."""
    profiles_path = tmp_path / "profiles.tsv"
    experiments_path = tmp_path / "experiments.tsv"
    _write_profiles(
        profiles_path,
        [
            ("bsdb:1", "bsdb:1-1", "inc", [0, 0, 0, 0, 0, 7, 100]),
            ("bsdb:1", "bsdb:1-1", "dec", [0, 0, 0, 0, 0, 7, 100]),
        ],
    )
    _write_experiments(
        experiments_path,
        [("bsdb:1", "bsdb:1-1", "Periodontitis", "Saliva", "Homo sapiens")],
    )
    profiles, _ = load_profiles(profiles_path, experiments_path, OXYGEN_GUILDS)

    shaded = build_matrix(profiles, max_conditions=10, max_body_sites=10, min_studies=1)
    assert "Periodontitis" in shaded.conditions

    unshaded = build_matrix(
        profiles,
        max_conditions=10,
        max_body_sites=10,
        min_studies=1,
        mark_no_signal=False,
    )
    assert unshaded.conditions == []
    assert "Periodontitis" in unshaded.dropped_empty_conditions


def test_subset_derives_the_other_axis_and_drops_blank_bands(corpus: tuple[Path, Path]) -> None:
    """A cut named on one axis takes the other axis from where values were drawn."""
    profiles, _ = load_profiles(*corpus, OXYGEN_GUILDS)
    matrix = build_matrix(profiles, max_conditions=10, max_body_sites=10, min_studies=1)

    by_condition = subset_matrix(matrix, conditions=["Asthma"])
    assert by_condition.conditions == ["Asthma"]
    # Only the sites Asthma was actually drawn in survive, and every column that
    # remains belongs to one of them.
    assert by_condition.body_sites
    assert set(by_condition.body_sites) <= set(matrix.body_sites)
    assert all(
        column.split("|", 1)[0] in by_condition.body_sites for column in by_condition.columns
    )
    assert any(value is not None for row in by_condition.values for value in row)

    # Naming the axis the other way round drops the conditions with nothing in it.
    by_site = subset_matrix(matrix, sites=by_condition.body_sites)
    assert by_site.body_sites == by_condition.body_sites
    assert "Asthma" in by_site.conditions
    assert all(any(v is not None for v in row) for row in by_site.values)

    with pytest.raises(ValueError, match="Nowhere"):
        subset_matrix(matrix, sites=["Nowhere"])


def _study_guild_profile(
    study: str, experiment: str, direction: str, aerobe: int, anaerobe: int, condition: str
) -> Profile:
    return Profile(
        study_id=study,
        experiment_id=experiment,
        inc_dec=direction,
        condition=condition,
        body_site="Feces",
        guild_counts={"aerobe": aerobe, "anaerobe": anaerobe},
        total=aerobe + anaerobe + 10,
    )


def test_rank_correlations_count_each_study_once_across_its_cells() -> None:
    """#165: a study reaching several condition cells is one observation, not several.

    The fixture deliberately does not mirror the decreased side, so the two guild
    contrasts are not perfect opposites and rho is a value rather than -1 by
    construction (#289). The expectation is recomputed here from the four
    study-level contrasts, so counting bsdb:1 once per condition - the bug this
    test exists to catch - gives six observations and a different rho, and fails.
    """
    # (aerobe_inc, anaerobe_inc, aerobe_dec, anaerobe_dec) per experiment.
    spans = {
        "A": (9, 1, 2, 6),
        "B": (7, 2, 3, 5),
        "C": (8, 3, 1, 7),
    }
    profiles: list[Profile] = []
    for index, (condition, (a_inc, n_inc, a_dec, n_dec)) in enumerate(spans.items()):
        profiles += [
            _study_guild_profile("bsdb:1", f"bsdb:1-{index}", "inc", a_inc, n_inc, condition),
            _study_guild_profile("bsdb:1", f"bsdb:1-{index}", "dec", a_dec, n_dec, condition),
        ]
    others = {
        "bsdb:2": (6, 3, 4, 4),
        "bsdb:3": (2, 7, 5, 2),
        # bsdb:4 shifts BOTH guilds up, so the two contrast series are not
        # perfectly inversely ranked and rho is -0.8 rather than a forced -1.
        "bsdb:4": (6, 6, 2, 2),
    }
    for study, (a_inc, n_inc, a_dec, n_dec) in others.items():
        profiles += [
            _study_guild_profile(study, f"{study}-1", "inc", a_inc, n_inc, "A"),
            _study_guild_profile(study, f"{study}-1", "dec", a_dec, n_dec, "A"),
        ]

    result = compute_correlation(profiles)

    # Four studies, however many cells each reaches.
    assert {n for _rho, _p, n in result.spearman.values()} == {4}

    # Recompute the expectation independently: one contrast per study, where a
    # study's contrast averages its experiments' fractions on each side.
    def share(aerobe: int, anaerobe: int, guild: str) -> float:
        count = aerobe if guild == "aerobe" else anaerobe
        return count / (aerobe + anaerobe + 10)

    def contrast(experiments: list[tuple[int, int, int, int]], guild: str) -> float:
        increased = fmean([share(a_inc, n_inc, guild) for a_inc, n_inc, _, _ in experiments])
        decreased = fmean([share(a_dec, n_dec, guild) for _, _, a_dec, n_dec in experiments])
        return increased - decreased

    by_study = [list(spans.values()), *([value] for value in others.values())]
    expected_rho, _ = spearmanr(
        [contrast(experiments, "aerobe") for experiments in by_study],
        [contrast(experiments, "anaerobe") for experiments in by_study],
    )

    rho, _p, n_studies = result.spearman["aerobe_contrast~anaerobe_contrast"]
    assert n_studies == 4
    assert rho == pytest.approx(float(expected_rho))
    # Not the degenerate -1 that a mirrored fixture would force regardless of code.
    assert -1.0 < rho < 1.0


def test_a_body_site_containing_the_column_separator_is_rejected() -> None:
    profiles = [
        Profile("bsdb:1", "bsdb:1-1", "inc", "A", "Plaque | saliva", {"aerobe": 1}, 2),
        Profile("bsdb:1", "bsdb:1-1", "dec", "A", "Plaque | saliva", {"aerobe": 1}, 2),
    ]
    with pytest.raises(ValueError, match="Plaque"):
        build_matrix(profiles, max_conditions=5, max_body_sites=5, min_studies=1)


def test_a_missing_viz_extra_names_the_install_command(monkeypatch: pytest.MonkeyPatch) -> None:
    """plotly is a 'viz' extra, so its absence must read as a setup step (#288)."""
    import importlib.util

    real_find_spec = importlib.util.find_spec

    def without_plotly(name: str, package: str | None = None) -> object | None:
        if name == "plotly":
            return None
        return real_find_spec(name, package)

    monkeypatch.setattr(importlib.util, "find_spec", without_plotly)
    with pytest.raises(ModuleNotFoundError, match="viz"):
        require_plotly()


def test_require_plotly_passes_when_the_extra_is_installed() -> None:
    """The guard is a probe, not a second import path: present means it does nothing."""
    require_plotly()


@pytest.mark.parametrize("site_headers", [False, True])
def test_the_plotting_area_starts_below_the_last_subtitle_line(
    tmp_path: Path,
    rendered_title_clearance: Callable[[Path], float],
    site_headers: bool,
) -> None:
    """Require breathing room and stable clearance after a two-line note (#839, #929)."""
    # The wrapped site name makes this fixture sensitive to reduced header spacing.
    matrix = build_matrix(
        [
            _guild_profile(str(index), f"Condition {index}", "Small intestine terminal ileum")
            for index in range(8)
        ],
        max_conditions=10,
        max_body_sites=10,
        min_studies=1,
        guild_order=("aerobe",),
    )
    assert len(matrix.conditions) >= 8, (
        "Rendered layout fixture needs at least eight condition rows"
    )
    clearances = []
    # A <br> inside one note entry reproduces the #849 line-counting regression.
    for name, note in (("short", None), ("long", "Population note.<br>A second population line.")):
        path = tmp_path / f"{name}.png"
        write_figure(
            matrix,
            path,
            "v1.3.1",
            guild_order=("aerobe",),
            note=note,
            separate_body_sites=site_headers,
        )
        clearances.append(rendered_title_clearance(path))
    assert min(clearances) >= 30, f"Title/header clearance below 30 layout px: {clearances}"
    assert abs(clearances[1] - clearances[0]) <= 4, (
        f"A two-line note changed title/header clearance by more than 4 layout px: {clearances}"
    )


def _fixed_total_profile(
    study: str, direction: str, aerobe: int, anaerobe: int, total: int = 20
) -> Profile:
    """One profile on a fixed denominator, so a guild fraction can be held constant.

    ``_study_guild_profile`` derives the total from the two counts, which makes every
    fraction move when either count does. Pinning the total is what lets a
    fixture vary one direction's series while holding another's flat.
    """
    return Profile(
        study_id=study,
        experiment_id=f"{study}-1",
        inc_dec=direction,
        condition="A",
        body_site="Feces",
        guild_counts={"aerobe": aerobe, "anaerobe": anaerobe},
        total=total,
    )


# (aerobe_inc, anaerobe_inc, aerobe_dec, anaerobe_dec) per study, on a total of 20.
# aerobe_inc is 4 in every study, so the aerobe increased series is constant at
# 0.20 while the other three series take four distinct values each.
_CONSTANT_SERIES_STUDIES = {
    "bsdb:1": (4, 1, 2, 7),
    "bsdb:2": (4, 2, 3, 5),
    "bsdb:3": (4, 3, 5, 4),
    "bsdb:4": (4, 6, 9, 1),
}
_VARYING_FIRST_STUDY = (1, 1, 2, 7)


def _constant_series_profiles(*, hold_aerobe_inc_constant: bool) -> list[Profile]:
    rows = dict(_CONSTANT_SERIES_STUDIES)
    if not hold_aerobe_inc_constant:
        # One study's increased aerobe count differs, so that series varies too.
        rows["bsdb:1"] = _VARYING_FIRST_STUDY
    profiles: list[Profile] = []
    for study, (a_inc, n_inc, a_dec, n_dec) in rows.items():
        profiles.append(_fixed_total_profile(study, "inc", a_inc, n_inc))
        profiles.append(_fixed_total_profile(study, "dec", a_dec, n_dec))
    return profiles


def test_a_constant_guild_series_is_omitted_from_the_correlations_not_written_as_nan() -> None:
    """Spearman against a flat series is undefined; the pair must be absent (#289).

    ``spearmanr`` returns nan for a constant input rather than raising, so
    dropping the guard writes nan into the correlation table and into the figure.
    The fixture carries four studies, so the pair is past the ``len(paired) < 3``
    sample-size guard: the control run below varies the same series and recovers
    every pair at the same n, which is what separates the two guards.
    """
    result = compute_correlation(_constant_series_profiles(hold_aerobe_inc_constant=True))

    # The studies are in the ranked population: none is dropped as unannotated.
    assert result.studies == 4
    assert result.excluded_no_annotation == 0

    # Exact key set, because "aerobe_inc" is a substring of "anaerobe_inc" and a
    # substring test would wrongly pass by also excluding the anaerobe pairs.
    assert set(result.spearman) == {
        "aerobe_dec~anaerobe_inc",
        "aerobe_dec~anaerobe_dec",
        "anaerobe_inc~anaerobe_dec",
        "aerobe_contrast~anaerobe_contrast",
    }
    assert {n for _rho, _p, n in result.spearman.values()} == {4}
    assert not [key for key, (rho, _p, _n) in result.spearman.items() if isnan(rho)]

    control = compute_correlation(_constant_series_profiles(hold_aerobe_inc_constant=False))

    # Same four studies, same denominators: only the constancy changed, and the
    # three aerobe_inc pairs come back with the same n.
    assert set(control.spearman) == set(result.spearman) | {
        "aerobe_inc~aerobe_dec",
        "aerobe_inc~anaerobe_inc",
        "aerobe_inc~anaerobe_dec",
    }
    assert {n for _rho, _p, n in control.spearman.values()} == {4}
    assert not [key for key, (rho, _p, _n) in control.spearman.items() if isnan(rho)]


def test_cell_status_separates_the_three_reasons_a_cell_is_blank() -> None:
    """A blank cell is one of three different claims, not one (#294).

    Both the figure writer and the integrated evidence table classify cells
    here, so the distinction has to hold in one place: a cell held back for too
    little support, a cell sampled in which nothing carried the annotation, and
    a cell nothing was sampled in.
    """
    assert cell_status(0.42, False, 7) == CELL_VALUE
    # A zero mean is a value, not an absence.
    assert cell_status(0.0, False, 7) == CELL_VALUE
    assert cell_status(None, True, 7) == CELL_NO_SIGNAL
    assert cell_status(None, False, 2) == CELL_BELOW_MIN_STUDIES
    assert cell_status(None, False, 0) == CELL_NOT_SAMPLED
    # The support state must come from the study count, not from the value.
    assert cell_status(None, False, 1) != cell_status(None, False, 0)


def _guild_profile(
    study: str,
    condition: str,
    site: str,
    *,
    signal: int = 5,
    experiment: str = "1",
    direction: str = "inc",
) -> Profile:
    return Profile(
        study_id=study,
        experiment_id=f"{study}-{experiment}",
        inc_dec=direction,
        condition=condition,
        body_site=site,
        guild_counts={"aerobe": signal},
        total=10,
    )


@pytest.mark.parametrize("mark_no_signal", [False, True])
def test_subsets_apply_the_combined_matrix_no_signal_policy(mark_no_signal: bool) -> None:
    """No-signal cells keep both rows and derived sites only while shaded."""
    matrix = build_matrix(
        [
            _guild_profile("1", "A", "Feces"),
            _guild_profile("2", "A", "Saliva", signal=0),
            _guild_profile("3", "B", "Saliva"),
            _guild_profile("4", "NoSig", "Saliva", signal=0),
            _guild_profile("5", "Outside", "Feces"),
        ],
        max_conditions=10,
        max_body_sites=10,
        min_studies=1,
        guild_order=("aerobe",),
        include_other_conditions=False,
        include_other_body_sites=False,
        mark_no_signal=mark_no_signal,
    )
    assert matrix.mark_no_signal is mark_no_signal

    by_condition = subset_matrix(matrix, conditions=["A"])
    assert by_condition.conditions == ["A"]
    assert by_condition.body_sites == (["Feces", "Saliva"] if mark_no_signal else ["Feces"])
    assert by_condition.mark_no_signal is mark_no_signal

    by_site = subset_matrix(matrix, sites=["Saliva"])
    assert set(by_site.conditions) == ({"A", "B", "NoSig"} if mark_no_signal else {"B"})
    assert by_site.mark_no_signal is mark_no_signal
    assert "Outside" not in by_site.conditions

    only_no_signal = subset_matrix(matrix, conditions=["NoSig"])
    assert only_no_signal.conditions == (["NoSig"] if mark_no_signal else [])
    assert only_no_signal.body_sites == (["Saliva"] if mark_no_signal else [])
    if mark_no_signal:
        assert only_no_signal.values == [[None, None]]
        assert only_no_signal.study_counts == [[1, 0]]
        assert only_no_signal.sampled_no_signal == [[True, False]]


def test_subset_preserves_cell_evidence_and_cluster_labels_when_sites_are_reordered() -> None:
    matrix = build_matrix(
        [
            _guild_profile("1", "A", "Feces", signal=2),
            _guild_profile("2", "A", "Feces", signal=6),
            _guild_profile("3", "A", "Saliva", signal=0),
            _guild_profile("4", "B", "Saliva", signal=8),
            _guild_profile("5", "C", "Feces", signal=9),
        ],
        max_conditions=10,
        max_body_sites=10,
        min_studies=1,
        guild_order=("aerobe",),
        include_other_conditions=False,
        include_other_body_sites=False,
    )
    cut = subset_matrix(matrix, conditions=["B", "A"], sites=["Saliva", "Feces"])
    assert cut.conditions == [name for name in matrix.conditions if name in {"A", "B"}]
    assert cut.body_sites == ["Saliva", "Feces"]
    assert cut.columns == [
        "Saliva|aerobe_inc",
        "Saliva|aerobe_dec",
        "Feces|aerobe_inc",
        "Feces|aerobe_dec",
    ]
    a_row = cut.conditions.index("A")
    b_row = cut.conditions.index("B")
    assert cut.values[a_row] == [None, None, pytest.approx(0.4), None]
    assert cut.study_counts[a_row] == [1, 0, 2, 0]
    assert cut.sampled_no_signal[a_row] == [True, False, False, False]
    assert cut.values[b_row] == [pytest.approx(0.8), None, None, None]
    assert cut.study_counts[b_row] == [1, 0, 0, 0]
    assert cut.sampled_no_signal[b_row] == [False, False, False, False]
    assert cut.row_clusters == [
        matrix.row_clusters[matrix.conditions.index(name)] for name in cut.conditions
    ]


def test_site_block_drops_are_distinct_from_study_masking() -> None:
    profiles = [
        _guild_profile(study, condition, site)
        for condition, site in [
            ("SiteDropped", "Rare"),
            ("Kept", "Shared"),
            ("Both", "Shared"),
            ("Both", "Elsewhere"),
        ]
        for study in ("1", "2")
    ]
    profiles.append(_guild_profile("3", "Masked", "Rare"))
    matrix = build_matrix(
        profiles,
        max_conditions=10,
        max_body_sites=10,
        min_studies=2,
        min_site_conditions=2,
        guild_order=("aerobe",),
        include_other_conditions=False,
        include_other_body_sites=False,
        mark_no_signal=False,
    )
    assert matrix.body_sites == ["Shared"]
    assert set(matrix.conditions) == {"Kept", "Both"}
    assert set(matrix.dropped_body_sites) == {"Rare", "Elsewhere"}
    assert matrix.dropped_site_conditions == ["SiteDropped"]
    assert matrix.dropped_empty_conditions == ["Masked"]


def test_excluded_sites_require_real_per_cell_study_support() -> None:
    """Site exclusions must not manufacture support by pooling sites or repeats."""
    matrix = build_matrix(
        [
            _guild_profile("1", "Supported", "Saliva"),
            _guild_profile("2", "Supported", "Saliva"),
            _guild_profile("1", "Repeated", "Saliva", experiment="1"),
            _guild_profile("1", "Repeated", "Saliva", experiment="2"),
            _guild_profile("1", "SplitSites", "Saliva"),
            _guild_profile("2", "SplitSites", "Feces"),
            _guild_profile("1", "SplitDirections", "Saliva", direction="inc"),
            _guild_profile("2", "SplitDirections", "Saliva", direction="dec"),
        ],
        max_conditions=10,
        max_body_sites=0,
        min_studies=2,
        guild_order=("aerobe",),
        include_other_conditions=False,
        include_other_body_sites=False,
        mark_no_signal=False,
    )
    assert matrix.conditions == []
    assert matrix.body_sites == []
    assert matrix.dropped_site_conditions == ["Supported"]
    assert set(matrix.dropped_empty_conditions) == {"Repeated", "SplitSites", "SplitDirections"}


def test_site_cap_does_not_drop_a_condition_drawn_in_a_retained_site() -> None:
    matrix = build_matrix(
        [
            _guild_profile("1", "Both", "Feces"),
            _guild_profile("2", "Both", "Saliva"),
            _guild_profile("3", "SiteDropped", "Saliva"),
            _guild_profile("4", "Anchor", "Feces"),
            _guild_profile("5", "Anchor", "Feces"),
        ],
        max_conditions=10,
        max_body_sites=1,
        min_studies=1,
        guild_order=("aerobe",),
        include_other_conditions=False,
        include_other_body_sites=False,
        mark_no_signal=False,
    )
    assert matrix.body_sites == ["Feces"]
    assert set(matrix.conditions) == {"Both", "Anchor"}
    assert matrix.dropped_site_conditions == ["SiteDropped"]
    assert matrix.dropped_empty_conditions == []


@pytest.mark.parametrize("mark_no_signal", [False, True])
@pytest.mark.parametrize("max_body_sites", [0, 10])
def test_site_drop_attribution_honors_no_signal_shading(
    mark_no_signal: bool, max_body_sites: int
) -> None:
    matrix = build_matrix(
        [_guild_profile("1", "NoSig", "Saliva", signal=0)],
        max_conditions=10,
        max_body_sites=max_body_sites,
        min_studies=3,
        min_site_conditions=2,
        guild_order=("aerobe",),
        include_other_conditions=False,
        include_other_body_sites=False,
        mark_no_signal=mark_no_signal,
    )
    assert matrix.conditions == []
    assert matrix.dropped_site_conditions == (["NoSig"] if mark_no_signal else [])
    if mark_no_signal:
        assert "NoSig" not in matrix.dropped_empty_conditions


def test_dropped_other_site_preserves_condition_pool_membership() -> None:
    matrix = build_matrix(
        [
            *[_guild_profile(str(study), "Anchor", "Feces") for study in range(4)],
            _guild_profile("4", "PooledA", "Feces"),
            _guild_profile("5", "PooledB", "Feces"),
        ],
        max_conditions=1,
        max_body_sites=0,
        min_studies=2,
        min_site_conditions=3,
        guild_order=("aerobe",),
        mark_no_signal=False,
    )
    assert matrix.conditions == []
    assert matrix.dropped_body_sites == [OTHER_LABEL]
    assert set(matrix.dropped_site_conditions) == {"Anchor", OTHER_LABEL}
    assert matrix.dropped_empty_conditions == []
    condition_pool = next(pool for pool in matrix.pooled if pool.kind == "condition")
    assert condition_pool.bucketed
    assert condition_pool.pooled_names == ("PooledA", "PooledB")


@pytest.mark.parametrize("bad_count", ["", "-1", "1.5", "NaN", "not-a-count", "3_0"])
def test_method_loader_rejects_invalid_counts_before_filtering(
    tmp_path: Path, bad_count: str
) -> None:
    profiles_path = tmp_path / "profiles.tsv"
    experiments_path = tmp_path / "experiments.tsv"
    _write_profiles(
        profiles_path,
        [
            ("bsdb:1", "bsdb:1-1", "inc", [6, 0, 0, 1, 0, 3, 100]),
            ("bsdb:1", "bsdb:1-1", "dec", [1, 0, 0, 6, 0, 3, 100]),
            ("bsdb:4", "bsdb:4-1", "inc", [9, 0, 0, 1, 0, bad_count, 100]),
        ],
    )
    _write_experiments(
        experiments_path,
        [
            ("bsdb:1", "bsdb:1-1", "Colorectal cancer", "Feces", "Homo sapiens"),
            ("bsdb:4", "bsdb:4-1", "Colitis", "Caecum", "Mus musculus"),
        ],
    )

    with pytest.raises(ValueError, match="count"):
        load_profiles(profiles_path, experiments_path, OXYGEN_GUILDS)

    if bad_count != "":
        with pytest.raises(ValueError, match="count"):
            load_profiles(
                profiles_path,
                experiments_path,
                OXYGEN_GUILDS,
                reject_missing_counts=False,
            )


def test_method_loader_rejects_blank_by_default_and_allows_sparse_zero_override(
    tmp_path: Path,
) -> None:
    profiles_path = tmp_path / "profiles.tsv"
    experiments_path = tmp_path / "experiments.tsv"
    _write_profiles(
        profiles_path,
        [
            ("bsdb:1", "bsdb:1-1", "inc", [6, 0, 0, 1, 0, "", 100]),
            ("bsdb:1", "bsdb:1-1", "dec", [1, 0, 0, 6, 0, 3, 100]),
        ],
    )
    _write_experiments(
        experiments_path,
        [("bsdb:1", "bsdb:1-1", "Colorectal cancer", "Feces", "Homo sapiens")],
    )

    # One label, one statement that missing is not zero, and the file:line:column
    # of the blank cell (#1300).
    with pytest.raises(
        ValueError,
        match=(
            rf"^Missing trait count in {re.escape(str(profiles_path))}:2, "
            r"column 'METPO:1000615'; missing is not an observed zero$"
        ),
    ):
        load_profiles(profiles_path, experiments_path, OXYGEN_GUILDS)

    profiles, _ = load_profiles(
        profiles_path,
        experiments_path,
        OXYGEN_GUILDS,
        reject_missing_counts=False,
    )
    inc = next(profile for profile in profiles if profile.inc_dec == "inc")
    assert inc.total == 7


def test_method_loader_relabels_only_missing_count_errors_when_the_path_names_the_marker(
    tmp_path: Path,
) -> None:
    # The relabel keys on the start of the reader's message, not on the phrase
    # appearing anywhere in it, so a path carrying the phrase changes nothing (#1300).
    profiles_path = tmp_path / "Missing count in here.tsv"
    experiments_path = tmp_path / "experiments.tsv"
    _write_experiments(
        experiments_path,
        [("bsdb:1", "bsdb:1-1", "Colorectal cancer", "Feces", "Homo sapiens")],
    )
    _write_profiles(
        profiles_path,
        [
            ("bsdb:1", "bsdb:1-1", "inc", [6, 0, 0, 1, 0, 3, 100]),
            ("bsdb:1", "bsdb:1-1", "dec", [1, 0, 0, 6, 0, 3, 100]),
            ("bsdb:1", "bsdb:1-1", "up", [1, 0, 0, 6, 0, 3, 100]),
        ],
    )

    path = re.escape(str(profiles_path))
    with pytest.raises(
        ValueError, match=rf"^Unknown profile direction in {path}:4: 'up'; "
    ) as unknown:
        load_profiles(profiles_path, experiments_path, OXYGEN_GUILDS)
    assert "Missing trait count" not in str(unknown.value)

    _write_profiles(
        profiles_path,
        [
            ("bsdb:1", "bsdb:1-1", "inc", [6, 0, 0, 1, 0, "", 100]),
            ("bsdb:1", "bsdb:1-1", "dec", [1, 0, 0, 6, 0, 3, 100]),
        ],
    )
    with pytest.raises(
        ValueError,
        match=(
            rf"^Missing trait count in {path}:2, column 'METPO:1000615'; "
            r"missing is not an observed zero$"
        ),
    ):
        load_profiles(profiles_path, experiments_path, OXYGEN_GUILDS)


@pytest.mark.parametrize("strict", [False, True])
def test_method_loader_rejects_duplicate_profile_keys(tmp_path: Path, strict: bool) -> None:
    profiles_path = tmp_path / "profiles.tsv"
    experiments_path = tmp_path / "experiments.tsv"
    _write_profiles(
        profiles_path,
        [
            ("bsdb:1", "bsdb:1-1", "inc", [6, 0, 0, 1, 0, 3, 100]),
            ("bsdb:1", "bsdb:1-1", "inc", [6, 0, 0, 1, 0, 3, 100]),
            ("bsdb:1", "bsdb:1-1", "dec", [1, 0, 0, 6, 0, 3, 100]),
        ],
    )
    _write_experiments(
        experiments_path,
        [("bsdb:1", "bsdb:1-1", "Colorectal cancer", "Feces", "Homo sapiens")],
    )

    # Anchored, so the reader's label must survive the loader unchanged (#1295).
    duplicate = re.escape(repr(("bsdb:1", "bsdb:1-1", "inc")))
    with pytest.raises(
        ValueError,
        match=rf"^Duplicate key {duplicate} in {re.escape(str(profiles_path))}:3; ",
    ):
        load_profiles(
            profiles_path,
            experiments_path,
            OXYGEN_GUILDS,
            reject_missing_counts=strict,
        )


@pytest.mark.parametrize(
    ("study_id", "experiment_id"),
    [
        pytest.param("bsdb:1", "bsdb:1-1", id="kept_human_row"),
        pytest.param("bsdb:4", "bsdb:4-1", id="non_human_row"),
        pytest.param("bsdb:9", "bsdb:9-1", id="unmatched_row"),
    ],
)
@pytest.mark.parametrize("strict", [False, True])
def test_method_loader_rejects_unknown_direction(
    tmp_path: Path, strict: bool, study_id: str, experiment_id: str
) -> None:
    # The direction is checked before the metadata join and the host filter, so
    # an unknown direction fails even on a row those filters would drop (#1295).
    profiles_path = tmp_path / "profiles.tsv"
    experiments_path = tmp_path / "experiments.tsv"
    _write_profiles(
        profiles_path,
        [
            ("bsdb:1", "bsdb:1-1", "inc", [6, 0, 0, 1, 0, 3, 100]),
            ("bsdb:1", "bsdb:1-1", "dec", [1, 0, 0, 6, 0, 3, 100]),
            (study_id, experiment_id, "up", [1, 0, 0, 6, 0, 3, 100]),
        ],
    )
    _write_experiments(
        experiments_path,
        [
            ("bsdb:1", "bsdb:1-1", "Colorectal cancer", "Feces", "Homo sapiens"),
            ("bsdb:4", "bsdb:4-1", "Colitis", "Caecum", "Mus musculus"),
        ],
    )

    with pytest.raises(
        ValueError,
        match=(
            rf"^Unknown profile direction in {re.escape(str(profiles_path))}:4: "
            r"'up'; expected one of "
        ),
    ):
        load_profiles(
            profiles_path,
            experiments_path,
            OXYGEN_GUILDS,
            reject_missing_counts=strict,
        )


@pytest.mark.parametrize("bad_count", ["-1", "3_0"])
@pytest.mark.parametrize("strict", [False, True])
def test_method_loader_rejects_negative_and_underscore_counts_on_kept_human_rows(
    tmp_path: Path, bad_count: str, strict: bool
) -> None:
    profiles_path = tmp_path / "profiles.tsv"
    experiments_path = tmp_path / "experiments.tsv"
    _write_profiles(
        profiles_path,
        [
            ("bsdb:1", "bsdb:1-1", "inc", [6, 0, 0, 1, 0, bad_count, 100]),
            ("bsdb:1", "bsdb:1-1", "dec", [1, 0, 0, 6, 0, 3, 100]),
        ],
    )
    _write_experiments(
        experiments_path,
        [("bsdb:1", "bsdb:1-1", "Colorectal cancer", "Feces", "Homo sapiens")],
    )

    with pytest.raises(ValueError, match="count"):
        load_profiles(
            profiles_path,
            experiments_path,
            OXYGEN_GUILDS,
            reject_missing_counts=strict,
        )


def test_method_loader_validates_yet_excludes_metadata_columns_from_total(
    tmp_path: Path,
) -> None:
    profiles_path = tmp_path / "profiles_with_taxa_count.tsv"
    experiments_path = tmp_path / "experiments.tsv"
    header = [
        "study_id",
        "experiment_id",
        "inc_dec",
        "taxa_count",
        AEROBE_A,
        AEROBE_B,
        AEROBE_C,
        ANAEROBE_A,
        ANAEROBE_B,
        "METPO:1000615",
        MEDIA_COLUMN,
    ]

    def write(inc_taxa_count: str) -> None:
        rows = [
            "\t".join(header),
            "\t".join(
                ["bsdb:1", "bsdb:1-1", "inc", inc_taxa_count, "6", "0", "0", "1", "0", "3", "100"]
            ),
            "\t".join(["bsdb:1", "bsdb:1-1", "dec", "40", "1", "0", "0", "6", "0", "3", "100"]),
        ]
        profiles_path.write_text("\n".join(rows) + "\n", encoding="utf-8")

    write("50")
    _write_experiments(
        experiments_path,
        [("bsdb:1", "bsdb:1-1", "Colorectal cancer", "Feces", "Homo sapiens")],
    )

    profiles, _ = load_profiles(profiles_path, experiments_path, OXYGEN_GUILDS)
    inc = next(profile for profile in profiles if profile.inc_dec == "inc")
    assert inc.total == 10

    # Excluded from the total is not exempt from validation: a malformed or
    # negative metadata cell fails in either mode, and a blank one fails by
    # default (#1295).
    location = rf"{re.escape(str(profiles_path))}:2, column 'taxa_count'"
    for bad in ("-5", "4.5", "fifty"):
        write(bad)
        for strict in (True, False):
            with pytest.raises(
                ValueError,
                match=(
                    rf"^Expected a finite nonnegative integer count in {location}; "
                    rf"got {re.escape(repr(bad))}$"
                ),
            ):
                load_profiles(
                    profiles_path, experiments_path, OXYGEN_GUILDS, reject_missing_counts=strict
                )
    write("")
    with pytest.raises(ValueError, match=rf"^Missing trait count in {location}; "):
        load_profiles(profiles_path, experiments_path, OXYGEN_GUILDS)


def test_method_loader_and_column_resolver_read_tsv_gz_like_tsv(
    tmp_path: Path, corpus: tuple[Path, Path]
) -> None:
    # The same matrix, gzip-compressed, loads to the same profiles, filter
    # counts and resolved columns (#1295).
    plain, experiments = corpus
    packed = tmp_path / "profiles.tsv.gz"
    packed.write_bytes(gzip.compress(plain.read_bytes()))

    plain_profiles, plain_stats = load_profiles(plain, experiments, OXYGEN_GUILDS)
    packed_profiles, packed_stats = load_profiles(packed, experiments, OXYGEN_GUILDS)
    assert plain_stats.kept_profiles == 4
    assert packed_profiles == plain_profiles
    assert packed_stats == plain_stats

    plain_columns = resolve_trait_columns(plain, OXYGEN_GUILDS)
    assert set(plain_columns) == set(OXYGEN_GUILDS)
    assert resolve_trait_columns(packed, OXYGEN_GUILDS) == plain_columns
