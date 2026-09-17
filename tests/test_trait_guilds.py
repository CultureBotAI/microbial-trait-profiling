"""Corpus filters, guild aggregation, and the guild-map statistics."""

import math
import struct
from collections.abc import Callable
from pathlib import Path

import pytest
from scipy.stats import binomtest

from bsdbng.guild_maps import (
    CELL_SIZE_PX,
    MIN_FIGURE_WIDTH_PX,
    MIN_PLOT_HEIGHT_PX,
    OTHER_LABEL,
    GuildMatrix,
    _figure_height,
    _figure_width,
    _panel_figure_width,
    _panel_plot_width,
    build_matrix,
    compute_correlation,
    study_direction_means,
    write_body_site_panels,
    write_figure,
    write_matrix_outputs,
)
from bsdbng.trait_guilds import (
    GRAM_GUILDS,
    OXYGEN_GUILDS,
    Profile,
    load_profiles,
    read_bugsigdb_vintage,
)

AEROBE_A, AEROBE_B, AEROBE_C = OXYGEN_GUILDS["aerobe"]
ANAEROBE_A, ANAEROBE_B = OXYGEN_GUILDS["anaerobe"]
MEDIA_COLUMN = "METPO:2000517|mediadive.medium:J14"


@pytest.mark.parametrize("row_count", [0, 1, 40])
def test_figure_height_keeps_plot_area_outside_margins(row_count: int) -> None:
    top_margin = 445
    bottom_margin = 300
    plot_height = (
        _figure_height(row_count, top_margin=top_margin, bottom_margin=bottom_margin)
        - top_margin
        - bottom_margin
    )

    assert plot_height >= MIN_PLOT_HEIGHT_PX
    assert plot_height >= CELL_SIZE_PX * max(1, row_count)


@pytest.mark.parametrize("columns", [4, 60])
def test_figure_width_keeps_margins_outside_cell_area(columns: int) -> None:
    left_margin = 300
    right_margin = 180
    figure_width = _figure_width(columns, left_margin=left_margin, right_margin=right_margin)
    cell_area_width = CELL_SIZE_PX * columns
    effective_right_margin = figure_width - left_margin - cell_area_width

    assert figure_width >= MIN_FIGURE_WIDTH_PX
    assert effective_right_margin >= right_margin


def test_panel_width_accounts_for_subplot_spacing() -> None:
    left_margin = 300
    right_margin = 180
    columns_per_panel = [4, 12]
    horizontal_spacing = 0.06
    plot_width = _panel_plot_width(columns_per_panel, horizontal_spacing)
    figure_width = _panel_figure_width(
        columns_per_panel,
        horizontal_spacing=horizontal_spacing,
        left_margin=left_margin,
        right_margin=right_margin,
    )
    effective_right_margin = figure_width - left_margin - plot_width
    cell_domains_width = plot_width * (1.0 - horizontal_spacing * (len(columns_per_panel) - 1))

    assert figure_width >= MIN_FIGURE_WIDTH_PX
    assert effective_right_margin >= right_margin
    for column_count in columns_per_panel:
        panel_share = column_count / sum(columns_per_panel)
        assert cell_domains_width * panel_share / column_count >= CELL_SIZE_PX


def _one_row_matrix() -> GuildMatrix:
    columns = [
        "Feces|aerobe_inc",
        "Feces|aerobe_dec",
        "Feces|anaerobe_inc",
        "Feces|anaerobe_dec",
    ]
    return GuildMatrix(
        conditions=["One condition"],
        body_sites=["Feces"],
        columns=columns,
        values=[[0.2, 0.1, 0.3, 0.4]],
        study_counts=[[2, 2, 2, 2]],
        row_clusters=[1],
        pooled=[],
        dropped_empty_conditions=[],
        dropped_body_sites=[],
        sampled_no_signal=[[False] * len(columns)],
    )


@pytest.mark.parametrize(
    ("writer", "expected_height"),
    [
        (write_figure, 445 + 300 + MIN_PLOT_HEIGHT_PX),
        (write_body_site_panels, 425 + 260 + MIN_PLOT_HEIGHT_PX),
    ],
)
def test_one_row_figures_render_at_the_size_floors(
    tmp_path: Path,
    writer: Callable[..., Path],
    expected_height: int,
) -> None:
    """Exercise Kaleido: short matrices used to fail during axis scaling."""
    output = tmp_path / f"{writer.__name__}.png"
    kwargs = {"sites": ["Feces"]} if writer is write_body_site_panels else {}

    writer(_one_row_matrix(), output, None, **kwargs)

    png = output.read_bytes()
    assert png.startswith(b"\x89PNG\r\n\x1a\n")
    width, height = struct.unpack(">II", png[16:24])
    assert (width, height) == (MIN_FIGURE_WIDTH_PX * 2, expected_height * 2)


def _write_profiles(path: Path, rows: list[tuple[str, str, str, list[int]]]) -> None:
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
    matrix = build_matrix(profiles, max_conditions=10, max_body_sites=10, min_studies=1)
    result = compute_correlation(profiles, matrix)

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
    # The paired odds ratio is the discordant ratio; with no reverse-discordant
    # studies it is infinite.
    assert result.odds_ratio == math.inf


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

    write_matrix_outputs(matrix, tmp_path, prefix="table1_traits")
    status = (tmp_path / "table1_traits_condition_body_site_cell_status.tsv").read_text(
        encoding="utf-8"
    )
    assert "no_annotation_signal" in status
    assert "no_guild_signal" not in status


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
