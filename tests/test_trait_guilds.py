"""Corpus filters, guild aggregation, and the guild-map statistics."""

import math
from pathlib import Path

import pytest
from scipy.stats import binomtest

from bsdbng.count_inputs import COVERAGE_COLUMNS
from bsdbng.guild_maps import OTHER_LABEL, build_matrix, compute_correlation, study_direction_means
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


def _write_coverage(path: Path, keys: list[tuple[str, str, str]]) -> None:
    lines = ["\t".join(COVERAGE_COLUMNS)]
    for key in keys:
        lines.append(
            "\t".join(
                [
                    *key,
                    "1",
                    "1",
                    "0",
                    "1",
                    "0",
                    "complete_annotation_coverage",
                    "[]",
                    "[]",
                    "columns:6:sha256:test",
                ]
            )
        )
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


@pytest.mark.parametrize("bad_count", ["", "-1", "1.5", "NaN", "not-a-count"])
def test_method_loader_rejects_invalid_counts_before_filtering(
    tmp_path: Path, corpus: tuple[Path, Path], bad_count: str
) -> None:
    profiles_path, experiments_path = corpus
    text = profiles_path.read_text(encoding="utf-8")
    profiles_path = tmp_path / f"bad-{bad_count or 'blank'}.tsv"
    profiles_path.write_text(text.replace("\t6\t0\t0", f"\t{bad_count}\t0\t0", 1), encoding="utf-8")
    with pytest.raises(ValueError, match="count"):
        load_profiles(profiles_path, experiments_path, OXYGEN_GUILDS)


def test_method_loader_joins_coverage_without_turning_it_into_zero(tmp_path: Path) -> None:
    profiles = tmp_path / "profiles.tsv"
    experiments = tmp_path / "experiments.tsv"
    coverage = tmp_path / "coverage.tsv"
    keys = [("s1", "e1", "inc"), ("s1", "e1", "dec")]
    _write_profiles(
        profiles,
        [
            (*keys[0], [0, 0, 0, 0, 0, 0, 0]),
            (*keys[1], [0, 0, 0, 0, 0, 0, 0]),
        ],
    )
    _write_experiments(experiments, [("s1", "e1", "X", "Feces", "Homo sapiens")])
    _write_coverage(coverage, keys)

    loaded, stats = load_profiles(
        profiles,
        experiments,
        OXYGEN_GUILDS,
        coverage_path=coverage,
        require_coverage=True,
    )

    assert len(loaded) == 2
    assert all(profile.total == 0 for profile in loaded)
    assert all(profile.source_coverage is not None for profile in loaded)
    assert all(profile.source_coverage.is_complete for profile in loaded if profile.source_coverage)
    assert stats.coverage_rows == 2
    assert stats.incomplete_coverage_profiles == 0
    assert stats.coverage_not_provided_profiles == 0


def test_method_loader_rejects_profile_coverage_key_mismatch(tmp_path: Path) -> None:
    profiles = tmp_path / "profiles.tsv"
    experiments = tmp_path / "experiments.tsv"
    coverage = tmp_path / "coverage.tsv"
    keys = [("s1", "e1", "inc"), ("s1", "e1", "dec")]
    _write_profiles(
        profiles,
        [
            (*keys[0], [1, 0, 0, 0, 0, 0, 0]),
            (*keys[1], [1, 0, 0, 0, 0, 0, 0]),
        ],
    )
    _write_experiments(experiments, [("s1", "e1", "X", "Feces", "Homo sapiens")])
    _write_coverage(coverage, [keys[0]])
    with pytest.raises(ValueError, match="Profile/coverage key mismatch"):
        load_profiles(profiles, experiments, OXYGEN_GUILDS, coverage_path=coverage)


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
