"""Corpus filters, guild aggregation, and the guild-map statistics."""

import math
from pathlib import Path

import pytest
from scipy.stats import binomtest

from bsdbng.guild_maps import (
    OTHER_LABEL,
    build_matrix,
    compute_correlation,
    study_direction_means,
    subset_matrix,
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
