"""Every guild-map summary table must be independent of the hash seed (#2322)."""

from __future__ import annotations

import csv
import itertools
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from mtp import guild_maps
from mtp.trait_guilds import Profile

EXPECTED_FILES = {
    "guild_condition_body_site_matrix.tsv",
    "guild_condition_body_site_study_counts.tsv",
    "guild_condition_body_site_cell_status.tsv",
    "guild_condition_row_clusters.tsv",
    "guild_aerobe_anaerobe_correlation.tsv",
}
CONDITIONS = ("condition_a", "condition_b")
SITES = ("Feces", "Oral")
DIRECTION_KEYS = ("aerobe_inc", "aerobe_dec", "anaerobe_inc", "anaerobe_dec")


def write_guild_fixture(output: Path) -> None:
    """Compute real matrix/correlation summaries from paired, unequal profile counts."""
    profiles = []
    # The first cell's study means average to 0.2000005, a six-decimal boundary.
    # Two unequal experiments per study and cell retain within-study averaging.
    # An unordered study reduction can therefore change the serialized table,
    # rather than disappearing beneath the writer's six-decimal precision.
    for study_index, count in enumerate((200000, 400000, 600003)):
        for condition_index, condition in enumerate(CONDITIONS):
            for site_index, site in enumerate(SITES):
                region = 2 * condition_index + site_index
                for experiment, offset in enumerate((-1, 1)):
                    aerobe_inc = count + 50000 * region + offset
                    aerobe_dec = (
                        aerobe_inc + 150000
                        if study_index == 1
                        else count // 2 + 30000 + 50000 * region - offset
                    )
                    anaerobe_inc = (800000, 550000, 300000)[study_index] + 25000 * region + offset
                    anaerobe_dec = anaerobe_inc + (-120000 if study_index == 1 else 120000)
                    for direction, aerobe, anaerobe in (
                        ("inc", aerobe_inc, anaerobe_inc),
                        ("dec", aerobe_dec, anaerobe_dec),
                    ):
                        profiles.append(
                            Profile(
                                study_id=f"study-{study_index}",
                                experiment_id=f"{condition}-{site}-{experiment}",
                                inc_dec=direction,
                                condition=condition,
                                body_site=site,
                                guild_counts={"aerobe": aerobe, "anaerobe": anaerobe},
                                total=2000000,
                            )
                        )
    assert len(profiles) == 48
    matrix = guild_maps.build_matrix(
        profiles,
        max_conditions=2,
        max_body_sites=2,
        min_studies=3,
        n_clusters=2,
        include_other_conditions=False,
        include_other_body_sites=False,
    )
    correlation = guild_maps.compute_correlation(profiles)
    written = guild_maps.write_outputs(matrix, correlation, output, "synthetic-2322")
    assert {path.name for path in written} == EXPECTED_FILES


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows, f"{path.name} must contain data rows"
    return rows


def _check_population(root: Path) -> None:
    """Require all configured matrix cells and every correlation comparison."""
    columns = [f"{site}|{direction}" for site in SITES for direction in DIRECTION_KEYS]
    for suffix in ("matrix", "study_counts", "cell_status"):
        rows = _rows(root / f"guild_condition_body_site_{suffix}.tsv")
        assert len(rows) == len(CONDITIONS)
        assert {row["condition"] for row in rows} == set(CONDITIONS)
        assert list(rows[0]) == ["condition", *columns]
        for row in rows:
            for column in columns:
                value = row[column]
                if suffix == "matrix":
                    assert value
                    assert math.isfinite(float(value))
                    assert 0 < float(value) < 1
                elif suffix == "study_counts":
                    assert value == "3"
                else:
                    assert value == "value"
    clusters = _rows(root / "guild_condition_row_clusters.tsv")
    assert len(clusters) == len(CONDITIONS)
    assert {row["condition"] for row in clusters} == set(CONDITIONS)
    assert [int(row["row_order"]) for row in clusters] == list(range(1, len(CONDITIONS) + 1))
    assert all(int(row["cluster"]) > 0 for row in clusters)
    correlations = _rows(root / "guild_aerobe_anaerobe_correlation.tsv")
    expected_pairs = {
        f"{left}~{right}" for left, right in itertools.combinations(DIRECTION_KEYS, 2)
    }
    expected_pairs.add("aerobe_contrast~anaerobe_contrast")
    expected_statistics = {
        *(("studies_shifted", key) for key in DIRECTION_KEYS),
        ("discordant_studies", "aerobe_up_anaerobe_down"),
        ("discordant_studies", "anaerobe_up_aerobe_down"),
        ("mcnemar_odds_ratio", "aerobe_vs_anaerobe_paired_signs (tested)"),
        *(("spearman_rho_across_studies", pair) for pair in expected_pairs),
    }
    assert len(correlations) == len(expected_statistics)
    assert {(row["statistic"], row["comparison"]) for row in correlations} == expected_statistics
    for row in correlations:
        assert row["value"]
        assert math.isfinite(float(row["value"]))
        assert row["n"] == "3"
        assert row["bugsigdb_vintage"] == "synthetic-2322"
        if row["statistic"] in {"mcnemar_odds_ratio", "spearman_rho_across_studies"}:
            assert row["p_value"]
            assert 0 <= float(row["p_value"]) <= 1
            assert row["p_value_type"] == "nominal"
            assert row["n_tests"] == "8"


@pytest.fixture(scope="module")
def seeded_guild_outputs(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Path]:
    root = tmp_path_factory.mktemp("guild-map-hash-seeds")
    outputs = (root / "seed-1", root / "seed-2")
    for seed, output in zip((1, 2), outputs, strict=True):
        result = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), str(output)],
            env={**os.environ, "PYTHONHASHSEED": str(seed)},
            capture_output=True,
            text=True,
            check=False,
            timeout=30,
        )
        assert result.returncode == 0, result.stdout + result.stderr
        assert "GUILD_MAPS_FILE=" in result.stdout
        print(f"seed={seed}: {result.stdout.strip()}", flush=True)
    return outputs


def test_guild_summary_tables_are_hash_seed_independent(
    seeded_guild_outputs: tuple[Path, Path],
) -> None:
    first, second = seeded_guild_outputs
    for root in (first, second):
        files = {path.relative_to(root).as_posix() for path in root.rglob("*") if path.is_file()}
        assert files == EXPECTED_FILES
        _check_population(root)
    differences = []
    for filename in sorted(EXPECTED_FILES):
        left, right = first / filename, second / filename
        if left.read_bytes() != right.read_bytes():
            differing_cells = [
                {"row": index, "column": key, "seed_1": a[key], "seed_2": b[key]}
                for index, (a, b) in enumerate(zip(_rows(left), _rows(right), strict=True), start=2)
                for key in a
                if a[key] != b[key]
            ]
            differences.append(f"{filename}: {differing_cells[:10]}")
    assert not differences, "\n".join(differences)


@pytest.mark.parametrize("filename", sorted(EXPECTED_FILES))
@pytest.mark.parametrize("damage", ["missing", "header-only"])
def test_incomplete_guild_outputs_fail_even_in_both_seeds(
    seeded_guild_outputs: tuple[Path, Path], tmp_path: Path, filename: str, damage: str
) -> None:
    """Equal file sets/bytes cannot make missing or unpopulated tables pass."""
    copies = (tmp_path / "seed-1", tmp_path / "seed-2")
    for source, copy in zip(seeded_guild_outputs, copies, strict=True):
        shutil.copytree(source, copy)
        path = copy / filename
        if damage == "missing":
            path.unlink()
        else:
            path.write_bytes(path.read_bytes().splitlines(keepends=True)[0])
    with pytest.raises(AssertionError):
        test_guild_summary_tables_are_hash_seed_independent(copies)


if __name__ == "__main__":
    print(f"GUILD_MAPS_FILE={guild_maps.__file__}", flush=True)
    write_guild_fixture(Path(sys.argv[1]))
