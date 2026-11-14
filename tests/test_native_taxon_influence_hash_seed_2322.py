"""Native summary table payloads must be independent of Python's hash seed (#2322)."""

from __future__ import annotations

import csv
import gzip
import io
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

import mtp.native_taxon_influence as native
from mtp.contribution_inputs import Base, Key, Pair, digest
from mtp.taxon_contributions import Record, write_table

COLUMNS = tuple(f"METPO:fixture{index}" for index in range(8))
DISPLAY_TRAITS = COLUMNS[:6]
SCENARIO_REMOVALS = {
    "exclude_genus_NCBITaxon:901": ("NCBITaxon:901", {COLUMNS[0]: 2, COLUMNS[1]: 1}),
    "exclude_family_NCBITaxon:902": ("NCBITaxon:902", {COLUMNS[2]: 1, COLUMNS[3]: 3}),
}
SCENARIOS = ("baseline", *SCENARIO_REMOVALS)
SCOPE_CONDITIONS = {"pooled": "*", "condition_a": "condition_a"}
SCOPE_PROFILES = {"pooled": 24, "condition_a": 12}
REFERENCES = ("refitted", "fixed_baseline_sd")
BACKGROUNDS = ("fixed_baseline", "recomputed_remaining")
NATIVE_FILES = {
    "native_composition.tsv.gz",
    "native_composition_comparison.tsv",
    "native_ora_fold_comparison.tsv",
    "native_ora_families.tsv",
    "native_summary.md",
    *(f"{scenario.replace(':', '_')}/native_ora_profiles.tsv.gz" for scenario in SCENARIOS),
}
INPUT_FILES = {
    "scenarios.tsv",
    *(
        f"{scenario.replace(':', '_')}/{filename}"
        for scenario in SCENARIOS
        for filename in ("removed_source_counts.tsv.gz", "profile_counts.tsv.gz")
    ),
}


def _payload(path: Path) -> bytes:
    data = path.read_bytes()
    return gzip.decompress(data) if path.suffix == ".gz" else data


def _table(path: Path) -> list[dict[str, str]]:
    return list(csv.DictReader(io.StringIO(_payload(path).decode("utf-8")), delimiter="\t"))


def _base() -> Base:
    """Six studies and two experiments each give unequal native and study-weighted vectors."""
    profiles: dict[Key, dict[str, int]] = {}
    metadata: dict[Pair, dict[str, str]] = {}
    reported: dict[Key, dict[str, set[str]]] = {}
    for study in range(6):
        for experiment in range(2):
            pair = (f"study-{study}", f"experiment-{experiment}")
            metadata[pair] = {
                "condition": "condition_a" if study % 2 == 0 else "condition_b",
                "body_site": "site",
            }
            for direction_index, direction in enumerate(("inc", "dec")):
                key = (*pair, direction)
                profiles[key] = {
                    trait: 5
                    + (study * 11 + experiment * 17 + direction_index * 23 + index * (3 + study))
                    % (41 + index * 13)
                    for index, trait in enumerate(COLUMNS)
                }
                reported[key] = {
                    taxon: {f"signature-{study}-{experiment}-{direction}"}
                    for taxon, present in (
                        ("NCBITaxon:901", (study + experiment + direction_index) % 2 == 0),
                        ("NCBITaxon:902", (study + 2 * experiment + direction_index) % 3 != 0),
                    )
                    if present
                }
    return Base(
        list(COLUMNS),
        {f"display-{index}": trait for index, trait in enumerate(DISPLAY_TRAITS)},
        profiles,
        metadata,
        reported,
        {},
        {},
        {},
        [],
    )


def write_native_fixture(output: Path) -> None:
    """Call the real native exporter with a small, internally checked sparse-removal audit."""
    base = _base()
    output.mkdir()
    write_table(
        output / "scenarios.tsv",
        ["scenario", "status"],
        [
            *({"scenario": scenario, "status": "run"} for scenario in SCENARIO_REMOVALS),
            {"scenario": "exclude_genus_not_selected", "status": "outside_rule_declared_panel"},
        ],
    )
    for scenario in SCENARIOS:
        directory = output / scenario.replace(":", "_")
        directory.mkdir()
        target, removal = SCENARIO_REMOVALS.get(scenario, ("", {}))
        write_table(
            directory / "removed_source_counts.tsv.gz",
            ["reported_taxon", "trait_key", "removed_count"],
            [
                {"reported_taxon": target, "trait_key": trait, "removed_count": count}
                for trait, count in removal.items()
            ],
        )
        entries: list[Record] = []
        for key, original in base.profiles.items():
            remaining = {
                trait: count - (removal.get(trait, 0) if target in base.reported[key] else 0)
                for trait, count in original.items()
            }
            assert min(remaining.values()) > 0
            entries.extend(
                {
                    "study_id": key[0],
                    "experiment_id": key[1],
                    "direction": key[2],
                    "trait_key": trait,
                    "original_count": original[trait],
                    "remaining_count": remaining[trait],
                    "original_total": sum(original.values()),
                    "remaining_total": sum(remaining.values()),
                }
                for trait in DISPLAY_TRAITS
            )
        write_table(directory / "profile_counts.tsv.gz", list(entries[0]), entries)
    before = {
        path.relative_to(output).as_posix(): digest(path)
        for path in output.rglob("*")
        if path.is_file()
    }
    assert set(before) == INPUT_FILES
    result = native.export_native(
        base,
        {
            "protocol": {
                "scopes": [
                    {"id": scope, "condition": condition, "body_site": "*"}
                    for scope, condition in SCOPE_CONDITIONS.items()
                ]
            }
        },
        {
            "version": "190-native-v3",
            "composition": {
                "minimum_profile_support": 1,
                "maximum_profile_support_fraction": 1.0,
                "minimum_pairs": 3,
                "transform_multiplier": 1000,
                "family_size_field": "family_selected_traits",
                "untestable_correction_p": 1.0,
            },
            "ora": {"backgrounds": list(BACKGROUNDS)},
        },
        output,
    )
    assert result["scenarios"] == len(SCENARIOS)
    assert result["scopes"] == len(SCOPE_CONDITIONS)
    assert set(result["output_sha256"]) == NATIVE_FILES
    assert set(result["input_sha256"]) == INPUT_FILES
    assert all(digest(output / relative) == checksum for relative, checksum in before.items())


@pytest.fixture(scope="module")
def seeded_outputs(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Path]:
    root = tmp_path_factory.mktemp("native-hash-seeds")
    outputs = (root / "seed-1", root / "seed-2")
    for seed, output in zip((1, 2), outputs, strict=True):
        result = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), str(output)],
            env={**os.environ, "PYTHONHASHSEED": str(seed)},
            capture_output=True,
            text=True,
            check=False,
            timeout=60,
        )
        assert result.returncode == 0, result.stdout + result.stderr
        identity = next(
            line for line in result.stdout.splitlines() if line.startswith("NATIVE_MODULE_FILE=")
        )
        print(f"seed={seed}: {identity}", flush=True)
    return outputs


def _assert_populated_native_outputs(root: Path) -> None:
    actual = {path.relative_to(root).as_posix() for path in root.rglob("*") if path.is_file()}
    assert actual == INPUT_FILES | NATIVE_FILES
    for relative in NATIVE_FILES:
        assert _payload(root / relative), relative
        if relative.endswith((".tsv", ".tsv.gz")):
            assert _table(root / relative), relative

    stats = _table(root / "native_composition.tsv.gz")
    expected_stats = {
        (scenario, scope, trait)
        for scenario in SCENARIOS
        for scope in SCOPE_CONDITIONS
        for trait in COLUMNS
    }
    assert len(stats) == len(expected_stats)
    assert {(row["scenario"], row["scope"], row["trait_key"]) for row in stats} == expected_stats
    assert all(row["effect"] and row["reference_sd"] for row in stats)
    assert any(row["status"] == "tested" and row["p"] and row["q"] for row in stats)

    for filename, dimension, values in (
        ("native_composition_comparison.tsv", "reference", REFERENCES),
        ("native_ora_fold_comparison.tsv", "background", BACKGROUNDS),
    ):
        records = _table(root / filename)
        expected = {
            (scenario, scope, trait, value)
            for scenario in SCENARIOS
            for scope in SCOPE_CONDITIONS
            for trait in DISPLAY_TRAITS
            for value in values
        }
        assert len(records) == len(expected), filename
        assert {
            (row["scenario"], row["scope"], row["trait_key"], row[dimension]) for row in records
        } == expected, filename
        assert all(row["before_effect"] and row["after_effect"] for row in records), filename
        assert all(int(row["before_studies"]) >= 3 for row in records), filename
        assert any(
            float(row["matched_change_effect"]) != 0
            for row in records
            if row["scenario"] != "baseline"
        ), filename

    families = _table(root / "native_ora_families.tsv")
    expected_families = {
        (scenario, scope, background)
        for scenario in SCENARIOS
        for scope in SCOPE_CONDITIONS
        for background in BACKGROUNDS
    }
    assert len(families) == len(expected_families)
    assert {
        (row["scenario"], row["scope"], row["background"]) for row in families
    } == expected_families
    assert all(
        int(row[f"{tail}_tests"]) == SCOPE_PROFILES[row["scope"]] * len(COLUMNS)
        for row in families
        for tail in ("enrichment", "depletion")
    )

    base = _base()
    for scenario in SCENARIOS:
        relative = f"{scenario.replace(':', '_')}/native_ora_profiles.tsv.gz"
        records = _table(root / relative)
        expected_cells = {
            (scope, background, *key, trait)
            for scope, condition in SCOPE_CONDITIONS.items()
            for background in BACKGROUNDS
            for key in base.profiles
            if condition == "*" or base.metadata[key[:2]]["condition"] == condition
            for trait in DISPLAY_TRAITS
        }
        assert len(records) == len(expected_cells), relative
        assert {
            (
                row["scope"],
                row["background"],
                row["study_id"],
                row["experiment_id"],
                row["direction"],
                row["trait_key"],
            )
            for row in records
        } == expected_cells, relative
        assert all(row["fold"] and row["enrichment_q"] and row["depletion_q"] for row in records), (
            relative
        )

    summary = (root / "native_summary.md").read_text(encoding="utf-8")
    assert "2 rule-declared exclusion scenarios" in summary
    assert all(f"| display-{index} |" in summary for index in range(len(DISPLAY_TRAITS)))


def test_every_native_output_is_hash_seed_independent(
    seeded_outputs: tuple[Path, Path],
) -> None:
    first, second = seeded_outputs
    for root in seeded_outputs:
        _assert_populated_native_outputs(root)
    differences: list[str] = []
    for relative in sorted(NATIVE_FILES):
        if _payload(first / relative) == _payload(second / relative):
            continue
        if relative.endswith((".tsv", ".tsv.gz")):
            left_rows, right_rows = (_table(root / relative) for root in seeded_outputs)
            cells = [
                {"row": index, "column": column, "seed_1": left[column], "seed_2": right[column]}
                for index, (left, right) in enumerate(zip(left_rows, right_rows, strict=True), 2)
                for column in left
                if left[column] != right[column]
            ]
            differences.append(f"{relative}: {cells[:10]}")
        else:
            differences.append(f"{relative}: Markdown bytes differ")
    assert not differences, "\n".join(differences)


@pytest.mark.parametrize("relative", sorted(NATIVE_FILES))
@pytest.mark.parametrize("corruption", ["missing", "empty"])
def test_missing_or_empty_native_output_fails(
    seeded_outputs: tuple[Path, Path], tmp_path: Path, relative: str, corruption: str
) -> None:
    """Both seeds losing the same export must fail independently of byte comparison."""
    copies = (tmp_path / "seed-1", tmp_path / "seed-2")
    for source, copy in zip(seeded_outputs, copies, strict=True):
        shutil.copytree(source, copy)
        path = copy / relative
        if corruption == "missing":
            path.unlink()
        else:
            data = _payload(path).splitlines(keepends=True)[0] if ".tsv" in relative else b""
            path.write_bytes(gzip.compress(data, mtime=0) if path.suffix == ".gz" else data)
    with pytest.raises(AssertionError):
        test_every_native_output_is_hash_seed_independent(copies)


if __name__ == "__main__":
    print(f"NATIVE_MODULE_FILE={native.__file__}", flush=True)
    write_native_fixture(Path(sys.argv[1]))
