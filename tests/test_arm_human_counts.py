"""Human-only legacy inputs retain values, repeats, unpaired rows, and zero profiles (#658)."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from mtp.paths import file_sha256
from mtp.workflow.arms.shared import select_human_count_profiles

HEADER = "study_id\texperiment_id\tinc_dec\tMETPO:1\tmediadive.medium:1\n"


def _sources(tmp_path: Path) -> tuple[Path, Path, Path, Path]:
    profiles = tmp_path / "all.tsv"
    profiles.write_text(
        HEADER
        + "human\te1\tinc\t0\t0\n"
        + "mouse\te1\tinc\t99\t88\n"
        + "human\te1\tinc\t1.0\t2\n"
        + "human\te2\tdec\t5\t6\n"
        + "unknown\te1\tdec\t7\t8\n"
    )
    experiments = tmp_path / "experiments.tsv"
    experiments.write_text(
        "study_id\texperiment_id\thost_species\n"
        "human\te1\tHomo sapiens\n"
        "mouse\te1\tMus musculus\n"
        "human\te2\tHomo sapiens\n"
        "unknown\te1\t\n"
    )
    return profiles, experiments, tmp_path / "human.tsv", tmp_path / "summary.json"


def test_exact_host_filter_preserves_all_human_values_and_order(tmp_path: Path) -> None:
    profiles, experiments, output, summary = _sources(tmp_path)
    select_human_count_profiles(profiles, experiments, output, summary)
    assert output.read_text() == (
        HEADER + "human\te1\tinc\t0\t0\n" + "human\te1\tinc\t1.0\t2\n" + "human\te2\tdec\t5\t6\n"
    )
    record = json.loads(summary.read_text())
    assert record["input_rows"] == 5
    assert record["retained_rows"] == 3
    assert record["excluded_rows"] == 2
    assert record["excluded_rows_by_host"] == {"Mus musculus": 1, "<unspecified>": 1}
    assert record["output"]["sha256"] == file_sha256(output)
    assert record["profiles"]["sha256"] == file_sha256(profiles)
    assert record["experiments"]["sha256"] == file_sha256(experiments)
    with output.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 3
    assert rows[0]["METPO:1"] == "0"
    assert rows[1]["mediadive.medium:1"] == "2"


@pytest.mark.parametrize("corruption", ["missing", "duplicate", "malformed", "nonhuman"])
def test_unaccountable_or_empty_cohort_fails_before_writes(tmp_path: Path, corruption: str) -> None:
    profiles, experiments, output, summary = _sources(tmp_path)
    if corruption == "missing":
        experiments.write_text(experiments.read_text().replace("human\te2\tHomo sapiens\n", ""))
    elif corruption == "duplicate":
        experiments.write_text(experiments.read_text() + "human\te1\tMus musculus\n")
    elif corruption == "malformed":
        profiles.write_text(profiles.read_text() + "human\te1\tinc\t1\n")
    else:
        experiments.write_text(experiments.read_text().replace("Homo sapiens", "Mus musculus"))
    with pytest.raises(ValueError, match=r"metadata|identity|header|Homo sapiens"):
        select_human_count_profiles(profiles, experiments, output, summary)
    assert not output.exists()
    assert not summary.exists()


def test_previous_output_is_never_replaced(tmp_path: Path) -> None:
    profiles, experiments, output, summary = _sources(tmp_path)
    output.write_text("previous result\n")
    with pytest.raises(FileExistsError):
        select_human_count_profiles(profiles, experiments, output, summary)
    assert output.read_text() == "previous result\n"
