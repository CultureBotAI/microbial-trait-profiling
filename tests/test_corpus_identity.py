"""Which corpus a run read, recorded so a number can be traced to it (#222)."""

import json
from pathlib import Path

import pytest

from mtp.corpus_identity import (
    UNKNOWN_RELEASE,
    corpus_identity,
    measure_matrix,
    read_kg_inputs,
    summary_lines,
    write_identity,
)


def _matrix(path: Path, *, mediadive: int, retained: int, rows: int) -> Path:
    header = ["study_id", "experiment_id", "inc_dec"]
    header += [f"mediadive.medium:{i}" for i in range(mediadive)]
    header += [f"METPO:{i}" for i in range(retained)]
    lines = ["\t".join(header)]
    counts = ["0"] * (mediadive + retained)
    for row in range(rows):
        lines.append("\t".join([f"bsdb:{row}", f"bsdb:{row}-1", "inc", *counts]))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def test_the_three_column_groups_account_for_the_whole_header(tmp_path: Path) -> None:
    identity = measure_matrix(_matrix(tmp_path / "m.tsv", mediadive=4, retained=7, rows=3))

    assert identity.identifier_columns == 3
    assert identity.mediadive_columns == 4
    assert identity.retained_trait_columns == 7
    assert identity.header_columns == 14
    assert identity.profile_rows == 3


def test_two_matrices_differing_only_in_mediadive_columns_are_distinguishable(
    tmp_path: Path,
) -> None:
    """The 2,290-vs-2,333 case: same retained count, different matrix (#222)."""
    first = measure_matrix(_matrix(tmp_path / "a.tsv", mediadive=4, retained=7, rows=2))
    second = measure_matrix(_matrix(tmp_path / "b.tsv", mediadive=6, retained=7, rows=2))

    assert first.retained_trait_columns == second.retained_trait_columns
    # Identical retained counts, so only the full record tells them apart.
    assert first.mediadive_columns != second.mediadive_columns
    assert first.header_columns != second.header_columns
    assert first.sha256 != second.sha256


def test_a_header_that_does_not_add_up_raises(tmp_path: Path) -> None:
    identity = measure_matrix(_matrix(tmp_path / "m.tsv", mediadive=2, retained=2, rows=1))
    broken = type(identity)(**{**identity.__dict__, "retained_trait_columns": 99})

    with pytest.raises(ValueError, match="not the"):
        broken.check()


def test_a_missing_release_tag_is_recorded_rather_than_omitted(tmp_path: Path) -> None:
    """ "Not known" and "not applicable" differ, so the record says which (#222)."""
    sidecar = tmp_path / "provenance.json"
    sidecar.write_text(
        json.dumps({"files": {"merged-kg_nodes.tsv": {"md5": "abc", "size": 1}}}),
        encoding="utf-8",
    )

    release, inputs = read_kg_inputs(sidecar)

    assert release == UNKNOWN_RELEASE
    assert [item.name for item in inputs] == ["merged-kg_nodes.tsv"]
    assert inputs[0].declared_md5 == "abc"
    # The sidecar names a file that is not beside it, and the record says so.
    assert inputs[0].present is False


def test_a_declared_release_tag_is_used_when_present(tmp_path: Path) -> None:
    sidecar = tmp_path / "provenance.json"
    sidecar.write_text(json.dumps({"release": "kg-microbe-2026-04", "files": {}}), encoding="utf-8")

    release, _ = read_kg_inputs(sidecar)

    assert release == "kg-microbe-2026-04"


def test_an_absent_or_malformed_sidecar_does_not_fail_the_run(tmp_path: Path) -> None:
    """A recipe must still record the matrix when the KG sidecar is unusable."""
    absent, inputs = read_kg_inputs(tmp_path / "nope.json")
    assert UNKNOWN_RELEASE in absent
    assert inputs == []

    broken = tmp_path / "broken.json"
    broken.write_text("not json", encoding="utf-8")
    release, inputs = read_kg_inputs(broken)
    assert UNKNOWN_RELEASE in release
    assert inputs == []


def test_the_written_record_names_both_sides_of_the_join(tmp_path: Path) -> None:
    matrix = _matrix(tmp_path / "m.tsv", mediadive=1, retained=2, rows=1)
    sidecar = tmp_path / "provenance.json"
    sidecar.write_text(
        json.dumps({"files": {"merged-kg_nodes.tsv": {"md5": "abc"}}}), encoding="utf-8"
    )

    identity = corpus_identity(matrix, bugsigdb_vintage="v", kg_provenance_path=sidecar)
    path = write_identity(identity, tmp_path / "out", code_commit="deadbeef")
    record = json.loads(path.read_text(encoding="utf-8"))

    assert record["matrix"]["sha256"] == identity.matrix.sha256
    assert record["bugsigdb_vintage"] == "v"
    assert record["kg_microbe"]["inputs"][0]["declared_md5"] == "abc"
    assert record["code_commit"] == "deadbeef"
    # A vintage alone does not identify a matrix, and the record says so.
    assert "join" in record["note"]


def test_the_printed_summary_states_the_matrix_and_both_inputs(tmp_path: Path) -> None:
    matrix = _matrix(tmp_path / "m.tsv", mediadive=1, retained=2, rows=1)
    sidecar = tmp_path / "provenance.json"
    sidecar.write_text(
        json.dumps({"files": {"merged-kg_nodes.tsv": {"md5": "abc"}}}), encoding="utf-8"
    )

    text = "\n".join(
        summary_lines(corpus_identity(matrix, bugsigdb_vintage="v", kg_provenance_path=sidecar))
    )

    assert "sha256" in text
    assert "1 MediaDive" in text
    assert "2 retained" in text
    assert "merged-kg_nodes.tsv" in text
