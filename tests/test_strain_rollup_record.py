"""The strain rollup records of the older pipeline's tables and the canonical check (#847)."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest

from mtp.strain_rollup_record import (
    StrainRollupRecordError,
    check_committed,
    content_sha256,
    main,
    record_committed,
    record_path,
    recorded_strain_rollup,
    require_declared_strain_rollup,
    require_strain_rollup,
    restore_committed,
    shared_strain_rollup,
    write_strain_rollup_record,
)
from mtp.taxa_trait_matrix import StrainRollup

OFF_TABLE = "study_id\texperiment_id\tinc_dec\tMETPO:1\nbsdb:1\tbsdb:1-1\tinc\t2\n"
PHENOTYPE_TABLE = (
    "study_id\texperiment_id\tinc_dec\tMETPO:1\tMETPO:2\nbsdb:1\tbsdb:1-1\tinc\t2\t1\n"
)


def _table(path: Path, content: str, scope: StrainRollup | None) -> Path:
    path.write_text(content, encoding="utf-8")
    if scope is not None:
        write_strain_rollup_record(path, scope)
    return path


def _committed(directory: Path, content: str = OFF_TABLE, scope: StrainRollup = "off") -> Path:
    """data/derived with a committed canonical table of ``scope`` and no local tables."""
    directory.mkdir(exist_ok=True)
    committed = directory / "study_trait_profiles_wide.tsv.gz"
    committed.write_bytes(gzip.compress(content.encode("utf-8"), mtime=0))
    write_strain_rollup_record(committed, scope)
    return directory / "study_trait_profiles_wide.tsv"


def test_a_record_binds_the_scope_to_the_table_content(tmp_path: Path) -> None:
    table = _table(tmp_path / "taxa_trait_matrix.tsv", PHENOTYPE_TABLE, "phenotype")

    assert record_path(table).name == "taxa_trait_matrix.tsv.strain_rollup.json"
    assert json.loads(record_path(table).read_text(encoding="utf-8")) == {
        "strain_rollup": "phenotype",
        "content_sha256": content_sha256(table),
    }
    assert recorded_strain_rollup(table) == require_strain_rollup(table) == "phenotype"


def test_a_gzip_table_is_bound_to_its_decompressed_content(tmp_path: Path) -> None:
    table = _committed(tmp_path / "derived")
    committed = table.with_name(table.name + ".gz")
    _table(table, OFF_TABLE, None)

    assert content_sha256(committed) == content_sha256(table)
    assert recorded_strain_rollup(committed) == "off"


def test_a_table_without_a_record_has_no_recorded_scope(tmp_path: Path) -> None:
    table = _table(tmp_path / "taxa_trait_matrix.tsv", OFF_TABLE, None)

    assert recorded_strain_rollup(table) is None
    with pytest.raises(StrainRollupRecordError, match="records no strain rollup scope"):
        require_strain_rollup(table)


def test_a_table_replaced_without_its_record_is_refused(tmp_path: Path) -> None:
    """The committed off table decompressed over a phenotype table keeps the phenotype
    record beside it; reading that pair under phenotype would be wrong."""
    table = _table(tmp_path / "study_trait_profiles_wide.tsv", PHENOTYPE_TABLE, "phenotype")
    table.write_text(OFF_TABLE, encoding="utf-8")

    with pytest.raises(StrainRollupRecordError, match="not the table its record"):
        recorded_strain_rollup(table)


@pytest.mark.parametrize(
    "record",
    ["not json", "[]", '{"content_sha256": "x"}', '{"strain_rollup": "strains"}'],
)
def test_an_invalid_record_is_refused(tmp_path: Path, record: str) -> None:
    table = _table(tmp_path / "taxa_trait_matrix.tsv", OFF_TABLE, None)
    record_path(table).write_text(record, encoding="utf-8")

    with pytest.raises(StrainRollupRecordError, match="not a valid strain rollup record"):
        recorded_strain_rollup(table)


def test_tables_read_together_must_not_record_different_scopes(tmp_path: Path) -> None:
    matrix = _table(tmp_path / "matrix.tsv", PHENOTYPE_TABLE, "phenotype")
    profiles = _table(tmp_path / "profiles.tsv", OFF_TABLE, "phenotype")
    subset = _table(tmp_path / "subset.tsv", OFF_TABLE, None)
    stale = _table(tmp_path / "other.tsv", OFF_TABLE, "off")

    assert shared_strain_rollup(matrix, profiles, subset) == "phenotype"
    assert shared_strain_rollup(subset) is None
    with pytest.raises(StrainRollupRecordError, match="different strain rollup scopes"):
        shared_strain_rollup(matrix, stale)


def test_a_step_refuses_recorded_tables_of_another_scope(tmp_path: Path) -> None:
    tables = [_table(tmp_path / "matrix.tsv", PHENOTYPE_TABLE, "phenotype")]

    require_declared_strain_rollup(tables, "phenotype")
    require_declared_strain_rollup([_table(tmp_path / "core.tsv", OFF_TABLE, None)], "off")
    with pytest.raises(StrainRollupRecordError, match="not the declared off"):
        require_declared_strain_rollup(tables, "off")


def test_the_committed_table_restored_with_its_record_passes_the_check(tmp_path: Path) -> None:
    table = _committed(tmp_path / "derived")
    matrix = _table(table.with_name("taxa_trait_matrix.tsv"), OFF_TABLE, "off")

    assert restore_committed(table, matrix) == "off"
    assert table.read_text(encoding="utf-8") == OFF_TABLE
    assert recorded_strain_rollup(table) == "off"
    assert check_committed(table, matrix) == []


def test_restoring_next_to_a_matrix_of_another_scope_is_refused(tmp_path: Path) -> None:
    """The old hint decompressed the committed off table next to a phenotype matrix."""
    table = _committed(tmp_path / "derived")
    matrix = _table(table.with_name("taxa_trait_matrix.tsv"), PHENOTYPE_TABLE, "phenotype")

    with pytest.raises(StrainRollupRecordError, match="off profile table next to the phenotype"):
        restore_committed(table, matrix)
    assert not table.exists()
    unrecorded = _table(table.with_name("unrecorded_matrix.tsv"), OFF_TABLE, None)
    with pytest.raises(StrainRollupRecordError, match="records no strain rollup scope"):
        restore_committed(table, unrecorded)


def test_the_check_names_both_scopes_and_offers_no_mixed_restore(tmp_path: Path) -> None:
    """A phenotype rebuild against the committed off table: rebuild or refresh, never
    restore beside the phenotype matrix (#847)."""
    table = _committed(tmp_path / "derived")
    matrix = _table(table.with_name("taxa_trait_matrix.tsv"), PHENOTYPE_TABLE, "phenotype")
    _table(table, PHENOTYPE_TABLE, "phenotype")

    messages = check_committed(table, matrix)

    assert messages[0] == f"{table} differs from the committed {table}.gz (#656)."
    assert (
        "built under the phenotype strain rollup scope, the committed table under off"
        in (messages[1])
    )
    assert messages[2] == "  Unintended: run just study-trait-profiles"
    assert "off profile table next to the phenotype matrix" in messages[3]
    assert messages[4].endswith(f"commit {table}.gz and {table}.gz.strain_rollup.json")
    assert not any("gzip -dc" in line or "-restore" in line for line in messages)


def test_the_check_offers_the_restore_only_when_the_scopes_match(tmp_path: Path) -> None:
    table = _committed(tmp_path / "derived")
    matrix = _table(table.with_name("taxa_trait_matrix.tsv"), PHENOTYPE_TABLE, "off")

    assert check_committed(table, matrix) == [
        f"{table} is missing or empty; run just study-trait-profiles, or "
        "just study-trait-profiles-restore"
    ]
    _table(table, PHENOTYPE_TABLE, "off")
    assert check_committed(table, matrix)[1] == (
        "  Unintended: run just study-trait-profiles, or just study-trait-profiles-restore"
    )


def test_the_check_refuses_a_local_table_without_a_valid_record(tmp_path: Path) -> None:
    table = _committed(tmp_path / "derived")
    matrix = table.with_name("taxa_trait_matrix.tsv")
    _table(table, OFF_TABLE, None)

    assert "records no strain rollup scope" in check_committed(table, matrix)[0]
    write_strain_rollup_record(table, "phenotype")
    table.write_text(OFF_TABLE + "bsdb:2\tbsdb:2-1\tdec\t1\n", encoding="utf-8")
    assert "not the table its record" in check_committed(table, matrix)[0]


def test_the_check_reports_a_committed_table_or_record_that_is_unusable(tmp_path: Path) -> None:
    table = _committed(tmp_path / "derived")
    matrix = table.with_name("taxa_trait_matrix.tsv")
    committed = table.with_name(table.name + ".gz")
    record_path(committed).unlink()

    assert "records no strain rollup scope" in check_committed(table, matrix)[0]
    committed.write_text("not gzip", encoding="utf-8")
    assert "is missing, empty or not gzip" in check_committed(table, matrix)[0]


def test_the_refresh_gives_the_committed_table_the_local_record(tmp_path: Path) -> None:
    table = _committed(tmp_path / "derived")
    committed = table.with_name(table.name + ".gz")
    _table(table, PHENOTYPE_TABLE, "phenotype")

    with pytest.raises(StrainRollupRecordError, match="does not decompress to"):
        record_committed(table)
    committed.write_bytes(gzip.compress(PHENOTYPE_TABLE.encode("utf-8"), mtime=0))
    record_committed(table)
    assert recorded_strain_rollup(committed) == "phenotype"
    assert check_committed(table, table.with_name("taxa_trait_matrix.tsv")) == []


def test_the_cli_exit_codes(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    table = _committed(tmp_path / "derived")
    matrix = table.with_name("taxa_trait_matrix.tsv")

    assert main(["check", str(table), "--matrix", str(matrix)]) == 1
    assert "is missing or empty" in capsys.readouterr().err
    assert main(["restore", str(table), "--matrix", str(matrix)]) == 0
    assert main(["check", str(table), "--matrix", str(matrix)]) == 0
    assert main(["show", str(table)]) == 0
    assert capsys.readouterr().out.splitlines()[-1] == "off"
    assert main(["show", str(matrix)]) == 1
