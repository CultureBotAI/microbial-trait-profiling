"""The strain rollup records of the older pipeline's tables and the canonical check (#847)."""

from __future__ import annotations

import gzip
import json
import os
import shlex
import subprocess
from pathlib import Path

import pytest

from mtp import strain_rollup_record
from mtp.paths import REPO_ROOT
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
from mtp.workflow.arms.runtime import INPUT_NAMES

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


def test_tables_read_together_record_one_scope_or_none(tmp_path: Path) -> None:
    """A table without a record is never read with recorded ones (#939)."""
    matrix = _table(tmp_path / "matrix.tsv", PHENOTYPE_TABLE, "phenotype")
    profiles = _table(tmp_path / "profiles.tsv", OFF_TABLE, "phenotype")
    unrecorded = _table(tmp_path / "unrecorded.tsv", OFF_TABLE, None)
    core = _table(tmp_path / "core.tsv", PHENOTYPE_TABLE, None)
    other = _table(tmp_path / "other.tsv", OFF_TABLE, "off")

    assert shared_strain_rollup(matrix, profiles) == "phenotype"
    assert shared_strain_rollup(unrecorded, core) is None
    with pytest.raises(StrainRollupRecordError, match=r"unrecorded\.tsv \(no record\)"):
        shared_strain_rollup(matrix, unrecorded)
    with pytest.raises(StrainRollupRecordError, match="record different strain rollup scopes"):
        shared_strain_rollup(matrix, other)
    other.write_text(PHENOTYPE_TABLE, encoding="utf-8")
    with pytest.raises(StrainRollupRecordError, match="rebuild them with the steps that wrote"):
        shared_strain_rollup(matrix, other)


def test_tables_in_data_derived_must_all_have_records(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The older pipeline records its tables, so two of them without records at their
    paths in data/derived, whose scopes the old restore hint may have mixed, do not pass
    together; elsewhere, such as in a core run, they are not compared (#945)."""
    derived = tmp_path / "data" / "derived"
    derived.mkdir(parents=True)
    monkeypatch.setattr(strain_rollup_record, "DERIVED_DATA_DIR", derived)
    matrix = _table(derived / "taxa_trait_matrix.tsv", PHENOTYPE_TABLE, None)
    profiles = _table(derived / "study_trait_profiles_wide.tsv", OFF_TABLE, None)

    with pytest.raises(
        StrainRollupRecordError,
        match="older pipeline's tables in data/derived must; rebuild them: just "
        "taxa-trait-matrix, then just study-trait-profiles",
    ):
        shared_strain_rollup(matrix, profiles)
    for table in (matrix, profiles):
        write_strain_rollup_record(table, "phenotype")
    assert shared_strain_rollup(matrix, profiles) == "phenotype"
    core = [_table(tmp_path / name, OFF_TABLE, None) for name in ("profiles.tsv", "matrix.tsv")]
    assert shared_strain_rollup(*core) is None


def test_subsets_and_arm_inputs_in_data_derived_are_not_the_pipeline_tables(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A core consumer's subsets and an arm's prepared inputs may lie in data/derived;
    only the older pipeline's own tables there must have records (#978)."""
    derived = tmp_path / "data" / "derived"
    monkeypatch.setattr(strain_rollup_record, "DERIVED_DATA_DIR", derived)
    inputs = derived / "rerun" / "inputs"
    inputs.mkdir(parents=True)
    arm = [
        _table(inputs / INPUT_NAMES[product], OFF_TABLE, None)
        for product in ("ora_profiles", "taxa_matrix")
    ]
    subset = _table(derived / "rerun" / "study_trait_profiles_pooled.tsv", OFF_TABLE, None)

    assert shared_strain_rollup(*arm) is None
    assert shared_strain_rollup(subset, arm[1]) is None


def test_the_rebuild_advice_names_the_recipe_of_each_table(tmp_path: Path) -> None:
    """The ORA table is written by its own recipe, not by just study-trait-profiles (#978)."""
    ora = _table(tmp_path / "study_trait_profiles_ora_wide.tsv", OFF_TABLE, None)
    matrix = _table(tmp_path / "taxa_trait_matrix.tsv", PHENOTYPE_TABLE, "phenotype")

    with pytest.raises(
        StrainRollupRecordError,
        match="rebuild them: just taxa-trait-matrix, then just study-trait-profiles-ora",
    ):
        shared_strain_rollup(ora, matrix)


def test_a_link_to_a_recorded_table_reads_its_record(tmp_path: Path) -> None:
    table = _table(tmp_path / "taxa_trait_matrix.tsv", PHENOTYPE_TABLE, "phenotype")
    link = tmp_path / "elsewhere" / "matrix.tsv"
    link.parent.mkdir()
    link.symlink_to(table)

    assert record_path(link) == record_path(table)
    assert recorded_strain_rollup(link) == "phenotype"


def test_a_step_refuses_recorded_tables_of_another_scope(tmp_path: Path) -> None:
    tables = [_table(tmp_path / "matrix.tsv", PHENOTYPE_TABLE, "phenotype")]

    require_declared_strain_rollup(tables, "phenotype")
    require_declared_strain_rollup([_table(tmp_path / "core.tsv", OFF_TABLE, None)], "off")
    with pytest.raises(StrainRollupRecordError, match="not the declared off; declare phenotype"):
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


def test_a_restore_refuses_a_committed_table_its_record_does_not_describe(
    tmp_path: Path,
) -> None:
    """The restored bytes are checked against the committed record before they replace
    the local table, and nothing partial is left behind (#939)."""
    table = _committed(tmp_path / "derived")
    committed = table.with_name(table.name + ".gz")
    _table(table, PHENOTYPE_TABLE, None)
    committed.write_bytes(gzip.compress(PHENOTYPE_TABLE.encode("utf-8"), mtime=0))

    with pytest.raises(StrainRollupRecordError, match="not the table its record describes"):
        restore_committed(table, table.with_name("taxa_trait_matrix.tsv"))
    assert table.read_text(encoding="utf-8") == PHENOTYPE_TABLE
    assert sorted(path.name for path in table.parent.iterdir()) == [
        table.name,
        committed.name,
        f"{committed.name}.strain_rollup.json",
    ]


def test_a_scope_change_is_refreshed_or_rebuilt_under_the_committed_scope(
    tmp_path: Path,
) -> None:
    """A phenotype rebuild against the committed off table cannot pass by rebuilding
    again: the committed copy is refreshed, or the off table is rebuilt elsewhere; it is
    never restored beside the phenotype matrix (#847, #939)."""
    table = _committed(tmp_path / "derived")
    matrix = _table(table.with_name("taxa_trait_matrix.tsv"), PHENOTYPE_TABLE, "phenotype")
    _table(table, PHENOTYPE_TABLE, "phenotype")

    assert check_committed(table, matrix) == [
        f"{table} differs from the committed {table}.gz (#656).",
        "  It was built under the phenotype strain rollup scope, the committed table under off "
        "(#788, #847); a rebuild gives a phenotype table again, because data/derived holds "
        "only the primary phenotype matrix.",
        "  If the phenotype table is intended: just study-trait-profiles-refresh-gz, then "
        f"commit {table}.gz and {table}.gz.strain_rollup.json",
        "  To rebuild the committed off table instead, outside data/derived: "
        "just study-trait-profiles-reproduce-committed <directory>",
    ]


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
    _table(matrix, PHENOTYPE_TABLE, "phenotype")
    table.unlink()
    assert check_committed(table, matrix) == [
        f"{table} is missing or empty; run just study-trait-profiles",
        "  Not offering the committed table: restoring it would put the off profile table "
        f"next to the phenotype matrix {matrix}.",
    ]


def test_the_check_advises_by_what_a_table_without_a_valid_record_holds(
    tmp_path: Path,
) -> None:
    """The committed content without its record is restored, unless the matrix records
    another scope; another table is rebuilt (#939)."""
    table = _committed(tmp_path / "derived")
    matrix = table.with_name("taxa_trait_matrix.tsv")
    record = f"{table.name}.strain_rollup.json"
    _table(table, OFF_TABLE, None)

    assert check_committed(table, matrix) == [
        f"{table} has the committed content, but it has no strain rollup record {record}; "
        "restore it with its record: run just study-trait-profiles-restore (#847)"
    ]
    _table(matrix, PHENOTYPE_TABLE, "phenotype")
    messages = check_committed(table, matrix)
    assert messages[0].endswith("; run just study-trait-profiles (#847)")
    assert "off profile table next to the phenotype matrix" in messages[1]
    write_strain_rollup_record(table, "phenotype")
    table.write_text(OFF_TABLE + "bsdb:2\tbsdb:2-1\tdec\t1\n", encoding="utf-8")
    assert check_committed(table, matrix)[0] == (
        f"{table} is not the committed table: its record {record} is invalid or describes "
        "other content; run just study-trait-profiles (#847)"
    )


def test_the_check_reports_a_committed_table_or_record_that_is_unusable(tmp_path: Path) -> None:
    table = _committed(tmp_path / "derived")
    matrix = table.with_name("taxa_trait_matrix.tsv")
    committed = table.with_name(table.name + ".gz")
    record_path(committed).unlink()

    assert "records no strain rollup scope" in check_committed(table, matrix)[0]
    committed.write_text("not gzip", encoding="utf-8")
    assert "is not gzip" in check_committed(table, matrix)[0]
    committed.write_bytes(b"")
    assert check_committed(table, matrix)[0].startswith(f"{committed} is missing or empty;")


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


def test_a_stale_committed_record_offers_the_refresh(tmp_path: Path) -> None:
    """When the committed table was regenerated without its record, as in the rerun
    PR's merge order, refreshing records it; a checkout restores the stale pair (#978)."""
    table = _committed(tmp_path / "derived")
    committed = table.with_name(table.name + ".gz")
    committed.write_bytes(gzip.compress(PHENOTYPE_TABLE.encode("utf-8"), mtime=0))

    messages = check_committed(table, table.with_name("taxa_trait_matrix.tsv"))
    assert "not the table its record" in messages[0]
    assert messages[2] == (
        "  Or, if the committed table was regenerated from the local one, record it: "
        "just study-trait-profiles-refresh-gz"
    )


def test_an_off_table_after_a_phenotype_commit_is_rebuilt_in_data_derived(
    tmp_path: Path,
) -> None:
    """Once the committed table is phenotype, an off local table is rebuilt with the
    canonical recipe, not reproduced elsewhere (#978)."""
    table = _committed(tmp_path / "derived", PHENOTYPE_TABLE, "phenotype")
    matrix = _table(table.with_name("taxa_trait_matrix.tsv"), PHENOTYPE_TABLE, "phenotype")
    _table(table, OFF_TABLE, "off")

    assert check_committed(table, matrix)[1:] == [
        "  It was built under the off strain rollup scope, the committed table under "
        "phenotype (#788, #847).",
        "  Unintended: run just study-trait-profiles, which builds the primary phenotype "
        "table, or just study-trait-profiles-restore",
    ]


@pytest.mark.parametrize("name", ["out dir", 'quoted "dir', "dollar $HOME", "tick `id`"])
def test_reproduce_committed_keeps_its_directory_argument_literal(
    tmp_path: Path, name: str
) -> None:
    """The recipe passes its directory through unexpanded (#975). A stub uv records the
    builds instead of running them, so the final comparison fails as expected."""
    stubs = tmp_path / "bin"
    stubs.mkdir()
    calls = tmp_path / "uv_calls.txt"
    stub = stubs / "uv"
    stub.write_text(
        "#!/usr/bin/env bash\n"
        f"printf '%s\\n' \"$@\" >> {shlex.quote(str(calls))}\n"
        'if [ "$5" = show ]; then echo off; fi\n',
        encoding="utf-8",
    )
    stub.chmod(0o755)
    directory = tmp_path / name
    environment = {**os.environ, "PATH": f"{stubs}{os.pathsep}{os.environ['PATH']}"}

    completed = subprocess.run(
        ["just", "study-trait-profiles-reproduce-committed", str(directory)],
        cwd=REPO_ROOT,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 1, completed.stderr
    assert directory.is_dir()
    assert f"{directory}/taxa_trait_matrix.tsv" in calls.read_text(encoding="utf-8").splitlines()
    assert f"{directory}/study_trait_profiles_wide.tsv differs from the committed off" in (
        completed.stderr
    )
