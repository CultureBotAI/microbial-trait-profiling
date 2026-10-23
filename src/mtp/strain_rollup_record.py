"""The strain rollup scope an older-pipeline table was built under (#847).

The core workflow records its scope in each run (#788). The matrix and profile tables of
the older ``data/derived`` pipeline carry theirs in a record beside the table,
``<table name>.strain_rollup.json``, bound to the table's content. A table replaced
without its record, such as the committed profile table decompressed over a local one,
is therefore refused rather than read under the scope of the table it replaced.

A table without a record, such as a core product, a stratum subset or a table built
before #847, has no recorded scope to check. Readers that combine tables refuse only
recorded scopes that disagree.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import shutil
import sys
import zlib
from collections.abc import Sequence
from pathlib import Path

from mtp.paths import DERIVED_DATA_DIR
from mtp.taxa_trait_matrix import StrainRollup, parse_strain_rollup

RECORD_SUFFIX = ".strain_rollup.json"
DEFAULT_MATRIX_PATH = DERIVED_DATA_DIR / "taxa_trait_matrix.tsv"
DEFAULT_PROFILE_PATH = DERIVED_DATA_DIR / "study_trait_profiles_wide.tsv"
# What a provenance record says for a table built without a record
UNRECORDED = "unrecorded: the table has no strain rollup record (#847)"


class StrainRollupRecordError(ValueError):
    """A table's strain rollup record is missing, unreadable, stale or contradicted."""


def record_path(table: Path) -> Path:
    """The record beside a table, named after its whole file name."""
    return table.with_name(table.name + RECORD_SUFFIX)


def content_sha256(table: Path) -> str:
    """The SHA-256 of a table's content; a gzip table is hashed decompressed."""
    digest = hashlib.sha256()
    with gzip.open(table, "rb") if table.suffix == ".gz" else table.open("rb") as handle:
        while chunk := handle.read(1 << 20):
            digest.update(chunk)
    return digest.hexdigest()


def write_strain_rollup_record(table: Path, strain_rollup: StrainRollup) -> Path:
    """Record the scope a table was just built under, bound to its content."""
    path = record_path(table)
    record = {"strain_rollup": strain_rollup, "content_sha256": content_sha256(table)}
    path.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    return path


def recorded_strain_rollup(table: Path) -> StrainRollup | None:
    """The scope a table's record states, or None when the table has no record.

    A record that is unreadable, names no valid scope or describes other content is
    refused: the table was replaced without its record.
    """
    path = record_path(table)
    if not path.is_file():
        return None
    try:
        record = json.loads(path.read_text(encoding="utf-8"))
        scope = parse_strain_rollup(record["strain_rollup"])
        recorded_sha = record["content_sha256"]
    except (ValueError, KeyError, TypeError) as error:
        raise StrainRollupRecordError(
            f"{path} is not a valid strain rollup record: {error}"
        ) from error
    if recorded_sha != content_sha256(table):
        raise StrainRollupRecordError(
            f"{table} is not the table its record {path.name} describes, so its strain "
            "rollup scope is unknown; rebuild it, or restore the committed profile table "
            "with `just study-trait-profiles-restore` (#847)"
        )
    return scope


def require_strain_rollup(table: Path) -> StrainRollup:
    """The scope a table records; refused when it has no record."""
    scope = recorded_strain_rollup(table)
    if scope is None:
        raise StrainRollupRecordError(
            f"{table} records no strain rollup scope ({record_path(table).name} is missing); "
            "rebuild it with the recipe that writes it, such as `just taxa-trait-matrix` "
            "(#847)"
        )
    return scope


def shared_strain_rollup(*tables: Path) -> StrainRollup | None:
    """The scope the recorded tables share, or None when none has a record.

    Tables without a record are not compared; two recorded scopes that differ are refused.
    """
    recorded = {table: recorded_strain_rollup(table) for table in tables}
    scopes = {scope for scope in recorded.values() if scope is not None}
    if len(scopes) > 1:
        listed = ", ".join(f"{table} {scope}" for table, scope in recorded.items() if scope)
        raise StrainRollupRecordError(
            f"These tables were built under different strain rollup scopes: {listed} (#847)"
        )
    return next(iter(scopes), None)


def require_declared_strain_rollup(tables: Sequence[Path], declared: StrainRollup) -> None:
    """Refuse recorded tables built under another scope than the one a step declares."""
    recorded = shared_strain_rollup(*tables)
    if recorded is not None and recorded != declared:
        names = ", ".join(str(table) for table in tables)
        raise StrainRollupRecordError(
            f"{names} were built under the {recorded} strain rollup scope, not the "
            f"declared {declared} (#847)"
        )


def committed_path(table: Path) -> Path:
    """The committed gzip copy of a canonical table (#656, #667)."""
    return table.with_name(table.name + ".gz")


def _restore_blocked(scope: StrainRollup, matrix: Path) -> str | None:
    """Why restoring a committed profile table of ``scope`` would mix scopes, if it would."""
    if not matrix.exists():
        return None
    try:
        matrix_scope = recorded_strain_rollup(matrix)
    except StrainRollupRecordError:
        return f"the matrix {matrix} is not the table its strain rollup record describes"
    if matrix_scope is None:
        return f"the matrix {matrix} records no strain rollup scope"
    if matrix_scope != scope:
        return (
            f"restoring it would put the {scope} profile table next to the {matrix_scope} "
            f"matrix {matrix}"
        )
    return None


def check_committed(table: Path, matrix: Path) -> list[str]:
    """Why ``table`` is not the committed canonical table; empty when it is (#656, #847)."""
    committed = committed_path(table)
    committed_record = record_path(committed)
    try:
        committed_sha = content_sha256(committed)
    except (OSError, EOFError, zlib.error) as error:
        return [
            f"{committed} is missing, empty or not gzip ({error}); restore it: "
            f"git checkout -- {committed} {committed_record} (#667)"
        ]
    try:
        committed_scope = require_strain_rollup(committed)
    except StrainRollupRecordError as error:
        return [f"{error}", f"  Restore it: git checkout -- {committed} {committed_record}"]
    blocked = _restore_blocked(committed_scope, matrix)
    rebuild = "run just study-trait-profiles" + (
        "" if blocked else ", or just study-trait-profiles-restore"
    )
    unrestorable = [f"  Not offering the committed table: {blocked}."] if blocked else []
    if not table.is_file() or table.stat().st_size == 0:
        return [f"{table} is missing or empty; {rebuild}", *unrestorable]
    try:
        scope = require_strain_rollup(table)
    except StrainRollupRecordError as error:
        return [f"{error}", f"  Unintended: {rebuild}", *unrestorable]
    same_content = content_sha256(table) == committed_sha
    if same_content and scope == committed_scope:
        return []
    if same_content:
        return [
            f"{table} has the committed content but records the {scope} strain rollup "
            f"scope, the committed record {committed_scope}; rebuild it (#847)"
        ]
    messages = [f"{table} differs from the committed {committed} (#656)."]
    if scope != committed_scope:
        messages.append(
            f"  It was built under the {scope} strain rollup scope, the committed table "
            f"under {committed_scope} (#788, #847)."
        )
    return [
        *messages,
        f"  Unintended: {rebuild}",
        *unrestorable,
        "  After an intended input change, including a new strain rollup scope: "
        f"just study-trait-profiles-refresh-gz, then commit {committed} and "
        f"{committed_record}",
    ]


def restore_committed(table: Path, matrix: Path) -> StrainRollup:
    """Decompress the committed table over ``table`` and give it the committed record.

    Refused when the local matrix records another scope or none, so a restore never puts
    a profile table next to a matrix of another scope (#847).
    """
    committed = committed_path(table)
    scope = require_strain_rollup(committed)
    blocked = _restore_blocked(scope, matrix)
    if blocked:
        raise StrainRollupRecordError(
            f"Not restoring {committed}: {blocked}; run just study-trait-profiles (#847)"
        )
    partial = table.with_name(f".{table.name}.restoring")
    with gzip.open(committed, "rb") as source, partial.open("wb") as target:
        shutil.copyfileobj(source, target)
    os.replace(partial, table)
    write_strain_rollup_record(table, scope)
    return scope


def record_committed(table: Path) -> Path:
    """Give a freshly compressed committed table the local table's record (#667, #847)."""
    committed = committed_path(table)
    scope = require_strain_rollup(table)
    if content_sha256(committed) != content_sha256(table):
        raise StrainRollupRecordError(f"{committed} does not decompress to {table}")
    return write_strain_rollup_record(committed, scope)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    show = commands.add_parser("show", help="Print the scope a table records; fail without one.")
    show.add_argument("table", type=Path)
    for name, summary in (
        ("check", "Check the profile table against its committed copy and record."),
        ("restore", "Restore the committed profile table and its record."),
        ("record-committed", "Record the scope of a freshly compressed committed table."),
    ):
        command = commands.add_parser(name, help=summary)
        command.add_argument("table", type=Path, nargs="?", default=DEFAULT_PROFILE_PATH)
        command.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX_PATH)
    args = parser.parse_args(argv)
    try:
        if args.command == "show":
            print(require_strain_rollup(args.table))
        elif args.command == "check":
            problems = check_committed(args.table, args.matrix)
            for line in problems:
                print(line, file=sys.stderr)
            return 1 if problems else 0
        elif args.command == "restore":
            scope = restore_committed(args.table, args.matrix)
            print(f"Restored {args.table} ({scope}) from {committed_path(args.table)}")
        else:
            print(f"Wrote {record_committed(args.table)}")
    except (StrainRollupRecordError, OSError, EOFError, zlib.error) as error:
        print(error, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
