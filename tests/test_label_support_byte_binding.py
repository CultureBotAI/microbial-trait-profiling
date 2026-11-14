"""A label-support parse cannot consume bytes absent from its input receipt (#1357)."""

from __future__ import annotations

import builtins
import hashlib
import io
import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pytest

from mtp import label_support_audit
from test_label_support_audit import inputs as inputs
from test_label_support_audit import run_audit, write_facultative_metadata

ROLES = ("matrix", "profiles", "taxa", "edges", "nodes", "madin", "columns")


@dataclass(frozen=True)
class ReadTrace:
    inventory_digests: list[str]
    parser_stream_modes: list[str]


def _contrasting_bytes(
    role: str, inputs: dict[str, Path], columns: Path
) -> tuple[Path, bytes, bytes]:
    """Keep the fixture reconcilable so only the consumed-byte binding can reject B."""
    target = columns if role == "columns" else inputs[role]
    original = target.read_bytes()
    if role == "profiles":
        # A and B both parse to integer 1, with different same-size byte spellings.
        assert b"\t1\t" in original
        original = original.replace(b"\t1\t", b"\t01\t", 1)
        target.write_bytes(original)
        replacement = original.replace(b"\t01\t", b"\t1 \t", 1)
    else:
        before, after = {
            "matrix": (b"Examplegenus 10", b"Examplegenus 90"),
            "taxa": (b"increased", b"INCREASED"),
            "edges": (b"infores:madin_etal", b"infores:other_etal"),
            "nodes": (b"A alpha", b"A omega"),
            "madin": (b"gold", b"GOLD"),
            "columns": (b"METPO:1000605", b"METPO:1000606"),
        }[role]
        assert before in original
        replacement = original.replace(before, after, 1)
        if role == "columns":
            assert b"METPO:1000612" in replacement
            replacement = replacement.replace(b"METPO:1000612", b"METPO:1000611", 1)
    assert replacement != original
    assert len(replacement) == len(original) < 65536
    return target, original, replacement


def _swap_parser_streams(
    monkeypatch: pytest.MonkeyPatch, target: Path, original: bytes, replacement: bytes
) -> ReadTrace:
    """Both inventories see real A; only the intervening parser opens receive B.

    This crosses the public file-opening boundary, so it also exercises a new binary
    stream parser without depending on the name or internals of its future helper.
    The real file remains A throughout; no restore timing or filesystem race is needed.
    """
    inventory = label_support_audit._inventory
    path_open = Path.open
    builtin_open = builtins.open
    selected = target.absolute()
    expected_digest = hashlib.sha256(original).hexdigest()
    inventory_depth = 0
    armed = False
    trace = ReadTrace([], [])

    def is_target(file: object) -> bool:
        return (
            isinstance(file, (str, bytes, os.PathLike))
            and Path(os.fsdecode(file)).absolute() == selected
        )

    def guarded_inventory(path: Path, *, allow_absent: bool = False) -> dict[str, Any]:
        nonlocal inventory_depth, armed
        inventory_depth += 1
        try:
            record = inventory(path, allow_absent=allow_absent)
        finally:
            inventory_depth -= 1
        if is_target(path):
            assert record["sha256"] == expected_digest
            assert record["size_bytes"] == len(original)
            trace.inventory_digests.append(record["sha256"])
            armed = True
        return record

    def parser_stream(file: object, args: tuple[Any, ...], kwargs: dict[str, Any]) -> Any:
        if not armed or inventory_depth or not is_target(file):
            return None
        mode = args[0] if args else kwargs.get("mode", "r")
        assert isinstance(mode, str)
        assert not any(flag in mode for flag in "wax+")
        trace.parser_stream_modes.append(mode)
        if "b" in mode:
            return io.BytesIO(replacement)
        encoding = kwargs.get("encoding", args[2] if len(args) > 2 else "utf-8") or "utf-8"
        errors = kwargs.get("errors", args[3] if len(args) > 3 else "strict") or "strict"
        newline = kwargs.get("newline", args[4] if len(args) > 4 else None)
        return io.StringIO(replacement.decode(encoding, errors), newline=newline)

    def opened_path(path: Path, *args: Any, **kwargs: Any) -> Any:
        altered = parser_stream(path, args, kwargs)
        return altered if altered is not None else path_open(path, *args, **kwargs)

    def opened_builtin(file: Any, *args: Any, **kwargs: Any) -> Any:
        altered = parser_stream(file, args, kwargs)
        return altered if altered is not None else builtin_open(file, *args, **kwargs)

    monkeypatch.setattr(label_support_audit, "_inventory", guarded_inventory)
    monkeypatch.setattr(Path, "open", opened_path)
    monkeypatch.setattr(builtins, "open", opened_builtin)
    return trace


@pytest.mark.parametrize("role", ROLES)
def test_audit_refuses_unrecorded_parser_bytes_after_inventory(
    inputs: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch, role: str
) -> None:
    columns = tmp_path / "selected_columns.tsv"
    write_facultative_metadata(columns, union=True)
    target, original, replacement = _contrasting_bytes(role, inputs, columns)
    output = tmp_path / "unbound-byte-audit"
    with monkeypatch.context() as patch:
        trace = _swap_parser_streams(patch, target, original, replacement)
        with pytest.raises(ValueError, match="Audit input bytes differ"):
            run_audit(inputs, output, columns=columns)
    assert trace.inventory_digests
    assert trace.parser_stream_modes, f"No parser read reached the {role} swap boundary"
    assert target.read_bytes() == original
    assert not output.exists()


@pytest.mark.parametrize("role", ROLES)
def test_audit_consumes_one_bound_binary_parser_stream_per_input(
    inputs: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch, role: str
) -> None:
    """Opening a checked stream is insufficient if a reader then reopens its path."""
    columns = tmp_path / "selected_columns.tsv"
    write_facultative_metadata(columns, union=True)
    target, original, _ = _contrasting_bytes(role, inputs, columns)
    output = tmp_path / "bound-byte-audit"
    with monkeypatch.context() as patch:
        trace = _swap_parser_streams(patch, target, original, original)
        run_audit(inputs, output, columns=columns)
    # Preliminary/terminal inventory and the early Madin pin check are excluded
    # by the guard. Both the header and rows must use this one binary source.
    assert trace.parser_stream_modes == ["rb"]
    assert len(trace.inventory_digests) == 2
    provenance = json.loads((output / "provenance.json").read_text(encoding="utf-8"))
    key = {"edges": "kg_edges", "nodes": "kg_nodes", "columns": "column_metadata"}.get(role, role)
    assert provenance["inputs"][key]["sha256"] == hashlib.sha256(original).hexdigest()
    assert provenance["inputs"][key]["size_bytes"] == len(original)
    assert len(provenance["outputs"]) == 5
    assert all((output / name).is_file() for name in provenance["outputs"])
    assert target.read_bytes() == original
