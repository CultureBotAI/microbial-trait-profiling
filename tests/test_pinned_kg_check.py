"""`just verify-pinned-kg` checks data/kg-microbe against the protocol's md5s (#695)."""

import hashlib
import json
from pathlib import Path

import pytest

from mtp import contribution_inputs

EDGES = b"id\tsubject\tpredicate\tobject\n"
NODES = b"id\tname\n"


def md5(data: bytes) -> str:
    return hashlib.md5(data, usedforsecurity=False).hexdigest()


def pinned_kg(tmp_path: Path, edges: bytes = EDGES, nodes: bytes | None = NODES) -> list[str]:
    """A KG folder and a protocol pinning EDGES and NODES; return the CLI arguments."""
    kg_dir = tmp_path / "kg-microbe"
    kg_dir.mkdir()
    (kg_dir / "merged-kg_edges.tsv").write_bytes(edges)
    if nodes is not None:
        (kg_dir / "merged-kg_nodes.tsv").write_bytes(nodes)
    protocol = tmp_path / "protocol.json"
    pinned = {"merged-kg_edges.tsv": md5(EDGES), "merged-kg_nodes.tsv": md5(NODES)}
    protocol.write_text(json.dumps({"kg_md5": pinned}), encoding="utf-8")
    return ["--kg-dir", str(kg_dir), "--protocol", str(protocol)]


def test_pinned_kg_passes(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    contribution_inputs.main(pinned_kg(tmp_path))

    assert "match the md5s in protocol.json" in capsys.readouterr().out


def test_changed_kg_file_is_refused_and_names_the_forced_refetch(tmp_path: Path) -> None:
    """`just download-kg` compares sizes only, so the remedy is the forced download (#703)."""
    with pytest.raises(
        SystemExit,
        match=r"Pinned KG mismatch .*merged-kg_edges\.tsv.*; refetch the KG with "
        r"`just download-kg-force`$",
    ):
        contribution_inputs.main(pinned_kg(tmp_path, edges=b"another release\n"))


def test_missing_kg_file_is_refused(tmp_path: Path) -> None:
    with pytest.raises(SystemExit, match=r"lacks the pinned KG files merged-kg_nodes\.tsv"):
        contribution_inputs.main(pinned_kg(tmp_path, nodes=None))
