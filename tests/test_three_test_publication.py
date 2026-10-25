"""Publication permissions without changing the caller's process mask (#1081)."""

from __future__ import annotations

import importlib
import multiprocessing
import os
from multiprocessing.connection import Connection
from unittest.mock import Mock

import pytest


def _observe_requested_mask(mask: int, connection: Connection) -> None:
    """Only this isolated test child sets a known mask; its parent stays untouched."""
    publication = importlib.import_module("mtp.three_test_table.publication")
    os.umask(mask)
    mode = publication.directory_mode()
    unchanged = os.umask(mask)
    connection.send_bytes(f"{mode},{unchanged}".encode("ascii"))
    connection.close()


def _silent_child(connection: Connection) -> None:
    connection.close()


def _stalled_child(connection: Connection) -> None:
    del connection
    multiprocessing.Event().wait(60)


def _invalid_child(connection: Connection) -> None:
    connection.send_bytes(b"512")
    connection.close()


def _oversized_child(connection: Connection) -> None:
    connection.send_bytes(b"1" * 32)
    connection.close()


def _failed_child(connection: Connection) -> None:
    connection.send_bytes(b"18")
    connection.close()
    raise SystemExit(1)


@pytest.mark.parametrize("mask", [0o002, 0o022, 0o027, 0o077])
def test_directory_mode_reads_inherited_mask_without_changing_it(mask: int) -> None:
    importlib.import_module("mtp.three_test_table.publication")
    context = multiprocessing.get_context("spawn")
    receive, send = context.Pipe(duplex=False)
    child = context.Process(target=_observe_requested_mask, args=(mask, send))
    child.start()
    send.close()
    try:
        assert receive.poll(20), "permission observer child did not respond"
        result = receive.recv_bytes(64).decode("ascii")
        child.join(5)
        assert child.exitcode == 0
        assert result == f"{0o777 & ~mask},{mask}"
    finally:
        receive.close()
        if child.is_alive():
            child.kill()
            child.join(5)
        child.close()


def test_directory_mode_never_calls_parent_umask(monkeypatch: pytest.MonkeyPatch) -> None:
    publication = importlib.import_module("mtp.three_test_table.publication")
    forbidden = Mock(side_effect=AssertionError("parent process umask must remain untouched"))
    monkeypatch.setattr(os, "umask", forbidden)
    result = publication.directory_mode()
    assert 0 <= result <= 0o777
    forbidden.assert_not_called()


@pytest.mark.parametrize("failure", ["closed_pipe", "timeout", "invalid", "oversized", "failed"])
def test_mask_query_failure_is_bounded_and_leaves_no_child(
    monkeypatch: pytest.MonkeyPatch, failure: str
) -> None:
    publication = importlib.import_module("mtp.three_test_table.publication")
    worker = {
        "closed_pipe": _silent_child,
        "timeout": _stalled_child,
        "invalid": _invalid_child,
        "oversized": _oversized_child,
        "failed": _failed_child,
    }[failure]
    monkeypatch.setattr(publication, "_read_umask", worker)
    before = {process.pid for process in multiprocessing.active_children()}
    with pytest.raises(RuntimeError, match="umask"):
        publication.directory_mode(timeout=0.5 if failure == "timeout" else 5)
    assert {process.pid for process in multiprocessing.active_children()} == before
