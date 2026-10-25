"""Read inherited publication permissions without changing the caller's umask (#1081)."""

from __future__ import annotations

import multiprocessing
import os
from multiprocessing.connection import Connection
from multiprocessing.process import BaseProcess


def _read_umask(connection: Connection) -> None:
    """Only the fresh child changes its own mask; the parent mask is never touched."""
    try:
        connection.send_bytes(str(os.umask(0o777)).encode("ascii"))
    finally:
        connection.close()


def _close_child(child: BaseProcess) -> None:
    """Reap this owned child on every exit, escalating only if it does not stop."""
    if child.pid is not None:
        if child.is_alive():
            child.terminate()
        child.join(2)
        if child.is_alive():
            child.kill()
            child.join(2)
        if child.is_alive():
            raise RuntimeError("Could not stop the umask worker")
    child.close()


def directory_mode(*, timeout: float = 10.0) -> int:
    """Return 0777 masked by the inherited umask, with a bounded isolated query.

    Python exposes no portable read-only umask getter. A set/restore pair in
    the parent races with concurrent file creation. Spawn instead of forking
    the already imported numerical runtime, and exchange only bounded bytes.
    """
    if timeout <= 0:
        raise ValueError("The umask query timeout must be positive")
    context = multiprocessing.get_context("spawn")
    receive, send = context.Pipe(duplex=False)
    child = context.Process(target=_read_umask, args=(send,))
    try:
        child.start()
        send.close()
        if not receive.poll(timeout):
            raise RuntimeError("The umask worker timed out")
        try:
            encoded = receive.recv_bytes(16)
        except (EOFError, OSError) as error:
            raise RuntimeError("The umask worker returned no valid response") from error
        child.join(timeout)
        if child.exitcode != 0:
            raise RuntimeError("The umask worker did not exit successfully")
        if not encoded.isdigit() or not 0 <= int(encoded) <= 0o777:
            raise RuntimeError("The umask worker returned an invalid mask")
        return 0o777 & ~int(encoded)
    finally:
        send.close()
        receive.close()
        _close_child(child)
