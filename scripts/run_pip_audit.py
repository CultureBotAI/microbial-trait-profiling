"""Run pip-audit with bounded retries on transient network errors (#1113).

Transient vulnerability-feed disconnects (such as PyPI or OSV closing a TLS
connection mid-read, connection timeouts, or HTTP 408/429/5xx responses) are
retried with exponential backoff up to a bounded attempt count.

Vulnerability findings and non-network failures (such as malformed versions,
HTTP 4xx client errors, or invalid CLI/source arguments) fail immediately on
the first attempt without retry, and persistent network failures exit with code
2 and an explicit diagnostic distinguishing them from vulnerability findings.
"""

from __future__ import annotations

import contextlib
import http.client
import importlib
import io
import logging
import sys
import time
from collections.abc import Callable, Iterator, Sequence
from typing import Any, TextIO

DEFAULT_MAX_ATTEMPTS = 4
DEFAULT_INITIAL_BACKOFF_SECONDS = 2.0
DEFAULT_BACKOFF_FACTOR = 2.0
NETWORK_FAILURE_EXIT_CODE = 2
TRANSIENT_HTTP_STATUS_CODES = frozenset({408, 429})


class _UnmatchedVulnServiceConnectionError(Exception):
    """Sentinel class replacing pip_audit._cli.VulnServiceConnectionError.

    pip_audit._cli.audit() catches VulnServiceConnectionError and calls
    sys.exit(1). Temporarily rebinding that name inside pip_audit._cli lets
    pip_audit._service.interface.ConnectionError propagate as an exception to
    our retry loop instead of terminating via SystemExit(1).
    """


def _iter_exception_chain(exc: BaseException) -> Iterator[BaseException]:
    """Yield ``exc`` and its linked ``__cause__`` / ``__context__`` chain."""
    seen: set[int] = set()
    current: BaseException | None = exc
    while current is not None and id(current) not in seen:
        seen.add(id(current))
        yield current
        current = current.__cause__ if current.__cause__ is not None else current.__context__


def _is_transient_http_status(status_code: object) -> bool:
    if not isinstance(status_code, int):
        return False
    return status_code in TRANSIENT_HTTP_STATUS_CODES or 500 <= status_code < 600


def is_transient_network_error(exc: BaseException) -> bool:
    """Return True when ``exc`` was caused by a transient network/HTTP failure."""
    requests_exceptions: Any = importlib.import_module("requests.exceptions")
    urllib3_exceptions: Any = importlib.import_module("urllib3.exceptions")
    pip_audit_interface: Any = importlib.import_module("pip_audit._service.interface")

    vuln_conn_error_cls = pip_audit_interface.ConnectionError
    http_error_cls = requests_exceptions.HTTPError
    transient_request_classes = (
        vuln_conn_error_cls,
        requests_exceptions.ConnectionError,
        requests_exceptions.Timeout,
        requests_exceptions.ChunkedEncodingError,
        requests_exceptions.TooManyRedirects,
        requests_exceptions.RetryError,
        urllib3_exceptions.HTTPError,
        http.client.HTTPException,
        ConnectionError,
        TimeoutError,
    )

    for item in _iter_exception_chain(exc):
        if isinstance(item, http_error_cls):
            response = getattr(item, "response", None)
            if response is None:
                return True
            status_code = getattr(response, "status_code", None)
            return _is_transient_http_status(status_code)
        if isinstance(item, transient_request_classes):
            return True
    return False


def _format_exception_detail(exc: BaseException) -> str:
    parts: list[str] = []
    for item in _iter_exception_chain(exc):
        text = str(item).strip()
        label = f"{type(item).__name__}: {text}" if text else type(item).__name__
        if label not in parts:
            parts.append(label)
    return " <- ".join(parts)


def _is_vulnerability_finding_output(stderr_text: str) -> bool:
    return "Found " in stderr_text and " known vulnerabilit" in stderr_text


def _normalize_exit_code(code: object) -> int:
    if code is None:
        return 0
    if isinstance(code, int):
        return code
    return 1


def _invoke_pip_audit_cli(argv: Sequence[str]) -> None:
    cli_mod: Any = importlib.import_module("pip_audit._cli")
    pip_audit_interface: Any = importlib.import_module("pip_audit._service.interface")
    saved_argv = sys.argv[:]
    sys.argv = ["pip-audit", *argv]
    cli_mod.VulnServiceConnectionError = _UnmatchedVulnServiceConnectionError
    try:
        raise SystemExit(cli_mod.audit())
    finally:
        cli_mod.VulnServiceConnectionError = pip_audit_interface.ConnectionError
        sys.argv = saved_argv


@contextlib.contextmanager
def _capture_attempt_streams() -> Iterator[tuple[io.StringIO, io.StringIO]]:
    buf_out = io.StringIO()
    buf_err = io.StringIO()
    pip_audit_logger = logging.getLogger("pip_audit")
    temp_handler = logging.StreamHandler(buf_err)
    saved_handlers = list(pip_audit_logger.handlers)
    saved_propagate = pip_audit_logger.propagate
    pip_audit_logger.handlers = [temp_handler]
    pip_audit_logger.propagate = False
    try:
        with contextlib.redirect_stdout(buf_out), contextlib.redirect_stderr(buf_err):
            yield buf_out, buf_err
    finally:
        temp_handler.flush()
        pip_audit_logger.handlers = saved_handlers
        pip_audit_logger.propagate = saved_propagate


def run_pip_audit(
    argv: Sequence[str],
    *,
    max_attempts: int = DEFAULT_MAX_ATTEMPTS,
    initial_backoff: float = DEFAULT_INITIAL_BACKOFF_SECONDS,
    backoff_factor: float = DEFAULT_BACKOFF_FACTOR,
    audit_fn: Callable[[Sequence[str]], int | str | None] | None = None,
    sleep_fn: Callable[[float], None] = time.sleep,
    stdout: TextIO | None = None,
    stderr: TextIO | None = None,
) -> int:
    """Execute ``pip-audit`` with bounded retries on transient network errors."""
    if max_attempts < 1:
        raise ValueError(f"max_attempts must be >= 1, got {max_attempts}")
    if initial_backoff < 0.0:
        raise ValueError(f"initial_backoff must be >= 0.0, got {initial_backoff}")
    if backoff_factor < 1.0:
        raise ValueError(f"backoff_factor must be >= 1.0, got {backoff_factor}")

    runner = _invoke_pip_audit_cli if audit_fn is None else audit_fn
    out_stream = sys.stdout if stdout is None else stdout
    err_stream = sys.stderr if stderr is None else stderr

    for attempt in range(1, max_attempts + 1):
        try:
            with _capture_attempt_streams() as (attempt_out, attempt_err):
                result = runner(argv)
                if result is not None:
                    raise SystemExit(result)
        except SystemExit as exit_exc:
            code = _normalize_exit_code(exit_exc.code)
            out_text = attempt_out.getvalue()
            err_text = attempt_err.getvalue()
            if out_text:
                out_stream.write(out_text)
            if err_text:
                err_stream.write(err_text)
            if code == 0:
                return 0
            if _is_vulnerability_finding_output(err_text):
                err_stream.write(
                    "[pip-audit] VULNERABILITY FINDING: pip-audit reported known "
                    "vulnerabilities; failing immediately without retry.\n"
                )
                return code
            return code
        except Exception as exc:
            if not is_transient_network_error(exc):
                out_text = attempt_out.getvalue()
                err_text = attempt_err.getvalue()
                if out_text:
                    out_stream.write(out_text)
                if err_text:
                    err_stream.write(err_text)
                raise

            detail = _format_exception_detail(exc)
            if attempt < max_attempts:
                delay = initial_backoff * (backoff_factor ** (attempt - 1))
                err_stream.write(
                    f"[pip-audit] NETWORK ERROR (attempt {attempt}/{max_attempts}): "
                    f"{detail}. Retrying in {delay:.1f}s...\n"
                )
                sleep_fn(delay)
                continue

            err_stream.write(
                f"[pip-audit] NETWORK FAILURE: pip-audit failed after {max_attempts} "
                f"attempts due to transient network error ({detail}); "
                "not a vulnerability finding.\n"
            )
            return NETWORK_FAILURE_EXIT_CODE
        else:
            out_text = attempt_out.getvalue()
            err_text = attempt_err.getvalue()
            if out_text:
                out_stream.write(out_text)
            if err_text:
                err_stream.write(err_text)
            return 0

    return NETWORK_FAILURE_EXIT_CODE


def main(argv: Sequence[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    return run_pip_audit(args)


if __name__ == "__main__":
    raise SystemExit(main())
