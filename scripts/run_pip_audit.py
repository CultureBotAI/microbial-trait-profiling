"""Run pip-audit with typed, bounded transient-network retries (#1113, #1185).

Only a completed successful audit exits zero. Findings and usage errors retain
pip-audit's nonzero codes; recognized network failures use 3, distinct from
argparse's 2. Certificate/hostname verification, proxy, redirect and unclassified
HTTP failures are not retried, even beneath a broader connection-error wrapper.
Typed TLS EOF and zero-return closures remain retryable.

The serial in-process adapter retains exception types instead of classifying
stderr text. It checks the installed private API before capture and restores
its actual prior bindings. Console output streams live, preserving retry
attempts and logging prefixes. See docs/pip_audit_retry.md for the subprocess
tradeoff and the process-global, single-threaded scope of this adapter.
"""

from __future__ import annotations

import contextlib
import http.client
import importlib
import importlib.util
import io
import logging
import ssl
import sys
import time
from collections.abc import Callable, Iterator, Sequence
from functools import partial
from importlib.metadata import PackageNotFoundError, version
from typing import Any, TextIO

DEFAULT_MAX_ATTEMPTS = 4
DEFAULT_INITIAL_BACKOFF_SECONDS = 2.0
DEFAULT_BACKOFF_FACTOR = 2.0
NETWORK_FAILURE_EXIT_CODE = 3
TRANSIENT_HTTP_STATUS_CODES = frozenset({408, 429})


class PipAuditCompatibilityError(RuntimeError):
    """The installed pip-audit no longer exposes the adapter's required API."""


class PipAuditNotInstalledError(RuntimeError):
    """The current environment lacks pip-audit, rather than an adapter API symbol."""


class _UnmatchedVulnServiceConnectionError(Exception):
    """Leave the real service error uncaught by pip-audit's CLI exit handler."""


def _compatibility_error(detail: str) -> PipAuditCompatibilityError:
    try:
        installed = version("pip-audit")
    except PackageNotFoundError:
        installed = "(distribution metadata unavailable)"
    return PipAuditCompatibilityError(f"pip-audit {installed}: incompatible adapter API: {detail}")


def _import_api(name: str) -> Any:
    try:
        return importlib.import_module(name)
    except ImportError as error:
        raise _compatibility_error(f"cannot import {name}: {error}") from error


def _exception_class(module: Any, name: str) -> type[Exception]:
    value = getattr(module, name, None)
    if not isinstance(value, type) or not issubclass(value, Exception):
        raise _compatibility_error(f"{module.__name__}.{name} is missing or not an exception class")
    return value


def _prepare_pip_audit() -> tuple[Any, type[Exception]]:
    if importlib.util.find_spec("pip_audit") is None:
        raise PipAuditNotInstalledError(
            "pip-audit is not installed in this environment (run just sync)"
        )
    # Import outside capture: basicConfig must not retain attempt 1's buffer.
    cli_mod = _import_api("pip_audit._cli")
    service_mod = _import_api("pip_audit._service.interface")
    connection_error = _exception_class(service_mod, "ConnectionError")
    _exception_class(cli_mod, "VulnServiceConnectionError")
    if not callable(getattr(cli_mod, "audit", None)):
        raise _compatibility_error("pip_audit._cli.audit is missing or not callable")
    return cli_mod, connection_error


def _iter_exception_chain(exc: BaseException) -> Iterator[BaseException]:
    """Visit typed causes, contexts, reasons and exception-valued args once each."""
    seen: set[int] = set()
    pending = [exc]
    while pending:
        current = pending.pop()
        if id(current) in seen:
            continue
        seen.add(id(current))
        yield current
        cause = current.__cause__ if current.__cause__ is not None else current.__context__
        for linked in (getattr(current, "reason", None), cause, *current.args):
            if isinstance(linked, BaseException):
                pending.append(linked)


def _is_transient_http_status(status_code: object) -> bool:
    if not isinstance(status_code, int) or isinstance(status_code, bool):
        return False
    return status_code in TRANSIENT_HTTP_STATUS_CODES or 500 <= status_code < 600


def is_transient_network_error(
    exc: BaseException, *, service_connection_error: type[Exception] | None = None
) -> bool:
    """Retry known transient causes only; permanent causes override broad wrappers."""
    requests_exceptions: Any = importlib.import_module("requests.exceptions")
    urllib3_exceptions: Any = importlib.import_module("urllib3.exceptions")
    if service_connection_error is None:
        service_connection_error = _exception_class(
            _import_api("pip_audit._service.interface"), "ConnectionError"
        )
    chain = list(_iter_exception_chain(exc))
    permanent_classes = (
        ssl.SSLCertVerificationError,
        ssl.CertificateError,
        requests_exceptions.ProxyError,
        requests_exceptions.TooManyRedirects,
        urllib3_exceptions.ProxyError,
    )
    if any(isinstance(item, permanent_classes) for item in chain):
        return False
    http_errors = [item for item in chain if isinstance(item, requests_exceptions.HTTPError)]
    if any(
        not _is_transient_http_status(getattr(getattr(item, "response", None), "status_code", None))
        for item in http_errors
    ):
        return False
    transient_classes = (
        service_connection_error,
        requests_exceptions.ConnectionError,
        requests_exceptions.Timeout,
        requests_exceptions.ChunkedEncodingError,
        urllib3_exceptions.ProtocolError,
        urllib3_exceptions.TimeoutError,
        ssl.SSLEOFError,
        ssl.SSLZeroReturnError,
        http.client.BadStatusLine,
        http.client.IncompleteRead,
        ConnectionError,
        TimeoutError,
    )
    return bool(http_errors) or any(isinstance(item, transient_classes) for item in chain)


def _is_recognized_network_error(exc: BaseException) -> bool:
    """Label typed network refusals without treating arbitrary service errors as network."""
    requests_exceptions: Any = importlib.import_module("requests.exceptions")
    urllib3_exceptions: Any = importlib.import_module("urllib3.exceptions")
    network_classes = (
        requests_exceptions.ConnectionError,
        requests_exceptions.Timeout,
        requests_exceptions.ChunkedEncodingError,
        requests_exceptions.HTTPError,
        requests_exceptions.TooManyRedirects,
        urllib3_exceptions.HTTPError,
        ssl.SSLError,
        ssl.CertificateError,
        http.client.BadStatusLine,
        http.client.IncompleteRead,
        ConnectionError,
        TimeoutError,
    )
    return any(isinstance(item, network_classes) for item in _iter_exception_chain(exc))


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


def _invoke_pip_audit_cli(argv: Sequence[str], *, cli_mod: Any = None) -> None:
    if cli_mod is None:
        cli_mod, _ = _prepare_pip_audit()
    saved_error = _exception_class(cli_mod, "VulnServiceConnectionError")
    saved_argv = sys.argv
    try:
        sys.argv = ["pip-audit", *argv]
        cli_mod.VulnServiceConnectionError = _UnmatchedVulnServiceConnectionError
        raise SystemExit(cli_mod.audit())
    finally:
        cli_mod.VulnServiceConnectionError = saved_error
        sys.argv = saved_argv


class _AttemptStream(io.TextIOBase):
    """Forward output immediately, optionally retaining stderr for finding labels."""

    def __init__(self, target: TextIO, capture: io.StringIO | None = None) -> None:
        self.target = target
        self.capture = capture

    def write(self, text: str) -> int:
        if self.capture is not None:
            self.capture.write(text)
        written = self.target.write(text)
        self.target.flush()
        return written

    def flush(self) -> None:
        self.target.flush()

    def isatty(self) -> bool:
        return self.target.isatty()

    def fileno(self) -> int:
        return self.target.fileno()


@contextlib.contextmanager
def _capture_attempt_streams(stdout: TextIO, stderr: TextIO) -> Iterator[io.StringIO]:
    captured_error = io.StringIO()
    live_out = _AttemptStream(stdout)
    live_err = _AttemptStream(stderr, captured_error)
    original_out, original_err = sys.stdout, sys.stderr
    root = logging.getLogger()
    package = logging.getLogger("pip_audit")
    restored_streams: list[tuple[Any, Any]] = []
    added_handlers: list[tuple[logging.Logger, logging.Handler]] = []
    console_loggers: set[logging.Logger] = set()
    seen_handlers: set[int] = set()
    try:
        for logger in (root, package):
            for handler in logger.handlers:
                if not isinstance(handler, logging.StreamHandler):
                    continue
                if id(handler) in seen_handlers:
                    console_loggers.add(logger)
                    continue
                if handler.stream is original_err or handler.stream is stderr:
                    destination = live_err
                elif handler.stream is original_out or handler.stream is stdout:
                    destination = live_out
                else:
                    continue
                console_loggers.add(logger)
                if id(handler) not in seen_handlers:
                    previous = handler.setStream(destination)
                    restored_streams.append((handler, previous))
                    seen_handlers.add(id(handler))
        if root not in console_loggers:
            handler = logging.StreamHandler(live_err)
            handler.setFormatter(logging.Formatter(logging.BASIC_FORMAT))
            if package in console_loggers:
                # Keep capture/file handlers, without duplicating a package console line.
                handler.addFilter(
                    lambda record: (
                        not (record.name == "pip_audit" or record.name.startswith("pip_audit."))
                    )
                )
            root.addHandler(handler)
            added_handlers.append((root, handler))
        if not package.propagate and package not in console_loggers:
            handler = logging.StreamHandler(live_err)
            handler.setFormatter(logging.Formatter(logging.BASIC_FORMAT))
            package.addHandler(handler)
            added_handlers.append((package, handler))
        with contextlib.redirect_stdout(live_out), contextlib.redirect_stderr(live_err):
            yield captured_error
    finally:
        for logger, handler in reversed(added_handlers):
            logger.removeHandler(handler)
            handler.flush()
            handler.close()
        for handler, previous in reversed(restored_streams):
            handler.setStream(previous)


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
    """Execute serial attempts without converting an incomplete audit into success."""
    if max_attempts < 1:
        raise ValueError(f"max_attempts must be >= 1, got {max_attempts}")
    if initial_backoff < 0.0:
        raise ValueError(f"initial_backoff must be >= 0.0, got {initial_backoff}")
    if backoff_factor < 1.0:
        raise ValueError(f"backoff_factor must be >= 1.0, got {backoff_factor}")

    service_connection_error = None
    runner: Callable[[Sequence[str]], int | str | None]
    if audit_fn is None:
        cli_mod, service_connection_error = _prepare_pip_audit()
        runner = partial(_invoke_pip_audit_cli, cli_mod=cli_mod)
    else:
        runner = audit_fn
    out_stream = sys.stdout if stdout is None else stdout
    err_stream = sys.stderr if stderr is None else stderr

    for attempt in range(1, max_attempts + 1):
        try:
            with _capture_attempt_streams(out_stream, err_stream) as attempt_err:
                result = runner(argv)
                if result is not None:
                    raise SystemExit(result)
        except SystemExit as exit_exc:
            code = _normalize_exit_code(exit_exc.code)
            if exit_exc.code is not None and not isinstance(exit_exc.code, int):
                err_stream.write(f"{exit_exc.code}\n")
            if code and _is_vulnerability_finding_output(attempt_err.getvalue()):
                err_stream.write(
                    "[pip-audit] VULNERABILITY FINDING: pip-audit reported known "
                    "vulnerabilities; failing immediately without retry.\n"
                )
            return code
        except Exception as exc:
            if not is_transient_network_error(
                exc, service_connection_error=service_connection_error
            ):
                if _is_recognized_network_error(exc):
                    err_stream.write(
                        "[pip-audit] NETWORK ERROR (not retried): "
                        f"{_format_exception_detail(exc)}\n"
                    )
                    err_stream.flush()
                    return NETWORK_FAILURE_EXIT_CODE
                raise
            detail = _format_exception_detail(exc)
            if attempt < max_attempts:
                delay = initial_backoff * (backoff_factor ** (attempt - 1))
                err_stream.write(
                    f"[pip-audit] NETWORK ERROR (attempt {attempt}/{max_attempts}): "
                    f"{detail}. Retrying in {delay:.1f}s...\n"
                )
                err_stream.flush()
                sleep_fn(delay)
                continue
            err_stream.write(
                f"[pip-audit] NETWORK FAILURE: pip-audit failed after {max_attempts} "
                f"attempts due to transient network error ({detail}); "
                "not a vulnerability finding; audit did not complete.\n"
            )
            return NETWORK_FAILURE_EXIT_CODE
        else:
            return 0
    return NETWORK_FAILURE_EXIT_CODE


def main(argv: Sequence[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    try:
        return run_pip_audit(args)
    except PipAuditNotInstalledError as error:
        print(f"[pip-audit] INSTALLATION ERROR: {error}", file=sys.stderr)
        return 1
    except PipAuditCompatibilityError as error:
        print(f"[pip-audit] COMPATIBILITY ERROR: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
