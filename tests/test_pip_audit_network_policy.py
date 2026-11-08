"""Retry policy excludes permanent causes even inside connection wrappers (#1163)."""

from __future__ import annotations

import importlib
import io
from typing import Any

import pytest

from test_pip_audit_retry import _load_script


def _failure(kind: str) -> BaseException:
    requests: Any = importlib.import_module("requests")
    urllib3: Any = importlib.import_module("urllib3.exceptions")
    if kind == "ssl":
        error = requests.exceptions.SSLError("certificate validation failed")
    elif kind == "proxy":
        error = requests.exceptions.ProxyError("proxy configuration refused")
    elif kind == "redirects":
        error = requests.exceptions.TooManyRedirects("redirect loop")
    elif kind == "urllib3_ssl":
        error = urllib3.SSLError("certificate validation failed")
    elif kind == "urllib3_proxy":
        error = urllib3.ProxyError("proxy refused", OSError("proxy connection refused"))
    elif kind == "retry_ssl_reason":
        error = urllib3.MaxRetryError(None, "/feed", urllib3.SSLError("certificate failed"))
    elif kind == "retry_proxy_reason":
        error = urllib3.MaxRetryError(
            None, "/feed", urllib3.ProxyError("proxy refused", OSError("proxy refused"))
        )
    elif kind == "http_without_status":
        error = requests.exceptions.HTTPError("unknown HTTP result")
    else:
        response = requests.Response()
        response.status_code = int(kind.removeprefix("http_"))
        error = requests.exceptions.HTTPError("HTTP request failed", response=response)
    assert isinstance(error, BaseException)
    return error


@pytest.mark.parametrize(
    "kind",
    [
        "ssl",
        "proxy",
        "redirects",
        "urllib3_ssl",
        "urllib3_proxy",
        "retry_ssl_reason",
        "retry_proxy_reason",
        "http_without_status",
        "http_400",
        "http_403",
    ],
)
@pytest.mark.parametrize("wrapped", [False, True])
def test_permanent_causes_are_not_retried_even_when_connection_wrapper_comes_first(
    kind: str, wrapped: bool
) -> None:
    wrapper = _load_script("run_pip_audit")
    error = _failure(kind)
    if wrapped:
        interface: Any = importlib.import_module("pip_audit._service.interface")
        outer = interface.ConnectionError("feed connection failed")
        outer.__cause__ = error
        error = outer
    calls: list[list[str]] = []
    sleeps: list[float] = []
    output, diagnostic = io.StringIO(), io.StringIO()

    def audit(argv: list[str]) -> None:
        calls.append(argv)
        print("attempt evidence")
        raise error

    with pytest.raises(type(error)):
        wrapper.run_pip_audit(
            [],
            max_attempts=3,
            audit_fn=audit,
            sleep_fn=sleeps.append,
            stdout=output,
            stderr=diagnostic,
        )
    assert calls == [[]]
    assert sleeps == []
    assert output.getvalue() == "attempt evidence\n"
    assert "NETWORK ERROR" not in diagnostic.getvalue()
    assert "NETWORK FAILURE" not in diagnostic.getvalue()


@pytest.mark.parametrize("status", [408, 429, 500, 503, 599])
def test_recognized_transient_http_status_retries_without_losing_attempt_output(
    status: int,
) -> None:
    wrapper = _load_script("run_pip_audit")
    error = _failure(f"http_{status}")
    calls = 0
    sleeps: list[float] = []
    output, diagnostic = io.StringIO(), io.StringIO()

    def audit(_argv: list[str]) -> int:
        nonlocal calls
        calls += 1
        print(f"attempt {calls}")
        if calls == 1:
            raise error
        return 0

    assert (
        wrapper.run_pip_audit(
            [],
            max_attempts=3,
            audit_fn=audit,
            sleep_fn=sleeps.append,
            stdout=output,
            stderr=diagnostic,
        )
        == 0
    )
    assert calls == 2
    assert sleeps == [2.0]
    assert output.getvalue() == "attempt 1\nattempt 2\n"


def test_exception_reason_cycles_are_bounded_and_permanent_causes_take_precedence() -> None:
    wrapper = _load_script("run_pip_audit")
    urllib3: Any = importlib.import_module("urllib3.exceptions")
    root = urllib3.MaxRetryError(None, "/feed", None)
    nested = urllib3.MaxRetryError(None, "/feed", root)
    root.reason = nested
    nested.__cause__ = urllib3.SSLError("certificate failed")
    assert wrapper.is_transient_network_error(root) is False


def test_network_exhaustion_exit_code_is_distinct_from_findings_and_usage() -> None:
    wrapper = _load_script("run_pip_audit")
    assert wrapper.NETWORK_FAILURE_EXIT_CODE == 3
