"""Retry policy excludes permanent causes even inside connection wrappers (#1163)."""

from __future__ import annotations

import importlib
import io
import ssl
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest

from test_pip_audit_retry import (
    REPO_ROOT,
    _load_script,
    _pip_audit_types,
    _write_pinned_requirements,
)


def _reported_exit(call: Callable[[], int]) -> int:
    """A recognized network refusal must return a label/code, never escape as a traceback."""
    try:
        return call()
    except Exception as error:
        raise AssertionError(
            "Recognized network error escaped instead of returning a labelled exit"
        ) from error


def _tls_transport_wrapper(error: BaseException) -> BaseException:
    """Use real exception-valued args/reason nesting, without a cause shortcut."""
    requests: Any = importlib.import_module("requests")
    urllib3: Any = importlib.import_module("urllib3.exceptions")
    wrapped = requests.exceptions.SSLError(
        urllib3.MaxRetryError(None, "/feed", urllib3.SSLError(error))
    )
    assert isinstance(wrapped, BaseException)
    assert wrapped.__cause__ is None
    assert wrapped.__context__ is None
    return wrapped


def _failure(kind: str) -> BaseException:
    requests: Any = importlib.import_module("requests")
    urllib3: Any = importlib.import_module("urllib3.exceptions")
    if kind == "ssl":
        error = requests.exceptions.SSLError(ssl.SSLCertVerificationError(1, "issuer refused"))
    elif kind == "certificate":
        error = ssl.SSLCertVerificationError(1, "issuer refused")
    elif kind == "hostname":
        error = ssl.CertificateError("hostname does not match")
    elif kind == "nested_certificate":
        error = _tls_transport_wrapper(ssl.SSLCertVerificationError(1, "issuer refused"))
    elif kind == "proxy":
        error = requests.exceptions.ProxyError("proxy configuration refused")
    elif kind == "redirects":
        error = requests.exceptions.TooManyRedirects("redirect loop")
    elif kind == "urllib3_ssl":
        error = urllib3.SSLError(ssl.SSLCertVerificationError(1, "issuer refused"))
    elif kind == "urllib3_proxy":
        error = urllib3.ProxyError("proxy refused", OSError("proxy connection refused"))
    elif kind == "retry_ssl_reason":
        error = urllib3.MaxRetryError(
            None, "/feed", urllib3.SSLError(ssl.SSLCertVerificationError(1, "issuer refused"))
        )
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
        "certificate",
        "hostname",
        "nested_certificate",
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

    assert wrapper.is_transient_network_error(error) is False
    exit_code = _reported_exit(
        lambda: wrapper.run_pip_audit(
            [],
            max_attempts=3,
            audit_fn=audit,
            sleep_fn=sleeps.append,
            stdout=output,
            stderr=diagnostic,
        )
    )
    assert exit_code == wrapper.NETWORK_FAILURE_EXIT_CODE
    assert exit_code not in {0, 1, 2}
    assert diagnostic.getvalue().startswith("[pip-audit] NETWORK ERROR (not retried): ")
    assert type(error).__name__ in diagnostic.getvalue()
    assert "Traceback" not in diagnostic.getvalue()
    assert "VULNERABILITY FINDING" not in diagnostic.getvalue()
    assert calls == [[]]
    assert sleeps == []
    assert output.getvalue() == "attempt evidence\n"
    assert "NETWORK FAILURE" not in diagnostic.getvalue()


@pytest.mark.parametrize("error_class", [ssl.SSLEOFError, ssl.SSLZeroReturnError])
@pytest.mark.parametrize("wrapped", [False, True])
def test_typed_tls_closure_retries_without_losing_attempt_output(
    error_class: type[ssl.SSLError], wrapped: bool
) -> None:
    wrapper = _load_script("run_pip_audit")
    error: BaseException = error_class(8, "peer closed TLS while reading")
    if wrapped:
        error = _tls_transport_wrapper(error)
    assert wrapper.is_transient_network_error(error) is True
    calls: list[list[str]] = []
    sleeps: list[float] = []
    output, diagnostic = io.StringIO(), io.StringIO()

    def audit(argv: list[str]) -> int:
        calls.append(argv)
        print(f"attempt {len(calls)}")
        print(f"stderr {len(calls)}", file=sys.stderr)
        if len(calls) == 1:
            raise error
        return 0

    assert (
        wrapper.run_pip_audit(
            ["--synthetic"],
            max_attempts=2,
            audit_fn=audit,
            sleep_fn=sleeps.append,
            stdout=output,
            stderr=diagnostic,
        )
        == 0
    )
    assert calls == [["--synthetic"], ["--synthetic"]]
    assert sleeps == [2.0]
    assert output.getvalue() == "attempt 1\nattempt 2\n"
    text = diagnostic.getvalue()
    assert (
        text.index("stderr 1") < text.index("NETWORK ERROR (attempt 1/2)") < text.index("stderr 2")
    )
    assert error_class.__name__ in text
    assert "not retried" not in text


@pytest.mark.parametrize("kind", ["certificate", "eof", "zero_return"])
def test_real_cli_typed_tls_policy_preserves_bindings_and_outputs(
    kind: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    wrapper = _load_script("run_pip_audit")
    _, pypi_service, _, _, _ = _pip_audit_types()
    cli: Any = importlib.import_module("pip_audit._cli")
    original_alias, original_argv = cli.VulnServiceConnectionError, sys.argv
    error_class = {
        "certificate": ssl.SSLCertVerificationError,
        "eof": ssl.SSLEOFError,
        "zero_return": ssl.SSLZeroReturnError,
    }[kind]
    failure = _tls_transport_wrapper(error_class(1, "controlled typed TLS fault"))
    calls = 0
    sleeps: list[float] = []
    output, diagnostic = io.StringIO(), io.StringIO()

    def query(_self: Any, spec: Any) -> tuple[Any, list[Any]]:
        nonlocal calls
        calls += 1
        print(f"query {calls}")
        if calls == 1:
            raise failure
        return spec, []

    monkeypatch.setattr(pypi_service, "query", query)
    result = _reported_exit(
        lambda: wrapper.run_pip_audit(
            [
                "-r",
                str(_write_pinned_requirements(tmp_path)),
                "--disable-pip",
                "--no-deps",
                "--progress-spinner",
                "off",
                "--cache-dir",
                str(tmp_path / "cache"),
            ],
            max_attempts=2,
            sleep_fn=sleeps.append,
            stdout=output,
            stderr=diagnostic,
        )
    )
    assert cli.VulnServiceConnectionError is original_alias
    assert sys.argv is original_argv
    assert "query 1" in output.getvalue()
    if kind == "certificate":
        assert result == wrapper.NETWORK_FAILURE_EXIT_CODE
        assert calls == 1
        assert sleeps == []
        assert "NETWORK ERROR (not retried):" in diagnostic.getvalue()
    else:
        assert result == 0
        assert calls == 2
        assert sleeps == [2.0]
        assert "query 2" in output.getvalue()
        assert "NETWORK ERROR (attempt 1/2)" in diagnostic.getvalue()
    assert "Traceback" not in diagnostic.getvalue()


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
    nested.__cause__ = urllib3.SSLError(ssl.SSLCertVerificationError(1, "issuer refused"))
    assert wrapper.is_transient_network_error(root) is False


def test_exception_argument_cycles_are_bounded_and_typed_certificate_wins() -> None:
    wrapper = _load_script("run_pip_audit")
    requests: Any = importlib.import_module("requests")
    outer = requests.exceptions.ConnectionError("feed unavailable")
    nested = RuntimeError(outer, ssl.SSLCertVerificationError(1, "issuer refused"))
    outer.args = (nested,)
    assert wrapper.is_transient_network_error(outer) is False


def test_certificate_words_do_not_classify_an_unknown_programming_error() -> None:
    wrapper = _load_script("run_pip_audit")
    calls: list[str] = []
    sleeps: list[float] = []
    diagnostic = io.StringIO()

    def audit(_argv: list[str]) -> None:
        calls.append("audit")
        raise ValueError("certificate proxy redirect TLS EOF: malformed application value")

    with pytest.raises(ValueError, match="malformed application value"):
        wrapper.run_pip_audit(
            [],
            audit_fn=audit,
            max_attempts=2,
            sleep_fn=sleeps.append,
            stderr=diagnostic,
        )
    assert calls == ["audit"]
    assert sleeps == []
    assert "NETWORK" not in diagnostic.getvalue()


def test_bare_unclassified_ssl_error_is_not_guessed_transient_or_certificate() -> None:
    wrapper = _load_script("run_pip_audit")
    error = ssl.SSLError(1, "certificate failed is only untyped diagnostic text")
    assert wrapper.is_transient_network_error(error) is False
    calls: list[str] = []
    sleeps: list[float] = []
    diagnostic = io.StringIO()

    def audit(_argv: list[str]) -> None:
        calls.append("audit")
        raise error

    code = _reported_exit(
        lambda: wrapper.run_pip_audit(
            [],
            audit_fn=audit,
            sleep_fn=sleeps.append,
            stderr=diagnostic,
        )
    )
    assert code == wrapper.NETWORK_FAILURE_EXIT_CODE
    assert calls == ["audit"]
    assert sleeps == []
    assert "NETWORK ERROR (not retried):" in diagnostic.getvalue()
    assert "SSLError:" in diagnostic.getvalue()


def test_retry_policy_documentation_names_typed_tls_and_missing_installation() -> None:
    text = " ".join((REPO_ROOT / "docs/pip_audit_retry.md").read_text().split())
    for name in (
        "SSLCertVerificationError",
        "CertificateError",
        "SSLEOFError",
        "SSLZeroReturnError",
    ):
        assert name in text
    assert "[pip-audit] NETWORK ERROR (not retried):" in text
    assert "not installed" in text
    assert "just sync" in text
    assert "TLS/SSL, proxy and redirect errors do not retry" not in text
    assert "exception-valued" in text
    assert "unclassified bare" in text


def test_network_exhaustion_exit_code_is_distinct_from_findings_and_usage() -> None:
    wrapper = _load_script("run_pip_audit")
    assert wrapper.NETWORK_FAILURE_EXIT_CODE == 3
