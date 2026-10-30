"""Tests for transient-network retry and failure classification in pip-audit (#1113)."""

from __future__ import annotations

import importlib
import importlib.util
import io
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
from packaging.version import InvalidVersion, Version

REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_script(name: str) -> ModuleType:
    path = REPO_ROOT / "scripts" / f"{name}.py"
    if not path.is_file():
        raise FileNotFoundError(f"Missing script: {path}")
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def _pip_audit_types() -> tuple[Any, Any, Any, Any, Any]:
    requests_mod = importlib.import_module("requests")
    interface_mod = importlib.import_module("pip_audit._service.interface")
    pypi_mod = importlib.import_module("pip_audit._service.pypi")
    return (
        requests_mod,
        pypi_mod.PyPIService,
        interface_mod.ConnectionError,
        interface_mod.ServiceError,
        interface_mod.VulnerabilityResult,
    )


def _write_pinned_requirements(tmp_path: Path) -> Path:
    req_file = tmp_path / "requirements.txt"
    req_file.write_text("demo-pkg==1.0.0\n", encoding="utf-8")
    return req_file


def test_retries_transient_network_errors_with_backoff_then_succeeds(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    run_pip_audit = _load_script("run_pip_audit")
    requests_mod, pypi_service_cls, vuln_conn_error_cls, _, _ = _pip_audit_types()
    req_file = _write_pinned_requirements(tmp_path)
    cache_dir = tmp_path / "cache"
    attempts = 0
    sleeps: list[float] = []

    def flaky_query(_self: Any, spec: Any) -> tuple[Any, list[Any]]:
        nonlocal attempts
        attempts += 1
        if attempts == 1:
            raise requests_mod.exceptions.ConnectionError(
                "('Connection aborted.', "
                "RemoteDisconnected('Remote end closed connection without response'))"
            )
        if attempts == 2:
            raise vuln_conn_error_cls("Could not connect to PyPI's vulnerability feed")
        return spec, []

    monkeypatch.setattr(pypi_service_cls, "query", flaky_query)
    stdout = io.StringIO()
    stderr = io.StringIO()

    exit_code = run_pip_audit.run_pip_audit(
        [
            "-r",
            str(req_file),
            "--disable-pip",
            "--no-deps",
            "--progress-spinner",
            "off",
            "--cache-dir",
            str(cache_dir),
        ],
        max_attempts=4,
        initial_backoff=2.0,
        backoff_factor=2.0,
        sleep_fn=sleeps.append,
        stdout=stdout,
        stderr=stderr,
    )

    assert exit_code == 0
    assert attempts == 3
    assert sleeps == [2.0, 4.0]
    err_text = stderr.getvalue()
    assert "NETWORK ERROR (attempt 1/4)" in err_text
    assert "RemoteDisconnected" in err_text
    assert "NETWORK ERROR (attempt 2/4)" in err_text
    assert "Could not connect to PyPI's vulnerability feed" in err_text
    assert "No known vulnerabilities found" in err_text


def test_retries_transient_http_5xx_service_error_then_succeeds(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    run_pip_audit = _load_script("run_pip_audit")
    requests_mod, pypi_service_cls, _, service_error_cls, _ = _pip_audit_types()
    req_file = _write_pinned_requirements(tmp_path)
    cache_dir = tmp_path / "cache"
    attempts = 0
    sleeps: list[float] = []

    def flaky_503_query(_self: Any, spec: Any) -> tuple[Any, list[Any]]:
        nonlocal attempts
        attempts += 1
        if attempts == 1:
            response = requests_mod.Response()
            response.status_code = 503
            http_error = requests_mod.HTTPError(
                "503 Server Error: Service Unavailable", response=response
            )
            raise service_error_cls from http_error
        return spec, []

    monkeypatch.setattr(pypi_service_cls, "query", flaky_503_query)
    stdout = io.StringIO()
    stderr = io.StringIO()

    exit_code = run_pip_audit.run_pip_audit(
        [
            "-r",
            str(req_file),
            "--disable-pip",
            "--no-deps",
            "--progress-spinner",
            "off",
            "--cache-dir",
            str(cache_dir),
        ],
        max_attempts=3,
        initial_backoff=1.5,
        backoff_factor=2.0,
        sleep_fn=sleeps.append,
        stdout=stdout,
        stderr=stderr,
    )

    assert exit_code == 0
    assert attempts == 2
    assert sleeps == [1.5]
    assert "NETWORK ERROR (attempt 1/3)" in stderr.getvalue()


def test_fails_immediately_on_vulnerability_finding_without_retry(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    run_pip_audit = _load_script("run_pip_audit")
    _, pypi_service_cls, _, _, vuln_result_cls = _pip_audit_types()
    req_file = _write_pinned_requirements(tmp_path)
    cache_dir = tmp_path / "cache"
    attempts = 0
    sleeps: list[float] = []

    def vulnerable_query(_self: Any, spec: Any) -> tuple[Any, list[Any]]:
        nonlocal attempts
        attempts += 1
        return (
            spec,
            [
                vuln_result_cls.create(
                    ids=["PYSEC-2026-9999"],
                    description="Synthetic test vulnerability",
                    fix_versions=[Version("1.0.1")],
                    published=None,
                )
            ],
        )

    monkeypatch.setattr(pypi_service_cls, "query", vulnerable_query)
    stdout = io.StringIO()
    stderr = io.StringIO()

    exit_code = run_pip_audit.run_pip_audit(
        [
            "-r",
            str(req_file),
            "--disable-pip",
            "--no-deps",
            "--progress-spinner",
            "off",
            "--cache-dir",
            str(cache_dir),
        ],
        max_attempts=4,
        initial_backoff=2.0,
        backoff_factor=2.0,
        sleep_fn=sleeps.append,
        stdout=stdout,
        stderr=stderr,
    )

    assert exit_code == 1
    assert attempts == 1
    assert sleeps == []
    assert "PYSEC-2026-9999" in stdout.getvalue()
    err_text = stderr.getvalue()
    assert "Found 1 known vulnerability in 1 package" in err_text
    assert "VULNERABILITY FINDING" in err_text
    assert "NETWORK FAILURE" not in err_text
    assert "NETWORK ERROR" not in err_text


def test_reports_persistent_network_failure_distinctly_from_vulnerability_finding(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    run_pip_audit = _load_script("run_pip_audit")
    requests_mod, pypi_service_cls, _, _, _ = _pip_audit_types()
    req_file = _write_pinned_requirements(tmp_path)
    cache_dir = tmp_path / "cache"
    attempts = 0
    sleeps: list[float] = []

    def always_disconnected(_self: Any, _spec: Any) -> tuple[Any, list[Any]]:
        nonlocal attempts
        attempts += 1
        raise requests_mod.exceptions.ConnectionError(
            "('Connection aborted.', "
            "RemoteDisconnected('Remote end closed connection without response'))"
        )

    monkeypatch.setattr(pypi_service_cls, "query", always_disconnected)
    stdout = io.StringIO()
    stderr = io.StringIO()

    exit_code = run_pip_audit.run_pip_audit(
        [
            "-r",
            str(req_file),
            "--disable-pip",
            "--no-deps",
            "--progress-spinner",
            "off",
            "--cache-dir",
            str(cache_dir),
        ],
        max_attempts=3,
        initial_backoff=2.0,
        backoff_factor=2.0,
        sleep_fn=sleeps.append,
        stdout=stdout,
        stderr=stderr,
    )

    assert exit_code == run_pip_audit.NETWORK_FAILURE_EXIT_CODE
    assert exit_code != 1
    assert attempts == 3
    assert sleeps == [2.0, 4.0]
    err_text = stderr.getvalue()
    assert "NETWORK FAILURE" in err_text
    assert "after 3 attempts" in err_text
    assert "not a vulnerability finding" in err_text
    assert "VULNERABILITY FINDING" not in err_text
    assert "Found " not in err_text


def test_does_not_retry_non_network_service_or_source_errors(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    run_pip_audit = _load_script("run_pip_audit")
    requests_mod, pypi_service_cls, _, service_error_cls, _ = _pip_audit_types()
    req_file = _write_pinned_requirements(tmp_path)
    cache_dir = tmp_path / "cache"
    attempts = 0
    sleeps: list[float] = []

    def malformed_version_query(_self: Any, _spec: Any) -> tuple[Any, list[Any]]:
        nonlocal attempts
        attempts += 1
        try:
            Version("not-a-version")
        except InvalidVersion as exc:
            raise service_error_cls("Received malformed version from PyPI") from exc
        raise AssertionError("unreachable")

    monkeypatch.setattr(pypi_service_cls, "query", malformed_version_query)
    stdout = io.StringIO()
    stderr = io.StringIO()

    with pytest.raises(service_error_cls, match="Received malformed version from PyPI"):
        run_pip_audit.run_pip_audit(
            [
                "-r",
                str(req_file),
                "--disable-pip",
                "--no-deps",
                "--progress-spinner",
                "off",
                "--cache-dir",
                str(cache_dir),
            ],
            max_attempts=4,
            sleep_fn=sleeps.append,
            stdout=stdout,
            stderr=stderr,
        )
    assert attempts == 1
    assert sleeps == []

    http_400_attempts = 0

    def http_400_query(_self: Any, _spec: Any) -> tuple[Any, list[Any]]:
        nonlocal http_400_attempts
        http_400_attempts += 1
        response = requests_mod.Response()
        response.status_code = 400
        raise service_error_cls from requests_mod.HTTPError(
            "400 Client Error: Bad Request", response=response
        )

    monkeypatch.setattr(pypi_service_cls, "query", http_400_query)
    with pytest.raises(service_error_cls):
        run_pip_audit.run_pip_audit(
            [
                "-r",
                str(req_file),
                "--disable-pip",
                "--no-deps",
                "--progress-spinner",
                "off",
                "--cache-dir",
                str(cache_dir),
            ],
            max_attempts=4,
            sleep_fn=sleeps.append,
            stdout=io.StringIO(),
            stderr=io.StringIO(),
        )
    assert http_400_attempts == 1
    assert sleeps == []

    missing_req = tmp_path / "missing.txt"
    missing_stdout = io.StringIO()
    missing_stderr = io.StringIO()
    exit_code = run_pip_audit.run_pip_audit(
        ["-r", str(missing_req), "--progress-spinner", "off", "--cache-dir", str(cache_dir)],
        max_attempts=4,
        sleep_fn=sleeps.append,
        stdout=missing_stdout,
        stderr=missing_stderr,
    )
    assert exit_code == 1
    assert sleeps == []
    assert "NETWORK FAILURE" not in missing_stderr.getvalue()
    assert "VULNERABILITY FINDING" not in missing_stderr.getvalue()


def test_nonzero_return_from_audit_function_gives_nonzero_exit(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    run_pip_audit = _load_script("run_pip_audit")
    cli_mod: Any = importlib.import_module("pip_audit._cli")
    cache_dir = tmp_path / "cache"

    stub_attempts = 0
    sleeps: list[float] = []

    def returning_audit_fn(_argv: Any) -> int:
        nonlocal stub_attempts
        stub_attempts += 1
        return 1

    stub_exit = run_pip_audit.run_pip_audit(
        ["--progress-spinner", "off", "--cache-dir", str(cache_dir)],
        max_attempts=4,
        audit_fn=returning_audit_fn,
        sleep_fn=sleeps.append,
        stdout=io.StringIO(),
        stderr=io.StringIO(),
    )
    assert stub_exit == 1
    assert stub_attempts == 1
    assert sleeps == []

    cli_attempts = 0

    def cli_audit_returning_one() -> int:
        nonlocal cli_attempts
        cli_attempts += 1
        return 1

    monkeypatch.setattr(cli_mod, "audit", cli_audit_returning_one)
    cli_exit = run_pip_audit.run_pip_audit(
        ["--progress-spinner", "off", "--cache-dir", str(cache_dir)],
        max_attempts=4,
        sleep_fn=sleeps.append,
        stdout=io.StringIO(),
        stderr=io.StringIO(),
    )
    assert cli_exit == 1
    assert cli_attempts == 1
    assert sleeps == []


def test_restores_vuln_service_connection_error_after_normal_and_exception_runs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    run_pip_audit = _load_script("run_pip_audit")
    _, pypi_service_cls, vuln_conn_error_cls, service_error_cls, _ = _pip_audit_types()
    cli_mod: Any = importlib.import_module("pip_audit._cli")
    interface_mod: Any = importlib.import_module("pip_audit._service.interface")
    req_file = _write_pinned_requirements(tmp_path)
    cache_dir = tmp_path / "cache"
    argv = [
        "-r",
        str(req_file),
        "--disable-pip",
        "--no-deps",
        "--progress-spinner",
        "off",
        "--cache-dir",
        str(cache_dir),
    ]

    observed_during_run: list[Any] = []

    def clean_query(_self: Any, spec: Any) -> tuple[Any, list[Any]]:
        observed_during_run.append(cli_mod.VulnServiceConnectionError)
        return spec, []

    monkeypatch.setattr(pypi_service_cls, "query", clean_query)
    cli_mod.VulnServiceConnectionError = interface_mod.ConnectionError
    exit_code = run_pip_audit.run_pip_audit(
        argv,
        max_attempts=2,
        stdout=io.StringIO(),
        stderr=io.StringIO(),
    )
    assert exit_code == 0
    assert observed_during_run == [run_pip_audit._UnmatchedVulnServiceConnectionError]
    assert cli_mod.VulnServiceConnectionError is interface_mod.ConnectionError
    assert cli_mod.VulnServiceConnectionError is vuln_conn_error_cls

    def exploding_query(_self: Any, _spec: Any) -> tuple[Any, list[Any]]:
        observed_during_run.append(cli_mod.VulnServiceConnectionError)
        raise service_error_cls("Non-network failure")

    monkeypatch.setattr(pypi_service_cls, "query", exploding_query)
    cli_mod.VulnServiceConnectionError = interface_mod.ConnectionError
    with pytest.raises(service_error_cls, match="Non-network failure"):
        run_pip_audit.run_pip_audit(
            argv,
            max_attempts=2,
            stdout=io.StringIO(),
            stderr=io.StringIO(),
        )
    assert cli_mod.VulnServiceConnectionError is interface_mod.ConnectionError
    assert cli_mod.VulnServiceConnectionError is vuln_conn_error_cls

    # Also verify recovery when VulnServiceConnectionError was dirty prior to invocation.
    cli_mod.VulnServiceConnectionError = Exception
    with pytest.raises(service_error_cls, match="Non-network failure"):
        run_pip_audit.run_pip_audit(
            argv,
            max_attempts=2,
            stdout=io.StringIO(),
            stderr=io.StringIO(),
        )
    assert cli_mod.VulnServiceConnectionError is interface_mod.ConnectionError


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"max_attempts": 0}, "max_attempts"),
        ({"initial_backoff": -0.1}, "initial_backoff"),
        ({"backoff_factor": 0.5}, "backoff_factor"),
    ],
)
def test_rejects_invalid_retry_parameters(kwargs: dict[str, Any], match: str) -> None:
    run_pip_audit = _load_script("run_pip_audit")
    with pytest.raises(ValueError, match=match):
        run_pip_audit.run_pip_audit([], **kwargs)


def test_justfile_and_pre_commit_invoke_retry_wrapper_without_loosening_audit_flags() -> None:
    justfile_text = (REPO_ROOT / "Justfile").read_text(encoding="utf-8")
    assert (
        'uv run python scripts/run_pip_audit.py -r "$reqfile" --disable-pip \\\n'
        '\t  --cache-dir "$tmpdir/pip-audit-cache"'
    ) in justfile_text
    assert (
        "uv run --locked --group qa python scripts/run_pip_audit.py "
        '-r "$audit_tmp/requirements-linux.txt" \\\n'
        '\t  --disable-pip --no-deps --cache-dir "$audit_tmp/pip-audit-cache"'
    ) in justfile_text

    precommit_text = (REPO_ROOT / ".pre-commit-config.yaml").read_text(encoding="utf-8")
    assert "entry: uv run python scripts/run_pip_audit.py" in precommit_text
