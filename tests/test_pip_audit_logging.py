"""Console capture preserves shared handlers, prefixes and nested state (#1163)."""

from __future__ import annotations

import io
import logging
import sys
from collections.abc import Sequence
from contextlib import closing
from pathlib import Path

import pytest

from test_pip_audit_retry import _load_script


def _configure_loggers(
    monkeypatch: pytest.MonkeyPatch,
    root_handlers: list[logging.Handler],
    package_handlers: list[logging.Handler],
) -> tuple[logging.Logger, logging.Logger]:
    root = logging.getLogger()
    package = logging.getLogger("pip_audit")
    monkeypatch.setattr(root, "handlers", root_handlers)
    monkeypatch.setattr(package, "handlers", package_handlers)
    monkeypatch.setattr(root, "propagate", True)
    monkeypatch.setattr(package, "propagate", False)
    for logger in (root, package):
        monkeypatch.setattr(logger, "level", logging.WARNING)
        monkeypatch.setattr(logger, "disabled", False)
        # Isolate cached level decisions along with each logger's configured state.
        monkeypatch.setattr(logger, "_cache", {})
    return root, package


def test_shared_root_and_package_console_handler_emits_each_warning_once(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    wrapper = _load_script("run_pip_audit")
    stderr = io.StringIO()
    handler = logging.StreamHandler(stderr)
    formatter = logging.Formatter("CUSTOM %(levelname)s %(name)s | %(message)s")
    handler.setFormatter(formatter)
    root_handlers: list[logging.Handler] = [handler]
    package_handlers: list[logging.Handler] = [handler]
    root, package = _configure_loggers(monkeypatch, root_handlers, package_handlers)

    def audit(_argv: Sequence[str]) -> int:
        package.warning("shared-handler-marker")
        return 0

    assert (
        wrapper.run_pip_audit(
            [], audit_fn=audit, max_attempts=1, stdout=io.StringIO(), stderr=stderr
        )
        == 0
    )
    assert stderr.getvalue() == "CUSTOM WARNING pip_audit | shared-handler-marker\n"
    assert handler.stream is stderr
    assert handler.formatter is formatter
    assert root.handlers is root_handlers
    assert package.handlers is package_handlers
    assert root.handlers == [handler]
    assert package.handlers == [handler]
    assert package.propagate is False


def test_preexisting_root_console_retargets_to_explicit_stderr_and_restores_stream(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    wrapper = _load_script("run_pip_audit")
    previous_stderr = sys.stderr
    supplied_stderr = io.StringIO()
    handler = logging.StreamHandler(previous_stderr)
    formatter = logging.Formatter("ORIGINAL %(levelname)s %(name)s | %(message)s")
    handler.setFormatter(formatter)
    root_handlers: list[logging.Handler] = [handler]
    package_handlers: list[logging.Handler] = []
    root, package = _configure_loggers(monkeypatch, root_handlers, package_handlers)

    def audit(_argv: Sequence[str]) -> int:
        root.warning("explicit-destination-marker")
        return 0

    assert (
        wrapper.run_pip_audit(
            [], audit_fn=audit, max_attempts=1, stdout=io.StringIO(), stderr=supplied_stderr
        )
        == 0
    )
    assert supplied_stderr.getvalue() == "ORIGINAL WARNING root | explicit-destination-marker\n"
    assert handler.stream is previous_stderr
    assert handler.formatter is formatter
    assert sys.stderr is previous_stderr
    assert root.handlers is root_handlers
    assert package.handlers is package_handlers
    assert root.handlers == [handler]
    assert package.handlers == []
    assert package.propagate is False


@pytest.mark.parametrize("outcome", ["success", "system_exit", "ordinary_error"])
def test_custom_prefix_and_all_logger_state_survive_every_exit(
    monkeypatch: pytest.MonkeyPatch, outcome: str
) -> None:
    wrapper = _load_script("run_pip_audit")
    stdout, stderr = io.StringIO(), io.StringIO()
    root_handler, package_handler = logging.StreamHandler(stdout), logging.StreamHandler(stderr)
    root_formatter = logging.Formatter("ROOT %(levelname)s | %(message)s")
    package_formatter = logging.Formatter("AUDIT %(levelname)s %(name)s :: %(message)s")
    root_handler.setFormatter(root_formatter)
    package_handler.setFormatter(package_formatter)
    root_handlers: list[logging.Handler] = [root_handler]
    package_handlers: list[logging.Handler] = [package_handler]
    root, package = _configure_loggers(monkeypatch, root_handlers, package_handlers)
    previous_stdout, previous_stderr = sys.stdout, sys.stderr

    def audit(_argv: Sequence[str]) -> int:
        root.warning("root-during")
        package.warning("package-during")
        if outcome == "system_exit":
            raise SystemExit(7)
        if outcome == "ordinary_error":
            raise ValueError("Synthetic permanent logging-path error")
        return 0

    if outcome == "ordinary_error":
        with pytest.raises(ValueError, match="Synthetic permanent logging-path error"):
            wrapper.run_pip_audit([], audit_fn=audit, max_attempts=1, stdout=stdout, stderr=stderr)
    else:
        expected = 7 if outcome == "system_exit" else 0
        assert (
            wrapper.run_pip_audit([], audit_fn=audit, max_attempts=1, stdout=stdout, stderr=stderr)
            == expected
        )

    assert sys.stdout is previous_stdout
    assert sys.stderr is previous_stderr
    assert root.handlers is root_handlers
    assert package.handlers is package_handlers
    assert root.handlers == [root_handler]
    assert package.handlers == [package_handler]
    assert root_handler.stream is stdout
    assert package_handler.stream is stderr
    assert root_handler.formatter is root_formatter
    assert package_handler.formatter is package_formatter
    assert root.level == logging.WARNING
    assert package.level == logging.WARNING
    assert root.disabled is False
    assert package.disabled is False
    assert root.propagate is True
    assert package.propagate is False
    package.warning("package-after")
    assert stdout.getvalue() == "ROOT WARNING | root-during\n"
    assert stderr.getvalue() == (
        "AUDIT WARNING pip_audit :: package-during\nAUDIT WARNING pip_audit :: package-after\n"
    )


class _RecordHandler(logging.Handler):
    def __init__(self) -> None:
        super().__init__()
        self.records: list[logging.LogRecord] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.records.append(record)


def test_fallback_console_preserves_file_and_nonstream_capture_handlers(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    wrapper = _load_script("run_pip_audit")
    stderr = io.StringIO()
    path = tmp_path / "retained-log.txt"
    with (
        closing(logging.FileHandler(path, encoding="utf-8")) as file_handler,
        closing(_RecordHandler()) as capture_handler,
    ):
        file_handler.setFormatter(logging.Formatter("FILE %(name)s | %(message)s"))
        root_handlers: list[logging.Handler] = [file_handler]
        package_handlers: list[logging.Handler] = [capture_handler]
        root, package = _configure_loggers(monkeypatch, root_handlers, package_handlers)
        previous_file_stream = file_handler.stream
        assert previous_file_stream is not None

        def audit(_argv: Sequence[str]) -> int:
            root.warning("root-fallback-marker")
            package.warning("package-fallback-marker")
            return 0

        assert (
            wrapper.run_pip_audit(
                [], audit_fn=audit, max_attempts=1, stdout=io.StringIO(), stderr=stderr
            )
            == 0
        )
        assert stderr.getvalue() == (
            "WARNING:root:root-fallback-marker\nWARNING:pip_audit:package-fallback-marker\n"
        )
        assert path.read_text(encoding="utf-8") == "FILE root | root-fallback-marker\n"
        assert [record.getMessage() for record in capture_handler.records] == [
            "package-fallback-marker"
        ]
        assert root.handlers is root_handlers
        assert package.handlers is package_handlers
        assert root.handlers == [file_handler]
        assert package.handlers == [capture_handler]
        assert file_handler.stream is previous_file_stream
        assert not previous_file_stream.closed
        assert package.propagate is False


def test_nested_wrapper_restores_outer_shared_handler_and_live_streams(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    wrapper = _load_script("run_pip_audit")
    stdout, stderr = io.StringIO(), io.StringIO()
    handler = logging.StreamHandler(stderr)
    formatter = logging.Formatter("NEST %(name)s | %(message)s")
    handler.setFormatter(formatter)
    root_handlers: list[logging.Handler] = [handler]
    package_handlers: list[logging.Handler] = [handler]
    root, package = _configure_loggers(monkeypatch, root_handlers, package_handlers)
    previous_stdout, previous_stderr = sys.stdout, sys.stderr

    def inner(_argv: Sequence[str]) -> int:
        package.warning("inner")
        return 0

    def outer(_argv: Sequence[str]) -> int:
        outer_handler_stream = handler.stream
        outer_stdout, outer_stderr = sys.stdout, sys.stderr
        package.warning("outer-before")
        assert (
            wrapper.run_pip_audit([], audit_fn=inner, max_attempts=1, stdout=stdout, stderr=stderr)
            == 0
        )
        assert handler.stream is outer_handler_stream
        assert sys.stdout is outer_stdout
        assert sys.stderr is outer_stderr
        assert root.handlers is root_handlers
        assert package.handlers is package_handlers
        package.warning("outer-after")
        return 0

    assert (
        wrapper.run_pip_audit([], audit_fn=outer, max_attempts=1, stdout=stdout, stderr=stderr) == 0
    )
    assert stderr.getvalue() == (
        "NEST pip_audit | outer-before\nNEST pip_audit | inner\nNEST pip_audit | outer-after\n"
    )
    assert handler.stream is stderr
    assert handler.formatter is formatter
    assert sys.stdout is previous_stdout
    assert sys.stderr is previous_stderr
    assert root.handlers is root_handlers
    assert package.handlers is package_handlers
    assert root.handlers == [handler]
    assert package.handlers == [handler]
    assert package.propagate is False
