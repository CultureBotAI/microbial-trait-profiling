"""The in-process pip-audit adapter restores bindings and refuses drift (#1185)."""

from __future__ import annotations

import importlib
import importlib.metadata
import importlib.util
import io
import sys
import tomllib
from importlib.machinery import ModuleSpec
from types import ModuleType
from typing import Any

import pytest

from test_pip_audit_retry import REPO_ROOT, _load_script


class PreviousConnectionError(Exception):
    """A caller-owned prior binding, deliberately different from pip-audit's default."""


def _installed_api() -> tuple[Any, Any]:
    return (
        importlib.import_module("pip_audit._cli"),
        importlib.import_module("pip_audit._service.interface"),
    )


@pytest.mark.parametrize("outcome", ["success", "system_exit", "service_error"])
def test_exact_prior_class_and_argv_binding_survive_every_adapter_outcome(
    monkeypatch: pytest.MonkeyPatch, outcome: str
) -> None:
    wrapper = _load_script("run_pip_audit")
    cli, interface = _installed_api()
    original_argv = ["caller-owned-program", "keep-this-argument"]
    monkeypatch.setattr(sys, "argv", original_argv)
    monkeypatch.setattr(cli, "VulnServiceConnectionError", PreviousConnectionError)
    seen: list[list[str]] = []
    sleeps: list[float] = []

    def audit() -> int:
        assert cli.VulnServiceConnectionError is wrapper._UnmatchedVulnServiceConnectionError
        assert sys.argv is not original_argv
        seen.append(sys.argv)
        if outcome == "system_exit":
            raise SystemExit(7)
        if outcome == "service_error":
            raise interface.ServiceError("Synthetic permanent audit failure")
        return 0

    monkeypatch.setattr(cli, "audit", audit)
    kwargs = {
        "max_attempts": 2,
        "sleep_fn": sleeps.append,
        "stdout": io.StringIO(),
        "stderr": io.StringIO(),
    }
    if outcome == "service_error":
        with pytest.raises(interface.ServiceError, match="Synthetic permanent audit failure"):
            wrapper.run_pip_audit(["--synthetic-case"], **kwargs)
    else:
        expected = 7 if outcome == "system_exit" else 0
        assert wrapper.run_pip_audit(["--synthetic-case"], **kwargs) == expected

    assert seen == [["pip-audit", "--synthetic-case"]]
    assert sleeps == []
    assert cli.VulnServiceConnectionError is PreviousConnectionError
    assert sys.argv is original_argv
    assert original_argv == ["caller-owned-program", "keep-this-argument"]


def test_nested_adapter_call_preserves_outer_binding_and_argv_identity(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    wrapper = _load_script("run_pip_audit")
    cli, _ = _installed_api()
    original_argv = ["caller-owned-program"]
    monkeypatch.setattr(sys, "argv", original_argv)
    monkeypatch.setattr(cli, "VulnServiceConnectionError", PreviousConnectionError)
    calls: list[str] = []

    def audit() -> int:
        calls.append(sys.argv[1])
        assert cli.VulnServiceConnectionError is wrapper._UnmatchedVulnServiceConnectionError
        if sys.argv[1] == "outer":
            outer_argv = sys.argv
            outer_binding = cli.VulnServiceConnectionError
            assert (
                wrapper.run_pip_audit(
                    ["inner"], max_attempts=1, stdout=io.StringIO(), stderr=io.StringIO()
                )
                == 0
            )
            assert cli.VulnServiceConnectionError is outer_binding
            assert sys.argv is outer_argv
            assert sys.argv == ["pip-audit", "outer"]
        else:
            assert sys.argv == ["pip-audit", "inner"]
        return 0

    monkeypatch.setattr(cli, "audit", audit)
    assert (
        wrapper.run_pip_audit(["outer"], max_attempts=1, stdout=io.StringIO(), stderr=io.StringIO())
        == 0
    )
    assert calls == ["outer", "inner"]
    assert cli.VulnServiceConnectionError is PreviousConnectionError
    assert sys.argv is original_argv


def test_missing_cli_exception_alias_refuses_before_audit_and_is_not_recreated(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    wrapper = _load_script("run_pip_audit")
    cli, _ = _installed_api()
    installed_version = importlib.metadata.version("pip-audit")
    original_argv = sys.argv
    calls: list[str] = []

    def audit() -> int:
        calls.append("audit")
        return 0

    monkeypatch.setattr(cli, "audit", audit)
    monkeypatch.delattr(cli, "VulnServiceConnectionError")
    with pytest.raises(RuntimeError, match=r"pip.audit") as refused:
        wrapper.run_pip_audit([], max_attempts=1, stdout=io.StringIO(), stderr=io.StringIO())

    assert installed_version in str(refused.value)
    assert "VulnServiceConnectionError" in str(refused.value)
    assert "_cli" in str(refused.value)
    assert calls == []
    assert not hasattr(cli, "VulnServiceConnectionError")
    assert sys.argv is original_argv


def test_main_reports_missing_cli_alias_without_recreating_it_or_auditing(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    wrapper = _load_script("run_pip_audit")
    cli, _ = _installed_api()
    original_argv = sys.argv
    calls: list[str] = []

    def audit() -> int:
        calls.append("audit")
        return 0

    monkeypatch.setattr(cli, "audit", audit)
    monkeypatch.delattr(cli, "VulnServiceConnectionError")
    try:
        exit_code = wrapper.main(["--synthetic-missing-alias"])
    except RuntimeError as error:
        raise AssertionError("The CLI must report adapter drift rather than raise it") from error
    assert exit_code == 1
    output = capsys.readouterr()
    assert output.out == ""
    assert output.err.startswith("[pip-audit] COMPATIBILITY ERROR: ")
    assert importlib.metadata.version("pip-audit") in output.err
    assert "VulnServiceConnectionError" in output.err
    assert "Traceback" not in output.err
    assert calls == []
    assert not hasattr(cli, "VulnServiceConnectionError")
    assert sys.argv is original_argv


def test_main_missing_package_names_installation_and_sync_before_private_imports(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    wrapper = _load_script("run_pip_audit")
    original_spec, original_import = importlib.util.find_spec, importlib.import_module
    original_argv = sys.argv
    probes: list[str] = []
    private_imports: list[str] = []

    def absent_package(name: str, package: str | None = None) -> ModuleSpec | None:
        if name == "pip_audit":
            probes.append(name)
            return None
        return original_spec(name, package)

    def unavailable_import(name: str, package: str | None = None) -> ModuleType:
        if name.startswith("pip_audit"):
            private_imports.append(name)
            raise ModuleNotFoundError("Controlled missing package", name="pip_audit")
        return original_import(name, package)

    monkeypatch.setattr(importlib.util, "find_spec", absent_package)
    monkeypatch.setattr(importlib, "import_module", unavailable_import)
    assert wrapper.main(["--synthetic-missing-package"]) != 0
    output = capsys.readouterr()
    assert output.out == ""
    assert "pip-audit is not installed in this environment" in output.err
    assert "just sync" in output.err
    assert "incompatible adapter API" not in output.err
    assert "COMPATIBILITY ERROR" not in output.err
    assert "Traceback" not in output.err
    assert probes == ["pip_audit"]
    assert private_imports == []
    assert sys.argv is original_argv


@pytest.mark.parametrize("invalid_alias", [None, 3, object()])
def test_nonexception_cli_alias_is_a_versioned_compatibility_refusal(
    monkeypatch: pytest.MonkeyPatch, invalid_alias: object
) -> None:
    wrapper = _load_script("run_pip_audit")
    cli, _ = _installed_api()
    installed_version = importlib.metadata.version("pip-audit")
    calls: list[str] = []

    def audit() -> int:
        calls.append("audit")
        return 0

    monkeypatch.setattr(cli, "audit", audit)
    monkeypatch.setattr(cli, "VulnServiceConnectionError", invalid_alias)
    with pytest.raises(RuntimeError, match=r"pip.audit") as refused:
        wrapper.run_pip_audit([], max_attempts=1, stdout=io.StringIO(), stderr=io.StringIO())

    assert installed_version in str(refused.value)
    assert "VulnServiceConnectionError" in str(refused.value)
    assert calls == []
    assert cli.VulnServiceConnectionError is invalid_alias


def test_missing_interface_module_has_a_versioned_preflight_refusal(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    wrapper = _load_script("run_pip_audit")
    cli, _ = _installed_api()
    installed_version = importlib.metadata.version("pip-audit")
    original_import = importlib.import_module
    original_argv = sys.argv
    calls: list[str] = []
    monkeypatch.setattr(cli, "VulnServiceConnectionError", PreviousConnectionError)

    def audit() -> int:
        calls.append("audit")
        return 0

    def missing_interface(name: str, package: str | None = None) -> ModuleType:
        if name == "pip_audit._service.interface":
            raise ModuleNotFoundError("Synthetic interface relocation", name=name)
        return original_import(name, package)

    monkeypatch.setattr(cli, "audit", audit)
    monkeypatch.setattr(importlib, "import_module", missing_interface)
    with pytest.raises(RuntimeError, match=r"pip.audit") as refused:
        wrapper.run_pip_audit([], max_attempts=1, stdout=io.StringIO(), stderr=io.StringIO())

    assert installed_version in str(refused.value)
    assert "pip_audit._service.interface" in str(refused.value)
    assert calls == []
    assert cli.VulnServiceConnectionError is PreviousConnectionError
    assert sys.argv is original_argv


def test_locked_dependency_exposes_the_private_names_required_by_the_adapter() -> None:
    lock = tomllib.loads((REPO_ROOT / "uv.lock").read_text(encoding="utf-8"))
    locked = [package["version"] for package in lock["package"] if package["name"] == "pip-audit"]
    assert len(locked) == 1
    assert importlib.metadata.version("pip-audit") == locked[0]
    cli, interface = _installed_api()
    assert callable(cli.audit)
    assert isinstance(interface.ConnectionError, type)
    assert issubclass(interface.ConnectionError, Exception)
    assert isinstance(interface.ServiceError, type)
    assert issubclass(interface.ServiceError, Exception)
    assert cli.VulnServiceConnectionError is interface.ConnectionError
