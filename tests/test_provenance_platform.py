"""Portable account lookup does not require POSIX or a configured home (#1777)."""

from __future__ import annotations

import builtins
import importlib.util
import os
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

from test_provenance_paths import api as api


def test_module_import_and_writer_work_without_pwd(
    tmp_path: Path,
    api: ModuleType,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    original_import = builtins.__import__

    def without_pwd(name: str, *args: Any, **kwargs: Any) -> ModuleType:
        if name == "pwd":
            raise ImportError("No password database on this platform")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", without_pwd)
    spec = importlib.util.spec_from_file_location("provenance_without_pwd", api.__file__)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    monkeypatch.setitem(sys.modules, spec.name, module)
    try:
        spec.loader.exec_module(module)
        paths = module.PortablePaths.create(tmp_path / "record.json", [tmp_path / "source.tsv"])
        assert paths.encode(tmp_path / "source.tsv") == "source.tsv"
    except ImportError as error:
        raise AssertionError("Portable writers must work without the pwd module") from error


@pytest.mark.parametrize("operation", ["create", "encode"])
def test_unknown_home_is_an_actionable_privacy_refusal(
    tmp_path: Path,
    api: ModuleType,
    monkeypatch: pytest.MonkeyPatch,
    operation: str,
) -> None:
    def unavailable_home() -> Path:
        raise RuntimeError("Could not determine home directory.")

    monkeypatch.setattr(Path, "home", unavailable_home)

    def exercise() -> None:
        if operation == "create":
            api.PortablePaths.create(tmp_path / "record.json", [tmp_path / "source.tsv"])
        else:
            api.PortablePaths(tmp_path / "record.json", tmp_path).encode(tmp_path / "source.tsv")

    try:
        with pytest.raises(
            ValueError, match=r"^Cannot determine the home directory for provenance privacy checks"
        ):
            exercise()
    except RuntimeError as error:
        raise AssertionError("An unconfigured home must produce a controlled ValueError") from error


def test_a_missing_password_entry_does_not_block_a_configured_writer(
    tmp_path: Path,
    api: ModuleType,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    password_module = ModuleType("pwd")

    def missing_entry(uid: int) -> None:
        raise KeyError(uid)

    monkeypatch.setattr(password_module, "getpwuid", missing_entry, raising=False)
    monkeypatch.setitem(sys.modules, "pwd", password_module)
    monkeypatch.setattr(os, "getuid", lambda: 987, raising=False)
    try:
        assert api._password_entry() is None
        paths = api.PortablePaths.create(tmp_path / "record.json", [tmp_path / "input.tsv"])
        assert paths.encode(tmp_path / "input.tsv") == "input.tsv"
    except KeyError as error:
        raise AssertionError("A numeric UID without a password entry must be supported") from error
