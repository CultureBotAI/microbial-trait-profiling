"""Madin placement refuses checkout identity aliases before writes or HTTP (#1178)."""

from __future__ import annotations

import errno
import os
from pathlib import Path

import httpx
import pytest

from mtp import madin_source
from mtp.madin_source import MadinPinError, download_madin_table
from test_madin_source import PAYLOAD, serve
from test_madin_source import pinned_payload as pinned_payload


def _invoke(destination: Path, interface: str) -> None:
    if interface == "cli":
        madin_source.main(["--dest-dir", str(destination)])
    else:
        download_madin_table(destination)


def _alias_pair(root: Path, spelling: str) -> tuple[Path, Path]:
    checkout = root / "Checkout"
    checkout.mkdir()
    alias = {
        "lexical": root / "unrelated-tree",
        "case": root / "cHECKOUT",
        "firmlink": root / "volume-view" / "Checkout",
    }[spelling]
    alias.mkdir(parents=True, exist_ok=True)
    return checkout, alias


def _bind_directory_identity(monkeypatch: pytest.MonkeyPatch, alias: Path, checkout: Path) -> None:
    """Model a directory alias that resolve() does not rewrite, on any host FS."""
    original_stat = Path.stat
    checkout_stat = checkout.stat()

    def aliased_stat(path: Path, *, follow_symlinks: bool = True) -> os.stat_result:
        if path == alias and follow_symlinks:
            return checkout_stat
        return original_stat(path, follow_symlinks=follow_symlinks)

    monkeypatch.setattr(Path, "stat", aliased_stat)


@pytest.mark.usefixtures("pinned_payload")
@pytest.mark.parametrize("interface", ["api", "cli"])
@pytest.mark.parametrize("spelling", ["lexical", "case", "firmlink"])
def test_madin_refuses_checkout_alias_identity_with_missing_output_suffix(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, interface: str, spelling: str
) -> None:
    checkout, alias = _alias_pair(tmp_path, spelling)
    monkeypatch.setattr(madin_source, "REPO_ROOT", checkout)
    destination = alias / "new" / "nested" / "madin"
    before = set(tmp_path.rglob("*"))
    requests: list[httpx.Request] = []
    serve(monkeypatch, PAYLOAD, requests)
    expected_error = MadinPinError if interface == "api" else SystemExit

    with monkeypatch.context() as identity_patch:
        _bind_directory_identity(identity_patch, alias, checkout)
        with pytest.raises(expected_error, match=r"(?i)repository|working.tree"):
            _invoke(destination, interface)

    assert requests == []
    assert set(tmp_path.rglob("*")) == before
    assert not destination.exists()


@pytest.mark.usefixtures("pinned_payload")
@pytest.mark.parametrize("interface", ["api", "cli"])
def test_madin_existing_target_alias_checks_directory_ancestry_before_reuse(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, interface: str
) -> None:
    checkout, alias = _alias_pair(tmp_path, "lexical")
    monkeypatch.setattr(madin_source, "REPO_ROOT", checkout)
    source = alias / "retained-source.csv"
    source.write_bytes(PAYLOAD)
    # The regular file has its own identity; containment is in its directory ancestry.
    assert (source.stat().st_dev, source.stat().st_ino) != (
        checkout.stat().st_dev,
        checkout.stat().st_ino,
    )
    destination = tmp_path / "external-data"
    destination.mkdir()
    target = destination / madin_source.MADIN_PIN.filename
    target.symlink_to(source)
    before = set(tmp_path.rglob("*"))
    requests: list[httpx.Request] = []
    serve(monkeypatch, PAYLOAD, requests)
    expected_error = MadinPinError if interface == "api" else SystemExit

    with monkeypatch.context() as identity_patch:
        _bind_directory_identity(identity_patch, alias, checkout)
        with pytest.raises(expected_error, match=r"(?i)repository|working.tree"):
            _invoke(destination, interface)

    assert requests == []
    assert set(tmp_path.rglob("*")) == before
    assert source.read_bytes() == PAYLOAD
    assert target.is_symlink()
    assert not (destination / madin_source.PROVENANCE_FILE).exists()


@pytest.mark.usefixtures("pinned_payload")
@pytest.mark.parametrize("interface", ["api", "cli"])
@pytest.mark.parametrize("failed_path", ["checkout", "candidate_ancestor"])
def test_madin_stat_permission_failure_refuses_before_writes_or_http(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    interface: str,
    failed_path: str,
) -> None:
    checkout, alias = _alias_pair(tmp_path, "lexical")
    monkeypatch.setattr(madin_source, "REPO_ROOT", checkout)
    destination = alias / "new" / "madin"
    denied = checkout if failed_path == "checkout" else alias
    before = set(tmp_path.rglob("*"))
    requests: list[httpx.Request] = []
    serve(monkeypatch, PAYLOAD, requests)
    original_stat = Path.stat
    expected_error = MadinPinError if interface == "api" else SystemExit

    def denied_stat(path: Path, *, follow_symlinks: bool = True) -> os.stat_result:
        if path == denied:
            raise PermissionError(errno.EACCES, "Permission denied", str(path))
        return original_stat(path, follow_symlinks=follow_symlinks)

    with monkeypatch.context() as identity_patch:
        identity_patch.setattr(Path, "stat", denied_stat)
        if interface == "cli" and failed_path == "candidate_ancestor":
            with pytest.raises(SystemExit) as error:
                _invoke(destination, interface)
            assert error.value.code == 2
            stderr = capsys.readouterr().err
            assert "Cannot inspect output path" in stderr
            assert "Permission denied" in stderr
        else:
            with pytest.raises(
                expected_error, match=r"(?i)inspect|identity|permission|stat|filesystem"
            ):
                _invoke(destination, interface)

    assert requests == []
    assert set(tmp_path.rglob("*")) == before
    assert not destination.exists()


@pytest.mark.usefixtures("pinned_payload")
@pytest.mark.parametrize("interface", ["api", "cli"])
def test_madin_safe_external_reuse_retains_pinned_bytes_without_http(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, interface: str
) -> None:
    checkout = tmp_path / "checkout"
    checkout.mkdir()
    monkeypatch.setattr(madin_source, "REPO_ROOT", checkout)
    destination = tmp_path / "external-data"
    destination.mkdir()
    target = destination / madin_source.MADIN_PIN.filename
    target.write_bytes(PAYLOAD)
    before = target.stat()
    requests: list[httpx.Request] = []
    serve(monkeypatch, PAYLOAD, requests)

    _invoke(destination, interface)

    assert requests == []
    assert target.read_bytes() == PAYLOAD
    assert (target.stat().st_dev, target.stat().st_ino) == (before.st_dev, before.st_ino)
    assert (destination / madin_source.PROVENANCE_FILE).is_file()
    assert list(checkout.iterdir()) == []
