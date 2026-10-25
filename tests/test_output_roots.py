"""Run and archive outputs a record cites must not live in temporary storage (#1049)."""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import pytest

from mtp.workflow.output_roots import temporary_root

MACOS = pytest.mark.skipif(sys.platform != "darwin", reason="macOS links /tmp to /private/tmp")


@pytest.fixture
def fixed_process_temporary(monkeypatch: pytest.MonkeyPatch) -> None:
    """Pin $TMPDIR and Python's temporary directory to roots no real path uses."""
    monkeypatch.setenv("TMPDIR", "/mtp-1049-environment-tmp")
    monkeypatch.setattr(tempfile, "tempdir", "/mtp-1049-python-tmp")


@pytest.mark.usefixtures("fixed_process_temporary")
@pytest.mark.parametrize(
    "output",
    [
        "/private/tmp/core",
        "/var/folders/ab/T/core",
        "/private/var/folders/ab/T/core",
        "/private/tmp",
    ],
)
def test_fixed_temporary_roots_are_refused(output: str) -> None:
    assert temporary_root(Path(output)) is not None


@MACOS
@pytest.mark.usefixtures("fixed_process_temporary")
def test_macos_tmp_is_refused_through_its_link() -> None:
    """/tmp is not a fixed root (#1051); on macOS it resolves to /private/tmp."""
    assert temporary_root(Path("/tmp/core")) == Path("/private/tmp")


@MACOS
def test_a_symlink_to_tmp_is_refused_on_macos(tmp_path: Path) -> None:
    if tmp_path.resolve().is_relative_to("/private/tmp"):
        # Without $TMPDIR, pytest's scratch is under /private/tmp itself (#1086)
        pytest.skip("the link's own location would explain the refusal")
    link = tmp_path / "link"
    link.symlink_to("/tmp", target_is_directory=True)
    # With only the fixed /private/tmp root, the link's own location cannot be the reason
    assert temporary_root(tmp_path / "core", [Path("/private/tmp")]) is None
    assert temporary_root(link / "core", [Path("/private/tmp")]) == Path("/private/tmp")


def test_the_process_temporary_directory_is_refused_on_any_system() -> None:
    assert temporary_root(Path(tempfile.gettempdir()) / "core") is not None


@pytest.mark.usefixtures("fixed_process_temporary")
@pytest.mark.parametrize(
    "output", ["/mtp-1049-environment-tmp/core", "/mtp-1049-python-tmp/nested/core"]
)
def test_the_process_temporary_directories_are_refused(output: str) -> None:
    """$TMPDIR and tempfile.gettempdir() count even when they are not a fixed root."""
    assert temporary_root(Path(output)) is not None


@pytest.mark.usefixtures("fixed_process_temporary")
@pytest.mark.parametrize(
    "output", ["/srv/mtp-runs/core", "/tmpfiles/core", "/var/foldersx/core", "/private/var/db"]
)
def test_durable_paths_are_accepted(output: str) -> None:
    """A shared prefix that is not a parent directory is not temporary storage."""
    assert temporary_root(Path(output)) is None


def test_a_symlink_into_temporary_storage_is_refused(tmp_path: Path) -> None:
    volatile, durable = tmp_path / "volatile", tmp_path / "durable"
    volatile.mkdir()
    durable.mkdir()
    (durable / "link").symlink_to(volatile, target_is_directory=True)

    assert temporary_root(durable / "link" / "core", [volatile]) == volatile.resolve()
    assert temporary_root(durable / "core", [volatile]) is None


def test_a_symlinked_temporary_root_still_holds_its_target(tmp_path: Path) -> None:
    """macOS's /tmp and /var are links; a root is compared in both spellings."""
    volatile, alias = tmp_path / "volatile", tmp_path / "alias"
    volatile.mkdir()
    alias.symlink_to(volatile, target_is_directory=True)

    assert temporary_root(volatile / "core", [alias]) is not None
    assert temporary_root(alias / "core", [volatile]) is not None


def test_a_path_written_under_temporary_storage_is_refused_wherever_it_leads(
    tmp_path: Path,
) -> None:
    """The written spelling is checked too, a conservative over-refusal.

    On macOS /tmp is not a written root (#1051), so /tmp/<link to durable> is accepted.
    """
    volatile, durable = tmp_path / "volatile", tmp_path / "durable"
    volatile.mkdir()
    durable.mkdir()
    (volatile / "link").symlink_to(durable, target_is_directory=True)

    assert temporary_root(volatile / "link" / "core", [volatile]) is not None


def test_a_relative_path_is_judged_where_it_resolves(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    volatile = tmp_path / "volatile"
    volatile.mkdir()
    monkeypatch.chdir(volatile)

    assert temporary_root(Path("core"), [volatile]) is not None
    assert temporary_root(Path("../durable/core"), [volatile]) is None
