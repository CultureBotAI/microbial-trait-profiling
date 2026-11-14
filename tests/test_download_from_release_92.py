"""Tests for downloading BugSigDB raw CSVs from a GitHub release tarball (#92)."""

from __future__ import annotations

import hashlib
import io
import json
import os
import shutil
import tarfile
from pathlib import Path

import pytest

from mtp import download as dl
from mtp.download import (
    DEFAULT_RELEASE_REPO,
    EXPORT_FILES,
    FULL_DUMP_NAME,
    PROVENANCE_FILE,
    RELEASE_ARCHIVE_PATTERN,
    RELEASE_EXPECTED_CSVS,
    SOURCE_RELEASE,
    IntegrityError,
    SourceResolutionError,
    download_from_release,
    main,
)
from mtp.paths import REPO_ROOT

SNAPSHOT_HEADER = (
    b"# BugSigDB 2026-03-30_00:00_UTC, License: Creative Commons Attribution 4.0 "
    b"International, URL: https://bugsigdb.org\n"
)


def _write_tar_gz(
    archive_path: Path,
    members: dict[str, bytes],
    *,
    include_dot_dir: bool = True,
    symlink_member: tuple[str, str] | None = None,
) -> Path:
    archive_path.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive_path, mode="w:gz") as tar:
        if include_dot_dir:
            dot_info = tarfile.TarInfo(name=".")
            dot_info.type = tarfile.DIRTYPE
            dot_info.mode = 0o755
            tar.addfile(dot_info)
        for name, payload in members.items():
            info = tarfile.TarInfo(name=name)
            info.size = len(payload)
            info.mode = 0o644
            tar.addfile(info, io.BytesIO(payload))
        if symlink_member is not None:
            link_name, link_target = symlink_member
            link_info = tarfile.TarInfo(name=link_name)
            link_info.type = tarfile.SYMTYPE
            link_info.linkname = link_target
            tar.addfile(link_info)
    return archive_path


def _default_csv_members(*, prefix: str = "./") -> dict[str, bytes]:
    return {
        f"{prefix}full_dump.csv": SNAPSHOT_HEADER + b"BSDB ID,Study\nbsdb:1/1/1,StudyOne\n",
        f"{prefix}studies.csv": b"Study,State\nbsdb:1,Complete\n",
        f"{prefix}experiments.csv": b"Study,Experiment\nbsdb:1,bsdb:1/1\n",
        f"{prefix}signatures.csv": b"Study,Experiment,Signature\nbsdb:1,bsdb:1/1,bsdb:1/1/1\n",
    }


def test_download_from_release_extracts_csvs_and_writes_provenance(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """download_from_release extracts expected CSVs and records release provenance (#92)."""
    calls: list[tuple[str | None, str, str]] = []
    members = _default_csv_members(prefix="./")

    def fake_fetch(dest_dir: Path, *, tag: str | None, repo: str) -> Path:
        calls.append((tag, repo, RELEASE_ARCHIVE_PATTERN))
        return _write_tar_gz(dest_dir / "bsdbng-raw-data-v2026.03.30.tar.gz", members)

    monkeypatch.setattr(dl, "_fetch_release_archive", fake_fetch)

    out_dir = tmp_path / "raw"
    extracted = download_from_release(tag="latest", output_dir=out_dir)

    assert calls == [(None, DEFAULT_RELEASE_REPO, "bsdbng-raw-data-*.tar.gz")]
    assert [p.name for p in extracted] == [*EXPORT_FILES, FULL_DUMP_NAME]
    for name in RELEASE_EXPECTED_CSVS:
        assert (out_dir / name).read_bytes() == members[f"./{name}"]

    prov_path = out_dir / PROVENANCE_FILE
    assert prov_path.is_file()
    prov = json.loads(prov_path.read_text(encoding="utf-8"))
    assert prov["source"] == SOURCE_RELEASE
    assert prov["release_tag"] == "v2026.03.30"
    assert prov["release_repo"] == DEFAULT_RELEASE_REPO
    assert prov["archive_name"] == "bsdbng-raw-data-v2026.03.30.tar.gz"
    assert len(prov["archive_sha256"]) == 64
    assert prov["bugsigdb_release"]["release_tag"] == "v2026.03.30"
    assert prov["bugsigdb_release"].get("bugsigdb_snapshot") == "2026-03-30_00:00_UTC"
    assert prov[FULL_DUMP_NAME]["source"] == SOURCE_RELEASE
    assert prov[FULL_DUMP_NAME].get("bugsigdb_snapshot") == "2026-03-30_00:00_UTC"
    assert prov[FULL_DUMP_NAME]["sha256"] == hashlib.sha256(members["./full_dump.csv"]).hexdigest()


def test_download_from_release_passes_explicit_tag_and_verifies_pinned_match(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Explicit tag is forwarded and matching pinned full_dump.csv is marked verified (#92)."""
    members = _default_csv_members(prefix="")
    full_dump_bytes = members[FULL_DUMP_NAME]
    monkeypatch.setattr(dl, "PINNED_FULL_DUMP_BYTES", len(full_dump_bytes))
    monkeypatch.setattr(
        dl,
        "PINNED_FULL_DUMP_MD5",
        hashlib.md5(full_dump_bytes, usedforsecurity=False).hexdigest(),
    )
    observed: list[tuple[str | None, str]] = []

    def fake_fetch(dest_dir: Path, *, tag: str | None, repo: str) -> Path:
        observed.append((tag, repo))
        return _write_tar_gz(
            dest_dir / "bsdbng-raw-data-v2026.04.01.tar.gz",
            members,
            include_dot_dir=False,
        )

    monkeypatch.setattr(dl, "_fetch_release_archive", fake_fetch)

    out_dir = tmp_path / "raw_pinned"
    assert (
        main(
            [
                "--from-release",
                "v2026.04.01",
                "--repo",
                "example/repo",
                "--output-dir",
                str(out_dir),
            ]
        )
        == 0
    )
    assert observed == [("v2026.04.01", "example/repo")]
    prov = json.loads((out_dir / PROVENANCE_FILE).read_text(encoding="utf-8"))
    assert prov["release_tag"] == "v2026.04.01"
    assert prov["release_repo"] == "example/repo"
    assert prov[FULL_DUMP_NAME]["pinned"] is True
    assert prov[FULL_DUMP_NAME]["integrity_verified"] is True


def test_download_from_release_rejects_corrupt_pinned_size_dump(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A full_dump.csv matching PINNED_FULL_DUMP_BYTES with wrong MD5 raises IntegrityError."""
    members = _default_csv_members(prefix="")
    monkeypatch.setattr(dl, "PINNED_FULL_DUMP_BYTES", len(members[FULL_DUMP_NAME]))
    monkeypatch.setattr(dl, "PINNED_FULL_DUMP_MD5", "0" * 32)

    monkeypatch.setattr(
        dl,
        "_fetch_release_archive",
        lambda dest_dir, **_kw: _write_tar_gz(
            dest_dir / "bsdbng-raw-data-v2026.03.30.tar.gz", members
        ),
    )

    out_dir = tmp_path / "corrupt"
    msg = ""
    try:
        download_from_release(tag="v2026.03.30", output_dir=out_dir)
    except IntegrityError as exc:
        msg = str(exc)
    assert "md5" in msg
    assert not (out_dir / FULL_DUMP_NAME).exists()


@pytest.mark.parametrize(
    "bad_member",
    [
        "../evil.csv",
        "/etc/passwd",
        "subdir/full_dump.csv",
        "unexpected.csv",
    ],
)
def test_download_from_release_rejects_unsafe_or_unexpected_tar_members(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, bad_member: str
) -> None:
    """Path traversal, subdirectories, and unexpected files in release archive fail closed."""
    members = _default_csv_members(prefix="")
    members[bad_member] = b"bad"

    monkeypatch.setattr(
        dl,
        "_fetch_release_archive",
        lambda dest_dir, **_kw: _write_tar_gz(
            dest_dir / "bsdbng-raw-data-v2026.03.30.tar.gz", members
        ),
    )

    out_dir = tmp_path / "unsafe"
    msg = ""
    try:
        download_from_release(tag="v2026.03.30", output_dir=out_dir)
    except IntegrityError as exc:
        msg = str(exc)
    assert msg != ""
    assert not (out_dir / FULL_DUMP_NAME).exists()


def test_download_from_release_rejects_symlink_tar_member(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Symlink members in release archive are rejected before extraction."""
    members = _default_csv_members(prefix="")

    monkeypatch.setattr(
        dl,
        "_fetch_release_archive",
        lambda dest_dir, **_kw: _write_tar_gz(
            dest_dir / "bsdbng-raw-data-v2026.03.30.tar.gz",
            members,
            symlink_member=("studies.csv", "/etc/passwd"),
        ),
    )

    out_dir = tmp_path / "symlink"
    msg = ""
    try:
        download_from_release(tag="v2026.03.30", output_dir=out_dir)
    except IntegrityError as exc:
        msg = str(exc)
    assert "non-regular" in msg or "Duplicate" in msg
    assert not (out_dir / FULL_DUMP_NAME).exists()


def test_download_from_release_rejects_missing_full_dump_and_bad_args(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Missing full_dump.csv or invalid repo/tag strings fail closed."""
    members = {"studies.csv": b"Study,State\nbsdb:1,Complete\n"}
    monkeypatch.setattr(
        dl,
        "_fetch_release_archive",
        lambda dest_dir, **_kw: _write_tar_gz(
            dest_dir / "bsdbng-raw-data-v2026.03.30.tar.gz", members
        ),
    )
    msg_missing = ""
    try:
        download_from_release(tag="v2026.03.30", output_dir=tmp_path / "missing_dump")
    except IntegrityError as exc:
        msg_missing = str(exc)
    assert FULL_DUMP_NAME in msg_missing

    valid_members = _default_csv_members(prefix="")
    monkeypatch.setattr(
        dl,
        "_fetch_release_archive",
        lambda dest_dir, **_kw: _write_tar_gz(
            dest_dir / "bsdbng-raw-data-v2026.03.30.tar.gz", valid_members
        ),
    )

    for bad_tag in ("--flag", "v1/../../bad", "v1 2"):
        msg_tag = ""
        try:
            download_from_release(tag=bad_tag, output_dir=tmp_path / "bad_tag")
        except ValueError as exc:
            msg_tag = str(exc)
        assert msg_tag != ""

    for bad_repo in ("", "invalid-repo", "-flag/repo", "owner/repo/extra"):
        msg_repo = ""
        try:
            download_from_release(tag="latest", output_dir=tmp_path / "bad_repo", repo=bad_repo)
        except ValueError as exc:
            msg_repo = str(exc)
        assert msg_repo != ""


def test_fetch_release_archive_handles_missing_gh_and_nonzero_exit(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """_fetch_release_archive fails with SourceResolutionError when gh is missing or fails."""
    monkeypatch.setattr(shutil, "which", lambda _name: None)
    msg_no_gh = ""
    try:
        dl._fetch_release_archive(tmp_path, tag=None, repo=DEFAULT_RELEASE_REPO)
    except SourceResolutionError as exc:
        msg_no_gh = str(exc)
    assert "gh" in msg_no_gh

    monkeypatch.setattr(shutil, "which", lambda _name: "/usr/bin/gh")
    monkeypatch.setattr(os, "posix_spawn", lambda *_a, **_kw: 999)
    monkeypatch.setattr(os, "waitpid", lambda _pid, _opts: (999, 256))
    monkeypatch.setattr(os, "waitstatus_to_exitcode", lambda _status: 1)
    msg_exit = ""
    try:
        dl._fetch_release_archive(tmp_path, tag="v2026.03.30", repo=DEFAULT_RELEASE_REPO)
    except SourceResolutionError as exc:
        msg_exit = str(exc)
    assert "exit code 1" in msg_exit


def test_justfile_declares_download_from_release_recipe() -> None:
    """Justfile defines download-from-release tag='latest' calling python -m mtp.download (#92)."""
    justfile_text = (REPO_ROOT / "Justfile").read_text(encoding="utf-8")
    assert 'download-from-release tag="latest":' in justfile_text
    assert 'uv run python -m mtp.download --from-release "{{tag}}"' in justfile_text
