"""`just download-madin` fetches only the pinned Madin et al. table (#1048)."""

from __future__ import annotations

import hashlib
import json
import os
import re
from collections.abc import Iterator
from dataclasses import replace
from pathlib import Path
from typing import Any

import httpx
import pytest

from mtp import madin_source
from mtp.madin_source import MADIN_PIN, MadinPinError, download_madin_table, verify_madin_table
from mtp.workflow import output_roots

PAYLOAD = b'"tax_id","species_tax_id","data_source"\n11,11,"gold"\n'


@pytest.fixture
def pinned_payload(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """PAYLOAD stands in for the pinned bytes; repository, commit and path stay real."""
    monkeypatch.setattr(
        madin_source,
        "MADIN_PIN",
        replace(MADIN_PIN, sha256=hashlib.sha256(PAYLOAD).hexdigest(), size_bytes=len(PAYLOAD)),
    )
    monkeypatch.setattr(
        output_roots, "system_temporary_directories", lambda: [tmp_path / "volatile"]
    )


def serve(
    monkeypatch: pytest.MonkeyPatch, body: bytes, requests: list[httpx.Request] | None = None
) -> None:
    """Route the download's httpx.Client through a MockTransport that returns *body*."""

    def handler(request: httpx.Request) -> httpx.Response:
        if requests is None:
            pytest.fail(f"unexpected request to {request.url}")
        requests.append(request)
        return httpx.Response(200, content=body)

    real_client = httpx.Client

    def client(**kwargs: Any) -> httpx.Client:
        return real_client(transport=httpx.MockTransport(handler), **kwargs)

    monkeypatch.setattr(httpx, "Client", client)


def test_pin_is_fixed_to_a_commit() -> None:
    assert re.fullmatch(r"[0-9a-f]{40}", MADIN_PIN.commit)
    # A typo in the pinned hash must fail here, not only when a real download is refused.
    assert MADIN_PIN.sha256 == "6100fc11679cd9facae65cad37a2969978bc64c3e57984d896f8147f78db6440"
    assert MADIN_PIN.url == (
        "https://raw.githubusercontent.com/bacteria-archaea-traits/bacteria-archaea-traits/"
        "396dcd4d98ca20b613ba280f5514611651b6fb91/output/condensed_traits_NCBI.csv"
    )
    assert MADIN_PIN.size_bytes == 46362149


@pytest.mark.usefixtures("pinned_payload")
def test_download_keeps_pinned_bytes_and_records_them(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    requests: list[httpx.Request] = []
    serve(monkeypatch, PAYLOAD, requests)
    (tmp_path / madin_source.PROVENANCE_FILE).write_text('{"full_dump.csv": {"pinned": true}}')

    path = download_madin_table(tmp_path)

    assert path == tmp_path / "condensed_traits_NCBI.csv"
    assert path.read_bytes() == PAYLOAD
    [request] = requests
    assert str(request.url) == madin_source.MADIN_PIN.url
    assert not {"authorization", "cookie", "from"} & set(request.headers.keys())
    assert request.headers["user-agent"].startswith("python-httpx/")
    provenance = json.loads((tmp_path / madin_source.PROVENANCE_FILE).read_text())
    assert provenance["full_dump.csv"] == {"pinned": True}
    record = provenance["condensed_traits_NCBI.csv"]
    assert record["commit"] == MADIN_PIN.commit
    assert record["sha256"] == hashlib.sha256(PAYLOAD).hexdigest()
    assert record["integrity_verified"] is True


@pytest.mark.usefixtures("pinned_payload")
@pytest.mark.parametrize(
    ("body", "message"),
    [(PAYLOAD + b"12,12,patric\n", "bytes, expected"), (PAYLOAD.upper(), "has sha256")],
    ids=["size", "sha256"],
)
def test_download_refuses_bytes_that_differ_from_the_pin(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, body: bytes, message: str
) -> None:
    serve(monkeypatch, body, [])

    with pytest.raises(MadinPinError, match=message):
        download_madin_table(tmp_path)

    assert list(tmp_path.iterdir()) == []


@pytest.mark.usefixtures("pinned_payload")
def test_existing_copy_is_checked_and_never_overwritten(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    serve(monkeypatch, PAYLOAD)
    target = tmp_path / "condensed_traits_NCBI.csv"
    target.write_bytes(PAYLOAD)
    assert download_madin_table(tmp_path) == target

    target.write_bytes(PAYLOAD.upper())
    with pytest.raises(MadinPinError, match="remove it to fetch the pinned file again"):
        download_madin_table(tmp_path)
    assert target.read_bytes() == PAYLOAD.upper()


@pytest.mark.usefixtures("pinned_payload")
def test_cli_exits_on_a_mismatch(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    serve(monkeypatch, PAYLOAD.upper(), [])

    with pytest.raises(SystemExit, match="has sha256"):
        madin_source.main(["--dest-dir", str(tmp_path / "madin")])


@pytest.mark.usefixtures("pinned_payload")
def test_verify_returns_the_pin_for_matching_bytes(tmp_path: Path) -> None:
    path = tmp_path / "table.csv"
    path.write_bytes(PAYLOAD)
    assert verify_madin_table(path) == madin_source.MADIN_PIN.descriptor()


@pytest.mark.usefixtures("pinned_payload")
@pytest.mark.parametrize(
    "location",
    [
        "relative_repo",
        "absolute_repo",
        "repo_alias",
        "repo_outgoing_alias",
        "temporary",
        "temporary_alias",
        "target_repo_alias",
        "target_temporary_alias",
    ],
)
def test_unsafe_destination_is_refused_before_writes_or_network(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, location: str
) -> None:
    checkout = tmp_path / "checkout"
    volatile = tmp_path / "volatile"
    checkout.mkdir()
    volatile.mkdir()
    monkeypatch.chdir(checkout)
    monkeypatch.setattr(madin_source, "REPO_ROOT", checkout, raising=False)
    destination = tmp_path / "data" / "new" / "madin"
    if location == "relative_repo":
        destination = Path("new") / "madin"
    elif location == "absolute_repo":
        destination = checkout / "new" / "madin"
    elif location == "repo_outgoing_alias":
        outside = tmp_path / "outside"
        outside.mkdir()
        alias = checkout / "outward"
        alias.symlink_to(outside, target_is_directory=True)
        destination = alias / "new" / "madin"
    elif location in {"repo_alias", "temporary_alias"}:
        alias = tmp_path / "alias"
        alias.symlink_to(
            checkout if location == "repo_alias" else volatile, target_is_directory=True
        )
        destination = alias / "new" / "madin"
    elif location == "temporary":
        destination = volatile / "new" / "madin"
    else:
        destination.mkdir(parents=True)
        source = (checkout if location == "target_repo_alias" else volatile) / "source.csv"
        source.write_bytes(PAYLOAD)
        (destination / MADIN_PIN.filename).symlink_to(source)
    before = set(tmp_path.rglob("*"))
    requests: list[httpx.Request] = []
    serve(monkeypatch, PAYLOAD, requests)

    with pytest.raises(MadinPinError, match=r"(?i)repository|working.tree|temporary|durable"):
        download_madin_table(destination)

    assert requests == []
    assert set(tmp_path.rglob("*")) == before


@pytest.mark.usefixtures("pinned_payload")
def test_relative_path_outside_checkout_is_allowed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    checkout = tmp_path / "checkout"
    checkout.mkdir()
    monkeypatch.chdir(checkout)
    monkeypatch.setattr(madin_source, "REPO_ROOT", checkout)
    requests: list[httpx.Request] = []
    serve(monkeypatch, PAYLOAD, requests)

    target = download_madin_table(Path("..") / "data" / "madin")

    assert target == tmp_path / "data" / "madin" / MADIN_PIN.filename
    assert target.read_bytes() == PAYLOAD
    assert len(requests) == 1
    assert list(checkout.iterdir()) == []


@pytest.mark.usefixtures("pinned_payload")
def test_cli_uses_durable_guard_before_creating_ancestors(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    requests: list[httpx.Request] = []
    serve(monkeypatch, PAYLOAD, requests)
    destination = tmp_path / "volatile" / "new" / "madin"
    with pytest.raises(SystemExit):
        madin_source.main(["--dest-dir", str(destination)])
    assert requests == []
    assert not destination.parent.exists()


@pytest.mark.usefixtures("pinned_payload")
def test_expanded_destination_is_used_for_the_actual_download(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    home = tmp_path / "user-home"
    home.mkdir()
    work = tmp_path / "work"
    work.mkdir()
    original_expanduser = os.path.expanduser

    def expanduser(path: str) -> str:
        return str(home) if path == "~" else original_expanduser(path)

    monkeypatch.setattr(os.path, "expanduser", expanduser)
    monkeypatch.chdir(work)
    requests: list[httpx.Request] = []
    serve(monkeypatch, PAYLOAD, requests)

    path = download_madin_table(Path("~") / "madin")

    assert path == home / "madin" / MADIN_PIN.filename
    assert path.read_bytes() == PAYLOAD
    assert not (work / "~").exists()
    assert len(requests) == 1


class OverflowStream(httpx.SyncByteStream):
    """Observe partial bytes and make reading beyond the first excess chunk fatal."""

    def __init__(self, destination: Path, *, forbid_further_read: bool) -> None:
        self.destination = destination
        self.forbid_further_read = forbid_further_read
        self.delivered = 0
        self.closed = False
        self.partial_sizes: list[int] = []

    def __iter__(self) -> Iterator[bytes]:
        for chunk in (PAYLOAD, b"x"):
            self.delivered += 1
            yield chunk
        if self.forbid_further_read:
            pytest.fail("Download consumed the response after exceeding the pinned size")

    def close(self) -> None:
        self.closed = True
        self.partial_sizes = [path.stat().st_size for path in self.destination.glob("*.part")]


def serve_stream(monkeypatch: pytest.MonkeyPatch, stream: OverflowStream) -> None:
    def handler(_request: httpx.Request) -> httpx.Response:
        return httpx.Response(200, stream=stream)

    real_client = httpx.Client

    def client(**kwargs: Any) -> httpx.Client:
        return real_client(transport=httpx.MockTransport(handler), **kwargs)

    monkeypatch.setattr(httpx, "Client", client)


@pytest.mark.usefixtures("pinned_payload")
@pytest.mark.parametrize(
    "forbid_further_read", [False, True], ids=["no-overflow-write", "abort-early"]
)
def test_overflow_is_stopped_before_write_or_further_read(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, forbid_further_read: bool
) -> None:
    stream = OverflowStream(tmp_path, forbid_further_read=forbid_further_read)
    serve_stream(monkeypatch, stream)
    deleted_partial_sizes: list[int] = []
    original_unlink = os.unlink

    def observe_unlink(path: str | Path, *, dir_fd: int | None = None) -> None:
        candidate = Path(path)
        if candidate.suffix == ".part" and candidate.is_file():
            deleted_partial_sizes.append(candidate.stat().st_size)
        original_unlink(path, dir_fd=dir_fd)

    monkeypatch.setattr(os, "unlink", observe_unlink)

    with pytest.raises(MadinPinError, match=r"(?i)bytes|size|exceed"):
        download_madin_table(tmp_path)

    assert stream.delivered == 2
    assert stream.closed
    assert all(size <= len(PAYLOAD) for size in stream.partial_sizes)
    assert deleted_partial_sizes
    assert all(size <= len(PAYLOAD) for size in deleted_partial_sizes)
    assert list(tmp_path.iterdir()) == []


@pytest.mark.usefixtures("pinned_payload")
@pytest.mark.parametrize("existing_record", [False, True])
def test_verified_reuse_writes_provenance_without_download(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, existing_record: bool
) -> None:
    target = tmp_path / MADIN_PIN.filename
    target.write_bytes(PAYLOAD)
    record_path = tmp_path / madin_source.PROVENANCE_FILE
    if existing_record:
        record_path.write_text(
            json.dumps(
                {
                    "another.csv": {"note": "preserve this public record"},
                    MADIN_PIN.filename: {
                        **madin_source.MADIN_PIN.descriptor(),
                        "downloaded_at": "2026-01-01T00:00:00+00:00",
                    },
                }
            )
        )
    serve(monkeypatch, PAYLOAD)

    assert download_madin_table(tmp_path) == target

    records = json.loads(record_path.read_bytes())
    record = records[MADIN_PIN.filename]
    for key, value in madin_source.MADIN_PIN.descriptor().items():
        assert record[key] == value
    assert record["integrity_verified"] is True
    assert record["retrieval"] == "verified_reuse"
    assert record["verified_at"]
    if existing_record:
        assert record["downloaded_at"] == "2026-01-01T00:00:00+00:00"
        assert records["another.csv"] == {"note": "preserve this public record"}
    else:
        assert "downloaded_at" not in record
    assert target.read_bytes() == PAYLOAD


@pytest.mark.usefixtures("pinned_payload")
@pytest.mark.parametrize("existing_table", [False, True])
@pytest.mark.parametrize(
    "record_bytes",
    [
        b"{invalid-json",
        b"[]",
        b'{"other.csv": null}',
        b'{"condensed_traits_NCBI.csv": []}',
        b'{"condensed_traits_NCBI.csv": {"commit": "contradictory"}}',
    ],
    ids=["invalid-json", "root-shape", "entry-shape", "madin-shape", "contradictory-pin"],
)
def test_bad_provenance_is_rejected_before_network_or_publication(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, record_bytes: bytes, existing_table: bool
) -> None:
    record_path = tmp_path / madin_source.PROVENANCE_FILE
    record_path.write_bytes(record_bytes)
    target = tmp_path / MADIN_PIN.filename
    if existing_table:
        target.write_bytes(PAYLOAD)
    before = {path.name: path.read_bytes() for path in tmp_path.iterdir()}
    requests: list[httpx.Request] = []
    serve(monkeypatch, PAYLOAD, requests)

    with pytest.raises(MadinPinError, match=r"(?i)provenance|record|pin|JSON"):
        download_madin_table(tmp_path)

    assert requests == []
    assert {path.name: path.read_bytes() for path in tmp_path.iterdir()} == before


@pytest.mark.usefixtures("pinned_payload")
@pytest.mark.parametrize("existing_table", [False, True])
def test_atomic_record_failure_preserves_old_record_and_cleans_staging(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, existing_table: bool
) -> None:
    record_path = tmp_path / madin_source.PROVENANCE_FILE
    old_record = b'{"another.csv": {"public": true}}\n'
    record_path.write_bytes(old_record)
    target = tmp_path / MADIN_PIN.filename
    if existing_table:
        target.write_bytes(PAYLOAD)
    serve(monkeypatch, PAYLOAD, [] if not existing_table else None)
    original_replace = os.replace

    def fail_record_replacement(source: str | Path, destination: str | Path) -> None:
        if Path(destination) == record_path:
            raise OSError("injected atomic record publication failure")
        original_replace(source, destination)

    monkeypatch.setattr(os, "replace", fail_record_replacement)
    with pytest.raises(OSError, match="injected atomic record"):
        download_madin_table(tmp_path)

    assert record_path.read_bytes() == old_record
    assert {path.name for path in tmp_path.iterdir()} <= {record_path.name, target.name}
    if target.exists():
        assert target.read_bytes() == PAYLOAD


@pytest.mark.usefixtures("pinned_payload")
def test_wrong_size_is_rejected_without_hashing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    target = tmp_path / "oversized.csv"
    target.write_bytes(PAYLOAD + b"x")

    def forbid_hash(_path: Path) -> str:
        pytest.fail("A known size mismatch must be rejected before hashing")

    monkeypatch.setattr(madin_source, "file_sha256", forbid_hash)
    with pytest.raises(MadinPinError, match="bytes, expected"):
        verify_madin_table(target)


@pytest.mark.usefixtures("pinned_payload")
@pytest.mark.parametrize("through_download", [False, True])
def test_directory_table_is_a_pin_error(tmp_path: Path, through_download: bool) -> None:
    target = tmp_path / MADIN_PIN.filename
    target.mkdir()
    operation = download_madin_table if through_download else verify_madin_table
    argument = tmp_path if through_download else target
    with pytest.raises(MadinPinError, match=r"(?i)file|directory"):
        operation(argument)
    assert list(tmp_path.iterdir()) == [target]
