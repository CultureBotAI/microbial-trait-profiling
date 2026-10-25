"""`just download-madin` fetches only the pinned Madin et al. table (#1048)."""

from __future__ import annotations

import hashlib
import json
import re
from dataclasses import replace
from pathlib import Path
from typing import Any

import httpx
import pytest

from mtp import madin_source
from mtp.madin_source import MADIN_PIN, MadinPinError, download_madin_table, verify_madin_table

PAYLOAD = b'"tax_id","species_tax_id","data_source"\n11,11,"gold"\n'


@pytest.fixture
def pinned_payload(monkeypatch: pytest.MonkeyPatch) -> None:
    """PAYLOAD stands in for the pinned bytes; repository, commit and path stay real."""
    monkeypatch.setattr(
        madin_source,
        "MADIN_PIN",
        replace(MADIN_PIN, sha256=hashlib.sha256(PAYLOAD).hexdigest(), size_bytes=len(PAYLOAD)),
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
    assert re.fullmatch(r"[0-9a-f]{64}", MADIN_PIN.sha256)
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
