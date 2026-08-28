import json
from pathlib import Path

import httpx
import pytest

from bsdbng import download as dl
from bsdbng.download import (
    BUGSIGDB_EXPORTS_DOI,
    BUGSIGDB_EXPORTS_VERSION,
    DEVEL_FULL_DUMP_URL,
    EXPORT_FILES,
    FULL_DUMP_NAME,
    PINNED_FULL_DUMP_URL,
    PROVENANCE_FILE,
    SOURCE_ENV_VAR,
    SourceResolutionError,
    describe_source,
    download_exports,
    parse_snapshot,
    release_descriptor,
    resolve_source,
)

SNAPSHOT_HEADER = (
    b"# BugSigDB 2026-04-24_00:41_UTC, License: Creative Commons Attribution 4.0 "
    b"International, URL: https://bugsigdb.org\n"
)


def _seed_raw_dir(tmp_path: Path) -> None:
    for name in EXPORT_FILES:
        (tmp_path / name).write_text("existing")
    (tmp_path / FULL_DUMP_NAME).write_text("existing")


def _mock_client(monkeypatch: pytest.MonkeyPatch, handler) -> None:
    """Route download_exports' httpx.Client through a MockTransport.

    Bind the real class first: dl.httpx *is* the httpx module, so patching its
    Client attribute is global, and a lambda that called httpx.Client would
    recurse into itself.
    """
    real_client = httpx.Client
    monkeypatch.setattr(
        dl.httpx,
        "Client",
        lambda **_kw: real_client(transport=httpx.MockTransport(handler)),
    )


def _skip_integrity(monkeypatch: pytest.MonkeyPatch) -> None:
    """Stub the md5/size check for tests that exercise other logic.

    These tests feed synthetic payloads through the pinned path to check
    provenance, caching and outage behaviour. The integrity check itself has
    dedicated tests below.
    """
    monkeypatch.setattr(dl, "verify_pinned_dump", lambda _content: None)


def test_download_skips_existing_files(tmp_path: Path) -> None:
    """When all target files exist, download_exports should not overwrite them."""
    _seed_raw_dir(tmp_path)

    result = download_exports(tmp_path, source="live")

    assert len(result) == 4
    # Files should still have original content (not re-fetched)
    assert (tmp_path / "studies.csv").read_text() == "existing"


def test_provenance_written(tmp_path: Path) -> None:
    """Provenance file should exist after download, even when all files are cached."""
    _seed_raw_dir(tmp_path)

    download_exports(tmp_path, source="live")

    prov_path = tmp_path / PROVENANCE_FILE
    assert prov_path.exists()
    data = json.loads(prov_path.read_text())
    assert isinstance(data, dict)


# --- source resolution (#115) ---------------------------------------------


def test_default_source_is_the_pinned_release(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv(SOURCE_ENV_VAR, raising=False)
    assert resolve_source() == "pinned"


def test_environment_overrides_the_default(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv(SOURCE_ENV_VAR, "LIVE")
    assert resolve_source() == "live"
    # An explicit argument still wins over the environment.
    assert resolve_source("pinned") == "pinned"


def test_unknown_source_is_rejected() -> None:
    with pytest.raises(SourceResolutionError, match="Unknown BugSigDB source"):
        resolve_source("v1.3.0")


def test_pinned_descriptor_carries_a_citable_identity() -> None:
    release = release_descriptor("pinned")
    assert release["pinned"] is True
    assert release["version"] == BUGSIGDB_EXPORTS_VERSION
    assert release["doi"] == BUGSIGDB_EXPORTS_DOI
    assert release["url"] == PINNED_FULL_DUMP_URL
    assert release["tag_sha"]


def test_live_descriptor_claims_nothing() -> None:
    release = release_descriptor("live")
    assert release["pinned"] is False
    assert release["version"] is None
    assert release["doi"] is None
    assert release["url"] == DEVEL_FULL_DUMP_URL


def test_describe_source_names_the_release_and_the_unpinned_remainder() -> None:
    summary = describe_source(release_descriptor("pinned"))
    assert BUGSIGDB_EXPORTS_VERSION in summary
    assert BUGSIGDB_EXPORTS_DOI in summary
    assert "unpinned" in summary


def test_live_source_is_labelled_as_not_pinned() -> None:
    assert "NOT pinned" in describe_source(release_descriptor("live"))


# --- snapshot header ------------------------------------------------------


def test_parse_snapshot_reads_the_leading_comment() -> None:
    assert parse_snapshot(SNAPSHOT_HEADER + b"BSDB ID,Study\n") == "2026-04-24_00:41_UTC"


def test_parse_snapshot_returns_none_without_a_header() -> None:
    assert parse_snapshot(b"BSDB ID,Study\nbsdb:1/1/1,x\n") is None


# --- provenance integrity -------------------------------------------------


def test_pinned_run_records_the_release_and_flags_the_wiki_csvs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A pinned run must record the release and mark only full_dump.csv as pinned."""
    requested: list[str] = []

    def handler(request: httpx.Request) -> httpx.Response:
        requested.append(str(request.url))
        return httpx.Response(200, content=SNAPSHOT_HEADER, headers={"ETag": '"abc"'})

    _skip_integrity(monkeypatch)
    _mock_client(monkeypatch, handler)

    download_exports(tmp_path, source="pinned")

    prov = json.loads((tmp_path / PROVENANCE_FILE).read_text())
    assert prov["bugsigdb_release"]["doi"] == BUGSIGDB_EXPORTS_DOI
    assert prov["bugsigdb_release"]["version"] == BUGSIGDB_EXPORTS_VERSION
    assert prov[FULL_DUMP_NAME]["pinned"] is True
    assert prov[FULL_DUMP_NAME]["url"] == PINNED_FULL_DUMP_URL
    for name in EXPORT_FILES:
        assert prov[name]["pinned"] is False, f"{name} is a live wiki table, not a release asset"
    assert PINNED_FULL_DUMP_URL in requested
    assert DEVEL_FULL_DUMP_URL not in requested


def test_snapshot_timestamp_is_recorded_even_for_a_live_run(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A live run cannot be repeated, but it must still be identifiable."""

    def handler(_request: httpx.Request) -> httpx.Response:
        return httpx.Response(200, content=SNAPSHOT_HEADER)

    _mock_client(monkeypatch, handler)

    download_exports(tmp_path, source="live")

    prov = json.loads((tmp_path / PROVENANCE_FILE).read_text())
    assert prov[FULL_DUMP_NAME]["bugsigdb_snapshot"] == "2026-04-24_00:41_UTC"
    assert prov["bugsigdb_release"]["bugsigdb_snapshot"] == "2026-04-24_00:41_UTC"
    assert prov["bugsigdb_release"]["pinned"] is False


def test_switching_source_does_not_reuse_the_other_sources_etag(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A cached ETag belongs to the URL it came from.

    Sent to a different source it is meaningless, and a bogus 304 would leave
    the previous vintage's bytes on disk while the sidecar claimed the new one.
    """
    (tmp_path / FULL_DUMP_NAME).write_text("devel-vintage")
    for name in EXPORT_FILES:
        (tmp_path / name).write_text("existing")
    (tmp_path / PROVENANCE_FILE).write_text(
        json.dumps(
            {
                FULL_DUMP_NAME: {
                    "downloaded_at": "2026-01-01T00:00:00+00:00",
                    "url": DEVEL_FULL_DUMP_URL,
                    "etag": '"devel-etag"',
                }
            }
        )
    )

    conditional_headers: list[str | None] = []

    def handler(request: httpx.Request) -> httpx.Response:
        if str(request.url) == PINNED_FULL_DUMP_URL:
            conditional_headers.append(request.headers.get("If-None-Match"))
            return httpx.Response(200, content=b"pinned-vintage", headers={"ETag": '"p"'})
        return httpx.Response(200, content=b"wiki")

    _skip_integrity(monkeypatch)
    _mock_client(monkeypatch, handler)

    download_exports(tmp_path, source="pinned")

    assert conditional_headers == [None], "the devel ETag must not be sent to the pinned tag"
    assert (tmp_path / FULL_DUMP_NAME).read_bytes() == b"pinned-vintage"


def test_unprovenanced_full_dump_is_refetched_when_pinning(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Bytes of unknown origin must not be claimed as the pinned release."""
    _seed_raw_dir(tmp_path)  # files present, no provenance at all

    _skip_integrity(monkeypatch)
    _mock_client(monkeypatch, lambda _r: httpx.Response(200, content=b"pinned-vintage"))

    download_exports(tmp_path, source="pinned")

    assert (tmp_path / FULL_DUMP_NAME).read_bytes() == b"pinned-vintage"
    # The unpinned wiki CSVs claim nothing, so the historical skip still applies.
    assert (tmp_path / "studies.csv").read_text() == "existing"


def test_wiki_csv_outage_does_not_abort_a_pinned_run(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The wiki CSVs are not ingest inputs, so their outage must not be fatal."""

    def handler(request: httpx.Request) -> httpx.Response:
        if str(request.url) == PINNED_FULL_DUMP_URL:
            return httpx.Response(200, content=SNAPSHOT_HEADER)
        return httpx.Response(404)

    _skip_integrity(monkeypatch)
    _mock_client(monkeypatch, handler)

    result = download_exports(tmp_path, source="pinned")

    assert [p.name for p in result] == [FULL_DUMP_NAME]
    prov = json.loads((tmp_path / PROVENANCE_FILE).read_text())
    assert "last_error" in prov["studies.csv"]
    assert prov["studies.csv"]["attempted_url"].endswith("studies.csv")
    assert prov[FULL_DUMP_NAME]["pinned"] is True


def test_full_dump_outage_is_fatal(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Losing the ingest input must fail loudly rather than silently."""

    def handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(404 if str(request.url) == PINNED_FULL_DUMP_URL else 200)

    _mock_client(monkeypatch, handler)

    with pytest.raises(httpx.HTTPStatusError):
        download_exports(tmp_path, source="pinned")


# --- pin strength (#151, #152) --------------------------------------------


def test_pinned_url_uses_the_commit_sha_not_the_mutable_tag_name() -> None:
    """Tags can be re-pointed; a SHA cannot. See #151."""
    assert dl.BUGSIGDB_EXPORTS_TAG_SHA in PINNED_FULL_DUMP_URL
    assert f"/{BUGSIGDB_EXPORTS_VERSION}/" not in PINNED_FULL_DUMP_URL


def test_verify_pinned_dump_rejects_a_truncated_download() -> None:
    with pytest.raises(dl.IntegrityError, match="bytes, expected"):
        dl.verify_pinned_dump(b"truncated")


def test_verify_pinned_dump_rejects_a_wrong_payload_of_the_right_size() -> None:
    with pytest.raises(dl.IntegrityError, match="md5"):
        dl.verify_pinned_dump(b"x" * dl.PINNED_FULL_DUMP_BYTES)


def test_a_corrupt_pinned_download_is_never_written_to_disk(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Bad bytes must not land on disk and then be recorded as pinned."""

    def handler(request: httpx.Request) -> httpx.Response:
        if str(request.url) == PINNED_FULL_DUMP_URL:
            return httpx.Response(200, content=b"<html>proxy error</html>")
        return httpx.Response(200, content=b"wiki")

    _mock_client(monkeypatch, handler)

    with pytest.raises(dl.IntegrityError):
        download_exports(tmp_path, source="pinned")

    assert not (tmp_path / FULL_DUMP_NAME).exists()


def test_live_downloads_are_not_integrity_checked(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """devel has no published checksum, so the check must not fire there."""
    _mock_client(monkeypatch, lambda _r: httpx.Response(200, content=b"whatever devel says"))

    download_exports(tmp_path, source="live")

    assert (tmp_path / FULL_DUMP_NAME).read_bytes() == b"whatever devel says"
