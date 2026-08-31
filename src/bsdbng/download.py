"""Download BugSigDB raw CSV exports."""

from __future__ import annotations

import json
import time
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import httpx

from bsdbng.paths import RAW_DATA_DIR

EXPORT_BASE = "https://bugsigdb.org/w/images/csv_reports"
EXPORT_FILES = ("studies.csv", "experiments.csv", "signatures.csv")
FULL_DUMP_URL = "https://raw.githubusercontent.com/waldronlab/BugSigDBExports/devel/full_dump.csv"
ZENODO_RECORD_ID = "19740674"
ZENODO_API_BASE = "https://zenodo.org/api/records"
PROVENANCE_FILE = "download_provenance.json"

MAX_RETRIES = 5
INITIAL_BACKOFF_SECONDS = 2.0


def _fetch_with_retry(
    client: httpx.Client,
    url: str,
    headers: dict[str, str] | None = None,
) -> httpx.Response:
    """Fetch a URL with exponential backoff on 429 and 5xx responses."""
    delay = INITIAL_BACKOFF_SECONDS
    for attempt in range(MAX_RETRIES):
        resp = client.get(url, headers=headers or {})
        is_retryable = resp.status_code == 429 or resp.status_code >= 500
        if is_retryable and attempt < MAX_RETRIES - 1:
            retry_after = resp.headers.get("Retry-After")
            wait = float(retry_after) if retry_after and retry_after.isdigit() else delay
            time.sleep(wait)
            delay *= 2
            continue
        # 304 Not Modified is a success — don't raise
        if resp.status_code != 304:
            resp.raise_for_status()
        return resp
    resp.raise_for_status()
    return resp  # unreachable, but keeps mypy happy


def download_exports(raw_dir: Path | None = None, *, force: bool = False) -> list[Path]:
    """Download the three BugSigDB CSV exports and full_dump.csv into *raw_dir*.

    Uses ETag and Last-Modified headers for conditional requests when
    provenance data is available. Skips download on 304 Not Modified.
    Falls back to skipping files that already exist when no provenance.
    Retries with exponential backoff on rate-limit (429) or server errors.
    Writes a provenance JSON sidecar recording download timestamps and cache headers.
    """
    raw_dir = raw_dir or RAW_DATA_DIR
    raw_dir.mkdir(parents=True, exist_ok=True)

    provenance_path = raw_dir / PROVENANCE_FILE
    provenance: dict[str, Any] = {}
    if provenance_path.exists():
        provenance = json.loads(provenance_path.read_text())

    urls: dict[str, str] = {name: f"{EXPORT_BASE}/{name}" for name in EXPORT_FILES}
    urls["full_dump.csv"] = FULL_DUMP_URL

    downloaded: list[Path] = []
    with httpx.Client(follow_redirects=True, timeout=120) as client:
        for name, url in urls.items():
            dest = raw_dir / name
            file_prov = provenance.get(name, {})
            if isinstance(file_prov, str):
                # Migrate old format: bare timestamp string → dict
                file_prov = {"downloaded_at": file_prov}

            # Build conditional request headers from cached provenance
            cond_headers: dict[str, str] = {}
            if dest.exists() and not force:
                etag = file_prov.get("etag")
                last_modified = file_prov.get("last_modified")
                if etag:
                    cond_headers["If-None-Match"] = etag
                if last_modified:
                    cond_headers["If-Modified-Since"] = last_modified
                if not cond_headers:
                    # No cache headers available — skip like before
                    downloaded.append(dest)
                    continue

            resp = _fetch_with_retry(client, url, headers=cond_headers)

            if resp.status_code == 304:
                # Upstream hasn't changed
                downloaded.append(dest)
                continue

            dest.write_bytes(resp.content)

            # Record provenance with cache headers for next conditional request
            file_prov = {
                "downloaded_at": datetime.now(tz=UTC).isoformat(),
                "etag": resp.headers.get("ETag"),
                "last_modified": resp.headers.get("Last-Modified"),
                "content_length": resp.headers.get("Content-Length"),
            }
            provenance[name] = file_prov
            downloaded.append(dest)

    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n")
    return downloaded


def download_zenodo(
    record_id: str | None = None,
    raw_dir: Path | None = None,
    *,
    force: bool = False,
    pattern: str = "full_dump.csv",
) -> list[Path]:
    """Download files matching *pattern* from Zenodo *record_id* into *raw_dir*.

    Defaults to record 19740674 (BugSigDBExports v1.3.1) and full_dump.csv.
    Uses provenance JSON to skip unchanged files via ETag/Last-Modified.
    """
    record_id = record_id or ZENODO_RECORD_ID
    raw_dir = raw_dir or RAW_DATA_DIR
    raw_dir.mkdir(parents=True, exist_ok=True)

    provenance_path = raw_dir / PROVENANCE_FILE
    provenance: dict[str, Any] = {}
    if provenance_path.exists():
        provenance = json.loads(provenance_path.read_text())

    api_url = f"{ZENODO_API_BASE}/{record_id}"
    downloaded: list[Path] = []

    with httpx.Client(follow_redirects=True, timeout=120) as client:
        resp = _fetch_with_retry(client, api_url)
        record = resp.json()
        files = record.get("files", [])

        for file_info in files:
            name = file_info.get("key") or file_info.get("filename")
            if not name or name != pattern:
                continue

            url = file_info["links"]["self"]
            dest = raw_dir / name
            file_prov = provenance.get(name, {})
            if isinstance(file_prov, str):
                file_prov = {"downloaded_at": file_prov}

            cond_headers: dict[str, str] = {}
            if dest.exists() and not force:
                etag = file_prov.get("etag")
                last_modified = file_prov.get("last_modified")
                if etag:
                    cond_headers["If-None-Match"] = etag
                if last_modified:
                    cond_headers["If-Modified-Since"] = last_modified
                if not cond_headers:
                    downloaded.append(dest)
                    continue

            file_resp = _fetch_with_retry(client, url, headers=cond_headers)
            if file_resp.status_code == 304:
                downloaded.append(dest)
                continue

            dest.write_bytes(file_resp.content)
            provenance[name] = {
                "downloaded_at": datetime.now(tz=UTC).isoformat(),
                "etag": file_resp.headers.get("ETag"),
                "last_modified": file_resp.headers.get("Last-Modified"),
                "content_length": file_resp.headers.get("Content-Length"),
                "zenodo_record": record_id,
            }
            downloaded.append(dest)

    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n")
    return downloaded
