"""Download BugSigDB CSV exports."""

from __future__ import annotations

import json
import time
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import httpx

from bsdbng.paths import RAW_DATA_DIR

EXPORT_BASE = "https://bugsigdb-csv.s3.us-east-va.perf.cloud.ovh.us"
EXPORT_FILES = ("studies.csv", "experiments.csv", "signatures.csv")
BUGSIGDB_EXPORTS_VERSION = "v1.3.1"
BUGSIGDB_EXPORTS_DOI = "10.5281/zenodo.19740674"
# Public Git revision identifying the pinned export, not a credential.
BUGSIGDB_EXPORTS_TAG_SHA = "5815c2638d6b56134ac8db19edc2bb3c4b17f1ca"  # pragma: allowlist secret
FULL_DUMP_URL = (
    "https://raw.githubusercontent.com/waldronlab/BugSigDBExports/"
    f"{BUGSIGDB_EXPORTS_VERSION}/full_dump.csv"
)
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
    """Download BugSigDB raw CSV exports and pinned release full dump into *raw_dir*.

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
                "url": url,
                "etag": resp.headers.get("ETag"),
                "last_modified": resp.headers.get("Last-Modified"),
                "content_length": resp.headers.get("Content-Length"),
                "downloaded_bytes": len(resp.content),
            }
            if name == "full_dump.csv":
                file_prov.update(
                    {
                        "source": "waldronlab/BugSigDBExports",
                        "version": BUGSIGDB_EXPORTS_VERSION,
                        "doi": BUGSIGDB_EXPORTS_DOI,
                        "tag_sha": BUGSIGDB_EXPORTS_TAG_SHA,
                    }
                )
            provenance[name] = file_prov
            downloaded.append(dest)

    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n")
    return downloaded
