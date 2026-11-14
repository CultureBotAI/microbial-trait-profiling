"""Download BugSigDB CSV exports.

Provenance model
----------------
The ingest (``mtp.ingest``) reads **only** ``full_dump.csv`` — the separate-CSV
path is a stub that raises ``NotImplementedError``. So pinning ``full_dump.csv``
pins the ingest, and it is the one file this module fetches from a citable
release by default.

The three MediaWiki exports (``studies.csv``, ``experiments.csv``,
``signatures.csv``) are live wiki tables with no pinned equivalent — the
BugSigDBExports release does not contain them. They are still fetched, because
``mtp.stats`` reports on them and downstream analyses join against
``experiments.csv``, but they are **not pinned** and the provenance sidecar says
so per file. Do not read "the ingest is pinned" as "everything in data/raw is
pinned".

Every dump also carries its own BugSigDB snapshot timestamp on line 1::

    # BugSigDB 2026-04-24_00:41_UTC, License: ..., URL: https://bugsigdb.org

``ingest._read_csv`` strips ``#`` lines, so nothing else reads it. It is parsed
into the sidecar here, which means the exact vintage is recorded even for an
unpinned ``live`` run: such a run cannot be repeated, but it can at least be
identified after the fact.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import time
from datetime import UTC, datetime
from pathlib import Path, PurePosixPath
from typing import Any

import httpx

from mtp.paths import RAW_DATA_DIR

# BugSigDB moved its CSV exports off the wiki; bugsigdb.org/w/images/csv_reports
# now 301s to a /W/ path that 404s.
EXPORT_BASE = "https://bugsigdb-csv.s3.us-east-va.perf.cloud.ovh.us"
EXPORT_FILES = ("studies.csv", "experiments.csv", "signatures.csv")
FULL_DUMP_NAME = "full_dump.csv"

# v1.3.1 is a REVIEWED-ONLY export: every row carries State=Complete and a named
# Reviewer, which v1.3.0 and devel do not. It ingests fewer studies than v1.3.0
# despite being a year newer. See issue #130 before changing this pin.
BUGSIGDB_EXPORTS_VERSION = "v1.3.1"
BUGSIGDB_EXPORTS_DOI = "10.5281/zenodo.19740674"
# Fetch by COMMIT SHA, not by tag name (#151). Tags are mutable: if the release
# were re-pointed, a tag-name URL would silently serve different bytes under the
# same version and DOI -- the moving-target failure this pinning exists to
# remove. The version and DOI stay as metadata; the SHA is what is fetched.
BUGSIGDB_EXPORTS_TAG_SHA = "5815c2638d6b56134ac8db19edc2bb3c4b17f1ca"
# md5 and size of full_dump.csv as published in the Zenodo record. The file
# served from the git tag is byte-identical to the Zenodo asset, so one pair of
# constants validates both access paths (#152).
PINNED_FULL_DUMP_MD5 = "00df5a86c684eac19fd767c029ff0130"
PINNED_FULL_DUMP_BYTES = 16241290
PINNED_FULL_DUMP_URL = (
    "https://raw.githubusercontent.com/waldronlab/BugSigDBExports/"
    f"{BUGSIGDB_EXPORTS_TAG_SHA}/full_dump.csv"
)
# Retained for the opt-out path only. `devel` is a moving branch: two runs on
# different days can differ, with no way to say which vintage was used.
DEVEL_FULL_DUMP_URL = (
    "https://raw.githubusercontent.com/waldronlab/BugSigDBExports/devel/full_dump.csv"
)
FULL_DUMP_URL = PINNED_FULL_DUMP_URL  # back-compat for importers

PROVENANCE_FILE = "download_provenance.json"

SOURCE_PINNED = "pinned"
SOURCE_LIVE = "live"
SOURCE_RELEASE = "release"
SOURCE_ENV_VAR = "BSDBNG_BUGSIGDB_SOURCE"
DEFAULT_SOURCE = SOURCE_PINNED

DEFAULT_RELEASE_REPO = "CultureBotAI/microbial-trait-profiling"
RELEASE_ARCHIVE_PATTERN = "bsdbng-raw-data-*.tar.gz"
RELEASE_EXPECTED_CSVS = (*EXPORT_FILES, FULL_DUMP_NAME)
_REPO_SLUG_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*/[A-Za-z0-9][A-Za-z0-9_.-]*$")
_TAG_SLUG_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
_ARCHIVE_TAG_RE = re.compile(r"^bsdbng-raw-data-(.+)\.tar\.gz$")

SNAPSHOT_RE = re.compile(rb"^#\s*BugSigDB\s+(\S+?),")

MAX_RETRIES = 5
INITIAL_BACKOFF_SECONDS = 2.0


class SourceResolutionError(RuntimeError):
    """Raised when the requested ingest source cannot be resolved."""


class IntegrityError(RuntimeError):
    """Raised when downloaded bytes do not match the pinned release."""


def verify_pinned_dump(content: bytes) -> None:
    """Check downloaded bytes against the release's published md5 and size.

    The pinned path's whole value is that two runs read the same bytes. Without
    this a truncated response or a proxy error page would be written to disk and
    then recorded as ``"pinned": true`` with a DOI attached (#152).
    """
    if len(content) != PINNED_FULL_DUMP_BYTES:
        msg = (
            f"{FULL_DUMP_NAME} from {BUGSIGDB_EXPORTS_VERSION} is "
            f"{len(content):,} bytes, expected {PINNED_FULL_DUMP_BYTES:,}."
        )
        raise IntegrityError(msg)
    digest = hashlib.md5(content, usedforsecurity=False).hexdigest()
    if digest != PINNED_FULL_DUMP_MD5:
        msg = (
            f"{FULL_DUMP_NAME} from {BUGSIGDB_EXPORTS_VERSION} has md5 {digest}, "
            f"expected {PINNED_FULL_DUMP_MD5} (Zenodo {BUGSIGDB_EXPORTS_DOI})."
        )
        raise IntegrityError(msg)


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


def resolve_source(source: str | None = None) -> str:
    """Resolve the requested ingest source, honouring the environment override."""
    resolved = (source or os.environ.get(SOURCE_ENV_VAR) or DEFAULT_SOURCE).strip().lower()
    if resolved not in (SOURCE_PINNED, SOURCE_LIVE):
        msg = (
            f"Unknown BugSigDB source {resolved!r}. "
            f"Expected {SOURCE_PINNED!r} (release {BUGSIGDB_EXPORTS_VERSION}) "
            f"or {SOURCE_LIVE!r} (moving devel branch)."
        )
        raise SourceResolutionError(msg)
    return resolved


def parse_snapshot(content: bytes) -> str | None:
    """Read the BugSigDB snapshot timestamp from a dump's leading comment line."""
    first, _, _ = content.partition(b"\n")
    match = SNAPSHOT_RE.match(first)
    return match.group(1).decode("utf-8", "replace") if match else None


def release_descriptor(source: str) -> dict[str, Any]:
    """Describe the resolved full-dump source, recorded verbatim in the sidecar."""
    if source == SOURCE_PINNED:
        return {
            "source": "waldronlab/BugSigDBExports",
            "pinned": True,
            "version": BUGSIGDB_EXPORTS_VERSION,
            "doi": BUGSIGDB_EXPORTS_DOI,
            "tag_sha": BUGSIGDB_EXPORTS_TAG_SHA,
            "url": PINNED_FULL_DUMP_URL,
        }
    return {
        "source": "waldronlab/BugSigDBExports @ devel (moving branch, not a release)",
        "pinned": False,
        "version": None,
        "doi": None,
        "tag_sha": None,
        "url": DEVEL_FULL_DUMP_URL,
    }


def describe_source(release: dict[str, Any]) -> str:
    """One-line, log-friendly summary of what the ingest is about to read."""
    if release["pinned"]:
        return (
            f"BugSigDB source: pinned release {release['version']} ({release['doi']}) — "
            f"{FULL_DUMP_NAME} only; the three MediaWiki CSVs remain live and unpinned."
        )
    return (
        "BugSigDB source: LIVE devel branch (NOT pinned) — re-running at a "
        "different time may silently yield different input data."
    )


def download_exports(
    raw_dir: Path | None = None,
    *,
    force: bool = False,
    source: str | None = None,
) -> list[Path]:
    """Download the three BugSigDB CSV exports and full_dump.csv into *raw_dir*.

    ``source`` selects where ``full_dump.csv`` comes from: ``"pinned"`` (default)
    uses BugSigDBExports ``v1.3.1``, ``"live"`` uses the moving devel branch. The
    environment variable ``BSDBNG_BUGSIGDB_SOURCE`` overrides the default when
    ``source`` is not passed explicitly.

    Uses ETag and Last-Modified headers for conditional requests when
    provenance data is available. Skips download on 304 Not Modified.
    Verifies pinned cached bytes before reuse; corrupt caches raise
    ``IntegrityError`` and can be replaced with a verified ``force=True`` fetch.
    Unpinned files without provenance retain skip-if-present behavior.
    Retries with exponential backoff on rate-limit (429) or server errors.
    Writes a provenance JSON sidecar recording the resolved release, the
    BugSigDB snapshot timestamp, download timestamps and cache headers.
    """
    raw_dir = raw_dir or RAW_DATA_DIR
    raw_dir.mkdir(parents=True, exist_ok=True)

    resolved_source = resolve_source(source)
    release = release_descriptor(resolved_source)
    print(describe_source(release))

    provenance_path = raw_dir / PROVENANCE_FILE
    provenance: dict[str, Any] = {}
    if provenance_path.exists():
        provenance = json.loads(provenance_path.read_text())

    provenance["bugsigdb_release"] = {
        **release,
        "resolved_at": datetime.now(tz=UTC).isoformat(),
    }

    urls: dict[str, str] = {name: f"{EXPORT_BASE}/{name}" for name in EXPORT_FILES}
    urls[FULL_DUMP_NAME] = release["url"]

    downloaded: list[Path] = []
    with httpx.Client(follow_redirects=True, timeout=120) as client:
        for name, url in urls.items():
            dest = raw_dir / name
            file_prov = provenance.get(name, {})
            if isinstance(file_prov, str):
                # Migrate old format: bare timestamp string → dict
                file_prov = {"downloaded_at": file_prov}

            # A cached ETag belongs to the URL it came from. Sending the devel
            # branch's ETag to the pinned tag (or vice versa) is meaningless,
            # and a bogus 304 would leave the previous source's bytes on disk
            # while the sidecar claimed the new release.
            cached_url = file_prov.get("url")
            is_pinned_target = name == FULL_DUMP_NAME and resolved_source == SOURCE_PINNED
            # Pre-existing bytes with no recorded source cannot be *claimed* as
            # the pinned release, so re-fetch rather than assert a vintage we
            # cannot substantiate. For the unpinned wiki CSVs nothing is being
            # claimed, so the historical skip-if-present behaviour stands.
            reusable = cached_url == url or (cached_url is None and not is_pinned_target)

            # Build conditional request headers from cached provenance
            cond_headers: dict[str, str] = {}
            if dest.exists() and not force and reusable:
                if is_pinned_target:
                    # A matching URL or ETag identifies the remote object, not
                    # the bytes currently on disk. Never reuse a corrupt cache.
                    verify_pinned_dump(dest.read_bytes())
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

            try:
                resp = _fetch_with_retry(client, url, headers=cond_headers)
            except httpx.HTTPError as exc:
                # full_dump.csv is the ingest input — losing it is fatal. The
                # three wiki CSVs are not read by the ingest at all, so an
                # outage there must not block a pinned run. The export host has
                # gone down for hours at a time before, taking the whole
                # download with it including files that fetched fine.
                if name == FULL_DUMP_NAME:
                    raise
                provenance[name] = {
                    **file_prov,
                    "last_error": str(exc),
                    "last_error_at": datetime.now(tz=UTC).isoformat(),
                    "attempted_url": url,
                }
                state = "keeping the existing copy" if dest.exists() else "no local copy"
                print(f"WARNING: {name} could not be fetched ({exc}) — {state}.")
                if dest.exists():
                    downloaded.append(dest)
                continue

            if resp.status_code == 304:
                # Upstream hasn't changed
                if is_pinned_target:
                    if not dest.is_file():
                        raise IntegrityError(
                            f"{FULL_DUMP_NAME} returned HTTP 304 without a local cached file."
                        )
                    verify_pinned_dump(dest.read_bytes())
                downloaded.append(dest)
                continue

            if is_pinned_target:
                verify_pinned_dump(resp.content)

            dest.write_bytes(resp.content)

            # Record provenance with cache headers for next conditional request
            file_prov = {
                "downloaded_at": datetime.now(tz=UTC).isoformat(),
                "url": url,
                "pinned": is_pinned_target,
                "etag": resp.headers.get("ETag"),
                "last_modified": resp.headers.get("Last-Modified"),
                "content_length": resp.headers.get("Content-Length"),
                "downloaded_bytes": len(resp.content),
            }
            if is_pinned_target:
                file_prov["md5"] = PINNED_FULL_DUMP_MD5
                file_prov["integrity_verified"] = True
            snapshot = parse_snapshot(resp.content)
            if snapshot:
                # The vintage the data itself declares, as opposed to the vintage
                # the URL implies. Present for pinned and live runs alike.
                file_prov["bugsigdb_snapshot"] = snapshot
            if name == FULL_DUMP_NAME:
                file_prov.update(
                    {
                        "source": release["source"],
                        "version": release["version"],
                        "doi": release["doi"],
                        "tag_sha": release["tag_sha"],
                    }
                )
                if snapshot:
                    provenance["bugsigdb_release"]["bugsigdb_snapshot"] = snapshot
            provenance[name] = file_prov
            downloaded.append(dest)

    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n")
    return downloaded


def _validate_release_repo(repo: str) -> str:
    cleaned = repo.strip()
    if not _REPO_SLUG_RE.match(cleaned):
        raise ValueError(f"Invalid GitHub repository slug {repo!r}; expected 'owner/repo'")
    return cleaned


def _validate_release_tag(tag: str | None) -> str | None:
    if tag is None:
        return None
    cleaned = tag.strip()
    if not cleaned or cleaned.lower() == "latest":
        return None
    if ".." in cleaned or not _TAG_SLUG_RE.match(cleaned):
        raise ValueError(f"Invalid release tag {tag!r}")
    return cleaned


def _fetch_release_archive(dest_dir: Path, *, tag: str | None, repo: str) -> Path:
    import shutil

    gh_bin = shutil.which("gh")
    if not gh_bin:
        raise SourceResolutionError(
            "gh CLI not found on PATH; install GitHub CLI to run download-from-release."
        )
    cmd = [gh_bin, "release", "download"]
    if tag is not None:
        cmd.append(tag)
    cmd.extend(
        [
            "--repo",
            repo,
            "--dir",
            str(dest_dir),
            "--pattern",
            RELEASE_ARCHIVE_PATTERN,
        ]
    )
    pid = os.posix_spawn(gh_bin, cmd, os.environ)
    _, wait_status = os.waitpid(pid, 0)
    exit_code = os.waitstatus_to_exitcode(wait_status)
    if exit_code != 0:
        raise SourceResolutionError(
            f"gh release download failed with exit code {exit_code} "
            f"for repo={repo!r} tag={tag or 'latest'!r}."
        )
    matches = sorted(dest_dir.glob(RELEASE_ARCHIVE_PATTERN))
    if not matches:
        raise SourceResolutionError(
            f"No release archive matching {RELEASE_ARCHIVE_PATTERN!r} downloaded "
            f"from {repo!r} ({tag or 'latest'})."
        )
    if len(matches) > 1:
        raise IntegrityError(
            f"Multiple release archives matched {RELEASE_ARCHIVE_PATTERN!r}: "
            f"{[m.name for m in matches]}"
        )
    return matches[0]


def _infer_release_tag(archive_name: str, requested_tag: str | None) -> str:
    if requested_tag is not None:
        return requested_tag
    match = _ARCHIVE_TAG_RE.match(archive_name)
    return match.group(1) if match else "latest"


def _read_release_archive_members(archive_path: Path) -> dict[str, bytes]:
    import tarfile

    expected_set = set(RELEASE_EXPECTED_CSVS)
    ignorable_set = {PROVENANCE_FILE, ".gitkeep"}
    extracted: dict[str, bytes] = {}
    with tarfile.open(archive_path, mode="r:gz") as tar:
        for member in tar.getmembers():
            raw_name = member.name
            if not raw_name or raw_name.startswith("/") or "\\" in raw_name:
                raise IntegrityError(f"Unsafe path in release archive: {raw_name!r}")
            parts = PurePosixPath(raw_name).parts
            if ".." in parts:
                raise IntegrityError(f"Unsafe path traversal in release archive: {raw_name!r}")
            if (
                member.issym()
                or member.islnk()
                or member.isdev()
                or member.isfifo()
                or member.ischr()
                or member.isblk()
            ):
                raise IntegrityError(
                    f"Refusing non-regular tar member {raw_name!r} in release archive"
                )
            if member.isdir():
                if parts != ():
                    raise IntegrityError(
                        f"Unexpected directory member {raw_name!r} in release archive"
                    )
                continue
            if not member.isfile():
                raise IntegrityError(
                    f"Refusing non-regular tar member {raw_name!r} in release archive"
                )
            if len(parts) != 1:
                raise IntegrityError(f"Unexpected nested member {raw_name!r} in release archive")
            filename = parts[0]
            if filename in ignorable_set:
                continue
            if filename not in expected_set:
                raise IntegrityError(
                    f"Unexpected file {raw_name!r} in release archive; "
                    f"expected only {sorted(expected_set)}"
                )
            if filename in extracted:
                raise IntegrityError(f"Duplicate member {filename!r} in release archive")
            file_obj = tar.extractfile(member)
            if file_obj is None:
                raise IntegrityError(f"Unable to read tar member {raw_name!r}")
            with file_obj:
                extracted[filename] = file_obj.read()
    if FULL_DUMP_NAME not in extracted or not extracted[FULL_DUMP_NAME]:
        raise IntegrityError(
            f"Release archive {archive_path.name!r} is missing non-empty {FULL_DUMP_NAME}"
        )
    return extracted


def download_from_release(
    tag: str | None = None,
    output_dir: Path | None = None,
    repo: str = DEFAULT_RELEASE_REPO,
) -> list[Path]:
    """Download and extract BugSigDB CSVs from a GitHub release tarball (#92).

    Resolves ``bsdbng-raw-data-*.tar.gz`` via ``gh release download`` into a
    temporary directory, validates every archive member against path traversal
    and unexpected entries, verifies ``full_dump.csv`` against the pinned digest
    when its size or checksum matches the pinned release, writes the CSV files
    into ``output_dir``, and records ``download_provenance.json`` with
    ``"source": "release"``.
    """
    validated_repo = _validate_release_repo(repo)
    resolved_tag = _validate_release_tag(tag)
    target_dir = output_dir if output_dir is not None else RAW_DATA_DIR
    import tempfile

    with tempfile.TemporaryDirectory(prefix="mtp-release-download-") as tmp_str:
        tmp_dir = Path(tmp_str)
        archive_path = _fetch_release_archive(tmp_dir, tag=resolved_tag, repo=validated_repo)
        archive_bytes = archive_path.read_bytes()
        archive_sha256 = hashlib.sha256(archive_bytes).hexdigest()
        archive_name = archive_path.name
        release_tag = _infer_release_tag(archive_name, resolved_tag)
        extracted_payloads = _read_release_archive_members(archive_path)

    full_dump_bytes = extracted_payloads[FULL_DUMP_NAME]
    full_dump_md5 = hashlib.md5(full_dump_bytes, usedforsecurity=False).hexdigest()
    pinned_match = (
        len(full_dump_bytes) == PINNED_FULL_DUMP_BYTES and full_dump_md5 == PINNED_FULL_DUMP_MD5
    )
    if len(full_dump_bytes) == PINNED_FULL_DUMP_BYTES or full_dump_md5 == PINNED_FULL_DUMP_MD5:
        verify_pinned_dump(full_dump_bytes)

    target_dir.mkdir(parents=True, exist_ok=True)
    now_iso = datetime.now(tz=UTC).isoformat()
    full_dump_snapshot = parse_snapshot(full_dump_bytes)

    provenance: dict[str, Any] = {
        "source": SOURCE_RELEASE,
        "release_tag": release_tag,
        "release_repo": validated_repo,
        "archive_name": archive_name,
        "archive_sha256": archive_sha256,
        "downloaded_at": now_iso,
        "bugsigdb_release": {
            "source": f"github-release:{validated_repo}@{release_tag}",
            "pinned": pinned_match,
            "version": BUGSIGDB_EXPORTS_VERSION if pinned_match else release_tag,
            "doi": BUGSIGDB_EXPORTS_DOI if pinned_match else None,
            "tag_sha": BUGSIGDB_EXPORTS_TAG_SHA if pinned_match else None,
            "release_tag": release_tag,
            "release_repo": validated_repo,
            "archive_name": archive_name,
            "archive_sha256": archive_sha256,
            "resolved_at": now_iso,
            **({"bugsigdb_snapshot": full_dump_snapshot} if full_dump_snapshot is not None else {}),
        },
    }

    written_paths: list[Path] = []
    for name in RELEASE_EXPECTED_CSVS:
        if name not in extracted_payloads:
            continue
        payload = extracted_payloads[name]
        dest = target_dir / name
        dest.write_bytes(payload)
        written_paths.append(dest)
        file_md5 = hashlib.md5(payload, usedforsecurity=False).hexdigest()
        file_sha256 = hashlib.sha256(payload).hexdigest()
        is_pinned_file = name == FULL_DUMP_NAME and pinned_match
        file_prov: dict[str, Any] = {
            "downloaded_at": now_iso,
            "source": SOURCE_RELEASE,
            "release_tag": release_tag,
            "release_repo": validated_repo,
            "archive_name": archive_name,
            "archive_sha256": archive_sha256,
            "pinned": is_pinned_file,
            "downloaded_bytes": len(payload),
            "md5": file_md5,
            "sha256": file_sha256,
        }
        snapshot = parse_snapshot(payload)
        if snapshot is not None:
            file_prov["bugsigdb_snapshot"] = snapshot
        if is_pinned_file:
            file_prov["integrity_verified"] = True
        provenance[name] = file_prov

    (target_dir / PROVENANCE_FILE).write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
    )
    return written_paths


def build_parser() -> Any:
    """Create CLI argument parser for ``python -m mtp.download`` (#92)."""
    import argparse

    parser = argparse.ArgumentParser(
        prog="python -m mtp.download",
        description="Download BugSigDB CSV exports or restore them from a GitHub release tarball.",
    )
    parser.add_argument(
        "--from-release",
        nargs="?",
        const="latest",
        default=None,
        metavar="TAG",
        help="Download raw CSVs from a GitHub release tarball (default tag: latest).",
    )
    parser.add_argument(
        "--repo",
        default=DEFAULT_RELEASE_REPO,
        help="GitHub owner/repo for --from-release.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=RAW_DATA_DIR,
        help="Destination directory for downloaded CSV files (default: data/raw).",
    )
    parser.add_argument(
        "--source",
        default=None,
        help="Ingest source for direct download ('pinned' or 'live').",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force re-download even when cached files exist.",
    )
    return parser


def main(argv: list[str] | tuple[str, ...] | None = None) -> int:
    """CLI entry point for ``python -m mtp.download``."""
    parser = build_parser()
    args = parser.parse_args(list(argv) if argv is not None else None)
    if args.from_release is not None:
        download_from_release(
            tag=args.from_release,
            output_dir=args.output_dir,
            repo=args.repo,
        )
    else:
        download_exports(
            raw_dir=args.output_dir,
            force=bool(args.force),
            source=args.source,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
