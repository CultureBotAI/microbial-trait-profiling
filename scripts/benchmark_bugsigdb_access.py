#!/usr/bin/env python3
"""Benchmark the non-R BugSigDB data access methods.

Measures four access paths (raw Help:Export CSVs, BugSigDBExports full_dump,
BugSigDBExports GMT bundle, Zenodo pinned release). For each, records the
file count, declared total size, downloaded bytes, total elapsed seconds,
and throughput.

Writes one TSV per method plus a single Markdown summary to
``docs/bugsigdb-access-benchmark.md``.

Usage:
    uv run python scripts/benchmark_bugsigdb_access.py
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from datetime import UTC, datetime

import httpx

from mtp.paths import REPO_ROOT

BENCHMARK_DIR = REPO_ROOT / "data" / "benchmarks"
REPORT_PATH = REPO_ROOT / "docs" / "bugsigdb-access-benchmark.md"
USER_AGENT = "bsdbng-benchmark/0.1"
BUGSIGDB_EXPORTS_VERSION = "v1.3.1"
BUGSIGDB_EXPORTS_DOI = "10.5281/zenodo.19740674"
ZENODO_RECORD_ID = "19740674"

RAW_CSV_URLS = [
    "https://bugsigdb-csv.s3.us-east-va.perf.cloud.ovh.us/studies.csv",
    "https://bugsigdb-csv.s3.us-east-va.perf.cloud.ovh.us/experiments.csv",
    "https://bugsigdb-csv.s3.us-east-va.perf.cloud.ovh.us/signatures.csv",
]
FULL_DUMP_URL = (
    "https://raw.githubusercontent.com/waldronlab/BugSigDBExports/"
    f"{BUGSIGDB_EXPORTS_VERSION}/full_dump.csv"
)
BUGSIGDB_EXPORTS_CONTENTS_API = (
    "https://api.github.com/repos/waldronlab/BugSigDBExports/contents"
    f"?ref={BUGSIGDB_EXPORTS_VERSION}"
)
ZENODO_API = f"https://zenodo.org/api/records/{ZENODO_RECORD_ID}"


@dataclass
class MethodResult:
    name: str
    files: int
    declared_bytes: int
    downloaded_bytes: int
    seconds: float

    @property
    def throughput_mib_per_s(self) -> float:
        if self.seconds <= 0:
            return 0.0
        return self.downloaded_bytes / self.seconds / (1024 * 1024)


def _download_bytes(url: str, client: httpx.Client) -> tuple[int, float]:
    """Stream a URL; return (bytes, seconds) without persisting to disk."""
    start = time.perf_counter()
    size = 0
    with client.stream("GET", url, follow_redirects=True) as response:
        response.raise_for_status()
        for chunk in response.iter_bytes():
            size += len(chunk)
    return size, time.perf_counter() - start


def _list_gh_contents(client: httpx.Client, api_url: str) -> list[dict]:
    """Fetch a GitHub contents listing, following Link header pagination."""
    entries: list[dict] = []
    separator = "&" if "?" in api_url else "?"
    url: str | None = f"{api_url}{separator}per_page=100"
    while url:
        response = client.get(url, follow_redirects=True)
        response.raise_for_status()
        entries.extend(response.json())
        next_link = response.links.get("next")
        url = next_link["url"] if next_link else None
    return entries


def _list_zenodo_files(client: httpx.Client, api_url: str) -> list[dict]:
    response = client.get(api_url, follow_redirects=True)
    response.raise_for_status()
    return response.json()["files"]


def benchmark_raw_csv(client: httpx.Client) -> MethodResult:
    start = time.perf_counter()
    total_bytes = 0
    for url in RAW_CSV_URLS:
        size, _seconds = _download_bytes(url, client)
        total_bytes += size
    elapsed = time.perf_counter() - start
    return MethodResult(
        name="raw_csv_exports",
        files=len(RAW_CSV_URLS),
        declared_bytes=total_bytes,
        downloaded_bytes=total_bytes,
        seconds=elapsed,
    )


def benchmark_full_dump(client: httpx.Client) -> MethodResult:
    start = time.perf_counter()
    size, _seconds = _download_bytes(FULL_DUMP_URL, client)
    return MethodResult(
        name=f"bugsigdbexports_{BUGSIGDB_EXPORTS_VERSION.replace('.', '_')}_full_dump",
        files=1,
        declared_bytes=size,
        downloaded_bytes=size,
        seconds=time.perf_counter() - start,
    )


def benchmark_gmt_bundle(client: httpx.Client) -> MethodResult:
    start = time.perf_counter()
    entries = [
        e
        for e in _list_gh_contents(client, BUGSIGDB_EXPORTS_CONTENTS_API)
        if e.get("type") == "file" and e.get("name", "").endswith(".gmt")
    ]
    declared_total = sum(int(e.get("size", 0)) for e in entries)
    total_bytes = 0
    for entry in entries:
        size, _seconds = _download_bytes(entry["download_url"], client)
        total_bytes += size
    return MethodResult(
        name=f"bugsigdbexports_{BUGSIGDB_EXPORTS_VERSION.replace('.', '_')}_gmt_bundle",
        files=len(entries),
        declared_bytes=declared_total,
        downloaded_bytes=total_bytes,
        seconds=time.perf_counter() - start,
    )


def benchmark_zenodo_release(client: httpx.Client) -> MethodResult:
    start = time.perf_counter()
    files = _list_zenodo_files(client, ZENODO_API)
    declared_total = sum(int(f["size"]) for f in files)
    total_bytes = 0
    for f in files:
        size, _seconds = _download_bytes(f["links"]["self"], client)
        total_bytes += size
    return MethodResult(
        name=f"zenodo_release_{BUGSIGDB_EXPORTS_VERSION.replace('.', '_')}",
        files=len(files),
        declared_bytes=declared_total,
        downloaded_bytes=total_bytes,
        seconds=time.perf_counter() - start,
    )


def _write_tsv(result: MethodResult) -> None:
    row = "\t".join(
        [
            result.name,
            str(result.files),
            str(result.declared_bytes),
            str(result.downloaded_bytes),
            f"{result.seconds:.6f}",
        ]
    )
    (BENCHMARK_DIR / f"{result.name}.tsv").write_text(row + "\n", encoding="utf-8")


def _fmt_mib(size: int) -> str:
    return f"{size / (1024 * 1024):.2f} MiB"


def _write_report(results: list[MethodResult]) -> None:
    timestamp = datetime.now(UTC).replace(microsecond=0).isoformat()
    lines = [
        "# BugSigDB Access Benchmark",
        "",
        "This file records a single-run benchmark of the non-R BugSigDB data access methods.",
        "",
        f"Generated at: `{timestamp}`",
        "",
        "## Methodology",
        "",
        "- one end-to-end run per access method",
        "- files downloaded sequentially",
        "- bytes read fully but not persisted to disk",
        "- timings include manifest resolution for methods that require API discovery",
        "- host environment: local workstation run from `bsdbng`",
        "",
        "## Results",
        "",
        "| Method | Files | Declared size | Downloaded size | Total seconds | Throughput |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for r in results:
        lines.append(
            f"| {r.name} | {r.files} | {_fmt_mib(r.declared_bytes)} | "
            f"{_fmt_mib(r.downloaded_bytes)} | {r.seconds:.2f} | "
            f"{r.throughput_mib_per_s:.2f} MiB/s |"
        )
    lines.extend(
        [
            "",
            "## Notes",
            "",
            "- `raw_csv_exports` is the live source-faithful path and remains "
            "the preferred surface once separate-table ingestion is implemented.",
            f"- `bugsigdbexports_{BUGSIGDB_EXPORTS_VERSION.replace('.', '_')}_full_dump` "
            "is the citable convenience single-file path currently used by bsdbng.",
            f"- `bugsigdbexports_{BUGSIGDB_EXPORTS_VERSION.replace('.', '_')}_gmt_bundle` "
            "is the heaviest path and is optimized "
            "for set-based analysis, not canonical ingestion.",
            f"- `zenodo_release_{BUGSIGDB_EXPORTS_VERSION.replace('.', '_')}` "
            f"is the reproducible pinned-release path ({BUGSIGDB_EXPORTS_DOI}).",
            "",
        ]
    )
    REPORT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    BENCHMARK_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)

    results: list[MethodResult] = []
    with httpx.Client(headers={"User-Agent": USER_AGENT}, timeout=httpx.Timeout(60.0)) as client:
        for label, fn in [
            ("raw-csv", benchmark_raw_csv),
            ("full-dump", benchmark_full_dump),
            ("gmt-bundle", benchmark_gmt_bundle),
            ("zenodo-release", benchmark_zenodo_release),
        ]:
            print(f"benchmarking {label}...", flush=True)
            result = fn(client)
            _write_tsv(result)
            results.append(result)
            print(
                f"  {result.name}: {result.files} files, "
                f"{_fmt_mib(result.downloaded_bytes)} in {result.seconds:.2f}s "
                f"({result.throughput_mib_per_s:.2f} MiB/s)"
            )

    _write_report(results)
    print(f"\nWrote {REPORT_PATH}")


if __name__ == "__main__":
    main()
