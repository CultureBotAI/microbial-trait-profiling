#!/usr/bin/env python3
"""Download the merged KG-Microbe nodes and edges TSVs.

Fetches two files from a pinned public Google Drive folder (Marcin
Joachimiak's published KG-Microbe build) into ``data/kg-microbe/``:

- ``merged-kg_nodes.tsv``
- ``merged-kg_edges.tsv``

Also writes ``data/kg-microbe/provenance.json`` recording file IDs, sizes,
md5 hashes, and the UTC fetch time, so it's obvious when the upstream
files have been refreshed.

The file IDs are defined as module-level constants. Swap them when a newer
build is published to a different Drive folder.

Usage:
    uv run --extra kg python scripts/download_kg_microbe.py
    uv run --extra kg python scripts/download_kg_microbe.py --force

The ``kg`` extra provides ``gdown`` (pure Python, no system dependencies).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from dataclasses import dataclass
from datetime import UTC, datetime

import gdown

from mtp.paths import KG_MICROBE_DIR

OUT_DIR = KG_MICROBE_DIR
PROVENANCE_PATH = OUT_DIR / "provenance.json"


@dataclass(frozen=True)
class DriveFile:
    name: str
    file_id: str


# Drive file IDs from Marcin Joachimiak's published 20260120 build folder.
# Not secrets; the detect-secrets pragma silences a base64-pattern false positive.
_NODES_ID = "1lDVjE67E82LMUdNtqWGvyczrHvNu3BdY"  # pragma: allowlist secret
_EDGES_ID = "1iIVd-FEleeRIv7CPS0l1F732p3J9JO5d"  # pragma: allowlist secret

FILES = [
    DriveFile("merged-kg_nodes.tsv", _NODES_ID),
    DriveFile("merged-kg_edges.tsv", _EDGES_ID),
]


def _md5(path) -> str:
    h = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _load_provenance() -> dict:
    if PROVENANCE_PATH.exists():
        return json.loads(PROVENANCE_PATH.read_text(encoding="utf-8"))
    return {"files": {}}


def _should_skip(f: DriveFile, provenance: dict) -> bool:
    dest = OUT_DIR / f.name
    record = provenance.get("files", {}).get(f.name)
    if not dest.exists() or not record:
        return False
    if record.get("file_id") != f.file_id:
        return False
    return record.get("size") == dest.stat().st_size


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-download even if the local files already match the recorded fingerprint",
    )
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    provenance = _load_provenance()
    provenance["files"] = provenance.get("files", {})

    for f in FILES:
        dest = OUT_DIR / f.name

        if not args.force and _should_skip(f, provenance):
            print(f"skip {f.name} (fingerprint matches provenance.json)")
            continue

        print(f"download {f.name} from Drive id {f.file_id}")
        t0 = time.time()
        result = gdown.download(id=f.file_id, output=str(dest), quiet=False)
        if result is None:
            print(
                f"ERROR: gdown failed for {f.name}. The file may have been "
                "un-shared, moved, or need quota approval.",
                file=sys.stderr,
            )
            return 1
        elapsed = time.time() - t0

        size = dest.stat().st_size
        md5 = _md5(dest)
        provenance["files"][f.name] = {
            "file_id": f.file_id,
            "size": size,
            "md5": md5,
            "fetched_at": datetime.now(UTC).replace(microsecond=0).isoformat(),
            "elapsed_seconds": round(elapsed, 2),
        }
        print(f"  {size:,} bytes  md5={md5}  in {elapsed:.1f}s")

    PROVENANCE_PATH.write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"\nWrote {PROVENANCE_PATH}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
