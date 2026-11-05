"""Inspect static correction anchors for #851 (2026-11-04).

Input: repository Python sources and the focused registry guard.
Output: JSON discovery evidence at an explicitly supplied path; no analysis runs.
"""

from __future__ import annotations

import argparse
import json
import runpy
from collections.abc import Iterator
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    guard = runpy.run_path(str(root / "tests/test_multiple_testing_registry.py"))
    source_count = 0

    def counted_sources() -> Iterator[tuple[str, str]]:
        nonlocal source_count
        for path, source in guard["repository_sources"](root):
            source_count += 1
            yield path, source

    sites, _ = guard["discover"](counted_sources())
    records = [
        {
            "path": path,
            "function": function,
            "procedures": sorted(site["procedures"]),
            "calls": site["calls"],
            "definition": site["definition"],
        }
        for (path, function), site in sorted(sites.items())
    ]
    args.output.write_text(
        json.dumps({"source_files": source_count, "anchors": records}, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"{source_count} source files; {len(records)} correction anchors; {args.output}")


if __name__ == "__main__":
    main()
