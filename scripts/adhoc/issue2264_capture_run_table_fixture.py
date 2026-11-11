#!/usr/bin/env python3
"""Copy the run tables that NCASM plate 21 reads into a test fixture.

Purpose: let ``tests/test_ncasm_tufte_plates.py`` draw the table-driven plate from real
run output without the run directory, and pin what that output said.
Inputs: RUN_ROOT, a run root holding ``comparison/`` and ``arms/`` (the fixture in the
repository was taken from ``F-511c891/runs/tier12-off``).
Outputs: under OUTPUT_DIR, the two per-trait tables copied byte for byte at their
run-relative paths, ``comparison/comparison_run.json`` reduced to the seven fields the
loader reads (the full record names local absolute paths), and ``README.md`` listing
every file with the SHA-256 of its source.
Issue #2264; 2026-11-11.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
from pathlib import Path

VERBATIM = (
    "comparison/arm1_pooled_stats/trait_profile_paired_contrast_trait_stats_pooled.tsv",
    "arms/noel/trait_direction_screen/trait_direction_screen.tsv",
)
RUN_RECORD = "comparison/comparison_run.json"
RUN_FIELDS = (
    "id",
    "contrast_population",
    "direction_basis",
    "status",
    "core_run_id",
    "git_commit",
    "git_dirty",
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("run_root", type=Path)
    parser.add_argument("output_dir", type=Path, help="fixture folder to create (must not exist)")
    args = parser.parse_args(argv)
    if args.output_dir.exists():
        print(f"{args.output_dir} exists; refusing to replace it", file=sys.stderr)
        return 1

    lines = [
        "# Run tables for the NCASM plate tests",
        "",
        "Copied by `scripts/adhoc/issue2264_capture_run_table_fixture.py` from one analysis",
        "run (see `comparison/comparison_run.json` for its population and code revision).",
        "They are real run output, kept so the table-driven plates can be drawn and checked",
        "without the run directory. Do not edit them; capture a new fixture instead.",
        "",
        "| File | Source SHA-256 | Copy |",
        "| --- | --- | --- |",
    ]
    for relative in VERBATIM:
        source = args.run_root / relative
        target = args.output_dir / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        digest = _sha256(source)
        if _sha256(target) != digest:
            print(f"{relative}: the copy differs from its source", file=sys.stderr)
            return 1
        lines.append(f"| `{relative}` | `{digest}` | byte for byte |")

    source = args.run_root / RUN_RECORD
    record = json.loads(source.read_text(encoding="utf-8"))
    reduced = {field: record[field] for field in RUN_FIELDS}
    target = args.output_dir / RUN_RECORD
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(reduced, indent=2) + "\n", encoding="utf-8")
    lines.append(
        f"| `{RUN_RECORD}` | `{_sha256(source)}` | reduced to {len(RUN_FIELDS)} fields: "
        f"{', '.join(RUN_FIELDS)} |"
    )
    (args.output_dir / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {len(VERBATIM) + 2} files under {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
