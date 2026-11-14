#!/usr/bin/env python3
"""Compare two folders of plate PNGs pixel by pixel.

Purpose: show whether the port in ``src/mtp/ncasm_tufte_plates/`` draws the same pixels
as the original slide-wizard builders (rendered by
``issue2264_render_original_builders.sh``), and how both differ from the stored
2026-10-22 renders.
Inputs: two directories; PNG files are matched by name.
Outputs: one line per file on stdout (shape, differing pixels, largest channel
difference) and an optional TSV. Exit status 1 when any matched pair differs or a file
has no partner, unless ``--report-only`` is given.
Issue #2264; 2026-11-11.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
from PIL import Image


def _load(path: Path) -> np.ndarray:
    with Image.open(path) as image:
        return np.asarray(image.convert("RGBA"), dtype=np.int16)


def compare(left: Path, right: Path) -> list[dict[str, str]]:
    """Rows describing every PNG name found in either directory."""
    names = sorted({p.name for p in left.glob("*.png")} | {p.name for p in right.glob("*.png")})
    rows: list[dict[str, str]] = []
    for name in names:
        left_path, right_path = left / name, right / name
        if not left_path.is_file() or not right_path.is_file():
            missing = "left" if not left_path.is_file() else "right"
            rows.append({"file": name, "status": f"missing_{missing}"})
            continue
        a, b = _load(left_path), _load(right_path)
        if a.shape != b.shape:
            rows.append(
                {
                    "file": name,
                    "status": "shape_differs",
                    "left_shape": "x".join(map(str, a.shape)),
                    "right_shape": "x".join(map(str, b.shape)),
                }
            )
            continue
        difference = np.abs(a - b)
        differing = int(np.count_nonzero(difference.max(axis=2)))
        rows.append(
            {
                "file": name,
                "status": "identical" if differing == 0 else "differs",
                "left_shape": "x".join(map(str, a.shape)),
                "right_shape": "x".join(map(str, b.shape)),
                "differing_pixels": str(differing),
                "pixels": str(a.shape[0] * a.shape[1]),
                "max_channel_difference": str(int(difference.max())),
            }
        )
    return rows


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("left", type=Path)
    parser.add_argument("right", type=Path)
    parser.add_argument("--tsv", type=Path, help="also write the rows to this TSV")
    parser.add_argument("--report-only", action="store_true", help="always exit 0")
    args = parser.parse_args(argv)

    rows = compare(args.left, args.right)
    if not rows:
        print("no PNG files found", file=sys.stderr)
        return 1
    fields = [
        "file",
        "status",
        "left_shape",
        "right_shape",
        "differing_pixels",
        "pixels",
        "max_channel_difference",
    ]
    for row in rows:
        print("\t".join(row.get(field, "") for field in fields))
    if args.tsv is not None:
        with args.tsv.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
    identical = sum(row["status"] == "identical" for row in rows)
    print(f"{identical} of {len(rows)} PNG files identical", file=sys.stderr)
    if args.report_only:
        return 0
    return 0 if identical == len(rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
