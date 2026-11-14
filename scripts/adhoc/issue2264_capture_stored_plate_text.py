#!/usr/bin/env python3
"""Record the text drawn on the stored 2026-10-22 renders of NCASM plates I-XX.

Purpose: a reference for ``tests/test_ncasm_tufte_plates.py``, which checks that the
ported builders still draw this text. It is a record of what the plates said, not
evidence that any number on them is right (#776-#779).
Inputs: the folder holding the stored ``plate_01_*.svg`` ... ``plate_20_*.svg`` (the
slide-wizard deck folder's ``tufte_plates/``). Those SVGs keep text as text.
Outputs: a JSON file mapping each plate stem to its text nodes in drawing order, with
the SHA-256 of the SVG it was read from.
Issue #2264; 2026-11-11.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

from defusedxml import ElementTree


def text_nodes(svg: Path) -> list[str]:
    """The text of every ``<text>`` element of ``svg``, in document order."""
    root = ElementTree.parse(svg).getroot()
    return ["".join(node.itertext()) for node in root.iter("{http://www.w3.org/2000/svg}text")]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("stored_dir", type=Path, help="folder with the stored plate SVGs")
    parser.add_argument("output", type=Path, help="JSON file to write (must not exist)")
    args = parser.parse_args(argv)

    if args.output.exists():
        print(f"{args.output} exists; refusing to replace it", file=sys.stderr)
        return 1
    record: dict[str, dict[str, object]] = {}
    for number in range(1, 21):
        matches = sorted(args.stored_dir.glob(f"plate_{number:02d}_*.svg"))
        if len(matches) != 1:
            print(f"plate {number}: expected one SVG, found {len(matches)}", file=sys.stderr)
            return 1
        svg = matches[0]
        nodes = text_nodes(svg)
        if not nodes:
            print(f"{svg.name}: no text nodes; was it saved with text as paths?", file=sys.stderr)
            return 1
        record[svg.stem] = {
            "svg_sha256": hashlib.sha256(svg.read_bytes()).hexdigest(),
            "text": nodes,
        }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(record, indent=1, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(f"{len(record)} plates, {sum(len(v['text']) for v in record.values())} text nodes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
