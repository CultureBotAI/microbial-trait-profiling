#!/usr/bin/env python3
"""Capture the reviewed historical audit-number fixture for #1709, 2026-11-06.

Inputs: the Git source blob at retirement commit 92539eb and the committed numeric
extractors in tests/test_cellulose_historical_numeric_snapshot.py.
Output: exclusively create the requested JSON fixture; never update an existing
snapshot. Construct six artist trees and three SVG/16-slide strings; do not save
or rasterize figures, import analysis outputs, or regenerate scientific data.
Run only in the scheduled serial lane. The default source is deliberately fixed:
changing current producer code cannot silently rebaseline historical numbers.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from pathlib import Path
from tempfile import TemporaryDirectory
from types import ModuleType
from typing import cast

import matplotlib.pyplot as plt
from git import Repo
from matplotlib.figure import Figure

ROOT = Path(__file__).resolve().parents[2]
SOURCE = "src/mtp/cellulose_decr/tufte_plates.py"
REVIEWED_COMMIT = "92539eb"
EXTRACTORS = ROOT / "tests/test_cellulose_historical_numeric_snapshot.py"
DEFAULT_OUTPUT = ROOT / "tests/fixtures/cellulose_historical_numeric_snapshot.json"
MAX_SOURCE_BYTES = 512 * 1024
Snapshot = dict[str, dict[str, list[str]]]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def load_extractors() -> ModuleType:
    """Load reusable test helpers without invoking fixtures or running tests."""
    spec = importlib.util.spec_from_file_location("_cellulose_numeric_extractors", EXTRACTORS)
    if spec is None or spec.loader is None:
        raise ValueError("Numeric extractor module specification or loader is missing")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    require(
        Path(module.tufte_plates.__file__).resolve() == (ROOT / SOURCE).resolve(),
        "Set PYTHONPATH to this checkout's src directory; imported producer differs",
    )
    return module


@contextmanager
def reviewed_module() -> Iterator[ModuleType]:
    """Import a bounded copy of the immutable Git blob in configured temporary storage."""
    repository = Repo(ROOT)
    try:
        blob = repository.commit(REVIEWED_COMMIT).tree / SOURCE
        payload = blob.data_stream.read(MAX_SOURCE_BYTES + 1)
    finally:
        repository.close()
    require(0 < len(payload) <= MAX_SOURCE_BYTES, "Historical source exceeds bounded size")
    with TemporaryDirectory(prefix="cellulose-numeric-") as temporary:
        filename = Path(temporary) / "tufte_plates.py"
        filename.write_bytes(payload)
        spec = importlib.util.spec_from_file_location("_cellulose_historical_92539eb", filename)
        if spec is None or spec.loader is None:
            raise ValueError("Historical module specification or loader is missing")
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        try:
            spec.loader.exec_module(module)
            yield module
        finally:
            sys.modules.pop(spec.name, None)


def capture(module: ModuleType, extractors: ModuleType) -> Snapshot:
    figure_tokens = cast(Callable[[Figure], list[str]], extractors.figure_numeric_tokens)
    svg_tokens = cast(Callable[[str], list[str]], extractors.svg_numeric_tokens)
    slide_tokens = cast(Callable[[str], dict[str, list[str]]], extractors.html_numeric_tokens)
    builders: dict[str, list[str]] = {}
    for name in extractors.BUILDERS:
        builder = cast(Callable[[], Figure], getattr(module, name))
        figure = builder()
        try:
            builders[name] = figure_tokens(figure)
        finally:
            plt.close(figure)
    svgs = {
        name: svg_tokens(cast(Callable[[], str], getattr(module, name))())
        for name in extractors.SVG_RENDERERS
    }
    slides = slide_tokens(module.render_slide_deck_html())
    require(len(builders) == 6, "Expected six historical builders")
    require(len(svgs) == 3, "Expected three standalone historical SVGs")
    require(set(slides) == set(extractors.SLIDES), "Expected all sixteen historical slides")
    snapshot = {"builders": builders, "svgs": svgs, "slides": slides}
    for group, outputs in snapshot.items():
        require(all(outputs.values()), f"Empty numeric tokens in {group}")
    return snapshot


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    require(not args.output.exists(), "Snapshot already exists; automatic rebaselining is refused")
    require(not args.output.is_symlink(), "Aliased snapshot output is refused")
    extractors = load_extractors()
    with reviewed_module() as baseline:
        snapshot = capture(baseline, extractors)
    payload = {
        "schema_version": 1,
        "baseline_commit": REVIEWED_COMMIT,
        "source": SOURCE,
        "purpose": (
            "Retain the reviewed historical audit display; these unverified numbers "
            "are not current analysis results or a scientific correctness claim."
        ),
        "extraction": (
            "Ordered visible Figure text/active-axis labels, SVG text nodes, and "
            "each HTML slide's text; Unicode normalization preserves signed numeric tokens."
        ),
        "numeric_tokens": snapshot,
    }
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    written = json.loads(args.output.read_text(encoding="utf-8"))
    require(written == payload, "Written snapshot differs from captured baseline")
    print(
        json.dumps(
            {
                "baseline_commit": REVIEWED_COMMIT,
                "output_bytes": args.output.stat().st_size,
                "outputs_per_group": {group: len(outputs) for group, outputs in snapshot.items()},
                "numeric_tokens_per_group": {
                    group: sum(len(tokens) for tokens in outputs.values())
                    for group, outputs in snapshot.items()
                },
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
