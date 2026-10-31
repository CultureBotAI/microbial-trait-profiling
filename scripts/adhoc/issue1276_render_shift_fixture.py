"""Export a synthetic six-trait layout fixture for #1276/#1280 (2026-10-30).

Inputs: fixed, explicitly synthetic display values and the current figure builder.
Outputs: PNG/PDF/SVG, Plotly JSON and an export/hash record in a fresh durable directory.
No research data is read and no statistical inference or scientific validation runs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import asdict
from pathlib import Path

from git import Repo

from mtp import table1_trait_maps
from mtp.paths import REPO_ROOT
from mtp.workflow.output_roots import require_durable_output_path


def fixture_result() -> table1_trait_maps.ShiftResult:
    """Fixed plotting inputs; probability/interval values are illustrative only."""
    traits = table1_trait_maps.TRAIT_ORDER
    return table1_trait_maps.ShiftResult(
        shifts=[
            table1_trait_maps.TraitShift(
                trait=trait,
                studies_up=500,
                studies_down=400,
                studies_tied=22,
                share_up=500 / 900,
                ci_low=0.52,
                ci_high=0.59,
                p_value=0.001,
                q_value=0.002,
            )
            for trait in traits
        ],
        spearman={
            f"{left}~{right}": (0.3, 0.01, 900)
            for index, left in enumerate(traits)
            for right in traits[index + 1 :]
        },
        studies=922,
    )


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    output = require_durable_output_path(args.output_dir.expanduser(), allowed=False)
    if output.is_relative_to(REPO_ROOT.resolve()):
        parser.error("Fixture exports must be outside the source repository")
    if args.output_dir.is_symlink() or output.exists():
        parser.error("Fixture export directory must be fresh and nonsymlink")
    sources = (Path(__file__).resolve(), Path(table1_trait_maps.__file__).resolve())
    source_hashes = {str(path.relative_to(REPO_ROOT)): sha256(path) for path in sources}
    result = fixture_result()
    figure = table1_trait_maps.build_shift_figure(
        result, "synthetic fixture, 2026-10-30", table1_trait_maps.TRAIT_ORDER
    )
    figure.layout.title.text = "Synthetic fixture: " + str(figure.layout.title.text)
    output.mkdir(parents=True, exist_ok=False)
    for extension in ("png", "pdf", "svg"):
        figure.write_image(str(output / f"shift-fixture.{extension}"), scale=2)
    figure.write_json(output / "shift-fixture.plotly.json")
    exports = {}
    for path in sorted(output.iterdir()):
        size = path.stat().st_size
        if not size:
            raise ValueError(f"Empty fixture export: {path.name}")
        exports[path.name] = {"sha256": sha256(path), "size_bytes": size}
    if source_hashes != {str(path.relative_to(REPO_ROOT)): sha256(path) for path in sources}:
        raise ValueError("Fixture producer source changed during rendering")
    repository = Repo(REPO_ROOT)
    record = {
        "status": "exported_synthetic_fixture_only",
        "warning": "Illustrative values only; not research results or inferential evidence.",
        "code_commit": repository.head.commit.hexsha,
        "code_dirty": repository.is_dirty(untracked_files=True),
        "source_sha256": source_hashes,
        "fixture": asdict(result),
        "width_layout_px": figure.layout.width,
        "height_layout_px": figure.layout.height,
        "export_scale": 2,
        "exports": exports,
    }
    with (output / "fixture-export.json").open("x", encoding="utf-8") as handle:
        json.dump(record, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(output / "fixture-export.json")


if __name__ == "__main__":
    main()
