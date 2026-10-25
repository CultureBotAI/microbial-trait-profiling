"""Exercise synthetic publication from the actual CLI entry point (#1088; 2026-10-25).

Inputs: --output-root, inherited umask, and the repository's tiny source fixture.
Outputs: a real three-test table bundle and result.json under the fresh output root.
The fixture uses synthetic source adapters and stubs execution verification;
native calculations and CLI publication run normally. A multiprocessing re-import
fails with a deterministic marker before heavy imports, so the regression never
depends on host load or a timing threshold.
"""

from __future__ import annotations

import argparse
import csv
import importlib
import json
import os
import stat
import sys
from pathlib import Path

if __name__ == "__mp_main__":
    raise RuntimeError("ISSUE1088_FORBIDDEN_MAIN_REIMPORT")

# Import the same heavy CLI/stage graph in the real __main__ process. A spawned
# permission worker would re-execute this module under __mp_main__ and fail above.
cli = importlib.import_module("mtp.three_test_table.__main__")
stage = importlib.import_module("mtp.three_test_table.stage")
pytest = importlib.import_module("pytest")
ROOT = Path(__file__).resolve().parents[2]


def _forbid_parent_umask(mask: int) -> int:
    del mask
    raise AssertionError("Publication must not query or change the parent process umask")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    root = stage.require_durable(args.output_root)
    root.mkdir(parents=True, exist_ok=False)
    sys.path.insert(0, str(ROOT / "tests"))
    fixture = importlib.import_module("test_three_test_stage")
    parent = root / "published"
    parent.mkdir()
    previous_mode = stat.S_IMODE(parent.stat().st_mode)
    output = parent / "output"
    with pytest.MonkeyPatch.context() as monkeypatch:
        sources = fixture.sources_for_stage(root, monkeypatch)
        monkeypatch.setattr(os, "umask", _forbid_parent_umask)
        exit_code = cli.main(
            [
                "run",
                "--comparison-run",
                str(sources.comparison.manifest_path),
                "--output",
                str(output),
                "--scratch-dir",
                str(root / "scratch"),
            ]
        )
        with (output / "per_trait_evidence.tsv").open(encoding="utf-8", newline="") as handle:
            records = list(csv.DictReader(handle, delimiter="\t"))
        report = {
            "cli_exit_code": exit_code,
            "output_mode": stat.S_IMODE(output.stat().st_mode),
            "parent_mode_before": previous_mode,
            "parent_mode_after": stat.S_IMODE(parent.stat().st_mode),
            "siblings": sorted(path.name for path in parent.iterdir()),
            "trait_count": len(records),
            "expected_trait_count": len(sources.features["count_profiles"]),
            "unique_traits": len({row["trait"] for row in records}),
            "family_sizes": sorted({row["composition_family_size"] for row in records}),
        }
    (root / "result.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
