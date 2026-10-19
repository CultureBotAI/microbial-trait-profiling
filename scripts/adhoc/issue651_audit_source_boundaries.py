#!/usr/bin/env python3
"""Audit source reconciliation on an existing bundle (#651, 2026-10-19).

Input: a core run manifest and its available raw/Study/export artifact paths.
Output: JSON reconciliation counts and elapsed seconds. This checks these two
scientific boundaries only; it does not assert full workflow-manifest validity.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from time import perf_counter

from mtp.workflow.models import WorkflowRun
from mtp.workflow.runtime import artifact_path
from mtp.workflow.source_contracts import validate_export_sources, validate_study_sources


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    args = parser.parse_args()
    run = WorkflowRun.model_validate_json(args.manifest.read_text(encoding="utf-8"))
    run.manifest_path = str(args.manifest.resolve())
    paths = {item.product_id: artifact_path(run, item.product_id) for item in run.artifacts or []}
    started = perf_counter()
    ingest = validate_study_sources(paths["raw_dump"], paths["studies"], paths["ingest_log"])
    ingest_seconds = perf_counter() - started
    started = perf_counter()
    exports = validate_export_sources(paths["studies"], paths)
    print(
        json.dumps(
            {
                "run_id": run.run_id,
                "scope": "raw-to-Study and Study-to-export reconciliation only",
                "ingest": {**ingest, "elapsed_seconds": round(ingest_seconds, 3)},
                "exports": {**exports, "elapsed_seconds": round(perf_counter() - started, 3)},
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
