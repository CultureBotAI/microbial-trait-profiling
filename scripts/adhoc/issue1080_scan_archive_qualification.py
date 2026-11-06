"""Prepare and compare synthetic v2 scan relocation evidence (#1080, 2026-10-28).

Inputs: the committed selection-scan fixture builder and the real recipe outputs
made by the companion shell driver. Outputs: synthetic inputs, bound producer
inventory and strict parity JSON under a fresh durable qualification directory.
No historical F inputs are read. Only the report generation timestamp may differ.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

from git import Repo

from mtp.paths import REPO_ROOT, file_sha256
from mtp.workflow.output_roots import require_durable_output

FIXTURE = "tests/test_cellulose_decr_selection_scan.py"
BOUND = (
    FIXTURE,
    "scripts/adhoc/issue1080_scan_archive_qualification.py",
    "scripts/adhoc/issue1080_scan_archive_qualification.sh",
    "src/mtp/cellulose_decr/selection_scan.py",
    "src/mtp/cellulose_decr/selection_scan_scopes.py",
    "methods/cellulose_decr_selection_scan_protocol_v2.json",
    "methods/cellulose_decr_selection_scan.md",
    "Justfile",
    "uv.lock",
)
SCOPES = ("phenotype", "off")


def producer() -> dict[str, Any]:
    repository = Repo(REPO_ROOT)
    if repository.is_dirty(untracked_files=True):
        raise ValueError("Qualification requires a clean committed producer")
    return {
        "commit": repository.head.commit.hexsha,
        "files": {name: file_sha256(REPO_ROOT / name) for name in BOUND},
    }


def write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def prepare(root: Path) -> None:
    identity = producer()
    root.mkdir(parents=True, exist_ok=False)
    inputs = root / "original" / "inputs"
    inputs.mkdir(parents=True)
    fixture_path = REPO_ROOT / FIXTURE
    spec = importlib.util.spec_from_file_location("issue1080_fixture", fixture_path)
    if spec is None or spec.loader is None:
        raise ValueError("Cannot load the committed synthetic fixture builder")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    bindings = module._write_core_inputs(inputs)
    identity["synthetic_inputs"] = {
        name: {"path": str(path), "sha256": file_sha256(path)} for name, path in bindings.items()
    }
    write_json(root / "producer.json", identity)


def single_manifest(root: Path) -> Path:
    manifests = list(root.glob("run-*/run.json"))
    if len(manifests) != 1:
        raise ValueError(f"Expected exactly one core manifest in {root}")
    return manifests[0]


def compare(root: Path) -> None:
    identity = json.loads((root / "producer.json").read_bytes())
    current = producer()
    if {key: identity[key] for key in current} != current:
        raise ValueError("Qualification producer changed during execution")
    if (root / "original").exists() or (root / "archives").exists():
        raise ValueError("Original logical inputs and original archives must be unavailable")
    before = root / "unmoved"
    after = root / "relocated"
    table = "selection_scan_scopes.tsv"
    report = "selection_scan_scopes.json"
    if (before / table).read_bytes() != (after / table).read_bytes():
        raise ValueError("Relocation changed the paired trait table")
    reports = [json.loads((directory / report).read_bytes()) for directory in (before, after)]
    timestamps = [value["provenance"].pop("generated_at_utc") for value in reports]
    if reports[0] != reports[1]:
        raise ValueError("Relocation changed report content beyond generation time")
    paired = reports[1]
    if paired["protocol"] != "cellulose-decr-selection-scan-v2":
        raise ValueError("The qualification must exercise protocol v2")
    if paired["follows_protocol"] is not False:
        raise ValueError("The small synthetic null is not the production protocol run")
    if [(item["scope"], item["role"]) for item in paired["scopes"]] != [
        ("phenotype", "primary"),
        ("off", "declared_sensitivity"),
    ]:
        raise ValueError("Scope roles or ordering changed")
    records = {}
    for scope, summary in zip(SCOPES, paired["scopes"], strict=True):
        counts = summary["counts"]
        if (
            counts["eligible_conditions"],
            counts["triples"],
            counts["traits"],
            counts["pool_studies"],
        ) != (5, 10, 2, 18):
            raise ValueError(f"Unexpected synthetic scan population for {scope}: {counts}")
        archived = root / "restored" / scope / "archive.json"
        resolved = Path((root / f"resolved-{scope}.txt").read_text().strip())
        if not resolved.is_relative_to(root / "restored" / scope / "data"):
            raise ValueError("Evidence resolution escaped the restored archive")
        source = root / "retired-original" / f"scan-{scope}" / "core_consumption.json"
        if source.read_bytes() != resolved.read_bytes():
            raise ValueError("Archive restore changed the consumer record")
        manifest = single_manifest(root / "retired-original" / f"core-{scope}")
        run = json.loads(manifest.read_bytes())
        if run["git_commit"] != current["commit"] or run["git_dirty"] is not False:
            raise ValueError("Synthetic core has the wrong producer identity")
        records[scope] = {
            "archive_sha256": file_sha256(archived),
            "core_manifest_sha256": file_sha256(manifest),
            "core_run_id": run["run_id"],
            "consumer_sha256": file_sha256(resolved),
            "resolved_consumer": str(resolved),
            "counts": counts,
        }
    write_json(
        root / "parity.json",
        {
            "passed": True,
            "synthetic_only": True,
            "follows_production_permutation_protocol": False,
            "producer": identity,
            "scopes": records,
            "table_sha256": file_sha256(before / table),
            "report_sha256": {
                "unmoved": file_sha256(before / report),
                "relocated": file_sha256(after / report),
            },
            "allowed_difference": "provenance.generated_at_utc only",
            "generation_timestamps": timestamps,
        },
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("operation", choices=("prepare", "compare"))
    parser.add_argument("root", type=Path)
    args = parser.parse_args()
    if not args.root.is_absolute():
        parser.error("Qualification requires an absolute output root")
    require_durable_output(parser, args.root, allowed=False)
    root = args.root.expanduser().resolve()
    if root.is_relative_to(REPO_ROOT):
        parser.error("Qualification evidence must be outside the working tree")
    (prepare if args.operation == "prepare" else compare)(root)


if __name__ == "__main__":
    main()
