"""Compare #1117 output to a pinned archived numerical reference (2026-10-28).

Inputs: explicit old/new consumer and core pins, a declared trait list and
trial/parent Git SHAs. Outputs: durable parity JSON with same-buffer identities.
The old reference proves archived table bytes, not a clean historical producer.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
from typing import Any

from git import Blob, Commit, Repo

from mtp.trait_lists import parse_trait_list
from mtp.workflow.output_roots import require_durable_output_path

TABLES = (
    "matrix_trait_overlap.tsv",
    "tier12_coverage_by_rank.tsv",
    "filter_cascade.tsv",
    "parkinsons_stool_study_overlap_summary.tsv",
    "parkinsons_stool_study_pairs.tsv",
)
ADDITIONS = (
    "tier12_input_occurrence_rows",
    "tier12_excluded_occurrence_rows",
    "tier12_excluded_profile_keys",
)
PD_ADDITIONS = (
    "input_occurrence_rows",
    "excluded_outside_active_profile_rows",
    "excluded_outside_active_profile_keys",
)
SCRIPTS = (
    "src/mtp/matrix_overlap_coverage.py",
    "src/mtp/trait_lists.py",
    "src/mtp/matrix_overlap_inputs.py",
    "src/mtp/human_host_corpus.py",
    "src/mtp/workflow/coverage.py",
    "src/mtp/workflow/output_roots.py",
)
METHOD = "methods/matrix_coverage_and_pd_overlap.md"


def identity(raw: bytes) -> dict[str, Any]:
    return {"sha256": hashlib.sha256(raw).hexdigest(), "size_bytes": len(raw)}


def table(raw: bytes) -> list[dict[str, str]]:
    reader = csv.DictReader(io.StringIO(raw.decode()), delimiter="\t")
    names = reader.fieldnames
    if not names or len(names) != len(set(names)):
        raise ValueError("Invalid TSV header")
    result = list(reader)
    if any(set(row) != set(names) or None in row.values() for row in result):
        raise ValueError("Malformed TSV row")
    return result


def bundle(
    path: Path, expected: str, snapshots: dict[Path, bytes]
) -> tuple[dict[str, bytes], dict[str, Any]]:
    raw = path.read_bytes()
    if identity(raw)["sha256"] != expected:
        raise ValueError("Consumer manifest pin differs")
    record = json.loads(raw)
    if record["consumer"] != "matrix-overlap-coverage":
        raise ValueError("Wrong reference consumer")
    snapshots[path] = raw
    files: dict[str, bytes] = {}
    inputs: dict[str, Any] = {path.name: {"path": str(path), **identity(raw)}}
    for name in (*TABLES, "summary.json"):
        matches = [entry for entry in record["outputs"] if Path(entry["path"]).name == name]
        if len(matches) != 1:
            raise ValueError(f"Absent or ambiguous consumer output: {name}")
        entry = matches[0]
        local = path.parent / name
        data = local.read_bytes()
        if identity(data) != {key: entry[key] for key in ("sha256", "size_bytes")}:
            raise ValueError(f"Output differs from bound bytes: {name}")
        files[name] = data
        snapshots[local] = data
        inputs[name] = {
            "path": str(local),
            "recorded_path": entry["path"],
            "relocation": "unique_basename_colocated_with_explicit_consumer",
            **identity(data),
        }
    return files, inputs


def without(row: dict[str, Any], names: tuple[str, ...]) -> dict[str, Any]:
    return {key: value for key, value in row.items() if key not in names}


def pinned_core(path: Path, expected: str, snapshots: dict[Path, bytes]) -> dict[str, Any]:
    """Bind core identity and declared scope without replaying scientific verification."""
    raw = path.read_bytes()
    if identity(raw)["sha256"] != expected:
        raise ValueError("Core manifest differs from its explicit pin")
    record: dict[str, Any] = json.loads(raw)
    if record.get("status") != "completed" or record.get("git_dirty") is not False:
        raise ValueError("Pinned core must declare a completed, clean source run")
    parameters = {item["name"]: item["value"] for item in record["parameters"]}
    if len(parameters) != len(record["parameters"]):
        raise ValueError("Pinned core contains duplicate parameters")
    for name, required in {
        "contrast_population": "disease_vs_healthy_tier1_2_subject_controls",
        "direction_basis": "disease_vs_comparison",
    }.items():
        if record.get(name) != required or parameters.get(name, required) != required:
            raise ValueError(f"Pinned core has an incompatible {name}")
    if parameters.get("strain_rollup") not in {"off", "phenotype", "full"}:
        raise ValueError("Pinned core lacks an explicit strain-rollup scope")
    snapshots[path] = raw
    return {
        "core_run_id": record["run_id"],
        "workflow_id": record["definition"]["id"],
        "strain_rollup": parameters["strain_rollup"],
        "contrast_population": record["contrast_population"],
        "core_git_commit_sha": record["git_commit"],
    }


def bind_core_record(record: dict[str, Any], core: dict[str, Any], *, provenance: bool) -> None:
    fields = core if provenance else without(core, ("core_git_commit_sha",))
    if any(record.get(key) != value for key, value in fields.items()):
        raise ValueError("Summary or provenance identity/scope differs from the pinned core")


def git_digest(trial: Commit, relative: str) -> str:
    try:
        blob = trial.tree / relative
    except KeyError as error:
        raise ValueError(
            f"Producer script or method absent from trial Git tree: {relative}"
        ) from error
    if not isinstance(blob, Blob):
        raise ValueError(f"Producer script or method is not a Git blob: {relative}")
    raw: bytes = blob.data_stream.read()
    return hashlib.sha256(raw).hexdigest()


def bind_producer(provenance: dict[str, Any], trial: Commit) -> None:
    scripts = {relative: {"sha256": git_digest(trial, relative)} for relative in SCRIPTS}
    if provenance.get("scripts") != scripts:
        raise ValueError("Producer code differs from the exact script inventory")
    if provenance.get("script_sha256") != scripts[SCRIPTS[0]]["sha256"]:
        raise ValueError("Producer primary script digest alias differs")
    methods = {METHOD: {"sha256": git_digest(trial, METHOD)}}
    if provenance.get("methods") != methods:
        raise ValueError("Producer methods inventory differs from the trial Git blob")
    if provenance.get("methods_sha256") != methods[METHOD]["sha256"]:
        raise ValueError("Producer methods digest alias differs")


def compare(args: argparse.Namespace) -> dict[str, Any]:
    output = require_durable_output_path(args.output_dir, allowed=False)
    old_manifest = args.reference_consumer.resolve(strict=True)
    new_manifest = args.current_consumer.resolve(strict=True)
    core_manifest = args.core_manifest.resolve(strict=True)
    for root in (old_manifest.parent, new_manifest.parent, core_manifest.parent):
        if output.is_relative_to(root) or root.is_relative_to(output):
            raise ValueError("Parity output must be disjoint from consumer and core bundles")
    if output.exists():
        raise ValueError("Parity output must be fresh")
    repo = Repo(Path(__file__).resolve().parents[2])
    trial = repo.commit(args.trial_commit)
    expected_parents = [
        repo.commit(args.issue_commit).hexsha,
        repo.commit(args.compatibility_commit).hexsha,
    ]
    if [parent.hexsha for parent in trial.parents] != expected_parents:
        raise ValueError("Trial does not have the declared ordered parents")
    if repo.head.commit != trial or repo.is_dirty(untracked_files=True):
        raise ValueError("Run this helper from the clean, pinned trial checkout")
    snapshots: dict[Path, bytes] = {}
    core = pinned_core(core_manifest, args.core_sha256, snapshots)
    old, old_inputs = bundle(old_manifest, args.reference_sha256, snapshots)
    new, new_inputs = bundle(new_manifest, args.current_sha256, snapshots)
    for manifest in (old_manifest, new_manifest):
        consumer = json.loads(snapshots[manifest])
        if (
            consumer.get("run_id") != core["core_run_id"]
            or consumer.get("run_manifest_sha256") != args.core_sha256
        ):
            raise ValueError("Consumer does not bind the same explicitly pinned core")
    trait_path = args.trait_list.resolve(strict=True)
    trait_bytes = trait_path.read_bytes()
    traits = parse_trait_list(trait_bytes.decode())
    snapshots[trait_path] = trait_bytes
    summary_old, summary_new = json.loads(old["summary.json"]), json.loads(new["summary.json"])
    for summary in (summary_old, summary_new):
        bind_core_record(summary, core, provenance=False)
    declared = {
        "trait_list_file": args.trait_list.name,
        "trait_list_sha256": identity(trait_bytes)["sha256"],
        "trait_count": len(traits),
        "trait_ids": list(traits),
        "source_name": args.trait_list.name,
        "sha256": identity(trait_bytes)["sha256"],
    }
    pair_count = len(traits) * (len(traits) - 1) // 2
    if len(traits) < 2 or summary_new["declared_traits"] != {**declared, "pair_count": pair_count}:
        raise ValueError("Current declared trait metadata differs from the supplied list")
    old_traits = summary_old["declared_traits"]["trait_ids"]
    if not set(traits) <= set(old_traits):
        raise ValueError("New traits have no archived matrix reference")
    results = []
    for name in TABLES:
        before, after = table(old[name]), table(new[name])
        expected = before
        checked = after
        if name == "matrix_trait_overlap.tsv":
            expected = [
                row for row in before if row["trait_a_id"] in traits and row["trait_b_id"] in traits
            ]
            if len(after) != pair_count:
                raise ValueError("Declared trait pair count differs from the matrix overlap rows")
        elif name == "parkinsons_stool_study_overlap_summary.tsv":
            checked = [without(row, PD_ADDITIONS) for row in after]
        if checked != expected:
            raise ValueError(f"Unexpected numerical, label, key or ordering difference: {name}")
        results.append(
            {
                "table": name,
                "reference_rows": len(before),
                "current_rows": len(after),
                "exact_projected_row_equality": True,
                "byte_identical": old[name] == new[name],
            }
        )
    projected = dict(summary_new)
    projected["declared_traits"] = summary_old["declared_traits"]
    projected["coverage_and_cascade"] = without(summary_new["coverage_and_cascade"], ADDITIONS)
    projected["parkinsons_stool_overlap"] = [
        without(row, PD_ADDITIONS) for row in summary_new["parkinsons_stool_overlap"]
    ]
    if projected != summary_old:
        raise ValueError("Unexpected summary difference after declared list/exclusion additions")
    provenance_path = new_manifest.parent / "provenance.json"
    provenance_bytes = provenance_path.read_bytes()
    provenance = json.loads(provenance_bytes)
    consumer = json.loads(snapshots[new_manifest])
    matching = [
        entry for entry in consumer["outputs"] if Path(entry["path"]).name == "provenance.json"
    ]
    if len(matching) != 1 or identity(provenance_bytes) != {
        key: matching[0][key] for key in ("sha256", "size_bytes")
    }:
        raise ValueError("Current provenance not bound by consumer")
    if provenance["code_commit"] != trial.hexsha:
        raise ValueError("Current sidecar does not name the clean trial producer")
    snapshots[provenance_path] = provenance_bytes
    bind_core_record(provenance, core, provenance=True)
    if provenance.get("declared_traits") != declared:
        raise ValueError(
            "Current provenance declared trait metadata differs from the supplied list"
        )
    bind_producer(provenance, trial)
    helper = Path(__file__).resolve()
    snapshots[helper] = helper.read_bytes()
    payload = {
        "issue": 1117,
        "producer": {"trial": trial.hexsha, "parents": expected_parents, "tree": trial.tree.hexsha},
        "pinned_core": {
            "path": str(core_manifest),
            **identity(snapshots[core_manifest]),
            **core,
            "direction_basis": "disease_vs_comparison",
            "scientific_reverification_performed_by_this_helper": False,
        },
        "reference": old_inputs,
        "current": new_inputs,
        "current_provenance": identity(provenance_bytes),
        "trait_list": {
            "path": str(trait_path),
            "requested_path": str(args.trait_list),
            **identity(trait_bytes),
        },
        "helper": identity(snapshots[helper]),
        "removed_traits": [trait for trait in old_traits if trait not in traits],
        "tables": results,
        "summary_unchanged_after_declared_additions": True,
        "new_exclusion_counters": {
            key: summary_new["coverage_and_cascade"][key] for key in ADDITIONS
        },
        "unexpected_changes": [],
        "historical_clean_producer_verified": False,
        "limitation": "Archived consumer-bound numerical reference only; historical sidecar "
        "producer was not independently bound. Trial retains historical F core policy, "
        "not canonical-union P outputs. This helper checks core byte identity and declared "
        "scope; scientific core verification belongs to the real producer command.",
    }
    for path, raw in snapshots.items():
        if path.read_bytes() != raw:
            raise ValueError(f"Input changed during comparison: {path}")
    if repo.head.commit != trial or repo.is_dirty(untracked_files=True):
        raise ValueError("Pinned trial checkout changed or became dirty during comparison")
    output.mkdir(parents=True, exist_ok=False)
    report = (json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n").encode()
    (output / "parity.json").write_bytes(report)
    return {"output": str(output / "parity.json"), **identity(report), "tables": results}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in (
        "reference-consumer",
        "current-consumer",
        "core-manifest",
        "trait-list",
        "output-dir",
    ):
        parser.add_argument(f"--{name}", type=Path, required=True)
    for name in (
        "reference-sha256",
        "current-sha256",
        "core-sha256",
        "trial-commit",
        "issue-commit",
        "compatibility-commit",
    ):
        parser.add_argument(f"--{name}", required=True)
    print(json.dumps(compare(parser.parse_args()), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
