"""Audit list/catalog compatibility for the off quadrant (#1115, 2026-10-28).

Inputs: an explicitly pinned completed F parent, the current declared trait list,
and optional reference files whose byte identities are recorded without verification.
Outputs: a fresh JSON membership report under a caller-supplied durable location.
This is a read-only byte/identity audit, not full source verification or a rerun.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import stat
from collections import deque
from pathlib import Path
from typing import Any

OXYGEN_TRAITS = ("METPO:1000608", "METPO:1000607")
LEDGER_KEYS = {
    "composition_pair_delta_z.tsv": ("study_id", "experiment_id"),
    "composition_study_delta_z.tsv": ("study_id",),
}
_MAX_SYMLINK_EXPANSIONS = 40


def identity(raw: bytes) -> dict[str, str | int]:
    return {"sha256": hashlib.sha256(raw).hexdigest(), "size_bytes": len(raw)}


def _spellings(path: Path) -> list[Path]:
    original = path if path.is_absolute() else Path.cwd() / path
    spellings = [original]
    visited = {original}
    pending = deque(original.parts[1:])
    resolved = Path(original.anchor)
    expansions = 0
    while pending:
        component = pending.popleft()
        if component == "..":
            resolved = resolved.parent
            continue
        candidate = resolved / component
        try:
            mode = candidate.lstat().st_mode
        except FileNotFoundError:
            resolved = candidate
            continue
        if not stat.S_ISLNK(mode):
            resolved = candidate
            continue
        expansions += 1
        if expansions > _MAX_SYMLINK_EXPANSIONS:
            raise ValueError("Symlink expansion limit reached while protecting input bundles")
        spellings.append(candidate.joinpath(*pending))
        target = candidate.readlink()
        expanded = target if target.is_absolute() else candidate.parent / target
        spelling = expanded.joinpath(*pending)
        if spelling in visited:
            raise ValueError("Symlink cycle while protecting input bundles")
        visited.add(spelling)
        spellings.append(spelling)
        if target.is_absolute():
            resolved = Path(target.anchor)
            target_parts = target.parts[1:]
        else:
            resolved = candidate.parent
            target_parts = target.parts
        pending.extendleft(reversed(target_parts))
    spellings.append(resolved)
    return list(dict.fromkeys(spellings))


def _expanded_spellings(path: Path) -> list[Path]:
    candidates = [*_spellings(path), Path(os.path.abspath(path))]
    result: list[Path] = []
    for candidate in candidates:
        for item in (*_spellings(candidate), Path(os.path.abspath(candidate)), candidate.resolve()):
            if item not in result:
                result.append(item)
    return result


def validate_output_path(
    output: Path,
    parent: Path,
    trait_list: Path,
    reference_files: list[Path],
) -> None:
    output_spellings = _expanded_spellings(output)
    protected_roots = [
        spelling if spelling.is_dir() else spelling.parent
        for bundle_input in (parent, *reference_files)
        for spelling in _expanded_spellings(bundle_input)
    ]
    protected = [
        *(spelling for root in protected_roots for spelling in _expanded_spellings(root)),
        *_expanded_spellings(trait_list),
    ]
    if any(
        candidate == source or candidate.is_relative_to(source) or source.is_relative_to(candidate)
        for candidate in output_spellings
        for source in protected
    ):
        raise ValueError(
            "Report output must be outside completed input bundles and disjoint from audit inputs"
        )


def audit(
    parent: Path,
    parent_sha256: str,
    trait_list: Path,
    reference_files: list[Path],
    output: Path | None = None,
) -> dict[str, Any]:
    if output is not None:
        validate_output_path(output, parent, trait_list, reference_files)
    parent_bytes = parent.read_bytes()
    if hashlib.sha256(parent_bytes).hexdigest() != parent_sha256:
        raise ValueError("Parent manifest differs from the caller's explicit pin")
    record: dict[str, Any] = json.loads(parent_bytes)
    if record["strain_rollup"] != "off":
        raise ValueError("This audit requires the off sensitivity parent")
    retained: list[str] = record["retained_traits"]
    if not retained or len(retained) != len(set(retained)):
        raise ValueError("Parent retained catalog must be nonempty and unique")
    snapshots = {parent: parent_bytes}
    outputs = {}
    for name, keys in LEDGER_KEYS.items():
        path = parent.parent / name
        raw = path.read_bytes()
        if identity(raw) != record["outputs"][name]:
            raise ValueError(f"Parent output inventory mismatch: {name}")
        header = next(csv.reader(io.StringIO(raw.decode("utf-8")), delimiter="\t"))
        if header != [*keys, *retained]:
            raise ValueError(f"Parent output header differs from retained catalog: {name}")
        snapshots[path] = raw
        outputs[name] = identity(raw)
    evidence_path = parent.parent / "per_trait_evidence.tsv"
    evidence_bytes = evidence_path.read_bytes()
    if identity(evidence_bytes) != record["outputs"][evidence_path.name]:
        raise ValueError("Parent evidence inventory mismatch")
    evidence = list(csv.DictReader(io.StringIO(evidence_bytes.decode("utf-8")), delimiter="\t"))
    if [row["trait"] for row in evidence] != retained:
        raise ValueError("Evidence vocabulary/order differs from retained catalog")
    snapshots[evidence_path] = evidence_bytes
    outputs[evidence_path.name] = identity(evidence_bytes)
    list_bytes = trait_list.read_bytes()
    declared = [
        value
        for line in list_bytes.decode("utf-8").splitlines()
        if (value := line.strip()) and not value.startswith("#")
    ]
    if not declared or len(declared) != len(set(declared)):
        raise ValueError("Declared list must be nonempty and unique")
    snapshots[trait_list] = list_bytes
    missing = [trait for trait in declared if trait not in retained]
    oxygen_missing = [trait for trait in OXYGEN_TRAITS if trait not in retained]
    references = []
    for path in reference_files:
        raw = path.read_bytes()
        snapshots[path] = raw
        references.append({"path": str(path.resolve()), **identity(raw)})
    for path, raw in snapshots.items():
        if path.read_bytes() != raw:
            raise ValueError(f"Input changed during the read-only audit: {path.name}")
    return {
        "issue": 1115,
        "audit_kind": "input_compatibility_only_not_full_verification",
        "helper_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "parent_manifest_path": str(parent.resolve()),
        "parent_manifest": identity(parent_bytes),
        "parent_run_id": record["run_id"],
        "parent_code_commit": record["code_commit"],
        "strain_rollup": record["strain_rollup"],
        "contrast_population": record["contrast_population"],
        "direction_basis": record["direction_basis"],
        "parent_outputs": outputs,
        "parent_output_paths": {name: str((parent.parent / name).resolve()) for name in outputs},
        "reference_file_identities_only": references,
        "retained_trait_count": len(retained),
        "trait_list": {
            "filename": trait_list.name,
            "path": str(trait_list.resolve()),
            **identity(list_bytes),
            "traits": declared,
        },
        "missing_declared_traits": missing,
        "missing_oxygen_traits": oxygen_missing,
        "legacy_blocking_trait": {
            "trait": "METPO:1000605",
            "in_current_list": "METPO:1000605" in declared,
            "in_off_catalog": "METPO:1000605" in retained,
        },
        "list_membership_blocks_off_quadrant": bool(missing or oxygen_missing),
        "full_upstream_verification_performed": False,
        "analysis_rerun_performed": False,
        "input_bytes_unchanged": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parent-manifest", type=Path, required=True)
    parser.add_argument("--parent-sha256", required=True)
    parser.add_argument("--trait-list", type=Path, required=True)
    parser.add_argument("--reference-file", type=Path, action="append", default=[])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    validate_output_path(args.output, args.parent_manifest, args.trait_list, args.reference_file)
    result = audit(
        args.parent_manifest,
        args.parent_sha256,
        args.trait_list,
        args.reference_file,
        output=args.output,
    )
    validate_output_path(args.output, args.parent_manifest, args.trait_list, args.reference_file)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2)
        handle.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
