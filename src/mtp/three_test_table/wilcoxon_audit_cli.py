"""Publish a pinned-parent Wilcoxon sensitivity audit without rerunning sources (#1073)."""

from __future__ import annotations

import argparse
import json
import re
import shutil
import uuid
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from mtp.multiple_testing import canonical_sha256
from mtp.paths import file_sha256
from mtp.three_test_table import stage
from mtp.three_test_table.wilcoxon_audit import audit_parent
from mtp.workflow.output_roots import temporary_root

VERSION = "three-test-wilcoxon-audit-v1"
MANIFEST = "wilcoxon_audit_run.json"
TABLE = "wilcoxon_audit.tsv"
METHODS = "methods/three_test_wilcoxon_audit.md"
_PARENT_MANIFEST = "three_test_run.json"
_PARENT_TABLES = (
    "per_trait_evidence.tsv",
    "composition_pair_delta_z.tsv",
    "composition_study_delta_z.tsv",
)


def _manifest_aliases(supplied: Path, target: Path) -> set[Path]:
    """Keep every file link in this pinned manifest's chain as an input identity."""
    aliases = {supplied, target}
    visited: set[Path] = set()
    current = supplied
    # resolve(strict=True) already refused loops; bound traversal also handles races.
    for _ in range(40):
        physical = current.parent.resolve(strict=True) / current.name
        aliases.add(physical)
        if physical in visited:
            raise ValueError("Parent manifest alias cycle")
        visited.add(physical)
        if not physical.is_symlink():
            if physical.resolve(strict=True) != target:
                raise ValueError("Parent manifest alias target changed during preflight")
            return aliases
        link = physical.readlink()
        current = link if link.is_absolute() else physical.parent / link
        aliases.add(current)
    raise ValueError("Parent manifest alias expansion limit exceeded")


def _parent_paths(parent: Path) -> tuple[Path, Path, set[Path], set[Path]]:
    supplied = (parent / _PARENT_MANIFEST if parent.is_dir() else parent).absolute()
    manifest = supplied.resolve(strict=True)
    root = manifest.parent
    paths = {supplied, manifest, *(root / name for name in _PARENT_TABLES)}
    if supplied.parent.resolve(strict=True) == root:
        paths.update(supplied.parent / name for name in _PARENT_TABLES)
    return supplied, root, paths, _manifest_aliases(supplied, manifest)


def _recheck_aliases(targets: dict[Path, Path], expected_sha256: str) -> None:
    for path, target in targets.items():
        if path.resolve(strict=True) != target or file_sha256(path) != expected_sha256:
            raise ValueError("Parent source manifest alias changed after audit")


def _producer(code_commit: str | None) -> dict[str, Any]:
    producer = dict(stage._producer(code_commit))
    inventory = dict(producer["code_inventory"])
    inventory[METHODS] = file_sha256(Path(__file__).resolve().parents[3] / METHODS)
    producer["code_inventory"] = inventory
    producer["implementation_sha256"] = canonical_sha256(inventory)
    producer["audit_methods_sha256"] = inventory[METHODS]
    return producer


def _input_inventory(
    hashes: dict[Path, str], targets: dict[Path, Path], root: Path
) -> dict[str, dict[str, str | int]]:
    if set(hashes) != set(targets):
        raise ValueError("Audit input inventory differs from the protected parent files")
    inventory: dict[str, dict[str, str | int]] = {}
    for path, expected in hashes.items():
        target = path.resolve(strict=True)
        if target != targets[path] or file_sha256(path) != expected:
            raise ValueError("Parent source input changed after audit")
        if target.parent != root:
            raise ValueError("Audit input escapes the pinned parent root")
        name = target.name
        item: dict[str, str | int] = {"sha256": expected, "size_bytes": path.stat().st_size}
        if name in inventory and inventory[name] != item:
            raise ValueError("Conflicting audit input identities")
        inventory[name] = item
    if set(inventory) != {_PARENT_MANIFEST, *_PARENT_TABLES}:
        raise ValueError("Audit input inventory is not the complete pinned parent")
    return inventory


def run_audit(
    parent: Path,
    output: Path,
    *,
    expected_sha256: str,
    scope: str,
    code_commit: str | None = None,
) -> dict[str, Any]:
    """Read pinned historical outputs and atomically publish one fresh audit bundle."""
    if re.fullmatch(r"[0-9a-f]{64}", expected_sha256) is None:
        raise ValueError("A lowercase SHA256 parent pin is required")
    if scope not in {"off", "phenotype"}:
        raise ValueError("Audit scope must be off or phenotype")
    if temporary_root(output) is not None:
        raise ValueError("Use a durable directory for audit output")
    if output.exists() or output.is_symlink():
        raise FileExistsError("Choose a fresh audit output directory")
    output = output.resolve()
    manifest, root, paths, aliases = _parent_paths(parent)
    protected = stage._manifest_roots(list(paths | aliases), multiplicity=None, label_support=None)
    stage._disjoint(output, protected)
    targets = {path: path.resolve(strict=True) for path in paths}
    alias_targets = {path: path.resolve(strict=True) for path in aliases}
    if set(alias_targets.values()) != {targets[manifest]}:
        raise ValueError("Parent manifest alias target changed during preflight")
    producer = _producer(code_commit)
    result = audit_parent(manifest, expected_sha256, scope=scope)
    _recheck_aliases(alias_targets, expected_sha256)
    inputs = _input_inventory(result.input_hashes, targets, root)
    if inputs[_PARENT_MANIFEST]["sha256"] != expected_sha256:
        raise ValueError("Audit parent input differs from the requested SHA256 pin")
    if _producer(code_commit) != producer:
        raise ValueError("Audit producer changed during calculation")
    record: dict[str, Any] = {
        "version": VERSION,
        "run_id": f"wilcoxon-audit-{uuid.uuid4().hex}",
        "recorded_at": datetime.now(UTC).isoformat(),
        **producer,
        "parent": {"manifest": _PARENT_MANIFEST, "sha256": expected_sha256, "scope": scope},
        "inputs": inputs,
        "summary": result.summary,
    }
    # Refuse malformed audit metadata before creating any requested output ancestors.
    json.dumps(record, allow_nan=False)
    stage._disjoint(output, protected)
    ancestor = output.parent
    while not ancestor.exists() and not ancestor.is_symlink():
        ancestor = ancestor.parent
    if not ancestor.is_dir():
        raise NotADirectoryError("Audit output ancestor is not a directory")
    staging = ancestor / f".{output.name}-staging-{uuid.uuid4().hex}"
    stage._disjoint(staging, protected)
    staging.mkdir()
    published = False
    try:
        result.table.to_csv(staging / TABLE, sep="\t", index=False, lineterminator="\n")
        record["outputs"] = {TABLE: stage._inventory(staging / TABLE)}
        _recheck_aliases(alias_targets, expected_sha256)
        if _input_inventory(result.input_hashes, targets, root) != inputs:
            raise ValueError("Audit input inventory changed before publication")
        if _producer(code_commit) != producer:
            raise ValueError("Audit producer changed before publication")
        (staging / MANIFEST).write_text(
            json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8"
        )
        stage._disjoint(output, protected)
        output.parent.mkdir(parents=True, exist_ok=True)
        stage._publish_directory(staging, output)
        published = True
    finally:
        if not published:
            shutil.rmtree(staging)
    return record


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        epilog="Requires nonempty pinned parent pair and study ledgers; empty-population audits "
        "are outside this audit's input scope.",
    )
    parser.add_argument("--parent", type=Path, required=True)
    parser.add_argument("--parent-sha256", required=True)
    parser.add_argument("--scope", choices=("off", "phenotype"), required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--code-commit")
    args = parser.parse_args(argv)
    record = run_audit(
        args.parent,
        args.output,
        expected_sha256=args.parent_sha256,
        scope=args.scope,
        code_commit=args.code_commit,
    )
    print(f"{record['run_id']}: {record['summary']['family_size']} traits; {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
