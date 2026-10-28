"""Read-only full historical core verification for #1129 (2026-10-27).

Inputs: explicitly pinned core manifest, optionally an explicitly pinned evidence
archive. Outputs: a fresh durable JSON record with full artifact inventories
checked before/after, verifier result, producer identity and timing. No core,
profile, arm or downstream result is regenerated or modified.
"""

from __future__ import annotations

import argparse
import json
import time
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from mtp.datamodel.workflow_pydantic import WorkflowRun
from mtp.paths import REPO_ROOT, file_sha256
from mtp.workflow.archive import ArchiveResolver
from mtp.workflow.locations import LocationResolver, locate
from mtp.workflow.operations import OPERATIONS
from mtp.workflow.output_roots import require_durable_output
from mtp.workflow.runtime import _code_identity, artifact_path, inventory, repository_state
from mtp.workflow.trait_policy import resolve_trait_policy
from mtp.workflow.verification import verify_core_run


def _inventories(run: WorkflowRun, resolver: LocationResolver | None) -> dict[str, Any]:
    kinds = {item.id: str(item.storage_kind) for item in run.definition.products}
    result = {}
    for artifact in run.artifacts or []:
        path = artifact_path(run, artifact.product_id, resolver)
        digest, size, members = inventory(path, kinds[artifact.product_id])
        if (digest, size, members) != (
            artifact.sha256,
            artifact.size_bytes,
            artifact.inventory,
        ):
            raise ValueError(f"Historical artifact inventory changed: {artifact.product_id}")
        result[artifact.product_id] = {
            "path": str(path),
            "sha256": digest,
            "size_bytes": size,
            "members": [member.model_dump(mode="json") for member in members],
        }
    return result


def _producer() -> dict[str, Any]:
    commit, dirty = repository_state(REPO_ROOT)
    code_hash, entries = _code_identity(REPO_ROOT, OPERATIONS)
    return {
        "commit": commit,
        "dirty": dirty,
        "code_sha256": code_hash,
        "code_inventory": [entry.model_dump(mode="json") for entry in entries],
        "helper_sha256": file_sha256(Path(__file__)),
        "lock_sha256": file_sha256(REPO_ROOT / "uv.lock"),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--core", type=Path, required=True)
    parser.add_argument("--core-sha256", required=True)
    parser.add_argument("--archive", type=Path)
    parser.add_argument("--archive-sha256")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if (args.archive is None) != (args.archive_sha256 is None):
        parser.error("--archive and --archive-sha256 must be supplied together")
    require_durable_output(parser, args.output, allowed=False)
    output = args.output.resolve()
    if output.exists() or args.output.is_symlink() or output.is_relative_to(REPO_ROOT):
        parser.error("Choose a fresh output file outside the source repository")
    resolver: LocationResolver | None = None
    if args.archive is not None:
        if file_sha256(args.archive) != args.archive_sha256:
            raise ValueError("Archive manifest pin mismatch")
        resolver = ArchiveResolver(args.archive).resolve
    manifest = locate(args.core.absolute(), resolver).resolve(strict=True)
    if file_sha256(manifest) != args.core_sha256:
        raise ValueError("Core manifest pin mismatch")
    run = WorkflowRun.model_validate_json(manifest.read_bytes())
    run.manifest_path = str(manifest)
    protected = [manifest.parent]
    if args.archive is not None:
        protected.append(args.archive.resolve().parent)
    for product in run.definition.products:
        path = artifact_path(run, product.id, resolver).resolve()
        protected.append(path if str(product.storage_kind) == "collection" else path.parent)
    if any(output.is_relative_to(root) or root.is_relative_to(output.parent) for root in protected):
        raise ValueError("Review output overlaps a historical source tree")
    producer = _producer()
    started = datetime.now(UTC)
    timer = time.monotonic()
    before = _inventories(run, resolver)
    print(f"{started.isoformat()} source inventories verified; starting full core verification")
    error: str | None = None
    try:
        verified = verify_core_run(args.core.absolute(), resolver=resolver)
        if verified.model_dump(mode="json") != run.model_dump(mode="json"):
            raise ValueError("Verification changed the historical core record")
    except (ValueError, OSError) as failure:
        error = f"{type(failure).__name__}: {failure}"
    after = _inventories(run, resolver)
    if before != after or file_sha256(manifest) != args.core_sha256:
        raise ValueError("Historical inputs changed during verification")
    if args.archive is not None and file_sha256(args.archive) != args.archive_sha256:
        raise ValueError("Archive manifest changed during verification")
    if producer != _producer():
        raise ValueError("Verifier implementation changed during verification")
    policy = resolve_trait_policy(run)
    result = {
        "issue": 1129,
        "analysis": "read_only_full_core_verification_no_replay",
        "started_at": started.isoformat(),
        "finished_at": datetime.now(UTC).isoformat(),
        "elapsed_seconds": time.monotonic() - timer,
        "status": "passed" if error is None else "failed",
        "error": error,
        "core_manifest_sha256": args.core_sha256,
        "core_run_id": run.run_id,
        "core_producer_commit": run.git_commit,
        "archive_manifest_sha256": args.archive_sha256,
        "resolved_trait_key_policy": policy.policy,
        "include_policy_metadata": policy.include_policy_metadata,
        "producer": producer,
        "before_artifact_inventories": before,
        "after_artifact_inventories": after,
        "historical_bytes_unchanged": True,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8") as handle:
        json.dump(result, handle, sort_keys=True, indent=2)
        handle.write("\n")
    print(f"{result['finished_at']} {result['status']}: {error or run.run_id}")
    return int(error is not None)


if __name__ == "__main__":
    raise SystemExit(main())
