"""Preserve and verify retired composition PR evidence for #207 and #528.

Purpose: record final archived source hashes without relabeling historical results.
Inputs: git-archive tar files for PRs #197/#526 and their expected commit IDs.
Output: deterministic JSON manifest with corrected LF source hashes and inventory.
Date: 2026-10-16. Run with `just audit-retired-composition-prs`.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import tarfile
from pathlib import Path
from typing import Any

AUDIT_DIRECTORY = "trait-visualization/denominator-harmonization/"


def read_file(archive: tarfile.TarFile, name: str) -> bytes:
    """Read a regular archived file without extracting paths to the filesystem."""
    member = archive.getmember(name)
    if not member.isfile():
        raise ValueError(f"Expected a regular archived file: {name}")
    stream = archive.extractfile(member)
    if stream is None:
        raise ValueError(f"Cannot read archived file: {name}")
    with stream:
        return stream.read()


def source_hashes(
    archive: tarfile.TarFile, recorded: dict[str, str]
) -> tuple[dict[str, str], dict[str, dict[str, str]]]:
    """Measure LF-normalized source bytes and retain every historical discrepancy."""
    measured = {}
    corrections = {}
    for name, previous in sorted(recorded.items()):
        digest = hashlib.sha256(read_file(archive, name).replace(b"\r\n", b"\n")).hexdigest()
        measured[name] = digest
        if digest != previous:
            corrections[name] = {"historical_sha256": previous, "archived_source_sha256": digest}
    return measured, corrections


def inspect_archive(path: Path, expected_commit: str, pr: int) -> dict[str, Any]:
    """Require the pinned git-archive identity and describe the preserved evidence."""
    with path.open("rb") as stream:
        archive_sha256 = hashlib.file_digest(stream, "sha256").hexdigest()
    with tarfile.open(path, "r:") as archive:
        commit = archive.pax_headers.get("comment")
        if commit != expected_commit:
            raise ValueError(f"PR #{pr} archive commit {commit!r} differs from {expected_commit}")
        result: dict[str, Any] = {
            "git_commit": commit,
            "archive_tag": f"archive/pr-{pr}-20261016",
            "git_archive_sha256": archive_sha256,
        }
        if pr == 197:
            provenance = json.loads(read_file(archive, AUDIT_DIRECTORY + "provenance.json"))
            measured, corrections = source_hashes(archive, provenance["implementation_sha256"])
            replay = json.loads(
                read_file(archive, AUDIT_DIRECTORY + "raw_input_replay_provenance.json")
            )
            replay_measured, replay_corrections = source_hashes(archive, replay["source_sha256"])
            members = sorted(
                (
                    member
                    for member in archive.getmembers()
                    if member.isfile() and member.name.startswith(AUDIT_DIRECTORY)
                ),
                key=lambda member: member.name,
            )
            inventory = hashlib.sha256()
            for member in members:
                inventory.update(member.name.encode("utf-8") + b"\0")
                inventory.update(hashlib.sha256(read_file(archive, member.name)).digest())
            result.update(
                {
                    "historical_results_directory": AUDIT_DIRECTORY,
                    "historical_results_file_count": len(members),
                    "historical_results_bytes": sum(member.size for member in members),
                    "historical_results_inventory_sha256": inventory.hexdigest(),
                    "archived_implementation_sha256_lf": measured,
                    "historical_implementation_hash_corrections": corrections,
                    "archived_replay_source_sha256_lf": replay_measured,
                    "historical_replay_hash_corrections": replay_corrections,
                }
            )
        elif pr == 526:
            name = "tests/test_composition_denominator_callers.py"
            result["archived_caller_tests_sha256_lf"] = hashlib.sha256(
                read_file(archive, name).replace(b"\r\n", b"\n")
            ).hexdigest()
        return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pr197", type=Path, required=True)
    parser.add_argument("--pr197-commit", required=True)
    parser.add_argument("--pr526", type=Path, required=True)
    parser.add_argument("--pr526-commit", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    manifest = {
        "audit_date": "2026-10-16",
        "scope": (
            "Final archived source identity and historical evidence preservation; "
            "not a rerun, validation of historical numerical claims, or adoption of "
            "those results into current main."
        ),
        "source_hash_normalization": "CRLF to LF before SHA256",
        "inventory_hash_encoding": "sorted path UTF-8 + NUL + raw file SHA256 digest",
        "archives": {
            "197": inspect_archive(args.pr197, args.pr197_commit, 197),
            "526": inspect_archive(args.pr526, args.pr526_commit, 526),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Wrote verified archive manifest: {args.output}")


if __name__ == "__main__":
    main()
