"""Audit two uv lock files without modifying them (#1286, 2026-10-30).

Inputs: --before and --after TOML files. Output: a fresh --output JSON report.
Only the root revision and package sdist/wheel upload-time fields are reported
separately as format metadata. All other parsed fields and array order remain
in the normalized unified diff. This is not a resolver or security audit.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import tomllib
from collections.abc import Sequence
from copy import deepcopy
from difflib import unified_diff
from pathlib import Path
from typing import Any


def identity(raw: bytes) -> dict[str, str | int]:
    return {"sha256": hashlib.sha256(raw).hexdigest(), "size_bytes": len(raw)}


def normalize(record: dict[str, Any]) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Remove only explicitly located format fields; retain dictionary/list order."""
    normalized = deepcopy(record)
    metadata: list[dict[str, Any]] = []
    if "revision" in normalized:
        metadata.append({"path": ["revision"], "value": normalized.pop("revision")})
    for index, package in enumerate(normalized.get("package", [])):
        sdist = package.get("sdist")
        if isinstance(sdist, dict) and "upload-time" in sdist:
            metadata.append(
                {
                    "path": ["package", index, "sdist", "upload-time"],
                    "value": sdist.pop("upload-time"),
                }
            )
        for wheel_index, wheel in enumerate(package.get("wheels", [])):
            if "upload-time" in wheel:
                metadata.append(
                    {
                        "path": ["package", index, "wheels", wheel_index, "upload-time"],
                        "value": wheel.pop("upload-time"),
                    }
                )
    return normalized, metadata


def package_versions(record: dict[str, Any]) -> dict[str, list[str | None]]:
    """Keep all entries in source order, including multiple versions of one name."""
    result: dict[str, list[str | None]] = {}
    packages = record.get("package", [])
    if not isinstance(packages, list):
        raise ValueError("Lock package inventory must be an array")
    for package in packages:
        if not isinstance(package, dict) or not isinstance(package.get("name"), str):
            raise ValueError("Each lock package must declare a name")
        version = package.get("version")
        if version is not None and not isinstance(version, str):
            raise ValueError("A package version must be a string")
        result.setdefault(package["name"], []).append(version)
    return result


def compare_locks(before: Path, after: Path) -> dict[str, Any]:
    """Hash and parse the same input buffers; report every non-format difference."""
    snapshots: dict[str, dict[str, Any]] = {}
    normalized_text: dict[str, str] = {}
    versions: dict[str, dict[str, list[str | None]]] = {}
    metadata: dict[str, list[dict[str, Any]]] = {}
    for role, path in (("before", before), ("after", after)):
        raw = path.read_bytes()
        record = tomllib.loads(raw.decode("utf-8"))
        versions[role] = package_versions(record)
        normalized, metadata[role] = normalize(record)
        # No sorting: array order and parsed mapping order must remain visible.
        text = json.dumps(normalized, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
        normalized_text[role] = text
        snapshots[role] = {
            "path": str(path.resolve()),
            **identity(raw),
            "normalized_sha256": identity(text.encode("utf-8"))["sha256"],
        }
    names = dict.fromkeys((*versions["before"], *versions["after"]))
    changes = {
        name: {"before": versions["before"].get(name, []), "after": versions["after"].get(name, [])}
        for name in names
        if versions["before"].get(name, []) != versions["after"].get(name, [])
    }
    diff = "".join(
        unified_diff(
            normalized_text["before"].splitlines(keepends=True),
            normalized_text["after"].splitlines(keepends=True),
            fromfile="before.normalized.json",
            tofile="after.normalized.json",
        )
    )
    return {
        "version": "uv-lock-audit-v1",
        "verification_scope": "Parsed lock comparison only; no resolution or security verification",
        **snapshots,
        "substantive_equal": not diff,
        "package_version_changes": changes,
        "format_metadata": metadata,
        "normalized_diff": diff,
        "helper": identity(Path(__file__).read_bytes()),
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--before", type=Path, required=True)
    parser.add_argument("--after", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    result = compare_locks(args.before, args.after)
    content = json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    # Exclusive creation also refuses existing or dangling symlinks.
    with args.output.open("x", encoding="utf-8") as handle:
        handle.write(content)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
