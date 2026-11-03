"""Inventory PR #1219's deferred deletions and first-party references (#1344).

Date: 2026-11-02. Inputs: repository and a saved Git deleted-path list.
Outputs: fresh TSV inventories plus a JSON search-scope/identity record.
Recreate the 522-path input list (#1358) with:

    git diff --diff-filter=D --name-only ab9156f 1eda796

Save that command's stdout to a new file and pass it as --deleted-list.
Reads ignored files without consulting Git ignore rules. Deletes nothing;
does not import repository modules, execute commands, or parse XML/binary assets.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

ROOTS = ("docs", "methods", "tests", "scripts", "src", ".github", "schema", "workflows")
SKIP_DIRS = {
    ".git",
    ".venv",
    "venv",
    "env",
    "node_modules",
    "vendor",
    "__pycache__",
    ".mypy_cache",
    ".pytest_cache",
    ".ruff_cache",
    "site-packages",
}
TEXT_SUFFIXES = {
    ".md",
    ".rst",
    ".txt",
    ".py",
    ".sh",
    ".bash",
    ".json",
    ".jsonl",
    ".yaml",
    ".yml",
    ".toml",
    ".ini",
    ".cfg",
    ".tsv",
    ".csv",
    ".html",
    ".js",
    ".css",
}
DOCUMENT_SUFFIXES = {".md", ".rst", ".txt"}
REPORT_PREFIX = "issue1344_deferred_generated_outputs"
HELPER = "scripts/adhoc/issue1344_inventory_deferred_deletions.py"
MAX_BYTES = 10_000_000


def deleted_paths(filename: Path) -> tuple[list[str], str]:
    raw = filename.read_bytes()
    values = raw.decode("utf-8").splitlines()
    if not values or len(values) != len(set(values)):
        raise ValueError("Deleted-path list must be nonempty and unique")
    for value in values:
        relative = Path(value)
        if relative.is_absolute() or ".." in relative.parts or value != relative.as_posix():
            raise ValueError(f"Invalid repository-relative deletion: {value}")
        if not value.startswith(("data/analysis/", "trait-visualization/")):
            raise ValueError(f"Unexpected non-output deletion requires review: {value}")
    return sorted(values), hashlib.sha256(raw).hexdigest()


def source_files(repository: Path) -> tuple[list[Path], list[dict[str, str]]]:
    """Enumerate first-party text directly, including hidden and ignored files."""
    found: set[Path] = set()
    skipped: list[dict[str, str]] = []
    roots = [(name, False) for name in ROOTS]
    roots.extend([("trait-visualization", True), ("data/analysis", True)])
    for name, documents_only in roots:
        root = repository / name
        if not root.exists():
            skipped.append({"path": name, "reason": "root_absent"})
            continue
        if root.is_symlink():
            skipped.append({"path": name, "reason": "symlink_root"})
            continue
        for directory, children, filenames in os.walk(root, followlinks=False):
            parent = Path(directory)
            retained = []
            for child in sorted(children):
                candidate = parent / child
                if child in SKIP_DIRS or candidate.is_symlink():
                    skipped.append(
                        {
                            "path": candidate.relative_to(repository).as_posix(),
                            "reason": "excluded_directory_or_symlink",
                        }
                    )
                else:
                    retained.append(child)
            children[:] = retained
            for filename in sorted(filenames):
                candidate = parent / filename
                relative = candidate.relative_to(repository).as_posix()
                if relative == HELPER or filename.startswith(REPORT_PREFIX):
                    continue
                suffixes = DOCUMENT_SUFFIXES if documents_only else TEXT_SUFFIXES
                if candidate.suffix.lower() not in suffixes and filename != "Justfile":
                    continue
                found.add(candidate)
    for candidate in repository.iterdir():
        if candidate.is_file() and (
            candidate.suffix.lower() in DOCUMENT_SUFFIXES or candidate.name == "Justfile"
        ):
            found.add(candidate)
    return sorted(found), skipped


def token_map(paths: list[str], *, directories: bool = False) -> dict[str, set[str]]:
    result: dict[str, set[str]] = defaultdict(set)
    for value in paths:
        if directories:
            for parent in Path(value).parents:
                if parent.as_posix() != ".":
                    result[parent.as_posix()].add(parent.as_posix())
                    result[parent.name].add(parent.as_posix())
        else:
            result[value].add(value)
            result[Path(value).name].add(value)
    return dict(result)


def matcher(tokens: dict[str, set[str]]) -> re.Pattern[str]:
    return re.compile(
        "|".join(re.escape(token) for token in sorted(tokens, key=lambda s: (-len(s), s)))
    )


def write_tsv(filename: Path, columns: list[str], rows: list[dict[str, Any]]) -> None:
    with filename.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def grouped_references(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Store each source/target/token once while retaining every matched line."""
    groups: dict[tuple[tuple[str, Any], ...], list[int]] = defaultdict(list)
    for row in rows:
        key = tuple((name, value) for name, value in row.items() if name != "line")
        groups[key].append(row["line"])
    return [
        {**dict(key), "line_numbers": ",".join(map(str, sorted(lines)))}
        for key, lines in groups.items()
    ]


def inventory(repository: Path, paths: list[str]) -> dict[str, Any]:
    candidates, skipped = source_files(repository)
    path_tokens, directory_tokens = token_map(paths), token_map(paths, directories=True)
    file_pattern, directory_pattern = matcher(path_tokens), matcher(directory_tokens)
    references: list[dict[str, Any]] = []
    directories: list[dict[str, Any]] = []
    scanned: list[dict[str, Any]] = []
    deleted = set(paths)
    counts: dict[str, Counter[str]] = defaultdict(Counter)
    for candidate in candidates:
        relative = candidate.relative_to(repository).as_posix()
        if candidate.is_symlink():
            skipped.append({"path": relative, "reason": "symlink_file"})
            continue
        before = candidate.stat()
        if before.st_size > MAX_BYTES:
            skipped.append({"path": relative, "reason": "over_10_MB"})
            continue
        raw = candidate.read_bytes()
        after = candidate.stat()
        if (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns) != (
            after.st_dev,
            after.st_ino,
            after.st_size,
            after.st_mtime_ns,
        ):
            raise ValueError(f"Source changed during scan: {relative}")
        try:
            content = raw.decode("utf-8")
        except UnicodeDecodeError:
            skipped.append({"path": relative, "reason": "not_UTF8"})
            continue
        if "\x00" in content:
            skipped.append({"path": relative, "reason": "binary_NUL"})
            continue
        scanned.append(
            {"path": relative, "sha256": hashlib.sha256(raw).hexdigest(), "size_bytes": len(raw)}
        )
        for number, line in enumerate(content.splitlines(), 1):
            for token in sorted(set(file_pattern.findall(line))):
                for target in sorted(path_tokens[token]):
                    kind = "literal_path" if token == target else "basename"
                    references.append(
                        {
                            "deleted_path": target,
                            "source_path": relative,
                            "line": number,
                            "match_kind": kind,
                            "matched_token": token,
                            "candidate_paths_for_token": len(path_tokens[token]),
                            "source_is_deletion_candidate": relative in deleted,
                        }
                    )
                    counts[target][kind] += 1
                    if relative not in deleted:
                        counts[target]["surviving_source_hits"] += 1
            for token in sorted(set(directory_pattern.findall(line))):
                for target in sorted(directory_tokens[token]):
                    directories.append(
                        {
                            "directory": target,
                            "source_path": relative,
                            "line": number,
                            "match_kind": "directory_path"
                            if token == target
                            else "directory_basename",
                            "matched_token": token,
                            "candidate_directories_for_token": len(directory_tokens[token]),
                            "source_is_deletion_candidate": relative in deleted,
                        }
                    )
    files = []
    for target in paths:
        candidate = repository / target
        state = (
            "symlink" if candidate.is_symlink() else "present" if candidate.is_file() else "absent"
        )
        files.append(
            {
                "deleted_path": target,
                "directory": Path(target).parent.as_posix(),
                "suffix": Path(target).suffix,
                "checkout_state": state,
                "literal_path_hits": counts[target]["literal_path"],
                "basename_hits": counts[target]["basename"],
                "surviving_source_hits": counts[target]["surviving_source_hits"],
                "decision": "deferred_no_deletion",
            }
        )
    digest = hashlib.sha256(
        json.dumps(scanned, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return {
        "files": files,
        "references": references,
        "directories": directories,
        "scanned": scanned,
        "skipped": skipped,
        "scanned_inventory_sha256": digest,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository", type=Path, required=True)
    parser.add_argument("--deleted-list", type=Path, required=True)
    parser.add_argument("--source-base", required=True)
    parser.add_argument("--source-head", required=True)
    parser.add_argument("--review-base", required=True)
    parser.add_argument("--output", type=Path, required=True, help="Fresh report directory")
    args = parser.parse_args()
    repository = args.repository.resolve(strict=True)
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError("Report directory must be fresh")
    paths, input_sha = deleted_paths(args.deleted_list)
    result = inventory(repository, paths)
    args.output.mkdir(parents=True, exist_ok=False)
    write_tsv(args.output / "files.tsv", list(result["files"][0]), result["files"])
    for name in ("references", "directories"):
        if not result[name]:
            raise ValueError(f"Unexpected empty {name} inventory; inspect scan scope")
        grouped = grouped_references(result[name])
        write_tsv(args.output / f"{name}.tsv", list(grouped[0]), grouped)
    metadata = {
        "issue": 1344,
        "source_pr": 1219,
        "source_base": args.source_base,
        "source_head": args.source_head,
        "review_base": args.review_base,
        "deleted_list_sha256": input_sha,
        "deleted_count": len(paths),
        "helper_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "ignored_files_included": True,
        "first_party_roots": list(ROOTS),
        "documents_only_roots": ["trait-visualization", "data/analysis"],
        "root_files": "Justfile and Markdown/reStructuredText/text files",
        "excluded_directory_names": sorted(SKIP_DIRS),
        "max_file_bytes": MAX_BYTES,
        "scanned_file_count": len(result["scanned"]),
        "scanned_inventory_sha256": result["scanned_inventory_sha256"],
        "reference_count": len(result["references"]),
        "directory_reference_count": len(result["directories"]),
        "reference_storage": "Grouped by source/target/token; every matched line is retained",
        "scanned_files": result["scanned"],
        "skipped": result["skipped"],
        "limitations": [
            "Lexical matches do not establish runtime consumption or regeneration equivalence.",
            "Basename and directory matches can be ambiguous; inspect the cited source lines.",
            "No-hit paths are not proven unused: constructed names, binary decks "
            "and external consumers remain.",
            "The review-base names the starting checkout; scanned file hashes "
            "bind its observed working bytes.",
        ],
    }
    with (args.output / "scan.json").open("x", encoding="utf-8") as stream:
        json.dump(metadata, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(
        json.dumps(
            {
                key: metadata[key]
                for key in (
                    "deleted_count",
                    "scanned_file_count",
                    "reference_count",
                    "directory_reference_count",
                )
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
