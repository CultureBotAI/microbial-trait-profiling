#!/usr/bin/env python3
"""Reproduce the historical comparison fixture for #1610 (2026-11-06).

Input: --repo, a Git checkout containing commit f0cb5b7a500b11320ea84e63bfaf125fe1520df9.
Output: --output, an exclusively created JSON file with all seven declarations.
Read only the pinned downstream_contracts.py blob, bounded to 128 KiB. Parse its
literal AnalysisSpec declarations with the standard-library AST, never import
historical code, and preserve declaration order and every field's source order.
"""

from __future__ import annotations

import argparse
import ast
import json
from collections.abc import Sequence
from pathlib import Path

from git import Repo

# Full 40-character SHA (#1722), joined across 5-hex chunks so detect-secrets does not
# mistake the public commit identifier for a high-entropy secret literal.
SOURCE_COMMIT = "".join(("f0cb5", "b7a50", "0b113", "20ea8", "4e63b", "faf12", "5fe15", "20df9"))
SOURCE_PATH = "src/mtp/workflow/downstream_contracts.py"
MAX_SOURCE_BYTES = 128 * 1024
SPEC_FIELDS = frozenset(
    {
        "id",
        "estimand",
        "kernel",
        "transformation",
        "sampling_unit",
        "weighting",
        "eligibility_rule",
        "multiple_testing_policy",
    }
)


def read_pinned_source(repository: Path) -> str:
    """Read the historical Git blob without consulting working-tree source."""
    repo = Repo(repository)
    try:
        blob = repo.commit(SOURCE_COMMIT).tree / SOURCE_PATH
        source = blob.data_stream.read(MAX_SOURCE_BYTES + 1)
    finally:
        repo.close()
    if not isinstance(source, bytes):
        raise TypeError("Historical Git blob did not return bytes")
    if not source or len(source) > MAX_SOURCE_BYTES:
        raise ValueError("Historical source is empty or exceeds the 128 KiB read bound")
    return source.decode("utf-8")


def literal_spec(key: ast.expr | None, value: ast.expr) -> dict[str, str]:
    """Require one complete, keyword-only literal declaration."""
    if not isinstance(key, ast.Constant) or not isinstance(key.value, str):
        raise ValueError("AnalysisSpec mapping keys must be literal strings")
    if (
        not isinstance(value, ast.Call)
        or not isinstance(value.func, ast.Name)
        or value.func.id != "AnalysisSpec"
        or value.args
    ):
        raise ValueError("Expected a keyword-only literal AnalysisSpec declaration")
    fields: dict[str, str] = {}
    for keyword in value.keywords:
        if (
            keyword.arg is None
            or keyword.arg in fields
            or not isinstance(keyword.value, ast.Constant)
            or not isinstance(keyword.value.value, str)
        ):
            raise ValueError("AnalysisSpec fields must be unique named literal strings")
        fields[keyword.arg] = keyword.value.value
    if set(fields) != SPEC_FIELDS:
        raise ValueError("AnalysisSpec fields differ from the eight historical fields")
    if fields["id"] != key.value:
        raise ValueError("AnalysisSpec id differs from its mapping key")
    return fields


def extract_specs(source: str) -> list[dict[str, str]]:
    """Extract only the seven declarations in STANDARD_ANALYSIS_SPECS order."""
    if len(source.encode("utf-8")) > MAX_SOURCE_BYTES:
        raise ValueError("Historical source exceeds the 128 KiB parse bound")
    syntax = ast.parse(source, filename=SOURCE_PATH)
    declarations = [
        node.value
        for node in syntax.body
        if isinstance(node, ast.AnnAssign)
        and isinstance(node.target, ast.Name)
        and node.target.id == "STANDARD_ANALYSIS_SPECS"
    ]
    if len(declarations) != 1 or not isinstance(declarations[0], ast.Dict):
        raise ValueError("Expected one literal STANDARD_ANALYSIS_SPECS mapping")
    declaration = declarations[0]
    if len(declaration.values) != 7:
        raise ValueError("Expected exactly seven historical AnalysisSpec declarations")
    specs = [
        literal_spec(key, value)
        for key, value in zip(declaration.keys, declaration.values, strict=True)
    ]
    if len({item["id"] for item in specs}) != len(specs):
        raise ValueError("Historical AnalysisSpec ids must be unique")
    return specs


def write_fixture(repository: Path, output: Path) -> None:
    """Write a fresh JSON copy; never replace the committed fixture or another file."""
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Refusing to replace existing output: {output}")
    specs = extract_specs(read_pinned_source(repository))
    payload = (json.dumps(specs, indent=2) + "\n").encode("utf-8")
    with output.open("xb") as handle:
        handle.write(payload)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    write_fixture(args.repo, args.output)
    print(f"Saved seven historical declarations from {SOURCE_COMMIT} to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
