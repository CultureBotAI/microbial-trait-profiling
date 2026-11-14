"""Harden authoritative LinkML pattern endings for #1316 (2026-11-08).

Inputs: the five committed schema files in a selected repository root.
Outputs: those same YAML files, with pattern bodies preserved and actual-end
assertions appended. Run `just gen-pydantic` separately; never patch its products.
"""

from __future__ import annotations

import argparse
import copy
import re
from pathlib import Path
from typing import Any

import yaml

SCHEMAS = ("mtp", "bsdbng", "workflow", "core_products", "multiple_testing")
END_ASSERTION = r"(?![\s\S])$"


def hardened(document: Any) -> tuple[Any, int]:
    """Return an otherwise identical schema with every old pattern hardened."""
    if isinstance(document, dict):
        result = copy.deepcopy(document)
        count = 0
        for key, value in document.items():
            if key == "pattern" and isinstance(value, str):
                if not value.endswith("$") or value.endswith(END_ASSERTION):
                    raise ValueError(f"Expected an unhardened dollar-ending pattern: {value!r}")
                result[key] = value[:-1] + END_ASSERTION
                count += 1
            else:
                result[key], changed = hardened(value)
                count += changed
        return result, count
    if isinstance(document, list):
        pairs = [hardened(value) for value in document]
        return [value for value, _ in pairs], sum(count for _, count in pairs)
    return document, 0


def rewrite(text: str) -> tuple[str, int]:
    """Retain all original layout outside the pattern value itself."""
    before = yaml.safe_load(text)
    expected, count = hardened(before)
    lines: list[str] = []
    replacements = 0
    for line in text.splitlines(keepends=True):
        match = re.fullmatch(r"(\s*pattern:\s*)([^\n]+)(\n?)", line)
        if match is None:
            lines.append(line)
            continue
        prefix, raw, newline = match.groups()
        value = yaml.safe_load("pattern: " + raw)["pattern"]
        if not isinstance(value, str):
            raise ValueError("Expected a scalar string pattern")
        new_value = (value[:-1] + END_ASSERTION).replace("'", "''")
        lines.append(f"{prefix}'{new_value}'{newline}")
        replacements += 1
    result = "".join(lines)
    if replacements != count or yaml.safe_load(result) != expected:
        raise ValueError("Schema rewrite differs from the declared pattern-only change")
    return result, count


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    proposed = []
    for name in SCHEMAS:
        path = args.root / "schema" / f"{name}.yaml"
        result, count = rewrite(path.read_text(encoding="utf-8"))
        proposed.append((path, result, count))
    if sum(count for _, _, count in proposed) != 86:
        raise ValueError("Expected the reviewed 86 pattern declarations")
    for path, result, count in proposed:
        path.write_text(result, encoding="utf-8")
        print(f"{path.name}: {count} actual-end patterns")


if __name__ == "__main__":
    main()
