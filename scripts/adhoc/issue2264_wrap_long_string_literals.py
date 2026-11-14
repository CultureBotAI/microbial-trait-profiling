#!/usr/bin/env python3
"""Split over-long string literals into adjacent literals without changing their value.

Purpose: the NCASM plate builders ported under ``src/mtp/ncasm_tufte_plates/`` carry
several hundred text lines longer than the 100 columns ruff allows (E501). This rewrites
each such literal as adjacent literals on separate lines, which Python joins back into
the same string, and writes en dashes as ``\\u2013`` (RUF001). A file is rewritten only
when ``ast.dump`` of the result equals ``ast.dump`` of the input, so no string value,
name or call can change.
Inputs: Python files. Outputs: the same files, rewritten in place; run ``ruff format``
afterwards. Docstrings and comments are left for hand editing.
Issue #2264; 2026-11-11.
"""

from __future__ import annotations

import argparse
import ast
import io
import re
import sys
import tokenize
from pathlib import Path

LIMIT = 100
UNIT = re.compile(
    r"\\(?:u[0-9a-fA-F]{4}|U[0-9a-fA-F]{8}|x[0-9a-fA-F]{2}|N\{[^}]*\}|[0-7]{1,3}|.)"
    r"|\{\{|\}\}|.",
    re.DOTALL,
)
LITERAL = re.compile(r'(?P<prefix>[fF]?)"(?P<body>(?:[^"\\]|\\.)*)"\Z', re.DOTALL)


def _pieces(prefix: str, body: str, width: int) -> list[str] | None:
    """The literal's source split into chunks of at most ``width`` characters.

    Breaks fall after a space where one is available, never inside an escape sequence
    and never inside the braces of an f-string field."""
    units = UNIT.findall(body)
    chunks: list[str] = []
    current: list[str] = []
    length = 0
    last_break = -1  # index in ``current`` after which a break is preferred
    depth = 0
    for unit in units:
        if length + len(unit) > width and current:
            cut = last_break + 1 if last_break >= 0 else len(current)
            if depth > 0 and last_break < 0:
                return None
            chunks.append("".join(current[:cut]))
            current = current[cut:]
            length = sum(len(part) for part in current)
            last_break = -1
        current.append(unit)
        length += len(unit)
        if prefix and unit == "{":
            depth += 1
        elif prefix and unit == "}":
            depth -= 1
        elif unit == " " and depth == 0:
            last_break = len(current) - 1
    if current:
        chunks.append("".join(current))
    if any(len(chunk) > width for chunk in chunks):
        return None
    return chunks


def rewrite(source: str) -> str:
    """``source`` with every over-long single-line string literal split."""
    lines = source.split("\n")
    depth = 0
    candidates: dict[int, tuple[int, int, int]] = {}  # row -> (start col, end col, depth)
    for token in tokenize.generate_tokens(io.StringIO(source).readline):
        if token.type == tokenize.OP:
            if token.string in "([{":
                depth += 1
            elif token.string in ")]}":
                depth -= 1
        if token.type != tokenize.STRING or token.start[0] != token.end[0]:
            continue
        row = token.start[0]
        if len(lines[row - 1]) <= LIMIT:
            continue
        known = candidates.get(row)
        if known is None or token.end[1] - token.start[1] > known[1] - known[0]:
            candidates[row] = (token.start[1], token.end[1], depth)

    for row in sorted(candidates, reverse=True):
        start, end, nesting = candidates[row]
        line = lines[row - 1]
        match = LITERAL.match(line[start:end])
        if match is None or nesting == 0:
            continue
        prefix, body = match["prefix"], match["body"]
        indent = len(line) - len(line.lstrip())
        before, after = line[:start], line[end:]
        pure = not before.strip() and after.strip() in ("", ",")
        piece_indent = indent if pure else indent + 4
        width = LIMIT - piece_indent - len(prefix) - 2 - (1 if pure else 0)
        chunks = _pieces(prefix, body, width)
        if chunks is None or len(chunks) < 2:
            continue
        pad = " " * piece_indent
        literals = [f'{pad}{prefix}"{chunk}"' for chunk in chunks]
        if pure:
            literals[-1] += after
            lines[row - 1 : row] = literals
        else:
            lines[row - 1 : row] = [f"{before}(", *literals, f"{' ' * indent}){after}"]
    return "\n".join(lines)


def _escape_en_dashes(source: str) -> str:
    """En dashes inside ordinary (non-raw) string literals written as ``\\u2013``."""
    out: list[str] = []
    position = (1, 0)
    lines = source.split("\n")

    def between(begin: tuple[int, int], finish: tuple[int, int]) -> str:
        if begin[0] == finish[0]:
            return lines[begin[0] - 1][begin[1] : finish[1]]
        parts = [lines[begin[0] - 1][begin[1] :]]
        parts.extend(lines[begin[0] : finish[0] - 1])
        parts.append(lines[finish[0] - 1][: finish[1]])
        return "\n".join(parts)

    for token in tokenize.generate_tokens(io.StringIO(source).readline):
        out.append(between(position, token.start))
        text = between(token.start, token.end)
        raw = text[: len(text) - len(text.lstrip("rRbBuUfF"))].lower()
        if token.type == tokenize.STRING and "r" not in raw and "b" not in raw:
            text = text.replace("\u2013", "\\u2013")
        out.append(text)
        position = token.end
    return "".join(out)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("files", nargs="+", type=Path)
    args = parser.parse_args(argv)
    status = 0
    for path in args.files:
        source = path.read_text(encoding="utf-8")
        result = rewrite(_escape_en_dashes(source))
        if ast.dump(ast.parse(result)) != ast.dump(ast.parse(source)):
            print(f"{path}: the rewrite would change the syntax tree; left unchanged")
            status = 1
            continue
        if result != source:
            path.write_text(result, encoding="utf-8")
        long_before = sum(len(line) > LIMIT for line in source.split("\n"))
        long_after = sum(len(line) > LIMIT for line in result.split("\n"))
        print(f"{path}: lines over {LIMIT} columns {long_before} -> {long_after}")
    return status


if __name__ == "__main__":
    sys.exit(main())
