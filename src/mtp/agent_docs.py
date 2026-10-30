"""Repository instruction links, adapted from Jed Dongjin Kim-Ozaeta's PR #245 (#1262)."""

from __future__ import annotations

import re


def _escaped(text: str, position: int) -> bool:
    start = position
    while start and text[start - 1] == "\\":
        start -= 1
    return (position - start) % 2 == 1


def _outside_code(content: str) -> str:
    """Mask code without changing offsets, line endings, or the eventual output."""
    masked = list(content)

    def hide(start: int, end: int) -> None:
        for index in range(start, end):
            if masked[index] not in "\r\n":
                masked[index] = " "

    fence = ""
    position = 0
    reference_continuation = False
    for line in content.splitlines(keepends=True):
        marker = re.match(r" {0,3}(`{3,}|~{3,})(.*)", line.rstrip("\r\n"))
        if fence:
            hide(position, position + len(line))
            if re.fullmatch(
                r" {0,3}" + re.escape(fence[0]) + "{" + str(len(fence)) + r",}[ \t]*",
                line.rstrip("\r\n"),
            ):
                fence = ""
        elif marker and not (marker[1][0] == "`" and "`" in marker[2]):
            fence = marker[1]
            hide(position, position + len(line))
        elif line.startswith(("    ", "\t")) and not reference_continuation:
            hide(position, position + len(line))
        reference_continuation = not fence and bool(
            re.fullmatch(r" {0,3}\[[^\]\r\n]+\]:[ \t]*", line.rstrip("\r\n"))
        )
        position += len(line)
    code = "".join(masked)
    ticks = list(re.finditer(r"`+", code))
    index = 0
    while index < len(ticks):
        opening = ticks[index]
        if not _escaped(code, opening.start()):
            boundary = re.search(r"\n[ \t\r]*\n", code[opening.end() :])
            paragraph_end = opening.end() + boundary.start() if boundary else len(code)
            closing = next(
                (
                    item
                    for item in range(index + 1, len(ticks))
                    if ticks[item].start() < paragraph_end and ticks[item][0] == opening[0]
                ),
                None,
            )
            if closing is not None:
                hide(opening.start(), ticks[closing].end())
                index = closing + 1
                continue
        index += 1
    return "".join(masked)


def _space_end(text: str, position: int) -> int | None:
    start = position
    while position < len(text) and text[position].isspace():
        position += 1
    return None if re.search(r"\n[ \t\r]*\n", text[start:position]) else position


def _destination(text: str, start: int) -> tuple[int, int, int] | None:
    """Return destination bounds and the first following character, preserving <>."""
    whitespace_end = _space_end(text, start)
    if whitespace_end is None:
        return None
    start = whitespace_end
    if start == len(text):
        return None
    angle = text[start] == "<"
    start += int(angle)
    position, depth = start, 0
    while position < len(text):
        char = text[position]
        if char == "\\" and position + 1 < len(text):
            position += 2
            continue
        if angle:
            if char == ">":
                return start, position, position + 1
            if char in "\r\n<":
                return None
        elif char.isspace() or (char == ")" and depth == 0):
            break
        elif char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
        elif char == "<":
            return None
        position += 1
    if angle or depth or start == position:
        return None
    return start, position, position


def _title_end(text: str, position: int) -> int | None:
    if position == len(text) or text[position] not in "\"'(":
        return None
    end = ")" if text[position] == "(" else text[position]
    position += 1
    while position < len(text) and (text[position] != end or _escaped(text, position)):
        position += 1
    return position + 1 if position < len(text) else None


def _inline_end(text: str, position: int) -> int | None:
    """Require the enclosing closing parenthesis, allowing an optional link title."""
    whitespace_end = _space_end(text, position)
    if whitespace_end is None:
        return None
    position = whitespace_end
    title_end = _title_end(text, position)
    if title_end is not None:
        whitespace_end = _space_end(text, title_end)
        if whitespace_end is None:
            return None
        position = whitespace_end
    return position + 1 if position < len(text) and text[position] == ")" else None


def _line_end(text: str, position: int) -> int:
    end = text.find("\n", position)
    return end if end != -1 else len(text)


def _reference_end(text: str, position: int) -> int | None:
    """Validate the rest of a definition; titles may continue on the next line."""
    line_end = _line_end(text, position)
    on_next_line = not text[position:line_end].strip()
    if on_next_line:
        position = min(line_end + 1, len(text))
    while position < len(text) and text[position] in " \t\r":
        position += 1
    title_end = _title_end(text, position)
    if title_end is None:
        return line_end if on_next_line else None
    end = _line_end(text, title_end)
    if text[title_end:end].strip():
        return line_end if on_next_line else None
    return end


def _targets(content: str, visible: str) -> list[tuple[int, int]]:
    targets: list[tuple[int, int]] = []
    masked = list(visible)
    for reference in re.finditer(r"(?m)^ {0,3}\[[^\]\r\n]+\]:[ \t]*", visible):
        destination = _destination(content, reference.end())
        if destination is not None:
            end = _reference_end(content, destination[2])
            if end is not None:
                targets.append(destination[:2])
                # Link-like text inside a definition's title is literal.
                for index in range(reference.start(), end):
                    if masked[index] not in "\r\n":
                        masked[index] = " "
    visible = "".join(masked)
    brackets: list[int] = []
    position = 0
    while position < len(visible):
        char = visible[position]
        if char == "[" and not _escaped(visible, position):
            brackets.append(position)
        elif char == "]" and brackets and not _escaped(visible, position):
            brackets.pop()
            if visible[position + 1 : position + 2] == "(":
                destination = _destination(content, position + 2)
                if destination is not None:
                    start, end, after = destination
                    closing = _inline_end(content, after)
                    if closing is not None:
                        targets.append((start, end))
                        position = closing
                        continue
        position += 1
    return targets


def links_relative_to_github(content: str) -> str:
    """Rebase AGENTS-root Markdown destinations into `.github`, leaving code intact.

    Inline links/images and reference definitions retain their labels, titles,
    angle brackets and exact surrounding bytes. Schemes, root URLs and local
    fragments/queries remain unchanged. Parent-relative paths gain one parent
    too: this transforms source text once, rather than treating it as idempotent.
    """
    rewritten = content
    for start, end in sorted(set(_targets(content, _outside_code(content))), reverse=True):
        target = content[start:end]
        if (
            target
            and not target.startswith(("#", "/", "?"))
            and not re.match(r"[A-Za-z][A-Za-z0-9+.-]*:", target)
        ):
            rewritten = rewritten[:start] + "../" + target.removeprefix("./") + rewritten[end:]
    return rewritten
