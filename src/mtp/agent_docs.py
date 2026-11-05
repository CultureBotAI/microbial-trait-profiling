"""Repository instruction links, adapted from Jed Dongjin Kim-Ozaeta's PR #245 (#1262)."""

from __future__ import annotations

import re

_NEWLINE = r"(?:\r\n|\n|\r(?!\n))"
_BLANK_LINE = rf"{_NEWLINE}[ \t]*{_NEWLINE}"
_REF_LABEL = r"(\[(?:\\[^\r\n]|[^\\\[\]\r\n])+\]):[ \t]*"


def _escaped(text: str, position: int) -> bool:
    start = position
    while start and text[start - 1] == "\\":
        start -= 1
    return (position - start) % 2 == 1


def _advance_indent(text: str, index: int, column: int) -> tuple[int, int]:
    while index < len(text):
        if text[index] == " ":
            column += 1
        elif text[index] == "\t":
            column += 4 - (column % 4)
        else:
            break
        index += 1
    return index, column


def _outside_code(content: str) -> str:
    """Mask code without changing offsets, line endings, or the eventual output."""
    masked = list(content)

    def hide(start: int, end: int) -> None:
        for index in range(start, end):
            if masked[index] not in "\r\n":
                masked[index] = " "

    fence = ""
    html_comment = False
    position = 0
    reference_end_offset = 0
    in_paragraph = False
    quote_depth = 0
    list_cols: list[int] = []
    block_boundaries: set[int] = set()
    for line in re.findall(r"[^\r\n]*(?:\r\n|\n|\r)|[^\r\n]+", content):
        text = line.rstrip("\r\n")
        if fence:
            hide(position, position + len(line))
            if re.fullmatch(
                r" {0,3}" + re.escape(fence[0]) + "{" + str(len(fence)) + r",}[ \t]*",
                text,
            ):
                fence = ""
            position += len(line)
            continue
        if html_comment:
            hide(position, position + len(line))
            html_comment = "-->" not in text
            position += len(line)
            continue
        index, column = _advance_indent(text, 0, 0)
        base_col = 0
        line_quotes = 0
        while column - base_col <= 3 and index < len(text) and text[index] == ">":
            line_quotes += 1
            index += 1
            column += 1
            if index < len(text) and text[index] == " ":
                index += 1
                column += 1
            elif index < len(text) and text[index] == "\t":
                column += 4 - (column % 4)
                index += 1
            base_col = column
            index, column = _advance_indent(text, index, column)
        if line_quotes != quote_depth:
            list_cols = []
        kept_cols = (
            [col for col in list_cols if col <= column] if index < len(text) else list(list_cols)
        )
        dedented = len(kept_cols) < len(list_cols)
        if kept_cols:
            base_col = kept_cols[-1]
        if line_quotes > quote_depth:
            block_boundaries.add(position)
            in_paragraph = False
        was_in_paragraph = in_paragraph
        same_container = (line_quotes == quote_depth) and not dedented
        container_paragraph = in_paragraph and same_container
        rest = text[index:]
        indent = column - base_col
        setext = container_paragraph and bool(re.fullmatch(r"(?:=+|-+)[ \t]*", rest))
        thematic = not setext and bool(
            re.fullmatch(r"(?:(?:\*[ \t]*){3,}|(?:-[ \t]*){3,}|(?:_[ \t]*){3,})", rest)
        )
        had_item = False
        if indent <= 3 and not setext and not thematic and position >= reference_end_offset:
            item = (
                re.match(r"([-+*]|1[.)])(?=[ \t]+\S)", rest)
                if container_paragraph
                else re.match(r"([-+*]|[0-9]{1,9}[.)])(?=[ \t]|$)", rest)
            )
            if item:
                had_item = True
                block_boundaries.add(position)
                in_paragraph = False
                container_paragraph = False
                list_cols = kept_cols
                marker_end = index + len(item[1])
                marker_col = column + len(item[1])
                if marker_end == len(text):
                    base_col = marker_col + 1
                    index, column = marker_end, marker_col
                else:
                    after_idx, after_col = _advance_indent(text, marker_end, marker_col)
                    base_col = after_col if after_col - marker_col <= 4 else marker_col + 1
                    index, column = after_idx, after_col
                list_cols.append(base_col)
                kept_cols = list(list_cols)
                rest = text[index:]
                indent = column - base_col
                thematic = bool(
                    re.fullmatch(r"(?:(?:\*[ \t]*){3,}|(?:-[ \t]*){3,}|(?:_[ \t]*){3,})", rest)
                )
        if line_quotes > 0 or had_item:
            hide(position, position + index)
        marker = re.match(r"(`{3,}|~{3,})(.*)", rest) if indent <= 3 else None
        comment_start = indent <= 3 and rest.startswith("<!--")
        heading = (
            re.match(r"#{1,6}(?:[ \t]|$)", rest)
            if indent <= 3 and position >= reference_end_offset
            else None
        )
        thematic = indent <= 3 and thematic and position >= reference_end_offset
        setext = indent <= 3 and setext and position >= reference_end_offset
        if marker and not (marker[1][0] == "`" and "`" in marker[2]):
            block_boundaries.add(position)
            fence = marker[1]
            in_paragraph = False
            list_cols = kept_cols
            hide(position, position + len(line))
        elif comment_start:
            block_boundaries.add(position)
            in_paragraph = False
            list_cols = kept_cols
            hide(position, position + len(line))
            html_comment = "-->" not in rest
        elif heading or thematic or setext:
            block_boundaries.add(position)
            in_paragraph = False
            list_cols = kept_cols
            if heading:
                # An ATX heading is one line; its code spans cannot continue
                # into the following paragraph. List-item continuations can.
                block_boundaries.add(position + len(line))
        elif not rest or position < reference_end_offset:
            in_paragraph = False
        elif indent >= 4 and not in_paragraph:
            list_cols = kept_cols
            hide(position, position + len(line))
        else:
            ref_match = re.match(_REF_LABEL, rest)
            if ref_match:
                if indent <= 3 and not in_paragraph:
                    dest = _destination(content, position + index + ref_match.end())
                    ref_end = _reference_end(content, dest[2]) if dest is not None else None
                    if ref_end is not None:
                        reference_end_offset = ref_end
                        list_cols = kept_cols
                    else:
                        masked[position + index + ref_match.end(1)] = " "
                        in_paragraph = True
                else:
                    masked[position + index + ref_match.end(1)] = " "
                    in_paragraph = True
            else:
                in_paragraph = True
            if not was_in_paragraph or not dedented:
                list_cols = kept_cols
        quote_depth = line_quotes
        position += len(line)
    code = "".join(masked)
    tokens = list(re.finditer(r"`+|<!--", code))
    index = 0
    while index < len(tokens):
        opening = tokens[index]
        if not _escaped(code, opening.start()):
            boundary = re.search(_BLANK_LINE, code[opening.end() :])
            paragraph_end = opening.end() + boundary.start() if boundary else len(code)
            block_end = min(
                (offset for offset in block_boundaries if offset > opening.start()),
                default=len(code),
            )
            limit = min(paragraph_end, block_end)
            if opening[0] == "<!--":
                comment_end = code.find("-->", opening.end(), limit)
                if comment_end != -1:
                    hide(opening.start(), comment_end + 3)
                    while index < len(tokens) and tokens[index].start() < comment_end + 3:
                        index += 1
                    continue
            else:
                closing = next(
                    (
                        item
                        for item in range(index + 1, len(tokens))
                        if tokens[item].start() < limit and tokens[item][0] == opening[0]
                    ),
                    None,
                )
                if closing is not None:
                    hide(opening.start(), tokens[closing].end())
                    index = closing + 1
                    continue
        index += 1
    return "".join(masked)


def _space_end(text: str, position: int) -> int | None:
    start = position
    while position < len(text) and text[position] in " \t\r\n":
        position += 1
    return None if re.search(_BLANK_LINE, text[start:position]) else position


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
        elif char in " \t\r\n" or (char == ")" and depth == 0):
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
    opening = text[position]
    end = ")" if opening == "(" else opening
    start = position + 1
    position = start
    while position < len(text) and (text[position] != end or _escaped(text, position)):
        if opening == "(" and text[position] == "(" and not _escaped(text, position):
            return None
        position += 1
    if position == len(text) or re.search(_BLANK_LINE, text[start:position]):
        return None
    return position + 1


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
    match = re.search(r"\r\n|\n|\r", text[position:])
    return position + match.start() if match else len(text)


def _reference_end(text: str, position: int) -> int | None:
    """Validate the rest of a definition; titles may continue on the next line."""
    line_end = _line_end(text, position)
    on_next_line = not text[position:line_end].strip()
    if on_next_line:
        step = 2 if text[line_end : line_end + 2] == "\r\n" else 1
        position = min(line_end + step, len(text))
    while position < len(text) and text[position] in " \t":
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
    for reference in re.finditer(rf"(?m)^[ \t]*{_REF_LABEL}", visible):
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
    brackets: list[tuple[bool, bool]] = []
    position = 0
    while position < len(visible):
        char = visible[position]
        if char == "[" and not _escaped(visible, position):
            is_image = (
                position > 0
                and visible[position - 1] == "!"
                and not _escaped(visible, position - 1)
            )
            brackets.append((is_image, True))
        elif char == "]" and brackets and not _escaped(visible, position):
            is_image, active = brackets.pop()
            if active and visible[position + 1 : position + 2] == "(":
                destination = _destination(content, position + 2)
                if destination is not None:
                    start, end, after = destination
                    closing = _inline_end(content, after)
                    if closing is not None:
                        targets.append((start, end))
                        if not is_image:
                            brackets = [(img, False if not img else act) for img, act in brackets]
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
