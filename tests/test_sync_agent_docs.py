"""Copilot instructions retain root links and literal code (#1262, adapted from #245)."""

from __future__ import annotations

import importlib.util
from pathlib import Path
from types import ModuleType

import pytest

from mtp.agent_docs import _destination, links_relative_to_github


def test_copilot_links_are_relative_to_github_directory() -> None:
    """Port Jed's #245 links while preserving parent targets from the new directory."""
    source = (
        "[guide](methods/guide.md) [local](./docs/local.md) "
        "[parent](../shared.md) [web](https://example.org/x) [heading](#scope)"
    )
    assert links_relative_to_github(source) == (
        "[guide](../methods/guide.md) [local](../docs/local.md) "
        "[parent](../../shared.md) [web](https://example.org/x) [heading](#scope)"
    )


@pytest.mark.parametrize(
    ("source", "expected"),
    [
        (
            '[guide](<methods/a guide.md#scope> "Guide")',
            '[guide](<../methods/a guide.md#scope> "Guide")',
        ),
        ("[guide](docs/a(b).md)", "[guide](../docs/a(b).md)"),
        (r"[guide](docs/a\(b\).md)", r"[guide](../docs/a\(b\).md)"),
        ("[guide](docs/a.md 'A title')", "[guide](../docs/a.md 'A title')"),
        ("[guide](docs/a.md (A title))", "[guide](../docs/a.md (A title))"),
        ("[guide](\n docs/a.md\n)", "[guide](\n ../docs/a.md\n)"),
        ("[guide](\n <docs/a guide.md>\n)", "[guide](\n <../docs/a guide.md>\n)"),
        ("![figure](images/a.svg)", "![figure](../images/a.svg)"),
        ("[`guide.md`](docs/a.md)", "[`guide.md`](../docs/a.md)"),
        ("[nested [label]](docs/a.md)", "[nested [label]](../docs/a.md)"),
        ("[guide][ref]\n\n[ref]: docs/a.md#scope\n", "[guide][ref]\n\n[ref]: ../docs/a.md#scope\n"),
        (
            '[ref]: <docs/a guide.md> "Title"\n[ref][]',
            '[ref]: <../docs/a guide.md> "Title"\n[ref][]',
        ),
        ("  [ref]: ./methods/a.md\n[ref]", "  [ref]: ../methods/a.md\n[ref]"),
        ("[ref]:\n  docs/a.md\n[ref]", "[ref]:\n  ../docs/a.md\n[ref]"),
        ("[ref]:\n    <docs/a guide.md>\n[ref]", "[ref]:\n    <../docs/a guide.md>\n[ref]"),
        (
            '[ref]: docs/a.md "[literal](docs/title.md)"\n[ref]',
            '[ref]: ../docs/a.md "[literal](docs/title.md)"\n[ref]',
        ),
        (
            '[ref]: docs/a.md\n  "[literal](docs/title.md)"\n[ref]',
            '[ref]: ../docs/a.md\n  "[literal](docs/title.md)"\n[ref]',
        ),
        (
            "[a](docs/a.md)\r\n\r\n[b]: methods/b.md\r\n",
            "[a](../docs/a.md)\r\n\r\n[b]: ../methods/b.md\r\n",
        ),
        (r"[guide](docs/a\)b.md)", r"[guide](../docs/a\)b.md)"),
        (r"[guide](<docs/a\>b.md>)", r"[guide](<../docs/a\>b.md>)"),
        ('[guide](docs/a.md "title" )', '[guide](../docs/a.md "title" )'),
        (
            r'[guide](docs/a.md "title \" [literal](docs/b.md)")',
            r'[guide](../docs/a.md "title \" [literal](docs/b.md)")',
        ),
        (
            '[guide](docs/a.md "[literal](docs/b.md)")',
            '[guide](../docs/a.md "[literal](docs/b.md)")',
        ),
        ("[ref]: docs/a(b).md\n", "[ref]: ../docs/a(b).md\n"),
        ("[ref]: docs/a.md", "[ref]: ../docs/a.md"),
        (
            '[ref]: docs/a.md "[literal](docs/b.md)"',
            '[ref]: ../docs/a.md "[literal](docs/b.md)"',
        ),
        (
            '[ref]: docs/a.md\n"quoted" [guide](docs/b.md)\n',
            '[ref]: ../docs/a.md\n"quoted" [guide](../docs/b.md)\n',
        ),
        (
            "[ref]: docs/a.md\n    [literal](docs/b.md)\n",
            "[ref]: ../docs/a.md\n    [literal](docs/b.md)\n",
        ),
        (r"\`not code [guide](docs/a.md) `later`", r"\`not code [guide](../docs/a.md) `later`"),
        ("` `` ` [guide](docs/a.md) ``", "` `` ` [guide](../docs/a.md) ``"),
        ("```info`string\n[guide](docs/a.md)\n", "```info`string\n[guide](../docs/a.md)\n"),
        ("    ```\n[guide](docs/a.md)\n", "    ```\n[guide](../docs/a.md)\n"),
        (
            "```\n[literal](docs/a.md)\n````\n[guide](docs/b.md)\n",
            "```\n[literal](docs/a.md)\n````\n[guide](../docs/b.md)\n",
        ),
    ],
)
def test_markdown_destination_forms_preserve_other_bytes(source: str, expected: str) -> None:
    assert links_relative_to_github(source) == expected


@pytest.mark.parametrize(
    "source",
    [
        "`[literal](docs/a.md)`",
        r"\\`[literal](docs/a.md)`",
        "``a ` [literal](docs/a.md)``",
        "`multi\n[literal](docs/a.md)\nline`",
        "[ref `code\rspan`]: docs/a.md\n",
        "```markdown\n[literal](docs/a.md)\n[ref]: docs/a.md\n```\n",
        "```\n[literal](docs/a.md)\n\n[also literal](docs/b.md)\n```\n",
        "```\n[literal](docs/a.md)\n~~~\n[still literal](docs/b.md)\n```\n",
        "```\n[literal](docs/a.md)\n    ```\n[still literal](docs/b.md)\n```\n",
        "```\n[literal](docs/a.md)\n``` trailing\n[still literal](docs/b.md)\n```\n",
        "~~~info`string\n[literal](docs/a.md)\n~~~\n",
        "~~~ [literal](docs/a.md)\n~~~\n",
        "~~~~\n[literal](docs/a.md)\n~~~\n[still literal](docs/b.md)\n~~~~\n",
        "    [literal](docs/a.md)\n\t[also literal](docs/b.md)\n",
        "    [ref]:\n    [literal](docs/a.md)\n",
        "[prev]:\n    [ref]: docs/a.md\n",
        "para\n    - x\n2. `code\n3. [literal](docs/a.md)\nend`",
        "1234567890. `code\n2. [literal](docs/a.md)\nend`",
        r"\[literal](docs/a.md)",
        r"[literal\](docs/a.md)",
        "not a link ](docs/a.md)",
        "[guide](<>)",
        "[incomplete](docs/a.md",
        "[incomplete](\n",
        "[invalid](docs/a(b)",
        "[invalid](docs/a<b.md)",
        "[invalid](<docs/a<b.md>)",
        "[invalid](<docs/a\nb.md>)",
        "[invalid](<docs/a\rb.md>)",
        "[invalid](\n\n docs/a.md)",
        "[invalid](docs/a.md\n\n)",
        '[invalid](docs/a.md "title"\n\n)',
        '[invalid](docs/a.md "title"junk)',
        '[invalid](docs/a.md "title"',
        '[web](<https://example.org/a(b)> "Title")',
        "[web](https://example.org/x) [mail](mailto:info@example.org)",
        (
            "[upper](HTTPS://example.org/x) [digit](h323:host) "
            "[plus](svn+ssh://host/x) [dot](z39.50r:host) [dash](custom-scheme:path)"
        ),
        "[root](/docs/a.md) [remote](//example.org/a) [fragment](#a) [query](?view=raw)",
        "[web]: https://example.org/x\n[anchor]: #scope\n[web] [anchor]",
        "[ref]: <docs/a.md",
        "[ref]: docs/a(b.md\n",
        '[ref]: docs/a.md "unclosed title',
        "[ref]: docs/a.md x",
        '[ref]: docs/a.md "Title" junk\n',
        "[ref]: docs/a.md stray words\n[ref]",
    ],
)
def test_code_nonlinks_and_nonlocal_targets_are_unchanged(source: str) -> None:
    assert links_relative_to_github(source) == source


def test_links_after_code_blocks_are_rebased() -> None:
    source = (
        "```md\n[literal](docs/literal.md)\n```\n"
        "[guide](docs/a.md) and `[literal](docs/literal.md)` then [next](methods/b.md)\n"
        "~~~\n[ref]: docs/literal.md\n~~~\n[ref]: docs/c.md\n"
    )
    expected = (
        source.replace("(docs/a.md)", "(../docs/a.md)")
        .replace("(methods/b.md)", "(../methods/b.md)")
        .replace("[ref]: docs/c.md", "[ref]: ../docs/c.md")
    )
    assert links_relative_to_github(source) == expected


@pytest.mark.parametrize("paragraph", ["\n\n", "\n \t\n", "\r\n\r\n"])
def test_unmatched_code_delimiter_cannot_hide_a_later_paragraph(paragraph: str) -> None:
    """#1267: CommonMark code spans cannot join separate block paragraphs."""
    source = f"`unclosed{paragraph}[guide](methods/guide.md){paragraph}`later`"
    expected = source.replace("(methods/guide.md)", "(../methods/guide.md)")
    assert links_relative_to_github(source) == expected


@pytest.mark.parametrize("newline", ["\n", "\r\n", "\r"], ids=["lf", "crlf", "cr"])
@pytest.mark.parametrize("boundary", ["heading", "list", "thematic", "fence"])
def test_unmatched_code_delimiter_stops_at_recognized_blocks(boundary: str, newline: str) -> None:
    """#1281: the next block's real link must not become a previous code span."""
    blocks = {
        "heading": ["# Heading [guide](docs/guide.md)"],
        "list": ["- [guide](docs/guide.md)"],
        "thematic": ["***", "[guide](docs/guide.md)"],
        "fence": ["~~~markdown", "[literal](docs/literal.md)", "~~~", "[guide](docs/guide.md)"],
    }
    source = newline.join(["`unclosed", *blocks[boundary], "`later`"])
    expected = source.replace("[guide](docs/guide.md)", "[guide](../docs/guide.md)")
    assert links_relative_to_github(source) == expected


def test_heading_code_delimiter_stops_at_the_heading_end() -> None:
    source = "# `unclosed\n[guide](docs/guide.md)\n`later`"
    expected = source.replace("(docs/guide.md)", "(../docs/guide.md)")
    assert links_relative_to_github(source) == expected


def test_heading_start_boundary_stops_prior_unmatched_code_delimiter() -> None:
    """#1312: removing heading from the start-boundary check must fail."""
    source = "`unclosed\n# Heading [guide](docs/guide.md) `code`"
    expected = source.replace("(docs/guide.md)", "(../docs/guide.md)")
    assert links_relative_to_github(source) == expected


def test_multiline_code_span_immediately_after_heading() -> None:
    """#1312: heading-end boundary at opening.start() must not cut the next span."""
    source = "# Heading\n`code\n[literal](docs/literal.md)\nend`"
    assert links_relative_to_github(source) == source


@pytest.mark.parametrize(
    "source",
    [
        "`unclosed\n* [guide](docs/guide.md)\n`later`",
        "`unclosed\n+ [guide](docs/guide.md)\n`later`",
        "`unclosed\n1. [guide](docs/guide.md)\n`later`",
        "`unclosed\n1) [guide](docs/guide.md)\n`later`",
        "`unclosed\n___\n[guide](docs/guide.md)\n`later`",
        "- `unclosed\n---\n[guide](docs/guide.md)\n`later`",
        "`unclosed\n###### Heading [guide](docs/guide.md)\n`later`",
        "`unclosed\n######\n[guide](docs/guide.md)\n`later`",
        "> `unclosed\n-\n[guide](docs/guide.md)\n`later`",
    ],
    ids=[
        "star-item",
        "plus-item",
        "ordered-dot",
        "ordered-paren",
        "underscore-break",
        "dash-break",
        "six-hash-heading",
        "empty-six-hash-heading",
        "empty-item-after-blockquote",
    ],
)
def test_list_and_thematic_marker_variants_stop_unmatched_code_delimiter(source: str) -> None:
    """#1312: each list and thematic marker branch must be exercised."""
    expected = source.replace("(docs/guide.md)", "(../docs/guide.md)")
    assert links_relative_to_github(source) == expected


@pytest.mark.parametrize("marker", ["2.", "10)", "2025."])
def test_non_interrupting_ordered_markers_preserve_multiline_code_spans(marker: str) -> None:
    """#1311: ordered markers other than 1 cannot interrupt a paragraph."""
    source = f"`code\n{marker} [literal](docs/literal.md)\nend`"
    assert links_relative_to_github(source) == source


@pytest.mark.parametrize("marker", ["1.", "1)", "*", "+"])
def test_empty_list_markers_preserve_multiline_code_spans(marker: str) -> None:
    """#1311: empty list markers cannot interrupt a paragraph."""
    source = f"`code\n{marker}\n[literal](docs/literal.md)\nend`"
    assert links_relative_to_github(source) == source


def test_ordered_list_context_distinguishes_paragraph_from_sibling() -> None:
    """#1311, #1312: non-1 ordered markers are boundaries in lists, after blanks, and fences."""
    source = (
        "Paragraph\n\n"
        "2. `unclosed\n"
        "3. [guide](docs/guide.md)\n"
        "`later`\n"
        "~~~\n"
        "code\n"
        "~~~\n"
        "2. `unclosed\n"
        "3. [next](methods/next.md)\n"
        "`later`\n\n"
        "`code\n"
        "2. [literal](docs/literal.md)\n"
        "end`"
    )
    expected = source.replace("(docs/guide.md)", "(../docs/guide.md)").replace(
        "(methods/next.md)", "(../methods/next.md)"
    )
    assert links_relative_to_github(source) == expected


@pytest.mark.parametrize("underline", ["===", "--", "=", "-"])
def test_setext_heading_underlines_stop_unmatched_code_delimiter(underline: str) -> None:
    """#1311: a setext heading underline is a block boundary."""
    source = f"Heading `unclosed\n{underline}\n[guide](docs/guide.md)\n`later`"
    expected = source.replace("(docs/guide.md)", "(../docs/guide.md)")
    assert links_relative_to_github(source) == expected


@pytest.mark.parametrize(
    "blockquote",
    ["> [guide](docs/guide.md)", "> quote\n[guide](docs/guide.md)"],
    ids=["direct", "lazy-continuation"],
)
def test_blockquote_stops_unmatched_code_delimiter(blockquote: str) -> None:
    """#1311: a blockquote interrupting a paragraph is a block boundary."""
    source = f"`unclosed\n{blockquote}\n`later`"
    expected = source.replace("(docs/guide.md)", "(../docs/guide.md)")
    assert links_relative_to_github(source) == expected


@pytest.mark.parametrize(
    "line",
    [
        "#not-a-heading",
        "####### not-a-heading",
        "-not-a-list",
        "**",
        "*** not-a-break",
        "[ref]:\n    # not-a-heading",
        "[ref]:\n    1. not-a-list",
        "[ref]:\n    > not-a-quote",
        "[ref]:\n    ===",
        "[ref]:\n    ***",
        "[ref]:\n    - not-a-list",
    ],
)
def test_nonblock_lines_preserve_multiline_code_spans(line: str) -> None:
    source = f"`code\n{line}\n[literal](docs/literal.md)\nend`"
    assert links_relative_to_github(source) == source


def test_destination_helper_rejects_empty_bare_and_bounds_trailing_backslash() -> None:
    """#1283: direct unit checks for _destination edge conditions masked by caller."""
    assert _destination("[guide]()", 8) is None
    assert _destination("docs/a.md\\", 0) == (0, 10, 10)


@pytest.mark.parametrize("continuation", ["  continuation", "==="])
def test_same_list_item_preserves_multiline_code_span(continuation: str) -> None:
    source = f"- `code\n{continuation}\n[literal](docs/literal.md)\n  end`"
    assert links_relative_to_github(source) == source


@pytest.mark.parametrize("paragraph", ["\n\n", "\n \t\n", "\r\n\r\n"])
def test_reference_destination_cannot_cross_a_blank_paragraph(paragraph: str) -> None:
    """#1268: retain both newlines when checking a reference destination boundary."""
    source = f"[ref]:{paragraph}  docs/a.md\n"
    assert links_relative_to_github(source) == source


@pytest.mark.parametrize(
    "source",
    [
        "- item\n    - [guide](docs/a.md)\n",
        "- item\n\t- [guide](docs/a.md)\n",
        "Paragraph\n    [guide](docs/a.md)\n",
        "Paragraph\n\t[guide](docs/a.md)\n",
        "1.  item\n\n    [guide](docs/a.md)\n",
        "1.\titem\n\n\t[guide](docs/a.md)\n",
        "    - x\n2. `code\n3. [guide](docs/a.md)\nend`",
    ],
    ids=[
        "nested-list-4-spaces",
        "nested-list-tab",
        "lazy-continuation-4-spaces",
        "lazy-continuation-tab",
        "item-paragraph-4-spaces",
        "item-paragraph-tab",
        "list-after-indented-code",
    ],
)
def test_indented_nested_lists_lazy_continuations_and_item_paragraphs_rebase_links(
    source: str,
) -> None:
    """#1282: 4-space or tab indentation in list items and lazy continuations is not code."""
    expected = source.replace("(docs/a.md)", "(../docs/a.md)")
    assert links_relative_to_github(source) == expected


@pytest.mark.parametrize(
    ("source", "expected"),
    [
        ("> [ref]: docs/a.md\n> [ref]\n", "> [ref]: ../docs/a.md\n> [ref]\n"),
        ("- [ref]: docs/a.md\n  [ref]\n", "- [ref]: ../docs/a.md\n  [ref]\n"),
        (
            r"[a\]b]: docs/a.md" + "\n" + r"[a\]b]" + "\n",
            r"[a\]b]: ../docs/a.md" + "\n" + r"[a\]b]" + "\n",
        ),
    ],
    ids=["blockquote-ref", "list-item-ref", "escaped-bracket-label"],
)
def test_reference_definitions_in_containers_and_with_escaped_brackets_rebase_links(
    source: str, expected: str
) -> None:
    """#1282: reference definitions inside > or - containers and \\[ / \\] labels rebase."""
    assert links_relative_to_github(source) == expected


@pytest.mark.parametrize(
    "separator",
    ["\u2028", "\x0c", "\x85"],
    ids=["line-separator-u2028", "form-feed-x0c", "nel-x85"],
)
def test_non_commonmark_line_separators_before_fence_markers_do_not_hide_links(
    separator: str,
) -> None:
    """#1282: CommonMark does not split lines on U+2028, form feed, or NEL."""
    source = f"text{separator}```\n[guide](docs/a.md)\n"
    expected = source.replace("(docs/a.md)", "(../docs/a.md)")
    assert links_relative_to_github(source) == expected


@pytest.mark.parametrize(
    "source",
    [
        "Paragraph\n[ref]: docs/a.md\n[ref]\n",
        "> Paragraph\n> [ref]: docs/a.md\n> [ref]\n",
        "- Paragraph\n  [ref]: docs/a.md\n  [ref]\n",
    ],
    ids=["top-level-paragraph", "blockquote-paragraph", "list-item-paragraph"],
)
def test_reference_like_line_interrupting_paragraph_is_unchanged(source: str) -> None:
    """#1282: a reference definition cannot interrupt a paragraph."""
    assert links_relative_to_github(source) == source


@pytest.mark.parametrize(
    ("source", "expected"),
    [
        ("[[inner](docs/b.md)](docs/a.md)", "[[inner](../docs/b.md)](docs/a.md)"),
        (
            "[outer [inner](docs/b.md) text](docs/a.md)",
            "[outer [inner](../docs/b.md) text](docs/a.md)",
        ),
        ("[![alt](images/b.svg)](docs/a.md)", "[![alt](../images/b.svg)](../docs/a.md)"),
        ("![[inner](docs/b.md)](images/a.svg)", "![[inner](../docs/b.md)](../images/a.svg)"),
    ],
    ids=["nested-link", "outer-text-nested-link", "image-in-link", "link-in-image"],
)
def test_nested_link_rebases_only_inner_target_while_preserving_image_combinations(
    source: str, expected: str
) -> None:
    """#1282: CommonMark forbids links inside links; outer (a) is literal text."""
    assert links_relative_to_github(source) == expected


@pytest.mark.parametrize(
    "source",
    [
        "<!-- [literal](docs/a.md) -->\n",
        "<!--\n[literal](docs/a.md)\n[ref]: docs/b.md\n-->\n",
        "Paragraph <!-- [literal](docs/a.md) --> text\n",
    ],
    ids=["single-line-block-comment", "multiline-block-comment", "inline-comment"],
)
def test_html_comments_are_unchanged(source: str) -> None:
    """#1282: links inside HTML comments are not rebased."""
    assert links_relative_to_github(source) == source


@pytest.mark.parametrize(
    ("source", "expected"),
    [
        ("[guide](docs/a.md (title (nested)))\n", "[guide](docs/a.md (title (nested)))\n"),
        ("[ref]: docs/a.md (title (nested)\n[ref]\n", "[ref]: docs/a.md (title (nested)\n[ref]\n"),
        (
            r"[guide](docs/a.md (title \(escaped\)))",
            r"[guide](../docs/a.md (title \(escaped\)))",
        ),
    ],
    ids=["inline-nested-parens", "reference-nested-parens", "inline-escaped-parens"],
)
def test_parenthesised_titles_reject_unescaped_nested_parentheses(
    source: str, expected: str
) -> None:
    """#1282: unescaped parentheses inside (...) link titles are invalid in CommonMark."""
    assert links_relative_to_github(source) == expected


@pytest.mark.parametrize(
    ("source", "expected"),
    [
        ("[invalid](\r\r docs/a.md)", "[invalid](\r\r docs/a.md)"),
        ("[invalid](docs/a.md\r\r)", "[invalid](docs/a.md\r\r)"),
        ('[invalid](docs/a.md "title"\r\r)', '[invalid](docs/a.md "title"\r\r)'),
        ("[ref]:\r\r  docs/a.md\n[ref]\n", "[ref]:\r\r  docs/a.md\n[ref]\n"),
        (
            "`unclosed\r\r[guide](docs/a.md)\r\r`later`",
            "`unclosed\r\r[guide](../docs/a.md)\r\r`later`",
        ),
    ],
    ids=[
        "inline-before-dest-cr-blank",
        "inline-after-dest-cr-blank",
        "inline-after-title-cr-blank",
        "ref-dest-cr-blank",
        "code-span-cr-blank",
    ],
)
def test_cr_only_blank_lines_break_paragraphs_and_reject_multiline_destinations(
    source: str, expected: str
) -> None:
    """#1282: blank lines ending in CR only are paragraph breaks in CommonMark."""
    assert links_relative_to_github(source) == expected


def load_sync_script() -> ModuleType:
    path = Path(__file__).resolve().parents[1] / "scripts" / "sync_agent_docs.py"
    spec = importlib.util.spec_from_file_location("sync_agent_docs_test", path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_existing_copilot_generation_is_unchanged() -> None:
    """#1281 must not require new generated bytes for the current AGENTS source."""
    module = load_sync_script()
    agents = module.AGENTS_MD.read_text(encoding="utf-8")
    expected = module.COPILOT_HEADER + links_relative_to_github(agents) + module.PROJECT_CONVENTIONS
    assert module.COPILOT_MD.read_text(encoding="utf-8") == expected


def test_sync_entrypoint_rebases_links_and_preserves_source(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source = tmp_path / "AGENTS.md"
    output = tmp_path / "copilot-instructions.md"
    raw = "[guide](<methods/guide.md>)\n`[example](docs/example.md)`\n\n[ref]: docs/ref.md\n"
    source.write_text(raw, encoding="utf-8")
    module = load_sync_script()
    monkeypatch.setattr(module, "AGENTS_MD", source)
    monkeypatch.setattr(module, "COPILOT_MD", output)
    module.main()
    assert output.read_text(encoding="utf-8") == (
        module.COPILOT_HEADER
        + "[guide](<../methods/guide.md>)\n`[example](docs/example.md)`\n\n[ref]: ../docs/ref.md\n"
        + module.PROJECT_CONVENTIONS
    )
    assert source.read_text(encoding="utf-8") == raw
