"""Copilot instructions retain root links and literal code (#1262, adapted from #245)."""

from __future__ import annotations

import importlib.util
import re
from pathlib import Path
from types import ModuleType
from typing import Any, SupportsIndex

import pytest
from markdown_it import MarkdownIt

from mtp import agent_docs
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


def commonmark_targets(source: str) -> tuple[list[str], list[str]]:
    """Collect rendered links and definitions independently of the rewriter."""
    environment: dict[str, Any] = {}
    tokens = MarkdownIt("commonmark").parse(source, environment)
    links: list[str] = []
    for token in tokens:
        for child in token.children or []:
            if child.type == "link_open":
                target = child.attrGet("href")
                assert isinstance(target, str)
                links.append(target)
    definitions = sorted(
        reference["href"] for reference in environment.get("references", {}).values()
    )
    return links, definitions


def assert_commonmark_rewrite(source: str, expected: str) -> None:
    """Pin literal output and independently confirm each intended destination."""
    original_links, original_definitions = commonmark_targets(source)
    assert commonmark_targets(expected) == (
        ["../" + target for target in original_links],
        ["../" + target for target in original_definitions],
    )
    assert links_relative_to_github(source) == expected


@pytest.mark.parametrize(
    ("source", "target"),
    [
        ("> ```\n> code\n> ```\n\nSee [y](docs/b.md).\n", "docs/b.md"),
        (
            "> Example:\n>\n> ```bash\n> just check\n> ```\n\n## Next\n\nSee [y](docs/b.md).\n",
            "docs/b.md",
        ),
        ("- a\n  - b\n    ```bash\n    code\n    ```\n- [y](docs/b.md)\n", "docs/b.md"),
        (
            "1. Run:\n\n    ```bash\n    just check\n    ```\n\n2. Read [guide](docs/guide.md).\n",
            "docs/guide.md",
        ),
        ("> - item\n>   ```\n>   [x](docs/a.md)\n>   ```\n> [y](docs/b.md)\n", "docs/b.md"),
    ],
    ids=["quoted", "quoted-shell-example", "nested-list", "ordered-list", "quoted-list"],
)
def test_review_container_fences_close_before_following_links(source: str, target: str) -> None:
    """#1373: opening and closing fences use the same container coordinates."""
    assert_commonmark_rewrite(source, source.replace(f"({target})", f"(../{target})"))


@pytest.mark.parametrize(
    ("source", "expected"),
    [
        (
            "> <!-- note\n\nSee [guide](docs/guide.md).\n",
            "> <!-- note\n\nSee [guide](../docs/guide.md).\n",
        ),
        (
            "> ```\n\nSee [guide](docs/guide.md).\n",
            "> ```\n\nSee [guide](../docs/guide.md).\n",
        ),
        ("- <!--\n[a](docs/a.md)\n", "- <!--\n[a](../docs/a.md)\n"),
        ("- ```\n[a](docs/a.md)\n", "- ```\n[a](../docs/a.md)\n"),
        ("- ```\n  code\n```\n[x](docs/a.md)\n", "- ```\n  code\n```\n[x](docs/a.md)\n"),
        (
            "- ```\n\n  [literal](docs/a.md)\n  ```\n\n[guide](docs/b.md)\n",
            "- ```\n\n  [literal](docs/a.md)\n  ```\n\n[guide](../docs/b.md)\n",
        ),
    ],
    ids=[
        "quote-comment-ends",
        "quote-fence-ends",
        "item-comment-ends",
        "item-fence-ends",
        "dedented-fence-opens-new-block",
        "blank-preserves-item-fence",
    ],
)
def test_review_fences_and_comments_end_with_their_container(source: str, expected: str) -> None:
    """#1374: container exit ends its block, but list-item blank lines do not."""
    assert_commonmark_rewrite(source, expected)


@pytest.mark.parametrize(
    ("source", "expected"),
    [
        (
            "- [Required]:\n```bash\n[x](docs/a.md)\n```\n",
            "- [Required]:\n```bash\n[x](docs/a.md)\n```\n",
        ),
        ("> [Note]:\n>\n> [x](docs/a.md)\n", "> [Note]:\n>\n> [x](../docs/a.md)\n"),
        ("- [Note]:\n---\n", "- [Note]:\n---\n"),
        (
            "- [Note]:\n# Heading\n\n[x](docs/a.md)\n",
            "- [Note]:\n# Heading\n\n[x](../docs/a.md)\n",
        ),
        (
            "- [Note]:\n<!-- comment -->\n\n[x](docs/a.md)\n",
            "- [Note]:\n<!-- comment -->\n\n[x](../docs/a.md)\n",
        ),
    ],
    ids=["fence", "empty-quote-line", "thematic-break", "heading", "html-comment"],
)
def test_review_reference_destination_cannot_consume_another_block(
    source: str, expected: str
) -> None:
    """#1375: a new block cannot supply an incomplete definition's destination."""
    assert_commonmark_rewrite(source, expected)


def test_review_lazy_quote_continuation_preserves_paragraph_state() -> None:
    """#1377: an explicit > after a lazy line continues the same paragraph."""
    source = "> Para\nlazy\n> [r]: docs/a.md\n> [r]\n"
    assert commonmark_targets(source) == ([], [])
    assert_commonmark_rewrite(source, source)


@pytest.mark.parametrize("spaces", [5, 6, 8])
def test_review_five_or_more_spaces_after_list_marker_are_indented_code(spaces: int) -> None:
    """#1378: pin the item-5plus-spaces-rule-dropped mutation."""
    source = "-" + " " * spaces + "[x](docs/a.md)\n"
    assert commonmark_targets(source) == ([], [])
    assert_commonmark_rewrite(source, source)


def test_review_entering_quote_starts_a_new_reference_paragraph() -> None:
    """#1378: pin the bq-entry-no-paragraph-reset mutation."""
    source = "Text\n> [r]: docs/a.md\n\n[r]\n"
    expected = "Text\n> [r]: ../docs/a.md\n\n[r]\n"
    assert commonmark_targets(source) == (["docs/a.md"], ["docs/a.md"])
    assert_commonmark_rewrite(source, expected)


@pytest.mark.parametrize(
    ("source", "expected"),
    [
        ("   >\t[](d)", "   >\t[](../d)"),
        ("0) x\n)\n0) [r]:d", "0) x\n)\n0) [r]:../d"),
        ("  -\n\t[](>)", "  -\n\t[](../>)"),
    ],
    ids=["quote-tab-content", "list-after-lazy-line", "empty-list-item-content-column"],
)
def test_review_surviving_mutant_witnesses_preserve_commonmark(source: str, expected: str) -> None:
    """#1378: pin reviewer minima where the existing branch matches CommonMark."""
    assert_commonmark_rewrite(source, expected)


@pytest.mark.parametrize(
    ("source", "expected"),
    [
        (">\n`\n>`[](d)`", ">\n`\n>`[](d)`"),
        (">\r`\n>[](d)`", ">\r`\n>[](../d)`"),
        ("-\tx\n>     [b]:d", "-\tx\n>     [b]:d"),
    ],
    ids=["new-quote-code-span", "new-quote-link", "new-quote-drops-list-indent"],
)
def test_review_quote_state_reentry_after_container_exit(source: str, expected: str) -> None:
    """#1377: only an actual lazy paragraph retains old container state."""
    assert_commonmark_rewrite(source, expected)


@pytest.mark.parametrize(
    ("source", "expected"),
    [
        (
            "   > - ```\n>   [literal](d)\n>   ```\n> [x](e)\n",
            "   > - ```\n>   [literal](d)\n>   ```\n> [x](../e)\n",
        ),
        (
            "   > - <!--\n>   [literal](d)\n>   -->\n> [x](e)\n",
            "   > - <!--\n>   [literal](d)\n>   -->\n> [x](../e)\n",
        ),
        ("- [r]:\n  <div>\n\n[r]\n", "- [r]:\n  <div>\n\n[r]\n"),
        (
            "- [r]: d\n  'title\noutside [x](e)'\n\n[r]\n",
            "- [r]: ../d\n  'title\noutside [x](e)'\n\n[r]\n",
        ),
        ("> [r]: d\n'[x](e)'\n[r]\n", "> [r]: ../d\n'[x](e)'\n[r]\n"),
        (
            "  <docs/auto.md>\n  >  `\r1. `\n >\ttext [b](docs/b.md) `x\n"
            "> [r]: docs/r.md 'title'\n",
            "  <docs/auto.md>\n  >  `\r1. `\n >\ttext [b](../docs/b.md) `x\n"
            "> [r]: docs/r.md 'title'\n",
        ),
        (">- '\n    ~~~\n \t[](d)", ">- '\n    ~~~\n \t[](../d)"),
        (">- )\n    1. [r]:'", ">- )\n    1. [r]:'"),
    ],
    ids=[
        "fence-variable-quote-prefix-indent",
        "comment-variable-quote-prefix-indent",
        "definition-destination-html-block",
        "list-reference-title-lazy-continuation",
        "quote-reference-title-lazy-continuation",
        "unquoted-list-ends-old-quote",
        "lazy-quote-line-is-not-fence",
        "lazy-quote-line-is-not-list-reference",
    ],
)
def test_review_container_coordinates_and_valid_definition_titles(
    source: str, expected: str
) -> None:
    """#1373-#1377: preserve actual containers and valid lazy reference titles."""
    assert_commonmark_rewrite(source, expected)


@pytest.mark.parametrize(
    ("destination", "expected_destination"),
    [
        ("<script>", "<script>"),
        ("<DIV>", "<DIV>"),
        ("<?pi?>", "<?pi?>"),
        ("<![CDATA[x]]>", "<![CDATA[x]]>"),
        ("<!DOCTYPE html>", "<!DOCTYPE html>"),
        ("<foo>", "<../foo>"),
        ("<docs/a>", "<../docs/a>"),
        ("<docs/a guide.md>", "<../docs/a guide.md>"),
    ],
    ids=[
        "script-html-block",
        "uppercase-div-html-block",
        "processing-instruction-block",
        "cdata-block",
        "doctype-block",
        "generic-tag-angle-destination",
        "path-angle-destination",
        "spaced-path-angle-destination",
    ],
)
def test_review_definition_continuation_distinguishes_html_blocks_and_angle_destinations(
    destination: str, expected_destination: str
) -> None:
    """#1375: only paragraph-interrupting HTML blocks forbid a destination."""
    source = f"- [r]:\n  {destination}\n\n[r]\n"
    expected = f"- [r]:\n  {expected_destination}\n\n[r]\n"
    assert_commonmark_rewrite(source, expected)


@pytest.mark.parametrize(
    ("source", "expected"),
    [
        (">`\n2. `[](d)`", ">`\n2. `[](d)`"),
        (">`\n2. [](d)`", ">`\n2. [](../d)`"),
        (">)\n  1.\n    [](d)", ">)\n  1.\n    [](d)"),
    ],
    ids=["non-one-ordered-code-span", "non-one-ordered-link", "empty-ordered-item-code"],
)
def test_review_unquoted_ordered_item_starts_after_quote_exit(source: str, expected: str) -> None:
    """#1377: a missing quote marker can start a new list, not a lazy paragraph."""
    assert_commonmark_rewrite(source, expected)


class MeasuredText(str):
    """Count copied characters while preserving tracking across string operations."""

    work: dict[str, int]

    def __new__(cls, value: str, work: dict[str, int]) -> MeasuredText:
        text = super().__new__(cls, value)
        text.work = work
        return text

    def __getitem__(self, key: SupportsIndex | slice) -> str:
        result = super().__getitem__(key)
        if isinstance(key, slice):
            self.work["copied"] += len(result)
            return MeasuredText(result, self.work)
        return result

    def __add__(self, other: str) -> str:
        result = super().__add__(other)
        self.work["copied"] += len(result)
        return MeasuredText(result, self.work)

    def __radd__(self, other: str) -> str:
        result = str.__add__(other, self)
        self.work["copied"] += len(result)
        return MeasuredText(result, self.work)


@pytest.mark.parametrize("count", [64, 128])
@pytest.mark.parametrize("prefix", ["", "> ", "- "])
def test_reference_definition_work_is_linear(
    count: int, prefix: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    """#1431: count actual parser operations and copies, not wall-clock timings."""
    source = "".join(f"{prefix}[r{i:04d}]: docs/{i:04d}.md\n" for i in range(count))
    work = {"line_ends": 0, "copied": 0}
    original = agent_docs._line_end

    def measured_line_end(*args: Any, **kwargs: Any) -> int:
        work["line_ends"] += 1
        return original(*args, **kwargs)

    monkeypatch.setattr(agent_docs, "_line_end", measured_line_end)
    assert links_relative_to_github(MeasuredText(source, work)) == source.replace(
        ": docs/", ": ../docs/"
    )
    assert work["line_ends"] <= 12 * count, work
    assert work["copied"] <= 20 * len(source), work


@pytest.mark.parametrize("count", [64, 128])
def test_definition_followed_by_long_paragraph_has_linear_copy_work(count: int) -> None:
    """#1431: a single forward walk must not copy a shrinking suffix at each line."""
    source = "[r]: docs/r.md\n" + "ordinary text line\n" * count
    work = {"copied": 0}
    assert links_relative_to_github(MeasuredText(source, work)) == source.replace(
        ": docs/", ": ../docs/"
    )
    assert work["copied"] <= 20 * len(source), work


@pytest.mark.parametrize("newline", ["\n", "\r\n"])
@pytest.mark.parametrize(
    "definition",
    [
        '[r]: docs/r.md\n  "next-line title"\n',
        '[r]: docs/r.md "multi\nline title"\n',
        '[r]: docs/r.md\n  "unfinished title\n',
        '[r]: docs/r.md\n  "title [literal](docs/literal.md)"\n',
    ],
)
def test_definition_bound_reuse_preserves_titles(definition: str, newline: str) -> None:
    """#1431: caching cannot truncate valid titles or change an unfinished fallback."""
    source = ("[first]: docs/first.md\n" + definition + "\n[r]\n").replace("\n", newline)
    expected = source.replace(": docs/", ": ../docs/")
    assert_commonmark_rewrite(source, expected)


def test_definition_bound_cache_distinguishes_pending_destination() -> None:
    """A cached continuation cannot bypass a pending destination's container guard."""
    source = '- [first]: docs/first.md\n  [pending]:\n"outside"\n'
    expected = source.replace(": docs/first.md", ": ../docs/first.md")
    assert links_relative_to_github(source) == expected


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_contributor_docs_agree_on_package_paths_and_check_before_commit() -> None:
    """#1341, #463: contributor docs state mtp package paths and check-before-commit."""
    repo_root = REPO_ROOT
    pyproject = (repo_root / "pyproject.toml").read_text(encoding="utf-8")
    agents = (repo_root / "AGENTS.md").read_text(encoding="utf-8")
    copilot = (repo_root / ".github" / "copilot-instructions.md").read_text(encoding="utf-8")
    contributing = (repo_root / "CONTRIBUTING.md").read_text(encoding="utf-8")
    readme = (repo_root / "README.md").read_text(encoding="utf-8")

    assert 'name = "mtp"' in pyproject

    for name, doc in (
        ("AGENTS.md", agents),
        (".github/copilot-instructions.md", copilot),
        ("CONTRIBUTING.md", contributing),
        ("README.md", readme),
    ):
        assert "before committing" in doc, f"{name} must state the check-before-commit rule"
        assert "before requesting review" not in doc, f"{name} has conflicting check timing"
        assert "before opening a PR" not in doc, f"{name} has conflicting check timing"

    assert "3. **Run `just check`** — all checks must pass before committing" in agents
    assert "4. **Commit, push, and open a PR.**" in agents
    assert "3. **Run `just check`** — all checks must pass before committing" in copilot
    assert "4. **Commit, push, and open a PR.**" in copilot
    assert "4. Run `just check` — **all checks must pass before committing**" in readme
    assert "### 4. Run checks before committing" in contributing
    assert "**All checks must pass before committing.**" in contributing

    assert "- Source code goes in `src/mtp/` (with compatibility shims in `src/bsdbng/`)" in (
        contributing
    )
    assert "- Source code goes in `src/bsdbng/`" not in contributing

    assert "The Python package remains `bsdbng`" not in readme
    assert "The Python package is `mtp` (`src/mtp/`, `schema/mtp.yaml`)" in readme
    assert "compatibility shims" in readme
    assert "schema/mtp.yaml" in readme
    assert "src/mtp/" in readme
    assert "datamodel/mtp_pydantic.py" in readme
    assert "1900+" not in readme
    assert "1,900+" not in readme
    assert "1,245 reviewed studies" not in readme
    assert "1,245 studies in pinned v1.3.1" not in readme
    assert "parse into per-study YAML for the pinned BugSigDB release (~22s)" in readme
    assert "load study YAMLs into MongoDB (~45s; for the pinned BugSigDB release)" in readme

    for name, doc in (
        ("AGENTS.md", agents),
        (".github/copilot-instructions.md", copilot),
        ("CONTRIBUTING.md", contributing),
    ):
        flat = " ".join(doc.split())
        assert "Do not use `gh issue develop`" in flat, (
            f"{name} missing gh issue develop prohibition"
        )
        assert re.search(r"(?<!Do not )(?<!do not )\b[Uu]se `gh issue develop`", flat) is None, (
            f"{name} permits gh issue develop"
        )
        assert "gh pr view N --json closingIssuesReferences" in flat, (
            f"{name} missing closingIssuesReferences command"
        )
        assert "which lists issues only, never PRs, or the PR's Development sidebar" in flat, (
            f"{name} missing closingIssuesReferences issues-only scope"
        )
        assert (
            "a link that comes only from a closing keyword goes away once the keyword is removed "
            "from the PR body"
        ) in flat, f"{name} missing closing-keyword removal rule"
        assert (
            "whereas a Development-sidebar or branch link persists and must be removed in the "
            "PR's Development sidebar in the web UI (there is no `gh` command, so an agent must "
            "ask a maintainer)"
        ) in flat, f"{name} missing Development-sidebar removal rule"
        assert (
            "A closing keyword in a commit message (and so in the squash message) also closes "
            "the issue on merge, appears in neither check, and is not removed by editing the PR "
            'body: write `Refs #<number>` in commits, and never a negated keyword such as "does '
            'not close #<number>" (which still closes).'
        ) in flat, f"{name} missing commit-message and negated-keyword rule"
        assert "Remove any link to an issue the PR only partly addresses." in flat, (
            f"{name} missing remove-partial-link rule"
        )
        assert "Keep any link" not in flat, f"{name} reverses remove-partial-link rule"
        for sentence in re.split(r"(?<=[.!?])\s+", flat):
            if ("sidebar" in sentence or "branch link" in sentence) and "goes away" in sentence:
                assert "only from a closing keyword goes away" in sentence, (
                    f"{name} claims sidebar/branch link goes away: {sentence!r}"
                )
                assert "whereas" in sentence, (
                    f"{name} claims sidebar/branch link goes away: {sentence!r}"
                )
