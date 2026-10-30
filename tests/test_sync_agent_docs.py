"""Copilot instructions retain root links and literal code (#1262, adapted from #245)."""

from __future__ import annotations

import importlib.util
from pathlib import Path
from types import ModuleType

import pytest

from mtp.agent_docs import links_relative_to_github


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
        ("[guide][ref]\n[ref]: docs/a.md#scope\n", "[guide][ref]\n[ref]: ../docs/a.md#scope\n"),
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
            "[a](docs/a.md)\r\n[b]: methods/b.md\r\n",
            "[a](../docs/a.md)\r\n[b]: ../methods/b.md\r\n",
        ),
    ],
)
def test_markdown_destination_forms_preserve_other_bytes(source: str, expected: str) -> None:
    assert links_relative_to_github(source) == expected


@pytest.mark.parametrize(
    "source",
    [
        "`[literal](docs/a.md)`",
        "``a ` [literal](docs/a.md)``",
        "`multi\n[literal](docs/a.md)\nline`",
        "```markdown\n[literal](docs/a.md)\n[ref]: docs/a.md\n```\n",
        "~~~~\n[literal](docs/a.md)\n~~~\n[still literal](docs/b.md)\n~~~~\n",
        "    [literal](docs/a.md)\n\t[also literal](docs/b.md)\n",
        r"\[literal](docs/a.md)",
        r"[literal\](docs/a.md)",
        "not a link ](docs/a.md)",
        "[incomplete](docs/a.md",
        "[invalid](\n\n docs/a.md)",
        '[web](<https://example.org/a(b)> "Title")',
        "[web](https://example.org/x) [mail](mailto:info@example.org)",
        "[root](/docs/a.md) [remote](//example.org/a) [fragment](#a) [query](?view=raw)",
        "[web]: https://example.org/x\n[anchor]: #scope\n[web] [anchor]",
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


@pytest.mark.parametrize("paragraph", ["\n\n", "\n \t\n", "\r\n\r\n"])
def test_reference_destination_cannot_cross_a_blank_paragraph(paragraph: str) -> None:
    """#1268: retain both newlines when checking a reference destination boundary."""
    source = f"[ref]:{paragraph}  docs/a.md\n"
    assert links_relative_to_github(source) == source


def load_sync_script() -> ModuleType:
    path = Path(__file__).resolve().parents[1] / "scripts" / "sync_agent_docs.py"
    spec = importlib.util.spec_from_file_location("sync_agent_docs_test", path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_sync_entrypoint_rebases_links_and_preserves_source(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source = tmp_path / "AGENTS.md"
    output = tmp_path / "copilot-instructions.md"
    raw = "[guide](<methods/guide.md>)\n`[example](docs/example.md)`\n[ref]: docs/ref.md\n"
    source.write_text(raw, encoding="utf-8")
    module = load_sync_script()
    monkeypatch.setattr(module, "AGENTS_MD", source)
    monkeypatch.setattr(module, "COPILOT_MD", output)
    module.main()
    assert output.read_text(encoding="utf-8") == (
        module.COPILOT_HEADER
        + "[guide](<../methods/guide.md>)\n`[example](docs/example.md)`\n[ref]: ../docs/ref.md\n"
        + module.PROJECT_CONVENTIONS
    )
    assert source.read_text(encoding="utf-8") == raw
