"""Pin CR references, bare control rejection and definition bounds (#1380/#1381/#1475)."""

import pytest

from mtp.agent_docs import _destination, links_relative_to_github


@pytest.mark.parametrize("newline", ["\r", "\n", "\r\n"])
def test_reference_line_endings(newline: str) -> None:
    source = newline.join(["Intro", "", "[r]: docs/a.md", "[r]", ""])
    assert links_relative_to_github(source) == source.replace("docs/a.md", "../docs/a.md")


@pytest.mark.parametrize(
    "control", [chr(n) for n in range(32) if chr(n) not in " \t\r\n"] + ["\x7f"]
)
@pytest.mark.parametrize("prefix", ["", "\\"])
def test_bare_destination_controls(control: str, prefix: str) -> None:
    target = f"docs/a{prefix}{control}.md"
    assert _destination(target, 0) is None
    for source in (f"[x]({target})\n", f"[r]: {target}\n"):
        assert links_relative_to_github(source) == source


@pytest.mark.parametrize("source", ["[x](docs/a\u00a0.md)\n", "[r]: docs/a\u00a0.md\n"])
def test_nbsp_is_allowed(source: str) -> None:
    assert links_relative_to_github(source) == source.replace("docs/", "../docs/")


def test_reference_title_limit() -> None:
    source = '[guide]: docs/guide.md "Guide\n> quoted [x](docs/x.md)"\n\n[guide]\n'
    assert links_relative_to_github(source) == source.replace("docs/x.md", "../docs/x.md")


def test_outside_code_destination_limit() -> None:
    source = "- [note]:\n[ref]: docs/a.md\n\n[ref]\n"
    assert links_relative_to_github(source) == source
