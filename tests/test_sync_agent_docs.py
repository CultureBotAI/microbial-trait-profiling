"""Generation of nested tool instructions from the root AGENTS.md."""

from bsdbng.agent_docs import links_relative_to_github


def test_copilot_links_are_relative_to_github_directory() -> None:
    source = (
        "[guide](methods/guide.md) [local](./docs/local.md) "
        "[parent](../shared.md) [web](https://example.org/x) [heading](#scope)"
    )

    assert links_relative_to_github(source) == (
        "[guide](../methods/guide.md) [local](../docs/local.md) "
        "[parent](../shared.md) [web](https://example.org/x) [heading](#scope)"
    )
