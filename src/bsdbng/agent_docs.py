"""Helpers for generating tool-specific repository instructions."""

from urllib.parse import urlsplit


def links_relative_to_github(content: str) -> str:
    """Rebase repository-relative Markdown links from the root into `.github`."""
    chunks = content.split("](")
    rewritten = [chunks[0]]
    for chunk in chunks[1:]:
        target, separator, remainder = chunk.partition(")")
        if (
            separator
            and target
            and not urlsplit(target).scheme
            and not target.startswith(("#", "/", "../"))
        ):
            target = "../" + target.removeprefix("./")
        rewritten.append("](" + target + separator + remainder)
    return "".join(rewritten)
