#!/usr/bin/env python3
"""Validate that a pull request body references a repository issue (#12).

Enforces the project's issue -> branch -> PR contract by requiring every
non-exempt pull request body to include at least one explicit issue reference:
- Closing keywords (`Closes #N`, `Fixes #N`, `Resolves #N`, etc.) when the PR
  fully resolves the issue.
- Reference keywords (`Refs #N`, `Ref #N`, `References #N`) when the PR partially
  addresses or relates to the issue (#1317).

HTML comments (`<!-- ... -->`), fenced code blocks, and inline code spans are
stripped before matching so template placeholders or code examples cannot satisfy
the check. Negated phrases (such as `does not close #12` or `without closing #12`)
are detected and rejected unless accompanied by a non-negated reference such as
`Refs #12`.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

DEFAULT_REPOSITORY = "CultureBotAI/microbial-trait-profiling"
ISSUE_NUMBER_PATTERN = r"[1-9][0-9]*"

CLOSING_KEYWORDS: tuple[str, ...] = (
    "close",
    "closes",
    "closed",
    "fix",
    "fixes",
    "fixed",
    "resolve",
    "resolves",
    "resolved",
)
REFERENCE_KEYWORDS: tuple[str, ...] = (
    "ref",
    "refs",
    "reference",
    "references",
)
GERUND_KEYWORDS: tuple[str, ...] = (
    "closing",
    "fixing",
    "resolving",
    "referencing",
)

EXEMPT_BOT_ACTORS: frozenset[str] = frozenset(
    {
        "dependabot[bot]",
        "dependabot",
        "renovate[bot]",
        "renovate",
        "github-actions[bot]",
    }
)

_HTML_COMMENT_RE = re.compile(r"<!--.*?-->", re.DOTALL)
_FENCED_BACKTICK_RE = re.compile(r"```.*?```", re.DOTALL)
_FENCED_TILDE_RE = re.compile(r"~~~.*?~~~", re.DOTALL)
_INLINE_CODE_RE = re.compile(r"`[^`\n]+`")

_NEGATION_PREFIX_RE = re.compile(
    r"(?:\b(?:does|do|did|will|would|should|can|could|must)\s+not"
    r"|\b(?:doesn['\u2019]t|don['\u2019]t|didn['\u2019]t|won['\u2019]t"
    r"|wouldn['\u2019]t|shouldn['\u2019]t|can['\u2019]t|cannot|mustn['\u2019]t)"
    r"|\b(?:not|never|without|neither|nor))\s+$",
    re.IGNORECASE,
)

_CLAUSE_SPLIT_RE = re.compile(r"(?:[;\n\r]+|\.(?=\s|$))")


@dataclass(frozen=True)
class ValidationResult:
    """Structured outcome of validating a pull request body."""

    ok: bool
    exempt: bool
    closing_issues: tuple[int, ...]
    referenced_issues: tuple[int, ...]
    negated_issues: tuple[int, ...]
    message: str


def strip_non_prose_regions(body: str) -> str:
    """Remove HTML comments, fenced code blocks, and inline code spans."""
    cleaned = _HTML_COMMENT_RE.sub(" ", body)
    cleaned = _FENCED_BACKTICK_RE.sub(" ", cleaned)
    cleaned = _FENCED_TILDE_RE.sub(" ", cleaned)
    cleaned = _INLINE_CODE_RE.sub(" ", cleaned)
    return cleaned


def is_exempt_actor(actor: str | None, head_ref: str | None = None) -> bool:
    """Return True when the PR author is an exempt automation bot account."""
    del head_ref
    if not actor:
        return False
    return actor.strip().lower() in EXEMPT_BOT_ACTORS


def _build_target_pattern(repository: str, *, with_groups: bool = False) -> str:
    escaped_repo = re.escape(repository.strip())
    if with_groups:
        return (
            r"(?:"
            rf"#(?P<num_short>{ISSUE_NUMBER_PATTERN})"
            rf"|{escaped_repo}#(?P<num_repo>{ISSUE_NUMBER_PATTERN})"
            rf"|https?://github\.com/{escaped_repo}/issues/(?P<num_url>{ISSUE_NUMBER_PATTERN})"
            r")(?![0-9A-Za-z_-])"
        )
    return (
        r"(?:"
        rf"#{ISSUE_NUMBER_PATTERN}"
        rf"|{escaped_repo}#{ISSUE_NUMBER_PATTERN}"
        rf"|https?://github\.com/{escaped_repo}/issues/{ISSUE_NUMBER_PATTERN}"
        r")(?![0-9A-Za-z_-])"
    )


def _extract_target_numbers(
    target_blob: str,
    target_re: re.Pattern[str],
) -> list[int]:
    numbers: list[int] = []
    for match in target_re.finditer(target_blob):
        raw_num = match.group("num_short") or match.group("num_repo") or match.group("num_url")
        if raw_num is not None:
            value = int(raw_num)
            if value not in numbers:
                numbers.append(value)
    return numbers


def _dedup_preserve_order(items: Sequence[int]) -> tuple[int, ...]:
    seen: set[int] = set()
    ordered: list[int] = []
    for item in items:
        if item not in seen:
            seen.add(item)
            ordered.append(item)
    return tuple(ordered)


def validate_pr_issue_reference(
    body: str | None,
    *,
    actor: str | None = None,
    head_ref: str | None = None,
    repository: str = DEFAULT_REPOSITORY,
) -> ValidationResult:
    """Validate that `body` contains at least one non-negated issue reference."""
    if is_exempt_actor(actor, head_ref):
        return ValidationResult(
            ok=True,
            exempt=True,
            closing_issues=(),
            referenced_issues=(),
            negated_issues=(),
            message=f"Exempt automation actor ({actor}).",
        )

    cleaned = strip_non_prose_regions(body or "")
    target_pattern = _build_target_pattern(repository, with_groups=False)
    target_re = re.compile(_build_target_pattern(repository, with_groups=True), re.IGNORECASE)

    all_keywords = (*CLOSING_KEYWORDS, *REFERENCE_KEYWORDS, *GERUND_KEYWORDS)
    keyword_group = "|".join(re.escape(word) for word in all_keywords)
    list_continuation = rf"(?:\s*(?:,|\band\b|\bor\b)\s*{target_pattern})*"
    statement_re = re.compile(
        rf"\b(?P<keyword>{keyword_group})\b\s*:?\s*(?P<targets>{target_pattern}{list_continuation})",
        re.IGNORECASE,
    )

    closing_found: list[int] = []
    referenced_found: list[int] = []
    negated_found: list[int] = []

    for clause in _CLAUSE_SPLIT_RE.split(cleaned):
        clause_negated = False
        last_end = 0
        for match in statement_re.finditer(clause):
            prefix = clause[last_end : match.start()]
            between = clause[last_end : match.start()].strip().lower()
            is_directly_negated = _NEGATION_PREFIX_RE.search(prefix) is not None
            is_chained_negated = clause_negated and between in {"or", "and", ",", "nor"}
            keyword = match.group("keyword").lower()
            numbers = _extract_target_numbers(match.group("targets"), target_re)

            if is_directly_negated or is_chained_negated or keyword in GERUND_KEYWORDS:
                if is_directly_negated or is_chained_negated:
                    clause_negated = True
                    negated_found.extend(numbers)
                last_end = match.end()
                continue

            clause_negated = False
            if keyword in CLOSING_KEYWORDS:
                closing_found.extend(numbers)
            elif keyword in REFERENCE_KEYWORDS:
                referenced_found.extend(numbers)
            last_end = match.end()

    closing_tuple = _dedup_preserve_order(closing_found)
    referenced_tuple = _dedup_preserve_order(referenced_found)
    negated_tuple = _dedup_preserve_order(negated_found)

    if closing_tuple or referenced_tuple:
        parts: list[str] = []
        if closing_tuple:
            parts.append("closing=" + ",".join(f"#{n}" for n in closing_tuple))
        if referenced_tuple:
            parts.append("refs=" + ",".join(f"#{n}" for n in referenced_tuple))
        return ValidationResult(
            ok=True,
            exempt=False,
            closing_issues=closing_tuple,
            referenced_issues=referenced_tuple,
            negated_issues=negated_tuple,
            message=f"Valid PR issue reference ({'; '.join(parts)}).",
        )

    if negated_tuple:
        negated_str = ", ".join(f"#{n}" for n in negated_tuple)
        return ValidationResult(
            ok=False,
            exempt=False,
            closing_issues=(),
            referenced_issues=(),
            negated_issues=negated_tuple,
            message=(
                f"PR body contains only negated issue reference(s) ({negated_str}). "
                "Include an affirmative reference such as 'Closes #123' (when fully "
                "resolving an issue) or 'Refs #123' (for partial or related work)."
            ),
        )

    return ValidationResult(
        ok=False,
        exempt=False,
        closing_issues=(),
        referenced_issues=(),
        negated_issues=(),
        message=(
            "PR body must include an explicit issue reference such as "
            "'Closes #123', 'Fixes #123', 'Resolves #123', or 'Refs #123'."
        ),
    )


def load_pull_request_event(
    event_path: Path,
) -> tuple[str | None, str | None, str | None, str]:
    """Extract `(body, actor, head_ref, repository)` from a GitHub event JSON file."""
    payload = json.loads(event_path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or not isinstance(payload.get("pull_request"), dict):
        raise ValueError(
            f"GitHub event payload at {event_path} does not contain a 'pull_request' object."
        )
    pr = payload["pull_request"]
    raw_body = pr.get("body")
    body = str(raw_body) if isinstance(raw_body, str) else None

    user = pr.get("user")
    sender = payload.get("sender")
    actor: str | None = None
    if isinstance(user, dict) and isinstance(user.get("login"), str):
        actor = user["login"]
    elif isinstance(sender, dict) and isinstance(sender.get("login"), str):
        actor = sender["login"]

    head = pr.get("head")
    raw_ref = head.get("ref") if isinstance(head, dict) else None
    head_ref: str | None = raw_ref if isinstance(raw_ref, str) else None

    repo_obj = payload.get("repository")
    raw_repo = repo_obj.get("full_name") if isinstance(repo_obj, dict) else None
    repository: str = raw_repo if isinstance(raw_repo, str) else DEFAULT_REPOSITORY
    return body, actor, head_ref, repository


def build_parser() -> argparse.ArgumentParser:
    """Build the CLI argument parser."""
    parser = argparse.ArgumentParser(
        description="Check that a PR body contains a required issue reference (#12)."
    )
    parser.add_argument(
        "--event-path",
        type=Path,
        default=None,
        help="Path to the GitHub Actions pull_request event JSON (defaults to $GITHUB_EVENT_PATH).",
    )
    parser.add_argument(
        "--body",
        type=str,
        default=None,
        help="Explicit pull request body text (overrides --event-path when provided).",
    )
    parser.add_argument(
        "--actor",
        type=str,
        default=None,
        help="Pull request author login (optional override).",
    )
    parser.add_argument(
        "--head-ref",
        type=str,
        default=None,
        help="Pull request head branch name (optional override).",
    )
    parser.add_argument(
        "--repository",
        type=str,
        default=DEFAULT_REPOSITORY,
        help=f"Repository full name (default: {DEFAULT_REPOSITORY}).",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """CLI entry point returning 0 on valid PR body and 1 on failure."""
    parser = build_parser()
    args = parser.parse_args(argv)

    body: str | None = args.body
    actor: str | None = args.actor
    head_ref: str | None = args.head_ref
    repository: str = args.repository

    event_path = args.event_path
    if event_path is None and body is None:
        env_event = os.environ.get("GITHUB_EVENT_PATH")
        if env_event:
            event_path = Path(env_event)

    if event_path is not None and body is None:
        try:
            ev_body, ev_actor, ev_head_ref, ev_repo = load_pull_request_event(event_path)
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            print(f"Error reading pull_request event payload: {exc}", file=sys.stderr)
            return 1
        body = ev_body
        if actor is None:
            actor = ev_actor
        if head_ref is None:
            head_ref = ev_head_ref
        if args.repository == DEFAULT_REPOSITORY:
            repository = ev_repo

    result = validate_pr_issue_reference(
        body,
        actor=actor,
        head_ref=head_ref,
        repository=repository,
    )
    if result.ok:
        print(result.message)
        return 0
    print(result.message, file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
