#!/usr/bin/env python3
"""Validate that a pull request head branch uses an issue-linked name (#11).

Enforces the project's issue -> branch -> PR convention on pull request head
branches:
- `<issue-number>-<slug>` (standing convention in `CONTRIBUTING.md`, e.g.
  `11-branch-name-issue-prefix`, `30-aggregate-tables`)
- `issue-<issue-number>-<slug>` (proposed convention in #11, e.g.
  `issue-12-bootstrap-ingest`, `issue-27-tighten-schema`)

where `<issue-number>` is a positive integer (`[1-9][0-9]*`) and `<slug>` is a
non-empty lowercase alphanumeric slug (`[a-z0-9]+(?:[._-][a-z0-9]+)*`).

Automated bot accounts (`dependabot[bot]`, `renovate[bot]`, `github-actions[bot]`,
`copilot-swe-agent[bot]`) and automated tool branch prefixes (`dependabot/`,
`renovate/`, `copilot/`) are explicitly exempted.
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

ISSUE_NUMBER_PATTERN = r"[1-9][0-9]*"
SLUG_PATTERN = r"[a-z0-9]+(?:[._-][a-z0-9]+)*"
BRANCH_PATTERN = re.compile(
    rf"^(?:issue-)?(?P<issue_number>{ISSUE_NUMBER_PATTERN})-(?P<slug>{SLUG_PATTERN})$"
)

EXEMPT_BOT_ACTORS: frozenset[str] = frozenset(
    {
        "dependabot[bot]",
        "dependabot",
        "renovate[bot]",
        "renovate",
        "github-actions[bot]",
        "copilot-swe-agent[bot]",
        "copilot",
    }
)
EXEMPT_BRANCH_PREFIXES: tuple[str, ...] = (
    "dependabot/",
    "renovate/",
    "copilot/",
)


@dataclass(frozen=True)
class BranchValidationResult:
    """Structured outcome of validating a pull request head branch name."""

    ok: bool
    exempt: bool
    branch: str
    issue_number: int | None
    slug: str | None
    message: str


def normalize_branch_ref(head_ref: str | None) -> str:
    """Strip whitespace and an optional leading `refs/heads/` prefix."""
    raw = (head_ref or "").strip()
    if raw.startswith("refs/heads/"):
        return raw[len("refs/heads/") :]
    return raw


def is_exempt_branch(head_ref: str | None, *, actor: str | None = None) -> bool:
    """Return True when the actor or head branch prefix is an exempt automation bot."""
    normalized = normalize_branch_ref(head_ref)
    if any(normalized.startswith(prefix) for prefix in EXEMPT_BRANCH_PREFIXES):
        return True
    return bool(actor and actor.strip().lower() in EXEMPT_BOT_ACTORS)


def validate_pr_branch_name(
    head_ref: str | None,
    *,
    actor: str | None = None,
) -> BranchValidationResult:
    """Validate that `head_ref` follows the issue-prefixed branch naming contract."""
    normalized = normalize_branch_ref(head_ref)
    if is_exempt_branch(normalized, actor=actor):
        return BranchValidationResult(
            ok=True,
            exempt=True,
            branch=normalized,
            issue_number=None,
            slug=None,
            message=f"Exempt automation branch or actor ({normalized or actor}).",
        )

    if not normalized:
        return BranchValidationResult(
            ok=False,
            exempt=False,
            branch="",
            issue_number=None,
            slug=None,
            message=(
                "PR head branch name is empty. Use '<issue-number>-<slug>' "
                "(e.g. '30-aggregate-tables') or 'issue-<issue-number>-<slug>' "
                "(e.g. 'issue-12-bootstrap-ingest')."
            ),
        )

    match = BRANCH_PATTERN.match(normalized)
    if match is not None:
        issue_number = int(match.group("issue_number"))
        slug = match.group("slug")
        return BranchValidationResult(
            ok=True,
            exempt=False,
            branch=normalized,
            issue_number=issue_number,
            slug=slug,
            message=(
                f"Valid issue-linked PR branch '{normalized}' "
                f"(issue #{issue_number}, slug '{slug}')."
            ),
        )

    return BranchValidationResult(
        ok=False,
        exempt=False,
        branch=normalized,
        issue_number=None,
        slug=None,
        message=(
            f"Invalid PR head branch name '{normalized}'. Branch names must match "
            "'<issue-number>-<slug>' (e.g. '30-aggregate-tables') or "
            "'issue-<issue-number>-<slug>' (e.g. 'issue-12-bootstrap-ingest') "
            "with a positive integer issue number and lowercase alphanumeric slug."
        ),
    )


def load_pull_request_event(event_path: Path) -> tuple[str | None, str | None]:
    """Extract `(head_ref, actor)` from a GitHub Actions `pull_request` event JSON."""
    payload = json.loads(event_path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or not isinstance(payload.get("pull_request"), dict):
        raise ValueError(
            f"GitHub event payload at {event_path} does not contain a 'pull_request' object."
        )
    pr = payload["pull_request"]
    head = pr.get("head")
    raw_ref = head.get("ref") if isinstance(head, dict) else None
    head_ref: str | None = raw_ref if isinstance(raw_ref, str) else None

    user = pr.get("user")
    sender = payload.get("sender")
    actor: str | None = None
    if isinstance(user, dict) and isinstance(user.get("login"), str):
        actor = user["login"]
    elif isinstance(sender, dict) and isinstance(sender.get("login"), str):
        actor = sender["login"]

    return head_ref, actor


def build_parser() -> argparse.ArgumentParser:
    """Build the CLI argument parser."""
    parser = argparse.ArgumentParser(
        description="Check that a PR head branch follows issue-based naming (#11)."
    )
    parser.add_argument(
        "--event-path",
        type=Path,
        default=None,
        help="Path to the GitHub Actions pull_request event JSON (defaults to $GITHUB_EVENT_PATH).",
    )
    parser.add_argument(
        "--head-ref",
        type=str,
        default=None,
        help="Explicit pull request head branch name (overrides --event-path when provided).",
    )
    parser.add_argument(
        "--actor",
        type=str,
        default=None,
        help="Pull request author login (optional override).",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """CLI entry point returning 0 on valid branch name and 1 on failure."""
    parser = build_parser()
    args = parser.parse_args(argv)

    head_ref: str | None = args.head_ref
    actor: str | None = args.actor

    event_path = args.event_path
    if event_path is None and head_ref is None:
        env_event = os.environ.get("GITHUB_EVENT_PATH")
        if env_event:
            event_path = Path(env_event)

    if event_path is not None and head_ref is None:
        try:
            ev_head_ref, ev_actor = load_pull_request_event(event_path)
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            print(f"Error reading pull_request event payload: {exc}", file=sys.stderr)
            return 1
        head_ref = ev_head_ref
        if actor is None:
            actor = ev_actor

    result = validate_pr_branch_name(head_ref, actor=actor)
    if result.ok:
        print(result.message)
        return 0
    print(result.message, file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
