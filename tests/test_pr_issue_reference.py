"""Tests for PR body issue-reference enforcement (#12)."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from types import ModuleType

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "scripts" / "check_pr_issue_reference.py"
WORKFLOW_PATH = REPO_ROOT / ".github" / "workflows" / "pr-issue-ref.yaml"


def _load_validator_module() -> ModuleType:
    assert SCRIPT_PATH.is_file(), f"Missing validator script: {SCRIPT_PATH}"
    spec = importlib.util.spec_from_file_location("check_pr_issue_reference", SCRIPT_PATH)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize(
    ("body", "expected_closing", "expected_referenced"),
    [
        ("Closes #12", (12,), ()),
        ("closes: #123", (123,), ()),
        ("Closed #45", (45,), ()),
        ("Fix #7", (7,), ()),
        ("Fixes #89", (89,), ()),
        ("Fixed #90", (90,), ()),
        ("Resolve #11", (11,), ()),
        ("Resolves #1487", (1487,), ()),
        ("Resolved #2048", (2048,), ()),
        ("Refs #12", (), (12,)),
        ("Ref: #1487", (), (1487,)),
        ("Reference #568", (), (568,)),
        ("References #744", (), (744,)),
        ("Closes #12, #15", (12, 15), ()),
        ("Refs #1146, #1148", (), (1146, 1148)),
        ("Summary text.\n\nCloses #12\nRefs #1487\n", (12,), (1487,)),
        (
            "Fixes CultureBotAI/microbial-trait-profiling#42",
            (42,),
            (),
        ),
        (
            "Refs https://github.com/CultureBotAI/microbial-trait-profiling/issues/965",
            (),
            (965,),
        ),
        (
            "Closes https://github.com/culturebotai/microbial-trait-profiling/issues/1209",
            (1209,),
            (),
        ),
    ],
)
def test_valid_issue_references_pass(
    body: str,
    expected_closing: tuple[int, ...],
    expected_referenced: tuple[int, ...],
) -> None:
    mod = _load_validator_module()
    result = mod.validate_pr_issue_reference(body)
    assert result.ok is True
    assert result.exempt is False
    assert result.closing_issues == expected_closing
    assert result.referenced_issues == expected_referenced
    assert result.negated_issues == ()
    assert "Valid PR issue reference" in result.message


@pytest.mark.parametrize(
    "body",
    [
        None,
        "",
        "   \n\t  ",
        "## Summary\n\nJust a summary with no issue link.",
        "Mentions bare #12 without a linking keyword.",
        "Closes #0",
        "Refs #0",
        "Closes #012",
        "Refs #12abc",
        "Closes #12-slug",
        "<!-- Closes #12 -->",
        '## Issue\n\n<!-- Link the issue this PR addresses: "Closes #N" or "Refs #N" -->\n',
        "```markdown\nCloses #12\n```",
        "~~~\nRefs #12\n~~~",
        "Use `Closes #12` or `Refs #12` in the PR body.",
        "Closes other-org/other-repo#12",
        "Refs https://github.com/other-org/other-repo/issues/12",
        "Closes https://github.com/CultureBotAI/microbial-trait-profiling/pull/12",
    ],
)
def test_missing_commented_code_or_malformed_references_fail(body: str | None) -> None:
    mod = _load_validator_module()
    result = mod.validate_pr_issue_reference(body)
    assert result.ok is False
    assert result.exempt is False
    assert result.closing_issues == ()
    assert result.referenced_issues == ()
    assert "Closes #" in result.message
    assert "Refs #" in result.message


@pytest.mark.parametrize(
    ("body", "expected_negated"),
    [
        ("This PR does not close #12.", (12,)),
        ("Do not close #12 yet.", (12,)),
        ("Did not fix #45.", (45,)),
        ("Will not resolve #88 in this PR.", (88,)),
        ("This doesn't close #12.", (12,)),
        ("Don't fix #12 here.", (12,)),
        ("Won't resolve #99.", (99,)),
        ("Cannot close #12 until part 2 lands.", (12,)),
        ("Should not close #12.", (12,)),
        ("Not closing #12; not fixes #12; never closes #12.", (12,)),
        ("Without closing #12 or resolving #13.", (12, 13)),
        ("Does not ref #12.", (12,)),
    ],
)
def test_negated_closing_phrases_are_rejected(
    body: str,
    expected_negated: tuple[int, ...],
) -> None:
    mod = _load_validator_module()
    result = mod.validate_pr_issue_reference(body)
    assert result.ok is False
    assert result.closing_issues == ()
    assert result.referenced_issues == ()
    assert result.negated_issues == expected_negated
    assert "negated" in result.message.lower()


def test_negated_closing_phrase_plus_explicit_refs_passes() -> None:
    """A PR body that explains 'does not close #1487' and writes 'Refs #1487' is valid."""
    mod = _load_validator_module()
    body = (
        "This PR covers only the diagnostic export markers and does not close #1487.\n\n"
        "Refs #1487\n"
    )
    result = mod.validate_pr_issue_reference(body)
    assert result.ok is True
    assert result.closing_issues == ()
    assert result.referenced_issues == (1487,)
    assert result.negated_issues == (1487,)


@pytest.mark.parametrize(
    ("actor", "head_ref", "expected_exempt"),
    [
        ("dependabot[bot]", "dependabot/pip/ruff-0.8.0", True),
        ("dependabot", "dependabot/github_actions/actions-checkout-7", True),
        ("renovate[bot]", "renovate/configure", True),
        ("github-actions[bot]", "automated-update", True),
        ("realmarcin", "12-pr-body-issue-reference", False),
        ("turbomam", "dependabot/pip/fake-human-branch", False),
        (None, "12-pr-body-issue-reference", False),
    ],
)
def test_bot_actor_exemption_contract(
    actor: str | None,
    head_ref: str | None,
    expected_exempt: bool,
) -> None:
    mod = _load_validator_module()
    assert mod.is_exempt_actor(actor, head_ref) is expected_exempt
    result = mod.validate_pr_issue_reference("", actor=actor, head_ref=head_ref)
    assert result.ok is expected_exempt
    assert result.exempt is expected_exempt


def test_cli_with_event_path_and_flags(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    mod = _load_validator_module()
    valid_event = tmp_path / "valid_event.json"
    valid_event.write_text(
        json.dumps(
            {
                "pull_request": {
                    "number": 2140,
                    "body": "## Summary\nAdd schema version.\n\nRefs #1487\n",
                    "user": {"login": "realmarcin"},
                    "head": {"ref": "1487-diagnostic-export-schema-versions"},
                },
                "repository": {"full_name": "CultureBotAI/microbial-trait-profiling"},
            }
        ),
        encoding="utf-8",
    )
    rc = mod.main(["--event-path", str(valid_event)])
    captured = capsys.readouterr()
    assert rc == 0
    assert "#1487" in captured.out
    assert captured.err == ""

    null_body_event = tmp_path / "null_body_event.json"
    null_body_event.write_text(
        json.dumps(
            {
                "pull_request": {
                    "number": 2145,
                    "body": None,
                    "user": {"login": "realmarcin"},
                    "head": {"ref": "12-pr-body-issue-reference"},
                },
                "repository": {"full_name": "CultureBotAI/microbial-trait-profiling"},
            }
        ),
        encoding="utf-8",
    )
    rc_fail = mod.main(["--event-path", str(null_body_event)])
    captured_fail = capsys.readouterr()
    assert rc_fail == 1
    assert "Closes #" in captured_fail.err
    assert "Refs #" in captured_fail.err

    bot_event = tmp_path / "bot_event.json"
    bot_event.write_text(
        json.dumps(
            {
                "pull_request": {
                    "number": 2146,
                    "body": "Bumps actions/checkout from 6 to 7.",
                    "user": {"login": "dependabot[bot]"},
                    "head": {"ref": "dependabot/github_actions/actions-checkout-7"},
                },
                "repository": {"full_name": "CultureBotAI/microbial-trait-profiling"},
            }
        ),
        encoding="utf-8",
    )
    rc_bot = mod.main(["--event-path", str(bot_event)])
    captured_bot = capsys.readouterr()
    assert rc_bot == 0
    assert "Exempt" in captured_bot.out

    assert mod.main(["--body", "Closes #12"]) == 0
    capsys.readouterr()
    assert mod.main(["--body", "No issue reference here"]) == 1
    capsys.readouterr()


def test_cli_rejects_malformed_event_payload(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    mod = _load_validator_module()
    bad_event = tmp_path / "bad_event.json"
    bad_event.write_text(json.dumps({"action": "opened"}), encoding="utf-8")
    rc = mod.main(["--event-path", str(bad_event)])
    captured = capsys.readouterr()
    assert rc == 1
    assert "pull_request" in captured.err


def test_workflow_file_contract() -> None:
    assert WORKFLOW_PATH.is_file(), f"Missing workflow file: {WORKFLOW_PATH}"
    raw = WORKFLOW_PATH.read_text(encoding="utf-8")
    parsed = yaml.safe_load(raw)
    assert isinstance(parsed, dict)
    assert parsed.get("name") == "PR Issue Reference"
    on_block = parsed.get("on") or parsed.get(True)
    assert isinstance(on_block, dict)
    assert "pull_request" in on_block
    pr_trigger = on_block["pull_request"]
    assert isinstance(pr_trigger, dict)
    assert pr_trigger.get("types") == ["opened", "edited", "synchronize", "reopened"]
    assert parsed.get("permissions") == {"contents": "read"}
    jobs = parsed.get("jobs")
    assert isinstance(jobs, dict)
    assert "issue-reference" in jobs
    job = jobs["issue-reference"]
    assert job.get("runs-on") == "ubuntu-latest"
    steps = job.get("steps")
    assert isinstance(steps, list)
    assert len(steps) >= 2
    assert any(
        isinstance(step, dict)
        and "scripts/check_pr_issue_reference.py" in str(step.get("run", ""))
        and "GITHUB_EVENT_PATH" in str(step.get("run", ""))
        for step in steps
    )
