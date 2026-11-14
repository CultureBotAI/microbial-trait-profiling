"""Tests for PR head branch-name enforcement (#11)."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from types import ModuleType

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "scripts" / "check_pr_branch_name.py"
WORKFLOW_PATH = REPO_ROOT / ".github" / "workflows" / "pr-branch-name.yaml"


def _load_validator_module() -> ModuleType:
    assert SCRIPT_PATH.is_file(), f"Missing validator script: {SCRIPT_PATH}"
    spec = importlib.util.spec_from_file_location("check_pr_branch_name", SCRIPT_PATH)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize(
    ("head_ref", "expected_issue", "expected_slug"),
    [
        ("11-branch-name-issue-prefix", 11, "branch-name-issue-prefix"),
        ("12-pr-body-issue-reference", 12, "pr-body-issue-reference"),
        ("30-aggregate-tables", 30, "aggregate-tables"),
        ("1487-diagnostic-export-schema-versions", 1487, "diagnostic-export-schema-versions"),
        ("issue-12-bootstrap-ingest", 12, "bootstrap-ingest"),
        ("issue-27-tighten-schema", 27, "tighten-schema"),
        ("568-three_arm_ora.v2", 568, "three_arm_ora.v2"),
        ("refs/heads/744-validation-feasibility-docs", 744, "validation-feasibility-docs"),
    ],
)
def test_valid_issue_branch_names_pass(
    head_ref: str,
    expected_issue: int,
    expected_slug: str,
) -> None:
    mod = _load_validator_module()
    result = mod.validate_pr_branch_name(head_ref)
    assert result.ok is True
    assert result.exempt is False
    assert result.issue_number == expected_issue
    assert result.slug == expected_slug
    assert f"#{expected_issue}" in result.message


@pytest.mark.parametrize(
    "head_ref",
    [
        None,
        "",
        "   ",
        "main",
        "master",
        "develop",
        "patch-1",
        "feature/add-new-table",
        "fix-calibration",
        "0-zero-issue",
        "012-leading-zero",
        "issue-0-zero-issue",
        "issue-012-leading-zero",
        "12",
        "12-",
        "issue-12",
        "issue-12-",
        "12--double-dash",
        "12-trailing-dash-",
        "12-UpperCaseSlug",
        "12-slug/with/slash",
        "12-slug with spaces",
    ],
)
def test_invalid_branch_names_fail(head_ref: str | None) -> None:
    mod = _load_validator_module()
    result = mod.validate_pr_branch_name(head_ref)
    assert result.ok is False
    assert result.exempt is False
    assert result.issue_number is None
    assert result.slug is None
    assert "<issue-number>-<slug>" in result.message


@pytest.mark.parametrize(
    ("head_ref", "actor", "expected_exempt"),
    [
        ("dependabot/pip/ruff-0.8.0", "dependabot[bot]", True),
        ("dependabot/github_actions/actions-checkout-7", "realmarcin", True),
        ("renovate/configure", "renovate[bot]", True),
        ("copilot/fix-12", "copilot-swe-agent[bot]", True),
        ("automated-sync", "github-actions[bot]", True),
        ("feature-branch", "realmarcin", False),
        ("patch-1", "turbomam", False),
    ],
)
def test_bot_and_prefix_exemptions(
    head_ref: str,
    actor: str | None,
    expected_exempt: bool,
) -> None:
    mod = _load_validator_module()
    assert mod.is_exempt_branch(head_ref, actor=actor) is expected_exempt
    result = mod.validate_pr_branch_name(head_ref, actor=actor)
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
                    "user": {"login": "realmarcin"},
                    "head": {"ref": "1487-diagnostic-export-schema-versions"},
                }
            }
        ),
        encoding="utf-8",
    )
    rc = mod.main(["--event-path", str(valid_event)])
    captured = capsys.readouterr()
    assert rc == 0
    assert "#1487" in captured.out
    assert captured.err == ""

    bad_branch_event = tmp_path / "bad_branch_event.json"
    bad_branch_event.write_text(
        json.dumps(
            {
                "pull_request": {
                    "number": 2145,
                    "user": {"login": "realmarcin"},
                    "head": {"ref": "feature-no-issue-number"},
                }
            }
        ),
        encoding="utf-8",
    )
    rc_fail = mod.main(["--event-path", str(bad_branch_event)])
    captured_fail = capsys.readouterr()
    assert rc_fail == 1
    assert "<issue-number>-<slug>" in captured_fail.err

    bot_event = tmp_path / "bot_event.json"
    bot_event.write_text(
        json.dumps(
            {
                "pull_request": {
                    "number": 2146,
                    "user": {"login": "dependabot[bot]"},
                    "head": {"ref": "dependabot/github_actions/actions-checkout-7"},
                }
            }
        ),
        encoding="utf-8",
    )
    rc_bot = mod.main(["--event-path", str(bot_event)])
    captured_bot = capsys.readouterr()
    assert rc_bot == 0
    assert "Exempt" in captured_bot.out

    assert mod.main(["--head-ref", "11-branch-name-issue-prefix"]) == 0
    capsys.readouterr()
    assert mod.main(["--head-ref", "feature-branch"]) == 1
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
    assert parsed.get("name") == "PR Branch Name"
    on_block = parsed.get("on") or parsed.get(True)
    assert isinstance(on_block, dict)
    assert "pull_request" in on_block
    pr_trigger = on_block["pull_request"]
    assert isinstance(pr_trigger, dict)
    assert pr_trigger.get("types") == ["opened", "edited", "synchronize", "reopened"]
    assert parsed.get("permissions") == {"contents": "read"}
    jobs = parsed.get("jobs")
    assert isinstance(jobs, dict)
    assert "branch-name" in jobs
    job = jobs["branch-name"]
    assert job.get("runs-on") == "ubuntu-latest"
    steps = job.get("steps")
    assert isinstance(steps, list)
    assert len(steps) >= 2
    assert any(
        isinstance(step, dict)
        and "scripts/check_pr_branch_name.py" in str(step.get("run", ""))
        and "GITHUB_EVENT_PATH" in str(step.get("run", ""))
        for step in steps
    )
