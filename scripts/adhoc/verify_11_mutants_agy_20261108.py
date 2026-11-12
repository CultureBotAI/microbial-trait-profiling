#!/usr/bin/env python3
"""Isolated scratch-copy mutant verifier for Issue #11 (PR head branch-name check).

Purpose:
  Verify that every core branch and contract in `scripts/check_pr_branch_name.py`
  and `.github/workflows/pr-branch-name.yaml` is killed with `AssertionError` by
  `tests/test_pr_branch_name.py` in an isolated scratch copy (never mutating
  the live worktree), and verify import/file origin inside pytest on the unmutated
  control.

Inputs:
  - scripts/check_pr_branch_name.py
  - .github/workflows/pr-branch-name.yaml
  - tests/test_pr_branch_name.py

Outputs:
  Console table of control + 12 mutant runs; exits 0 iff control passes and all
  12 mutants fail with AssertionError (survived=0, other_fails=0).

Issue / PR:
  Issue #11

Date:
  2026-11-08
"""

from __future__ import annotations

import contextlib
import hashlib
import io
import shutil
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import pytest

WORKTREE_ROOT = Path(__file__).resolve().parents[2]


@dataclass(frozen=True)
class Mutant:
    mutant_id: str
    rel_path: str
    description: str
    target: str
    replacement: str


MUTANTS: tuple[Mutant, ...] = (
    Mutant(
        mutant_id="M01_allow_zero_issue_number",
        rel_path="scripts/check_pr_branch_name.py",
        description="Allow #0 / leading zeros in branch issue prefix",
        target='ISSUE_NUMBER_PATTERN = r"[1-9][0-9]*"',
        replacement='ISSUE_NUMBER_PATTERN = r"[0-9]+"',
    ),
    Mutant(
        mutant_id="M02_allow_uppercase_or_double_dash_slug",
        rel_path="scripts/check_pr_branch_name.py",
        description="Relax SLUG_PATTERN to allow uppercase or double dashes",
        target='SLUG_PATTERN = r"[a-z0-9]+(?:[._-][a-z0-9]+)*"',
        replacement='SLUG_PATTERN = r"[A-Za-z0-9._-]+"',
    ),
    Mutant(
        mutant_id="M03_drop_optional_issue_prefix",
        rel_path="scripts/check_pr_branch_name.py",
        description="Reject 'issue-<number>-<slug>' branch convention",
        target=(
            '    rf"^(?:issue-)?(?P<issue_number>{ISSUE_NUMBER_PATTERN})-(?P<slug>{SLUG_PATTERN})$"'
        ),
        replacement=('    rf"^(?P<issue_number>{ISSUE_NUMBER_PATTERN})-(?P<slug>{SLUG_PATTERN})$"'),
    ),
    Mutant(
        mutant_id="M04_require_issue_literal_prefix",
        rel_path="scripts/check_pr_branch_name.py",
        description="Reject '<number>-<slug>' standing convention by requiring 'issue-'",
        target=(
            '    rf"^(?:issue-)?(?P<issue_number>{ISSUE_NUMBER_PATTERN})-(?P<slug>{SLUG_PATTERN})$"'
        ),
        replacement=(
            '    rf"^issue-(?P<issue_number>{ISSUE_NUMBER_PATTERN})-(?P<slug>{SLUG_PATTERN})$"'
        ),
    ),
    Mutant(
        mutant_id="M05_skip_refs_heads_normalization",
        rel_path="scripts/check_pr_branch_name.py",
        description="Do not strip leading refs/heads/ prefix in normalize_branch_ref",
        target='    if raw.startswith("refs/heads/"):\n        return raw[len("refs/heads/") :]',
        replacement="    if False:\n        return raw",
    ),
    Mutant(
        mutant_id="M06_disable_branch_prefix_exemption",
        rel_path="scripts/check_pr_branch_name.py",
        description="Disable dependabot/renovate/copilot branch prefix exemption",
        target=(
            "    if any(normalized.startswith(prefix) for prefix in EXEMPT_BRANCH_PREFIXES):\n"
            "        return True"
        ),
        replacement="    if False:\n        return True",
    ),
    Mutant(
        mutant_id="M07_disable_actor_exemption",
        rel_path="scripts/check_pr_branch_name.py",
        description="Disable EXEMPT_BOT_ACTORS exemption in is_exempt_branch",
        target="    return bool(actor and actor.strip().lower() in EXEMPT_BOT_ACTORS)",
        replacement="    return False",
    ),
    Mutant(
        mutant_id="M08_exempt_all_branches",
        rel_path="scripts/check_pr_branch_name.py",
        description="Return True unconditionally from is_exempt_branch",
        target="    return bool(actor and actor.strip().lower() in EXEMPT_BOT_ACTORS)",
        replacement="    return True",
    ),
    Mutant(
        mutant_id="M09_accept_empty_branch_name",
        rel_path="scripts/check_pr_branch_name.py",
        description="Return ok=True when branch name is empty",
        target=(
            "    if not normalized:\n        return BranchValidationResult(\n            ok=False,"
        ),
        replacement=(
            "    if not normalized:\n        return BranchValidationResult(\n            ok=True,"
        ),
    ),
    Mutant(
        mutant_id="M10_main_returns_zero_on_invalid",
        rel_path="scripts/check_pr_branch_name.py",
        description="Return exit code 0 from main() even when branch validation fails",
        target="    print(result.message, file=sys.stderr)\n    return 1",
        replacement="    print(result.message, file=sys.stderr)\n    return 0",
    ),
    Mutant(
        mutant_id="M11_workflow_omit_synchronize_trigger",
        rel_path=".github/workflows/pr-branch-name.yaml",
        description="Omit 'synchronize' from pull_request.types in workflow YAML",
        target="    types:\n      - opened\n      - edited\n      - synchronize\n      - reopened",
        replacement="    types:\n      - opened\n      - edited\n      - reopened",
    ),
    Mutant(
        mutant_id="M12_workflow_broaden_permissions",
        rel_path=".github/workflows/pr-branch-name.yaml",
        description="Change workflow permissions from contents: read to contents: write",
        target="permissions:\n  contents: read",
        replacement="permissions:\n  contents: write",
    ),
)

REL_FILES: tuple[str, ...] = (
    "scripts/check_pr_branch_name.py",
    ".github/workflows/pr-branch-name.yaml",
    "tests/test_pr_branch_name.py",
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _stage_scratch(scratch_root: Path) -> None:
    for rel in REL_FILES:
        src = WORKTREE_ROOT / rel
        dst = scratch_root / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


def _run_pytest(scratch_root: Path) -> tuple[int, str]:
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        code = int(
            pytest.main(
                [
                    "-q",
                    "--tb=native",
                    str(scratch_root / "tests" / "test_pr_branch_name.py"),
                ]
            )
        )
    return code, buf.getvalue()


def main() -> int:
    worktree_shas = {rel: _sha256(WORKTREE_ROOT / rel) for rel in REL_FILES}

    with tempfile.TemporaryDirectory(prefix="verify_11_mutants_") as tmp_dir:
        scratch_root = Path(tmp_dir).resolve()
        _stage_scratch(scratch_root)

        origin_test = scratch_root / "tests" / "test_pr_branch_name.py"
        original_test_text = origin_test.read_text(encoding="utf-8")
        origin_assertion = (
            "\n\ndef test_scratch_origin_binding() -> None:\n"
            f"    assert str(SCRIPT_PATH).startswith({str(scratch_root)!r})\n"
            f"    assert str(WORKFLOW_PATH).startswith({str(scratch_root)!r})\n"
        )
        origin_test.write_text(original_test_text + origin_assertion, encoding="utf-8")

        control_rc, control_out = _run_pytest(scratch_root)
        if control_rc != 0:
            print("[ERROR] Unmutated control failed:", file=sys.stderr)
            print(control_out, file=sys.stderr)
            return 2
        print(
            f"{'CONTROL':45s} | {'PASSED (rc=0)':24s} | Unmutated scratch copy with origin binding"
        )

        killed = 0
        survived = 0
        other_fails = 0

        for mutant in MUTANTS:
            target_file = scratch_root / mutant.rel_path
            base_text = (WORKTREE_ROOT / mutant.rel_path).read_text(encoding="utf-8")
            occurrences = base_text.count(mutant.target)
            if occurrences != 1:
                print(
                    f"[ERROR] {mutant.mutant_id}: target matched {occurrences} "
                    f"times in {mutant.rel_path}",
                    file=sys.stderr,
                )
                return 3

            try:
                target_file.write_text(
                    base_text.replace(mutant.target, mutant.replacement, 1),
                    encoding="utf-8",
                )
                rc, combined = _run_pytest(scratch_root)
                if rc != 0 and "AssertionError" in combined:
                    status = "KILLED (AssertionError)"
                    killed += 1
                elif rc == 0:
                    status = "SURVIVED"
                    survived += 1
                else:
                    status = f"OTHER_FAIL (rc={rc})"
                    other_fails += 1
                print(f"{mutant.mutant_id:45s} | {status:24s} | {mutant.description}")
            finally:
                target_file.write_text(base_text, encoding="utf-8")

    for rel, expected_sha in worktree_shas.items():
        actual_sha = _sha256(WORKTREE_ROOT / rel)
        if actual_sha != expected_sha:
            print(f"[ERROR] Worktree file modified during mutant run: {rel}", file=sys.stderr)
            return 4

    print(
        f"\nSummary: total={len(MUTANTS)} killed={killed} survived={survived} "
        f"other_fails={other_fails}"
    )
    return 0 if (killed == len(MUTANTS) and survived == 0 and other_fails == 0) else 1


if __name__ == "__main__":
    raise SystemExit(main())
