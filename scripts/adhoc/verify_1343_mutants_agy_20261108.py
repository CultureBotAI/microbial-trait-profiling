"""Scratch-copy mutant verifier for #1343 (closing-link docs and zero-signal wording).

Purpose:
    Verify that tests/test_sync_agent_docs.py and tests/test_biological_interpretation.py
    kill mutations to AGENTS.md, CONTRIBUTING.md, .github/copilot-instructions.md, and
    methods/guild_condition_body_site_heatmap.md in an isolated temporary scratch copy
    without mutating the live worktree.
Inputs:
    AGENTS.md
    CONTRIBUTING.md
    .github/copilot-instructions.md
    methods/guild_condition_body_site_heatmap.md
    tests/test_sync_agent_docs.py
    tests/test_biological_interpretation.py
Outputs:
    Per-mutant pass/fail summary on stdout; exits non-zero if baseline fails,
    origin verification fails, or any mutant survives.
Issue:
    #1343
Date:
    2026-11-08
"""

from __future__ import annotations

import contextlib
import io
import multiprocessing
import re
import shutil
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
STAGED_RELS: tuple[Path, ...] = (
    Path("pyproject.toml"),
    Path("AGENTS.md"),
    Path("CONTRIBUTING.md"),
    Path("README.md"),
    Path("MANUSCRIPT_AUDIT.md"),
    Path("docs/table1_trait_semantics_audit.md"),
    Path(".github/PULL_REQUEST_TEMPLATE.md"),
    Path(".github/copilot-instructions.md"),
    Path("scripts/sync_agent_docs.py"),
    Path("tests/test_sync_agent_docs.py"),
    Path("tests/test_biological_interpretation.py"),
)
AG = Path("AGENTS.md")
CP = Path(".github/copilot-instructions.md")
CO = Path("CONTRIBUTING.md")
HEAT = Path("methods/guild_condition_body_site_heatmap.md")

SIDEBAR_AG = "   a maintainer). A closing keyword in a commit message (and so, by default,\n"
SIDEBAR_CO = "must ask a maintainer). A closing keyword in a\n"
ALSO = "A Development-sidebar link also goes away when the PR body is edited. "


@dataclass(frozen=True)
class Mutant:
    """Source-text mutation across one or more staged files."""

    name: str
    edits: tuple[tuple[Path, str, str], ...]


MUTANTS: tuple[Mutant, ...] = (
    Mutant(
        name="M1_agents_missing_gh_pr_view_command_1343",
        edits=(
            (
                AG,
                "(`gh pr view N --json closingIssuesReferences`,\n"
                "   which lists issues only, never PRs, or the PR's Development sidebar)",
                "(the PR's Development sidebar)",
            ),
        ),
    ),
    Mutant(
        name="M2_agents_missing_sidebar_removal_rule_1343",
        edits=(
            (
                AG,
                "a link that comes only from a closing keyword goes away once the\n"
                "   keyword is removed from the PR body, whereas a Development-sidebar or\n"
                "   branch link persists and must be removed in the PR's Development\n"
                "   sidebar in the web UI (no `gh` subcommand does this, and only a linked\n"
                "   branch can be unlinked through the GraphQL API, so an agent must ask\n"
                "   a maintainer).",
                "a link can outlast a later edit of the body.",
            ),
        ),
    ),
    Mutant(
        name="M3_contributing_missing_gh_pr_view_command_1343",
        edits=(
            (
                CO,
                "(`gh pr view N --json closingIssuesReferences`,\n"
                "which lists issues only, never PRs, or the PR's Development sidebar)",
                "(the PR's Development sidebar)",
            ),
        ),
    ),
    Mutant(
        name="M4_copilot_out_of_sync_1343",
        edits=(
            (
                CP,
                "(`gh pr view N --json closingIssuesReferences`,\n"
                "   which lists issues only, never PRs, or the PR's Development sidebar)",
                "(the PR's Development sidebar)",
            ),
        ),
    ),
    Mutant(
        name="M5_guild_heatmap_revert_to_sampled_but_1343",
        edits=(
            (
                HEAT,
                '"sampled, with no recorded guild\nannotation signal"',
                '"sampled, but no recorded guild\nannotation signal"',
            ),
        ),
    ),
    Mutant(
        name="M6_sidebar_link_also_goes_away_on_body_edit_2097",
        edits=(
            (AG, SIDEBAR_AG, SIDEBAR_AG.replace("a maintainer). ", f"a maintainer). {ALSO}")),
            (CP, SIDEBAR_AG, SIDEBAR_AG.replace("a maintainer). ", f"a maintainer). {ALSO}")),
            (CO, SIDEBAR_CO, SIDEBAR_CO.replace("a maintainer). ", f"a maintainer). {ALSO}")),
        ),
    ),
    Mutant(
        name="M7_keep_links_to_partly_addressed_issues_2097",
        edits=(
            (
                AG,
                "Remove any link to an issue the PR only partly",
                "Keep any link to an issue the PR only partly",
            ),
            (
                CP,
                "Remove any link to an issue the PR only partly",
                "Keep any link to an issue the PR only partly",
            ),
            (
                CO,
                "Remove any link to an issue the PR\nonly partly",
                "Keep any link to an issue the PR\nonly partly",
            ),
        ),
    ),
    Mutant(
        name="M8_use_gh_issue_develop_2097",
        edits=(
            (
                AG,
                "   Do not use `gh issue develop`. It records",
                "   Use `gh issue develop`. It records",
            ),
            (
                CP,
                "   Do not use `gh issue develop`. It records",
                "   Use `gh issue develop`. It records",
            ),
            (CO, "Do not use `gh issue develop`: its", "Use `gh issue develop`: its"),
        ),
    ),
    Mutant(
        name="M9_heatmap_zero_signal_means_absent_2097",
        edits=(
            (
                HEAT,
                'annotation signal" distinguishes mapped zero signal from unsampled cells.',
                'annotation signal" shows the guild is absent at those sites.',
            ),
        ),
    ),
    Mutant(
        name="M10_missing_commit_message_and_negated_keyword_rule_2096",
        edits=(
            (
                AG,
                "A closing keyword in a commit message (and so, by default,\n"
                "   in the squash message) also closes the issue on merge, appears in neither\n"
                "   check, and is not removed by editing the PR body: write `Refs #<number>`\n"
                '   in commits, and never a negated keyword such as "does not close #<number>"\n'
                "   (which still closes). ",
                "",
            ),
            (
                CP,
                "A closing keyword in a commit message (and so, by default,\n"
                "   in the squash message) also closes the issue on merge, appears in neither\n"
                "   check, and is not removed by editing the PR body: write `Refs #<number>`\n"
                '   in commits, and never a negated keyword such as "does not close #<number>"\n'
                "   (which still closes). ",
                "",
            ),
            (
                CO,
                "A closing keyword in a\n"
                "commit message (and so, by default, in the squash message) also closes the "
                "issue on merge,\n"
                "appears in neither check, and is not removed by editing the PR body: write\n"
                '`Refs #<number>` in commits, and never a negated keyword such as "does not\n'
                'close #<number>" (which still closes). ',
                "",
            ),
        ),
    ),
)


class _OriginCheckPlugin:
    def __init__(self, scratch_root: Path) -> None:
        self.scratch_root = scratch_root.resolve()
        self.origin_ok = False
        self.passed = 0
        self.assert_fails = 0
        self.other_fails = 0

    def pytest_runtest_logreport(self, report: Any) -> None:
        if report.when == "call" and report.passed:
            self.passed += 1
        elif report.failed:
            text = str(getattr(report, "longreprtext", ""))
            if re.search(r"AssertionError|^E\s+assert\b", text, re.MULTILINE):
                self.assert_fails += 1
            else:
                self.other_fails += 1

    def pytest_sessionfinish(self, session: Any) -> None:
        items = getattr(session, "items", [])
        if not items:
            return
        for item in items:
            item_path = Path(str(item.fspath)).resolve()
            if not item_path.is_relative_to(self.scratch_root):
                return
            mod = getattr(item, "module", None)
            if mod is not None and hasattr(mod, "REPO_ROOT"):
                root_path = Path(mod.REPO_ROOT).resolve()
                if not root_path.is_relative_to(self.scratch_root):
                    return
        self.origin_ok = True


def _stage_scratch_tree(scratch_root: Path) -> None:
    shutil.copytree(REPO_ROOT / "methods", scratch_root / "methods")
    for rel in STAGED_RELS:
        src = REPO_ROOT / rel
        dst = scratch_root / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


def _child_worker(scratch_root_str: str, queue: Any) -> None:
    scratch_root = Path(scratch_root_str).resolve()
    sys.path.insert(0, str(REPO_ROOT / "src"))
    plugin = _OriginCheckPlugin(scratch_root)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        rc = int(
            pytest.main(
                [
                    str(scratch_root / "tests/test_biological_interpretation.py"),
                    str(scratch_root / "tests/test_sync_agent_docs.py"),
                    "-q",
                    "--tb=short",
                ],
                plugins=[plugin],
            )
        )
    out = buf.getvalue()
    if plugin.origin_ok:
        out += (
            f"\nPYTEST_ORIGIN_OK passed={plugin.passed} "
            f"assert_fails={plugin.assert_fails} other_fails={plugin.other_fails}\n"
        )
    queue.put((rc, out, plugin.passed, plugin.assert_fails, plugin.other_fails))


def _run_pytest(scratch_root: Path) -> tuple[int, str, int, int, int]:
    ctx = multiprocessing.get_context("spawn")
    queue = ctx.Queue()
    proc = ctx.Process(target=_child_worker, args=(str(scratch_root), queue))
    proc.start()
    proc.join(timeout=60)
    if proc.is_alive():
        proc.terminate()
        proc.join()
        return (124, "Child timed out", 0, 0, 1)
    if queue.empty():
        return (proc.exitcode or 1, "Child exited without result", 0, 0, 1)
    rc, out, passed, assert_fails, other_fails = queue.get()
    return (int(rc), str(out), int(passed), int(assert_fails), int(other_fails))


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="verify_1343_") as tmp_dir:
        scratch_root = Path(tmp_dir)
        _stage_scratch_tree(scratch_root)

        base_rc, baseline_out, passed, assert_fails, other_fails = _run_pytest(scratch_root)
        if (
            base_rc != 0
            or "PYTEST_ORIGIN_OK" not in baseline_out
            or assert_fails != 0
            or other_fails != 0
        ):
            print("Baseline unmutated scratch run failed:", file=sys.stderr)
            print(baseline_out, file=sys.stderr)
            return 1
        print(
            f"Baseline unmutated scratch run passed "
            f"(PYTEST_ORIGIN_OK, passed={passed}, assert_fails=0, other_fails=0)."
        )

        for mutant in MUTANTS:
            originals: dict[Path, str] = {}
            for rel_path, old, new in mutant.edits:
                target_file = scratch_root / rel_path
                current = target_file.read_text(encoding="utf-8")
                originals.setdefault(target_file, current)
                count = current.count(old)
                if count != 1:
                    print(
                        f"Mutant {mutant.name} target occurrence count={count} in {rel_path} "
                        "(expected 1)",
                        file=sys.stderr,
                    )
                    return 1
                target_file.write_text(current.replace(old, new, 1), encoding="utf-8")
            try:
                rc, combined, _, m_assert_fails, m_other_fails = _run_pytest(scratch_root)
            finally:
                for target_file, original_text in originals.items():
                    target_file.write_text(original_text, encoding="utf-8")

            killed = (
                rc != 0
                and m_assert_fails > 0
                and m_other_fails == 0
                and "PYTEST_ORIGIN_OK" in combined
            )
            if not killed:
                print(
                    f"SURVIVED or invalid failure for {mutant.name} "
                    f"(rc={rc}, assert_fails={m_assert_fails}, other_fails={m_other_fails}):",
                    file=sys.stderr,
                )
                print(combined, file=sys.stderr)
                return 1
            print(
                f"KILLED {mutant.name} "
                f"(rc={rc}, assert_fails={m_assert_fails}, other_fails=0, PYTEST_ORIGIN_OK)"
            )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
