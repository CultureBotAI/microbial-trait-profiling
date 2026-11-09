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
    Path(".github/copilot-instructions.md"),
    Path("scripts/sync_agent_docs.py"),
    Path("methods/guild_condition_body_site_heatmap.md"),
    Path("tests/test_sync_agent_docs.py"),
    Path("tests/test_biological_interpretation.py"),
)


@dataclass(frozen=True)
class Mutant:
    """Single source-text mutation on a staged file."""

    name: str
    rel_path: Path
    old: str
    new: str


MUTANTS: tuple[Mutant, ...] = (
    Mutant(
        name="M1_agents_missing_gh_pr_view_command_1343",
        rel_path=Path("AGENTS.md"),
        old="(`gh pr view N --json closingIssuesReferences`,\n   or its Development sidebar)",
        new="(its Development sidebar, or `closingIssuesReferences` in the GraphQL API)",
    ),
    Mutant(
        name="M2_agents_missing_sidebar_removal_rule_1343",
        rel_path=Path("AGENTS.md"),
        old=(
            "a link that comes only from a closing\n"
            "   keyword goes away when the PR body is edited, whereas a\n"
            "   Development-sidebar or branch link persists and must be removed in\n"
            "   the Development sidebar."
        ),
        new="a link can outlast a later edit of the body.",
    ),
    Mutant(
        name="M3_contributing_missing_gh_pr_view_command_1343",
        rel_path=Path("CONTRIBUTING.md"),
        old="(`gh pr view N --json closingIssuesReferences`,\nor the PR's Development sidebar)",
        new="(the PR's Development sidebar)",
    ),
    Mutant(
        name="M4_copilot_out_of_sync_1343",
        rel_path=Path(".github/copilot-instructions.md"),
        old="(`gh pr view N --json closingIssuesReferences`,\n   or its Development sidebar)",
        new="(its Development sidebar, or `closingIssuesReferences` in the GraphQL API)",
    ),
    Mutant(
        name="M5_guild_heatmap_revert_to_sampled_but_1343",
        rel_path=Path("methods/guild_condition_body_site_heatmap.md"),
        old='"sampled, with no recorded guild\nannotation signal"',
        new='"sampled, but no recorded guild\nannotation signal"',
    ),
)


class _OriginCheckPlugin:
    def __init__(self, scratch_root: Path) -> None:
        self.scratch_root = scratch_root.resolve()
        self.origin_ok = False

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
                    "-k",
                    (
                        "test_guild_methods_distinguish_annotation_zero_from_biological_absence or "
                        "test_contributor_docs_agree_on_package_paths_and_check_before_commit or "
                        "test_existing_copilot_generation_is_unchanged"
                    ),
                    "-q",
                    "--tb=short",
                ],
                plugins=[plugin],
            )
        )
    out = buf.getvalue()
    if plugin.origin_ok:
        out += "\nPYTEST_ORIGIN_OK\n"
    queue.put((rc, out))


def _run_pytest(scratch_root: Path) -> tuple[int, str]:
    ctx = multiprocessing.get_context("spawn")
    queue = ctx.Queue()
    proc = ctx.Process(target=_child_worker, args=(str(scratch_root), queue))
    proc.start()
    proc.join(timeout=60)
    if proc.is_alive():
        proc.terminate()
        proc.join()
        return (124, "Child timed out")
    if queue.empty():
        return (proc.exitcode or 1, "Child exited without result")
    rc, out = queue.get()
    return (int(rc), str(out))


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="verify_1343_") as tmp_dir:
        scratch_root = Path(tmp_dir)
        _stage_scratch_tree(scratch_root)

        base_rc, baseline_out = _run_pytest(scratch_root)
        if base_rc != 0 or "PYTEST_ORIGIN_OK" not in baseline_out:
            print("Baseline unmutated scratch run failed:", file=sys.stderr)
            print(baseline_out, file=sys.stderr)
            return 1
        print("Baseline unmutated scratch run passed (PYTEST_ORIGIN_OK).")

        for mutant in MUTANTS:
            target_file = scratch_root / mutant.rel_path
            original = (REPO_ROOT / mutant.rel_path).read_text(encoding="utf-8")
            count = original.count(mutant.old)
            if count != 1:
                print(
                    f"Mutant {mutant.name} target occurrence count={count} (expected 1)",
                    file=sys.stderr,
                )
                return 1
            mutated = original.replace(mutant.old, mutant.new, 1)
            target_file.write_text(mutated, encoding="utf-8")
            try:
                rc, combined = _run_pytest(scratch_root)
            finally:
                target_file.write_text(original, encoding="utf-8")

            killed = rc != 0 and "AssertionError" in combined and "PYTEST_ORIGIN_OK" in combined
            if not killed:
                print(f"SURVIVED or invalid failure for {mutant.name}:", file=sys.stderr)
                print(combined, file=sys.stderr)
                return 1
            print(f"KILLED {mutant.name} (rc={rc}, PYTEST_ORIGIN_OK)")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
