"""Scratch-copy mutant verifier for #1314 (historical audit rows C04 and P09).

Purpose:
    Verify that tests/test_biological_interpretation.py kills mutations to
    methods/biological_interpretation_audit.md for C04, P09, Status checkpoint,
    and Adoption record in an isolated temporary scratch copy without mutating
    the live worktree.
Inputs:
    methods/biological_interpretation_audit.md
    tests/test_biological_interpretation.py
Outputs:
    Per-mutant pass/fail summary on stdout; exits non-zero if baseline fails,
    origin verification fails, or any mutant survives.
Issue:
    #1314
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
TARGET_REL = Path("methods/biological_interpretation_audit.md")
TEST_REL = Path("tests/test_biological_interpretation.py")


@dataclass(frozen=True)
class Mutant:
    """Single source-text mutation on methods/biological_interpretation_audit.md."""

    name: str
    old: str
    new: str


MUTANTS: tuple[Mutant, ...] = (
    Mutant(
        name="M1_restore_c04_open_work_wording_1314",
        old=(
            "At the original audit, #165/#182 and #189 tracked the statistical-dependence, "
            "calibration and trait-semantics work; see the "
            "[status checkpoint](#status-checkpoint-2026-10-30) for their subsequent closures "
            "(closing them does not establish semantic equivalence or validate a mechanism)."
        ),
        new="#165/#182 and #189 retain their work.",
    ),
    Mutant(
        name="M2_restore_p09_open_work_wording_1314",
        old=(
            "in the historical PR #173 review, #131/#132/#182 were separate statistical "
            "requirements (subsequently closed; see the "
            "[status checkpoint](#status-checkpoint-2026-10-30); their closure does not "
            "validate the biological explanation)."
        ),
        new="#131/#132/#182 remain independent requirements.",
    ),
    Mutant(
        name="M3_drop_c04_non_validation_caveat_1314",
        old="closing them does not establish semantic equivalence or validate a mechanism",
        new="closing them establishes semantic equivalence",
    ),
    Mutant(
        name="M4_drop_p09_non_validation_caveat_1314",
        old="their closure does not validate the biological explanation",
        new="their closure validates the biological explanation",
    ),
    Mutant(
        name="M5_omit_131_132_165_from_status_checkpoint_1314",
        old=(
            "| [#131](https://github.com/CultureBotAI/microbial-trait-profiling/issues/131), "
            "enrichment sampling unit | Closed | 2026-10-07 |\n"
            "| [#132](https://github.com/CultureBotAI/microbial-trait-profiling/issues/132), "
            "depletion structural zeros | Closed | 2026-10-07 |\n"
            "| [#165](https://github.com/CultureBotAI/microbial-trait-profiling/issues/165), "
            "guild/Table 1 dependence and ties | Closed | 2026-10-13 |\n"
        ),
        new="",
    ),
    Mutant(
        name="M6_omit_131_132_165_from_adoption_record_1314",
        old="for #131, #132, #165, #194/#190, #182 and #188-#191.",
        new="for #194/#190, #182 and #188-#191.",
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
            if mod is not None and hasattr(mod, "BIOLOGICAL_INTERPRETATION_AUDIT_MD"):
                target_path = Path(mod.BIOLOGICAL_INTERPRETATION_AUDIT_MD).resolve()
                if not target_path.is_relative_to(self.scratch_root):
                    return
        self.origin_ok = True


def _stage_scratch_tree(scratch_root: Path) -> None:
    for rel in (TARGET_REL, TEST_REL):
        src = REPO_ROOT / rel
        dst = scratch_root / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


def _child_worker(scratch_root_str: str, queue: Any) -> None:
    scratch_root = Path(scratch_root_str).resolve()
    plugin = _OriginCheckPlugin(scratch_root)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        rc = int(
            pytest.main(
                [
                    str(scratch_root / TEST_REL),
                    "-k",
                    (
                        "test_audit_separates_dated_tracker_closures_from_biological_validation or "
                        "test_c04_and_p09_date_scope_closed_issues_without_claiming_biological_validation"
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
    source_text = (REPO_ROOT / TARGET_REL).read_text(encoding="utf-8")

    with tempfile.TemporaryDirectory(prefix="verify_1314_") as tmp_dir:
        scratch_root = Path(tmp_dir)
        _stage_scratch_tree(scratch_root)

        base_rc, baseline_out = _run_pytest(scratch_root)
        if base_rc != 0 or "PYTEST_ORIGIN_OK" not in baseline_out:
            print("Baseline unmutated scratch run failed:", file=sys.stderr)
            print(baseline_out, file=sys.stderr)
            return 1
        print("Baseline unmutated scratch run passed (PYTEST_ORIGIN_OK).")

        for mutant in MUTANTS:
            count = source_text.count(mutant.old)
            if count != 1:
                print(
                    f"Mutant {mutant.name} target occurrence count={count} (expected 1)",
                    file=sys.stderr,
                )
                return 1
            mutated = source_text.replace(mutant.old, mutant.new, 1)
            (scratch_root / TARGET_REL).write_text(mutated, encoding="utf-8")

            rc, combined = _run_pytest(scratch_root)
            killed = rc != 0 and "AssertionError" in combined and "PYTEST_ORIGIN_OK" in combined
            if not killed:
                print(f"SURVIVED or invalid failure for {mutant.name}:", file=sys.stderr)
                print(combined, file=sys.stderr)
                return 1
            print(f"KILLED {mutant.name} (rc={rc}, PYTEST_ORIGIN_OK)")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
