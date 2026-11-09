"""Scratch-copy mutant verifier for #1314, #2084, #2085, #2087 (historical audit rows).

Purpose:
    Verify that running the full `tests/test_biological_interpretation.py` suite
    kills all 13 mutations (including the 5 meaning-reversing mutations from
    #2084 and the PR #173 intro / P10 scoping mutations from #2085) on
    `methods/biological_interpretation_audit.md` in an isolated temporary
    scratch copy without mutating the live worktree, and that every failure is
    an assertion failure (#2087).
Inputs:
    methods/biological_interpretation_audit.md
    tests/test_biological_interpretation.py
Outputs:
    Per-mutant pass/fail summary on stdout; exits non-zero if baseline fails,
    origin verification fails, any failure is not an AssertionError, or any
    mutant survives.
Issue:
    #1314, #2084, #2085, #2087
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
TARGET_REL = Path("methods/biological_interpretation_audit.md")
TEST_REL = Path("tests/test_biological_interpretation.py")

SCRATCH_FILES: tuple[Path, ...] = (
    Path("AGENTS.md"),
    Path("CONTRIBUTING.md"),
    Path("MANUSCRIPT_AUDIT.md"),
    Path("README.md"),
    Path(".github/PULL_REQUEST_TEMPLATE.md"),
    Path("docs/table1_trait_semantics_audit.md"),
    Path("methods/README.md"),
    Path("methods/biological_interpretation.md"),
    Path("methods/biological_interpretation_audit.md"),
    Path("methods/delta_z_derived.md"),
    Path("methods/guild_condition_body_site_heatmap.md"),
    Path("methods/table1_trait_condition_body_site_maps.md"),
    Path("methods/three_arm_statistics.md"),
    Path("methods/trait_key_policy.md"),
    Path("tests/test_biological_interpretation.py"),
)


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
    Mutant(
        name="M7_p09_subsequently_closed_becomes_still_open_2084",
        old="(subsequently closed; see the",
        new="(still open; see the",
    ),
    Mutant(
        name="M8_c04_work_remains_open_2084",
        old="calibration and trait-semantics work; see the",
        new="calibration and trait-semantics work, which remains open; see the",
    ),
    Mutant(
        name="M9_c04_appended_remain_open_2084",
        old="or validate a mechanism). |",
        new="or validate a mechanism). #165, #182 and #189 remain open. |",
    ),
    Mutant(
        name="M10_adoption_closures_do_not_supersede_2084",
        old="The dated tracker closures above supersede",
        new="The dated tracker closures above do not supersede",
    ),
    Mutant(
        name="M11_p02_rewritten_as_open_requirements_2084",
        old="#131/#182 limitations remain.",
        new="#131/#132/#182 are still open requirements.",
    ),
    Mutant(
        name="M12_overlay_intro_drops_2026_09_14_review_state_2085",
        old=(
            "Issue references in the\n"
            '"Owner / remaining work" column record the 2026-09-14 review state, and the\n'
            "[status checkpoint](#status-checkpoint-2026-10-30) dates their closures."
        ),
        new="The Owner / remaining work column lists open items.",
    ),
    Mutant(
        name="M13_p10_restores_closed_173_2085",
        old="reconcile the generated legend and actual mark provenance in #172 before",
        new="reconcile the generated legend and actual mark provenance in #172/#173 before",
    ),
)


class _OriginAndAssertionPlugin:
    def __init__(self, scratch_root: Path) -> None:
        self.scratch_root = scratch_root.resolve()
        self.origin_ok = False
        self.passed = 0
        self.assertion_failures = 0
        self.other_failures = 0

    def pytest_runtest_logreport(self, report: Any) -> None:
        if report.when == "call":
            if report.passed:
                self.passed += 1
            elif report.failed:
                longrepr = str(getattr(report, "longreprtext", "") or "")
                if re.search(r"AssertionError|^E\s+assert\b", longrepr, re.MULTILINE):
                    self.assertion_failures += 1
                else:
                    self.other_failures += 1
        elif report.failed:
            self.other_failures += 1

    def pytest_sessionfinish(self, session: Any) -> None:
        items = getattr(session, "items", [])
        if not items:
            return
        for item in items:
            item_path = Path(str(item.fspath)).resolve()
            if not item_path.is_relative_to(self.scratch_root):
                return
            mod = getattr(item, "module", None)
            if mod is None:
                return
            if Path(getattr(mod, "__file__", "")).resolve() != (self.scratch_root / TEST_REL):
                return
            if Path(getattr(mod, "BIOLOGICAL_INTERPRETATION_AUDIT_MD", "")).resolve() != (
                self.scratch_root / TARGET_REL
            ):
                return
        self.origin_ok = True


def _stage_scratch_tree(scratch_root: Path) -> None:
    for rel in SCRATCH_FILES:
        src = REPO_ROOT / rel
        dst = scratch_root / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


def _child_worker(scratch_root_str: str, queue: Any) -> None:
    scratch_root = Path(scratch_root_str).resolve()
    src_dir = str(REPO_ROOT / "src")
    if src_dir not in sys.path:
        sys.path.insert(0, src_dir)
    plugin = _OriginAndAssertionPlugin(scratch_root)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        rc = int(
            pytest.main(
                [
                    "-o",
                    "addopts=",
                    str(scratch_root / TEST_REL),
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
            f"assert_fails={plugin.assertion_failures} other_fails={plugin.other_failures}\n"
        )
    queue.put((rc, plugin.passed, plugin.assertion_failures, plugin.other_failures, out))


def _run_pytest(scratch_root: Path) -> tuple[int, int, int, int, str]:
    ctx = multiprocessing.get_context("spawn")
    queue = ctx.Queue()
    proc = ctx.Process(target=_child_worker, args=(str(scratch_root), queue))
    proc.start()
    proc.join(timeout=60)
    if proc.is_alive():
        proc.terminate()
        proc.join()
        return (124, 0, 0, 1, "Child timed out")
    if queue.empty():
        return (proc.exitcode or 1, 0, 0, 1, "Child exited without result")
    rc, passed, assert_fails, other_fails, out = queue.get()
    return (int(rc), int(passed), int(assert_fails), int(other_fails), str(out))


def main() -> int:
    source_text = (REPO_ROOT / TARGET_REL).read_text(encoding="utf-8")

    with tempfile.TemporaryDirectory(prefix="verify_1314_") as tmp_dir:
        scratch_root = Path(tmp_dir)
        _stage_scratch_tree(scratch_root)

        base_rc, base_passed, base_af, base_of, baseline_out = _run_pytest(scratch_root)
        if (
            base_rc != 0
            or base_passed < 30
            or base_af != 0
            or base_of != 0
            or "PYTEST_ORIGIN_OK" not in baseline_out
        ):
            print("Baseline unmutated scratch run failed:", file=sys.stderr)
            print(baseline_out, file=sys.stderr)
            return 1
        print(
            "Baseline unmutated scratch run passed "
            f"(PYTEST_ORIGIN_OK, passed={base_passed}, "
            f"assert_fails={base_af}, other_fails={base_of})."
        )

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

            rc, passed, assert_fails, other_fails, combined = _run_pytest(scratch_root)
            killed = (
                rc != 0
                and assert_fails >= 1
                and other_fails == 0
                and "PYTEST_ORIGIN_OK" in combined
            )
            if not killed:
                print(f"SURVIVED or invalid failure for {mutant.name}:", file=sys.stderr)
                print(combined, file=sys.stderr)
                return 1
            print(
                f"KILLED {mutant.name} (rc={rc}, passed={passed}, "
                f"assert_fails={assert_fails}, other_fails={other_fails}, PYTEST_ORIGIN_OK)"
            )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
