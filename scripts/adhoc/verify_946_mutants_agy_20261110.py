"""Isolated scratch-copy mutation verifier for issue #946.

Purpose:
    Prove that tests/test_study_trait_profile_heatmap.py kills >= 6 distinct
    mutations of scripts/create_study_trait_profile_heatmap.py with
    AssertionError while keeping the live worktree untouched.
Inputs:
    scripts/create_study_trait_profile_heatmap.py
    tests/test_study_trait_profile_heatmap.py
Outputs:
    Per-mutant pass/kill report on stdout; exits 0 iff control passes and all
    mutants fail with AssertionError.
Serves:
    Issue #946
Date:
    2026-11-10
"""

from __future__ import annotations

import hashlib
import shutil
import sys
import tempfile
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
TARGET_REL = Path("scripts/create_study_trait_profile_heatmap.py")
TEST_REL = Path("tests/test_study_trait_profile_heatmap.py")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _replace_nth(text: str, old: str, new: str, occurrence: int) -> str:
    parts = text.split(old)
    if len(parts) <= occurrence:
        raise ValueError(f"Expected occurrence {occurrence} of {old!r}")
    before = old.join(parts[:occurrence])
    after = old.join(parts[occurrence:])
    return f"{before}{new}{after}"


MUTANTS: tuple[tuple[str, str, str, int], ...] = (
    (
        "M01_cell_border_linewidth_0_05",
        "CELL_BORDER_LINEWIDTH = 0.25",
        "CELL_BORDER_LINEWIDTH = 0.05",
        1,
    ),
    (
        "M02_save_heatmap_linewidth_0_05",
        "        linewidths=CELL_BORDER_LINEWIDTH,",
        "        linewidths=0.05,",
        1,
    ),
    (
        "M03_save_paired_heatmap_linewidth_0_05",
        "        linewidths=CELL_BORDER_LINEWIDTH,",
        "        linewidths=0.05,",
        2,
    ),
    (
        "M04_palette_uses_set2",
        "selected_colors: list[str | tuple[float, float, float]] = list(rotated[: len(ordered)])",
        (
            "selected_colors: list[str | tuple[float, float, float]] = "
            'list(sns.color_palette("Set2", n_colors=len(ordered)))'
        ),
        1,
    ),
    (
        "M05_palette_collides_with_other_grey",
        '    "#D55E00",',
        '    "#b3b3b3",',
        1,
    ),
    (
        "M06_palette_deuteranopia_collision",
        '    "#D55E00",',
        '    "#786418",',
        1,
    ),
    (
        "M07_condition_rotation_zero",
        '    "condition": 2,',
        '    "condition": 0,',
        1,
    ),
    (
        "M08_other_color_swapped",
        '    colors["Other"] = "#b0b0b0"',
        '    colors["Other"] = "#e5e7eb"',
        1,
    ),
)


class _AssertionCapturePlugin:
    def __init__(self) -> None:
        self.failed_reports: list[str] = []

    def pytest_runtest_logreport(self, report: pytest.TestReport) -> None:
        if report.when in ("call", "setup") and report.failed:
            self.failed_reports.append(str(report.longrepr))


def _stage_scratch(scratch: Path) -> None:
    for rel in (TARGET_REL, TEST_REL):
        src_path = REPO_ROOT / rel
        dst_path = scratch / rel
        dst_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src_path, dst_path)
    conftest_path = scratch / "tests" / "conftest.py"
    conftest_path.write_text(
        "from pathlib import Path\n"
        "import tests.test_study_trait_profile_heatmap as _test_mod\n\n"
        "def pytest_sessionfinish(session, exitstatus):\n"
        "    scratch_root = Path(__file__).resolve().parents[1]\n"
        "    origin = Path(_test_mod.heatmap.__file__).resolve()\n"
        "    if not origin.is_relative_to(scratch_root):\n"
        "        raise RuntimeError(f'Loaded from outside scratch: {origin}')\n",
        encoding="utf-8",
    )


def _run_pytest(scratch: Path) -> tuple[int, list[str]]:
    saved_path = list(sys.path)
    transient_names = {
        "create_study_trait_profile_heatmap",
        "test_study_trait_profile_heatmap",
        "conftest",
        "tests",
        "tests.test_study_trait_profile_heatmap",
        "tests.conftest",
    }
    saved_modules = {name: mod for name, mod in sys.modules.items() if name in transient_names}
    for name in transient_names:
        sys.modules.pop(name, None)
    sys.path.insert(0, str(scratch))
    sys.path.insert(1, str(REPO_ROOT / "src"))
    plugin = _AssertionCapturePlugin()
    try:
        exit_code = int(
            pytest.main(
                [
                    "-q",
                    "--tb=long",
                    "-o",
                    f"cache_dir={scratch / '.pytest_cache'}",
                    str(scratch / TEST_REL),
                    "-k",
                    (
                        "test_rendered_heatmaps_enforce_minimum_cell_border_linewidth "
                        "or test_palette_for_categories_enforces_accessible_pairwise_"
                        "color_separation "
                        "or test_annotation_track_palettes_use_distinct_rotations"
                    ),
                ],
                plugins=[plugin],
            )
        )
    finally:
        for name in transient_names:
            sys.modules.pop(name, None)
        sys.modules.update(saved_modules)
        sys.path[:] = saved_path
    return exit_code, plugin.failed_reports


def main() -> int:
    target_path = REPO_ROOT / TARGET_REL
    test_path = REPO_ROOT / TEST_REL
    orig_target_sha = _sha256(target_path)
    orig_test_sha = _sha256(test_path)
    original_source = target_path.read_text(encoding="utf-8")

    with tempfile.TemporaryDirectory(prefix="verify_946_") as tmp_dir:
        base = Path(tmp_dir)
        control_dir = base / "control"
        _stage_scratch(control_dir)
        rc, failures = _run_pytest(control_dir)
        if rc != 0 or failures:
            print(f"CONTROL FAILED (rc={rc}):\n" + "\n".join(failures))
            return 1
        print("CONTROL: rc=0 PYTEST_ORIGIN_OK (9 passed)")

        killed = 0
        for name, old, new, occurrence in MUTANTS:
            mutated_source = _replace_nth(original_source, old, new, occurrence)
            if mutated_source == original_source:
                print(f"MUTANT {name}: replacement did not change source")
                return 1
            mut_dir = base / name
            _stage_scratch(mut_dir)
            (mut_dir / TARGET_REL).write_text(mutated_source, encoding="utf-8")
            m_rc, m_failures = _run_pytest(mut_dir)
            joined = "\n".join(m_failures)
            has_assertion = "AssertionError" in joined
            print(f"{name}: rc={m_rc} assertion={has_assertion}")
            if m_rc == 0 or not has_assertion:
                print(f"SURVIVED OR NON-ASSERTION FAILURE:\n{joined}")
                return 1
            killed += 1

    if _sha256(target_path) != orig_target_sha or _sha256(test_path) != orig_test_sha:
        print("ERROR: Canonical worktree files were modified during verification!")
        return 1

    print(f"All {killed}/{len(MUTANTS)} mutants killed by AssertionError; worktree untouched.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
