#!/usr/bin/env python3
"""Scratch-copy mutation verifier for #1921, #1922, #2005, and #2006.

Purpose:
  Verify that the slide stroke-floor and provenance inventory assertions in
  tests/test_cellulose_slide_figures.py (#1921) and the multi-condition
  palette/strip/legend assertions in tests/test_enrichment_infographic_followup.py
  (#1922, #2005, #2006) kill source-level mutations in
  src/mtp/cellulose_decr/slide_figures.py, src/mtp/cellulose_decr/figures.py,
  and scripts/analyze_enrichment_infographic_followup.py via AssertionError
  (rc=1, no KeyError) when executed against an isolated scratch copy with
  verified import origin inside pytest.

Inputs:
  - src/mtp/cellulose_decr/slide_figures.py
  - src/mtp/cellulose_decr/figures.py
  - scripts/analyze_enrichment_infographic_followup.py
  - tests/test_cellulose_slide_figures.py
  - tests/test_enrichment_infographic_followup.py

Outputs:
  - stdout summary of CONTROL (rc=0) and 12 killed mutants (rc=1, AssertionError,
    PYTEST_ORIGIN_OK).

Issue: #1921, #1922, #2005, #2006
Date: 2026-11-08
"""

from __future__ import annotations

import contextlib
import io
import multiprocessing
import os
import shutil
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SLIDE_FIGURES_REL = Path("src/mtp/cellulose_decr/slide_figures.py")
FIGURES_REL = Path("src/mtp/cellulose_decr/figures.py")
FOLLOWUP_REL = Path("scripts/analyze_enrichment_infographic_followup.py")
TEST_SLIDE_REL = Path("tests/test_cellulose_slide_figures.py")
TEST_FOLLOWUP_REL = Path("tests/test_enrichment_infographic_followup.py")

# pytest.raises reports a missing exception as "Failed: DID NOT RAISE", a test assertion.
ASSERTION_MARKERS = ("AssertionError", "Failed: DID NOT RAISE")
ERROR_MARKERS = (
    "KeyError",
    "ImportError",
    "ModuleNotFoundError",
    "SyntaxError",
    "IndentationError",
    "NameError",
    "TypeError",
    "AttributeError",
    "INTERNALERROR",
    "ERROR collecting",
)


@dataclass(frozen=True)
class Mutant:
    """One source-level mutation for #1921 or #1922."""

    name: str
    target_rel: Path
    old: str
    new: str
    test_rel: Path
    test_filter: str


MUTANTS: tuple[Mutant, ...] = (
    Mutant(
        name="M1_axis_strokes_below_floor_X4_1921",
        target_rel=SLIDE_FIGURES_REL,
        old=(
            "    ax.tick_params(labelsize=12, width=0.8)\n"
            "    for spine in ax.spines.values():\n"
            "        spine.set_linewidth(0.8)"
        ),
        new=(
            "    ax.tick_params(labelsize=12, width=0.3)\n"
            "    for spine in ax.spines.values():\n"
            "        spine.set_linewidth(0.3)"
        ),
        test_rel=TEST_SLIDE_REL,
        test_filter=(
            "test_every_slide_page_meets_native_slide_size_font_and_canvas_bounds and synthetic"
        ),
    ),
    Mutant(
        name="M2_spine_strokes_below_floor_1921",
        target_rel=SLIDE_FIGURES_REL,
        old="        spine.set_linewidth(0.8)",
        new="        spine.set_linewidth(0.3)",
        test_rel=TEST_SLIDE_REL,
        test_filter=(
            "test_every_slide_page_meets_native_slide_size_font_and_canvas_bounds and synthetic"
        ),
    ),
    Mutant(
        name="M3_tick_strokes_below_floor_1921",
        target_rel=SLIDE_FIGURES_REL,
        old="    ax.tick_params(labelsize=12, width=0.8)",
        new="    ax.tick_params(labelsize=12, width=0.3)",
        test_rel=TEST_SLIDE_REL,
        test_filter=(
            "test_every_slide_page_meets_native_slide_size_font_and_canvas_bounds and synthetic"
        ),
    ),
    Mutant(
        name="M4_provenance_drops_slide_inventory_X6_1921",
        target_rel=FIGURES_REL,
        old='        provenance["slide_pages"] = slide_inventory\n',
        new="",
        test_rel=TEST_SLIDE_REL,
        test_filter="test_slide_exports_are_opt_in_and_registered_with_captions_and_hashes",
    ),
    Mutant(
        name="M5_heatmap_palette_reverts_equal_lightness_1922",
        target_rel=FOLLOWUP_REL,
        old=('    "#c2410c",\n    "#60a5fa",\n    "#eab308",'),
        new=('    "#b91c1c",\n    "#2563eb",\n    "#047857",'),
        test_rel=TEST_FOLLOWUP_REL,
        test_filter="test_companion_heatmap_multi_condition_palette_and_legend_clear_dendrogram",
    ),
    Mutant(
        name="M6_heatmap_legend_overlaps_dendrogram_1922",
        target_rel=FOLLOWUP_REL,
        old='        loc="center right",\n        bbox_to_anchor=(-0.04, 0.5),',
        new='        loc="center",\n        bbox_to_anchor=(-0.08, 0.5),',
        test_rel=TEST_FOLLOWUP_REL,
        test_filter="test_companion_heatmap_multi_condition_palette_and_legend_clear_dendrogram",
    ),
    Mutant(
        name="M7_sixth_colour_reverts_to_7c3aed_2005",
        target_rel=FOLLOWUP_REL,
        old=('    "#ab7578",\n'),
        new=('    "#7c3aed",\n'),
        test_rel=TEST_FOLLOWUP_REL,
        test_filter="test_companion_heatmap_multi_condition_palette_and_legend_clear_dendrogram",
    ),
    Mutant(
        # #2336: the inline capacity check is now a shared helper call. Remove that
        # backstop and wrap the module-level palette index: the same ninth-condition
        # silent colour reuse must still fail the direct-render test.
        name="M8_guard_removed_wrap_restored_2005",
        target_rel=FOLLOWUP_REL,
        old=(
            "    validate_companion_heatmap_condition_count(len(TARGET_CONDITIONS))\n"
            "    condition_palette = {\n"
            "        condition: COMPANION_HEATMAP_PALETTE_COLORS[index]\n"
            "        for index, condition in enumerate(TARGET_CONDITIONS)\n"
            "    }"
        ),
        new=(
            "    condition_palette = {\n"
            "        condition: COMPANION_HEATMAP_PALETTE_COLORS[\n"
            "            index % len(COMPANION_HEATMAP_PALETTE_COLORS)\n"
            "        ]\n"
            "        for index, condition in enumerate(TARGET_CONDITIONS)\n"
            "    }"
        ),
        test_rel=TEST_FOLLOWUP_REL,
        test_filter="test_companion_heatmap_refuses_more_conditions_than_distinct_colours",
    ),
    Mutant(
        name="M9_fourth_colour_equals_trait_2005",
        target_rel=FOLLOWUP_REL,
        old=('    "#3b0764",\n'),
        new=('    "#155e75",\n'),
        test_rel=TEST_FOLLOWUP_REL,
        test_filter="test_companion_heatmap_multi_condition_palette_and_legend_clear_dendrogram",
    ),
    Mutant(
        name="M10_legend_over_heatmap_2006",
        target_rel=FOLLOWUP_REL,
        old='        loc="center right",\n        bbox_to_anchor=(-0.04, 0.5),',
        new='        loc="center right",\n        bbox_to_anchor=(4.0, 0.5),',
        test_rel=TEST_FOLLOWUP_REL,
        test_filter="test_companion_heatmap_multi_condition_palette_and_legend_clear_dendrogram",
    ),
    Mutant(
        # #2336: strict lookup replaced fallback mapping. Keep row validation, but
        # replace the strip lookup with one colour so it disagrees with the legend.
        name="M11_strip_ignores_palette_2006",
        target_rel=FOLLOWUP_REL,
        old=("    row_colors = row_conditions.map(lambda condition: condition_palette[condition])"),
        new=('    row_colors = row_conditions.map(lambda condition: "#155e75")'),
        test_rel=TEST_FOLLOWUP_REL,
        test_filter="test_companion_heatmap_multi_condition_palette_and_legend_clear_dendrogram",
    ),
    Mutant(
        name="M12_cyan_ramp_palette_2006",
        target_rel=FOLLOWUP_REL,
        old=('    "#c2410c",\n    "#60a5fa",\n    "#eab308",'),
        new=('    "#083344",\n    "#0891b2",\n    "#a5f3fc",'),
        test_rel=TEST_FOLLOWUP_REL,
        test_filter="test_companion_heatmap_multi_condition_palette_and_legend_clear_dendrogram",
    ),
)


def _setup_scratch(tmp_dir: Path) -> None:
    """Populate a minimal isolated source + test tree under tmp_dir."""
    shutil.copytree(REPO_ROOT / "src/mtp", tmp_dir / "src/mtp")
    (tmp_dir / "scripts").mkdir(parents=True, exist_ok=True)
    shutil.copy2(REPO_ROOT / FOLLOWUP_REL, tmp_dir / FOLLOWUP_REL)
    (tmp_dir / "tests").mkdir(parents=True, exist_ok=True)
    for rel in (
        TEST_SLIDE_REL,
        Path("tests/test_cellulose_nominal_labels.py"),
        TEST_FOLLOWUP_REL,
    ):
        shutil.copy2(REPO_ROOT / rel, tmp_dir / rel)
    (tmp_dir / "data/analysis").mkdir(parents=True, exist_ok=True)
    os.symlink(
        REPO_ROOT / "data/analysis/cellulose_decr",
        tmp_dir / "data/analysis/cellulose_decr",
    )
    expected_sf = (tmp_dir / SLIDE_FIGURES_REL).resolve()
    expected_fg = (tmp_dir / FIGURES_REL).resolve()
    expected_fu = (tmp_dir / FOLLOWUP_REL).resolve()
    conftest = tmp_dir / "tests/conftest.py"
    conftest.write_text(
        "import sys\n"
        "from pathlib import Path\n"
        "import mtp.cellulose_decr.figures as _fg\n"
        "import mtp.cellulose_decr.slide_figures as _sf\n\n"
        "def pytest_runtest_setup(item):\n"
        "    assert Path(_sf.__file__).resolve() == Path("
        f"{str(expected_sf)!r}).resolve()\n"
        "    assert Path(_fg.__file__).resolve() == Path("
        f"{str(expected_fg)!r}).resolve()\n"
        "    fu = sys.modules.get('enrichment_infographic_followup')\n"
        "    if fu is not None:\n"
        "        assert Path(fu.__file__).resolve() == Path("
        f"{str(expected_fu)!r}).resolve()\n"
        "    print('PYTEST_ORIGIN_OK', flush=True)\n",
        encoding="utf-8",
    )


def _worker(tmp_dir_str: str, test_rel_str: str, test_filter: str, out_file_str: str) -> None:
    """Run pytest in a clean child process rooted at tmp_dir."""
    tmp_dir = Path(tmp_dir_str)
    os.chdir(tmp_dir)
    sys.path.insert(0, str(tmp_dir / "tests"))
    sys.path.insert(0, str(tmp_dir / "src"))
    for name in [
        k
        for k in list(sys.modules)
        if k == "mtp"
        or k.startswith(("mtp.", "bsdbng", "enrichment_infographic_followup", "test_cellulose"))
    ]:
        sys.modules.pop(name, None)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        rc = int(
            pytest.main(
                [
                    "-q",
                    "-s",
                    "--tb=long",
                    "--rootdir",
                    str(tmp_dir),
                    str(tmp_dir / test_rel_str),
                    "-k",
                    test_filter,
                ]
            )
        )
    Path(out_file_str).write_text(f"RC={rc}\n{buf.getvalue()}", encoding="utf-8")


def _run_case(tmp_dir: Path, test_rel: Path, test_filter: str) -> tuple[int, str]:
    out_file = tmp_dir / "pytest_out.txt"
    ctx = multiprocessing.get_context("spawn")
    proc = ctx.Process(
        target=_worker,
        args=(str(tmp_dir), str(test_rel), test_filter, str(out_file)),
    )
    proc.start()
    proc.join()
    if proc.exitcode != 0 or not out_file.exists():
        raise RuntimeError(f"worker process failed with exitcode={proc.exitcode}")
    raw = out_file.read_text(encoding="utf-8")
    first, _, rest = raw.partition("\n")
    return int(first.removeprefix("RC=")), rest


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="mtp_1922_control_") as tmp:
        tmp_dir = Path(tmp)
        _setup_scratch(tmp_dir)
        for test_rel, test_filter in (
            (
                TEST_SLIDE_REL,
                (
                    "test_every_slide_page_meets_native_slide_size_font_and_canvas_bounds "
                    "and synthetic or "
                    "test_slide_exports_are_opt_in_and_registered_with_captions_and_hashes"
                ),
            ),
            (
                TEST_FOLLOWUP_REL,
                "test_companion_heatmap_multi_condition_palette_and_legend_clear_dendrogram",
            ),
        ):
            rc, out = _run_case(tmp_dir, test_rel, test_filter)
            if rc != 0 or "PYTEST_ORIGIN_OK" not in out:
                print(f"FAIL [CONTROL {test_rel}]: rc={rc}\n{out}")
                return 1
        print("CONTROL: rc=0 PYTEST_ORIGIN_OK")

    for mutant in MUTANTS:
        with tempfile.TemporaryDirectory(prefix=f"mtp_1922_{mutant.name}_") as tmp:
            tmp_dir = Path(tmp)
            _setup_scratch(tmp_dir)
            target_path = tmp_dir / mutant.target_rel
            source = target_path.read_text(encoding="utf-8")
            count = source.count(mutant.old)
            if count != 1:
                print(f"FAIL [{mutant.name}]: expected 1 occurrence, found {count}")
                return 1
            target_path.write_text(source.replace(mutant.old, mutant.new, 1), encoding="utf-8")
            rc, out = _run_case(tmp_dir, mutant.test_rel, mutant.test_filter)
            if rc != 1:
                print(f"FAIL [{mutant.name}]: expected rc=1, got rc={rc}\n{out}")
                return 1
            if "PYTEST_ORIGIN_OK" not in out:
                print(f"FAIL [{mutant.name}]: missing PYTEST_ORIGIN_OK\n{out}")
                return 1
            if any(err in out for err in ERROR_MARKERS):
                print(f"FAIL [{mutant.name}]: non-assertion error marker in output\n{out}")
                return 1
            if not any(marker in out for marker in ASSERTION_MARKERS):
                print(f"FAIL [{mutant.name}]: no assertion marker in output\n{out}")
                return 1
            print(f"KILLED {mutant.name}: rc=1 AssertionError PYTEST_ORIGIN_OK")

    print(f"ALL {len(MUTANTS)} #1921/#1922 MUTANTS KILLED CLEANLY.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
