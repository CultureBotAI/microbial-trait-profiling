"""Require every real infographic mutation anchor to match current source (#2336)."""

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
VERIFIERS = (
    ("verify_2035_2047_mutants_agy_20261111.py", "rel_path"),
    ("verify_1922_mutants_agy_20261108.py", "target_rel"),
)
# #2350: a stale anchor must be repaired, not removed with its mutant.
EXPECTED_MUTANTS = {
    "verify_2035_2047_mutants_agy_20261111.py": (
        "M1_restore_unknown_grey_fallback_2035",
        "M2_allow_unlisted_non_unknown_condition_2035",
        "M3_skip_unlisted_check_on_dominant_condition_column_2035",
        "M4_skip_main_paired_top_early_guard_2047",
        "M5_skip_main_target_conditions_early_guard_2047",
        "M6_mkdir_before_condition_count_guard_in_main_2047",
        "M7_widen_condition_count_limit_to_nine_2047",
        "M8_remove_heatmap_backstop_call_2047",
    ),
    "verify_1922_mutants_agy_20261108.py": (
        "M1_axis_strokes_below_floor_X4_1921",
        "M2_spine_strokes_below_floor_1921",
        "M3_tick_strokes_below_floor_1921",
        "M4_provenance_drops_slide_inventory_X6_1921",
        "M5_heatmap_palette_reverts_equal_lightness_1922",
        "M6_heatmap_legend_overlaps_dendrogram_1922",
        "M7_sixth_colour_reverts_to_7c3aed_2005",
        "M8_guard_removed_wrap_restored_2005",
        "M9_fourth_colour_equals_trait_2005",
        "M10_legend_over_heatmap_2006",
        "M11_strip_ignores_palette_2006",
        "M12_cyan_ramp_palette_2006",
    ),
}


def _anchor_cases() -> list[tuple[str, Path, str]]:
    """Read the verifiers' specifications without invoking their mutation runners."""
    cases: list[tuple[str, Path, str]] = []
    for filename, target_attribute in VERIFIERS:
        path = ROOT / "scripts/adhoc" / filename
        spec = importlib.util.spec_from_file_location(f"anchor2336_{path.stem}", path)
        assert spec is not None
        assert spec.loader is not None
        verifier = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = verifier
        spec.loader.exec_module(verifier)
        mutants = verifier.MUTANTS
        assert tuple(mutant.name for mutant in mutants) == EXPECTED_MUTANTS[filename], (
            f"{filename}: the mutant list changed; repair anchors instead of dropping mutants"
        )
        for mutant in mutants:
            name = mutant.name
            target = getattr(mutant, target_attribute)
            anchor = mutant.old
            assert isinstance(name, str)
            assert isinstance(target, Path)
            assert isinstance(anchor, str)
            assert anchor, f"{filename}:{name} must define a nonempty source anchor"
            cases.append((f"{path.stem}:{name}", ROOT / target, anchor))
    return cases


ANCHOR_CASES = _anchor_cases()


@pytest.mark.parametrize(
    ("case_name", "target_path", "anchor"),
    ANCHOR_CASES,
    ids=[case[0] for case in ANCHOR_CASES],
)
def test_real_mutation_anchor_matches_current_source_once(
    case_name: str, target_path: Path, anchor: str
) -> None:
    source = target_path.read_text(encoding="utf-8")
    count = source.count(anchor)
    assert count == 1, (
        f"{case_name}: expected exactly one anchor in {target_path.relative_to(ROOT)}, "
        f"found {count}"
    )
