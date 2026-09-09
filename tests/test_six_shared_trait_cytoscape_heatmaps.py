import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace

import pandas as pd

MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "create_six_shared_trait_cytoscape_heatmaps.py"
)
SPEC = importlib.util.spec_from_file_location(
    "create_six_shared_trait_cytoscape_heatmaps",
    MODULE_PATH,
)
assert SPEC is not None
assert SPEC.loader is not None
heatmaps = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = heatmaps
SPEC.loader.exec_module(heatmaps)


def test_prepare_pattern_cells_masks_by_fisher_p_value_and_support() -> None:
    view = SimpleNamespace(
        view_key="increased_grouped",
        view_label="Increased, grouped conditions",
        focus_direction="increased",
        condition_mode="grouped condition labels",
    )
    traits = pd.DataFrame(
        {
            "trait_id": ["trait1"],
            "trait_label": ["trait one"],
            "trait_rank": [0],
        }
    )
    edges = pd.DataFrame(
        [
            {
                "condition": "Condition A",
                "body_site_group": "Feces",
                "trait_id": "trait1",
                "trait_label": "trait one",
                "focus_direction": "increased",
                "study_count": 5,
                "focus_profile_count": 20,
                "focus_enriched_profile_count": 10,
                "focus_enriched_study_count": 4,
                "focus_no_enriched_trait_profile_count": 10,
                "comparator_profile_count": 20,
                "comparator_enriched_profile_count": 1,
                "focus_enriched_mean_guild_z": 1.2,
                "focus_enriched_minus_comparator_guild_z": 1.1,
                "rank_score": 1.0,
            },
            {
                "condition": "Condition B",
                "body_site_group": "Feces",
                "trait_id": "trait1",
                "trait_label": "trait one",
                "focus_direction": "increased",
                "study_count": 1,
                "focus_profile_count": 20,
                "focus_enriched_profile_count": 10,
                "focus_enriched_study_count": 1,
                "focus_no_enriched_trait_profile_count": 10,
                "comparator_profile_count": 20,
                "comparator_enriched_profile_count": 1,
                "focus_enriched_mean_guild_z": 1.2,
                "focus_enriched_minus_comparator_guild_z": 1.1,
                "rank_score": 1.0,
            },
        ]
    )

    cells = heatmaps.prepare_pattern_cells(edges, view, traits)

    assert bool(cells.loc[0, "displayed_in_heatmap"]) is True
    assert cells.loc[0, "fisher_p_value"] <= heatmaps.P_VALUE_CUTOFF
    assert bool(cells.loc[1, "displayed_in_heatmap"]) is False
    assert cells.loc[1, "blank_reason"] == "low positive-z profile or study support"


def test_prepare_contrast_cells_uses_signed_increased_minus_decreased_fraction() -> None:
    mode = heatmaps.ContrastMode(
        key="grouped",
        label="Grouped conditions",
        condition_mode="grouped condition labels",
        increased_view_key="increased_grouped",
        decreased_view_key="decreased_grouped",
    )
    traits = pd.DataFrame(
        {
            "trait_id": ["trait1"],
            "trait_label": ["trait one"],
            "trait_rank": [0],
        }
    )
    increased_edges = pd.DataFrame(
        [
            {
                "condition": "Condition A",
                "body_site_group": "Feces",
                "trait_id": "trait1",
                "trait_label": "trait one",
                "focus_direction": "increased",
                "study_count": 5,
                "focus_profile_count": 20,
                "focus_enriched_profile_count": 10,
                "comparator_profile_count": 20,
                "comparator_enriched_profile_count": 1,
            }
        ]
    )
    decreased_edges = pd.DataFrame(
        [
            {
                "condition": "Condition A",
                "body_site_group": "Feces",
                "trait_id": "trait1",
                "trait_label": "trait one",
                "focus_direction": "decreased",
                "study_count": 5,
                "focus_profile_count": 20,
                "focus_enriched_profile_count": 1,
                "comparator_profile_count": 20,
                "comparator_enriched_profile_count": 10,
            }
        ]
    )

    cells = heatmaps.prepare_contrast_cells(increased_edges, decreased_edges, mode, traits)

    assert bool(cells.loc[0, "displayed_in_heatmap"]) is True
    assert cells.loc[0, "increased_minus_decreased_positive_fraction"] == 0.45
    assert cells.loc[0, "fisher_p_value"] <= heatmaps.P_VALUE_CUTOFF


def test_select_pattern_rows_prioritizes_significant_trait_recurrence() -> None:
    cells = pd.DataFrame(
        [
            {
                "condition_body_site_pair": "A | Feces",
                "displayed_in_heatmap": True,
                "focus_positive_profile_count": 3,
                "strength": 1.0,
            },
            {
                "condition_body_site_pair": "A | Feces",
                "displayed_in_heatmap": True,
                "focus_positive_profile_count": 3,
                "strength": 1.0,
            },
            {
                "condition_body_site_pair": "B | Feces",
                "displayed_in_heatmap": True,
                "focus_positive_profile_count": 20,
                "strength": 20.0,
            },
        ]
    )

    assert heatmaps.select_pattern_rows(cells, max_rows=1) == ["A | Feces"]
