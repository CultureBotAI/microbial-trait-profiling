"""Aerobe and anaerobe v2 condition filters and the paired-d minimum (#645)."""

import importlib.util
import math
import sys
from pathlib import Path
from types import ModuleType

import pandas as pd
import pytest

SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))


def load_script(name: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(
        f"guild_filters_{name}", SCRIPTS_DIR / f"{name}.py"
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    # Dataclasses in the scripts look their module up in sys.modules.
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


aerobe = load_script("investigate_facultative_aerobe_expansion")
anaerobe = load_script("investigate_gut_anaerobe_depletion")
guild_maps = load_script("create_condition_body_site_guild_maps")
workflow_diagram = load_script("create_profile_enrichment_guild_assignment_diagram")
oxygen_switches = load_script("investigate_taxa_oxygen_switches")


@pytest.mark.parametrize(
    "label",
    [
        "Time",
        "Timepoint",
        "Clinical treatment",
        "Health study participation",
        "Glucose",
        "Socioeconomic status",
        "Tumor stage",
    ],
)
def test_broad_covariate_labels_are_excluded_across_guild_scripts(label: str) -> None:
    assert aerobe.is_excluded_condition_label(label)
    assert anaerobe.is_excluded_condition_label(label)
    assert workflow_diagram.is_excluded_condition_label(label)
    for config in (guild_maps.FACULTATIVE_CONFIG, guild_maps.ANAEROBE_CONFIG):
        assert guild_maps.is_excluded_condition_label(label, config.excluded_condition_labels)


def test_excluded_condition_label_sets_are_aligned() -> None:
    expected = set(aerobe.EXCLUDED_CONDITION_SUMMARY_LABELS)
    assert set(anaerobe.EXCLUDED_CONDITION_SUMMARY_LABELS) == expected
    assert set(workflow_diagram.EXCLUDED_CONDITION_LABELS) == expected
    assert set(guild_maps.FACULTATIVE_CONFIG.excluded_condition_labels) == expected
    assert set(guild_maps.ANAEROBE_CONFIG.excluded_condition_labels) == expected


def test_disease_labels_stay_in_aerobe_condition_summaries() -> None:
    assert not aerobe.is_excluded_condition_label("Parkinson's disease")
    assert not anaerobe.is_excluded_condition_label("Parkinson's disease")
    excluded = guild_maps.FACULTATIVE_CONFIG.excluded_condition_labels
    assert not guild_maps.is_excluded_condition_label("Parkinson's disease", excluded)


def test_select_top_conditions_requires_at_least_three_comparator_profiles() -> None:
    aerobe_summary = pd.DataFrame(
        [
            {
                "condition": "Sparse comparator condition",
                "increased_profile_count": 5,
                "decreased_profile_count": 2,
                "study_count": 4,
                "increased_minus_decreased_guild_z": 1.2,
                "increased_mean_guild_z": 0.8,
            },
            {
                "condition": "Balanced comparator condition",
                "increased_profile_count": 5,
                "decreased_profile_count": 3,
                "study_count": 4,
                "increased_minus_decreased_guild_z": 0.9,
                "increased_mean_guild_z": 0.7,
            },
        ]
    )
    aerobe_top = aerobe.select_top_conditions(aerobe_summary, limit=10)
    assert list(aerobe_top["condition"]) == ["Balanced comparator condition"]

    anaerobe_summary = pd.DataFrame(
        [
            {
                "condition": "Acute graft vs. host disease",
                "decreased_profile_count": 6,
                "increased_profile_count": 2,
                "study_count": 3,
                "decreased_minus_increased_guild_z": 1.3,
                "decreased_mean_guild_z": 0.5,
            },
            {
                "condition": "Crohn's disease",
                "decreased_profile_count": 6,
                "increased_profile_count": 3,
                "study_count": 3,
                "decreased_minus_increased_guild_z": 0.9,
                "decreased_mean_guild_z": 0.4,
            },
        ]
    )
    anaerobe_top = anaerobe.select_top_conditions(anaerobe_summary, limit=10)
    assert list(anaerobe_top["condition"]) == ["Crohn's disease"]


@pytest.mark.parametrize("module", [aerobe, anaerobe], ids=["aerobe", "anaerobe"])
def test_target_taxon_stacked_bars_have_distinct_colors_and_hatches(
    module: ModuleType,
) -> None:
    colors = module.TARGET_TAXON_COLORS
    hatches = module.TARGET_TAXON_HATCHES
    assert set(colors) == set(hatches)
    assert len(set(colors.values())) == len(colors)
    assert len(set(hatches.values())) == len(hatches)


def paired_profiles(pair_count: int) -> pd.DataFrame:
    """`pair_count` experiments, each with one increased and one decreased profile."""
    rows = []
    for index in range(pair_count):
        experiment = f"s{index}-1"
        rows.append((f"s{index}", experiment, "inc", 1.0 + index * 0.1))
        rows.append((f"s{index}", experiment, "dec", -0.5 + index * 0.05))
    return pd.DataFrame(rows, columns=["study_id", "experiment_id", "inc_dec", "guild_z_mean"])


@pytest.mark.parametrize("module", [aerobe, anaerobe], ids=["aerobe", "anaerobe"])
def test_paired_d_needs_five_matched_pairs(module: ModuleType) -> None:
    four = module.paired_guild_enrichment_stats(paired_profiles(4))
    five = module.paired_guild_enrichment_stats(paired_profiles(5))

    assert four["paired_profile_pair_count"] == 4
    assert math.isnan(four["paired_standardized_mean_difference"])
    assert math.isnan(four["paired_t_test_p_value"])
    assert five["paired_profile_pair_count"] == 5
    assert math.isfinite(five["paired_standardized_mean_difference"])
    assert math.isfinite(five["paired_t_test_p_value"])


def test_oxygen_switch_summary_uses_side_specific_study_counts() -> None:
    taxon_condition = pd.DataFrame(
        [
            {
                "taxon_id": "NCBITaxon:1",
                "taxon_name": "Example taxon",
                "taxonomic_rank": "species",
                "oxygen_class": "anaerobic_only",
                "condition": "Alcoholic hepatitis",
                "dominant_direction": "increased",
                "profile_count": 10,
                "study_count": 7,
                "increased_profile_count": 6,
                "decreased_profile_count": 4,
                "increased_study_count": 4,
                "decreased_study_count": 4,
                "aerobic_traits": "",
                "anaerobic_traits": "anaerobic",
                "body_site_displays": "Gut / feces",
            },
            {
                "taxon_id": "NCBITaxon:1",
                "taxon_name": "Example taxon",
                "taxonomic_rank": "species",
                "oxygen_class": "anaerobic_only",
                "condition": "Single-study side condition",
                "dominant_direction": "increased",
                "profile_count": 6,
                "study_count": 3,
                "increased_profile_count": 4,
                "decreased_profile_count": 2,
                "increased_study_count": 1,
                "decreased_study_count": 2,
                "aerobic_traits": "",
                "anaerobic_traits": "anaerobic",
                "body_site_displays": "Gut / feces",
            },
            {
                "taxon_id": "NCBITaxon:1",
                "taxon_name": "Example taxon",
                "taxonomic_rank": "species",
                "oxygen_class": "anaerobic_only",
                "condition": "COVID-19",
                "dominant_direction": "decreased",
                "profile_count": 21,
                "study_count": 14,
                "increased_profile_count": 6,
                "decreased_profile_count": 15,
                "increased_study_count": 5,
                "decreased_study_count": 11,
                "aerobic_traits": "",
                "anaerobic_traits": "anaerobic",
                "body_site_displays": "Gut / feces",
            },
        ]
    )
    switches = oxygen_switches.summarize_taxon_switches(taxon_condition)
    assert len(switches) == 1
    row = switches.iloc[0]
    assert row["increased_dominant_condition_count"] == 1
    assert row["top_increased_conditions"] == "Alcoholic hepatitis (6 profiles, 4 studies)"
    assert row["top_decreased_conditions"] == "COVID-19 (15 profiles, 11 studies)"
