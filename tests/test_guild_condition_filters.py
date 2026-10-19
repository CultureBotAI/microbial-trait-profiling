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


@pytest.mark.parametrize("label", ["Time", "Timepoint"])
def test_time_labels_are_excluded_from_aerobe_condition_summaries(label: str) -> None:
    assert aerobe.is_excluded_condition_label(label)
    excluded = guild_maps.FACULTATIVE_CONFIG.excluded_condition_labels
    assert guild_maps.is_excluded_condition_label(label, excluded)


def test_disease_labels_stay_in_aerobe_condition_summaries() -> None:
    assert not aerobe.is_excluded_condition_label("Parkinson's disease")
    excluded = guild_maps.FACULTATIVE_CONFIG.excluded_condition_labels
    assert not guild_maps.is_excluded_condition_label("Parkinson's disease", excluded)


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
