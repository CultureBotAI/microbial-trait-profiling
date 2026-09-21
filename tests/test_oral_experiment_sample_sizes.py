"""Guard sample-size denominators and membership in the saved oral analysis."""

import importlib
import sys
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest

with patch.object(sys, "path", [str(Path(__file__).resolve().parents[1] / "scripts"), *sys.path]):
    sizes = importlib.import_module("plot_oral_experiment_sample_sizes")


def test_dump_reader_preserves_url_fragments_and_multiline_fields(tmp_path: Path) -> None:
    path = tmp_path / "full_dump.csv"
    path.write_text(
        "# BugSigDB snapshot\nBSDB ID,URL,Description\n"
        'bsdb:1/1/1,https://example.org/#section,"first line\nsecond line"\n'
    )
    frame = sizes.read_raw_dump(path)
    assert len(frame) == 1
    assert frame.loc[0, "URL"] == "https://example.org/#section"
    assert frame.loc[0, "Description"] == "first line\nsecond line"


def test_oral_membership_deduplicates_directions_and_keeps_disease_only_rows() -> None:
    metadata = pd.DataFrame(
        {
            "study_id": ["bsdb:1"] * 4,
            "experiment_id": ["bsdb:1-1", "bsdb:1-1", "bsdb:1-2", "bsdb:1-3"],
            "inc_dec": ["inc", "dec", "inc", "inc"],
            "condition": ["Gingivitis", "Gingivitis", "Colorectal cancer", "Periodontitis"],
            "body_site": ["Saliva", "Saliva", "Feces", "Feces"],
            "body_site_group": ["Saliva", "Saliva", "Feces", "Feces"],
            "host_species": ["Homo sapiens"] * 4,
        }
    )
    cells = pd.DataFrame(
        {
            "view_key": ["individual_conditions_direction_contrast"] * 3,
            "condition": ["Gingivitis", "Colorectal cancer", "Periodontitis"],
            "body_site_group": ["Oral cavity / dental", "Feces", "Feces"],
            "trait_id": ["trait"] * 3,
            "displayed_in_heatmap": [True] * 3,
            "selected_after_bh": [False] * 3,
            "fisher_p_value": [1.0] * 3,
            "fisher_q_value": [1.0] * 3,
        }
    )
    cells["condition_body_site_pair"] = cells["condition"] + " | " + cells["body_site_group"]
    result = sizes.oral_experiments(metadata, cells)
    assert result["experiment_id"].tolist() == ["bsdb:1-1", "bsdb:1-3"]
    assert result["profile_directions"].tolist() == ["dec;inc", "inc"]
    with pytest.raises(ValueError, match="missing from profile metadata"):
        sizes.oral_experiments(metadata.iloc[:2], cells)


def test_signature_deduplication_missing_conflicts_and_reference_counting() -> None:
    experiments = pd.DataFrame(
        {"study_id": ["bsdb:1"] * 6, "experiment_id": [f"bsdb:1-{i}" for i in range(1, 7)]}
    )
    raw = pd.DataFrame(
        {
            "BSDB ID": [
                "bsdb:1/1/1",
                "bsdb:1/1/2",
                "bsdb:1/2/1",
                "bsdb:1/3/1",
                "bsdb:1/3/2",
                "bsdb:1/4/1",
                "bsdb:1/5/1",
            ],
            "Group 0 sample size": ["20", "20.0", "NA", "10", "12", "100", "-1"],
            "Group 1 sample size": ["30", "30", "80", "20", "20", "100", "20"],
            "Group 0 name": ["controls"] * 7,
            "Group 1 name": ["cases"] * 7,
            "Study design": ["case-control"] * 7,
        }
    )
    audit = sizes.sample_size_audit(experiments, raw).set_index("experiment_id")
    assert audit.loc["bsdb:1-1", "total_n"] == 50  # two signatures, one experiment
    assert audit.loc["bsdb:1-2", "size_status"] == "group_0:missing"
    assert audit.loc["bsdb:1-3", "size_status"] == "group_0:conflicting"
    assert audit.loc["bsdb:1-5", "size_status"] == "group_0:invalid"
    assert audit.loc["bsdb:1-6", "size_status"] == "group_0:unmatched;group_1:unmatched"
    assert audit["total_n"].isna().sum() == 4
    summary = sizes.summarize(audit)
    assert summary["experiments_with_sizes"] == 2
    assert summary["median_n"] == 125
    assert summary["below_200"] == 1
    assert summary["equal_200"] == 1  # external reference is not an observation
    assert summary["oral_studies"] == 1  # repeated experiments do not become extra studies


@pytest.mark.parametrize("value", ["0", "-3", "12.5", "inf", "100 people", "10-20"])
def test_ambiguous_or_invalid_counts_are_not_plotted(value: str) -> None:
    assert sizes.parse_size(value) == (None, "invalid")
