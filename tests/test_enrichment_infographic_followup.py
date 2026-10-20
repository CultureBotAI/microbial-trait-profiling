"""The enrichment follow-up names target conditions it cannot summarize (#586)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest

MODULE_PATH = (
    Path(__file__).resolve().parents[1] / "scripts" / "analyze_enrichment_infographic_followup.py"
)
SPEC = importlib.util.spec_from_file_location("enrichment_infographic_followup", MODULE_PATH)
assert SPEC is not None
assert SPEC.loader is not None
followup = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = followup
SPEC.loader.exec_module(followup)


def test_targets_below_the_profile_minimum_are_reported_with_their_counts() -> None:
    present, rare, absent = followup.TARGET_CONDITIONS
    metadata = pd.DataFrame(
        {
            "profile_id": ["p1", "p2", "p3"],
            "condition": [present, present, rare],
            "study_id": ["s1", "s1", "s2"],
            "experiment_id": ["e1", "e2", "e3"],
            "profile_trait_count": [5, 3, 4],
            "total_trait_signal": [10, 6, 8],
        }
    )

    depth = followup.summarize_condition_depth(metadata, min_profiles=2)

    assert depth["condition"].tolist() == [present]
    assert followup.dropped_target_conditions(metadata, depth) == {rare: 1, absent: 0}


def _profiles(rows: list[tuple[str, str, str, str]]) -> pd.DataFrame:
    """Profiles from (profile_id, study_id, experiment_id, condition) rows."""
    frame = pd.DataFrame(rows, columns=["profile_id", "study_id", "experiment_id", "condition"])
    frame["inc_dec"] = frame["profile_id"].str.rsplit("|", n=1).str[-1]
    frame["direction_label"] = frame["inc_dec"].map({"inc": "increased", "dec": "decreased"})
    frame["body_site"] = "Feces"
    frame["location"] = "Japan"
    frame["host_species"] = "Homo sapiens"
    frame["title"] = "A study"
    frame["experiment_label"] = "Experiment 1: control vs case"
    frame["profile_trait_count"] = 3
    frame["total_trait_signal"] = 6.0
    return frame.set_index("profile_id", drop=False)


def test_title_table_lists_corpus_experiments_and_marks_missing_titles(tmp_path: Path) -> None:
    target = followup.TARGET_CONDITIONS[0]
    metadata = _profiles(
        [
            ("s1|s1-1|inc", "s1", "s1-1", target),
            ("s1|s1-1|dec", "s1", "s1-1", target),
            ("s2|s2-1|inc", "s2", "s2-1", target),
        ]
    )
    experiments = pd.DataFrame(
        {
            "study_id": ["s1", "s2", "s3"],
            "experiment_id": ["s1-1", "s2-1", "s3-1"],
            "condition": [target] * 3,
            "efo_id": ["MONDO:1"] * 3,
            "experiment_label": ["Experiment 1: control vs case"] * 3,
            "group_1_definition": ["cases"] * 3,
            "sequencing_type": ["16S"] * 3,
        }
    )
    studies = pd.DataFrame(
        {"study_id": ["s1"], "title": ["Study one"], "publication_year": ["2020"]}
    ).assign(doi="", url="")

    table = followup.write_target_title_table(
        metadata, experiments, studies, tmp_path / "titles.tsv"
    ).set_index("experiment_id")

    # s3-1 carries the label in the dump but has no profile in this corpus.
    assert table.index.tolist() == ["s1-1", "s2-1"]
    assert table.loc["s1-1", "profile_count"] == 2
    assert table.loc["s1-1", "title"] == "Study one"
    assert table.loc["s2-1", "title"] == followup.MISSING_TITLE
    assert table.loc["s2-1", "publication_year"] == "unknown"


def test_cross_target_check_leaves_out_targets_it_never_queried(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    queried, unqueried = "Queried condition", "Rare condition"
    monkeypatch.setattr(followup, "TARGET_CONDITIONS", (queried, unqueried))
    metadata = _profiles(
        [
            ("s1|s1-1|inc", "s1", "s1-1", queried),
            ("s2|s2-1|inc", "s2", "s2-1", queried),
            ("s3|s3-1|inc", "s3", "s3-1", unqueried),
            ("s4|s4-1|inc", "s4", "s4-1", "Other condition"),
        ]
    )
    values = pd.DataFrame(
        [[3, 1, 0], [2, 2, 0], [0, 1, 4], [1, 0, 2]],
        index=metadata.index,
        columns=["METPO:1", "METPO:2", "METPO:3"],
    )

    neighbors, crosscheck, _ = followup.compute_profile_neighbors(
        values, metadata, top_n=5, min_condition_profiles=2
    )

    assert set(neighbors["query_condition"]) == {queried}
    assert crosscheck["query_condition"].tolist() == [queried]
