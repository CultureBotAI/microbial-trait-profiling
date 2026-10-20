"""Preserve missing PR #526 caller acceptance coverage for issue #528.

Synthetic profiles exercise eligibility, trait selection and stale-cache rejection
without downloading data or rendering figures. Reviewed 2026-10-16.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType, SimpleNamespace

import pandas as pd
import pytest

MEDIA = "METPO:2000517|MediaDive.Medium:42"
TRAITS = pd.Index(["METPO:1", "METPO:2", "METPO:3"])
SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"


def load_script_module(name: str, monkeypatch: pytest.MonkeyPatch) -> ModuleType:
    monkeypatch.syspath_prepend(str(SCRIPTS_DIR))
    module_name = f"denominator_callers_{name}"
    spec = importlib.util.spec_from_file_location(module_name, SCRIPTS_DIR / f"{name}.py")
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    monkeypatch.setitem(sys.modules, module_name, module)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def counts() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "METPO:1": [5, 1, 0, 4, 0],
            "METPO:2": [1, 3, 2, 0, 0],
            "METPO:3": [0, 2, 6, 1, 0],
            MEDIA: [4000, 0, 7, 90, 10000],
        },
        index=["s1|e1|inc", "s1|e1|dec", "s2|e2|inc", "s2|e2|dec", "s3|e3|inc"],
    )


@pytest.fixture
def profile_path(tmp_path: Path, counts: pd.DataFrame) -> Path:
    frame = counts.reset_index(names="key")
    keys = frame.pop("key").str.split("|", expand=True)
    for position, column in enumerate(["study_id", "experiment_id", "inc_dec"]):
        frame.insert(position, column, keys[position])
    path = tmp_path / "profiles.tsv"
    frame.to_csv(path, sep="\t", index=False)
    return path


def trait_labels(columns: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        dict.fromkeys(
            [
                "trait_label",
                "display_label",
                "column_label",
                "object_name",
                "object_category",
                "display_category",
            ],
            columns,
        ),
        index=columns,
    )


@pytest.mark.parametrize(
    "script",
    [
        "create_study_trait_profile_heatmap",
        "create_condition_trait_module_heatmap",
        "create_trait_profile_infographic",
    ],
)
def test_trait_selection_ignores_media_signal_and_never_selects_media(
    script: str, counts: pd.DataFrame, monkeypatch: pytest.MonkeyPatch
) -> None:
    module = load_script_module(script, monkeypatch)
    if hasattr(module, "REQUIRED_PAIRED_TRAITS"):
        monkeypatch.setattr(module, "REQUIRED_PAIRED_TRAITS", {})
    selections = []
    for frame in [counts, counts.drop(columns=MEDIA)]:
        if script == "create_trait_profile_infographic":
            selection = module.select_stratum_traits(
                frame,
                trait_labels(frame.columns),
                max_traits=10,
                min_profile_support=1,
                max_profile_support_fraction=1.0,
            )
        else:
            _, selection = module.select_traits(
                frame, trait_labels(frame.columns), 10, 1, 1.0, False
            )
        assert set(selection["trait_id"]) == set(TRAITS)
        selections.append(selection)
    pd.testing.assert_frame_equal(selections[0], selections[1])
    if script != "create_trait_profile_infographic":
        with pytest.raises(ValueError, match="MediaDive"):
            module.select_traits(counts, trait_labels(counts.columns), 10, 1, 1.0, True)


def test_study_profile_eligibility_ignores_media_only_signal(
    profile_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    module = load_script_module("create_study_trait_profile_heatmap", monkeypatch)
    matrix, metadata, columns = module.load_profile_matrix(profile_path)

    assert columns == TRAITS.tolist()
    assert metadata["total_trait_signal"].tolist() == [6, 6, 8, 5, 0]
    assert metadata["trait_count"].tolist() == [2, 3, 2, 2, 0]
    selected = module.select_profiles(matrix, metadata, 0, 2, 5)
    assert selected.tolist() == ["s1|e1|inc", "s1|e1|dec", "s2|e2|inc", "s2|e2|dec"]
    assert module.select_profiles(matrix, metadata, 1, 2, 5).tolist() == ["s2|e2|inc"]


def test_condition_module_eligibility_and_totals_ignore_media_only_signal(
    profile_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    module = load_script_module("create_condition_trait_module_heatmap", monkeypatch)
    matrix, metadata = module.load_aggregated_study_profiles(profile_path)
    assert metadata["total_trait_signal"].to_dict() == {"s1": 12, "s2": 13, "s3": 0}
    assert metadata["trait_count"].to_dict() == {"s1": 3, "s2": 3, "s3": 0}
    assert matrix.columns.tolist() == TRAITS.tolist()
    metadata["condition"] = ["c1", "c2", "c3"]
    metadata["body_site"] = "Feces"
    metadata["location"] = "USA"
    metadata["host_species"] = "Homo sapiens"

    # A caller may supply its own matrix; group filtering must enforce the same rule.
    matrix[MEDIA] = [4000, 97, 10000]
    groups, group_metadata = module.build_condition_modules(
        matrix, metadata, ["condition"], 1, 1, 1, 1, True, 10
    )
    assert group_metadata["total_trait_signal"].to_dict() == {"c2": 13, "c1": 12}
    assert group_metadata["trait_count"].to_dict() == {"c2": 3, "c1": 3}
    assert groups.columns.tolist() == TRAITS.tolist()


@pytest.mark.parametrize(
    "script",
    [
        "create_six_shared_trait_sankeys",
        "create_pd_cellulose_analysis",
        "create_pd_significant_condition_taxon_trait_network",
    ],
)
def test_consumers_reject_media_inclusive_cached_totals(
    script: str,
    profile_path: Path,
    counts: pd.DataFrame,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    module = load_script_module(script, monkeypatch)
    monkeypatch.setattr(module, "PROFILE_MATRIX_PATH", profile_path)
    metadata = pd.DataFrame({"total_trait_signal": [4006, 6, 15, 95, 10000]}, index=counts.index)

    def read_values() -> pd.DataFrame:
        if script == "create_six_shared_trait_sankeys":
            traits = tuple(SimpleNamespace(trait_id=trait) for trait in TRAITS)
            values = module.load_raw_trait_counts(traits, metadata)
        elif script == "create_pd_cellulose_analysis":
            values = module.build_zscore_matrix(metadata, TRAITS.tolist())
        else:
            values = module.build_zscore_matrix(metadata, pd.DataFrame({"trait_id": TRAITS}))
        assert isinstance(values, pd.DataFrame)
        return values

    with pytest.raises(ValueError, match="Stale or incompatible non-MediaDive denominator cache"):
        read_values()
    metadata["total_trait_signal"] = [6, 6, 8, 5, 0]
    values = read_values()
    assert values.index.tolist() == counts.index.tolist()
    assert values.columns.tolist() == TRAITS.tolist()


def test_followup_profile_metadata_excludes_media_from_counts_and_totals(
    profile_path: Path, counts: pd.DataFrame, monkeypatch: pytest.MonkeyPatch
) -> None:
    module = load_script_module("analyze_enrichment_infographic_followup", monkeypatch)
    experiments = pd.DataFrame(
        {
            "study_id": ["s1", "s2", "s3"],
            "experiment_id": ["e1", "e2", "e3"],
            "condition": ["export label"] * 3,
        }
    )
    studies = pd.DataFrame({"study_id": ["s1", "s2", "s3"]})
    labels = pd.DataFrame(
        {"condition": ["c1", "c1", "c2", "c2", "c3"]},
        index=counts.index,
    ).reindex(columns=module.PROFILE_LABEL_COLUMNS)
    values, metadata, columns = module.load_profiles(profile_path, experiments, studies, labels)

    assert metadata["condition"].tolist() == ["c1", "c1", "c2", "c2", "c3"]
    assert columns == TRAITS.tolist()
    assert metadata.index.tolist() == counts.index.tolist()
    assert metadata["total_trait_signal"].tolist() == [6, 6, 8, 5, 0]
    assert metadata["profile_trait_count"].tolist() == [2, 3, 2, 2, 0]
    pd.testing.assert_frame_equal(values, counts.drop(columns=MEDIA).astype("float32"))
