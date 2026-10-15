"""Composition callers must divide by non-MediaDive totals (#174, #514).

Each test adds a large MediaDive column to a synthetic profile matrix; the
fractions or z-scores of the real traits must not change.
"""

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import pandas as pd
import pytest

from mtp.composition_denominators import retained_totals

MEDIA = "mediadive.medium:42"
TRAITS = pd.Index(["METPO:1", "METPO:2", "METPO:3"])


def load_script_module(script_name: str) -> ModuleType:
    module_path = Path(__file__).resolve().parents[1] / "scripts" / script_name
    name = script_name.removesuffix(".py")
    spec = importlib.util.spec_from_file_location(name, module_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def counts(with_media: bool) -> pd.DataFrame:
    frame = pd.DataFrame(
        {"METPO:1": [5, 1, 0, 4], "METPO:2": [1, 3, 2, 0], "METPO:3": [0, 2, 6, 1]},
        index=["s1|e1|inc", "s1|e1|dec", "s2|e2|inc", "s2|e2|dec"],
    )
    if with_media:
        frame[MEDIA] = [40, 0, 7, 90]
    return frame


def test_retained_totals_ignores_mediadive() -> None:
    assert retained_totals(counts(True)).tolist() == [6, 6, 8, 5]


@pytest.mark.parametrize(
    ("script", "function"),
    [
        ("create_study_trait_profile_heatmap.py", "build_plot_matrix"),
        ("create_condition_trait_module_heatmap.py", "build_plot_matrix"),
        ("create_trait_profile_infographic.py", "build_zscore_matrix"),
    ],
)
def test_plot_values_do_not_depend_on_mediadive(script: str, function: str) -> None:
    build = getattr(load_script_module(script), function)
    pd.testing.assert_frame_equal(build(counts(True), TRAITS), build(counts(False), TRAITS))


def test_study_heatmap_metadata_total_excludes_mediadive(tmp_path: Path) -> None:
    heatmap = load_script_module("create_study_trait_profile_heatmap.py")
    frame = counts(True).reset_index(names="key")
    keys = frame.pop("key").str.split("|", expand=True)
    frame.insert(0, "inc_dec", keys[2])
    frame.insert(0, "experiment_id", keys[1])
    frame.insert(0, "study_id", keys[0])
    path = tmp_path / "profiles.tsv"
    frame.to_csv(path, sep="\t", index=False)

    _, metadata, _ = heatmap.load_profile_matrix(path)

    assert metadata["total_trait_signal"].tolist() == [6, 6, 8, 5]


def test_cached_totals_with_mediadive_are_rejected(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1] / "scripts"))
    sankeys = load_script_module("create_six_shared_trait_sankeys.py")
    frame = counts(True).reset_index(names="key")
    keys = frame.pop("key").str.split("|", expand=True)
    for position, column in enumerate(["study_id", "experiment_id", "inc_dec"]):
        frame.insert(position, column, keys[position])
    matrix_path = tmp_path / "profiles.tsv"
    frame.to_csv(matrix_path, sep="\t", index=False)
    metadata = frame[["study_id", "experiment_id", "inc_dec"]].copy()
    metadata.index = counts(True).index
    metadata["condition"] = "c"
    metadata["body_site_group"] = "Feces"
    metadata_path = tmp_path / "metadata.tsv"

    monkeypatch.setattr(sankeys, "PROFILE_MATRIX_PATH", matrix_path)
    monkeypatch.setattr(sankeys, "PROFILE_METADATA_PATH", metadata_path)

    metadata["total_trait_signal"] = counts(True).sum(axis=1)  # stale: includes MediaDive
    metadata.to_csv(metadata_path, sep="\t")
    with pytest.raises(ValueError, match="non-MediaDive denominator"):
        sankeys.load_profile_metadata()

    metadata["total_trait_signal"] = retained_totals(counts(True))
    metadata.to_csv(metadata_path, sep="\t")
    assert len(sankeys.load_profile_metadata()) == 4


def trait_labels(columns: pd.Index) -> pd.DataFrame:
    labels = pd.DataFrame(index=columns)
    for name in [
        "trait_label",
        "display_label",
        "column_label",
        "object_name",
        "object_category",
        "display_category",
    ]:
        labels[name] = [str(column) for column in columns]
    return labels


@pytest.mark.parametrize(
    ("script", "call"),
    [
        (
            "create_study_trait_profile_heatmap.py",
            lambda module, frame: module.select_traits(
                frame, trait_labels(frame.columns), 10, 1, 1.0, False
            )[1],
        ),
        (
            "create_condition_trait_module_heatmap.py",
            lambda module, frame: module.select_traits(
                frame, trait_labels(frame.columns), 10, 1, 1.0, False
            )[1],
        ),
        (
            "create_trait_profile_infographic.py",
            lambda module, frame: module.select_stratum_traits(
                frame,
                trait_labels(frame.columns),
                max_traits=10,
                min_profile_support=1,
                max_profile_support_fraction=1.0,
            ).set_index("trait_id", drop=False),
        ),
    ],
)
def test_trait_selection_does_not_depend_on_mediadive(
    script: str, call: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    module = load_script_module(script)
    if hasattr(module, "REQUIRED_PAIRED_TRAITS"):
        monkeypatch.setattr(module, "REQUIRED_PAIRED_TRAITS", {})
    assert callable(call)
    with_media = call(module, counts(True)).loc[TRAITS, "composition_variance"]
    without_media = call(module, counts(False)).loc[TRAITS, "composition_variance"]
    pd.testing.assert_series_equal(with_media, without_media)


def test_condition_module_totals_exclude_mediadive(tmp_path: Path) -> None:
    heatmap = load_script_module("create_condition_trait_module_heatmap.py")
    frame = counts(True).reset_index(names="key")
    keys = frame.pop("key").str.split("|", expand=True)
    for position, column in enumerate(["study_id", "experiment_id", "inc_dec"]):
        frame.insert(position, column, keys[position])
    path = tmp_path / "profiles.tsv"
    frame.to_csv(path, sep="\t", index=False)

    study_matrix, study_metadata = heatmap.load_aggregated_study_profiles(path)
    assert study_metadata["total_trait_signal"].to_dict() == {"s1": 12, "s2": 13}
    assert study_metadata["trait_count"].to_dict() == {"s1": 3, "s2": 3}

    for column, value in [
        ("condition", "c"),
        ("body_site", "Feces"),
        ("location", "x"),
        ("host_species", "Homo sapiens"),
    ]:
        study_metadata[column] = value
    _, groups = heatmap.build_condition_modules(
        study_matrix,
        study_metadata.assign(condition=study_metadata.index.map({"s1": "c1", "s2": "c2"})),
        ["condition"],
        1,
        1,
        0,
        0,
        True,
        10,
    )
    assert groups["total_trait_signal"].to_dict() == {"c1": 12, "c2": 13}


def test_guild_relative_signal_does_not_depend_on_mediadive(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    followup = load_script_module("analyze_enrichment_infographic_followup.py")
    monkeypatch.setattr(followup, "GUILD_TRAITS", {"guild": {"METPO:1": "trait one"}})
    monkeypatch.setattr(followup, "TARGET_CONDITIONS", ("c",))
    metadata = pd.DataFrame(
        {
            "condition": "c",
            "study_id": ["s1", "s1", "s2", "s2"],
            "experiment_id": ["e1", "e1", "e2", "e2"],
            "inc_dec": ["inc", "dec", "inc", "dec"],
        },
        index=counts(True).index,
    )
    pd.testing.assert_frame_equal(
        followup.summarize_guild_connections(counts(True), metadata),
        followup.summarize_guild_connections(counts(False), metadata),
    )


@pytest.mark.parametrize(
    "script",
    ["create_pd_cellulose_analysis.py", "create_pd_significant_condition_taxon_trait_network.py"],
)
def test_other_consumers_reject_mediadive_cached_totals(
    script: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1] / "scripts"))
    module = load_script_module(script)
    frame = counts(True).reset_index(names="key")
    keys = frame.pop("key").str.split("|", expand=True)
    for position, column in enumerate(["study_id", "experiment_id", "inc_dec"]):
        frame.insert(position, column, keys[position])
    matrix_path = tmp_path / "profiles.tsv"
    frame.to_csv(matrix_path, sep="\t", index=False)
    monkeypatch.setattr(module, "PROFILE_MATRIX_PATH", matrix_path)
    metadata = frame[["study_id", "experiment_id", "inc_dec"]].copy()
    metadata.index = counts(True).index
    metadata["total_trait_signal"] = counts(True).sum(axis=1)  # stale: includes MediaDive

    with pytest.raises(ValueError, match="non-MediaDive denominator"):
        if script == "create_pd_cellulose_analysis.py":
            metadata["condition"] = "c"
            metadata["body_site_group"] = "Feces"
            metadata_path = tmp_path / "metadata.tsv"
            metadata.to_csv(metadata_path, sep="\t")
            monkeypatch.setattr(module, "PROFILE_METADATA_PATH", metadata_path)
            module.load_profile_metadata()
        else:
            module.build_zscore_matrix(metadata, pd.DataFrame({"trait_id": list(TRAITS)}))
