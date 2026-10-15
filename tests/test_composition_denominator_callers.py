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
