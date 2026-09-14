import importlib.util
from pathlib import Path
from types import ModuleType
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

from bsdbng.trait_composition import retained_trait_counts, trait_fractions


def script(name: str) -> ModuleType:
    path = Path(__file__).resolve().parents[1] / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


composition = script("create_trait_profile_infographic")
heatmap = script("create_study_trait_profile_heatmap")


def corpus() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    ids = [f"s{s}|e|{direction}" for s in range(4) for direction in ["inc", "dec"]]
    counts = pd.DataFrame(
        {
            "A": [10, 2, 0, 0, 3, 5, 8, 1],
            "B": [90, 4, 0, 1, 5, 3, 1, 8],
            "C": [0, 3, 0, 2, 1, 0, 3, 0],
        },
        index=ids,
    )
    labels = pd.DataFrame(
        {
            name: list("ABC")
            for name in [
                "trait_label",
                "display_label",
                "column_label",
                "object_name",
                "object_category",
                "display_category",
            ]
        },
        index=list("ABC"),
    )
    pairs = pd.DataFrame(
        {
            "study_id": [f"s{s}" for s in range(4)],
            "increased_profile_ids": ids[::2],
            "decreased_profile_ids": ids[1::2],
            "body_site_group": ["Feces", "Feces", "Oral", "Oral"],
        }
    )
    return counts, labels, pairs


@pytest.mark.parametrize("scale", [1, 1000])
def test_media_changes_do_not_change_selection_zscores_or_paired_tests(
    scale: int, tmp_path: Path
) -> None:
    counts, labels, pairs = corpus()
    with_media = counts.assign(
        **{
            "METPO:2000517|mediadive.medium:1": np.array([900, 0, 70, 5, 99, 2, 40, 800]) * scale,
            "METPO:2000517|MediaDive.Medium:2": np.arange(8) * scale,
        }
    )
    pd.testing.assert_frame_equal(trait_fractions(counts), trait_fractions(with_media))
    expected = composition.select_stratum_traits(
        counts,
        labels,
        min_profile_support=1,
        max_profile_support_fraction=1.0,
    )
    actual = composition.select_stratum_traits(
        with_media,
        labels,
        min_profile_support=1,
        max_profile_support_fraction=1.0,
    )
    pd.testing.assert_frame_equal(expected, actual)
    pd.testing.assert_frame_equal(
        composition.build_zscore_matrix(counts, expected.index),
        composition.build_zscore_matrix(with_media, actual.index),
    )
    with patch.object(heatmap, "REQUIRED_PAIRED_TRAITS", {}):
        old_selected, old_metadata = heatmap.select_traits(counts, labels, 2, 1, 1.0, False)
        new_selected, new_metadata = heatmap.select_traits(with_media, labels, 2, 1, 1.0, False)
    pd.testing.assert_index_equal(old_selected, new_selected)
    pd.testing.assert_frame_equal(old_metadata, new_metadata)
    pd.testing.assert_frame_equal(
        heatmap.build_plot_matrix(counts, old_selected),
        heatmap.build_plot_matrix(with_media, new_selected),
    )
    for data, output in [(counts, tmp_path / "without"), (with_media, tmp_path / "with")]:
        composition.write_stratified_trait_stats(
            data,
            pairs,
            labels,
            output_dir=output,
            min_pairs=2,
            min_profile_support=1,
            max_profile_support_fraction=1.0,
        )
    for token in ["feces", "feces_excluded", "oral"]:
        filename = f"trait_profile_paired_contrast_trait_stats_{token}.tsv"
        pd.testing.assert_frame_equal(
            pd.read_csv(tmp_path / "without" / filename, sep="\t"),
            pd.read_csv(tmp_path / "with" / filename, sep="\t"),
        )
    kwargs = {
        "stratum_kind": "pooled",
        "stratum": "pooled",
        "min_profile_support": 1,
        "max_profile_support_fraction": 1.0,
    }
    pd.testing.assert_frame_equal(
        composition.compute_stratum_trait_stats(counts, pairs, labels, **kwargs),
        composition.compute_stratum_trait_stats(with_media, pairs, labels, **kwargs),
    )


def test_denominator_includes_unselected_traits_and_retains_zero_signal_profiles(
    tmp_path: Path,
) -> None:
    frame = pd.DataFrame(
        {"A": [10, 0, 0, 20], "B": [90, 0, 0, 80], "mediadive.medium:1": [900, 0, 700, 5]}
    )
    fractions = trait_fractions(frame)
    assert fractions["A"].tolist() == [0.1, 0.0, 0.0, 0.2]
    assert len(fractions) == 4
    assert fractions.iloc[1].eq(0).all()
    assert fractions.iloc[2].eq(0).all()
    expected_log = np.log1p(np.array([0.1, 0.0, 0.0, 0.2]) * 1000)
    expected_z = (expected_log - expected_log.mean()) / expected_log.std(ddof=0)
    actual = composition.build_zscore_matrix(frame, pd.Index(["A"]))
    np.testing.assert_allclose(actual["A"], expected_z)
    assert actual.iloc[1, 0] == actual.iloc[2, 0] < 0
    assert np.isfinite(actual.to_numpy()).all()

    path = tmp_path / "profiles.tsv"
    frame.assign(study_id=["s1", "s2", "s3", "s4"], experiment_id="e", inc_dec="inc").to_csv(
        path,
        sep="\t",
        index=False,
    )
    loaded, metadata, columns = heatmap.load_profile_matrix(path)
    assert columns == ["A", "B"]
    assert metadata.total_trait_signal.tolist() == [100, 0, 0, 100]
    assert metadata.trait_count.tolist() == [2, 0, 0, 2]
    assert len(loaded) == 4


def test_all_media_matrix_and_constant_traits_are_finite() -> None:
    only_media = pd.DataFrame({"mediadive.medium:1": [10, 0]})
    assert retained_trait_counts(only_media).shape == (2, 0)
    assert trait_fractions(only_media).shape == (2, 0)
    empty_trait = only_media.assign(A=0)
    assert composition.build_zscore_matrix(empty_trait, pd.Index(["A"])).eq(0).all().all()


def test_legacy_media_opt_in_is_rejected() -> None:
    counts, labels, _ = corpus()
    with pytest.raises(ValueError, match="MediaDive traits cannot be included"):
        heatmap.select_traits(counts, labels, 2, 1, 1.0, True)
