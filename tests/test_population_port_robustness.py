"""Output-preserving robustness from Jed's #1219, reviewed in #1354."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType
from typing import Any
from unittest.mock import Mock

import matplotlib
import pandas as pd
import pytest

from mtp import native_composition

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"


def load_script(name: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("name", ["create_clustered_heatmap", "create_condition_taxon_heatmap"])
def test_standalone_taxon_heatmaps_explicitly_select_headless_backend(
    name: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    # Observe the requested backend without disturbing caller figures/backend state.
    choose_backend = Mock()
    monkeypatch.setattr(matplotlib, "use", choose_backend)
    load_script(name)
    choose_backend.assert_called_once_with("Agg")


@pytest.mark.parametrize(
    ("name", "builder", "renderer", "metadata_name"),
    [
        ("create_clustered_heatmap", "build_binary_matrix", "save_clustered_heatmap", "inc_dec"),
        (
            "create_condition_taxon_heatmap",
            "build_condition_genus_matrix",
            "save_condition_heatmap",
            "study_count",
        ),
    ],
)
def test_matrix_export_creates_only_its_parent_and_preserves_table_bytes(
    name: str,
    builder: str,
    renderer: str,
    metadata_name: str,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    module = load_script(name)
    matrix = pd.DataFrame({"Genus A": [1, 0], "Genus B": [0, 1]}, index=["a", "b"])
    metadata = pd.Series(
        ["inc", "dec"] if metadata_name == "inc_dec" else [2, 3], index=matrix.index
    )
    before = matrix.copy(deep=True)
    monkeypatch.setattr(module, builder, Mock(return_value=(matrix, metadata)))
    render = Mock()
    monkeypatch.setattr(module, renderer, render)
    existing = tmp_path / "existing.tsv"
    fresh = tmp_path / "new" / "nested" / "matrix.tsv"
    image = tmp_path / "separate-images" / "heatmap.png"
    for destination in (existing, fresh):
        monkeypatch.setattr(
            sys,
            "argv",
            [name, "--matrix-output", str(destination), "--output", str(image)],
        )
        module.main()
    expected = before.copy()
    expected.insert(0, metadata_name, metadata)
    assert fresh.read_bytes() == existing.read_bytes() == expected.to_csv(sep="\t").encode()
    pd.testing.assert_frame_equal(matrix, before)
    assert render.call_count == 2
    assert not image.parent.exists()


def labels_for(columns: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        dict.fromkeys(
            (
                "trait_label",
                "display_label",
                "column_label",
                "object_name",
                "object_category",
                "display_category",
            ),
            columns,
        ),
        index=columns,
    )


def test_native_selection_owns_joined_metadata_without_changing_output_bytes(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    counts = pd.DataFrame({"METPO:A": [1, 2, 0, 0], "METPO:B": [0, 0, 2, 1]})
    labels = labels_for(counts.columns)
    expected = native_composition.select_stratum_traits(counts, labels, min_profile_support=1)
    joins: list[tuple[pd.DataFrame, pd.DataFrame]] = []
    original_join = pd.DataFrame.join

    def observe_join(frame: pd.DataFrame, *args: Any, **kwargs: Any) -> pd.DataFrame:
        joined = original_join(frame, *args, **kwargs)
        joins.append((joined, joined.copy(deep=True)))
        return joined

    monkeypatch.setattr(pd.DataFrame, "join", observe_join)
    observed = native_composition.select_stratum_traits(counts, labels, min_profile_support=1)
    assert observed.to_csv(sep="\t").encode() == expected.to_csv(sep="\t").encode()
    assert len(joins) == 1
    pd.testing.assert_frame_equal(*joins[0])


def test_missing_optional_powerpoint_does_not_fail_regeneration(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    module = load_script("regenerate_trait_profile_investigation")
    monkeypatch.setattr(module, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(module, "PPTX_PATH", tmp_path / "absent.pptx")
    module.refresh_powerpoint()
    assert "skipping PowerPoint refresh" in capsys.readouterr().out
    assert not list(tmp_path.iterdir())


def test_present_powerpoint_still_requires_its_source_images(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    module = load_script("regenerate_trait_profile_investigation")
    deck = tmp_path / "existing.pptx"
    deck.write_bytes(b"existing deck must not change")
    monkeypatch.setattr(module, "PPTX_PATH", deck)
    monkeypatch.setattr(module, "SLIDE_IMAGE_SOURCES", {"image": tmp_path / "absent.png"})
    with pytest.raises(FileNotFoundError, match="Missing slide image source"):
        module.refresh_powerpoint()
    assert deck.read_bytes() == b"existing deck must not change"


@pytest.mark.parametrize("metadata_present", [False, True])
def test_investigation_prepares_dependencies_without_changing_existing_inputs(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    metadata_present: bool,
) -> None:
    module = load_script("regenerate_trait_profile_investigation")
    metadata = tmp_path / "trait-visualization" / "study_trait_profile_study_metadata.tsv"
    if metadata_present:
        metadata.parent.mkdir()
        metadata.write_bytes(b"existing recorded population and directions\n")
    monkeypatch.setattr(module, "REPO_ROOT", tmp_path)
    calls: list[tuple[str, ...]] = []
    monkeypatch.setattr(module, "run_script", lambda *args: calls.append(args))
    monkeypatch.setattr(module, "refresh_powerpoint", Mock())
    module.main()
    preparation = (
        "create_study_trait_profile_heatmap.py",
        "--direction-basis",
        "group_1_vs_group_0",
    )
    if metadata_present:
        assert preparation not in calls
        assert metadata.read_bytes() == b"existing recorded population and directions\n"
    else:
        assert calls[0] == preparation
    additional = [
        ("create_six_shared_trait_sankeys.py",),
        ("create_six_shared_trait_sankeys.py", "--combined-high-support-only"),
        ("create_six_shared_trait_sankeys.py", "--high-support-layouts-only"),
        ("investigate_taxa_oxygen_switches.py",),
    ]
    assert [call for call in calls if call in additional] == additional
    assert calls.index(additional[0]) > calls.index(("create_trait_profile_infographic.py",))
    assert calls.index(additional[-1]) < calls.index(
        ("analyze_parkinsons_disease_trait_profiles.py",)
    )
    assert not any("--contrast-population" in call for call in calls)
    assert not any(call[0] == "create_taxa_coverage_histogram.py" for call in calls)
    prior_order = [
        "investigate_facultative_aerobe_expansion.py",
        "investigate_facultative_aerobe_condition_followup.py",
        "investigate_gut_anaerobe_depletion.py",
        "create_condition_body_site_guild_maps.py",
        "create_all_bugsigdb_condition_body_site_maps.py",
        "create_profile_enrichment_guild_assignment_diagram.py",
        "create_trait_profile_infographic.py",
        "analyze_parkinsons_disease_trait_profiles.py",
    ]
    assert [call[0] for call in calls if call[0] in prior_order] == prior_order


def test_anaerobic_cap_exemption_would_change_default_selected_table(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Characterize the deferred scientific change; do not adopt its new default."""
    module = load_script("create_study_trait_profile_heatmap")
    counts = pd.DataFrame(
        {
            "METPO:1000615": [10, 1] * 10,
            "METPO:1000698": [2, 8] * 10,
            "METPO:1000703": [5, 3] * 10,
            "METPO:1000603": [9, 1] * 10,
            "GO:0030245": [4, 0] * 10,
        }
    )
    labels = labels_for(counts.columns)
    options = {
        "max_traits": 180,
        "min_profile_support": 6,
        "max_profile_support_fraction": 0.95,
        "include_media_traits": False,
    }
    old_ids, old_table = module.select_traits(counts, labels, **options)
    assert "METPO:1000603" not in old_ids
    # In this fixture all eligible traits fit. Extending the existing exemption
    # reproduces the source PR's added anaerobic eligibility without changing code.
    monkeypatch.setattr(
        module,
        "REQUIRED_PAIRED_TRAITS",
        {
            **module.REQUIRED_PAIRED_TRAITS,
            "METPO:1000603": "anaerobic",
        },
    )
    new_ids, new_table = module.select_traits(counts, labels, **options)
    assert set(new_ids) - set(old_ids) == {"METPO:1000603"}
    assert new_table.loc["METPO:1000603", "profile_support"] == 20
    assert new_table.loc["METPO:1000603", "composition_variance"] > 0
    assert new_table.to_csv(sep="\t").encode() != old_table.to_csv(sep="\t").encode()
