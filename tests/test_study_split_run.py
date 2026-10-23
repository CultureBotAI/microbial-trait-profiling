"""Output isolation and replay boundary regressions for the split sidecar (#854)."""

from pathlib import Path
from types import SimpleNamespace

import pytest

from mtp.cellulose_decr import study_split_run as runner
from mtp.cellulose_decr.study_split import Protocol
from mtp.multiplicity import sources
from mtp.paths import file_sha256


def test_consumed_scanner_protocol_is_bound_in_provenance() -> None:
    assert "methods/cellulose_decr_selection_scan_protocol.json" in runner.IMPLEMENTATIONS


def test_active_scientific_helper_changes_provenance(monkeypatch: pytest.MonkeyPatch) -> None:
    before = runner._provenance({}, Protocol())
    original_hash = file_sha256

    def changed_helper(path: Path) -> str:
        if path.name == "table1_trait_maps.py":
            return "f" * 64
        return original_hash(path)

    monkeypatch.setattr(runner, "file_sha256", changed_helper)
    assert runner._provenance({}, Protocol()) != before


@pytest.mark.parametrize("location", ["child", "same", "ancestor", "symlink"])
def test_outputs_cannot_modify_an_upstream_bundle(tmp_path: Path, location: str) -> None:
    source = tmp_path / "verified" / "arms" / "noel"
    source.mkdir(parents=True)
    manifest = source / "arm_run.json"
    manifest.write_text("{}\n")
    outputs = {
        "child": source / "sidecar",
        "same": source,
        "ancestor": source.parent,
        "symlink": tmp_path / "alias",
    }
    outputs["symlink"].symlink_to(source, target_is_directory=True)
    with pytest.raises(ValueError, match="immutable upstream bundle"):
        runner.preflight_bundle_roots(outputs[location], {"arm:noel:manifest": manifest})
    assert manifest.read_text() == "{}\n"
    assert not (source / "sidecar").exists()


def test_sibling_sidecar_is_allowed_without_writing(tmp_path: Path) -> None:
    output = tmp_path / "sidecars" / "splits"
    runner.preflight_bundle_roots(
        output,
        {
            "comparison:manifest": tmp_path / "comparison" / "comparison_run.json",
            "core:manifest": tmp_path / "core" / "run.json",
            "arm:noel:manifest": tmp_path / "arms" / "noel" / "arm_run.json",
        },
    )
    assert not output.exists()


def test_directory_input_uses_the_verified_manifest_path(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    comparison_dir = tmp_path / "comparison"
    manifest = comparison_dir / "comparison_run.json"
    upstream = SimpleNamespace(
        core=SimpleNamespace(
            contrast_population="disease_vs_healthy_tier1_2_subject_controls",
            direction_basis="disease_vs_comparison",
        ),
        scope="off",
        comparison=SimpleNamespace(
            manifest_path=str(manifest),
            core_manifest_path=str(tmp_path / "core" / "run.json"),
            noel_manifest_path=str(tmp_path / "noel" / "arm_run.json"),
            jed_manifest_path=str(tmp_path / "jed" / "arm_run.json"),
            ilyas_manifest_path=str(tmp_path / "ilyas" / "arm_run.json"),
        ),
        product=lambda name: tmp_path / "core" / name,
    )
    monkeypatch.setattr(sources, "open_verified_sources", lambda _path: upstream)
    _, paths = runner._inputs(comparison_dir)
    assert paths["comparison:manifest"] == manifest
    runner.preflight_bundle_roots(tmp_path / "sidecar", paths)


@pytest.mark.parametrize(
    ("population", "direction", "scope"),
    [
        ("all_paired", "group_1_vs_group_0", "off"),
        ("disease_vs_healthy_tier1_2_subject_controls", "group_1_vs_group_0", "off"),
        ("disease_vs_healthy_tier1_2_subject_controls", "disease_vs_comparison", "full"),
    ],
)
def test_source_boundary_refuses_other_populations_directions_and_scopes(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    population: str,
    direction: str,
    scope: str,
) -> None:
    upstream = SimpleNamespace(
        core=SimpleNamespace(contrast_population=population, direction_basis=direction),
        scope=scope,
    )
    monkeypatch.setattr(sources, "open_verified_sources", lambda _path: upstream)
    with pytest.raises(ValueError, match=r"Tier 1\+2|phenotype sensitivity"):
        runner.execute(tmp_path / "comparison_run.json", tmp_path / "out", Protocol())
    assert not (tmp_path / "out").exists()
