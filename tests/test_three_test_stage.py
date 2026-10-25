"""Run identity, non-destructive output and replay guards for #1044."""

from __future__ import annotations

import json
import os
import sys
import tempfile
from collections.abc import Iterator
from dataclasses import replace
from pathlib import Path

import pandas as pd
import pytest

from mtp.datamodel.multiple_testing_pydantic import PostArmMultiplicityRun
from mtp.datamodel.workflow_pydantic import ContrastPopulation, DirectionBasis
from mtp.multiplicity.sources import VerifiedSources
from mtp.multiplicity.stage import run_post_arm, verify_post_arm
from mtp.paths import file_sha256
from mtp.three_test_table import stage
from mtp.three_test_table.carriers import carrier_table as carrier_reader
from multiplicity_fixture import fixture_sources, write_table


@pytest.fixture
def tmp_path() -> Iterator[Path]:
    """Keep positive-stage artifacts durable independently of pytest's temporary root (#1058)."""
    root = Path(__file__).resolve().parents[1] / ".tmp" / "three-test-stage-tests"
    stage.require_durable(root).mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="case-", dir=root) as temporary:
        path = Path(temporary)
        with stage.durable_scratch(path / "default-scratch"):
            yield path


def sources_for_stage(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> VerifiedSources:
    sources = fixture_sources(tmp_path / "inputs", monkeypatch)
    sources.core.git_commit = "a" * 40
    sources.comparison.id = "fixture-comparison"
    sources.comparison.git_commit = "b" * 40
    for name, arm in sources.arms.items():
        arm.id = f"fixture-{name}"
        arm.git_commit = "c" * 40
    monkeypatch.setattr(stage, "open_verified_sources", lambda _path: sources)
    return sources


def test_run_persists_complete_table_and_replays_without_machine_paths(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    comparison = Path(sources.comparison.manifest_path)
    output = tmp_path / "evidence"
    before = file_sha256(sources.product("count_profiles"))
    record = stage.run_table(comparison, output, scratch=tmp_path / "scratch")
    table = pd.read_csv(output / "per_trait_evidence.tsv", sep="\t")
    assert list(table.trait) == list(sources.features["count_profiles"])
    assert len(table) == 7
    assert table.ora_calibration_status.eq(
        "uncalibrated (methods/calibration.md excludes native ORA)"
    ).all()
    assert table.carrier_status.eq("not_supplied").all()
    assert record["core_run_id"] == "fixture-core"
    assert record["comparison_run_id"] == "fixture-comparison"
    assert record["source_run_ids"] == {
        "noel": "fixture-noel",
        "jed": "fixture-jed",
        "ilyas": "fixture-ilyas",
    }
    assert len(record["code_commit"]) == 40
    assert record["primary_ranking"] == "undecided"
    assert str(tmp_path) not in (output / "three_test_run.json").read_text()
    assert file_sha256(sources.product("count_profiles")) == before
    assert stage.verify_table(output, comparison, scratch=tmp_path / "scratch") == record


def test_existing_output_and_source_ancestor_are_rejected(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    comparison = Path(sources.comparison.manifest_path)
    existing = tmp_path / "existing"
    existing.mkdir()
    (existing / "keep.txt").write_text("keep")
    with pytest.raises(FileExistsError):
        stage.run_table(comparison, existing, scratch=tmp_path / "scratch")
    assert (existing / "keep.txt").read_text() == "keep"
    with pytest.raises(ValueError, match="disjoint"):
        stage.run_table(comparison, comparison.parent / "nested", scratch=tmp_path / "scratch")


@pytest.mark.parametrize("path", ["/tmp/mtp-table", "/private/tmp/mtp-table"])
def test_system_temporary_output_rejected_before_source_read(
    path: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(tempfile, "gettempdir", lambda: "/tmp")
    with pytest.raises(ValueError, match="durable"):
        stage.run_table(Path("absent.json"), Path(path), scratch=Path(path) / "scratch")


def test_falsely_requested_code_commit_is_rejected(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    with pytest.raises(ValueError, match="commit"):
        stage.run_table(
            Path(sources.comparison.manifest_path),
            tmp_path / "output",
            scratch=tmp_path / "scratch",
            code_commit="0" * 40,
        )
    assert not (tmp_path / "output" / "three_test_run.json").exists()


def test_verification_rejects_table_tampering_even_if_local_hash_is_updated(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    comparison = Path(sources.comparison.manifest_path)
    output = tmp_path / "output"
    stage.run_table(comparison, output, scratch=tmp_path / "scratch")
    table_path = output / "per_trait_evidence.tsv"
    table = pd.read_csv(table_path, sep="\t")
    table.loc[0, "study_mean_delta_z"] = 200.0
    table.to_csv(table_path, sep="\t", index=False)
    manifest = output / "three_test_run.json"
    record = json.loads(manifest.read_text())
    record["outputs"]["per_trait_evidence.tsv"]["sha256"] = file_sha256(table_path)
    record["outputs"]["per_trait_evidence.tsv"]["size_bytes"] = table_path.stat().st_size
    manifest.write_text(json.dumps(record))
    with pytest.raises(ValueError, match="replay"):
        stage.verify_table(output, comparison, scratch=tmp_path / "scratch")


def test_source_tampering_is_not_promoted_to_a_new_snapshot(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    table = sources.product("count_profiles")
    table.write_text(table.read_text().replace("\t2\t", "\t9\t", 1))
    with pytest.raises(ValueError, match=r"(recorded|hash|Source)"):
        stage.run_table(
            Path(sources.comparison.manifest_path),
            tmp_path / "output",
            scratch=tmp_path / "scratch",
        )
    assert not (tmp_path / "output" / "three_test_run.json").exists()


def test_validated_population_and_direction_strings_publish_successfully(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    # Assignment validation applies the same use_enum_values model configuration
    # as reading a manifest; model_construct in the base fixture bypasses it.
    sources.comparison.contrast_population = ContrastPopulation.all_paired
    sources.comparison.direction_basis = DirectionBasis.group_1_vs_group_0
    assert type(sources.comparison.contrast_population) is str
    assert type(sources.comparison.direction_basis) is str
    output = tmp_path / "output"
    record = stage.run_table(
        Path(sources.comparison.manifest_path), output, scratch=tmp_path / "scratch"
    )
    assert record["contrast_population"] == "all_paired"
    assert record["direction_basis"] == "group_1_vs_group_0"
    table = pd.read_csv(output / "per_trait_evidence.tsv", sep="\t")
    assert table.contrast_population.eq("all_paired").all()
    assert table.direction_basis.eq("group_1_vs_group_0").all()


@pytest.mark.parametrize("prefix", ["", "biolink:has_phenotype|"])
@pytest.mark.parametrize(
    ("identifier", "correct", "incorrect"),
    [
        ("METPO:1000603", "anaerobic", "obligately anaerobic"),
        ("METPO:1000607", "obligately anaerobic", "anaerobic"),
    ],
)
def test_oxygen_identity_pin_applies_to_bare_and_qualified_traits(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    prefix: str,
    identifier: str,
    correct: str,
    incorrect: str,
) -> None:
    sources = fixture_sources(tmp_path / "inputs", monkeypatch)
    feature = prefix + identifier
    selected = replace(sources, features={**sources.features, "count_profiles": (feature,)})
    metadata = sources.product("feature_metadata")
    write_table(metadata, [{"column_id": feature, "object_name": correct, "column_label": ""}])
    assert stage._labels(selected).to_dict("records") == [
        {"trait": feature, "trait_label": correct}
    ]
    write_table(metadata, [{"column_id": feature, "object_name": incorrect, "column_label": ""}])
    with pytest.raises(ValueError, match=identifier):
        stage._labels(selected)


@pytest.mark.parametrize(
    "location", ["core", "comparison", "noel", "jed", "ilyas", "output", "multiplicity", "label"]
)
def test_verify_rejects_overlapping_scratch_before_creating_any_directory(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, location: str
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    comparison = Path(sources.comparison.manifest_path)
    output = tmp_path / "output"
    stage.run_table(comparison, output, scratch=tmp_path / "initial-scratch")
    options: dict[str, Path] = {}
    if location == "core":
        protected = Path(sources.core.manifest_path).parent
    elif location == "comparison":
        protected = comparison.parent
    elif location in sources.arms:
        protected = Path(sources.arms[location].manifest_path).parent
    elif location == "output":
        protected = output
    else:
        protected = tmp_path / location
        protected.mkdir()
        options["multiplicity" if location == "multiplicity" else "label_support"] = protected
    before = {path.relative_to(protected): path.is_dir() for path in protected.rglob("*")}
    scratch = protected / "must-not-be-created"
    with pytest.raises(ValueError, match="disjoint"):
        stage.verify_table(output, comparison, scratch=scratch, **options)
    after = {path.relative_to(protected): path.is_dir() for path in protected.rglob("*")}
    assert after == before
    assert not scratch.exists()


@pytest.mark.parametrize("mutation", ["artifact", "manifest", "extra", "remove", "symlink"])
def test_optional_multiplicity_change_after_verification_is_rejected(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, mutation: str
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    comparison = Path(sources.comparison.manifest_path)
    sidecar = tmp_path / "multiplicity"
    original = run_post_arm(comparison, sidecar)
    target = sidecar / original.families[0].corrected.identity

    def verified_then_changed(path: Path) -> PostArmMultiplicityRun:
        record = verify_post_arm(path)
        if mutation == "artifact":
            target.write_bytes(target.read_bytes() + b"\n")
        elif mutation == "manifest":
            manifest = sidecar / "multiplicity_run.json"
            manifest.write_bytes(manifest.read_bytes() + b"\n")
        elif mutation == "extra":
            (sidecar / "unexpected.txt").write_text("unexpected")
        elif mutation == "remove":
            target.unlink()
        else:
            (sidecar / "unexpected-link").symlink_to(target)
        return record

    monkeypatch.setattr(stage, "verify_post_arm", verified_then_changed)
    output = tmp_path / "output"
    with pytest.raises(ValueError, match="Multiplicity"):
        stage.run_table(comparison, output, scratch=tmp_path / "scratch", multiplicity=sidecar)
    assert not (output / stage.MANIFEST).exists()


def test_verified_optional_multiplicity_is_bound_and_replays(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    comparison = Path(sources.comparison.manifest_path)
    sidecar = tmp_path / "multiplicity"
    multiplicity = run_post_arm(comparison, sidecar)
    output = tmp_path / "output"
    record = stage.run_table(comparison, output, scratch=tmp_path / "scratch", multiplicity=sidecar)
    manifest = sidecar / "multiplicity_run.json"
    identity = f"{manifest.name}@sha256:{file_sha256(manifest)}"
    assert record["inputs"][identity]["sha256"] == file_sha256(manifest)
    for family in multiplicity.families:
        for artifact in (family.declared, family.raw, family.corrected):
            identity = f"{Path(artifact.identity).name}@sha256:{artifact.sha256}"
            assert record["inputs"][identity]["sha256"] == artifact.sha256
    assert (
        stage.verify_table(output, comparison, scratch=tmp_path / "scratch", multiplicity=sidecar)
        == record
    )


def test_optional_inventory_change_during_adaptation_is_rejected_before_publication(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    comparison = Path(sources.comparison.manifest_path)
    sidecar = tmp_path / "multiplicity"
    run_post_arm(comparison, sidecar)

    def change_after_initial_binding(
        bound: VerifiedSources, audit: Path | None
    ) -> tuple[pd.DataFrame, tuple[Path, ...]]:
        result = carrier_reader(bound, audit)
        (sidecar / "late-extra.txt").write_text("unexpected")
        return result

    monkeypatch.setattr(stage, "carrier_table", change_after_initial_binding)
    output = tmp_path / "output"
    with pytest.raises(ValueError, match="Multiplicity input inventory"):
        stage.run_table(comparison, output, scratch=tmp_path / "scratch", multiplicity=sidecar)
    assert not (output / stage.MANIFEST).exists()


@pytest.mark.parametrize("source", ["process", "environment"])
def test_declared_temporary_roots_are_rejected_as_written_and_resolved(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, source: str
) -> None:
    process = tmp_path / "process-temporary"
    environment = tmp_path / "environment-temporary"
    monkeypatch.setattr(tempfile, "gettempdir", lambda: str(process))
    monkeypatch.setenv("TMPDIR", str(environment))
    temporary = process if source == "process" else environment
    for candidate in (temporary, temporary / "output"):
        with pytest.raises(ValueError, match="durable"):
            stage.require_durable(candidate)
    sibling = tmp_path / f"{temporary.name}-sibling" / "output"
    assert stage.require_durable(sibling) == sibling.resolve()


@pytest.mark.parametrize("direction", ["into_temporary", "out_of_temporary"])
def test_temporary_path_aliases_cannot_hide_written_or_resolved_roots(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, direction: str
) -> None:
    temporary = tmp_path / "declared-temporary"
    durable = tmp_path / "durable"
    temporary.mkdir()
    durable.mkdir()
    monkeypatch.setattr(tempfile, "gettempdir", lambda: str(temporary))
    monkeypatch.setenv("TMPDIR", str(temporary))
    if direction == "into_temporary":
        alias = durable / "link"
        alias.symlink_to(temporary, target_is_directory=True)
    else:
        alias = temporary / "link"
        alias.symlink_to(durable, target_is_directory=True)
    with pytest.raises(ValueError, match="durable"):
        stage.require_durable(alias / "output")


@pytest.mark.skipif(sys.platform != "linux", reason="Documents the Linux-only #1051 gap")
@pytest.mark.parametrize("path", ["/tmp/mtp-table", "/var/tmp/mtp-table"])
def test_linux_root_not_declared_by_process_or_environment_remains_a_known_gap(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, path: str
) -> None:
    custom = tmp_path / "custom-temporary"
    monkeypatch.setattr(tempfile, "gettempdir", lambda: str(custom))
    monkeypatch.setenv("TMPDIR", str(custom))
    # No directory is created: this states the accepted compatibility limitation.
    assert stage.require_durable(Path(path)) == Path(path).resolve()


def test_verify_nested_replay_preserves_inherited_temporary_environment(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    comparison = Path(sources.comparison.manifest_path)
    output = tmp_path / "output"
    original = stage.run_table(comparison, output, scratch=tmp_path / "scratch")
    inherited = os.environ["TMPDIR"]
    assert stage.verify_table(output, comparison, scratch=tmp_path / "scratch") == original
    assert os.environ["TMPDIR"] == inherited
    assert tempfile.gettempdir() == inherited
