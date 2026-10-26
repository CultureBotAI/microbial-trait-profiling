"""Fail-closed publication tests for the two-parent comparison (#1041)."""

from __future__ import annotations

import json
import shutil
import tempfile
from collections.abc import Iterator
from pathlib import Path
from types import SimpleNamespace
from typing import Any
from unittest.mock import Mock

import pandas as pd
import pytest

from mtp.paths import file_sha256
from mtp.rollup_comparison import stage
from mtp.three_test_table.stage import _inventory, _publish_directory


@pytest.fixture
def tmp_path() -> Iterator[Path]:
    root = Path(__file__).resolve().parents[1] / ".tmp" / "rollup-stage-tests"
    root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="case-", dir=root) as directory:
        yield Path(directory)


def parents(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Any:
    bundles = []
    comparisons = []
    for scope in ("off", "phenotype"):
        root = tmp_path / scope
        root.mkdir()
        source = tmp_path / f"{scope}-source"
        source.mkdir()
        records = {
            name: source / f"{name}.json" for name in ("core", "comparison", "noel", "jed", "ilyas")
        }
        for path in records.values():
            path.write_text("{}\n", encoding="utf-8")
        records["comparison"].write_text(
            json.dumps(
                {
                    "manifest_path": str(records["comparison"]),
                    **{
                        f"{name}_manifest_path": str(path)
                        for name, path in records.items()
                        if name != "comparison"
                    },
                }
            ),
            encoding="utf-8",
        )
        evidence = root / "per_trait_evidence.tsv"
        evidence.write_text("trait\tstudy_t_q_by\ntrait:a\t0.012000\n", encoding="utf-8")
        record = {
            "run_id": f"parent-{scope}",
            "code_commit": "a" * 40,
            "code_dirty": False,
            "implementation_sha256": "b" * 64,
            "environment": {"python": "fixture"},
            "strain_rollup": scope,
            "scope_role": "primary" if scope == "phenotype" else "sensitivity",
            "contrast_population": "disease_vs_healthy_tier1_2_subject_controls",
            "direction_basis": "disease_vs_comparison",
            "family_size": 1,
            "retained_traits": ["trait:a"],
            "source_commits": {"core": "c" * 40},
            "core_run_id": f"core-{scope}",
            "comparison_run_id": f"comparison-{scope}",
            "source_run_ids": {"noel": f"noel-{scope}"},
            "analysis_bases": {},
            "inputs": {},
        }
        manifest = root / "three_test_run.json"
        manifest.write_text(json.dumps(record), encoding="utf-8")
        snapshots = {path: _inventory(path) for path in (manifest, evidence, *records.values())}
        sources = SimpleNamespace(
            core=SimpleNamespace(manifest_path=str(records["core"])),
            comparison=SimpleNamespace(manifest_path=str(records["comparison"])),
            arms={
                name: SimpleNamespace(manifest_path=str(records[name]))
                for name in ("noel", "jed", "ilyas")
            },
        )
        bundles.append(
            SimpleNamespace(
                manifest=manifest,
                root=root,
                record=record,
                sources=sources,
                snapshots=snapshots,
                evidence=pd.DataFrame({"trait": ["trait:a"]}),
                pair_delta_z=pd.DataFrame(
                    {"trait:a": [1.0]}, index=pd.MultiIndex.from_tuples([("001", "01")])
                ),
                study_delta_z=pd.DataFrame({"trait:a": [1.0]}, index=["001"]),
            )
        )
        comparisons.append(records["comparison"])
    lookup = {bundle.manifest: bundle for bundle in bundles}
    loader = Mock(side_effect=lambda manifest, _comparison, **_kwargs: lookup[manifest])
    monkeypatch.setattr(stage, "load_completed_table", loader)

    def recheck(bundle: Any) -> None:
        for path, expected in bundle.snapshots.items():
            if _inventory(path) != expected:
                raise ValueError("Parent identity changed")

    monkeypatch.setattr(stage, "recheck_completed_table", recheck)
    comparable = Mock()
    monkeypatch.setattr(stage, "verify_comparable", comparable)
    producer = {"code_commit": "d" * 40, "code_dirty": False, "implementation_sha256": "e" * 64}
    monkeypatch.setattr(stage, "_producer", lambda _commit: producer.copy())
    tables = SimpleNamespace(
        shared_traits=pd.DataFrame({"trait": ["trait:a"], "off_study_t_q_by": [0.012]}),
        trait_catalog=pd.DataFrame({"trait": ["trait:a"], "membership": ["shared"]}),
        rank_transitions=pd.DataFrame(
            {"trait": ["trait:a"], "view": ["noel_by_0_05"], "membership_state": ["retained"]}
        ),
    )
    builder = Mock(return_value=tables)
    monkeypatch.setattr(stage, "build_tables", builder)
    return SimpleNamespace(
        bundles=bundles,
        comparisons=comparisons,
        loader=loader,
        comparable=comparable,
        builder=builder,
        producer=producer,
    )


def run(fixture: Any, output: Path) -> dict[str, Any]:
    off, phenotype = fixture.bundles
    return stage.run_comparison(
        off.manifest,
        fixture.comparisons[0],
        phenotype.manifest,
        fixture.comparisons[1],
        output,
        off_sha256=file_sha256(off.manifest),
        phenotype_sha256=file_sha256(phenotype.manifest),
    )


def test_publishes_exact_evidence_copies_with_both_parent_pins(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fixture = parents(tmp_path, monkeypatch)
    output = tmp_path / "output"
    record = run(fixture, output)
    assert record == json.loads((output / stage.MANIFEST).read_text())
    assert record["primary_scope"] == "phenotype"
    assert record["sensitivity_scope"] == "off"
    for scope, bundle in zip(("off", "phenotype"), fixture.bundles, strict=True):
        assert record["parents"][scope]["manifest_sha256"] == file_sha256(bundle.manifest)
        assert record["parents"][scope]["run_id"] == bundle.record["run_id"]
        copied = output / f"{scope}_per_trait_evidence.tsv"
        assert copied.read_bytes() == (bundle.root / "per_trait_evidence.tsv").read_bytes()
        assert record["outputs"][copied.name] == _inventory(copied)
    fixture.comparable.assert_called_once_with(*fixture.bundles)
    assert set(record["outputs"]) == {
        "off_per_trait_evidence.tsv",
        "phenotype_per_trait_evidence.tsv",
        "shared_traits.tsv",
        "trait_catalog.tsv",
        "rank_transitions.tsv",
    }
    assert str(tmp_path) not in json.dumps(record)


@pytest.mark.parametrize("scope", [0, 1])
def test_protects_both_sources_before_reader_or_probe(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scope: int
) -> None:
    fixture = parents(tmp_path, monkeypatch)
    probe = Mock()
    monkeypatch.setattr(stage, "_probe_publication", probe)
    with pytest.raises(ValueError, match="disjoint"):
        run(fixture, fixture.comparisons[scope].parent / "output")
    fixture.loader.assert_not_called()
    probe.assert_not_called()


def test_incompatible_parents_leave_no_output_or_probe(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fixture = parents(tmp_path, monkeypatch)
    fixture.comparable.side_effect = ValueError("different cohort")
    probe = Mock()
    monkeypatch.setattr(stage, "_probe_publication", probe)
    with pytest.raises(ValueError, match="different cohort"):
        run(fixture, tmp_path / "output")
    probe.assert_not_called()
    assert not (tmp_path / "output").exists()


@pytest.mark.parametrize("scope", [0, 1])
def test_late_parent_change_rejects_and_cleans_own_staging(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scope: int
) -> None:
    fixture = parents(tmp_path, monkeypatch)
    original = shutil.copyfile

    def mutate(source: Path, destination: Path) -> Any:
        result = original(source, destination)
        if Path(source).parent == fixture.bundles[scope].root:
            Path(source).write_text("changed\n", encoding="utf-8")
        return result

    monkeypatch.setattr(shutil, "copyfile", mutate)
    with pytest.raises(ValueError, match="identity"):
        run(fixture, tmp_path / "output")
    assert not (tmp_path / "output").exists()
    assert not list(tmp_path.glob(".output-staging-*"))


def test_corrupt_copy_cannot_pass_restored_source_recheck(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fixture = parents(tmp_path, monkeypatch)
    original = shutil.copyfile

    def corrupt(source: Path, destination: Path) -> Any:
        result = original(source, destination)
        Path(destination).write_text("corrupt\n", encoding="utf-8")
        return result

    monkeypatch.setattr(shutil, "copyfile", corrupt)
    with pytest.raises(ValueError, match="Copied parent evidence"):
        run(fixture, tmp_path / "output")
    assert not (tmp_path / "output").exists()


def test_existing_output_is_preserved(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    fixture = parents(tmp_path, monkeypatch)
    output = tmp_path / "output"
    output.mkdir()
    sentinel = output / "owner.txt"
    sentinel.write_text("keep", encoding="utf-8")
    with pytest.raises(FileExistsError):
        run(fixture, output)
    fixture.loader.assert_not_called()
    assert sentinel.read_text() == "keep"


def test_competing_destination_is_not_overwritten(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fixture = parents(tmp_path, monkeypatch)
    original = _publish_directory
    output = tmp_path / "output"

    def compete(staging: Path, destination: Path) -> None:
        destination.mkdir()
        (destination / "owner.txt").write_text("keep", encoding="utf-8")
        original(staging, destination)

    monkeypatch.setattr(stage, "_publish_directory", compete)
    with pytest.raises(FileExistsError):
        run(fixture, output)
    assert (output / "owner.txt").read_text() == "keep"
    assert not list(tmp_path.glob(".output-staging-*"))


def test_source_alias_target_directory_is_protected(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fixture = parents(tmp_path, monkeypatch)
    alias_root = tmp_path / "aliases"
    alias_root.mkdir()
    alias = alias_root / "phenotype.json"
    alias.symlink_to(fixture.comparisons[1])
    fixture.comparisons[1] = alias
    probe = Mock()
    monkeypatch.setattr(stage, "_probe_publication", probe)
    with pytest.raises(ValueError, match="disjoint"):
        run(fixture, alias_root / "output")
    fixture.loader.assert_not_called()
    probe.assert_not_called()


def test_producer_change_before_publication_is_rejected(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fixture = parents(tmp_path, monkeypatch)
    producer = Mock(
        side_effect=[fixture.producer, fixture.producer, {**fixture.producer, "code_dirty": True}]
    )
    monkeypatch.setattr(stage, "_producer", producer)
    with pytest.raises(ValueError, match="producer changed"):
        run(fixture, tmp_path / "output")
    assert not (tmp_path / "output").exists()
    assert not list(tmp_path.glob(".output-staging-*"))


@pytest.mark.parametrize("scope", [0, 1])
def test_new_validation_inputs_are_protected_before_probe(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scope: int
) -> None:
    fixture = parents(tmp_path, monkeypatch)
    external = tmp_path / "external-identity"
    external.mkdir()
    identity = external / "source_identity.json"
    identity.write_text("{}\n", encoding="utf-8")

    def discover(*_bundles: Any) -> None:
        fixture.bundles[scope].snapshots[identity] = _inventory(identity)

    fixture.comparable.side_effect = discover
    probe = Mock()
    monkeypatch.setattr(stage, "_probe_publication", probe)
    with pytest.raises(ValueError, match="disjoint"):
        run(fixture, external / "output")
    probe.assert_not_called()
    assert not (external / "output").exists()


def test_scope_selection_discloses_both_observed_results(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fixture = parents(tmp_path, monkeypatch)
    record = run(fixture, tmp_path / "output")
    assert record["scope_selection"] == (
        "phenotype_primary_selected_after_off_and_phenotype_results_were_seen"
    )
