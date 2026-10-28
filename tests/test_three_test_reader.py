"""Read completed delta-z bundles without replay or changing their producer (#1052)."""

from __future__ import annotations

import json
import tempfile
from collections.abc import Iterator
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any
from unittest.mock import Mock

import pandas as pd
import pytest
from git import Actor, Repo

from mtp.datamodel.workflow_pydantic import ComparisonRun, ExecutionStatus
from mtp.delta_z_derived import stage as derived_stage
from mtp.multiple_testing import canonical_sha256
from mtp.multiplicity import sources as source_module
from mtp.multiplicity.sources import VerifiedSources
from mtp.paths import file_sha256
from mtp.three_test_table import reader, stage
from mtp.three_test_table.signed_rank import (
    LEGACY_WILCOXON_CONTRACT,
    WILCOXON_CONTRACT,
)
from mtp.workflow.runtime import artifact_path
from multiplicity_fixture import fixture_sources


@pytest.fixture
def tmp_path() -> Iterator[Path]:
    root = Path(__file__).resolve().parents[1] / ".tmp" / "three-test-reader-tests"
    root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="case-", dir=root) as directory:
        yield Path(directory)


@dataclass
class ParentBundle:
    root: Path
    comparison: Path
    record: dict[str, Any]
    sources: VerifiedSources
    repository: Path

    @property
    def manifest(self) -> Path:
        return self.root / stage.MANIFEST

    def save(self) -> str:
        self.manifest.write_text(json.dumps(self.record, sort_keys=True) + "\n", encoding="utf-8")
        return file_sha256(self.manifest)

    def refresh(self, name: str) -> str:
        path = self.root / name
        self.record["outputs"][name] = {
            "sha256": file_sha256(path),
            "size_bytes": path.stat().st_size,
        }
        return self.save()

    def load(self) -> reader.CompletedTable:
        return reader.load_completed_table(
            self.manifest, self.comparison, expected_sha256=file_sha256(self.manifest)
        )


@pytest.fixture
def parent_bundle(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> ParentBundle:
    sources = fixture_sources(tmp_path / "inputs", monkeypatch)
    sources.core.git_commit = "a" * 40
    sources.comparison.id = "fixture-comparison"
    sources.comparison.git_commit = "a" * 40
    sources.comparison.status = ExecutionStatus.completed
    sources.comparison.started_at = datetime(2026, 10, 25, tzinfo=UTC)
    sources.comparison.git_dirty = False
    sources.comparison.analysis_specs = []
    sources.comparison.eligibility_ledgers = []
    sources.comparison.calibration_verdicts = []
    checked = ComparisonRun.model_validate(sources.comparison.model_dump())
    comparison = Path(checked.manifest_path)
    comparison.write_text(checked.model_dump_json(), encoding="utf-8")
    for name, arm in sources.arms.items():
        arm.id = f"fixture-{name}"
        arm.git_commit = "a" * 40
    monkeypatch.setattr(stage, "open_verified_sources", lambda _path: sources)
    monkeypatch.setattr(reader, "open_verified_sources", lambda _path: sources)
    repository_root = tmp_path / "producer"
    with Repo.init(repository_root) as repository:
        names = [
            "src/mtp/three_test_table/stage.py",
            stage.METHODS,
            "uv.lock",
            "pyproject.toml",
            "Justfile",
        ]
        for name in names:
            path = repository_root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("synthetic producer\n", encoding="utf-8")
        repository.index.add(names)
        actor = Actor("fixture", "fixture@example.invalid")
        repository.index.commit("Synthetic parent producer", author=actor, committer=actor)
        monkeypatch.setattr(stage, "__file__", str(repository_root / names[0]))
        monkeypatch.setattr(
            reader, "__file__", str(repository_root / "src/mtp/three_test_table/reader.py")
        )
        output = tmp_path / "parent"
        record = stage.run_table(comparison, output, scratch=tmp_path / "scratch")
        (repository_root / "consumer.txt").write_text("new consumer\n", encoding="utf-8")
        repository.index.add(["consumer.txt"])
        repository.index.commit("Later consumer", author=actor, committer=actor)
        assert repository.head.commit.hexsha != record["code_commit"]
    return ParentBundle(output, comparison, record, sources, repository_root)


def test_reads_older_clean_producer_without_replay_or_recomputing_delta_z(
    parent_bundle: ParentBundle, monkeypatch: pytest.MonkeyPatch
) -> None:
    forbidden = Mock(side_effect=AssertionError("reader must not recompute or replay"))
    monkeypatch.setattr(stage, "run_table", forbidden)
    monkeypatch.setattr(stage, "verify_table", forbidden)
    monkeypatch.setattr(stage, "composition_table", forbidden)
    loaded = parent_bundle.load()
    assert loaded.record == parent_bundle.record
    assert loaded.manifest == parent_bundle.manifest.resolve()
    assert loaded.root == parent_bundle.root.resolve()
    assert loaded.pair_delta_z.index.names == ["study_id", "experiment_id"]
    assert loaded.study_delta_z.index.name == "study_id"
    assert loaded.pair_delta_z.shape == (8, 7)
    assert loaded.study_delta_z.shape == (8, 7)
    assert loaded.evidence.trait.tolist() == list(parent_bundle.sources.features["count_profiles"])
    reader.recheck_completed_table(loaded)
    forbidden.assert_not_called()


def test_reader_preserves_round_trip_float_values(tmp_path: Path) -> None:
    path = tmp_path / "round-trip.tsv"
    expected = 0.1 + 0.2
    pd.DataFrame({"study_id": ["study"], "trait": [expected]}).to_csv(path, sep="\t", index=False)
    frame = reader._read_table(path, ["study_id"], expected=reader._snapshot(path))
    assert float(frame.trait.iloc[0]) == expected


def test_manifest_pin_is_checked_before_source_verification(
    parent_bundle: ParentBundle, monkeypatch: pytest.MonkeyPatch
) -> None:
    forbidden = Mock(side_effect=AssertionError("source verification ran before manifest pin"))
    monkeypatch.setattr(reader, "open_verified_sources", forbidden)
    with pytest.raises(ValueError, match=r"manifest.*(hash|pin)|SHA"):
        reader.load_completed_table(
            parent_bundle.manifest, parent_bundle.comparison, expected_sha256="0" * 64
        )
    forbidden.assert_not_called()


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("version", "unsupported"),
        ("code_dirty", True),
        ("code_commit", "f" * 40),
        ("family_size", 6),
        ("strain_rollup", "phenotype"),
        ("scope_role", "primary"),
        ("contrast_population", "case_vs_control_subject"),
        ("direction_basis", "disease_vs_comparison"),
        ("core_run_id", "other-core"),
        ("comparison_run_id", "other-comparison"),
        ("comparison_manifest_sha256", "0" * 64),
        ("ora_declaration", "undeclared"),
        ("primary_ranking", "posthoc"),
        ("environment", {}),
    ],
)
def test_rejects_rebound_contract_or_identity(
    parent_bundle: ParentBundle, field: str, value: object
) -> None:
    parent_bundle.record[field] = value
    parent_bundle.save()
    with pytest.raises(ValueError, match=r"contract|producer|identity|inventory|environment"):
        parent_bundle.load()


@pytest.mark.parametrize("field", ["source_commits", "source_run_ids", "analysis_bases"])
def test_rejects_source_record_mismatch(parent_bundle: ParentBundle, field: str) -> None:
    parent_bundle.record[field] = {}
    parent_bundle.save()
    with pytest.raises(ValueError, match=r"source|identity|basis"):
        parent_bundle.load()


def test_rejects_rehashed_false_code_inventory(parent_bundle: ParentBundle) -> None:
    inventory = parent_bundle.record["code_inventory"]
    inventory["src/mtp/three_test_table/stage.py"] = "f" * 64
    parent_bundle.record["implementation_sha256"] = canonical_sha256(inventory)
    parent_bundle.save()
    with pytest.raises(ValueError, match=r"producer|inventory|Git"):
        parent_bundle.load()


def test_rejects_incomplete_producer_inventory(parent_bundle: ParentBundle) -> None:
    inventory = parent_bundle.record["code_inventory"]
    del inventory["src/mtp/three_test_table/stage.py"]
    parent_bundle.record["implementation_sha256"] = canonical_sha256(inventory)
    parent_bundle.save()
    with pytest.raises(ValueError, match=r"producer|inventory|Git"):
        parent_bundle.load()


@pytest.mark.parametrize("mutation", ["extra", "missing", "symlink", "unbound", "traversal"])
def test_rejects_output_inventory_changes(parent_bundle: ParentBundle, mutation: str) -> None:
    path = parent_bundle.root / "composition_pair_delta_z.tsv"
    if mutation == "extra":
        (parent_bundle.root / "extra.txt").write_text("unexpected", encoding="utf-8")
    elif mutation == "missing":
        path.unlink()
    elif mutation == "symlink":
        target = parent_bundle.root.parent / "other.tsv"
        path.rename(target)
        path.symlink_to(target)
    elif mutation == "unbound":
        path.write_text(path.read_text(encoding="utf-8") + "\n", encoding="utf-8")
    else:
        item = parent_bundle.record["outputs"].pop(path.name)
        parent_bundle.record["outputs"]["../composition_pair_delta_z.tsv"] = item
        parent_bundle.save()
    with pytest.raises(ValueError, match=r"output|inventory|bundle|hash"):
        parent_bundle.load()


@pytest.mark.parametrize("mutation", ["duplicate", "nonfinite", "missing_pair", "study_mean"])
def test_rejects_rehashed_ledger_corruption(parent_bundle: ParentBundle, mutation: str) -> None:
    name = (
        "composition_study_delta_z.tsv"
        if mutation == "study_mean"
        else "composition_pair_delta_z.tsv"
    )
    path = parent_bundle.root / name
    frame = pd.read_csv(path, sep="\t", keep_default_na=False)
    if mutation == "duplicate":
        frame = pd.concat([frame, frame.iloc[[0]]], ignore_index=True)
    elif mutation == "nonfinite":
        frame.iloc[0, 2] = float("inf")
    elif mutation == "missing_pair":
        frame = frame.iloc[1:]
    else:
        frame.iloc[0, 1] = float(str(frame.iloc[0, 1])) + 0.25
    frame.to_csv(path, sep="\t", index=False)
    parent_bundle.refresh(name)
    with pytest.raises(ValueError, match=r"ledger|identity|population|finite|mean"):
        parent_bundle.load()


def test_rejects_duplicate_raw_tsv_headers(parent_bundle: ParentBundle) -> None:
    path = parent_bundle.root / "composition_pair_delta_z.tsv"
    lines = path.read_text(encoding="utf-8").splitlines()
    columns = lines[0].split("\t")
    columns[-1] = columns[-2]
    lines[0] = "\t".join(columns)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    parent_bundle.refresh(path.name)
    with pytest.raises(ValueError, match=r"column|catalog|header"):
        parent_bundle.load()


@pytest.mark.parametrize(
    ("column", "value"),
    [
        ("trait_label", "invented label"),
        ("noel_studies_up", 123),
        ("noel_q_by_nominal", 0.12345),
        ("composition_pairs", 9),
        ("experiment_mean_delta_z", 999.0),
        ("study_mean_delta_z", 999.0),
        ("experiment_mean_delta_z_rank_descending", 999),
        ("composition_contract", "unsupported"),
    ],
)
def test_rejects_rehashed_summary_corruption(
    parent_bundle: ParentBundle, column: str, value: str | int | float
) -> None:
    path = parent_bundle.root / "per_trait_evidence.tsv"
    frame = pd.read_csv(path, sep="\t", keep_default_na=False)
    frame[column] = frame[column].astype(object)
    frame.loc[0, column] = value
    frame.to_csv(path, sep="\t", index=False)
    parent_bundle.refresh(path.name)
    with pytest.raises(ValueError, match=r"summary|label|Noel|contract|mean|count"):
        parent_bundle.load()


def test_source_verification_failure_is_not_swallowed(
    parent_bundle: ParentBundle, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        reader, "open_verified_sources", Mock(side_effect=ValueError("invalid core execution"))
    )
    with pytest.raises(ValueError, match="invalid core execution"):
        parent_bundle.load()


def test_rejects_fractional_noel_count_within_float_tolerance(parent_bundle: ParentBundle) -> None:
    path = parent_bundle.root / "per_trait_evidence.tsv"
    frame = pd.read_csv(path, sep="\t", keep_default_na=False, float_precision="round_trip")
    frame["noel_studies_up"] = frame["noel_studies_up"].astype(float)
    frame.loc[0, "noel_studies_up"] = float(str(frame.loc[0, "noel_studies_up"])) + 1e-13
    frame.to_csv(path, sep="\t", index=False)
    parent_bundle.refresh(path.name)
    with pytest.raises(ValueError, match=r"Noel.*count"):
        parent_bundle.load()


def test_missing_bound_parent_input_is_rejected(parent_bundle: ParentBundle) -> None:
    identity = next(
        key
        for key in parent_bundle.record["inputs"]
        if key.startswith("trait_direction_screen.tsv@")
    )
    del parent_bundle.record["inputs"][identity]
    parent_bundle.save()
    with pytest.raises(ValueError, match=r"input|bound"):
        parent_bundle.load()


def test_recheck_detects_changes_after_loading(parent_bundle: ParentBundle) -> None:
    loaded = parent_bundle.load()
    path = parent_bundle.root / "per_trait_evidence.tsv"
    path.write_text(path.read_text(encoding="utf-8") + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match=r"output|hash|changed"):
        reader.recheck_completed_table(loaded)


def test_late_input_mutation_during_verification_is_rejected(
    parent_bundle: ParentBundle, monkeypatch: pytest.MonkeyPatch
) -> None:
    def mutate(_comparison: Path) -> VerifiedSources:
        parent_bundle.manifest.write_text("{}\n", encoding="utf-8")
        return parent_bundle.sources

    monkeypatch.setattr(reader, "open_verified_sources", mutate)
    with pytest.raises(ValueError, match=r"manifest|changed|pin"):
        parent_bundle.load()


def external_experiments_alias(
    parent: ParentBundle, monkeypatch: pytest.MonkeyPatch, *, requested_only: bool = False
) -> tuple[Path, Path]:
    """Retain real reader/source checks while placing one bound product behind an alias (#1094)."""
    products = {
        artifact.product_id: parent.sources.product(artifact.product_id)
        for artifact in parent.sources.core.artifacts or []
    }
    original = products["experiments"]
    alias = parent.root.parent / "external-alias" / original.name
    target = parent.root.parent / "external-target" / original.name
    alias.parent.mkdir()
    target.parent.mkdir()
    original.rename(target)
    alias.symlink_to(target)
    products["experiments"] = target if requested_only else alias
    for artifact in parent.sources.core.artifacts or []:
        artifact.path = str(products[artifact.product_id])
        if artifact.product_id == "experiments":
            artifact.requested_path = str(alias)
    core_manifest = Path(parent.sources.core.manifest_path)
    core_payload = json.loads(core_manifest.read_text(encoding="utf-8"))
    core_payload["artifacts"] = [
        artifact.model_dump(mode="json") for artifact in parent.sources.core.artifacts or []
    ]
    core_manifest.write_text(json.dumps(core_payload), encoding="utf-8")
    core_sha = file_sha256(core_manifest)
    parent.sources.comparison.core_manifest_sha256 = core_sha
    for arm in parent.sources.arms.values():
        arm.core_manifest_sha256 = core_sha
    for key in list(parent.record["inputs"]):
        if key.startswith(f"{core_manifest.name}@sha256:"):
            del parent.record["inputs"][key]
    parent.record["inputs"][f"{core_manifest.name}@sha256:{core_sha}"] = {
        "sha256": core_sha,
        "size_bytes": core_manifest.stat().st_size,
    }
    parent.save()
    monkeypatch.setattr(source_module, "artifact_path", artifact_path)
    return alias, target


@pytest.mark.parametrize("contract", [LEGACY_WILCOXON_CONTRACT, WILCOXON_CONTRACT])
def test_accepts_pinned_v1_and_v2_wilcoxon_contracts(
    parent_bundle: ParentBundle, contract: str
) -> None:
    path = parent_bundle.root / "per_trait_evidence.tsv"
    frame = pd.read_csv(path, sep="\t", keep_default_na=False)
    frame["study_wilcoxon_contract"] = contract
    frame.to_csv(path, sep="\t", index=False)
    parent_bundle.refresh(path.name)
    loaded = parent_bundle.load()
    assert set(loaded.evidence["study_wilcoxon_contract"]) == {contract}


@pytest.mark.parametrize("mutation", ["unknown", "mixed"])
def test_rejects_unknown_or_mixed_wilcoxon_contracts(
    parent_bundle: ParentBundle, mutation: str
) -> None:
    path = parent_bundle.root / "per_trait_evidence.tsv"
    frame = pd.read_csv(path, sep="\t", keep_default_na=False)
    frame["study_wilcoxon_contract"] = WILCOXON_CONTRACT
    frame.loc[0, "study_wilcoxon_contract"] = (
        "study-mean-delta-z-wilcoxon-v99" if mutation == "unknown" else LEGACY_WILCOXON_CONTRACT
    )
    frame.to_csv(path, sep="\t", index=False)
    parent_bundle.refresh(path.name)
    with pytest.raises(ValueError, match="study_wilcoxon_contract"):
        parent_bundle.load()


@pytest.mark.parametrize("location", ["alias", "target"])
@pytest.mark.parametrize("requested_only", [False, True])
def test_real_reader_protects_external_product_locations_before_stage_probe(
    parent_bundle: ParentBundle,
    monkeypatch: pytest.MonkeyPatch,
    location: str,
    requested_only: bool,
) -> None:
    alias, target = external_experiments_alias(
        parent_bundle, monkeypatch, requested_only=requested_only
    )
    output = (alias if location == "alias" else target).parent / "derived"
    original_loader = reader.load_completed_table
    loaded_bundles: list[reader.CompletedTable] = []

    def load(parent: Path, comparison: Path, *, expected_sha256: str) -> reader.CompletedTable:
        bundle = original_loader(parent, comparison, expected_sha256=expected_sha256)
        loaded_bundles.append(bundle)
        return bundle

    probe = Mock(side_effect=AssertionError("probe reached a protected source location"))
    monkeypatch.setattr(derived_stage, "load_completed_table", load)
    monkeypatch.setattr(derived_stage, "_probe_publication", probe)
    monkeypatch.setattr(derived_stage, "_producer", lambda _commit: {"code_dirty": False})
    before = {path: path.read_bytes() for path in (alias, target)}
    with pytest.raises(ValueError, match="disjoint"):
        derived_stage.run_derived(
            parent_bundle.manifest,
            parent_bundle.comparison,
            output,
            expected_sha256=file_sha256(parent_bundle.manifest),
            summary_only=True,
        )
    assert len(loaded_bundles) == 1
    assert alias in loaded_bundles[0].snapshots
    assert target in loaded_bundles[0].snapshots
    assert {path: path.read_bytes() for path in (alias, target)} == before
    assert not output.exists()
    assert not list(output.parent.glob(".derived-probe-*"))
    probe.assert_not_called()


@pytest.mark.parametrize("requested_only", [False, True])
def test_reader_recheck_rejects_identical_bytes_at_retargeted_product_alias(
    parent_bundle: ParentBundle, monkeypatch: pytest.MonkeyPatch, requested_only: bool
) -> None:
    alias, target = external_experiments_alias(
        parent_bundle, monkeypatch, requested_only=requested_only
    )
    bundle = parent_bundle.load()
    replacement = target.parent / "same-bytes.tsv"
    replacement.write_bytes(target.read_bytes())
    alias.unlink()
    alias.symlink_to(replacement)
    assert file_sha256(replacement) == file_sha256(target)
    with pytest.raises(ValueError, match=r"alias|target|changed"):
        reader.recheck_completed_table(bundle)


@pytest.mark.parametrize("argument", ["manifest", "comparison"])
def test_reader_retains_supplied_manifest_alias_and_rechecks_its_target(
    parent_bundle: ParentBundle, argument: str
) -> None:
    target = parent_bundle.manifest if argument == "manifest" else parent_bundle.comparison
    alias = parent_bundle.root.parent / f"{argument}-alias" / target.name
    alias.parent.mkdir()
    alias.symlink_to(target)
    bundle = reader.load_completed_table(
        alias if argument == "manifest" else parent_bundle.manifest,
        alias if argument == "comparison" else parent_bundle.comparison,
        expected_sha256=file_sha256(parent_bundle.manifest),
    )
    assert alias in bundle.snapshots
    assert target in bundle.snapshots
    replacement = alias.parent / "replacement.json"
    replacement.write_bytes(target.read_bytes())
    alias.unlink()
    alias.symlink_to(replacement)
    with pytest.raises(ValueError, match=r"alias|target|changed"):
        reader.recheck_completed_table(bundle)


def test_supplied_parent_alias_preserves_symlink_dotdot_resolution(
    parent_bundle: ParentBundle,
) -> None:
    route = parent_bundle.root.parent / "route"
    target = parent_bundle.root.parent / "bridge-target"
    route.mkdir()
    target.mkdir()
    link = route / "link"
    link.symlink_to(target, target_is_directory=True)
    supplied = link / ".." / parent_bundle.root.name / stage.MANIFEST
    assert supplied.resolve() == parent_bundle.manifest
    bundle = reader.load_completed_table(
        supplied,
        parent_bundle.comparison,
        expected_sha256=file_sha256(parent_bundle.manifest),
    )
    assert bundle.manifest == parent_bundle.manifest
    assert supplied in bundle.snapshots
    reader.recheck_completed_table(bundle)


def test_reader_parses_the_hashed_evidence_bytes_during_transient_path_change(
    parent_bundle: ParentBundle, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The filesystem can return to its original hash after parsing false evidence."""
    path = parent_bundle.root / "per_trait_evidence.tsv"
    original = path.read_bytes()
    rows = original.decode("utf-8").splitlines()
    transient = (
        "\n".join(
            f"{line}\t{'unbound_read_marker' if index == 0 else 'injected'}"
            for index, line in enumerate(rows)
        )
        + "\n"
    ).encode("utf-8")
    original_reader = pd.read_csv

    def transient_reader(source: Any, *args: Any, **kwargs: Any) -> pd.DataFrame:
        if isinstance(source, Path) and source == path:
            path.write_bytes(transient)
            try:
                result: pd.DataFrame = original_reader(source, *args, **kwargs)
                return result
            finally:
                path.write_bytes(original)
        result = original_reader(source, *args, **kwargs)
        return result

    monkeypatch.setattr(pd, "read_csv", transient_reader)
    loaded = parent_bundle.load()
    assert path.read_bytes() == original
    assert "unbound_read_marker" not in loaded.evidence.columns


@pytest.mark.parametrize("name", reader._OUTPUTS)
def test_reader_checks_the_parsed_buffer_against_the_bound_output_identity(
    parent_bundle: ParentBundle, monkeypatch: pytest.MonkeyPatch, name: str
) -> None:
    """A changed read buffer is not validated by hashing the restored on-disk file."""
    path = parent_bundle.root / name
    original_reader = Path.read_bytes

    def transient_bytes(source: Path) -> bytes:
        raw = original_reader(source)
        return raw + b"\n" if source == path else raw

    monkeypatch.setattr(Path, "read_bytes", transient_bytes)
    with pytest.raises(ValueError, match=r"Completed-table.*(hash|identity|size)"):
        parent_bundle.load()
