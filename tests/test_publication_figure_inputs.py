"""#1046: pinned synthetic figure inputs, historical producer and read-race checks."""

from __future__ import annotations

import json
import shutil
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from unittest.mock import Mock

import pandas as pd
import pytest
from git import Actor, Repo

from mtp.multiple_testing import canonical_sha256
from mtp.multiplicity.sources import VerifiedSources
from mtp.paths import file_sha256
from mtp.publication_figures import inputs
from mtp.three_test_table.reader import CompletedTable

FA = "METPO:1000608"
OA = "METPO:1000607"
PAIR = "composition_pair_delta_z.tsv"
STUDY = "composition_study_delta_z.tsv"
MANIFEST = "delta_z_derived_run.json"
METHODS = "methods/delta_z_derived.md"
ORA_COUNTS = {
    "ora_background_k": 2,
    "ora_direction_studies_up": 2,
    "ora_direction_studies_down": 1,
    "ora_direction_studies_tied": 1,
    "ora_direction_studies_scored": 3,
    "ora_direction_studies_eligible": 4,
    "ora_direction_studies_excluded": 0,
    "ora_direction_studies_total": 4,
    "ora_direction_pairs_eligible": 5,
    "ora_direction_pairs_excluded": 0,
    "ora_direction_pairs_total": 5,
}


def identity(path: Path) -> dict[str, Any]:
    return {"sha256": file_sha256(path), "size_bytes": path.stat().st_size}


@dataclass
class Fixture:
    parent: CompletedTable
    comparison: Path
    root: Path
    record: dict[str, Any]
    loader: Mock
    validator: Mock
    metadata: Path
    experiments: Path

    @property
    def manifest(self) -> Path:
        return self.root / MANIFEST

    def seal(self) -> str:
        self.record["outputs"] = {
            path.name: identity(path) for path in sorted(self.root.glob("*.tsv"))
        }
        self.manifest.write_text(json.dumps(self.record, sort_keys=True), encoding="utf-8")
        return file_sha256(self.manifest)

    def load(self, **kwargs: Any) -> inputs.FigureInputs:
        options = {
            "three_test_sha256": file_sha256(self.parent.manifest),
            "derived_sha256": file_sha256(self.manifest),
            **kwargs,
        }
        return inputs.load_inputs(self.parent.manifest, self.comparison, self.root, **options)

    def rebound_ora(self, changes: dict[str, Any]) -> None:
        """Rebind an inherited row and its copied summary, avoiding a hash-only rejection."""
        evidence = self.parent.evidence
        for name, value in changes.items():
            evidence[name] = evidence[name].astype(object)
            evidence.at[0, name] = value
        evidence_path = self.parent.root / "per_trait_evidence.tsv"
        evidence.to_csv(evidence_path, sep="\t", index=False)
        parsed = pd.read_csv(
            evidence_path,
            sep="\t",
            keep_default_na=False,
            float_precision="round_trip",
            dtype={"trait": str},
        )
        for column in parsed:
            evidence[column] = parsed[column]
        summary_path = self.root / "trait_summary.tsv"
        summary = pd.read_csv(
            summary_path,
            sep="\t",
            keep_default_na=False,
            float_precision="round_trip",
            dtype={"trait": str},
        )
        for column in evidence:
            summary[column] = evidence[column]
        summary.to_csv(summary_path, sep="\t", index=False)
        self.parent.record["outputs"][evidence_path.name] = identity(evidence_path)
        self.parent.manifest.write_text(json.dumps(self.parent.record), encoding="utf-8")
        for path in (evidence_path, self.parent.manifest):
            self.parent.snapshots[path] = identity(path)
        self.record["parent_manifest_sha256"] = file_sha256(self.parent.manifest)
        self.record["inputs"] = {
            f"{path.name}@sha256:{item['sha256']}": item
            for path, item in self.parent.snapshots.items()
        }
        self.seal()


@pytest.fixture
def fixture(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Fixture:
    """Stub only the separately tested completed-parent verification boundary."""
    root, parent_root, source = (tmp_path / name for name in ("derived", "parent", "source"))
    for path in (root, parent_root, source):
        path.mkdir()
    keys = [("001", "01"), ("001", "02"), ("002", "01"), ("003", "01"), ("004", "01")]
    pairs = pd.DataFrame(
        {FA: [0.1 + 0.2, -1.0, 2.0, -2.0, 0.0], OA: [-1.0, 2.0, 3.0, -4.0, 1.0]},
        index=pd.MultiIndex.from_tuples(keys, names=["study_id", "experiment_id"]),
    )
    studies = pairs.groupby(level="study_id").mean()
    labels = {FA: "Synthetic facultative annotation", OA: "Synthetic obligate annotation"}
    evidence = pd.DataFrame(
        {
            "trait": [FA, OA],
            "trait_label": list(labels.values()),
            "experiment_mean_delta_z": pairs.mean().to_numpy(),
            "study_mean_delta_z": studies.mean().to_numpy(),
        }
    )
    for name, count in ORA_COUNTS.items():
        evidence[name] = count
    evidence["ora_direction_status"] = "tested"
    for name in ("ora_enriched_inc", "ora_enriched_dec", "ora_depleted_inc", "ora_depleted_dec"):
        evidence[name] = 1
    for name in ("ora_direction_p", "ora_direction_bh_q", "ora_direction_by_q"):
        evidence[name] = 0.5
    pairs.to_csv(parent_root / PAIR, sep="\t")
    studies.to_csv(parent_root / STUDY, sep="\t")
    evidence.to_csv(parent_root / "per_trait_evidence.tsv", sep="\t", index=False)
    metadata = source / "feature_metadata.tsv"
    pd.DataFrame(
        {
            "column_id": [FA, OA],
            "object_name": [labels[FA], ""],
            "column_label": ["Unused fallback", labels[OA]],
        }
    ).to_csv(metadata, sep="\t", index=False)
    experiments = source / "experiments.tsv"
    frame = pd.DataFrame(keys, columns=["study_id", "experiment_id"])
    frame["body_site"] = ["Feces", "Feces,Blood", "Feces", "Saliva", "Feces"]
    frame["condition"] = ["NA", "NA", "Disease B", "Disease C", "Disease D"]
    frame.to_csv(experiments, sep="\t", index=False)
    comparison = source / "comparison_run.json"
    comparison.write_text("{}\n", encoding="utf-8")
    parent_record: dict[str, Any] = {
        "run_id": "three-test-fixture",
        "code_commit": "a" * 40,
        "code_dirty": False,
        "core_run_id": "core-fixture",
        "comparison_run_id": "comparison-fixture",
        "source_run_ids": {name: f"arm-{name}" for name in ("noel", "jed", "ilyas")},
        "source_commits": dict.fromkeys(("core", "comparison", "noel", "jed", "ilyas"), "b" * 40),
        "strain_rollup": "phenotype",
        "scope_role": "primary",
        "contrast_population": "disease_vs_healthy_tier1_2_subject_controls",
        "direction_basis": "disease_vs_comparison",
        "analysis_bases": {"count_profiles": {"id": "basis"}},
        "family_size": 2,
        "retained_traits": [FA, OA],
        "outputs": {
            name: identity(parent_root / name) for name in (PAIR, STUDY, "per_trait_evidence.tsv")
        },
        "inputs": {
            f"{path.name}@sha256:{file_sha256(path)}": identity(path)
            for path in (metadata, experiments, comparison)
        },
    }
    parent_manifest = parent_root / "three_test_run.json"
    parent_manifest.write_text(json.dumps(parent_record), encoding="utf-8")
    snapshots = {
        path: identity(path)
        for path in (
            parent_manifest,
            parent_root / PAIR,
            parent_root / STUDY,
            parent_root / "per_trait_evidence.tsv",
            metadata,
            experiments,
            comparison,
        )
    }
    sources = Mock(spec=VerifiedSources)
    sources.product.side_effect = {
        "feature_metadata": metadata,
        "experiments": experiments,
    }.__getitem__
    sources.expected_digest.side_effect = lambda path: snapshots[path]["sha256"]
    parent = CompletedTable(
        parent_manifest,
        parent_root,
        parent_record,
        sources,
        evidence,
        pairs,
        studies,
        snapshots,
        {path: path for path in snapshots},
    )
    loader = Mock(return_value=parent)
    monkeypatch.setattr(inputs, "load_completed_table", loader)
    validator = Mock()
    monkeypatch.setattr(inputs, "validate_evidence", validator, raising=False)

    def recheck(bundle: CompletedTable) -> None:
        for path, expected in bundle.snapshots.items():
            if identity(path) != expected or path.resolve() != bundle.resolved_targets[path]:
                raise ValueError("Parent input changed")

    monkeypatch.setattr(inputs, "recheck_completed_table", recheck)
    repository_root = tmp_path / "producer"
    names = [
        "src/mtp/delta_z_derived/stage.py",
        "methods/three_test_table.md",
        METHODS,
        "uv.lock",
        "pyproject.toml",
        "Justfile",
    ]
    with Repo.init(repository_root) as repository:
        for name in names:
            path = repository_root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(f"Synthetic producer: {name}\n", encoding="utf-8")
        repository.index.add(names)
        actor = Actor("fixture", "fixture@example.invalid")
        commit = repository.index.commit(
            "Synthetic historical producer", author=actor, committer=actor
        )
        inventory = {name: file_sha256(repository_root / name) for name in names}
        (repository_root / "src/mtp/publication_figures").mkdir()
        monkeypatch.setattr(
            inputs, "__file__", str(repository_root / "src/mtp/publication_figures/inputs.py")
        )
        (repository_root / "later.txt").write_text("later consumer\n", encoding="utf-8")
        repository.index.add(["later.txt"])
        repository.index.commit("Later consumer", author=actor, committer=actor)
    for name in (PAIR, STUDY):
        shutil.copyfile(parent_root / name, root / name)
    summary = evidence.copy()
    summary["sdz_up"] = studies.gt(1e-12).sum().to_numpy()
    summary["sdz_down"] = studies.lt(-1e-12).sum().to_numpy()
    summary["sdz_active_studies"] = summary.sdz_up + summary.sdz_down
    summary["sdz_tied"] = len(studies) - summary.sdz_active_studies
    summary["sdz_studies_total"] = len(studies)
    summary["sdz_active_threshold"] = 1e-12
    summary["sdz_active_at_least_50"] = False
    summary["sdz_quantity"] = "sign_of_study_mean_delta_z"
    summary["noel_quantity"] = "native_raw_fraction_direction"
    summary.to_csv(root / "trait_summary.tsv", sep="\t", index=False)
    quadrant = frame.copy()
    quadrant["fa_trait"], quadrant["oa_trait"] = FA, OA
    quadrant["fa_delta_z"], quadrant["oa_delta_z"] = pairs[FA].to_numpy(), pairs[OA].to_numpy()
    quadrant["quadrant"] = ["fa_up_oa_down", "fa_down_oa_up", "both_up", "both_down", "axis"]
    quadrant["is_stool"] = frame.body_site.eq("Feces")
    quadrant.to_csv(root / "quadrant_pairs.tsv", sep="\t", index=False)
    summary_rows = []
    for population, count, study_count, q2, down in (("all", 5, 4, 1, 1), ("stool", 3, 3, 0, 0)):
        summary_rows.append(
            {
                "population": population,
                "fa_trait": FA,
                "oa_trait": OA,
                "fa_label": labels[FA],
                "oa_label": labels[OA],
                "fa_up_oa_down": 1,
                "fa_down_oa_up": q2,
                "both_up": 1,
                "both_down": down,
                "axis": 1,
                "n_comparisons": count,
                "n_studies": study_count,
                "sign_test_n": 1 + q2,
                "sign_threshold": 0.0,
                "test_unit": "comparison",
                "p_adjustment": "none_nominal",
                "stool_field": "experiments.body_site",
                "stool_value": "Feces",
                "sign_test_status": "tested",
                "sign_test_p": 0.5,
                "pearson_status": "tested",
                "pearson_r": 0.2,
                "pearson_p": 0.5,
            }
        )
    pd.DataFrame(summary_rows).to_csv(root / "quadrant_summary.tsv", sep="\t", index=False)
    for name in (
        "stool_disease_shifts.tsv",
        "stool_study_sensitivity.tsv",
        "site_disease_shifts.tsv",
    ):
        pd.DataFrame({"trait": [FA], "mean_delta_z": [0.5]}).to_csv(
            root / name, sep="\t", index=False
        )
    record = {
        name: deepcopy(parent_record[name])
        for name in (
            "core_run_id",
            "comparison_run_id",
            "source_run_ids",
            "source_commits",
            "strain_rollup",
            "scope_role",
            "contrast_population",
            "direction_basis",
            "analysis_bases",
            "family_size",
            "retained_traits",
        )
    }
    record.update(
        version="delta-z-derived-v1",
        run_id="derived-fixture",
        code_commit=commit.hexsha,
        code_dirty=False,
        code_inventory=inventory,
        implementation_sha256=canonical_sha256(inventory),
        methods_sha256=inventory[METHODS],
        lock_sha256=inventory["uv.lock"],
        environment=dict.fromkeys(("python", "numpy", "pandas", "scipy"), "1.2.3"),
        parent_run_id=parent_record["run_id"],
        parent_manifest_sha256=file_sha256(parent_manifest),
        parent_producer={name: parent_record[name] for name in ("code_commit", "code_dirty")},
        parent_input_inventory=parent_record["inputs"],
        summary_only=False,
        inputs={f"{path.name}@sha256:{item['sha256']}": item for path, item in snapshots.items()},
        counts={
            "traits": 2,
            "comparisons": 5,
            "studies": 4,
            "traits_with_at_least_50_active_studies": 0,
        },
        policies={
            "row2_x_trait": FA,
            "row2_y_trait": OA,
            "row2_zero_policy": "exact_zero",
            "stool_metadata_field": "body_site",
            "stool_metadata_value": "Feces",
            "row1_activity_epsilon": 1e-12,
            "row1_minimum_active_studies": 50,
        },
    )
    result = Fixture(parent, comparison, root, record, loader, validator, metadata, experiments)
    result.seal()
    return result


def test_reads_complete_historical_bundle_without_replay(
    fixture: Fixture, monkeypatch: pytest.MonkeyPatch
) -> None:
    from mtp.delta_z_derived import stage

    forbidden = Mock(side_effect=AssertionError("No inference or Delta-z replay"))
    monkeypatch.setattr(stage, "run_derived", forbidden)
    monkeypatch.setattr(stage, "build_tables", forbidden)
    loaded = fixture.load()
    assert loaded.parent is fixture.parent
    assert len(loaded.tables) == 8
    assert loaded.labels[FA] == "Synthetic facultative annotation"
    assert loaded.labels[OA] == "Synthetic obligate annotation"
    pairs = loaded.tables["quadrant_pairs.tsv"]
    assert len(pairs) == 5
    assert pairs.study_id.tolist() == ["001", "001", "002", "003", "004"]
    assert float(pairs.fa_delta_z.iloc[0]) == 0.1 + 0.2
    assert pairs.quadrant.tolist()[-1] == "axis"
    assert len(pairs.loc[pairs.is_stool]) == 3
    inputs.recheck_inputs(loaded)
    fixture.loader.assert_called_once_with(
        fixture.parent.manifest,
        fixture.comparison,
        expected_sha256=file_sha256(fixture.parent.manifest),
    )
    forbidden.assert_not_called()


def test_derived_pin_refused_before_parent_verification(fixture: Fixture) -> None:
    with pytest.raises(ValueError, match=r"pin|SHA|hash"):
        fixture.load(derived_sha256="0" * 64)
    fixture.loader.assert_not_called()


def test_parent_verification_failure_propagates(fixture: Fixture) -> None:
    fixture.loader.side_effect = ValueError("Parent source pin refused")
    with pytest.raises(ValueError, match="Parent source pin"):
        fixture.load()


def test_uses_native_persisted_evidence_validator(fixture: Fixture) -> None:
    fixture.load()
    fixture.validator.assert_called_once_with(fixture.parent)


def test_native_persisted_evidence_rejection_propagates(fixture: Fixture) -> None:
    fixture.validator.side_effect = ValueError("Inherited t statistic/CI disagrees")
    with pytest.raises(ValueError, match="Inherited t statistic/CI"):
        fixture.load()


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("version", "other"),
        ("summary_only", True),
        ("code_dirty", True),
        ("code_commit", "f" * 40),
        ("parent_run_id", "other"),
        ("parent_manifest_sha256", "0" * 64),
        ("parent_producer", {}),
        ("core_run_id", "other"),
        ("comparison_run_id", "other"),
        ("source_run_ids", {}),
        ("source_commits", {}),
        ("strain_rollup", "off"),
        ("scope_role", "sensitivity"),
        ("contrast_population", "all_paired"),
        ("direction_basis", "group_1_vs_group_0"),
        ("analysis_bases", {}),
        ("family_size", 3),
        ("retained_traits", [OA, FA]),
        ("parent_input_inventory", {}),
        ("inputs", {}),
        ("environment", {}),
        ("lock_sha256", "0" * 64),
    ],
)
def test_rejects_identity_or_contract_changes(fixture: Fixture, field: str, value: Any) -> None:
    fixture.record[field] = value
    fixture.seal()
    with pytest.raises(ValueError, match=r"Derived|derived|producer"):
        fixture.load()


@pytest.mark.parametrize("mutation", ["false", "missing"])
def test_producer_inventory_matches_historical_git(fixture: Fixture, mutation: str) -> None:
    if mutation == "false":
        fixture.record["code_inventory"][METHODS] = "f" * 64
    else:
        fixture.record["code_inventory"].pop(METHODS)
    fixture.record["implementation_sha256"] = canonical_sha256(fixture.record["code_inventory"])
    fixture.seal()
    with pytest.raises(ValueError, match=r"producer|inventory|Git"):
        fixture.load()


@pytest.mark.parametrize("mutation", ["missing", "extra", "symlink", "unbound", "traversal"])
def test_exact_derived_inventory(fixture: Fixture, mutation: str) -> None:
    path = fixture.root / "quadrant_pairs.tsv"
    if mutation == "missing":
        path.unlink()
    elif mutation == "extra":
        (fixture.root / "extra.txt").write_text("unlisted", encoding="utf-8")
    elif mutation == "symlink":
        target = fixture.root.parent / "linked.tsv"
        path.rename(target)
        path.symlink_to(target)
    elif mutation == "unbound":
        path.write_bytes(path.read_bytes() + b"\n")
    else:
        fixture.record["outputs"]["../escape.tsv"] = fixture.record["outputs"].pop(path.name)
        fixture.manifest.write_text(json.dumps(fixture.record), encoding="utf-8")
    with pytest.raises(ValueError, match=r"inventory|bundle|hash|identity|tables"):
        fixture.load()


@pytest.mark.parametrize(
    "mutation", ["dropped", "duplicate", "changed", "reversed", "axis", "stool", "key"]
)
def test_rehashed_quadrant_corruption_refused(fixture: Fixture, mutation: str) -> None:
    path = fixture.root / "quadrant_pairs.tsv"
    frame = pd.read_csv(
        path,
        sep="\t",
        dtype={"study_id": str, "experiment_id": str},
        keep_default_na=False,
        float_precision="round_trip",
    )
    if mutation == "dropped":
        frame = frame.iloc[:-1]
    elif mutation == "duplicate":
        frame = pd.concat([frame, frame.iloc[[0]]], ignore_index=True)
    elif mutation == "changed":
        frame.loc[0, "fa_delta_z"] = 99.0
    elif mutation == "reversed":
        frame["fa_trait"], frame["oa_trait"] = OA, FA
    elif mutation == "axis":
        frame.loc[4, "quadrant"] = "both_up"
    elif mutation == "stool":
        frame.loc[1, "is_stool"] = True
    else:
        frame.loc[0, "experiment_id"] = "other"
    frame.to_csv(path, sep="\t", index=False)
    fixture.seal()
    with pytest.raises(ValueError, match=r"quadrant|Quadrant|pair|Pair|stool"):
        fixture.load()


@pytest.mark.parametrize(
    ("name", "column"),
    [
        ("quadrant_summary.tsv", "n_comparisons"),
        ("quadrant_summary.tsv", "fa_up_oa_down"),
        ("trait_summary.tsv", "study_mean_delta_z"),
        ("trait_summary.tsv", "sdz_active_studies"),
    ],
)
def test_rehashed_summaries_refused(fixture: Fixture, name: str, column: str) -> None:
    path = fixture.root / name
    frame = pd.read_csv(path, sep="\t", keep_default_na=False, float_precision="round_trip")
    frame.loc[0, column] = 99
    frame.to_csv(path, sep="\t", index=False)
    fixture.seal()
    with pytest.raises(ValueError, match=r"summary|activity"):
        fixture.load()


def test_rehashed_copied_parent_ledger_refused(fixture: Fixture) -> None:
    path = fixture.root / PAIR
    path.write_bytes(path.read_bytes() + b"\n")
    fixture.seal()
    with pytest.raises(ValueError, match=r"parent|copied|ledger"):
        fixture.load()


def test_changed_feature_metadata_refused(fixture: Fixture) -> None:
    fixture.metadata.write_text("column_id\tobject_name\tcolumn_label\n", encoding="utf-8")
    with pytest.raises(ValueError, match=r"hash|identity|changed|metadata"):
        fixture.load()


def test_parses_checked_bytes_not_a_second_path_read(
    fixture: Fixture, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = fixture.root / "quadrant_pairs.tsv"
    original = path.read_bytes()
    changed = original.replace(b"0.30000000000000004", b"9.30000000000000004")
    assert changed != original
    original_read = Path.read_bytes
    calls = 0

    def read(candidate: Path) -> bytes:
        nonlocal calls
        if candidate == path:
            calls += 1
            return changed if calls == 1 else original
        return original_read(candidate)

    monkeypatch.setattr(Path, "read_bytes", read)
    with pytest.raises(ValueError, match=r"hash|identity"):
        fixture.load()


def test_aliases_preserved_and_identical_byte_retarget_refused(fixture: Fixture) -> None:
    alias = fixture.root.parent / "derived-link"
    alias.symlink_to(fixture.root, target_is_directory=True)
    loaded = inputs.load_inputs(
        fixture.parent.manifest,
        fixture.comparison,
        alias,
        three_test_sha256=file_sha256(fixture.parent.manifest),
        derived_sha256=file_sha256(fixture.manifest),
    )
    assert alias / MANIFEST in loaded.snapshots
    assert fixture.manifest in loaded.snapshots
    alternate = fixture.root.parent / "identical-derived"
    shutil.copytree(fixture.root, alternate)
    alias.unlink()
    alias.symlink_to(alternate, target_is_directory=True)
    with pytest.raises(ValueError, match=r"alias|target|changed"):
        inputs.recheck_inputs(loaded)


def test_recheck_refuses_late_output_and_parent_changes(fixture: Fixture) -> None:
    loaded = fixture.load()
    target = fixture.root / "quadrant_summary.tsv"
    target.write_bytes(target.read_bytes() + b"\n")
    with pytest.raises(ValueError, match=r"hash|identity|changed"):
        inputs.recheck_inputs(loaded)


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("ora_direction_studies_up", -1),
        ("ora_direction_studies_up", 0.5),
        ("ora_direction_studies_up", True),
        ("ora_direction_studies_up", None),
        ("ora_direction_studies_eligible", 1),
        ("ora_direction_studies_scored", 4),
        ("ora_direction_studies_excluded", 1),
        ("ora_direction_studies_total", 5),
        ("ora_direction_pairs_total", 6),
        ("ora_direction_pairs_excluded", 1),
        ("ora_direction_pairs_eligible", 1),
        ("ora_background_k", 0),
        ("ora_direction_status", "unknown"),
        ("ora_direction_status", "all_tied"),
        ("ora_direction_status", "no_eligible_pairs"),
        ("ora_direction_status", "no_background_taxa"),
    ],
)
def test_inherited_ora_incoherence_refused(fixture: Fixture, field: str, value: Any) -> None:
    fixture.rebound_ora({field: value})
    with pytest.raises(ValueError, match="ORA"):
        fixture.load()


@pytest.mark.parametrize(
    "status", ["tested", "all_tied", "no_eligible_pairs", "no_background_taxa"]
)
def test_inherited_ora_valid_status_and_missing_fields_preserved(
    fixture: Fixture, status: str
) -> None:
    changes: dict[str, Any] = {"ora_direction_status": status}
    if status != "tested":
        changes.update(
            ora_direction_studies_up=0, ora_direction_studies_down=0, ora_direction_studies_scored=0
        )
        for name in ("ora_direction_p", "ora_direction_bh_q", "ora_direction_by_q"):
            changes[name] = ""
    if status == "all_tied":
        changes["ora_direction_studies_tied"] = 4
    if status in {"no_eligible_pairs", "no_background_taxa"}:
        changes.update(
            ora_direction_studies_tied=0,
            ora_direction_studies_eligible=0,
            ora_direction_studies_excluded=4,
            ora_direction_pairs_eligible=0,
            ora_direction_pairs_excluded=5,
        )
    if status == "no_background_taxa":
        changes["ora_background_k"] = 0
        for name in (
            "ora_enriched_inc",
            "ora_enriched_dec",
            "ora_depleted_inc",
            "ora_depleted_dec",
        ):
            changes[name] = ""
    fixture.rebound_ora(changes)
    loaded = fixture.load()
    row = loaded.tables["trait_summary.tsv"].iloc[0]
    assert row.ora_direction_status == status
    if status != "tested":
        for name in ("ora_direction_p", "ora_direction_bh_q", "ora_direction_by_q"):
            assert row[name] == ""
    if status == "no_background_taxa":
        for name in (
            "ora_enriched_inc",
            "ora_enriched_dec",
            "ora_depleted_inc",
            "ora_depleted_dec",
        ):
            assert row[name] == ""
