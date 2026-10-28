"""Publication and provenance regressions for the Delta-z consumer (#1052)."""

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

from mtp.delta_z_derived import stage
from mtp.paths import file_sha256
from mtp.three_test_table.stage import _publish_directory

FA = "METPO:1000608"
OA = "METPO:1000607"
PAIR_FILE = "composition_pair_delta_z.tsv"
STUDY_FILE = "composition_study_delta_z.tsv"


@pytest.fixture
def tmp_path() -> Iterator[Path]:
    """Use durable synthetic artifacts independently of the process temporary root."""
    root = Path(__file__).resolve().parents[1] / ".tmp" / "delta-z-stage-tests"
    root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="case-", dir=root) as directory:
        yield Path(directory)


def synthetic_parent(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Any:
    """Stub only the separately tested completed-table reader and pure calculations."""
    parent = tmp_path / "parent"
    parent.mkdir()
    source = tmp_path / "sources"
    source.mkdir()
    records: dict[str, Path] = {}
    for name in ("core", "comparison", "noel", "jed", "ilyas"):
        directory = source / name
        directory.mkdir()
        records[name] = directory / "run.json"
        records[name].write_text("{}\n", encoding="utf-8")
    comparison_record = {
        "manifest_path": str(records["comparison"]),
        "core_manifest_path": str(records["core"]),
        **{f"{name}_manifest_path": str(records[name]) for name in ("noel", "jed", "ilyas")},
    }
    records["comparison"].write_text(json.dumps(comparison_record), encoding="utf-8")
    experiments = records["core"].parent / "experiments.tsv"
    experiments.write_text(
        "study_id\texperiment_id\tbody_site\tcondition\n"
        "001\t01\tFeces\tNA\n002\t02\tSaliva\tCondition B\n",
        encoding="utf-8",
    )
    pairs = pd.DataFrame(
        {FA: [1.0, -2.0], OA: [-1.0, 2.0]},
        index=pd.MultiIndex.from_tuples(
            [("001", "01"), ("002", "02")], names=["study_id", "experiment_id"]
        ),
    )
    studies = pairs.groupby(level="study_id").mean()
    evidence = pd.DataFrame(
        {"trait": [FA, OA], "trait_label": ["facultatively aerobic", "obligately anaerobic"]}
    )
    # Noncanonical float spelling makes a dataframe round-trip observably different.
    (parent / PAIR_FILE).write_text(
        f"study_id\texperiment_id\t{FA}\t{OA}\n001\t01\t1.0000\t-1.0000\n"
        "002\t02\t-2.0000\t2.0000\n",
        encoding="utf-8",
    )
    (parent / STUDY_FILE).write_text(
        f"study_id\t{FA}\t{OA}\n001\t1.0000\t-1.0000\n002\t-2.0000\t2.0000\n",
        encoding="utf-8",
    )
    evidence.to_csv(parent / "per_trait_evidence.tsv", sep="\t", index=False)
    record = {
        "version": "three-test-per-trait-v1",
        "run_id": "parent-fixture",
        "code_commit": "a" * 40,
        "code_dirty": False,
        "core_run_id": "core-fixture",
        "comparison_run_id": "comparison-fixture",
        "source_run_ids": {name: f"{name}-fixture" for name in ("noel", "jed", "ilyas")},
        "source_commits": dict.fromkeys(records, "b" * 40),
        "strain_rollup": "phenotype",
        "scope_role": "primary",
        "contrast_population": "disease_vs_healthy_tier1_2_subject_controls",
        "direction_basis": "disease_vs_comparison",
        "family_size": 2,
        "retained_traits": [FA, OA],
        "inputs": {"fixture_input": {"sha256": "c" * 64, "size_bytes": 1}},
        "analysis_bases": {"count_profiles": {"id": "fixture-basis"}},
    }
    manifest = parent / "three_test_run.json"
    manifest.write_text(json.dumps(record), encoding="utf-8")
    paths = [manifest, parent / PAIR_FILE, parent / STUDY_FILE, *records.values(), experiments]
    snapshots = {
        path: {"sha256": file_sha256(path), "size_bytes": path.stat().st_size} for path in paths
    }
    sources = SimpleNamespace(
        core=SimpleNamespace(
            manifest_path=str(records["core"]),
            definition=SimpleNamespace(products=[]),
            artifacts=[],
        ),
        comparison=SimpleNamespace(manifest_path=str(records["comparison"])),
        arms={
            name: SimpleNamespace(manifest_path=str(records[name]))
            for name in ("noel", "jed", "ilyas")
        },
        product=lambda name: experiments if name == "experiments" else None,
    )
    bundle = SimpleNamespace(
        manifest=manifest,
        root=parent,
        record=record,
        sources=sources,
        evidence=evidence,
        pair_delta_z=pairs,
        study_delta_z=studies,
        snapshots=snapshots,
    )
    loader = Mock(return_value=bundle)
    monkeypatch.setattr(stage, "load_completed_table", loader)

    def recheck(_bundle: Any) -> None:
        for path, expected in snapshots.items():
            if file_sha256(path) != expected["sha256"]:
                raise ValueError("Parent source identity changed")

    monkeypatch.setattr(stage, "recheck_completed_table", recheck)
    producer = {
        "code_commit": "d" * 40,
        "code_dirty": False,
        "code_inventory": {"src/mtp/delta_z_derived/stage.py": "e" * 64},
        "implementation_sha256": "e" * 64,
        "methods_sha256": "f" * 64,
        "lock_sha256": "0" * 64,
        "environment": {"python": "fixture"},
    }
    monkeypatch.setattr(stage, "_producer", lambda _commit: producer.copy())
    summary = pd.DataFrame(
        {
            "trait": [FA, OA],
            "sdz_up": [1, 1],
            "sdz_down": [1, 1],
            "sdz_active_studies": [2, 2],
            "sdz_active_at_least_50": [False, False],
        }
    )
    derived = SimpleNamespace(
        trait_summary=summary,
        quadrant_pairs=pd.DataFrame({"study_id": ["001"], "quadrant": ["Q4"]}),
        quadrant_summary=pd.DataFrame({"population": ["all"], "n_comparisons": [2]}),
        stool_disease_shifts=pd.DataFrame({"condition": ["NA"], "trait": [FA]}),
        stool_study_sensitivity=pd.DataFrame({"condition": ["NA"], "trait": [FA]}),
        site_disease_shifts=pd.DataFrame({"body_site": ["Feces"], "trait": [FA]}),
    )
    builder = Mock(return_value=derived)
    monkeypatch.setattr(stage, "build_tables", builder)
    monkeypatch.setattr(stage, "trait_summary", Mock(return_value=summary))
    trait_list = tmp_path / "traits.txt"
    trait_list.write_text(f"# declared order\n{OA}\n{FA}\n", encoding="utf-8")
    return SimpleNamespace(
        bundle=bundle,
        comparison=records["comparison"],
        experiments=experiments,
        trait_list=trait_list,
        loader=loader,
        builder=builder,
        derived=derived,
        producer=producer,
        digest=file_sha256(manifest),
    )


def run(parent: Any, output: Path, **kwargs: Any) -> dict[str, Any]:
    options = {"expected_sha256": parent.digest, "trait_list": parent.trait_list, **kwargs}
    return stage.run_derived(parent.bundle.manifest, parent.comparison, output, **options)


def test_exports_parent_bytes_and_binds_input_and_output_identity(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    parent = synthetic_parent(tmp_path, monkeypatch)
    output = tmp_path / "output"
    record = run(parent, output)
    assert record == json.loads((output / "delta_z_derived_run.json").read_text())
    assert record["parent_run_id"] == "parent-fixture"
    assert record["parent_manifest_sha256"] == parent.digest
    for field in (
        "source_run_ids",
        "source_commits",
        "core_run_id",
        "comparison_run_id",
        "strain_rollup",
        "scope_role",
        "contrast_population",
        "direction_basis",
    ):
        assert record[field] == parent.bundle.record[field]
    assert record["code_commit"] == "d" * 40
    assert record["code_dirty"] is False
    assert record["counts"] == {
        "traits": 2,
        "comparisons": 2,
        "studies": 2,
        "traits_with_at_least_50_active_studies": 0,
    }
    assert record["trait_list"]["trait_ids"] == [OA, FA]
    assert record["trait_list"]["trait_list_sha256"] == file_sha256(parent.trait_list)
    assert set(record["outputs"]) == {
        PAIR_FILE,
        STUDY_FILE,
        "trait_summary.tsv",
        "quadrant_pairs.tsv",
        "quadrant_summary.tsv",
        "stool_disease_shifts.tsv",
        "stool_study_sensitivity.tsv",
        "site_disease_shifts.tsv",
    }
    for name, identity in record["outputs"].items():
        assert identity == {
            "sha256": file_sha256(output / name),
            "size_bytes": (output / name).stat().st_size,
        }
    for filename in (PAIR_FILE, STUDY_FILE):
        assert (output / filename).read_bytes() == (parent.bundle.root / filename).read_bytes()
    assert not list(tmp_path.glob(".output-*"))
    parent.loader.assert_called_once_with(
        parent.bundle.manifest, parent.comparison, expected_sha256=parent.digest
    )
    args, kwargs = parent.builder.call_args
    assert args[3].iloc[0].to_dict() == {
        "study_id": "001",
        "experiment_id": "01",
        "body_site": "Feces",
        "condition": "NA",
    }
    assert list(args[4]) == [OA, FA]
    assert kwargs == {"top_diseases": 7}


def test_summary_only_uses_study_table_without_headline_or_full_analysis(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    parent = synthetic_parent(tmp_path, monkeypatch)
    parent.bundle.record["strain_rollup"] = "off"
    parent.builder.side_effect = AssertionError("Full analysis was invoked")
    record = run(parent, tmp_path / "output", summary_only=True, trait_list=None)
    assert record["summary_only"] is True
    assert record["trait_list"] is None
    assert set(record["outputs"]) == {PAIR_FILE, STUDY_FILE, "trait_summary.tsv"}


@pytest.mark.parametrize("summary_only", [False, True])
def test_fresh_nested_output_creates_parents_after_validation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, summary_only: bool
) -> None:
    """Publication accepts durable destinations with missing ancestors (#1158)."""
    parent = synthetic_parent(tmp_path, monkeypatch)
    output = tmp_path / "new" / "nested" / "output"
    assert not (tmp_path / "new").exists()
    options = {"summary_only": True, "trait_list": None} if summary_only else {}
    record = run(parent, output, **options)
    assert json.loads((output / "delta_z_derived_run.json").read_text()) == record
    assert {path.name for path in output.parent.iterdir()} == {"output"}
    for name, expected in record["outputs"].items():
        assert file_sha256(output / name) == expected["sha256"]


@pytest.mark.parametrize("failure_point", ["load", "late_source_recheck"])
def test_source_failure_leaves_missing_output_ancestors_absent(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, failure_point: str
) -> None:
    parent = synthetic_parent(tmp_path, monkeypatch)
    output = tmp_path / "new" / "nested" / "output"
    if failure_point == "load":
        parent.loader.side_effect = ValueError("source verification failed")
    else:

        def change_source(*_args: Any, **_kwargs: Any) -> Any:
            parent.experiments.write_text(parent.experiments.read_text() + "\n")
            return parent.derived

        parent.builder.side_effect = change_source
    with pytest.raises(ValueError, match=r"source|identity"):
        run(parent, output)
    assert not (tmp_path / "new").exists()
    assert not list(tmp_path.glob(".output-*"))


def test_full_run_requires_declared_trait_list_before_writing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    parent = synthetic_parent(tmp_path, monkeypatch)
    with pytest.raises(ValueError, match="trait list"):
        run(parent, tmp_path / "new" / "output", trait_list=None)
    assert not (tmp_path / "new").exists()


@pytest.mark.parametrize("kind", ["directory", "file", "dangling_symlink"])
def test_existing_destination_is_untouched_before_loading(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, kind: str
) -> None:
    parent = synthetic_parent(tmp_path, monkeypatch)
    output = tmp_path / "output"
    if kind == "directory":
        output.mkdir()
    elif kind == "file":
        output.write_text("unowned")
    else:
        output.symlink_to(tmp_path / "unowned-missing")
    original = output.lstat()
    with pytest.raises(FileExistsError):
        run(parent, output)
    parent.loader.assert_not_called()
    assert output.lstat().st_ino == original.st_ino


def test_shared_temporary_guard_precedes_loading(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    parent = synthetic_parent(tmp_path, monkeypatch)
    output = tmp_path / "temporary" / "output"
    monkeypatch.setattr(stage, "temporary_root", lambda _path: tmp_path / "temporary")
    with pytest.raises(ValueError, match=r"durable|temporary"):
        run(parent, output)
    parent.loader.assert_not_called()
    assert not output.parent.exists()


@pytest.mark.parametrize("location", ["parent", "core", "comparison", "noel"])
def test_declared_sources_are_protected_before_any_probe_or_loading(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, location: str
) -> None:
    parent = synthetic_parent(tmp_path, monkeypatch)
    root = parent.bundle.root if location == "parent" else tmp_path / "sources" / location
    output = root / "new" / "output"
    with pytest.raises(ValueError, match=r"disjoint|source"):
        run(parent, output)
    parent.loader.assert_not_called()
    assert not output.parent.exists()


@pytest.mark.parametrize("which", ["parent", "comparison"])
def test_file_alias_and_target_are_protected(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, which: str
) -> None:
    parent = synthetic_parent(tmp_path, monkeypatch)
    aliases = tmp_path / "aliases"
    aliases.mkdir()
    actual = parent.bundle.manifest if which == "parent" else parent.comparison
    alias = aliases / actual.name
    alias.symlink_to(actual)
    output = aliases / "new" / "output"
    with pytest.raises(ValueError, match=r"disjoint|source"):
        stage.run_derived(
            alias if which == "parent" else parent.bundle.manifest,
            alias if which == "comparison" else parent.comparison,
            output,
            expected_sha256=parent.digest,
            trait_list=parent.trait_list,
        )
    assert not output.parent.exists()
    parent.loader.assert_not_called()


@pytest.mark.parametrize("location", ["alias", "target"])
def test_verified_source_file_alias_parents_are_protected_before_probe(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, location: str
) -> None:
    parent = synthetic_parent(tmp_path, monkeypatch)
    target = tmp_path / "external-source" / "source.tsv"
    target.parent.mkdir()
    target.write_text("bound source bytes\n", encoding="utf-8")
    alias = tmp_path / "source-alias" / "source.tsv"
    alias.parent.mkdir()
    alias.symlink_to(target)
    parent.bundle.snapshots[alias] = {
        "sha256": file_sha256(target),
        "size_bytes": target.stat().st_size,
    }
    output = (alias.parent if location == "alias" else target.parent) / "new" / "output"
    with pytest.raises(ValueError, match=r"disjoint|source"):
        run(parent, output)
    parent.loader.assert_called_once()
    parent.builder.assert_not_called()
    assert not output.parent.exists()


@pytest.mark.parametrize("changed", ["parent", "experiments", "trait_list"])
def test_mutation_during_calculation_refuses_publication(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, changed: str
) -> None:
    parent = synthetic_parent(tmp_path, monkeypatch)
    path = {
        "parent": parent.bundle.manifest,
        "experiments": parent.experiments,
        "trait_list": parent.trait_list,
    }[changed]

    def mutate(*_args: Any, **_kwargs: Any) -> Any:
        path.write_text(path.read_text() + "\n", encoding="utf-8")
        return parent.derived

    parent.builder.side_effect = mutate
    with pytest.raises(ValueError, match=r"changed|identity"):
        run(parent, tmp_path / "output")
    assert not (tmp_path / "output").exists()
    assert not list(tmp_path.glob(".output-*"))


def test_changed_producer_refuses_publication(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    parent = synthetic_parent(tmp_path, monkeypatch)
    monkeypatch.setattr(
        stage,
        "_producer",
        Mock(side_effect=[parent.producer, {**parent.producer, "code_dirty": True}]),
    )
    with pytest.raises(ValueError, match=r"producer|changed"):
        run(parent, tmp_path / "output")
    assert not (tmp_path / "output").exists()
    assert not list(tmp_path.glob(".output-*"))


def test_mutation_while_writing_is_rechecked_and_owned_staging_is_removed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    parent = synthetic_parent(tmp_path, monkeypatch)
    original = pd.DataFrame.to_csv

    def change_source(self: pd.DataFrame, *args: Any, **kwargs: Any) -> Any:
        result = original(self, *args, **kwargs)
        parent.experiments.write_text(parent.experiments.read_text() + "\n")
        return result

    monkeypatch.setattr(pd.DataFrame, "to_csv", change_source)
    with pytest.raises(ValueError, match=r"changed|identity"):
        run(parent, tmp_path / "output")
    assert not (tmp_path / "output").exists()
    assert not list(tmp_path.glob(".output-*"))


def test_competing_destination_is_preserved_and_staging_removed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    parent = synthetic_parent(tmp_path, monkeypatch)
    output = tmp_path / "output"
    original = _publish_directory

    def race(source: Path, destination: Path) -> None:
        destination.mkdir()
        (destination / "keep.txt").write_text("unowned competitor")
        original(source, destination)

    monkeypatch.setattr(stage, "_publish_directory", race)
    with pytest.raises(FileExistsError):
        run(parent, output)
    assert (output / "keep.txt").read_text() == "unowned competitor"
    assert not list(tmp_path.glob(".output-staging-*"))


@pytest.mark.parametrize("filename", [PAIR_FILE, STUDY_FILE])
def test_copied_ledger_must_match_pinned_parent_after_transient_source_change(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, filename: str
) -> None:
    """A source restored after copying must not hide divergent published bytes (#1095)."""
    parent = synthetic_parent(tmp_path, monkeypatch)
    original_copy = shutil.copyfile
    input_path = parent.bundle.root / filename
    original_bytes = input_path.read_bytes()
    expected = parent.bundle.snapshots[input_path]

    def copy_transient(source: Path, destination: Path) -> Path:
        if source.name != filename:
            return original_copy(source, destination)
        source.write_bytes(original_bytes.replace(b"1.0000", b"9.0000"))
        try:
            return original_copy(source, destination)
        finally:
            source.write_bytes(original_bytes)

    monkeypatch.setattr(shutil, "copyfile", copy_transient)
    with pytest.raises(ValueError, match=r"[Cc]opied|ledger|parent|identity"):
        run(parent, tmp_path / "output")
    assert input_path.read_bytes() == original_bytes
    assert file_sha256(input_path) == expected["sha256"]
    assert not (tmp_path / "output").exists()
    assert not list(tmp_path.glob(".output-staging-*"))


def test_capability_failure_precedes_computation_and_leaves_no_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    parent = synthetic_parent(tmp_path, monkeypatch)
    monkeypatch.setattr(stage, "_probe_publication", Mock(side_effect=OSError("unsupported")))
    with pytest.raises(OSError, match="unsupported"):
        run(parent, tmp_path / "output")
    parent.builder.assert_not_called()
    assert not (tmp_path / "output").exists()


def test_experiments_parse_uses_the_bytes_bound_to_the_record(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """An ABA path change must never publish a body site under another payload's hash."""
    parent = synthetic_parent(tmp_path, monkeypatch)
    original = parent.experiments.read_bytes()
    original_reader = pd.read_csv

    def transient_reader(source: Any, *args: Any, **kwargs: Any) -> pd.DataFrame:
        if isinstance(source, Path) and source == parent.experiments:
            source.write_bytes(original.replace(b"Feces", b"Blood"))
            try:
                result: pd.DataFrame = original_reader(source, *args, **kwargs)
                return result
            finally:
                source.write_bytes(original)
        result = original_reader(source, *args, **kwargs)
        return result

    monkeypatch.setattr(pd, "read_csv", transient_reader)
    record = run(parent, tmp_path / "output")
    assert parent.experiments.read_bytes() == original
    assert parent.builder.call_args.args[3].loc[0, "body_site"] == "Feces"
    assert (
        record["inputs"][f"experiments.tsv@sha256:{file_sha256(parent.experiments)}"]
        == (parent.bundle.snapshots[parent.experiments])
    )


def test_experiments_rejects_changed_bytes_even_when_path_inventory_stays_original(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The returned byte buffer itself must match the verified source identity."""
    parent = synthetic_parent(tmp_path, monkeypatch)
    original_reader = Path.read_bytes

    def transient_bytes(path: Path) -> bytes:
        raw = original_reader(path)
        return raw.replace(b"Feces", b"Blood") if path == parent.experiments else raw

    monkeypatch.setattr(Path, "read_bytes", transient_bytes)
    with pytest.raises(ValueError, match=r"Experiments.*identity|input.*identity"):
        run(parent, tmp_path / "output")
    parent.builder.assert_not_called()
    assert not (tmp_path / "output").exists()


@pytest.mark.parametrize(("active_counts", "expected"), [([49, 50], 1), ([50, 51], 2)])
def test_stage_counts_the_actual_inclusive_activity_flags(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    active_counts: list[int],
    expected: int,
) -> None:
    """Exercise the real row-1 calculation and inspect its exported stage count."""
    from mtp.delta_z_derived.tables import trait_summary

    parent = synthetic_parent(tmp_path, monkeypatch)
    studies = pd.DataFrame(
        {
            trait: [1.0] * n + [0.0] * (51 - n)
            for trait, n in zip([FA, OA], active_counts, strict=True)
        },
        index=pd.Index([f"study{i}" for i in range(51)], name="study_id"),
    )
    parent.bundle.study_delta_z = studies
    parent.bundle.evidence = pd.DataFrame(
        {
            "trait": [FA, OA],
            "trait_label": ["verified FA", "verified OA"],
            "study_mean_delta_z": studies.mean().to_numpy(),
            "experiment_mean_delta_z": studies.mean().to_numpy(),
            "noel_studies_up": active_counts,
            "noel_studies_down": [0, 0],
            "noel_studies_tied": [51 - n for n in active_counts],
        }
    )
    monkeypatch.setattr(stage, "trait_summary", trait_summary)
    output = tmp_path / "output"
    record = run(parent, output, summary_only=True, trait_list=None)
    table = pd.read_csv(output / "trait_summary.tsv", sep="\t")
    assert table["sdz_active_studies"].tolist() == active_counts
    assert table["sdz_active_at_least_50"].tolist() == (
        [False, True] if expected == 1 else [True, True]
    )
    assert record["counts"]["traits_with_at_least_50_active_studies"] == expected
    assert int(table["sdz_active_at_least_50"].sum()) == expected
