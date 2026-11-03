"""Verified-core coverage uses its own population, directions and KG (#1344)."""

from __future__ import annotations

import importlib
import json
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any, cast
from unittest.mock import Mock

import pandas as pd
import pytest

from mtp.corpus_identity import CorpusIdentity, measure_matrix
from mtp.paths import file_sha256
from mtp.workflow.consumers import CoreBundle
from test_taxa_coverage_histogram import coverage
from test_taxa_coverage_histogram import inputs as inputs

TIER12 = "disease_vs_healthy_tier1_2_subject_controls"


@pytest.fixture
def core(
    inputs: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> SimpleNamespace:
    """Mock the independently tested native verifier, keep all coverage reads real."""
    module = importlib.import_module("mtp.taxa_coverage_inputs")
    manifest = tmp_path / "core" / "run.json"
    manifest.parent.mkdir()
    manifest.write_text('{"fixture": "verified boundary"}\n')
    raw = tmp_path / "raw.csv"
    raw.write_text("# BugSigDB coverage fixture, License: CC0\nStudy\nfixture\n")
    nodes = tmp_path / "nodes.tsv"
    nodes.write_text("id\tname\nNCBITaxon:2\tBacteria\n")
    identity_path = tmp_path / "source_identity.json"
    identity_path.write_text('{"fixture": "core KG identities"}\n')
    products = {
        "reported_taxa": inputs["taxa"],
        "taxa_matrix": inputs["matrix"],
        "count_profiles": inputs["profiles"],
        "experiments": inputs["experiments"],
        "raw_dump": raw,
        "kg_edges": inputs["taxonomy_edges"],
        "kg_nodes": nodes,
        "source_identity": identity_path,
    }
    run = SimpleNamespace(
        run_id="synthetic-core",
        contrast_population="all_paired",
        direction_basis="group_1_vs_group_0",
        parameters=[SimpleNamespace(name="strain_rollup", value="phenotype")],
    )
    bundle = cast(
        CoreBundle, SimpleNamespace(manifest=manifest, run=run, path=products.__getitem__)
    )
    native = Mock(return_value=bundle)
    monkeypatch.setattr(module, "open_core_bundle", native)
    identity = CorpusIdentity(
        matrix=measure_matrix(inputs["profiles"]),
        bugsigdb_vintage="fixture",
        kg_release_tag=None,
        kg_provenance_path=str(identity_path),
        kg_provenance_source="verified core execution input bytes; KG release identity unverified",
    )
    monkeypatch.setattr(
        coverage, "core_corpus_identity", Mock(return_value=identity), raising=False
    )
    preflight = Mock()
    monkeypatch.setattr(coverage, "preflight_consumer_outputs", preflight, raising=False)
    writes: list[dict[str, Any]] = []

    def record(_bundle: CoreBundle, destination: Path, **kwargs: Any) -> Path:
        assert all(path.is_file() for path in kwargs["outputs"])
        writes.append(kwargs)
        destination.write_text('{"fixture": "consumer publication boundary"}\n')
        return destination

    monkeypatch.setattr(coverage, "write_consumer_record", record, raising=False)
    return SimpleNamespace(
        module=module,
        manifest=manifest,
        bundle=bundle,
        products=products,
        native=native,
        writes=writes,
        preflight=preflight,
    )


def invoke(core: SimpleNamespace, output: Path, **kwargs: Any) -> None:
    coverage.run_coverage(
        core_run=core.manifest, strain_rollup=None, output_dir=output, tables_only=True, **kwargs
    )


def test_full_core_retains_unpaired_zero_profiles_and_binds_outputs(
    core: SimpleNamespace, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def forbidden(*_args: Any, **_kwargs: Any) -> None:
        raise AssertionError("Core coverage must not read a checkout KG provenance sidecar")

    monkeypatch.setattr(coverage, "corpus_identity", forbidden)
    output = tmp_path / "full-core"
    invoke(core, output)
    result = pd.read_csv(output / "profile_trait_mapping_summary.tsv", sep="\t")
    assert list(zip(result.study_id, result.experiment_id, result.inc_dec, strict=True)) == [
        ("s1", "e1", "inc"),
        ("s1", "e1", "dec"),
        ("s2", "e2", "inc"),
    ]
    assert result.mapped_trait_count.tolist() == [1, 0, 1]
    record = json.loads((output / "provenance.json").read_text())
    assert record["population"]["contrast_filter"] == "all_paired"
    assert record["population"]["direction_basis"] == "group_1_vs_group_0"
    assert record["strain_rollup"] == "phenotype"
    assert record["strain_rollup_evidence"]["status"] == "verified_core"
    assert record["taxonomy_diagnosis"]["matrix_kg_identity"] == "verified_core"
    assert "kg_provenance" not in record["inputs"]
    for name, product in (
        ("taxonomy_edges", "kg_edges"),
        ("kg_nodes", "kg_nodes"),
        ("source_identity", "source_identity"),
        ("profiles", "count_profiles"),
    ):
        assert record["inputs"][name]["sha256"] == file_sha256(core.products[product])
    core.native.assert_called_once_with(core.manifest)
    core.preflight.assert_called_once()
    assert len(core.writes) == 1
    assert set(core.writes[0]["additional_products"]) == {"reported_taxa", "kg_edges", "kg_nodes"}
    assert {item.name for item in core.writes[0]["outputs"]} == {
        "taxa_coverage_by_rank.tsv",
        "profile_coverage_by_rank.tsv",
        "profile_trait_mapping_summary.tsv",
        "taxa_kgmicrobe_coverage_details.tsv",
        "top_zero_trait_bugsigdb_taxa.tsv",
        "README.md",
        "provenance.json",
    }


def test_oriented_core_never_reclassifies_or_flips_again(
    core: SimpleNamespace, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from mtp import contrast_orientation

    core.bundle.run.contrast_population = TIER12
    core.bundle.run.direction_basis = "disease_vs_comparison"
    profiles = pd.read_csv(core.products["count_profiles"], sep="\t")
    profiles = profiles.loc[profiles.study_id.eq("s1")].copy()
    profiles["inc_dec"] = profiles.inc_dec.map({"inc": "dec", "dec": "inc"})
    profiles.to_csv(core.products["count_profiles"], sep="\t", index=False)
    taxa = pd.read_csv(core.products["reported_taxa"], sep="\t")
    taxa = taxa.loc[taxa.study_id.eq("s1")].copy()
    taxa["direction"] = taxa.direction.map({"increased": "decreased", "decreased": "increased"})
    taxa.to_csv(core.products["reported_taxa"], sep="\t", index=False)

    def reject_reclassification(*_args: Any, **_kwargs: Any) -> None:
        raise AssertionError("Verified-core directions are already oriented")

    monkeypatch.setattr(contrast_orientation, "read_experiment_records", reject_reclassification)
    output = tmp_path / "oriented-core"
    invoke(core, output)
    result = pd.read_csv(output / "profile_trait_mapping_summary.tsv", sep="\t")
    assert result.inc_dec.tolist() == ["dec", "inc"]
    assert result.mapped_trait_count.tolist() == [1, 0]
    record = json.loads((output / "provenance.json").read_text())
    assert record["population"]["contrast_filter"] == TIER12
    assert record["population"]["direction_basis"] == "disease_vs_comparison"
    assert (
        "No contrast-class, control-tier, pairing or direction-orientation restriction"
        not in (output / "README.md").read_text()
    )


@pytest.mark.parametrize(
    ("argument", "value"),
    [
        ("strain_rollup", "off"),
        ("contrast_population", TIER12),
        ("direction_basis", "disease_vs_comparison"),
    ],
)
def test_conflicting_core_declarations_refuse_before_output(
    core: SimpleNamespace, tmp_path: Path, argument: str, value: str
) -> None:
    output = tmp_path / "conflict"
    arguments = {"strain_rollup": None, argument: value}
    with pytest.raises(ValueError, match=r"conflict|differs"):
        coverage.run_coverage(
            core_run=core.manifest, output_dir=output, tables_only=True, **arguments
        )
    assert not output.exists()
    assert not core.writes


@pytest.mark.parametrize("argument", ["profiles", "taxonomy_edges", "kg_provenance"])
def test_foreign_core_input_refuses_before_output(
    core: SimpleNamespace, tmp_path: Path, argument: str
) -> None:
    foreign = tmp_path / "foreign.tsv"
    foreign.write_text("different input\n")
    output = tmp_path / "conflict"
    with pytest.raises(ValueError, match=r"conflict|core"):
        invoke(core, output, **{argument: foreign})
    assert not output.exists()


def test_core_verification_failure_propagates_before_data_reads(
    core: SimpleNamespace, tmp_path: Path
) -> None:
    core.native.side_effect = ValueError("native core verification failed")
    output = tmp_path / "invalid-core"
    with pytest.raises(ValueError, match="native core verification failed"):
        invoke(core, output)
    assert not output.exists()


@pytest.mark.parametrize("role", ["kg_nodes", "source_identity"])
def test_core_extra_inputs_are_rechecked_before_consumer_publication(
    core: SimpleNamespace, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, role: str
) -> None:
    original = coverage.write_readme

    def mutate(*args: Any, **kwargs: Any) -> None:
        original(*args, **kwargs)
        with core.products[role].open("a") as handle:
            handle.write("changed\n")

    monkeypatch.setattr(coverage, "write_readme", mutate)
    output = tmp_path / "changed"
    with pytest.raises(ValueError, match="Input changed"):
        invoke(core, output)
    assert not core.writes
    assert not (output / "core_consumption.json").exists()


def test_matching_core_declarations_are_accepted(core: SimpleNamespace, tmp_path: Path) -> None:
    coverage.run_coverage(
        core_run=core.manifest,
        output_dir=tmp_path / "matched",
        tables_only=True,
        strain_rollup="phenotype",
        contrast_population="all_paired",
        direction_basis="group_1_vs_group_0",
        profiles=core.products["count_profiles"],
    )
    assert len(core.writes) == 1


def test_wrong_producer_commit_refuses_before_output(core: SimpleNamespace, tmp_path: Path) -> None:
    output = tmp_path / "wrong-producer"
    with pytest.raises(ValueError, match="code commit"):
        invoke(core, output, code_commit="wrong-commit")
    assert not output.exists()


def test_core_cli_needs_no_standalone_scope_or_defaults(
    core: SimpleNamespace, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    output = tmp_path / "core-cli"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            str(coverage.SCRIPT) if hasattr(coverage, "SCRIPT") else "coverage",
            "--core-run",
            str(core.manifest),
            "--output-dir",
            str(output),
            "--tables-only",
        ],
    )
    coverage.main()
    assert (output / "core_consumption.json").is_file()


def test_explicit_as_supplied_preserves_standalone_tables(
    inputs: dict[str, Path], tmp_path: Path
) -> None:
    old, explicit = tmp_path / "default", tmp_path / "as-supplied"
    coverage.run_coverage(output_dir=old, strain_rollup="off", tables_only=True, **inputs)
    coverage.run_coverage(
        output_dir=explicit,
        strain_rollup="off",
        tables_only=True,
        contrast_population="as_supplied",
        **inputs,
    )
    for name in (
        "taxa_coverage_by_rank.tsv",
        "profile_coverage_by_rank.tsv",
        "profile_trait_mapping_summary.tsv",
        "taxa_kgmicrobe_coverage_details.tsv",
        "top_zero_trait_bugsigdb_taxa.tsv",
    ):
        assert (old / name).read_bytes() == (explicit / name).read_bytes()


def test_as_supplied_records_explicit_oriented_input_basis_without_flipping(
    inputs: dict[str, Path], tmp_path: Path
) -> None:
    output = tmp_path / "already-oriented"
    coverage.run_coverage(
        output_dir=output,
        strain_rollup="off",
        tables_only=True,
        contrast_population="as_supplied",
        direction_basis="disease_vs_comparison",
        **inputs,
    )
    result = pd.read_csv(output / "profile_trait_mapping_summary.tsv", sep="\t")
    assert result.inc_dec.tolist() == ["inc", "dec", "inc"]
    assert result.mapped_trait_count.tolist() == [1, 0, 1]
    record = json.loads((output / "provenance.json").read_text())
    assert record["population"]["direction_basis"] == "disease_vs_comparison"
    assert "Input direction basis: `disease_vs_comparison`" in (output / "README.md").read_text()


@pytest.mark.parametrize("mode", ["core", "explicit-selection"])
def test_selected_reads_preserve_numeric_and_zero_padded_identifiers(
    core: SimpleNamespace, inputs: dict[str, Path], tmp_path: Path, mode: str
) -> None:
    for name in ("taxa", "profiles", "experiments"):
        frame = pd.read_csv(inputs[name], sep="\t")
        frame["study_id"] = frame.study_id.map({"s1": "001", "s2": "2", "s3": "003", "s4": "4"})
        frame["experiment_id"] = frame.experiment_id.map(
            {"e1": "01", "e2": "2", "e3": "03", "e4": "4"}
        )
        frame.to_csv(inputs[name], sep="\t", index=False)
    output = tmp_path / mode
    if mode == "core":
        invoke(core, output)
    else:
        coverage.run_coverage(
            output_dir=output,
            strain_rollup="off",
            tables_only=True,
            contrast_population="all_paired",
            **inputs,
        )
    result = pd.read_csv(
        output / "profile_trait_mapping_summary.tsv",
        sep="\t",
        dtype={"study_id": str, "experiment_id": str},
    )
    assert list(zip(result.study_id, result.experiment_id, result.inc_dec, strict=True)) == [
        ("001", "01", "inc"),
        ("001", "01", "dec"),
        ("2", "2", "inc"),
    ]
    assert result.mapped_trait_count.tolist() == [1, 0, 1]
    assert result.total_unique_taxa.tolist() == [2, 1, 1]
