"""#1046 slice 2: immutable scan evidence and complete candidate distributions."""

from __future__ import annotations

import gzip
import hashlib
import json
from dataclasses import dataclass, replace
from itertools import combinations
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from matplotlib.text import Text

from mtp.cellulose_decr import selection_scan as native
from mtp.datamodel.workflow_pydantic import ConsumerRecord
from mtp.publication_figures import scan
from mtp.publication_figures.inputs import _table
from mtp.publication_figures.plots import validate_authored_wording
from mtp.workflow.trait_policy import HISTORICAL_MATRIX_POLICIES, MATRIX_SOURCE

TRAIT = "METPO:2000103|GO:0030245"
LABEL = "capable of: cellulose catabolic process GO:0030245 [METPO:2000103]"
CONDITIONS = ["Colorectal cancer", "HIV infection", "Parkinson's disease", "Z condition"]


def identity(path: Path) -> dict[str, Any]:
    raw = path.read_bytes()
    return {"sha256": hashlib.sha256(raw).hexdigest(), "size_bytes": len(raw)}


@dataclass
class Fixture:
    root: Path
    core: Path
    metadata: Path
    summary: dict[str, Any]
    provenance: dict[str, Any]
    triples: pd.DataFrame
    traits: pd.DataFrame
    distribution: pd.DataFrame
    consumer: dict[str, Any]

    @property
    def manifest(self) -> Path:
        return self.root / native.CONSUMER_RECORD

    def seal(self) -> str:
        for name, frame in (("traits", self.traits), ("distribution", self.distribution)):
            frame.to_csv(self.root / native.OUTPUTS[name], sep="\t", index=False)
        self.root.joinpath(native.OUTPUTS["triples"]).write_bytes(
            gzip.compress(self.triples.to_csv(sep="\t", index=False).encode(), mtime=0)
        )
        self.summary["traits"] = self.traits.to_dict("records")
        for name, value in (("summary", self.summary), ("provenance", self.provenance)):
            self.root.joinpath(native.OUTPUTS[name]).write_text(json.dumps(value), encoding="utf-8")
        self.consumer["outputs"] = [
            {"path": str(path), **identity(path)}
            for path in sorted(self.root.iterdir())
            if path != self.manifest
        ]
        self.manifest.write_text(json.dumps(self.consumer), encoding="utf-8")
        return str(identity(self.manifest)["sha256"])

    def load(self) -> scan.ScanFigureInputs:
        return scan.load_scan_inputs(self.root, manifest_sha256=self.seal())


@pytest.fixture
def fixture(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Fixture:
    """Only the separately tested native source-verification boundary is stubbed."""
    root, core = tmp_path / "scan", tmp_path / "core"
    root.mkdir()
    core.mkdir()
    artifacts, inputs = [], {}
    for product in native.CORE_PRODUCTS:
        path = core / f"{product}.tsv"
        path.write_text("column_id\tcolumn_label\n" + TRAIT + "\t" + LABEL + "\n", encoding="utf-8")
        item = identity(path)
        artifacts.append(
            {
                "id": product,
                "product_id": product,
                "run_id": "core-fixture",
                "path": str(path),
                **item,
                "inventory": [],
                "validations": [],
            }
        )
        inputs[product] = {"path": path.name, "sha256": item["sha256"], "bytes": item["size_bytes"]}
    core_path = core / "run.json"
    core_record = {
        "run_id": "core-fixture",
        "manifest_path": str(core_path),
        "definition": {"id": "fixture", "version": "1", "products": [], "steps": []},
        "status": "completed",
        "started_at": "2026-10-29T00:00:00Z",
        "spec_path": "fixture.yaml",
        "spec_sha256": "a" * 64,
        "definition_sha256": "a" * 64,
        "execution_fingerprint": "a" * 64,
        "code_sha256": "a" * 64,
        "code_inventory": [
            {"path": MATRIX_SOURCE, "sha256": digest, "size_bytes": size}
            for (digest, size), policy in HISTORICAL_MATRIX_POLICIES.items()
            if policy == "raw_object_v1"
        ],
        "lock_sha256": "a" * 64,
        "git_commit": "a" * 40,
        "git_dirty": False,
        "runtime_json": "{}",
        "artifacts": artifacts,
        "parameters": [{"name": "strain_rollup", "value": "phenotype", "source": "override"}],
        "contrast_population": native.TIER1_2,
        "direction_basis": "disease_vs_comparison",
    }
    core_path.write_text(json.dumps(core_record), encoding="utf-8")
    protocol = native.load_protocol(native.PROTOCOLS[native.VERSION_V1])
    rollup = {
        "scope": "phenotype",
        "role": "declared_sensitivity",
        "declared_primary": "off",
        "declared_sensitivity": ["phenotype"],
    }
    summary: dict[str, Any] = {
        "protocol": protocol.version,
        "strain_rollup": rollup,
        "core_run": {
            "run_id": "core-fixture",
            "contrast_population": native.TIER1_2,
            "direction_basis": "disease_vs_comparison",
        },
        "lead_trait": TRAIT,
        "lead_triple": {
            "conditions": list(protocol.lead_triple),
            "eligible": dict.fromkeys(protocol.lead_triple, True),
        },
        "counts": {
            "triples": 4,
            "traits": 1,
            "conditions_with_stool_comparisons": 4,
            "eligible_conditions": 4,
            "study_disjoint_triples": 4,
            "untestable_triples": 1,
        },
        "rule": {
            "min_studies": 3,
            "top_n": 10,
            "always_included": TRAIT,
            "permutations": protocol.permutations,
            "seed": protocol.seed,
            "protocol_permutations": protocol.permutations,
            "protocol_seed": protocol.seed,
            "distance_tolerance": protocol.tolerance,
            "primary_statistic": native.MEAN_DISTANCE,
            "secondary_statistic": native.CORE,
            "secondary_rank": "u6_z",
        },
    }
    provenance = {
        "protocol": native.protocol_record(protocol.version),
        "strain_rollup": rollup,
        "core_run": {
            **summary["core_run"],
            "manifest_sha256": identity(core_path)["sha256"],
            "git_commit": "a" * 40,
            "git_dirty": False,
        },
        "code_commit": "b" * 40,
        "inputs": inputs,
        "permutations": protocol.permutations,
        "seed": protocol.seed,
        "follows_protocol": True,
    }
    rows = []
    for i, conditions in enumerate(combinations(CONDITIONS, 3)):
        tested = i < 3
        rows.append(
            {
                "trait_id": TRAIT,
                **dict(zip(("condition_a", "condition_b", "condition_c"), conditions, strict=True)),
                "n_a": 3,
                "n_b": 3,
                "n_c": 3,
                "u6_units_a": 3,
                "u6_units_b": 3,
                "u6_units_c": 3,
                "shared_studies": 0,
                "three_way_core": 1,
                "three_way_union": 5,
                native.MEAN_DISTANCE: 0.5,
                "u6_p_mean_distance": [0.05, 0.05, 0.9, ""][i],
                "u6_p_core": [0.8, 0.1, 0.8, ""][i],
                "rank_mean_distance": [1, 1, 3, ""][i],
                "rank_core": [2, 1, 2, ""][i],
                "u6_draws": 100000 if tested else 0,
                "u6_draws_without_distance": 0,
                "u6_draws_source": "shared" if tested else "none",
                "test_status": "tested" if tested else "pool_too_small",
                "is_lead_triple": i == 0,
            }
        )
    triples = pd.DataFrame(rows)
    trait: dict[str, Any] = {
        "trait_id": TRAIT,
        "trait_label": LABEL,
        "is_lead_trait": True,
        "lead_in_scan": True,
        "lead_condition_a": CONDITIONS[0],
        "lead_condition_b": CONDITIONS[1],
        "lead_condition_c": CONDITIONS[2],
        "lead_test_status": "tested",
        "lead_u6_p_mean_distance": 0.05,
        "lead_u6_p_core": 0.8,
        "lead_mean_pairwise_jaccard_distance": 0.5,
        "lead_three_way_core": 1,
    }
    histogram = []
    for statistic, bounds in ((native.MEAN_DISTANCE, (1, 2)), (native.CORE, (2, 3))):
        for reference in native.REFERENCE_SETS:
            prefix = f"{statistic}_{reference}_"
            trait.update(
                {
                    prefix + "triples": 3,
                    prefix + "rank_min": bounds[0],
                    prefix + "rank_max": bounds[1],
                    prefix + "share_at_most": bounds[1] / 3,
                }
            )
            for index in range(20):
                count = (
                    (2 if index == 0 else 1 if index == 17 else 0)
                    if statistic == native.MEAN_DISTANCE
                    else (1 if index == 1 else 2 if index == 15 else 0)
                )
                histogram.append(
                    {
                        "trait_id": TRAIT,
                        "statistic": statistic,
                        "reference_set": reference,
                        "p_above": index / 20,
                        "p_at_most": (index + 1) / 20,
                        "triples": count,
                    }
                )
    pd.DataFrame(
        {
            "condition": CONDITIONS,
            "stool_comparisons": [3] * 4,
            "stool_studies": [3] * 4,
            "u6_units": [3] * 4,
            "eligible": [True] * 4,
        }
    ).to_csv(root / native.OUTPUTS["conditions"], sep="\t", index=False)
    coverage = root / "core_consumption.coverage.tsv"
    coverage.write_text("fixture\n", encoding="utf-8")
    consumer = {
        "id": "scan-fixture",
        "run_id": "core-fixture",
        "run_manifest_path": str(core_path),
        "run_manifest_sha256": identity(core_path)["sha256"],
        "consumer": native.CONSUMER,
        "product_ids": list(native.CORE_PRODUCTS),
        "artifact_ids": list(native.CORE_PRODUCTS),
        "feature_policy": "exclude_mediadive",
        "cohort": "Homo sapiens",
        "created_at": "2026-10-29T00:00:00Z",
        "coverage_path": str(coverage),
        "feature_filter_id": "sha256:" + "a" * 64,
        "coverage_scope": "input_cohort_before_analysis_specific_exclusions",
        "outputs": [],
    }
    monkeypatch.setattr(
        scan,
        "verify_consumer_record",
        lambda path: ConsumerRecord.model_validate_json(path.read_bytes()),
    )
    return Fixture(
        root,
        core_path,
        core / "feature_metadata.tsv",
        summary,
        provenance,
        triples,
        pd.DataFrame([trait]),
        pd.DataFrame(histogram),
        consumer,
    )


def test_historical_v1_preserves_role_labels_and_every_candidate(fixture: Fixture) -> None:
    bundle = fixture.load()
    assert bundle.context["protocol_version"] == native.VERSION_V1
    assert bundle.context["scope_role"] == "declared_sensitivity"
    assert bundle.context["historical"] is True
    assert bundle.context["report_beside"] == []
    audit = scan.scan_audit(bundle)
    assert len(audit) == 8
    assert set(audit["trait_label"]) == {LABEL}
    assert set(audit["rank_min"].dropna()) == {1, 2, 3}
    pages = list(scan.scan_pages(bundle))
    assert len(pages) == 1
    assert pages[0].mark_count == 6
    assert len(pages[0].row_ids) == len(set(pages[0].row_ids)) == 6
    assert pages[0].details["candidate_triples"] == 4
    assert pages[0].details["tested_triples"] == 3
    assert pages[0].details["untestable_triples"] == 1
    authored = " ".join(item.get_text() for item in pages[0].figure.findobj(Text))
    assert "Candidate-triple" in authored
    assert "nominal" in authored
    assert "Rank 1\u20132 of 3" in authored
    assert "Rank 2\u20133 of 3" in authored
    assert "permutation-null histogram" not in authored


def test_v2_requires_explicit_companion_context(fixture: Fixture) -> None:
    v2_fixture(fixture)
    bundle = fixture.load()
    assert bundle.context["report_beside"] == ["off"]
    assert bundle.context["scope_role"] == "primary"
    assert bundle.context["historical"] is False


@pytest.mark.parametrize(
    "change",
    [
        "missing_row",
        "duplicate",
        "lead",
        "rank",
        "negative_support",
        "p",
        "bin",
        "edge",
        "label",
        "core",
        "protocol",
        "count",
    ],
)
def test_rebound_scientific_inconsistency_is_refused(fixture: Fixture, change: str) -> None:
    fixture.load()
    if change == "missing_row":
        fixture.triples = fixture.triples.iloc[:-1]
    elif change == "duplicate":
        fixture.triples.iloc[1] = fixture.triples.iloc[0]
    elif change == "lead":
        fixture.triples.loc[1, "is_lead_triple"] = True
    elif change == "rank":
        fixture.triples.loc[0, "rank_mean_distance"] = 2
    elif change == "negative_support":
        fixture.triples.loc[0, "n_a"] = -1
    elif change == "p":
        fixture.triples.loc[0, "u6_p_core"] = 1.1
    elif change == "bin":
        fixture.distribution.loc[0, "triples"] = 1
    elif change == "edge":
        fixture.distribution.loc[0, "p_at_most"] = 0.04
    elif change == "label":
        fixture.traits.loc[0, "trait_label"] = "Unverified physiology"
    elif change == "core":
        fixture.provenance["core_run"]["git_commit"] = "c" * 40
    elif change == "protocol":
        fixture.provenance["protocol"]["sha256"] = "f" * 64
    else:
        fixture.summary["counts"]["triples"] = 3
    with pytest.raises(ValueError, match=r"(?i)scan|candidate|callout|recorded|protocol"):
        fixture.load()


def test_wrong_manifest_pin_and_native_failure_are_refused(
    fixture: Fixture, monkeypatch: pytest.MonkeyPatch
) -> None:
    fixture.seal()
    with pytest.raises(ValueError, match=r"SHA256|pin"):
        scan.load_scan_inputs(fixture.root, manifest_sha256="f" * 64)

    def refuse(_path: Path) -> ConsumerRecord:
        raise ValueError("native source invalid")

    monkeypatch.setattr(scan, "verify_consumer_record", refuse)
    with pytest.raises(ValueError, match="native source invalid"):
        fixture.load()


@pytest.mark.parametrize("name", ["summary", "triples", "provenance"])
def test_changed_after_native_verification_fails(
    fixture: Fixture, monkeypatch: pytest.MonkeyPatch, name: str
) -> None:
    fixture.seal()

    def mutate(path: Path) -> ConsumerRecord:
        record = ConsumerRecord.model_validate_json(path.read_bytes())
        target = fixture.root / native.OUTPUTS[name]
        target.write_bytes(target.read_bytes() + b" ")
        return record

    monkeypatch.setattr(scan, "verify_consumer_record", mutate)
    with pytest.raises(ValueError, match=r"hash|size|identity|changed"):
        scan.load_scan_inputs(
            fixture.root, manifest_sha256=str(identity(fixture.manifest)["sha256"])
        )


def test_recheck_detects_bound_metadata_change(fixture: Fixture, tmp_path: Path) -> None:
    alias = tmp_path / "scan-alias"
    alias.symlink_to(fixture.root, target_is_directory=True)
    pin = fixture.seal()
    bundle = scan.load_scan_inputs(alias, manifest_sha256=pin)
    scan.recheck_scan_inputs(bundle)
    fixture.metadata.write_bytes(fixture.metadata.read_bytes() + b"\n")
    with pytest.raises(ValueError, match=r"hash|changed"):
        scan.recheck_scan_inputs(bundle)


def test_authored_scan_words_pass_existing_guard(
    fixture: Fixture, monkeypatch: pytest.MonkeyPatch
) -> None:
    bundle = fixture.load()
    original = validate_authored_wording
    called: list[str] = []

    def checked(text: str) -> str:
        called.append(text)
        return original(text)

    monkeypatch.setattr(scan, "validate_authored_wording", checked)
    list(scan.scan_pages(bundle))
    assert any("Candidate-triple" in text for text in called)
    assert any("nominal" in text for text in called)


def test_recheck_refuses_identical_byte_alias_retarget(fixture: Fixture, tmp_path: Path) -> None:
    import shutil

    alias = tmp_path / "logical-scan"
    alias.symlink_to(fixture.root, target_is_directory=True)
    pin = fixture.seal()
    bundle = scan.load_scan_inputs(alias, manifest_sha256=pin)
    other = tmp_path / "same-bytes"
    shutil.copytree(fixture.root, other)
    alias.unlink()
    alias.symlink_to(other, target_is_directory=True)
    with pytest.raises(ValueError, match="alias target changed"):
        scan.recheck_scan_inputs(bundle)


def test_reader_uses_same_hashed_bytes_for_provenance(
    fixture: Fixture, monkeypatch: pytest.MonkeyPatch
) -> None:
    pin = fixture.seal()
    path = fixture.root / native.OUTPUTS["provenance"]
    original = path.read_bytes()
    changed_record = json.loads(original)
    changed_record["code_commit"] = "c" * 40
    changed = json.dumps(changed_record).encode()
    read = Path.read_bytes
    calls = 0

    def alternating(current: Path) -> bytes:
        nonlocal calls
        if current == path:
            calls += 1
            return changed if calls == 2 else original
        return read(current)

    monkeypatch.setattr(Path, "read_bytes", alternating)
    # A second identity read is allowed, but it must fail closed; a second unbound
    # parse may never silently replace the producer attached to the hashed bytes.
    with pytest.raises(ValueError, match=r"hash|size|identity|changed"):
        scan.load_scan_inputs(fixture.root, manifest_sha256=pin)


def test_plot_contains_each_candidate_twice_without_thinning(fixture: Fixture) -> None:
    bundle = fixture.load()
    page = next(scan.scan_pages(bundle))
    assert len(page.figure.axes) == 2
    for axis in page.figure.axes:
        rugs = [line for line in axis.lines if line.get_marker() == "|"]
        assert len(rugs) == 1
        assert np.asarray(rugs[0].get_xdata()).size == 3
    assert sum(page.details["panels"][0]["bin_counts"]) == 3
    assert sum(page.details["panels"][1]["bin_counts"]) == 3


def v2_fixture(fixture: Fixture) -> None:
    """Populate the native v2 comparability declarations, without invoking the scan."""
    protocol = native.load_protocol()
    rollup = native.scope_record("phenotype", native.RECORDED, protocol)
    fixture.summary.update(protocol=protocol.version, strain_rollup=rollup)
    fixture.provenance.update(
        protocol=native.protocol_record(protocol.version),
        strain_rollup=rollup,
        environment={"python": "fixture", "numpy": "fixture"},
        methods={"path": "methods/cellulose_decr_selection_scan.md", "sha256": "d" * 64},
        chunk=32,
        scripts={"src/mtp/cellulose_decr/selection_scan.py": "e" * 64},
    )
    fixture.provenance["core_run"].update(
        parameters={"strain_rollup": "phenotype", "source_mode": "pinned"},
        sources={"raw_dump": "f" * 64},
    )


def test_native_pairing_identity_and_exact_core_pin_are_exposed(fixture: Fixture) -> None:
    v2_fixture(fixture)
    bundle = fixture.load()
    assert bundle.core_manifest_sha256 == identity(fixture.core)["sha256"]
    assert bundle.native_run.core_run_id == "core-fixture"
    assert bundle.native_run.protocol.version == native.VERSION
    assert bundle.pairing_identity == {
        "protocol": fixture.provenance["protocol"],
        "code_commit": "b" * 40,
        "scripts": fixture.provenance["scripts"],
        "environment": fixture.provenance["environment"],
        "methods": fixture.provenance["methods"],
        "chunk": 32,
        "rule": fixture.summary["rule"],
        "core population": (native.TIER1_2, "disease_vs_comparison"),
        "core git_commit": ("a" * 40, False),
        "core parameters": {"source_mode": "pinned"},
        "core sources": {"raw_dump": "f" * 64},
    }


def test_historical_scan_needs_no_v2_comparability_declarations(fixture: Fixture) -> None:
    bundle = fixture.load()
    assert bundle.pairing_identity is None
    assert bundle.core_manifest_sha256 == fixture.consumer["run_manifest_sha256"]


@pytest.mark.parametrize("field", ["environment", "methods", "chunk"])
def test_missing_native_pairing_metadata_refuses_v2(fixture: Fixture, field: str) -> None:
    v2_fixture(fixture)
    fixture.provenance.pop(field)
    with pytest.raises(ValueError, match=f"missing {field}"):
        fixture.load()


def test_pairing_uses_native_contract_on_already_loaded_records(
    fixture: Fixture, monkeypatch: pytest.MonkeyPatch
) -> None:
    from mtp.cellulose_decr import selection_scan_scopes

    v2_fixture(fixture)
    original = selection_scan_scopes.comparable
    calls: list[selection_scan_scopes.ScanRun] = []

    def capture(run: selection_scan_scopes.ScanRun) -> dict[str, Any]:
        calls.append(run)
        return original(run)

    monkeypatch.setattr(selection_scan_scopes, "comparable", capture)
    fixture.load()
    assert len(calls) == 1
    assert calls[0].core_run_id == "core-fixture"
    assert calls[0].summary == fixture.summary
    assert calls[0].provenance == fixture.provenance


def native_consumer_products(fixture: Fixture) -> None:
    """Match the producer's default nine references, excluding scan-only reported_taxa.

    The native verifier boundary is still stubbed by the synthetic fixture. The
    additional default products are real core artifact records here, while all
    five scientific scan input bindings remain unchanged.
    """
    core = json.loads(fixture.core.read_bytes())
    defaults = [
        "count_profiles",
        "count_coverage",
        "feature_catalog",
        "feature_metadata",
        "experiments",
        "taxa_matrix",
        "raw_dump",
        "studies",
        "source_identity",
    ]
    existing = {item["product_id"] for item in core["artifacts"]}
    for name in defaults:
        if name not in existing:
            path = fixture.core.parent / f"{name}.tsv"
            path.write_text("synthetic bound source\n", encoding="utf-8")
            core["artifacts"].append(
                {
                    "id": name,
                    "product_id": name,
                    "run_id": "core-fixture",
                    "path": str(path),
                    **identity(path),
                    "inventory": [],
                    "validations": [],
                }
            )
    fixture.consumer["product_ids"] = defaults
    fixture.consumer["artifact_ids"] = defaults.copy()
    rebind_core(fixture, core)


def rebind_core(fixture: Fixture, core: dict[str, Any]) -> None:
    fixture.core.write_text(json.dumps(core), encoding="utf-8")
    digest = identity(fixture.core)["sha256"]
    fixture.consumer["run_manifest_sha256"] = digest
    fixture.provenance["core_run"]["manifest_sha256"] = digest


@pytest.mark.parametrize("version", [native.VERSION_V1, native.VERSION])
def test_native_default_consumer_binds_reported_taxa_through_verified_core(
    fixture: Fixture, version: str
) -> None:
    if version == native.VERSION:
        v2_fixture(fixture)
    native_consumer_products(fixture)
    assert "reported_taxa" not in fixture.consumer["product_ids"]
    bundle = fixture.load()
    source = fixture.core.parent / fixture.provenance["inputs"]["reported_taxa"]["path"]
    assert bundle.snapshots[source] == identity(source)
    assert bundle.targets[source] == source.resolve()
    assert source.parent in bundle.source_roots
    assert len(scan.scan_audit(bundle)) == 8
    assert bundle.context["protocol_version"] == version


@pytest.mark.parametrize(
    "change",
    ["missing", "duplicate", "wrong_path", "wrong_sha", "wrong_size", "explicit_reference"],
)
def test_scan_only_core_product_still_requires_unique_matching_bound_artifact(
    fixture: Fixture, change: str
) -> None:
    native_consumer_products(fixture)
    core = json.loads(fixture.core.read_bytes())
    item = next(row for row in core["artifacts"] if row["product_id"] == "reported_taxa")
    if change == "missing":
        core["artifacts"].remove(item)
    elif change == "duplicate":
        core["artifacts"].append({**item, "id": "another-reported-taxa"})
    elif change == "wrong_path":
        fixture.provenance["inputs"]["reported_taxa"]["path"] = "experiments.tsv"
    elif change == "wrong_sha":
        fixture.provenance["inputs"]["reported_taxa"]["sha256"] = "f" * 64
    elif change == "wrong_size":
        fixture.provenance["inputs"]["reported_taxa"]["bytes"] += 1
    else:
        fixture.consumer["product_ids"].append("reported_taxa")
        fixture.consumer["artifact_ids"].append("experiments")
    rebind_core(fixture, core)
    with pytest.raises(ValueError, match=r"Scan (core product|input|explicit product)"):
        fixture.load()


def test_late_scan_only_source_mutation_is_refused(fixture: Fixture) -> None:
    native_consumer_products(fixture)
    bundle = fixture.load()
    source = fixture.core.parent / "reported_taxa.tsv"
    source.write_bytes(source.read_bytes() + b"changed\n")
    with pytest.raises(ValueError, match="Scan input hash changed"):
        scan.recheck_scan_inputs(bundle)


def set_metadata(fixture: Fixture, object_name: str) -> None:
    frame = pd.read_csv(fixture.metadata, sep="\t", keep_default_na=False)
    frame["object_name"] = object_name
    frame.to_csv(fixture.metadata, sep="\t", index=False)
    core = json.loads(fixture.core.read_bytes())
    item = next(row for row in core["artifacts"] if row["product_id"] == "feature_metadata")
    item.update(identity(fixture.metadata))
    fixture.provenance["inputs"]["feature_metadata"].update(
        sha256=item["sha256"], bytes=item["size_bytes"]
    )
    rebind_core(fixture, core)


@pytest.mark.parametrize("object_name", ["cellulose catabolic process", ""])
def test_display_label_uses_object_name_with_native_column_label_check(
    fixture: Fixture, object_name: str
) -> None:
    set_metadata(fixture, object_name)
    bundle = fixture.load()
    expected = object_name or LABEL
    assert bundle.labels[TRAIT] == expected
    assert set(scan.scan_audit(bundle)["trait_label"]) == {expected}
    text = " ".join(item.get_text() for item in next(scan.scan_pages(bundle)).figure.findobj(Text))
    assert expected in text
    if object_name:
        assert LABEL not in text
    fixture.traits.loc[0, "trait_label"] = "A different native label"
    with pytest.raises(ValueError, match="label differs from bound metadata"):
        fixture.load()


@pytest.mark.parametrize("status", [native.NO_TAXA, native.POOL_TOO_SMALL])
def test_audit_keeps_untestable_candidates_without_plotted_marks(
    fixture: Fixture, status: str
) -> None:
    fixture.triples.loc[3, "test_status"] = status
    if status == native.NO_TAXA:
        fixture.triples.loc[3, ["n_a", "n_b", "n_c", native.CORE, "three_way_union"]] = 0
    bundle = fixture.load()
    audit = scan.scan_audit(bundle)
    assert len(audit) == 8  # Four triples, two statistics, including untested rows.
    untested = audit.loc[audit["test_status"].eq(status)]
    assert len(untested) == 2
    assert not untested["plotted"].any()
    assert untested["p_value"].isna().all()
    assert untested["rank_min"].isna().all()
    page = next(scan.scan_pages(bundle))
    assert page.mark_count == 6
    assert set(page.row_ids) == set(audit.loc[audit["plotted"], "row_id"])
    assert not set(untested["row_id"]) & set(page.row_ids)
    assert page.details["candidate_triples"] == 4
    assert page.details["tested_triples"] == 3
    assert page.details["untestable_triples"] == 1


def test_row_ids_are_unique_across_scope_pages(fixture: Fixture) -> None:
    phenotype = fixture.load()
    off = replace(phenotype, context={**phenotype.context, "strain_rollup": "off"})
    audits = pd.concat([scan.scan_audit(phenotype), scan.scan_audit(off)], ignore_index=True)
    assert not audits["row_id"].duplicated().any()
    ids = [
        row_id for bundle in (phenotype, off) for row_id in next(scan.scan_pages(bundle)).row_ids
    ]
    assert len(ids) == len(set(ids)) == 12


@pytest.mark.parametrize("policy", ["raw_object_v1", "facultative_union_v1"])
def test_scan_resolves_and_displays_bound_policy_once(
    fixture: Fixture, monkeypatch: pytest.MonkeyPatch, policy: str
) -> None:
    from mtp.workflow.trait_policy import resolve_trait_policy

    core = json.loads(fixture.core.read_bytes())
    digest, size = next(key for key, value in HISTORICAL_MATRIX_POLICIES.items() if value == policy)
    core["code_inventory"] = [{"path": MATRIX_SOURCE, "sha256": digest, "size_bytes": size}]
    rebind_core(fixture, core)
    calls: list[str] = []

    def resolved(run: Any) -> Any:
        calls.append(run.run_id)
        return resolve_trait_policy(run)

    monkeypatch.setattr(scan, "resolve_trait_policy", resolved, raising=False)
    bundle = fixture.load()
    assert calls == ["core-fixture"]
    assert bundle.context["trait_key_policy"] == policy
    page = next(scan.scan_pages(bundle))
    assert page.details["trait_key_policy"] == policy
    text = " ".join(item.get_text() for item in page.figure.findobj(Text))
    assert f"trait keys: {policy}" in text


def test_scan_refuses_unknown_historical_policy(fixture: Fixture) -> None:
    core = json.loads(fixture.core.read_bytes())
    core["code_inventory"][0]["sha256"] = "0" * 64
    rebind_core(fixture, core)
    with pytest.raises(ValueError, match="Unknown historical matrix implementation"):
        fixture.load()


def test_scan_explains_rank_null_and_readable_recorded_role(fixture: Fixture) -> None:
    page = next(scan.scan_pages(fixture.load()))
    text = " ".join(" ".join(item.get_text() for item in page.figure.findobj(Text)).split())
    assert "rank 1 = smallest nominal p-value" in text
    assert "U6: whole-study-signature null" in text
    assert "distinct studies" in text
    assert "declared sensitivity" in text
    assert "declared_sensitivity" not in text


def test_callout_conditions_come_from_record_not_fixed_abbreviations(fixture: Fixture) -> None:
    bundle = fixture.load()
    summary = {**bundle.summary, "lead_triple": {"conditions": ["Alpha", "Beta", "Gamma"]}}
    page = next(scan.scan_pages(replace(bundle, summary=summary)))
    text = " ".join(item.get_text() for item in page.figure.findobj(Text))
    assert "Alpha / Beta / Gamma" in text
    assert "PD / HIV / CRC" not in text
    assert "Parkinson's disease" not in text


def test_numeric_callout_equality_does_not_depend_on_text_spelling(fixture: Fixture) -> None:
    fixture.triples[native.MEAN_DISTANCE] = fixture.triples[native.MEAN_DISTANCE].astype(object)
    fixture.triples[native.CORE] = fixture.triples[native.CORE].astype(object)
    fixture.triples.loc[0, native.MEAN_DISTANCE] = "1"
    fixture.traits.loc[0, "lead_" + native.MEAN_DISTANCE] = 1.0
    fixture.triples.loc[0, native.CORE] = "1.0"
    fixture.traits.loc[0, "lead_" + native.CORE] = 1
    fixture.load()


def test_late_focal_rows_have_consistent_whole_table_parsing(fixture: Fixture) -> None:
    """#1240: blanks after the first 32,768 rows must not create mixed Python types."""
    bundle = fixture.load()
    before = pd.concat([fixture.triples.iloc[[0]]] * 33000, ignore_index=True)
    before["trait_id"] = "other-trait"
    frame = pd.concat([before, fixture.triples], ignore_index=True)
    # Force the native-size parser chunk boundary (at least 33 columns).
    for index in range(12):
        frame[f"extra_{index}"] = 0
    parsed = _table(frame.to_csv(sep="\t", index=False).encode())
    column = parsed["u6_p_mean_distance"]
    assert {type(value) for value in column} == {str}
    focal = parsed.loc[parsed["trait_id"].eq(TRAIT)].copy()
    assert len(focal) == 4
    assert len(parsed) > 32768
    reloaded = replace(bundle, tables={**bundle.tables, "triples": focal})
    scan._validate(reloaded, native.load_protocol(native.PROTOCOLS[native.VERSION_V1]))
    assert scan.scan_audit(reloaded)["p_value"].dropna().tolist() == [
        0.05,
        0.05,
        0.9,
        0.8,
        0.1,
        0.8,
    ]
