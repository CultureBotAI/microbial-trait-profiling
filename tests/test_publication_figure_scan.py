"""#1046 slice 2: immutable scan evidence and complete candidate distributions."""

from __future__ import annotations

import gzip
import hashlib
import json
from dataclasses import dataclass
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
from mtp.publication_figures.plots import validate_authored_wording

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
        "code_inventory": [],
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
    assert len(audit) == 6
    assert set(audit["trait_label"]) == {LABEL}
    assert set(audit["rank_min"]) == {1, 2, 3}
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
