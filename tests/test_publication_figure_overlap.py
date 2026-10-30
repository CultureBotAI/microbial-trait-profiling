"""Declared carrier-overlap figures bind saved evidence without selecting pairs (#1046)."""

from __future__ import annotations

import hashlib
import io
import json
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Any
from unittest.mock import Mock

import pandas as pd
import pytest
from git import Actor, Repo
from matplotlib.text import Text

from mtp.datamodel.workflow_pydantic import ConsumerRecord, WorkflowRun
from mtp.publication_figures import overlap
from mtp.workflow.trait_policy import HISTORICAL_MATRIX_POLICIES, MATRIX_SOURCE

SCRIPTS = (
    "src/mtp/matrix_overlap_coverage.py",
    "src/mtp/trait_lists.py",
    "src/mtp/matrix_overlap_inputs.py",
    "src/mtp/human_host_corpus.py",
    "src/mtp/workflow/coverage.py",
    "src/mtp/workflow/output_roots.py",
)
METHOD = "methods/matrix_coverage_and_pd_overlap.md"
TABLE = "matrix_trait_overlap.tsv"
FILES = (
    TABLE,
    "tier12_coverage_by_rank.tsv",
    "filter_cascade.tsv",
    "parkinsons_stool_study_overlap_summary.tsv",
    "parkinsons_stool_study_pairs.tsv",
    "summary.json",
    "provenance.json",
    "core_consumption.coverage.tsv",
)


def identity(path: Path) -> dict[str, Any]:
    raw = path.read_bytes()
    return {"sha256": hashlib.sha256(raw).hexdigest(), "size_bytes": len(raw)}


def historical_implementation(policy: str) -> dict[str, Any]:
    digest, size = next(key for key, value in HISTORICAL_MATRIX_POLICIES.items() if value == policy)
    return {"path": MATRIX_SOURCE, "sha256": digest, "size_bytes": size}


@dataclass
class Fixture:
    root: Path
    traits: Path
    metadata: Path
    core: Path
    provenance: dict[str, Any]
    summary: dict[str, Any]
    consumer: dict[str, Any]
    frame: pd.DataFrame
    verifier: Mock

    @property
    def manifest(self) -> Path:
        return self.root / "core_consumption.json"

    def seal(self) -> str:
        self.frame.to_csv(self.root / TABLE, sep="\t", index=False, float_format="%.12g")
        (self.root / "summary.json").write_text(json.dumps(self.summary), encoding="utf-8")
        self.provenance["outputs"] = {
            name: {"relative_path": name, "sha256": identity(self.root / name)["sha256"]}
            for name in FILES
            if name not in {"provenance.json", "core_consumption.coverage.tsv"}
        }
        (self.root / "provenance.json").write_text(json.dumps(self.provenance), encoding="utf-8")
        self.consumer["outputs"] = [
            {"path": str(self.root / name), **identity(self.root / name)} for name in FILES
        ]
        self.manifest.write_text(json.dumps(self.consumer), encoding="utf-8")
        return str(identity(self.manifest)["sha256"])

    def load(self) -> overlap.OverlapInputs:
        return overlap.load_overlap_inputs(
            self.manifest,
            manifest_sha256=str(identity(self.manifest)["sha256"]),
            trait_list=self.traits,
        )


@pytest.fixture
def fixture(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Fixture:
    root, core_root = tmp_path / "overlap", tmp_path / "core"
    root.mkdir()
    core_root.mkdir()
    ids = [f"METPO:{1000000 + index}" for index in range(9)]
    labels = [f"Verified annotation {index}" for index in range(9)]
    traits = tmp_path / "declared.txt"
    traits.write_text("# ordered nine\n" + "\n".join(ids) + "\n", encoding="utf-8")
    metadata = core_root / "taxa_trait_matrix_columns.tsv"
    pd.DataFrame({"column_id": ids, "object_name": labels}).to_csv(metadata, sep="\t", index=False)
    core_record: dict[str, Any] = {
        "run_id": "core-fixture",
        "status": "completed",
        "git_commit": "a" * 40,
        "git_dirty": False,
        "definition": {"id": "bugsigdb-core-subset", "version": "1", "products": [], "steps": []},
        "started_at": "2026-10-29T00:00:00Z",
        "manifest_path": str(core_root / "run.json"),
        "spec_path": "synthetic.yaml",
        "spec_sha256": "b" * 64,
        "definition_sha256": "b" * 64,
        "execution_fingerprint": "b" * 64,
        "code_sha256": "b" * 64,
        "code_inventory": [historical_implementation("raw_object_v1")],
        "lock_sha256": "b" * 64,
        "runtime_json": "{}",
        "contrast_population": "disease_vs_healthy_tier1_2_subject_controls",
        "direction_basis": "disease_vs_comparison",
        "parameters": [{"name": "strain_rollup", "value": "phenotype", "source": "override"}],
        "artifacts": [
            {
                "id": "feature_metadata",
                "product_id": "feature_metadata",
                "run_id": "core-fixture",
                "path": str(metadata),
                "inventory": [],
                "validations": [],
                **identity(metadata),
            }
        ],
    }
    WorkflowRun.model_validate(core_record)
    core = core_root / "run.json"
    core.write_text(json.dumps(core_record), encoding="utf-8")
    declared = {
        "trait_list_file": traits.name,
        "trait_list_sha256": identity(traits)["sha256"],
        "trait_count": 9,
        "trait_ids": ids,
        "source_name": traits.name,
        "sha256": identity(traits)["sha256"],
    }
    common = {
        "schema_version": "1.0.0",
        "consumer": "matrix-overlap-coverage",
        "core_run_id": core_record["run_id"],
        "workflow_id": core_record["definition"]["id"],
        "strain_rollup": "phenotype",
        "contrast_population": core_record["contrast_population"],
    }
    repo_root = tmp_path / "producer"
    with Repo.init(repo_root) as repository:
        for name in (*SCRIPTS, METHOD):
            path = repo_root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(f"Synthetic historical file {name}\n", encoding="utf-8")
        repository.index.add([*SCRIPTS, METHOD])
        actor = Actor("Fixture", "fixture@example.invalid")
        commit = repository.index.commit("Synthetic producer", author=actor, committer=actor)
    monkeypatch.setattr(overlap, "__file__", str(repo_root / "src/mtp/overlap.py"))
    provenance = {
        **common,
        "code_commit": commit.hexsha,
        "core_git_commit_sha": core_record["git_commit"],
        "declared_traits": declared.copy(),
        "scripts": {name: {"sha256": identity(repo_root / name)["sha256"]} for name in SCRIPTS},
        "methods": {METHOD: {"sha256": identity(repo_root / METHOD)["sha256"]}},
        "script_sha256": identity(repo_root / SCRIPTS[0])["sha256"],
        "methods_sha256": identity(repo_root / METHOD)["sha256"],
        "consumed_products": {
            "feature_metadata": {
                "relative_path": metadata.name,
                "sha256": identity(metadata)["sha256"],
            }
        },
    }
    summary = {
        **common,
        "declared_traits": {**declared, "pair_count": 36},
        "coverage_and_cascade": {"matrix_taxa_total": 10},
    }
    carriers = [{1, 2}, {2, 3}, set(), {4}, {1, 2, 3, 4}, {1}, {3}, set(), {2, 4}]
    rows = []
    for a, b in combinations(range(9), 2):
        intersection, union = len(carriers[a] & carriers[b]), len(carriers[a] | carriers[b])
        fraction = intersection / union if union else None
        rows.append(
            {
                "trait_a_id": ids[a],
                "trait_a_label": labels[a],
                "trait_b_id": ids[b],
                "trait_b_label": labels[b],
                "matrix_taxa_total": 10,
                "trait_a_taxa": len(carriers[a]),
                "trait_b_taxa": len(carriers[b]),
                "intersection_taxa": intersection,
                "union_taxa": union,
                "jaccard_fraction": fraction,
                "jaccard_pct": 100 * fraction if fraction is not None else None,
            }
        )
    for name in FILES:
        (root / name).write_text("field\nvalue\n", encoding="utf-8")
    consumer = {
        "id": "consumer:fixture",
        "run_id": "core-fixture",
        "run_manifest_path": str(core),
        "run_manifest_sha256": identity(core)["sha256"],
        "consumer": "matrix-overlap-coverage",
        "product_ids": ["feature_metadata"],
        "artifact_ids": ["feature_metadata"],
        "feature_policy": "exclude_mediadive",
        "cohort": "all",
        "created_at": "2026-10-29T00:00:00Z",
        "coverage_path": str(root / "core_consumption.coverage.tsv"),
        "feature_filter_id": "sha256:" + "b" * 64,
        "coverage_scope": "input_cohort_before_analysis_specific_exclusions",
    }
    verifier = Mock(side_effect=lambda path: ConsumerRecord.model_validate_json(path.read_bytes()))
    monkeypatch.setattr(overlap, "verify_consumer_record", verifier)
    result = Fixture(
        root, traits, metadata, core, provenance, summary, consumer, pd.DataFrame(rows), verifier
    )
    result.seal()
    return result


def test_complete_saved_pairs_labels_context_and_undefined_audit(fixture: Fixture) -> None:
    bundle = fixture.load()
    fixture.verifier.assert_called_once_with(fixture.manifest)
    audit = overlap.overlap_audit(bundle)
    assert len(audit) == 36
    assert (
        audit[["trait_a_id", "trait_b_id"]].values.tolist()
        == fixture.frame[["trait_a_id", "trait_b_id"]].values.tolist()
    )
    assert (audit.overlap_status == "undefined_empty_union").sum() == 1
    assert ((audit.jaccard_fraction == 0) & (audit.overlap_status == "available")).any()
    assert bundle.context == {
        "core_run_id": "core-fixture",
        "strain_rollup": "phenotype",
        "contrast_population": "disease_vs_healthy_tier1_2_subject_controls",
        "direction_basis": "disease_vs_comparison",
        "source_commit": "a" * 40,
        "consumer_commit": fixture.provenance["code_commit"],
        "trait_key_policy": "raw_object_v1",
    }
    assert bundle.record == fixture.provenance
    assert {fixture.root, fixture.core.parent}.issubset(bundle.source_roots)


@pytest.mark.parametrize("policy", ["raw_object_v1", "facultative_union_v1"])
def test_historical_policy_uses_audited_producer_identity(fixture: Fixture, policy: str) -> None:
    record = json.loads(fixture.core.read_bytes())
    assert "trait_key_policy" not in record
    record["code_inventory"] = [historical_implementation(policy)]
    fixture.core.write_text(json.dumps(record))
    fixture.consumer["run_manifest_sha256"] = identity(fixture.core)["sha256"]
    fixture.seal()
    assert fixture.load().context["trait_key_policy"] == policy


def test_unknown_historical_policy_is_not_inferred_from_labels(fixture: Fixture) -> None:
    record = json.loads(fixture.core.read_bytes())
    record["code_inventory"][0]["sha256"] = "c" * 64
    fixture.core.write_text(json.dumps(record))
    fixture.consumer["run_manifest_sha256"] = identity(fixture.core)["sha256"]
    fixture.seal()
    with pytest.raises(ValueError, match="Unknown historical matrix implementation"):
        fixture.load()


@pytest.mark.parametrize(
    "corruption",
    [
        "missing_pair",
        "duplicate_pair",
        "reverse_pair",
        "label",
        "negative_count",
        "noninteger_count",
        "union",
        "fraction",
        "percent",
        "nan",
        "empty_union_value",
    ],
)
def test_rebound_saved_table_contract_is_checked(fixture: Fixture, corruption: str) -> None:
    if corruption == "missing_pair":
        fixture.frame = fixture.frame.iloc[:-1].copy()
    elif corruption == "duplicate_pair":
        fixture.frame.iloc[-1] = fixture.frame.iloc[0]
    elif corruption == "reverse_pair":
        fixture.frame.loc[0, ["trait_a_id", "trait_b_id"]] = fixture.frame.loc[
            0, ["trait_b_id", "trait_a_id"]
        ].to_numpy()
    elif corruption == "label":
        fixture.frame.loc[0, "trait_a_label"] = "invented annotation"
    elif corruption == "empty_union_value":
        fixture.frame.loc[fixture.frame.union_taxa == 0, "jaccard_fraction"] = 0.0
    else:
        changes: dict[str, tuple[str, str | float]] = {
            "negative_count": ("trait_a_taxa", -1),
            "noninteger_count": ("trait_a_taxa", 1.5),
            "union": ("union_taxa", 9),
            "fraction": ("jaccard_fraction", 0.7),
            "percent": ("jaccard_pct", 88.0),
            "nan": ("jaccard_fraction", "NaN"),
        }
        field, value = changes[corruption]
        fixture.frame[field] = fixture.frame[field].astype(object)
        fixture.frame.loc[0, field] = value
    fixture.seal()
    message = {
        "missing_pair": "pair inventory",
        "duplicate_pair": "pair inventory",
        "reverse_pair": "pair inventory",
        "label": "label differs",
        "negative_count": "nonnegative integers",
        "noninteger_count": "nonnegative integers",
        "union": "count identities",
        "fraction": "Saved Jaccard",
        "percent": "Saved Jaccard",
        "nan": "Nonfinite overlap",
        "empty_union_value": "Empty union",
    }[corruption]
    with pytest.raises(ValueError, match=message):
        fixture.load()


@pytest.mark.parametrize(
    "corruption",
    [
        "list_bytes",
        "ten_traits",
        "summary_order",
        "source_scope",
        "producer_inventory",
        "output_inventory",
        "unbound_table",
    ],
)
def test_list_source_producer_and_output_bindings(fixture: Fixture, corruption: str) -> None:
    if corruption == "list_bytes":
        fixture.traits.write_text(fixture.traits.read_text() + "# changed bytes\n")
    elif corruption == "ten_traits":
        fixture.traits.write_text(fixture.traits.read_text() + "METPO:9999999\n")
    elif corruption == "summary_order":
        fixture.summary["declared_traits"]["trait_ids"] = list(
            reversed(fixture.summary["declared_traits"]["trait_ids"])
        )
        fixture.seal()
    elif corruption == "source_scope":
        fixture.provenance["strain_rollup"] = "off"
        fixture.seal()
    elif corruption == "producer_inventory":
        del fixture.provenance["scripts"][SCRIPTS[-1]]
        fixture.seal()
    elif corruption == "output_inventory":
        fixture.consumer["outputs"] = fixture.consumer["outputs"][:-1]
        fixture.manifest.write_text(json.dumps(fixture.consumer))
    else:
        (fixture.root / TABLE).write_text(
            (fixture.root / TABLE).read_text().replace("0.333333333333", "0.533333333333")
        )
    message = {
        "list_bytes": "Declared trait list",
        "ten_traits": "explicit nine-trait list",
        "summary_order": "Declared trait list",
        "source_scope": "record scope",
        "producer_inventory": "producer inventory",
        "output_inventory": "output inventory",
        "unbound_table": "hash or size differs",
    }[corruption]
    with pytest.raises(ValueError, match=message):
        fixture.load()


def test_parsed_matrix_uses_the_pinned_buffer(
    fixture: Fixture, monkeypatch: pytest.MonkeyPatch
) -> None:
    original = Path.read_bytes
    count = 0
    table = fixture.root / TABLE
    raw = original(table)
    changed = raw.replace(b"0.333333333333", b"0.533333333333")
    assert changed != raw

    def read(path: Path) -> bytes:
        nonlocal count
        if path == table:
            count += 1
            if count == 1:
                return changed
        return original(path)

    monkeypatch.setattr(Path, "read_bytes", read)
    with pytest.raises(ValueError, match=r"hash|identity"):
        fixture.load()


def test_consumer_native_parse_must_equal_pinned_buffer(fixture: Fixture) -> None:
    verified = ConsumerRecord.model_validate_json(fixture.manifest.read_bytes())
    fixture.verifier.side_effect = None
    fixture.verifier.return_value = verified.model_copy(update={"id": "consumer:other"})
    with pytest.raises(ValueError, match="Consumer native parse differs from pinned buffer"):
        fixture.load()


def test_recheck_detects_later_change(fixture: Fixture) -> None:
    bundle = fixture.load()
    fixture.metadata.write_text(fixture.metadata.read_text() + "extra\n")
    with pytest.raises(ValueError, match="hash or size differs"):
        overlap.recheck_overlap_inputs(bundle)


@pytest.mark.parametrize("boundary", ["negative_zero", "above_one"])
def test_ratio_tolerance_does_not_relax_domain(fixture: Fixture, boundary: str) -> None:
    if boundary == "negative_zero":
        fixture.frame.loc[1, "jaccard_fraction"] = -1e-13
    else:
        fixture.frame.loc[0, "intersection_taxa"] = 2
        fixture.frame.loc[0, "union_taxa"] = 2
        fixture.frame["jaccard_fraction"] = fixture.frame.jaccard_fraction.astype(object)
        fixture.frame.loc[0, "jaccard_fraction"] = "1.000000000001"
        fixture.frame.loc[0, "jaccard_pct"] = 100.0
    fixture.seal()
    with pytest.raises(ValueError, match="outside its domain"):
        fixture.load()


def test_parser_receives_the_exact_hashed_matrix_buffer(
    fixture: Fixture, monkeypatch: pytest.MonkeyPatch
) -> None:
    original = pd.read_csv
    expected = (fixture.root / TABLE).read_bytes()
    seen: list[bytes] = []

    def read_csv(buffer: Any, *args: Any, **kwargs: Any) -> pd.DataFrame:
        assert isinstance(buffer, io.BytesIO), "Figure parser reopened a pathname"
        payload = buffer.getvalue()
        if payload.startswith(b"trait_a_id\t"):
            seen.append(payload)
        result = original(buffer, *args, **kwargs)
        assert isinstance(result, pd.DataFrame)
        return result

    monkeypatch.setattr(pd, "read_csv", read_csv)
    fixture.load()
    assert seen == [expected]


def test_recheck_rejects_identical_byte_alias_retarget(fixture: Fixture, tmp_path: Path) -> None:
    alias = tmp_path / "consumer-alias.json"
    alias.symlink_to(fixture.manifest)
    bundle = overlap.load_overlap_inputs(
        alias, manifest_sha256=str(identity(fixture.manifest)["sha256"]), trait_list=fixture.traits
    )
    replacement = tmp_path / "same-consumer.json"
    replacement.write_bytes(fixture.manifest.read_bytes())
    alias.unlink()
    alias.symlink_to(replacement)
    with pytest.raises(ValueError, match="alias changed"):
        overlap.recheck_overlap_inputs(bundle)


def test_complete_render_and_authored_wording(fixture: Fixture) -> None:
    bundle = fixture.load()
    pages = list(overlap.overlap_pages(bundle))
    assert len(pages) == 1
    page = pages[0]
    page.figure.canvas.draw()
    assert page.mark_count == len(page.row_ids) == len(set(page.row_ids)) == 36
    assert len(page.figure.axes[0].patches) == 36
    assert page.figure.get_figwidth() == pytest.approx(7.2)
    texts = [item.get_text() for item in page.figure.findobj(Text) if item.get_visible()]
    joined = " ".join(texts)
    assert not re.search(r"\b(redox|coupled|enriched|depleted|abundance|top 10)\b", joined, re.I)
    assert "annotation" in joined.lower()
    for trait, label in bundle.labels.items():
        assert trait in joined.replace("\n", "")
        assert label in joined.replace("\n", " ")
    assert page.details["pair_count"] == 36
    assert page.details["undefined_pairs"] == 1
    assert page.details["trait_ids"] == list(bundle.trait_ids)


def test_display_distinguishes_small_positive_zero_and_missing_overlap(fixture: Fixture) -> None:
    """#1243: positive overlaps must remain visibly distinct from a true zero."""
    carriers = [
        {1},
        set(range(1, 2312)),
        set(range(1, 2002)),
        set(range(1, 2001)),
        set(range(1, 2000)),
        {3000},
        set(),
        set(),
        {1, 2},
    ]
    for index, (a, b) in enumerate(combinations(carriers, 2)):
        intersection, union = len(a & b), len(a | b)
        fraction = intersection / union if union else None
        fixture.frame.loc[
            index,
            [
                "matrix_taxa_total",
                "trait_a_taxa",
                "trait_b_taxa",
                "intersection_taxa",
                "union_taxa",
                "jaccard_fraction",
                "jaccard_pct",
            ],
        ] = [
            3000,
            len(a),
            len(b),
            intersection,
            union,
            fraction,
            100 * fraction if fraction is not None else None,
        ]
    fixture.summary["coverage_and_cascade"]["matrix_taxa_total"] = 3000
    fixture.seal()
    bundle = fixture.load()
    original = overlap.overlap_audit(bundle).copy(deep=True)
    (page,) = overlap.overlap_pages(bundle)
    displayed = {item.get_position(): item.get_text() for item in page.figure.axes[0].texts}
    assert displayed[(1.5, 0.5)] == "<0.1"  # One shared carrier out of 2,311.
    assert displayed[(2.5, 0.5)] == "<0.1"  # Immediately below 0.05 percent.
    assert displayed[(3.5, 0.5)] == "0.1"  # Exactly 0.05 percent.
    assert displayed[(4.5, 0.5)] == "0.1"  # Immediately above 0.05 percent.
    assert displayed[(5.5, 0.5)] == "0.0"  # Disjoint nonempty carrier sets.
    assert displayed[(7.5, 6.5)] == "—"  # Two empty sets: undefined Jaccard.
    assert page.mark_count == 36
    pd.testing.assert_frame_equal(overlap.overlap_audit(bundle), original)


@pytest.mark.parametrize("matrix_taxa", [10, 12345])
def test_matrix_taxon_denominator_is_printed_from_saved_values(
    fixture: Fixture, matrix_taxa: int
) -> None:
    """#1244: label the validated universe; never substitute the F reference count."""
    fixture.frame["matrix_taxa_total"] = matrix_taxa
    fixture.summary["coverage_and_cascade"]["matrix_taxa_total"] = matrix_taxa
    fixture.seal()
    (page,) = overlap.overlap_pages(fixture.load())
    text = " ".join(item.get_text() for item in page.figure.findobj(Text))
    assert f"all {matrix_taxa:,} matrix taxa" in text
    assert page.details["matrix_taxa_total"] == matrix_taxa


@pytest.mark.parametrize("file_format", ["pdf", "svg"])
def test_vector_exports_do_not_rasterize_the_colorbar(fixture: Fixture, file_format: str) -> None:
    """#1228: inspect real exports; a low-resolution embedded scale must fail."""
    (page,) = overlap.overlap_pages(fixture.load())
    exported = io.BytesIO()
    page.figure.savefig(exported, format=file_format)
    raw = exported.getvalue()
    if file_format == "pdf":
        assert raw.startswith(b"%PDF")
        assert b"/Subtype /Image" not in raw, "PDF unexpectedly embeds a raster image"
    else:
        document = ET.fromstring(raw)
        assert document.tag == "{http://www.w3.org/2000/svg}svg"
        assert not list(document.iter("{http://www.w3.org/2000/svg}image")), (
            "SVG unexpectedly embeds a raster image"
        )


def test_verified_label_is_not_rewritten_by_wording_guard(fixture: Fixture) -> None:
    frame = pd.read_csv(fixture.metadata, sep="\t", dtype=str)
    frame.loc[0, "object_name"] = "redox coupled annotation"
    frame.to_csv(fixture.metadata, sep="\t", index=False)
    core = json.loads(fixture.core.read_bytes())
    core["artifacts"][0].update(identity(fixture.metadata))
    fixture.core.write_text(json.dumps(core))
    fixture.consumer["run_manifest_sha256"] = identity(fixture.core)["sha256"]
    fixture.provenance["consumed_products"]["feature_metadata"]["sha256"] = identity(
        fixture.metadata
    )["sha256"]
    fixture.frame.loc[fixture.frame.trait_a_id == frame.loc[0, "column_id"], "trait_a_label"] = (
        frame.loc[0, "object_name"]
    )
    fixture.seal()
    page = next(overlap.overlap_pages(fixture.load()))
    assert any("redox coupled annotation" in item.get_text() for item in page.figure.findobj(Text))
