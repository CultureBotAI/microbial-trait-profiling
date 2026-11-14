"""Explicit contribution keys preserve historical raw and current union inputs (#1354)."""

from __future__ import annotations

import csv
import gzip
import importlib.util
import json
import sys
from collections.abc import Callable
from pathlib import Path
from types import ModuleType, SimpleNamespace
from typing import cast

import pandas as pd
import pytest

from mtp import influence_extension as extension
from mtp import taxon_contributions, taxon_influence
from mtp.contribution_inputs import (
    Base,
    digest,
    load_graph,
    resolve_inputs_trait_key_policy,
)
from mtp.datamodel.workflow_pydantic import FileInventoryEntry, WorkflowRun
from mtp.taxa_trait_matrix import StrainRollup
from mtp.trait_key_policy import TraitKeyPolicy
from mtp.workflow.trait_policy import HISTORICAL_MATRIX_POLICIES, MATRIX_SOURCE

MEMBERS = {"METPO:1000605", "METPO:1000608", "METPO:1000612"}
NEGATED = "METPO:2000021|METPO:1000605"


@pytest.fixture
def graph_files(tmp_path: Path) -> tuple[Path, Path]:
    edges, nodes = tmp_path / "edges.tsv", tmp_path / "nodes.tsv"
    edges.write_text(
        "id\tsubject\tpredicate\tobject\tprimary_knowledge_source\n"
        "raw605\tNCBITaxon:1\tbiolink:has_phenotype\tMETPO:1000605\tinfores:bacdive\n"
        "raw608\tNCBITaxon:2\tbiolink:has_phenotype\tMETPO:1000608\tinfores:madin_etal\n"
        "raw612\tNCBITaxon:3\tbiolink:has_phenotype\tMETPO:1000612\tinfores:bactotraits\n"
        "negative\tNCBITaxon:1\tMETPO:2000021\tMETPO:1000605\tinfores:bacdive\n"
        "parent\tkgmicrobe.strain:4\tbiolink:subclass_of\tNCBITaxon:4\tinfores:bacdive\n"
        "strain605\tkgmicrobe.strain:4\tbiolink:has_phenotype\tMETPO:1000605\tinfores:bacdive\n",
        encoding="utf-8",
    )
    nodes.write_text(
        "id\tname\nNCBITaxon:1\tSpecies one\nNCBITaxon:2\tSpecies two\n"
        "NCBITaxon:3\tSpecies three\nNCBITaxon:4\tSpecies four\n",
        encoding="utf-8",
    )
    return edges, nodes


@pytest.mark.parametrize("scope", ["off", "phenotype", "full"])
def test_raw_keys_preserve_objects_predicates_and_strain_scope(
    graph_files: tuple[Path, Path], scope: StrainRollup
) -> None:
    columns = MEMBERS | {NEGATED}
    graph = load_graph(
        *graph_files, columns, columns, strain_rollup=scope, trait_key_policy="raw_object_v1"
    )
    expected = {
        "NCBITaxon:1": {"METPO:1000605", NEGATED},
        "NCBITaxon:2": {"METPO:1000608"},
        "NCBITaxon:3": {"METPO:1000612"},
    }
    if scope != "off":
        expected["NCBITaxon:4"] = {"METPO:1000605"}
    if scope == "full":
        # The established full scope removes even directly asserted negative keys.
        expected["NCBITaxon:1"].remove(NEGATED)
    assert graph.traits == expected
    assert {
        (taxon, key): [(a.edge_id, a.object_id, a.primary_source) for a in evidence]
        for (taxon, key), evidence in graph.assertions.items()
    } == {
        ("NCBITaxon:1", "METPO:1000605"): [("raw605", "METPO:1000605", "infores:bacdive")],
        ("NCBITaxon:2", "METPO:1000608"): [("raw608", "METPO:1000608", "infores:madin_etal")],
        ("NCBITaxon:3", "METPO:1000612"): [("raw612", "METPO:1000612", "infores:bactotraits")],
        ("NCBITaxon:1", NEGATED): [("negative", "METPO:1000605", "infores:bacdive")],
        **(
            {("NCBITaxon:4", "METPO:1000605"): [("strain605", "METPO:1000605", "infores:bacdive")]}
            if scope != "off"
            else {}
        ),
    }


def test_raw_policy_is_not_inferred_from_the_retained_column_vocabulary(
    graph_files: tuple[Path, Path],
) -> None:
    columns = {"METPO:1000608"}
    raw = load_graph(*graph_files, columns, columns, trait_key_policy="raw_object_v1")
    union = load_graph(*graph_files, columns, columns, trait_key_policy="facultative_union_v1")
    assert raw.traits == {"NCBITaxon:2": columns}
    assert set(raw.assertions) == {("NCBITaxon:2", "METPO:1000608")}
    assert union.traits == {f"NCBITaxon:{i}": columns for i in range(1, 5)}
    assert load_graph(*graph_files, columns, columns) == union
    # Selecting a different policy does not change the next call's default.
    assert load_graph(*graph_files, columns, columns, trait_key_policy="raw_object_v1") == raw


def test_unknown_policy_is_refused_before_reading_inputs(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="trait_key_policy"):
        load_graph(
            tmp_path / "absent_edges.tsv",
            tmp_path / "absent_nodes.tsv",
            set(),
            set(),
            trait_key_policy=cast(TraitKeyPolicy, "unknown"),
        )


@pytest.mark.parametrize("policy", ["raw_object_v1", "facultative_union_v1"])
def test_extension_forwards_the_verified_core_policy(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, policy: TraitKeyPolicy
) -> None:
    # Core execution verification is a separate expensive boundary. The real policy
    # resolver still runs here against a validated model and its audited producer ID.
    fixture = Path(__file__).parent / "fixtures/trait_key_policy/historical_raw_core/core/run.json"
    run = WorkflowRun.model_validate_json(fixture.read_bytes())
    identity = next(key for key, value in HISTORICAL_MATRIX_POLICIES.items() if value == policy)
    run.code_inventory = [
        FileInventoryEntry(path=MATRIX_SOURCE, sha256=identity[0], size_bytes=identity[1])
    ]
    source = tmp_path / "source.tsv"
    source.write_text("fixture input\n", encoding="utf-8")
    bundle = SimpleNamespace(run=run, manifest=fixture, path=lambda _role: source)
    monkeypatch.setattr(extension, "open_core_bundle", lambda _path: bundle)
    monkeypatch.setattr(extension, "preflight_consumer_outputs", lambda *_args, **_kwargs: None)
    base = Base(sorted(MEMBERS), {}, {}, {}, {}, {}, {}, {}, [])
    monkeypatch.setattr(extension, "load_base", lambda *_args: base)

    class GraphReached(Exception):
        pass

    def inspect_graph(*_args: object, **kwargs: object) -> None:
        assert kwargs.get("trait_key_policy") == policy
        raise GraphReached

    monkeypatch.setattr(extension, "load_graph", inspect_graph)
    with pytest.raises(GraphReached):
        extension.run_extension(fixture, tmp_path / "output")
    assert not (tmp_path / "output").exists()


def test_extension_inventory_binds_trait_policy_implementations() -> None:
    inventory = extension.implementation_hashes()
    assert {"trait_key_policy.py", "workflow/trait_policy.py"}.issubset(inventory)


TRAIT_608 = "METPO:1000608"
COLUMNS_SIX = [
    TRAIT_608,
    "GO:0030245",
    "METPO:1000668",
    "METPO:1000683",
    "pathways:arsenate_reduction",
    "METPO:1000623",
]


def _write_core_manifest(path: Path, policy: TraitKeyPolicy) -> Path:
    fixture = Path(__file__).parent / "fixtures/trait_key_policy/historical_raw_core/core/run.json"
    run = WorkflowRun.model_validate_json(fixture.read_bytes())
    identity = next(key for key, value in HISTORICAL_MATRIX_POLICIES.items() if value == policy)
    run.code_inventory = [
        FileInventoryEntry(path=MATRIX_SOURCE, sha256=identity[0], size_bytes=identity[1])
    ]
    path.write_text(json.dumps(run.model_dump(mode="json"), indent=2) + "\n", encoding="utf-8")
    return path


def _write_synthetic_contribution_inputs(
    root: Path,
    policy: TraitKeyPolicy,
    *,
    with_run_manifest: bool = True,
    with_column_metadata: bool = False,
) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    edges = root / "merged-kg_edges.tsv"
    nodes = root / "merged-kg_nodes.tsv"
    edges.write_text(
        "id\tsubject\tpredicate\tobject\n"
        "a\tNCBITaxon:1\tbiolink:has_phenotype\tMETPO:1000605\n"
        "b\tNCBITaxon:2\tbiolink:has_phenotype\tMETPO:1000608\n"
        "c\tNCBITaxon:3\tbiolink:has_phenotype\tMETPO:1000612\n"
        "f\tNCBITaxon:5\tbiolink:has_phenotype\tMETPO:1000609\n"
        "d\tkgmicrobe.strain:4\tbiolink:subclass_of\tNCBITaxon:4\n"
        "e\tkgmicrobe.strain:4\tbiolink:has_phenotype\tMETPO:1000605\n"
        + "".join(
            f"p{i}\tNCBITaxon:{i}\tbiolink:subclass_of\tNCBITaxon:9\n" for i in (1, 2, 3, 4, 5)
        ),
        encoding="utf-8",
    )
    nodes.write_text("id\tname\nNCBITaxon:9\tGenus\n", encoding="utf-8")
    count = 1 if policy == "raw_object_v1" else 4
    extra_cols = ["METPO:1000605", "METPO:1000609"] if policy == "raw_object_v1" else []
    all_cols = [*COLUMNS_SIX, *extra_cols]
    extra_vals = "\t2\t1" if policy == "raw_object_v1" else ""
    values = f"{count}\t0\t0\t0\t0\t0{extra_vals}\n"
    (root / "study_trait_profiles_wide.tsv").write_text(
        "study_id\texperiment_id\tinc_dec\t"
        + "\t".join(all_cols)
        + "\n"
        + "S\tE\tinc\t"
        + values
        + "S\tE\tdec\t"
        + values,
        encoding="utf-8",
    )
    (root / "taxa_trait_matrix.tsv").write_text(
        "taxon_id\ttaxon_name\ttaxonomic_rank\t"
        + "\t".join(all_cols)
        + "\n"
        + "NCBITaxon:9\tGenus\tgenus\t"
        + values,
        encoding="utf-8",
    )
    (root / "experiments.tsv").write_text(
        "study_id\texperiment_id\tcondition\tbody_site\thost_species\n"
        "S\tE\tParkinson disease\tFeces\tHomo sapiens\n",
        encoding="utf-8",
    )
    (root / "bugsigdb_taxa.tsv").write_text(
        "study_id\texperiment_id\tdirection\ttaxon_id\ttaxon_name\ttaxonomic_rank\tsignature_id\n"
        "S\tE\tincreased\tNCBITaxon:9\tGenus\tgenus\tinc\n"
        "S\tE\tdecreased\tNCBITaxon:9\tGenus\tgenus\tdec\n",
        encoding="utf-8",
    )
    if with_run_manifest:
        _write_core_manifest(root / "run.json", policy)
    if with_column_metadata:
        (root / "taxa_trait_matrix_columns.tsv").write_text(
            f"column_id\tobject_id\ttrait_key_policy\n{TRAIT_608}\t{TRAIT_608}\t{policy}\n",
            encoding="utf-8",
        )
    protocol = root / "protocol.json"
    protocol.write_text(
        json.dumps({"kg_md5": {p.name: digest(p, "md5") for p in (edges, nodes)}}),
        encoding="utf-8",
    )
    return protocol


@pytest.fixture
def pd_contributors_mod() -> ModuleType:
    path = Path(__file__).resolve().parents[1] / "scripts/create_pd_trait_contributors.py"
    sys.path.insert(0, str(path.parent))
    try:
        spec = importlib.util.spec_from_file_location("test_pd_contributors_1359", path)
        assert spec is not None
        assert spec.loader is not None
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    finally:
        sys.path.pop(0)


def _assert_raises_value_error(fn: Callable[[], object], match: str) -> None:
    message = ""
    try:
        fn()
    except ValueError as exc:
        message = str(exc)
    else:
        raise AssertionError(f"Expected ValueError matching {match!r}")
    assert match in message, f"Expected {match!r} in {message!r}"


@pytest.mark.parametrize("policy", ["raw_object_v1", "facultative_union_v1"])
def test_resolve_inputs_trait_key_policy_inherits_and_validates_core_and_column_metadata(
    tmp_path: Path,
    policy: TraitKeyPolicy,
) -> None:
    other: TraitKeyPolicy = "facultative_union_v1" if policy == "raw_object_v1" else "raw_object_v1"
    manifest_dir = tmp_path / f"manifest_{policy}"
    _write_synthetic_contribution_inputs(
        manifest_dir, policy, with_run_manifest=True, with_column_metadata=False
    )
    assert resolve_inputs_trait_key_policy(manifest_dir) == policy
    assert resolve_inputs_trait_key_policy(manifest_dir, policy) == policy
    _assert_raises_value_error(
        lambda: resolve_inputs_trait_key_policy(manifest_dir, other),
        "contradicts recorded input policy",
    )

    col_dir = tmp_path / f"columns_{policy}"
    _write_synthetic_contribution_inputs(
        col_dir, policy, with_run_manifest=False, with_column_metadata=True
    )
    assert resolve_inputs_trait_key_policy(col_dir) == policy
    _assert_raises_value_error(
        lambda: resolve_inputs_trait_key_policy(col_dir, other),
        "contradicts recorded input policy",
    )

    conflict_manifest = _write_core_manifest(tmp_path / f"conflict_{policy}.json", other)
    _assert_raises_value_error(
        lambda: resolve_inputs_trait_key_policy(col_dir, core_manifest=conflict_manifest),
        "contradicts",
    )


def test_taxon_contributions_cli_inherits_core_manifest_policy_and_rejects_mismatch(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    inputs = tmp_path / "raw_core_inputs"
    protocol = _write_synthetic_contribution_inputs(
        inputs, "raw_object_v1", with_run_manifest=True, with_column_metadata=True
    )
    output = tmp_path / "ledger_output"
    monkeypatch.setattr(taxon_influence, "run_influence", lambda *_args: {})
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "ledger",
            "--kg-dir",
            str(inputs),
            "--inputs",
            str(inputs),
            "--protocol",
            str(protocol),
            "--output",
            str(output),
            "--strain-rollup",
            "phenotype",
        ],
    )
    try:
        taxon_contributions.main()
    except ValueError as exc:
        raise AssertionError(f"Unexpected reconciliation failure: {exc}") from exc
    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["completed"] is True
    assert manifest["trait_key_policy"] == "raw_object_v1"
    assert manifest["mismatches"] == 0
    with gzip.open(output / "source_contributions.tsv.gz", "rt", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["source_taxon"] == "NCBITaxon:2"
    assert rows[0]["trait_key"] == TRAIT_608

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "ledger",
            "--kg-dir",
            str(inputs),
            "--inputs",
            str(inputs),
            "--protocol",
            str(protocol),
            "--output",
            str(tmp_path / "bad_output"),
            "--strain-rollup",
            "phenotype",
            "--trait-key-policy",
            "facultative_union_v1",
        ],
    )
    _assert_raises_value_error(
        taxon_contributions.main,
        "contradicts recorded input policy",
    )


def test_pd_contributors_inherits_core_policy_and_rejects_policy_mismatch(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    pd_contributors_mod: ModuleType,
) -> None:
    inputs = tmp_path / "pd_raw_core"
    protocol = _write_synthetic_contribution_inputs(
        inputs, "raw_object_v1", with_run_manifest=True, with_column_metadata=True
    )
    try:
        _, graph, ledger = pd_contributors_mod.load_contribution_ledger(
            inputs,
            inputs,
            protocol,
            "phenotype",
        )
    except ValueError as exc:
        raise AssertionError(f"Unexpected PD ledger failure: {exc}") from exc
    assert ledger.mismatches == []
    assert ledger.reconstructed["NCBITaxon:9"][TRAIT_608] == 1
    assert ledger.reconstructed["NCBITaxon:9"]["METPO:1000605"] == 2
    assert ledger.reconstructed["NCBITaxon:9"]["METPO:1000609"] == 1
    assert graph.traits["NCBITaxon:5"] == {"METPO:1000609"}
    assert "METPO:1000608" not in graph.traits["NCBITaxon:1"]

    _assert_raises_value_error(
        lambda: pd_contributors_mod.load_contribution_ledger(
            inputs,
            inputs,
            protocol,
            "phenotype",
            trait_key_policy="facultative_union_v1",
        ),
        "contradicts recorded input policy",
    )

    output = tmp_path / "pd_output"
    output.mkdir()
    metadata = pd.DataFrame(
        [
            {
                "study_id": "S",
                "experiment_id": "E",
                "inc_dec": direction,
                "profile_id": f"S_E_{direction}",
                "condition": pd_contributors_mod.PD_CONDITION,
                "body_site_group": pd_contributors_mod.FECES,
            }
            for direction in ("inc", "dec")
        ]
    )
    traits = pd.DataFrame([{"trait_key": TRAIT_608, "trait_slug": "facultatively_aerobic"}])
    cells = pd.DataFrame(columns=["rank", "taxon", "trait_slug"])
    selected = pd.DataFrame(
        columns=[
            "scope",
            "rank",
            "taxon",
            "name",
            "display_name",
            "pd_feces_profiles",
            "other_feces_profiles",
            "pd_other_feces_counts",
            "pd_feces_studies",
            "other_feces_studies",
            "pd_specificity_fraction",
            "pd_vs_other_feces_log2_ratio",
            "pd_feces_total_profiles",
            "other_feces_total_profiles",
            "pd_feces_prevalence",
            "other_feces_prevalence",
            "pd_vs_other_feces_log2_prevalence_ratio",
            "pd_minus_other_feces_prevalence_difference",
            "pd_other_feces_display",
        ]
    )
    for name, value in {
        "REPO_ROOT": tmp_path,
        "KG_DIR": inputs,
        "DERIVED_DIR": inputs,
        "KG_PROTOCOL_PATH": protocol,
        "OUTPUT_DIR": output,
    }.items():
        monkeypatch.setattr(pd_contributors_mod, name, value)
    monkeypatch.setattr(pd_contributors_mod, "load_profile_metadata", lambda: metadata.copy())
    monkeypatch.setattr(
        pd_contributors_mod, "require_current_trait_taxa_evidence", lambda *_a: None
    )
    monkeypatch.setattr(pd_contributors_mod, "read_tsv", lambda *_a: pd.DataFrame())
    monkeypatch.setattr(pd_contributors_mod, "trait_directions", lambda *_a: pd.DataFrame())
    monkeypatch.setattr(pd_contributors_mod, "check_trait_intervals", lambda *_a: None)
    monkeypatch.setattr(pd_contributors_mod, "trait_frame", lambda *_a: traits.copy())
    monkeypatch.setattr(
        pd_contributors_mod,
        "build_contributor_tables",
        lambda _c, cov, _t, _n: (cells.copy(), cov.copy(), selected.copy()),
    )
    monkeypatch.setattr(pd_contributors_mod, "augment_selected_taxa", lambda frame, _m: frame)
    monkeypatch.setattr(
        pd_contributors_mod, "build_direction_differences", lambda *_a: pd.DataFrame()
    )
    monkeypatch.setattr(
        pd_contributors_mod, "select_large_difference_taxa", lambda *_a: selected.copy()
    )
    monkeypatch.setattr(
        pd_contributors_mod,
        "build_trait_cocarriage_tables",
        lambda *_a: (pd.DataFrame(), pd.DataFrame(), pd.DataFrame()),
    )
    monkeypatch.setattr(pd_contributors_mod, "plot_contributor_heatmap", lambda *_a: None)
    monkeypatch.setattr(pd_contributors_mod, "plot_trait_cocarriage_network", lambda *_a: None)
    monkeypatch.setattr(pd_contributors_mod, "write_note", lambda *_a: None)
    monkeypatch.setattr(sys, "argv", ["pd"])

    try:
        pd_contributors_mod.main()
    except ValueError as exc:
        raise AssertionError(f"Unexpected PD main failure: {exc}") from exc
    pd_manifest = json.loads(
        (output / f"{pd_contributors_mod.STEM}_manifest.json").read_text(encoding="utf-8")
    )
    assert pd_manifest["trait_key_policy"] == "raw_object_v1"
