"""Explicit contribution keys preserve historical raw and current union inputs (#1354)."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from typing import cast

import pytest

from mtp import influence_extension as extension
from mtp.contribution_inputs import Base, load_graph
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
