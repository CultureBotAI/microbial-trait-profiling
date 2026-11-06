"""Execution-bound trait policies and non-bypass historical compatibility (#1129)."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from mtp.workflow.operations import OPERATIONS, validate_execution_product
from mtp.workflow.runtime import artifact_path, execute_workflow, validate_run
from mtp.workflow.verification import verify_core_run

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "workflow"


@pytest.fixture
def policy_inputs(tmp_path: Path) -> dict[str, Path]:
    """Keep three disjoint oxygen sources so raw and union cannot be identical."""
    edge_fields = ["id", "subject", "predicate", "object", "primary_knowledge_source"]
    with (FIXTURE / "merged-kg_edges.tsv").open(encoding="utf-8", newline="") as handle:
        edges = [
            row for row in csv.DictReader(handle, delimiter="\t") if row["id"] != "annotation3"
        ]
    for row in edges:
        if row["id"] == "annotation2":
            row["object"] = "METPO:1000612"
            row["primary_knowledge_source"] = "infores:bactotraits"
    strain = "kgmicrobe.strain:policy_fixture"
    for edge_id, predicate, obj in (
        ("strain_parent", "biolink:subclass_of", "NCBITaxon:103"),
        ("strain_605", "biolink:has_phenotype", "METPO:1000605"),
    ):
        edges.append(
            dict(
                zip(edge_fields, [edge_id, strain, predicate, obj, "infores:bacdive"], strict=True)
            )
        )
    edges_path = tmp_path / "policy-edges.tsv"
    with edges_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=edge_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(edges)
    with (FIXTURE / "merged-kg_nodes.tsv").open(encoding="utf-8", newline="") as handle:
        nodes = list(csv.DictReader(handle, delimiter="\t"))
    for node in nodes:
        if node["id"] == "METPO:1000608":
            node["name"] = "facultatively aerobic"
    nodes.extend(
        [
            {
                "id": "METPO:1000605",
                "category": "biolink:PhenotypicQuality",
                "name": "facultatively anaerobic",
            },
            {
                "id": "METPO:1000612",
                "category": "biolink:PhenotypicQuality",
                "name": "facultative oxygen preference",
            },
            {
                "id": strain,
                "category": "biolink:OrganismTaxon",
                "name": "Examplegenus gamma fixture strain",
            },
        ]
    )
    nodes_path = tmp_path / "policy-nodes.tsv"
    with nodes_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["id", "category", "name"], delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(nodes)
    return {"raw_dump": FIXTURE / "full_dump.csv", "kg_edges": edges_path, "kg_nodes": nodes_path}


@pytest.mark.parametrize("spec_name", ["bugsigdb_core.yaml", "bugsigdb_core_subset.yaml"])
@pytest.mark.parametrize("policy", ["raw_object_v1", "facultative_union_v1"])
def test_execution_binds_policy_to_run_steps_metadata_and_complete_verification(
    tmp_path: Path, spec_name: str, policy: str, policy_inputs: dict[str, Path]
) -> None:
    run = execute_workflow(
        ROOT / "workflows" / spec_name,
        policy_inputs,
        tmp_path / "runs",
        OPERATIONS,
        parameter_overrides={
            "source_mode": "fixture",
            "trait_key_policy": policy,
            "strain_rollup": "phenotype",
        },
        validator=validate_execution_product,
        repository_root=ROOT,
    )
    assert run.trait_key_policy == policy
    assert {item.name: item.value for item in run.parameters or []}["trait_key_policy"] == policy
    with artifact_path(run, "taxa_matrix").open(encoding="utf-8", newline="") as handle:
        matrix = {row["taxon_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    members = ["METPO:1000605", "METPO:1000608", "METPO:1000612"]
    expected_members = members if policy == "raw_object_v1" else ["METPO:1000608"]
    genus = matrix["NCBITaxon:100"]
    assert set(genus).intersection(members) == set(expected_members)
    if policy == "raw_object_v1":
        # Three disjoint sources: BacDive strain605, direct608, BactoTraits612.
        assert [int(genus[member]) for member in members] == [1, 1, 1]
        assert [int(matrix["NCBITaxon:101"][member]) for member in members] == [0, 1, 0]
        assert [int(matrix["NCBITaxon:103"][member]) for member in members] == [1, 0, 0]
    else:
        assert int(genus["METPO:1000608"]) == 3
        assert int(matrix["NCBITaxon:101"]["METPO:1000608"]) == 1
        assert int(matrix["NCBITaxon:103"]["METPO:1000608"]) == 1
    for product in ("count_profiles", "ora_profiles"):
        with artifact_path(run, product).open(encoding="utf-8", newline="") as handle:
            header = next(csv.reader(handle, delimiter="\t"))
        assert set(header).intersection(members) == set(expected_members)
    for step in run.steps or []:
        if step.step_id in {"taxa_matrix", "feature_catalog"}:
            assert {item.name: item.value for item in step.parameters or []}[
                "trait_key_policy"
            ] == policy
    expected = ["METPO:1000605", "METPO:1000608", "METPO:1000612"]
    for product in ("feature_metadata", "feature_catalog"):
        with artifact_path(run, product).open(encoding="utf-8", newline="") as handle:
            records = list(csv.DictReader(handle, delimiter="\t"))
        focal = next(row for row in records if row["column_id"] == "METPO:1000608")
        assert focal["object_name"] == "facultatively aerobic"
        assert focal["column_label"] == "facultatively aerobic METPO:1000608"
        assert focal["trait_key_policy"] == policy
        assert json.loads(focal["trait_member_ids_json"]) == (
            expected if policy == "facultative_union_v1" else ["METPO:1000608"]
        )
        assert {row["column_id"] for row in records}.intersection(members) == set(expected_members)
        for row in records:
            assert row["trait_key_policy"] == policy
            if row["column_id"] in expected_members:
                assert json.loads(row["trait_member_ids_json"]) == (
                    expected if policy == "facultative_union_v1" else [row["column_id"]]
                )
    verified = verify_core_run(Path(run.manifest_path))
    assert verified == run
    original = Path(run.manifest_path).read_bytes()
    raw = json.loads(original)
    raw["trait_key_policy"] = (
        "raw_object_v1" if policy == "facultative_union_v1" else "facultative_union_v1"
    )
    Path(run.manifest_path).write_text(json.dumps(raw), encoding="utf-8")
    with pytest.raises(ValueError, match=r"(?i)policy"):
        validate_run(Path(run.manifest_path))
    Path(run.manifest_path).write_bytes(original)
    assert verify_core_run(Path(run.manifest_path)) == run


@pytest.mark.parametrize(
    ("policy", "expected_label"),
    [
        ("facultative_union_v1", "facultative (aerobe or anaerobe) (METPO:1000608)"),
        ("raw_object_v1", "facultatively aerobic (METPO:1000608)"),
    ],
)
def test_table1_trait_maps_core_run_binds_feature_metadata_and_renders_policy_label(
    tmp_path: Path,
    policy_inputs: dict[str, Path],
    policy: str,
    expected_label: str,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """table1_trait_maps --core-run binds feature_metadata under steps/taxa_matrix/derived."""
    from mtp import table1_trait_maps

    run = execute_workflow(
        ROOT / "workflows" / "bugsigdb_core.yaml",
        policy_inputs,
        tmp_path / "runs",
        OPERATIONS,
        parameter_overrides={
            "source_mode": "fixture",
            "trait_key_policy": policy,
            "strain_rollup": "phenotype",
        },
        validator=validate_execution_product,
        repository_root=ROOT,
    )
    manifest = Path(run.manifest_path)
    out_dir = tmp_path / "maps_out"
    monkeypatch.setattr(
        "plotly.graph_objs.Figure.write_image",
        lambda _self, path, **_kw: Path(path).write_bytes(b"png"),
    )
    monkeypatch.setattr(
        table1_trait_maps, "TABLE1_TRAITS", {"facultatively_aerobic": ("METPO:1000608",)}
    )
    monkeypatch.setattr(table1_trait_maps, "TRAIT_ORDER", ("facultatively_aerobic",))
    monkeypatch.setattr(table1_trait_maps, "FOCUS_PANELS", ())
    monkeypatch.setattr(
        "sys.argv",
        [
            "table1_trait_maps",
            "--core-run",
            str(manifest),
            "--output-dir",
            str(out_dir),
            "--min-studies",
            "1",
            "--min-site-conditions",
            "1",
        ],
    )
    table1_trait_maps.main()
    shift_tsv = out_dir / "table1_traits_direction_shift.tsv"
    with shift_tsv.open(encoding="utf-8", newline="") as handle:
        shift_rows = list(csv.DictReader(handle, delimiter="\t"))
    fac_row = next(row for row in shift_rows if row["trait"] == "facultatively_aerobic")
    assert fac_row["trait_label"] == expected_label

    monkeypatch.setattr(
        "sys.argv",
        [
            "table1_trait_maps",
            "--core-run",
            str(manifest),
            "--column-labels",
            str(tmp_path / "other_columns.tsv"),
            "--output-dir",
            str(tmp_path / "maps_conflict"),
        ],
    )
    with pytest.raises(ValueError, match="--core-run conflicts with explicit --column-labels"):
        table1_trait_maps.main()
