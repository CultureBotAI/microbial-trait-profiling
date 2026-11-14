"""Audit KG readers disclose and honor raw versus union keys (#1136)."""

from __future__ import annotations

import csv
import importlib.util
import json
import sys
from pathlib import Path
from types import ModuleType
from typing import cast

import pytest

from mtp.paths import REPO_ROOT
from mtp.trait_key_policy import TraitKeyPolicy

COVERAGE_SCRIPT = "issue774_strain_rollup_coverage.py"
RECOUNT_SCRIPT = "recount_table1_audit_kg_columns_807.py"


def load_script(filename: str) -> ModuleType:
    path = REPO_ROOT / "scripts" / "adhoc" / filename
    spec = importlib.util.spec_from_file_location(f"policy_audit_{path.stem}", path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def kg_inputs(tmp_path: Path) -> dict[str, Path]:
    """One direct carrier per raw member, plus a strain-only facultative carrier."""
    paths = {name: tmp_path / f"{name}.tsv" for name in ("edges", "nodes", "reported_taxa")}
    paths["reported_taxa"].write_text(
        "direction\ttaxon_id\ttaxon_name\ttaxonomic_rank\n"
        + "".join(f"increased\tNCBITaxon:{n}\tSpecies {n}\tspecies\n" for n in range(1, 5)),
        encoding="utf-8",
    )
    edges = [
        ("raw605", "NCBITaxon:1", "biolink:has_phenotype", "METPO:1000605"),
        ("raw608", "NCBITaxon:2", "biolink:has_phenotype", "METPO:1000608"),
        ("raw612", "NCBITaxon:3", "biolink:has_phenotype", "METPO:1000612"),
        ("strain_link", "kgmicrobe.strain:4", "biolink:subclass_of", "NCBITaxon:4"),
        ("strain605", "kgmicrobe.strain:4", "biolink:has_phenotype", "METPO:1000605"),
        ("coccus", "NCBITaxon:1", "biolink:has_phenotype", "METPO:1000668"),
        ("sphere", "NCBITaxon:1", "biolink:has_phenotype", "METPO:1000683"),
    ]
    # The recount's historical source table uses primary_knowledge_source at column 10.
    header = [
        "id",
        "subject",
        "predicate",
        "object",
        *(f"aux_{n}" for n in range(6)),
        "primary_knowledge_source",
    ]
    paths["edges"].write_text(
        "\t".join(header)
        + "\n"
        + "".join(
            "\t".join((*edge, *("" for _ in range(6)), "infores:fixture")) + "\n" for edge in edges
        ),
        encoding="utf-8",
    )
    names = {
        "METPO:1000605": "facultatively anaerobic",
        "METPO:1000608": "facultatively aerobic",
        "METPO:1000612": "facultative aerobe anaerobe",
        "METPO:1000668": "coccus shaped",
        "METPO:1000683": "sphere shaped",
    }
    paths["nodes"].write_text(
        "id\tcategory\tname\n"
        + "".join(f"{key}\tbiolink:PhenotypicQuality\t{name}\n" for key, name in names.items()),
        encoding="utf-8",
    )
    return paths


@pytest.mark.parametrize(
    ("policy", "off", "phenotype"),
    [("raw_object_v1", 1, 1), ("facultative_union_v1", 3, 4)],
)
def test_coverage_api_honors_and_records_policy(
    kg_inputs: dict[str, Path], policy: TraitKeyPolicy, off: int, phenotype: int
) -> None:
    module = load_script(COVERAGE_SCRIPT)
    rows, summary = module.audit(
        kg_inputs["reported_taxa"],
        kg_inputs["edges"],
        kg_inputs["nodes"],
        None,
        trait_key_policy=policy,
    )
    assert summary["trait_key_policy"] == policy
    assert all(row["trait_key_policy"] == policy for row in rows)
    facultative = {row["trait_key"]: row for row in rows}["METPO:1000608"]
    assert facultative["ncbitaxon_nodes_off"] == off
    assert facultative["reported_taxa_off"] == off
    assert facultative["ncbitaxon_nodes_phenotype"] == phenotype
    assert facultative["reported_taxa_phenotype"] == phenotype
    keys = {row["trait_key"] for row in rows}
    if policy == "raw_object_v1":
        assert {"METPO:1000605", "METPO:1000612"} <= keys
    else:
        assert not {"METPO:1000605", "METPO:1000612"} & keys


@pytest.mark.parametrize(
    ("policy", "off", "phenotype"),
    [("raw_object_v1", 1, 1), ("facultative_union_v1", 3, 4)],
)
def test_recount_api_separates_source_objects_from_mapped_counts(
    kg_inputs: dict[str, Path],
    policy: TraitKeyPolicy,
    off: int,
    phenotype: int,
    capsys: pytest.CaptureFixture[str],
) -> None:
    module = load_script(RECOUNT_SCRIPT)
    module.recount(kg_inputs["edges"], kg_inputs["reported_taxa"], trait_key_policy=policy)
    output = capsys.readouterr().out
    assert f"trait_key_policy={policy}" in output
    assert "Source edge counts: raw KG objects (not policy-mapped)" in output
    assert "Coverage counts: trait-key-policy-mapped" in output
    source_line = next(line for line in output.splitlines() if line.startswith("  facultatively"))
    # Only one raw 1000608 edge exists under either policy; 1000605/1000612 are distinct sources.
    assert source_line.split()[-1] == "1"
    mapped_line = next(line for line in output.splitlines() if line.startswith("facultatively"))
    assert f"direct NCBITaxon nodes={off} phenotype nodes={phenotype}" in mapped_line
    assert f"reported taxa bearing each trait, off: [{off}, 0, 1, 1, 0, 0]" in output
    assert f"reported taxa bearing each trait, phenotype: [{phenotype}, 0, 1, 1, 0, 0]" in output


@pytest.mark.parametrize("policy", ["raw_object_v1", "facultative_union_v1"])
def test_coverage_cli_records_consumed_policy(
    kg_inputs: dict[str, Path],
    tmp_path: Path,
    policy: TraitKeyPolicy,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    module = load_script(COVERAGE_SCRIPT)
    output = tmp_path / "coverage"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            COVERAGE_SCRIPT,
            "--reported-taxa",
            str(kg_inputs["reported_taxa"]),
            "--edges",
            str(kg_inputs["edges"]),
            "--nodes",
            str(kg_inputs["nodes"]),
            "--output-dir",
            str(output),
            "--trait-key-policy",
            policy,
        ],
    )
    module.main()
    summary = json.loads((output / "summary.json").read_text(encoding="utf-8"))
    assert summary["trait_key_policy"] == policy
    assert json.loads(capsys.readouterr().out)["trait_key_policy"] == policy
    expected = 1 if policy == "raw_object_v1" else 4
    assert (
        summary["scopes"]["phenotype"]["focal_traits"]["METPO:1000608"]["reported_taxa"] == expected
    )
    assert (output / "trait_coverage.tsv").stat().st_size > 0
    with (output / "trait_coverage.tsv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows
    assert all(row["trait_key_policy"] == policy for row in rows)


@pytest.mark.parametrize("policy", ["raw_object_v1", "facultative_union_v1"])
def test_recount_cli_records_consumed_policy(
    kg_inputs: dict[str, Path],
    policy: TraitKeyPolicy,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    module = load_script(RECOUNT_SCRIPT)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            RECOUNT_SCRIPT,
            "--reported-taxa",
            str(kg_inputs["reported_taxa"]),
            "--edges",
            str(kg_inputs["edges"]),
            "--trait-key-policy",
            policy,
        ],
    )
    module.main()
    output = capsys.readouterr().out
    assert f"trait_key_policy={policy}" in output
    expected = 1 if policy == "raw_object_v1" else 4
    assert f"reported taxa bearing each trait, phenotype: [{expected}, 0, 1, 1, 0, 0]" in output


@pytest.mark.parametrize("script", [COVERAGE_SCRIPT, RECOUNT_SCRIPT])
def test_unknown_policy_is_refused_before_audit_inputs_are_read(
    script: str, tmp_path: Path
) -> None:
    module = load_script(script)
    absent = tmp_path / "absent.tsv"
    call = module.audit if script == COVERAGE_SCRIPT else module.recount
    args = (absent, absent, absent, None) if script == COVERAGE_SCRIPT else (absent, absent)
    with pytest.raises(ValueError, match="trait_key_policy"):
        call(*args, trait_key_policy=cast(TraitKeyPolicy, "unknown"))


@pytest.mark.parametrize("script", [COVERAGE_SCRIPT, RECOUNT_SCRIPT])
def test_default_api_policy_remains_current_union_and_is_disclosed(
    script: str, kg_inputs: dict[str, Path], capsys: pytest.CaptureFixture[str]
) -> None:
    module = load_script(script)
    if script == COVERAGE_SCRIPT:
        _, summary = module.audit(
            kg_inputs["reported_taxa"], kg_inputs["edges"], kg_inputs["nodes"], None
        )
        assert summary["trait_key_policy"] == "facultative_union_v1"
    else:
        module.recount(kg_inputs["edges"], kg_inputs["reported_taxa"])
        assert "trait_key_policy=facultative_union_v1" in capsys.readouterr().out
