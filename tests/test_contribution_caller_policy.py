"""Raw and union core products remain distinct at contribution entrypoints (#1359)."""

from __future__ import annotations

import csv
import importlib.util
import json
import sys
from pathlib import Path
from types import ModuleType

import pandas as pd
import pytest

from mtp import taxon_contributions, taxon_influence
from mtp.contribution_inputs import digest
from mtp.trait_key_policy import TraitKeyPolicy

TRAIT = "METPO:1000608"
COLUMNS = [
    TRAIT,
    "GO:0030245",
    "METPO:1000668",
    "METPO:1000683",
    "pathways:arsenate_reduction",
    "METPO:1000623",
]


@pytest.fixture
def contributors() -> ModuleType:
    path = Path(__file__).resolve().parents[1] / "scripts/create_pd_trait_contributors.py"
    sys.path.insert(0, str(path.parent))
    try:
        spec = importlib.util.spec_from_file_location("policy_pd_contributors", path)
        assert spec is not None
        assert spec.loader is not None
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    finally:
        sys.path.pop(0)


def write_inputs(root: Path, policy: TraitKeyPolicy) -> Path:
    """Hand-authored counts: literal608 has one source; union has four incl. strain605."""
    root.mkdir()
    edges = root / "merged-kg_edges.tsv"
    nodes = root / "merged-kg_nodes.tsv"
    edges.write_text(
        "id\tsubject\tpredicate\tobject\n"
        "a\tNCBITaxon:1\tbiolink:has_phenotype\tMETPO:1000605\n"
        "b\tNCBITaxon:2\tbiolink:has_phenotype\tMETPO:1000608\n"
        "c\tNCBITaxon:3\tbiolink:has_phenotype\tMETPO:1000612\n"
        "d\tkgmicrobe.strain:4\tbiolink:subclass_of\tNCBITaxon:4\n"
        "e\tkgmicrobe.strain:4\tbiolink:has_phenotype\tMETPO:1000605\n"
        + "".join(f"p{i}\tNCBITaxon:{i}\tbiolink:subclass_of\tNCBITaxon:9\n" for i in range(1, 5)),
        encoding="utf-8",
    )
    nodes.write_text("id\tname\nNCBITaxon:9\tGenus\n", encoding="utf-8")
    count = 1 if policy == "raw_object_v1" else 4
    values = f"{count}\t0\t0\t0\t0\t0\n"
    (root / "study_trait_profiles_wide.tsv").write_text(
        "study_id\texperiment_id\tinc_dec\t"
        + "\t".join(COLUMNS)
        + "\n"
        + "S\tE\tinc\t"
        + values
        + "S\tE\tdec\t"
        + values,
        encoding="utf-8",
    )
    (root / "taxa_trait_matrix.tsv").write_text(
        "taxon_id\ttaxon_name\ttaxonomic_rank\t"
        + "\t".join(COLUMNS)
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
    protocol = root / "protocol.json"
    protocol.write_text(
        json.dumps({"kg_md5": {p.name: digest(p, "md5") for p in (edges, nodes)}}), encoding="utf-8"
    )
    return protocol


@pytest.mark.parametrize(
    ("policy", "expected"), [("raw_object_v1", 1), ("facultative_union_v1", 4)]
)
def test_real_ledger_cli_records_and_reconstructs_policy(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    policy: TraitKeyPolicy,
    expected: int,
) -> None:
    inputs = tmp_path / policy
    protocol = write_inputs(inputs, policy)
    output = tmp_path / "output"
    # Influence statistics are a separate consumer; input parsing, KG verification,
    # graph reconstruction, reconciliation and all ledger exports remain real.
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
            "--trait-key-policy",
            policy,
        ],
    )
    taxon_contributions.main()
    manifest = json.loads((output / "manifest.json").read_text())
    assert manifest["completed"] is True
    assert manifest["trait_key_policy"] == policy
    assert manifest["mismatches"] == 0
    with (output / "reconciliation.tsv").open() as stream:
        assert list(csv.DictReader(stream, delimiter="\t")) == []
    import gzip

    with gzip.open(output / "source_contributions.tsv.gz", "rt") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert len(rows) == expected
    assert {row["trait_key"] for row in rows} == {TRAIT}


@pytest.mark.parametrize(
    ("policy", "expected"), [("raw_object_v1", 1), ("facultative_union_v1", 4)]
)
def test_pd_ledger_reconstructs_selected_policy_and_refuses_the_other(
    tmp_path: Path,
    contributors: ModuleType,
    policy: TraitKeyPolicy,
    expected: int,
) -> None:
    inputs = tmp_path / policy
    protocol = write_inputs(inputs, policy)
    _, graph, ledger = contributors.load_contribution_ledger(
        inputs,
        inputs,
        protocol,
        "phenotype",
        trait_key_policy=policy,
    )
    assert ledger.mismatches == []
    assert ledger.reconstructed["NCBITaxon:9"][TRAIT] == expected
    assert len(graph.traits) == expected
    wrong = "facultative_union_v1" if policy == "raw_object_v1" else "raw_object_v1"
    with pytest.raises(ValueError, match="counts differ"):
        contributors.load_contribution_ledger(
            inputs,
            inputs,
            protocol,
            "phenotype",
            trait_key_policy=wrong,
        )


@pytest.mark.parametrize("policy", ["raw_object_v1", "facultative_union_v1"])
def test_pd_cli_forwards_policy_before_analysis(
    monkeypatch: pytest.MonkeyPatch,
    contributors: ModuleType,
    policy: TraitKeyPolicy,
) -> None:
    class Reached(Exception):
        pass

    monkeypatch.setattr(contributors, "load_profile_metadata", lambda: None)
    monkeypatch.setattr(contributors, "require_current_trait_taxa_evidence", lambda *_a: None)
    monkeypatch.setattr(contributors, "read_tsv", lambda *_a: None)
    monkeypatch.setattr(contributors, "trait_directions", lambda *_a: None)
    monkeypatch.setattr(contributors, "check_trait_intervals", lambda *_a: None)

    def inspect(*_args: object, **kwargs: object) -> None:
        assert kwargs["trait_key_policy"] == policy
        raise Reached

    monkeypatch.setattr(contributors, "load_contribution_ledger", inspect)
    monkeypatch.setattr(sys, "argv", ["pd", "--trait-key-policy", policy])
    with pytest.raises(Reached):
        contributors.main()


@pytest.mark.parametrize(
    ("policy", "expected"), [("raw_object_v1", 1), ("facultative_union_v1", 4)]
)
def test_pd_main_persists_consumed_policy_and_input_binding(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    contributors: ModuleType,
    policy: TraitKeyPolicy,
    expected: int,
) -> None:
    """Keep the real ledger and manifest path; downstream tables and plots are out of scope."""
    inputs = tmp_path / policy
    protocol = write_inputs(inputs, policy)
    output = tmp_path / "pd-output"
    output.mkdir()
    metadata = pd.DataFrame(
        [
            {
                "study_id": "S",
                "experiment_id": "E",
                "inc_dec": direction,
                "profile_id": f"S_E_{direction}",
                "condition": contributors.PD_CONDITION,
                "body_site_group": contributors.FECES,
            }
            for direction in ("inc", "dec")
        ]
    )
    traits = pd.DataFrame([{"trait_key": TRAIT, "trait_slug": "facultatively_aerobic"}])
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
    coverage_seen: list[pd.DataFrame] = []

    def fixture_tables(
        counts: pd.DataFrame,
        coverage_counts: pd.DataFrame,
        selected_traits: pd.DataFrame,
        names: dict[str, str],
    ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        # The synthetic genus has no species/family resolution. Its four cells
        # (two ranks by two directions) must nevertheless carry the real sources.
        assert counts.empty
        assert set(coverage_counts["trait_key"]) == {TRAIT}
        assert set(coverage_counts["status"]) == {"unresolved_rank"}
        assert set(zip(coverage_counts["rank"], coverage_counts["direction"], strict=True)) == {
            (rank, direction) for rank in ("species", "family") for direction in ("inc", "dec")
        }
        assert coverage_counts["count"].tolist() == [expected] * 4
        pd.testing.assert_frame_equal(selected_traits, traits)
        assert names["NCBITaxon:9"] == "Genus"
        coverage_seen.append(coverage_counts.copy())
        return cells.copy(), coverage_counts.copy(), selected.copy()

    for name, value in {
        "REPO_ROOT": tmp_path,
        "KG_DIR": inputs,
        "DERIVED_DIR": inputs,
        "KG_PROTOCOL_PATH": protocol,
        "OUTPUT_DIR": output,
    }.items():
        monkeypatch.setattr(contributors, name, value)
    # Prerequisite figures are unrelated to policy loading. The KG reader,
    # reconciliation, profile selection, contribution counting, input hashing,
    # argparse and manifest serialization below all remain real.
    monkeypatch.setattr(contributors, "load_profile_metadata", lambda: metadata.copy())
    monkeypatch.setattr(contributors, "require_current_trait_taxa_evidence", lambda *_a: None)
    monkeypatch.setattr(contributors, "read_tsv", lambda *_a: pd.DataFrame())
    monkeypatch.setattr(contributors, "trait_directions", lambda *_a: pd.DataFrame())
    monkeypatch.setattr(contributors, "check_trait_intervals", lambda *_a: None)
    monkeypatch.setattr(contributors, "trait_frame", lambda *_a: traits.copy())
    monkeypatch.setattr(contributors, "build_contributor_tables", fixture_tables)
    monkeypatch.setattr(contributors, "augment_selected_taxa", lambda frame, _metadata: frame)
    monkeypatch.setattr(contributors, "build_direction_differences", lambda *_a: pd.DataFrame())
    monkeypatch.setattr(contributors, "select_large_difference_taxa", lambda *_a: selected.copy())
    monkeypatch.setattr(
        contributors,
        "build_trait_cocarriage_tables",
        lambda *_a: (pd.DataFrame(), pd.DataFrame(), pd.DataFrame()),
    )
    monkeypatch.setattr(contributors, "plot_contributor_heatmap", lambda *_a: None)
    monkeypatch.setattr(contributors, "plot_trait_cocarriage_network", lambda *_a: None)
    monkeypatch.setattr(contributors, "write_note", lambda *_a: None)
    monkeypatch.setattr(sys, "argv", ["pd", "--trait-key-policy", policy])

    contributors.main()

    assert len(coverage_seen) == 1
    manifest = json.loads((output / f"{contributors.STEM}_manifest.json").read_text())
    assert manifest["trait_key_policy"] == policy
    assert manifest["strain_rollup"] == "phenotype"
    assert manifest["pd_feces_profiles"] == 2
    assert manifest["ledger_profiles"] == 2
    assert manifest["ledger_increased_profiles"] == 1
    assert manifest["ledger_decreased_profiles"] == 1
    assert manifest["ledger_studies"] == 1
    assert manifest["traits"] == traits.to_dict(orient="records")
    assert manifest["kg_md5_verified_against"] == protocol.relative_to(tmp_path).as_posix()
    assert {name: "".join(chunks) for name, chunks in manifest["input_sha256_chunks"].items()} == {
        name: digest(inputs / name, "sha256") for name in contributors.LEDGER_INPUTS
    }


@pytest.mark.parametrize("policy", ["raw_object_v1", "facultative_union_v1"])
def test_taxon_contributions_main_forwards_explicit_trait_key_policy_to_load_graph(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    policy: TraitKeyPolicy,
) -> None:
    class Reached(Exception):
        pass

    inputs = tmp_path / policy
    protocol = write_inputs(inputs, policy)
    output = tmp_path / "output"
    captured_policies: list[object] = []

    def inspect_load_graph(*_args: object, **kwargs: object) -> None:
        captured_policies.append(kwargs.get("trait_key_policy"))
        raise Reached

    monkeypatch.setattr(taxon_contributions, "load_graph", inspect_load_graph)
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
            "--trait-key-policy",
            policy,
        ],
    )
    with pytest.raises(Reached):
        taxon_contributions.main()
    captured_policy = captured_policies[0] if len(captured_policies) == 1 else None
    assert captured_policy == policy
