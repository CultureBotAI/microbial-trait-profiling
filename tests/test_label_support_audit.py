"""Synthetic source reconciliation and scientific-estimand regression tests for #773."""

from __future__ import annotations

import ast
import importlib
import json
from pathlib import Path

import pandas as pd
import pytest

from mtp.label_support_audit import (
    AUDIT_TRAITS,
    MADIN_SOURCE,
    Scope,
    bind_occurrences,
    load_kg_provenance_graph,
    load_madin_subsources,
    main,
    metadata_source_status,
)
from mtp.paths import REPO_ROOT, file_sha256

CELL = "METPO:2000103|GO:0030245"
ANA = "METPO:1000603"
FAC = "METPO:1000608"
COCCUS = "METPO:1000668"


@pytest.fixture
def inputs(tmp_path: Path) -> dict[str, Path]:
    """Two synthetic genera with known carrier sources; three pairs in two studies."""
    paths = {
        name: tmp_path / f"{name}.tsv" for name in ("taxa", "matrix", "profiles", "edges", "nodes")
    }
    paths["madin"] = tmp_path / "madin.csv"
    records = []
    for child, parent in ((11, 10), (12, 10), (21, 20), (22, 20)):
        records.append(
            (f"NCBITaxon:{child}", "biolink:subclass_of", f"NCBITaxon:{parent}", "infores:ncbi")
        )
    for taxon, trait, source in (
        (11, CELL, MADIN_SOURCE),
        (11, ANA, MADIN_SOURCE),
        (12, ANA, MADIN_SOURCE),
        (21, CELL, MADIN_SOURCE),
        (21, FAC, MADIN_SOURCE),
        (21, COCCUS, MADIN_SOURCE),
        (22, CELL, MADIN_SOURCE),
        (22, ANA, "infores:bactotraits"),
    ):
        predicate, obj = trait.split("|", 1) if "|" in trait else ("biolink:has_phenotype", trait)
        records.append((f"NCBITaxon:{taxon}", predicate, obj, source))
    records.extend(
        [
            ("kgmicrobe.strain:fixture", "biolink:subclass_of", "NCBITaxon:11", "infores:bacdive"),
            ("kgmicrobe.strain:fixture", "biolink:has_phenotype", COCCUS, "infores:bacdive"),
        ]
    )
    edges = pd.DataFrame(
        records, columns=["subject", "predicate", "object", "primary_knowledge_source"]
    )
    edges.insert(0, "id", [f"edge-{index}" for index in range(len(edges))])
    edges.to_csv(paths["edges"], sep="\t", index=False)
    pd.DataFrame(
        {
            "id": [f"NCBITaxon:{i}" for i in (10, 11, 12, 20, 21, 22)],
            "name": ["Examplegenus A", "A alpha", "A beta", "Examplegenus B", "B alpha", "B beta"],
        }
    ).to_csv(paths["nodes"], sep="\t", index=False)
    matrix_rows = []
    for taxon, counts in ((10, {CELL: 1, ANA: 2}), (20, {CELL: 2, ANA: 1, FAC: 1, COCCUS: 1})):
        matrix_rows.append(
            {
                "taxon_id": f"NCBITaxon:{taxon}",
                "taxon_name": f"Examplegenus {taxon}",
                "taxonomic_rank": "genus",
                **dict.fromkeys(AUDIT_TRAITS, 0),
                **counts,
            }
        )
    matrix = pd.DataFrame(matrix_rows)
    matrix.to_csv(paths["matrix"], sep="\t", index=False)
    count_lookup = matrix.set_index("taxon_id")
    taxa, profiles = [], []
    for study, experiment, inc in (("s1", "e1", 20), ("s1", "e2", 10), ("s2", "e1", 20)):
        for direction, taxon in (("inc", inc), ("dec", 10)):
            tid = f"NCBITaxon:{taxon}"
            taxa.append(
                {
                    "study_id": study,
                    "experiment_id": experiment,
                    "direction": "increased" if direction == "inc" else "decreased",
                    "taxon_id": tid,
                }
            )
            profiles.append(
                {
                    "study_id": study,
                    "experiment_id": experiment,
                    "inc_dec": direction,
                    **{feature: count_lookup.loc[tid, feature] for feature in AUDIT_TRAITS},
                }
            )
    pd.DataFrame(taxa).to_csv(paths["taxa"], sep="\t", index=False)
    pd.DataFrame(profiles).to_csv(paths["profiles"], sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "tax_id": 11,
                "data_source": "gold",
                "pathways": "cellulose_degradation",
                "metabolism": "anaerobic",
            },
            {"tax_id": 12, "data_source": "faprotax", "pathways": "", "metabolism": "anaerobic"},
            {
                "tax_id": 21,
                "data_source": "faprotax",
                "pathways": "cellulose_degradation",
                "metabolism": "facultative",
                "cell_shape": "coccus",
            },
            {"tax_id": 22, "data_source": "patric", "pathways": "cellulose_degradation"},
        ]
    ).to_csv(paths["madin"], index=False)
    return paths


def run_audit(inputs: dict[str, Path], output: Path, *, madin: bool = True) -> None:
    arguments = [
        argument
        for name in ("taxa", "matrix", "profiles", "edges", "nodes")
        for argument in (f"--{name}", str(inputs[name]))
    ]
    if madin:
        arguments.extend(["--madin-csv", str(inputs["madin"])])
    assert main([*arguments, "--output-dir", str(output), "--strain-rollup", "off"]) == 0


@pytest.mark.parametrize("scope", ["off", "phenotype"])
@pytest.mark.parametrize("subject", ["NCBITaxon:11", "kgmicrobe.strain:fixture"])
def test_audit_graph_guards_assertion_negation_even_for_excluded_strains(
    inputs: dict[str, Path], scope: Scope, subject: str
) -> None:
    edges = pd.read_csv(inputs["edges"], sep="\t", dtype=str).fillna("")
    edges["meta"] = ""
    edges.loc[len(edges)] = [
        "unreviewed-polarity",
        subject,
        "biolink:has_phenotype",
        ANA,
        "fixture",
        '{"is_negated": true}',
    ]
    edges.to_csv(inputs["edges"], sep="\t", index=False)
    with pytest.raises(ValueError, match=r"unreviewed-polarity.*negat.*review"):
        load_kg_provenance_graph(inputs["edges"], inputs["nodes"], strain_rollup=scope)


@pytest.mark.parametrize("scope", ["off", "phenotype"])
def test_audit_graph_guards_observed_external_predicate_metadata(
    inputs: dict[str, Path], scope: Scope
) -> None:
    with inputs["edges"].open("a", encoding="utf-8") as handle:
        handle.write(f"external\tNCBITaxon:11\texample:negative_relation\t{ANA}\tfixture\n")
    with inputs["nodes"].open("a", encoding="utf-8") as handle:
        handle.write("example:negative_relation\tdoes not exhibit\n")
    with pytest.raises(ValueError, match=r"external.*negat.*review"):
        load_kg_provenance_graph(inputs["edges"], inputs["nodes"], strain_rollup=scope)


def test_sidecar_records_sources_reconciliations_and_study_means(
    inputs: dict[str, Path], tmp_path: Path
) -> None:
    output = tmp_path / "audit"
    run_audit(inputs, output)
    support = pd.read_csv(output / "taxon_trait_label_support.tsv", sep="\t")
    focal = support[(support["taxon_id"] == "NCBITaxon:10") & (support["trait_key"] == CELL)].iloc[
        0
    ]
    assert focal["carrier_count"] == 1
    assert not focal["carrier_at_reported_taxon"]
    assert focal["annotated_subtree_taxa"] == 2
    assert focal["annotated_subtree_fraction"] == pytest.approx(0.5)
    assert focal["source_status"] == "gold_patric_only"
    sensitivity = pd.read_csv(output / "support_threshold_sensitivity.tsv", sep="\t")
    expected_focal = {
        CELL,
        FAC,
        COCCUS,
        "METPO:1000683",
        "METPO:1000623",
        "METPO:2000103|pathways:arsenate_reduction",
    }
    assert set(sensitivity["trait_key"]) == expected_focal
    assert len(sensitivity) == 42
    assert sensitivity.groupby("trait_key")["rule_id"].nunique().eq(7).all()
    row = sensitivity[
        (sensitivity["trait_key"] == CELL)
        & (sensitivity["rule_id"] == "carriers_per_reported_taxon")
    ].iloc[0]
    assert row["n_experiments"] == 3
    assert row["n_studies"] == 2
    assert row["study_mean_increased_minus_decreased"] == pytest.approx(0.75)
    assert row["study_sign_p_value"] == pytest.approx(0.5)
    assert row["effect_units"] == "carriers_per_reported_taxon"
    native = sensitivity[
        (sensitivity["trait_key"] == CELL)
        & (sensitivity["rule_id"] == "native_annotation_fraction")
    ].iloc[0]
    assert native["study_mean_increased_minus_decreased"] == pytest.approx(0.05)
    ties = sensitivity[sensitivity["trait_key"] == "METPO:1000623"]
    assert ties["study_sign_p_value"].isna().all()
    census = pd.read_csv(output / "coasserted_label_census.tsv", sep="\t")
    dual = census[
        (census["label_pair"] == "aerobe_guild_and_anaerobe_guild") & (census["rank"] == "genus")
    ].iloc[0]
    assert dual["taxa_with_either"] == 2
    assert dual["taxa_with_both"] == 1
    provenance = json.loads((output / "provenance.json").read_text())
    assert provenance["scope"]["profile_counts_reconciled"] is True
    assert provenance["inputs"]["madin"]["sha256"] == file_sha256(inputs["madin"])
    for name, entry in provenance["outputs"].items():
        assert entry["sha256"] == file_sha256(output / name)
    assert "table1_trait_semantic_audit.tsv" not in provenance["outputs"]
    assert "/Users/" not in (output / "provenance.json").read_text()


def test_optional_madin_is_unknown_not_negative(inputs: dict[str, Path], tmp_path: Path) -> None:
    run_audit(inputs, tmp_path / "audit", madin=False)
    support = pd.read_csv(tmp_path / "audit" / "taxon_trait_label_support.tsv", sep="\t")
    assert support.loc[support["trait_key"] == CELL, "source_status"].eq("unknown").all()
    with pytest.raises(FileNotFoundError):
        load_madin_subsources(tmp_path / "absent.csv")
    assert metadata_source_status({MADIN_SOURCE}, set(), 1) == "unknown"
    assert metadata_source_status({MADIN_SOURCE}, {"gold"}, 0) == "gold_patric_only"
    assert metadata_source_status({MADIN_SOURCE}, {"gold"}, 1) == "unknown"
    assert metadata_source_status({MADIN_SOURCE}, {"faprotax"}, 0) == "includes_other_sources"
    sensitivity = pd.read_csv(tmp_path / "audit" / "support_threshold_sensitivity.tsv", sep="\t")
    cell = sensitivity.loc[sensitivity["trait_key"] == CELL].set_index("rule_id")
    assert cell.loc[
        "exclude_known_gold_patric_only", "study_mean_increased_minus_decreased"
    ] == pytest.approx(0.75)
    assert (
        cell.loc["exclude_known_gold_patric_only", "study_mean_increased_minus_decreased"]
        == (cell.loc["carriers_per_reported_taxon", "study_mean_increased_minus_decreased"])
    )


def test_empty_profile_columns_may_be_dropped(inputs: dict[str, Path], tmp_path: Path) -> None:
    profiles = pd.read_csv(inputs["profiles"], sep="\t")
    empty = [key for key in profiles.columns[3:] if not profiles[key].gt(0).any()]
    profiles.drop(columns=empty).to_csv(inputs["profiles"], sep="\t", index=False)
    run_audit(inputs, tmp_path / "audit")


def test_sensitivity_source_counts_use_only_bound_reported_taxa(
    inputs: dict[str, Path], tmp_path: Path
) -> None:
    matrix = pd.read_csv(inputs["matrix"], sep="\t")
    unused = matrix.iloc[0].copy()
    unused["taxon_id"] = "NCBITaxon:11"
    unused["taxon_name"] = "A alpha"
    unused["taxonomic_rank"] = "species"
    unused[ANA] = 1
    pd.concat([matrix, unused.to_frame().T]).to_csv(inputs["matrix"], sep="\t", index=False)
    output = tmp_path / "audit"
    run_audit(inputs, output, madin=False)
    sensitivity = pd.read_csv(output / "support_threshold_sensitivity.tsv", sep="\t")
    rows = sensitivity.loc[sensitivity["trait_key"] == CELL]
    assert rows["source_unknown_taxa"].eq(2).all()
    census = pd.read_csv(output / "coasserted_label_census.tsv", sep="\t")
    assert census.loc[census["rank"] == "ALL", "unique_reported_taxa"].eq(2).all()
    assert census.loc[census["rank"] == "species", "unique_reported_taxa"].eq(0).all()


def test_named_kg_columns_and_shared_phenotype_scope(inputs: dict[str, Path]) -> None:
    edges = pd.read_csv(inputs["edges"], sep="\t")
    edges.insert(4, "unused_evidence_field", "ignored")
    edges.to_csv(inputs["edges"], sep="\t", index=False)
    graph = load_kg_provenance_graph(inputs["edges"], inputs["nodes"])
    assert graph["sources"][CELL]["NCBITaxon:11"] == {MADIN_SOURCE}
    assert COCCUS not in graph["traits"]["NCBITaxon:11"]
    phenotype = load_kg_provenance_graph(
        inputs["edges"], inputs["nodes"], strain_rollup="phenotype"
    )
    assert COCCUS in phenotype["traits"]["NCBITaxon:11"]
    assert phenotype["sources"][COCCUS]["NCBITaxon:11"] == {"infores:bacdive"}


def test_known_bare_cr_kg_header_and_agent_record_are_reconstructed(
    inputs: dict[str, Path], tmp_path: Path
) -> None:
    records = pd.read_csv(inputs["edges"], sep="\t").to_dict("records")
    header = (
        "id\tsubject\tpredicate\tobject\trelation\tagent_type\tagent_type\r"
        "\tknowledge_level\tknowledge_source\tmeta\tprimary_knowledge_source\n"
    )
    lines = []
    for index, row in enumerate(records):
        agent = "manual_validation_of_automated_agent\r" if index == 4 else "manual_agent"
        lines.append(
            "\t".join(
                [
                    row["id"],
                    row["subject"],
                    row["predicate"],
                    row["object"],
                    "",
                    "automated_agent",
                    agent,
                    "knowledge_assertion",
                    "",
                    "",
                    row["primary_knowledge_source"],
                ]
            )
            + "\n"
        )
    inputs["edges"].write_bytes((header + "".join(lines)).encode())
    output = tmp_path / "audit"
    run_audit(inputs, output)
    provenance = json.loads((output / "provenance.json").read_text())
    serialization = provenance["kg_edge_serialization"]
    assert serialization["records"] == len(records)
    assert serialization["header_continuations_reconstructed"] == 1
    assert serialization["agent_field_continuations_reconstructed"] == 1
    evidence = pd.read_csv(output / "carrier_source_evidence.tsv", sep="\t")
    assert evidence.loc[evidence["trait_key"] == CELL, "madin_resolution"].eq("resolved").all()


def test_pipe_joined_primary_sources_resolve_madin(inputs: dict[str, Path], tmp_path: Path) -> None:
    edges = pd.read_csv(inputs["edges"], sep="\t")
    selected = (edges["subject"] == "NCBITaxon:11") & (edges["object"] == "GO:0030245")
    edges.loc[selected, "primary_knowledge_source"] = f"{MADIN_SOURCE}|infores:bactotraits"
    edges.to_csv(inputs["edges"], sep="\t", index=False)
    output = tmp_path / "audit"
    run_audit(inputs, output)
    evidence = pd.read_csv(output / "carrier_source_evidence.tsv", sep="\t")
    row = evidence.loc[
        (evidence["reported_taxon_id"] == "NCBITaxon:10") & (evidence["trait_key"] == CELL)
    ].iloc[0]
    assert row["madin_resolution"] == "resolved"
    assert row["madin_subsources"] == "gold"
    assert set(row["primary_sources"].split(";")) == {MADIN_SOURCE, "infores:bactotraits"}


def test_all_support_rules_share_zero_total_filtered_pairs(
    inputs: dict[str, Path], tmp_path: Path
) -> None:
    matrix = pd.read_csv(inputs["matrix"], sep="\t")
    unannotated = matrix.iloc[0].copy()
    unannotated["taxon_id"] = "NCBITaxon:99"
    unannotated["taxon_name"] = "Unannotated fixture taxon"
    unannotated.loc[list(AUDIT_TRAITS)] = 0
    excluded_only = unannotated.copy()
    excluded_only["taxon_id"] = "NCBITaxon:30"
    excluded_only["taxon_name"] = "Excluded-pair fixture taxon"
    excluded_only[CELL] = 1
    pd.concat([matrix, unannotated.to_frame().T, excluded_only.to_frame().T]).to_csv(
        inputs["matrix"], sep="\t", index=False
    )
    with inputs["edges"].open("a", encoding="utf-8") as handle:
        handle.write("excluded-only\tNCBITaxon:30\tMETPO:2000103\tGO:0030245\t\n")
    with inputs["nodes"].open("a", encoding="utf-8") as handle:
        handle.write("NCBITaxon:30\tExcluded-pair fixture taxon\n")
    profiles = pd.read_csv(inputs["profiles"], sep="\t")
    occurrences = pd.read_csv(inputs["taxa"], sep="\t")
    profile_rows, occurrence_rows = [], []
    for study, experiment, both_zero in (("s3", "e1", True), ("s1", "e3", False)):
        for direction in ("inc", "dec"):
            zero = both_zero or direction == "dec"
            row = profiles.iloc[0].copy()
            row.loc[["study_id", "experiment_id", "inc_dec"]] = [study, experiment, direction]
            row.loc[list(AUDIT_TRAITS)] = 0
            if not zero:
                row[CELL] = 1
            profile_rows.append(row)
            occurrence_rows.append(
                {
                    "study_id": study,
                    "experiment_id": experiment,
                    "direction": "increased" if direction == "inc" else "decreased",
                    "taxon_id": "NCBITaxon:99" if zero else "NCBITaxon:30",
                }
            )
    pd.concat([profiles, pd.DataFrame(profile_rows)]).to_csv(
        inputs["profiles"], sep="\t", index=False
    )
    pd.concat([occurrences, pd.DataFrame(occurrence_rows)]).to_csv(
        inputs["taxa"], sep="\t", index=False
    )
    output = tmp_path / "audit"
    run_audit(inputs, output)
    sensitivity = pd.read_csv(output / "support_threshold_sensitivity.tsv", sep="\t")
    assert sensitivity["n_experiments"].eq(3).all()
    assert sensitivity["n_studies"].eq(2).all()
    assert sensitivity["n_excluded_zero_annotation_pairs"].eq(2).all()
    assert sensitivity["n_input_paired_experiments"].eq(5).all()
    cell = sensitivity.loc[sensitivity["trait_key"] == CELL].set_index("rule_id")
    assert cell["source_unknown_taxa"].eq(0).all()
    assert cell.loc["carriers_per_reported_taxon", "study_mean_increased_minus_decreased"] == (
        pytest.approx(0.75)
    )


def test_legacy_label_support_module_is_canonical_alias() -> None:
    assert importlib.import_module("bsdbng.label_support_audit") is importlib.import_module(
        "mtp.label_support_audit"
    )
    canonical = ast.parse((REPO_ROOT / "src/mtp/label_support_audit.py").read_text())
    shim = ast.parse((REPO_ROOT / "src/bsdbng/label_support_audit.py").read_text())
    declarations: set[str] = set()
    for node in canonical.body:
        names: list[str] = []
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            names = [node.name]
        elif isinstance(node, ast.Assign):
            names = [target.id for target in node.targets if isinstance(target, ast.Name)]
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            names = [node.target.id]
        declarations.update(name for name in names if not name.startswith("_"))
    exports = next(
        ast.literal_eval(node.value)
        for node in shim.body
        if isinstance(node, ast.Assign)
        and any(isinstance(target, ast.Name) and target.id == "__all__" for target in node.targets)
    )
    imported = {
        name.asname or name.name
        for node in shim.body
        if isinstance(node, ast.ImportFrom) and node.module == "mtp.label_support_audit"
        for name in node.names
    }
    assert declarations <= set(exports) <= imported
    assert len(exports) == len(set(exports))


@pytest.mark.parametrize(
    "problem",
    ["missing_side", "profile_counts", "carrier_counts", "missing_feature", "negative_count"],
)
def test_mismatched_sources_fail_before_writing(
    inputs: dict[str, Path], tmp_path: Path, problem: str
) -> None:
    if problem == "missing_side":
        table = pd.read_csv(inputs["taxa"], sep="\t").iloc[1:]
        table.to_csv(inputs["taxa"], sep="\t", index=False)
    elif problem == "carrier_counts":
        edges = pd.read_csv(inputs["edges"], sep="\t")
        edges = edges[~((edges["subject"] == "NCBITaxon:11") & (edges["object"] == "GO:0030245"))]
        edges.to_csv(inputs["edges"], sep="\t", index=False)
    else:
        profiles = pd.read_csv(inputs["profiles"], sep="\t")
        if problem == "missing_feature":
            profiles = profiles.drop(columns=ANA)
        else:
            profiles.loc[0, CELL] = -1 if problem == "negative_count" else 9
        profiles.to_csv(inputs["profiles"], sep="\t", index=False)
    with pytest.raises(ValueError, match=r"missing|disagree|integer"):
        run_audit(inputs, tmp_path / "audit")
    assert not (tmp_path / "audit").exists()


def test_extra_occurrences_are_scoped_to_profiles(inputs: dict[str, Path]) -> None:
    taxa = pd.read_csv(inputs["taxa"], sep="\t")
    extra = taxa.iloc[:1].assign(study_id="outside_population")
    pd.concat([taxa, extra]).to_csv(inputs["taxa"], sep="\t", index=False)
    bound, facts = bind_occurrences(
        inputs["taxa"],
        pd.read_csv(inputs["matrix"], sep="\t"),
        pd.read_csv(inputs["profiles"], sep="\t"),
    )
    assert len(bound) == 6
    assert facts["input_occurrence_rows"] == 7


def test_sidecar_never_overwrites_existing_results(inputs: dict[str, Path], tmp_path: Path) -> None:
    output = tmp_path / "audit"
    output.mkdir()
    previous = output / "provenance.json"
    previous.write_text("existing provenance\n")
    with pytest.raises(FileExistsError):
        run_audit(inputs, output)
    assert previous.read_text() == "existing provenance\n"
