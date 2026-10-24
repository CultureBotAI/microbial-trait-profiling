"""Count KG phenotype co-assertions without rebuilding the pipeline (#846)."""

from __future__ import annotations

import argparse
import csv
import io
import json
import platform
from collections import defaultdict
from collections.abc import Sequence
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from mtp import taxa_trait_matrix as matrix_api
from mtp import trait_guilds
from mtp.paths import REPO_ROOT, file_sha256
from mtp.workflow import kg_inputs
from mtp.workflow.runtime import repository_state


@dataclass(frozen=True)
class Pair:
    """Literal classes; overlap does not establish physiological contradiction."""

    id: str
    left: tuple[str, ...]
    right: tuple[str, ...]


PAIRS = (
    Pair("aerobic_class_and_anaerobic_class", ("METPO:1000602",), ("METPO:1000603",)),
    Pair("obligately_aerobic_and_anaerobic", ("METPO:1000606",), ("METPO:1000603",)),
    Pair("mesophilic_and_psychrophilic", ("METPO:1000615",), ("METPO:1000614",)),
    Pair("facultatively_anaerobic_and_obligately_aerobic", ("METPO:1000605",), ("METPO:1000606",)),
    Pair(
        "aerobe_guild_and_anaerobe_guild",
        trait_guilds.OXYGEN_GUILDS["aerobe"],
        trait_guilds.OXYGEN_GUILDS["anaerobe"],
    ),
)
TARGET_KEYS = frozenset(key for pair in PAIRS for key in (*pair.left, *pair.right))
SUMMARY_FIELDS = (
    "scope",
    "pair_id",
    "left_keys",
    "right_keys",
    "left_taxon_nodes",
    "right_taxon_nodes",
    "either_taxon_nodes",
    "coasserted_taxon_nodes",
    "direct_taxon_has_both_nodes",
    "same_strain_has_both_nodes",
    "requires_union_across_subjects_nodes",
    "coasserted_strain_nodes",
    "unlinked_coasserted_strain_nodes",
)
TAXON_FIELDS = (
    "scope",
    "pair_id",
    "taxon_id",
    "taxon_name",
    "node_metadata_status",
    "left_keys",
    "right_keys",
    "left_direct_keys",
    "right_direct_keys",
    "left_strain_ids",
    "right_strain_ids",
    "same_strain_ids",
    "direct_taxon_has_both",
    "same_strain_has_both",
    "requires_union_across_subjects",
)
STRAIN_FIELDS = (
    "scope",
    "pair_id",
    "strain_id",
    "left_keys",
    "right_keys",
    "mapped_taxon_ids",
    "mapped_taxon_count",
)
EVIDENCE_FIELDS = (
    "scope",
    "pair_id",
    "taxon_id",
    "side",
    "trait_key",
    "assertion_subject",
    "origin",
    "edge_id",
    "predicate",
    "object",
    "primary_source",
)


@dataclass(frozen=True)
class Assertion:
    subject: str
    predicate: str
    obj: str
    key: str
    edge_id: str
    sources: tuple[str, ...]


@dataclass
class Graph:
    taxon_traits: dict[str, set[str]]
    strain_traits: dict[str, set[str]]
    strain_taxa: dict[str, set[str]]
    assertions: dict[str, list[Assertion]]
    serialization: dict[str, object]


def read_graph(edges: Path) -> Graph:
    """Retain target assertions and strain links through the shared KG reader."""
    taxa: defaultdict[str, set[str]] = defaultdict(set)
    strains: defaultdict[str, set[str]] = defaultdict(set)
    links: defaultdict[str, set[str]] = defaultdict(set)
    assertions: defaultdict[str, list[Assertion]] = defaultdict(list)
    serialization: dict[str, object] = {}
    for row in kg_inputs.iter_kg_edges(edges, serialization):
        if "primary_knowledge_source" not in row:
            raise ValueError("KG edges require primary_knowledge_source for source evidence")
        subject, predicate, obj = (row[key].strip() for key in ("subject", "predicate", "object"))
        strain = subject.startswith(matrix_api.STRAIN_PREFIX)
        if (
            strain
            and predicate == matrix_api.SUBCLASS_OF
            and obj.startswith(matrix_api.TAXON_PREFIX)
        ):
            links[subject].add(obj)
            continue
        key = matrix_api.make_trait_key(predicate, obj)
        if key not in TARGET_KEYS:
            continue
        if subject.startswith(matrix_api.TAXON_PREFIX):
            taxa[subject].add(key)
        elif strain:
            strains[subject].add(key)
        else:
            continue
        sources = tuple(
            sorted(
                {
                    part.strip()
                    for part in row["primary_knowledge_source"].split("|")
                    if part.strip()
                }
            )
        ) or ("unknown",)
        assertions[subject].append(
            Assertion(subject, predicate, obj, key, row.get("id", ""), sources)
        )
    if not serialization["records"]:
        raise ValueError("KG edge export has no records")
    return Graph(dict(taxa), dict(strains), dict(links), dict(assertions), serialization)


def read_names(path: Path, wanted: set[str]) -> dict[str, str]:
    """Read named columns; missing taxon names stay missing rather than invented."""
    names = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"id", "name"} <= set(reader.fieldnames or ()):
            raise ValueError("KG nodes require id and name columns")
        for row in reader:
            if row["id"] in wanted and row["name"].strip():
                names[row["id"]] = row["name"].strip()
    if missing := TARGET_KEYS - names.keys():
        raise ValueError(f"Missing audited class metadata: {sorted(missing)}")
    return names


def _joined(values: Sequence[str] | set[str]) -> str:
    return json.dumps(sorted(values), ensure_ascii=False)


def _evidence(
    graph: Graph,
    scope: str,
    pair: Pair,
    taxon: str,
    strains: set[str],
    admitted: dict[str, set[str]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    subjects = ([taxon] if taxon else []) + sorted(strains)
    for subject in subjects:
        allowed = (
            graph.taxon_traits.get(subject, set())
            if subject == taxon
            else admitted.get(subject, set())
        )
        for assertion in graph.assertions.get(subject, ()):
            if assertion.key not in allowed:
                continue
            for side, keys in (("left", pair.left), ("right", pair.right)):
                if assertion.key in keys:
                    rows.extend(
                        {
                            "scope": scope,
                            "pair_id": pair.id,
                            "taxon_id": taxon,
                            "side": side,
                            "trait_key": assertion.key,
                            "assertion_subject": subject,
                            "origin": "direct_taxon"
                            if subject == taxon
                            else ("mapped_strain" if taxon else "unmapped_strain"),
                            "edge_id": assertion.edge_id,
                            "predicate": assertion.predicate,
                            "object": assertion.obj,
                            "primary_source": source,
                        }
                        for source in assertion.sources
                    )
    return rows


def build_tables(
    graph: Graph, names: dict[str, str]
) -> dict[str, tuple[tuple[str, ...], list[dict[str, Any]]]]:
    """Count distinct carrier nodes and retain original strain/source support."""
    summary_rows, taxon_rows, strain_rows, evidence_rows = [], [], [], []
    taxon_strains: defaultdict[str, set[str]] = defaultdict(set)
    for strain, taxa in graph.strain_taxa.items():
        for taxon in taxa:
            taxon_strains[taxon].add(strain)
    for scope in matrix_api.STRAIN_ROLLUPS:
        scoped = matrix_api.scoped_taxon_traits(
            graph.taxon_traits, graph.strain_taxa, graph.strain_traits, scope
        )
        admitted = {
            strain: {key for key in keys if matrix_api.strain_trait_rolls_up(key, scope)}
            for strain, keys in graph.strain_traits.items()
        }
        for pair in PAIRS:
            left = {node for node, keys in scoped.items() if keys.intersection(pair.left)}
            right = {node for node, keys in scoped.items() if keys.intersection(pair.right)}
            local_taxa = []
            for taxon in sorted(left & right):
                direct = graph.taxon_traits.get(taxon, set())
                direct_left, direct_right = (
                    direct.intersection(pair.left),
                    direct.intersection(pair.right),
                )
                supporters = taxon_strains[taxon]
                left_strains = {
                    strain
                    for strain in supporters
                    if admitted.get(strain, set()).intersection(pair.left)
                }
                right_strains = {
                    strain
                    for strain in supporters
                    if admitted.get(strain, set()).intersection(pair.right)
                }
                same_strains = left_strains & right_strains
                record: dict[str, Any] = {
                    "scope": scope,
                    "pair_id": pair.id,
                    "taxon_id": taxon,
                    "taxon_name": names.get(taxon, ""),
                    "node_metadata_status": "named" if taxon in names else "missing_name",
                    "left_keys": _joined(scoped[taxon].intersection(pair.left)),
                    "right_keys": _joined(scoped[taxon].intersection(pair.right)),
                    "left_direct_keys": _joined(direct_left),
                    "right_direct_keys": _joined(direct_right),
                    "left_strain_ids": _joined(left_strains),
                    "right_strain_ids": _joined(right_strains),
                    "same_strain_ids": _joined(same_strains),
                    "direct_taxon_has_both": bool(direct_left and direct_right),
                    "same_strain_has_both": bool(same_strains),
                    "requires_union_across_subjects": not bool(
                        (direct_left and direct_right) or same_strains
                    ),
                }
                local_taxa.append(record)
                evidence_rows.extend(_evidence(graph, scope, pair, taxon, supporters, admitted))
            taxon_rows.extend(local_taxa)
            coasserted = sorted(
                strain
                for strain, keys in admitted.items()
                if keys.intersection(pair.left) and keys.intersection(pair.right)
            )
            for strain in coasserted:
                mapped = graph.strain_taxa.get(strain, set())
                strain_rows.append(
                    {
                        "scope": scope,
                        "pair_id": pair.id,
                        "strain_id": strain,
                        "left_keys": _joined(admitted[strain].intersection(pair.left)),
                        "right_keys": _joined(admitted[strain].intersection(pair.right)),
                        "mapped_taxon_ids": _joined(mapped),
                        "mapped_taxon_count": len(mapped),
                    }
                )
                if not mapped:
                    evidence_rows.extend(_evidence(graph, scope, pair, "", {strain}, admitted))
            summary_rows.append(
                {
                    "scope": scope,
                    "pair_id": pair.id,
                    "left_keys": _joined(pair.left),
                    "right_keys": _joined(pair.right),
                    "left_taxon_nodes": len(left),
                    "right_taxon_nodes": len(right),
                    "either_taxon_nodes": len(left | right),
                    "coasserted_taxon_nodes": len(left & right),
                    **{
                        f"{field}_nodes": sum(row[field] for row in local_taxa)
                        for field in (
                            "direct_taxon_has_both",
                            "same_strain_has_both",
                            "requires_union_across_subjects",
                        )
                    },
                    "coasserted_strain_nodes": len(coasserted),
                    "unlinked_coasserted_strain_nodes": sum(
                        not graph.strain_taxa.get(strain) for strain in coasserted
                    ),
                }
            )
    return {
        "coassertion_summary.tsv": (SUMMARY_FIELDS, summary_rows),
        "coasserted_taxon_nodes.tsv": (TAXON_FIELDS, taxon_rows),
        "coasserted_strain_nodes.tsv": (STRAIN_FIELDS, strain_rows),
        "coassertion_evidence.tsv": (EVIDENCE_FIELDS, evidence_rows),
    }


def _tsv(fields: tuple[str, ...], rows: list[dict[str, Any]]) -> bytes:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue().encode("utf-8")


def _report(summary: list[dict[str, Any]]) -> bytes:
    lines = [
        "# KG phenotype co-assertion audit",
        "",
        "These are counts of annotations at distinct NCBITaxon carrier nodes, before ancestor",
        "closure. They are not counts of BugSigDB taxa, experiments, patients or matrix cells.",
        "This is not a measurement of physiology or a classification of contradictory organisms.",
        "Broad oxygen labels can be compatible; strains and source conditions can differ.",
        "",
        "| Pair | Scope | Co-asserted taxon nodes | Both direct | "
        "Both on one strain | Requires union |",
        "| --- | --- | ---: | ---: | ---: | ---: |",
    ]
    lines.extend(
        f"| {row['pair_id']} | {row['scope']} | {row['coasserted_taxon_nodes']} | "
        f"{row['direct_taxon_has_both_nodes']} | {row['same_strain_has_both_nodes']} | "
        f"{row['requires_union_across_subjects_nodes']} |"
        for row in summary
    )
    lines.extend(
        [
            "",
            "Direct-both and same-strain-both flags can overlap; do not sum those columns.",
            "Requires-union means no original taxon or strain subject supplies both sides.",
            "Unlinked strains remain in the strain table and do not add taxon nodes.",
            "The full scope is an annotation audit only; no full-scope pipeline was run.",
            "",
            "No oxygen-guild result or sensitivity effect was computed here. See",
            "methods/coassertion_audit.md for the proposed exclusion sensitivity "
            "and owner decision.",
            "Final run-specific counts require this audit against the pinned KG.",
            "",
        ]
    )
    return "\n".join(lines).encode("utf-8")


def _identities(edges: Path, nodes: Path) -> dict[str, Any]:
    implementation = {
        path.relative_to(REPO_ROOT).as_posix(): file_sha256(path)
        for path in sorted((REPO_ROOT / "src" / "mtp").rglob("*.py"))
    }
    implementation["methods/coassertion_audit.md"] = file_sha256(
        REPO_ROOT / "methods/coassertion_audit.md"
    )
    implementation["uv.lock"] = file_sha256(REPO_ROOT / "uv.lock")
    return {
        "inputs": {
            label: {
                "filename": path.name,
                "sha256": file_sha256(path),
                "bytes": path.stat().st_size,
            }
            for label, path in (("kg_edges", edges), ("kg_nodes", nodes))
        },
        "implementation": implementation,
    }


def run_audit(edges: Path, nodes: Path, output: Path) -> Path:
    """Write a new independent audit with exact input, implementation and output hashes."""
    if output.is_symlink() or (output.exists() and (not output.is_dir() or any(output.iterdir()))):
        raise ValueError("Audit output must be a new or empty directory without a symlink alias")
    if any(path.resolve().is_relative_to(output.resolve()) for path in (edges, nodes)):
        raise ValueError("Audit output must not contain its KG inputs")
    identities = _identities(edges, nodes)
    graph = read_graph(edges)
    wanted = (
        set(graph.taxon_traits)
        | {node for nodes_ in graph.strain_taxa.values() for node in nodes_}
        | set(TARGET_KEYS)
    )
    names = read_names(nodes, wanted)
    tables = build_tables(graph, names)
    documents = {name: _tsv(fields, rows) for name, (fields, rows) in tables.items()}
    documents["README.md"] = _report(tables["coassertion_summary.tsv"][1])
    if _identities(edges, nodes) != identities:
        raise ValueError("Audit inputs or implementation changed during execution")
    commit, dirty = repository_state(REPO_ROOT)
    output.mkdir(parents=True, exist_ok=True)
    for name, content in documents.items():
        (output / name).write_bytes(content)
    record = {
        "analysis_kind": "kg_annotation_census_not_a_core_or_arm_rerun",
        "protocol": "scope-coassertion-census-v1",
        "recorded_at": datetime.now(UTC).isoformat(),
        **identities,
        "code_commit": commit,
        "code_dirty_at_execution": dirty,
        "python": platform.python_version(),
        "scope_policy": "mtp.taxa_trait_matrix.scoped_taxon_traits",
        "pairs": [asdict(pair) for pair in PAIRS],
        "class_labels": {key: names[key] for key in sorted(TARGET_KEYS)},
        "kg_edge_serialization": graph.serialization,
        "outputs": {name: file_sha256(output / name) for name in sorted(documents)},
        "sensitivity_status": "proposed_not_executed; owner decision and verified rerun required",
    }
    manifest = output / "provenance.json"
    manifest.write_bytes(
        (
            json.dumps(record, sort_keys=True, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
        ).encode("utf-8")
    )
    return manifest


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--edges", type=Path, required=True)
    parser.add_argument("--nodes", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    arguments = parser.parse_args(argv)
    print(run_audit(arguments.edges, arguments.nodes, arguments.output_dir))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
