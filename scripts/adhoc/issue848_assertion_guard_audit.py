"""Read-only #848 admission-encoding audit, 2026-10-22.

Inputs: KG edge/node TSV paths. Output: JSON parser/guard counts and input SHA-256s.
No matrix, analysis, arm or scientific output is generated.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

from mtp.paths import file_sha256
from mtp.taxa_trait_matrix import (
    HAS_PHENOTYPE,
    NEGATED_PREDICATES,
    PHENOTYPE_EXCLUDED_TRAITS,
    assertion_trait_key,
    read_assertion_term_metadata,
)
from mtp.workflow.kg_inputs import iter_kg_edges


def audit(edges: Path, nodes: Path) -> dict[str, object]:
    """Stream complete records, validating the same candidate assertions as the matrix."""
    before = {"edges": file_sha256(edges), "nodes": file_sha256(nodes)}
    terms = read_assertion_term_metadata(nodes, edges)
    parser: dict[str, object] = {}
    counts: Counter[str] = Counter({"guard_errors": 0})
    predicates: Counter[str] = Counter()
    changed_keys: Counter[tuple[str, str, str]] = Counter()
    errors: list[dict[str, str]] = []
    for record in iter_kg_edges(edges, parser):
        subject, predicate, obj = (
            record[key].strip() for key in ("subject", "predicate", "object")
        )
        if not subject.startswith(("NCBITaxon:", "kgmicrobe.strain:")) or not (
            predicate.startswith("METPO:") or obj.startswith("METPO:")
        ):
            continue
        counts["candidate_assertions"] += 1
        predicates[predicate] += 1
        counts["nonempty_edge_meta"] += bool(record.get("meta", "").strip())
        counts["metpo_object_without_has_phenotype"] += (
            obj.startswith("METPO:") and predicate != HAS_PHENOTYPE
        )
        changes_key = (
            obj.startswith("METPO:")
            and predicate != HAS_PHENOTYPE
            and predicate not in NEGATED_PREDICATES
        )
        counts["metpo_object_previously_collapsed_under_other_predicate"] += changes_key
        if changes_key:
            kind = "strain" if subject.startswith("kgmicrobe.strain:") else "taxon"
            changed_keys[(kind, predicate, obj)] += 1
        counts["direct_handling_class_assertions"] += (
            subject.startswith("NCBITaxon:") and obj in PHENOTYPE_EXCLUDED_TRAITS
        )
        try:
            assertion_trait_key(record, terms)
        except ValueError as error:
            counts["guard_errors"] += 1
            if len(errors) < 10:
                errors.append(
                    {"edge_id": record["id"], "predicate": predicate, "error": str(error)}
                )
        else:
            counts["guard_accepted"] += 1
    after = {"edges": file_sha256(edges), "nodes": file_sha256(nodes)}
    if after != before:
        raise ValueError("Source files changed while the read-only audit ran")
    return {
        "issue": 848,
        "scope": "source_encoding_only_no_matrix_regeneration",
        "input_sha256": before,
        "parser": parser,
        "guard_counts": dict(sorted(counts.items())),
        "candidate_predicates": dict(sorted(predicates.items())),
        "changed_key_assertion_inventory": [
            {"subject_kind": kind, "predicate": predicate, "object": obj, "assertions": count}
            for (kind, predicate, obj), count in sorted(changed_keys.items())
        ],
        "first_guard_errors": errors,
        "passed": not errors,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--edges", required=True, type=Path)
    parser.add_argument("--nodes", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = audit(args.edges, args.nodes)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
