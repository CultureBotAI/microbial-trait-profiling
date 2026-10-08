"""Audit Table 1 trait semantics, ontology overlap, and source evidence.

Exports data/derived/table1_trait_semantic_audit.tsv and verifies provenance against
data/kg-microbe/provenance.json.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from bsdbng.trait_semantics import (
    DEFAULT_OUTPUT_TSV,
    CoccusSphereComparison,
    build_semantic_audit_records,
    verify_provenance_checksums,
    write_semantic_audit_tsv,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Audit Table 1 trait semantics and export audit TSV."
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT_TSV,
        help=f"Output TSV file path (default: {DEFAULT_OUTPUT_TSV})",
    )
    parser.add_argument(
        "--strict-provenance",
        action="store_true",
        default=True,
        help="Enforce exact match against pinned KG-Microbe checksums in provenance.json",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    if args.strict_provenance:
        if not verify_provenance_checksums():
            sys.stderr.write(
                "Error: data/kg-microbe/provenance.json does not match pinned KG checksums.\n"
            )
            return 1
        print("KG-Microbe provenance checksums verified against pinned hashes.")

    records = build_semantic_audit_records()
    output_path = write_semantic_audit_tsv(records, args.output)
    print(f"Exported Table 1 semantic audit table with {len(records)} traits to: {output_path}")

    # Display summary
    print("\n" + "=" * 80)
    print("TABLE 1 TRAIT SEMANTIC AUDIT SUMMARY")
    print("=" * 80)
    for r in records:
        print(f"\n[{r.trait_key}] {r.feature_name}")
        print(f"  Ontology: {r.ontology_id} ('{r.ontology_label}')")
        print(f"  Sources:  {r.source_knowledge_source}")
        print(f"  Evidence: {r.evidence_type} | Status: {r.evidence_status}")
        print(
            f"  Support:  {r.direct_kg_taxa_count} direct KG taxa -> "
            f"{r.bugsigdb_matrix_taxa_count} BugSigDB taxa -> "
            f"{r.profile_support_count} profiles ({r.profile_support_pct:.1f}%)"
        )
        print(f"  Decision: {r.decision}")
        print(f"  Approved: {r.approved_wording}")
        print(f"  Prohibited: {r.unsupported_interpretations}")

    comp = CoccusSphereComparison()
    print("\n" + "=" * 80)
    print("COCCUS VS SPHERE COMPARISON SUMMARY")
    print("=" * 80)
    print(f"Ontology: {comp.ontology_relationship}")
    print(
        f"Direct KG Taxa: Coccus={comp.coccus_direct_taxa}, Sphere={comp.sphere_direct_taxa}, "
        f"Intersection={comp.direct_taxa_intersection} (Jaccard={comp.direct_taxa_jaccard:.4f})"
    )
    print(
        f"BugSigDB Matrix Taxa: Coccus={comp.coccus_bugsigdb_taxa}, "
        f"Sphere={comp.sphere_bugsigdb_taxa}, "
        f"Intersection={comp.bugsigdb_taxa_intersection} "
        f"(Jaccard={comp.bugsigdb_taxa_jaccard:.4f})"
    )
    print(
        f"Profile Support: Coccus={comp.coccus_profile_nonzero_count}, "
        f"Sphere={comp.sphere_profile_nonzero_count}"
    )
    print(
        f"Correlations: Pearson r={comp.pearson_correlation:.4f}, "
        f"Spearman rho={comp.spearman_correlation:.4f}"
    )
    print(f"Decision: {comp.decision}")
    print(f"Rationale: {comp.rationale}")
    print("=" * 80)

    return 0


if __name__ == "__main__":
    sys.exit(main())
