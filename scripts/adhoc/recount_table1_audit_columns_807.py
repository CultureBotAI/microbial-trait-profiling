"""Recount the Table 1 audit's matrix and profile columns for #807 (2026-10-23).

Purpose: recompute the numbers in docs/table1_trait_semantics_audit.md that come from a
taxa-trait matrix and a profile table. For each of the six Table 1 traits, the script
gives the matrix taxa that bear the trait, as a count and as a percentage of all matrix
taxa, and the profiles with a nonzero count. It also gives the coccus/sphere overlap in
the matrix, and the Pearson and Spearman correlations of the coccus and sphere counts
across every profile (#807, #1004).

Inputs:
- --matrix: a taxa_trait_matrix.tsv. Expected sha256 prefix 45df0fe7 for the final
  rerun's phenotype-scope matrix and f6de0576 for its off-scope matrix. Neither is
  committed; data/analysis/cellulose_decr/cellulose_decr_provenance.json pins 45df0fe7.
- --profiles: a study_trait_profiles_wide table (.tsv or .tsv.gz). The default is the
  committed data/derived/study_trait_profiles_wide.tsv.gz, which decompresses to sha256
  43b2e230 (phenotype). The off table, 519a095b, is that file as committed at 511c891.
- --benchmark: the committed data/analysis/cellulose_decr/trait_convergence_benchmark.tsv
  (sha256 ad18059e). Its matrix_taxon_count cross-checks the phenotype matrix.

Output: plain-text counts on stdout.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import pandas as pd
from scipy import stats

SIX = {
    "METPO:1000608": "facultatively aerobic",
    "METPO:2000103|GO:0030245": "cellulose catabolic process",
    "METPO:1000668": "coccus shaped",
    "METPO:1000683": "sphere shaped",
    "METPO:2000103|pathways:arsenate_reduction": "arsenate reduction",
    "METPO:1000623": "moderately halophilic",
}
COCCUS, SPHERE = "METPO:1000668", "METPO:1000683"


def recount(matrix_path: Path, profiles_path: Path, benchmark_path: Path) -> None:
    matrix = pd.read_csv(matrix_path, sep="\t", usecols=["taxon_id", *SIX], index_col="taxon_id")
    n_taxa = len(matrix)
    with benchmark_path.open(encoding="utf-8", newline="") as handle:
        benchmark = {r["trait_id"]: r for r in csv.DictReader(handle, delimiter="\t")}
    print(f"matrix taxa (rows): {n_taxa}")
    for trait, label in SIX.items():
        k = int((matrix[trait] > 0).sum())
        committed = benchmark.get(trait, {}).get("matrix_taxon_count", "not in benchmark")
        print(
            f"  {label:28s} K={k:5d} ({100 * k / n_taxa:.1f}%)  "
            f"benchmark matrix_taxon_count={committed}"
        )
    coccus = set(matrix.index[matrix[COCCUS] > 0])
    sphere = set(matrix.index[matrix[SPHERE] > 0])
    shared, union = len(coccus & sphere), len(coccus | sphere)
    print(
        f"coccus K={len(coccus)} sphere K={len(sphere)} intersection={shared} union={union} "
        f"sphere-in-coccus={100 * shared / len(sphere):.1f}% "
        f"coccus-in-sphere={100 * shared / len(coccus):.1f}% Jaccard={shared / union:.4f}"
    )

    profiles = pd.read_csv(profiles_path, sep="\t", usecols=["study_id", *SIX])
    n_profiles = len(profiles)
    print(f"profiles: {n_profiles}")
    for trait, label in SIX.items():
        nonzero = int((profiles[trait] > 0).sum())
        print(f"  {label:28s} nonzero={nonzero:5d} ({100 * nonzero / n_profiles:.1f}%)")
    pearson = stats.pearsonr(profiles[COCCUS], profiles[SPHERE])
    spearman = stats.spearmanr(profiles[COCCUS], profiles[SPHERE])
    print(
        f"coccus vs sphere profile counts: Pearson r={pearson[0]:.3f} (p={pearson[1]:.3g}); "
        f"Spearman rho={spearman[0]:.3f} (p={spearman[1]:.3g})"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument(
        "--profiles", type=Path, default=Path("data/derived/study_trait_profiles_wide.tsv.gz")
    )
    parser.add_argument(
        "--benchmark",
        type=Path,
        default=Path("data/analysis/cellulose_decr/trait_convergence_benchmark.tsv"),
    )
    args = parser.parse_args()
    recount(args.matrix, args.profiles, args.benchmark)


if __name__ == "__main__":
    main()
