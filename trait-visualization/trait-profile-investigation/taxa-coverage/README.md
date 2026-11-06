# Taxa coverage and profile annotation mapping

Population: Homo sapiens; sequencing types 16S, WMS.
No contrast-class, control-tier, pairing or direction-orientation restriction is applied.
Declared upstream strain-rollup scope: `phenotype`.
This required caller declaration is recorded as `declared_unverified`.
Input hashes and count reconciliation bind the consumed tables but
do not verify the upstream rollup scope; no scope is inferred from feature counts.
Profiles are study / experiment / curated direction; taxa and signature mentions
are restricted to those same profiles. The matrix is restricted to the profile vocabulary,
including MediaDive columns if present. Source hashes and exclusions are in `provenance.json`.
Core counts exclude every MediaDive-qualified column, case-insensitively. The
all-annotation counts and coverage figures retain those columns; the rows below
show how many taxa and profiles have only MediaDive annotations.

| Metric | Count |
| --- | ---: |
| Total profiles | 5629 |
| Profiles with mapped traits | 5538 |
| Zero-trait profiles | 91 |
| Unique taxa | 4750 |
| Trait-bearing taxa | 3274 |
| Taxa with non-MediaDive traits | 3272 |
| Taxa with MediaDive annotations only | 2 |
| Profiles with non-MediaDive traits | 5538 |
| Profiles with MediaDive annotations only | 0 |
| Zero-trait taxa | 1476 |
| Taxon signature mentions | 46865 |
| Zero-trait taxa recurring in ≥20 profiles | 12 |

Zero-trait profiles remain in the coverage tables. Only the mapped-profile histogram
and mapped-profile richness boxplot omit them; this audit does not declare eligibility
for downstream statistical tests. No recorded annotation is not biological absence.

| Rank | Unique taxa | Zero-trait taxa | Trait-bearing (%) | Signature mentions | Profiles with rank |
| --- | ---: | ---: | ---: | ---: | ---: |
| Species | 3165 | 1246 | 60.63 | 15507 | 2546 |
| Genus | 1086 | 147 | 86.46 | 20021 | 4158 |
| Family | 253 | 39 | 84.58 | 5840 | 2024 |
| Order | 132 | 25 | 81.06 | 2117 | 1062 |
| Class | 64 | 10 | 84.38 | 1471 | 872 |
| Phylum | 49 | 9 | 81.63 | 1908 | 1159 |
| Superkingdom | 1 | 0 | 100.00 | 1 | 1 |

| Distribution | Median | Q1 | Q3 | Mean |
| --- | ---: | ---: | ---: | ---: |
| Traits per mapped profile | 156 | 60 | 387 | 290.912 |
| Mapped taxa per mapped profile | 4 | 2 | 9 | 7.76887 |
| Traits per reported taxon | 8 | 0 | 20 | 21.6143 |
| Traits per trait-bearing taxon | 14 | 7 | 32 | 31.3586 |

Annotation counts can include traits aggregated from descendants; they do not measure
physiology, patient-level abundance or activity. Suggested zero-trait diagnoses are
heuristics for review, not verified causes or permission to propagate annotations downward.

Regenerate from this repository's root with the same input paths and scope:

```sh
just taxa-coverage --strain-rollup phenotype --taxa data/derived/bugsigdb_taxa.tsv --matrix data/derived/taxa_trait_matrix.tsv --profiles data/derived/study_trait_profiles_wide.tsv --experiments data/derived/experiments.tsv --output-dir trait-visualization/trait-profile-investigation/taxa-coverage --kg-provenance data/kg-microbe/provenance.json --raw-dump data/raw/full_dump.csv
```
The five TSVs, this README, and generated figures are inventoried by SHA-256 in `provenance.json`.
`profile_coverage_by_rank.tsv` uses `pct_of_mapped_profiles`, whose denominator is the
observed count of mapped profiles. Rank tables omit empty ranks and include
an other/unrecorded category when present in the selected corpus.
The columns consumed by the slide and plate builders remain available in the rank and taxon detail TSVs.
