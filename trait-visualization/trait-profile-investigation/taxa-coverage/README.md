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
| Total profiles | 1870 |
| Profiles with mapped traits | 1845 |
| Zero-trait profiles | 25 |
| Unique taxa | 3221 |
| Trait-bearing taxa | 2363 |
| Taxa with non-MediaDive traits | 2363 |
| Taxa with MediaDive annotations only | 0 |
| Profiles with non-MediaDive traits | 1845 |
| Profiles with MediaDive annotations only | 0 |
| Zero-trait taxa | 858 |
| Taxon signature mentions | 19566 |
| Zero-trait taxa recurring in ≥20 profiles | 3 |

Zero-trait profiles remain in the coverage tables. Only the mapped-profile histogram
and mapped-profile richness boxplot omit them; this audit does not declare eligibility
for downstream statistical tests. No recorded annotation is not biological absence.

| Rank | Unique taxa | Zero-trait taxa | Trait-bearing (%) | Signature mentions | Profiles with rank |
| --- | ---: | ---: | ---: | ---: | ---: |
| Species | 1983 | 696 | 64.90 | 6078 | 893 |
| Genus | 842 | 103 | 87.77 | 8556 | 1454 |
| Family | 202 | 26 | 87.13 | 2392 | 783 |
| Order | 99 | 19 | 80.81 | 964 | 444 |
| Class | 50 | 8 | 84.00 | 683 | 390 |
| Phylum | 45 | 6 | 86.67 | 893 | 515 |

| Distribution | Median | Q1 | Q3 | Mean |
| --- | ---: | ---: | ---: | ---: |
| Traits per mapped profile | 186 | 79 | 530 | 358.723 |
| Mapped taxa per mapped profile | 6 | 3 | 12 | 9.70894 |
| Traits per reported taxon | 10 | 0 | 24 | 25.8802 |
| Traits per trait-bearing taxon | 15 | 8 | 36 | 35.2772 |

Annotation counts can include traits aggregated from descendants; they do not measure
physiology, patient-level abundance or activity. Suggested zero-trait diagnoses are
heuristics for review, not verified causes or permission to propagate annotations downward.

Regenerate from this repository's root with the same input paths and scope:

```sh
just taxa-coverage --strain-rollup phenotype --contrast-population disease_vs_healthy_tier1_2_subject_controls --taxa data/derived/bugsigdb_taxa.tsv --matrix data/derived/taxa_trait_matrix.tsv --profiles data/derived/study_trait_profiles_wide.tsv --experiments data/derived/experiments.tsv --output-dir trait-visualization/trait-profile-investigation/taxa-coverage --kg-provenance data/kg-microbe/provenance.json --raw-dump data/raw/full_dump.csv
```
The five TSVs, this README, and generated figures are inventoried by SHA-256 in `provenance.json`.
`profile_coverage_by_rank.tsv` uses `pct_of_mapped_profiles`, whose denominator is the
observed count of mapped profiles. Rank tables omit empty ranks and include
an other/unrecorded category when present in the selected corpus.
The columns consumed by the slide and plate builders remain available in the rank and taxon detail TSVs.
