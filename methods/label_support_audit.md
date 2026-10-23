# Label-support sidecar

`just audit-label-support` audits the provenance and taxonomic support behind
selected trait annotations (#773). It writes independent tables and never
changes the feature catalog, trait semantics table, matrix, profiles, or arm
calls. These are annotation-support sensitivities, not measurements of microbial
physiology or independent replication of the original findings.

## Run against matching inputs

Use the reported-taxon, taxon-matrix, and wide-profile tables from one run, with
the exact KG snapshot used to produce them. Supply the run's strain scope:

```bash
just audit-label-support \
  --taxa RUN/bugsigdb_taxa.tsv \
  --matrix RUN/taxa_trait_matrix.tsv \
  --profiles RUN/study_trait_profiles_wide.tsv \
  --edges KG/merged-kg_edges.tsv \
  --nodes KG/merged-kg_nodes.tsv \
  --strain-rollup off \
  --output-dir RESULTS/label-support-off
```

Replace these example paths with the actual products. The destination must be
empty or new. No other pipeline step is triggered. Optionally supply
`--madin-csv SOURCES/madin.csv`; it is hashed but not downloaded or pinned.
Omitting it leaves Madin subsources unknown. A supplied but unreadable or
malformed table fails instead of silently falling back to unknown.

The `phenotype` scope uses the public `scoped_taxon_traits` and
`strain_trait_rolls_up` functions introduced by PR #828. Until that dependency
is present, `--strain-rollup phenotype` fails explicitly; no private duplicate
of the strain admission rules is used. `off` works independently. `full` is not
offered by this sidecar.

The profile keys select the audit population. Extra occurrence rows are
excluded; each supplied profile side must have reported taxa, and summing the
corresponding unique taxon-matrix rows must exactly reproduce every supplied
profile count. Every retained feature with positive counts among those taxa
must be in the profiles; columns empty in that population may be dropped. This prevents
an accidentally truncated denominator. All audited matrix cells are also
reconciled against distinct KG carrier taxa. Scope, snapshot, and orientation
mismatches fail before outputs are written. Reversed profiles require matching
oriented occurrence directions; the audit does not independently classify
contrasts or certify a host/population selection. Counts alone cannot resolve
an orientation mismatch when both sides happen to have identical counts.

## Products and interpretation

| Product | Meaning |
| --- | --- |
| `taxon_trait_label_support.tsv` | Positive audited matrix cells: distinct carrier taxa, support at the reported node versus descendants, fraction of annotated taxa in the subtree, source status, and bound occurrence counts. |
| `focal_taxa_trait_support.tsv` | The same audit for explicitly listed focal taxa, including absence from the supplied matrix. |
| `carrier_source_evidence.tsv` | Named KG source fields and optional Madin subsources for each carrier supporting an audited matrix cell. |
| `support_threshold_sensitivity.tsv` | Six focal traits under seven declared support rules, with effect units, paired/study counts, ties, and nominal study sign tests. |
| `coasserted_label_census.tsv` | Co-assertion by species, genus, family and other/unknown rank, among bound reported taxa and their occurrences; missing features remain unavailable. |
| `provenance.json` | SHA256 and size of all inputs, all five output tables, producer and matrix policy code, environment lock, scope reconciliation, source rules, and limitations. |

A carrier is a distinct NCBITaxon node at or below the reported taxon, including
the reported taxon itself when it carries an annotation under the selected scope.
`carrier_at_reported_taxon` describes that node's position, not the original KG
assertion subject: admitted strain assertions can also create a carrier there.
Multiple source records on
one carrier remain one matrix count. Therefore `single_carrier_flag` does not
assert that there was only one experimental record. The subtree denominator is
all taxa bearing any non-MediaDive KG annotation admitted by the selected scope
in that subtree; it is not all extant
species, all KG descendants, or the number of sampled organisms. Phenotype
strain assertions map to their NCBITaxon carrier and follow #828's deduplication.

`gold_patric_only` requires all primary support to be Madin, at least one
resolved subsource for every Madin carrier, and all resolved subsources to be
GOLD/PATRIC. Missing provenance is `unknown`, not evidence against metadata
support. `includes_other_sources` means another declared source exists; it
does not certify experimental validation. Madin mappings are explicit
field/value rules recorded in provenance. Taxon IDs and exact normalized
tokens are used, not name heuristics. Unsupported biological classifications
of carriers are not inferred.

Co-asserted oxygen or morphology labels can refer to different descendants,
compatible terms, or inconsistent source assertions. The census is not a claim
that individual cells express mutually exclusive physiological states.

## Sensitivity estimands

For each paired experiment, each rule computes increased minus decreased.
All eligible experiment differences, including zeros, are averaged within a
study; studies then receive equal weight. The two-sided sign test uses only
non-tied study means, with a `1e-12` absolute tie tolerance. No non-tied studies
means an unavailable test, not p = 1. These selected traits and rules are
exploratory: p-values are nominal, no adjusted calls are made, and the primary
arm families and calls remain unchanged.

| Rule | Profile score and denominator |
| --- | --- |
| `native_annotation_fraction` | Focal count divided by all retained non-MediaDive annotation counts. A pair with a zero annotation total is excluded and counted. |
| `carriers_per_reported_taxon` | Sum of focal carrier counts divided by unique reported taxa on that side. This can exceed one and is not the native fraction. |
| `binary_presence` | Fraction of reported taxa with at least one carrier. |
| `min_two_carriers` | Fraction of reported taxa with at least two carriers. |
| `min_subtree_fraction_0.10` | Fraction of reported taxa whose annotated-subtree support is at least 0.10. |
| `subtree_fraction` | Mean annotated-subtree support fraction across reported taxa. |
| `exclude_known_gold_patric_only` | Carrier count per reported taxon after removing only confidently resolved GOLD/PATRIC-only cells. Unknown support is retained and counted. |

The support table audits all supplied matrix taxa; its occurrence counts identify
which are used by this population. Sensitivity source-uncertainty counts include
only taxa bound to the supplied profiles, even when the matrix covers a larger
population. Unmapped reported taxa remain in the reported-taxon denominator with zero
mapped support; they are counted explicitly in provenance. Missing profile
sides are never replaced by zero. These different denominators answer different
questions, so changes in effect size across rules cannot alone establish that
one rule is biologically correct.

Hashes and exact count reconciliation bind this sidecar's inputs and outputs.
They do not constitute a verified core execution record or independently prove
the original corpus selection. No primary result or talk-ready number is
regenerated merely by adding this implementation.
