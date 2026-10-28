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
`--madin-csv MADIN/condensed_traits_NCBI.csv`, the pinned Madin et al. table
described below; a file with any other bytes fails before outputs are written.
Omitting it leaves Madin subsources unknown. A supplied but unreadable or
malformed table fails instead of silently falling back to unknown.

The `phenotype` scope uses the public `scoped_taxon_traits` and
`strain_trait_rolls_up` functions introduced by PR #828. Both `off` and
`phenotype` use the merged matrix scope policy. `full` is not offered by this
sidecar.

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

## Pinned Madin et al. table

KG-Microbe attributes every trait it takes from the Madin et al. (2020)
compilation to one source, `infores:madin_etal`. To resolve that source to the
compiled datasets, the audit reads the compilation's record-level table, pinned
in `mtp.madin_source` (#1048):

| Field | Value |
| --- | --- |
| Repository | `bacteria-archaea-traits/bacteria-archaea-traits` on GitHub |
| Commit | `396dcd4d98ca20b613ba280f5514611651b6fb91` (2020-10-16) |
| File | `output/condensed_traits_NCBI.csv` |
| URL | `https://raw.githubusercontent.com/bacteria-archaea-traits/bacteria-archaea-traits/396dcd4d98ca20b613ba280f5514611651b6fb91/output/condensed_traits_NCBI.csv` |
| SHA256 | `6100fc11679cd9facae65cad37a2969978bc64c3e57984d896f8147f78db6440` |
| Size | 46,362,149 bytes |

```bash
just download-madin ../mtp-data/madin
```

This fetches the table anonymously into a durable data directory outside the
working tree. Choose an absolute path or an outside-checkout relative path,
such as the sibling directory in this example. Temporary destinations and
paths inside the working tree are refused before a request or write; the
checks include resolved symlink destinations. The downloader records the
operation in `download_provenance.json` beside the table.

A response whose size or SHA256 differs from the pin is discarded and the
recipe fails. The stream stops before writing bytes beyond the pinned size;
local copies are checked for regular-file type and size before hashing.
An existing verified copy is reused without a request or overwrite, and
reuse also writes a verification record. A corrupt existing provenance
record is refused before any download or placement. Provenance is written
atomically, preserving unrelated entries. Reuse records verification time
without inventing a download time for a copy supplied by hand.

Passing the file as `--madin-csv` records
`madin_status: pinned_verified` and the pin in `provenance.json`. Without
`--madin-csv`, the status stays `not_supplied_source_support_unknown`, the pin
is still recorded, and Madin carriers keep `unknown` source status.

The pinned file has one row per source record, with the NCBI `tax_id` and
`data_source` columns the audit needs. The species-level
`condensed_species_NCBI.csv` merges records and has neither column, so it
cannot resolve subsources. KG-Microbe's own Madin download also fetches
`condensed_traits_NCBI.csv`, from the default branch.

How the pin was confirmed. #773 calls the table it traced the "2021-02-24
release"; no such release exists:

- The upstream repository has one release and one tag, v1.0.0 (2019-09-26).
  "2021-02-24" is the file modification time of the local copy the #773
  tracing read, that is, when it was downloaded.
- That copy's git blob hash, `05a10e84ea015f9ee2b99ea16aae6d26c557ce3b`, equals
  the blob of `output/condensed_traits_NCBI.csv` at commit `396dcd4`
  (2020-10-16). The repository's commit history shows `396dcd4` as the last
  commit to change the file, and the default branch (`master`) has pointed at
  it since.
- The file fetched from the pinned URL has the SHA256 above. So does a 2026
  local copy of KG-Microbe's raw Madin download. The `download.yaml` of
  `Knowledge-Graph-Hub/kg-microbe` writes that download, from the default
  branch, to `data/raw/madin_etal.csv`, which the repository gitignores; it is
  not a file in the repository.

Two earlier tables have other bytes, so neither can substitute for the pin:

| Table | Size | SHA256 |
| --- | --- | --- |
| v1.0.0 (tag commit `86dbbbcea2cfd7a3be7003fc556a3ec21b969186`) | 41,870,010 bytes | `92edd6e59687cdb832dfc926bf1e9263978b43120ccea664e42347890e0b6b20` |
| `jmadin/bacteria-archaea-traits` fork head (`e9304397f430b2966834b979769da5491da22a69`, 2020-02-18) | 45,677,024 bytes | `196421061bb2e15f43a39799e9cbf2701cec36bac8b237718aa91832d1baeb24` |

Both predate the pinned commit; neither is a newer version to upgrade to.

Limits of the pin:

- Madin labels do not establish experimental validation. `data_source` names
  the compiled dataset a record came from. Those datasets range from
  genome-project metadata to curated literature assignments and predictions.
  A verified pin fixes which bytes were read, not whether, how or in which
  strain a trait was measured.
- The KG export does not record which Madin bytes it ingested. Because the
  default branch has not changed the file since the pinned commit, a KG build
  that fetched it afterwards read these bytes. That is inferred from the
  repository history, not checked against the KG.
- Records reach carriers by `tax_id` alone. A Madin carrier with no row
  matching the trait's declared field and value stays `unknown`.
- Supplying the table changes only the Madin-derived fields of these sidecars:
  - `carrier_source_evidence.tsv`: `madin_subsources` and `madin_resolution`;
  - the two support tables: `madin_subsources`, `unresolved_madin_carriers` and
    `source_status`;
  - `support_threshold_sensitivity.tsv`: the `exclude_known_gold_patric_only`
    rows and the `source_unknown_taxa` column;
  - `provenance.json`: `madin_status`, `source_status_counts` and the Madin
    input record.

  The KG, matrix, profiles and arm calls are unchanged.

## Products and interpretation

| Product | Meaning |
| --- | --- |
| `taxon_trait_label_support.tsv` | Positive audited matrix cells: distinct carrier taxa, support at the reported node versus descendants, fraction of annotated taxa in the subtree, source status, and bound occurrence counts. |
| `focal_taxa_trait_support.tsv` | The same audit for explicitly listed focal taxa, including absence from the supplied matrix. |
| `carrier_source_evidence.tsv` | Named KG source fields and optional Madin subsources for each carrier supporting an audited matrix cell. |
| `support_threshold_sensitivity.tsv` | Six focal traits under seven declared support rules, with effect units, paired/study counts, ties, and nominal study sign tests. |
| `coasserted_label_census.tsv` | Co-assertion by species, genus, family and other/unknown rank, among bound reported taxa and their occurrences; missing features remain unavailable. |
| `provenance.json` | SHA256 and size of all inputs, all five output tables, producer, matrix policy and Madin verification code, environment lock, scope reconciliation, source rules, the Madin pin and status, and limitations. |

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
All seven rules use the same paired population: both sides must have positive
total retained non-MediaDive annotation counts. The input paired count and
excluded zero-total paired count are reported for every rule. This prevents a
change in the eligible population from being mistaken for a support-rule effect.
All eligible experiment differences, including zeros, are averaged within a
study; studies then receive equal weight. The two-sided sign test uses only
non-tied study means, with a `1e-12` absolute tie tolerance. No non-tied studies
means an unavailable test, not p = 1. These selected traits and rules are
exploratory: p-values are nominal, no adjusted calls are made, and the primary
arm families and calls remain unchanged.

| Rule | Profile score and denominator |
| --- | --- |
| `native_annotation_fraction` | Focal count divided by all retained non-MediaDive annotation counts, on the common positive-total paired population. |
| `carriers_per_reported_taxon` | Sum of focal carrier counts divided by unique reported taxa on that side. This can exceed one and is not the native fraction. |
| `binary_presence` | Fraction of reported taxa with at least one carrier. |
| `min_two_carriers` | Fraction of reported taxa with at least two carriers. |
| `min_subtree_fraction_0.10` | Fraction of reported taxa whose annotated-subtree support is at least 0.10. |
| `subtree_fraction` | Mean annotated-subtree support fraction across reported taxa. |
| `exclude_known_gold_patric_only` | Carrier count per reported taxon after removing only confidently resolved GOLD/PATRIC-only cells. Unknown support is retained and counted. |

The support table audits all supplied matrix taxa; its occurrence counts identify
which are used by this population. Sensitivity source-uncertainty counts include
only taxa in the common positive-total paired population, even when the matrix covers a larger
population. Unmapped reported taxa remain in the reported-taxon denominator with zero
mapped support; they are counted explicitly in provenance. Missing profile
sides are never replaced by zero. These different denominators answer different
questions, so changes in effect size across rules cannot alone establish that
one rule is biologically correct.

Hashes and exact count reconciliation bind this sidecar's inputs and outputs.
They do not constitute a verified core execution record or independently prove
the original corpus selection. No primary result or talk-ready number is
regenerated merely by adding this implementation.

The KG edge reader is the workflow's versioned `iter_kg_edges` contract. It
reconstructs the known bare-carriage-return header and agent-field encodings;
unknown malformed records fail. The sidecar records repair counts and hashes
the parser and shared trait-definition code. Pipe-joined primary sources are
split before Madin lookup, preserving each source and its resolution status.
Table 1 trait membership and the oxygen guild pair use `trait_guilds` definitions.

The broader node-level co-assertion census, mesophilic/psychrophilic and
facultative-anaerobe/obligate-aerobe pairs, and the off/phenotype/full scope
comparison belong to #846. This implementation continues to support off and
phenotype only. Feature-catalog and canonical semantics-audit integration remain
open under #773; these sidecars do not complete those items.
