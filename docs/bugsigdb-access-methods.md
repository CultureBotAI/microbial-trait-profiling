# BugSigDB Data Access Methods

This note documents the available ways to access BugSigDB data, what we
learned about each, and the criteria for choosing among them in `microbial-trait-profiling`.

## BugSigDB is a structured database

BugSigDB is a MediaWiki-based structured database. The data is entered
through forms with distinct Study, Experiment, and Signature entities.
The CSV exports are a flattened representation of this structure — the
hierarchy is encoded in the `BSDB ID` column (e.g. `bsdb:10202341/1/1`
= study/experiment/signature) and repeated metadata columns.

Our ingestion pipeline is essentially unflattening the CSV back into the
structure that BugSigDB already stores internally. This means we are
working around CSV encoding artifacts (comma-separated ontology IDs,
pipe/semicolon-encoded MetaPhlAn lineage paths) that wouldn't exist if
we had access to the structured data directly.

## Study count comparison (2026-04-17)

| Source | Studies | Notes |
|--------|---------|-------|
| BugSigDB main page (live) | 2,010 | Authoritative count |
| studies.csv (Help:Export) | 2,009 | 1,993 Complete, 10 Incomplete, 6 blank State |
| full_dump.csv (BugSigDBExports) | 1,949 | 60 studies lost to upstream filtering |
| This project's output | 1,923 | After our filtering (NA direction/taxa) |

The 60-study gap between studies.csv and full_dump.csv is caused by
BugSigDBExports' `dump_release.R` filtering. We don't control that filter.

## Available methods

### 1. Raw CSV exports from BugSigDB

BugSigDB exposes three raw export tables on `Help:Export`:

- `studies.csv` — 2,009 studies (13 columns including State)
- `experiments.csv` — 8,091 rows (experiment-level metadata)
- `signatures.csv` — 13,023 rows (taxa, direction, curator info)

These are the closest public exports to the current operational state of the
database. They are the least transformed input surface available to us.

**Important:** these are separate tables that must be joined to produce
the full study→experiment→signature→taxa hierarchy. The join keys are
encoded in the BSDB ID column and the Study/Experiment page name columns.

### 2. BugSigDBExports repository

The `waldronlab/BugSigDBExports` repository is a derived export layer built
from the three raw BugSigDB CSV tables.

Its `dump_release.R` script:

1. obtains and merges the study, experiment, and signature tables from
   `Help:Export`
2. **filters incomplete records** (this is where 60 studies are lost)
3. adds signature IDs
4. writes `full_dump.csv` (51 columns, pre-merged flat table)
5. writes GMT files for multiple identifier schemes and taxonomic levels

`full_dump.csv` is convenient (pre-joined, one file) but inherits upstream
filtering decisions we don't control.

### 3. Stable snapshot releases

BugSigDB publishes stable snapshot releases on Zenodo. Useful when
reproducibility matters more than freshness.

### 4. `bugsigdbr` R package

The Bioconductor R package for programmatic access to BugSigDB. May
provide structured access to the data without CSV flattening artifacts.
Tracked in issue #62.

## CSV encoding artifacts we work around

These are properties of the CSV flattening, not of the underlying data:

| Artifact | Where | What we do |
|----------|-------|------------|
| Pipe-separated lineage paths | NCBI Taxonomy IDs column | Split on `;` for taxa, take leaf of `\|`-separated lineage |
| MetaPhlAn rank prefixes | MetaPhlAn taxon names column | Parse `k__`, `g__`, `s__` etc. for rank |
| Comma-separated ontology IDs | UBERON ID, EFO ID columns | Split into multivalued lists |
| Mixed ontology prefixes in "EFO ID" | EFO ID column | 20+ prefixes (EFO, MONDO, HP, CHEBI, GO...), renamed to condition_ontology_id |
| Repeated study/experiment metadata | Every row | Verified: fields are consistent within each study/experiment group — first row is representative |

## Column coverage (51 columns in full_dump.csv)

- **45 captured** in `schema/bsdbng.yaml` (study, experiment, signature, taxon levels)
- **6 not captured** (curatorial metadata): Signature page name, Curated date, Curator, Revision editor, State, Reviewer
- **State column** should be used to filter: only ingest Complete records
  (not yet implemented — tracked in issue #64)

## Current approach in microbial-trait-profiling

`microbial-trait-profiling` currently uses `full_dump.csv` (from BugSigDBExports) as the
primary input, with the three raw CSVs downloaded but unused. The merge
from separate CSVs is a `NotImplementedError` stub.

**Planned change (issue #64):** switch to using the three separate CSVs
as primary input, implementing the join in Python, so we control filtering
and don't lose 60 studies to upstream decisions.

## Selection criteria

- source fidelity
- reproducibility
- freshness
- developer time
- machine-time performance
- good citizenship toward upstream infrastructure
- fit for purpose

## Sources

- BugSigDB `Help:Export`: https://bugsigdb.org/Help:Export
- BugSigDBExports README:
  https://github.com/waldronlab/BugSigDBExports/blob/devel/README.md
- BugSigDB paper:
  https://doi.org/10.1038/s41587-023-01872-y
  https://pmc.ncbi.nlm.nih.gov/articles/PMC11098749/
- Stable releases:
  https://zenodo.org/records/15272273
  https://zenodo.org/records/13997429
  https://zenodo.org/records/10407666
