# BugSigDB Data Access Methods

This note documents the available ways to access BugSigDB data, what we
learned about each, and the criteria for choosing among them in `microbial-trait-profiling`
(including the `clustered-heatmap-generation` workflow).

## Current Decision

Use `waldronlab/BugSigDBExports` tag `v1.3.1` as the citable BugSigDB input for
regenerated figures and derived analysis artifacts.

- version: `v1.3.1`
- DOI: `10.5281/zenodo.19740674`
- tag SHA: `5815c2638d6b56134ac8db19edc2bb3c4b17f1ca`
- branch input file: `data/raw/full_dump.csv`
- branch corpus after ingest: 1,245 usable studies
- human-host figure corpus: 5,829 profiles from 1,080 studies

The live raw exports are downloaded for comparison and future ingestion work,
but the branch does not use them as the primary input until separate-table
joining is implemented in Python.

## Count Comparison

Counts were refreshed on 2026-08-26.

| Source | Studies | Rows | Notes |
| --- | ---: | ---: | --- |
| raw `studies.csv` | 2,191 | 2,191 | 2,169 complete, 22 incomplete |
| raw `experiments.csv` | n/a | 9,140 | 7,210 Homo sapiens experiment rows |
| raw `signatures.csv` | n/a | 14,474 | 14,468 complete, 4 incomplete, 2 blank state |
| `BugSigDBExports` `v1.3.1` `full_dump.csv` | 1,246 parsed IDs | 7,425 | complete rows only |
| bsdbng YAML output | 1,245 | n/a | one parsed study skipped after ingestion filtering |

## Available Methods

### 1. Raw CSV Exports From BugSigDB

BugSigDB exposes three raw export tables on `Help:Export`:

- `studies.csv`
- `experiments.csv`
- `signatures.csv`

These are the closest public exports to the live operational state of the
database. They are separate tables, so they must be joined to recover the
Study -> Experiment -> Signature -> Taxa hierarchy. The join keys are encoded
in `BSDB ID` plus the study and experiment page-name columns.

bsdbng downloads these files from the current S3-backed `Help:Export` links,
but `ingest.py` still raises `NotImplementedError` for separate raw-table
merging when `full_dump.csv` is absent.

### 2. BugSigDBExports Tagged Full Dump

`waldronlab/BugSigDBExports` is a derived export layer built from the BugSigDB
tables. Its `full_dump.csv` is pre-joined and convenient for bsdbng's current
ingestion code.

For reproducible branch artifacts, use the `v1.3.1` tag rather than `devel`:

- `v1.3.1` is reviewed/complete enough for citation and has a DOI.
- `devel` is fresher but moving and uncitable.
- `v1.3.0` should not be used for new work; it is older and was only named in
  earlier scripts/issues because it was current when that code was written.

### 3. Zenodo Release

Zenodo is the authoritative citation target for the pinned release. The branch
stores the DOI and tag SHA in `data/raw/download_provenance.json` when
`download_exports(force=True)` refreshes `full_dump.csv`.

### 4. `bugsigdbr` R Package

The Bioconductor `bugsigdbr` package may eventually provide a better structured
access path. That is separate from the current Python regeneration workflow.

## Regeneration Workflow

The intended branch refresh command is:

```bash
just refresh-v131-trait-visualization
```

The same steps, one recipe at a time (`download-kg` fetches the KG-Microbe
edges and nodes that `taxa-trait-matrix` reads):

```bash
just download
just ingest
just study-tsvs
just taxa-trait-matrix
just study-trait-profiles
just trait-heatmaps
just trait-profile-investigation
just benchmark-bugsigdb-access
```

The visible figure outputs default to Homo sapiens experiments when
`data/raw/full_dump.csv` is available.

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
- BugSigDB main page: https://bugsigdb.org/Main_Page
- BugSigDBExports `v1.3.1`: https://github.com/waldronlab/BugSigDBExports/tree/v1.3.1
- Zenodo record: https://zenodo.org/records/19740674
- BugSigDB paper: https://doi.org/10.1038/s41587-023-01872-y
