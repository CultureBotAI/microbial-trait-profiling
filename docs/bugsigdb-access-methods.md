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
- human-host figure corpus: 5,629 profiles from 1,047 studies, 16S and shotgun
  (`WMS`) experiments only (#742; 5,829 from 1,080 before)

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
Study -> Experiment -> Signature -> Taxa hierarchy. Join the raw tables using
their study and experiment page references; `BSDB ID` is constructed later by
the flat-file exporter, rather than supplied in these raw headers.

bsdbng downloads these files from the current S3-backed `Help:Export` links,
but `ingest.py` still raises `NotImplementedError` for separate raw-table
merging when `full_dump.csv` is absent.

### 2. BugSigDBExports Tagged Full Dump

`waldronlab/BugSigDBExports` is a derived export layer built from the BugSigDB
tables. Its `full_dump.csv` is pre-joined and convenient for bsdbng's current
ingestion code.

For reproducible branch artifacts, use the `v1.3.1` tag rather than `devel`:

- `v1.3.1` is reviewed/complete enough for citation and has a DOI.
- `devel` is moving; cite a resolved commit and file checksum, not the branch name alone.
- `v1.3.0` should not be used for new work; it is older and was only named in
  earlier scripts/issues because it was current when that code was written.

### 3. Zenodo Release

Zenodo is the authoritative citation target for the pinned release. The branch
stores the DOI and tag SHA in `data/raw/download_provenance.json` when
`download_exports(force=True)` refreshes `full_dump.csv`.

### 4. `bugsigdbr` R Package

The Bioconductor `bugsigdbr` interfaces remain **untested** in this comparison
(#62): no R installation, API execution, or R timing was performed. Reading the
exporter's R source explains its transformations; it does not validate the R
client as an alternative ingestion interface.

## Decision matrix (#87)

Evidence checked on 2026-11-10: bounded CSV headers, a GMT record, the pinned
exporter, the upstream workflow, and release metadata. This is a documentary
comparison, not a new download-speed benchmark or a row-level completeness audit.
The [existing timings](bugsigdb-access-benchmark.md) remain a single historical run.

| Access path | Fidelity | Freshness cadence | Field completeness | Endpoint stability |
| --- | --- | --- | --- | --- |
| `Help:Export` raw `studies.csv` / `experiments.csv` / `signatures.csv` | Closest exported representation of authoritative MediaWiki tables; requires an explicit relational join and curation-state policy, not a revision-history archive. | Live database exports refreshed hourly, not an instantaneous query of every edit [S1]. | Separate study, experiment and signature fields, including Abstract and curation fields; inventory below [S2]. | `Special:Ask` queries depend on wiki properties/query configuration; `Help:Export` is the discovery page. Its S3 URLs can change: the current host differs from the pinned exporter's URLs [S1, S3]. |
| `BugSigDBExports/full_dump.csv` | Derived flat join; filters Complete records, synchronizes parent records, removes Abstract, and constructs a signature-level BSDB ID [S3]. | `devel` workflow scheduled hourly (minute 55 UTC); not a service-level guarantee. Failed jobs or unchanged output need not yield a new data commit [S6]. | Joined study/experiment/signature fields; retained signature curation fields do not imply preservation of study/experiment review metadata [S3, S4]. | Moving `devel` raw paths are convenient discovery, not a fixed dataset. Resolve a commit or release and verify a file checksum. |
| `BugSigDBExports/*.gmt` | Derived signature sets of microbial taxa (GMT's gene-set format); rank/identifier selection changes membership [S3, S5]. | Produced by the same scheduled export as the flat CSV; a chosen commit is fixed. | Set name, description slot and taxon members; full study and experiment metadata are stripped. Use a same-revision CSV for context [S5]. | Branch filenames can move; pin commit plus the chosen rank/ID/exactness filename. The extension alone does not specify the set universe. |
| Zenodo pinned release `10.5281/zenodo.10627578` | Release-time snapshot of derived exports, not a live MediaWiki mirror [S7]. | Fixed v1.2.1 release, issued 2024-02-06; immutable file version, with no automatic freshness [S8]. | Contains full_dump.csv, GMT variants and README; its CSV lacks Data transformation, present in the later v1.3.1 header [S7, S9]. | Permanent version DOI identifies this release. Resolve file links/checksums from its record; distinguish a version DOI from the all-versions concept DOI. |

**Release identity:** #87's DOI is **v1.2.1**, not the v1.3.1 DOI in the
current project decision. Its CSV declares export vintage `2024-02-04_18:58_UTC`,
distinct from the release date. The v1.3.1 CSV declares `2026-04-24_00:41_UTC`.
The older release is a historical-reproduction option, not a proposed replacement
for the project's current pin. Metadata and header inspection do not establish
equal row counts, equal content, or equal taxonomic coverage across releases [S4, S7–S9].

## Field-by-field completeness

The table inventories all headers observed in the three raw CSVs on 2026-11-10
and the pinned v1.3.1 CSV. “Retained” means a column exists, **not** that every
row has a value. This does not project today's raw schema backward onto v1.2.1.
Raw headers are mutable evidence [S2]; the flat-file header and transformation
code are commit-pinned [S3, S4].

| Raw table / fields | v1.3.1 full_dump.csv treatment | GMT / v1.2.1 distinction |
| --- | --- | --- |
| Study: `Study page name` | Renamed `Study`; repeated across joined experiment/signature rows. | May contribute to encoded set identity; no separate study table. |
| Study: `Study design`, `PMID`, `DOI`, `URL`, `Authors list`, `Title`, `Journal`, `Year`, `Keywords` | Retained. | These are not separate GMT columns; present in the v1.2.1 CSV header. |
| Study: `Abstract` | Explicitly removed by exporter. | No abstract field in either inspected flat CSV or GMT. |
| Study: `Curation method`, `Curation protocol`, `submittedBy`, `submissionDate`, `Review state`, `Quality control tag` | In current raw header, absent from the pinned v1.3.1 header. This comparison alone does not identify when they were added. | Also absent from the inspected v1.2.1 header and GMT. |
| Study and experiment: `State`, `Reviewer` | State used to filter Complete records; these levels' State/Reviewer columns removed before joining. | Cannot recover their distinct values from GMT or the flat CSV. |
| Experiment: `Experiment page name`, `Study` | First renamed `Experiment`; second joins to the study record. | Set name can encode experiment identity; not separate GMT columns. |
| Experiment: `Location of subjects`, `Host species`, `Body site`, `UBERON ID`, `Condition`, `EFO ID` | Retained; exporter attempts case normalization for labels sharing ontology IDs. | Context is not a structured GMT metadata table. |
| Experiment: `Group 0 name`, `Group 1 name`, `Group 1 definition`, `Group 0 sample size`, `Group 1 sample size` | Retained. | Names may be encoded in a set name; definitions and sizes require the CSV. |
| Experiment: `Antibiotics exclusion`, `Sequencing type`, `16S variable region`, `Sequencing platform` | Retained. | Not structured GMT fields. |
| Experiment: `Data transformation` | Present in raw and v1.3.1. | **absent from v1.2.1** CSV; not a GMT field. |
| Experiment: `Statistical test`, `Significance threshold`, `MHT correction`, `LDA Score above`, `Matched on`, `Confounders controlled for` | Retained. | Statistical-method and adjustment context require the CSV, not GMT. |
| Experiment: `Pielou`, `Shannon`, `Chao1`, `Simpson`, `Inverse Simpson`, `Richness` | Retained. | Not GMT metadata columns. |
| Signature: `Signature page name`, `Experiment`, `Study` | Signature name retained; Study/Experiment keys used for joining. | Set name carries signature identity; CSV needed to reconstruct explicit hierarchy. |
| Signature: `Source`, `Curated date`, `Curator`, `Revision editor`, `Description`, `State`, `Reviewer` | Retained from the signature level, after Complete filtering. | Description slot in inspected GMT is `NA`; it is not the signature's Description column. |
| Signature: `Abundance in Group 1`, `MetaPhlAn taxon names`, `NCBI Taxonomy IDs` | Retained; direction remains relative to the named comparison. | Direction can be encoded in set name; members depend on taxonomic rank, identifier type and exactness selection. |
| Derived: `BSDB ID` | Constructed from study/experiment/signature identifiers; not an observed column of the three raw headers. | Encoded in set name, rather than a separate column. |

The v1.2.1 header otherwise matches the enumerated flat-file fields; this is a
header comparison only [S9]. GMT should not be used to infer missing patient,
cohort or statistical metadata from a taxon list. Source direction and taxon
membership do not measure microbial activity.

## Use-case verdicts

- **Canonical ingestion / reproducible figures:** retain the project's pinned
  v1.3.1 CSV and recorded checksum. It fits the existing Python join contract;
  changing access paths would also change the population/filtering contract.
- **Live curation audit:** choose the three raw exports and save retrieval time,
  headers and hashes. Join and validate hierarchy explicitly; preserve excluded
  records and their State values. Separate-table ingestion remains tracked in #64.
- **Signature-set analysis:** use a pinned GMT variant when set membership is the
  required unit; use its same-revision CSV for study design, comparison direction
  and statistical context. Do not treat alternative GMT rank projections as
  independent studies.
- **Historical reproduction:** use the exact version DOI, such as v1.2.1
  `10.5281/zenodo.10627578`, and validate the record's file checksum. Its age is
  appropriate for reproducing that release, not for a current curation census.
- **R/API evaluation:** `bugsigdbr` remains untested under #62. Neither speed nor
  fidelity parity with the Python path is established by this record. Live
  `Special:Ask` query correctness/limits were not benchmarked here.

## Evidence for the decision record

- **S1:** [Help:Export, revision 259555](https://bugsigdb.org/w/index.php?title=Help:Export&oldid=259555), including its hourly-refresh statement and download links.
- **S2:** Header-only inspection of the linked [studies.csv](https://bugsigdb-csv.s3.us-east-va.perf.cloud.ovh.us/studies.csv), [experiments.csv](https://bugsigdb-csv.s3.us-east-va.perf.cloud.ovh.us/experiments.csv), and [signatures.csv](https://bugsigdb-csv.s3.us-east-va.perf.cloud.ovh.us/signatures.csv); field inventory above records the observation, not an immutable upstream schema.
- **S3:** [v1.3.1 exporter source](https://github.com/waldronlab/BugSigDBExports/blob/5815c2638d6b56134ac8db19edc2bb3c4b17f1ca/inst/scripts/dump_release.R): readFiles, addID, validation and GMT writing. Source inspected, R not executed.
- **S4:** [v1.3.1 CSV header](https://raw.githubusercontent.com/waldronlab/BugSigDBExports/5815c2638d6b56134ac8db19edc2bb3c4b17f1ca/full_dump.csv), first two lines only.
- **S5:** [Pinned mixed-rank NCBI GMT](https://raw.githubusercontent.com/waldronlab/BugSigDBExports/5815c2638d6b56134ac8db19edc2bb3c4b17f1ca/bugsigdb_signatures_mixed_ncbi.gmt), header and first record only; the description slot is NA and the remaining tokens are taxon IDs.
- **S6:** [Export workflow at inspected devel commit](https://github.com/waldronlab/BugSigDBExports/blob/2523676c1bb667e524dad5a8b76bc0a97bbb4bbf/.github/workflows/export-bugsigdb.yml), scheduled hourly; cadence is configuration, not measured uptime.
- **S7:** [Zenodo record metadata and file inventory](https://zenodo.org/api/records/10627578); v1.2.1 full_dump.csv lists MD5 `eb639bc00414bc6db02351c5e72b7da3`. Full release payload was not downloaded or rehashed in this review.
- **S8:** [DataCite DOI metadata](https://api.datacite.org/dois/10.5281/zenodo.10627578), title/version and issued date; concept DOI is `10.5281/zenodo.5606165`.
- **S9:** [v1.2.1 CSV header](https://zenodo.org/records/10627578/files/full_dump.csv), first two lines only.

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

Historical ingestion-schema accounting (not the access-path column inventory above):

- **45 captured** in `schema/bsdbng.yaml` (study, experiment, signature, taxon levels)
- **6 not captured** (curatorial metadata): Signature page name, Curated date, Curator, Revision editor, State, Reviewer
- **State column** should be used to filter: only ingest Complete records
  (not yet implemented — tracked in issue #64)

## Current approach in microbial-trait-profiling

`microbial-trait-profiling` currently uses `full_dump.csv` (from BugSigDBExports) as the
primary input, with the three raw CSVs downloaded but unused. The merge
from separate CSVs is a `NotImplementedError` stub.

**Planned change (issue #64):** switch to using the three separate CSVs
as primary input, implementing the join in Python, so filtering decisions and
excluded records are auditable. Historical count differences do not establish
how many additional usable studies a current join would retain.

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
