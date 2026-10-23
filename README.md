# microbial-trait-profiling

`microbial-trait-profiling` integrates curated microbial taxon signatures with
trait annotations for comparative analysis across studies and biological contexts.

It downloads [BugSigDB](https://bugsigdb.org) data, parses it into
per-study YAML files validated against a [LinkML](https://linkml.io) schema,
and aggregates KG-Microbe trait annotations for paired comparisons and
enrichment analyses.

The Python package remains `bsdbng`; existing imports, schema identifiers,
environment variables, and release asset names retain that name for compatibility.

All biological claims, figure captions, and pipeline interpretations follow the
project-wide biological interpretation guide in [`methods/biological_interpretation.md`](methods/biological_interpretation.md)
and trait semantics in [`docs/table1_trait_semantics_audit.md`](docs/table1_trait_semantics_audit.md). In short:
results quantify reported taxon annotation patterns across literature cohorts,
not direct in vivo physiological or metabolic measurements.

## Quick start

### Prerequisites

You need two tools installed. Everything else is managed automatically — with
one exception: the figure recipes (`just guild-maps`, `just table1-trait-maps`)
render PNG files through a headless Chrome, fetched version-pinned by
`just get-chrome` (see `methods/README.md`).

| Tool | What it does | Install |
|------|-------------|---------|
| **[uv](https://docs.astral.sh/uv/)** | Python environment and dependency management | `curl -LsSf https://astral.sh/uv/install.sh \| sh` |
| **[just](https://just.systems)** | Command runner (like make, but simpler) | `cargo install just` or `brew install just` or [other options](https://just.systems/man/en/packages.html) |

> **Windows users**: both `uv` and `just` work on Windows.
> See [uv Windows install](https://docs.astral.sh/uv/getting-started/installation/#windows)
> and [just Windows install](https://just.systems/man/en/packages.html).
> If you use WSL, the Linux instructions apply.

### Clone and run

```bash
git clone https://github.com/CultureBotAI/microbial-trait-profiling.git
cd microbial-trait-profiling
uv sync --all-extras --group dev --group qa
just check       # run all checks (should pass on a fresh clone)
just download    # download BugSigDB data (~28 MB, ~6s)
just ingest      # parse into per-study YAML (~22s, 1900+ studies)
```

That's it. After `just ingest`, look at `data/studies/` for the output.

### BugSigDB data provenance

`just download` pins the ingest input to a citable release rather than a moving branch.

| field | value |
|---|---|
| source | [waldronlab/BugSigDBExports](https://github.com/waldronlab/BugSigDBExports) release **v1.3.1** |
| version DOI | [`10.5281/zenodo.19740674`](https://doi.org/10.5281/zenodo.19740674) |
| pinned by | commit SHA `5815c2638d6b56134ac8db19edc2bb3c4b17f1ca` — not the tag name, which is mutable |
| integrity | md5 `00df5a86c684eac19fd767c029ff0130`, 16,241,290 bytes, checked on every download |
| release published | 2026-04-24 |
| pinned file | `full_dump.csv` |

**Why v1.3.1 and not v1.3.0.** v1.3.1 is a *reviewed-only* export: every row carries
`State = Complete` and a named `Reviewer`, which v1.3.0 and `devel` do not. It therefore
ingests *fewer* studies than v1.3.0 despite being a year newer (1,245 vs 1,368). v1.3.0 is
unfiltered, 16 months older, and contains four studies present in no current source. See
issue #130 before changing this pin.

**What is and is not pinned.** The ingest reads only `full_dump.csv`, so pinning that file
pins the ingest. The three MediaWiki exports — `studies.csv`, `experiments.csv`,
`signatures.csv` — are live wiki tables with **no pinned equivalent** (the release does not
contain them). They are still downloaded, because `just stats` reports on them and
downstream analyses join against `experiments.csv`, but they are **not pinned**: two runs on
different days can fetch different vintages. Do not read "the ingest is pinned" as
"everything in `data/raw/` is pinned". Because they are not ingest inputs, an outage on them
is a warning rather than a failure — a dead export URL should not take down a pinned run.

**Every dump declares its own vintage.** Line 1 of `full_dump.csv` is a comment the CSV
reader strips:

```
# BugSigDB 2026-04-24_00:41_UTC, License: ..., URL: https://bugsigdb.org
```

That timestamp is parsed into `data/raw/download_provenance.json` alongside the resolved
release, DOI, tag SHA, per-file source URL, ETag and a per-file `pinned` flag — and the
resolved source is printed at the start of every run. The sidecar is **committed**, so the
vintage of a checkout is recoverable from version control rather than only from the machine
that ran the download. It is recorded for unpinned runs too: a `live` run cannot be repeated,
but it can at least be identified afterwards.

To fetch the moving `devel` branch instead — for a deliberate freshness check, never for a
reproducible run:

```bash
BSDBNG_BUGSIGDB_SOURCE=live just download
```

### If `just download` is unreachable: use the release raw-data tarball

Each tagged release attaches a `bsdbng-raw-data-v<tag>.tar.gz` asset containing:
`studies.csv`, `experiments.csv`, `signatures.csv`, `full_dump.csv`, and
`download_provenance.json`.

If BugSigDB's export URLs are down or returning errors, you can unpack this asset
into `data/raw/` and skip `just download`:

```bash
mkdir -p release-assets
gh release download v0.4.0 --repo CultureBotAI/microbial-trait-profiling \
  --pattern 'bsdbng-raw-data-*.tar.gz' --dir release-assets
tar xzf release-assets/bsdbng-raw-data-v0.4.0.tar.gz -C data/raw/
just ingest   # skip `just download`
```

Motivation: BugSigDB's CSV export URLs returned 404s for several hours on
2026-04-20 due to a `/w/` -> `/W/` redirect that went nowhere. The release
tarball was the fastest recovery path.

### Explore the data

Install [yq](https://github.com/mikefarah/yq) (Mike Farah's Go version,
not the Python one) to query YAML from the command line:

```bash
# Find a study by PMID
grep -l "pmid: 10202341" data/studies/*.yaml

# Show taxa in a study
yq '.experiments[].signatures[].taxa[].taxon_name' data/studies/bsdb_10202341.yaml

# Count studies
ls data/studies/*.yaml | wc -l
```

## Developer commands

All commands go through the Justfile. Run `just --list` to see them all.

| Command | What it does |
|---------|-------------|
| `just check` | Run all checks (ruff, mypy, pytest, linkml, bandit, deptry, pip-audit) |
| `just precommit` | Run the full pre-commit stack, including the installed-environment dependency audit |
| `just lock --upgrade-package PACKAGE` | Update a dependency with the declared uv version; run `just sync` first |
| `just secrets-scan` / `just secrets-audit` | Print candidate secret findings / review the committed baseline |
| `just download` | Download BugSigDB CSVs (with caching and retry) |
| `just ingest` | Parse CSVs into per-study YAML |
| `just taxa-trait-matrix` | Build the taxon x trait matrix from BugSigDB taxa and KG-Microbe |
| `just study-trait-profiles` | Sum taxon trait vectors by study, experiment, and increased/decreased direction, for human hosts, without traits no human profile carries: the canonical table committed as `data/derived/study_trait_profiles_wide.tsv.gz` (`just study-trait-profiles-check` compares them and their strain rollup records) |
| `just table1-compositional-stats` | Paired log-ratio tests for all 15 Table 1 trait pairs; `--body-site LABEL` for one exact stratum |
| `just pipeline` | Run download + ingest |
| `just clean-all` | Remove all generated artifacts |
| `just format` | Auto-fix lint and formatting |
| `just test` | Run tests only |
| `just lint` | Run linter only |

Every pipeline step is timed and logged to `logs/pipeline.log`.

The audited `.secrets.baseline` contains only the public BugSigDB source revision
and dump checksum in `src/bsdbng/download.py`. It records those exact values as
false positives; it does not exempt the file or new values from detection. Review
candidate findings before updating the baseline. The lock and secret-audit CLIs
are declared in the QA group at the same versions as their pre-commit hooks.

### Downloading KG-Microbe for trait enrichment

The `kg` optional extra (installed by `uv sync --all-extras`) pulls in
[`gdown`](https://github.com/wkentaro/gdown), a pure-Python Google Drive
downloader. Use it to fetch the merged KG-Microbe nodes and edges from
[Marcin Joachimiak's published build](https://drive.google.com/drive/folders/1-oxGDvJS5oWlwhT1bV_ns2KMeumVLnTa)
(January 2026 snapshot, file IDs pinned in `scripts/download_kg_microbe.py`):

```bash
just download-kg          # skips files whose provenance fingerprint matches
just download-kg-force    # redownload unconditionally
```

Outputs land in `data/kg-microbe/` (gitignored TSVs, versioned
`provenance.json` with file IDs, sizes, md5s, and fetch time). When a
newer build lands in a new Drive folder, update the file IDs at the top
of `scripts/download_kg_microbe.py` and rerun.

### Querying studies with linkml-store over MongoDB

The `mongo` optional extra (installed by `uv sync --all-extras`) pulls in
[`linkml-store`](https://github.com/linkml/linkml-store), `pymongo`, and
`python-dotenv`. Use it to query the ingested studies three different ways:
exact-match, trigram similarity, or OpenAI embedding semantic search.

**Prereqs:**

- MongoDB running locally on `mongodb://localhost:27017` (override with
  `BSDBNG_MONGO_URI=...`).
- `data/studies/` populated (run `just pipeline` first).
- For embedding search only: `OPENAI_API_KEY` in `.env` (see `.env.example`).

**One-time setup:**

```bash
just mongo-load         # load 1,900+ study YAMLs into MongoDB (~45s)
just index-trigram      # build trigram index (free, ~seconds)
just index-embeddings   # build OpenAI embedding index (~minutes, $)
```

**Queries:**

```bash
just query-field experiments.condition "Colorectal cancer"
just search-trigram "gut bacteria in elderly patients"
just search-embeddings "autoimmune conditions linked to oral bacteria"

# side-by-side demo: exact-match returns nothing, embeddings return
# five semantically-spot-on studies (Alzheimer's, Parkinson's, ...)
just demo-embeddings
```

The Justfile also ships `yq-*` recipes (direct yq over YAML files) and
`mongo-native-*` recipes (mongosh aggregations) for contrast with the
linkml-store approach.

## Contributing

### Workflow

1. Pick an [open issue](https://github.com/CultureBotAI/microbial-trait-profiling/issues)
2. Create a linked branch: `gh issue develop <number> --repo CultureBotAI/microbial-trait-profiling --checkout`
3. Make your changes
4. Run `just check` — **all checks must pass before opening a PR**
5. Commit and push
6. Open a PR that references the issue (e.g., "Closes #30")

### Using LLMs

LLM-generated code is welcome — the guardrails don't care who wrote the code.
But the code must pass `just check` before committing. The workflow is:

1. Generate code (Claude Code, Copilot, ChatGPT, whatever)
2. Run `just check`
3. Fix any failures
4. Commit

Do not commit code that hasn't passed `just check`. Do not skip checks.
If a check seems wrong, open an issue to discuss changing the check.

See [AGENTS.md](AGENTS.md) for the full policy on agent-assisted development,
including what to do when your agent works around a missing tool.

### File types

This project uses a minimal set of file types:

- **CSV** — raw BugSigDB input (downloaded, not committed)
- **YAML** — per-study output and schema definition
- **Python** — source code and tests
- **JSON** — provenance and ingest logs (generated, not committed)

No DuckDB, no Parquet, no Jupyter notebooks in this repo.

## Project structure

```
microbial-trait-profiling/
  schema/bsdbng.yaml              # LinkML schema (source of truth)
  src/bsdbng/
    download.py                    # BugSigDB CSV download with caching
    ingest.py                      # CSV → per-study YAML pipeline
    datamodel/bsdbng_pydantic.py   # generated Pydantic models (do not edit)
    paths.py                       # project path constants
  tests/                           # pytest tests
  data/
    raw/                           # downloaded CSVs (gitignored)
    studies/                       # per-study YAML output (gitignored)
    derived/                       # aggregate artifacts (gitignored)
  Justfile                         # developer commands
  pyproject.toml                   # project config
  logs/                            # pipeline timing logs (gitignored)
```

## Schema

The [execution-checked core workflow](docs/workflow.md) adds a tracked LinkML
specification for raw snapshot inputs, study collections, flattened tables,
trait matrices, distinct count/ORA profiles, coverage, and actual execution
records. Start with `just workflow-check` and `just core-fixture`; use
`just core-validate RUN` to recheck a completed bundle before downstream analysis.

The schema defines four entity types, all inheriting from `NamedThing`:

- **Study** — one published microbiome study (tree root)
- **Experiment** — a contrasted-group comparison within a study
- **Signature** — a set of taxa with a shared direction (increased/decreased)
- **Taxon** — a microbial organism identified by NCBITaxon CURIE

See `schema/bsdbng.yaml` for the full definition.

## Data quality

The ingest pipeline logs every data drop and transformation assumption to
`data/derived/ingest_log.json`. See [#44](https://github.com/CultureBotAI/microbial-trait-profiling/issues/44)
for the full audit of what gets discarded and why.

## Background

- [Project origins](docs/origin-story.md) — how this repo relates to the
  KG-Microbe-search prototype and how to talk about your contributions
- [BugSigDB data access methods](docs/bugsigdb-access-methods.md)
- [BugSigDB live counts](docs/bugsigdb-live-counts.md)
