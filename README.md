# bsdbng

`bsdbng` is a clean-room BugSigDB ingestion and normalization project.

It downloads [BugSigDB](https://bugsigdb.org) data, parses it into
per-study YAML files validated against a [LinkML](https://linkml.io) schema,
and provides a foundation for trait aggregation and phenotype cross-referencing.

Trait results describe KG-Microbe annotation representation among taxa in curated
BugSigDB signatures. They do not directly measure patient-level trait abundance,
pathway activity, environmental exposure or causation. Result descriptions
are governed by the [biological interpretation guide](methods/biological_interpretation.md);
the [claim audit](methods/biological_interpretation_audit.md) records corrections,
proposed revisions, and outstanding evidence needs across current and proposed outputs.

## Quick start

### Prerequisites

You need two tools installed. Everything else is managed automatically.

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
git clone https://github.com/turbomam/bsdbng.git
cd bsdbng
uv sync --all-extras --group dev --group qa
just check       # run all checks (should pass on a fresh clone)
just download    # download BugSigDB data (~28 MB, ~6s)
just ingest      # parse into per-study YAML (~22s, 1900+ studies)
```

That's it. After `just ingest`, look at `data/studies/` for the output.

### If `just download` is unreachable: use the release raw-data tarball

Each tagged release attaches a `bsdbng-raw-data-v<tag>.tar.gz` asset containing:
`studies.csv`, `experiments.csv`, `signatures.csv`, `full_dump.csv`, and
`download_provenance.json`.

If BugSigDB's export URLs are down or returning errors, you can unpack this asset
into `data/raw/` and skip `just download`:

```bash
mkdir -p release-assets
gh release download v0.4.0 --repo turbomam/bsdbng \
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
| `just download` | Download BugSigDB CSVs (with caching and retry) |
| `just ingest` | Parse CSVs into per-study YAML |
| `just taxa-trait-matrix` | Build the taxon x trait matrix from BugSigDB taxa and KG-Microbe |
| `just study-trait-profiles` | Sum taxon trait vectors by study, experiment, and increased/decreased direction |
| `just pipeline` | Run download + ingest |
| `just clean-all` | Remove all generated artifacts |
| `just format` | Auto-fix lint and formatting |
| `just test` | Run tests only |
| `just lint` | Run linter only |

Every pipeline step is timed and logged to `logs/pipeline.log`.

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

1. Pick an [open issue](https://github.com/turbomam/bsdbng/issues)
2. Create a branch: `git checkout -b issue-<number>-<short-slug>`
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
bsdbng/
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

The schema defines four entity types, all inheriting from `NamedThing`:

- **Study** — one published microbiome study (tree root)
- **Experiment** — a contrasted-group comparison within a study
- **Signature** — a set of taxa with a shared direction (increased/decreased)
- **Taxon** — a microbial organism identified by NCBITaxon CURIE

See `schema/bsdbng.yaml` for the full definition.

## Data quality

The ingest pipeline logs every data drop and transformation assumption to
`data/derived/ingest_log.json`. See [#44](https://github.com/turbomam/bsdbng/issues/44)
for the full audit of what gets discarded and why.

## Background

- [Where bsdbng came from](docs/origin-story.md) — how this repo relates to the
  KG-Microbe-search prototype and how to talk about your contributions
- [BugSigDB data access methods](docs/bugsigdb-access-methods.md)
- [BugSigDB live counts](docs/bugsigdb-live-counts.md)
