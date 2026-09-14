# Project origins

## The short version

`microbial-trait-profiling`, originally named `bsdbng`, is a clean-room
reimplementation of work that was prototyped in
[CultureBotAI/KG-Microbe-search](https://github.com/CultureBotAI/KG-Microbe-search).
The prototyping was done by undergraduate researchers. Their exploratory
work proved the concept and revealed the data structure; this project packages
those findings into a pipeline with a formal schema and automated validation.

## What the prototyping proved

Each finding below links to the specific prototype artifact that
demonstrated it, and the current project's code that builds on it.

### 1. BugSigDB data can be structured into per-study records

The prototype organized studies into per-study directories with nested
experiments, signatures, and taxa. This hierarchy was discovered through
hands-on data exploration.

- Prototype: [colorectal-cancer-studies/](https://github.com/CultureBotAI/KG-Microbe-search/tree/main/colorectal-cancer-studies) —
  per-study folders with `analysis_results.yaml`, metadata, and taxa files
- This project: `schema/bsdbng.yaml` — the Study → Experiment → Signature → Taxon
  hierarchy formalized as a LinkML schema

### 2. LinkML schemas can validate BugSigDB data

A schema was written and used with `linkml validate` to check study
outputs against a formal specification. Valid and invalid example files
were created to test the schema.

- Prototype: [schema/bugsigdb_kgmicrobe_schema.yaml](https://github.com/CultureBotAI/KG-Microbe-search/blob/main/schema/bugsigdb_kgmicrobe_schema.yaml) —
  the first BugSigDB LinkML schema
- Prototype: [schema/valid example output yaml/](https://github.com/CultureBotAI/KG-Microbe-search/tree/main/schema/valid%20example%20output%20yaml) and
  [schema/invalid example output yaml/](https://github.com/CultureBotAI/KG-Microbe-search/tree/main/schema/invalid%20example%20output%20yaml) —
  test data files
- This project: `schema/bsdbng.yaml` — refined schema with NamedThing inheritance,
  id_prefixes, per-class patterns, and Pydantic validation on write

### 3. Trait profile tables can be generated from structured study data

A script was built to produce wide and long TSV formats summarizing
taxa and their abundance directions across studies.

- Prototype: [Studies/reshape_trait_profiles.py](https://github.com/CultureBotAI/KG-Microbe-search/blob/main/Studies/reshape_trait_profiles.py)
  and the universal generator on branch
  [updated-noel-tsv](https://github.com/CultureBotAI/KG-Microbe-search/tree/updated-noel-tsv)
- This project: issue [#30](https://github.com/CultureBotAI/microbial-trait-profiling/issues/30) —
  the same capability, reading from schema-validated YAML

### 4. KG-Microbe phenotype data can be cross-referenced with BugSigDB taxa

An analysis pipeline looked up taxa in KG-Microbe's knowledge graph
to find phenotypic traits associated with differentially abundant organisms.

- Prototype: [scripts/analyze_taxa.py](https://github.com/CultureBotAI/KG-Microbe-search/blob/main/scripts/analyze_taxa.py) —
  taxa lookup against KG-Microbe edges.tsv and nodes.tsv
- Prototype: [bugsigdb-analysis-prompt-template.md](https://github.com/CultureBotAI/KG-Microbe-search/blob/main/bugsigdb-analysis-prompt-template.md) —
  LLM prompt for BugSigDB taxa extraction and KG-Microbe analysis
- This project: issue [#31](https://github.com/CultureBotAI/microbial-trait-profiling/issues/31) —
  deterministic version of the same cross-referencing

### 5. BugSigDB uses MetaPhlAn lineage encoding

Pipes separate taxonomic levels within one taxon, semicolons separate
distinct taxa. This is not documented in BugSigDB's own exports and was
learned through data exploration.

- Prototype: discovered through analysis of the CSV data and the
  [MetaPhlAn taxon names column](https://github.com/CultureBotAI/KG-Microbe-search/blob/main/schema/YAML_CONVERSION_GUIDE.md)
- This project: `src/bsdbng/ingest.py`, function `_parse_taxa` — implements
  correct semicolon/pipe parsing with the `METAPHLAN_RANK_MAP`

### 6. 1,743 BugSigDB studies were successfully analyzed

The prototype validated the approach at scale using LLM-driven analysis.

- Prototype: [PR #41](https://github.com/CultureBotAI/KG-Microbe-search/pull/41) —
  1,743 study analysis YAML files validated against the schema
- This project: `just pipeline` produces 1,923 studies deterministically in ~22s

## Archived prototype branches

The prototype branches in KG-Microbe-search were tagged and deleted
after their code PRs were merged. The tags preserve the exact commit
history. To reconstitute any branch:

```bash
# List archive tags
git tag -l 'archive/*'

# Inspect a tag
git log --oneline archive/updated-noel-tsv

# Recreate a branch from a tag
git checkout -b restored-noel-tsv archive/updated-noel-tsv

# Browse files at that point without creating a branch
git show archive/updated-noel-tsv:generate_trait_tsv_from_yaml.py
```

| Tag | Original branch | What it contains |
|-----|----------------|------------------|
| `archive/content-pr41-yaml-outputs` | `content-pr41-bugsigdb-analysis-yaml` | ~100 generated study YAML files (PR #46) |
| `archive/content-pr42-curated-fixes` | `content-pr42-curated-study-fixes` | Curated epilepsy/skin-wound study updates (PR #47) |
| `archive/updated-noel-tsv` | `updated-noel-tsv` | YAML-based trait profile TSV generator (PR #51) |
| `archive/noel-trait-tsv-orphan` | `noel_trait_tsv` | Early exploratory work: trait profiles, GTDB/NCBI tables |

## What this project adds

`microbial-trait-profiling` takes the proven concepts and packages them with:

- **A single LinkML schema** as the source of truth, with generated Pydantic
  models enforcing validation on every write
- **Automated checks** (ruff, mypy strict, pytest, linkml lint/validate,
  bandit, deptry, pip-audit) running in CI
- **Deterministic ingestion** — the same input always produces the same
  output, with no LLM in the critical path
- **Comprehensive logging** of every data drop and transformation assumption
- **Correct taxon parsing** based on the MetaPhlAn lineage encoding
  discovered during prototyping
- **HTTP caching** with ETag/If-Modified-Since for efficient re-downloads
- **A reproducible command surface** via Justfile
