set shell := ["bash", "-cu"]

schema := "schema/bsdbng.yaml"
lintconfig := "schema/.linkmllint.yaml"
logdir := "logs"

default: check

# --- timed wrapper (logs wall-clock seconds per step) ---

[private]
timed step +cmd:
	mkdir -p {{logdir}}
	start=$(date +%s) && \
	{{cmd}} && \
	elapsed=$(($(date +%s) - start)) && \
	echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) {{step}} ${elapsed}s" | tee -a {{logdir}}/pipeline.log

gen-pydantic:
	mkdir -p src/bsdbng/datamodel
	uv run gen-pydantic {{schema}} > src/bsdbng/datamodel/bsdbng_pydantic.py

gen-derived: gen-pydantic

check-generated:
	tmpfile="$(mktemp)" && \
	trap 'rm -f "$tmpfile"' EXIT && \
	uv run gen-pydantic {{schema}} > "$tmpfile" && \
	diff -u src/bsdbng/datamodel/bsdbng_pydantic.py "$tmpfile"

sync:
	uv sync --all-extras --group dev --group qa

lint:
	uv run ruff check src tests scripts
	uv run ruff format --check src tests scripts

format:
	uv run ruff check --fix src tests scripts
	uv run ruff format src tests scripts

trait-profile-investigation:
	uv run python scripts/regenerate_trait_profile_investigation.py

# Reported sample sizes of the experiments behind the saved oral heatmap.
oral-experiment-sizes raw_dump="data/raw/full_dump.csv":
    uv run python scripts/plot_oral_experiment_sample_sizes.py --raw-dump "{{raw_dump}}"

# Regenerate six-trait pattern and full-table signed direction heatmaps.
six-shared-trait-heatmaps:
	uv run python scripts/create_six_shared_trait_cytoscape_heatmaps.py

# Focused view of saved Jed contrasts; preserves their original Fisher/BH inference.
oral-six-trait-heatmap:
	uv run python scripts/create_six_shared_trait_cytoscape_heatmaps.py --oral-only --source-commit "$(git log -1 --format=%H -- trait-visualization/trait-profile-investigation/six-shared-traits/six_shared_trait_condition_body_site_direction_contrasts.tsv)"

# Trace PD lollipop traits to observed taxa and profile evidence.
pd-trait-taxa:
	uv run python scripts/lookup_pd_trait_taxa.py

# Plot species-level carrier counts for the four positive-side PD traits.
pd-trait-taxa-heatmap:
	uv run python scripts/create_pd_four_trait_taxon_heatmap.py

# Overlay increased (solid) and decreased (dashed) high-support Sankey links.
six-shared-trait-combined-sankey:
	uv run python scripts/create_six_shared_trait_sankeys.py --combined-high-support-only

# Redraw the saved increased, decreased and combined high-support selections.
six-shared-trait-high-support-sankeys:
	uv run python scripts/create_six_shared_trait_sankeys.py --high-support-layouts-only

# Flatten per-study YAML into data/derived/*.tsv.
study-tsvs:
	uv run python -m bsdbng.study_tsv

# Build data/derived/taxa_trait_matrix.tsv from BugSigDB taxa and KG-Microbe.
# Requires data/kg-microbe/merged-kg_{edges,nodes}.tsv; run `just download-kg` first.
taxa-trait-matrix:
	uv run python scripts/create_taxa_trait_matrix.py

# Sum taxon trait vectors by study, experiment, and increased/decreased profile.
study-trait-profiles:
	uv run python scripts/create_study_trait_profiles.py

trait-heatmaps:
	uv run python scripts/create_clustered_heatmap.py
	uv run python scripts/create_condition_taxon_heatmap.py
	uv run python scripts/create_taxa_guild_heatmap.py
	uv run python scripts/create_study_trait_profile_heatmap.py
	uv run python scripts/create_condition_trait_module_heatmap.py

trait-visualization:
	just download-kg
	just pipeline
	just study-tsvs
	just taxa-trait-matrix
	just study-trait-profiles
	just trait-heatmaps
	just trait-profile-investigation

refresh-v131-trait-visualization:
	just clean-all
	just trait-visualization

typecheck:
	uv run mypy src tests

test:
	uv run pytest

test-cov:
	uv run pytest --cov=bsdbng --cov-report=term-missing

install-hooks:
	uv run pre-commit install

# Benchmark the four non-R BugSigDB data access methods (raw Help:Export CSVs,
# BugSigDBExports full_dump, GMT bundle, Zenodo pinned release). Writes per-
# method TSVs to data/benchmarks/ and a summary to docs/bugsigdb-access-benchmark.md.
benchmark-bugsigdb-access:
	uv run python scripts/benchmark_bugsigdb_access.py

# Download the merged KG-Microbe nodes and edges TSVs from a pinned Google Drive
# folder (Marcin Joachimiak's published build) into data/kg-microbe/. Skips files
# whose provenance fingerprint still matches. Re-download with `just download-kg-force`.
download-kg:
	uv run --extra kg python scripts/download_kg_microbe.py

download-kg-force:
	uv run --extra kg python scripts/download_kg_microbe.py --force

security:
	uv run bandit -c pyproject.toml -r src scripts
	uv run deptry src
	uv run pip-audit

precommit:
	uv run pre-commit run --all-files

validate-schema:
	uv run linkml validate {{schema}}

lint-schema:
	uv run linkml lint --config {{lintconfig}} --validate --ignore-warnings {{schema}}
	uv run linkml lint --config {{lintconfig}} --ignore-warnings {{schema}}

check-schema: validate-schema lint-schema

sync-agent-docs:
	uv run python scripts/sync_agent_docs.py

check-agent-docs:
	uv run python scripts/sync_agent_docs.py
	git diff --exit-code .github/copilot-instructions.md

check: check-generated check-agent-docs lint typecheck test check-schema security

# --- pipeline steps (timed, logged) ---

download:
	just timed download "uv run python -c 'from bsdbng.download import download_exports; download_exports()'"

ingest:
	just timed ingest "uv run python -c 'from bsdbng.ingest import ingest; ingest()'"

validate-studies:
	just timed validate-studies "uv run python -c 'from bsdbng.validate import validate_all_studies; validate_all_studies()'"

validate-studies-linkml:
	just timed validate-studies-linkml "uv run linkml validate -s {{schema}} -C Study data/studies/*.yaml"

stats:
	@uv run python -c "from bsdbng.stats import print_stats; print_stats()"

pipeline: download ingest

# --- linkml-store over MongoDB (localhost:27017/bsdbng.studies) ---
# Run once to populate the collection:
mongo-load:
	just timed mongo-load "uv run python scripts/linkml_store_load.py"

# Wipe the studies collection AND all its indexes, then reload:
mongo-reload:
	just timed mongo-reload "uv run python scripts/linkml_store_load.py --drop"

# Build a trigram (substring) index -- no API key, ~seconds:
index-trigram:
	just timed index-trigram "uv run python scripts/linkml_store_index.py --kind simple"

# Build an OpenAI embedding index -- requires OPENAI_API_KEY, ~minutes:
index-embeddings:
	just timed index-embeddings "uv run python scripts/linkml_store_index.py --kind llm"

# Rebuild an index (drop + recreate) -- use after re-loading data:
rebuild-trigram:
	uv run python scripts/linkml_store_index.py --kind simple --rebuild

rebuild-embeddings:
	uv run python scripts/linkml_store_index.py --kind llm --rebuild

# --- linkml-store queries ---

# Exact match on a field path (nested paths use dots):
#   just query-field experiments.condition "Colorectal cancer"
#   just query-field experiments.condition_ontology_id "EFO:0005842"
query-field FIELD VALUE:
	uv run python scripts/linkml_store_query.py --field {{FIELD}} --value "{{VALUE}}"

# Trigram substring search:
#   just search-trigram "gut bacteria in elderly patients"
search-trigram QUERY:
	uv run python scripts/linkml_store_query.py --search "{{QUERY}}"

# Embedding semantic search:
#   just search-embeddings "autoimmune conditions linked to oral bacteria"
search-embeddings QUERY:
	uv run python scripts/linkml_store_query.py --embed "{{QUERY}}"

# Canonical demo of what embeddings buy you: a query whose top-5 hits are
# impossible to retrieve with exact-match (no study contains the phrase
# "diseases of aging" -- the matches are Alzheimer's, Parkinson's,
# age-related sepsis, aging-mouse obesity, saliva-age signatures).
demo-embeddings:
	@echo "--- baseline: exact match on 'diseases of aging' (returns nothing) ---"
	@just query-field title "diseases of aging" || true
	@echo
	@echo "--- embeddings: semantically related studies ---"
	@just search-embeddings "diseases of aging"

# --- yq: direct queries over YAML files (no DB, no index) ---

# Find all studies whose title matches a substring (case-insensitive):
#   just yq-title-search colorectal
yq-title-search SUBSTRING:
	@SUB="{{SUBSTRING}}" yq -r 'select((.title // "") | test("(?i)" + strenv(SUB))) | [.id, .title] | @tsv' \
	  data/studies/*.yaml | grep -v '^$' | head -20

# Find all studies where any experiment has a given condition (exact match):
#   just yq-condition "Parkinson's disease"
yq-condition CONDITION:
	@COND="{{CONDITION}}" yq -r 'select(.experiments[].condition == strenv(COND)) | [.id, .title] | @tsv' \
	  data/studies/*.yaml | grep -v '^$' | head -20

# --- native MongoDB (mongosh) ---

# Count studies by host species via aggregation pipeline:
mongo-native-hosts:
	@mongosh --quiet bsdbng --eval '\
	  db.studies.aggregate([ \
	    {$unwind: "$experiments"}, \
	    {$group: {_id: "$experiments.host_species", n: {$sum: 1}}}, \
	    {$sort: {n: -1}}, \
	    {$limit: 10} \
	  ]).forEach(d => print(d._id + "\t" + d.n))'

# Fetch one study by id:
#   just mongo-native-one bsdb:31182740
mongo-native-one ID:
	@mongosh --quiet bsdbng --eval 'printjson(db.studies.findOne({id: "{{ID}}"}, {id:1, title:1, "experiments.condition":1}))'

# --- cleanup ---

clean-raw:
	rm -f data/raw/*.csv data/raw/download_provenance.json

clean-studies:
	rm -f data/studies/*.yaml

clean-derived:
	find data/derived -type f ! -name .gitkeep -delete

clean-logs:
	rm -f {{logdir}}/pipeline.log

clean-all: clean-raw clean-studies clean-derived clean-logs
