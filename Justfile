set shell := ["bash", "-cu"]

# Every recipe draws with matplotlib's headless Agg backend. The macOS backend lays a figure
# out at the display's pixel ratio, so the same command wrote different bytes in different
# sessions (#719).
export MPLBACKEND := "Agg"

schema := "schema/mtp.yaml"
legacy_schema := "schema/bsdbng.yaml"
lintconfig := "schema/.linkmllint.yaml"
logdir := "logs"
p_threshold := "0.001"
min_fold := "1.5"
code_paths := "src scripts methods schema Justfile pyproject.toml uv.lock"

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
	mkdir -p src/mtp/datamodel src/bsdbng/datamodel
	uv run gen-pydantic {{schema}} > src/mtp/datamodel/mtp_pydantic.py
	uv run gen-pydantic {{legacy_schema}} > src/bsdbng/datamodel/bsdbng_pydantic.py
	uv run gen-pydantic schema/workflow.yaml > src/mtp/datamodel/workflow_pydantic.py
	uv run gen-pydantic schema/core_products.yaml > src/mtp/datamodel/core_products_pydantic.py
	uv run gen-pydantic schema/multiple_testing.yaml > src/mtp/datamodel/multiple_testing_pydantic.py

gen-derived: gen-pydantic

check-generated:
	tmpfile="$(mktemp)" && \
	trap 'rm -f "$tmpfile"' EXIT && \
	uv run gen-pydantic {{schema}} > "$tmpfile" && \
	diff -u src/mtp/datamodel/mtp_pydantic.py "$tmpfile" && \
	uv run gen-pydantic {{legacy_schema}} > "$tmpfile" && \
	diff -u src/bsdbng/datamodel/bsdbng_pydantic.py "$tmpfile" && \
	uv run gen-pydantic schema/workflow.yaml > "$tmpfile" && \
	diff -u src/mtp/datamodel/workflow_pydantic.py "$tmpfile" && \
	uv run gen-pydantic schema/core_products.yaml > "$tmpfile" && \
	diff -u src/mtp/datamodel/core_products_pydantic.py "$tmpfile" && \
	uv run gen-pydantic schema/multiple_testing.yaml > "$tmpfile" && \
	diff -u src/mtp/datamodel/multiple_testing_pydantic.py "$tmpfile"

# Validate the executable specification and its graph without reading large inputs.
workflow-check spec="workflows/bugsigdb_core.yaml":
	uv run python -m mtp.workflow check "{{spec}}"
	uv run linkml validate -s schema/workflow.yaml -C WorkflowDefinition "{{spec}}"
	uv run linkml lint --config {{lintconfig}} --ignore-warnings schema/workflow.yaml
	uv run linkml lint --config {{lintconfig}} --ignore-warnings schema/core_products.yaml
	uv run linkml lint --config {{lintconfig}} --ignore-warnings schema/multiple_testing.yaml
	uv run linkml validate -s schema/multiple_testing.yaml -C MultipleTestingPolicy methods/multiple_testing_policy_v1.yaml
	uv run linkml validate -s schema/multiple_testing.yaml -C HypothesisFamily tests/fixtures/multiple_testing_family.yaml

workflow-describe:
	uv run python -m mtp.workflow describe workflows/bugsigdb_core.yaml > docs/workflow-products.md
	uv run python -m mtp.workflow describe workflows/bugsigdb_core_subset.yaml > docs/workflow-subset-products.md

workflow-subset-check:
	uv run python -m mtp.workflow check workflows/bugsigdb_core_subset.yaml
	uv run linkml validate -s schema/workflow.yaml -C WorkflowDefinition workflows/bugsigdb_core_subset.yaml

# Explicit bindings, e.g. --input raw_dump=/path/full_dump.csv --output-root /path/runs.
[positional-arguments]
core-run *ARGS:
	uv run python -m mtp.workflow run "$@"

core-validate run:
	uv run python -m mtp.workflow validate "{{run}}"

[positional-arguments]
core-lineage run product *ARGS:
	uv run python -m mtp.workflow lineage "$@"

# A network-free execution from raw fixture inputs through both profile units.
core-fixture output=".tmp/core-fixture":
	uv run python -m mtp.workflow fixture --output-root "{{output}}"

# Historical claims remain explicitly unverified; no execution history is invented.
core-import-legacy source output:
	uv run python -m mtp.workflow import-legacy "{{source}}" "{{output}}"

core-compare run baseline_dir output:
	uv run python -m mtp.workflow compare-legacy "{{run}}" "{{baseline_dir}}" "{{output}}"

core-screen run output:
	uv run python -m mtp.trait_screen --core-run "{{run}}" --output-dir "{{output}}"

core-composition run output:
	uv run python -m mtp.trait_composition --core-run "{{run}}" --output-dir "{{output}}"

# After study-trait-profiles and trait-heatmaps: pooled/site tables with hashed provenance.
[positional-arguments]
trait-profile-site-stats *ARGS:
	uv run python scripts/create_trait_profile_infographic.py --stratified-stats-only "$@"

# Synthetic argument-guard probe using a copied producer (#950).
[positional-arguments]
probe-site-stats-guard *ARGS:
	uv run python scripts/adhoc/issue928_guard_isolation.py "$@"

# Complete exports carry the verified input fingerprints needed for integration.
[positional-arguments]
core-ora run output tail="upper":
	uv run python scripts/trait_enrichment.py --core-run "$1" \
	  --output "$2/report.txt" --tsv-output "$2/results.tsv" \
	  --all-results-output "$2/all_results.tsv" --tail "$3" \
	  --p-threshold 0.001 --min-fold 1.5 --max-fold 0.5 --fdr-scope global

[positional-arguments]
audit-ora-refresh manifest output:
	uv run python scripts/adhoc/issue720_audit_ora_refresh.py --manifest "$1" --output "$2"

# Bind all anatomical partitions and summaries to the same verified ORA core.
[positional-arguments]
core-strata run output *ARGS:
	uv run python scripts/run_stratified_enrichment.py --core-run "$1" --output-dir "$2" "${@:3}"

[positional-arguments]
audit-core-strata historical current:
	uv run python scripts/adhoc/issue727_compare_core_strata.py "$1" "$2"

core-verify-consumer record:
	uv run python -m mtp.workflow verify-consumer "{{record}}"

# Retain original records and declared dependencies; --retain adds explicit supporting evidence.
# A destination under the system temporary directory needs --allow-temporary-output (#1049).
[positional-arguments]
evidence-archive output *ARGS:
	uv run python -m mtp.workflow archive --output "$1" "${@:2}"

evidence-verify archive:
	uv run python -m mtp.workflow verify-archive "{{archive}}"

evidence-resolve archive original:
	uv run python -m mtp.workflow archive-resolve "{{archive}}" "{{original}}"

# Its strain rollup scope is the primary phenotype unless --strain-rollup off or full asks
# for the sensitivity or exploratory core (#788).
# Build one shared raw-to-core execution before rerunning any owner analysis.
# The rerun recipes refuse an output root under the system temporary directory unless
# --allow-temporary-output marks a throwaway run (#1049).
[positional-arguments]
rerun-core output raw edges nodes *ARGS:
	uv run python -m mtp.workflow prepare-core "$@"

# Standard ORA, depletion, rankings, strata, calibration and optional figures.
[positional-arguments]
rerun-ilyas run output *ARGS:
	uv run python -m mtp.workflow.arms ilyas "$@"

# Standard directional screen, calibration, orientation and evidence products.
[positional-arguments]
rerun-noel run output *ARGS:
	uv run python -m mtp.workflow.arms noel "$@"

# Standard paired composition, sensitivity, calibration, audits and optional figures.
[positional-arguments]
rerun-jed run output *ARGS:
	uv run python -m mtp.workflow.arms jed "$@"

[positional-arguments]
rerun-verify record:
	uv run python -m mtp.workflow.arms verify "$@"

# Execute complete cross-arm comparison consuming verified core, Noel, Jed, and Ilyas manifests.
[positional-arguments]
rerun-compare core_run noel_run jed_run ilyas_run output *ARGS:
	uv run python -m mtp.workflow.arms compare "$@"

[positional-arguments]
rerun-verify-comparison record:
	uv run python -m mtp.workflow.arms verify-comparison "$@"

# Additive BY/BH sidecars from a verified completed comparison; native tables unchanged.
[positional-arguments]
rerun-multiplicity comparison output *ARGS:
	uv run python -m mtp.multiplicity run --comparison-run "$1" --output "$2" "${@:3}"

[positional-arguments]
rerun-verify-multiplicity record:
	uv run python -m mtp.multiplicity verify --run "$1"

# Additive #771 lineage extension and complete six-fraction study LOSO.
[positional-arguments]
table1-influence-extension *ARGS:
	uv run python -m mtp.influence_extension "$@"

# Synthetic checks for the separately declared #771 sensitivity protocol.
test-influence-extension:
	uv run pytest tests/test_influence_extension.py tests/test_influence_extension_core.py tests/test_influence_extension_joint.py -v

# Isolated logical-guard probes; each invocation requires a fresh output directory (#969).
[positional-arguments]
probe-influence-extension-guards mutation output:
	uv run python scripts/adhoc/issue969_probe_influence_guards.py "$1" --output-dir "$2"

# Isolated original/current arm-parser regression evidence (#973).
[positional-arguments]
probe-arm-parser-guards *ARGS:
	uv run python scripts/adhoc/issue973_probe_parser_guards.py "$@"

sync:
	uv sync --all-extras --group dev --group qa

# Isolated copies prove actual plotted-label and metadata propagation (#955, #962).
[positional-arguments]
probe-direction-forwarding *ARGS:
	uv run python scripts/adhoc/issue955_probe_direction_forwarding.py "$@"

# Use the same declared uv version as the lock hook; run `just sync` first.
[positional-arguments]
lock *ARGS:
	uv run --frozen --group qa uv lock "$@"

# Print candidates by default; --baseline updates entries for explicit audit (#537).
[positional-arguments]
secrets-scan *ARGS:
	uv run --group qa detect-secrets scan "$@"

secrets-audit:
	uv run --group qa detect-secrets audit .secrets.baseline

lint:
	uv run ruff check src tests scripts
	uv run ruff format --check src tests scripts

format:
	uv run ruff check --fix src tests scripts
	uv run ruff format src tests scripts

trait-profile-investigation:
	uv run python scripts/regenerate_trait_profile_investigation.py

# Native pooled statistics only, with explicit zero-total eligibility (#551).
[positional-arguments]
trait-profile-pooled-stats *ARGS:
	uv run python scripts/create_trait_profile_infographic.py --stats-only --output-dir data/derived "$@"

# Regenerate six-trait pattern and full-table signed direction heatmaps.
six-shared-trait-heatmaps:
	uv run python scripts/create_six_shared_trait_cytoscape_heatmaps.py

# Signed increased-vs-decreased contrasts, as in the six-trait heatmaps, for a declared
# TRAIT_LIST on one run's PROFILES, profile METADATA (its total_trait_signal must match the
# profiles) and trait COLUMNS, into OUTPUT. Writes edge tables, the contrast table and
# trait_list_provenance.json only; no figures, and the six-trait directory is untouched (#1054).
[positional-arguments]
trait-list-direction-contrasts trait_list profiles metadata columns output *ARGS:
	uv run python scripts/create_six_shared_trait_sankeys.py --code-commit "$(git rev-parse HEAD)$(git diff --quiet HEAD -- {{code_paths}} && git diff --cached --quiet -- {{code_paths}} || echo +dirty)" --trait-list "$1" --profiles "$2" --profile-metadata "$3" --columns "$4" --output-dir "$5" "${@:6}"

# Check regenerated clustering and oral exports against a saved contrast table.
audit-six-trait-diagnostics baseline:
	uv run python scripts/adhoc/audit_clustering_stack.py --baseline-contrasts "{{baseline}}"

# Preserve final source identity separately from historical result provenance (#207).
# Fetch the named archive tags first; this recipe reads local git objects only.
audit-retired-composition-prs:
	mkdir -p .tmp/retired-composition-prs
	git archive --format=tar archive/pr-197-20261016 > .tmp/retired-composition-prs/pr197.tar
	git archive --format=tar archive/pr-526-20261016 > .tmp/retired-composition-prs/pr526.tar
	uv run python scripts/adhoc/audit_retired_composition_prs.py \
	  --pr197 .tmp/retired-composition-prs/pr197.tar \
	  --pr197-commit 865b15fcdc2be0bc9f0be976a37dbd0a2ed9f33d \
	  --pr526 .tmp/retired-composition-prs/pr526.tar \
	  --pr526-commit cf4679c5a8becc638f88d749386a40fa3599aa24 \
	  --output docs/reviews/2026-10-16-composition-archive-manifest.json

# Focused view of saved Jed contrasts; preserves their original Fisher/BH inference.
oral-six-trait-heatmap:
	source_commit="uncommitted" && \
	if git diff --quiet HEAD -- trait-visualization/trait-profile-investigation/six-shared-traits/six_shared_trait_condition_body_site_direction_contrasts.tsv trait-visualization/trait-profile-investigation/six-shared-traits/six_shared_trait_selection.tsv; then \
	  source_commit="$(git rev-parse HEAD)"; \
	fi && \
	uv run python scripts/create_six_shared_trait_cytoscape_heatmaps.py --oral-only --source-commit "$source_commit"

# PD fecal cellulose tables, figures and interpretations.md from the six-trait profile
# z-scores and the individual-condition Sankey edge tables that `just six-shared-trait-sankeys`
# writes. It stops when those tables no longer match the profile inputs (#687).
pd-cellulose-analysis:
	uv run python scripts/create_pd_cellulose_analysis.py

# Trace all six PD lollipop traits, each with its observed direction, to observed taxa and
# profile evidence, and record the PD fecal profile set and the hashes of the inputs read and
# the tables written, which the readers check (#584, #693, #703). After a data/derived change,
# run `just trait-heatmaps` and `just pd-cellulose-analysis` first (PD-cellulose/README.md).
pd-trait-taxa:
	uv run python scripts/lookup_pd_trait_taxa.py

# Plot thresholded species-only and rank-grouped species/genus/family PD influence contrasts.
pd-trait-taxa-heatmap:
	uv run python scripts/create_pd_six_trait_taxon_heatmap.py
	uv run python scripts/create_pd_exploratory_taxa_influence_heatmap.py

# The ledger reconciles data/derived under the primary phenotype strain rollup scope and
# records it in the manifest (#788).
# Build PD trait contributor tables from the Table 1 contribution ledger (#190) on the pinned
# KG in data/kg-microbe, then draw the dot heatmap and co-carriage network (#584).
pd-trait-contributors:
	uv run python scripts/create_pd_trait_contributors.py

# Test all eligible fecal condition-taxon and condition-trait pairs, apply global
# BH correction, and redraw the FDR-tiered evidence landscape.
pd-significant-ecological-network:
	uv run python scripts/create_pd_significant_condition_taxon_trait_network.py

# Build the ecological-reporter tables from the profile inputs (paired within-study
# contrasts of the six-trait z-scores, and the PD genus/species taxon panel), then draw
# the combined figure for the conditions with a descriptive 95% interval excluding zero
# for at least one of the six traits.
pd-ecological-reporter-combined:
	uv run python scripts/create_pd_ecological_reporter_combined.py

# Every six-shared-trait table and Sankey figure from the profile inputs, including
# six_shared_trait_summary.tsv and the individual-condition edge tables that
# `just pd-cellulose-analysis` reads. It removes the shared high-support node tables, which
# the two recipes below write again, so run them after it.
six-shared-trait-sankeys:
	uv run python scripts/create_six_shared_trait_sankeys.py

# Overlay increased (solid) and decreased (dashed) high-support Sankey links.
six-shared-trait-combined-sankey:
	uv run python scripts/create_six_shared_trait_sankeys.py --combined-high-support-only

# Redraw the saved increased, decreased and combined high-support selections.
six-shared-trait-high-support-sankeys:
	uv run python scripts/create_six_shared_trait_sankeys.py --high-support-layouts-only

trait-heatmaps:
	uv run python scripts/create_clustered_heatmap.py
	uv run python scripts/create_condition_taxon_heatmap.py
	uv run python scripts/create_taxa_guild_heatmap.py --direction-basis group_1_vs_group_0
	uv run python scripts/create_study_trait_profile_heatmap.py --direction-basis group_1_vs_group_0
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

# Canonical exports and legacy import-star contracts (#937); no scientific inputs.
test-shim-exports:
	uv run pytest tests/test_mtp_imports.py

[positional-arguments]
audit-shim-exports *ARGS:
	uv run python scripts/adhoc/issue937_audit_exports.py "$@"

# Additive policy and full-family BH/BY contract; no scientific rerun inputs.
test-multiple-testing:
	uv run pytest tests/test_multiple_testing.py

test-multiplicity-adapters:
	uv run pytest tests/test_multiplicity_adapters.py tests/test_multiplicity_basis_integration.py tests/test_multiplicity_regression_guards.py tests/test_multiplicity_cross_product.py tests/test_multiplicity_reason_vocabulary.py

# Analyzer, ranking exporter, and complete-universe density regression checks (#664).
test-stratification:
	uv run pytest tests/test_stratification.py tests/test_stratification_core.py tests/test_enrichment_pipeline.py tests/test_corpus_identity.py

test-cov:
	uv run pytest --cov=mtp --cov=bsdbng --cov-report=term-missing

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

# Check data/kg-microbe against the md5s pinned in methods/table1_taxon_influence_protocol.json,
# which the Table 1 influence audit and `just pd-trait-contributors` read. Refetch a changed
# file with `just download-kg-force`: `just download-kg` compares sizes only.
verify-pinned-kg:
	uv run python -m mtp.contribution_inputs

# Check that every identifier the guild, gram-stain and Table 1 figures group exists
# in data/kg-microbe/merged-kg_nodes.tsv with its expected label (run `just download-kg` first).
verify-trait-terms:
	uv run python scripts/verify_trait_terms.py

# Requires BugSigDB CSVs from `just download`; `just ingest` creates the data/studies/*.yaml inputs.
# Flatten per-study YAML into data/derived/*.tsv, including bugsigdb_taxa.tsv.
study-tsvs: ingest
	uv run python -m mtp.study_tsv

# Classify what every BugSigDB experiment compares (#188): disease vs healthy (with
# its orientation), disease vs disease, treatment/exposure, time point, other, or
# unresolved. Writes a versioned manifest and compares the six Table 1 directions and
# 15 log-ratios on the oriented disease-vs-healthy subset with the historical corpus,
# to data/derived/contrast_orientation/.
[positional-arguments]
contrast-orientation *ARGS:
	uv run python -m mtp.contrast_orientation --code-commit "$(git rev-parse HEAD)$(git diff --quiet HEAD -- {{code_paths}} && git diff --cached --quiet -- {{code_paths}} || echo +dirty)" "$@"

# Redraw the stratified hand-review sample behind methods/contrast_orientation_review.tsv.
# Writes an unfilled template to /dev/stdout's sibling path given by SAMPLE; the seed and
# the number per class are the only things that decide which experiments are drawn, so the
# same corpus and seed always give the same cases. Verdicts are then filled in by hand and
# the file committed; `just contrast-orientation` joins it back as validation_status.
contrast-review-sample SAMPLE="methods/contrast_orientation_review_sample.tsv" SEED="188" PER_CLASS="10":
	uv run python -m mtp.contrast_orientation --write-review-sample {{SAMPLE}} --review-seed {{SEED}} --review-per-class {{PER_CLASS}}

# Rebuild the data/analysis/cellulose_decr/ tables (methods/cellulose_decr.md). Stops when
# an input's hash differs from the committed provenance or the contrast manifest came from
# another classifier version; --input-root points at a checkout holding the pinned inputs,
# and --allow-input-drift regenerates every table from changed ones (#562). The pins are
# the v6 manifest and the phenotype-scope inputs, adopted with --repin in the final rerun
# at 511c891. A checkout with another classifier version refuses until it reruns
# `just contrast-orientation` and then --repin; the tables reproduce at the provenance's
# code_commit.
[positional-arguments]
cellulose-decr *ARGS:
	uv run python -m mtp.cellulose_decr.analysis --code-commit "$(git rev-parse HEAD)$(git diff --quiet HEAD -- {{code_paths}} && git diff --cached --quiet -- {{code_paths}} || echo +dirty)" "$@"

# Draw the cellulose_decr figures from its tables. The figures are regenerated, not
# committed (CONTRIBUTING.md, "Generated outputs on main"; #560). Their provenance records
# this checkout's commit, with +dirty when its code or tables are modified; the deck links
# the figures' code and tables at that commit (#556). That commit describes this
# checkout's tables only, so --analysis-dir is refused here (#609).
[positional-arguments]
cellulose-decr-figures *ARGS:
	#!/usr/bin/env bash
	set -euo pipefail
	for arg in "$@"; do
		case "$arg" in
			--analysis-dir* | --code-commit*)
				echo "the recipe records this checkout and its commit; drop $arg" >&2
				exit 1
				;;
		esac
	done
	commit="$(git rev-parse HEAD)$(git diff --quiet HEAD -- {{code_paths}} data/analysis/cellulose_decr && git diff --cached --quiet -- {{code_paths}} data/analysis/cellulose_decr || echo +dirty)"
	uv run python -m mtp.cellulose_decr.figures --code-commit "$commit" "$@"

# The #835 selection scan (methods/cellulose_decr_selection_scan.md): cellulose_decr's
# overlap and U6 null for every triple of eligible conditions and the top decreasing
# traits, on a verified Tier 1+2 core RUN (its run.json), into a new or empty OUTPUT
# directory. Pass --plan to print the counts and a timed estimate without running the null.
# Protocol v2 (#1050) makes phenotype primary and off the declared sensitivity.
# Scan a core of each scope, then report both with cellulose-decr-selection-scan-scopes.
[positional-arguments]
cellulose-decr-selection-scan run output *ARGS:
	uv run python -m mtp.cellulose_decr.selection_scan --code-commit "$(git rev-parse HEAD)$(git diff --quiet HEAD -- {{code_paths}} && git diff --cached --quiet -- {{code_paths}} || echo +dirty)" --core-run "$1" --output-dir "$2" "${@:3}"

# The selection scan's two strain-rollup scopes side by side (#1050): give the output
# directories of `just cellulose-decr-selection-scan` on a Tier 1+2 x phenotype core and
# on a Tier 1+2 x off core, in either order, and a new or empty OUTPUT directory.
# Both scans are verified and written together, so neither scope is reported alone.
[positional-arguments]
cellulose-decr-selection-scan-scopes scan_a scan_b output:
	uv run python -m mtp.cellulose_decr.selection_scan_scopes --code-commit "$(git rev-parse HEAD)$(git diff --quiet HEAD -- {{code_paths}} && git diff --cached --quiet -- {{code_paths}} || echo +dirty)" --scan "$1" --scan "$2" --output-dir "$3"

# The two methods plates and seven Tufte plates of the selected-figure deck, drawn
# from the rerun tables into this checkout's data/derived/tufte_plates/ (#573). Plates IV
# and V also quote the cellulose_decr results, so run `just cellulose-decr-figures` first.
# Pass --input-root to read another checkout's gitignored tables, and --overwrite to
# replace existing plates; they have no other copy, so back them up first. The plates
# record this checkout's commit (+dirty when its code is modified), at which the deck
# links their code (#613), and each plate's inputs with their hashes, which the deck
# checks (#637).
[positional-arguments]
figure-review-plates *ARGS:
	#!/usr/bin/env bash
	set -euo pipefail
	for arg in "$@"; do
		case "$arg" in
			--code-commit*)
				echo "the recipe records this checkout's commit; drop $arg" >&2
				exit 1
				;;
		esac
	done
	commit="$(git rev-parse HEAD)$(git diff --quiet HEAD -- {{code_paths}} && git diff --cached --quiet -- {{code_paths}} || echo +dirty)"
	uv run python scripts/generate_tufte_plates.py --code-commit "$commit" "$@"

deck_slides := "data/derived/figure_review_2026-09-21/slides"
deck_provenance := "data/derived/figure_review_2026-09-21/provenance"
deck_schema := "data/derived/figure_review_2026-09-21/provenance/figure_provenance_schema.yaml"

# The selected-figure deck: the HTML, slide manifest and assets in {{deck_slides}}, the
# figure provenance YAML in {{deck_provenance}} (validated against its LinkML schema),
# the PDF, and deck_build_status.json recording that validation and printing succeeded.
# Many inputs are gitignored local outputs (earlier deck manifests, data/derived figures),
# so a fresh clone cannot build it. The build refuses to replace an existing deck without
# --overwrite (back up slides/ and provenance/ first) (#573, #592). The slide text reads
# its numbers from the rerun's tables (#579). The plates and cellulose_decr figures come
# from this checkout, where their recipes draw them; the rest from --input-root. Slides
# whose code was never merged link it at its pull request (#614). The build checks every
# committed file a slide links or quotes against the revision the slide cites, from
# `git ls-tree` listings written here, and refuses when one differs; --allow-unverified
# builds a marked review deck instead (#637). Choose other folders with the variables, e.g.
# `just deck_slides=/tmp/s deck_provenance=/tmp/p figure-review-deck`, so that every
# step uses them; --output-dir and --provenance-dir are refused here. A worktree has no
# data/derived, so there also set deck_schema= to the schema beside the input root's YAML.
[positional-arguments]
figure-review-deck *ARGS:
	#!/usr/bin/env bash
	set -euo pipefail
	overwrite=0
	for arg in "$@"; do
		case "$arg" in
			-h | --help)
				uv run python scripts/build_slide_deck.py --help
				exit 0
				;;
			--output-dir* | --provenance-dir*)
				echo "set deck_slides= and deck_provenance= instead of $arg" >&2
				exit 1
				;;
			--git-trees* | --print-refs)
				echo "the recipe writes the git listings itself; drop $arg" >&2
				exit 1
				;;
			--overwrite)
				overwrite=1
				;;
		esac
	done
	if [ ! -f "{{deck_schema}}" ]; then
		echo "no LinkML schema at {{deck_schema}}; set deck_schema=" >&2
		exit 1
	fi
	slides="{{deck_slides}}"
	provenance="{{deck_provenance}}"
	# The status describes the deck it was written for and is rewritten only if every step
	# passes. With --overwrite it is removed first, so a build that fails partway cannot
	# leave it beside a half-rebuilt deck; a refusal after that leaves no status at all,
	# and only a build refused without --overwrite keeps it.
	if [ "$overwrite" = 1 ]; then
		rm -f "$slides/deck_build_status.json"
	fi
	# The committed files at each revision the slides cite, for the build to check (#637).
	# The unmerged pull requests' commits are fetched when this clone lacks them (#614); a
	# fetch slower than 1 kB/s for 30 s gives up, and the build then refuses that ref (#686).
	trees="$(mktemp -d)"
	trap 'rm -rf "$trees"' EXIT
	for ref in $(uv run python scripts/build_slide_deck.py --print-refs); do
		git cat-file -e "$ref^{commit}" 2>/dev/null || git -c http.lowSpeedLimit=1000 -c http.lowSpeedTime=30 fetch --quiet origin "$ref" || true
		git ls-tree -r -z --full-tree "$ref" > "$trees/$ref.tree" || rm -f "$trees/$ref.tree"
	done
	uv run python scripts/build_slide_deck.py --output-dir "$slides" --provenance-dir "$provenance" --git-trees "$trees" "$@"
	rm -f "$slides/deck_build_status.json"
	uv run linkml-validate -s "{{deck_schema}}" "$provenance/figure_provenance_and_results.yaml"
	# A review build prints under its own name and never replaces the deck's PDFs (#637).
	if grep -q '"review_build": true' "$slides/deck_validation_v6.json"; then
		review=true
		pdf="microbial-trait-profiling_v6_review.pdf"
		just html-to-pdf "$slides/selected_figures_v6.html" "$slides/$pdf"
	else
		review=false
		pdf="microbial-trait-profiling_v6.pdf"
		just html-to-pdf "$slides/selected_figures_v6.html" "$slides/$pdf"
		cp "$slides/$pdf" "$slides/microbial-trait-profiling.pdf"
		cp "$slides/$pdf" "$slides/selected_figures.pdf"
	fi
	pdf_bytes="$(wc -c < "$slides/$pdf" | tr -d ' ')"
	printf '{"provenance_linkml_validated": true, "review_build": %s, "pdf": "%s", "pdf_bytes": %s, "finished_utc": "%s"}\n' "$review" "$pdf" "$pdf_bytes" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$slides/deck_build_status.json"

# Print an HTML page to PDF with headless Chrome: $CHROME if set (it must be executable),
# else Playwright's Chromium, Google Chrome, or a chromium / google-chrome on PATH.
[positional-arguments]
html-to-pdf HTML PDF:
	#!/usr/bin/env bash
	set -euo pipefail
	html="$1"
	pdf="$2"
	if [ ! -f "$html" ]; then
		echo "no such page: $html" >&2
		exit 1
	fi
	if [ -n "${CHROME:-}" ]; then
		if [ ! -x "$CHROME" ]; then
			echo "CHROME is not an executable: $CHROME" >&2
			exit 1
		fi
		chrome="$CHROME"
	else
		shopt -s nullglob
		candidates=(
			"$HOME"/Library/Caches/ms-playwright/chromium*/chrome-mac*/Chromium.app/Contents/MacOS/Chromium
			"$HOME"/Library/Caches/ms-playwright/chromium_headless_shell*/chrome-headless-shell-mac*/chrome-headless-shell
			"$HOME"/.cache/ms-playwright/chromium*/chrome-linux*/chrome
			"$HOME"/.cache/ms-playwright/chromium_headless_shell*/chrome-headless-shell-linux*/chrome-headless-shell
			"/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
		)
		for name in chromium chromium-browser google-chrome google-chrome-stable; do
			if command -v "$name" > /dev/null; then
				candidates+=("$(command -v "$name")")
			fi
		done
		chrome=""
		for candidate in "${candidates[@]}"; do
			if [ -x "$candidate" ]; then
				chrome="$candidate"
				break
			fi
		done
		if [ -z "$chrome" ]; then
			echo "no headless Chrome found; set CHROME to its binary" >&2
			exit 1
		fi
	fi
	page="$(python3 -c 'import pathlib, sys; print(pathlib.Path(sys.argv[1]).resolve().as_uri())' "$html")"
	rm -f "$pdf"
	"$chrome" --headless --disable-gpu --no-pdf-header-footer --print-to-pdf-no-header --print-to-pdf="$pdf" "$page"
	test -s "$pdf"

#   - data/derived/bugsigdb_taxa.tsv (columns: direction, taxon_id, taxon_name,
#     taxonomic_rank), generated by `just study-tsvs` from per-study YAML.
#   - data/kg-microbe/merged-kg_{edges,nodes}.tsv (run `just download-kg`).
# Its BacDive strain rollup scope is the primary phenotype, as in the core workflow (#788).
# Build data/derived/taxa_trait_matrix.tsv from BugSigDB taxa and KG-Microbe traits.
taxa-trait-matrix: download-kg study-tsvs
	uv run python scripts/create_taxa_trait_matrix.py

# Takes explicit inputs: --reported-taxa, --edges, --nodes and --output-dir, and optionally
# a contrast-selected --study-dir for the signature taxa that gain a trait label.
# Compare trait coverage across the BacDive strain rollup scopes off, phenotype and full (#788).
[positional-arguments]
audit-strain-rollup *ARGS:
	uv run python scripts/adhoc/issue774_strain_rollup_coverage.py "$@"

# Sidecar carrier/source audit; explicit matching inputs, no primary output changes (#773).
[positional-arguments]
audit-label-support *ARGS:
	uv run python -m mtp.label_support_audit "$@"

# dest is a durable data directory, not a temporary one. The request is anonymous, and a file
# whose size or sha256 differs from the pin is refused. Pass the file to
# `just audit-label-support --madin-csv`.
# Download the Madin et al. condensed trait table pinned in mtp.madin_source (#1048).
download-madin dest:
	uv run python -m mtp.madin_source --dest-dir "{{dest}}"

# The analyses read human hosts only, without the trait columns that are zero in every
# human profile, as committed in data/derived/study_trait_profiles_wide.tsv.gz (#656).
# They also read 16S and shotgun (WMS) experiments only, the values of
# mtp.human_host_corpus.CORPUS_SEQUENCING_TYPES (#742).
# Build the canonical study trait profile table (human hosts, 16S and WMS experiments,
# non-empty trait columns).
study-trait-profiles: taxa-trait-matrix
	uv run python scripts/create_study_trait_profiles.py --host-species "Homo sapiens" --sequencing-types 16S WMS --drop-empty-traits

# Run after rebuilding data/derived: canonical human 16S/WMS coverage, README and hashes.
# Explicit paths are supported; this recipe never rebuilds upstream tables (#820/#821).
# Declare the upstream scope explicitly, e.g. just taxa-coverage --strain-rollup phenotype.
[positional-arguments]
taxa-coverage *ARGS:
	uv run python scripts/create_taxa_coverage_histogram.py "$@"

# After an intended input change, such as a BugSigDB release or a new KG pin, refresh the
# committed copy with `just study-trait-profiles-refresh-gz` and commit it (#667).
# Check that the profile table is the committed canonical one (#656).
study-trait-profiles-check:
	#!/usr/bin/env bash
	set -euo pipefail
	table=data/derived/study_trait_profiles_wide.tsv
	committed="$table.gz"
	if [ ! -s "$committed" ] || ! gzip -t "$committed" 2>/dev/null; then
		echo "$committed is missing, empty or not gzip; restore it: git checkout -- $committed (#667)" >&2
		exit 1
	fi
	if [ ! -s "$table" ]; then
		echo "$table is missing or empty; run just study-trait-profiles, or gzip -dc $committed > $table" >&2
		exit 1
	fi
	if ! gzip -dc "$committed" | cmp -s - "$table"; then
		echo "$table differs from the committed $committed (#656)." >&2
		echo "  Unintended: run just study-trait-profiles, or gzip -dc $committed > $table" >&2
		echo "  After an intended input change: just study-trait-profiles-refresh-gz, then commit $committed" >&2
		exit 1
	fi

# Replace the committed canonical profile table with the current one, then commit it (#667).
study-trait-profiles-refresh-gz:
	gzip -9 -n -c data/derived/study_trait_profiles_wide.tsv > data/derived/study_trait_profiles_wide.tsv.gz

# The chain builds data/derived under the primary phenotype strain rollup scope, so the
# ledger reconciles that scope (#788).
# Source-count reconciliation and rule-declared genus/family influence (#190).
# The default entry point declares the complete clean-clone input chain.
table1-taxon-influence: study-trait-profiles
	just table1-taxon-influence-from "data/kg-microbe" "data/derived" "data/derived/table1_taxon_influence" phenotype

# strain_rollup must be the scope the inputs' matrix was built with: off, phenotype or full.
# A core's derived products carry the scope its run recorded (#788).
# Use already-materialized, independently pinned inputs without rebuilding the defaults.
table1-taxon-influence-from kg_dir inputs="data/derived" output="data/derived/table1_taxon_influence" strain_rollup="phenotype":
	uv run python -m mtp.taxon_contributions --kg-dir "{{kg_dir}}" --inputs "{{inputs}}" --output "{{output}}" --strain-rollup "{{strain_rollup}}"
	just table1-native-influence "{{inputs}}" "{{output}}"

# Requires an unchanged, completed contribution ledger; no KG reload is needed.
table1-native-influence inputs="data/derived" output="data/derived/table1_taxon_influence":
	dirty_flag="" && \
	if ! (git diff --quiet HEAD -- {{code_paths}} && git diff --cached --quiet -- {{code_paths}}); then dirty_flag="--code-dirty"; fi && \
	uv run python -m mtp.native_taxon_influence --inputs "{{inputs}}" --output "{{output}}" --code-commit "$(git rev-parse HEAD)" $dirty_flag
	just table1-native-summary "{{output}}"

table1-native-summary output="data/derived/table1_taxon_influence":
	if ! (git diff --quiet HEAD -- {{code_paths}} && git diff --cached --quiet -- {{code_paths}}); then echo "Refusing to publish from a dirty tree" >&2; exit 1; fi
	uv run python -m mtp.native_influence_publication --source "{{output}}" --destination "{{output}}/review"

study-trait-profiles-ora:
	uv run python scripts/create_study_trait_profiles.py --binary --include-taxa-count --host-species "Homo sapiens" --sequencing-types 16S WMS

# Trait enrichment analysis (ORA) with hypergeometric test and Benjamini-Hochberg FDR.
analyze: study-trait-profiles-ora
	just timed analyze "uv run python scripts/trait_enrichment.py \
	  --input data/derived/study_trait_profiles_ora_wide.tsv \
	  --taxa-matrix data/derived/taxa_trait_matrix.tsv \
	  --output data/derived/trait_enrichment_report.txt \
	  --tsv-output data/derived/trait_enrichment_results.tsv \
	  --all-results-output data/derived/trait_enrichment_all_results.tsv \
	  --p-threshold {{p_threshold}} \
	  --min-fold {{min_fold}} \
	  --fdr-scope global"
	just timed export-rankings "uv run python scripts/export_trait_enrichment_rankings.py \
	  --input data/derived/trait_enrichment_all_results.tsv \
	  --all-results \
	  --fdr-scope global \
	  --overall-output data/derived/trait_enrichment_traits_overall_ranked.tsv \
	  --pair-output data/derived/trait_enrichment_trait_study_pairs_ranked.tsv"

# Hypergeometric depletion analysis.
analyze-depletion: study-trait-profiles-ora
	just timed analyze-depletion "uv run python scripts/trait_depletion.py \
	  --input data/derived/study_trait_profiles_ora_wide.tsv \
	  --taxa-matrix data/derived/taxa_trait_matrix.tsv \
	  --output data/derived/trait_depletion_report.txt \
	  --tsv-output data/derived/trait_depletion_results.tsv \
	  --all-results-output data/derived/trait_depletion_all_results.tsv \
	  --p-threshold {{p_threshold}} \
	  --max-fold 0.5 \
	  --fdr-scope global"
	just timed export-depletion-rankings "uv run python scripts/export_trait_enrichment_rankings.py \
	  --input data/derived/trait_depletion_all_results.tsv \
	  --all-results \
	  --fdr-scope global \
	  --overall-output data/derived/trait_depletion_traits_overall_ranked.tsv \
	  --pair-output data/derived/trait_depletion_trait_study_pairs_ranked.tsv"

# Deliver full rankings and export summaries for both enrichment and depletion (#254, #349).
export-rankings: analyze analyze-depletion
	@echo "Complete exports and ranking summaries generated."

# Run the ORA harness benchmark check suite over synthetic null and mixed scenarios (Refs #263, #343, #348, #351).
calibrate replicates="10" seed="42" scenario="all": study-trait-profiles-ora
	just calibrate-from data/derived/study_trait_profiles_ora_wide.tsv data/derived/taxa_trait_matrix.tsv data/derived/ora_calibration "{{replicates}}" "{{seed}}" "{{scenario}}"

# Calibrate existing, explicitly selected products without rebuilding another corpus (#551).
[positional-arguments]
calibrate-from profiles taxa output replicates="10" seed="42" scenario="all":
	uv run python -m mtp.ora_calibration.harness \
	  --profiles "$1" \
	  --taxa-matrix "$2" \
	  --output-dir "$3" \
	  --replicates "$4" \
	  --seed "$5" \
	  --scenario "$6" \
	  --alpha 0.05 \
	  --fdr-q 0.05 \
	  --q-thresholds 0.05 0.001 \
	  --signal-strength 0.6 \
	  --n-spike-traits 3 \
	  --test-type enrichment \
	  --min-fold 1.5 \
	  --max-fold 0.5

# Run paired composition arm calibration against the study-flip complete null (#263).
calibrate-composition replicates="100" seed="182":
	just timed calibrate-composition "uv run python -m mtp.composition_calibration.composition \
	  --profiles data/derived/study_trait_profiles_wide.tsv \
	  --columns data/derived/taxa_trait_matrix_columns.tsv \
	  --output-dir data/derived/calibration/composition \
	  --replicates {{replicates}} \
	  --seed {{seed}}"

# Run ORA background stratification and annotation-density sensitivity evaluation (#118).
# Run against explicit existing products without rebuilding source data (#664).
[positional-arguments]
analyze-strata-inputs *ARGS:
	uv run python scripts/run_stratified_enrichment.py "$@"

# Compare legacy selected denominators with complete-universe outputs (#664).
[positional-arguments]
audit-stratification-universe before after:
	uv run python scripts/adhoc/issue664_compare_stratification.py "$@"

analyze-strata: study-trait-profiles-ora
	just timed analyze-strata "uv run python scripts/run_stratified_enrichment.py \
	  --profiles data/derived/study_trait_profiles_ora_wide.tsv \
	  --experiments data/derived/experiments.tsv \
	  --taxa-matrix data/derived/taxa_trait_matrix.tsv \
	  --output-dir data/derived \
	  --p-threshold {{p_threshold}} \
	  --min-fold {{min_fold}}"

# Run site-stratified enrichment for an anatomical site (e.g. Feces, Saliva) (#118).
analyze-site SITE: study-trait-profiles-ora
	just timed analyze-site "uv run python scripts/run_stratified_enrichment.py \
	  --profiles data/derived/study_trait_profiles_ora_wide.tsv \
	  --experiments data/derived/experiments.tsv \
	  --taxa-matrix data/derived/taxa_trait_matrix.tsv \
	  --output-dir data/derived/sites/{{SITE}} \
	  --site {{SITE}} \
	  --p-threshold {{p_threshold}} \
	  --min-fold {{min_fold}}"

# Predeclared sensitivity scenarios for the 15 log-ratios (#183): zero exclusions,
# population-matched fractions, common-support references, fecal strata and
# leave-one-study-out influence. Needs data/derived/bugsigdb_taxa.tsv (`just study-tsvs`).
# Scenarios are fixed in methods/table1_log_ratio_sensitivity_protocol.json.
[positional-arguments]
table1-log-ratio-sensitivity *ARGS:
	uv run python -m mtp.log_ratio_sensitivity --code-commit "$(git rev-parse HEAD)$(git diff --quiet HEAD -- {{code_paths}} && git diff --cached --quiet -- {{code_paths}} || echo +dirty)" "$@"

# PNG export (guild-maps, table1-trait-maps) runs through kaleido 1.x, which
# drives a headless Chrome. This downloads (~170 MB) the version-pinned Chrome
# for Testing build into choreographer's data directory; kaleido prefers that
# copy over a system Chrome, so run it once for reproducible figure renders.
get-chrome:
	uv run kaleido_get_chrome

# Requires data/derived/study_trait_profiles_wide.tsv (`just study-trait-profiles`)
# and data/derived/experiments.tsv (`just study-tsvs`). Human samples only, matched
# inc/dec pairs only, MediaDive traits excluded. Writes the heatmap plus its value
# matrix, study-count matrix, and correlation table to data/derived/guild_maps/.
# Needs a local Chrome for the PNG (`just get-chrome`).
guild-maps:
	just timed guild-maps "uv run python scripts/create_guild_maps.py --code-commit \"$(git rev-parse HEAD)$(git diff --quiet HEAD -- {{code_paths}} && git diff --cached --quiet -- {{code_paths}} || echo +dirty)\""

# Same inputs and same corpus filters as guild-maps, but over the six individual
# traits of manuscript Table 1 instead of the two pooled oxygen guilds. Writes the
# combined condition x body-site heatmap, the paired directional-shift figure, and
# their backing tables to data/derived/table1_traits/.
# Needs a local Chrome for the PNG files (`just get-chrome`).
table1-trait-maps:
	just timed table1-trait-maps "uv run python scripts/create_table1_trait_maps.py --code-commit \"$(git rev-parse HEAD)$(git diff --quiet HEAD -- {{code_paths}} && git diff --cached --quiet -- {{code_paths}} || echo +dirty)\""

# Generate modern publication-ready figures inspired by the project slide deck (#347, #354).
# Renders Figures 1-5 with zero synthetic data and dynamic M universe calculations.
# Writes PNG figures to data/derived/deck_figures/.
# Needs a local Chrome for the PNGs (`just get-chrome`).
deck-figures: study-trait-profiles-ora
	just timed deck-figures "uv run python scripts/generate_deck_visualizations.py"

# Study-level paired log-ratios for all 15 Table 1 trait pairs; writes tables,
# coverage/provenance and a Markdown report. No figure renderer is needed.
[positional-arguments]
table1-compositional-stats *ARGS:
	uv run python -m mtp.trait_composition --code-commit "$(git rev-parse HEAD)$(git diff --quiet HEAD -- {{code_paths}} && git diff --cached --quiet -- {{code_paths}} || echo +dirty)" "$@"

# The study sign test over every retained non-MediaDive trait, not just the six
# of Table 1 (#255), with the significant set reported at both the nominal and
# the effective family size (#265). Candidates are declared in
# methods/trait_screen_manifest.json. No figure renderer is needed.
[positional-arguments]
trait-direction-screen *ARGS:
	just timed trait-direction-screen "uv run python -m mtp.trait_screen --code-commit \"$(git rev-parse HEAD)$(git diff --quiet HEAD -- {{code_paths}} && git diff --cached --quiet -- {{code_paths}} || echo +dirty)\" $*"

# Audit trait correlation, annotation density adjustment, and near-twin collinearity.
[positional-arguments]
audit-trait-density *ARGS:
	just timed audit-trait-density "uv run python scripts/audit_trait_density_correlation.py --code-commit \"$(git rev-parse HEAD)$(git diff --quiet HEAD -- {{code_paths}} && git diff --cached --quiet -- {{code_paths}} || echo +dirty)\" $*"

# Density-adjusted paired effects for a declared TRAIT_LIST (one identifier per line) on one
# run's PROFILES and COLUMNS, into OUTPUT; PC1 is still fit on every non-MediaDive trait and
# trait_list_provenance.json records the list's SHA-256 (#1054). Without a list: audit-trait-density.
[positional-arguments]
density-adjustment-trait-list trait_list profiles columns output *ARGS:
	uv run python -m mtp.density_adjustment --code-commit "$(git rev-parse HEAD)$(git diff --quiet HEAD -- {{code_paths}} && git diff --cached --quiet -- {{code_paths}} || echo +dirty)" --trait-list "$1" --profiles "$2" --columns "$3" --output-dir "$4" "${@:5}"

# Audit trait coupling sensitivity via LOTO, LCO, and three-arm evaluation (#260).
[positional-arguments]
audit-trait-coupling *ARGS:
	just timed audit-trait-coupling "uv run python scripts/audit_trait_coupling.py --code-commit \"$(git rev-parse HEAD)$(git diff --quiet HEAD -- {{code_paths}} && git diff --cached --quiet -- {{code_paths}} || echo +dirty)\" $*"

# Host marker evaluation can omit greenlet on arm64; see security-linux and CONTRIBUTING.md (#888).
# Audit source, dependencies, secrets, and host-selected locked packages.
security:
	uv run bandit -c pyproject.toml -r src scripts
	uv run deptry src
	uv run pre-commit run detect-secrets --all-files
	tmpdir="$(mktemp -d)" && \
	trap 'rm -rf "$tmpdir"' EXIT && \
	reqfile="$tmpdir/requirements.txt" && \
	uv export --quiet --frozen --all-extras --group dev --group qa --format requirements-txt \
	  --no-emit-project --no-emit-workspace \
	  --output-file "$reqfile" && \
	uv run pip-audit -r "$reqfile" --disable-pip \
	  --cache-dir "$tmpdir/pip-audit-cache"

# Audit locked packages selected for Linux x86_64 / Python 3.11, on any host.
security-linux:
	audit_tmp="$(mktemp -d)" && \
	trap 'rm -rf "$audit_tmp"' EXIT && \
	uv run --locked --group qa uv export --quiet --frozen --all-extras --group dev --group qa \
	  --format requirements-txt --no-emit-project --no-emit-workspace \
	  --output-file "$audit_tmp/requirements.txt" && \
	uv run --locked --group qa uv pip compile "$audit_tmp/requirements.txt" \
	  --python-platform x86_64-unknown-linux-gnu --python-version 3.11 --no-deps \
	  --no-header --no-annotate --output-file "$audit_tmp/requirements-linux.txt" && \
	uv run --locked --group qa pip-audit -r "$audit_tmp/requirements-linux.txt" \
	  --disable-pip --no-deps --cache-dir "$audit_tmp/pip-audit-cache"

precommit:
	uv run pre-commit run --all-files

validate-schema:
	uv run linkml validate {{schema}}
	uv run linkml validate {{legacy_schema}}

lint-schema:
	uv run linkml lint --config {{lintconfig}} --ignore-warnings {{schema}}
	uv run linkml lint --config {{lintconfig}} --ignore-warnings {{legacy_schema}}

check-schema: validate-schema lint-schema

sync-agent-docs:
	uv run python scripts/sync_agent_docs.py

check-agent-docs:
	uv run python scripts/sync_agent_docs.py
	git diff --exit-code .github/copilot-instructions.md

check: check-generated check-agent-docs lint typecheck test check-schema workflow-check workflow-subset-check security

# --- pipeline steps (timed, logged) ---

download:
	just timed download "uv run python -c 'from mtp.download import download_exports; download_exports()'"

ingest:
	just timed ingest "uv run python -c 'from mtp.ingest import ingest; ingest()'"

validate-studies:
	just timed validate-studies "uv run python -c 'from mtp.validate import validate_all_studies; validate_all_studies()'"

validate-studies-linkml:
	just timed validate-studies-linkml "uv run linkml validate -s {{schema}} -C Study data/studies/*.yaml"

# Inventory the measured evidence that could independently test the trait signals
# (#191): what each study's curated metadata mentions, which experiments are
# stratified by a measured marker, which inferred-function traits exist, and which
# external cohorts overlap the discovery corpus. Acquires no new data.
[positional-arguments]
validation-evidence *ARGS:
	uv run python -m mtp.validation_evidence "$@"

stats:
	@uv run python -c "from mtp.stats import print_stats; print_stats()"

# Empirical calibration of this repository's directional tests (#182): rejection
# rates when their own nulls are true, power under a declared injected effect, and
# the cost of selecting traits on the data. Scenarios, tolerances and the replicate
# budget are fixed in methods/calibration_protocol.json. About twenty-five minutes at
# the declared budget of 1400 null and 400 power replicates; `just calibration-report
# --replicates 5 --power-replicates 5` is the quick pass, and it reports every judged
# criterion as insufficient_budget because five replicates can certify nothing.
[positional-arguments]
calibration-report *ARGS:
	uv run python -m mtp.calibration --code-commit "$(git rev-parse HEAD)$(git diff --quiet HEAD -- {{code_paths}} && git diff --cached --quiet -- {{code_paths}} || echo +dirty)" "$@"

# Compare the calibration guild's eligible population with integrated production (#730).
[positional-arguments]
audit-calibration-eligibility profiles experiments output:
	uv run python scripts/adhoc/issue730_audit_calibration_eligibility.py --profiles "$1" --experiments "$2" --output "$3"

# Preserve scientific results while refreshing applicable assessment evidence (#729).
[positional-arguments]
audit-assessment-refresh *ARGS:
	uv run python scripts/adhoc/issue729_audit_assessment_refresh.py "$@"

# Just the shared flip-assignment artefact of #263: which studies each replicate of
# complete_null_study_flip exchanges, for another arm to rerun its own pipeline on
# the identical null. Seconds, not minutes.
[positional-arguments]
calibration-flip-assignments *ARGS:
	uv run python -m mtp.calibration --flip-assignments-only --code-commit "$(git rev-parse HEAD)$(git diff --quiet HEAD -- {{code_paths}} && git diff --cached --quiet -- {{code_paths}} || echo +dirty)" "$@"

pipeline: download ingest

# Additive binary-presence and equal-study depth sensitivity of a verified core (#772).
[positional-arguments]
presence-depth-sensitivity *ARGS:
	uv run python -m mtp.presence_depth_sensitivity "$@"

# The integrated evidence table (#178): one row per method result, each keeping its
# own family, units and states. Reads the other analyses' status files from
# data/derived/ when they exist; anything absent is recorded as pending, never as a
# negative result. Nothing is pooled across methods.
[positional-arguments]
integrated-evidence *ARGS:
	uv run python -m mtp.integrated_evidence --code-commit "$(git rev-parse HEAD)$(git diff --quiet HEAD -- {{code_paths}} && git diff --cached --quiet -- {{code_paths}} || echo +dirty)" "$@"
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

# The widened directional screen of #255 run under the directional arm's complete
# study-flip null (#263): how many traits each correction calls when nothing is there,
# beside the observed calls. Replicate i is replicate i of the shared flip artefact.
[positional-arguments]
calibration-screen *ARGS:
	uv run python -m mtp.screen_calibration --code-commit "$(git rev-parse HEAD)$(git diff --quiet HEAD -- {{code_paths}} && git diff --cached --quiet -- {{code_paths}} || echo +dirty)" "$@"

# Additive training/held-out study splits; native arm outputs remain unchanged (#854).
[positional-arguments]
selection-study-splits COMPARISON_RUN OUTPUT *ARGS:
	uv run python -m mtp.cellulose_decr.study_split_run --comparison-run "$1" --output-dir "$2" "${@:3}"

# Audit literal phenotype co-assertions for all strain scopes without rerunning arms (#846).
[positional-arguments]
audit-coassertions *ARGS:
	uv run python -m mtp.coassertion_audit "$@"
