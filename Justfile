set shell := ["bash", "-cu"]

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
	diff -u src/mtp/datamodel/core_products_pydantic.py "$tmpfile"

# Validate the executable specification and its graph without reading large inputs.
workflow-check spec="workflows/bugsigdb_core.yaml":
	uv run python -m mtp.workflow check "{{spec}}"
	uv run linkml validate -s schema/workflow.yaml -C WorkflowDefinition "{{spec}}"
	uv run linkml lint --config {{lintconfig}} --ignore-warnings schema/workflow.yaml
	uv run linkml lint --config {{lintconfig}} --ignore-warnings schema/core_products.yaml

workflow-describe:
	uv run python -m mtp.workflow describe workflows/bugsigdb_core.yaml > docs/workflow-products.md

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

core-ora run output:
	uv run python scripts/trait_enrichment.py --core-run "{{run}}" \
	  --output "{{output}}/report.txt" --tsv-output "{{output}}/results.tsv"

core-verify-consumer record:
	uv run python -m mtp.workflow verify-consumer "{{record}}"

sync:
	uv sync --all-extras --group dev --group qa

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

# Regenerate six-trait pattern and full-table signed direction heatmaps.
six-shared-trait-heatmaps:
	uv run python scripts/create_six_shared_trait_cytoscape_heatmaps.py

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

# Trace PD lollipop traits to observed taxa and profile evidence.
pd-trait-taxa:
	uv run python scripts/lookup_pd_trait_taxa.py

# Plot thresholded species-only and rank-grouped species/genus/family PD influence contrasts.
pd-trait-taxa-heatmap:
	uv run python scripts/create_pd_six_trait_taxon_heatmap.py
	uv run python scripts/create_pd_exploratory_taxa_influence_heatmap.py

# Refresh PD trait contributor dot heatmap with specificity and co-carriage checks.
pd-trait-contributors:
	uv run python scripts/create_pd_trait_contributors.py

# Test all eligible fecal condition-taxon and condition-trait pairs, apply global
# BH correction, and redraw the FDR-tiered evidence landscape.
pd-significant-ecological-network:
	uv run python scripts/create_pd_significant_condition_taxon_trait_network.py

# Redraw the combined ecological-reporter figure after retaining conditions with
# a descriptive 95% interval excluding zero for at least one of the six traits.
pd-ecological-reporter-combined:
	uv run python scripts/create_pd_ecological_reporter_combined.py

# Overlay increased (solid) and decreased (dashed) high-support Sankey links.
six-shared-trait-combined-sankey:
	uv run python scripts/create_six_shared_trait_sankeys.py --combined-high-support-only

# Redraw the saved increased, decreased and combined high-support selections.
six-shared-trait-high-support-sankeys:
	uv run python scripts/create_six_shared_trait_sankeys.py --high-support-layouts-only

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
# and --allow-input-drift regenerates every table from changed ones (#562). The pins are a
# v2 manifest and the metadata table before #555, so until #580 re-pins them a checkout
# with a newer classifier refuses; the tables reproduce at the provenance's code_commit.
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
# Many inputs are gitignored local outputs (earlier deck manifests, svg_sources/,
# data/derived figures), so a fresh clone cannot build it. The build refuses to replace
# an existing deck without --overwrite (back up slides/ and provenance/ first) and
# refuses contrast-orientation inputs other than the version its typed numbers describe
# (#573, #579, #592). The plates and cellulose_decr figures come from this checkout, where
# their recipes draw them; the rest from --input-root. The build checks every committed
# file a slide links or quotes against the revision the slide cites, from `git ls-tree`
# listings written here, and refuses when one differs; --allow-unverified builds a marked
# review deck instead (#637). Choose other folders with the variables, e.g.
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
	# The committed files at each revision the slides cite, for the build to check (#637)
	trees="$(mktemp -d)"
	trap 'rm -rf "$trees"' EXIT
	for ref in $(uv run python scripts/build_slide_deck.py --print-refs); do
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
# Build data/derived/taxa_trait_matrix.tsv from BugSigDB taxa and KG-Microbe traits.
taxa-trait-matrix: download-kg study-tsvs
	uv run python scripts/create_taxa_trait_matrix.py

# The canonical profile table the analyses read: human hosts only, without the trait
# columns that are zero in every human profile. It reproduces the committed
# data/derived/study_trait_profiles_wide.tsv.gz, which `just study-trait-profiles-check`
# compares it with (#656).
study-trait-profiles: taxa-trait-matrix
	uv run python scripts/create_study_trait_profiles.py --host-species "Homo sapiens" --drop-empty-traits

study-trait-profiles-check:
	gzip -dc data/derived/study_trait_profiles_wide.tsv.gz | cmp -s - data/derived/study_trait_profiles_wide.tsv || { echo "data/derived/study_trait_profiles_wide.tsv differs from the committed .gz; run just study-trait-profiles (#656)" >&2; exit 1; }

# Source-count reconciliation and rule-declared genus/family influence (#190).
# The default entry point declares the complete clean-clone input chain.
table1-taxon-influence: study-trait-profiles
	just table1-taxon-influence-from "data/kg-microbe" "data/derived" "data/derived/table1_taxon_influence"

# Use already-materialized, independently pinned inputs without rebuilding the defaults.
table1-taxon-influence-from kg_dir inputs="data/derived" output="data/derived/table1_taxon_influence":
	uv run python -m mtp.taxon_contributions --kg-dir "{{kg_dir}}" --inputs "{{inputs}}" --output "{{output}}"
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
	uv run python scripts/create_study_trait_profiles.py --binary --include-taxa-count --host-species "Homo sapiens"

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
	just timed calibrate "uv run python -m mtp.ora_calibration.harness \
	  --profiles data/derived/study_trait_profiles_ora_wide.tsv \
	  --taxa-matrix data/derived/taxa_trait_matrix.tsv \
	  --output-dir data/derived/ora_calibration \
	  --replicates {{replicates}} \
	  --seed {{seed}} \
	  --scenario {{scenario}} \
	  --alpha 0.05 \
	  --fdr-q 0.05 \
	  --signal-strength 0.6 \
	  --n-spike-traits 3 \
	  --test-type enrichment \
	  --min-fold 1.5 \
	  --max-fold 0.67"

# Run paired composition arm calibration against the study-flip complete null (#263).
calibrate-composition replicates="100" seed="182":
	just timed calibrate-composition "uv run python -m mtp.composition_calibration.composition \
	  --profiles data/derived/study_trait_profiles_wide.tsv \
	  --columns data/derived/taxa_trait_matrix_columns.tsv \
	  --output-dir data/derived/calibration/composition \
	  --replicates {{replicates}} \
	  --seed {{seed}}"

# Run ORA background stratification and annotation-density sensitivity evaluation (#118).
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

# Audit trait coupling sensitivity via LOTO, LCO, and three-arm evaluation (#260).
[positional-arguments]
audit-trait-coupling *ARGS:
	just timed audit-trait-coupling "uv run python scripts/audit_trait_coupling.py --code-commit \"$(git rev-parse HEAD)$(git diff --quiet HEAD -- {{code_paths}} && git diff --cached --quiet -- {{code_paths}} || echo +dirty)\" $*"

security:
	uv run bandit -c pyproject.toml -r src scripts
	uv run deptry src
	uv run pre-commit run detect-secrets --all-files
	tmpdir="$(mktemp -d)" && \
	trap 'rm -rf "$tmpdir"' EXIT && \
	reqfile="$tmpdir/requirements.txt" && \
	uv export --quiet --frozen --no-dev --format requirements-txt \
	  --no-emit-project --no-emit-workspace \
	  --output-file "$reqfile" && \
	uv run pip-audit -r "$reqfile" --disable-pip \
	  --cache-dir "$tmpdir/pip-audit-cache"

precommit:
	uv run pre-commit run --all-files

validate-schema:
	uv run linkml validate {{schema}}
	uv run linkml validate {{legacy_schema}}

lint-schema:
	uv run linkml lint --config {{lintconfig}} --validate --ignore-warnings {{schema}}
	uv run linkml lint --config {{lintconfig}} --ignore-warnings {{schema}}
	uv run linkml lint --config {{lintconfig}} --validate --ignore-warnings {{legacy_schema}}
	uv run linkml lint --config {{lintconfig}} --ignore-warnings {{legacy_schema}}

check-schema: validate-schema lint-schema

sync-agent-docs:
	uv run python scripts/sync_agent_docs.py

check-agent-docs:
	uv run python scripts/sync_agent_docs.py
	git diff --exit-code .github/copilot-instructions.md

check: check-generated check-agent-docs lint typecheck test check-schema workflow-check security

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

# Just the shared flip-assignment artefact of #263: which studies each replicate of
# complete_null_study_flip exchanges, for another arm to rerun its own pipeline on
# the identical null. Seconds, not minutes.
[positional-arguments]
calibration-flip-assignments *ARGS:
	uv run python -m mtp.calibration --flip-assignments-only --code-commit "$(git rev-parse HEAD)$(git diff --quiet HEAD -- {{code_paths}} && git diff --cached --quiet -- {{code_paths}} || echo +dirty)" "$@"

pipeline: download ingest

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
