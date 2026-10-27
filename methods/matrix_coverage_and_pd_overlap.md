# Matrix Trait Overlap, Tier 1+2 Coverage Cascade, and Parkinson's Stool Study Overlap

This note documents the definitions, data sources, and exact denominator semantics implemented in [`src/mtp/matrix_overlap_coverage.py`](../src/mtp/matrix_overlap_coverage.py) (`just matrix-overlap-coverage`, issue `#1053`). All interpretations follow [`methods/biological_interpretation.md`](biological_interpretation.md): KG-Microbe trait columns represent curated organism-level annotations rather than direct in vivo physiological measurements, and BugSigDB signatures record reported taxon abundance shifts across published comparisons.

---

## 1. Pairwise Trait Overlap in the Taxon × Trait Matrix (`matrix_trait_overlap.tsv`)

`mtp.matrix_overlap_coverage` requires an explicit newline-delimited trait list via `--trait-list` (parsed and validated by [`src/mtp/trait_lists.py`](../src/mtp/trait_lists.py)). For any declared pair of retained (non-MediaDive) trait columns $(a, b)$ in `taxa_trait_matrix.tsv` (with $N_{\text{matrix}} = 5,829$ taxon rows in the `tier12-phenotype` core run):

- Let $T_a = \{i : M_{i,a} > 0\}$ and $T_b = \{i : M_{i,b} > 0\}$ denote the sets of matrix taxa carrying a positive annotation count for traits $a$ and $b$.
- The pairwise Jaccard overlap is:
  $$J(a, b) = \frac{|T_a \cap T_b|}{|T_a \cup T_b|}$$

### Focal Pairs in `matrix_trait_overlap.tsv` (`tier12-phenotype`, `methods/trait_lists/headline_ten_phenotype.txt`, $N_{\text{matrix}} = 5,829$)

| Trait A | Trait B | $|T_a|$ | $|T_b|$ | $|T_a \cap T_b|$ | $|T_a \cup T_b|$ | Jaccard (%) |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| `METPO:1000605` (`facultatively anaerobic`) | `METPO:1000608` (`facultatively aerobic`) | 712 | 920 | 535 | 1,097 | 48.77% (48.8%) |
| `METPO:2000103\|GO:0030245` (`cellulose catabolic process`) | `METPO:2000103\|GO:0045493` (`xylan catabolic process`) | 189 | 210 | 109 | 290 | 37.59% (37.6%) |
| `METPO:1000603` (`anaerobic`) | `METPO:1000608` (`facultatively aerobic`) | 1,884 | 920 | 493 | 2,311 | 21.33% (21.3%) |

Full 45-pair output for the 10-trait list [`methods/trait_lists/headline_ten_phenotype.txt`](trait_lists/headline_ten_phenotype.txt) is written to `matrix_trait_overlap.tsv`.

---

## 2. Tier 1+2 Taxon Coverage by Rank (`tier12_coverage_by_rank.tsv`)

### Why `just taxa-coverage` (`scripts/create_taxa_coverage_histogram.py`) Computes a Different Quantity

[`scripts/create_taxa_coverage_histogram.py`](../scripts/create_taxa_coverage_histogram.py) (`just taxa-coverage`) computes a different upstream table set (`taxa_coverage_by_rank.tsv`, `taxa_kgmicrobe_coverage_details.tsv`):

1. **Corpus scope**: `select_corpus` filters `experiments.tsv`, `bugsigdb_taxa.tsv`, and `study_trait_profiles_wide.tsv` to the full human 16S/WMS population (`host_species == "Homo sapiens"` and `sequencing_type in {"16S", "WMS"}`, $1,047$ studies and $3,372$ experiments, with `contrast_filter: "none"`) and restricts `taxa_trait_matrix.tsv` to the reported taxa in that population, rather than restricting to the $1,870$ active Tier 1+2 paired signatures (`disease_vs_healthy_tier1_2_subject_controls`, $548$ studies, $935$ experiments, $3,221$ unique NCBI taxa).
2. **Trait universe**: `load_and_process_data` classifies primary coverage (`has_traits`, `coverage_category`) using `trait_count > 0` across all profile-vocabulary trait columns in `study_trait_profiles_wide.tsv` — including the `METPO:2000517|mediadive.medium:*` growth-medium columns when present — while recording `core_trait_count` separately, whereas `mtp.matrix_overlap_coverage` defines labeled/unlabeled taxa strictly over the $344$ retained non-MediaDive trait columns.
3. **Occurrence vs deduplicated-mention and study-rarity breakdown**: `scripts/create_taxa_coverage_histogram.py` tabulates raw signature occurrences (`signature_occurrences`) and rank/profile coverage across the unselected human 16S/WMS corpus, whereas `mtp.matrix_overlap_coverage` deduplicates `(study_id, experiment_id, inc_dec, taxon_id)` profile mentions and tabulates rank-stratified taxon/mention coverage plus distinct study-count rarity ($n_{\text{studies}} = 1$, $n_{\text{studies}} \le 2$, median $n_{\text{studies}}$) for labeled versus unlabeled taxa on the active Tier 1+2 corpus in `tier12_coverage_by_rank.tsv`.

### Exact Tier 1+2 Coverage Definitions (`mtp.matrix_overlap_coverage`)

In `mtp.matrix_overlap_coverage` (`tier12_coverage_by_rank.tsv` and `summary.json`):

- **Retained trait set ($344$ columns)**: All columns in `taxa_trait_matrix_columns.tsv` (`feature_metadata`) present in `taxa_trait_matrix.tsv` except `METPO:2000517|mediadive.medium:*` growth-medium terms. A taxon is **labeled** if and only if it has a positive value in at least one of those $344$ retained columns.
- **Active Tier 1+2 mentions ($19,388$ distinct `(profile, taxon_id)` mentions)**: `steps/export_tables/derived/bugsigdb_taxa.tsv` (`reported_taxa`) preserves repeated within-signature occurrences (`coverage_policy: preserve_repeated_occurrences`, $19,566$ raw rows across the $1,870$ active Tier 1+2 profiles). Deduplicating by `(study_id, experiment_id, inc_dec, taxon_id)` yields $19,388$ distinct `(profile, taxon_id)` mentions across $3,221$ unique NCBI taxa.
- **Key Tier 1+2 `tier12-phenotype` figures (`tier12_coverage_by_rank.tsv`)**:
  - **All ranks (`rank == "all"`, $3,221$ unique taxa, $19,388$ mentions)**: $2,363$ labeled taxa ($73.36\%$) and $858$ unlabeled taxa ($26.64\%$); $17,913$ labeled mentions ($92.39\%$) and $1,475$ unlabeled mentions ($7.61\%$).
  - **Genus rank (`rank == "genus"`, $842$ unique genera, $8,495$ mentions)**: $739$ labeled genera ($87.77\%$) and $103$ unlabeled genera ($12.23\%$); $8,324$ of $8,495$ genus mentions ($97.99\%$, i.e. $98.0\%$) carry at least one retained trait annotation ($171$ unlabeled mentions, $2.01\%$).
  - **Species rank (`rank == "species"`, $1,983$ unique species, $6,040$ mentions)**: $1,287$ labeled species ($64.90\%$) and $696$ unlabeled species ($35.10\%$, i.e. $35.1\%$); $4,858$ labeled mentions ($80.43\%$) and $1,182$ unlabeled mentions ($19.57\%$).
  - **Single-study rarity by labeling status (`rank == "all"`)**: Among the $858$ unlabeled taxa, $717$ ($83.57\%$, i.e. $84\%$) appear in only $1$ study and $805$ ($93.82\%$, i.e. $94\%$) appear in $\le 2$ studies (median $1$ study). Among the $2,363$ labeled taxa, $1,124$ ($47.57\%$, i.e. $48\%$) appear in $1$ study and $1,465$ ($62.00\%$, i.e. $62\%$) appear in $\le 2$ studies (median $2$ studies).

---

## 3. Two-Axis Filter Cascade (`filter_cascade.tsv`)

### Trait Axis

1. **KG-Microbe trait catalog (`kg_trait_catalog`)**: $2,677$ total trait columns in `feature_catalog.tsv` (`feature_catalog`) and `taxa_trait_matrix_columns.tsv` (`feature_metadata`), of which $2,612$ have positive counts in the Tier 1+2 `count_profiles` table (`study_trait_profiles_wide.tsv`) after `--drop-empty-traits` drops $65$ zero-sum columns on the Tier 1+2 subset.
2. **Retained non-MediaDive traits (`retained_non_mediadive_traits`)**: $344$ traits retained after dropping $2,333$ `METPO:2000517|mediadive.medium:*` growth-medium columns:
   - $162$ METPO phenotypes (`METPO:*` without compound predicate pipes)
   - $102$ carbon-source traits (`METPO:2000006|*`, comprising $101$ `CHEBI:*` substrates plus `METPO:2000006|carbon_substrates:N-acetylgulcosamine`)
   - $68$ pathway traits (`METPO:2000103|pathways:*`)
   - $12$ GO biological process traits (`METPO:2000103|GO:*`)

### Study and Signature Axis (Distinguishing the Upstream Populations in `filter_cascade.tsv`)

Across the workflow, four nested study/signature steps are recorded in `filter_cascade.tsv`:

- **All-hosts BugSigDB export (`all_hosts_bugsigdb_export`)**: $1,245$ studies, $4,435$ experiments, and $7,329$ non-empty directional signatures in `steps/source_export_tables/source/derived/experiments.tsv` (`source_experiments`), spanning all host species and sequencing types.
- **Canonical human 16S/WMS profile corpus (`full_human_16s_wms`)**: $1,047$ studies ($198$ dropped), $3,372$ experiments, and $5,629$ non-empty directional signatures ($1,700$ dropped) matching `Host species == "Homo sapiens"` and `Sequencing type in {"16S", "WMS"}`.
- **Paired human 16S/WMS corpus (`paired_human_16s_wms`)**: $932$ studies ($115$ dropped), $2,257$ paired experiments, and $4,514$ directional signatures ($1,115$ unpaired signatures dropped).
- **Tier 1+2 disease vs healthy subject controls (`tier12_subject_controls`)**: $548$ studies ($384$ dropped), $935$ paired comparisons, and $1,870$ directional signatures ($2,644$ dropped) across $211$ conditions, including $355$ stool studies and $595$ stool comparisons (`Body site == "Feces"`).

---

## 4. Parkinson's Disease Stool Study Taxon Overlap (`parkinsons_stool_study_overlap_summary.tsv`)

Experiments are matched when `Body site == "Feces"` (exact compartment match; multi-site strings such as `"Feces,Rectum"` are excluded) and `Condition` lists `"Parkinson's disease"` as an exact comma-separated condition token (`any(tok == "parkinson's disease" for tok in cond_tokens)`, which retains multi-condition curated rows such as `"Constipation,Parkinson's disease"` while excluding unrelated substrings such as `"Atypical parkinsonism"`). Within each scope, only the selected Parkinson's stool experiments of study $s$ contribute to its reported NCBI taxon set $U_s$. Across all $\binom{K}{2}$ unordered study pairs $(s_1, s_2)$:

$$J(s_1, s_2) = \frac{|U_{s_1} \cap U_{s_2}|}{|U_{s_1} \cup U_{s_2}|}$$

Four explicit scopes are emitted in `parkinsons_stool_study_overlap_summary.tsv` and `parkinsons_stool_study_pairs.tsv`:

- **Tier 1+2 paired (`tier12_paired`)**: $K = 36$ studies ($70$ comparisons, $\binom{36}{2} = 630$ pairs; $388$ nonzero, $242$ zero), mean pairwise taxon Jaccard $= 0.053302$ ($5.33\%$, i.e. $5.3\%$), median $= 0.038099$ ($3.81\%$, i.e. $3.8\%$).
- **Full export (`full_export`, comma-separated token match across all hosts in `source_experiments`)**: $K = 59$ studies ($157$ comparisons, $\binom{59}{2} = 1,711$ pairs; $947$ nonzero, $764$ zero), mean pairwise taxon Jaccard $= 0.041030$ ($4.10\%$, i.e. $4.1\%$), median $= 0.024390$ ($2.44\%$, i.e. $2.4\%$).
- **Human-host export (`export_human_hosts`, `Host species == "Homo sapiens"` in `source_experiments`)**: $K = 56$ studies ($149$ comparisons, $\binom{56}{2} = 1,540$ pairs; $860$ nonzero, $680$ zero), mean pairwise taxon Jaccard $= 0.041823$ ($4.18\%$), median $= 0.024390$ ($2.44\%$).
- **Full export exact condition (`full_export_exact_condition`, `Condition == "Parkinson's disease"` in `source_experiments`)**: $K = 58$ studies ($156$ comparisons, $\binom{58}{2} = 1,653$ pairs; $918$ nonzero, $735$ zero), mean pairwise taxon Jaccard $= 0.041188$ ($4.12\%$), median $= 0.024390$ ($2.44\%$).

---

## 5. Provenance and Core-Run Guards

`mtp.matrix_overlap_coverage` enforces four provenance invariants before writing any output files:

1. **Tier 1+2 population guard**: `--core-run` must be a verified core bundle with `contrast_population == "disease_vs_healthy_tier1_2_subject_controls"` and `direction_basis == "disease_vs_comparison"`. Passing an `all_paired` core run or a tampered manifest raises `ValueError`.
2. **Durable output root guard (`#1049`, `#1051`)**: Temporary output directories matched by [`src/mtp/workflow/output_roots.py`](../src/mtp/workflow/output_roots.py) (`/tmp`, `/private/tmp`, `/var/tmp`, `/private/var/tmp`, `/var/folders`, `/private/var/folders`, and `$TMPDIR` / `tempfile.gettempdir()`) are refused unless `--allow-temporary-output` is explicitly passed for ephemeral test runs.
3. **Clean output directory guard**: `--output-dir` must be a new or empty directory so stale tables from earlier runs cannot survive alongside newly emitted tables.
4. **Core-bundle non-overlap preflight and bound provenance**: `preflight_consumer_outputs` runs before `output_dir.mkdir(...)` so an output path inside the core run bundle is rejected without mutating the core run directory. `provenance.json` records the sidecar's own `code_commit`, `script_sha256` (`src/mtp/matrix_overlap_coverage.py` and `src/mtp/trait_lists.py`), `methods_sha256` (`methods/matrix_coverage_and_pd_overlap.md`), `core_git_commit_sha`, and run-relative `steps/...` paths for all consumed core products, and is bound alongside the six primary outputs inside `core_consumption.json`.
