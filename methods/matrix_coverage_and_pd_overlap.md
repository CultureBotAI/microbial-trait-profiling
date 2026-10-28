# Matrix Trait Overlap, Tier 1+2 Coverage Cascade, and Parkinson's Stool Study Overlap

This note documents the definitions, data sources, and exact denominator semantics implemented in [`src/mtp/matrix_overlap_coverage.py`](../src/mtp/matrix_overlap_coverage.py) (`just matrix-overlap-coverage`, issue `#1053`). All interpretations follow [`methods/biological_interpretation.md`](biological_interpretation.md): KG-Microbe trait columns represent curated organism-level annotations rather than direct in vivo physiological measurements, and BugSigDB signatures record reported taxon abundance shifts across published comparisons.

The numeric tables below are historical references from the F Tier 1+2
`phenotype` core, with separate raw-object trait keys and the then-current
ten-trait list. They are not results from a rerun with the canonical facultative
union. The current [declared list](trait_lists/headline_ten_phenotype.txt) contains
nine traits despite its historical filename and requests 36 unordered trait
pairs; the archived ten-trait list requested 45. Trait-key policy and strain
rollup are separate choices: `phenotype` alone does not specify the trait-key
policy.

---

## 1. Pairwise Trait Overlap in the Taxon × Trait Matrix (`matrix_trait_overlap.tsv`)

`mtp.matrix_overlap_coverage` requires an explicit newline-delimited trait list via `--trait-list` (parsed and validated by [`src/mtp/trait_lists.py`](../src/mtp/trait_lists.py)). For any declared pair of retained (non-MediaDive) trait columns $(a, b)$ in `taxa_trait_matrix.tsv`:

- Let $T_a = \{i : M_{i,a} > 0\}$ and $T_b = \{i : M_{i,b} > 0\}$ denote the sets of matrix taxa carrying a positive annotation count for traits $a$ and $b$.
- The pairwise Jaccard overlap is:
  $$J(a, b) = \frac{|T_a \cap T_b|}{|T_a \cup T_b|}$$

An empty union has undefined Jaccard, recorded as JSON `null` or a blank TSV
cell. A nonempty union with no shared taxa has Jaccard zero. Matrix cells must be
finite, nonnegative numbers; malformed values are refused rather than converted
to absence of annotation.

### Historical Focal Pairs (`tier12-phenotype`, raw-object ten-trait list, $N_{\text{matrix}} = 5,829$)

| Trait A | Trait B | $|T_a|$ | $|T_b|$ | $|T_a \cap T_b|$ | $|T_a \cup T_b|$ | Jaccard (%) |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| `METPO:1000605` (`facultatively anaerobic`) | `METPO:1000608` (`facultatively aerobic`) | 712 | 920 | 535 | 1,097 | 48.77% (48.8%) |
| `METPO:2000103\|GO:0030245` (`cellulose catabolic process`) | `METPO:2000103\|GO:0045493` (`xylan catabolic process`) | 189 | 210 | 109 | 290 | 37.59% (37.6%) |
| `METPO:1000603` (`anaerobic`) | `METPO:1000608` (`facultatively aerobic`) | 1,884 | 920 | 493 | 2,311 | 21.33% (21.3%) |

The historical `matrix_trait_overlap.tsv` contains all 45 pairs for the archived
ten-trait list. New output uses the explicitly supplied list and the verified
core's matrix; it does not reconstruct the historical ten-trait family from the
current nine-trait file.

---

## 2. Tier 1+2 Taxon Coverage by Rank (`tier12_coverage_by_rank.tsv`)

### Why `just taxa-coverage` (`scripts/create_taxa_coverage_histogram.py`) Computes a Different Quantity

[`scripts/create_taxa_coverage_histogram.py`](../scripts/create_taxa_coverage_histogram.py) (`just taxa-coverage`) computes a different upstream table set (`taxa_coverage_by_rank.tsv`, `taxa_kgmicrobe_coverage_details.tsv`):

1. **Corpus scope**: `select_corpus` filters `experiments.tsv`, `bugsigdb_taxa.tsv`, and `study_trait_profiles_wide.tsv` to the full human 16S/WMS population (`host_species == "Homo sapiens"` and `sequencing_type in {"16S", "WMS"}`, $1,047$ studies and $3,372$ experiments, with `contrast_filter: "none"`) and restricts `taxa_trait_matrix.tsv` to the reported taxa in that population, rather than restricting to the $1,870$ active Tier 1+2 paired signatures (`disease_vs_healthy_tier1_2_subject_controls`, $548$ studies, $935$ experiments, $3,221$ unique NCBI taxa).
2. **Trait universe**: `load_and_process_data` classifies primary coverage (`has_traits`, `coverage_category`) using `trait_count > 0` across all profile-vocabulary trait columns in `study_trait_profiles_wide.tsv` — including the `METPO:2000517|mediadive.medium:*` growth-medium columns when present — while recording `core_trait_count` separately, whereas `mtp.matrix_overlap_coverage` defines labeled/unlabeled taxa strictly over the $344$ retained non-MediaDive trait columns.
3. **Occurrence vs deduplicated-mention and study-rarity breakdown**: `scripts/create_taxa_coverage_histogram.py` tabulates raw signature occurrences (`signature_occurrences`) and rank/profile coverage across the unselected human 16S/WMS corpus, whereas `mtp.matrix_overlap_coverage` deduplicates `(study_id, experiment_id, inc_dec, taxon_id)` profile mentions and tabulates rank-stratified taxon/mention coverage plus distinct study-count rarity ($n_{\text{studies}} = 1$, $n_{\text{studies}} \le 2$, median $n_{\text{studies}}$) for labeled versus unlabeled taxa on the active Tier 1+2 corpus in `tier12_coverage_by_rank.tsv`.

### Exact Tier 1+2 Coverage Definitions (`mtp.matrix_overlap_coverage`)

In `mtp.matrix_overlap_coverage` (`tier12_coverage_by_rank.tsv` and `summary.json`):

- **Retained trait set (historically $344$ columns)**: Feature metadata identifies the matrix columns; every declared feature must be present. The shared [`feature_class`](../src/mtp/workflow/coverage.py) classifier excludes any feature ID containing `mediadive.medium` after lowercasing, rather than testing only one predicate prefix. A taxon is **labeled** if and only if it has a positive value in at least one retained column. The historical excluded terms used the `METPO:2000517|mediadive.medium:*` prefix.
- **Active Tier 1+2 mentions ($19,388$ distinct `(profile, taxon_id)` mentions)**: `steps/export_tables/derived/bugsigdb_taxa.tsv` (`reported_taxa`) preserves repeated within-signature occurrences (`coverage_policy: preserve_repeated_occurrences`, $19,566$ raw rows across the $1,870$ active Tier 1+2 profiles). Deduplicating by `(study_id, experiment_id, inc_dec, taxon_id)` yields $19,388$ distinct `(profile, taxon_id)` mentions across $3,221$ unique NCBI taxa.
- **Key Tier 1+2 `tier12-phenotype` figures (`tier12_coverage_by_rank.tsv`)**:
  - **All ranks (`rank == "all"`, $3,221$ unique taxa, $19,388$ mentions)**: $2,363$ labeled taxa ($73.36\%$) and $858$ unlabeled taxa ($26.64\%$); $17,913$ labeled mentions ($92.39\%$) and $1,475$ unlabeled mentions ($7.61\%$).
  - **Genus rank (`rank == "genus"`, $842$ unique genera, $8,495$ mentions)**: $739$ labeled genera ($87.77\%$) and $103$ unlabeled genera ($12.23\%$); $8,324$ of $8,495$ genus mentions ($97.99\%$, i.e. $98.0\%$) carry at least one retained trait annotation ($171$ unlabeled mentions, $2.01\%$).
  - **Species rank (`rank == "species"`, $1,983$ unique species, $6,040$ mentions)**: $1,287$ labeled species ($64.90\%$) and $696$ unlabeled species ($35.10\%$, i.e. $35.1\%$); $4,858$ labeled mentions ($80.43\%$) and $1,182$ unlabeled mentions ($19.57\%$).
  - **Single-study rarity by labeling status (`rank == "all"`)**: Among the $858$ unlabeled taxa, $717$ ($83.57\%$, i.e. $84\%$) appear in only $1$ study and $805$ ($93.82\%$, i.e. $94\%$) appear in $\le 2$ studies (median $1$ study). Among the $2,363$ labeled taxa, $1,124$ ($47.57\%$, i.e. $48\%$) appear in $1$ study and $1,465$ ($62.00\%$, i.e. $62\%$) appear in $\le 2$ studies (median $2$ studies).

Profile membership uses the complete `(study_id, experiment_id, inc_dec)` key;
experiment membership uses `(study_id, experiment_id)`. Valid reported-taxon
rows outside the active profile population are excluded and counted before
within-profile deduplication. The coverage summary records
`tier12_input_occurrence_rows`, `tier12_excluded_occurrence_rows` and
`tier12_excluded_profile_keys`; `tier12_raw_occurrence_rows` counts the selected
occurrences before deduplication. These population exclusions are distinct from
malformed input: blank or duplicate identity keys, unknown directions, missing
active-experiment metadata and reported taxa with no matching experiment
metadata are refused. A repeated reported-taxon occurrence is allowed and is
deduplicated only for the stated mention statistic.

Ratios with zero denominators and medians of empty groups are undefined, emitted
as JSON `null` and blank TSV cells. Counts still describe the empty group; an
undefined summary is not reported as a biological zero.

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

The human profile population uses the shared
[`experiment_in_profile_population`](../src/mtp/human_host_corpus.py) predicate:
host `Homo sapiens` and sequencing type `16S` or `WMS`, with surrounding
whitespace stripped and case preserved. Source exports are required inputs;
selected Tier 1+2 tables cannot substitute for the full upstream exports. Stool
counts use the literal raw body-site value `Feces`. Case variants, surrounding
whitespace, aliases and multi-site strings do not qualify. No compartment
normalization is applied.

---

## 4. Parkinson's Disease Stool Study Taxon Overlap (`parkinsons_stool_study_overlap_summary.tsv`)

Experiments are matched when the raw `Body site == "Feces"` (literal string equality, with no stripping or case folding; multi-site strings such as `"Feces,Rectum"` are excluded) and `Condition` lists `"Parkinson's disease"` as an exact comma-separated condition token after stripping and lowercasing each token (`any(tok == "parkinson's disease" for tok in cond_tokens)`, which retains multi-condition curated rows such as `"Constipation,Parkinson's disease"` while excluding unrelated substrings such as `"Atypical parkinsonism"`). Within each scope, only the selected Parkinson's stool experiments of study $s$ contribute to its reported NCBI taxon set $U_s$. Across all $\binom{K}{2}$ unordered study pairs $(s_1, s_2)$:

$$J(s_1, s_2) = \frac{|U_{s_1} \cap U_{s_2}|}{|U_{s_1} \cup U_{s_2}|}$$

Four explicit scopes are emitted in `parkinsons_stool_study_overlap_summary.tsv` and `parkinsons_stool_study_pairs.tsv`:

The following values are historical F references. For `tier12_paired`, taxon
rows must additionally belong to the active profile keys. Its summary records
`input_occurrence_rows`, `excluded_outside_active_profile_rows` and
`excluded_outside_active_profile_keys` before disease/site selection. The human
export scope applies the shared host predicate with no sequencing restriction.
When fewer than two studies provide taxon sets, no pairwise mean, median,
minimum or maximum is estimated; these fields are JSON `null` or blank TSV
cells. Pair and study counts remain explicit. The overlaps are descriptive;
study pairs share studies and are not independent biological replicates.

- **Tier 1+2 paired (`tier12_paired`)**: $K = 36$ studies ($70$ comparisons, $\binom{36}{2} = 630$ pairs; $427$ nonzero, $203$ zero), mean pairwise taxon Jaccard $= 0.053302$ ($5.33\%$, i.e. $5.3\%$), median $= 0.038099$ ($3.81\%$, i.e. $3.8\%$).
- **Full export (`full_export`, comma-separated token match across all hosts in `source_experiments`)**: $K = 59$ studies ($157$ comparisons, $\binom{59}{2} = 1,711$ pairs; $1,004$ nonzero, $707$ zero), mean pairwise taxon Jaccard $= 0.041030$ ($4.10\%$, i.e. $4.1\%$), median $= 0.024390$ ($2.44\%$, i.e. $2.4\%$).
- **Human-host export (`export_human_hosts`, `Host species == "Homo sapiens"` in `source_experiments`)**: $K = 56$ studies ($149$ comparisons, $\binom{56}{2} = 1,540$ pairs; $901$ nonzero, $639$ zero), mean pairwise taxon Jaccard $= 0.041823$ ($4.18\%$), median $= 0.024390$ ($2.44\%$).
- **Full export exact condition (`full_export_exact_condition`, `Condition == "Parkinson's disease"` in `source_experiments`)**: $K = 58$ studies ($156$ comparisons, $\binom{58}{2} = 1,653$ pairs; $972$ nonzero, $681$ zero), mean pairwise taxon Jaccard $= 0.041188$ ($4.12\%$), median $= 0.024390$ ($2.44\%$).

The nonzero/zero pair counts were corrected in
[#1165](https://github.com/CultureBotAI/microbial-trait-profiling/issues/1165)
against the archived integer intersection ledger, using the
[reference support-count audit](../scripts/adhoc/issue1117_audit_reference_pd_counts.py).
That audit pins the consumer manifest, verifies the summary and pair-ledger
bytes against its output inventory, checks complete unique unordered pairs,
and reconciles study/pair counts and positive versus zero intersections. It
does not reconstruct taxon sets or recompute means, medians or scientific tests.

---

## 5. Provenance and Core-Run Guards

`mtp.matrix_overlap_coverage` directly reads `count_profiles`, `experiments`,
`reported_taxa`, `source_experiments`, `source_reported_taxa`, `taxa_matrix`,
`feature_metadata` and `feature_catalog`, plus the declared trait list.
`contrast_selection` is not directly read by this sidecar. Core verification
and its consumption/coverage record also bind the required upstream provenance
products; those bindings do not imply that every bound product is parsed by the
sidecar's numerical calculations.

The sidecar enforces four provenance invariants before writing any output files:

1. **Tier 1+2 population guard**: `--core-run` must be a verified core bundle whose authoritative `contrast_population` and `direction_basis` fields are `"disease_vs_healthy_tier1_2_subject_controls"` and `"disease_vs_comparison"`. Conflicting same-name parameter values are refused; parameters cannot override those fields. Passing an `all_paired` core run or a tampered manifest raises `ValueError`.
2. **Durable output root guard (`#1049`, `#1051`)**: Both the library entry point and CLI use the shared policy in [`src/mtp/workflow/output_roots.py`](../src/mtp/workflow/output_roots.py), including its fixed roots, process temporary directory, `$TMPDIR` and resolved symlink spellings. A matched output path is refused unless temporary output is explicitly allowed (`--allow-temporary-output` for ephemeral CLI runs).
3. **Clean output directory guard**: `--output-dir` must be a new or empty directory so stale tables from earlier runs cannot survive alongside newly emitted tables.
4. **Core-bundle non-overlap preflight and bound provenance**: `preflight_consumer_outputs` runs before `output_dir.mkdir(...)` so an output path inside the core run bundle is rejected without mutating the core run directory. `provenance.json` records the sidecar's own `code_commit`, `script_sha256` for the primary module, a per-file `scripts` inventory, `methods_sha256` (`methods/matrix_coverage_and_pd_overlap.md`), `core_git_commit_sha`, and run-relative `steps/...` paths for all consumed core products, and is bound alongside the six primary outputs inside `core_consumption.json`.

The `scripts` inventory hashes `matrix_overlap_coverage.py`, `trait_lists.py`,
`matrix_overlap_inputs.py`, `human_host_corpus.py`, `workflow/coverage.py` and
`workflow/output_roots.py`, all beneath `src/mtp/`.

### Historical numerical reference identity

The archived F reference used core run
`run-88a731a42ddf47fd95ac2277726a25d9` at source commit
`511c891229b659e207cfe0f4281a97c3fc039a80`. Its core manifest SHA256 is
`0dddc45684fb2628046922a553a30e49d9ccfbe13de6b2d7566e9c8d8ccd02e7`;
the historical ten-trait list SHA256 is
`6a0b16238b6609962fccf5532d5aa59dbda771e1c4d45878cc0b079d7275947f`.
The archived sidecar's consumption record binds its numerical output tables to
that core. It predates the current sidecar producer binding: its `git_commit_sha`
identifies F, it lacks a separate sidecar commit and script inventory, and its
`provenance.json` is not among its consumer-bound outputs. It is therefore an
archived numerical reference, not independent proof of a clean sidecar producer
or a run satisfying every current provenance check.
