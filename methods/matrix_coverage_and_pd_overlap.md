# Matrix Trait Overlap, Tier 1+2 Coverage Cascade, and Parkinson's Stool Study Overlap

This note documents the definitions, data sources, and exact denominator semantics implemented in [`src/mtp/matrix_overlap_coverage.py`](../src/mtp/matrix_overlap_coverage.py) (`just matrix-overlap-coverage`, issue `#1053`). All interpretations follow [`methods/biological_interpretation.md`](biological_interpretation.md): KG-Microbe trait columns represent curated organism-level annotations rather than direct in vivo physiological measurements, and BugSigDB signatures record reported taxon abundance shifts across published comparisons.

---

## 1. Pairwise Trait Overlap in the Taxon × Trait Matrix (`matrix_trait_overlap.tsv`)

For any declared pair of retained (non-MediaDive) trait columns $(a, b)$ in `taxa_trait_matrix.tsv` (with $N_{\text{matrix}} = 5,829$ taxon rows in the `tier12-phenotype` core run):

- Let $T_a = \{i : M_{i,a} > 0\}$ and $T_b = \{i : M_{i,b} > 0\}$ denote the sets of matrix taxa carrying a positive annotation count for traits $a$ and $b$.
- The pairwise Jaccard overlap is:
  $$J(a, b) = \frac{|T_a \cap T_b|}{|T_a \cup T_b|}$$

### Reference Pairs (`tier12-phenotype`, $N_{\text{matrix}} = 5,829$)

| Trait A | Trait B | $|T_a|$ | $|T_b|$ | $|T_a \cap T_b|$ | $|T_a \cup T_b|$ | Jaccard (%) |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| `METPO:1000608` (`facultative`) | `METPO:1000605` (`facultatively anaerobic`) | 920 | 712 | 535 | 1,097 | 48.77% (48.8%) |
| `METPO:2000103\|GO:0030245` (`cellulose catabolic process`) | `METPO:2000103\|GO:0045493` (`xylan catabolic process`) | 189 | 210 | 109 | 290 | 37.59% (37.6%) |
| `METPO:1000608` (`facultative`) | `METPO:1000603` (`anaerobic`) | 920 | 1,884 | 493 | 2,311 | 21.33% (21.3%) |
| `METPO:1000608` (`facultative`) | `METPO:1000607` (`obligately anaerobic`) | 920 | 414 | 96 | 1,238 | 7.75% (7.8%) |

Both `METPO:1000603` (`anaerobic`, broad parent term with $1,884$ matrix taxa) and `METPO:1000607` (`obligately anaerobic`, strict term with $414$ matrix taxa) are included in `DEFAULT_MATRIX_OVERLAP_TRAITS` so downstream consumers can distinguish general anaerobic annotations from strict obligate anaerobe annotations without ambiguity.

---

## 2. Tier 1+2 Taxon Coverage by Rank (`tier12_coverage_by_rank.tsv`)

### Why `just taxa-coverage` (`scripts/create_taxa_coverage_histogram.py`) Computes a Different Quantity

[`scripts/create_taxa_coverage_histogram.py`](../scripts/create_taxa_coverage_histogram.py) (`just taxa-coverage`) summarizes the static `taxa_trait_matrix.tsv` table directly:

1. **Taxon universe**: It iterates over all $5,829$ rows in `taxa_trait_matrix.tsv` (which spans all BugSigDB reported taxa across hosts and experimental designs), rather than the $3,221$ unique NCBI taxa reported in the $1,870$ active Tier 1+2 signatures (`disease_vs_healthy_tier1_2_subject_controls`).
2. **Trait universe**: Unless MediaDive columns are explicitly excluded, `taxa_trait_matrix.tsv` includes $2,333$ `METPO:2000517|mediadive.medium:*` growth-medium columns alongside the $344$ retained biological trait columns. A taxon with only MediaDive medium annotations is unannotated under the $344$-trait biological analysis set.
3. **Mention and study-frequency granularity**: `scripts/create_taxa_coverage_histogram.py` does not join `bugsigdb_taxa.tsv` to the active Tier 1+2 signature keys (`count_profiles`), does not compute signature-level mention coverage, and does not stratify labeled versus unlabeled taxa by the number of studies ($n_{\text{studies}} = 1$, $n_{\text{studies}} \le 2$, median $n_{\text{studies}}$) in which each taxon appears.

### Exact Tier 1+2 Coverage Definitions (`mtp.matrix_overlap_coverage`)

In `mtp.matrix_overlap_coverage`:

- **Retained trait set ($344$ columns)**: All columns in `taxa_trait_matrix_columns.tsv` (`feature_metadata`) present in `taxa_trait_matrix.tsv` except `METPO:2000517|mediadive.medium:*` growth-medium terms. A taxon is **labeled** if and only if it has a positive value in at least one of those $344$ retained columns.
- **Active Tier 1+2 mentions ($19,388$ distinct `(profile, taxon_id)` mentions)**: `steps/export_tables/derived/bugsigdb_taxa.tsv` (`reported_taxa`) preserves repeated within-signature occurrences (`coverage_policy: preserve_repeated_occurrences`, $19,566$ raw rows across the $1,870$ active Tier 1+2 profiles). Deduplicating by `(study_id, experiment_id, inc_dec, taxon_id)` yields $19,388$ distinct `(profile, taxon_id)` mentions across $3,221$ unique NCBI taxa.
- **Key Tier 1+2 `tier12-phenotype` figures**:
  - **All ranks ($3,221$ unique taxa, $19,388$ mentions)**: $2,363$ labeled taxa ($73.36\%$) and $858$ unlabeled taxa ($26.64\%$); $16,743$ labeled mentions ($86.36\%$) and $2,645$ unlabeled mentions ($13.64\%$).
  - **Genus rank ($842$ unique genera, $8,495$ mentions)**: $739$ labeled genera ($87.77\%$) and $103$ unlabeled genera ($12.23\%$); $8,324$ of $8,495$ genus mentions ($97.99\%$, i.e. $98.0\%$) carry at least one retained trait annotation.
  - **Species rank ($1,983$ unique species, $7,952$ mentions)**: $1,287$ labeled species ($64.90\%$) and $696$ unlabeled species ($35.10\%$, i.e. $35.1\%$).
  - **Single-study rarity by labeling status**: Among the $858$ unlabeled taxa, $717$ ($83.57\%$, i.e. $84\%$) appear in only $1$ study and $805$ ($93.82\%$, i.e. $94\%$) appear in $\le 2$ studies (median $1$ study). Among the $2,363$ labeled taxa, $1,124$ ($47.57\%$, i.e. $48\%$) appear in $1$ study and $1,465$ ($62.00\%$, i.e. $62\%$) appear in $\le 2$ studies (median $2$ studies).

---

## 3. Two-Axis Filter Cascade (`filter_cascade.tsv`)

### Trait Axis

1. **KG-Microbe trait catalog (`kg_trait_catalog`)**: $2,677$ total trait columns in `feature_catalog.tsv` and `taxa_trait_matrix_columns.tsv` ($2,636$ of which have positive counts in the Tier 1+2 profile table after $41$ canary-0 columns with zero hits on the subset).
2. **Retained non-MediaDive traits (`retained_non_mediadive_traits`)**: $344$ traits retained after dropping $2,333$ `METPO:2000517|mediadive.medium:*` growth-medium columns:
   - $162$ METPO phenotypes (`METPO:*` without compound predicate pipes)
   - $102$ carbon-source traits (`METPO:2000006|*`, comprising $101$ `CHEBI:*` substrates plus `METPO:2000006|carbon_substrates:N-acetylgulcosamine`)
   - $68$ pathway traits (`METPO:2000103|pathways:*`)
   - $12$ GO biological process traits (`METPO:2000103|GO:*`)

### Study and Signature Axis (Distinguishing the Three "Full Corpus" Populations)

Across the workflow, three distinct upstream populations exist and must not be conflated:

- **All-hosts BugSigDB export (`all_hosts_bugsigdb_export`)**: $1,245$ studies and $4,435$ experiments in `bugsigdb_experiments.tsv` (`source_experiments`), spanning all host species and sequencing types. (By comparison, the `full-phenotype` / `full-off` core workflow runs retain $1,212$ studies, $4,299$ experiments, and $7,111$ non-empty directional profiles across all hosts.)
- **Canonical human 16S/WMS profile corpus (`full_human_16s_wms`)**: $1,047$ studies, $3,221$ experiments, and $5,629$ non-empty directional taxon lists (`Host species == "Homo sapiens"` and `Sequencing type in {"16S", "WMS"}`).
- **Paired human 16S/WMS corpus (`paired_human_16s_wms`)**: $932$ studies, $2,257$ paired experiments, and $4,514$ directional lists ($115$ studies and $1,115$ unpaired lists dropped).
- **Tier 1+2 disease vs healthy subject controls (`tier12_subject_controls`)**: $548$ studies, $935$ paired comparisons, and $1,870$ directional lists across $211$ conditions ($384$ studies and $2,644$ lists dropped; includes $355$ stool studies and $595$ stool comparisons where `Body site == "Feces"`).

---

## 4. Parkinson's Disease Stool Study Taxon Overlap (`parkinsons_stool_study_overlap_summary.tsv`)

For each scope (`tier12_paired` and `full_export`), experiments matching `Condition` containing `Parkinson` (case-insensitive) and `Body site == "Feces"` are grouped by `study_id` (using the active Tier 1+2 paired subset for `tier12_paired` and `source_experiments` (`steps/source_export_tables/source/derived/experiments.tsv`) for `full_export`). Each study $s$ is represented by the union of all reported NCBI taxon IDs $U_s$ across its active experiments in that scope. Across all $\binom{K}{2}$ unordered study pairs $(s_1, s_2)$:

$$J(s_1, s_2) = \frac{|U_{s_1} \cap U_{s_2}|}{|U_{s_1} \cup U_{s_2}|}$$

- **Tier 1+2 paired (`tier12_paired`)**: $K = 36$ studies ($70$ comparisons, $\binom{36}{2} = 630$ pairs), mean pairwise taxon Jaccard $= 0.053302$ ($5.33\%$, i.e. $5.3\%$), median $= 0.038099$ ($3.81\%$, i.e. $3.8\%$).
- **Full export (`full_export`)**: $K = 59$ studies ($157$ comparisons, $\binom{59}{2} = 1,711$ pairs), mean pairwise taxon Jaccard $= 0.041030$ ($4.10\%$, i.e. $4.1\%$), median $= 0.024390$ ($2.44\%$, i.e. $2.4\%$).
