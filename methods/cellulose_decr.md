# Cross-condition overlap of stool-decreased cellulose-annotated taxa

This note describes how the tables and figures in
[`data/analysis/cellulose_decr/`](../data/analysis/cellulose_decr/README.md) are
produced (#552). They come from `just cellulose-decr` (tables) and
`just cellulose-decr-figures` (figures). The analysis asks two questions about
taxa annotated `METPO:2000103|GO:0030245` (cellulose catabolic process) in
KG-Microbe:

1. Within Parkinson's disease (PD), HIV infection and colorectal cancer (CRC), do
   curated decreased signatures carry a larger cellulose annotation share than the
   paired increased signatures?
2. Do the three conditions report the same cellulose-annotated taxa decreased in
   feces more often than chance? And against which chance model?

Both are questions about curated signature annotations. Neither measures
cellulose degradation, fiber intake, short-chain fatty acids or any other
physiology (see [`biological_interpretation.md`](biological_interpretation.md)).

## Inputs and population

| Input | Use |
| --- | --- |
| `trait-visualization/study_trait_profile_study_metadata.tsv` | Defines the 5,829 retained human profiles (one per study, experiment and direction) and their non-MediaDive annotation totals. Its condition and body-site columns are **not** used (see below). |
| `data/derived/experiments.tsv` | Condition, body site and host of each experiment. It matches `data/raw/full_dump.csv` for all 4,435 experiments. |
| `data/derived/study_trait_profiles_wide.tsv` | Annotation counts per profile, used for the cellulose annotation share. |
| `data/derived/bugsigdb_taxa.tsv` | Taxa reported in each signature. |
| `data/derived/taxa_trait_matrix.tsv` and `taxa_trait_matrix_columns.tsv` | Which of the 5,829 matrix taxa bear each trait, and the trait labels. |
| `data/kg-microbe/merged-kg_edges.tsv` | Which taxa carry the annotation edge themselves rather than through a descendant. |
| `data/derived/contrast_orientation/contrast_manifest.tsv` | Contrast class and orientation of every experiment ([`contrast_orientation.md`](contrast_orientation.md)). |
| `data/derived/studies.tsv`, `data/raw/full_dump.csv` | Study identifiers, authors, journals and locations for the evidence tables, and the ontology IDs of condition labels. |

The run's `cellulose_decr_provenance.json` records:

- the SHA-256 and size of every input;
- the contrast-classifier version;
- the code commit, seeds and package versions.

Numbers in this note come from BugSigDB vintage `2026-04-24_00:41_UTC` and
classifier `contrast-orientation-v2`.

**Pinned inputs.** The inputs are gitignored, and other work rebuilds them in
place. The provenance therefore keeps the reference hashes in `pinned_inputs`,
apart from the `inputs` a run actually read. A run stops in four cases:

- an input's hash differs from those pins;
- the pinned provenance is missing;
- the manifest's `classifier_version` differs from
  `mtp.contrast_orientation.CLASSIFIER_VERSION`;
- an input changes while the run reads it, which a second hash just before
  writing detects.

Two options change this:

- `--allow-input-drift` runs on changed inputs and lists them, but carries the
  old pins forward, so the next ordinary run still refuses them.
- `--repin` adopts the current inputs as the pins.

See #562 and #571. The pins describe the v2 contrast manifest and the profile metadata
table before #555 relabelled it per experiment. Until #580 re-pins the analysis on the
current manifest, a checkout with a newer classifier (v3 since #554) or the relabelled
metadata refuses to run, and the committed tables reproduce by checking out the
provenance's `code_commit` with the pinned inputs (#642).

**Condition and body site come from each experiment's own record**
(`experiments.tsv`). The October 16 outputs took them from the profile metadata
table, which then held one modal condition and body site per *study* (#555), so
1,112 of its 5,829 profiles carried a label that belongs to another experiment of
the same study. Switching to `experiments.tsv` changed 61 experiments in the
three conditions, among them saliva and nasal experiments that had been counted
as PD feces; `label_source_corrections.tsv` lists them.

`--label-source study_modal` recomputes the study-modal labels from
`full_dump.csv` (#594), whatever the metadata table now holds. It reproduces the
October 16 sets, exact tests, catalog counts and feces direction statistics. U5
(#567), U4's Monte Carlo p-values, the oral and other-site compartment rows and
two pairs of benchmark ranks differ.

**Condition labels.** A condition matches only its exact BugSigDB label
(primary). Other single labels for the same disease are listed in
`condition_label_inventory.tsv` with their status, so no exclusion is silent:

| Disease | Synonym labels |
| --- | --- |
| HIV | "HIV-1 infection" (EFO:0000180), "Human immunodeficiency virus" (NCBITaxon:12721) |
| CRC | colorectal carcinoma (EFO:1001951), colorectal adenocarcinoma (EFO:0000365), metastatic colorectal cancer (EFO:1001480), colon carcinoma (EFO:1001950) |

These synonyms are counted only in the `synonyms_merged` sensitivity, or
throughout with `--condition-labels synonyms`. Colorectal adenoma, colorectal
stromal tumour, HIV mother-to-child transmission, simian immunodeficiency virus
and multi-condition labels such as `Constipation,Parkinson's disease` are
different conditions and stay out of both.

**Body sites.** Feces is the exact body site `Feces`. Labels naming more than one
site, such as `Feces,Blood` or `Oral cavity,Feces`, are grouped with other non-GI.
Single sites are grouped by keyword into non-fecal GI (mucosa or lumen), oral or
oropharyngeal, and other non-GI. The grouping actually applied to the three
conditions is written to the summary JSON.

## Direction tests

For each profile, the cellulose annotation share is the cellulose count divided
by the sum of all retained non-MediaDive counts. For each experiment with both
signatures, the paired difference is the share in the increased signature minus
the share in the decreased signature. *Down* means the decreased signature had
the larger share. Ties are excluded.

- **Experiment level:** an exact two-sided sign test and a Wilcoxon signed-rank
  test over paired experiments. Experiments from one study are not independent,
  so these p-values are optimistic.
- **Study level:** one vote per study, the sign of the mean paired difference over
  its experiments. The tests are an exact two-sided sign test and a Wilcoxon test
  on the study means. A second rule votes by the majority of each study's
  experiment signs (`study_majority_*`). No study contributes to more than one of
  the three conditions.
- **Populations:**
  - all contrasts;
  - disease-versus-healthy contrasts, with reversed experiments relabelled so
    *decreased* means lower in disease (`disease_vs_healthy_oriented`);
  - the same without experiments flagged `weak_comparison_evidence`;
  - for HIV, the oriented population without experiments whose group names or
    descriptions mention MSM. That last filter is lexical. It does not adjust for
    sexual behaviour, which was not recorded for most studies.
- **Condition definitions.** `condition_definition_sensitivity.tsv` repeats the
  feces tests under these definitions:
  - the synonym labels;
  - the HIV– versus HIV+ experiment the classifier calls disease versus disease;
  - three ART-treated HIV+ versus uninfected experiments the classifier calls
    treatment exposure (#569);
  - without the one HIV study whose "feces" was rectal mucosa (#566);
  - without the CRC experiment the classifier counts as disease vs healthy
    although it compares cancerous with adjacent non-cancerous tissue of the
    same patients (#569).

  `condition_definition_interaction.tsv` repeats the Feces-vs-GI test under each
  of these definitions.
- **Compartments:**
  - For each compartment and condition, the table reports the smallest two-sided
    p a sign test can reach with the observed untied pairs, and the power to
    detect the pooled fecal share at alpha 0.05.
  - Fisher's exact test of direction by compartment compares feces with each
    other compartment, at experiment and study level. It reports the conditional
    maximum-likelihood odds ratio and its exact 95% interval.
  - Pooling the three conditions mixes condition with compartment, because they
    sample the compartments unevenly. The pooled rows therefore also give an
    exact conditional test of a common odds ratio across the three conditions,
    with its conditional MLE, exact interval and Mantel-Haenszel estimate.
  - As in Fisher's test, the p-value counts outcomes no more probable than the
    observed one, while the interval inverts two one-sided tails, so the two can
    disagree near 0.05.
  - Units are treated as independent, and a few studies sampled two compartments.
- **Multiplicity:** `significant_results_summary.tsv` adds Benjamini-Hochberg
  q-values over its 29 tests.

## Overlap tests

A condition's set holds the matrix taxa that bear the trait and that at least one
of its feces decreased signatures reports, at any rank. Overlap is measured three
ways:

- the Jaccard similarity J and distance d_J = 1 - J for each pair;
- their mean over the three pairs;
- the number of taxa in all three sets.

| Null | What is drawn | What it holds fixed |
| --- | --- | --- |
| U1 | Uniform draws from the 189 matrix taxa with the annotation | set sizes |
| U2 | Uniform draws from the annotated taxa ever reported in a feces signature | set sizes |
| U3 | Uniform draws from the annotated taxa ever reported decreased in feces | set sizes |
| U4 | A random trait placed on as many of the conditions' decreased feces taxa as the observed union | each condition's decreased feces taxa |
| U5 | Draws from U3, the taxa an observed set can contain, weighted by the number of other feces decreased studies reporting the taxon plus 0.5 | set sizes, reporting prevalence |
| U6 | Whole study signatures: 49, 19 and 31 distinct feces-decreased studies drawn from all 646, each giving one of its condition-level signatures | number of studies, taxa reported together |
| U6 size-matched | U6 draws whose three set sizes are each within max(2, 10%) of the observed | as U6, plus set sizes |
| Within-condition reference | Random subsets of each condition's own decreased feces taxa (any trait), of the observed set sizes | each condition's decreased feces taxa, set sizes |

The within-condition reference is not one of the six nulls. It shows how little
three sets of the observed sizes would overlap if each were drawn from its own
condition's decreased feces taxa, whatever their traits. It is written to the
`within_cond_random_subset_*` columns of `pairwise_jaccard_significance.tsv` and
to `null_models.within` in the summary JSON.

**Exact tests.** U1 to U3 have exact hypergeometric tails for each pair and an
exact convolution for the 3-way core.

**Monte Carlo tests.** Each Monte Carlo null makes 100,000 draws with a fixed
seed, and p = (1 + draws at least as extreme) / (1 + draws). The size-matched U6
keeps the 5,284 of its draws that match. An earlier U5 drew from U2, which
includes 27 taxa never reported decreased in feces. Those taxa can never appear
in an observed set, so that version overstated significance (#567).

**Which null to trust.** U1 to U4 treat taxa as independent draws, which they are
not. Genera, their families and their phyla are reported together, and some taxa
are reported in most feces studies whatever the condition. U6 keeps both features
and is the reference null for claims about overlap.
`condition_definition_overlap.tsv` repeats the overlap and U6 under each
condition definition.

**Sensitivity table.** It repeats the overlap under other set definitions:

- profiles with a positive cellulose z-score only;
- genus and species ranks only;
- taxa with more decreased than increased profiles or studies;
- all body sites;
- disease-versus-healthy contrasts;
- the synonym labels;
- without the rectal-mucosa study;
- every decreased feces taxon regardless of trait.

**Trait benchmark.** It ranks every matrix trait with at least 10 decreased feces
taxa in each condition by mean pairwise d_J. Its ranks describe this corpus and
are not a test; traits share taxa.

## Outputs

`just cellulose-decr` runs `mtp.cellulose_decr.analysis`
(`src/mtp/cellulose_decr/analysis.py`), which writes these committed tables to
`data/analysis/cellulose_decr/`:

| File | Contents |
| --- | --- |
| `cellulose_decr_summary.json` | Every result, definition, universe, null summary, label inventory and label-source comparison |
| `cellulose_decr_provenance.json` | Input hashes and pins, classifier version, code commit, seeds, package versions |
| `pairwise_jaccard_significance.tsv` | Pairwise, mean and 3-way overlap under U1–U6 and the within-condition reference |
| `sensitivity_definitions_jaccard.tsv` | Overlap under eleven set definitions |
| `compartment_jaccard_and_direction_stats.tsv` | Direction tests by population, compartment and condition, at experiment and study level |
| `compartment_direction_interaction.tsv` | Fisher tests of direction by compartment, crude and stratified |
| `condition_definition_sensitivity.tsv` | Feces direction tests under each condition definition |
| `condition_definition_overlap.tsv` | Overlap and U6 under each condition definition |
| `condition_definition_interaction.tsv` | Feces-vs-GI tests under each condition definition |
| `condition_label_inventory.tsv` | Every human condition label naming the three diseases, with its status |
| `cellulose_decreased_taxa_catalog.tsv` | The taxa in any condition's set, with counts and flags |
| `profile_and_study_pair_jaccard.tsv` | Overlap between individual profiles and studies |
| `trait_convergence_benchmark.tsv` | Every benchmarked trait with its label, predicate and rank |
| `significant_results_summary.tsv` | Headline tests with unit, population, p and BH q |
| `label_source_corrections.tsv` | Experiments whose study-level label differed from their own |
| `four_highlighted_pairs_taxon_drivers.tsv`, `four_highlighted_pairs_study_evidence.tsv`, `four_highlighted_pairs_summary.json` | Taxa and studies behind the four exploratory leads |

`just cellulose-decr-figures` runs `mtp.cellulose_decr.figures`
(`src/mtp/cellulose_decr/figures.py`), which draws the figures from those tables
only. The figures are written beside the tables and are gitignored, as
CONTRIBUTING.md requires for images outside `trait-visualization/` (#560).
Reruns on the same inputs give identical tables and figures; only the provenance
timestamp changes. Figure provenance is in `figures/figures_provenance.json`,
with the commit that drew the figures (`+dirty` when the code or tables were
modified), and `figures/figure_captions.json` holds the slide text for each
figure, worded from the same tables (#556).

```bash
just cellulose-decr
just cellulose-decr-figures
```

The recipe records `git rev-parse HEAD` as the code commit, with `+dirty` when
the analysis scripts or the repository's code paths (`src`, `scripts`,
`methods`, `schema`, `Justfile`, `pyproject.toml`, `uv.lock`) differ from it.

**Building the inputs.**

- BugSigDB: `just download`, `just ingest` and `just study-tsvs`.
- KG-Microbe: `just download-kg`, `just taxa-trait-matrix` and
  `just study-trait-profiles`.
- Contrast manifest: `just contrast-orientation`.

From a worktree, `--input-root <checkout>` points the analysis at a checkout
holding those inputs. For the October 16 labels, pass
`--label-source study_modal --output-dir <elsewhere>`, and for synonym labels
throughout, `--condition-labels synonyms --output-dir <elsewhere>`. Without
`--output-dir` either run would overwrite the committed tables. The code is the
`mtp.cellulose_decr` package (`src/mtp/cellulose_decr/`, #557), tested in
`tests/test_cellulose_decr.py` and `tests/test_cellulose_decr_figures.py`.

## Limits

- **Annotations are inherited.** A taxon bears an annotation when it or any
  descendant has the KG edge. None of the 17 taxa in all three sets has an edge
  of its own: genera, families and phyla alike inherit it from annotated species
  or strains.
- **Sets measure reported labels, not consistent loss.** A taxon enters a set once
  any study reports it decreased. Several core taxa are reported increased in more
  studies than decreased in some condition.
- **The HIV direction depends on definitions.** It changes with which labels and
  contrasts count as HIV disease-versus-healthy experiments: the study-level p
  runs from 0.002 to 0.18 across the definitions above.
- **Compartments are not comparable.** Sampling effort and contrast types differ
  greatly between compartments, so stronger results in feces do not show that the
  pattern is specific to feces.
- **Nothing here is confirmatory.** The three conditions and the trait are
  exploratory leads picked from this corpus (`EXPLORATORY_LEADS` in
  `src/mtp/deck_visualizations.py`).
