# Cellulose-annotated taxa reported decreased in feces in Parkinson's disease, HIV infection and colorectal cancer

This directory compares the taxa annotated `METPO:2000103|GO:0030245`
(**cellulose catabolic process**) that curated BugSigDB signatures report decreased
in human feces in Parkinson's disease (PD), HIV infection (HIV) and colorectal
cancer (CRC). The quantitative blocks below bind to the committed tables and JSON
in this directory. Their sources are checked section by section by
`tests/test_cellulose_readme_numbers_799.py`, including a census that refuses
unbound numeric claims. The historical notes describe changes; they no longer
quote old-run values that this committed bundle cannot establish.

<!-- numbers: producer -->
The main analysis tables record code commit `511c891` and classifier `contrast-orientation-v6`. The provenance pins 10 input files.
<!-- /numbers -->

The analysis runs through `just cellulose-decr`
([`src/mtp/cellulose_decr/analysis.py`](../../../src/mtp/cellulose_decr/analysis.py)).
Its method is described in
[`methods/cellulose_decr.md`](../../../methods/cellulose_decr.md).
The Tier 1–2 catalog also has its own entries in the committed summary; its
displayed counts are checked directly against that catalog.

These are annotation counts of curated signature taxa. They do not measure
cellulose degradation, fiber intake or short-chain fatty acids. A taxon counted
here was reported decreased in a study, not consistently lost
([`methods/biological_interpretation.md`](../../../methods/biological_interpretation.md)).
The trait and the conditions were chosen as an exploratory lead from this
same corpus, so every p-value here is nominal. Section 4 adds adjusted values.

## What changed on 2026-10-23 (final rerun, #799)

The tables were regenerated and re-pinned with `--repin`. The contrast manifest
uses the updated classifier, the profile counts use the primary `phenotype`
trait scope (#788), and the taxa-trait matrix, column table and profile metadata
were rebuilt. Sections 1–9 now quote those tables.

The CRC Feces-vs-GI study test no longer depends on whether the paired-tissue
experiment is counted; section 5 shows the current values. The benchmark's
membership and ranks changed; section 7 gives the committed ranking. The
disease-vs-healthy study-level PD q-value had also changed in the earlier #767
rerun without being reflected here. Section 4 now binds each q-value to its
named row in the significant-results table.

The earlier recount helper checked selected strings across several documents.
The dedicated README tests now bind each displayed result to its section and
source row. Unsupported before/after values have been removed rather than
silently attributed to the current rerun.

## What changed on 2026-10-22 (#781)

`condition_definition_sensitivity.tsv` and `condition_definition_interaction.tsv`
had the defect that #769 fixed in the overlap table: case-control definitions
were repeated on an all-contrasts population that did not apply their filters.
The analysis now tests those definitions on disease-versus-healthy contrasts
only. The exact-label rows remain, so removing copies does not change the
reported range over definitions. The figures use the disease-versus-healthy
rows for these comparisons.

## What changed on 2026-10-21 (#769)

`condition_definition_overlap.tsv` originally built every condition-definition
set from all contrasts. The case-control definitions now apply to
disease-versus-healthy contrasts only, and every row names its
`contrast_population`. Exact-label, synonym and rectal-mucosa definitions have
rows on each population, so comparisons stay within the same population.

## What changed on 2026-10-19 (#580, #672)

Disease-versus-healthy populations retain `disease_vs_healthy` contrasts whose
`contrast_type` is `case_vs_control_subject`. Tumour-vs-adjacent-tissue and
comorbidity contrasts no longer count in the primary reading. Tier 3 controls
still count; sensitivity scenarios set them aside.

The classifier, rather than hand-coded scenarios, makes the HIV status and
paired-tissue calls (#569, #554). Contrasts whose cases are named as treated,
restricted control tiers, paired tissue and comorbidity are separate sensitivity
definitions. Yachida's polyp-bearing controls are subject controls (#672).
The summary's `case_control_inventory` records these distinctions.

The earlier before/after direction counts and p-values are not preserved here as
current results. Their dated explanations remain historical; the current
direction and compartment tables below are the quantitative record.

## What changed on 2026-10-17 (#552)

The earlier display had no generator. The rebuilt analysis distinguished
experiment labels from study-modal labels: study-modal labels had put saliva
and nasal experiments into PD feces and moved other experiments out. The
comparison saved in the committed summary establishes the following values.

<!-- numbers: label_mode_comparison -->
The committed rerun compares study-modal and experiment labels on its own inputs; these are not preserved results from the earlier dated runs.

| Condition | Study-modal down : up : tied | Experiment down : up : tied |
| --- | --- | --- |
| PD | 56 : 26 : 19 | 57 : 24 : 16 |
| HIV | 19 : 6 : 3 | 18 : 7 : 3 |
| CRC | 30 : 7 : 10 | 32 : 6 : 11 |

The label-source corrections table records 61 experiments in the target conditions whose labels differ.
<!-- /numbers -->

The benchmark carries canonical trait labels from the matrix column table.
Earlier figures and slides had mislabeled acetogenesis, pleomorphic shaped,
spiral shaped, thermophilic and filament shaped. Spore shaped is not a matrix
column in this analysis.

The rebuild added a whole-signature null, study-level direction tests,
disease-versus-healthy subsets, compartment interactions and catalog flags.
It also added synonym-label sensitivity (#561), pinned-input refusal (#562),
direct-edge checks (#565), the rectal-mucosa sensitivity (#566), revised null
definitions and study-majority tests (#567), and reproducible figures (#560).

The old Monte Carlo estimates, compartment universes, benchmark ranks and
literature/catalog before-and-after counts are omitted because the current
bundle cannot establish each historical endpoint. The supported label-mode
comparison above is explicitly a comparison under the committed rerun.

## 1. Overlap of the taxon sets

A condition's set is every matrix taxon that bears the trait, at any rank, and
that any of the condition's feces decreased signatures reports.

<!-- numbers: overlap_table -->
| Pair | Set sizes | Shared | Union | J | d_J |
| --- | --- | --- | --- | --- | --- |
| PD vs. HIV | 38, 24 | 21 | 41 | 0.512 | 0.488 |
| PD vs. CRC | 38, 34 | 27 | 45 | 0.600 | 0.400 |
| HIV vs. CRC | 24, 34 | 18 | 40 | 0.450 | 0.550 |
| Mean of the pairs |  |  |  | 0.521 | 0.479 |
| All conditions |  | 17 | 47 |  |  |
<!-- /numbers -->

<!-- numbers: overlap_support -->
32 of the 47 catalog taxa occur in at least 2 condition sets.
<!-- /numbers -->

Sources:
[`pairwise_jaccard_significance.tsv`](pairwise_jaccard_significance.tsv) and
[`cellulose_decreased_taxa_catalog.tsv`](cellulose_decreased_taxa_catalog.tsv).

## 2. Is the overlap unusual?

<!-- numbers: null_draws -->
Monte Carlo nulls use 100,000 draws. The U6 size-window subset retains 5,557 draws under the recorded tolerance: `each set size within max(2, 10% of observed)`. Its set sizes need not equal the observed sizes. The displayed central intervals use the saved 2.5th and 97.5th percentiles. Values at the simulation floor are displayed as ≤ 1e-5, a rounded display convention rather than the exact finite-simulation minimum.
<!-- /numbers -->

<!-- numbers: null_table -->
| Null | Expected mean d_J (central interval) | p | Expected core (central interval) | p |
| --- | --- | --- | --- | --- |
| U1: uniform, 189 annotated matrix taxa | 0.909 (0.865–0.946) | ≤ 1e-5 | 0.87 | exact 4.8e-23 |
| U2: uniform, 135 annotated taxa seen in feces | 0.868 (0.817–0.911) | ≤ 1e-5 | 1.70 | exact 7.5e-18 |
| U3: uniform, 108 annotated taxa ever decreased in feces | 0.829 (0.772–0.878) | ≤ 1e-5 | 2.66 | exact 2.3e-14 |
| U4: random trait on 466 decreased feces taxa | 0.773 (0.668–0.873) | ≤ 1e-5 | 5.14 (2–9) | ≤ 1e-5 |
| U5: reporting-weighted U3 taxa | 0.586 (0.502–0.667) | 0.0086 | 13.7 (10–17) | 0.072 |
| U6: whole study signatures | 0.527 (0.420–0.634) | 0.19 | 18.3 (12–24) | 0.73 |
| U6: retained size-window draws | 0.514 (0.419–0.608) | 0.23 | 16.9 (13–21) | 0.59 |
| Within-condition reference, any trait | 0.971 (0.944–0.994) | ≤ 1e-5 | 0.162 (0–1) | ≤ 1e-5 |
<!-- /numbers -->

The mean-distance intervals and simulated-core intervals come from
`cellulose_decr_summary.json:null_models`. Uniform-null core expectations and
exact probabilities come from `pairwise_jaccard_significance.tsv`; they are not
the simulated-core estimates.

- **Independent draws (U1–U4).** These nulls treat taxa as independent, making
  the overlap look extreme. The within-condition reference draws from each
  condition's own decreased feces taxa, whatever their traits.
- **Weighted by reporting (U5).** Weighting taxa that can appear in a set by
  their reporting frequency leaves a weaker signal. The older sampling universe
  also admitted taxa never reported decreased in feces, inflating significance
  (#567).
- **Whole study signatures (U6).** Taxa reported together remain together.
  The expected core depends on the simulated set sizes, which is why the
  size-window comparison is reported separately.

<!-- numbers: null_units -->
U6 samples whole signatures: PD 47, HIV 19, CRC 31, from 629 feces studies with a decreased signature. Its pairwise p-values are PD vs. HIV: 0.21, PD vs. CRC: 0.098, HIV vs. CRC: 0.55.
<!-- /numbers -->

Most shared taxa are ones that feces studies often report decreased, whatever
the condition. These null tests do not establish condition-specific convergence.
The following all-contrasts sensitivity rows come from
`condition_definition_overlap.tsv`.

<!-- numbers: null_scenarios -->
| All-contrasts condition definition | PD / HIV / CRC | Core / union | U6 mean p | U6 core p |
| --- | --- | --- | --- | --- |
| Synonym labels counted under their disease | 38 / 27 / 37 | 20 / 50 | 0.14 | 0.54 |
| Exact labels without bsdb:24451087 (rectal mucosa curated as Feces) | 38 / 23 / 34 | 16 / 47 | 0.25 | 0.79 |
<!-- /numbers -->

## 3. Core taxa

The table gives stool studies reporting each taxon decreased / increased.
“Net decreased” means more decreased than increased studies in every target
condition (`is_strictly_net_decreased_all3`).

<!-- numbers: core_table -->
| Taxon | Rank | PD | HIV | CRC | Net decreased |
| --- | --- | --- | --- | --- | --- |
| Roseburia | genus | 20 / 4 | 2 / 0 | 6 / 0 | yes |
| Blautia | genus | 13 / 2 | 1 / 3 | 6 / 1 | no |
| Clostridium | genus | 7 / 3 | 5 / 2 | 6 / 3 | yes |
| Bacteroides | genus | 6 / 3 | 8 / 0 | 2 / 2 | no |
| Ruminococcus | genus | 6 / 2 | 4 / 0 | 4 / 3 | yes |
| Lachnospira | genus | 5 / 2 | 3 / 0 | 6 / 1 | yes |
| Anaerostipes | genus | 2 / 1 | 3 / 0 | 3 / 1 | yes |
| Butyrivibrio | genus | 2 / 2 | 1 / 2 | 1 / 0 | no |
| Eubacterium | genus | 1 / 1 | 1 / 1 | 3 / 2 | no |
| Lachnospiraceae | family | 20 / 7 | 5 / 3 | 9 / 5 | yes |
| Bacillota | phylum | 8 / 5 | 4 / 2 | 4 / 1 | yes |
| Clostridiaceae | family | 8 / 2 | 3 / 0 | 2 / 2 | no |
| Oscillospiraceae | family | 8 / 13 | 4 / 0 | 2 / 1 | no |
| Clostridia | class | 5 / 7 | 4 / 0 | 7 / 1 | no |
| Eubacteriales | order | 5 / 6 | 4 / 1 | 4 / 1 | no |
| Bacteroidia | class | 4 / 2 | 2 / 0 | 1 / 4 | no |
| Bacteroidota | phylum | 3 / 0 | 3 / 1 | 1 / 4 | no |
<!-- /numbers -->

<!-- numbers: core_annotations -->
The core contains 17 taxa: 9 genera and 8 family-to-phylum nodes. 0 core taxa have direct annotation edges. The KG has 186 NCBITaxon subjects carrying the edge. Among matrix trait-bearing taxa, 0 of 49 genera and 38 of 44 species have direct edges. The catalog union contains 7 taxa with direct edges, 7 of them species.
<!-- /numbers -->

A matrix taxon bears the trait when it or a descendant has the annotation edge,
so the core annotations are inherited. Examples include *Roseburia* from
*R. faecis*, *Lachnospira* from *L. multipara*, *Anaerostipes* from *A. hadrus*,
and *Bacteroides* from a *B. cellulosilyticus* strain. *Clostridium* and
*Ruminococcus* inherit annotations from several species.

<!-- numbers: core_net -->
Annotated matrix taxa at or below the named nodes are Lachnospiraceae: 17; Bacillota: 86. Net decreased across all conditions: 7 taxa by studies (5 genera), and 7 by profiles.
<!-- /numbers -->

The profile-based net-decreased definition includes *Bacteroides* instead of
*Anaerostipes*. `is_direct_or_genus_species` is a rank flag for the core, not
evidence of a direct annotation.

<!-- numbers: core_rectal -->
Eubacterium has 1 HIV decreased-study support. Excluding the rectal-mucosa study leaves 16 core taxa.
<!-- /numbers -->

The rectal-mucosa study is McHardy et al. (`bsdb:24451087`); BugSigDB curates its
sponge sample as feces (#566). The sensitivity changes the retained population,
not the biological meaning of the inherited annotation.

## 4. Direction within each condition

Each row compares paired experiments in feces. “Down” means the taxa annotated
with cellulose catabolic process made up a larger annotation share of the
decreased signature than of the increased signature. Study rows give each study
an equally weighted sign of its mean paired difference.

<!-- numbers: direction_table -->
| Condition | Contrasts | Experiments down : up : tied | Sign test p | Studies down : up : tied | Sign test p | Study Wilcoxon p |
| --- | --- | --- | --- | --- | --- | --- |
| PD | all | 57 : 24 : 16 | 3.2e-4 | 26 : 14 : 3 | 0.081 | 0.013 |
| PD | disease vs healthy | 46 : 15 : 11 | 8.8e-5 | 23 : 10 : 3 | 0.035 | 0.019 |
| HIV | all | 18 : 7 : 3 | 0.043 | 11 : 5 : 2 | 0.21 | 0.32 |
| HIV | disease vs healthy | 13 : 3 : 3 | 0.021 | 10 : 2 : 2 | 0.039 | 0.18 |
| CRC | all | 32 : 6 : 11 | 2.4e-5 | 19 : 2 : 7 | 2.2e-4 | 2.0e-4 |
| CRC | disease vs healthy | 22 : 5 : 8 | 0.0015 | 16 : 2 : 7 | 0.0013 | 0.0010 |
| Pooled | all | 107 : 37 : 30 | 4.4e-9 | 56 : 21 : 12 | 8.2e-5 | 5.6e-5 |
| Pooled | disease vs healthy | 81 : 23 : 22 | 9.3e-9 | 49 : 14 : 12 | 1.1e-5 | 8.8e-5 |
<!-- /numbers -->

<!-- numbers: direction_fdr -->
Benjamini–Hochberg correction covers the 29 rows in `significant_results_summary.tsv`; these are the study-direction q-values.

| Condition | All contrasts q | Disease vs healthy q |
| --- | --- | --- |
| PD | 0.10 | 0.054 |
| HIV | 0.23 | 0.056 |
| CRC | 8.0e-4 | 0.0035 |
| Pooled | 4.3e-4 | 1.1e-4 |
<!-- /numbers -->

Source: `compartment_jaccard_and_direction_stats.tsv` for the direction table;
`significant_results_summary.tsv` for the table-scoped BH correction. These are
the historical bundle's recorded tests, not a new multiplicity policy or rerun.

Collapsing to studies weakens the PD and HIV sign tests. The study-level
Wilcoxon and sign tests need not agree because they use different information.
Study majority voting is a separate collapse rule:

<!-- numbers: direction_majority -->
| Condition | Population | Study majority down : up : tied | p |
| --- | --- | --- | --- |
| PD | all_contrasts | 27 : 12 : 4 | 0.024 |
| PD | disease_vs_healthy_oriented | 23 : 9 : 4 | 0.020 |
| CRC | disease_vs_healthy_oriented | 14 : 2 : 9 | 0.0042 |
<!-- /numbers -->

Weak-comparison and MSM-name filters are separate population definitions:

<!-- numbers: direction_additional_filters -->
| HIV population | Experiments | p | Studies | p |
| --- | --- | --- | --- | --- |
| disease_vs_healthy_strong_evidence | 12 : 2 : 1 | 0.013 | 9 : 1 : 0 | 0.021 |
| disease_vs_healthy_oriented_excl_msm_mention | 12 : 3 : 3 | 0.035 | 10 : 2 : 2 | 0.039 |
<!-- /numbers -->

**The HIV result is sensitive to population definitions.** These rows bind to
`condition_definition_sensitivity.tsv`:

<!-- numbers: direction_hiv -->
| HIV disease-vs-healthy definition | Experiments | p | Studies | p |
| --- | --- | --- | --- | --- |
| Exact condition labels (primary analysis) | 13 : 3 : 3 | 0.021 | 10 : 2 : 2 | 0.039 |
| Synonym labels counted under their disease | 13 : 5 : 3 | 0.096 | 10 : 4 : 2 | 0.18 |
| Exact labels; contrasts whose cases are named as treated set aside | 12 : 1 : 3 | 0.0034 | 10 : 0 : 2 | 0.0020 |
| Exact labels; Tier 1-2 controls only (hospital, comorbid, tissue and treated-case contrasts set aside) | 12 : 1 : 3 | 0.0034 | 10 : 0 : 2 | 0.0020 |
| Exact labels; tumour-vs-adjacent-tissue (paired tissue) contrasts also counted | 13 : 3 : 3 | 0.021 | 10 : 2 : 2 | 0.039 |
| Exact labels; comorbidity contrasts also counted | 13 : 3 : 3 | 0.021 | 10 : 2 : 2 | 0.039 |
| Exact labels without bsdb:24451087 (rectal mucosa curated as Feces) | 12 : 3 : 3 | 0.035 | 9 : 2 : 2 | 0.065 |
<!-- /numbers -->

The named-treated-case flag reads group names, not actual treatment exposure
(#676). Some treated participants are described only in group definitions and
are therefore retained by that lexical scenario. The MSM-name filter likewise
cannot separate HIV from sexual behaviour. The literature discussion in
`FOUR_HIGHLIGHTED_PAIRS_LITERATURE_SUPPORT.md` supplies the biological context;
these tests do not identify its causal contribution.

Synonym labels change the retained studies and experiments:

<!-- numbers: direction_synonyms -->
| Synonym labels | Population | Experiments | p | Studies | p |
| --- | --- | --- | --- | --- | --- |
| HIV | all_contrasts | 19 : 10 : 3 | 0.14 | 12 : 7 : 2 | 0.36 |
| CRC | all_contrasts | 37 : 10 : 12 | 9.8e-5 | 22 : 3 : 7 | 1.6e-4 |
| CRC | disease_vs_healthy_oriented | 25 : 6 : 8 | 8.8e-4 | 18 : 3 : 7 | 0.0015 |

The pooled study-level sign test stays at p ≤ 1.21e-4 in every scenario.
<!-- /numbers -->

The former corpus-wide comparison cited an untracked
`data/derived/contrast_orientation/disease_control_comparison.tsv`. Its
quantitative claim is removed because this committed directory cannot verify
that external run. It must not be read as another result from the bound tables.

## 5. Outside feces

<!-- numbers: compartment_table -->
| Compartment | Decreased-signature studies PD / HIV / CRC | Pooled experiments down : up |
| --- | --- | --- |
| Feces | 47 / 19 / 31 | 107 : 37 |
| Non-fecal GI (mucosa/lumen) | 2 / 4 / 14 | 9 : 12 |
| Oral / Oropharyngeal | 2 / 2 / 4 | 3 : 5 |
| Other non-GI | 3 / 4 / 0 | 4 : 2 |
<!-- /numbers -->

<!-- numbers: compartment_power -->
Outside feces, individual-condition strata other than CRC GI have at most 5 untied experiments. CRC GI has 15, with 4 pointing down. The non-fecal pooled strata have 21, 8, 6 untied experiments, in the table order.
<!-- /numbers -->

Sparse strata have limited sign-test resolution. The interactions compare
direction between compartments rather than treating a weak within-compartment
test as evidence of no difference.

<!-- numbers: interaction_table -->
| Comparison | Unit | Feces | GI mucosa / lumen | Conditional OR (exact interval) | p |
| --- | --- | --- | --- | --- | --- |
| CRC, all_contrasts | experiment | 32 : 6 | 4 : 11 | 13.6 (2.9–81) | 1.2e-4 |
| CRC, all_contrasts | study | 19 : 2 | 3 : 8 | 21.6 (2.73–309) | 5.6e-4 |
| Pooled, stratified by condition, all_contrasts | experiment | 107 : 37 | 9 : 12 | 4.79 (1.52–15.8) | 0.0049 |
| Pooled, stratified by condition, all_contrasts | study | 56 : 21 | 5 : 9 | 7.87 (1.67–46.8) | 0.0047 |
| CRC, disease_vs_healthy_oriented | experiment | 22 : 5 | 2 : 5 | 10 (1.23–134) | 0.014 |
| CRC, disease_vs_healthy_oriented | study | 16 : 2 | 2 : 5 | 16.6 (1.57–302) | 0.0069 |
| Pooled, stratified by condition, disease_vs_healthy_oriented | study | 49 : 14 | 4 : 6 | 7.6 (1.31–52.2) | 0.010 |
| CRC, paired_tissue_counted, disease_vs_healthy_oriented | study | 16 : 2 | 3 : 4 | 9.35 (0.899–151) | 0.032 |
| CRC, control_tiers_1_2, disease_vs_healthy_oriented | study | 15 : 2 | 2 : 5 | 15.5 (1.46–284) | 0.0086 |
<!-- /numbers -->

Sources: `compartment_direction_interaction.tsv` and
`condition_definition_interaction.tsv`. Odds ratios are conditional
maximum-likelihood estimates with exact intervals. The TSVs also record sample
and Mantel–Haenszel estimates. The stratified rows test a common odds ratio
across conditions; crude pooling mixes condition with compartment.

Exact-test p-values count outcomes no more likely than the observed outcome,
whereas the intervals invert tail probabilities; their decisions can differ
near a chosen threshold. Paired-tissue and subject-control definitions remain
distinct, including the tumour/adjacent-tissue experiments `bsdb:22761885-2`
and `bsdb:25699023-2`.

<!-- numbers: interaction_shared_studies -->
The CRC disease-vs-healthy comparison includes 2 studies sampled in both compartments; the tests nevertheless treat units as independent.
<!-- /numbers -->

The direction differs between feces and CRC mucosa or lumen. With this sampling
effort and mixture of contrasts, that does not establish feces specificity.

## 6. Sensitivity to the set definition

<!-- numbers: sensitivity_table -->
| Definition | PD / HIV / CRC | Mean d_J | Core / union | U1 exact core p |
| --- | --- | --- | --- | --- |
| Feces dec, all cellulose-annotated ranks (primary) | 38 / 24 / 34 | 0.479 | 17 / 47 | 4.8e-23 |
| Feces dec, profiles with cellulose z > 0 only | 34 / 24 / 33 | 0.496 | 16 / 46 | 4.0e-22 |
| Feces dec, genus + species ranks only | 19 / 12 / 20 | 0.478 | 9 / 25 | 3.5e-12 |
| Feces net-decreased (dec profiles > inc profiles) | 23 / 19 / 28 | 0.641 | 7 / 40 | 5.6e-9 |
| Feces net-decreased (dec studies > inc studies) | 22 / 19 / 26 | 0.643 | 7 / 39 | 2.2e-9 |
| All human body sites dec, cellulose-annotated | 43 / 28 / 44 | 0.529 | 19 / 60 | 2.0e-21 |
| Baseline: all Feces dec matrix taxa (any trait) | 279 / 151 / 231 | 0.773 | 51 / 466 | 2.4e-103 |
| Feces dec, disease-vs-healthy contrasts only (oriented) | 33 / 21 / 31 | 0.512 | 14 / 43 | 8.7e-20 |
| Feces dec, disease-vs-healthy contrasts with strong comparison evidence | 33 / 19 / 31 | 0.508 | 14 / 42 | 8.9e-21 |
| Feces dec, without bsdb:24451087 (rectal mucosa curated as Feces) | 38 / 23 / 34 | 0.496 | 16 / 47 | 2.5e-21 |
| Feces dec, synonym condition labels counted under their disease | 38 / 27 / 37 | 0.456 | 20 / 50 | 1.8e-27 |
<!-- /numbers -->

Source: `sensitivity_definitions_jaccard.tsv`. The U1 probabilities assume
independent draws and therefore have the limitation described in section 2.
Restricting to genera and species leaves a similar mean distance.

## 7. Trait benchmark

<!-- numbers: benchmark_summary -->
The benchmark contains 186 traits, 146 excluding growth media. Cellulose ranks 2 among those traits and 23 overall. The non-media median mean d_J is 0.681. The smallest per-condition set size among all benchmarked traits is 10.
<!-- /numbers -->

<!-- numbers: benchmark_table -->
| Rank | Trait | ID | Mean d_J | Core / union |
| --- | --- | --- | --- | --- |
| 1 | Acetogenesis | `METPO:1000845` | 0.332 | 10 / 18 |
| 2 | cellulose catabolic process | `METPO:2000103\|GO:0030245` | 0.479 | 17 / 47 |
| 3 | chitin catabolic process | `METPO:2000103\|GO:0006032` | 0.510 | 8 / 24 |
| 4 | nitrate_reduction_to_ammonia | `METPO:2000103\|pathways:nitrate_reduction_to_ammonia` | 0.515 | 8 / 26 |
| 5 | aromatic_compound_degradation | `METPO:2000103\|pathways:aromatic_compound_degradation` | 0.517 | 13 / 39 |
| 6 | xylan catabolic process | `METPO:2000103\|GO:0045493` | 0.523 | 20 / 60 |
| 7 | spiral shaped | `METPO:1000684` | 0.525 | 10 / 32 |
| 8 | hydrogen_oxidation_dark | `METPO:2000103\|pathways:hydrogen_oxidation_dark` | 0.536 | 8 / 26 |
| 9 | sulfate_reduction | `METPO:2000103\|pathways:sulfate_reduction` | 0.542 | 6 / 19 |
| 10 | filament shaped | `METPO:1000674` | 0.548 | 11 / 37 |
<!-- /numbers -->

Source: `trait_convergence_benchmark.tsv`. Traits annotated on the same
frequently reported taxa rank alike here. The rank describes this corpus; it
does not test whether cellulose-annotated taxa are singled out.

## 8. Counting rules

The earlier catalog and literature file used different labels, not different
taxon-expansion rules. The present table reports the current committed bundle;
it does not mix those counts with unsupported earlier snapshots.

<!-- numbers: counting_table -->
| Count | Current pinned run |
| --- | --- |
| Studies with any decreased feces signature PD / HIV / CRC | 47 / 19 / 31 |
| Studies reporting an annotated taxon decreased PD / HIV / CRC | 39 / 16 / 22 |
| Roseburia decreased / increased studies PD; HIV; CRC | 20 / 4; 2 / 0; 6 / 0 |
<!-- /numbers -->

The rows distinguish any decreased feces signature (`n_studies_dec`) from a
signature containing an annotated taxon (`n_studies_a` in
`profile_and_study_pair_jaccard.tsv`). Taxon-specific counts come directly
from the catalog.

## 9. Files

| File | Contents |
| --- | --- |
| [`cellulose_decr_summary.json`](cellulose_decr_summary.json) | Results, definitions, universes, null summaries, label inventories and label-source comparison |
| [`cellulose_decr_provenance.json`](cellulose_decr_provenance.json) | Input hashes, classifier version, code commit, seeds and package versions |
| [`pairwise_jaccard_significance.tsv`](pairwise_jaccard_significance.tsv) | Pairwise, mean and core overlap under U1–U6 and the within-condition reference |
| [`sensitivity_definitions_jaccard.tsv`](sensitivity_definitions_jaccard.tsv) | Set-definition sensitivity |
| [`compartment_jaccard_and_direction_stats.tsv`](compartment_jaccard_and_direction_stats.tsv) | Direction tests by population, compartment and condition, with study collapses and power notes |
| [`compartment_direction_interaction.tsv`](compartment_direction_interaction.tsv) | Crude and stratified direction-by-compartment tests |
| [`condition_definition_sensitivity.tsv`](condition_definition_sensitivity.tsv), [`condition_definition_overlap.tsv`](condition_definition_overlap.tsv), [`condition_definition_interaction.tsv`](condition_definition_interaction.tsv) | Direction, overlap and compartment tests under condition definitions, with the case-control definitions restricted to disease-versus-healthy contrasts (#769, #781) |
| [`condition_label_inventory.tsv`](condition_label_inventory.tsv) | Human condition labels naming PD, HIV or CRC and their status |
| [`cellulose_decreased_taxa_catalog.tsv`](cellulose_decreased_taxa_catalog.tsv) | Taxon study/profile counts and annotation flags |
| [`cellulose_decreased_taxa_catalog_tier1_2.tsv`](cellulose_decreased_taxa_catalog_tier1_2.tsv) | Catalog on disease-versus-healthy contrasts with Tier 1–2 controls (#813) |
| [`profile_and_study_pair_jaccard.tsv`](profile_and_study_pair_jaccard.tsv) | Overlap between individual profiles and studies |
| [`trait_convergence_benchmark.tsv`](trait_convergence_benchmark.tsv) | Ranked traits, canonical labels and predicates |
| [`significant_results_summary.tsv`](significant_results_summary.tsv) | Headline tests, populations, p-values and BH q-values; see [`SIGNIFICANT_RESULTS_SUMMARY.md`](SIGNIFICANT_RESULTS_SUMMARY.md) |
| [`label_source_corrections.tsv`](label_source_corrections.tsv) | Experiments whose own labels differ from study-modal labels |
| [`four_highlighted_pairs_taxon_drivers.tsv`](four_highlighted_pairs_taxon_drivers.tsv), [`four_highlighted_pairs_study_evidence.tsv`](four_highlighted_pairs_study_evidence.tsv), [`four_highlighted_pairs_summary.json`](four_highlighted_pairs_summary.json) | Highlighted condition–trait pairs and their supporting taxa/studies |

<!-- numbers: tier_catalog -->
The committed Tier 1–2 catalog contains 41 taxa, 14 in the core, and 0 core taxa with a direct edge. These counts bind directly to the catalog and are cross-checked against its entries in the committed summary.
<!-- /numbers -->

Some outputs are written here but not committed: under CONTRIBUTING.md, images
stay out of git outside `trait-visualization/`.

- **Figures.** `just cellulose-decr-figures` writes presentation figures to
  `figures_graphical/` and statistical reference figures to `figures/`, as
  SVG, PDF and PNG. Their names describe the contents rather than a conclusion
  (#570): `fig1_overlap_and_null_models`,
  `fig2_shared_taxa_and_set_definitions`, `fig3_direction_by_body_site`,
  `fig4_trait_overlap_ranking`, `fig1_pairwise_and_3way_jaccard_nulls`,
  `fig2_direction_by_body_site_and_definition` and
  `fig3_core_taxa_and_trait_benchmark`.
- **Combined PDFs.** `cellulose_decr_graphical_figures.pdf` and
  `cellulose_decr_figures.pdf`.
- **Slide captions.** `figures/figure_captions.json` gives titles, takeaways,
  caveats and evidence tables. The selected-figure deck takes its cellulose
  slides from that file (#556).
- **Retired historical Tufte plates and deck.** The fixed study-modal display is
  retired because some numbers have no compatible tracked source or contradict
  the cited source
  ([#1572](https://github.com/CultureBotAI/microbial-trait-profiling/issues/1572)).
  It is excluded from default generation. For audit only,
  `just cellulose-decr-tufte --audit-retired` invokes
  `src/mtp/cellulose_decr/tufte_plates.py` to reproduce historical plates,
  supplemental SVGs, the combined `cellulose_decr_tufte_figures.pdf` and
  `cellulose_decr_slide_deck.html`. Every newly generated audit plate,
  standalone SVG and HTML slide carries the banner
  `RETIRED: historical, unverified numbers (#1572)`.
  Audit opt-in does not validate or recompute the historical numbers. Source
  snapshots and source-pin tests remain as evidence. Each slide has its own
  HTML banner even without referenced images. Inline fig8/fig9 omit the duplicate
  SVG banner band and keep their original size; referenced plates retain their
  own banners (#1738). A replacement is pending under
  [#1670](https://github.com/CultureBotAI/microbial-trait-profiling/issues/1670).
  The separate generator `scripts/generate_tufte_plates.py` is unaffected.
  Preview old unmarked outputs with `just cellulose-decr-cleanup-retired`;
  adding `--delete` removes only its declared historical files, preserving
  current figures.

Reruns of the current figures on the same tables give identical bytes.

## 10. Reproduce

```bash
just cellulose-decr                # tables, summary and provenance
just cellulose-decr-figures        # current table-driven figures only (gitignored)
```

The Monte Carlo nulls use fixed seeds, so reruns on the same inputs give identical
tables; the provenance timestamp changes.

**Inputs.** The pinned input inventory is in
`cellulose_decr_provenance.json`. Most inputs are gitignored. They are built by:

- `just download`, `just ingest` and `just study-tsvs` (BugSigDB);
- `just download-kg`, `just taxa-trait-matrix` and `just study-trait-profiles`
  (KG-Microbe);
- `just contrast-orientation` (contrast manifest).

**Drift check.** A run refuses changed input hashes, missing pinned provenance,
a different contrast classifier version, or inputs changed while being read.

From a worktree, use `--input-root <checkout>` to find the pinned inputs.
`--allow-input-drift` records changed inputs while retaining the old pins;
`--repin` adopts new pins. Neither option turns a changed input into the
historical run documented here.

The pins use phenotype-scope profiles and each experiment's own metadata labels
(#555). The committed tables reproduce at the provenance's `code_commit`
(#642). An earlier checkout may need `--pinned-provenance` pointing at the
committed provenance, or the commit that adds those pins (#674).

**Other labels.** Write alternative-label runs elsewhere, preserving the
committed tables:

- study-modal labels: `just cellulose-decr --label-source study_modal --output-dir <elsewhere>`;
- synonym labels: `just cellulose-decr --condition-labels synonyms --output-dir <elsewhere>`.
