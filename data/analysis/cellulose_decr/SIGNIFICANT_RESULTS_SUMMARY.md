# Results summary: feces-decreased `cellulose catabolic process` taxa in PD, HIV and CRC

This is a short account of which results hold up, which are fragile and which do
not hold. It covers `METPO:2000103|GO:0030245` (cellulose catabolic process) in
Parkinson's disease (PD), HIV infection (HIV) and colorectal cancer (CRC). The
machine-readable version is
[`significant_results_summary.tsv`](significant_results_summary.tsv). The full
account and definitions are in [`README.md`](README.md) and
[`methods/cellulose_decr.md`](../../../methods/cellulose_decr.md).

All results describe curated BugSigDB signature taxa and their KG-Microbe
annotations. None of them measures cellulose degradation or any other physiology
([`methods/biological_interpretation.md`](../../../methods/biological_interpretation.md)).
The trait and conditions were picked as exploratory leads from this corpus, so
no result here is confirmatory. The p-values are nominal; the TSV adds
Benjamini-Hochberg q-values over its 29 tests (`bh_q_within_table`). The numbers
are those of the tables regenerated in the final rerun at code commit 511c891
(#799).

## What holds up

1. **In PD and CRC, the decreased feces signature carries a larger cellulose
   annotation share than the paired increased signature.** The same holds
   across the three conditions pooled. In most studies, taxa annotated with
   cellulose catabolic process form a larger share of the decreased signature's
   annotations than of the increased signature's. The share counts annotations,
   not taxa.
   - **CRC:** 19 of 21 untied studies (p = 2.2e-4). Among disease-versus-healthy
     contrasts it is 16 of 18 (p = 0.0013, q = 0.0035), and 18 of 21 with synonym
     labels counted (p = 0.0015). Keeping only Tier 1–2 controls sets aside the
     Yachida 2019 study, whose controls had "a few polyps" (#672), giving 15 of 17
     (p = 0.0023).
   - **PD, disease vs healthy:** 23 of 33 studies (p = 0.035, q = 0.054). This q
     is just above 0.05; before the final rerun it was 0.037.
   - **PD, all contrasts:** the study-level sign test gives 26 of 40 (p = 0.081).
     The study-level Wilcoxon test gives p = 0.013, and a majority vote of
     experiment signs gives 27 of 39 (p = 0.024).
   - **Pooled:** 56 of 77 untied studies (p = 8.2e-5) and 49 of 63 among
     disease-versus-healthy contrasts (p = 1.1e-5). The pooled result stays at
     p ≤ 1.3e-4 under every label and contrast definition tried
     ([`condition_definition_sensitivity.tsv`](condition_definition_sensitivity.tsv)).

## What is fragile

1. **HIV.** Counting the exact label "HIV infection" and disease-versus-healthy
   contrasts, 10 of 12 untied studies point down (p = 0.039, q = 0.056). That is
   one of several defensible readings, and across them the study-level p runs from
   0.002 to 0.18:
   - Since contrast classifier v4, this reading counts the HIV– vs HIV+ experiment
     and the ART-treated HIV+ vs uninfected experiments as disease vs healthy
     (#569, #580).
   - The two studies that point up, `bsdb:38383483` and `bsdb:39880958`, compare
     people with HIV on antiretroviral therapy against people without HIV. Setting
     aside contrasts whose cases are named as treated gives 10 down and none up
     (p = 0.0020). Keeping only Tier 1–2 controls gives the same.
   - The treated flag reads group names only. `bsdb:32415070-1`, whose group
     definition names people on antiretroviral therapy, is not flagged and points
     down. The readings separate contrasts by their wording, not by treatment
     (#676).
   - Two studies that compare people with and without HIV are filed as "Human
     immunodeficiency virus". Counting them gives 10 down and 4 up (p = 0.18).

   Table A2 lists every reading (#561, #569, #580).
2. **The direction in CRC mucosa and lumen differs from feces, but the samples
   differ too.** In feces, 32 of 38 untied experiments point down; in mucosa or lumen,
   4 of 15 do. The conditional odds ratio is 13.6 (95% CI 2.9–81.0,
   p = 1.2e-4). However:
   - Only 11 of the 31 paired mucosa and lumen experiments compare patients with
     controls (`contrast_type` `case_vs_control_subject`). The rest compare tumour
     with adjacent tissue, cancer subgroups or treatments.
   - Among disease-versus-healthy studies, the test gives odds ratio 16.6
     (95% CI 1.57–302, p = 0.0069). Counting the tumour-vs-adjacent-tissue contrast
     `bsdb:25699023-2` gives 9.35 (0.90–151, p = 0.032). Keeping only Tier 1–2
     controls gives 15.5 (1.46–284, p = 0.0086). The p stays below 0.05 either
     way, on 7 untied GI studies (#569, #580, #672, #673).

## What does not hold

1. **The three conditions do not share detectably more cellulose-annotated taxa
   than random sets of feces studies.**
   - Drawing whole study signatures (U6) gives an expected mean pairwise Jaccard
     distance of 0.527, against 0.479 observed. That is slightly more overlap
     than expected, but not significantly (p = 0.19).
   - U6's larger expected core (18.3 taxa) comes from its larger simulated sets.
     With set sizes matched it is 16.9 against 17 observed (p = 0.59).
   - The very small p-values under U1–U4 come from treating co-reported taxa as
     independent draws. Weighting taxa by how often other studies report them
     (U5) already weakens them to p = 0.0086 (mean) and 0.072 (core).
2. **The overlap is not a "same bacteria lost" result.** A taxon enters a set
   once any study reports it decreased. The 17 taxa in all three sets have these
   features:
   - none has its own annotation edge; all of them, the nine genera included,
     inherit it from annotated species or strains;
   - only seven are reported decreased more often than increased in every
     condition;
   - *Eubacterium* is in the core only because of one HIV study that sampled
     rectal mucosa but is curated as feces (#566).
3. **Specificity to feces is not shown.**
   - Outside feces, single conditions have at most five untied pairs, too few for
     p < 0.05, except CRC mucosa and lumen.
   - The CRC contrasts there are mostly not patient-versus-control, and the
     disease-vs-healthy comparison has 7 untied GI studies (fragile #2).
4. **The trait ranking is descriptive.** Cellulose catabolic process ranks second
   of 146 non-MediaDive traits by cross-condition distance, after Acetogenesis
   (`METPO:1000845`). The ranking is not a test, and traits share taxa.

## Table A. Direction in feces

The table gives down : up : tied counts. "Down" means the decreased signature had
the larger cellulose annotation share. Sign tests are exact and two-sided; the
study-level Wilcoxon test is on study means. The q-values (BH over the 29 tests
of the TSV) are for the study sign test.

| Condition | Contrasts | Experiments | p | Studies | p | q | Study Wilcoxon p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| PD | all | 57 : 24 : 16 | 3.2e-4 | 26 : 14 : 3 | 0.081 | 0.10 | 0.013 |
| PD | disease vs healthy | 46 : 15 : 11 | 8.8e-5 | 23 : 10 : 3 | 0.035 | 0.054 | 0.019 |
| HIV | all | 18 : 7 : 3 | 0.043 | 11 : 5 : 2 | 0.21 | 0.23 | 0.32 |
| HIV | disease vs healthy | 13 : 3 : 3 | 0.021 | 10 : 2 : 2 | 0.039 | 0.056 | 0.18 |
| CRC | all | 32 : 6 : 11 | 2.4e-5 | 19 : 2 : 7 | 2.2e-4 | 8.0e-4 | 2.0e-4 |
| CRC | disease vs healthy | 22 : 5 : 8 | 0.0015 | 16 : 2 : 7 | 0.0013 | 0.0035 | 0.0010 |
| All three | all | 107 : 37 : 30 | 4.4e-9 | 56 : 21 : 12 | 8.2e-5 | 4.3e-4 | 5.6e-5 |
| All three | disease vs healthy | 81 : 23 : 22 | 9.3e-9 | 49 : 14 : 12 | 1.1e-5 | 1.1e-4 | 8.8e-5 |

Experiment-level p-values treat experiments from one study as independent, so
they are optimistic. The majority-vote study counts are in
`compartment_jaccard_and_direction_stats.tsv` (`study_majority_*`).

## Table A2. HIV feces, disease vs healthy, under each reading

| Reading | Experiments | p | Studies | p |
| --- | --- | --- | --- | --- |
| Exact label (primary) | 13 : 3 : 3 | 0.021 | 10 : 2 : 2 | 0.039 |
| + synonym labels | 13 : 5 : 3 | 0.096 | 10 : 4 : 2 | 0.18 |
| Contrasts whose cases are named as treated set aside | 12 : 1 : 3 | 0.0034 | 10 : 0 : 2 | 0.0020 |
| Tier 1–2 controls only | 12 : 1 : 3 | 0.0034 | 10 : 0 : 2 | 0.0020 |
| Paired-tissue contrasts also counted | 13 : 3 : 3 | 0.021 | 10 : 2 : 2 | 0.039 |
| Comorbidity contrasts also counted | 13 : 3 : 3 | 0.021 | 10 : 2 : 2 | 0.039 |
| Exact label, without the rectal-mucosa study | 12 : 3 : 3 | 0.035 | 9 : 2 : 2 | 0.065 |

## Table B. Overlap against the null models

The table gives p-values for the two summary endpoints. For the 3-way core,
U1–U3 give exact hypergeometric tails. All other p-values are Monte Carlo with
100,000 draws; the size-matched U6 keeps 5,557 of them.

| Null | Mean pairwise d_J (0.479) | Taxa in all three (17) |
| --- | --- | --- |
| U1–U3: uniform draws of taxa | ≤ 1e-5 | 4.8e-23 to 2.3e-14 |
| U4: random trait on the decreased feces taxa | ≤ 1e-5 | ≤ 1e-5 |
| U5: U3 taxa weighted by reporting prevalence | 0.0086 | 0.072 |
| U6: whole study signatures | 0.19 | 0.73 |
| U6, draws with matching set sizes | 0.23 | 0.59 |
| U6 with synonym labels counted (sets 38 / 27 / 37, core 20) | 0.14 | 0.54 |

## Table C. The 17 taxa in all three sets

| Group | Taxa | Mean d_J if only this kind of taxon is kept |
| --- | --- | --- |
| Genus rank (9) | *Roseburia*, *Blautia*, *Clostridium*, *Bacteroides*, *Lachnospira*, *Ruminococcus*, *Anaerostipes*, *Eubacterium*, *Butyrivibrio* | 0.478 with genus and species ranks only: 9 of 25 taxa in all three, 17.1 times the U1 expectation, exact p = 3.5e-12 |
| Family-to-phylum rank (8) | Lachnospiraceae, Clostridiaceae, Oscillospiraceae (families); Eubacteriales (order); Clostridia, Bacteroidia (classes); Bacillota, Bacteroidota (phyla) | |
| More decreased than increased studies in all three conditions (7) | *Roseburia*, *Clostridium*, *Lachnospira*, *Ruminococcus*, *Anaerostipes*, Lachnospiraceae, Bacillota | 0.643 if sets keep only taxa with more decreased than increased studies: 7 of 39 in all three, 23.0 times the U1 expectation, exact p = 2.2e-9 (counting profiles: 0.641, 7 of 40, p = 5.6e-9) |
| Reported increased at least as often as decreased in one condition | *Blautia* (HIV 1 / 3 studies), *Bacteroides* (CRC 2 / 2), *Eubacterium*, *Butyrivibrio*, Oscillospiraceae (PD 8 / 13), Clostridia (PD 5 / 7), Eubacteriales (PD 5 / 6), Clostridiaceae (CRC 2 / 2), Bacteroidia (CRC 1 / 4), Bacteroidota (CRC 1 / 4) | |

None of the 17 has its own annotation edge (`has_direct_trait_edge`). Every
annotation in the table is inherited from annotated descendants, genera included.
The U1 expectations assume independent draws; Table B shows how much that
overstates significance. Per-taxon counts are in
[`cellulose_decreased_taxa_catalog.tsv`](cellulose_decreased_taxa_catalog.tsv).
The flags are `has_direct_trait_edge`, `is_direct_or_genus_species`,
`is_strictly_net_decreased_all3`, `is_net_decreased_all3_by_profiles` and
`core_class`.

## Table D. Feces against GI mucosa or lumen

| Comparison | Units | Feces down : up | GI down : up | Conditional odds ratio (95% CI) | p |
| --- | --- | --- | --- | --- | --- |
| CRC | experiments | 32 : 6 | 4 : 11 | 13.6 (2.9–81.0) | 1.2e-4 |
| CRC | studies | 19 : 2 | 3 : 8 | 21.6 (2.7–309) | 5.6e-4 |
| CRC, disease vs healthy | studies | 16 : 2 | 2 : 5 | 16.6 (1.57–302) | 0.0069 |
| CRC, disease vs healthy, paired-tissue contrasts counted | studies | 16 : 2 | 3 : 4 | 9.35 (0.90–151) | 0.032 |
| All three, stratified by condition | experiments | 107 : 37 | 9 : 12 | 4.79 (1.52–15.8) | 0.0049 |
| All three, stratified by condition | studies | 56 : 21 | 5 : 9 | 7.87 (1.67–46.8) | 0.0047 |

The stratified rows test a common odds ratio across PD, HIV and CRC. Crude
pooling mixes condition with compartment, because pooled stool pairs are mostly
PD and pooled GI pairs mostly CRC. PD and HIV have one and five untied GI pairs,
too few to test. Oral and other non-GI strata are in
[`compartment_direction_interaction.tsv`](compartment_direction_interaction.tsv).
