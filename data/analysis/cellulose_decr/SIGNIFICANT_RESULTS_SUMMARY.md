# Results summary: feces-decreased `cellulose catabolic process` taxa in PD, HIV and CRC

A short account of which results hold up and which do not, for
`METPO:2000103|GO:0030245` (cellulose catabolic process) in Parkinson's disease
(PD), HIV infection (HIV) and colorectal cancer (CRC). The machine-readable
version is [`significant_results_summary.tsv`](significant_results_summary.tsv);
the full account and definitions are in [`README.md`](README.md) and
[`methods/cellulose_decr.md`](../../../methods/cellulose_decr.md).

All results describe curated BugSigDB signature taxa and their KG-Microbe
annotations, not measured cellulose degradation or any other physiology
([`methods/biological_interpretation.md`](../../../methods/biological_interpretation.md)).
The trait and conditions were picked as exploratory leads from this corpus, so
no result here is confirmatory.

## What holds up

1. **Within each condition, decreased feces signatures carry a larger cellulose
   annotation share than the paired increased signatures.** The result is
   clearest in CRC (19 of 21 untied studies, p = 2.2e-4) and across the three
   conditions pooled (58 of 78 untied studies, p = 2.0e-5). For PD and HIV it
   reaches p < 0.05 at the study level only among disease-versus-healthy
   contrasts (PD 24 of 34 studies, p = 0.024; HIV 9 of 9, p = 0.0039).
2. **The direction differs between feces and CRC mucosa or lumen**: 32 of 38
   untied experiments in feces against 4 of 15 in mucosa or lumen (conditional
   odds ratio 13.6, 95% CI 2.9–81.0, p = 1.2e-4; at the study level 21.6,
   2.7–309, p = 5.6e-4).

## What does not

1. **The three conditions do not share more cellulose-annotated taxa than random
   sets of feces studies.** Drawing whole study signatures (U6) gives an expected
   mean pairwise Jaccard distance of 0.527 against 0.479 observed (p = 0.19) and
   an expected 18.3 taxa in all three against 17 observed (p = 0.73). The very
   small p-values under U1–U4 come from treating co-reported taxa as independent
   draws.
2. **The overlap is not a "same bacteria lost" result.** A taxon enters a set
   once any study reports it decreased. Of the 17 taxa in all three sets, eight
   are family-to-phylum nodes with rolled-up annotations and only seven are
   reported decreased more often than increased in every condition.
3. **Specificity to feces is not shown.** Outside feces most strata have fewer
   than six untied pairs, so p < 0.05 is often unattainable, and the CRC mucosa
   and lumen experiments mix tumour-versus-adjacent and subgroup contrasts with
   patient-versus-control ones.
4. **The trait ranking is descriptive.** Cellulose catabolic process ranks second
   of 142 non-MediaDive traits by cross-condition distance, after Acetogenesis
   (`METPO:1000845`). The ranking is not a test and traits share taxa.

## Table A. Direction in feces

Down : up : tied counts; down means the decreased signature had the larger
cellulose annotation share. Sign tests are exact and two-sided; the study-level
Wilcoxon test is on study means.

| Condition | Contrasts | Experiments | p | Studies | p | Study Wilcoxon p |
| --- | --- | --- | --- | --- | --- | --- |
| PD | all | 58 : 24 : 17 | 2.2e-4 | 27 : 14 : 4 | 0.060 | 0.0093 |
| PD | disease vs healthy | 47 : 15 : 10 | 5.8e-5 | 24 : 10 : 4 | 0.024 | 0.015 |
| HIV | all | 18 : 7 : 3 | 0.043 | 12 : 4 : 2 | 0.077 | 0.30 |
| HIV | disease vs healthy | 11 : 1 : 3 | 0.0063 | 9 : 0 : 2 | 0.0039 | 0.0039 |
| CRC | all | 32 : 6 : 11 | 2.4e-5 | 19 : 2 : 7 | 2.2e-4 | 2.0e-4 |
| CRC | disease vs healthy | 21 : 4 : 8 | 9.1e-4 | 14 : 3 : 7 | 0.013 | 0.0021 |
| All three | all | 108 : 37 : 31 | 3.0e-9 | 58 : 20 : 13 | 2.0e-5 | 4.1e-5 |
| All three | disease vs healthy | 79 : 20 : 21 | 1.8e-9 | 47 : 13 : 13 | 1.2e-5 | 2.0e-5 |

Experiment-level p-values treat experiments from one study as independent and
are optimistic.

## Table B. Overlap against the null models

p-values for the two summary endpoints. U1–U3 give exact hypergeometric tails for
the 3-way core; all other p-values are Monte Carlo with 100,000 draws.

| Null | Mean pairwise d_J (0.479) | Taxa in all three (17) |
| --- | --- | --- |
| U1–U3: uniform draws of taxa | ≤ 1e-5 | 4.8e-23 to 2.3e-14 |
| U4: random trait on the decreased feces taxa | ≤ 1e-5 | ≤ 1e-5 |
| U5: taxa weighted by reporting prevalence | 0.0048 | 0.053 |
| U6: whole study signatures | 0.19 | 0.73 |
| U6, draws with matching set sizes | 0.24 | 0.58 |

## Table C. The 17 taxa in all three sets

| Group | Taxa | Mean d_J if only this kind of taxon is kept |
| --- | --- | --- |
| Genera (9) | *Roseburia*, *Blautia*, *Clostridium*, *Bacteroides*, *Lachnospira*, *Ruminococcus*, *Anaerostipes*, *Eubacterium*, *Butyrivibrio* | 0.478 with genus and species ranks only: 9 of 25 taxa in all three, 17.1 times the U1 expectation, exact p = 3.5e-12 |
| Higher-rank nodes with rolled-up annotations (8) | Lachnospiraceae, Clostridiaceae, Oscillospiraceae (families); Eubacteriales (order); Clostridia, Bacteroidia (classes); Bacillota, Bacteroidota (phyla) | |
| More decreased than increased studies in all three conditions (7) | *Roseburia*, *Clostridium*, *Lachnospira*, *Ruminococcus*, *Anaerostipes*, Lachnospiraceae, Bacillota | 0.643 if sets keep only taxa with more decreased than increased studies: 7 of 39 in all three, 23.0 times the U1 expectation, exact p = 2.2e-9 (counting profiles: 0.641, 7 of 40, p = 5.6e-9) |
| Reported increased at least as often as decreased in one condition | *Blautia* (HIV 1 / 3 studies), *Bacteroides* (CRC 2 / 2), *Eubacterium*, *Butyrivibrio*, Oscillospiraceae (PD 8 / 13), Clostridia (PD 5 / 7), Eubacteriales (PD 5 / 6), Clostridiaceae (CRC 2 / 2), Bacteroidia (CRC 1 / 4), Bacteroidota (CRC 1 / 4) | |

The U1 expectations in this table assume independent draws; Table B shows how
much that overstates significance. Per-taxon counts are in
[`cellulose_decreased_taxa_catalog.tsv`](cellulose_decreased_taxa_catalog.tsv)
(`is_direct_or_genus_species`, `is_strictly_net_decreased_all3`,
`is_net_decreased_all3_by_profiles`, `core_class`).

## Table D. Feces against GI mucosa or lumen

| Comparison | Units | Feces down : up | GI down : up | Conditional odds ratio (95% CI) | p |
| --- | --- | --- | --- | --- | --- |
| CRC | experiments | 32 : 6 | 4 : 11 | 13.6 (2.9–81.0) | 1.2e-4 |
| CRC | studies | 19 : 2 | 3 : 8 | 21.6 (2.7–309) | 5.6e-4 |
| All three | experiments | 108 : 37 | 9 : 12 | 3.85 (1.37–11.3) | 0.0048 |
| All three | studies | 58 : 20 | 5 : 9 | 5.11 (1.35–21.9) | 0.0096 |

PD and HIV have one and five untied GI pairs, too few to test. Oral and other
non-GI strata are in
[`compartment_direction_interaction.tsv`](compartment_direction_interaction.tsv).
