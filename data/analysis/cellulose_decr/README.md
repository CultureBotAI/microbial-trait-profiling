# Cellulose-annotated taxa reported decreased in feces in Parkinson's disease, HIV infection and colorectal cancer

This directory compares the taxa annotated `METPO:2000103|GO:0030245`
(**cellulose catabolic process**) that curated BugSigDB signatures report decreased
in human feces in Parkinson's disease (PD), HIV infection (HIV) and colorectal
cancer (CRC). Every number below comes from one run of
[`scripts/cellulose_decr_analysis.py`](scripts/cellulose_decr_analysis.py)
(provenance in [`cellulose_decr_provenance.json`](cellulose_decr_provenance.json));
the method is described in [`methods/cellulose_decr.md`](../../../methods/cellulose_decr.md).

These are annotation counts of curated signature taxa. They do not measure
cellulose degradation, fiber intake or short-chain fatty acids, and a taxon
counted here was reported decreased in at least one study, not consistently lost
([`methods/biological_interpretation.md`](../../../methods/biological_interpretation.md)).

## What changed on 2026-10-17 (#552)

The October 16 version of this directory had no generator. The rebuilt analysis
reproduces it exactly when given the same labels, and then corrects it:

- **Labels.** Each profile now takes the condition and body site of its own
  experiment. The October 16 run used one modal label per study, which put
  saliva and nasal experiments into "PD feces" and moved other experiments out.
  61 experiments in these three conditions changed
  ([`label_source_corrections.tsv`](label_source_corrections.tsv)).
- **Trait names.** The benchmark carries `trait_label` from
  `data/derived/taxa_trait_matrix_columns.tsv`. `METPO:1000845` is
  **Acetogenesis**, `METPO:1000679` **pleomorphic shaped**, `METPO:1000684`
  **spiral shaped**, `METPO:1000616` **thermophilic** and `METPO:1000674`
  **filament shaped**; earlier figures and slides labelled them as spore-forming,
  motile, rod-shaped, cell envelope or obligately anaerobic traits.
  `METPO:1000682` (spore shaped) is not a matrix column.
- **A whole-signature null (U6).** Drawing whole study signatures keeps taxa that
  studies report together, and it removes the evidence that the three conditions
  share more taxa than random sets of feces studies do (section 2).
- **Study-level and disease-versus-healthy tests** of direction (section 4), a
  compartment interaction test with power notes (section 5), catalog flags
  separating genera from higher-rank nodes and net-decreased taxa (section 3),
  reconciled counts (section 8) and a corrected literature file.

| Quantity | Oct 16 (study labels) | Now (experiment labels) |
| --- | --- | --- |
| Decreased cellulose-annotated taxa, PD / HIV / CRC | 40 / 24 / 34 | 38 / 24 / 34 |
| Mean pairwise d_J; taxa in all three / union | 0.4957; 17 / 49 | 0.4793; 17 / 47 |
| PD paired experiments, down : up : tied | 57 : 26 : 20 | 58 : 24 : 17 |
| HIV paired experiments; Wilcoxon p | 19 : 6 : 3; 0.037 | 18 : 7 : 3; 0.135 |
| CRC paired experiments | 30 : 7 : 10 | 32 : 6 : 11 |
| Non-MediaDive traits in the benchmark | 141 (reported as 90) | 142 |

## 1. Overlap of the three taxon sets

A condition's set is every matrix taxon annotated with the trait, at any rank,
that at least one of its feces decreased signatures reports.

| Pair | Set sizes | Shared | Union | J | d_J |
| --- | --- | --- | --- | --- | --- |
| PD vs. HIV | 38, 24 | 21 | 41 | 0.512 | 0.488 |
| PD vs. CRC | 38, 34 | 27 | 45 | 0.600 | 0.400 |
| HIV vs. CRC | 24, 34 | 18 | 40 | 0.450 | 0.550 |
| Mean of the three pairs | | | | 0.521 | **0.479** |
| All three | | **17** | **47** | | |

32 of the 47 taxa appear in at least two sets
([`pairwise_jaccard_significance.tsv`](pairwise_jaccard_significance.tsv),
[`cellulose_decreased_taxa_catalog.tsv`](cellulose_decreased_taxa_catalog.tsv)).

## 2. Is the overlap unusual?

Expected values under each null, with the range holding 95% of the draws.
Monte Carlo p-values use 100,000 draws; "≤ 1e-5" means no draw was as extreme.

| Null | Expected mean d_J (95% range) | p | Expected taxa in all three (95% range) | p |
| --- | --- | --- | --- | --- |
| U1: uniform, 189 annotated matrix taxa | 0.909 (0.865–0.946) | ≤ 1e-5 | 0.87; exact p = 4.8e-23 | |
| U2: uniform, 135 annotated taxa seen in feces | 0.868 (0.817–0.911) | ≤ 1e-5 | 1.70; exact p = 7.5e-18 | |
| U3: uniform, 108 annotated taxa ever decreased in feces | 0.829 (0.772–0.878) | ≤ 1e-5 | 2.66; exact p = 2.3e-14 | |
| U4: random trait on the 471 decreased feces taxa of the three conditions | 0.775 (0.669–0.875) | ≤ 1e-5 | 5.09 (1–9) | ≤ 1e-5 |
| U5: weighted by reporting in other feces studies | 0.596 (0.513–0.677) | 0.0048 | 13.4 (9–17) | 0.053 |
| **U6: whole study signatures** | **0.527 (0.420–0.636)** | **0.19** | **18.3 (12–24)** | **0.73** |
| U6, draws with matching set sizes (5,284) | 0.513 (0.420–0.606) | 0.24 | 16.9 (13–21) | 0.58 |

The nulls that treat taxa as independent draws (U1–U4) make the overlap look
extreme. Weighting taxa by how often other feces studies report them decreased
(U5) leaves a modest signal for the mean distance (p = 0.0048) and none for the
3-way core (p = 0.053). Drawing 49, 19 and 31 whole study signatures at random
from all 646 feces studies with a decreased signature (U6) produces as much
overlap as observed: the expected mean d_J is 0.527 against 0.479 observed, and
the expected core is larger than the observed 17. The pairwise U6 p-values are
0.21 (PD–HIV), 0.10 (PD–CRC) and 0.55 (HIV–CRC).

So the shared taxa are mostly taxa that many feces studies report decreased,
together, whatever the condition. The data give no evidence of a
condition-specific convergence beyond that.

## 3. The 17 taxa in all three sets

Stool studies reporting each taxon decreased / increased
([`cellulose_decreased_taxa_catalog.tsv`](cellulose_decreased_taxa_catalog.tsv)).
"Net decreased" means more decreased than increased studies in all three
conditions (`is_strictly_net_decreased_all3`).

| Taxon | Rank | PD | HIV | CRC | Net decreased |
| --- | --- | --- | --- | --- | --- |
| *Roseburia* | genus | 20 / 4 | 2 / 0 | 6 / 0 | yes |
| *Blautia* | genus | 13 / 2 | 1 / 3 | 6 / 1 | no |
| *Clostridium* | genus | 7 / 3 | 5 / 2 | 6 / 3 | yes |
| *Bacteroides* | genus | 6 / 3 | 8 / 0 | 2 / 2 | no |
| *Lachnospira* | genus | 5 / 2 | 3 / 0 | 6 / 1 | yes |
| *Ruminococcus* | genus | 6 / 2 | 4 / 0 | 4 / 3 | yes |
| *Anaerostipes* | genus | 2 / 1 | 3 / 0 | 3 / 1 | yes |
| *Eubacterium* | genus | 1 / 1 | 1 / 1 | 3 / 2 | no |
| *Butyrivibrio* | genus | 2 / 2 | 1 / 2 | 1 / 0 | no |
| Lachnospiraceae | family | 20 / 7 | 5 / 3 | 9 / 5 | yes |
| Bacillota | phylum | 8 / 5 | 4 / 2 | 4 / 1 | yes |
| Clostridia | class | 5 / 7 | 4 / 0 | 7 / 1 | no |
| Oscillospiraceae | family | 8 / 13 | 4 / 0 | 2 / 1 | no |
| Clostridiaceae | family | 8 / 2 | 3 / 0 | 2 / 2 | no |
| Eubacteriales | order | 5 / 6 | 4 / 1 | 4 / 1 | no |
| Bacteroidia | class | 5 / 2 | 2 / 0 | 1 / 4 | no |
| Bacteroidota | phylum | 3 / 0 | 3 / 1 | 1 / 4 | no |

Nine of the 17 are genera (`is_direct_or_genus_species`, which is rank-based;
the matrix does not record whether an annotation is direct). The other eight are
family, order, class and phylum nodes whose annotation is rolled up from
descendants; `matrix_trait_bearing_taxa_at_or_below` gives how many annotated
taxa sit beneath each one (17 under Lachnospiraceae, 86 under Bacillota).
Seven taxa are net decreased in all three conditions, five of them genera. Several
others are reported increased more often than decreased in one condition.

Counting decreased minus increased *profiles* instead of studies (sensitivity
definition 4, `is_net_decreased_all3_by_profiles`) also gives seven taxa, with
*Bacteroides* in place of *Anaerostipes*.

## 4. Direction within each condition

Paired experiments in feces: down means the decreased signature had the larger
cellulose annotation share. Study rows give each study one vote, the sign of its
mean paired difference
([`compartment_jaccard_and_direction_stats.tsv`](compartment_jaccard_and_direction_stats.tsv)).

| Condition | Contrasts | Experiments down : up : tied | Sign test p | Studies down : up : tied | Sign test p |
| --- | --- | --- | --- | --- | --- |
| PD | all | 58 : 24 : 17 | 2.2e-4 | 27 : 14 : 4 | 0.060 |
| PD | disease vs healthy | 47 : 15 : 10 | 5.8e-5 | 24 : 10 : 4 | 0.024 |
| HIV | all | 18 : 7 : 3 | 0.043 | 12 : 4 : 2 | 0.077 |
| HIV | disease vs healthy | 11 : 1 : 3 | 0.0063 | 9 : 0 : 2 | 0.0039 |
| CRC | all | 32 : 6 : 11 | 2.4e-5 | 19 : 2 : 7 | 2.2e-4 |
| CRC | disease vs healthy | 21 : 4 : 8 | 9.1e-4 | 14 : 3 : 7 | 0.013 |
| All three | all | 108 : 37 : 31 | 3.0e-9 | 58 : 20 : 13 | 2.0e-5 |
| All three | disease vs healthy | 79 : 20 : 21 | 1.8e-9 | 47 : 13 : 13 | 1.2e-5 |

Collapsing experiments to studies weakens PD and HIV, as expected when several
experiments come from one cohort. Keeping only disease-versus-healthy contrasts
strengthens every condition; no experiment in the three conditions was
classified as reversed. Dropping experiments with weak comparison evidence
changes only HIV (11 : 1 : 1 experiments, 9 : 0 : 0 studies). One HIV feces
experiment mentions MSM in its group names; dropping it leaves 10 : 1 : 3
(p = 0.012). That lexical filter cannot separate HIV from sexual behaviour, which
Noguera-Julian et al. (2016, PMID 27077120) found drives much of the reported
HIV-associated shift.

The corpus-wide disease-versus-healthy comparison in
`data/derived/contrast_orientation/disease_control_comparison.tsv`, which pools
all conditions and is a separate analysis, moves the same way: cellulose
annotation share is lower in the decreased signature in 441 of 774 untied studies
(p = 1.2e-4) overall and 278 of 457 (p = 4.2e-6) for oriented
disease-versus-healthy contrasts.

## 5. Outside feces

| Compartment | Decreased-signature studies, PD / HIV / CRC | Experiments down : up (all three) | Note |
| --- | --- | --- | --- |
| Feces | 49 / 19 / 31 | 108 : 37 | |
| Non-fecal GI (mucosa, lumen) | 2 / 4 / 14 | 9 : 12 | PD has 1 untied pair, HIV 5: p < 0.05 is unattainable |
| Oral / oropharyngeal | 2 / 3 / 4 | 3 : 5 | p < 0.05 unattainable in every condition |
| Other non-GI | 3 / 3 / 0 | 4 : 2 | no CRC study |

Outside feces there are too few untied pairs to test anything except CRC
mucosa and lumen, where the decreased signature had the larger share in 4 of 15
untied experiments. Fisher's exact test of direction by compartment
([`compartment_direction_interaction.tsv`](compartment_direction_interaction.tsv)):

| Comparison | Units | Feces | GI mucosa / lumen | Odds ratio (exact 95% CI) | p |
| --- | --- | --- | --- | --- | --- |
| CRC | experiments | 32 : 6 | 4 : 11 | 13.6 (2.9–81.0) | 1.2e-4 |
| CRC | studies | 19 : 2 | 3 : 8 | 21.6 (2.7–309) | 5.6e-4 |
| All three | experiments | 108 : 37 | 9 : 12 | 3.85 (1.37–11.3) | 0.0048 |
| All three | studies | 58 : 20 | 5 : 9 | 5.11 (1.35–21.9) | 0.0096 |
| CRC, disease vs healthy | experiments | 21 : 4 | 3 : 5 | 8.0 (1.1–75.9) | 0.020 |
| CRC, disease vs healthy | studies | 14 : 3 | 3 : 4 | 5.7 (0.6–65.5) | 0.13 |

The odds ratios are conditional maximum-likelihood estimates; the sample odds
ratios are in the TSV (14.7 for CRC experiments). Only 16 of the 37 CRC mucosa
and lumen experiments compare patients with controls; the rest compare tumour
with adjacent tissue, cancer subgroups or aspirin with placebo. Two CRC studies
sampled both compartments, and the units are treated as independent. The
direction differs between feces and CRC mucosa or lumen; with this sampling
effort and mix of contrasts, that does not show the pattern is specific to feces.

## 6. Sensitivity to the set definition

| Definition | PD / HIV / CRC | Mean d_J | All three / union | U1 exact p, 3-way |
| --- | --- | --- | --- | --- |
| Primary: all ranks | 38 / 24 / 34 | 0.479 | 17 / 47 | 4.8e-23 |
| Genus and species only | 19 / 12 / 20 | 0.478 | 9 / 25 | 3.5e-12 |
| Profiles with cellulose z > 0 | 34 / 24 / 33 | 0.496 | 16 / 46 | 4.0e-22 |
| Disease-vs-healthy contrasts | 33 / 19 / 29 | 0.517 | 13 / 41 | 5.8e-19 |
| All body sites | 43 / 28 / 44 | 0.529 | 19 / 60 | 2.0e-21 |
| Net decreased, profiles | 23 / 19 / 28 | 0.641 | 7 / 40 | 5.6e-9 |
| Net decreased, studies | 22 / 19 / 26 | 0.643 | 7 / 39 | 2.2e-9 |
| Every decreased feces taxon, any trait | 284 / 151 / 231 | 0.775 | 51 / 471 | 6.5e-103 |

The U1 p-values in this table assume independent draws and overstate
significance in the way section 2 shows
([`sensitivity_definitions_jaccard.tsv`](sensitivity_definitions_jaccard.tsv)).
Restricting to genera and species leaves the mean distance unchanged.

## 7. Trait benchmark

Of the 182 matrix traits with at least 10 decreased feces taxa in each condition,
142 are not MediaDive growth media. Ranked by mean pairwise d_J, cellulose
catabolic process is second of the 142 (23rd of all 182); the rank describes
this corpus and is not a test
([`trait_convergence_benchmark.tsv`](trait_convergence_benchmark.tsv), median
d_J 0.679).

| Rank | Trait | ID | Mean d_J | All three / union |
| --- | --- | --- | --- | --- |
| 1 | Acetogenesis | `METPO:1000845` | 0.332 | 10 / 18 |
| 2 | cellulose catabolic process | `METPO:2000103\|GO:0030245` | 0.479 | 17 / 47 |
| 3 | chitin catabolic process | `METPO:2000103\|GO:0006032` | 0.510 | 8 / 24 |
| 4 | nitrate_reduction_to_ammonia | `METPO:2000103\|pathways:nitrate_reduction_to_ammonia` | 0.515 | 8 / 26 |
| 5 | aromatic_compound_degradation | `METPO:2000103\|pathways:aromatic_compound_degradation` | 0.517 | 13 / 39 |
| 6 | thermophilic | `METPO:1000616` | 0.521 | 13 / 42 |
| 7 | xylan catabolic process | `METPO:2000103\|GO:0045493` | 0.523 | 20 / 60 |
| 8 | spiral shaped | `METPO:1000684` | 0.525 | 10 / 32 |
| 9 | hydrogen_oxidation_dark | `METPO:2000103\|pathways:hydrogen_oxidation_dark` | 0.536 | 8 / 26 |
| 10 | sulfate_reduction | `METPO:2000103\|pathways:sulfate_reduction` | 0.542 | 6 / 19 |

Traits annotated on the same frequently reported taxa will rank alike here, so
the ranking does not show that cellulose-annotated taxa are singled out.

## 8. Counting rules

The October 16 catalog and the four-pairs literature file disagreed because they
used different labels, not different expansion rules. With experiment labels
both give the same counts:

| Count | Oct 16 catalog | Oct 16 literature file | Now |
| --- | --- | --- | --- |
| Studies with any feces decreased signature, PD / HIV / CRC | 46 / 17 / 27 | not reported | 49 / 19 / 31 |
| Studies reporting at least one annotated taxon decreased | 38 / 15 / 20 | 40 / 16 / 22 | 40 / 16 / 22 |
| *Roseburia* decreased / increased studies, PD; HIV; CRC | 18/4; 2/0; 5/0 | 20/4; 2/0; 6/0 | 20/4; 2/0; 6/0 |

The first row counts studies with any decreased feces signature in the condition
(`n_studies_dec`), the second only those whose decreased signature contains an
annotated taxon (`n_studies_a` in
[`profile_and_study_pair_jaccard.tsv`](profile_and_study_pair_jaccard.tsv) and
`n_supporting_studies_in_target_dir` in the four-pairs summary).

## 9. Files

| File | Contents |
| --- | --- |
| [`cellulose_decr_summary.json`](cellulose_decr_summary.json) | Every result, definition, universe, null summary and label-source comparison |
| [`cellulose_decr_provenance.json`](cellulose_decr_provenance.json) | Input hashes, code commit, seeds, package versions |
| [`pairwise_jaccard_significance.tsv`](pairwise_jaccard_significance.tsv) | Pairwise, mean and 3-way overlap under U1–U6; `endpoint` separates the mean distance from the 3-way core |
| [`sensitivity_definitions_jaccard.tsv`](sensitivity_definitions_jaccard.tsv) | Overlap under nine set definitions |
| [`compartment_jaccard_and_direction_stats.tsv`](compartment_jaccard_and_direction_stats.tsv) | Direction tests by contrast population, compartment and condition, at experiment and study level, with power notes |
| [`compartment_direction_interaction.tsv`](compartment_direction_interaction.tsv) | Fisher tests of direction by compartment |
| [`cellulose_decreased_taxa_catalog.tsv`](cellulose_decreased_taxa_catalog.tsv) | The 47 taxa with study and profile counts and flags |
| [`profile_and_study_pair_jaccard.tsv`](profile_and_study_pair_jaccard.tsv) | Overlap between individual profiles and studies across conditions |
| [`trait_convergence_benchmark.tsv`](trait_convergence_benchmark.tsv) | All 182 benchmarked traits with `trait_label` and `trait_predicate` |
| [`significant_results_summary.tsv`](significant_results_summary.tsv) | Headline tests with unit, population and p ([`SIGNIFICANT_RESULTS_SUMMARY.md`](SIGNIFICANT_RESULTS_SUMMARY.md)) |
| [`label_source_corrections.tsv`](label_source_corrections.tsv) | Experiments whose study-level label differed from their own |
| [`four_highlighted_pairs_taxon_drivers.tsv`](four_highlighted_pairs_taxon_drivers.tsv), [`four_highlighted_pairs_study_evidence.tsv`](four_highlighted_pairs_study_evidence.tsv), [`four_highlighted_pairs_summary.json`](four_highlighted_pairs_summary.json) | Taxa and studies behind the four highlighted condition–trait pairs ([`FOUR_HIGHLIGHTED_PAIRS_LITERATURE_SUPPORT.md`](FOUR_HIGHLIGHTED_PAIRS_LITERATURE_SUPPORT.md)) |
| [`figures_graphical/`](figures_graphical/), [`cellulose_decr_graphical_figures.pdf`](cellulose_decr_graphical_figures.pdf) | Four presentation figures (SVG, PDF, PNG) |
| [`figures/`](figures/), [`cellulose_decr_figures.pdf`](cellulose_decr_figures.pdf) | Three statistical reference figures |

The Tufte plates in `figures_tufte/` and the `cellulose_decr_slide_deck` files
made earlier in this directory have no generator, predate these corrections and
are not maintained here (#550).

## 10. Reproduce

```bash
COMMIT="$(git rev-parse HEAD)$(git diff --quiet HEAD -- data/analysis/cellulose_decr/scripts || echo +dirty)"
uv run python data/analysis/cellulose_decr/scripts/cellulose_decr_analysis.py --code-commit "$COMMIT"
uv run python data/analysis/cellulose_decr/scripts/cellulose_decr_figures.py
```

`--label-source study_modal` reproduces the October 16 tables. The run takes
about 20 seconds; the Monte Carlo nulls use fixed seeds, so reruns are identical.
