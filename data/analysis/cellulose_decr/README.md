# Cellulose-annotated taxa reported decreased in feces in Parkinson's disease, HIV infection and colorectal cancer

This directory compares the taxa annotated `METPO:2000103|GO:0030245`
(**cellulose catabolic process**) that curated BugSigDB signatures report decreased
in human feces in Parkinson's disease (PD), HIV infection (HIV) and colorectal
cancer (CRC). Every number below comes from one run of
[`src/mtp/cellulose_decr/analysis.py`](../../../src/mtp/cellulose_decr/analysis.py) through
`just cellulose-decr` (provenance in
[`cellulose_decr_provenance.json`](cellulose_decr_provenance.json)). The one
exception is a corpus-wide count in section 4, quoted from another analysis's
output. The method is described in
[`methods/cellulose_decr.md`](../../../methods/cellulose_decr.md).

These are annotation counts of curated signature taxa. They do not measure
cellulose degradation, fiber intake or short-chain fatty acids. A taxon counted
here was reported decreased in at least one study, not consistently lost
([`methods/biological_interpretation.md`](../../../methods/biological_interpretation.md)).
The trait and the three conditions were chosen as an exploratory lead from this
same corpus, so every p-value here is nominal. Section 4 adds adjusted values.

## What changed on 2026-10-19 (#580)

The analysis now reads contrast classifier v4 and is re-pinned on it.

- **Disease vs healthy means patients against controls.** The disease-versus-healthy
  populations keep `disease_vs_healthy` contrasts whose `contrast_type` is
  `case_vs_control_subject`. Tumour-vs-adjacent-tissue and comorbidity contrasts no
  longer count
  ([`methods/contrast_orientation.md`](../../../methods/contrast_orientation.md)).
  Tier 3 contrasts, including those whose cases are named as treated, still count;
  the scenarios set them aside.
- **The classifier makes the HIV and CRC calls.**
  - Since v4 it counts the HIV– vs HIV+ and ART-treated HIV+ vs uninfected contrasts as
    disease vs healthy (#569).
  - Since v3 it has marked tumour-vs-adjacent tissue as paired tissue (#554).
  - The scenarios that hand-coded these calls are retired. The condition-definition
    tables now set aside the contrasts whose cases are named as treated, keep only
    Tier 1–2 controls, or count the paired-tissue or the comorbidity contrasts as well.
  - The summary's `case_control_inventory` lists the 21 experiments in the three
    conditions that one of these rules sets aside.
- **One CRC study is miscalled (#672).** The classifier calls Yachida et al. 2019
  (`bsdb:31171880-1` to `-4`) comorbidity contrasts because its controls are "normal
  and a few polyps".
  - The primary analysis therefore leaves out these four paired stool experiments, all
    of which point down.
  - Counting the comorbidity contrasts (`comorbidity_counted`, which also adds
    `bsdb:36959686-4`) gives CRC stool 22 : 5 : 8 experiments (p = 0.0015) and
    15 : 3 : 7 studies (p = 0.0075).
  - Without them it is 18 : 5 : 8 (p = 0.011) and 14 : 3 : 7 (p = 0.013).
- **What moved.** The all-contrast tests, sets, nulls, catalog and benchmark are
  unchanged. Only their Benjamini-Hochberg q-values move, because those are taken over
  every test in the summary table. Among disease-versus-healthy contrasts:
  - HIV stool studies go from 9 : 0 : 2 (p = 0.0039) to 10 : 2 : 2 (p = 0.039).
  - CRC stool experiments go from 21 : 4 : 8 (p = 9.1e-4) to 18 : 5 : 8 (p = 0.011),
    mostly because of the Yachida study above.
  - The CRC study majority vote goes from 14 : 2 : 8 (p = 0.0042) to 13 : 2 : 9
    (p = 0.0074).
  - PD ties go from 10 to 12: the two PD experiments whose cases are named as treated
    (`bsdb:31815177-8`, `-11`) are tied.
  - The CRC stool-vs-GI study test goes from p = 0.13 to 0.021 (section 5). It rests
    on one experiment either way: counting `bsdb:25699023-2` gives 0.13 again.

## What changed on 2026-10-17 (#552)

The October 16 version of this directory had no generator. Given the same labels
(`--label-source study_modal`), the rebuilt analysis reproduces its sets, exact
tests, catalog counts and feces direction statistics. Other outputs differ:
- U5 was redefined (#567).
- U4's Monte Carlo p-values moved slightly (HIV–CRC 0.00149 to 0.00185).
- The compartment rows outside feces changed, including their annotated-taxon
  universes (non-fecal GI 78 to 71, oral 85 to 80, other 102 to 112).
- The sensitivity table no longer gives an expected mean distance for its
  mean rows.
- Two pairs of benchmark ranks swap.
- The four-pairs tables differ. Their October 16 counts match the experiment
  labels instead (section 8).

The rebuilt analysis then makes these corrections:

- **Labels.** Each profile now takes the condition and body site of its own
  experiment. The October 16 run used one modal label per study, which put
  saliva and nasal experiments into "PD feces" and moved other experiments out.
  61 experiments in these three conditions changed
  ([`label_source_corrections.tsv`](label_source_corrections.tsv); the same
  problem affects other scripts, #555).
- **Trait names.** The benchmark carries `trait_label` from
  `data/derived/taxa_trait_matrix_columns.tsv`. The canonical labels are:
  - `METPO:1000845`: **Acetogenesis**
  - `METPO:1000679`: **pleomorphic shaped**
  - `METPO:1000684`: **spiral shaped**
  - `METPO:1000616`: **thermophilic**
  - `METPO:1000674`: **filament shaped**

  Earlier figures and slides gave them other names. `METPO:1000682` (spore shaped)
  is not a matrix column.
- **A whole-signature null (U6).** U6 draws whole study signatures, so taxa that a
  study reports together stay together. Under it, the three conditions do not
  share detectably more taxa than random sets of feces studies (section 2).
- **New tests and checks:**
  - study-level and disease-versus-healthy tests of direction (section 4);
  - a compartment interaction test with power notes (section 5);
  - catalog flags (section 3);
  - reconciled counts (section 8);
  - a corrected literature file.

The adversarial review of #558 then found further problems. These changes are fixed in that PR:

- **Synonym condition labels (#561).** The analysis keeps only the exact labels.
  Other labels for the same disease are now listed in
  [`condition_label_inventory.tsv`](condition_label_inventory.tsv):
  - HIV: "Human immunodeficiency virus" and "HIV-1 infection";
  - CRC: colorectal carcinoma, adenocarcinoma, metastatic colorectal cancer and
    colon carcinoma.

  Three tables repeat the main results with those labels counted:
  - [`condition_definition_sensitivity.tsv`](condition_definition_sensitivity.tsv):
    the feces direction tests (sign and Wilcoxon, all contrasts and disease vs
    healthy);
  - [`condition_definition_overlap.tsv`](condition_definition_overlap.tsv):
    the overlap;
  - [`condition_definition_interaction.tsv`](condition_definition_interaction.tsv):
    the Feces-vs-GI test.

  The HIV result depends on that choice (section 4).
- **Pinned inputs (#562).** A run stops when an input's hash differs from the
  committed provenance, or when the contrast manifest came from another
  classifier version. This matters because `data/derived/` is regenerated by
  other work.
- **No core taxon is directly annotated (#565).** All 17 inherit the annotation
  from annotated species or strains (section 3).
- **One HIV "feces" study sampled rectal mucosa (#566).** It is reported with
  and without that study.
- **Null and test details (#567):**
  - U5 now draws only from taxa that can appear in an observed set;
  - multi-site body sites are grouped with Other non-GI;
  - studies are also collapsed by majority vote;
  - the pooled compartment test is stratified by condition.
- **Figures are regenerated, not committed (#560).**

| Quantity | Oct 16 (study labels) | Now (experiment labels) |
| --- | --- | --- |
| Decreased cellulose-annotated taxa, PD / HIV / CRC | 40 / 24 / 34 | 38 / 24 / 34 |
| Mean pairwise d_J; taxa in all three / union | 0.4957; 17 / 49 | 0.4793; 17 / 47 |
| PD paired experiments, down : up : tied | 57 : 26 : 20 | 58 : 24 : 17 |
| HIV paired experiments; Wilcoxon p | 19 : 6 : 3; 0.037 | 18 : 7 : 3; 0.135 |
| CRC paired experiments | 30 : 7 : 10 | 32 : 6 : 11 |
| Non-MediaDive traits in the benchmark | 141 (reported as 90) | 142 |

## 1. Overlap of the three taxon sets

A condition's set is every matrix taxon that bears the trait, at any rank, and
that at least one of the condition's feces decreased signatures reports.

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

Each null gives its expected value, with the range that holds 95% of its draws.
Monte Carlo p-values come from 100,000 draws; the size-matched U6 keeps the
5,284 draws whose set sizes match the observed ones. "≤ 1e-5" means no draw was
as extreme.

| Null | Expected mean d_J (95% range) | p | Expected taxa in all three (95% range) | p |
| --- | --- | --- | --- | --- |
| U1: uniform, 189 annotated matrix taxa | 0.909 (0.865–0.946) | ≤ 1e-5 | 0.87 | exact 4.8e-23 |
| U2: uniform, 135 annotated taxa seen in feces | 0.868 (0.817–0.911) | ≤ 1e-5 | 1.70 | exact 7.5e-18 |
| U3: uniform, 108 annotated taxa ever decreased in feces | 0.829 (0.772–0.878) | ≤ 1e-5 | 2.66 | exact 2.3e-14 |
| U4: random trait on the 471 decreased feces taxa of the three conditions | 0.775 (0.669–0.875) | ≤ 1e-5 | 5.09 (1–9) | ≤ 1e-5 |
| U5: U3 taxa weighted by reporting in other feces studies | 0.585 (0.501–0.665) | 0.0097 | 13.7 (10–17) | 0.076 |
| **U6: whole study signatures** | **0.527 (0.420–0.636)** | **0.19** | **18.3 (12–24)** | **0.73** |
| U6, draws with matching set sizes (5,284) | 0.513 (0.420–0.606) | 0.24 | 16.9 (13–21) | 0.58 |
| Within-condition reference: random subsets of each condition's own decreased feces taxa, any trait | 0.972 (0.945–0.994) | ≤ 1e-5 | 0.16 (0–1) | ≤ 1e-5 |

The nulls show three things:

- **Independent draws (U1–U4).** These nulls treat taxa as independent, and they make the overlap look
  extreme. So does the within-condition reference, which is not one of the six
  nulls: it draws each set from its own condition's decreased feces taxa,
  whatever their traits.
- **Weighted by reporting (U5).** Weighting the 108 taxa that can appear in a set
  by how often other feces studies report them decreased leaves a weaker signal.
  The mean distance gives p = 0.0097 and the 3-way core p = 0.076. An earlier U5
  also drew 27 taxa never reported decreased in feces, which inflated its
  significance (#567).
- **Whole study signatures (U6).** This null draws 49, 19 and 31 whole
  signatures at random from all 646 feces studies with a decreased signature.
  - The observed mean d_J is 0.479 against 0.527 expected: slightly more overlap
    than expected, but not detectably so (p = 0.19).
  - U6's larger expected core (18.3) comes from its larger simulated sets. With
    set sizes matched, the expected core is 16.9 against 17 observed (p = 0.58).
  - The pairwise U6 p-values are 0.21 (PD–HIV), 0.10 (PD–CRC) and 0.55 (HIV–CRC).

So most of the shared taxa are ones that many feces studies report decreased,
whatever the condition. The data give no evidence of a condition-specific
convergence beyond that. With synonym labels counted, the sets are 38 / 27 / 37 and
the core is 20 of 50, and U6 still gives p = 0.15 for the mean distance and 0.54
for the core. Without the rectal-mucosa study the U6 p-values are 0.26 and 0.79
([`condition_definition_overlap.tsv`](condition_definition_overlap.tsv)).

## 3. The 17 taxa in all three sets

The table gives stool studies reporting each taxon decreased / increased
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

Nine of the 17 are genera and eight are family-to-phylum nodes. No taxon in the
core has an annotation edge of its own in the KG
(`has_direct_trait_edge`; 186 NCBITaxon subjects carry the edge). A matrix taxon
bears the trait when it or any descendant has the edge, so every core taxon
inherits it. Among the matrix taxa that bear the trait, none of the 49 genera has
its own edge, while 38 of the 44 species do. Examples of inherited annotations:

- *Roseburia* from *R. faecis*;
- *Lachnospira* from *L. multipara*;
- *Anaerostipes* from *A. hadrus*;
- *Bacteroides* from a *B. cellulosilyticus* strain;
- *Clostridium* and *Ruminococcus* from several annotated species.

`matrix_trait_bearing_taxa_at_or_below` counts the annotated taxa beneath each
node: 17 under Lachnospiraceae and 86 under Bacillota. Only seven taxa in the 47-taxon
union have their own edge, all of them species.
`is_direct_or_genus_species` is therefore a rank flag for every core taxon.

Seven taxa are net decreased in all three conditions, five of them genera.
Several others are reported increased more often than decreased in one condition.
Counting decreased minus increased *profiles* instead of studies (sensitivity
definition 4, `is_net_decreased_all3_by_profiles`) also gives seven taxa, with
*Bacteroides* in place of *Anaerostipes*.

*Eubacterium* is in the core only through one HIV study, McHardy et al. 2013
(`bsdb:24451087`). That study sampled rectal mucosa with a sponge, but BugSigDB
curates it as feces (#566). Without it the core has 16 taxa
(sensitivity definition 9).

## 4. Direction within each condition

Each row compares paired experiments in feces. "Down" means the taxa annotated
with cellulose catabolic process made up a larger annotation share of the
decreased signature than of the increased one. Study rows give each study one
vote: the sign of its mean paired difference
([`compartment_jaccard_and_direction_stats.tsv`](compartment_jaccard_and_direction_stats.tsv)).

| Condition | Contrasts | Experiments down : up : tied | Sign test p | Studies down : up : tied | Sign test p | Study Wilcoxon p |
| --- | --- | --- | --- | --- | --- | --- |
| PD | all | 58 : 24 : 17 | 2.2e-4 | 27 : 14 : 4 | 0.060 | 0.0093 |
| PD | disease vs healthy | 47 : 15 : 12 | 5.8e-5 | 24 : 10 : 4 | 0.024 | 0.015 |
| HIV | all | 18 : 7 : 3 | 0.043 | 12 : 4 : 2 | 0.077 | 0.30 |
| HIV | disease vs healthy | 13 : 3 : 3 | 0.021 | 10 : 2 : 2 | 0.039 | 0.18 |
| CRC | all | 32 : 6 : 11 | 2.4e-5 | 19 : 2 : 7 | 2.2e-4 | 2.0e-4 |
| CRC | disease vs healthy | 18 : 5 : 8 | 0.011 | 14 : 3 : 7 | 0.013 | 0.0026 |
| All three | all | 108 : 37 : 31 | 3.0e-9 | 58 : 20 : 13 | 2.0e-5 | 4.1e-5 |
| All three | disease vs healthy | 78 : 23 : 23 | 3.6e-8 | 48 : 15 : 13 | 3.8e-5 | 1.2e-4 |

[`significant_results_summary.tsv`](significant_results_summary.tsv) gives
Benjamini-Hochberg q-values over its 29 tests. The four disease-vs-healthy
study-level rows have q = 0.041 (PD), 0.057 (HIV), 0.026 (CRC) and 2.2e-4
(pooled). The all-contrasts study-level PD and HIV rows have q = 0.079 and 0.093.

**How the rows change:**

- **Collapsing to studies.** This weakens PD and HIV, as expected when several
  experiments come from one cohort. PD across all contrasts stays significant
  by the study-level Wilcoxon test (p = 0.0093) but not by the sign test
  (p = 0.060).
- **Majority vote.** Taking each study's majority of experiment signs instead of
  its mean difference gives PD 28 : 12 : 5 (p = 0.017). For disease vs healthy
  it gives PD 24 : 9 : 5 (p = 0.014) and CRC 13 : 2 : 9 (p = 0.0074).
- **Disease-vs-healthy contrasts only.**
  - This strengthens the PD and HIV sign tests. PD's study-level Wilcoxon
    weakens slightly, from 0.0093 to 0.015.
  - It weakens CRC, though CRC stays below 0.05: the study sign-test p moves
    from 2.2e-4 to 0.013.
  - No experiment in the three conditions is classified as reversed.
- **Weak comparison evidence.** Dropping these experiments changes only HIV:
  12 : 2 : 1 experiments (p = 0.013) and 9 : 1 : 0 studies (p = 0.021).

**The HIV result is fragile.** It depends on which experiments count as HIV
disease-versus-healthy contrasts
([`condition_definition_sensitivity.tsv`](condition_definition_sensitivity.tsv)):

| HIV feces, disease vs healthy | Experiments | p | Studies | p |
| --- | --- | --- | --- | --- |
| Exact label "HIV infection" (primary) | 13 : 3 : 3 | 0.021 | 10 : 2 : 2 | 0.039 |
| + "Human immunodeficiency virus", "HIV-1 infection" labels | 13 : 5 : 3 | 0.096 | 10 : 4 : 2 | 0.18 |
| Contrasts whose cases are named as treated set aside (`treated_disease_group`) | 12 : 1 : 3 | 0.0034 | 10 : 0 : 2 | 0.0020 |
| Tier 1–2 controls only (Tier 3 hospital, comorbid and treated-case contrasts set aside) | 12 : 1 : 3 | 0.0034 | 10 : 0 : 2 | 0.0020 |
| Paired-tissue contrasts also counted | 13 : 3 : 3 | 0.021 | 10 : 2 : 2 | 0.039 |
| Comorbidity contrasts also counted | 13 : 3 : 3 | 0.021 | 10 : 2 : 2 | 0.039 |
| Exact label, without the rectal-mucosa study (bsdb:24451087) | 12 : 3 : 3 | 0.035 | 9 : 2 : 2 | 0.065 |

The two "Human immunodeficiency virus" studies compare people with and without
HIV, and in both the decreased signature has the smaller cellulose share. Across
all contrasts, HIV goes from 18 : 7 : 3 (p = 0.043) to 19 : 10 : 3 (p = 0.14) with
the synonyms counted. PD is unaffected. CRC stays essentially unchanged with its
synonyms:
- all contrasts: 37 : 10 : 12 experiments (p = 9.8e-5) and 22 : 3 : 7 studies
  (p = 1.6e-4);
- disease vs healthy: 16 : 4 : 7 studies (p = 0.012).

The pooled study-level result stays at p ≤ 1.4e-4 in every scenario.

Since classifier v4 the primary reading counts the HIV– vs HIV+ experiment
(`bsdb:33102244-1`) and the ART-treated HIV+ vs uninfected experiments as disease vs
healthy (#569, #580). The only two HIV studies that point up, `bsdb:38383483` and
`bsdb:39880958`, compare people with HIV on antiretroviral therapy against people
without HIV. Setting aside the contrasts whose cases are named as treated removes both.
That flag reads group names only, so these readings separate contrasts by their
wording, not by treatment (#676). `bsdb:32415070-1`, whose group definition names
people on antiretroviral therapy, is not flagged; it stays in that reading and points
down.

One HIV feces experiment mentions MSM in its group names; dropping it leaves
12 : 3 : 3 (p = 0.035). That lexical filter cannot separate HIV from sexual
behaviour. Two studies bear on this:

- Noguera-Julian et al. (2016, PMID 27077120) found that much of the reported
  HIV-associated shift tracks sexual practice.
- Vujkovic-Cvijin et al. (2020, PMID 32415070, a corpus study,
  `bsdb:32415070`) compared antiretroviral-treated people with HIV with
  seronegative controls matched for sexual practice. They still reported
  depletion of Lachnospiraceae and Ruminococcaceae.

A separate analysis points the same way:
`data/derived/contrast_orientation/disease_control_comparison.tsv` pools all
conditions and is not an input to this analysis. In it, the decreased signature
has the larger cellulose annotation share in 441 of 774 untied studies (p = 1.2e-4).
The pooling takes the historical human paired population, which does not depend
on the contrast classifier. The oriented disease-versus-healthy subset of that
file does depend on the classifier version, so it is not quoted here.

## 5. Outside feces

| Compartment | Decreased-signature studies, PD / HIV / CRC | Experiments down : up (all three) | Note |
| --- | --- | --- | --- |
| Feces | 49 / 19 / 31 | 108 : 37 | |
| Non-fecal GI (mucosa, lumen) | 2 / 4 / 14 | 9 : 12 | PD has 1 untied pair, HIV 5: p < 0.05 is unattainable |
| Oral / oropharyngeal | 2 / 2 / 4 | 3 : 5 | p < 0.05 unattainable in every condition |
| Other non-GI, including multi-site labels | 3 / 4 / 0 | 4 : 2 | no CRC study |

Outside feces, single conditions have at most 5 untied pairs, so no test can
reach p < 0.05. The exception is CRC mucosa and lumen, with 15 untied experiments.
There, the decreased signature had the larger share in only 4. The pooled strata
have 21, 8 and 6 untied pairs and could each reach p < 0.05.

The table below gives Fisher's exact test of direction by compartment
([`compartment_direction_interaction.tsv`](compartment_direction_interaction.tsv)):

| Comparison | Units | Feces | GI mucosa / lumen | Odds ratio (exact 95% CI) | p |
| --- | --- | --- | --- | --- | --- |
| CRC | experiments | 32 : 6 | 4 : 11 | 13.6 (2.9–81.0) | 1.2e-4 |
| CRC | studies | 19 : 2 | 3 : 8 | 21.6 (2.7–309) | 5.6e-4 |
| All three, stratified by condition | experiments | 108 : 37 | 9 : 12 | 4.79 (1.52–15.8) | 0.0049 |
| All three, stratified by condition | studies | 58 : 20 | 5 : 9 | 8.46 (1.79–49.9) | 0.0021 |
| CRC, disease vs healthy | experiments | 18 : 5 | 2 : 5 | 8.2 (0.99–111) | 0.026 |
| CRC, disease vs healthy | studies | 14 : 3 | 2 : 5 | 10.2 (1.07–158) | 0.021 |
| CRC, disease vs healthy, paired-tissue contrasts counted | studies | 14 : 3 | 3 : 4 | 5.66 (0.61–65.5) | 0.13 |
| All three, stratified, disease vs healthy | studies | 48 : 15 | 4 : 6 | 5.79 (1.03–37.3) | 0.039 |

**How to read these odds ratios:**

- They are conditional maximum-likelihood estimates; the TSV also gives sample
  and Mantel-Haenszel estimates.
- The stratified rows test a common odds ratio across PD, HIV and CRC. Crude
  pooling, also in the TSV, mixes condition with compartment: pooled stool pairs
  are mostly PD and pooled GI pairs mostly CRC.
- As in Fisher's test, each p-value counts outcomes no more likely than the
  observed one, while each interval comes from two one-sided tails. The two can
  disagree near 0.05, as in the disease-vs-healthy CRC experiment row.

The units are treated as independent. Only 11 of the 31 paired CRC mucosa and
lumen experiments compare patients with controls (`contrast_type`
`case_vs_control_subject`). The rest compare tumour with adjacent or normal
tissue, adenoma with cancer, cancer subgroups, or aspirin with placebo. Two of
them compare paracancerous or adjacent non-cancerous with cancerous tissue of the
same patients (`bsdb:22761885-2`, `bsdb:25699023-2`). The classifier calls them
disease vs healthy but `paired_tissue_autologous` (since v3), so they do not count
here.
- Counting `bsdb:25699023-2` changes the disease-vs-healthy study-level test from
  p = 0.021 to p = 0.13
  ([`condition_definition_interaction.tsv`](condition_definition_interaction.tsv)).
- `bsdb:22761885-2` changes nothing.
- That test therefore rests on one experiment either way.
- Counting the comorbidity contrasts instead, among them the Yachida study (#672),
  gives p = 0.017.

Two CRC studies sampled both compartments. The direction differs between feces and CRC mucosa or lumen. With
this sampling effort and this mix of contrasts, that does not show the pattern is
specific to feces.

## 6. Sensitivity to the set definition

| Definition | PD / HIV / CRC | Mean d_J | All three / union | U1 exact p, 3-way |
| --- | --- | --- | --- | --- |
| Primary: all ranks | 38 / 24 / 34 | 0.479 | 17 / 47 | 4.8e-23 |
| Genus and species only | 19 / 12 / 20 | 0.478 | 9 / 25 | 3.5e-12 |
| Profiles with cellulose z > 0 | 34 / 24 / 33 | 0.496 | 16 / 46 | 4.0e-22 |
| Disease-vs-healthy contrasts | 33 / 21 / 30 | 0.503 | 14 / 42 | 4.8e-20 |
| Disease-vs-healthy contrasts, strong comparison evidence | 33 / 19 / 30 | 0.499 | 14 / 41 | 4.9e-21 |
| Synonym condition labels counted | 38 / 27 / 37 | 0.456 | 20 / 50 | 1.8e-27 |
| Without the rectal-mucosa study | 38 / 23 / 34 | 0.496 | 16 / 47 | 2.5e-21 |
| All body sites | 43 / 28 / 44 | 0.529 | 19 / 60 | 2.0e-21 |
| Net decreased, profiles | 23 / 19 / 28 | 0.641 | 7 / 40 | 5.6e-9 |
| Net decreased, studies | 22 / 19 / 26 | 0.643 | 7 / 39 | 2.2e-9 |
| Every decreased feces taxon, any trait (reference, not a null) | 284 / 151 / 231 | 0.775 | 51 / 471 | 6.5e-103 |

The U1 p-values in this table assume independent draws, so they overstate
significance in the way section 2 shows
([`sensitivity_definitions_jaccard.tsv`](sensitivity_definitions_jaccard.tsv)).
Restricting to genera and species leaves the mean distance unchanged.

## 7. Trait benchmark

182 matrix traits have at least 10 decreased feces taxa in each condition; 142 of
them are not MediaDive growth media. Ranked by mean pairwise d_J, cellulose
catabolic process is second of the 142 (23rd of all 182). The rank describes
this corpus and is not a test
([`trait_convergence_benchmark.tsv`](trait_convergence_benchmark.tsv); median
d_J 0.679; `percentile_excluding_growth_media` is the rank divided by 142).

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

Traits annotated on the same frequently reported taxa rank alike here. The
ranking therefore does not show that cellulose-annotated taxa are singled out.

## 8. Counting rules

The October 16 catalog and the four-pairs literature file disagreed because they
used different labels, not different expansion rules. With experiment labels,
both give the same counts:

| Count | Oct 16 catalog | Oct 16 literature file | Now |
| --- | --- | --- | --- |
| Studies with any feces decreased signature, PD / HIV / CRC | 46 / 17 / 27 | not reported | 49 / 19 / 31 |
| Studies reporting at least one annotated taxon decreased | 38 / 15 / 20 | 40 / 16 / 22 | 40 / 16 / 22 |
| *Roseburia* decreased / increased studies, PD; HIV; CRC | 18/4; 2/0; 5/0 | 20/4; 2/0; 6/0 | 20/4; 2/0; 6/0 |

The first row counts studies with any decreased feces signature in the condition
(`n_studies_dec`). The second counts only those whose decreased signature contains
an annotated taxon: `n_studies_a` in
[`profile_and_study_pair_jaccard.tsv`](profile_and_study_pair_jaccard.tsv), and
`n_supporting_studies_in_target_dir` in the four-pairs summary.

## 9. Files

| File | Contents |
| --- | --- |
| [`cellulose_decr_summary.json`](cellulose_decr_summary.json) | Every result, definition, universe, null summary, label inventory, case-control inventory and label-source comparison |
| [`cellulose_decr_provenance.json`](cellulose_decr_provenance.json) | Input hashes, contrast-classifier version, code commit, seeds, package versions |
| [`pairwise_jaccard_significance.tsv`](pairwise_jaccard_significance.tsv) | Pairwise, mean and 3-way overlap under U1–U6 and the within-condition reference (`within_cond_random_subset_*`); `endpoint` separates the mean distance from the 3-way core |
| [`sensitivity_definitions_jaccard.tsv`](sensitivity_definitions_jaccard.tsv) | Overlap under eleven set definitions |
| [`compartment_jaccard_and_direction_stats.tsv`](compartment_jaccard_and_direction_stats.tsv) | Direction tests by contrast population, compartment and condition, at experiment and study level (mean and majority vote), with power notes |
| [`compartment_direction_interaction.tsv`](compartment_direction_interaction.tsv) | Fisher tests of direction by compartment, crude and stratified by condition |
| [`condition_definition_sensitivity.tsv`](condition_definition_sensitivity.tsv), [`condition_definition_overlap.tsv`](condition_definition_overlap.tsv), [`condition_definition_interaction.tsv`](condition_definition_interaction.tsv) | Feces direction, cross-condition overlap (with U6) and the Feces-vs-GI test, under each condition definition: synonym labels, cases named as treated set aside, Tier 1–2 controls only, paired-tissue or comorbidity contrasts counted, and without the rectal-mucosa study |
| [`condition_label_inventory.tsv`](condition_label_inventory.tsv) | Every human condition label naming PD, HIV or CRC, with its status |
| [`cellulose_decreased_taxa_catalog.tsv`](cellulose_decreased_taxa_catalog.tsv) | The 47 taxa with study and profile counts and flags, including `has_direct_trait_edge` |
| [`profile_and_study_pair_jaccard.tsv`](profile_and_study_pair_jaccard.tsv) | Overlap between individual profiles and studies across conditions |
| [`trait_convergence_benchmark.tsv`](trait_convergence_benchmark.tsv) | All 182 benchmarked traits with `trait_label` and `trait_predicate` |
| [`significant_results_summary.tsv`](significant_results_summary.tsv) | Headline tests with unit, population, p and BH q ([`SIGNIFICANT_RESULTS_SUMMARY.md`](SIGNIFICANT_RESULTS_SUMMARY.md)) |
| [`label_source_corrections.tsv`](label_source_corrections.tsv) | Experiments whose study-level label differed from their own |
| [`four_highlighted_pairs_taxon_drivers.tsv`](four_highlighted_pairs_taxon_drivers.tsv), [`four_highlighted_pairs_study_evidence.tsv`](four_highlighted_pairs_study_evidence.tsv), [`four_highlighted_pairs_summary.json`](four_highlighted_pairs_summary.json) | Taxa and studies behind the four highlighted condition–trait pairs ([`FOUR_HIGHLIGHTED_PAIRS_LITERATURE_SUPPORT.md`](FOUR_HIGHLIGHTED_PAIRS_LITERATURE_SUPPORT.md)) |

Some outputs are written here but not committed: under CONTRIBUTING.md, images
stay out of git outside `trait-visualization/`.

- **Figures.** `just cellulose-decr-figures` writes four presentation figures to
  `figures_graphical/` and three statistical reference figures to `figures/`, as
  SVG, PDF and PNG. Their names say what each shows, not a conclusion (#570):
  `fig1_overlap_and_null_models`, `fig2_shared_taxa_and_set_definitions`,
  `fig3_direction_by_body_site` and `fig4_trait_overlap_ranking`; and
  `fig1_pairwise_and_3way_jaccard_nulls`, `fig2_direction_by_body_site_and_definition`
  and `fig3_core_taxa_and_trait_benchmark`.
- **Combined PDFs.** It also writes `cellulose_decr_graphical_figures.pdf` and
  `cellulose_decr_figures.pdf`.
- **Slide captions.** `figures/figure_captions.json` gives each figure a slide
  title, takeaway, caveat and evidence table, worded from the same tables as the
  figure, plus two statements other slides quote. The selected-figure deck
  (`scripts/build_slide_deck.py`) takes its cellulose slides from it, so a rerun
  updates the slides with the figures (#556).

Reruns on the same tables give identical bytes. The Tufte plates in
`figures_tufte/` and the `cellulose_decr_slide_deck` files made earlier in this
directory have no generator and predate these corrections. The deck no longer
uses them (#550, #556).

## 10. Reproduce

```bash
just cellulose-decr                # tables, summary and provenance
just cellulose-decr-figures        # figures (gitignored)
```

The run takes about a minute. The Monte Carlo nulls use fixed seeds, so reruns
on the same inputs give identical tables; only the provenance timestamp changes.

**Inputs.** Nine of the ten inputs are gitignored. They are built by:

- `just download`, `just ingest` and `just study-tsvs` (BugSigDB);
- `just download-kg`, `just taxa-trait-matrix` and `just study-trait-profiles`
  (KG-Microbe);
- `just contrast-orientation` (contrast manifest).

**Drift check.** The run stops in four cases:
- an input's hash differs from the `pinned_inputs` of
  [`cellulose_decr_provenance.json`](cellulose_decr_provenance.json);
- there is no pinned provenance to compare with;
- the contrast manifest came from another classifier version;
- an input changed while the run was reading it.

The options:
- From a worktree, use `--input-root <checkout>` to point at a checkout that
  holds the pinned inputs.
- `--allow-input-drift` runs on changed inputs and records them, but keeps the old
  pins, so the next ordinary run still refuses them.
- `--repin` adopts the current inputs as the new pins, for example after a
  BugSigDB release.

The pins are the v4 contrast manifest and the profile metadata table with each
experiment's own labels (#555), adopted with `--repin` in #580. The committed tables
reproduce from the pinned inputs at the provenance's `code_commit` (#642). A commit
made before a re-pin still carries the old pins, so there, pass
`--pinned-provenance` with this provenance file, or check out the commit that adds
it; its analysis code is the same (#674).

**Other labels.**
 Write these runs to another directory, so the committed
tables are not overwritten:
- the October 16 labels (what they reproduce is listed under "What changed"):
  `just cellulose-decr --label-source study_modal --output-dir <elsewhere>`;
- synonym labels throughout:
  `just cellulose-decr --condition-labels synonyms --output-dir <elsewhere>`.
