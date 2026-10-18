# What each experiment compares, and which way its direction points

Produced by `src/bsdbng/contrast_orientation.py`, driven by `just contrast-orientation`,
which writes `data/derived/contrast_orientation/`: a versioned manifest of every
experiment in the export (`contrast_manifest.tsv`), class counts within the analysed
corpus (`contrast_class_counts.tsv`), the historical-versus-oriented comparison
(`disease_control_comparison.tsv`) and a provenance record. It implements
[#188](https://github.com/CultureBotAI/microbial-trait-profiling/issues/188).

## The question

BugSigDB records each signature as "Abundance in Group 1" relative to Group 0, so
a profile labelled `inc` holds taxa more abundant in Group 1. The analyses keep that
curated direction. Reading `inc` as "higher in disease" is right only when Group 1
is a disease group and Group 0 a healthy or unaffected comparison, and until now
nothing checked that. The human and paired filters do not: a before/after
treatment experiment or a mild-versus-severe comparison passes them just as a
case-control experiment does.

## The manifest

`contrast_manifest.tsv` has one row per experiment in the pinned v1.3.1 export
(4,457 experiments, 2,335 of them in the analysed human paired corpus). Each row
carries the export's own fields verbatim - PMID, DOI, study design, condition and
EFO ID, body site and UBERON ID, host, both group names, the Group 1 definition,
both sample sizes and the description - with the curated direction convention, the
assigned `contrast_class`, `orientation`, `contrast_type`, and `control_tier`, the
matched evidence, flags, whether the experiment is in the analysed corpus, its
`validation_status` from the committed hand review, the classifier version
(`contrast-orientation-v3`) and the BugSigDB vintage. The export holds no Group 0
definition, so every row carries the flag `no_group_0_definition_in_export`; a Group
0 reading rests on its name alone.

Whitespace inside the curated free-text fields is collapsed before writing. 187 of
the 4,457 descriptions and definitions contain a line break; quoting them keeps the
file valid CSV but puts one experiment on several physical lines, which `grep`,
`cut` and `head` then read wrongly.

## Classes, orientation, and control-tier stratification

Classification is lexical and ordered, and ambiguity is resolved toward exclusion
from the disease-versus-healthy population rather than into it.

1. A missing group name: `unresolved`.
2. A group named as resolved disease (recovered, convalescent): `longitudinal_timepoint`
   with `contrast_type = resolved_disease`.
3. Group names that differ in time (baseline against follow-up, pre against post,
   day 4 against month 1, or two calendar years such as 2008 against 2010):
   `longitudinal_timepoint`. A time word shared by both groups ("at baseline",
   "10-day old") is a shared timepoint, not a contrast; "baseline" or
   "pre-treatment" / "before treatment" beside a healthy comparison only dates the
   disease sample (`disease_sample_dated_at_baseline`,including when the baseline
   disease group is Group 0 and the healthy control is Group 1). An age or duration
   band ("2-12 years" against "18-62 years", "Controls aged > 10 years", "PD duration
   less than 5 years" against "PD duration more than 5 years") is not a time point
   and is flagged `age_band_not_timepoint`. Clinical stages and birth categories
   starting with `pre-` ("Pre-DM", "Pre-diabetes", "Pre-hypertension", "Pre-eclampsia",
   "Pre-term", "Pre-cancerous", "Pre-malignant", "Pre-menopausal", "Pre-frail") are
   likewise excluded from time-point matching via `NON_TEMPORAL_PRE_POST`.
4. A treatment or exposure in either group name (placebo, probiotic, smokers, users,
   responders, antibiotic, levodopa, carbidopa, polluted region), an intervention arm
   opposite a bare control, a treatment or exposure condition (response to a drug,
   antimicrobial agent, caesarean section, breastfeeding, air/environmental pollution),
   or an intervention named in the Group 1 definition or description and matched
   approximately to a group name: `treatment_exposure`.
   - **Negated treatment words** ("without antiretroviral treatment", "treatment-naive")
     name untreated disease and do not count.
   - **Shared exposure descriptors** ("Nonsmoking, nondrinking, and dentate healthy
     controls" vs "Nonsmoking, nondrinking, and dentate pre-treatment OSCC patients")
     are cohort matching criteria rather than the contrast when both groups carry the
     exact same treatment/exposure tokens and one group is a comparison group; they
     are flagged `shared_treatment_term`.
   - **Treated disease groups opposite a healthy control** ("Healthy Controls" vs
     "Parkinson's disease Patients treated with entacapone") compare healthy controls
     against medicated cases rather than two arms of a trial; they are retained in
     `disease_vs_healthy` with flag `treated_disease_group` and stratified into
     `tier_3_hospital_comorbid_or_tissue`.
5. One group named by time only, a condition that is a variable (population, age,
   age at assessment, ethnic group, physical activity, pregnancy, parity, gravidity,
   geography, socio-economic status, any "measurement"), or two groups that are
   levels of one variable (low against high LDL): `longitudinal_timepoint` or `other`.
   A condition can hold several comma-separated values and counts as a variable only
   when every value is one, so "Body mass index,Physical activity" is a variable and
   "Coinfection,Age" is not.
6. **Exactly one group reads as a comparison: `disease_vs_healthy`.** Strong
   comparison words are healthy, control, HC, unaffected, uninfected, disease-free,
   full-term / term birth / term infants, normal controls/subjects/groups/samples,
   normal tissues/mucosa/adjacent tissue, and volunteers. Weak ones are normal
   weight, lean, normoglycaemic, normotensive, a test marker negative only in that
   group ("HPV-" against "HPV+"), a bare "negative" whose subject is the condition,
   absence of the condition itself ("non-IBS", "no asthma", "free from colorectal
   cancer"), and self-scoped disease negation where one group explicitly negates the
   other group's disease noun ("Non-diabetic" vs "Diabetic" under condition
   "Metabolic toxicity"). When both groups share a generic strong comparison token
   ("normal-weight controls" vs "obese controls") and only one carries specific
   condition-negative evidence, the shared token is stripped so the specific
   comparison group wins. The orientation is `native` when the comparison is Group 0.
   It is `reversed` only when Group 1 carries a strong comparison word and Group 0
   has disease evidence; otherwise the experiment is `unresolved` with the flag
   `possible_reversal_unverified`. A contrast resting on weak evidence only is
   flagged `weak_comparison_evidence`.

   **Control-tier stratification (`contrast_type` and `control_tier`).** Not every
   experiment in `disease_vs_healthy` compares diseased patients against healthy,
   disease-free individuals. Each `disease_vs_healthy` experiment is therefore
   stratified along two orthogonal columns:

   - `contrast_type`:
     - `case_vs_control_subject` (983 analysed experiments, 567 studies): subject-level
       case versus control comparisons.
     - `comorbidity_within_disease` (33 analysed experiments, 18 studies): both groups
       share a background disease or clinical state and differ in a second condition
       (for example, "Maintenance Hemodialysis without Constipation" vs "Maintenance
       Hemodialysis with Constipation", flagged `comorbidity_control`).
     - `paired_tissue_autologous` (22 analysed experiments, 13 studies): tumour or
       lesional tissue compared against adjacent normal tissue from the same patients
       (flagged `paired_tissue_comparison`).
   - `control_tier`:
     - `tier_1_healthy_normal` (821 analysed experiments, 491 studies): explicit
       healthy, normal, unaffected, uninfected or full-term subject controls
       (`case_vs_control_subject` with strong comparison evidence and no hospital,
       comorbidity or treatment confound).
     - `tier_2_condition_negative_screened` (147 analysed experiments, 76 studies):
       asymptomatic condition-negative or screened subject controls ("non-diabetic",
       "non-obese", "HPV-", "Normal weight"). Together, Tier 1 + Tier 2 comprise
       **968 subject-level Disease Case vs. Healthy/Normal Control experiments across
       561 studies**.
     - `tier_3_hospital_comorbid_or_tissue` (70 analysed experiments, 38 studies):
       the 22 `paired_tissue_autologous` experiments, the 33
       `comorbidity_within_disease` experiments, and 15 `case_vs_control_subject`
       experiments whose control group is a hospital, surgical, symptomatic or
       benign-disease patient control or whose disease group is confounded by named
       pharmacotherapy (`treated_disease_group`).

   **What the comparison reading requires.** Group names say both more and less than
   they appear to, so four conditions scope the reading to the actual comparison
   context. Each withdrawal is recorded in `evidence` and flagged
   `comparison_reading_withdrawn` rather than dropped silently.

   - **A negation must scope the condition, not a modifier of it.** The first content
     word after the negator must name the condition, and a negation that runs on into
     another pathological finding does not count. "non-IBS" and "no asthma at 5 years
     old" deny the condition; "non-exacerbated asthma", "non-cachectic cancer",
     "Mild/Non-severe Preeclampsia", "Non-inflammatory depression", "No liver
     metastasis" and "Non-obstruction tumour group" deny a modifier or a neighbouring
     finding, and those subjects still have the disease. Filler words between the
     negator and the condition are skipped, so "healthy controls without any form of
     Parkinson's disease" keeps its reading, and words that say *which* case is meant
     ("prevalent", "functional", "idiopathic") are filler while words that say *how
     bad* ("severe", "exacerbated", "ulcerative") are not.
   - **A bare "negative" counts only where its subject is the condition.** The
     subject is the whole run of words before the sign, so "group B Streptococcus
     culture-negative" denies the organism, but any disease word in that run
     disqualifies it, because a disease named beside the sign is one the group has:
     "Her2 negative tumor", "crohn disease ASCA negative patients", "Metaformin
     negative", "Parkinson's disease dementia negative" and "HIV infected cART
     negative" all describe someone who is ill. The same requirement applies to a
     marker sign, so "tumor HP-" against "tumor HP+" is H. pylori status within
     gastric cancer rather than cancer against health.
   - **A group that only resembles the comparison group is not it.** "HC-like IBS
     (Healthy Control like Irritable Bowel Syndrome)" is an irritable-bowel subgroup
     whose profile looks like a healthy one.
   - **A condition with several values needs a group that is well in general.**
     Denying one value leaves the others standing: "Epstein-Barr virus non-infected
     (EBV- GC)" lacks the virus and has the gastric cancer, and "young HIV-infected
     group without coinfection" lacks the coinfection and has the HIV.

   **An untreated arm is not a healthy control.** A group labelled only "Control"
   opposite an "Experimental group", "Trial group", "Intervention group" or "usage
   group" is the untreated arm of an intervention, and the contrast is
   `treatment_exposure`. The trailing "group" is required, so the disease
   "Experimental Autoimmune Encephalomyelitis" and the "Trial Making Test" are not
   read as arms.
7. Both groups read as identical comparisons differing by demographic cohort or body
   site ("Lean Ghanaians" vs "Lean Americans", "Control Adults Saliva" vs "Control
   Adults Subgingival plaque"): `other`. Otherwise, if both groups read as
   comparisons: `unresolved`.
8. Both groups read as disease groups (disease, patients, stage, severity, activity,
   subtype, or naming the condition, with possessive suffixes stripped before
   initialism extraction so "PD" matches "Parkinson's disease"): `disease_vs_disease`.
9. Anything else: `unresolved`.

**Analysis populations.** Five populations are reported side by side in
`disease_control_comparison.tsv`:
- `historical_human_paired` (2,335 experiments, 959 studies): all human paired
  experiments in the historical corpus.
- `disease_vs_healthy_oriented` (1,038 experiments, 587 studies): all experiments
  classed `disease_vs_healthy`, with the 1 reversed experiment relabelled so `inc`
  means higher in the disease group.
- `disease_vs_healthy_strong_evidence` (863 experiments, 512 studies):
  `disease_vs_healthy` experiments without `weak_comparison_evidence`.
- `disease_vs_healthy_tier1_healthy_normal` (821 experiments, 491 studies):
  `tier_1_healthy_normal` subject-level healthy/normal controls only.
- `disease_vs_healthy_tier1_2_subject_controls` (968 experiments, 561 studies):
  combined `tier_1_healthy_normal` and `tier_2_condition_negative_screened`
  subject-level controls, excluding all `tier_3_hospital_comorbid_or_tissue`
  contrasts.

## Which fields the classifier reads

The study design, condition and both group names drive the classification. The Group
1 definition is read for a missing-definition flag and, with the description, for
intervention evidence: a group label can name an intervention the term list cannot
match - "Metaformin positive" misspells metformin, "Meditation group" is a practice
no drug list holds - while the curated prose spells it out. The prose supplies the
vocabulary and the group names decide whether it applies, because each intervention
found in the prose must also be approximately present in a group name. Without that
restriction, prose describing a cohort's history ("had history of taking
anti-diabetic medication ... metformin" beside "Healthy controls") turns a
disease-versus-healthy contrast into a drug trial. Cancer therapy is excluded from
that vocabulary for the same reason: a case-control description routinely dates its
samples against chemotherapy or radiotherapy, and including those words moved eight
genuine healthy-control contrasts out of the population. Nothing else is read from
the definition or the description, and the evidence string names which field supplied
each match.

## Development and review

The rules were written, reviewed against samples of every class, and corrected. The
reviews were done by a coding agent reading the group names, definition, design and
condition of each sampled experiment, not by the author; they check the rules against
the export's own words, not against the source papers, and they require review by the
author. Earlier rounds (seeds 188, 2026 and 9131 at the v1 rules) are superseded by
the review below, which is the one that is committed and reproducible.

**The committed review.** `methods/contrast_orientation_review.tsv` holds the sample,
the verdict for each experiment and a note for every disagreement. It was drawn at
seed 188 (10 analysed experiments per v2 contrast class, 60 in total). Under
`contrast-orientation-v3`, 15 of the 17 v2 disagreements are resolved by the rule
refinements above, raising agreement on the 60 reviewed experiments from **43/60
(71.7%)** to **58/60 (96.7%)**:

| Original v2 sampled class | Sampled | Agreed (v2) | Agreed (v3) | Remaining v3 disagreements |
| --- | --- | --- | --- | --- |
| disease_vs_healthy | 10 | 10 | 10 | - |
| disease_vs_disease | 10 | 6 | 10 | - (all 4 resolved: paired tissue `30584008-4`, multi-valued variable `36381318-5`, self-scoped negation `39417117-26`, term birth `39569373-1`) |
| treatment_exposure | 10 | 8 | 10 | - (both resolved: treated disease group `31815177-11` moved to Tier 3 `disease_vs_healthy`; shared exposure `37542310-2` moved to Tier 1 `disease_vs_healthy`) |
| longitudinal_timepoint | 10 | 6 | 8 | 2 intentional exclusions (`33431578-2` recovered COVID-19 and `39642873-2` RSV convalescence), which carry `contrast_type = resolved_disease` under `longitudinal_timepoint` |
| other | 10 | 10 | 10 | - |
| unresolved | 10 | 3 | 10 | - (all 7 resolved: `Age at assessment`, calendar years `2008` vs `2010`, PD duration bands, lean Ghanaians vs Americans, pollution levels, levodopa regimens, and control adults by body site) |

## Retained, reversed, excluded and unresolved

Within the 2,335 analysed experiments (959 studies):

| Class | Orientation | Experiments | Studies |
| --- | --- | --- | --- |
| disease_vs_healthy | native | 1,037 | 587 |
| disease_vs_healthy | reversed | 1 | 1 |
| disease_vs_disease | - | 254 | 168 |
| treatment_exposure | - | 359 | 173 |
| longitudinal_timepoint | - | 200 | 105 |
| other | - | 321 | 123 |
| unresolved | - | 163 | 102 |

Studies appear in more than one class when their experiments differ. 1,038
experiments from 587 studies are retained in `disease_vs_healthy`, 1 of them
reversed; 1,134 are excluded as other classes and 163 are unresolved. Within the
1,038 `disease_vs_healthy` experiments:

| Contrast type | Control tier | Experiments | Studies |
| --- | --- | --- | --- |
| `case_vs_control_subject` | `tier_1_healthy_normal` | 821 | 491 |
| `case_vs_control_subject` | `tier_2_condition_negative_screened` | 147 | 76 |
| `case_vs_control_subject` | `tier_3_hospital_comorbid_or_tissue` | 15 | 11 |
| `comorbidity_within_disease` | `tier_3_hospital_comorbid_or_tissue` | 33 | 18 |
| `paired_tissue_autologous` | `tier_3_hospital_comorbid_or_tissue` | 22 | 13 |

## The oriented and control-tier subsets against the historical corpus

The same two analyses were run on all five populations in
`disease_control_comparison.tsv`: the six fraction directions (studies scored by the
sign of their mean increased-minus-decreased fraction, exact binomial test,
Benjamini-Hochberg over six, zero-total profiles at fraction 0 as in the historical
analysis) and the 15 paired log-ratios (positive-count complete cases,
Benjamini-Yekutieli over 15). In all four oriented populations "increased" means
higher in the disease group.

| Comparison | Historical (959 studies) | All `disease_vs_healthy` (587) | Tier 1 + 2 subject controls (561) | Tier 1 healthy/normal only (491) |
| --- | --- | --- | --- | --- |
| facultatively aerobic | 0.581 [0.549, 0.614], q 5.7e-06 | 0.582 [0.540, 0.623], q 3.9e-04 | 0.575 [0.531, 0.617], q 0.0015 | 0.570 [0.523, 0.615], q 0.0044 |
| cellulose catabolic process | 0.430 [0.395, 0.466], q 2.3e-04 | 0.392 [0.348, 0.438], q 1.8e-05 | 0.405 [0.360, 0.452], q 4.2e-04 | 0.388 [0.340, 0.438], q 5.8e-05 |
| coccus shaped | 0.545 [0.512, 0.578], q 0.0073 | 0.570 [0.528, 0.612], q 0.0016 | 0.563 [0.520, 0.606], q 0.0050 | 0.557 [0.511, 0.603], q 0.017 |
| sphere shaped | 0.554 [0.520, 0.588], q 0.0029 | 0.547 [0.502, 0.590], q 0.038 | 0.549 [0.504, 0.594], q 0.034 | 0.554 [0.506, 0.602], q 0.027 |
| arsenate reduction | 0.560 [0.520, 0.600], q 0.0038 | 0.577 [0.525, 0.628], q 0.0043 | 0.585 [0.531, 0.637], q 0.0028 | 0.589 [0.532, 0.643], q 0.0039 |
| moderately halophilic | 0.572 [0.538, 0.605], q 1.0e-04 | 0.579 [0.535, 0.622], q 7.8e-04 | 0.577 [0.532, 0.622], q 0.0015 | 0.577 [0.529, 0.624], q 0.0039 |
| facultatively aerobic / cellulose | 0.595 [0.549, 0.640], q 0.0011 | 0.635 [0.574, 0.692], q 5.4e-04 | 0.624 [0.562, 0.683], q 0.0039 | 0.626 [0.560, 0.689], q 0.0079 |
| cellulose / sphere shaped | 0.407 [0.361, 0.456], q 0.0026 | 0.393 [0.332, 0.456], q 0.013 | 0.405 [0.343, 0.470], q 0.046 | 0.407 [0.341, 0.476], q 0.097 (p 0.0078) |
| cellulose / arsenate reduction | 0.365 [0.304, 0.429], q 0.0011 | 0.344 [0.270, 0.425], q 0.0034 | 0.340 [0.263, 0.424], q 0.0039 | 0.344 [0.263, 0.431], q 0.011 |
| cellulose / moderately halophilic | 0.435 [0.386, 0.485], q 0.12 | 0.404 [0.341, 0.470], q 0.050 | 0.399 [0.334, 0.467], q 0.046 | 0.400 [0.332, 0.471], q 0.094 (p 0.0057) |

All 105 rows across the five populations, with study counts, are in
`disease_control_comparison.tsv`. Every one of the six fraction directions keeps its
direction and significance ($q \le 0.038$) across all five populations, including the
strictest `tier_1_healthy_normal` subset (491 studies) and the combined subject-level
control subset `tier_1_2_subject_controls` (561 studies). All four significant
log-ratios keep their direction and significance ($q \le 0.05$) in both
`disease_vs_healthy_oriented` (587 studies) and `disease_vs_healthy_tier1_2_subject_controls`
(561 studies), and the two strongest ratios (`facultatively_aerobic / cellulose` and
`cellulose / arsenate_reduction`) remain significant under Benjamini-Yekutieli even
in the smaller `tier_1_healthy_normal` subset ($q = 0.0079$ and $q = 0.011$, with
the other two at nominal $p = 0.0078$ and $p = 0.0057$).

## Limits

The classification reads curated group names, the Group 1 definition and the
description, not the source papers, and is only as good as those fields; the
validation sample is small. Excluded classes are not wrong data - a severity or
treatment contrast answers a different question - and their directions are not
reversed or discarded, only kept out of this population. The oriented subset overlaps
the historical corpus, so agreement between them is a robustness check, not
independent replication. Nothing here measures disease activity, oxygen exposure or
abundance in patients (#179), and none of the procedures has been calibrated (#182).

Known remaining limitations:

- **Resolved disease is kept under `longitudinal_timepoint` with
  `contrast_type = resolved_disease`.** Healthy controls against recovered COVID-19
  cases (33431578-2) and against RSV convalescence (39642873-2) carry
  `contrast_type = resolved_disease` so they are distinguishable from longitudinal
  follow-up contrasts while remaining excluded from `disease_vs_healthy`.
- **Sibling experiments of one study can be read inconsistently when their group
  labels differ in specificity.** 33294219-1 is `other` and 33294219-2 is
  `disease_vs_disease`, though both compare hypertensive Han against hypertensive
  Tibetan subjects at different altitudes.
- **"Normal (N-like group)" against "Case (OB-like group)" (25202708-1) stays
  `unresolved`.** The only rule that admits "Normal" followed by a parenthesis also
  admits "Normal (MN) group" against "Medium-dose (MMD) group" and "Model high-dose
  (MHD) group", two dose arms in a mouse study, so it would trade one recovered
  experiment for two new errors.

## Reproducing

`just contrast-orientation` on the pinned BugSigDBExports v1.3.1 inputs (export
timestamp 2026-04-24_00:41_UTC, `full_dump.csv` SHA-256
`af84bfca7bfe76913a20c1fedbaadbfad4e4bf78605a2cfd7d00eb0c0fdd5788`). The provenance
record gives the code revision, the SHA-256 of the export, profile matrix, experiment
table, classifier source and committed review, the class counts, the control-tier
counts and the review tally. `just contrast-review-sample` redraws the review sample
from its seed.

`tests/test_contrast_orientation.py` holds a regression fixture, with the curated
strings copied verbatim, for every contrast that review
[#291](https://github.com/CultureBotAI/microbial-trait-profiling/issues/291) found
wrongly admitted to or wrongly kept out of the population, alongside the genuine
comparisons the narrowed rules must not take away, as well as the v3 control-tier
stratification and review-disagreement fixes ([#554](https://github.com/CultureBotAI/microbial-trait-profiling/issues/554)).

The classifier, the reviews and this description were drafted with the assistance
of coding agents (Claude Code and Antigravity) and require review by the author
before manuscript use.
