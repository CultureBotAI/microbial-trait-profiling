# What each experiment compares, and which way its direction points

Produced by `src/mtp/contrast_orientation.py`, driven by `just contrast-orientation`,
which writes `data/derived/contrast_orientation/`: a versioned manifest of every
experiment in the export (`contrast_manifest.tsv`), class counts within the analysed
corpus (`contrast_class_counts.tsv`), the historical-versus-oriented comparison
(`disease_control_comparison.tsv`) and a provenance record. It implements
[#188](https://github.com/CultureBotAI/microbial-trait-profiling/issues/188),
[#291](https://github.com/CultureBotAI/microbial-trait-profiling/issues/291),
[#554](https://github.com/CultureBotAI/microbial-trait-profiling/issues/554),
[#569](https://github.com/CultureBotAI/microbial-trait-profiling/issues/569),
[#574](https://github.com/CultureBotAI/microbial-trait-profiling/issues/574),
[#581](https://github.com/CultureBotAI/microbial-trait-profiling/issues/581), and
[#793](https://github.com/CultureBotAI/microbial-trait-profiling/issues/793).

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
(4,457 experiments, 2,257 of them in the analysed human paired corpus across 932 studies). Each row
carries the export's own fields verbatim - PMID, DOI, study design, condition and
EFO ID, body site and UBERON ID, host, both group names, the Group 1 definition,
both sample sizes and the description - with the curated direction convention, the
assigned `contrast_class`, `orientation`, `contrast_type`, and `control_tier`, the
matched evidence, flags, whether the experiment is in the analysed corpus, its
`validation_status` from the committed hand review, the classifier version
(`contrast-orientation-v6`) and the BugSigDB vintage. The export holds no Group 0
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
2. A group named as resolved disease (recovered, convalescent, post-COVID): `other`
   when paired against a comparison group (for example, healthy controls against
   recovered patients) and `longitudinal_timepoint` otherwise (active disease against
   recovered patients), with `contrast_type = resolved_disease` in both cases.
3. Group names that differ in time (baseline against follow-up, pre against post,
   day 4 against month 1, or two calendar years such as 2008 against 2010):
   `longitudinal_timepoint`. A time word shared by both groups ("at baseline",
   "10-day old") is a shared timepoint, not a contrast; "baseline" or
   "pre-treatment" / "before treatment" beside a healthy comparison only dates the
   disease sample (`disease_sample_dated_at_baseline`, including when the baseline
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
   - **Treated disease groups opposite a comparison group** ("Healthy Controls" vs
     "Parkinson's disease Patients treated with entacapone", in either native or
     reversed group order) compare controls against medicated cases rather than two
     arms of a trial; they are retained in `disease_vs_healthy` with flag
     `treated_disease_group` and stratified into `tier_3_hospital_comorbid_or_tissue`.
     Since `contrast-orientation-v4` (#569) a native comparison group (`Group 0`) may
     rest on weak evidence too ("Lean" vs "Obese patients treated with liraglutide"),
     for every condition; a reversed comparison group (`Group 1`) requires strong
     comparison evidence (`strong_1`) so a weak Group 1 comparison does not receive
     `treated_disease_group` on a row that later lands in `unresolved` (#574, #630).
   - **People with HIV on antiretroviral therapy opposite an HIV-negative comparison
     group** (`contrast-orientation-v4`, #569) are the same design: "healthy control"
     vs "HAART-treated HIV+ patients", "Controls/ People living without HIV (PLWOH)"
     vs "People living with HIV on Antiretroviral Therapy (PLWH on ART)", and "HIV-" vs
     "People Living With HIV on Antiretroviral therapy (ART+ PLWH)". The therapy is
     what the cases are on, not what the two groups differ by, so these are
     `disease_vs_healthy` with flag `treated_disease_group`. The cue is "ART-",
     "cART-" or "HAART-treated", or "on", "taking", "receiving" or "under" ART, cART,
     HAART or antiretroviral therapy, and in native group order a weak comparison
     group ("HIV-") is enough. "PLWH", "PWH" or the virus spelled out names the disease
     group as "HIV" does (#627). Treated against untreated people who all have HIV
     ("HAART-naive HIV-1 patients" vs "HAART-treated HIV-1 patients") has no
     comparison group and stays `treatment_exposure`, as do the same people before and
     on therapy. A contrast with a time word in either group keeps the treatment
     reading, because clearing it would leave only the time reading (see Limits).

     The treatment is recorded as the existing flag and tier rather than a new manifest
     column. Tier 3 already holds treated cases opposite healthy controls, and keeping
     these with them means the manifest columns are unchanged and the HIV contrasts
     are stratified like any other treated-case contrast. A consumer that wants people
     with HIV against people without it regardless of therapy filters on
     `contrast_class` or `contrast_type` (see "Selecting patients against controls"); one
     that wants to set aside cases whose group name states a therapy drops rows flagged
     `treated_disease_group`. An unflagged group is not known to be untreated: its name
     just states no therapy.
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
   "Metabolic toxicity"). The marker's sign may be a hyphen, a minus sign (U+2212) or
   an en dash (U+2013, "HIV–"). A minus sign or en dash touching its marker also
   counts before a word ("HIV– individuals"); a hyphen there does not, because it is
   as often a separator ("Periodontitis- Shallow pocket", "BV- Positive") as a minus.
   A condition that writes an acronym in capitals is also denied when the acronym is
   spelled out ("without Human immunodeficiency virus (HIV)" under "HIV infection"), and
   where HIV is the condition the initialism PLWOH or PWOH (people living without
   HIV) is a comparison too (`contrast-orientation-v4`, #569). When both groups share
   a generic strong comparison token ("normal-weight controls" vs "obese controls")
   and only one carries specific condition-negative evidence, the shared token is
   stripped so the specific comparison group wins. The orientation is `native` when
   the comparison is Group 0.
   It is `reversed` only when Group 1 carries a strong comparison word and Group 0
   has disease evidence; otherwise the experiment is `unresolved` with the flag
   `possible_reversal_unverified`. A contrast resting on weak evidence only is
   flagged `weak_comparison_evidence`.

   **Control-tier stratification (`contrast_type` and `control_tier`).** Not every
   experiment in `disease_vs_healthy` compares diseased patients against healthy,
   disease-free individuals. Each `disease_vs_healthy` experiment is therefore
   stratified along two orthogonal columns:

   - `contrast_type`:
     - `case_vs_control_subject` (961 analysed experiments, 559 studies): subject-level
       case versus control comparisons.
     - `comorbidity_within_disease` (29 analysed experiments, 17 studies): both groups
       share a background disease or clinical state and differ in a second condition
       (for example, "Maintenance Hemodialysis without Constipation" vs "Maintenance
       Hemodialysis with Constipation", flagged `comorbidity_control`).
     - `paired_tissue_autologous` (16 analysed experiments, 11 studies): tumour or
       lesional tissue compared against adjacent normal tissue from the same patients
       (flagged `paired_tissue_comparison`).
   - `control_tier`:
     - `tier_1_healthy_normal` (792 analysed experiments, 478 studies): explicit
       healthy, normal, unaffected, uninfected or full-term subject controls
       (`case_vs_control_subject` with strong comparison evidence and no hospital,
       comorbidity or treatment confound).
     - `tier_2_condition_negative_screened` (143 analysed experiments, 76 studies):
       asymptomatic condition-negative or screened subject controls ("non-diabetic",
       "non-obese", "HPV-", "Normal weight"). Together, Tier 1 + Tier 2 comprise
       **935 subject-level Disease Case vs. Healthy/Normal Control experiments across
       548 studies**.
     - `tier_3_hospital_comorbid_or_tissue` (71 analysed experiments, 42 studies):
       the 16 `paired_tissue_autologous` experiments, the 29
       `comorbidity_within_disease` experiments, and 26 `case_vs_control_subject`
       experiments whose control group is a hospital, surgical, symptomatic,
       precursor-lesion or benign-disease patient control (flagged
       `hospital_or_symptomatic_control`: 15 experiments in 9 studies, including 4
       colorectal-cancer contrasts against controls with a few polyps, #672) or whose
       disease group is confounded by named pharmacotherapy (`treated_disease_group`:
       11 experiments in 6 studies, 8 of them people with HIV on antiretroviral
       therapy).

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
- `historical_human_paired` (2,257 experiments, 932 studies): all human paired
  experiments in the historical corpus.
- `disease_vs_healthy_oriented` (1,006 experiments, 574 studies): all experiments
  classed `disease_vs_healthy`, with the 1 reversed experiment relabelled so `inc`
  means higher in the disease group.
- `disease_vs_healthy_strong_evidence` (834 experiments, 498 studies):
  `disease_vs_healthy` experiments without `weak_comparison_evidence`.
- `disease_vs_healthy_tier1_healthy_normal` (792 experiments, 478 studies):
  `tier_1_healthy_normal` subject-level healthy/normal controls only.
- `disease_vs_healthy_tier1_2_subject_controls` (935 experiments, 548 studies):
  combined `tier_1_healthy_normal` and `tier_2_condition_negative_screened`
  subject-level controls, excluding all `tier_3_hospital_comorbid_or_tissue`
  contrasts.

## Selecting patients against controls

`contrast_class = disease_vs_healthy` says only that exactly one group reads as a
comparison. Of its 1,006 analysed experiments, 16 are `paired_tissue_autologous`:
tumour or lesional tissue against adjacent normal tissue from the same patients, such
as 25699023-2 ("adjacent non-cancerous tissues" against "cancerous tissues", colorectal
mucosa). Another 29 are `comorbidity_within_disease`. A consumer that filters on the
class alone counts these as patients against controls. To keep only comparisons of
people with the condition against people without it:

- filter the manifest on `contrast_type = case_vs_control_subject` (961 analysed
  experiments, 559 studies). This drops the paired-tissue and comorbidity contrasts and
  keeps the Tier 3 hospital, symptomatic, precursor-lesion and treated-case contrasts; or
- filter on `control_tier` in {`tier_1_healthy_normal`,
  `tier_2_condition_negative_screened`} (935 experiments, 548 studies), which drops all
  of Tier 3. `oriented_profiles(profiles, classes, control_tiers=...)` applies this
  filter in code, and it defines `disease_vs_healthy_tier1_2_subject_controls`.

Either way the 1 reversed experiment still has to be relabelled, which
`oriented_profiles` does. It filters only on control tiers and strong evidence, so a
`contrast_type` filter has to be applied to the classes before it is called. Filtering on `contrast_type` and on `flags` can be combined:
dropping rows flagged `treated_disease_group` from the `case_vs_control_subject` rows
keeps hospital and symptomatic controls but sets aside cases named as treated.

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
`contrast-orientation-v3`, all 17 v2 disagreements are resolved by the rule
refinements above, raising agreement on the 60 reviewed experiments from **43/60
(71.7%)** to **60/60 (100.0%)**:

| Original v2 sampled class | Sampled | Agreed (v2) | Agreed (v3) | Remaining v3 disagreements |
| --- | --- | --- | --- | --- |
| disease_vs_healthy | 10 | 10 | 10 | - |
| disease_vs_disease | 10 | 6 | 10 | - (all 4 resolved: paired tissue `30584008-4`, multi-valued variable `36381318-5`, self-scoped negation `39417117-26`, term birth `39569373-1`) |
| treatment_exposure | 10 | 8 | 10 | - (both resolved: treated disease group `31815177-11` moved to Tier 3 `disease_vs_healthy`; shared exposure `37542310-2` moved to Tier 1 `disease_vs_healthy`) |
| longitudinal_timepoint | 10 | 6 | 10 | - (all 4 resolved: `Pre-DM` `24013136-6`, `Pre-term` `25444008-2`, and cross-sectional healthy-vs-recovered `33431578-2` and `39642873-2` routed to `other` with `contrast_type = resolved_disease`) |
| other | 10 | 10 | 10 | - |
| unresolved | 10 | 3 | 10 | - (all 7 resolved: `Age at assessment`, calendar years `2008` vs `2010`, PD duration bands, lean Ghanaians vs Americans, pollution levels, levodopa regimens, and control adults by body site) |

Neither `contrast-orientation-v4` (#569) nor `contrast-orientation-v5` (#574) changes
the two HIV rows in the 60-experiment sample (38162648-3, "HIV+ MASLD-" against "HIV-
MASLD+", stays `unresolved`), while `contrast-orientation-v5` resolves the two
remaining v2/v3 cross-sectional recovered-versus-healthy rows (`33431578-2` and
`39642873-2`, routed from `longitudinal_timepoint` to `other` with
`contrast_type = resolved_disease`), bringing agreement on the 60 reviewed experiments
to **60/60 (100.0%)**. `contrast-orientation-v6` (#581) moves no reviewed experiment:
29776865-2 ("PD duration less than 5 years" against "PD duration more than 5 years")
stays `disease_vs_disease`, flagged `age_band_not_timepoint` instead of
`shared_timepoint`, and agreement on the pinned export stays 60/60.

## Retained, reversed, excluded and unresolved

Within the 2,257 analysed experiments (932 studies):

| Class | Orientation | Experiments | Studies |
| --- | --- | --- | --- |
| disease_vs_healthy | native | 1,005 | 574 |
| disease_vs_healthy | reversed | 1 | 1 |
| disease_vs_disease | - | 236 | 160 |
| treatment_exposure | - | 345 | 166 |
| longitudinal_timepoint | - | 189 | 97 |
| other | - | 320 | 128 |
| unresolved | - | 161 | 100 |

Studies appear in more than one class when their experiments differ. 1,006
experiments from 574 studies are retained in `disease_vs_healthy`, 1 of them
reversed; 1,090 are excluded as other classes and 161 are unresolved. Within the
1,006 `disease_vs_healthy` experiments:

| Contrast type | Control tier | Experiments | Studies |
| --- | --- | --- | --- |
| `case_vs_control_subject` | `tier_1_healthy_normal` | 792 | 478 |
| `case_vs_control_subject` | `tier_2_condition_negative_screened` | 143 | 76 |
| `case_vs_control_subject` | `tier_3_hospital_comorbid_or_tissue` | 26 | 15 |
| `comorbidity_within_disease` | `tier_3_hospital_comorbid_or_tissue` | 29 | 17 |
| `paired_tissue_autologous` | `tier_3_hospital_comorbid_or_tissue` | 16 | 11 |

## Classifier versions

`CLASSIFIER_VERSION` in `src/mtp/contrast_orientation.py` names the rules. Every
manifest row, the provenance record and each workflow contrast selection carry it, so two
records with the same version must classify every experiment alike. A rule change that
moves any experiment's class, orientation, contrast type, control tier or flags needs a
new version and a section below saying which experiments moved, measured by classifying
the same pinned export under the old and the new rules. Rewording an evidence string
alone needs no new version.

`tests/test_contrast_orientation.py` enforces this (#793). For each version,
`tests/fixtures/contrast_orientation_guard_digests.json` pins SHA-256 digests of a fixed
set of cases and of the classifier's decisions on them, and the test fails when the
decisions change but the version does not, or when the version has no section below.
The cases are the repository's committed tables of experiment strings: the 60
experiments of the committed review, from the condition and group names it records, and
the 52 cases of the #291 and #569 regression tables. The guard sees only these cases.
They tell each version from the next (v3, v4, v5 and v6 give four different digests),
and a guard pinned before #739 would have failed that change. But of the 321 single-term
deletions from the classifier's term lists, 138 move at least one experiment of the
pinned export and only 37 of those move a guard case. A rule change can therefore still
move experiments without failing the guard, and the comparison on the pinned export
remains the measure of what a version changed.

## Changes from `contrast-orientation-v3` to `-v4` (#569)

Both versions were run on the same pinned inputs (`full_dump.csv` as under
Reproducing; profile matrix SHA-256
`c07d589f6cc88724e43c683ebd132bcd5c7f316f36b95df3ff9acad24ebbf62f`; experiment table
SHA-256 `403afea98d9b0a5b56939f68bcc9aace56a8eda350232a586806fe7bf30ad230`) and their
manifests compared experiment by experiment on `contrast_class`, `orientation`,
`contrast_type` and `control_tier`. Twelve experiments changed, ten of them in the analysed corpus, and
every one compares people without HIV against people with HIV. No other experiment
changed, and no evidence or flag changed elsewhere.

| From (v3) | To (v4) | Analysed | Whole export |
| --- | --- | --- | --- |
| `treatment_exposure` | `disease_vs_healthy` native, `case_vs_control_subject`, Tier 3, `treated_disease_group` | 8 | 9 |
| `disease_vs_disease` | `disease_vs_healthy` native, `case_vs_control_subject`, Tier 2 | 2 | 3 |

| Experiment | Group 0 | Group 1 | v3 class | v4 tier |
| --- | --- | --- | --- | --- |
| 28118207-1 | HIV-uninfected MSM | ART-treated HIV infected MSM | `treatment_exposure` | Tier 3, treated |
| 29411528-5 | healthy control | HAART-treated HIV+ patients | `treatment_exposure` | Tier 3, treated |
| 38383483-1, -2, -3, -7, -8 (and -6, not analysed) | Controls/ People living without HIV (PLWOH) | People living with HIV on Antiretroviral Therapy (PLWH on ART) | `treatment_exposure` | Tier 3, treated |
| 39880958-1 | HIV- | People Living With HIV on Antiretroviral therapy (ART+ PLWH) | `treatment_exposure` | Tier 3, treated |
| 33102244-1 | HIV– individuals (en dash) | HIV+ individuals | `disease_vs_disease` | Tier 2 |
| 37721699-2 (and -1 at Visit B, not analysed) | Pregnant women ≥ 18 years of age without Human immunodeficiency virus (HIV) at Postpartum (PPt) visit | Pregnant women ≥ 18 years of age with Human immunodeficiency virus (HIV) at Postpartum (PPt) visit. | `disease_vs_disease` | Tier 2 |

`disease_vs_healthy` gains 10 analysed experiments and 5 studies, 2 of the experiments in
Tiers 1 + 2. Analysed `disease_vs_healthy` experiments whose condition names HIV go from
30 (20 studies) to 40 (25 studies). 25699023-2, the colorectal tumour-against-adjacent-tissue
contrast raised in the same issue, was already `paired_tissue_autologous` in Tier 3
under v3 and is unchanged; what #569 adds for it is the guidance in "Selecting
patients against controls" above.

## Changes from `contrast-orientation-v4` to `-v5` (#574, #630, #672)

Compared on the same pinned inputs, `contrast-orientation-v5` keeps the
`disease_vs_healthy` class membership, Tier 1+2 subject controls (`970` experiments,
`563` studies), and all 105 rows of `disease_control_comparison.tsv` identical to
`contrast-orientation-v4`:

- **Cross-sectional resolved-disease contrasts opposite a comparison group**: 16
  experiments in the export (9 of them in the analysed corpus, across 8 studies,
  including reviewed experiments `33431578-2` and `39642873-2`) move from
  `longitudinal_timepoint` to `other` with `contrast_type = resolved_disease`, while
  active-versus-recovered cohorts without a healthy comparison group stay
  `longitudinal_timepoint`.
- **Precursor lesions in control groups (`PRECURSOR_LESION_CONTROL`, #672)**: 4
  colorectal-cancer experiments in Yachida et al. 2019 (`31171880-1` to `-4`, all in
  the analysed corpus, 1 study, "healthy control (normal and a few polyps)" against
  CRC stages) move within `tier_3_hospital_comorbid_or_tissue` from
  `comorbidity_within_disease` (flagged `comorbidity_control`) to
  `case_vs_control_subject` (flagged `hospital_or_symptomatic_control`), raising
  `case_vs_control_subject` across Tiers 1–3 from 993 experiments (573 studies) to
  997 experiments (574 studies). Explicitly negated adenoma controls
  ("Adenoma-free (controls)", "Non-adenoma") remain Tier 1 or Tier 2.
- **Hospital and symptomatic control flag (`hospital_or_symptomatic_control`)**: 23
  experiments in the export (16 of them in the analysed corpus, across 10 studies,
  including the 4 Yachida et al. 2019 contrasts above) carry the explicit
  `hospital_or_symptomatic_control` flag in `flags` under `case_vs_control_subject`
  (`tier_3_hospital_comorbid_or_tissue`).
- **Symmetric reversed `treated_disease_group` rule and review-table validation**:
  adds the reversed `Group 0` branch of `treated_disease_group` (requiring `strong_1`)
  and enforces label and verdict consistency in `read_review` (60/60 `reviewed_agreed`).

## Changes from `contrast-orientation-v5` to `-v6` (#581, #739, #793)

`contrast-orientation-v6` names the rules that #739 (f52923a) introduced without bumping
the version. From f52923a until #793 the classifier still reported
`contrast-orientation-v5`, so a record labelled v5 can hold either set of rules. The
classifier source it records tells them apart: `analysis_source_sha256` in the
contrast-orientation provenance or `classifier.source_sha256` in a workflow contrast
selection is `a0a7819ff0d202f05bd1282c75abbd0e9be129815f9fbdb7fc3861cfc1c769ff` under the
v5 rules and `e7d5f4c9ecf37ce566f3ef763362d445d6a32e673d0ee2435d768739924df5d1` under the
v6 rules.

- **Duration bands**: `DURATION_BAND` masks treatment and disease duration qualifiers
  ("for more than 1 year", "duration less than 5 years") out of `_timepoints` and
  records `age_band_not_timepoint`.
- **Spaced `pre - <therapy>` and parenthesized `(ART) Naïve`**: `SPACED_BASELINE_PRE`
  reads "pre - antiretroviral therapy(ART)" as a baseline before therapy, and
  `NEGATED_TREATMENT` reads it and "Antiretroviral therapy (ART) Naïve" as untreated
  disease.

The effect was measured by comparing per-experiment classifications in the
`contrast_selection.json` of two Tier 1+2 core runs on the same pinned export
(`full_dump.csv` SHA-256 `af84bfca…`): run `443ffc89…` at 94213f5, under the v5 rules,
and run `4a1ebaf5…` at 8b280a3, under the v6 rules. Four of the 4,435 experiments
changed class, contrast type, control tier or orientation, all of them HIV infection
contrasts that v5 read as treatment contrasts. Classifying all 4,457 experiments of the
export under the two classifier sources finds the same four:

| Experiment | Group 0 | Group 1 | v5 class, type, tier, orientation | v6 class, type, tier, orientation |
| --- | --- | --- | --- | --- |
| 29808701-1 | Age-matched HIV-uninfected women | HIV-positive women who have been taking ART for more than 1 year | `treatment_exposure`, `treatment_exposure`, `not_applicable`, `not_applicable` | `disease_vs_healthy`, `case_vs_control_subject`, `tier_3_hospital_comorbid_or_tissue`, `native`, flagged `treated_disease_group` |
| 33324725-1 | Human immunodeficiency virus negative (HIV-)Uninfected Controls | Human immunodeficiency virus positive(HIV+) pre - antiretroviral therapy(ART) | `treatment_exposure`, `treatment_exposure`, `not_applicable`, `not_applicable` | `disease_vs_healthy`, `case_vs_control_subject`, `tier_1_healthy_normal`, `native` |
| 33324725-4 | Human immunodeficiency virus positive(HIV+) pre- antiretroviral therapy (ART) in Fiebig stage 1-2(F1-2) | Human immunodeficiency virus positive (HIV+) pre- antiretroviral therapy(ART) in Fiebig stage 3(F3) | `treatment_exposure`, `treatment_exposure`, `not_applicable`, `not_applicable` | `disease_vs_disease`, `disease_vs_disease`, `not_applicable`, `not_applicable` |
| 38310301-1 | Healthy controls (HC) | Antiretroviral therapy (ART) Naïve at baseline | `treatment_exposure`, `treatment_exposure`, `not_applicable`, `not_applicable` | `disease_vs_healthy`, `case_vs_control_subject`, `tier_1_healthy_normal`, `native` |

Each of the four has a signature in one direction only, so none is in the analysed
corpus, and the Tier 1+2 core selection excludes all four as `unpaired` under both
versions; the two runs select the same experiments apart from the 35 that #742's
sequencing-type filter removes. Five more experiments changed their flags but not their
class, type, tier or orientation: 29776865-2 (`age_band_not_timepoint` instead of
`shared_timepoint`), 31824239-1, -2 and -3 (`age_band_not_timepoint` added) and
38310301-2 (`negated_treatment_term` added). Four changed their evidence string only
(33324725-3, -5, -6 and -7, whose time word `pre` now reads `pre-therapy`). No other
experiment changed, so the class counts, control tiers and populations reported above
are the same under both versions.

## The oriented and control-tier subsets against the historical corpus

The same two analyses were run on all five populations in
`disease_control_comparison.tsv`: the six fraction directions (studies scored by the
sign of their mean increased-minus-decreased fraction, exact binomial test,
Benjamini-Hochberg over six, zero-total profiles at fraction 0 as in the historical
analysis) and the 15 paired log-ratios (positive-count complete cases,
Benjamini-Yekutieli over 15). In all four oriented populations "increased" means
higher in the disease group.

| Comparison | Historical (932 studies) | All `disease_vs_healthy` (574) | Tier 1 + 2 subject controls (548) | Tier 1 healthy/normal only (478) |
| --- | --- | --- | --- | --- |
| facultatively aerobic | 0.580 [0.547, 0.613], q 1.0e-05 | 0.580 [0.538, 0.622], q 5.6e-04 | 0.573 [0.529, 0.616], q 0.0020 | 0.567 [0.520, 0.613], q 0.0067 |
| cellulose catabolic process | 0.433 [0.397, 0.469], q 4.3e-04 | 0.396 [0.351, 0.442], q 4.2e-05 | 0.407 [0.360, 0.454], q 5.4e-04 | 0.391 [0.342, 0.441], q 9.3e-05 |
| coccus shaped | 0.545 [0.511, 0.578], q 0.0083 | 0.571 [0.529, 0.613], q 0.0014 | 0.563 [0.519, 0.606], q 0.0052 | 0.555 [0.508, 0.601], q 0.024 |
| sphere shaped | 0.556 [0.521, 0.590], q 0.0021 | 0.548 [0.503, 0.592], q 0.036 | 0.551 [0.506, 0.596], q 0.027 | 0.558 [0.509, 0.606], q 0.021 |
| arsenate reduction | 0.563 [0.522, 0.603], q 0.0027 | 0.581 [0.529, 0.633], q 0.0026 | 0.587 [0.533, 0.640], q 0.0020 | 0.590 [0.533, 0.646], q 0.0038 |
| moderately halophilic | 0.574 [0.540, 0.608], q 6.6e-05 | 0.583 [0.538, 0.626], q 5.6e-04 | 0.581 [0.536, 0.626], q 0.0013 | 0.581 [0.532, 0.628], q 0.0038 |
| facultatively aerobic / cellulose | 0.595 [0.548, 0.641], q 0.0012 | 0.629 [0.568, 0.688], q 0.0012 | 0.623 [0.560, 0.684], q 0.0044 | 0.624 [0.556, 0.688], q 0.011 |
| cellulose / sphere shaped | 0.403 [0.357, 0.452], q 0.0016 | 0.391 [0.330, 0.454], q 0.014 | 0.398 [0.335, 0.464], q 0.032 | 0.402 [0.335, 0.471], q 0.079 (p 0.0063) |
| cellulose / arsenate reduction | 0.363 [0.301, 0.428], q 0.0012 | 0.346 [0.271, 0.428], q 0.0052 | 0.338 [0.261, 0.422], q 0.0044 | 0.342 [0.260, 0.430], q 0.012 |
| cellulose / moderately halophilic | 0.431 [0.382, 0.481], q 0.094 | 0.401 [0.337, 0.468], q 0.045 | 0.393 [0.328, 0.461], q 0.032 | 0.395 [0.326, 0.467], q 0.079 (p 0.0048) |

All 105 rows across the five populations, with study counts, are in
`disease_control_comparison.tsv`. Every one of the six fraction directions keeps its
direction and significance ($q \le 0.036$) across all five populations, including the
strictest `tier_1_healthy_normal` subset (478 studies) and the combined subject-level
control subset `tier_1_2_subject_controls` (548 studies). The four log-ratios above
keep their direction in every population. All four have $q \le 0.05$ in both
`disease_vs_healthy_oriented` (574 studies) and
`disease_vs_healthy_tier1_2_subject_controls` (548 studies). The two strongest ratios
(`facultatively_aerobic / cellulose` and `cellulose / arsenate_reduction`) remain
significant under Benjamini-Yekutieli even in the smaller `tier_1_healthy_normal`
subset ($q = 0.011$ and $q = 0.012$, with the other two at nominal $p = 0.0063$ and
$p = 0.0048$).

Two log-ratios crossed $q = 0.05$ between `contrast-orientation-v3` and `-v4`, both by
small margins and both in populations the HIV experiments joined. All 10 joined
`disease_vs_healthy_oriented`; 7 joined `disease_vs_healthy_strong_evidence`, which grew
from 863 to 870 experiments, while 33102244-1, 37721699-2 and 39880958-1 rest on weak
evidence:
`cellulose / moderately_halophilic` in `disease_vs_healthy_oriented` went from
$q = 0.0498$ to $0.052$, and `cellulose / sphere_shaped` in
`disease_vs_healthy_strong_evidence` (514 studies) went from $q = 0.043$ to $0.054$. In
`disease_vs_healthy_tier1_2_subject_controls`, which the two Tier 2 HIV experiments
joined, both of these ratios went from $q = 0.046$ to $0.038$. The
`tier_1_healthy_normal` population is unchanged.

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

- **Resolved disease is excluded from `disease_vs_healthy` with
  `contrast_type = resolved_disease`.** Healthy controls against recovered COVID-19
  cases (33431578-2) and against RSV convalescence (39642873-2) are routed to `other`
  because no time contrast separates the two groups, while active disease against
  recovered/post-COVID cohorts (34100340-2, 37119437-10) stays `longitudinal_timepoint`.
- **Sibling experiments of one study can be read inconsistently when their group
  labels differ in specificity.** 33294219-1 is `other` and 33294219-2 is
  `disease_vs_disease`, though both compare hypertensive Han against hypertensive
  Tibetan subjects at different altitudes.
- **"Normal (N-like group)" against "Case (OB-like group)" (25202708-1) stays
  `unresolved`.** The only rule that admits "Normal" followed by a parenthesis also
  admits "Normal (MN) group" against "Medium-dose (MMD) group" and "Model high-dose
  (MHD) group", two dose arms in a mouse study, so it would trade one recovered
  experiment for two new errors.
- **Five HIV-negative against HIV-positive contrasts in the analysed corpus remain
  excluded because of their curated condition scope (#581, #625, #626).** 25057045-1
  and -2 ("HIV-uninfected controls" against "Chronic HIV-infected cases" on suppressive
  ART in `group_1_definition`, feces) and 37660163-13 (healthy controls against HCV/HIV
  coinfection) stay `treatment_exposure` because their curated condition is "Response
  to antiviral drug", which the classifier reads before the group names. Letting a
  disease-versus-control group-name pair override `TREATMENT_CONDITION` (#626) was
  tested across all 78 export contrasts with a treatment condition and a comparison
  group: doing so would falsely admit faecal-microbiota-transplant recipient trials
  whose group names mirror donor phenotypes (such as `38477534-4`, "Healthy controls -
  Fecal Microbiota Transplantation (HC-FMT)" against "Chronic Constipation - Fecal
  Microbiota Transplantation (CC-FMT)", `31391545-1`, "Healthy Control Mice (HC Mice)"
  against "Schizophrenia Mice (SCZ Mice)", and `40568683-2`, "HC-FMT" against
  "HGG-FMT") as well as post-intervention comparisons (`40180909-1`, "Healthy spouses
  after prebiotics" against "Parkinson's Disease patients after prebiotics"). Similarly,
  41011452-3 and -4 (HIV-negative against HIV-positive people with oral candidiasis or
  HPV infection) stay `disease_vs_disease` because the condition has two values and both
  cohorts share the oral candidiasis or HPV infection; only a general comparison token
  (`GENERAL_COMPARISON`) is allowed to clear a multi-valued condition. Outside the
  analysed corpus, 29808701-1 ("for more than 1 year"), 33324725-1 ("pre -
  antiretroviral therapy(ART)"), and 38310301-1 ("Antiretroviral therapy (ART) Naïve at
  baseline") are resolved in `contrast-orientation-v6` (#581), while 29411528-3 and
  33324725-2 (untreated vs ART-treated HIV+) and 28743816-1 ("HIV seronegative" against
  "HIV seroconverted", a prospective seroconversion cohort) remain excluded. A hyphen
  followed by a word ("HIV- women") is not read as a sign; no HIV comparison group in
  the export is written that way.
- **Tumour against normal tissue is recognised only when a group name says so.**
  "Normal group" against "Tumor group" in colon mucosa after radical resection
  (40013792-1, 19 against 19) and "Normal samples free from Colorectal Cancer" against
  "Egyptian / African American / European American Samples with Colorectal Cancer" in
  colorectal mucosa (37483698-1 to -3, equal group sizes in each) are read as
  `case_vs_control_subject`, Tier 1. The names do not say whether the normal tissue
  came from the same patients; the equal group sizes suggest it may have, which only
  the papers can settle.

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
stratification and review-disagreement fixes ([#554](https://github.com/CultureBotAI/microbial-trait-profiling/issues/554)),
the v4 HIV-status and treated-case readings with the HIV contrasts that must stay
`treatment_exposure` and the tumour-against-adjacent-tissue contrast that the tiers keep
out of subject-level populations ([#569](https://github.com/CultureBotAI/microbial-trait-profiling/issues/569)),
the v5 symmetric treated-disease rule, cross-sectional resolved-disease routing,
precursor-lesion control handling, `hospital_or_symptomatic_control` flag and
review-table invariants
([#574](https://github.com/CultureBotAI/microbial-trait-profiling/issues/574),
[#630](https://github.com/CultureBotAI/microbial-trait-profiling/issues/630),
[#672](https://github.com/CultureBotAI/microbial-trait-profiling/issues/672)), the v6
duration and pre-ART baseline parsing
([#581](https://github.com/CultureBotAI/microbial-trait-profiling/issues/581)) and the
condition-scope contrasts that stay excluded
([#626](https://github.com/CultureBotAI/microbial-trait-profiling/issues/626)). It also
holds the classifier-version guard described under "Classifier versions"
([#793](https://github.com/CultureBotAI/microbial-trait-profiling/issues/793)).

The classifier, the reviews and this description were drafted with the assistance
of coding agents (Claude Code and Antigravity) and require review by the author
before manuscript use.
