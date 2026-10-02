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
assigned contrast class and orientation, the matched evidence, flags, whether the
experiment is in the analysed corpus, its `validation_status` from the committed
hand review, the classifier version (`contrast-orientation-v2`) and the BugSigDB
vintage. The export holds no Group 0 definition, so every row carries the flag
`no_group_0_definition_in_export`; a Group 0 reading rests on its name alone.

Whitespace inside the curated free-text fields is collapsed before writing. 187 of
the 4,457 descriptions and definitions contain a line break; quoting them keeps the
file valid CSV but puts one experiment on several physical lines, which `grep`,
`cut` and `head` then read wrongly.

## Classes and the orientation rule

Classification is lexical and ordered, and ambiguity is resolved toward exclusion
from the disease-versus-healthy population rather than into it.

1. A missing group name: `unresolved`.
2. A group named as resolved disease (recovered, convalescent): `longitudinal_timepoint`.
3. Group names that differ in time (baseline against follow-up, pre against post,
   day 4 against month 1): `longitudinal_timepoint`. A time word shared by both
   groups ("at baseline") is a shared timepoint, not a contrast; "baseline" beside a
   healthy comparison only dates the disease sample. An age band is not a time point:
   "2-12 years" against "18-62 years" and "Controls aged > 10 years" say how old the
   subjects were, not when they were sampled, and are flagged
   `age_band_not_timepoint` (27 experiments).
4. A treatment or exposure in either group name (placebo, probiotic, smokers, users,
   responders, antibiotic), an intervention arm opposite a bare control, a treatment
   or exposure condition (response to a drug, antimicrobial agent, caesarean section,
   breastfeeding), or an intervention named in the Group 1 definition or description
   and matched approximately to a group name: `treatment_exposure`. A negated
   treatment word ("without antiretroviral treatment", "treatment-naive") names
   untreated disease and does not count.
5. One group named by time only, a condition that is a variable (population, age,
   ethnic group, pregnancy, any "measurement"), or two groups that are levels of one
   variable (low against high LDL): `longitudinal_timepoint` or `other`. A condition
   can hold several comma-separated values and counts as a variable only when every
   value is one, so "Age" is a variable and "Coinfection,Age" is not.
6. **Exactly one group reads as a comparison: `disease_vs_healthy`.** Strong
   comparison words are healthy, control, HC, unaffected, uninfected, disease-free,
   normal controls/subjects/groups/samples and volunteers. Weak ones are normal
   weight, lean, normoglycaemic, a test marker negative only in that group ("HPV-"
   against "HPV+"), a bare "negative" whose subject is the condition, and absence of
   the condition itself ("non-IBS", "no asthma", "free from colorectal cancer"). The
   orientation is `native` when the comparison is Group 0. It is `reversed` only when
   Group 1 carries a strong comparison word and Group 0 has disease evidence;
   otherwise the experiment is `unresolved` with the flag
   `possible_reversal_unverified`. A contrast resting on weak evidence only is
   flagged `weak_comparison_evidence`.

   **What the comparison reading requires.** Group names say both more and less than
   they appear to, so four conditions scope the reading to the actual comparison
   context. Each withdrawal is recorded in `evidence` and flagged
   `comparison_reading_withdrawn` rather than dropped silently (10 experiments in the
   corpus).

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
7. Both groups read as comparisons: `unresolved`.
8. Both groups read as disease groups (disease, patients, stage, severity, activity,
   subtype, or naming the condition): `disease_vs_disease`.
9. Anything else: `unresolved`.

**Analysis population.** Experiments classed `disease_vs_healthy`. For the one
reversed experiment the increased and decreased profiles are swapped, so that in
the whole population `inc` means higher in the disease group. The historical
results are not overwritten; the oriented population is reported beside them.

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
the verdict for each experiment and a note for every disagreement. It is drawn by

```sh
just contrast-review-sample methods/contrast_orientation_review.tsv 188 10
```

which samples 10 analysed experiments per contrast class. An experiment's place in
the draw is a digest of the seed, its class and its own identifier rather than a
shuffle, so the same seed always gives the same 60 experiments, an experiment's place
does not depend on what else is in the corpus, and asking for a larger sample extends
this one rather than replacing it. A class with fewer than ten experiments contributes
all of them. The verdicts are filled in by
hand and committed; `just contrast-orientation` joins the file back, writing each
experiment's `validation_status` into the manifest as `reviewed_agreed`,
`reviewed_disagreed` or `not_reviewed`, and records the file's SHA-256 and the tally
in the provenance. A sampled row with a blank verdict is an error rather than a
silent non-answer.

| Class | Sampled | Agreed | What the disagreements were |
| --- | --- | --- | --- |
| disease_vs_healthy | 10 | 10 | - |
| disease_vs_disease | 10 | 6 | normal against tumour tissue read as two diseases; a body-mass contrast under a condition only half of whose values are variables; "Non-diabetic" against "Diabetic" where the curated condition is "Metabolic toxicity", so the negation has no condition to reach; term against preterm birth, where the unaffected group matches no comparison word |
| treatment_exposure | 10 | 8 | healthy controls against Parkinson's patients "treated with entacapone"; healthy controls against OSCC patients where "Nonsmoking, nondrinking" is a matching criterion shared by both groups rather than the contrast |
| longitudinal_timepoint | 10 | 6 | "Pre-DM" against type 2 diabetes, where "pre" names a disease stage; full-term against pre-term infants, where "Pre-term" is a birth category; healthy controls against recovered COVID-19 cases and against RSV convalescence, excluded deliberately by the resolved-disease rule but under a class that describes neither contrast |
| other | 10 | 10 | - |
| unresolved | 10 | 3 | seven conservative misses: an age contrast whose condition is "Age at assessment"; one cohort sampled in 2008 and 2010; Parkinson's disease by duration; lean Ghanaians against lean Americans; two pollution levels; two levodopa regimens; control adults by body site |

43 of 60 agreed. The figure is deliberately less flattering than the "24 of 25" it
replaces, because the sample is stratified over all six classes rather than weighted
toward the ones that were already right, and because a disagreement is recorded
wherever a better class exists, not only where the experiment would reach the
analysis population. The distribution matters more than the total:

- **All ten sampled disease-versus-healthy experiments were correct, in class and in
  orientation.** That is the population the Table 1 analyses consume, and it is the
  part of the classification this change was meant to narrow.
- **`unresolved` is where the conservatism lives.** Seven of ten sampled experiments
  belong in some class. None of the seven belongs in disease-versus-healthy, so the
  conservatism costs coverage of the excluded classes rather than the analysis
  population.
- **Three of the four disease-versus-disease disagreements are experiments that
  belong in disease-versus-healthy**, so eligibility is now narrow enough to lose
  genuine comparisons as well as spurious ones. Each is a named case a later change
  can be measured against.
- With ten experiments per class the per-class rates are very imprecise. 10 of 10
  puts the exact 95% interval for the disease-versus-healthy agreement share at 0.69
  to 1.00, which bounds nothing useful; an earlier draw of the same size, discarded
  when the sampling method changed, gave 8 of 10 on the same rules. The value of the
  table is the named cases, not the rate.

## Retained, reversed, excluded and unresolved

Within the 2,335 analysed experiments (959 studies):

| Class | Orientation | Experiments | Studies |
| --- | --- | --- | --- |
| disease_vs_healthy | native | 1,003 | 569 |
| disease_vs_healthy | reversed | 1 | 1 |
| disease_vs_disease | - | 260 | 173 |
| treatment_exposure | - | 366 | 174 |
| longitudinal_timepoint | - | 200 | 105 |
| other | - | 275 | 105 |
| unresolved | - | 230 | 136 |

Studies appear in more than one class when their experiments differ. 1,004
experiments from 569 studies are retained, 1 of them reversed; 1,101 are excluded
as other classes and 230 are unresolved. Of the retained experiments 166 (83
studies) rest on weak comparison evidence only; 838 experiments from 498 studies
rest on strong evidence.

Against the v1 rules the retained population falls by 48 experiments and 20 studies
(1,052 from 589 to 1,004 from 569). 102 experiments change class in all. 66 leave
disease-versus-healthy: 30 for disease-versus-disease, 22 for `unresolved` where a
marker sign no longer reads as absence of the condition, 12 for `treatment_exposure`
and 2 for `other` (the pregnancy contrasts). 18 enter it, 15 from `unresolved` and 3
from disease-versus-disease. The remaining 18 move between the excluded classes.
Every analysed experiment still carries exactly one class, and the class counts still
sum to 2,335.

## The oriented subset against the historical corpus

The same two analyses were run on three populations, each with its own native
families: the six fraction directions (studies scored by the sign of their mean
increased-minus-decreased fraction, exact binomial test, Benjamini-Hochberg over
six, zero-total profiles at fraction 0 as in the historical analysis) and the 15
paired log-ratios (positive-count complete cases, Benjamini-Yekutieli over 15). The
share is of non-tied studies shifting toward the increased side, with an exact 95%
interval for that share only. In the oriented populations "increased" means higher
in the disease group.

| Comparison | Historical (959 studies) | Disease vs healthy (569) | Strong evidence only (498) |
| --- | --- | --- | --- |
| facultatively aerobic | 0.582 [0.550, 0.614], q 4.1e-06 | 0.584 [0.541, 0.626], q 2.1e-04 | 0.574 [0.528, 0.619], q 0.0030 |
| cellulose catabolic process | 0.432 [0.396, 0.467], q 3.1e-04 | 0.392 [0.347, 0.438], q 2.5e-05 | 0.379 [0.332, 0.428], q 7.9e-06 |
| coccus shaped | 0.545 [0.512, 0.578], q 0.0067 | 0.570 [0.527, 0.612], q 0.0019 | 0.566 [0.520, 0.611], q 0.0055 |
| sphere shaped | 0.553 [0.519, 0.587], q 0.0030 | 0.546 [0.501, 0.591], q 0.044 | 0.555 [0.507, 0.602], q 0.025 |
| arsenate reduction | 0.562 [0.521, 0.601], q 0.0030 | 0.582 [0.529, 0.634], q 0.0026 | 0.585 [0.529, 0.639], q 0.0039 |
| moderately halophilic | 0.574 [0.540, 0.608], q 5.8e-05 | 0.589 [0.544, 0.633], q 2.1e-04 | 0.585 [0.537, 0.632], q 0.0013 |
| facultatively aerobic / cellulose | 0.595 [0.549, 0.640], q 0.0011 | 0.633 [0.572, 0.691], q 8.2e-04 | 0.630 [0.566, 0.692], q 0.0035 |
| cellulose / sphere shaped | 0.407 [0.360, 0.456], q 0.0026 | 0.393 [0.331, 0.457], q 0.015 | 0.395 [0.330, 0.463], q 0.029 |
| cellulose / arsenate reduction | 0.365 [0.304, 0.429], q 0.0011 | 0.336 [0.261, 0.417], q 0.0015 | 0.336 [0.257, 0.421], q 0.0037 |
| cellulose / moderately halophilic | 0.435 [0.386, 0.485], q 0.12 | 0.392 [0.329, 0.458], q 0.016 | 0.391 [0.324, 0.461], q 0.029 |

The other eleven ratios have q >= 0.12 in every population; all 45 rows, with study
counts, are in `disease_control_comparison.tsv`.

Narrowing eligibility moved none of these. Every one of the six fraction directions
and all four significant log-ratios keep their direction and their significance at
q <= 0.05 in both oriented populations, and no further ratio becomes significant.
The shares move by at most 0.015 in the oriented population, which is well inside
the intervals. That the headline numbers barely move is not a reason the systematic
error did not matter: the 66 experiments that left the population included
intervention arms, exposure contrasts and within-disease marker comparisons, and the
claim that the population is disease against healthy is now supportable where before
it was not.

Every direction of the historical analysis is preserved in the disease-versus-healthy
subset and in its strong-evidence restriction, with wider intervals from the smaller
populations. Cellulose-associated annotations are relatively lower on the disease
side, a little more so than in the historical corpus, and the facultatively aerobic,
coccus shaped, sphere shaped, arsenate reduction and moderately halophilic
directions keep their increased-side shares. On this evidence the historical
directions are not an artefact of mixing treatment, time-point and within-disease
contrasts into the corpus, and the disease-versus-healthy subset is a defensible
population for disease-enriched or disease-depleted wording of these annotation
contrasts.

## Limits

The classification reads curated group names, the Group 1 definition and the
description, not the source papers, and is only as good as those fields; the
validation sample is small. Excluded classes are not wrong data - a severity or
treatment contrast answers a different question - and their directions are not
reversed or discarded, only kept out of this population. The oriented subset overlaps
the historical corpus, so agreement between them is a robustness check, not
independent replication. Nothing here measures disease activity, oxygen exposure or
abundance in patients (#179), and none of the procedures has been calibrated (#182).

Known remaining weaknesses, each named by a reviewed experiment rather than assumed:

- **Eligibility is now narrow enough to lose genuine comparisons.** Three sampled
  disease-versus-disease experiments belong in disease-versus-healthy: normal against
  tumour tissue (30584008-4), where a bare "normal" is a comparison word only when
  followed by controls, subjects, individuals, group or samples; "Non-diabetic"
  against "Diabetic" (39417117-26), where the curated condition is "Metabolic
  toxicity" so the negation has no condition to reach; and term against preterm birth
  (39569373-1), where the unaffected group matches no comparison word at all.
- **A shared exposure descriptor still reads as an exposure contrast.** Healthy
  controls against OSCC patients, both described as "Nonsmoking, nondrinking"
  (37542310-2), becomes `treatment_exposure` on a word both groups share.
- **A treated disease group is excluded with its healthy comparison.** Healthy
  controls against Parkinson's patients "treated with entacapone" (31815177-11) is
  `treatment_exposure`. The drug confounds the contrast, so the exclusion is
  defensible, but it is not the reason the rule fires.
- **"Pre" still names a disease stage as well as a time.** "Pre-DM" against type 2
  diabetes (24013136-6) and full-term against "Pre-term" infants (25444008-2) are
  `longitudinal_timepoint`. The first was recorded at the v1 rules and is unfixed.
- **Resolved disease is excluded under a class that describes it wrongly.** Healthy
  controls against recovered COVID-19 cases (33431578-2) and against RSV
  convalescence (39642873-2) are `longitudinal_timepoint` although no time separates
  the two groups.
- **`unresolved` is conservative by design and by accident.** Seven of ten sampled
  `unresolved` experiments belong in some class - an age contrast whose condition is
  "Age at assessment", one cohort sampled in two years, Parkinson's disease by
  duration, two lean populations, two pollution levels, two levodopa regimens, and
  control adults by body site. None belongs in disease-versus-healthy.
- **Sibling experiments of one study can be read inconsistently.** 33294219-1 is
  `other` and 33294219-2 is `disease_vs_disease`, though both compare hypertensive
  Han against hypertensive Tibetan subjects at different altitudes.
- **"Normal (N-like group)" against "Case (OB-like group)" (25202708-1) stays
  `unresolved`.** The reviewer asked for it. The only rule that admits "Normal"
  followed by a parenthesis also admits "Normal (MN) group" against "Medium-dose
  (MMD) group" and "Model high-dose (MHD) group", two dose arms in a mouse study, so
  it would trade one recovered experiment for two new errors. It is left out until
  the dose-arm case has a rule of its own.

## Reproducing

`just contrast-orientation` on the pinned BugSigDBExports v1.3.1 inputs (export
timestamp 2026-04-24_00:41_UTC, `full_dump.csv` SHA-256
`af84bfca7bfe76913a20c1fedbaadbfad4e4bf78605a2cfd7d00eb0c0fdd5788`). The provenance
record gives the code revision, the SHA-256 of the export, profile matrix, experiment
table, classifier source and committed review, the class counts and the review
tally. `just contrast-review-sample` redraws the review sample from its seed.

`tests/test_contrast_orientation.py` holds a regression fixture, with the curated
strings copied verbatim, for every contrast that review
[#291](https://github.com/CultureBotAI/microbial-trait-profiling/issues/291) found
wrongly admitted to or wrongly kept out of the population, alongside the genuine
comparisons the narrowed rules must not take away. Each fixture asserts the expected
class and orientation rather than merely that a wrong answer was avoided. The rest of
the suite covers the reversal rule, the strong-evidence population, both-groups
comparisons, unverified reversals, the randomized-trial branch, baseline dating,
variable conditions, age bands, the class counts, the comparison table, the review
round trip and the command-line entry point.

Every rule was checked by breaking it on purpose and confirming a test failed: 35
such mutations, all caught. Four tests that passed on mutated code were strengthened
(a tautological classifier-version assertion, a study count a fixture could not
distinguish from an experiment count, a review sample whose ordering was not pinned,
and a `validation_status` no test read back), and two pieces of production code that
no test and no corpus experiment exercised were removed rather than left in.

The classifier, the reviews and this description were drafted with the assistance
of a coding agent (Claude Code) and require review by the author before manuscript
use.
