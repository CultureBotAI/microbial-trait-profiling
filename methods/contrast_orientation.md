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
experiment is in the analysed corpus, the classifier version
(`contrast-orientation-v1`) and the BugSigDB vintage. The export holds no Group 0
definition, so every row carries the flag `no_group_0_definition_in_export`; a
Group 0 reading rests on its name alone.

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
   healthy comparison only dates the disease sample.
4. A treatment or exposure in either group name (placebo, probiotic, smokers, users,
   responders, antibiotic) or a treatment or exposure condition (response to a drug,
   antimicrobial agent, caesarean section, breastfeeding): `treatment_exposure`. A
   negated treatment word ("without antiretroviral treatment", "treatment-naive") names
   untreated disease and does not count.
5. One group named by time only, a condition that is a variable (population, age,
   ethnic group, any "measurement"), or two groups that are levels of one variable
   (low against high LDL): `longitudinal_timepoint` or `other`.
6. **Exactly one group reads as a comparison: `disease_vs_healthy`.** Strong
   comparison words are healthy, control, HC, unaffected, uninfected, disease-free
   and volunteers. Weak ones are normal weight, lean, normoglycaemic, negative, a
   test marker negative only in that group ("HPV-" against "HPV+"), and absence of
   the condition itself ("non-IBS", "no asthma", "free from colorectal cancer"). The
   orientation is `native` when the comparison is Group 0. It is `reversed` only when
   Group 1 carries a strong comparison word and Group 0 has disease evidence;
   otherwise the experiment is `unresolved` with the flag
   `possible_reversal_unverified`. A contrast resting on weak evidence only is
   flagged `weak_comparison_evidence`.
7. Both groups read as comparisons: `unresolved`.
8. Both groups read as disease groups (disease, patients, stage, severity, activity,
   subtype, or naming the condition): `disease_vs_disease`.
9. Anything else: `unresolved`.

**Analysis population.** Experiments classed `disease_vs_healthy`. For the one
reversed experiment the increased and decreased profiles are swapped, so that in
the whole population `inc` means higher in the disease group. The historical
results are not overwritten; the oriented population is reported beside them.

## Development and review

The rules were written, reviewed against samples of every class, corrected, and
then frozen before any trait result was computed on the oriented population. The
reviews were done by a coding agent reading the group names, definitions, design
and condition of each sampled experiment, not by the author; they check the rules
against the export's own words, not against the source papers.

**First review** (12 random experiments per class, seed 188, plus all 4 then
classed as reversed). All 12 native disease-versus-healthy and all 12 `other`
experiments were correct and about 11 of 12 disease-versus-disease. Every reversal
was wrong: "abnormal pH" against "normal pH", "Adults" against "normal weight
children", a calprotectin-positive against calprotectin-normal Parkinson's
comparison, and "Asymptomatic" against "No Growth". Other errors: the condition
"Non-alcoholic fatty liver disease" matched the alcohol-exposure rule; a time word
shared by both groups made a stable-against-progressed comparison longitudinal;
"HPV-/HIV-" against "HPV+/HIV-" read both groups as negative; `Control_Chinese`
was missed at the underscore; and absence-of-condition names such as "Non-OSA" were
not recognised. The rules were changed as described above.

**Second review** (12 per class, seed 2026). All 13 disease-versus-healthy
experiments sampled (including the one remaining reversal, polycystic ovary
syndrome against "Healthy control") and all 12 `other` were correct. Errors were
outside the analysis population or conservative: low-against-high LDL as
disease-versus-disease, an HIV-positive group "without antiretroviral treatment" as
treatment, age bands as time points, and a caesarean-section comparison as
longitudinal. Those four rules were fixed and the classifier frozen.

**Validation review of the frozen rules** (seed 9131, no changes made afterwards):

| Class | Sampled | Agreed | Disagreements |
| --- | --- | --- | --- |
| disease_vs_healthy | 25 | 24 | mild/non-severe against severe preeclampsia, read as absence of the condition |
| disease_vs_disease | 8 | 7 | normal against very poor postnatal growth, which is disease versus unaffected |
| treatment_exposure | 8 | 7 | meconium against transitional stool within caesarean births, a time contrast |
| longitudinal_timepoint | 8 | 5 | "Pre-DM" against type 2 diabetes, and "Pre - High risk" patients against rheumatoid arthritis and against healthy controls: "pre" named a disease stage, not a time |
| other | 8 | 8 | - |
| unresolved | 8 | 8 | none wrongly excluded from a class; two conservative misses ("neurotypical" against autistic, normal against cancer within H. pylori-positive) |

Orientation was correct in all 24 correctly classed disease-versus-healthy
experiments. The one disease-versus-healthy error is the known weakness of the
absence rule, which can admit a within-disease severity comparison; a sensitivity
population restricted to strong comparison evidence addresses it below. With 25
sampled experiments, 24 agreements bound the population's error rate only loosely
(exact 95% interval for the agreement share 0.80 to 1.00).

## Retained, reversed, excluded and unresolved

Within the 2,335 analysed experiments (959 studies):

| Class | Orientation | Experiments | Studies |
| --- | --- | --- | --- |
| disease_vs_healthy | native | 1,051 | 589 |
| disease_vs_healthy | reversed | 1 | 1 |
| disease_vs_disease | - | 225 | 154 |
| treatment_exposure | - | 354 | 169 |
| longitudinal_timepoint | - | 204 | 107 |
| other | - | 269 | 102 |
| unresolved | - | 231 | 145 |

Studies appear in more than one class when their experiments differ. 1,052
experiments from 589 studies are retained, 1 of them reversed; 1,052 are excluded
as other classes and 231 are unresolved. Of the retained experiments 216 (108
studies) rest on weak comparison evidence only; 836 experiments from 497 studies
rest on strong evidence.

## The oriented subset against the historical corpus

The same two analyses were run on three populations, each with its own native
families: the six fraction directions (studies scored by the sign of their mean
increased-minus-decreased fraction, exact binomial test, Benjamini-Hochberg over
six, zero-total profiles at fraction 0 as in the historical analysis) and the 15
paired log-ratios (positive-count complete cases, Benjamini-Yekutieli over 15). The
share is of non-tied studies shifting toward the increased side, with an exact 95%
interval for that share only. In the oriented populations "increased" means higher
in the disease group.

| Comparison | Historical (959 studies) | Disease vs healthy (589) | Strong evidence only (497) |
| --- | --- | --- | --- |
| facultatively aerobic | 0.582 [0.550, 0.614], q 4.1e-06 | 0.579 [0.537, 0.621], q 5.6e-04 | 0.566 [0.520, 0.611], q 0.0074 |
| cellulose catabolic process | 0.432 [0.396, 0.467], q 3.1e-04 | 0.401 [0.357, 0.447], q 1.1e-04 | 0.382 [0.334, 0.431], q 1.3e-05 |
| coccus shaped | 0.545 [0.512, 0.578], q 0.0067 | 0.572 [0.530, 0.614], q 0.0011 | 0.563 [0.517, 0.608], q 0.0083 |
| sphere shaped | 0.553 [0.519, 0.587], q 0.0030 | 0.548 [0.504, 0.592], q 0.031 | 0.555 [0.507, 0.602], q 0.025 |
| arsenate reduction | 0.562 [0.521, 0.601], q 0.0030 | 0.567 [0.515, 0.618], q 0.013 | 0.582 [0.526, 0.636], q 0.0074 |
| moderately halophilic | 0.574 [0.540, 0.608], q 5.8e-05 | 0.581 [0.537, 0.624], q 5.6e-04 | 0.579 [0.531, 0.626], q 0.0037 |
| facultatively aerobic / cellulose | 0.595 [0.549, 0.640], q 0.0011 | 0.620 [0.560, 0.677], q 0.0018 | 0.623 [0.559, 0.685], q 0.0041 |
| cellulose / sphere shaped | 0.407 [0.360, 0.456], q 0.0026 | 0.395 [0.335, 0.458], q 0.015 | 0.394 [0.329, 0.461], q 0.024 |
| cellulose / arsenate reduction | 0.365 [0.304, 0.429], q 0.0011 | 0.335 [0.262, 0.415], q 0.0018 | 0.331 [0.253, 0.417], q 0.0041 |
| cellulose / moderately halophilic | 0.435 [0.386, 0.485], q 0.12 | 0.394 [0.332, 0.459], q 0.015 | 0.389 [0.323, 0.459], q 0.024 |

The other eleven ratios have q >= 0.12 in every population; all 45 rows, with study
counts, are in `disease_control_comparison.tsv`.

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

The classification reads curated group names, not the source papers, and is only
as good as those names; the validation sample is small. Excluded classes are not
wrong data - a severity or treatment contrast answers a different question - and
their directions are not reversed or discarded, only kept out of this population.
The oriented subset overlaps the historical corpus, so agreement between them is a
robustness check, not independent replication. Nothing here measures disease
activity, oxygen exposure or abundance in patients (#179), and none of the
procedures has been calibrated (#182).

## Reproducing

`just contrast-orientation` on the pinned BugSigDBExports v1.3.1 inputs (export
timestamp 2026-04-24_00:41_UTC). The provenance record gives the code revision, the
SHA-256 of the export, profile matrix, experiment table and classifier source, and
the class counts. The fixtures in `tests/test_contrast_orientation.py` cover reversed
groups, treatment and time-point contrasts, disease-versus-disease comparisons,
missing definitions and names, the non-alcoholic fatty liver disease guard, test
markers and shared timepoints.

The classifier, the reviews and this description were drafted with the assistance
of a coding agent (Claude Code) and require review by the author before manuscript
use.
