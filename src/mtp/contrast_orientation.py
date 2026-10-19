"""What each BugSigDB experiment compares, and which way its direction points (#188).

BugSigDB records a signature as "Abundance in Group 1" relative to Group 0, so
``inc`` means higher in Group 1. The Table 1 analyses keep that curated direction.
Reading it as "higher in disease" is only right when Group 1 is a disease group and
Group 0 a healthy or unaffected comparison, and nothing upstream checks that.

This module classifies every experiment from the export's own fields into one of six
contrast classes, records the matched evidence, and derives an orientation only where
the evidence supports one. The study design, condition and both group names drive the
classification; the Group 1 definition and the description are read for intervention
evidence alone, because a group label such as "Experimental group" says nothing about
what was done to it while its definition does. Which field supplied each piece of
evidence is named in the evidence string.

* ``disease_vs_healthy`` - exactly one group reads as a healthy, unaffected or
  negative comparison. Orientation is ``native`` when that group is Group 0 and
  ``reversed`` when it is Group 1.
* ``disease_vs_disease`` - both groups read as disease groups (subtype, stage,
  severity, activity) and neither as a comparison group.
* ``treatment_exposure`` - an intervention, exposure or treatment response.
* ``longitudinal_timepoint`` - the same subjects at different times.
* ``other`` - a comparison of a non-disease variable (population, age, lifestyle).
* ``unresolved`` - nothing above is supported, or the evidence conflicts.

The rules are lexical and deliberately conservative. A contrast that could belong
to two classes is not forced into the disease-versus-healthy population. The
export holds no Group 0 definition, so a Group 0 reading rests on its name alone,
and that limitation is recorded rather than filled in.

Eligibility for ``disease_vs_healthy`` is scoped to the actual comparison context:
a negation must scope the condition itself rather than a modifier of it ("non-IBS",
not "non-exacerbated asthma"); a bare "negative" counts only where its subject is
the condition ("HIV negative", not "Her2 negative"); a comparison reading is
withdrawn from a group that also carries disease evidence of its own ("HC-like IBS",
"Non-obstruction tumour group"); and a bare "control" opposite an intervention arm
("Experimental group", "Trial group") is an untreated arm, not a healthy control.
"""

from __future__ import annotations

import argparse
import csv
import difflib
import hashlib
import json
import re
from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import asdict, dataclass, replace
from pathlib import Path

from mtp.paths import DERIVED_DATA_DIR, REPO_ROOT, file_sha256
from mtp.table1_trait_maps import compute_direction_shifts
from mtp.trait_composition import compute_log_ratio_shifts
from mtp.trait_guilds import (
    DEFAULT_EXPERIMENTS_PATH,
    DEFAULT_FULL_DUMP_PATH,
    DEFAULT_PROFILE_PATH,
    TABLE1_TRAITS,
    Profile,
    load_profiles,
    read_bugsigdb_vintage,
    write_stats_report,
)

CLASSIFIER_VERSION = "contrast-orientation-v4"
DEFAULT_OUTPUT_DIR = DERIVED_DATA_DIR / "contrast_orientation"
DIRECTION_CONVENTION = "inc = higher abundance in Group 1 than in Group 0 (BugSigDB)"

DISEASE_VS_HEALTHY = "disease_vs_healthy"
DISEASE_VS_DISEASE = "disease_vs_disease"
TREATMENT = "treatment_exposure"
LONGITUDINAL = "longitudinal_timepoint"
OTHER = "other"
UNRESOLVED = "unresolved"
CLASSES = (DISEASE_VS_HEALTHY, DISEASE_VS_DISEASE, TREATMENT, LONGITUDINAL, OTHER, UNRESOLVED)

NATIVE = "native"
REVERSED = "reversed"
NOT_APPLICABLE = "not_applicable"

CASE_VS_CONTROL_SUBJECT = "case_vs_control_subject"
COMORBIDITY_WITHIN_DISEASE = "comorbidity_within_disease"
PAIRED_TISSUE_AUTOLOGOUS = "paired_tissue_autologous"
RESOLVED_DISEASE = "resolved_disease"
CONTRAST_TYPES = (
    CASE_VS_CONTROL_SUBJECT,
    COMORBIDITY_WITHIN_DISEASE,
    PAIRED_TISSUE_AUTOLOGOUS,
    DISEASE_VS_DISEASE,
    TREATMENT,
    LONGITUDINAL,
    RESOLVED_DISEASE,
    OTHER,
    NOT_APPLICABLE,
)

TIER_1_HEALTHY_NORMAL = "tier_1_healthy_normal"
TIER_2_CONDITION_NEGATIVE_SCREENED = "tier_2_condition_negative_screened"
TIER_3_HOSPITAL_COMORBID_OR_TISSUE = "tier_3_hospital_comorbid_or_tissue"
CONTROL_TIERS = (
    TIER_1_HEALTHY_NORMAL,
    TIER_2_CONDITION_NEGATIVE_SCREENED,
    TIER_3_HOSPITAL_COMORBID_OR_TISSUE,
    NOT_APPLICABLE,
)


def _pattern(*terms: str) -> re.Pattern[str]:
    return re.compile(r"(?<![\w])(?:" + "|".join(terms) + r")(?![\w])", re.IGNORECASE)


#: A group that is plainly well, whatever the condition is. A multi-valued condition
#: needs one of these, because denying one of its values leaves the others standing.
GENERAL_COMPARISON = _pattern(
    r"healthy",
    r"controls?",
    r"[yo]controls?",
    r"hcs?",
    r"unaffected",
    r"non[- ]?diseased",
    r"disease[- ]free",
    (
        r"normal (?:controls?|subjects?|individuals?|groups?|samples?"
        r"|tissues?|mucosa|skin|adjacent\w*)"
    ),
    r"adjacent normal(?:\s+\w+)*",
    r"paracancerous(?:\s+tissues?)?",
    r"volunteers?",
    r"full[- ]?term(?:\s+infants?)?",
    r"term[- ]birth",
    r"typical(?:ly)? developing",
    r"neurotypical",
    r"well[- ]nourished",
)
#: A group that is plainly a healthy or unaffected comparison. Only these can
#: support a reversed orientation.
STRONG_COMPARISON = _pattern(
    r"healthy",
    r"controls?",
    r"[yo]controls?",
    r"hcs?",
    r"unaffected",
    r"non[- ]?diseased",
    r"disease[- ]free",
    r"neoplasia[- ]free",
    r"uninfected",
    r"non[- ]?infected",
    (
        r"normal (?:controls?|subjects?|individuals?|groups?|samples?"
        r"|tissues?|mucosa|skin|adjacent\w*)"
    ),
    r"adjacent normal(?:\s+\w+)*",
    r"adjacent non[- ]?cancerous(?:\s+tissues?)?",
    r"paracancerous(?:\s+tissues?)?",
    r"early stage normal",
    r"non[- ]?lesional(?:\s+skin)?",
    r"volunteers?",
    r"full[- ]?term(?:\s+infants?)?",
    r"term[- ]birth",
    r"typical(?:ly)? developing",
    r"neurotypical",
    r"well[- ]nourished",
)
#: A comparison group named by the absence of the condition's defining feature.
#: A bare "negative" is not here: it is a comparison only where its subject is the
#: condition, which :func:`_condition_negative` decides.
WEAK_COMPARISON = _pattern(
    r"normal[- ]weight",
    r"normal glucose tolerance",
    r"lean",
    r"normo\w+",
    r"nonobese",
)
#: "<subject> negative", with everything before the subject kept so that a negation
#: qualifying one disease inside another can be told apart from a negation of the
#: condition: "HIV negative" and "HIV-1-negative" deny HIV infection, while "Her2
#: negative", "ASCA negative" and "cART negative" name a marker or drug.
CONDITION_NEGATIVE = re.compile(
    "(?P<before>.*?)(?P<subject>[\\w'\\u2019.]+)[^\\w]+negative(?![\\w])", re.IGNORECASE
)
#: People living without HIV, as an initialism. It names no condition, so it counts as a
#: comparison only where HIV is the condition; elsewhere HIV status describes the cohort.
WITHOUT_HIV = _pattern(r"plwoh", r"pwoh")
#: A condition that is HIV itself: the acronym or the virus spelled out. Testing the text,
#: not the initialism, keeps "Helminth infection" or "HOMA-IR" out (#627).
HIV_CONDITION = re.compile(r"(?<![a-z])hiv(?![a-z])|human immunodeficiency virus", re.IGNORECASE)
#: Names for people with HIV that carry no "HIV" token: "PLWH", "PWH", or the virus
#: spelled out without its acronym (#627).
WITH_HIV = re.compile(
    r"(?<![a-z])(?:plwh|pwh)(?![a-z])|human immunodeficiency virus", re.IGNORECASE
)
#: Small Roman numerals, as in "Type II diabetes mellitus": numbers, not acronyms (#627).
ROMAN_NUMERAL = re.compile(r"X{0,3}(?:IX|IV|V?I{0,3})")
#: "Negative for COVID-19 viral RNA": the sign comes before its subject.
NEGATIVE_FOR = re.compile(
    "(?<![\\w])negative\\s+for\\s+((?:[\\w'\\u2019.]+[\\s-]*){1,4})", re.IGNORECASE
)
#: "HC-like IBS", "Healthy Control like Irritable Bowel Syndrome". A group that
#: resembles the comparison group is not the comparison group: this is an
#: irritable-bowel subgroup whose profile happens to look a healthy one.
RESEMBLANCE = re.compile(
    r"(?<![\w])(?:healthy|controls?|hcs?|normal)(?:[\s-]+(?:healthy|controls?|hcs?|normal))*"
    r"[\s-]*like(?![\w])",
    re.IGNORECASE,
)
#: Words that sit between a negator and the condition without changing what is denied:
#: "without any form of Parkinson's disease" denies Parkinson's disease.
NEGATION_FILLER = frozenset(
    {
        "a",
        "an",
        "the",
        "any",
        "all",
        "form",
        "forms",
        "of",
        "known",
        "clinical",
        "clinically",
        "diagnosed",
        "reported",
        "self",
        "history",
        "evidence",
        "signs",
        "symptoms",
        "prior",
        "previous",
        "current",
        "active",
        "overt",
        "apparent",
        # Words saying which case is meant, not how severe it is. A child "without
        # functional constipation" is not constipated, whereas someone with
        # "non-exacerbated asthma" still has asthma, so severity and activity words
        # ("severe", "exacerbated", "cachectic", "ulcerative") are deliberately absent.
        "prevalent",
        "incident",
        "functional",
        "idiopathic",
        "primary",
        "secondary",
        "definite",
        "probable",
        "confirmed",
        "suspected",
        "psychiatric",
    }
)
#: A pathological finding that a negation can modify instead of denying. "non-obstruction
#: tumour group" and "non-endometrioma cysts" still have a tumour and a cyst.
PATHOLOGY_NOUN = _pattern(
    r"cysts?",
    r"tumou?rs?",
    r"lesions?",
    r"polyps?",
    r"nodules?",
    r"ulcers?",
    r"plaques?",
    r"cancers?",
    r"carcinomas?",
    r"metasta\w*",
)
#: An arm label. "Experimental group" or "Trial group" opposite a bare "control" is an
#: untreated arm against an intervention arm, not a disease against a healthy
#: comparison. The trailing "group" is required: "Experimental Autoimmune
#: Encephalomyelitis" and "Trial Making Test" are a disease and a test, not arms.
ARM_GROUP = re.compile(
    r"(?:experimental|trial|intervention|usage)[ -]+groups?(?![\w])", re.IGNORECASE
)
#: Resolved disease, which is neither active disease nor a healthy comparison.
RESOLVED = _pattern(r"recovered", r"recovery", r"convalescen\w*", r"post[- ]?covid")
#: The same subjects sampled at different times. Compound "pre-\w+" / "post-\w+" precede
#: bare "pre" / "post" so that "pre-treatment" and "pre-DM" match as full tokens.
TIMEPOINT = _pattern(
    r"before",
    r"after",
    r"pre-\w+",
    r"post-\w+",
    r"pre",
    r"post",
    r"post\w*partum",
    r"baseline",
    r"endpoint",
    r"follow[- ]?up",
    r"(?:day|week|month|year)s?\s*\d+",
    r"\d+\s*-?\s*(?:day|week|month|year)s?(?![- ]?old)",
    r"(?:one|two|three|four|six|twelve)[- ](?:day|week|month|year)s?",
    r"visit \d+",
    r"time ?points?",
    r"t\d+",
    r"admission",
    r"discharge",
)
#: Pre-/post-/before/after words that name clinical stages, reproductive states or birth
#: cohorts rather than longitudinal sampling of the same subjects.
NON_TEMPORAL_PRE_POST = _pattern(
    r"born\s+(?:before|after)",
    r"pre[- ]?dm",
    r"pre[- ]?diabet\w*",
    r"pre[- ]?term",
    r"pre[- ]?eclamp\w*",
    r"pre[- ]?hypertens\w*",
    r"pre[- ]?malignan\w*",
    r"pre[- ]?neoplas\w*",
    r"pre[- ]?cancer\w*",
    r"pre[- ]?tumou?r\w*",
    r"pre[- ]?lesion\w*",
    r"(?:pre|post)[- ]?menopaus\w*",
    r"pre[- ]?pregnan\w*",
    r"pre[- ]?school\w*",
    r"pre[- ]?clinical\w*",
    r"pre[- ]?symptomatic",
    r"pre[- ]?manifest\w*",
    r"pre[- ]?frail\w*",
    r"pre[- ]?cachec\w*",
    r"pre[- ]?obese",
    r"pre[- ]?adolescen\w*",
    r"pre[- ]?pubert\w*",
    r"(?:pre|post)[- ]?natal\w*",
)
#: Four-digit calendar years, treated as time points only when both groups carry one.
CALENDAR_YEAR = re.compile(
    r"(?<![\w])(?:19[89]\d|20[0-2]\d)(?!\s*(?:mg|g|ml|ug|iu|kcal))(?![\w])",
    re.IGNORECASE,
)
#: An age band, which reads like a duration but dates nobody: "2-12 years" and
#: "aged > 10 years" are how old the subjects were, not when they were sampled.
AGE_BAND = re.compile(
    "\\d+\\s*[-\\u2013]\\s*\\d+\\s*(?:years?|yrs?|months?|mo)(?![\\w])"
    "|aged?\\s*[<>=\\u2264\\u2265\\s]*\\d+\\s*(?:years?|yrs?|months?|mo)?(?![\\w])",
    re.IGNORECASE,
)
#: Time words that, beside a healthy comparison, only date the disease sample.
BASELINE_ONLY = {
    "baseline",
    "admission",
    "pre-treatment",
    "pretreatment",
    "pre-therapy",
    "pre-operation",
    "pre-operative",
    "pre-surgery",
}
#: An intervention, exposure or response.
TREATMENT_TERMS = _pattern(
    r"placebo",
    r"treat\w*",
    r"therap\w*",
    r"antibiotics?",
    r"probiotics?",
    r"prebiotics?",
    r"synbiotics?",
    r"supplement\w*",
    r"drugs?",
    r"medication\w*",
    r"vaccin\w*",
    r"interventions?",
    r"exposed",
    r"exposure",
    r"pollut\w*",
    r"poluuted",
    r"users?",
    r"non[- ]?users?",
    r"smok\w*",
    r"non[- ]?smok\w*",
    r"drinkers?",
    r"non[- ]?drinkers?",
    r"diets?",
    r"dietary",
    r"fib(?:er|re)",
    r"responders?",
    r"non[- ]?responders?",
    r"surgery",
    r"transplant\w*",
    r"chemotherapy",
    r"radiotherapy",
    r"art",
    r"antiretroviral",
    r"metformin",
    r"statins?",
    r"ppis?",
    r"levodopa",
    r"carbidopa",
    r"entacapone",
    r"prednisolon\w*",
)
#: A disease group named by its disease, stage, subtype, severity or activity.
DISEASE_TERMS = _pattern(
    r"patients?",
    r"cases?",
    r"diseases?",
    r"disorders?",
    r"syndromes?",
    r"cancers?",
    r"cancerous",
    r"\w+omas?",
    r"\w+itis",
    r"tumou?rs?",
    r"tumou?ral",
    r"peritumou?ral",
    r"precancer\w*",
    r"polyps?",
    r"lesions?",
    r"lesional",
    r"infect\w*",
    r"stages?",
    r"grades?",
    r"subtypes?",
    r"types?",
    r"active",
    r"remission",
    r"relapse\w*",
    r"flares?",
    r"severe\w*",
    r"mild",
    r"moderate",
    r"acute",
    r"chronic",
    r"advanced",
    r"benign",
    r"malignant",
    r"metasta\w*",
    r"diagnosed",
    r"obes\w*",
    r"overweight",
    r"diabet\w*",
    r"pre[- ]?dm",
    r"pre[- ]?diabet\w*",
    r"pre[- ]?term",
    r"preterm",
    r"hypertensi\w*",
    r"positive",
)
#: Unambiguous disease nouns that can be denied directly in Group 0 when the opposite
#: group has the same disease and the study condition is a broader umbrella label.
SELF_SCOPED_DISEASE = _pattern(
    r"diabet\w*",
    r"obes\w*",
    r"hypertensi\w*",
    r"asthma\w*",
    r"cirrho\w*",
    r"colitis",
    r"constipat\w*",
    r"diarrh\w*",
    r"depress\w*",
    r"anemi\w*",
    r"pneumoni\w*",
    r"osteoporo\w*",
    r"periodontit\w*",
    r"psoria\w*",
    r"dermatit\w*",
)
#: Condition labels that are variables rather than diseases. Tested against each
#: comma-separated value of a multi-valued condition, so "Age" is recognised in
#: "Coinfection,Age"; the condition as a whole is a variable only when every value is.
NON_DISEASE_CONDITION = re.compile(
    r"measurement|status|lifestyle|population|ethnic|(?<![\w])age(?![\w])|^aging$"
    r"|(?<![\w])sex(?![\w])|participation|protocol|location|geograph|body mass index"
    r"|physical activity|socioeconomic|^pregnancy$|altitude|sampling",
    re.IGNORECASE,
)
#: Condition labels that are a treatment or exposure. "Non-alcoholic" is a disease
#: name, not alcohol exposure, hence the guard against a preceding hyphen.
TREATMENT_CONDITION = re.compile(
    r"(?<![\w-])(?:response to|antimicrobial agent|antibiotic\w*|smoking|alcohol\w*|diet"
    r"|drug\w*|vaccin\w*|surgery|transplant\w*|c(?:a)?esarean\w*|breast ?feeding|delivery"
    r"|air pollution|pollution)(?![\w])",
    re.IGNORECASE,
)
#: Demographic, geographic, or anatomical site labels that distinguish two healthy
#: comparison cohorts (such as "Healthy men" vs "Healthy women").
DEMOGRAPHIC_OR_SITE_TERMS = _pattern(
    r"men",
    r"women",
    r"males?",
    r"females?",
    r"boys?",
    r"girls?",
    r"infants?",
    r"children",
    r"adults?",
    r"elderly",
    r"young",
    r"old",
    r"urban",
    r"rural",
    r"high[- ]altitude",
    r"low[- ]altitude",
    r"altitude",
    r"americans?",
    r"ghanaians?",
    r"belgians?",
    r"norwegians?",
    r"norweigans?",
    r"chinese",
    r"japanese",
    r"indians?",
    r"europeans?",
    r"africans?",
    r"han",
    r"tibetans?",
    r"saliva",
    r"plaque",
    r"subgingival",
    r"supragingival",
    r"stool",
    r"feces",
    r"faeces",
    r"fecal",
    r"faecal",
    r"duoden\w*",
    r"ileum",
    r"ileal",
    r"colon\w*",
    r"rectal",
    r"rectum",
    r"mucosa\w*",
    r"lumen\w*",
    r"luminal",
    r"oral",
    r"nasal",
    r"skin",
    r"vaginal",
    r"blood",
    r"serum",
    r"plasma",
    r"urine",
)
#: Paired / autologous tissue comparison groups (normal adjacent tissue vs tumor).
PAIRED_TISSUE_CONTROL = _pattern(
    (
        r"normal (?:tissues?|mucosa|skin|adjacent\w*|colonic\w*"
        r"|gastric\w*|esophageal\w*|oral\w*|breast\w*|lung\w*)"
    ),
    r"adjacent normal(?:\s+\w+)*",
    r"normal adjacent(?:\s+\w+)*",
    r"adjacent non[- ]?cancerous(?:\s+tissues?)?",
    r"non[- ]?cancerous(?:\s+tissues?)",
    r"non[- ]?malignant(?:\s+tissues?)",
    r"non[- ]?tumou?r(?:al)?(?:\s+tissues?)?",
    r"paracancerous(?:\s+tissues?)?",
    r"peritumou?ral(?:\s+tissues?)?",
    r"early stage normal",
    r"non[- ]?lesional(?:\s+skin)?",
    r"unaffected mucosa",
    r"non[- ]?inflamed mucosa",
    r"healthy tissue",
)
PAIRED_TISSUE_DISEASE = _pattern(
    r"tumou?rs?",
    r"tumou?ral",
    r"peritumou?ral",
    r"cancers?",
    r"cancerous",
    r"precancer\w*",
    r"carcinomas?",
    r"adenocarcinomas?",
    r"adenomas?",
    r"polyps?",
    r"lesions?",
    r"lesional",
    r"malignant",
)
#: Disease terms that, when present un-negated in a comparison group's name and distinct
#: from the study condition, indicate a comorbidity contrast within a shared background
#: disease rather than a healthy-normal control.
COMORBIDITY_DISEASE_TERMS = _pattern(
    r"obes\w*",
    r"overweight",
    r"diabet\w*",
    r"mody\d*",
    r"hypertensi\w*",
    r"pneumoni\w*",
    r"copd",
    r"urti",
    r"cirrho\w*",
    r"periodontit\w*",
    r"ibs",
    r"irritable bowel",
    r"ibd",
    r"colitis",
    r"crohn\w*",
    r"glaucoma\w*",
    r"polyps?",
    r"adenomas?",
    r"atopic",
    r"stroke",
    r"covid(?:[- ]?19)?",
    r"sars[- ]?cov[- ]?2",
    r"hemodialysis",
    r"haemodialysis",
    r"transplant\w*",
    r"gvhd",
    r"graft[- ]versus[- ]host",
    r"urinary findings",
)
#: Hospital, surgical, symptomatic or patient-derived non-case control labels.
HOSPITAL_OR_SYMPTOMATIC_CONTROL = _pattern(
    r"hospital controls?",
    r"surgical controls?",
    r"symptomatic controls?",
    r"patient controls?",
    r"control patients?",
    r"(?<!non )(?<!non-)disease controls?",
    r"oncological controls?",
    r"non[- ]?ibd controls?",
    r"unrelated [\w\s-]+ conditions?",
    r"benign",
)
DIRECTLY_NEGATED_DISEASE = re.compile(
    r"\b(?:free (?:from|of)|without|no|non|(?<!positive and )negative(?:\s+for)?)[- ]+"
    r"(?:history of\s+|any\s+|acute\s+|chronic\s+)?[\w'\u2019.-]+"
    r"|\b[\w'\u2019.]+[- ]free\b",
    re.IGNORECASE,
)
TREATED_WITH_CONF = re.compile(r"\btreated\s+with\s+[\w-]+", re.IGNORECASE)
#: A disease group on antiretroviral therapy: "HAART-treated HIV+ patients", "People living
#: with HIV on Antiretroviral Therapy (PLWH on ART)". Opposite an HIV-negative comparison
#: group the therapy is what the cases are on, not what the contrast compares.
ON_ANTIRETROVIRAL = _pattern(
    r"(?:c|ha)?art[- ]treated",
    r"(?:on|taking|receiving|under)\s+(?:(?:c|ha)?art|(?:combination\s+)?antiretroviral\w*)",
)
#: Intervention words unambiguous enough to be trusted in free text, so that the Group
#: 1 definition and the description can settle what a bare arm label leaves open: a
#: misspelt "Metaformin" group whose description says metformin, an "Experimental
#: group" defined as acupoint massage therapy, a "Meditation group" of monks who
#: meditate. Deliberately narrower than :data:`TREATMENT_TERMS`, which holds words
#: ("cases", "users", "diet", "art") that are ordinary in a case-control description.
#: Cancer therapy is deliberately absent: a case-control description routinely dates its
#: samples against chemotherapy or radiotherapy ("before chemotherapy", "treatment-naive")
#: without the experiment being about either, and including them moved eight genuine
#: healthy-control contrasts out of the population.
INTERVENTION_EVIDENCE = _pattern(
    r"placebo",
    r"metformin",
    r"statins?",
    r"acupoints?",
    r"acupunctur\w*",
    r"massages?",
    r"meditat\w*",
    r"fa?ecal microbiota transplant\w*",
    r"intra\w*gastric administration",
)
#: A treatment word that is negated or pre-baseline names untreated disease, not a treatment.
NEGATED_TREATMENT = re.compile(
    r"(?:without|no|not|never|free\s+of|prior\s+to)\s+"
    r"(?:taking\s+|using\s+|receiving\s+|on\s+|received\s+|treated\s+with\s+|exposed\s+to\s+)?"
    r"(?:[\w-]+\s+){0,3}"
    r"(?:treat\w*|therap\w*|antibiotic\w*|antiretroviral\w*|medication\w*|drugs?|probiotic\w*"
    r"|prebiotic\w*|supplement\w*|ppis?|statins?|metformin|chemotherapy|radiotherapy|art\b)"
    r"|(?:treatment|therapy|drug|antibiotic|medication|art|c?art)[- ](?:naive|free)"
    r"|untreated|unmedicated|non[- ]?medicated|unexposed|non[- ]?exposed"
    r"|pre[- ]?(?:treat\w*|therap\w*|operat\w*|surg\w*)",
    re.IGNORECASE,
)
#: Both groups named as low/high levels of one variable.
LEVEL_SPLIT = _pattern(
    r"low\w*", r"high\w*", r"levels?", r"above", r"below", r"quartiles?", r"tertiles?"
)
MULTI_DESIGN_TREATMENT = re.compile(r"randomized controlled trial", re.IGNORECASE)
MULTI_DESIGN_TIME = re.compile(r"time series|longitudinal", re.IGNORECASE)
#: A sign counts as a test result where the name or a sub-name ends after it, with or
#: without a space before it ("HPV-", "HPV +"). A minus sign or an en dash that touches its
#: marker also counts before a word ("HIV\u2013 individuals"), because neither can be the
#: hyphen that is as often a separator there ("Periodontitis- Shallow pocket", "BV-
#: Positive") as a minus. The space before a sign is matched outside the sign's group, so
#: the group holds only the sign (#596).
MARKER = re.compile(
    "([a-z][a-z0-9.]*)"
    "(?:\\s+(?=[+\\-\u2212\u2013]\\s*(?:$|[/),;(])))?"
    "([+\\-\u2212\u2013](?=\\s*(?:$|[/),;(]))|[\u2212\u2013](?=\\s+\\w))",
    re.IGNORECASE,
)
#: Words naming the same disease category, so that a negation can deny a condition
#: stated in other words: "noncancer individuals" against oral squamous cell carcinoma,
#: or "Non-obese" against "Obesity".
CONDITION_SYNONYMS = (
    frozenset(
        {
            "cancer",
            "cancerous",
            "carcinoma",
            "neoplasm",
            "neoplasia",
            "tumour",
            "tumor",
            "adenocarcinoma",
            "sarcoma",
            "malignancy",
        }
    ),
    frozenset({"obese", "obesity"}),
    frozenset({"diabetes", "diabetic", "diabetics"}),
    frozenset({"hypertension", "hypertensive"}),
)
CONDITION_STOPWORDS = {"disease", "disorder", "syndrome", "infection", "and", "with", "of", "the"}


@dataclass(frozen=True)
class ExperimentRecord:
    """The export's fields for one experiment, verbatim."""

    study_id: str
    experiment_id: str
    pmid: str
    doi: str
    study_design: str
    condition: str
    efo_id: str
    body_site: str
    uberon_id: str
    host_species: str
    group_0_name: str
    group_1_name: str
    group_1_definition: str
    group_0_sample_size: str
    group_1_sample_size: str
    description: str


@dataclass(frozen=True)
class Classification:
    contrast_class: str
    orientation: str
    evidence: str
    flags: str
    contrast_type: str = NOT_APPLICABLE
    control_tier: str = NOT_APPLICABLE


def _clean(value: str) -> str:
    # Curated descriptions and definitions contain line breaks, which would split a
    # manifest row across physical lines: valid CSV, but grep, cut and head all read it
    # wrongly. Whitespace is collapsed so every cell holds one line.
    value = " ".join(value.split())
    return "" if value.upper() == "NA" else value


def _normalised(name: str) -> str:
    """Underscores separate words in group names ("Control_Chinese"); slashes stay tokens."""
    return re.sub(r"\s*/\s*", " / ", name.replace("_", " "))


def _matches(pattern: re.Pattern[str], text: str) -> list[str]:
    return sorted({match.group(0).lower() for match in pattern.finditer(text)})


def _markers(name: str) -> dict[str, str]:
    """Test markers and their sign, as in "HPV-", "HIV +" or "HIV\u2013 individuals"."""
    return {marker.lower(): "+" if sign == "+" else "-" for marker, sign in MARKER.findall(name)}


def _condition_words(condition: str) -> tuple[set[str], set[str]]:
    """Six-letter stems of the condition's words, and its initialisms.

    A condition already written as an acronym ("HIV infection") contributes that acronym
    as a stem, since it is too short for the six-letter rule and its initialism ("hi")
    describes nothing. Possessive suffixes ("Crohn's", "Parkinson's") are stripped
    before tokenising so they do not produce a spurious 's' initial.
    """
    stems: set[str] = set()
    initials: set[str] = set()
    for part in condition.split(","):
        cleaned = re.sub(r"['\u2019]s\b", "", part.lower())
        words = [w for w in re.findall(r"[a-z]+", cleaned) if w not in CONDITION_STOPWORDS]
        stems.update(w[:6] for w in words if len(w) >= 4)
        stems.update(w for w in words if 2 <= len(w) <= 3)
        full = [w for w in re.findall(r"[a-z]+", cleaned) if w not in {"and", "of", "the", "with"}]
        if len(full) >= 2:
            initials.add("".join(w[0] for w in full))
    return stems, initials


def _same_disease_category(word: str, condition: str) -> bool:
    """Do a word and the condition name the same disease category?

    BugSigDB states a condition as precisely as the paper does, so a group named by the
    absence of the broad category ("noncancer") denies a condition named as a specific
    instance of it ("Squamous cell carcinoma") even though the two share no letters.
    """
    words = set(re.findall(r"[a-z]+", condition.lower()))
    return any(word.lower() in group and words & group for group in CONDITION_SYNONYMS)


def _mentions_condition(text: str, condition: str) -> bool:
    stems, initials = _condition_words(condition)
    tokens = re.findall(r"[a-z0-9]+", text.lower())
    return (
        any(t[:6] in stems for t in tokens if len(t) >= 4)
        or any(t in stems for t in tokens if len(t) <= 3)
        or any(_same_disease_category(t, condition) for t in tokens)
        # Either direction: the condition's initialism can be the longer string, as
        # "hpv" against the initialism of "Human papilloma virus infection".
        or any(
            t.startswith(i) or (len(t) >= 3 and i.startswith(t))
            for t in tokens
            for i in initials
            if len(i) >= 2
        )
    )


def _names_condition(name: str, condition: str) -> bool:
    """Does a group name name the condition? For HIV, "PLWH", "PWH" and the spelled-out
    virus count too, which :func:`_mentions_condition` cannot see (#627)."""
    return _mentions_condition(name, condition) or bool(
        HIV_CONDITION.search(condition) and WITH_HIV.search(name)
    )


def _acronym_spelled_out(words: list[str], condition: str) -> list[str]:
    """Leading words that spell out an acronym of the condition, replaced by the acronym.

    A condition named by an acronym ("HIV infection") is denied as often in full as by
    the acronym: "without Human immunodeficiency virus (HIV)" lacks HIV. Only a word the
    condition writes in capitals is an acronym, so the "to" of "Response to antiviral
    drug" or the "Non" of "Non-alcoholic fatty liver disease" is never spelled out.
    """
    acronyms = {
        word.lower()
        for word in re.findall(r"(?<![\w])[A-Z]{2,5}(?![a-z])", condition)
        if not ROMAN_NUMERAL.fullmatch(word)
    }
    for acronym in sorted(acronyms, key=len, reverse=True):
        if "".join(word[0] for word in words[: len(acronym)]) == acronym:
            return [acronym, *words[len(acronym) :]]
    return words


def _negation_scopes_condition(phrase: str, condition: str) -> bool:
    """Does a negated phrase deny the condition itself rather than a modifier of it?

    "non-IBS" and "free from colorectal cancer" deny the condition. "non-exacerbated
    asthma", "non-obstruction tumour group" and "non-endometrioma cysts" deny a
    modifier or a neighbouring finding: the subjects still have asthma, a tumour and a
    cyst. The condition must therefore be named by the first content word of the
    phrase, and every word after it must either name the condition too or be a generic
    cohort noun.
    """

    def denies(word: str) -> bool:
        return _mentions_condition(word, condition) or _same_disease_category(word, condition)

    words = re.findall(r"[\w'\u2019.]+", phrase.lower())
    while words and words[0] in NEGATION_FILLER:
        words.pop(0)
    words = _acronym_spelled_out(words, condition)
    if not words or not denies(words[0]):
        return False
    rest = list(words[1:])
    # A condition stopword ("disease", "syndrome") continues the condition's name even
    # though it is too generic to be one of its stems - but only after the condition has
    # already been named, so "without infection" alone never denies an unrelated cancer.
    while rest and (denies(rest[0]) or rest[0] in CONDITION_STOPWORDS):
        rest.pop(0)
    # "non-obstruction tumour group": the negation qualifies the tumour, not the patient.
    return not (rest and _matches(PATHOLOGY_NOUN, rest[0]))


NEGATORS = re.compile(
    "(?<![\\w])(?:free (?:from|of)[- ]+|non[- ]*|(?:no|without)[- ]+)"
    "((?:[\\w'\u2019.]+[\\s-]*){1,6})",
    re.IGNORECASE,
)


def _absence_of_condition(name: str, condition: str) -> list[tuple[str, tuple[int, int]]]:
    """ "non-IBS", "no asthma", "without any form of Parkinson's" - absence of this condition.

    Each hit is returned with its span, so that a condition named only inside a negation
    is not later read as the group having it.
    """
    found = []
    # Separator-insensitive, so neither "Non-alcoholic" nor "nonalcoholic" reads as the
    # absence of non-alcoholic fatty liver disease.
    folded = re.sub(r"[\W_]+", "", condition.lower())
    for match in NEGATORS.finditer(name):
        # Either way round, because the match window can stop inside the condition's
        # name or run past its end: "non-alcoholic fatty liver" and "non-alcoholic fatty
        # liver disease patients" both name the condition rather than deny it.
        squashed = re.sub(r"[\W_]+", "", match.group(0).lower())
        if squashed and (folded.startswith(squashed) or squashed.startswith(folded)):
            continue
        if _negation_scopes_condition(match.group(1), condition):
            found.append((match.group(0).strip().lower(), match.span()))
    return found


def _condition_negative(name: str, condition: str) -> list[tuple[str, tuple[int, int]]]:
    """ "HIV negative" - the sign reads as absence only when its subject is the condition.

    What is denied is whatever the sign's subject is, and the subject can be a word or
    two back: "group B Streptococcus culture-negative" and "SARS Cov-2 Saliva Negative"
    deny the organism. So the whole run of words before the sign is the subject, and the
    condition must be named somewhere in it. Any disease word in that run disqualifies
    it, because a disease named alongside the sign is one the group has: "crohn disease
    ASCA negative patients", "Parkinson's disease dementia negative", "HIV infected cART
    negative" and "Type 2 diabetes mellitus negative and periodontitis positive" each
    qualify one disease inside another. "Her2 negative tumor" and "Metaformin negative"
    never name the condition in the first place.
    """
    found = [
        (match.group(0).strip().lower(), match.span())
        for match in NEGATIVE_FOR.finditer(name)
        if _negation_scopes_condition(match.group(1), condition)
    ]
    for match in CONDITION_NEGATIVE.finditer(name):
        subject = f"{match.group('before')} {match.group('subject')}"
        if not _mentions_condition(subject, condition) or _matches(DISEASE_TERMS, subject):
            continue
        found.append((match.group(0).strip().lower(), match.span()))
    return found


def _self_scoped_absence(name: str, other: str) -> list[tuple[str, tuple[int, int]]]:
    """ "non-diabetic" vs "diabetic" when the condition field names a broader category.

    When one group explicitly negates the exact disease noun present in the other group
    (and carries no other disease term outside that negation), it acts as the condition-
    negative comparison group even if the study's Condition field uses a broader EFO
    term (for example, "glucose metabolism disease" for a diabetic vs non-diabetic contrast).
    """
    other_clean = NEGATORS.sub(" ", other)
    other_diseases = set(_matches(SELF_SCOPED_DISEASE, other_clean))
    if not other_diseases:
        return []
    found = []
    for match in NEGATORS.finditer(name):
        phrase_words = re.findall(r"[a-z]+", match.group(1).lower())
        while phrase_words and phrase_words[0] in NEGATION_FILLER:
            phrase_words.pop(0)
        if phrase_words and phrase_words[0] in other_diseases:
            outside = name[: match.start()] + " " + name[match.end() :]
            if not _matches(DISEASE_TERMS, outside):
                found.append((match.group(0).strip().lower(), match.span()))
    return found


def _condition_values(condition: str) -> list[str]:
    """The condition's comma-separated values, which BugSigDB allows more than one of."""
    return [value for value in (part.strip() for part in condition.split(",")) if value]


def _variable_condition(condition: str) -> bool:
    """Is every value of the condition a variable rather than a disease?

    "Age" alone is a variable; "Coinfection,Age" is not, because coinfection is a
    disease, and an experiment that compares coinfected against uncoinfected subjects is
    not a comparison of ages.
    """
    values = _condition_values(condition)
    return bool(values) and all(NON_DISEASE_CONDITION.search(value) for value in values)


def _age_bands(name: str) -> list[str]:
    return sorted({match.group(0).strip().lower() for match in AGE_BAND.finditer(name)})


def _timepoints(name: str) -> set[str]:
    """Time words in a group name, excluding those inside an age band or clinical stage.

    "2-12 years" against "18-62 years" and "Controls aged > 10 years" look like
    durations to :data:`TIMEPOINT` but say how old the subjects were, and "pre-DM" or
    "pre-hypertension" name a clinical stage rather than a follow-up visit.
    """
    spans = [match.span() for match in AGE_BAND.finditer(name)] + [
        match.span() for match in NON_TEMPORAL_PRE_POST.finditer(name)
    ]
    return {
        match.group(0).lower()
        for match in TIMEPOINT.finditer(name)
        if not any(start <= match.start() and match.end() <= end for start, end in spans)
    }


def _stated_intervention(record: ExperimentRecord) -> list[str]:
    """Interventions that a group is named after, spelled correctly by the curator.

    A group label can name an intervention the term list cannot match - "Metaformin
    positive" misspells metformin, "Meditation group" is a practice no drug list holds -
    while the Group 1 definition and the description spell it out. So the definition and
    description supply the vocabulary and the group names decide whether it applies:
    each intervention found in the prose must also be approximately present in a group
    name. Without that, prose naming a cohort's history ("had history of taking
    anti-diabetic medication ... metformin" beside "Healthy controls") would turn a
    disease-versus-healthy contrast into a drug trial. Negated mentions are stripped
    first, so a description that rules out recent antibiotics is not read as a trial.
    """
    prose = NEGATED_TREATMENT.sub(
        " ", f"{_clean(record.group_1_definition)} | {_clean(record.description)}"
    )
    tokens = re.findall(r"[a-z]{4,}", f"{record.group_0_name} {record.group_1_name}".lower())
    return [
        term
        for term in _matches(INTERVENTION_EVIDENCE, prose)
        if difflib.get_close_matches(term, tokens, n=1, cutoff=0.85)
    ]


def _comparison(name: str, other: str, condition: str) -> tuple[list[str], list[str], list[str]]:
    """Strong and weak evidence that ``name`` is the comparison group, and what withdraws it.

    A group that only resembles the comparison group is not it: "HC-like IBS (Healthy
    Control like Irritable Bowel Syndrome)" is an irritable-bowel subgroup whose profile
    looks like a healthy one. A condition with several values needs a group that is well
    in general, since denying one value leaves the others standing. A withdrawn reading
    is returned as the third element so the reason is recorded rather than silently
    dropped.
    """
    resemblance = RESEMBLANCE.sub(" ", name)
    strong = _matches(STRONG_COMPARISON, resemblance)
    absence = (
        _absence_of_condition(name, condition)
        + _condition_negative(name, condition)
        + _self_scoped_absence(name, other)
    )
    weak = _matches(WEAK_COMPARISON, resemblance) + [text for text, _span in absence]
    if HIV_CONDITION.search(condition):
        weak += _matches(WITHOUT_HIV, name)
    mine, theirs = _markers(name), _markers(other)
    negative = [
        m
        for m, sign in mine.items()
        if sign == "-" and theirs.get(m) != "-" and _mentions_condition(m, condition)
    ]
    if negative:
        weak.append(f"test-negative:{','.join(sorted(negative))}")
    if not (strong or weak):
        resembled = _matches(STRONG_COMPARISON, name) + _matches(WEAK_COMPARISON, name)
        return [], [], resembled
    # A condition with several values needs a group that is well in general. "Epstein-Barr
    # virus non-infected (EBV- GC)" lacks the virus and has the gastric cancer; "young
    # HIV-infected group without coinfection" lacks the coinfection and has the HIV.
    if len(_condition_values(condition)) > 1 and not _matches(GENERAL_COMPARISON, resemblance):
        return [], [], strong + weak
    return strong, weak, []


def _disease_evidence(name: str, condition: str) -> list[str]:
    found = _matches(DISEASE_TERMS, name)
    if _mentions_condition(name, condition):
        found.append("names the condition")
    if any(sign == "+" for sign in _markers(name).values()):
        found.append("test-positive")
    return found


def _classify_control_tier(
    control_group_name: str,
    disease_group_name: str,
    condition: str,
    strong_control: list[str],
    flags: list[str],
) -> tuple[str, str]:
    """Stratify a ``disease_vs_healthy`` contrast into ``(contrast_type, control_tier)``."""
    # 1. Paired / autologous tissue comparison (e.g. tumour vs adjacent normal mucosa).
    if _matches(PAIRED_TISSUE_CONTROL, control_group_name) or (
        ("tissue" in control_group_name.lower() or "tissues" in control_group_name.lower())
        and _matches(PAIRED_TISSUE_DISEASE, disease_group_name)
    ):
        return PAIRED_TISSUE_AUTOLOGOUS, TIER_3_HOSPITAL_COMORBID_OR_TISSUE

    # 2. Background comorbidity in the control group (disease A + B vs disease A without B).
    if re.search(r"non[- ]+stroke\s+associated\s+pneumonia", control_group_name, re.IGNORECASE):
        return COMORBIDITY_WITHIN_DISEASE, TIER_3_HOSPITAL_COMORBID_OR_TISSUE

    ctrl_markers = _markers(control_group_name)
    if any(sign == "+" for sign in ctrl_markers.values()):
        return COMORBIDITY_WITHIN_DISEASE, TIER_3_HOSPITAL_COMORBID_OR_TISSUE

    residual = DIRECTLY_NEGATED_DISEASE.sub(" ", control_group_name)
    comorbid_hits = [
        d
        for d in _matches(COMORBIDITY_DISEASE_TERMS, residual)
        if not _mentions_condition(d, condition)
    ]
    if comorbid_hits:
        return COMORBIDITY_WITHIN_DISEASE, TIER_3_HOSPITAL_COMORBID_OR_TISSUE

    if re.search(
        r"\b(?:no|without|non)[- ]+cirrhosis\b", control_group_name, re.IGNORECASE
    ) and not _mentions_condition("cirrhosis", condition):
        return COMORBIDITY_WITHIN_DISEASE, TIER_3_HOSPITAL_COMORBID_OR_TISSUE

    # 3. Hospital, surgical, symptomatic or patient controls, or treated-disease confound.
    if (
        _matches(HOSPITAL_OR_SYMPTOMATIC_CONTROL, control_group_name)
        or "treated_disease_group" in flags
        or (not strong_control and re.search(r"\bpatients?\b", control_group_name, re.IGNORECASE))
        or re.search(r"\bhiv\s+negative\s+patients?\b", control_group_name, re.IGNORECASE)
    ):
        return CASE_VS_CONTROL_SUBJECT, TIER_3_HOSPITAL_COMORBID_OR_TISSUE

    # 4. Subject-level controls: Tier 1 (explicit healthy/normal) vs Tier 2 (screened negative).
    if strong_control:
        return CASE_VS_CONTROL_SUBJECT, TIER_1_HEALTHY_NORMAL
    return CASE_VS_CONTROL_SUBJECT, TIER_2_CONDITION_NEGATIVE_SCREENED


def classify(record: ExperimentRecord) -> Classification:
    """Assign one contrast class, orientation, contrast type and control tier."""
    g0 = _normalised(_clean(record.group_0_name))
    g1 = _normalised(_clean(record.group_1_name))
    design, condition = _clean(record.study_design), _clean(record.condition)
    flags: list[str] = ["no_group_0_definition_in_export"]
    if not _clean(record.group_1_definition):
        flags.append("no_group_1_definition")
    if MULTI_DESIGN_TIME.search(design):
        flags.append("design_longitudinal")
    if MULTI_DESIGN_TREATMENT.search(design):
        flags.append("design_rct")
    if "meta-analysis" in design.lower():
        flags.append("design_meta_analysis")

    def result(
        contrast_class: str,
        orientation: str,
        evidence: str,
        contrast_type: str | None = None,
        control_tier: str = NOT_APPLICABLE,
    ) -> Classification:
        ctype = (
            contrast_type
            if contrast_type is not None
            else (contrast_class if contrast_class != UNRESOLVED else NOT_APPLICABLE)
        )
        return Classification(
            contrast_class, orientation, evidence, ",".join(flags), ctype, control_tier
        )

    if not g0 or not g1:
        return result(UNRESOLVED, NOT_APPLICABLE, "a group name is missing")

    strong_0, weak_0, withdrawn_0 = _comparison(g0, g1, condition)
    strong_1, weak_1, withdrawn_1 = _comparison(g1, g0, condition)

    # When both groups carry the same generic label ("normal-weight controls" vs
    # "obese controls") and only one carries specific condition-negative evidence,
    # strip the shared strong token so the specific comparison group wins.
    shared_strong = set(strong_0) & set(strong_1)
    if shared_strong and ((weak_0 and not weak_1) or (weak_1 and not weak_0)):
        strong_0 = [s for s in strong_0 if s not in shared_strong]
        strong_1 = [s for s in strong_1 if s not in shared_strong]

    control_0, control_1 = strong_0 + weak_0, strong_1 + weak_1
    if withdrawn_0 or withdrawn_1:
        flags.append("comparison_reading_withdrawn")
    variable_condition = _variable_condition(condition)

    resolved = _matches(RESOLVED, f"{g0} | {g1}")
    if resolved:
        return result(
            LONGITUDINAL,
            NOT_APPLICABLE,
            f"a group is resolved disease: {resolved}",
            contrast_type=RESOLVED_DISEASE,
        )

    time_0, time_1 = _timepoints(g0), _timepoints(g1)
    years_0 = {m.group(0).lower() for m in CALENDAR_YEAR.finditer(g0)}
    years_1 = {m.group(0).lower() for m in CALENDAR_YEAR.finditer(g1)}
    if years_0 and years_1:
        time_0 = time_0 | years_0
        time_1 = time_1 | years_1
    if _age_bands(g0) or _age_bands(g1):
        flags.append("age_band_not_timepoint")
    one_sided_time = False
    if time_0 == time_1 and time_0:
        flags.append("shared_timepoint")
    elif time_0 and time_1:
        return result(
            LONGITUDINAL,
            NOT_APPLICABLE,
            f"group names differ in time: {sorted(time_0)} / {sorted(time_1)}",
        )
    elif time_0 or time_1:
        dated_only = (control_0 and time_1 <= BASELINE_ONLY and not time_0) or (
            control_1 and time_0 <= BASELINE_ONLY and not time_1
        )
        if dated_only:
            flags.append("disease_sample_dated_at_baseline")
        else:
            one_sided_time = True

    untreated_0 = NEGATED_TREATMENT.sub(" ", g0)
    untreated_1 = NEGATED_TREATMENT.sub(" ", g1)
    if f"{untreated_0} | {untreated_1}" != f"{g0} | {g1}":
        flags.append("negated_treatment_term")
    treat_0 = _matches(TREATMENT_TERMS, untreated_0)
    treat_1 = _matches(TREATMENT_TERMS, untreated_1)
    if treat_0 and treat_0 == treat_1 and (control_0 or control_1):
        flags.append("shared_treatment_term")
        treatment_terms: list[str] = []
    else:
        treatment_terms = sorted(set(treat_0) | set(treat_1))

    # A healthy control compared against a disease group that happens to be on therapy
    # ("Controls" vs "HIV Patients treated with Maraviroc", "HIV-" vs "People Living With
    # HIV on Antiretroviral therapy") is a confounded disease-vs-healthy contrast (Tier 3),
    # not a treatment-arm trial. Treated against untreated people who all have the disease
    # has no comparison group and stays a treatment contrast. A treated group also named by
    # time keeps the treatment reading: without it only the time word would be left, and
    # in "taking ART for more than 1 year" that word is a duration, not a later visit. The
    # guard reads a time word in either group. "on HAART" or "on cART" counts as therapy
    # even though TREATMENT_TERMS reads only a whole "ART" (#627).
    if (
        (treatment_terms or ON_ANTIRETROVIRAL.search(g1))
        and control_0
        and not control_1
        and not treat_0
        and not one_sided_time
        and _names_condition(g1, condition)
        and (TREATED_WITH_CONF.search(g1) or ON_ANTIRETROVIRAL.search(g1))
    ):
        flags.append("treated_disease_group")
        treatment_terms = []

    if treatment_terms:
        if control_0 or control_1:
            flags.append("comparison_group_present")
        return result(
            TREATMENT,
            NOT_APPLICABLE,
            f"group names mention treatment or exposure: {treatment_terms}",
        )
    arms = sorted(
        {match.group(0).lower() for name in (g0, g1) for match in ARM_GROUP.finditer(name)}
    )
    if arms:
        if control_0 or control_1:
            flags.append("comparison_group_present")
        return result(TREATMENT, NOT_APPLICABLE, f"a group is an intervention arm: {arms}")
    if TREATMENT_CONDITION.search(condition):
        return result(
            TREATMENT, NOT_APPLICABLE, f"condition is a treatment or exposure: {condition!r}"
        )
    stated = _stated_intervention(record)
    if stated:
        if control_0 or control_1:
            flags.append("comparison_group_present")
        return result(
            TREATMENT,
            NOT_APPLICABLE,
            f"the Group 1 definition or description names an intervention: {stated}",
        )
    if one_sided_time:
        return result(
            LONGITUDINAL,
            NOT_APPLICABLE,
            f"one group is named by time: {sorted(time_0)} / {sorted(time_1)}",
        )
    if variable_condition and not (control_0 or control_1):
        return result(OTHER, NOT_APPLICABLE, f"condition is a variable: {condition!r}")
    if "design_rct" in flags and not (control_0 or control_1):
        return result(
            TREATMENT, NOT_APPLICABLE, "randomized controlled trial without a comparison-group name"
        )
    if _matches(LEVEL_SPLIT, g0) and _matches(LEVEL_SPLIT, g1) and not (strong_0 or strong_1):
        return result(OTHER, NOT_APPLICABLE, "both groups are levels of one variable")

    if control_0 and control_1:
        if variable_condition:
            return result(
                OTHER,
                NOT_APPLICABLE,
                f"comparison group named, but the condition is a variable: {condition!r}",
            )
        demo_0 = set(_matches(DEMOGRAPHIC_OR_SITE_TERMS, g0))
        demo_1 = set(_matches(DEMOGRAPHIC_OR_SITE_TERMS, g1))
        if set(control_0) == set(control_1) and (demo_0 ^ demo_1):
            return result(
                OTHER,
                NOT_APPLICABLE,
                (
                    "both groups are comparisons differing by cohort or site: "
                    f"{sorted(demo_0 ^ demo_1)}"
                ),
            )
        return result(
            UNRESOLVED,
            NOT_APPLICABLE,
            f"both groups read as comparisons: {control_0} / {control_1}",
        )
    if control_0 or control_1:
        if variable_condition:
            return result(
                OTHER,
                NOT_APPLICABLE,
                f"comparison group named, but the condition is a variable: {condition!r}",
            )
        if control_0:
            if not strong_0:
                flags.append("weak_comparison_evidence")
            ctype, ctier = _classify_control_tier(g0, g1, condition, strong_0, flags)
            if ctype == PAIRED_TISSUE_AUTOLOGOUS:
                flags.append("paired_tissue_comparison")
            elif ctype == COMORBIDITY_WITHIN_DISEASE:
                flags.append("comorbidity_control")
            return result(
                DISEASE_VS_HEALTHY,
                NATIVE,
                f"Group 0 reads as a comparison: {control_0}",
                contrast_type=ctype,
                control_tier=ctier,
            )
        disease_0 = _disease_evidence(g0, condition)
        if strong_1 and disease_0:
            ctype, ctier = _classify_control_tier(g1, g0, condition, strong_1, flags)
            if ctype == PAIRED_TISSUE_AUTOLOGOUS:
                flags.append("paired_tissue_comparison")
            elif ctype == COMORBIDITY_WITHIN_DISEASE:
                flags.append("comorbidity_control")
            return result(
                DISEASE_VS_HEALTHY,
                REVERSED,
                f"Group 1 reads as a comparison: {strong_1}; Group 0 as disease: {disease_0}",
                contrast_type=ctype,
                control_tier=ctier,
            )
        flags.append("possible_reversal_unverified")
        return result(
            UNRESOLVED,
            NOT_APPLICABLE,
            f"Group 1 may be a comparison ({control_1}); the evidence does not support reversing",
        )

    if variable_condition:
        return result(OTHER, NOT_APPLICABLE, f"condition is a variable: {condition!r}")
    disease_0, disease_1 = _disease_evidence(g0, condition), _disease_evidence(g1, condition)
    if disease_0 and disease_1:
        return result(
            DISEASE_VS_DISEASE,
            NOT_APPLICABLE,
            f"both groups read as disease groups: {disease_0} / {disease_1}",
        )
    return result(UNRESOLVED, NOT_APPLICABLE, "no comparison, disease, time or treatment reading")


def read_experiment_records(full_dump_path: Path) -> dict[tuple[str, str], ExperimentRecord]:
    """One record per experiment from the export, keyed as the ingest keys them.

    A BugSigDB ID ``bsdb:11/1/1`` is study ``bsdb:11``, experiment ``bsdb:11-1``.
    Experiment-level fields repeat on every signature row, so the first row is kept
    and a later row that disagrees on them is an error.
    """
    records: dict[tuple[str, str], ExperimentRecord] = {}
    with full_dump_path.open(encoding="utf-8", newline="") as handle:
        first = handle.readline()
        if not first.startswith("#"):
            handle.seek(0)
        for row in csv.DictReader(handle):
            parts = row["BSDB ID"].removeprefix("bsdb:").split("/")
            if len(parts) != 3:
                continue
            study, experiment = f"bsdb:{parts[0]}", f"bsdb:{parts[0]}-{parts[1]}"
            record = ExperimentRecord(
                study_id=study,
                experiment_id=experiment,
                pmid=_clean(row["PMID"]),
                doi=_clean(row["DOI"]),
                study_design=_clean(row["Study design"]),
                condition=_clean(row["Condition"]),
                efo_id=_clean(row["EFO ID"]),
                body_site=_clean(row["Body site"]),
                uberon_id=_clean(row["UBERON ID"]),
                host_species=_clean(row["Host species"]),
                group_0_name=_clean(row["Group 0 name"]),
                group_1_name=_clean(row["Group 1 name"]),
                group_1_definition=_clean(row["Group 1 definition"]),
                group_0_sample_size=_clean(row["Group 0 sample size"]),
                group_1_sample_size=_clean(row["Group 1 sample size"]),
                description=_clean(row["Description"]),
            )
            key = (study, experiment)
            previous = records.get(key)
            if previous is None:
                records[key] = record
            elif (previous.group_0_name, previous.group_1_name, previous.condition) != (
                record.group_0_name,
                record.group_1_name,
                record.condition,
            ):
                raise ValueError(
                    f"experiment {key} has conflicting group metadata across signatures"
                )
    return records


def oriented_profiles(
    profiles: Sequence[Profile],
    classes: Mapping[tuple[str, str], Classification],
    *,
    strong_evidence_only: bool = False,
    control_tiers: Iterable[str] | None = None,
) -> list[Profile]:
    """Disease-versus-healthy profiles, with reversed experiments relabelled.

    After this, ``inc`` means higher in the disease group for every kept profile.
    The relabelling is explicit and only applied where the classification says
    Group 1 is the comparison group. When ``control_tiers`` is supplied, only
    experiments whose ``control_tier`` is in that set are kept.
    """
    allowed_tiers = set(control_tiers) if control_tiers is not None else None
    flip = {"inc": "dec", "dec": "inc"}
    kept: list[Profile] = []
    for profile in profiles:
        found = classes.get((profile.study_id, profile.experiment_id))
        if found is None or found.contrast_class != DISEASE_VS_HEALTHY:
            continue
        if strong_evidence_only and "weak_comparison_evidence" in found.flags.split(","):
            continue
        if allowed_tiers is not None and found.control_tier not in allowed_tiers:
            continue
        if found.orientation == REVERSED:
            kept.append(replace(profile, inc_dec=flip[profile.inc_dec]))
        else:
            kept.append(profile)
    return kept


#: The committed hand review. Tracked, because the agreement rate it supports is only
#: auditable if the cases behind it are, and because a later rule change has to be
#: measurable against the same experiments.
DEFAULT_REVIEW_PATH = REPO_ROOT / "methods" / "contrast_orientation_review.tsv"
REVIEW_FIELDS = (
    "study_id",
    "experiment_id",
    "pmid",
    "condition",
    "group_0_name",
    "group_1_name",
    "sampled_class",
    "sampled_orientation",
    "reviewer",
    "verdict",
    # The note is free text and is empty for most agreements, so it is not the last
    # column: a trailing empty field is trailing whitespace, which the pre-commit hook
    # strips, which would silently drop a column from every agreeing row.
    "note",
    "reviewed_class",
    "reviewed_orientation",
)
NOT_REVIEWED = "not_reviewed"


def review_sample(
    records: Mapping[tuple[str, str], ExperimentRecord],
    classes: Mapping[tuple[str, str], Classification],
    corpus: Iterable[tuple[str, str]],
    *,
    per_class: int,
    seed: int,
) -> list[dict[str, str]]:
    """Draw a reproducible sample of analysed experiments, stratified by contrast class.

    Each experiment's place in the draw is a digest of the seed, its class and its own
    identifier, so the same seed always gives the same experiments, a larger sample
    extends a smaller one instead of replacing it, and an experiment's place does not
    depend on what else is in the corpus. A class with fewer than ``per_class``
    analysed experiments contributes all of them.
    """
    by_class: dict[str, list[tuple[str, str]]] = {label: [] for label in CLASSES}
    for key in sorted(set(corpus)):
        by_class[classes[key].contrast_class].append(key)

    def place(label: str, key: tuple[str, str]) -> str:
        return hashlib.sha256(f"{seed}:{label}:{key[0]}/{key[1]}".encode()).hexdigest()

    drawn: list[tuple[str, str]] = []
    for label in CLASSES:
        keys = sorted(by_class[label], key=lambda key: place(label, key))
        drawn.extend(sorted(keys[:per_class]))
    rows = []
    for key in drawn:
        record, found = records[key], classes[key]
        rows.append(
            {
                "study_id": record.study_id,
                "experiment_id": record.experiment_id,
                "pmid": record.pmid,
                "condition": record.condition,
                "group_0_name": record.group_0_name,
                "group_1_name": record.group_1_name,
                "sampled_class": found.contrast_class,
                "sampled_orientation": found.orientation,
                "reviewer": "",
                "verdict": "",
                "note": "",
                "reviewed_class": "",
                "reviewed_orientation": "",
            }
        )
    return rows


def read_review(path: Path) -> dict[tuple[str, str], str]:
    """Each reviewed experiment's validation status, from the committed review table."""
    if not path.exists():
        return {}
    statuses: dict[tuple[str, str], str] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            verdict = row["verdict"].strip().lower()
            if verdict not in {"agree", "disagree"}:
                msg = (
                    f"{path}: experiment {row['experiment_id']} has verdict {row['verdict']!r}; "
                    "every sampled row must be reviewed as 'agree' or 'disagree'"
                )
                raise ValueError(msg)
            statuses[(row["study_id"], row["experiment_id"])] = f"reviewed_{verdict}d"
    return statuses


MANIFEST_FIELDS = (
    "study_id",
    "experiment_id",
    "pmid",
    "doi",
    "study_design",
    "condition",
    "efo_id",
    "body_site",
    "uberon_id",
    "host_species",
    "group_0_name",
    "group_1_name",
    "group_1_definition",
    "group_0_sample_size",
    "group_1_sample_size",
    "description",
    "direction_convention",
    "contrast_class",
    "orientation",
    "contrast_type",
    "control_tier",
    "evidence",
    "flags",
    "in_analysed_corpus",
    "validation_status",
    "classifier_version",
    "bugsigdb_vintage",
)


def write_manifest(
    path: Path,
    records: Mapping[tuple[str, str], ExperimentRecord],
    classes: Mapping[tuple[str, str], Classification],
    corpus: Iterable[tuple[str, str]],
    vintage: str | None,
    validation: Mapping[tuple[str, str], str] | None = None,
) -> Path:
    in_corpus = set(corpus)
    reviewed = validation or {}
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=MANIFEST_FIELDS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for key in sorted(records):
            found = classes[key]
            writer.writerow(
                {
                    **asdict(records[key]),
                    "direction_convention": DIRECTION_CONVENTION,
                    **asdict(found),
                    "in_analysed_corpus": "true" if key in in_corpus else "false",
                    "validation_status": reviewed.get(key, NOT_REVIEWED),
                    "classifier_version": CLASSIFIER_VERSION,
                    "bugsigdb_vintage": vintage or "",
                }
            )
    return path


def class_counts(
    classes: Mapping[tuple[str, str], Classification], corpus: Iterable[tuple[str, str]]
) -> list[dict[str, object]]:
    """Experiments and studies per class and orientation within the analysed corpus."""
    experiments: Counter[tuple[str, str]] = Counter()
    studies: dict[tuple[str, str], set[str]] = {}
    for key in corpus:
        found = classes[key]
        label = (found.contrast_class, found.orientation)
        experiments[label] += 1
        studies.setdefault(label, set()).add(key[0])
    return [
        {
            "contrast_class": label[0],
            "orientation": label[1],
            "experiments": experiments[label],
            "studies": len(studies[label]),
        }
        for label in sorted(experiments, key=lambda item: (CLASSES.index(item[0]), item[1]))
    ]


def control_tier_counts(
    classes: Mapping[tuple[str, str], Classification], corpus: Iterable[tuple[str, str]]
) -> list[dict[str, object]]:
    """Experiments and studies per contrast_type and control_tier within disease_vs_healthy."""
    experiments: Counter[tuple[str, str]] = Counter()
    studies: dict[tuple[str, str], set[str]] = {}
    for key in corpus:
        found = classes[key]
        if found.contrast_class != DISEASE_VS_HEALTHY:
            continue
        label = (found.contrast_type, found.control_tier)
        experiments[label] += 1
        studies.setdefault(label, set()).add(key[0])
    return [
        {
            "contrast_type": label[0],
            "control_tier": label[1],
            "experiments": experiments[label],
            "studies": len(studies[label]),
        }
        for label in sorted(
            experiments,
            key=lambda item: (CONTROL_TIERS.index(item[1]), item[0]),
        )
    ]


def comparison_rows(
    populations: Sequence[tuple[str, Sequence[Profile]]],
) -> list[dict[str, object]]:
    """Six fraction directions and 15 log-ratios on each population, side by side.

    Each population keeps the native families of the two analyses: Benjamini-Hochberg
    over six fractions and Benjamini-Yekutieli over 15 ratios. Nothing historical is
    overwritten.
    """
    rows: list[dict[str, object]] = []
    for population, profiles in populations:
        studies = len({p.study_id for p in profiles})
        experiments = len({(p.study_id, p.experiment_id) for p in profiles})
        shifts = compute_direction_shifts(list(profiles))
        for shift in shifts.shifts:
            rows.append(
                {
                    "population": population,
                    "comparison": "fraction",
                    "numerator": shift.trait,
                    "reference": "",
                    "population_studies": studies,
                    "population_experiments": experiments,
                    "studies_scored": shift.studies_scored,
                    "studies_up": shift.studies_up,
                    "studies_down": shift.studies_down,
                    "studies_tied": shift.studies_tied,
                    "share_up": shift.share_up,
                    "ci_low": shift.ci_low,
                    "ci_high": shift.ci_high,
                    "mean_effect": None,
                    "p_value": shift.p_value,
                    "q_value": shift.q_value,
                    "family": "Benjamini-Hochberg over 6 fractions",
                }
            )
        for summary in compute_log_ratio_shifts(list(profiles)).summaries:
            rows.append(
                {
                    "population": population,
                    "comparison": "log_ratio",
                    "numerator": summary.numerator,
                    "reference": summary.reference,
                    "population_studies": studies,
                    "population_experiments": experiments,
                    "studies_scored": summary.studies_scored,
                    "studies_up": summary.studies_up,
                    "studies_down": summary.studies_down,
                    "studies_tied": summary.studies_tied,
                    "share_up": summary.share_ratio_up,
                    "ci_low": summary.share_ci_low,
                    "ci_high": summary.share_ci_high,
                    "mean_effect": summary.mean_log2_ratio_change,
                    "p_value": summary.sign_p,
                    "q_value": summary.by_q,
                    "family": "Benjamini-Yekutieli over 15 ratios",
                }
            )
    return rows


def _format(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def _write_table(path: Path, rows: Sequence[Mapping[str, object]]) -> Path:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        if rows:
            writer.writerow(list(rows[0]))
            writer.writerows([_format(value) for value in row.values()] for row in rows)
    return path


_sha256 = file_sha256


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile-input", type=Path, default=DEFAULT_PROFILE_PATH)
    parser.add_argument("--experiments-input", type=Path, default=DEFAULT_EXPERIMENTS_PATH)
    parser.add_argument("--full-dump-input", type=Path, default=DEFAULT_FULL_DUMP_PATH)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--code-commit", help="Source revision, supplied by the Justfile recipe")
    parser.add_argument(
        "--review-input",
        type=Path,
        default=DEFAULT_REVIEW_PATH,
        help="Committed hand review, joined to the manifest as validation_status",
    )
    parser.add_argument(
        "--write-review-sample",
        type=Path,
        help="Draw a fresh review sample to this path instead of classifying the corpus",
    )
    parser.add_argument("--review-seed", type=int, default=188)
    parser.add_argument("--review-per-class", type=int, default=10)
    args = parser.parse_args(argv)

    vintage = read_bugsigdb_vintage(args.full_dump_input)
    records = read_experiment_records(args.full_dump_input)
    classes = {key: classify(record) for key, record in records.items()}
    profiles, stats = load_profiles(
        args.profile_input, args.experiments_input, TABLE1_TRAITS, reject_missing_counts=True
    )
    write_stats_report(stats, vintage)
    corpus = sorted({(p.study_id, p.experiment_id) for p in profiles})
    unmatched = [key for key in corpus if key not in classes]
    if unmatched:
        raise ValueError(f"{len(unmatched)} analysed experiments have no export record")

    if args.write_review_sample:
        rows = review_sample(
            records,
            classes,
            corpus,
            per_class=args.review_per_class,
            seed=args.review_seed,
        )
        args.write_review_sample.parent.mkdir(parents=True, exist_ok=True)
        with args.write_review_sample.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=REVIEW_FIELDS, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)
        print(f"{len(rows)} experiments sampled to {args.write_review_sample}")
        print(f"seed {args.review_seed}, {args.review_per_class} per class")
        return 0

    validation = read_review(args.review_input)
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    counts = class_counts(classes, corpus)
    tier_counts = control_tier_counts(classes, corpus)
    oriented = oriented_profiles(profiles, classes)
    written = [
        write_manifest(
            output / "contrast_manifest.tsv", records, classes, corpus, vintage, validation
        ),
        _write_table(output / "contrast_class_counts.tsv", counts),
        _write_table(
            output / "disease_control_comparison.tsv",
            comparison_rows(
                [
                    ("historical_human_paired", profiles),
                    ("disease_vs_healthy_oriented", oriented),
                    (
                        "disease_vs_healthy_strong_evidence",
                        oriented_profiles(profiles, classes, strong_evidence_only=True),
                    ),
                    (
                        "disease_vs_healthy_tier1_healthy_normal",
                        oriented_profiles(profiles, classes, control_tiers={TIER_1_HEALTHY_NORMAL}),
                    ),
                    (
                        "disease_vs_healthy_tier1_2_subject_controls",
                        oriented_profiles(
                            profiles,
                            classes,
                            control_tiers={
                                TIER_1_HEALTHY_NORMAL,
                                TIER_2_CONDITION_NEGATIVE_SCREENED,
                            },
                        ),
                    ),
                ]
            ),
        ),
    ]
    provenance = {
        "classifier_version": CLASSIFIER_VERSION,
        "direction_convention": DIRECTION_CONVENTION,
        "bugsigdb_vintage": vintage,
        "code_commit": args.code_commit,
        "full_dump_sha256": _sha256(args.full_dump_input),
        "profile_sha256": _sha256(args.profile_input),
        "experiments_sha256": _sha256(args.experiments_input),
        "analysis_source_sha256": _sha256(Path(__file__)),
        "analysed_experiments": len(corpus),
        "analysis_population": (
            "disease_vs_healthy experiments; reversed ones relabelled so inc = higher in disease"
        ),
        "class_counts": counts,
        "control_tier_counts": tier_counts,
        "review": {
            "path": str(args.review_input.relative_to(REPO_ROOT))
            if args.review_input.is_relative_to(REPO_ROOT)
            else str(args.review_input),
            "sha256": _sha256(args.review_input) if args.review_input.exists() else None,
            "reviewed_experiments": len(validation),
            "agreed": sum(status == "reviewed_agreed" for status in validation.values()),
            "disagreed": sum(status == "reviewed_disagreed" for status in validation.values()),
        },
        "limitations": [
            "lexical rules over group names; not a curated re-reading of each paper",
            "the export has no Group 0 definition",
            "validation status: stratified-sample review only; every other experiment is "
            "not_reviewed in the manifest",
        ],
    }
    provenance_path = output / "contrast_orientation_provenance.json"
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    written.append(provenance_path)
    for row in counts:
        print(row)
    for path in written:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
