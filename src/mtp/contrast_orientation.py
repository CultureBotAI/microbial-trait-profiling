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

from mtp.paths import DERIVED_DATA_DIR, REPO_ROOT
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

CLASSIFIER_VERSION = "contrast-orientation-v2"
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


def _pattern(*terms: str) -> re.Pattern[str]:
    return re.compile(r"(?<![\w])(?:" + "|".join(terms) + r")(?![\w])", re.IGNORECASE)


#: A group that is plainly well, whatever the condition is. A multi-valued condition
#: needs one of these, because denying one of its values leaves the others standing.
GENERAL_COMPARISON = _pattern(
    r"healthy",
    r"controls?",
    r"hcs?",
    r"unaffected",
    r"non[- ]?diseased",
    r"disease[- ]free",
    r"normal (?:controls?|subjects?|individuals?|groups?|samples?)",
    r"volunteers?",
)
#: A group that is plainly a healthy or unaffected comparison. Only these can
#: support a reversed orientation. The last two deny one named thing rather than
#: asserting general health.
STRONG_COMPARISON = _pattern(
    r"healthy",
    r"controls?",
    r"hcs?",
    r"unaffected",
    r"non[- ]?diseased",
    r"disease[- ]free",
    r"uninfected",
    r"non[- ]?infected",
    r"normal (?:controls?|subjects?|individuals?|groups?|samples?)",
    r"volunteers?",
)
#: A comparison group named by the absence of the condition's defining feature.
#: A bare "negative" is not here: it is a comparison only where its subject is the
#: condition, which :func:`_condition_negative` decides.
WEAK_COMPARISON = _pattern(
    r"normal[- ]weight",
    r"normal glucose tolerance",
    r"lean",
    r"normo\w+",
)
#: "<subject> negative", with everything before the subject kept so that a negation
#: qualifying one disease inside another can be told apart from a negation of the
#: condition: "HIV negative" and "HIV-1-negative" deny HIV infection, while "Her2
#: negative", "ASCA negative", "cART negative" and "Parkinson's disease dementia
#: negative" name a marker, a drug or a second disease in someone who is ill.
CONDITION_NEGATIVE = re.compile(
    "(?P<before>.*?)(?P<subject>[\\w'\\u2019.]+)[^\\w]+negative(?![\\w])", re.IGNORECASE
)
#: "Negative for COVID-19 viral RNA": the sign comes before its subject.
NEGATIVE_FOR = re.compile(
    "(?<![\\w])negative\\s+for\\s+((?:[\\w'\\u2019.]+[\\s-]*){1,4})", re.IGNORECASE
)
#: "HC-like IBS", "Healthy Control like Irritable Bowel Syndrome". A group that
#: resembles the comparison group is not the comparison group: this is an
#: irritable-bowel subgroup whose profile happens to look like a healthy one.
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
#: The same subjects sampled at different times.
TIMEPOINT = _pattern(
    r"before",
    r"after",
    r"pre",
    r"post",
    r"pre-\w+",
    r"post-\w+",
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
#: An age band, which reads like a duration but dates nobody: "2-12 years" and
#: "aged > 10 years" are how old the subjects were, not when they were sampled.
AGE_BAND = re.compile(
    "\\d+\\s*[-\\u2013]\\s*\\d+\\s*(?:years?|yrs?|months?|mo)(?![\\w])"
    "|aged?\\s*[<>=\\u2264\\u2265\\s]*\\d+\\s*(?:years?|yrs?|months?|mo)?(?![\\w])",
    re.IGNORECASE,
)
#: Time words that, beside a healthy comparison, only date the disease sample.
BASELINE_ONLY = {"baseline", "admission"}
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
)
#: A disease group named by its disease, stage, subtype, severity or activity.
DISEASE_TERMS = _pattern(
    r"patients?",
    r"cases?",
    r"diseases?",
    r"disorders?",
    r"syndromes?",
    r"cancers?",
    r"\w+omas?",
    r"\w+itis",
    r"tumou?rs?",
    r"polyps?",
    r"lesions?",
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
    r"hypertensi\w*",
    r"positive",
)
#: Condition labels that are variables rather than diseases. Tested against each
#: comma-separated value of a multi-valued condition, so "Age" is recognised in
#: "Coinfection,Age"; the condition as a whole is a variable only when every value is.
NON_DISEASE_CONDITION = re.compile(
    r"measurement|status|lifestyle|population|ethnic|^age$|^aging$|^sex$|participation"
    r"|protocol|location|geograph|body mass index|socioeconomic|^pregnancy$",
    re.IGNORECASE,
)
#: Condition labels that are a treatment or exposure. "Non-alcoholic" is a disease
#: name, not alcohol exposure, hence the guard against a preceding hyphen.
TREATMENT_CONDITION = re.compile(
    r"(?<![\w-])(?:response to|antimicrobial agent|antibiotic\w*|smoking|alcohol\w*|diet"
    r"|drug\w*|vaccin\w*|surgery|transplant\w*|c(?:a)?esarean\w*|breast ?feeding|delivery"
    r")(?![\w])",
    re.IGNORECASE,
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
#: A treatment word that is negated names untreated disease, not a treatment.
NEGATED_TREATMENT = re.compile(
    r"(?:without|no|not|never)\s+(?:[\w-]+\s+){0,2}(?:treat\w*|therap\w*|antibiotic\w*"
    r"|antiretroviral\w*|medication\w*)|(?:treatment|therapy|drug|antibiotic)[- ]naive|untreated",
    re.IGNORECASE,
)
#: Both groups named as low/high levels of one variable.
LEVEL_SPLIT = _pattern(
    r"low\w*", r"high\w*", r"levels?", r"above", r"below", r"quartiles?", r"tertiles?"
)
MULTI_DESIGN_TREATMENT = re.compile(r"randomized controlled trial", re.IGNORECASE)
MULTI_DESIGN_TIME = re.compile(r"time series|longitudinal", re.IGNORECASE)
#: A sign counts as a test result only where the name or a sub-name ends after it.
MARKER = re.compile("([a-z][a-z0-9.]*)\\s*([+\\-\u2212])(?=\\s*(?:$|[/),;(]))", re.IGNORECASE)
#: Words naming the same disease category, so that a negation can deny a condition
#: stated in other words: "noncancer individuals" against oral squamous cell carcinoma.
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
    """Test markers and their sign, as in "HPV-" or "HIV +"."""
    return {marker.lower(): "-" if sign != "+" else "+" for marker, sign in MARKER.findall(name)}


def _condition_words(condition: str) -> tuple[set[str], set[str]]:
    """Six-letter stems of the condition's words, and its initialisms.

    A condition already written as an acronym ("HIV infection") contributes that acronym
    as a stem, since it is too short for the six-letter rule and its initialism ("hi")
    describes nothing.
    """
    stems: set[str] = set()
    initials: set[str] = set()
    for part in condition.split(","):
        words = [w for w in re.findall(r"[a-z]+", part.lower()) if w not in CONDITION_STOPWORDS]
        stems.update(w[:6] for w in words if len(w) >= 4)
        stems.update(w for w in words if 2 <= len(w) <= 3)
        full = [
            w for w in re.findall(r"[a-z]+", part.lower()) if w not in {"and", "of", "the", "with"}
        ]
        if len(full) >= 2:
            initials.add("".join(w[0] for w in full))
    return stems, initials


def _mentions_condition(text: str, condition: str) -> bool:
    stems, initials = _condition_words(condition)
    tokens = re.findall(r"[a-z0-9]+", text.lower())
    return (
        any(t[:6] in stems for t in tokens if len(t) >= 4)
        or any(t in stems for t in tokens if len(t) <= 3)
        # Either direction: the condition's initialism can be the longer string, as
        # "hpv" against the initialism of "Human papilloma virus infection".
        or any(
            t.startswith(i) or (len(t) >= 3 and i.startswith(t))
            for t in tokens
            for i in initials
            if len(i) >= 2
        )
    )


def _same_disease_category(word: str, condition: str) -> bool:
    """Do a word and the condition name the same disease category?

    BugSigDB states a condition as precisely as the paper does, so a group named by the
    absence of the broad category ("noncancer") denies a condition named as a specific
    instance of it ("Squamous cell carcinoma") even though the two share no letters.
    """
    words = set(re.findall(r"[a-z]+", condition.lower()))
    return any(word.lower() in group and words & group for group in CONDITION_SYNONYMS)


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
    """Time words in a group name, excluding those inside an age band.

    "2-12 years" against "18-62 years" and "Controls aged > 10 years" look like
    durations to :data:`TIMEPOINT` but say how old the subjects were, so the same
    subjects are not being followed through time.
    """
    spans = [match.span() for match in AGE_BAND.finditer(name)]
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
    absence = _absence_of_condition(name, condition) + _condition_negative(name, condition)
    weak = _matches(WEAK_COMPARISON, resemblance) + [text for text, _span in absence]
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


def classify(record: ExperimentRecord) -> Classification:
    """Assign one contrast class and an orientation, with the evidence for both."""
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

    def result(contrast_class: str, orientation: str, evidence: str) -> Classification:
        return Classification(contrast_class, orientation, evidence, ",".join(flags))

    if not g0 or not g1:
        return result(UNRESOLVED, NOT_APPLICABLE, "a group name is missing")

    strong_0, weak_0, withdrawn_0 = _comparison(g0, g1, condition)
    strong_1, weak_1, withdrawn_1 = _comparison(g1, g0, condition)
    control_0, control_1 = strong_0 + weak_0, strong_1 + weak_1
    if withdrawn_0 or withdrawn_1:
        flags.append("comparison_reading_withdrawn")
    variable_condition = _variable_condition(condition)

    resolved = _matches(RESOLVED, f"{g0} | {g1}")
    if resolved:
        return result(LONGITUDINAL, NOT_APPLICABLE, f"a group is resolved disease: {resolved}")

    time_0, time_1 = _timepoints(g0), _timepoints(g1)
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

    untreated = NEGATED_TREATMENT.sub(" ", f"{g0} | {g1}")
    if untreated != f"{g0} | {g1}":
        flags.append("negated_treatment_term")
    treatment_terms = _matches(TREATMENT_TERMS, untreated)
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
            return result(DISEASE_VS_HEALTHY, NATIVE, f"Group 0 reads as a comparison: {control_0}")
        disease_0 = _disease_evidence(g0, condition)
        if strong_1 and disease_0:
            return result(
                DISEASE_VS_HEALTHY,
                REVERSED,
                f"Group 1 reads as a comparison: {strong_1}; Group 0 as disease: {disease_0}",
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
) -> list[Profile]:
    """Disease-versus-healthy profiles, with reversed experiments relabelled.

    After this, ``inc`` means higher in the disease group for every kept profile.
    The relabelling is explicit and only applied where the classification says
    Group 1 is the comparison group.
    """
    flip = {"inc": "dec", "dec": "inc"}
    kept: list[Profile] = []
    for profile in profiles:
        found = classes.get((profile.study_id, profile.experiment_id))
        if found is None or found.contrast_class != DISEASE_VS_HEALTHY:
            continue
        if strong_evidence_only and "weak_comparison_evidence" in found.flags.split(","):
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


def comparison_rows(
    populations: Sequence[tuple[str, Sequence[Profile]]],
) -> list[dict[str, object]]:
    """Six fraction directions and 15 log-ratios on both populations, side by side.

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


def _sha256(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


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
