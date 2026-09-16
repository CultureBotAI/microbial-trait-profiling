"""What each BugSigDB experiment compares, and which way its direction points (#188).

BugSigDB records a signature as "Abundance in Group 1" relative to Group 0, so
``inc`` means higher in Group 1. The Table 1 analyses keep that curated direction.
Reading it as "higher in disease" is only right when Group 1 is a disease group and
Group 0 a healthy or unaffected comparison, and nothing upstream checks that.

This module classifies every experiment from the export's own fields - study
design, condition, both group names, the Group 1 definition and the description -
into one of six contrast classes, records the matched evidence, and derives an
orientation only where the evidence supports one:

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
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import asdict, dataclass, replace
from pathlib import Path

from bsdbng.paths import DERIVED_DATA_DIR
from bsdbng.table1_trait_maps import compute_direction_shifts
from bsdbng.trait_composition import compute_log_ratio_shifts
from bsdbng.trait_guilds import (
    DEFAULT_EXPERIMENTS_PATH,
    DEFAULT_FULL_DUMP_PATH,
    DEFAULT_PROFILE_PATH,
    TABLE1_TRAITS,
    Profile,
    load_profiles,
    read_bugsigdb_vintage,
    write_stats_report,
)

CLASSIFIER_VERSION = "contrast-orientation-v1"
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


#: A group that is plainly a healthy or unaffected comparison. Only these can
#: support a reversed orientation.
STRONG_COMPARISON = _pattern(
    r"healthy",
    r"controls?",
    r"hcs?",
    r"unaffected",
    r"non[- ]?diseased",
    r"disease[- ]free",
    r"uninfected",
    r"non[- ]?infected",
    r"normal (?:controls?|subjects?|individuals?)",
    r"volunteers?",
)
#: A comparison group named by the absence of the condition's defining feature.
WEAK_COMPARISON = _pattern(
    r"normal[- ]weight",
    r"normal glucose tolerance",
    r"lean",
    r"normo\w+",
    r"negative",
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
#: Condition labels that are variables rather than diseases.
NON_DISEASE_CONDITION = re.compile(
    r"measurement|status|lifestyle|population|ethnic|^age$|^aging$|^sex$|participation"
    r"|protocol|location|geograph|body mass index|socioeconomic",
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
    value = value.strip()
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
    """Six-letter stems of the condition's words, and its initialisms."""
    stems: set[str] = set()
    initials: set[str] = set()
    for part in condition.split(","):
        words = [w for w in re.findall(r"[a-z]+", part.lower()) if w not in CONDITION_STOPWORDS]
        stems.update(w[:6] for w in words if len(w) >= 4)
        full = [
            w for w in re.findall(r"[a-z]+", part.lower()) if w not in {"and", "of", "the", "with"}
        ]
        if len(full) >= 2:
            initials.add("".join(w[0] for w in full))
    return stems, initials


def _mentions_condition(text: str, condition: str) -> bool:
    stems, initials = _condition_words(condition)
    tokens = re.findall(r"[a-z0-9]+", text.lower())
    return any(t[:6] in stems for t in tokens if len(t) >= 4) or any(
        t.startswith(i) for t in tokens for i in initials if len(i) >= 2
    )


def _absence_of_condition(name: str, condition: str) -> list[str]:
    """ "non-IBS", "no asthma", "without Parkinson's" - absence of this condition."""
    found = []
    folded = re.sub(r"[\W_]+", " ", condition.lower())
    for match in re.finditer(
        "(?<![\\w])free (?:from|of)\\s+((?:[\\w'\u2019.]+\\s*){1,3})", name, re.IGNORECASE
    ):
        if _mentions_condition(match.group(1), condition):
            found.append(match.group(0).strip().lower())
    for match in re.finditer(
        "(?<![\\w])(non|no|without)[- ]+((?:[\\w'\u2019.]+\\s*){1,3})", name, re.IGNORECASE
    ):
        phrase = re.sub(r"[\W_]+", " ", match.group(0).lower()).strip()
        # "Non-alcoholic fatty liver disease" names the condition rather than its absence.
        if " ".join(phrase.split()[:2]) in folded:
            continue
        if _mentions_condition(match.group(2), condition):
            found.append(match.group(0).strip().lower())
    return found


def _comparison(name: str, other: str, condition: str) -> tuple[list[str], list[str]]:
    """Strong and weak evidence that ``name`` is the comparison group."""
    strong = _matches(STRONG_COMPARISON, name)
    weak = _matches(WEAK_COMPARISON, name) + _absence_of_condition(name, condition)
    mine, theirs = _markers(name), _markers(other)
    negative = [m for m, sign in mine.items() if sign == "-" and theirs.get(m) != "-"]
    if negative:
        weak.append(f"test-negative:{','.join(sorted(negative))}")
    return strong, weak


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

    strong_0, weak_0 = _comparison(g0, g1, condition)
    strong_1, weak_1 = _comparison(g1, g0, condition)
    control_0, control_1 = strong_0 + weak_0, strong_1 + weak_1

    resolved = _matches(RESOLVED, f"{g0} | {g1}")
    if resolved:
        return result(LONGITUDINAL, NOT_APPLICABLE, f"a group is resolved disease: {resolved}")

    time_0, time_1 = set(_matches(TIMEPOINT, g0)), set(_matches(TIMEPOINT, g1))
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
    if TREATMENT_CONDITION.search(condition):
        return result(
            TREATMENT, NOT_APPLICABLE, f"condition is a treatment or exposure: {condition!r}"
        )
    if one_sided_time:
        return result(
            LONGITUDINAL,
            NOT_APPLICABLE,
            f"one group is named by time: {sorted(time_0)} / {sorted(time_1)}",
        )
    if NON_DISEASE_CONDITION.search(condition) and not (control_0 or control_1):
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
        if NON_DISEASE_CONDITION.search(condition):
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

    if NON_DISEASE_CONDITION.search(condition):
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
    "classifier_version",
    "bugsigdb_vintage",
)


def write_manifest(
    path: Path,
    records: Mapping[tuple[str, str], ExperimentRecord],
    classes: Mapping[tuple[str, str], Classification],
    corpus: Iterable[tuple[str, str]],
    vintage: str | None,
) -> Path:
    in_corpus = set(corpus)
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

    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    counts = class_counts(classes, corpus)
    oriented = oriented_profiles(profiles, classes)
    written = [
        write_manifest(output / "contrast_manifest.tsv", records, classes, corpus, vintage),
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
        "limitations": [
            "lexical rules over group names; not a curated re-reading of each paper",
            "the export has no Group 0 definition",
            "validation status: representative-sample review only (see methods)",
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
