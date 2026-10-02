"""What measured evidence could independently test the Table 1 trait signals (#191).

The profile matrix holds annotation counts inferred from taxa in curated signatures.
A redox or short-chain fatty acid (SCFA) interpretation of those counts needs
measurements the matrix does not contain: inflammation markers, nitrate or oxygen,
SCFA concentrations, antibiotic exposure, diet, transit. This module takes stock of
what is recorded, without acquiring new data:

* **Per study**, which kinds of measured evidence the export's curated text mentions
  (condition, group names and definition, description, matching and confounder
  fields, title, keywords), and which sequencing modality the study used. A mention
  is a pointer to the source paper, not data: the export carries no sample-level
  values, and "not mentioned" never means "not measured".
* **Marker-split candidates**: every experiment whose condition or group fields name
  a marker of one of those evidence types, each classified as a split on the
  marker's measured value, an administered intervention, or a contrast that merely
  carries the marker's name. Only the first kind links a signature to a measurement,
  and even then it is part of the same BugSigDB corpus and not independent of it.
* **Inferred-function traits**: KG-Microbe trait columns whose labels concern SCFA
  production or fermentation. They are annotations, not measurements, and share
  their taxa with every other trait.
* **An external cohort collection** (Sun et al. 2024, 36 Chinese fecal metagenome
  studies): candidate overlaps with the discovery corpus by author surname, year and
  location, for confirmation before any "independent" claim.
"""

from __future__ import annotations

import argparse
import csv
import difflib
import hashlib
import json
import re
from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path

from bsdbng.paths import DERIVED_DATA_DIR, KG_MICROBE_DIR, REPO_ROOT
from bsdbng.trait_guilds import (
    DEFAULT_EXPERIMENTS_PATH,
    DEFAULT_FULL_DUMP_PATH,
    DEFAULT_PROFILE_PATH,
    MEDIADIVE_MARKER,
    PROFILE_ID_COLUMNS,
    TABLE1_TRAITS,
    load_profiles,
    read_bugsigdb_vintage,
)

DEFAULT_OUTPUT_DIR = DERIVED_DATA_DIR / "validation_evidence"
DEFAULT_NODES_PATH = KG_MICROBE_DIR / "merged-kg_nodes.tsv"
EXTERNAL_COHORTS_PATH = REPO_ROOT / "methods" / "external_cohorts_sun2024.tsv"
#: The external sample table, pinned to a commit so the snapshot is reproducible.
SUN2024_SAMPLE_GROUP_URL = (
    "https://raw.githubusercontent.com/yexianingyue/GM_common_diseases/"
    "bd254f4e6f6857c1d94eb6faf7941a0cdf489764/data/sample.group"
)
SUN2024_SAMPLE_GROUP_SHA256 = "0b3a8fef3c95731b53e6505f3f6e8c954cfd86d2b37ff3ca59d3bb0cead373e2"
#: Supplementary Table 1 of the same paper, as text, for the source study titles.
SUN2024_SUPPLEMENT_URL = (
    "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11442664/supplementaryFiles"
)
#: Two titles name the same publication at or above this character-level similarity
#: once case, spacing and punctuation are removed, which absorbs PDF artefacts such
#: as "Ametagenome" for "A metagenome" without matching merely similar titles.
TITLE_AGREEMENT = 0.9
VERSION = "validation-evidence-v1"

#: The export fields whose curated text is searched, in reporting order.
TEXT_FIELDS = (
    "Condition",
    "Group 0 name",
    "Group 1 name",
    "Group 1 definition",
    "Description",
    "Matched on",
    "Confounders controlled for",
    "Title",
    "Keywords",
)


def _terms(*terms: str) -> re.Pattern[str]:
    return re.compile(r"(?<![\w])(?:" + "|".join(terms) + r")(?![\w])", re.IGNORECASE)


#: Kinds of measured evidence, each with the words that point to it.
EVIDENCE_TYPES: dict[str, re.Pattern[str]] = {
    "inflammation_marker": _terms(
        r"calprotectin",
        r"lactoferrin",
        r"crp",
        r"c-reactive protein",
        r"hs-?crp",
        r"il-?\d+",
        r"interleukin[- ]?\d*",
        "tnf(?:-?\u03b1|-?alpha)?",
        r"esr",
        r"erythrocyte sedimentation rate",
        r"inflammatory markers?",
        r"cytokines?",
    ),
    "disease_activity_or_severity": _terms(
        r"mayo score",
        r"cdai",
        r"harvey[- ]bradshaw",
        r"updrs",
        r"hoehn(?: and |-)yahr",
        r"sledai",
        r"das28",
        r"pasi",
        r"scorad",
        r"nihss",
        r"meld",
        r"child[- ]pugh",
        r"disease activity",
        r"severity",
    ),
    "scfa_metabolite": _terms(
        r"scfas?",
        r"short[- ]chain fatty acids?",
        r"butyrate",
        r"butyric acid",
        r"propionate",
        r"propionic acid",
        r"acetate",
        r"acetic acid",
        r"valerate",
        r"metabolom\w*",
    ),
    "nitrate_oxygen_redox": _terms(
        r"nitrate",
        r"nitrite",
        r"nitric oxide",
        r"oxygen",
        r"hypoxi\w*",
        r"oxidative stress",
        r"redox",
    ),
    "antibiotic_exposure": _terms(
        r"antibiotics? (?:use|users|exposure|treatment|therapy|intake)",
        r"antibiotic-treated",
        r"received antibiotics?",
        r"prior antibiotics?",
    ),
    "diet_fiber_salt": _terms(
        r"diet\w*",
        r"fib(?:er|re)",
        r"salt",
        r"sodium",
        r"food frequency",
        r"nutri\w*",
    ),
    "transit": _terms(r"transit", r"bristol(?: stool)?", r"stool consistency", r"stool frequency"),
    "clinical_covariates": _terms(
        r"age",
        r"sex",
        r"gender",
        r"bmi",
        r"body mass index",
        r"smoking",
        r"medication\w*",
    ),
}

#: Markers whose measurement bears on the redox and SCFA reading of Table 1, by the
#: evidence type each belongs to. Disease-activity scores and clinical covariates are
#: deliberately absent: the export also holds experiments split on measured grey-matter
#: volume or triglycerides, which are measured splits but not of a marker in question.
#: The inflammation panel includes the growth factors assayed on the same multiplex as
#: the cytokines (bsdb:39103909), and names inflammation itself but not the disease
#: names that contain the word, so "inflammatory bowel disease" is not a marker.
MARKER_TERMS: dict[str, re.Pattern[str]] = {
    "inflammation_marker": _terms(
        r"calprotectin",
        r"lactoferrin",
        r"c-reactive protein",
        r"hs-?crp",
        r"crp",
        r"interleukin[- ]?\d*",
        r"il-?\d+\w*",
        r"cytokines?",
        r"chemokines?",
        r"tumou?r necrosis factor",
        "tnf(?:-?\u03b1|-?alpha)?",
        r"inflammation",
        r"inflammator\w* (?:response|marker|protein|cytokine|factor|mediator)s?",
        r"colony[- ]stimulating factor",
        r"g-?csf",
        r"monocyte chemo(?:tactic|attractant) protein(?: 1)?",
        r"mcp-?1",
        r"macrophage inflammatory protein(?: 1 ?[ab](?:lpha|eta)?)?",
        r"mip-?1[ab]?",
        r"interferon gamma[- ]induced protein(?: 10)?",
        r"ip-?10",
        r"vascular endothelial growth factor",
        r"vegf",
        r"platelet[- ]derived growth factor(?: bb)?",
        r"pdgf(?:-?bb)?",
        r"fibroblast growth factor(?: basic)?",
        r"fgf",
        r"erythrocyte sedimentation rate",
        r"esr",
    ),
    "scfa_metabolite": _terms(
        r"scfas?",
        r"short[- ]chain fatty acids?",
        r"butyrate",
        r"butyric acid",
        r"propionate",
        r"propionic acid",
        r"acetate",
        r"acetic acid",
        r"valerate",
    ),
    "nitrate_oxygen_redox": _terms(
        r"nitrate",
        r"nitrite",
        r"nitric oxide",
        r"oxygen",
        r"normoxi\w*",
        r"hypoxi\w*",
        r"oxidative stress",
        r"redox",
        r"malondialdehyde",
        r"mda",
        r"superoxide dismutase",
        r"sod",
    ),
}

#: A group *label* naming a measured value: a direction or level word, or a threshold
#: with units. Group names are labels for the two arms, so a direction word in one is
#: about the cut. This pattern is deliberately not applied to the Group 1 definition,
#: which is prose: the definition of bsdb:39787116-2 says AECOPD involves "increased
#: breathlessness, cough, mucus production, and inflammation", and that describes the
#: disease, not a cut on a measured inflammation value.
MEASURED_LEVEL = re.compile(
    r"(?<![\w])(?:increased|decreased|elevated|reduced|higher|lower|high|low|levels?"
    r"|concentrations?|quantified|measured|titers?|titres?|above|below)(?![\w])"
    r"|[<>]\s*\d|\d\s*(?:mg|ug|µg|ng|pg|mmol|µmol|nmol)(?![\w])",
    re.IGNORECASE,
)
#: Wording that states outright that the groups were cut on something, which is what
#: the Group 1 definition has to say before its prose is read as describing a cut.
SPLIT_WORDING = re.compile(
    r"grouped\s+(?:according\s+to|by)|defined\s+by|stratified\s+by|based\s+on"
    r"|divided\s+(?:by|according\s+to)|cut-?off|threshold",
    re.IGNORECASE,
)
#: A group label naming the marker's absence, the other half of a presence/absence cut.
MARKER_ABSENT = re.compile(r"^\s*(?:no|without|absence of|non-?)(?![\w])", re.IGNORECASE)
#: An administered intervention. Searched over the two group names only - the labels
#: that say what the two arms are. Not over the condition, which names the clinical
#: setting rather than the contrast: bsdb:31065565-2 is conditioned on "Response to
#: allogeneic hematopoietic stem cell transplant" but its groups are low against high
#: measured fecal butyrate. Not over the Group 1 definition either, whose prose often
#: names a therapy only to record an eligibility rule. "Dose" is absent on purpose:
#: bsdb:37789856 pools the same dose arms on both sides of a measured
#: superoxide-dismutase cut, so the dose is not what that contrast is.
INTERVENTION = re.compile(
    r"(?<![\w])(?:placebo|pro-?biotics?|pre-?biotics?|synbiotics?|supplement\w*"
    r"|therapy|therapeutic|therapies|treated|treatment|pills?|tablets?|capsules?"
    r"|vaccin\w*|transplant\w*|intervention|administration|infusion"
    r"|antibiotics?|antimicrobials?)(?![\w])",
    re.IGNORECASE,
)
#: Inferred-function trait labels about SCFA production or fermentation.
FUNCTION_LABEL = re.compile(
    r"(?<![\w-])(?:acetate|propionate|butyrate|butanoate|lactate|succinate)(?![a-z])"
    r"|short[- ]chain fatty|fermentation",
    re.IGNORECASE,
)


@dataclass(frozen=True)
class StudyRecord:
    study_id: str
    pmid: str
    doi: str
    title: str
    authors: str
    year: str
    location: str
    host_species: str
    sequencing_types: str
    antibiotics_exclusion_field: str
    texts: dict[str, list[str]]


def _clean(value: str | None) -> str:
    value = (value or "").strip()
    return "" if value.upper() == "NA" else value


def read_studies(full_dump_path: Path) -> dict[str, StudyRecord]:
    """One record per BugSigDB study, pooling curated text across its rows."""
    texts: defaultdict[str, defaultdict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    first: dict[str, dict[str, str]] = {}
    sequencing: defaultdict[str, set[str]] = defaultdict(set)
    exclusion: defaultdict[str, set[str]] = defaultdict(set)
    hosts: defaultdict[str, set[str]] = defaultdict(set)
    with full_dump_path.open(encoding="utf-8", newline="") as handle:
        line = handle.readline()
        if not line.startswith("#"):
            handle.seek(0)
        for row in csv.DictReader(handle):
            study = "bsdb:" + row["BSDB ID"].removeprefix("bsdb:").split("/")[0]
            first.setdefault(study, row)
            for name in TEXT_FIELDS:
                value = _clean(row.get(name))
                if value:
                    texts[study][name].add(value)
            if _clean(row.get("Sequencing type")):
                sequencing[study].add(_clean(row.get("Sequencing type")))
            if _clean(row.get("Antibiotics exclusion")):
                exclusion[study].add(_clean(row.get("Antibiotics exclusion")))
            if _clean(row.get("Host species")):
                hosts[study].add(_clean(row.get("Host species")))
    return {
        study: StudyRecord(
            study_id=study,
            pmid=_clean(row["PMID"]),
            doi=_clean(row["DOI"]),
            title=_clean(row["Title"]),
            authors=_clean(row["Authors list"]),
            year=_clean(row["Year"]),
            location=_clean(row["Location of subjects"]),
            host_species=";".join(sorted(hosts[study])),
            sequencing_types=";".join(sorted(sequencing[study])),
            antibiotics_exclusion_field=";".join(sorted(exclusion[study])),
            texts={name: sorted(values) for name, values in texts[study].items()},
        )
        for study, row in first.items()
    }


@dataclass(frozen=True)
class Mention:
    status: str
    fields: str
    terms: str


def find_mentions(record: StudyRecord, pattern: re.Pattern[str]) -> Mention:
    fields: list[str] = []
    terms: set[str] = set()
    for name in TEXT_FIELDS:
        found = {
            m.group(0).lower()
            for text in record.texts.get(name, [])
            for m in pattern.finditer(text)
        }
        if found:
            fields.append(name)
            terms |= found
    if fields:
        return Mention("mentioned_in_curated_metadata", ";".join(fields), ";".join(sorted(terms)))
    return Mention("not_mentioned_in_curated_metadata", "", "")


def evidence_rows(
    studies: Mapping[str, StudyRecord], discovery: set[str]
) -> list[dict[str, object]]:
    """The per-study evidence-availability matrix for human studies."""
    rows: list[dict[str, object]] = []
    for study_id, record in sorted(studies.items()):
        # A study whose experiments span hosts is kept; the host set is reported.
        if "Homo sapiens" not in record.host_species:
            continue
        row: dict[str, object] = {
            "study_id": study_id,
            "pmid": record.pmid,
            "doi": record.doi,
            "year": record.year,
            "location": record.location,
            "in_discovery_corpus": study_id in discovery,
            "sequencing_types": record.sequencing_types,
            "metagenomic_reads_possible": "WMS" in record.sequencing_types,
            "antibiotics_exclusion_field": record.antibiotics_exclusion_field,
            "observation_level": "study-level curated text; no sample-level values in the export",
            "access": "source article and supplements; not assessed",
        }
        for name, pattern in EVIDENCE_TYPES.items():
            mention = find_mentions(record, pattern)
            row[f"{name}_status"] = mention.status
            row[f"{name}_fields"] = mention.fields
            row[f"{name}_terms"] = mention.terms
        rows.append(row)
    return rows


def summary_rows(rows: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    """Studies mentioning each evidence type, overall and within the discovery corpus."""
    summary: list[dict[str, object]] = []
    discovery = [row for row in rows if row["in_discovery_corpus"]]
    for name in [*EVIDENCE_TYPES, "metagenomic_reads_possible"]:

        def hit(row: Mapping[str, object], name: str = name) -> bool:
            if name == "metagenomic_reads_possible":
                return bool(row[name])
            return row[f"{name}_status"] == "mentioned_in_curated_metadata"

        summary.append(
            {
                "evidence": name,
                "human_studies": len(rows),
                "human_studies_mentioning": sum(hit(row) for row in rows),
                "discovery_studies": len(discovery),
                "discovery_studies_mentioning": sum(hit(row) for row in discovery),
            }
        )
    return summary


#: The groups were cut on the marker's measured value.
MARKER_STRATIFIED = "marker_stratified"
#: The two groups differ by an administered intervention, so the cut is a trial arm.
INTERVENTION_CONTRAST = "intervention_contrast"
#: A marker is named, but nothing says the groups were cut on its measured value.
MARKER_NAMED_NOT_MEASURED = "marker_named_groups_not_measured"
#: The group fields whose text decides what the two groups are.
GROUP_FIELDS = ("Group 0 name", "Group 1 name", "Group 1 definition")


@dataclass(frozen=True)
class MarkerSplit:
    classification: str
    evidence_type: str
    terms: str
    split_evidence_field: str
    split_evidence_term: str
    exclusion_reason: str


def classify_marker_split(
    condition: str, group_0_name: str, group_1_name: str, group_1_definition: str
) -> MarkerSplit | None:
    """Decide whether an experiment's groups were cut on a measured marker.

    The condition label decides nothing on its own, in either direction. It can
    name a measurement whose groups are not one ("C-reactive protein measurement"
    labels bsdb:37180439-9, healthy controls against liver-disease patients), and
    it can fail to name the marker of a genuine measured cut ("Increased
    inflammatory response" labels bsdb:35888992-6, whose groups are inflammatory
    bowel disease patients either side of fecal calprotectin 1000 mg/kg). So the
    marker, and the evidence that the cut was on its measured value, are read from
    the group names and the Group 1 definition as well, and an intervention named in
    a group label disqualifies the contrast however the condition is labelled.

    Returns ``None`` when no marker of an evidence type in :data:`MARKER_TERMS` is
    named anywhere, which is the common case.
    """
    group_texts = (group_0_name, group_1_name, group_1_definition)
    everywhere = " | ".join(text for text in (condition, *group_texts) if text)
    found = {
        name: {match.group(0).lower() for match in pattern.finditer(everywhere)}
        for name, pattern in MARKER_TERMS.items()
        if pattern.search(everywhere)
    }
    if not found:
        return None
    evidence_type = ";".join(sorted(found))
    terms = ";".join(sorted({term for terms in found.values() for term in terms}))
    arm = INTERVENTION.search(" | ".join(text for text in (group_0_name, group_1_name) if text))
    if arm:
        return MarkerSplit(
            INTERVENTION_CONTRAST,
            evidence_type,
            terms,
            "",
            "",
            f"a group label names an administered intervention ({arm.group(0).lower()})",
        )
    for field, text in zip(GROUP_FIELDS, group_texts, strict=True):
        if not text or not any(pattern.search(text) for pattern in MARKER_TERMS.values()):
            continue
        is_label = field != "Group 1 definition"
        stated = SPLIT_WORDING.search(text) or (is_label and MEASURED_LEVEL.search(text))
        if stated:
            return MarkerSplit(
                MARKER_STRATIFIED, evidence_type, terms, field, stated.group(0).lower(), ""
            )
        if is_label and MARKER_ABSENT.search(text):
            return MarkerSplit(
                MARKER_STRATIFIED, evidence_type, terms, field, "marker named as absent", ""
            )
    return MarkerSplit(
        MARKER_NAMED_NOT_MEASURED,
        evidence_type,
        terms,
        "",
        "",
        "the group fields carry no measured level, threshold or stratification wording",
    )


def marker_split_candidates(
    full_dump_path: Path, discovery_experiments: set[tuple[str, str]]
) -> list[dict[str, object]]:
    """Every experiment naming a marker, with how its groups were actually cut.

    Candidates that are not measured splits are kept in the table with the reason,
    so the inventory can be audited rather than only counted.
    """
    seen: dict[tuple[str, str], dict[str, object]] = {}
    with full_dump_path.open(encoding="utf-8", newline="") as handle:
        if not handle.readline().startswith("#"):
            handle.seek(0)
        for row in csv.DictReader(handle):
            condition = _clean(row["Condition"])
            groups = tuple(_clean(row[name]) for name in GROUP_FIELDS)
            split = classify_marker_split(condition, *groups)
            if split is None:
                continue
            parts = row["BSDB ID"].removeprefix("bsdb:").split("/")
            key = (f"bsdb:{parts[0]}", f"bsdb:{parts[0]}-{parts[1]}")
            seen.setdefault(
                key,
                {
                    "study_id": key[0],
                    "experiment_id": key[1],
                    "classification": split.classification,
                    "marker_evidence_type": split.evidence_type,
                    "marker_terms": split.terms,
                    "split_evidence_field": split.split_evidence_field,
                    "split_evidence_term": split.split_evidence_term,
                    "exclusion_reason": split.exclusion_reason,
                    "condition": condition,
                    "host_species": _clean(row["Host species"]),
                    "body_site": _clean(row["Body site"]),
                    "group_0_name": groups[0],
                    "group_1_name": groups[1],
                    "group_1_definition": groups[2],
                    "in_discovery_corpus": key in discovery_experiments,
                    "independence": (
                        "same BugSigDB corpus; not independent of the discovery analyses"
                    ),
                },
            )
    return [seen[key] for key in sorted(seen)]


def _classification_counts(markers: Sequence[Mapping[str, object]]) -> dict[str, dict[str, int]]:
    """Candidate counts per classification, overall and for human and discovery subsets."""
    counts: dict[str, dict[str, int]] = {}
    for row in markers:
        bucket = counts.setdefault(
            str(row["classification"]), {"experiments": 0, "human": 0, "in_discovery_corpus": 0}
        )
        bucket["experiments"] += 1
        bucket["human"] += "Homo sapiens" in str(row["host_species"])
        bucket["in_discovery_corpus"] += bool(row["in_discovery_corpus"])
    return dict(sorted(counts.items()))


def inferred_function_traits(profile_path: Path, nodes_path: Path) -> list[dict[str, object]]:
    """Retained trait columns whose KG-Microbe label concerns SCFA or fermentation.

    A column can be qualified (``METPO:2000006|CHEBI:17272``); the qualifier's label
    says what the annotation asserts, and "uses as carbon source" is consumption, not
    production. Both labels are reported so the direction of the claim is explicit.
    """
    with profile_path.open(encoding="utf-8", newline="") as handle:
        header = next(csv.reader(handle, delimiter="\t"))
    columns = [
        c for c in header if c not in PROFILE_ID_COLUMNS and MEDIADIVE_MARKER not in c.lower()
    ]
    wanted = {part for column in columns for part in column.split("|")}
    labels: dict[str, str] = {}
    if nodes_path.exists():
        with nodes_path.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row["id"] in wanted:
                    labels[row["id"]] = row["name"]
                    if len(labels) == len(wanted):
                        break
    rows: list[dict[str, object]] = []
    for column in columns:
        term = column.rsplit("|", 1)[-1]
        qualifier = column.split("|", 1)[0] if "|" in column else ""
        label = labels.get(term, "")
        if not FUNCTION_LABEL.search(label):
            continue
        rows.append(
            {
                "column": column,
                "term": term,
                "label": label,
                "qualifier": qualifier,
                "qualifier_label": labels.get(qualifier, ""),
                "evidence_kind": "taxonomically inferred annotation, not a measurement",
            }
        )
    return rows


@dataclass(frozen=True)
class ExternalCohort:
    """One source study of the external collection, as its own tables describe it."""

    study_id: str
    reference_title: str
    year: str
    disease_name: str
    abbreviation: str
    accession: str
    region: str
    antibiotic_usage: str
    probiotics_usage: str
    medication_usage: str
    projects: str
    samples: int
    disease_samples: int
    control_samples: int
    role: str
    samples_with_age: int
    samples_with_bmi: int


def _first(pattern: str, text: str) -> str:
    match = re.search(pattern, text, re.IGNORECASE)
    return match.group(1).strip() if match else ""


def parse_supplement(supplement_text: str) -> dict[str, dict[str, str]]:
    """Read study titles, diseases, accessions and usage flags from the table text.

    The supplement is a PDF whose rows wrap, so each field is taken only where it
    parses unambiguously and anything else is left blank rather than guessed. The
    study identifier itself carries the year, which is used in preference to a year
    read out of a wrapped row.
    """
    fields: defaultdict[str, dict[str, str]] = defaultdict(dict)
    heading = re.compile(r"^\s*([A-Z][A-Za-z]+?[A-Z]?_\d{4}[ab]?)\s{2,}(.*)$")
    numeric = re.compile(r"\d+\.\d{3}")
    section = "disease"
    current: str | None = None
    for line in supplement_text.splitlines():
        if line.strip().startswith("Study ID") and "Reference" in line:
            section, current = "reference", None
            continue
        if "Sequencing platform" in line:
            section, current = "accession", None
            continue
        if section == "reference" and (numeric.search(line) or "Comparison ID" in line):
            # A later supplementary table has started.
            section, current = "done", None
            continue
        match = heading.match(line)
        rest = match.group(2) if match else line
        if match:
            current = match.group(1)
            if section == "reference":
                fields[current]["reference_title"] = rest.strip()
                continue
            if section == "disease":
                parts = re.split(r"\s{2,}", rest.strip())
                if len(parts) >= 7:
                    fields[current]["disease_name"] = parts[0]
                    fields[current]["abbreviation"] = parts[2]
                continue
        if current is None:
            continue
        if section == "reference" and re.search(r"\S\s{4,}\S", rest.strip()):
            # Titles are single-spaced; wide gaps mean a table header.
            section, current = "done", None
            continue
        if section == "reference" and rest.strip() and not rest.strip().startswith("Supplementary"):
            fields[current]["reference_title"] = (
                fields[current].get("reference_title", "") + " " + rest.strip()
            ).strip()
        elif section == "accession":
            accession = _first(r"((?:PRJ[A-Z]{2}|ERP|CNP)\d+)", rest)
            if accession:
                fields[current]["accession"] = accession
            region = _first(r"\s{2,}([A-Z][a-z]+),?\s{2,}(?:PRJ|ERP|CNP|\()", line)
            if region:
                fields[current]["region"] = region
            usage = re.findall(r"\s{2,}(Y|N|no data)(?=\s{2,})", line)
            if len(usage) >= 3:
                fields[current]["antibiotic_usage"] = usage[0]
                fields[current]["probiotics_usage"] = usage[1]
                fields[current]["medication_usage"] = usage[2]
    return dict(fields)


def summarize_sun2024(sample_group_path: Path, supplement_path: Path) -> list[ExternalCohort]:
    """Join the sample table's counts to the supplement's study descriptions."""
    described = parse_supplement(supplement_path.read_text(encoding="utf-8", errors="replace"))
    groups: defaultdict[str, list[dict[str, str]]] = defaultdict(list)
    for row in csv.DictReader(sample_group_path.open(encoding="utf-8", newline=""), delimiter="\t"):
        groups[row["Project_1"]].append(row)
    cohorts = []
    for study_id in sorted(set(described) | set(groups)):
        rows = groups.get(study_id, [])
        counts = Counter(row["Group"] for row in rows)
        described_row = described.get(study_id, {})
        cohorts.append(
            ExternalCohort(
                study_id=study_id,
                reference_title=described_row.get("reference_title", ""),
                year=_first(r"_(\d{4})", study_id),
                disease_name=described_row.get("disease_name", ""),
                abbreviation=described_row.get("abbreviation", ""),
                accession=described_row.get("accession", ""),
                region=described_row.get("region", ""),
                antibiotic_usage=described_row.get("antibiotic_usage", ""),
                probiotics_usage=described_row.get("probiotics_usage", ""),
                medication_usage=described_row.get("medication_usage", ""),
                projects=";".join(sorted({row["Project"] for row in rows})),
                samples=len(rows),
                disease_samples=counts.get("Disease", 0),
                control_samples=counts.get("Control", 0),
                role=";".join(sorted({row["Project_2"] for row in rows})),
                samples_with_age=sum(_clean(row["Age"]) != "" for row in rows),
                samples_with_bmi=sum(_clean(row["BMI"]) != "" for row in rows),
            )
        )
    return cohorts


def read_external_cohorts(path: Path) -> list[ExternalCohort]:
    with path.open(encoding="utf-8", newline="") as handle:
        return [
            ExternalCohort(
                study_id=row["study_id"],
                reference_title=row["reference_title"],
                year=row["year"],
                disease_name=row["disease_name"],
                abbreviation=row["abbreviation"],
                accession=row["accession"],
                region=row["region"],
                antibiotic_usage=row["antibiotic_usage"],
                probiotics_usage=row["probiotics_usage"],
                medication_usage=row["medication_usage"],
                projects=row["projects"],
                samples=int(row["samples"]),
                disease_samples=int(row["disease_samples"]),
                control_samples=int(row["control_samples"]),
                role=row["role"],
                samples_with_age=int(row["samples_with_age"]),
                samples_with_bmi=int(row["samples_with_bmi"]),
            )
            for row in csv.DictReader(handle, delimiter="\t")
        ]


def _title_key(title: str) -> str:
    return re.sub(r"[^a-z0-9]", "", title.lower())


def title_agreement(left: str, right: str) -> float:
    """Character-level similarity of two titles after removing case, spaces and punctuation."""
    a, b = _title_key(left), _title_key(right)
    if not a or not b:
        return 0.0
    matcher = difflib.SequenceMatcher(None, a, b, autojunk=False)
    if matcher.real_quick_ratio() < TITLE_AGREEMENT or matcher.quick_ratio() < TITLE_AGREEMENT:
        return 0.0
    return matcher.ratio()


def external_overlap(
    cohorts: Iterable[ExternalCohort], studies: Mapping[str, StudyRecord], discovery: set[str]
) -> list[dict[str, object]]:
    """Match each external source study to BugSigDB studies by publication title.

    The same paper curated in BugSigDB is the same cohort, so a title match is an
    overlap. A study without a title in the external table (a cohort first
    published by the external paper itself) cannot be assessed this way. No match
    is not proof of independence: BugSigDB may hold the study under a different
    title, or a later paper may reuse the same participants.
    """
    rows: list[dict[str, object]] = []
    titles = [
        (record, text)
        for record in studies.values()
        for text in record.texts.get("Title", [])
        if text
    ]
    for cohort in cohorts:
        matches: dict[str, tuple[StudyRecord, float]] = {}
        if cohort.reference_title:
            for record, text in titles:
                score = title_agreement(cohort.reference_title, text)
                if score >= TITLE_AGREEMENT:
                    matches[record.study_id] = (record, score)
        basis = "title" if matches else ""
        found = [record for record, _score in matches.values()]
        rows.append(
            {
                "external_study": cohort.study_id,
                "reference_title": cohort.reference_title,
                "disease_name": cohort.disease_name,
                "year": cohort.year,
                "accession": cohort.accession,
                "samples": cohort.samples,
                "role_in_source": cohort.role,
                "match_basis": basis,
                "matched_bugsigdb_studies": ";".join(sorted(matches)),
                "matched_pmids": ";".join(sorted(r.pmid for r in found if r.pmid)),
                "best_title_agreement": max((score for _r, score in matches.values()), default=0.0),
                "matches_in_discovery_corpus": ";".join(
                    sorted(r.study_id for r in found if r.study_id in discovery)
                ),
                "overlap_status": (
                    "not_assessable_without_title"
                    if not cohort.reference_title
                    else "overlaps_discovery_corpus"
                    if any(r.study_id in discovery for r in found)
                    else "curated_in_bugsigdb_outside_discovery_corpus"
                    if found
                    else "no_title_match_found"
                ),
                "measured_inflammation_or_scfa_available": False,
                "sample_level_covariates": (
                    f"age {cohort.samples_with_age}/{cohort.samples} samples, "
                    f"BMI {cohort.samples_with_bmi}/{cohort.samples}"
                ),
            }
        )
    return rows


def _format(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value)


def _write(path: Path, rows: Sequence[Mapping[str, object]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        if rows:
            writer.writerow(list(rows[0]))
            writer.writerows([_format(v) for v in row.values()] for row in rows)
    return path


def _sha256(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def write_external_snapshot(
    sample_group_path: Path, supplement_path: Path, output: Path = EXTERNAL_COHORTS_PATH
) -> Path:
    """Regenerate the committed cohort snapshot from the pinned sources."""
    digest = _sha256(sample_group_path)
    if digest != SUN2024_SAMPLE_GROUP_SHA256:
        raise ValueError(f"{sample_group_path} has SHA-256 {digest}, not the pinned table")
    return _write(
        output, [asdict(c) for c in summarize_sun2024(sample_group_path, supplement_path)]
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile-input", type=Path, default=DEFAULT_PROFILE_PATH)
    parser.add_argument("--experiments-input", type=Path, default=DEFAULT_EXPERIMENTS_PATH)
    parser.add_argument("--full-dump-input", type=Path, default=DEFAULT_FULL_DUMP_PATH)
    parser.add_argument("--nodes-input", type=Path, default=DEFAULT_NODES_PATH)
    parser.add_argument("--external-cohorts", type=Path, default=EXTERNAL_COHORTS_PATH)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--snapshot-from",
        type=Path,
        help=f"Regenerate the external cohort snapshot from a copy of {SUN2024_SAMPLE_GROUP_URL}",
    )
    parser.add_argument(
        "--snapshot-supplement",
        type=Path,
        help=f"Supplementary Table 1 as text, from {SUN2024_SUPPLEMENT_URL}",
    )
    args = parser.parse_args(argv)

    if args.snapshot_from is not None:
        if args.snapshot_supplement is None:
            parser.error("--snapshot-from also needs --snapshot-supplement")
        print(
            write_external_snapshot(
                args.snapshot_from, args.snapshot_supplement, args.external_cohorts
            )
        )
        return 0

    profiles, stats = load_profiles(
        args.profile_input, args.experiments_input, TABLE1_TRAITS, reject_missing_counts=True
    )
    discovery = {p.study_id for p in profiles}
    discovery_experiments = {(p.study_id, p.experiment_id) for p in profiles}
    studies = read_studies(args.full_dump_input)
    evidence = evidence_rows(studies, discovery)
    summary = summary_rows(evidence)
    markers = marker_split_candidates(args.full_dump_input, discovery_experiments)
    functions = inferred_function_traits(args.profile_input, args.nodes_input)
    overlap = external_overlap(read_external_cohorts(args.external_cohorts), studies, discovery)

    out = args.output_dir
    written = [
        _write(out / "evidence_availability_by_study.tsv", evidence),
        _write(out / "evidence_availability_summary.tsv", summary),
        _write(out / "marker_split_candidates.tsv", markers),
        _write(out / "inferred_function_traits.tsv", functions),
        _write(out / "external_cohort_overlap.tsv", overlap),
    ]
    provenance = {
        "version": VERSION,
        "bugsigdb_vintage": read_bugsigdb_vintage(args.full_dump_input),
        "full_dump_sha256": _sha256(args.full_dump_input),
        "profile_sha256": _sha256(args.profile_input),
        "nodes_available": args.nodes_input.exists(),
        "external_cohorts": str(args.external_cohorts),
        "external_cohorts_sha256": _sha256(args.external_cohorts),
        "external_source": SUN2024_SAMPLE_GROUP_URL,
        "external_source_sha256": SUN2024_SAMPLE_GROUP_SHA256,
        "discovery_studies": len(discovery),
        "discovery_experiments": stats.kept_experiments,
        "marker_split_candidates": _classification_counts(markers),
        "notes": [
            "a mention in curated metadata points to a source paper; it is not data",
            "not mentioned does not mean not measured",
            "marker-stratified experiments are part of the discovery corpus",
            "external overlap candidates require manual confirmation",
        ],
    }
    provenance_path = out / "validation_evidence_provenance.json"
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    written.append(provenance_path)
    for row in summary:
        print(row)
    for path in written:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
