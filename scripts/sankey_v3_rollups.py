"""Shared display-label rollups for v3 condition/body-site Sankey maps."""

from __future__ import annotations


def clean_text(value: object) -> str:
    """Normalize missing or empty display values."""
    text = str(value or "").strip()
    if not text or text.lower() == "nan" or text.upper() == "NA":
        return "unknown"
    return text


def has_any(text: str, terms: tuple[str, ...]) -> bool:
    """Return whether any lowercase term appears in text."""
    return any(term in text for term in terms)


def condition_display_label(value: object) -> tuple[str, str]:
    """Return a readable v3 condition display label and the applied rule."""
    label = clean_text(value)
    lower = label.lower()
    if label == "unknown":
        return "unknown", "unknown condition"

    # Keep transplant/GVHD concepts unmerged per v3 requirements.
    if has_any(lower, ("transplant", "graft", "gvhd")):
        return label, "kept raw transplant/graft label"

    if has_any(
        lower,
        (
            "response to ",
            "antimicrobial agent",
            "antibiotic",
            "antiviral",
            "treatment",
            "diet",
            "high fat",
            "ketogenic",
            "breastfeeding duration",
            "air pollution",
            "drug",
            "clinical treatment",
        ),
    ):
        return label, "kept raw exposure/intervention label"

    if has_any(
        lower,
        (
            "measurement",
            "timepoint",
            "time",
            "age",
            "ethnic group",
            "socioeconomic status",
            "species design",
            "library preparation",
            "population",
            "exposure",
            "health trait",
            "microbiome",
            "physical activity",
            "body fat distribution",
            "viral load",
        ),
    ):
        return label, "kept raw metadata/cohort-factor label"

    if has_any(lower, ("crohn", "ulcerative colitis", "inflammatory bowel disease")):
        return "Inflammatory bowel disease", "ibd rollup"

    if has_any(lower, ("colorectal cancer", "colorectal carcinoma", "colorectal adenoma")):
        return "Colorectal neoplasia", "colorectal neoplasia rollup"

    if has_any(lower, ("constipation", "diarrhea", "irritable bowel syndrome")):
        return "GI symptoms / motility", "gi-symptom rollup"

    if has_any(
        lower,
        (
            "obesity",
            "diabetes",
            "prediabetes",
            "metabolic syndrome",
            "polycystic ovary syndrome",
            "non-alcoholic fatty liver disease",
        ),
    ):
        return "Metabolic / endocrine", "metabolic/endocrine rollup"

    if has_any(lower, ("autism", "attention deficit hyperactivity disorder")):
        return "Autism / neurodevelopmental", "autism/neurodevelopmental rollup"

    if has_any(lower, ("ischemic stroke", "stroke")):
        return "Stroke", "stroke rollup"

    if has_any(lower, ("major depressive disorder", "depressive disorder")):
        return "Depressive disorder", "depressive-disorder rollup"

    if has_any(lower, ("periodontitis", "gingivitis")):
        return "Periodontal disease", "periodontal-disease rollup"

    return label, "kept raw"


def body_site_display_label(body_site: object, body_site_group: object) -> tuple[str, str]:
    """Return a readable v3 body-site display label and the applied rule."""
    raw_site = clean_text(body_site)
    raw_group = clean_text(body_site_group)
    text = f"{raw_site} {raw_group}".lower()

    if raw_site == "unknown" and raw_group == "unknown":
        return "Unknown body site", "unknown body-site rollup"

    if has_any(text, ("feces", "faeces", "stool")):
        return "Feces", "feces/stool rollup"

    if has_any(text, ("skin", "cutaneous", "axilla", "interdigital")):
        return "Skin", "skin rollup"

    if has_any(
        text,
        (
            "nasal",
            "nasopharynx",
            "nasopharyngeal",
            "nose epithelium",
            "middle nasal meatus",
            "inferior nasal concha",
            "oropharynx",
            "oropharyngeal",
            "throat",
            "pharynx",
            "upper respiratory",
        ),
    ):
        return "Upper airway / nasal", "upper-airway/nasal rollup"

    if has_any(text, ("sputum", "lung", "bronchoalveolar", "lower respiratory")):
        return "Lower airway / lung", "lower-airway/lung rollup"

    if has_any(
        text,
        (
            "saliva",
            "mouth",
            "oral",
            "dental plaque",
            "subgingival",
            "supragingival",
            "gingival",
            "tongue",
            "buccal",
            "tooth",
            "teeth",
        ),
    ):
        return "Oral cavity / dental", "oral/dental rollup"

    if has_any(text, ("colon", "caecum", "cecum", "rectum", "ileum", "intestin", "gut")):
        return "Lower GI tract", "lower-gi rollup"

    if has_any(text, ("vagina", "vaginal", "cervix", "cervical", "uter", "endomet")):
        return "Vaginal / reproductive", "vaginal/reproductive rollup"

    if has_any(text, ("blood", "serum", "plasma")):
        return "Blood", "blood rollup"

    if has_any(text, ("urine", "urinary", "bladder")):
        return "Urinary tract", "urinary rollup"

    if raw_group.lower() not in {"other", "unknown"}:
        return raw_group, "kept curated body-site group"

    return "Other body sites", "other body-site rollup"
