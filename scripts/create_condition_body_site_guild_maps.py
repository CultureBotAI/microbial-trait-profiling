"""Create condition-to-body-site guild z-score maps."""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from textwrap import shorten, wrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, PathPatch
from matplotlib.path import Path as MplPath

from mtp.composition_denominators import (
    DENOMINATOR_DEFINITION,
    load_retained_profiles,
    trait_fractions,
)
from sankey_v3_rollups import (
    body_site_display_label,
    condition_display_label,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
VIS_DIR = REPO_ROOT / "trait-visualization"
INVESTIGATION_DIR = VIS_DIR / "trait-profile-investigation"

MATRIX_PATH = VIS_DIR / "study_trait_profile_clustered_heatmap_matrix.tsv"
METADATA_PATH = VIS_DIR / "study_trait_profile_study_metadata.tsv"
TRAIT_METADATA_PATH = VIS_DIR / "study_trait_profile_trait_metadata.tsv"
PROFILE_MATRIX_PATH = REPO_ROOT / "data/derived/study_trait_profiles_wide.tsv"

MIN_FOCUS_PROFILES = 3
MIN_STUDIES = 2
MAX_DISPLAY_EDGES = 18
MAX_DISPLAY_CONDITIONS = 12
MAX_LABEL_WIDTH = 36
MAX_OTHER_SITES_PER_CONDITION = 3
V3_MAX_DISPLAY_EDGES = 32
V3_MAX_DISPLAY_CONDITIONS = 12
V3_MIN_FOCUS_ENRICHED_PROFILES = 3
V3_BACKGROUND_EDGE_Y_OFFSET = -0.007
V3_GROUPED_MODE = "grouped"
V3_INDIVIDUAL_MODE = "individual"

SLATE_50 = "#f8fafc"
SLATE_200 = "#e2e8f0"
SLATE_400 = "#94a3b8"
SLATE_500 = "#64748b"
SLATE_700 = "#334155"
SLATE_900 = "#0f172a"
EDGE_NEUTRAL = "#475569"


def wrapped_plot_text(value: str, width: int = 150) -> str:
    """Wrap long figure annotations so they do not widen the saved canvas."""
    return "\n".join(wrap(value, width=width))


@dataclass(frozen=True)
class MapConfig:
    """Configuration for one condition/body-site guild map."""

    key: str
    title: str
    focus_label: str
    comparator_label: str
    focus_direction: str
    comparator_direction: str
    output_dir: Path
    output_stem: str
    excluded_condition_labels: frozenset[str]

    @property
    def guild_selection_path(self) -> Path:
        return self.output_dir / "guild_trait_selection.tsv"

    @property
    def edge_tsv_path(self) -> Path:
        return self.output_dir / f"{self.output_stem}_edges.tsv"

    @property
    def edge_v2_tsv_path(self) -> Path:
        return self.output_dir / f"{self.output_stem}_v2_edges.tsv"

    @property
    def edge_v3_tsv_path(self) -> Path:
        return self.output_dir / f"{self.output_stem}_v3_edges.tsv"

    @property
    def condition_rollup_v3_tsv_path(self) -> Path:
        return self.output_dir / f"{self.output_stem}_v3_condition_rollups.tsv"

    @property
    def body_site_rollup_v3_tsv_path(self) -> Path:
        return self.output_dir / f"{self.output_stem}_v3_body_site_rollups.tsv"

    @property
    def sampled_site_tsv_path(self) -> Path:
        return self.output_dir / f"{self.output_stem}_sampled_body_sites.tsv"

    @property
    def sampled_site_v2_tsv_path(self) -> Path:
        return self.output_dir / f"{self.output_stem}_v2_sampled_body_sites.tsv"

    @property
    def png_path(self) -> Path:
        return self.output_dir / f"{self.output_stem}.png"

    @property
    def png_v2_path(self) -> Path:
        return self.output_dir / f"{self.output_stem}_v2.png"

    @property
    def png_v3_path(self) -> Path:
        return self.output_dir / f"{self.output_stem}_v3.png"

    @property
    def edge_v3_individual_tsv_path(self) -> Path:
        return self.output_dir / f"{self.output_stem}_v3_individual_edges.tsv"

    @property
    def edge_v3_grouped_tsv_path(self) -> Path:
        return self.output_dir / f"{self.output_stem}_v3_grouped_edges.tsv"

    @property
    def condition_rollup_v3_grouped_tsv_path(self) -> Path:
        return self.output_dir / f"{self.output_stem}_v3_grouped_condition_rollups.tsv"

    @property
    def body_site_rollup_v3_grouped_tsv_path(self) -> Path:
        return self.output_dir / f"{self.output_stem}_v3_grouped_body_site_rollups.tsv"

    @property
    def png_v3_individual_path(self) -> Path:
        return self.output_dir / f"{self.output_stem}_v3_individual.png"

    @property
    def png_v3_grouped_path(self) -> Path:
        return self.output_dir / f"{self.output_stem}_v3_grouped.png"


FACULTATIVE_CONFIG = MapConfig(
    key="facultative_aerobe_increased",
    title="Condition-to-Body-Site Map: Aerobe v2 Trait Enrichment in Increased Profiles",
    focus_label="increased",
    comparator_label="decreased",
    focus_direction="inc",
    comparator_direction="dec",
    output_dir=INVESTIGATION_DIR / "aerobe-enrichment-v2",
    output_stem="condition_body_site_guild_map",
    excluded_condition_labels=frozenset(
        {
            "Body fat distribution",
            "Clinical treatment",
            "Exposure",
            "Glucose",
            "Health study participation",
            "Health trait",
            "Microbiome",
            "Physical activity",
            "Population",
            "Socioeconomic status",
            "Time",
            "Timepoint",
            "Tumor stage",
        }
    ),
)

ANAEROBE_CONFIG = MapConfig(
    key="anaerobe_decreased",
    title="Condition-to-Body-Site Map: Anaerobe Trait Enrichment in Decreased Profiles",
    focus_label="decreased",
    comparator_label="increased",
    focus_direction="dec",
    comparator_direction="inc",
    output_dir=INVESTIGATION_DIR / "anaerobe-depletion-v2",
    output_stem="condition_body_site_guild_map",
    excluded_condition_labels=frozenset(
        {
            "Body fat distribution",
            "Clinical treatment",
            "Exposure",
            "Glucose",
            "Health study participation",
            "Health trait",
            "Microbiome",
            "Physical activity",
            "Population",
            "Socioeconomic status",
            "Time",
            "Timepoint",
            "Tumor stage",
        }
    ),
)


def clean_text(value: object) -> str:
    """Normalize missing or empty display values."""
    text = str(value or "").strip()
    if not text or text.lower() == "nan" or text.upper() == "NA":
        return "unknown"
    return text


def compact_label(value: object, width: int = MAX_LABEL_WIDTH) -> str:
    """Return a compact label for node text."""
    return shorten(clean_text(value), width=width, placeholder="...")


def wrap_label(value: object, width: int = MAX_LABEL_WIDTH) -> str:
    """Wrap labels without letting them dominate the diagram."""
    label = compact_label(value, width=width)
    return "\n".join(wrap(label, width=max(12, width // 2))) or "unknown"


def is_excluded_condition_label(value: object, excluded: frozenset[str]) -> bool:
    """Return whether a condition label is too broad for display."""
    label = clean_text(value)
    lower_label = label.lower()
    broad_terms = ("exposure", "measurement", "population")
    return (
        label in excluded or label == "unknown" or any(term in lower_label for term in broad_terms)
    )


def is_display_body_site(value: object) -> bool:
    """Return whether a body-site group should be shown as an anatomical node."""
    label = clean_text(value).lower()
    return label not in {"unknown", "other"}


def join_unique(values: pd.Series) -> str:
    """Join unique non-empty values while preserving first-seen order."""
    seen: list[str] = []
    for value in values:
        text = clean_text(value)
        if text != "unknown" and text not in seen:
            seen.append(text)
    return "; ".join(seen)


def split_joined_values(value: object) -> list[str]:
    """Split a semicolon-delimited value emitted by join_unique."""
    text = clean_text(value)
    if text == "unknown":
        return []
    return [part.strip() for part in text.split(";") if part.strip()]


def read_tsv(path: Path, *, index_col: int | str | None = None) -> pd.DataFrame:
    """Read a TSV file with stable defaults."""
    return pd.read_csv(path, sep="\t", index_col=index_col)


def load_study_heatmap() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load profile-level z-score matrix and metadata with readable trait names."""
    matrix = read_tsv(MATRIX_PATH, index_col=0)
    matrix = matrix.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    metadata = read_tsv(METADATA_PATH, index_col=0)
    trait_metadata = read_tsv(TRAIT_METADATA_PATH, index_col=0)

    if len(matrix) != len(metadata):
        raise ValueError("Heatmap matrix and metadata have different profile row counts")
    if (
        "denominator_definition" not in metadata.columns
        or not metadata["denominator_definition"].eq(DENOMINATOR_DEFINITION).all()
    ):
        raise ValueError(
            "Heatmap metadata is missing the harmonized non-MediaDive denominator marker"
        )
    matrix_keys = matrix.index.to_series().map(lambda value: str(value).split(" | ", 1)[0])
    matrix_pair_keys = matrix.index.to_series().map(
        lambda value: "/".join(part.strip() for part in str(value).split(" | ")[:2])
    )
    expected_keys = (
        metadata[["study_id", "experiment_id"]].astype(str).agg("/".join, axis=1).to_numpy()
    )
    if not (
        matrix.index.equals(metadata.index)
        or np.array_equal(matrix_keys.to_numpy(), expected_keys)
        or np.array_equal(matrix_pair_keys.to_numpy(), expected_keys)
    ):
        raise ValueError("Heatmap matrix row keys do not match profile metadata")
    matrix.index = metadata.index

    if "trait_id" in trait_metadata.columns:
        trait_metadata = trait_metadata.set_index("trait_id", drop=False)
    column_to_id = {
        row["display_label"]: str(trait_id) for trait_id, row in trait_metadata.iterrows()
    }
    column_labels = {
        row["display_label"]: row["object_name"] for _, row in trait_metadata.iterrows()
    }
    if not set(matrix.columns).issubset(column_to_id):
        raise ValueError("Heatmap matrix has columns missing from trait metadata")
    counts = load_retained_profiles(PROFILE_MATRIX_PATH, metadata)
    trait_ids = [column_to_id[column] for column in matrix.columns]
    expected_comp = trait_fractions(counts).loc[metadata.index, trait_ids]
    expected_scaled = (expected_comp * 1000).apply(lambda column: column.map(math.log1p))
    expected_z = (
        (expected_scaled - expected_scaled.mean(axis=0))
        / expected_scaled.std(axis=0, ddof=0).replace(0, float("nan"))
    ).fillna(0.0)
    expected_z.columns = matrix.columns
    if not np.allclose(matrix.to_numpy(dtype=float), expected_z.to_numpy(dtype=float), atol=1e-9):
        raise ValueError(
            "Stale or incompatible study_trait_profile_clustered_heatmap_matrix.tsv; "
            "regenerate from the non-MediaDive raw profile matrix."
        )
    matrix = matrix.rename(columns=column_labels)
    return matrix, metadata


def load_guild_traits(config: MapConfig) -> list[str]:
    """Load the selected traits for one analysis."""
    selection = read_tsv(config.guild_selection_path)
    trait_column = (
        "profile_matrix_object_name"
        if "profile_matrix_object_name" in selection.columns
        else "trait_name"
    )
    traits = [clean_text(value) for value in selection[trait_column]]
    return [trait for trait in traits if trait != "unknown"]


def add_guild_scores(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    config: MapConfig,
) -> pd.DataFrame:
    """Attach row-wise guild z-score means to profile metadata."""
    guild_traits = load_guild_traits(config)
    available_traits = [trait for trait in guild_traits if trait in matrix.columns]
    if not available_traits:
        msg = f"No selected traits for {config.key} were found in the profile matrix."
        raise ValueError(msg)

    scored = metadata.copy()
    scored["profile_id"] = scored.index
    scored["guild_z_mean"] = matrix[available_traits].mean(axis=1)
    scored["guild_positive_trait_count"] = matrix[available_traits].gt(0.0).sum(axis=1)
    scored["has_guild_enriched_trait"] = scored["guild_positive_trait_count"] > 0
    scored["condition"] = scored["condition"].map(clean_text)
    scored["body_site"] = scored["body_site"].map(clean_text)
    scored["body_site_group"] = scored["body_site_group"].map(clean_text)
    return scored


def nan_mean(series: pd.Series) -> float:
    """Return a float mean, preserving NaN for empty groups."""
    if series.empty:
        return math.nan
    return float(series.mean())


def condition_body_site_source(scored: pd.DataFrame, config: MapConfig) -> pd.DataFrame:
    """Return condition/body-site profile rows eligible for map summarization."""
    condition_is_allowed = ~scored["condition"].map(
        lambda value: is_excluded_condition_label(value, config.excluded_condition_labels)
    )
    return scored[condition_is_allowed & scored["body_site_group"].map(is_display_body_site)].copy()


def summarize_sampled_body_sites(
    scored: pd.DataFrame,
    config: MapConfig,
) -> pd.DataFrame:
    """Summarize every sampled condition/body-site pair before display filtering."""
    source = condition_body_site_source(scored, config)

    records: list[dict[str, object]] = []
    for (condition, body_site_group), group in source.groupby(
        ["condition", "body_site_group"],
        sort=False,
        dropna=False,
    ):
        focus = group[group["inc_dec"] == config.focus_direction]
        comparator = group[group["inc_dec"] == config.comparator_direction]

        focus_mean = nan_mean(focus["guild_z_mean"])
        comparator_mean = nan_mean(comparator["guild_z_mean"])
        contrast = focus_mean - comparator_mean if np.isfinite(comparator_mean) else math.nan
        records.append(
            {
                "analysis": config.key,
                "condition": clean_text(condition),
                "body_site_group": clean_text(body_site_group),
                "profile_count": len(group),
                "study_count": group["study_id"].nunique(),
                "study_ids": join_unique(group["study_id"]),
                "focus_direction": config.focus_label,
                "comparator_direction": config.comparator_label,
                "focus_profile_count": len(focus),
                "focus_study_count": focus["study_id"].nunique(),
                "focus_study_ids": join_unique(focus["study_id"]),
                "comparator_profile_count": len(comparator),
                "comparator_study_count": comparator["study_id"].nunique(),
                "comparator_study_ids": join_unique(comparator["study_id"]),
                "increased_profile_count": int((group["inc_dec"] == "inc").sum()),
                "decreased_profile_count": int((group["inc_dec"] == "dec").sum()),
                "mean_guild_z": nan_mean(group["guild_z_mean"]),
                "increased_mean_guild_z": nan_mean(
                    group.loc[group["inc_dec"] == "inc", "guild_z_mean"]
                ),
                "decreased_mean_guild_z": nan_mean(
                    group.loc[group["inc_dec"] == "dec", "guild_z_mean"]
                ),
                "focus_mean_guild_z": focus_mean,
                "comparator_mean_guild_z": comparator_mean,
                "focus_minus_comparator_guild_z": contrast,
            }
        )

    if not records:
        return pd.DataFrame()

    return pd.DataFrame.from_records(records)


def summarize_condition_body_site_edges(sampled_sites: pd.DataFrame) -> pd.DataFrame:
    """Build display-ranking edges from sampled condition/body-site pairs."""
    if sampled_sites.empty:
        return pd.DataFrame()

    edges = sampled_sites[sampled_sites["focus_profile_count"] > 0].copy()
    if edges.empty:
        return pd.DataFrame()

    edges["rank_signal"] = edges["focus_minus_comparator_guild_z"].fillna(
        edges["focus_mean_guild_z"]
    )
    edges["rank_score"] = edges["rank_signal"] * np.log1p(edges["focus_profile_count"])
    edges["displayed_in_map"] = False
    return edges.sort_values(
        ["rank_score", "rank_signal", "focus_profile_count", "study_count"],
        ascending=[False, False, False, False],
    )


def select_display_edges(edges: pd.DataFrame) -> pd.DataFrame:
    """Choose a readable set of high-signal supported edges for the PNG."""
    candidates = edges[
        (edges["focus_profile_count"] >= MIN_FOCUS_PROFILES)
        & (edges["study_count"] >= MIN_STUDIES)
        & (edges["focus_mean_guild_z"] > 0)
    ].copy()
    preferred = candidates[candidates["focus_minus_comparator_guild_z"] > 0].copy()
    if len(preferred) >= 6:
        candidates = preferred
    if candidates.empty:
        candidates = edges.head(MAX_DISPLAY_EDGES).copy()

    candidates = candidates.sort_values(
        ["rank_score", "rank_signal", "focus_profile_count", "study_count"],
        ascending=[False, False, False, False],
    )
    per_body_site = candidates.groupby("body_site_group", sort=False).head(2)
    selected_indices: list[int] = []
    selected_conditions: set[str] = set()

    def add_candidate(index: int, row: pd.Series) -> None:
        if index in selected_indices or len(selected_indices) >= MAX_DISPLAY_EDGES:
            return
        condition = str(row["condition"])
        condition_limit_reached = len(selected_conditions) >= MAX_DISPLAY_CONDITIONS
        if condition not in selected_conditions and condition_limit_reached:
            return
        selected_indices.append(index)
        selected_conditions.add(condition)

    for index, row in per_body_site.iterrows():
        add_candidate(index, row)
    for index, row in candidates.iterrows():
        add_candidate(index, row)

    display = candidates.loc[selected_indices].sort_values(
        ["rank_score", "rank_signal", "focus_profile_count", "study_count"],
        ascending=[False, False, False, False],
    )
    return display.copy()


def select_uncapped_display_edges(edges: pd.DataFrame) -> pd.DataFrame:
    """Return every supported high-signal edge for the uncapped v2 PNG."""
    candidates = edges[
        (edges["focus_profile_count"] >= MIN_FOCUS_PROFILES)
        & (edges["study_count"] >= MIN_STUDIES)
        & (edges["focus_mean_guild_z"] > 0)
    ].copy()
    return candidates.sort_values(
        ["rank_score", "rank_signal", "focus_profile_count", "study_count"],
        ascending=[False, False, False, False],
    )


def has_any_label_term(text: str, terms: tuple[str, ...]) -> bool:
    """Return whether any lowercase term appears in normalized label text."""
    return any(term in text for term in terms)


def v3_fully_grouped_condition_display_label(value: object) -> tuple[str, str]:
    """Return fully grouped condition labels for grouped v3 guild maps."""
    label = clean_text(value)
    lower = label.lower()
    if label == "unknown":
        return "Unknown condition", "unknown condition group"

    if lower == "healthy":
        return "Healthy reference", "healthy-reference condition group"

    if has_any_label_term(lower, ("crohn", "ulcerative colitis", "inflammatory bowel disease")):
        return "Inflammatory bowel disease", "ibd condition group"

    if has_any_label_term(
        lower,
        ("colorectal cancer", "colorectal carcinoma", "colorectal adenoma"),
    ):
        return "Colorectal neoplasia", "colorectal neoplasia condition group"

    if has_any_label_term(
        lower,
        (
            "periodontitis",
            "gingivitis",
            "peri-implantitis",
            "oral mucositis",
            "dental",
        ),
    ):
        return "Oral / periodontal disease", "oral/periodontal condition group"

    if has_any_label_term(
        lower,
        (
            "constipation",
            "diarrhea",
            "irritable bowel syndrome",
            "celiac",
            "gastroesophageal",
            "enteropathy",
            "pancreatitis",
            "biliary",
            "liver disease",
            "gastroenteritis",
            "colitis",
            "gastrointestinal disease",
            "gastrointestinal tract",
            "intestinal obstruction",
            "intestinal polyp",
            "peptic esophagitis",
            "lactose intolerance",
            "necrotizing enterocolitis",
            "ulcerative proctosigmoiditis",
            "abnormal stool",
            "colorectal health",
        ),
    ):
        return "GI / digestive disease", "gi/digestive condition group"

    if has_any_label_term(
        lower,
        (
            "macular",
            "glaucoma",
            "dry eye",
            "ophthalmopathy",
            "keratoconjunctivitis",
            "blepharitis",
            "retinopathy",
            "uveitis",
            "myopia",
            "meibomian",
        ),
    ):
        return "Ocular / eye disease", "ocular/eye condition group"

    if has_any_label_term(
        lower,
        (
            "obesity",
            "diabetes",
            "prediabetes",
            "metabolic syndrome",
            "polycystic ovary syndrome",
            "non-alcoholic fatty liver disease",
            "insulin",
            "metabolic",
            "endocrine",
            "fasting blood glucose",
            "glucose",
            "cholesterol",
            "lipoprotein",
            "phenylketonuria",
            "hypothyroidism",
            "thyroiditis",
            "hypercholesterolemia",
            "malnutrition",
        ),
    ):
        return "Metabolic / endocrine", "metabolic/endocrine condition group"

    if has_any_label_term(
        lower,
        (
            "pregnancy",
            "preterm",
            "premature birth",
            "cesarean",
            "maternal",
            "breastfeeding",
            "neonatal",
            "infant",
            "gestational",
            "postpartum",
            "fetal growth restriction",
            "delivery method",
            "delivery measurement",
            "preeclampsia",
        ),
    ):
        return "Pregnancy / neonatal", "pregnancy/neonatal condition group"

    if has_any_label_term(
        lower,
        (
            "waist circumference",
            "body mass index",
            "body weight",
            "increased body weight",
            "body fat",
            "visceral adipose",
            "growth rate",
            "postnatal growth",
            "short stature",
            "hypertrophy",
            "cachexia",
            "sarcopenia",
            "body composition",
        ),
    ):
        return "Body composition / growth", "body-composition/growth condition group"

    if has_any_label_term(lower, ("ischemic stroke", "stroke")):
        return "Stroke", "stroke condition group"

    if has_any_label_term(
        lower,
        (
            "autism",
            "attention deficit",
            "parkinson",
            "alzheimer",
            "multiple sclerosis",
            "neuromyelitis optica",
            "dementia",
            "lewy body",
            "cognitive impairment",
            "memory impairment",
            "epilepsy",
            "essential tremor",
            "amyotrophic lateral sclerosis",
            "huntington",
            "prion disease",
            "restless legs",
            "spinal cord injury",
            "brain hypoxia",
            "cerebral small vessel",
            "cerebral infarction",
            "encephalitis",
            "down syndrome",
            "diabetic neuropathy",
            "sleep",
            "insomnia",
            "neurologic",
            "neurodegenerative",
        ),
    ):
        return "Neurologic / cognitive disease", "neurologic/cognitive condition group"

    if has_any_label_term(
        lower,
        (
            "major depressive disorder",
            "depressive disorder",
            "unipolar depression",
            "postpartum depression",
            "bipolar disorder",
            "schizophrenia",
            "psychosis",
            "obsessive-compulsive",
            "anxiety",
            "mood disorder",
            "stress-related disorder",
            "mental or behavioural disorder",
            "anorexia nervosa",
            "eating behaviour",
            "behavior or behavioral disorder",
            "social interaction",
            "childhood trauma",
            "chronic fatigue syndrome",
            "fatigue",
        ),
    ):
        return "Psychiatric / behavioral health", "psychiatric/behavioral condition group"

    if has_any_label_term(
        lower,
        (
            "covid",
            "sars",
            "tuberculosis",
            "pneumonia",
            "hiv",
            "virus",
            "viral",
            "infection",
            "infectious",
            "helicobacter",
            "chlamydia",
            "zika",
            "bacterial disease",
            "hepatitis",
            "helminth",
            "ascariasis",
            "schistosomiasis",
            "trichuriasis",
            "clonorchiasis",
            "coccidiosis",
            "gonorrhea",
            "shigellosis",
            "influenza",
            "sepsis",
        ),
    ):
        return "Infectious disease", "infectious-disease condition group"

    if has_any_label_term(
        lower,
        (
            "asthma",
            "chronic obstructive pulmonary",
            "allergic rhinitis",
            "respiratory",
            "lung",
            "airway",
            "bronchiectasis",
            "wheezing",
            "pulmonary function",
            "pulmonary fibrosis",
            "pulmonary arterial hypertension",
            "rhinitis",
            "rhinosinusitis",
        ),
    ):
        return "Respiratory / airway disease", "respiratory/airway condition group"

    if has_any_label_term(
        lower,
        (
            "cancer",
            "carcinoma",
            "leukemia",
            "lymphoma",
            "myeloma",
            "melanoma",
            "neoplasm",
            "tumor",
            "adenocarcinoma",
            "hepatocellular",
            "glioma",
            "leukoplakia",
            "carcinogenicity",
            "pituitary gland adenoma",
            "peutz-jeghers polyp",
            "adenomatous colon polyp",
        ),
    ):
        return "Cancer / neoplasia", "cancer/neoplasia condition group"

    if has_any_label_term(
        lower,
        ("kidney", "nephro", "urinary", "bladder", "renal", "urolithiasis"),
    ):
        return "Renal / urinary disease", "renal/urinary condition group"

    if has_any_label_term(
        lower,
        (
            "ovarian",
            "cervical",
            "endometrial",
            "infertility",
            "reproductive",
            "vagina",
            "menopause",
            "endometriosis",
            "adenomyosis",
            "dysmenorrhea",
            "menstrual",
            "menopausal",
            "vulvovaginitis",
            "bacterial vaginosis",
            "mastitis",
            "sperm",
            "contraception",
            "medroxyprogesterone",
            "estradiol",
            "sex hormone",
            "puberty",
            "hormone replacement",
            "cervix erosion",
            "spontaneous abortion",
            "male homosexuality",
        ),
    ):
        return "Reproductive / hormonal health", "reproductive/hormonal condition group"

    if has_any_label_term(
        lower,
        (
            "cardiovascular",
            "hypertension",
            "myocardial infarction",
            "heart failure",
            "coronary artery",
            "atrial fibrillation",
            "patent ductus arteriosus",
            "congenital heart",
            "carotid artery",
            "syncope",
        ),
    ):
        return "Cardiovascular disease", "cardiovascular condition group"

    if has_any_label_term(
        lower,
        (
            "eczema",
            "dermatitis",
            "acne",
            "rosacea",
            "hidradenitis",
            "lichen planopilaris",
            "skin wound",
            "pruritus",
        ),
    ):
        return "Skin / hair disease", "skin/hair condition group"

    if has_any_label_term(
        lower,
        (
            "anemia",
            "thrombocythemia",
            "aplastic",
            "blood disease",
            "hematopoietic",
            "immunodeficiency",
            "chronic granulomatous",
            "cd4-positive",
        ),
    ):
        return "Blood / immune deficiency", "blood/immune-deficiency condition group"

    if has_any_label_term(
        lower,
        (
            "rheumatoid",
            "lupus",
            "behcet",
            "psoriasis",
            "vitiligo",
            "atopic eczema",
            "allergy",
            "allergic",
            "autoimmune",
            "inflammatory",
            "graves",
            "stevens-johnson",
            "scleroderma",
            "sjogren",
            "oxidative stress",
            "purpura",
            "granulomatosis",
            "polyangiitis",
            "antiphospholipid",
            "sapho",
            "pemphigus",
        ),
    ):
        return "Autoimmune / inflammatory disease", "autoimmune/inflammatory condition group"

    if has_any_label_term(
        lower,
        (
            "arthritis",
            "osteoarthritis",
            "spondylitis",
            "osteoporosis",
            "dermatomyositis",
            "bone mineral",
            "myostatin",
            "muscle",
            "fibromyalgia",
        ),
    ):
        return "Musculoskeletal / connective tissue", "musculoskeletal/connective condition group"

    if has_any_label_term(
        lower,
        (
            "otitis",
            "halitosis",
            "bitter taste",
            "body odor",
            "oral lichen planus",
            "mucositis",
            "rhinitis",
            "sinusitis",
        ),
    ):
        return "ENT / oral mucosal disease", "ent/oral-mucosal condition group"

    if has_any_label_term(lower, ("transplant", "graft", "gvhd")):
        return "Transplant / graft response", "transplant/graft condition group"

    if has_any_label_term(
        lower,
        (
            "antimicrobial agent",
            "antibiotic",
            "antiviral",
            "azithromycin",
        ),
    ):
        return "Antibiotic / antimicrobial exposure", "antibiotic/antimicrobial condition group"

    if has_any_label_term(
        lower,
        (
            "drug",
            "vaccine",
            "treatment",
            "response to",
            "clinical treatment",
        ),
    ):
        return "Medication / clinical treatment", "medication/clinical-treatment condition group"

    if has_any_label_term(
        lower,
        (
            "diet",
            "ketogenic",
            "high fat",
            "traditional chinese medicine",
            "nutraceutical",
            "curcumin",
            "vitamin",
            "fasting",
            "beverage",
            "fat intake",
            "iron",
            "equol",
            "lactobacillus",
        ),
    ):
        return "Diet / nutrition / supplements", "diet/nutrition/supplement condition group"

    if has_any_label_term(
        lower,
        (
            "smoking",
            "alcohol",
            "tobacco",
            "nicotine",
            "cannabis",
            "methamphetamine",
            "opioid",
            "substance-related",
        ),
    ):
        return "Smoking / alcohol / drug use", "smoking/alcohol/drug-use condition group"

    if has_any_label_term(
        lower,
        (
            "air pollution",
            "exposure",
            "environmental",
            "pesticide",
            "insecticide",
            "titanium dioxide",
            "silicosis",
            "chemical sensitivity",
            "light",
            "temperature",
            "high altitude",
            "sun exposure",
            "seasonality",
            "mosquito",
            "acetaldehyde",
            "hypoxia",
        ),
    ):
        return "Environmental / chemical exposure", "environmental/chemical condition group"

    if has_any_label_term(
        lower,
        (
            "age",
            "biological sex",
            "birth measurement",
            "developmental stage",
            "brain growth",
            "puberty",
            "aging",
            "development or differentiation",
            "breed",
        ),
    ):
        return "Age / development / sex", "age/development/sex condition group"

    if has_any_label_term(
        lower,
        (
            "population",
            "ethnic group",
            "socioeconomic",
            "social deprivation",
            "place of residence",
            "restricted to specific location",
            "health study participation",
            "mixed sex population",
        ),
    ):
        return "Population / social covariates", "population/social-covariate condition group"

    if has_any_label_term(
        lower,
        ("physical activity", "exercise", "lifestyle measurement"),
    ):
        return "Activity / lifestyle score", "activity/lifestyle-score condition group"

    if has_any_label_term(
        lower,
        (
            "autoantibody",
            "antibody",
            "igg4",
            "immune biomarker",
        ),
    ):
        return "Immune biomarker", "immune-biomarker condition group"

    if has_any_label_term(
        lower,
        (
            "bilirubin",
            "alanine",
            "liver fibrosis",
            "cholestasis",
            "gallstones",
            "bile",
        ),
    ):
        return "Liver / biliary biomarker", "liver/biliary-biomarker condition group"

    if has_any_label_term(
        lower,
        (
            "oxalate",
            "nitrite",
            "methionine",
            "creatine",
            "metabolite",
            "vitamin k",
        ),
    ):
        return "Metabolite biomarker", "metabolite-biomarker condition group"

    if has_any_label_term(
        lower,
        (
            "grey matter volume",
            "hippocampal fissure",
            "brain growth measurement",
        ),
    ):
        return "Brain morphology measurement", "brain-morphology condition group"

    if has_any_label_term(
        lower,
        (
            "periodontal measurement",
            "number of teeth",
            "root caries",
            "periodontal pocket",
            "gingival disease",
        ),
    ):
        return "Dental / periodontal measurement", "dental/periodontal-measurement group"

    if has_any_label_term(
        lower,
        (
            "disease staging",
            "illness severity",
            "disease prognosis",
            "disease progression",
            "remission",
            "survival time",
            "health-related quality of life",
            "risk factor",
        ),
    ):
        return "Disease activity / severity", "disease-activity/severity condition group"

    if has_any_label_term(
        lower,
        (
            "microbiome",
            "seasonal gut microbiome",
            "abnormal stool composition",
        ),
    ):
        return "Microbiome phenotype", "microbiome-phenotype condition group"

    if has_any_label_term(
        lower,
        (
            "extraction protocol",
            "library preparation",
            "nucleic acid extraction",
            "gene expression protocol",
        ),
    ):
        return "Extraction / library protocol", "extraction/library-protocol condition group"

    if has_any_label_term(
        lower,
        (
            "sample collection protocol",
            "freeze dried specimen",
            "transport",
        ),
    ):
        return "Sample collection / handling", "sample-collection/handling condition group"

    if has_any_label_term(
        lower,
        (
            "sampling site",
            "sampling time",
            "timepoint",
            "time",
        ),
    ):
        return "Sampling site / time protocol", "sampling-site/time condition group"

    if has_any_label_term(lower, ("species design",)):
        return "Host species / model design", "host-species/model-design condition group"

    if has_any_label_term(
        lower,
        (
            "species design",
            "growth condition",
            "stimulus or stress design",
            "phenotype",
            "vesicle membrane",
            "root",
            "peripheral neuron",
            "health trait",
        ),
    ):
        return "Experimental phenotype labels", "experimental-phenotype-label condition group"

    return "Other specific clinical labels", "other-specific-clinical-label condition group"


def v3_condition_display_label(value: object, config: MapConfig) -> tuple[str, str]:
    """Return the grouped v3 condition label for this map."""
    if config.key in {FACULTATIVE_CONFIG.key, ANAEROBE_CONFIG.key}:
        return v3_fully_grouped_condition_display_label(value)
    return condition_display_label(value)


def v3_body_site_display_label(
    body_site: object,
    body_site_group: object,
    config: MapConfig,
) -> tuple[str, str]:
    """Return the grouped v3 body-site label for this map."""
    label, rule = body_site_display_label(body_site, body_site_group)
    has_guild_specific_body_site_rollups = config.key in {
        FACULTATIVE_CONFIG.key,
        ANAEROBE_CONFIG.key,
    }
    if has_guild_specific_body_site_rollups and label == "Unknown body site":
        return "Missing body-site metadata", "missing body-site metadata group"
    if has_guild_specific_body_site_rollups and label == "Other body sites":
        raw_site = clean_text(body_site)
        raw_group = clean_text(body_site_group)
        text = f"{raw_site} {raw_group}".lower()

        if has_any_label_term(text, ("naris", "tonsillar", "tonsil")):
            return "Upper airway / nasal", "upper-airway/nasal body-site group"

        if has_any_label_term(text, ("bronchus", "bronchial", "trachea", "tracheal")):
            return "Bronchial / tracheal", "bronchial/tracheal body-site group"

        if has_any_label_term(text, ("stomach", "gastric", "esophagus", "oesophagus")):
            return "Stomach / esophagus", "stomach/esophagus body-site group"

        if has_any_label_term(
            text,
            (
                "duodenum",
                "duodenal",
                "jejunum",
                "ileal",
                "ileum",
                "small intestine",
                "colorectal mucosa",
                "mesenteric lymph node",
                "digestive tract",
                "rumen",
                "proventriculus",
                "cloaca",
            ),
        ):
            return "Intestinal / GI mucosa", "intestinal/gi-mucosa body-site group"

        if has_any_label_term(text, ("conjunct", "tear", "eyelid", "meibum", "ocular")):
            return "Ocular / eye", "ocular/eye body-site group"

        if has_any_label_term(text, ("breast", "milk", "colostrum", "mammary")):
            return "Breast / milk", "breast/milk body-site group"

        if has_any_label_term(text, ("meconium", "placenta", "marsupium")):
            return "Birth-associated samples", "birth-associated body-site group"

        if has_any_label_term(text, ("bile", "biliary", "liver", "pancreas", "pancreatic")):
            return "Hepatobiliary / pancreas", "hepatobiliary/pancreas body-site group"

        if has_any_label_term(
            text,
            ("cheek", "gingiva", "gingival", "dental", "dentine", "tooth", "teeth"),
        ):
            return "Oral cavity / dental", "oral/dental body-site group"

        if has_any_label_term(text, ("semen", "ovary", "oviduct", "prostate", "genitalia")):
            return "Reproductive / genital", "reproductive/genital body-site group"

        if has_any_label_term(text, ("kidney", "renal", "ureter", "urethra")):
            return "Urinary tract", "urinary body-site group"

        if has_any_label_term(text, ("ear", "otic")):
            return "Ear", "ear body-site group"

        if has_any_label_term(
            text,
            (
                "body proper",
                "anatomical surface",
                "integument",
                "entire surface",
                "forehead",
                "exoskeleton",
                "insect",
                "tadpole",
                "instar",
                "pupa",
            ),
        ):
            return "External surface / non-human", "external-surface/non-human body-site group"

        if has_any_label_term(text, ("mucosa", "mucus")):
            return "Unlocalized mucosal samples", "unlocalized-mucosal body-site group"

        if has_any_label_term(text, ("thyroid",)):
            return "Thyroid tissue", "thyroid-tissue body-site group"

        if has_any_label_term(text, ("nerve root",)):
            return "Neural tissue", "neural-tissue body-site group"

        return "Miscellaneous sampled tissues", "miscellaneous sampled tissue body-site group"
    return label, rule


def add_v3_display_rollups(
    scored: pd.DataFrame,
    config: MapConfig,
    *,
    display_mode: str = V3_GROUPED_MODE,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Attach v3 labels and return source rows plus audit rollup tables."""
    source = scored.copy()
    source["raw_condition"] = source["condition"]
    if display_mode == V3_GROUPED_MODE:
        condition_pairs = source["condition"].map(
            lambda value: v3_condition_display_label(value, config)
        )
        source["condition"] = [label for label, _rule in condition_pairs]
        source["condition_rollup_rule"] = [rule for _label, rule in condition_pairs]
    elif display_mode == V3_INDIVIDUAL_MODE:
        source["condition"] = source["condition"].map(clean_text)
        source["condition_rollup_rule"] = "kept raw individual condition label"
    else:
        msg = f"Unknown v3 display mode: {display_mode}"
        raise ValueError(msg)

    source["raw_body_site"] = source["body_site"]
    source["raw_body_site_group"] = source["body_site_group"]
    if display_mode == V3_GROUPED_MODE:
        body_pairs = [
            v3_body_site_display_label(row.body_site, row.body_site_group, config)
            for row in source.itertuples()
        ]
        source["body_site_group"] = [label for label, _rule in body_pairs]
        source["body_site_rollup_rule"] = [rule for _label, rule in body_pairs]
    else:
        source["body_site_group"] = source["body_site"].map(
            lambda value: (
                "Unknown body site" if clean_text(value) == "unknown" else clean_text(value)
            )
        )
        source["body_site_rollup_rule"] = "kept raw individual body-site label"
    source = source[source["condition"] != "unknown"].copy()

    condition_rollups = (
        source.groupby(
            ["raw_condition", "condition", "condition_rollup_rule"],
            dropna=False,
            sort=False,
        )
        .agg(
            profile_count=("profile_id", "count"),
            study_count=("study_id", "nunique"),
            focus_profile_count=(
                "inc_dec",
                lambda values: int((values == config.focus_direction).sum()),
            ),
            selected_guild_enriched_profile_count=(
                "has_guild_enriched_trait",
                lambda values: int(values.sum()),
            ),
        )
        .reset_index()
        .sort_values(
            ["profile_count", "study_count", "raw_condition"],
            ascending=[False, False, True],
        )
    )
    body_site_rollups = (
        source.groupby(
            [
                "raw_body_site",
                "raw_body_site_group",
                "body_site_group",
                "body_site_rollup_rule",
            ],
            dropna=False,
            sort=False,
        )
        .agg(
            profile_count=("profile_id", "count"),
            study_count=("study_id", "nunique"),
            focus_profile_count=(
                "inc_dec",
                lambda values: int((values == config.focus_direction).sum()),
            ),
            selected_guild_enriched_profile_count=(
                "has_guild_enriched_trait",
                lambda values: int(values.sum()),
            ),
        )
        .reset_index()
        .sort_values(
            ["profile_count", "study_count", "raw_body_site"],
            ascending=[False, False, True],
        )
    )
    return source, condition_rollups, body_site_rollups


def summarize_v3_condition_body_site_edges(scored: pd.DataFrame, config: MapConfig) -> pd.DataFrame:
    """Build grouped v3 condition/body-site edges after display-label rollups."""
    return summarize_v3_condition_body_site_edges_for_mode(
        scored,
        config,
        display_mode=V3_GROUPED_MODE,
    )


def summarize_v3_condition_body_site_edges_for_mode(
    scored: pd.DataFrame,
    config: MapConfig,
    *,
    display_mode: str,
) -> pd.DataFrame:
    """Build v3 condition/body-site edges for one label display mode."""
    source, _, _ = add_v3_display_rollups(scored, config, display_mode=display_mode)
    records: list[dict[str, object]] = []
    for (condition, body_site_group), group in source.groupby(
        ["condition", "body_site_group"],
        sort=False,
        dropna=False,
    ):
        focus = group[group["inc_dec"] == config.focus_direction]
        comparator = group[group["inc_dec"] == config.comparator_direction]
        focus_enriched = focus[focus["has_guild_enriched_trait"]]
        background_no_enriched = group[~group["has_guild_enriched_trait"]]
        focus_no_enriched = focus[~focus["has_guild_enriched_trait"]]
        comparator_enriched = comparator[comparator["has_guild_enriched_trait"]]

        focus_mean = nan_mean(focus["guild_z_mean"])
        focus_enriched_mean = nan_mean(focus_enriched["guild_z_mean"])
        comparator_mean = nan_mean(comparator["guild_z_mean"])
        signal_focus_mean = focus_enriched_mean if np.isfinite(focus_enriched_mean) else focus_mean
        contrast = (
            signal_focus_mean - comparator_mean
            if np.isfinite(signal_focus_mean) and np.isfinite(comparator_mean)
            else signal_focus_mean
        )
        records.append(
            {
                "analysis": config.key,
                "condition": clean_text(condition),
                "body_site_group": clean_text(body_site_group),
                "raw_conditions": join_unique(group["raw_condition"]),
                "raw_body_sites": join_unique(group["raw_body_site"]),
                "raw_body_site_groups": join_unique(group["raw_body_site_group"]),
                "profile_count": len(group),
                "study_count": group["study_id"].nunique(),
                "study_ids": join_unique(group["study_id"]),
                "focus_direction": config.focus_label,
                "comparator_direction": config.comparator_label,
                "focus_profile_count": len(focus),
                "focus_study_count": focus["study_id"].nunique(),
                "focus_study_ids": join_unique(focus["study_id"]),
                "focus_enriched_profile_count": len(focus_enriched),
                "focus_enriched_study_count": focus_enriched["study_id"].nunique(),
                "focus_enriched_study_ids": join_unique(focus_enriched["study_id"]),
                "focus_no_enriched_trait_profile_count": len(focus_no_enriched),
                "comparator_profile_count": len(comparator),
                "comparator_study_count": comparator["study_id"].nunique(),
                "comparator_study_ids": join_unique(comparator["study_id"]),
                "comparator_enriched_profile_count": len(comparator_enriched),
                "background_no_enriched_profile_count": len(background_no_enriched),
                "background_no_enriched_study_count": background_no_enriched["study_id"].nunique(),
                "background_no_enriched_study_ids": join_unique(background_no_enriched["study_id"]),
                "increased_profile_count": int((group["inc_dec"] == "inc").sum()),
                "decreased_profile_count": int((group["inc_dec"] == "dec").sum()),
                "selected_guild_enriched_profile_count": int(
                    group["has_guild_enriched_trait"].sum()
                ),
                "mean_guild_z": nan_mean(group["guild_z_mean"]),
                "focus_mean_guild_z": focus_mean,
                "focus_enriched_mean_guild_z": focus_enriched_mean,
                "comparator_mean_guild_z": comparator_mean,
                "focus_enriched_minus_comparator_guild_z": contrast,
            }
        )

    if not records:
        return pd.DataFrame()

    edges = pd.DataFrame.from_records(records)
    edges = edges[edges["focus_enriched_profile_count"] > 0].copy()
    if edges.empty:
        return edges
    edges["rank_signal"] = edges["focus_enriched_minus_comparator_guild_z"].fillna(
        edges["focus_enriched_mean_guild_z"]
    )
    edges["rank_score"] = (
        edges["focus_enriched_profile_count"]
        * np.log1p(edges["focus_enriched_study_count"].clip(lower=1))
        * (1.0 + np.maximum(edges["rank_signal"], 0.0))
    )
    edges["displayed_in_map"] = False
    return edges.sort_values(
        [
            "focus_enriched_profile_count",
            "rank_score",
            "rank_signal",
            "study_count",
        ],
        ascending=[False, False, False, False],
    )


def select_v3_display_edges(edges: pd.DataFrame) -> pd.DataFrame:
    """Choose grouped v3 edges, preserving body-site diversity."""
    if edges.empty:
        return edges.copy()

    candidates = edges[
        (edges["focus_enriched_profile_count"] >= V3_MIN_FOCUS_ENRICHED_PROFILES)
        & (edges["focus_enriched_study_count"] >= MIN_STUDIES)
    ].copy()
    if candidates.empty:
        candidates = edges.head(V3_MAX_DISPLAY_EDGES).copy()

    candidates = candidates.sort_values(
        [
            "focus_enriched_profile_count",
            "rank_score",
            "rank_signal",
            "study_count",
        ],
        ascending=[False, False, False, False],
    )
    per_body_site = candidates.groupby("body_site_group", sort=False).head(2)
    selected_indices: list[int] = []
    selected_conditions: set[str] = set()

    def add_candidate(index: int, row: pd.Series) -> None:
        if index in selected_indices or len(selected_indices) >= V3_MAX_DISPLAY_EDGES:
            return
        condition = str(row["condition"])
        if (
            condition not in selected_conditions
            and len(selected_conditions) >= V3_MAX_DISPLAY_CONDITIONS
        ):
            return
        selected_indices.append(index)
        selected_conditions.add(condition)

    for index, row in per_body_site.iterrows():
        add_candidate(index, row)
    for index, row in candidates.iterrows():
        add_candidate(index, row)

    return candidates.loc[selected_indices].copy()


def scaled_line_width(values: pd.Series) -> pd.Series:
    """Scale profile counts to readable line widths."""
    counts = values.astype(float)
    if counts.max() == counts.min():
        return pd.Series(4.0, index=values.index)
    scaled = (np.sqrt(counts) - np.sqrt(counts.min())) / (
        np.sqrt(counts.max()) - np.sqrt(counts.min())
    )
    return 1.2 + scaled * 8.0


def scaled_line_width_from_range(
    values: pd.Series,
    *,
    min_count: float,
    max_count: float,
    min_width: float = 0.9,
    max_width: float = 9.5,
) -> pd.Series:
    """Scale profile counts against a shared count range."""
    counts = values.astype(float)
    widths = pd.Series(0.0, index=values.index, dtype=float)
    positive = counts > 0
    if not positive.any():
        return widths
    if max_count <= min_count:
        widths.loc[positive] = (min_width + max_width) / 2
        return widths
    scaled = (np.sqrt(counts[positive]) - np.sqrt(min_count)) / (
        np.sqrt(max_count) - np.sqrt(min_count)
    )
    widths.loc[positive] = min_width + scaled.clip(0.0, 1.0) * (max_width - min_width)
    return widths


def node_positions(
    labels: list[str],
    *,
    top: float = 0.88,
    bottom: float = 0.16,
) -> dict[str, float]:
    """Return evenly spaced vertical positions for node labels."""
    if len(labels) == 1:
        return {labels[0]: 0.52}
    positions = np.linspace(top, bottom, len(labels))
    return dict(zip(labels, positions, strict=True))


def union_joined_counts(group: pd.DataFrame, column: str, fallback_column: str) -> int:
    """Count unique joined identifiers, falling back to the existing count column."""
    values: set[str] = set()
    if column in group:
        for value in group[column]:
            values.update(split_joined_values(value))
    if values:
        return len(values)
    if fallback_column in group and not group[fallback_column].empty:
        return int(group[fallback_column].max())
    return 0


def weighted_node_stats(edges: pd.DataFrame, node_column: str) -> pd.DataFrame:
    """Summarize displayed node support and focus mean z-score."""
    rows: list[dict[str, object]] = []
    for label, group in edges.groupby(node_column, sort=False):
        weights = group["focus_profile_count"].to_numpy(dtype=float)
        focus_values = group["focus_mean_guild_z"].to_numpy(dtype=float)
        contrast_values = (
            group["focus_minus_comparator_guild_z"]
            .fillna(group["focus_mean_guild_z"])
            .to_numpy(dtype=float)
        )
        rows.append(
            {
                node_column: label,
                "focus_profile_count": int(group["focus_profile_count"].sum()),
                "study_count": union_joined_counts(group, "study_ids", "study_count"),
                "focus_mean_guild_z": float(np.average(focus_values, weights=weights)),
                "focus_minus_comparator_guild_z": float(
                    np.average(contrast_values, weights=weights)
                ),
            }
        )
    return pd.DataFrame.from_records(rows).set_index(node_column)


def weighted_v3_node_stats(edges: pd.DataFrame, node_column: str) -> pd.DataFrame:
    """Summarize v3 node support over displayed grouped edges."""
    rows: list[dict[str, object]] = []
    for label, group in edges.groupby(node_column, sort=False):
        enriched_count = int(group["focus_enriched_profile_count"].sum())
        background_count = int(group["background_no_enriched_profile_count"].sum())
        total_count = int(group["profile_count"].sum())
        weights = group["focus_enriched_profile_count"].to_numpy(dtype=float)
        values = (
            group["focus_enriched_minus_comparator_guild_z"]
            .fillna(group["focus_enriched_mean_guild_z"])
            .fillna(0.0)
            .to_numpy(dtype=float)
        )
        rows.append(
            {
                node_column: label,
                "focus_enriched_profile_count": enriched_count,
                "background_no_enriched_profile_count": background_count,
                "profile_count": total_count,
                "study_count": union_joined_counts(group, "study_ids", "study_count"),
                "focus_enriched_study_count": union_joined_counts(
                    group,
                    "focus_enriched_study_ids",
                    "focus_enriched_study_count",
                ),
                "focus_enriched_minus_comparator_guild_z": float(
                    np.average(values, weights=weights) if weights.sum() else np.nan
                ),
            }
        )
    return pd.DataFrame.from_records(rows).set_index(node_column)


def draw_curve(
    ax: plt.Axes,
    *,
    start: tuple[float, float],
    end: tuple[float, float],
    color: str,
    linewidth: float,
    alpha: float = 0.72,
    zorder: float = 1.0,
    linestyle: str | tuple[float, tuple[float, ...]] = "-",
    capstyle: str = "round",
    control_y_offset: float = 0.0,
    control_x_offset: float = 0.0,
) -> None:
    """Draw one left-to-right Bezier edge."""
    x0, y0 = start
    x1, y1 = end
    control_dx = (x1 - x0) * 0.46
    path = MplPath(
        [
            (x0, y0),
            (x0 + control_dx + control_x_offset, y0 + control_y_offset),
            (x1 - control_dx + control_x_offset, y1 + control_y_offset),
            (x1, y1),
        ],
        [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4],
    )
    patch = PathPatch(
        path,
        facecolor="none",
        edgecolor=color,
        linewidth=linewidth,
        alpha=alpha,
        linestyle=linestyle,
        capstyle=capstyle,
        joinstyle="round",
        zorder=zorder,
    )
    ax.add_patch(patch)


def draw_node(
    ax: plt.Axes,
    *,
    x: float,
    y: float,
    label: str,
    stats_row: pd.Series,
    align: str,
    note: str = "",
) -> None:
    """Draw a node marker and its text label."""
    ax.add_patch(
        Circle(
            (x, y),
            0.012,
            facecolor="white",
            edgecolor=SLATE_700,
            linewidth=1.0,
            zorder=5,
        )
    )
    text_x = x - 0.025 if align == "right" else x + 0.025
    text_lines = [
        wrap_label(label),
        (
            f"n={int(stats_row.focus_profile_count)}, "
            f"s={int(stats_row.study_count)}, "
            f"z={stats_row.focus_mean_guild_z:+.2f}"
        ),
    ]
    if note:
        text_lines.extend(wrap(note, width=42))
    ax.text(
        text_x,
        y,
        "\n".join(text_lines),
        ha=align,
        va="center",
        fontsize=7.7 if note else 8.6,
        color=SLATE_900,
        linespacing=1.08,
    )


def draw_node_v3(
    ax: plt.Axes,
    *,
    x: float,
    y: float,
    label: str,
    stats_row: pd.Series,
    align: str,
) -> None:
    """Draw a v3 node marker with enriched and background counts."""
    ax.add_patch(
        Circle(
            (x, y),
            0.012,
            facecolor="white",
            edgecolor=SLATE_700,
            linewidth=1.0,
            zorder=5,
        )
    )
    text_x = x - 0.025 if align == "right" else x + 0.025
    label_lines = [
        wrap_label(label),
        (
            f"n={int(stats_row.focus_enriched_profile_count)}, "
            f"s={int(stats_row.focus_enriched_study_count)}, "
            f"bg={int(stats_row.background_no_enriched_profile_count)}"
        ),
    ]
    ax.text(
        text_x,
        y,
        "\n".join(label_lines),
        ha=align,
        va="center",
        fontsize=8.0,
        color=SLATE_900,
        linespacing=1.08,
    )


def build_condition_sampled_site_notes(
    display_edges: pd.DataFrame,
    sampled_sites: pd.DataFrame,
) -> dict[str, str]:
    """List sampled body sites for displayed conditions that are not drawn as edges."""
    notes: dict[str, str] = {}
    if sampled_sites.empty:
        return notes

    displayed_by_condition = {
        condition: set(group["body_site_group"])
        for condition, group in display_edges.groupby("condition", sort=False)
    }
    for condition, displayed_sites in displayed_by_condition.items():
        condition_sites = sampled_sites[
            (sampled_sites["condition"] == condition)
            & ~sampled_sites["body_site_group"].isin(displayed_sites)
        ].copy()
        if condition_sites.empty:
            continue
        condition_sites = condition_sites.sort_values(
            ["profile_count", "study_count", "body_site_group"],
            ascending=[False, False, True],
        )
        entries = [
            (
                f"{compact_label(row.body_site_group, width=18)} "
                f"profiles={int(row.profile_count)} studies={int(row.study_count)}"
            )
            for row in condition_sites.head(MAX_OTHER_SITES_PER_CONDITION).itertuples()
        ]
        remaining = len(condition_sites) - len(entries)
        if remaining > 0:
            entries.append(f"+{remaining} more")
        notes[condition] = "other sampled sites: " + "; ".join(entries)
    return notes


def plot_condition_body_site_map(
    display_edges: pd.DataFrame,
    sampled_sites: pd.DataFrame,
    config: MapConfig,
    output_path: Path,
    *,
    title_suffix: str = "",
    selection_note: str | None = None,
    dense_layout: bool = False,
) -> None:
    """Render a condition-to-body-site bipartite guild map."""
    if display_edges.empty:
        msg = f"No display edges available for {config.key}."
        raise ValueError(msg)

    condition_order = (
        display_edges.groupby("condition", sort=False)["rank_score"]
        .max()
        .sort_values(ascending=False)
        .index.tolist()
    )
    body_site_order = (
        display_edges.groupby("body_site_group", sort=False)["focus_profile_count"]
        .sum()
        .sort_values(ascending=False)
        .index.tolist()
    )
    node_top = 0.98 if dense_layout else 0.88
    node_bottom = 0.05 if dense_layout else 0.16
    header_y = 0.985 if dense_layout else 0.93
    condition_y = node_positions(condition_order, top=node_top, bottom=node_bottom)
    body_site_y = node_positions(body_site_order, top=node_top, bottom=node_bottom)
    condition_stats = weighted_node_stats(display_edges, "condition")
    body_site_stats = weighted_node_stats(display_edges, "body_site_group")
    condition_notes = build_condition_sampled_site_notes(display_edges, sampled_sites)

    height = max(
        8.0,
        1.2 + 0.64 * max(len(condition_order), len(body_site_order)) + 0.12 * len(condition_notes),
    )
    fig, ax = plt.subplots(figsize=(13.5, height))
    if dense_layout:
        fig.subplots_adjust(left=0.05, right=0.90, top=0.91, bottom=0.07)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    values = display_edges["focus_minus_comparator_guild_z"].fillna(
        display_edges["focus_mean_guild_z"]
    )
    vmax = max(0.4, float(values.abs().max()))
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    cmap = plt.get_cmap("RdBu_r")
    widths = scaled_line_width(display_edges["focus_profile_count"])

    left_x = 0.34
    right_x = 0.68
    for edge, linewidth in zip(display_edges.itertuples(), widths, strict=True):
        signal = (
            edge.focus_minus_comparator_guild_z
            if np.isfinite(edge.focus_minus_comparator_guild_z)
            else edge.focus_mean_guild_z
        )
        draw_curve(
            ax,
            start=(left_x, condition_y[edge.condition]),
            end=(right_x, body_site_y[edge.body_site_group]),
            color=cmap(norm(signal)),
            linewidth=float(linewidth),
        )

    for condition in condition_order:
        draw_node(
            ax,
            x=left_x,
            y=condition_y[condition],
            label=condition,
            stats_row=condition_stats.loc[condition],
            align="right",
            note=condition_notes.get(condition, ""),
        )
    for body_site in body_site_order:
        draw_node(
            ax,
            x=right_x,
            y=body_site_y[body_site],
            label=body_site,
            stats_row=body_site_stats.loc[body_site],
            align="left",
        )

    fig.text(
        0.05,
        0.965,
        f"{config.title}{title_suffix}",
        ha="left",
        va="top",
        fontsize=17,
        fontweight="bold",
    )
    fig.text(
        0.05,
        0.925,
        wrapped_plot_text(
            "Edges connect BugSigDB conditions to anatomical body-site groups. "
            f"Edge width = n ({config.focus_label} profile count); color = "
            f"{config.focus_label} minus {config.comparator_label} mean trait "
            "z-score; node key: n = focus-direction profiles, s = unique "
            "studies, z = focus mean trait z-score."
        ),
        ha="left",
        va="top",
        fontsize=10,
        color=SLATE_700,
    )
    ax.text(0.21, header_y, "Conditions", ha="center", va="center", fontsize=12, color=SLATE_700)
    ax.text(
        0.79,
        header_y,
        "Body-site groups",
        ha="center",
        va="center",
        fontsize=12,
        color=SLATE_700,
    )

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.030, pad=0.02, shrink=0.55)
    cbar.set_label(
        f"{config.focus_label} - {config.comparator_label} mean trait z-score",
        fontsize=9,
    )
    cbar.ax.tick_params(labelsize=8)

    representative_counts = sorted(
        {
            int(display_edges["focus_profile_count"].min()),
            int(display_edges["focus_profile_count"].median()),
            int(display_edges["focus_profile_count"].max()),
        }
    )
    legend_lines = [
        Line2D(
            [0],
            [0],
            color=EDGE_NEUTRAL,
            linewidth=float(scaled_line_width(pd.Series([count, max(representative_counts)]))[0]),
            alpha=0.72,
            solid_capstyle="round",
            label=f"n={count} profiles",
        )
        for count in representative_counts
    ]
    legend = ax.legend(
        handles=legend_lines,
        title=f"Edge width: {config.focus_label} profiles (n)",
        loc="lower center",
        bbox_to_anchor=(0.50, -0.015),
        ncol=len(legend_lines),
        frameon=False,
        fontsize=8,
        title_fontsize=9,
    )
    ax.add_artist(legend)

    if selection_note is None:
        selection_note = (
            f"Displayed edges require at least {MIN_FOCUS_PROFILES} {config.focus_label} profiles, "
            f"{MIN_STUDIES} studies, positive focus mean z-score, and are ranked by "
            "contrast-weighted support. Condition notes list sampled body sites for that "
            "condition that were not drawn as enriched edges."
        )

    fig.text(
        0.05,
        0.035,
        selection_note,
        ha="left",
        va="bottom",
        fontsize=8.2,
        color=SLATE_500,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_condition_body_site_map_v3(
    display_edges: pd.DataFrame,
    config: MapConfig,
    output_path: Path,
    *,
    display_mode: str = V3_GROUPED_MODE,
) -> None:
    """Render a v3 condition-to-body-site guild map."""
    if display_edges.empty:
        msg = f"No v3 {display_mode} display edges available for {config.key}."
        raise ValueError(msg)
    if display_mode == V3_INDIVIDUAL_MODE:
        condition_header = "Individual condition labels"
        body_site_header = "Individual body sites"
        pair_description = "same exact condition/body-site pair"
        footer_note = (
            "V3 individual display keeps raw condition labels and raw body-site labels, "
            f"requires at least {V3_MIN_FOCUS_ENRICHED_PROFILES} colored profiles and "
            f"{MIN_STUDIES} colored-profile studies where possible, and preserves all "
            "individual edge details in the companion TSV."
        )
    elif display_mode == V3_GROUPED_MODE:
        condition_header = "Condition display groups"
        body_site_header = "Body-site display groups"
        pair_description = "same display condition/body-site pair"
        if config.key in {FACULTATIVE_CONFIG.key, ANAEROBE_CONFIG.key}:
            footer_note = (
                "V3 grouped display assigns every condition and body site to display "
                "groups, "
                f"requires at least {V3_MIN_FOCUS_ENRICHED_PROFILES} colored profiles "
                f"and {MIN_STUDIES} colored-profile studies where possible, and "
                "preserves all grouped edge and rollup details in companion TSVs."
            )
        else:
            footer_note = (
                "V3 grouped display uses disease/body-site rollups, keeps "
                "transplant/graft and metadata/exposure labels raw, "
                f"requires at least {V3_MIN_FOCUS_ENRICHED_PROFILES} colored profiles "
                f"and {MIN_STUDIES} colored-profile studies where possible, and "
                "preserves all grouped edge and rollup details in companion TSVs."
            )
    else:
        msg = f"Unknown v3 display mode: {display_mode}"
        raise ValueError(msg)

    condition_order_frame = (
        display_edges.groupby("condition", sort=False)
        .agg(
            focus_enriched_profile_count=("focus_enriched_profile_count", "sum"),
            rank_score=("rank_score", "max"),
        )
        .sort_values(
            ["focus_enriched_profile_count", "rank_score"],
            ascending=[False, False],
        )
    )
    condition_order = condition_order_frame.index.tolist()
    body_site_order_frame = (
        display_edges.groupby("body_site_group", sort=False)
        .agg(
            focus_enriched_profile_count=("focus_enriched_profile_count", "sum"),
            background_no_enriched_profile_count=(
                "background_no_enriched_profile_count",
                "sum",
            ),
        )
        .sort_values(
            ["focus_enriched_profile_count", "background_no_enriched_profile_count"],
            ascending=[False, False],
        )
    )
    body_site_order = body_site_order_frame.index.tolist()

    node_top = 0.90
    node_bottom = 0.06
    header_y = 0.965
    condition_y = node_positions(condition_order, top=node_top, bottom=node_bottom)
    body_site_y = node_positions(body_site_order, top=node_top, bottom=node_bottom)
    condition_stats = weighted_v3_node_stats(display_edges, "condition")
    body_site_stats = weighted_v3_node_stats(display_edges, "body_site_group")

    height = max(9.0, 1.5 + 0.66 * max(len(condition_order), len(body_site_order)))
    fig, ax = plt.subplots(figsize=(15.0, height))
    fig.subplots_adjust(left=0.05, right=0.90, top=0.90, bottom=0.08)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    color_values = display_edges["focus_enriched_minus_comparator_guild_z"].fillna(
        display_edges["focus_enriched_mean_guild_z"]
    )
    vmax = max(0.4, float(color_values.abs().max()))
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    cmap = plt.get_cmap("RdBu_r")

    all_counts = pd.concat(
        [
            display_edges["focus_enriched_profile_count"],
            display_edges["background_no_enriched_profile_count"],
        ],
        ignore_index=True,
    )
    positive_counts = all_counts[all_counts > 0]
    min_count = float(positive_counts.min())
    max_count = float(positive_counts.max())
    enriched_widths = scaled_line_width_from_range(
        display_edges["focus_enriched_profile_count"],
        min_count=min_count,
        max_count=max_count,
        min_width=1.1,
        max_width=9.5,
    )
    background_widths = scaled_line_width_from_range(
        display_edges["background_no_enriched_profile_count"],
        min_count=min_count,
        max_count=max_count,
        min_width=1.6,
        max_width=12.0,
    )

    left_x = 0.34
    right_x = 0.68
    for edge, linewidth in zip(display_edges.itertuples(), background_widths, strict=True):
        if linewidth <= 0:
            continue
        draw_curve(
            ax,
            start=(left_x, condition_y[edge.condition] + V3_BACKGROUND_EDGE_Y_OFFSET),
            end=(right_x, body_site_y[edge.body_site_group] + V3_BACKGROUND_EDGE_Y_OFFSET),
            color=SLATE_500,
            linewidth=float(linewidth),
            alpha=0.48,
            zorder=1.0,
        )

    for edge, linewidth in zip(display_edges.itertuples(), enriched_widths, strict=True):
        signal = (
            edge.focus_enriched_minus_comparator_guild_z
            if np.isfinite(edge.focus_enriched_minus_comparator_guild_z)
            else edge.focus_enriched_mean_guild_z
        )
        draw_curve(
            ax,
            start=(left_x, condition_y[edge.condition]),
            end=(right_x, body_site_y[edge.body_site_group]),
            color=cmap(norm(signal)),
            linewidth=float(linewidth),
            alpha=0.76,
            zorder=1.5,
        )

    for condition in condition_order:
        draw_node_v3(
            ax,
            x=left_x,
            y=condition_y[condition],
            label=condition,
            stats_row=condition_stats.loc[condition],
            align="right",
        )
    for body_site in body_site_order:
        draw_node_v3(
            ax,
            x=right_x,
            y=body_site_y[body_site],
            label=body_site,
            stats_row=body_site_stats.loc[body_site],
            align="left",
        )

    fig.text(
        0.05,
        0.965,
        config.title,
        ha="left",
        va="top",
        fontsize=17,
        fontweight="bold",
    )
    fig.text(
        0.05,
        0.925,
        wrapped_plot_text(
            f"Colored edges = {config.focus_label} profiles with at least one "
            "selected guild trait above zero. Parallel grey edges = profiles in "
            f"the {pair_description} with zero selected guild traits. Conditions "
            "are ordered by colored n. Node key: n = colored positive-trait "
            "profiles, s = unique studies for n, bg = background profiles with "
            "zero selected guild traits."
        ),
        ha="left",
        va="top",
        fontsize=10,
        color=SLATE_700,
    )
    ax.text(
        0.21,
        header_y,
        condition_header,
        ha="center",
        va="center",
        fontsize=12,
        color=SLATE_700,
    )
    ax.text(
        0.79,
        header_y,
        body_site_header,
        ha="center",
        va="center",
        fontsize=12,
        color=SLATE_700,
    )

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.030, pad=0.02, shrink=0.55)
    cbar.set_label(
        f"{config.focus_label} enriched - {config.comparator_label} mean trait z-score",
        fontsize=9,
    )
    cbar.ax.tick_params(labelsize=8)

    representative_counts = sorted(
        {
            int(display_edges["focus_enriched_profile_count"].min()),
            int(display_edges["focus_enriched_profile_count"].median()),
            int(display_edges["focus_enriched_profile_count"].max()),
        }
    )
    width_lookup = scaled_line_width_from_range(
        pd.Series(representative_counts, index=representative_counts),
        min_count=min_count,
        max_count=max_count,
        min_width=1.1,
        max_width=9.5,
    )
    legend_lines = [
        Line2D(
            [0],
            [0],
            color=EDGE_NEUTRAL,
            linewidth=float(width_lookup.loc[count]),
            alpha=0.76,
            solid_capstyle="round",
            label=f"n={count} profiles",
        )
        for count in representative_counts
    ]
    legend_lines.append(
        Line2D(
            [0],
            [0],
            color=SLATE_500,
            linewidth=6.5,
            alpha=0.48,
            solid_capstyle="round",
            label="bg profiles: no selected guild traits",
        )
    )
    legend = ax.legend(
        handles=legend_lines,
        title=f"Edge width: {config.focus_label} positive-trait profiles (n)",
        loc="lower center",
        bbox_to_anchor=(0.50, -0.02),
        ncol=min(len(legend_lines), 4),
        frameon=False,
        fontsize=8,
        title_fontsize=9,
    )
    ax.add_artist(legend)

    fig.text(
        0.05,
        0.035,
        footer_note,
        ha="left",
        va="bottom",
        fontsize=8.2,
        color=SLATE_500,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def add_displayed_flag(
    edges: pd.DataFrame,
    display_edges: pd.DataFrame,
) -> pd.DataFrame:
    """Return an edge table with a boolean flag for the plotted subset."""
    display_keys = set(
        zip(display_edges["condition"], display_edges["body_site_group"], strict=True)
    )
    output = edges.copy()
    output["displayed_in_map"] = [
        (condition, body_site) in display_keys
        for condition, body_site in zip(output["condition"], output["body_site_group"], strict=True)
    ]
    return output


def remove_legacy_single_v3_outputs(config: MapConfig) -> None:
    """Remove superseded single v3 files after writing split v3 outputs."""
    for path in [
        config.edge_v3_tsv_path,
        config.condition_rollup_v3_tsv_path,
        config.body_site_rollup_v3_tsv_path,
        config.png_v3_path,
    ]:
        if path.exists():
            path.unlink()


def write_map_outputs(matrix: pd.DataFrame, metadata: pd.DataFrame, config: MapConfig) -> None:
    """Write the edge table and PNG for one analysis."""
    scored = add_guild_scores(matrix, metadata, config)
    sampled_sites = summarize_sampled_body_sites(scored, config)
    edges = summarize_condition_body_site_edges(sampled_sites)
    display_edges = select_display_edges(edges)
    display_edges_v2 = select_uncapped_display_edges(edges)
    _v3_source, condition_rollups_v3, body_site_rollups_v3 = add_v3_display_rollups(
        scored,
        config,
        display_mode=V3_GROUPED_MODE,
    )
    edges_v3 = summarize_v3_condition_body_site_edges(scored, config)
    display_edges_v3 = select_v3_display_edges(edges_v3)
    edges_v3_individual = summarize_v3_condition_body_site_edges_for_mode(
        scored,
        config,
        display_mode=V3_INDIVIDUAL_MODE,
    )
    display_edges_v3_individual = select_v3_display_edges(edges_v3_individual)
    display_keys = set(
        zip(display_edges["condition"], display_edges["body_site_group"], strict=True)
    )
    display_v2_keys = set(
        zip(display_edges_v2["condition"], display_edges_v2["body_site_group"], strict=True)
    )

    edge_output = edges.copy()
    edge_output["displayed_in_map"] = [
        (condition, body_site) in display_keys
        for condition, body_site in zip(
            edge_output["condition"], edge_output["body_site_group"], strict=True
        )
    ]
    sampled_output = sampled_sites.copy()
    sampled_output["displayed_in_map"] = [
        (condition, body_site) in display_keys
        for condition, body_site in zip(
            sampled_output["condition"], sampled_output["body_site_group"], strict=True
        )
    ]
    displayed_conditions = set(display_edges["condition"])
    sampled_output["sampled_not_displayed_for_displayed_condition"] = [
        condition in displayed_conditions and not displayed
        for condition, displayed in zip(
            sampled_output["condition"], sampled_output["displayed_in_map"], strict=True
        )
    ]

    edge_v2_output = edges.copy()
    edge_v2_output["displayed_in_map"] = [
        (condition, body_site) in display_v2_keys
        for condition, body_site in zip(
            edge_v2_output["condition"], edge_v2_output["body_site_group"], strict=True
        )
    ]
    sampled_v2_output = sampled_sites.copy()
    sampled_v2_output["displayed_in_map"] = [
        (condition, body_site) in display_v2_keys
        for condition, body_site in zip(
            sampled_v2_output["condition"], sampled_v2_output["body_site_group"], strict=True
        )
    ]
    displayed_v2_conditions = set(display_edges_v2["condition"])
    sampled_v2_output["sampled_not_displayed_for_displayed_condition"] = [
        condition in displayed_v2_conditions and not displayed
        for condition, displayed in zip(
            sampled_v2_output["condition"], sampled_v2_output["displayed_in_map"], strict=True
        )
    ]

    edge_v3_grouped_output = add_displayed_flag(edges_v3, display_edges_v3)
    edge_v3_individual_output = add_displayed_flag(
        edges_v3_individual,
        display_edges_v3_individual,
    )

    edge_output.to_csv(config.edge_tsv_path, sep="\t", index=False)
    sampled_output.to_csv(config.sampled_site_tsv_path, sep="\t", index=False)
    edge_v2_output.to_csv(config.edge_v2_tsv_path, sep="\t", index=False)
    sampled_v2_output.to_csv(config.sampled_site_v2_tsv_path, sep="\t", index=False)
    plot_condition_body_site_map(display_edges, sampled_sites, config, config.png_path)
    plot_condition_body_site_map(
        display_edges_v2,
        sampled_sites,
        config,
        config.png_v2_path,
        title_suffix=" (v2 uncapped)",
        selection_note=(
            f"V2 uncapped display: edges require at least {MIN_FOCUS_PROFILES} "
            f"{config.focus_label} profiles, {MIN_STUDIES} studies, and positive focus "
            "mean z-score. Every edge meeting those criteria is drawn, with no edge, "
            "condition, body-site, or per-body-site display cap."
        ),
        dense_layout=True,
    )
    edge_v3_individual_output.to_csv(
        config.edge_v3_individual_tsv_path,
        sep="\t",
        index=False,
    )
    edge_v3_grouped_output.to_csv(config.edge_v3_grouped_tsv_path, sep="\t", index=False)
    condition_rollups_v3.to_csv(
        config.condition_rollup_v3_grouped_tsv_path,
        sep="\t",
        index=False,
    )
    body_site_rollups_v3.to_csv(
        config.body_site_rollup_v3_grouped_tsv_path,
        sep="\t",
        index=False,
    )
    plot_condition_body_site_map_v3(
        display_edges_v3_individual,
        config,
        config.png_v3_individual_path,
        display_mode=V3_INDIVIDUAL_MODE,
    )
    plot_condition_body_site_map_v3(
        display_edges_v3,
        config,
        config.png_v3_grouped_path,
        display_mode=V3_GROUPED_MODE,
    )
    remove_legacy_single_v3_outputs(config)
    print(f"Wrote: {config.edge_tsv_path}")
    print(f"Wrote: {config.sampled_site_tsv_path}")
    print(f"Wrote: {config.png_path}")
    print(f"Wrote: {config.edge_v2_tsv_path}")
    print(f"Wrote: {config.sampled_site_v2_tsv_path}")
    print(f"Wrote: {config.png_v2_path}")
    print(f"Wrote: {config.edge_v3_individual_tsv_path}")
    print(f"Wrote: {config.png_v3_individual_path}")
    print(f"Wrote: {config.edge_v3_grouped_tsv_path}")
    print(f"Wrote: {config.condition_rollup_v3_grouped_tsv_path}")
    print(f"Wrote: {config.body_site_rollup_v3_grouped_tsv_path}")
    print(f"Wrote: {config.png_v3_grouped_path}")


def main() -> None:
    matrix, metadata = load_study_heatmap()
    for config in [FACULTATIVE_CONFIG, ANAEROBE_CONFIG]:
        write_map_outputs(matrix, metadata, config)


if __name__ == "__main__":
    main()
