"""Rebuild the cross-disease ``cellulose_decr`` analysis from pinned inputs (#552).

Compares the taxa annotated ``METPO:2000103|GO:0030245`` (cellulose catabolic
process) that curated BugSigDB signatures report decreased in human feces in
Parkinson's disease, HIV infection and colorectal cancer, and writes every table,
the summary JSON and a provenance record to ``data/analysis/cellulose_decr/``.

The October 16 outputs had no generator (adversarial review R11). They took each
profile's condition and body site from
``trait-visualization/study_trait_profile_study_metadata.tsv``, which then held one
modal label per *study* (#555), so 1,112 of 5,829 profiles carried another
experiment's label. ``--label-source study_modal`` recomputes those labels from
``data/raw/full_dump.csv`` (#594) and reproduces the October 16 sets, exact tests and
feces direction statistics; the default, ``experiment``, takes labels from
``data/derived/experiments.tsv``, which matches ``data/raw/full_dump.csv`` for every
experiment.

Other work regenerates the inputs in place, so a run stops when an input's hash
differs from the committed provenance or the contrast manifest came from another
classifier version (#562). The pins describe a v2 contrast manifest and the metadata
table before #555 relabelled it, so until #580 re-pins them a checkout with a newer
classifier or metadata refuses to run; the committed tables reproduce at the
provenance's ``code_commit`` from the pinned inputs (#642).

The counts are annotation counts of curated signature taxa, not measurements of
patients, cells or pathway activity (methods/biological_interpretation.md).
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import itertools
import json
import math
import re
import sys
from collections import Counter, defaultdict
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt
import pandas as pd
import scipy
from scipy import optimize, special, stats

from mtp.composition_denominators import load_retained_profiles, trait_fractions
from mtp.contrast_orientation import CLASSIFIER_VERSION
from mtp.paths import REPO_ROOT, file_sha256

ANALYSIS_DIR = REPO_ROOT / "data" / "analysis" / "cellulose_decr"
PINNED_PROVENANCE = ANALYSIS_DIR / "cellulose_decr_provenance.json"

TRAIT_ID = "METPO:2000103|GO:0030245"
HALOPHILIC_ID = "METPO:1000623"
CONDITIONS = {"PD": "Parkinson's disease", "HIV": "HIV infection", "CRC": "Colorectal cancer"}
# Other BugSigDB labels for the same three diseases (#561). The primary analysis
# keeps the exact labels; the "synonyms merged" sensitivity and
# ``--condition-labels synonyms`` count these experiments under their disease.
CONDITION_SYNONYMS = {
    "PD": (),
    "HIV": ("HIV-1 infection", "Human immunodeficiency virus"),
    "CRC": (
        "Colorectal carcinoma",
        "Colorectal adenocarcinoma",
        "Metastatic colorectal cancer",
        "Colon carcinoma",
    ),
}
# Labels that name one of the diseases or the colorectum, listed in the summary with
# their status so no exclusion is silent. Adenoma, polyp, neoplasm, stromal tumour,
# mother-to-child transmission, simian virus and cell-line labels are different
# conditions and stay out.
NEAR_LABEL_PATTERN = re.compile(
    r"parkinson|\bhiv\b|immunodeficiency virus|colorectal|\bcolon\b|\brectal\b|sigmoid", re.I
)
# Curated "Feces" whose paper sampled something else (#566); excluded in a sensitivity.
SITE_MISCURATIONS = {
    "bsdb:24451087": (
        "McHardy et al. 2013: rectal mucosa collected with a sponge, curated as Feces"
    ),
}
# HIV experiments the contrast classifier does not call disease-vs-healthy although
# they compare people with and without HIV (#569); counted in sensitivities only.
HIV_RECLASSIFIED = {
    ("bsdb:33102244", "bsdb:33102244-1"): "HIV- vs HIV+ (classified disease_vs_disease)",
}
HIV_ON_ART = {
    ("bsdb:29411528", "bsdb:29411528-5"): "HAART-treated HIV+ vs controls (treatment_exposure)",
    ("bsdb:38383483", "bsdb:38383483-3"): "PLWH on ART vs PLWOH (treatment_exposure)",
    ("bsdb:39880958", "bsdb:39880958-1"): "PLWH on ART vs HIV- (treatment_exposure)",
}
# A within-patient tumour-vs-adjacent-tissue contrast the classifier calls
# disease-vs-healthy (#569); left out in a sensitivity.
CRC_WITHIN_PATIENT = {
    ("bsdb:25699023", "bsdb:25699023-2"): (
        "adjacent non-cancerous vs cancerous tissue (classified disease_vs_healthy)"
    ),
}
SHORT = list(CONDITIONS)
PAIRS = list(itertools.combinations(SHORT, 2))
FECES = "Feces"
NON_FECAL_GI = "Non-fecal GI (mucosa/lumen)"
ORAL = "Oral / Oropharyngeal"
OTHER = "Other non-GI"
COMPARTMENTS = (FECES, NON_FECAL_GI, ORAL, OTHER)
POOLED = "Pooled (PD + HIV + CRC)"
GENUS_SPECIES = ("genus", "species")

ALL_CONTRASTS = "all_contrasts"
ORIENTED = "disease_vs_healthy_oriented"
STRONG = "disease_vs_healthy_strong_evidence"
NO_MSM = "disease_vs_healthy_oriented_excl_msm_mention"
POPULATIONS = (ALL_CONTRASTS, ORIENTED, STRONG)

INPUTS = {
    "profile_metadata": "trait-visualization/study_trait_profile_study_metadata.tsv",
    "experiments": "data/derived/experiments.tsv",
    "studies": "data/derived/studies.tsv",
    "profile_counts": "data/derived/study_trait_profiles_wide.tsv",
    "signature_taxa": "data/derived/bugsigdb_taxa.tsv",
    "taxa_trait_matrix": "data/derived/taxa_trait_matrix.tsv",
    "trait_columns": "data/derived/taxa_trait_matrix_columns.tsv",
    "contrast_manifest": "data/derived/contrast_orientation/contrast_manifest.tsv",
    "full_dump": "data/raw/full_dump.csv",
    "kg_edges": "data/kg-microbe/merged-kg_edges.tsv",
}

# Body sites outside exact "Feces" are grouped by keyword. Multi-site labels such
# as "Feces,Blood" or "Oral cavity,Feces" are "Other non-GI"; for single sites,
# oral terms are checked first, then sites that are clearly not gastrointestinal,
# then GI terms, and anything else is "Other non-GI". The mapping actually applied
# to the three conditions is written to the summary.
ORAL_TERMS = (
    "saliva",
    "oral",
    "mouth",
    "tongue",
    "dental",
    "plaque",
    "oropharyn",
    "throat",
    "buccal",
    "gingiv",
    "tonsil",
    "pharyn",
)
NON_GI_TERMS = (
    "nasal",
    "nose",
    "skin",
    "vagin",
    "cervi",
    "uter",
    "semen",
    "penis",
    "blood",
    "lung",
    "sputum",
    "milk",
    "urine",
    "breast",
    "eye",
    "conjunctiv",
)
GI_TERMS = (
    "colon",
    "colorect",
    "rect",
    "duoden",
    "intestin",
    "ileum",
    "ileal",
    "jejun",
    "cecum",
    "caecum",
    "stomach",
    "gastric",
    "sigmoid",
    "mucosa",
    "bowel",
    "esophag",
    "oesophag",
    "gut",
)

MSM_PATTERN = re.compile(
    r"\bMSM\b|men who have sex with men|men having sex with men|homosexual|\bgay\b|bisexual",
    re.IGNORECASE,
)

FLOAT_FORMAT = "%.6g"
ALPHA = 0.05
BATCH = 5000

BoolArray = npt.NDArray[np.bool_]
FloatArray = npt.NDArray[np.float64]
IntArray = npt.NDArray[np.int64]


@dataclass(frozen=True)
class Corpus:
    """Profiles, their signature taxa and the annotation matrix for one label source."""

    profiles: pd.DataFrame
    taxa: pd.DataFrame
    matrix: pd.DataFrame
    trait_labels: pd.DataFrame
    manifest: pd.DataFrame
    label_source: str
    condition_labels: str = "exact"


SYNONYM_TO_CONDITION = {
    label: CONDITIONS[short] for short, labels in CONDITION_SYNONYMS.items() for label in labels
}


def merge_synonyms(conditions: pd.Series) -> pd.Series:
    """Condition labels with CONDITION_SYNONYMS replaced by their canonical label."""
    merged: pd.Series = conditions.replace(SYNONYM_TO_CONDITION)
    return merged


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0] if __doc__ else None)
    parser.add_argument(
        "--input-root",
        type=Path,
        default=REPO_ROOT,
        help="Checkout holding data/ and trait-visualization/ inputs (default: this one).",
    )
    parser.add_argument("--output-dir", type=Path, default=ANALYSIS_DIR)
    parser.add_argument(
        "--label-source",
        choices=("experiment", "study_modal"),
        default="experiment",
        help="Where each profile's condition and body site come from (see module docstring).",
    )
    parser.add_argument(
        "--condition-labels",
        choices=("exact", "synonyms"),
        default="exact",
        help="Keep the exact condition labels, or also count CONDITION_SYNONYMS (#561).",
    )
    parser.add_argument("--permutations", type=int, default=100_000)
    parser.add_argument("--seed", type=int, default=20261017)
    parser.add_argument(
        "--code-commit",
        default="unrecorded",
        help="Revision of this script, e.g. $(git rev-parse HEAD) with +dirty when modified.",
    )
    parser.add_argument(
        "--pinned-provenance",
        type=Path,
        default=PINNED_PROVENANCE,
        help="Provenance whose pinned input hashes this run must match (default: the committed "
        "one). It must exist unless --repin is given.",
    )
    pinning = parser.add_mutually_exclusive_group()
    pinning.add_argument(
        "--allow-input-drift",
        action="store_true",
        help="Run even when inputs differ from the pins; the pins are carried forward unchanged "
        "and the drifted inputs are recorded (#562).",
    )
    pinning.add_argument(
        "--repin",
        action="store_true",
        help="Pin this run's input hashes as the new reference, e.g. after a BugSigDB release.",
    )
    return parser.parse_args(argv)


# --------------------------------------------------------------------------- inputs


def read_tsv(path: Path, **kwargs: Any) -> pd.DataFrame:
    frame: pd.DataFrame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, **kwargs)
    return frame


def compartment_of(body_site: str) -> str:
    """Assign an experiment body-site label to one of the four compartments."""
    site = body_site.strip().lower()
    if site == FECES.lower():
        return FECES
    if "," in site:
        return OTHER
    if any(term in site for term in ORAL_TERMS):
        return ORAL
    if any(term in site for term in NON_GI_TERMS):
        return OTHER
    if any(term in site for term in GI_TERMS):
        return NON_FECAL_GI
    return OTHER


def clean_label(value: str | None) -> str:
    """A BugSigDB label, with blanks and NA read as 'unknown'."""
    text = (value or "").strip()
    if not text or text.upper() == "NA" or text.lower() == "nan":
        return "unknown"
    return text


def most_common_label(values: Iterable[str]) -> str:
    """The most common known label, ties broken alphabetically."""
    counts = Counter(label for label in map(clean_label, values) if label != "unknown")
    if not counts:
        return "unknown"
    return min(counts.items(), key=lambda item: (-item[1], item[0]))[0]


def study_modal_labels(dump_path: Path) -> pd.DataFrame:
    """Each study's most common condition and body site over its signatures.

    These are the labels every profile carried in study_trait_profile_study_metadata.tsv
    before #555, computed the same way (``most_common_label`` in
    scripts/create_study_trait_profile_heatmap.py). Computing them here keeps
    ``--label-source study_modal`` reproducing the October 16 outputs whatever that table
    now holds (#594).
    """
    values: dict[str, dict[str, list[str]]] = defaultdict(
        lambda: {"condition": [], "body_site": []}
    )
    with dump_path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(line for line in handle if not line.startswith("#")):
            bsdb_id = clean_label(row.get("BSDB ID"))
            if bsdb_id == "unknown":
                continue
            study = "bsdb:" + bsdb_id.removeprefix("bsdb:").split("/", 1)[0]
            values[study]["condition"].append(row.get("Condition") or "")
            values[study]["body_site"].append(row.get("Body site") or "")
    rows = [
        {
            "study_id": study,
            "condition": most_common_label(labels["condition"]),
            "body_site": most_common_label(labels["body_site"]),
        }
        for study, labels in values.items()
    ]
    return pd.DataFrame(rows).set_index("study_id")


def load_corpus(
    paths: Mapping[str, Path], label_source: str, condition_labels: str = "exact"
) -> Corpus:
    metadata = read_tsv(paths["profile_metadata"], index_col=0)
    metadata.index.name = "profile_id"
    metadata["total_trait_signal"] = pd.to_numeric(metadata["total_trait_signal"])
    keys = ["study_id", "experiment_id"]
    experiments = read_tsv(paths["experiments"])
    if label_source == "experiment":
        labelled = (
            metadata[[*keys, "inc_dec"]]
            .reset_index()
            .merge(
                experiments[[*keys, "condition", "body_site", "host_species"]],
                on=keys,
                how="left",
                validate="many_to_one",
                indicator=True,
            )
            .set_index("profile_id")
        )
        unmatched = labelled["_merge"].ne("both")
        if unmatched.any():
            raise ValueError(f"{int(unmatched.sum())} profiles have no experiments.tsv row")
        hosts = set(labelled["host_species"])
        if hosts != {"Homo sapiens"}:
            raise ValueError(f"expected human profiles only, found hosts {sorted(hosts)}")
        profiles = labelled[[*keys, "inc_dec", "condition", "body_site"]].copy()
    else:
        modal = study_modal_labels(paths["full_dump"])
        profiles = metadata[[*keys, "inc_dec"]].join(modal, on="study_id")
        missing = profiles["condition"].isna()
        if missing.any():
            raise ValueError(f"{int(missing.sum())} profiles have no study in full_dump.csv")
    if condition_labels == "synonyms":
        profiles["condition"] = merge_synonyms(profiles["condition"])
    profiles["compartment"] = profiles["body_site"].map(compartment_of)

    # Annotation share over all retained non-MediaDive columns, aligned by key, and
    # the per-trait z-score of log1p(1000 * share) over all retained profiles, as in
    # scripts/create_pd_cellulose_analysis.py.
    counts = load_retained_profiles(paths["profile_counts"], metadata)
    share = trait_fractions(counts)[TRAIT_ID].astype(float)
    scaled = (share * 1000).map(math.log1p)
    profiles["frac"] = share
    profiles["z"] = (scaled - scaled.mean()) / scaled.std(ddof=0)

    matrix = pd.read_csv(paths["taxa_trait_matrix"], sep="\t", index_col=0)
    taxa = read_tsv(paths["signature_taxa"])
    direction = taxa["direction"].map({"increased": "inc", "decreased": "dec"})
    if direction.isna().any():
        raise ValueError(
            f"unexpected directions {sorted(set(taxa.loc[direction.isna(), 'direction']))}"
        )
    taxa["inc_dec"] = direction
    taxa["profile_id"] = taxa["study_id"] + "|" + taxa["experiment_id"] + "|" + taxa["inc_dec"]
    taxa = taxa[taxa["profile_id"].isin(profiles.index)]
    taxa = taxa.drop_duplicates(["profile_id", "taxon_id"])[["profile_id", "taxon_id"]]
    taxa = taxa[taxa["taxon_id"].isin(matrix.index)].reset_index(drop=True)
    labels = profiles.loc[taxa["profile_id"], [*keys, "inc_dec", "condition", "compartment"]]
    taxa = pd.concat([taxa, labels.reset_index(drop=True)], axis=1)

    trait_labels = read_tsv(paths["trait_columns"]).set_index("column_id")
    manifest = read_tsv(paths["contrast_manifest"]).set_index(keys)
    versions = set(manifest["classifier_version"])
    if versions != {CLASSIFIER_VERSION}:
        raise ValueError(
            f"contrast manifest written by {sorted(versions)}, but this checkout's "
            f"classifier is {CLASSIFIER_VERSION}; rerun `just contrast-orientation`"
        )
    return Corpus(profiles, taxa, matrix, trait_labels, manifest, label_source, condition_labels)


def msm_mentioned(manifest: pd.DataFrame) -> pd.Series:
    """Lexical flag: group names, definition or description mention MSM status."""
    text = manifest[["group_0_name", "group_1_name", "group_1_definition", "description"]]
    joined = text.agg(" | ".join, axis=1)
    return joined.map(lambda value: bool(MSM_PATTERN.search(value)))


def population(corpus: Corpus, name: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Profiles and taxa of one contrast population, with ``dir`` as the used direction.

    For the disease-versus-healthy populations, experiments classified ``reversed``
    are relabelled so ``dec`` always means lower in the disease group
    (mtp.contrast_orientation.oriented_profiles).
    """
    profiles = corpus.profiles
    if name == ALL_CONTRASTS:
        kept = profiles.assign(dir=profiles["inc_dec"])
    else:
        manifest = corpus.manifest
        index = pd.MultiIndex.from_frame(profiles[["study_id", "experiment_id"]])
        found = index.isin(manifest.index)
        if not found.all():
            raise ValueError(f"{int((~found).sum())} profiles have no contrast-manifest row")
        rows = manifest.reindex(index)
        keep = rows["contrast_class"].to_numpy() == "disease_vs_healthy"
        flags = rows["flags"].str.split(",")
        if name == STRONG:
            keep &= ~flags.map(lambda values: "weak_comparison_evidence" in values).to_numpy()
        if name == NO_MSM:
            keep &= ~msm_mentioned(rows).to_numpy()
        reversed_ = rows["orientation"].to_numpy() == "reversed"
        swapped = profiles["inc_dec"].map({"inc": "dec", "dec": "inc"})
        kept = profiles.assign(dir=np.where(reversed_, swapped, profiles["inc_dec"]))[keep]
    taxa = corpus.taxa[corpus.taxa["profile_id"].isin(kept.index)]
    taxa = taxa.assign(dir=taxa["profile_id"].map(kept["dir"]).to_numpy())
    return kept, taxa


def reclassified(
    corpus: Corpus, keys: Iterable[tuple[str, str]], contrast_class: str = "disease_vs_healthy"
) -> Corpus:
    """The corpus with the given experiments given another contrast class."""
    manifest = corpus.manifest.copy()
    chosen = list(keys)
    missing = [key for key in chosen if key not in manifest.index]
    if missing:
        raise ValueError(f"reclassified experiments missing from the manifest: {missing}")
    index = pd.MultiIndex.from_tuples(chosen, names=["study_id", "experiment_id"])
    manifest.loc[index, "contrast_class"] = contrast_class
    return dataclasses.replace(corpus, manifest=manifest)


def without_studies(corpus: Corpus, studies: Iterable[str]) -> Corpus:
    """The corpus without every profile of the given studies."""
    dropped = set(studies)
    profiles = corpus.profiles[~corpus.profiles["study_id"].isin(dropped)]
    taxa = corpus.taxa[corpus.taxa["profile_id"].isin(profiles.index)]
    return dataclasses.replace(corpus, profiles=profiles, taxa=taxa)


def with_synonyms(corpus: Corpus) -> Corpus:
    """The corpus with CONDITION_SYNONYMS counted under their disease."""
    profiles = corpus.profiles.assign(condition=merge_synonyms(corpus.profiles["condition"]))
    taxa = corpus.taxa.assign(condition=merge_synonyms(corpus.taxa["condition"]))
    return dataclasses.replace(corpus, profiles=profiles, taxa=taxa, condition_labels="synonyms")


# --------------------------------------------------------------------------- sets


def trait_taxa(matrix: pd.DataFrame, trait: str, ranks: Sequence[str] | None = None) -> set[str]:
    bearing = matrix[trait] > 0
    if ranks is not None:
        bearing &= matrix["taxonomic_rank"].isin(ranks)
    return set(matrix.index[bearing])


def reported(
    taxa: pd.DataFrame,
    condition: str | None,
    direction: str | None,
    compartment: str | None = FECES,
) -> pd.DataFrame:
    mask = pd.Series(True, index=taxa.index)
    if condition is not None:
        mask &= taxa["condition"].eq(condition)
    if direction is not None:
        mask &= taxa["dir"].eq(direction)
    if compartment is not None:
        mask &= taxa["compartment"].eq(compartment)
    return taxa[mask]


def condition_sets(
    taxa: pd.DataFrame, universe: set[str], compartment: str | None = FECES
) -> dict[str, set[str]]:
    return {
        short: set(reported(taxa, full, "dec", compartment)["taxon_id"]) & universe
        for short, full in CONDITIONS.items()
    }


def jaccard_distance(a: set[str], b: set[str]) -> float:
    union = len(a | b)
    return 1 - len(a & b) / union if union else float("nan")


def comparable_distance(a: set[str], b: set[str]) -> float:
    """Jaccard distance, undefined when either set is empty (no chance to overlap)."""
    return jaccard_distance(a, b) if a and b else float("nan")


def set_summary(sets: Mapping[str, set[str]]) -> dict[str, Any]:
    distances = {f"{a} vs. {b}": jaccard_distance(sets[a], sets[b]) for a, b in PAIRS}
    core = set.intersection(*sets.values())
    union = set.union(*sets.values())
    return {
        "sizes": {key: len(value) for key, value in sets.items()},
        "intersections": {f"{a} vs. {b}": len(sets[a] & sets[b]) for a, b in PAIRS},
        "unions": {f"{a} vs. {b}": len(sets[a] | sets[b]) for a, b in PAIRS},
        "jaccard_distance": distances,
        "mean_pairwise_jaccard_distance": float(np.mean(list(distances.values()))),
        "three_way_core": len(core),
        "three_way_union": len(union),
        "two_or_more_way": sum(
            1 for taxon in union if sum(taxon in value for value in sets.values()) >= 2
        ),
    }


# --------------------------------------------------------------------------- exact tests


def hypergeom_pair(n_universe: int, n_a: int, n_b: int, k: int) -> dict[str, float]:
    """Exact one-sided overlap test; d_J falls as k rises when set sizes are fixed."""
    support = np.arange(0, min(n_a, n_b) + 1)
    pmf = stats.hypergeom.pmf(support, n_universe, n_a, n_b)
    expected = n_a * n_b / n_universe
    distance = 1 - support / (n_a + n_b - support)
    rest = n_universe - n_a - n_b + k
    denominator = (n_a - k) * (n_b - k)
    return {
        "expected_intersection": expected,
        "expected_jaccard_distance": float((pmf * distance).sum()),
        "fold_enrichment": k / expected if expected else float("nan"),
        "odds_ratio": k * rest / denominator if denominator else float("inf"),
        "p": float(stats.hypergeom.sf(k - 1, n_universe, n_a, n_b)),
    }


def three_way_exact(n_universe: int, sizes: Sequence[int], k: int) -> dict[str, float]:
    """P(|A∩B∩C| >= k) by convolving the A∩B hypergeometric with a second draw."""
    n_a, n_b, n_c = sizes
    total = 0.0
    for x in range(min(n_a, n_b) + 1):
        weight = stats.hypergeom.pmf(x, n_universe, n_a, n_b)
        if weight > 0 and x >= k:
            total += weight * stats.hypergeom.sf(k - 1, n_universe, x, n_c)
    expected = n_a * n_b * n_c / n_universe**2
    return {
        "expected_3way_intersection": expected,
        "fold_3way_enrichment": k / expected,
        "p_3way_exact": total,
    }


# --------------------------------------------------------------------------- Monte Carlo


def choose(
    rng: np.random.Generator, n_iter: int, size: int, k: int, weights: FloatArray | None = None
) -> BoolArray:
    """``n_iter`` rows each marking ``k`` of ``size`` items, uniform or weighted.

    Weighted draws use Efraimidis-Spirakis keys ``log(u) / w``, which give the same
    distribution as successive weighted sampling without replacement.
    """
    noise = rng.random((n_iter, size))
    keys = noise if weights is None else np.log1p(-noise) / weights
    picked = np.argpartition(-keys, k - 1, axis=1)[:, :k] if k else np.zeros((n_iter, 0), int)
    mask = np.zeros((n_iter, size), dtype=bool)
    np.put_along_axis(mask, picked, True, axis=1)
    return mask


def mask_distance(a: BoolArray, b: BoolArray) -> FloatArray:
    union = (a | b).sum(axis=1)
    inter = (a & b).sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(union > 0, 1 - inter / np.maximum(union, 1), np.nan)


@dataclass
class NullDraws:
    """Per-iteration statistics accumulated over Monte Carlo batches."""

    pair_distance: dict[str, list[FloatArray]]
    pair_intersection: dict[str, list[IntArray]]
    mean_distance: list[FloatArray]
    core: list[IntArray]
    sizes: dict[str, list[IntArray]]

    @classmethod
    def empty(cls) -> NullDraws:
        names = [f"{a} vs. {b}" for a, b in PAIRS]
        return cls(
            {name: [] for name in names},
            {name: [] for name in names},
            [],
            [],
            {short: [] for short in SHORT},
        )

    def add(self, masks: Mapping[str, BoolArray]) -> None:
        distances = []
        for a, b in PAIRS:
            name = f"{a} vs. {b}"
            distance = mask_distance(masks[a], masks[b])
            distances.append(distance)
            self.pair_distance[name].append(distance)
            self.pair_intersection[name].append((masks[a] & masks[b]).sum(axis=1))
        self.mean_distance.append(np.nanmean(np.vstack(distances), axis=0))
        self.core.append((masks["PD"] & masks["HIV"] & masks["CRC"]).sum(axis=1))
        for short in SHORT:
            self.sizes[short].append(masks[short].sum(axis=1))

    def size_matched(self, observed: Mapping[str, Any], tolerance: float) -> tuple[NullDraws, int]:
        """Only the iterations whose set sizes are all within ``tolerance`` of observed."""
        sizes = {short: np.concatenate(values) for short, values in self.sizes.items()}
        keep = np.ones(next(iter(sizes.values())).size, dtype=bool)
        for short, values in sizes.items():
            target = observed["sizes"][short]
            keep &= np.abs(values - target) <= max(2, round(tolerance * target))
        matched = NullDraws.empty()
        for name in self.pair_distance:
            matched.pair_distance[name] = [np.concatenate(self.pair_distance[name])[keep]]
            matched.pair_intersection[name] = [np.concatenate(self.pair_intersection[name])[keep]]
        matched.mean_distance = [np.concatenate(self.mean_distance)[keep]]
        matched.core = [np.concatenate(self.core)[keep]]
        matched.sizes = {short: [values[keep]] for short, values in sizes.items()}
        return matched, int(keep.sum())

    def summary(self, observed: Mapping[str, Any]) -> dict[str, Any]:
        def p_low(draws: FloatArray, value: float) -> float:
            return float((1 + np.sum(draws <= value + 1e-12)) / (1 + draws.size))

        def p_high(draws: IntArray, value: int) -> float:
            return float((1 + np.sum(draws >= value)) / (1 + draws.size))

        result: dict[str, Any] = {"pairs": {}, "iterations": int(np.concatenate(self.core).size)}
        for name in self.pair_distance:
            distance = np.concatenate(self.pair_distance[name])
            result["pairs"][name] = {
                "expected_intersection": float(np.concatenate(self.pair_intersection[name]).mean()),
                "expected_jaccard_distance": float(np.nanmean(distance)),
                "p": p_low(distance[~np.isnan(distance)], observed["jaccard_distance"][name]),
            }
        mean_distance = np.concatenate(self.mean_distance)
        core = np.concatenate(self.core)
        result["mean_pairwise_jaccard_distance"] = {
            "expected": float(mean_distance.mean()),
            "sd": float(mean_distance.std(ddof=1)),
            "z": float(
                (observed["mean_pairwise_jaccard_distance"] - mean_distance.mean())
                / mean_distance.std(ddof=1)
            ),
            "p": p_low(mean_distance, observed["mean_pairwise_jaccard_distance"]),
            "p2_5": float(np.percentile(mean_distance, 2.5)),
            "p97_5": float(np.percentile(mean_distance, 97.5)),
        }
        result["three_way_core"] = {
            "expected": float(core.mean()),
            "p": p_high(core, observed["three_way_core"]),
            "p2_5": float(np.percentile(core, 2.5)),
            "p97_5": float(np.percentile(core, 97.5)),
        }
        result["synthetic_set_sizes"] = {
            short: {
                "mean": float(np.concatenate(values).mean()),
                "p2_5": float(np.percentile(np.concatenate(values), 2.5)),
                "p97_5": float(np.percentile(np.concatenate(values), 97.5)),
            }
            for short, values in self.sizes.items()
        }
        return result


def batches(total: int) -> Iterable[int]:
    done = 0
    while done < total:
        size = min(BATCH, total - done)
        done += size
        yield size


def uniform_null(
    rng: np.random.Generator, n_universe: int, sizes: Mapping[str, int], total: int
) -> NullDraws:
    draws = NullDraws.empty()
    for size in batches(total):
        draws.add({short: choose(rng, size, n_universe, sizes[short]) for short in SHORT})
    return draws


def weighted_null(
    rng: np.random.Generator, weights: FloatArray, sizes: Mapping[str, int], total: int
) -> NullDraws:
    draws = NullDraws.empty()
    for size in batches(total):
        draws.add(
            {short: choose(rng, size, weights.size, sizes[short], weights) for short in SHORT}
        )
    return draws


def u5_weights(feces_taxa: pd.DataFrame, cellulose: set[str]) -> FloatArray:
    """U5's weights over U3, the taxa an observed set can contain (#567).

    Each taxon weighs the number of studies outside the three conditions that report
    it decreased in feces, plus 0.5, so a taxon no other study reports can be drawn.
    """
    u3_taxa = sorted(set(feces_taxa.loc[feces_taxa["dir"].eq("dec"), "taxon_id"]) & cellulose)
    background = feces_taxa[
        feces_taxa["dir"].eq("dec") & ~feces_taxa["condition"].isin(CONDITIONS.values())
    ]
    background_studies = background.groupby("taxon_id")["study_id"].nunique()
    weights: FloatArray = np.asarray(
        background_studies.reindex(u3_taxa, fill_value=0), dtype=np.float64
    )
    return weights + 0.5


def shared_trait_null(
    rng: np.random.Generator, base: Mapping[str, BoolArray], trait_size: int, total: int
) -> NullDraws:
    """U4: a random trait on ``trait_size`` of the union of the baseline sets."""
    draws = NullDraws.empty()
    width = next(iter(base.values())).size
    for size in batches(total):
        trait = choose(rng, size, width, trait_size)
        draws.add({short: trait & base[short][None, :] for short in SHORT})
    return draws


def subset_null(
    rng: np.random.Generator, base: Mapping[str, BoolArray], sizes: Mapping[str, int], total: int
) -> NullDraws:
    """Random subsets of each condition's own baseline set, of the observed sizes."""
    draws = NullDraws.empty()
    width = next(iter(base.values())).size
    for size in batches(total):
        masks: dict[str, BoolArray] = {}
        for short in SHORT:
            members = np.flatnonzero(base[short])
            picked = choose(rng, size, members.size, sizes[short])
            mask = np.zeros((size, width), dtype=bool)
            mask[:, members] = picked
            masks[short] = mask
        draws.add(masks)
    return draws


def signature_null(
    rng: np.random.Generator,
    units: BoolArray,
    unit_study: IntArray,
    counts: Mapping[str, int],
    total: int,
) -> NullDraws:
    """U6: whole study signatures, keeping co-reported taxa together.

    Each iteration draws ``sum(counts)`` distinct studies, gives each one of its
    condition-level units at random, and splits them among the three synthetic
    conditions. As in the observed data, no study contributes to two conditions.
    """
    # Units of one study must be contiguous and the studies numbered 0..n-1, which the
    # searchsorted below relies on (#632).
    if np.any(np.diff(unit_study) < 0) or not np.array_equal(
        np.unique(unit_study), np.arange(int(unit_study.max()) + 1)
    ):
        raise ValueError("unit_study must be sorted and number the studies 0..n-1")
    draws = NullDraws.empty()
    n_studies = int(unit_study.max()) + 1
    first = np.searchsorted(unit_study, np.arange(n_studies))
    per_study = np.bincount(unit_study, minlength=n_studies)
    for size in batches(total):
        unit_of_study = first + (rng.random((size, n_studies)) * per_study).astype(np.int64)
        order = np.argsort(rng.random((size, n_studies)), axis=1)
        start = 0
        masks: dict[str, BoolArray] = {}
        for short in SHORT:
            studies = order[:, start : start + counts[short]]
            chosen = np.take_along_axis(unit_of_study, studies, axis=1)
            masks[short] = np.asarray(units[chosen].any(axis=1), dtype=bool)
            start += counts[short]
        draws.add(masks)
    return draws


# --------------------------------------------------------------------------- direction tests


def binomial_power(n: int, share: float, alpha: float = ALPHA) -> float:
    """Power of the exact two-sided sign test with ``n`` untied pairs at a true share."""
    if n == 0:
        return float("nan")
    rejecting = [k for k in range(n + 1) if stats.binomtest(k, n, 0.5).pvalue <= alpha]
    return float(sum(stats.binom.pmf(k, n, share) for k in rejecting))


def sign_tests(deltas: pd.Series) -> dict[str, float | int]:
    """Counts and tests for ``inc - dec`` deltas; ``down`` means dec > inc."""
    down = int((deltas < 0).sum())
    up = int((deltas > 0).sum())
    tied = int((deltas == 0).sum())
    untied = down + up
    nonzero = deltas[deltas != 0].to_numpy(dtype=float)
    nan = float("nan")
    return {
        "down": down,
        "up": up,
        "tied": tied,
        "share_down_untied": down / untied if untied else nan,
        "mean_delta": float(deltas.mean()) if len(deltas) else nan,
        "binom_p_down_one_sided": (
            float(stats.binomtest(down, untied, 0.5, alternative="greater").pvalue)
            if untied
            else nan
        ),
        "binom_p_two_sided": float(stats.binomtest(down, untied, 0.5).pvalue) if untied else nan,
        "wilcoxon_p_down_one_sided": (
            float(stats.wilcoxon(nonzero, alternative="less").pvalue) if nonzero.size else nan
        ),
        "wilcoxon_p_two_sided": float(stats.wilcoxon(nonzero).pvalue) if nonzero.size else nan,
    }


def paired_deltas(profiles: pd.DataFrame) -> pd.Series:
    """``frac(inc) - frac(dec)`` for experiments with both directions, indexed by key."""
    counts = profiles.groupby(["study_id", "experiment_id", "dir"]).size()
    if (counts > 1).any():
        raise ValueError("an experiment has two profiles in the same direction")
    wide = profiles.pivot_table(
        index=["study_id", "experiment_id"], columns="dir", values="frac", aggfunc="first"
    )
    if not {"inc", "dec"} <= set(wide.columns):
        keys = pd.MultiIndex.from_arrays([[], []], names=["study_id", "experiment_id"])
        return pd.Series(dtype=float, index=keys)
    wide = wide.dropna(subset=["inc", "dec"])
    return wide["inc"] - wide["dec"]


def study_votes(deltas: pd.Series, rule: str = "mean") -> pd.Series:
    """One signed value per study: its mean delta, or its majority of experiment signs.

    Under ``majority`` a study votes -1 (down) when more of its untied experiments
    are down than up, +1 in the reverse case, and 0 on a tie.
    """
    by_study = deltas.groupby(level="study_id")
    if rule == "mean":
        means: pd.Series = by_study.mean()
        return means
    votes: pd.Series = by_study.agg(
        lambda values: int(np.sign((values > 0).sum() - (values < 0).sum()))
    )
    return votes.astype(float)


def direction_row(
    profiles: pd.DataFrame, taxa: pd.DataFrame, cellulose: set[str], share: float
) -> dict[str, Any]:
    dec = profiles[profiles["dir"].eq("dec")]
    inc = profiles[profiles["dir"].eq("inc")]
    deltas = paired_deltas(profiles)
    experiment = sign_tests(deltas)
    study = sign_tests(study_votes(deltas))
    vote = sign_tests(study_votes(deltas, "majority"))
    untied = int(experiment["down"]) + int(experiment["up"])
    nan = float("nan")
    mwu_one = (
        float(stats.mannwhitneyu(dec["frac"], inc["frac"], alternative="greater").pvalue)
        if len(dec) and len(inc)
        else nan
    )
    mwu_two = (
        float(stats.mannwhitneyu(dec["frac"], inc["frac"], alternative="two-sided").pvalue)
        if len(dec) and len(inc)
        else nan
    )
    power = binomial_power(untied, share)
    min_p = min(1.0, 2 * 0.5**untied) if untied else nan
    if not untied:
        note = "no untied paired experiments"
    elif min_p > ALPHA:
        pairs = "pair" if untied == 1 else "pairs"
        note = f"p < {ALPHA} unattainable with {untied} untied {pairs}"
    elif power < 0.8:
        note = f"power {power:.2f} at the pooled feces share"
    else:
        note = ""
    dec_taxa = set(taxa.loc[taxa["dir"].eq("dec"), "taxon_id"]) & cellulose
    return {
        "n_studies_dec": dec["study_id"].nunique(),
        "n_profiles_dec": len(dec),
        "n_studies_inc": inc["study_id"].nunique(),
        "n_profiles_inc": len(inc),
        "mean_cellulose_frac_dec": dec["frac"].mean() if len(dec) else nan,
        "mean_cellulose_frac_inc": inc["frac"].mean() if len(inc) else nan,
        "mean_cellulose_z_dec": dec["z"].mean() if len(dec) else nan,
        "mean_cellulose_z_inc": inc["z"].mean() if len(inc) else nan,
        "uniq_dec_cellulose_taxa": len(dec_taxa),
        "paired_experiments": len(deltas),
        "paired_exp_down_dec_gt_inc": experiment["down"],
        "paired_exp_up_inc_gt_dec": experiment["up"],
        "paired_exp_tied": experiment["tied"],
        "share_down_untied": experiment["share_down_untied"],
        "mean_paired_delta_inc_minus_dec": experiment["mean_delta"],
        "binom_p_down_one_sided": experiment["binom_p_down_one_sided"],
        "binom_p_two_sided": experiment["binom_p_two_sided"],
        "wilcoxon_p_down_one_sided": experiment["wilcoxon_p_down_one_sided"],
        "wilcoxon_p_two_sided": experiment["wilcoxon_p_two_sided"],
        "mwu_p_dec_gt_inc_one_sided": mwu_one,
        "mwu_p_two_sided": mwu_two,
        "paired_studies": int(study["down"]) + int(study["up"]) + int(study["tied"]),
        "study_down": study["down"],
        "study_up": study["up"],
        "study_tied": study["tied"],
        "study_share_down_untied": study["share_down_untied"],
        "study_mean_delta_inc_minus_dec": study["mean_delta"],
        "study_binom_p_two_sided": study["binom_p_two_sided"],
        "study_wilcoxon_p_two_sided": study["wilcoxon_p_two_sided"],
        "study_majority_down": vote["down"],
        "study_majority_up": vote["up"],
        "study_majority_tied": vote["tied"],
        "study_majority_binom_p_two_sided": vote["binom_p_two_sided"],
        "min_attainable_binom_p_two_sided": min_p,
        "binom_power_at_pooled_feces_share": power,
        "power_note": note,
    }


def compartment_table(corpus: Corpus, cellulose: set[str]) -> pd.DataFrame:
    rows = []
    for name in (*POPULATIONS, NO_MSM):
        profiles, taxa = population(corpus, name)
        target = profiles[profiles["condition"].isin(CONDITIONS.values())]
        feces = target[target["compartment"].eq(FECES)]
        feces_deltas = paired_deltas(feces)
        untied = feces_deltas[feces_deltas != 0]
        share = float((untied < 0).mean()) if len(untied) else 0.5
        compartments = (FECES,) if name == NO_MSM else COMPARTMENTS
        for compartment in compartments:
            sets = condition_sets(taxa, cellulose, compartment)
            distances = [comparable_distance(sets[a], sets[b]) for a, b in PAIRS]
            defined = [value for value in distances if not math.isnan(value)]
            core = set.intersection(*sets.values())
            union = set.union(*sets.values())
            in_compartment = set(taxa.loc[taxa["compartment"].eq(compartment), "taxon_id"])
            shared = {
                "compartment_pairs_with_both_sets": len(defined),
                "compartment_mean_pairwise_jaccard_distance": (
                    float(np.mean(defined)) if defined else float("nan")
                ),
                "compartment_3way_shared_taxa": len(core),
                "compartment_union_taxa": len(union),
                # Cellulose taxa reported in this compartment, any condition or direction.
                "compartment_cellulose_universe_N": len(in_compartment & cellulose),
            }
            strata = [(short, [full]) for short, full in CONDITIONS.items()]
            strata.append((POOLED, list(CONDITIONS.values())))
            for short, members in strata:
                stratum = target[
                    target["compartment"].eq(compartment) & target["condition"].isin(members)
                ]
                stratum_taxa = taxa[
                    taxa["compartment"].eq(compartment) & taxa["condition"].isin(members)
                ]
                row: dict[str, Any] = {
                    "contrast_population": name,
                    "compartment": compartment,
                    "condition": POOLED if short == POOLED else CONDITIONS[short],
                    "short_condition": "Pooled" if short == POOLED else short,
                }
                row.update(direction_row(stratum, stratum_taxa, cellulose, share))
                row.update(shared)
                rows.append(row)
    return pd.DataFrame(rows)


def log_choose(n: int, k: IntArray) -> FloatArray:
    result: FloatArray = (
        special.gammaln(n + 1) - special.gammaln(k + 1) - special.gammaln(n - k + 1)
    )
    return result


def stratified_exact(tables: Sequence[Sequence[Sequence[int]]]) -> dict[str, float]:
    """Exact conditional inference on a common odds ratio over 2 x 2 strata.

    Given every stratum's margins, the top-left counts are independent Fisher
    noncentral hypergeometric variables, and their sum is sufficient for a common
    odds ratio. Returns the Mantel-Haenszel estimate, the conditional maximum-
    likelihood estimate with its exact central 95% interval, and the two-sided
    p-value at an odds ratio of 1 (outcomes no more probable than the observed
    one, as in ``scipy.stats.fisher_exact``). With one stratum the estimate and
    interval equal ``scipy.stats.contingency.odds_ratio``.
    """
    strata: list[tuple[int, int, int, int, int]] = []
    observed = 0
    mh_numerator = mh_denominator = 0.0
    for table in tables:
        (a, b), (c, d) = table
        total = a + b + c + d
        if total == 0:
            continue
        mh_numerator += a * d / total
        mh_denominator += b * c / total
        first_row, first_column = a + b, a + c
        low = max(0, first_row - (total - first_column))
        high = min(first_row, first_column)
        if low < high:
            strata.append((total, first_column, first_row, low, high))
            observed += a
    nan = float("nan")
    if mh_denominator:
        mantel_haenszel = mh_numerator / mh_denominator
    else:
        mantel_haenszel = math.inf if mh_numerator else nan
    if not strata:
        return {
            "mantel_haenszel_or": mantel_haenszel,
            "conditional_or": nan,
            "low": nan,
            "high": nan,
            "p": nan,
        }

    def distribution(log_odds: float) -> tuple[FloatArray, int]:
        dist: FloatArray = np.array([1.0])
        offset = 0
        for total, first_column, first_row, low, high in strata:
            support = np.arange(low, high + 1)
            log_weight = (
                log_choose(first_column, support)
                + log_choose(total - first_column, first_row - support)
                + support * log_odds
            )
            weight = np.exp(log_weight - log_weight.max())
            dist = np.convolve(dist, weight / weight.sum())
            offset += low
        return dist, offset

    null, offset = distribution(0.0)
    smallest, largest = offset, offset + null.size - 1
    p_observed = null[observed - offset]
    p_value = float(min(1.0, null[null <= p_observed * (1 + 1e-7)].sum()))

    def mean(log_odds: float) -> float:
        dist, start = distribution(log_odds)
        return float((np.arange(start, start + dist.size) * dist).sum())

    def upper_tail(log_odds: float) -> float:
        dist, start = distribution(log_odds)
        return float(dist[observed - start :].sum())

    def lower_tail(log_odds: float) -> float:
        dist, start = distribution(log_odds)
        return float(dist[: observed - start + 1].sum())

    def solve(function: Callable[[float], float], target: float) -> float:
        return float(np.exp(optimize.brentq(lambda x: function(x) - target, -50, 50)))

    if observed == smallest:
        estimate, low_ci, high_ci = 0.0, 0.0, solve(lower_tail, 0.025)
    elif observed == largest:
        estimate, low_ci, high_ci = math.inf, solve(upper_tail, 0.025), math.inf
    else:
        estimate = solve(mean, observed)
        low_ci, high_ci = solve(upper_tail, 0.025), solve(lower_tail, 0.025)
    return {
        "mantel_haenszel_or": mantel_haenszel,
        "conditional_or": estimate,
        "low": low_ci,
        "high": high_ci,
        "p": p_value,
    }


def direction_counts(deltas: pd.Series, unit: str) -> list[int]:
    """[down, up] untied units of a paired-delta series, per experiment or per study."""
    values = deltas if unit == "experiment" else study_votes(deltas)
    return [int((values < 0).sum()), int((values > 0).sum())]


POOLED_STRATIFIED = "Pooled, stratified by condition"


def interaction_table(corpus: Corpus) -> pd.DataFrame:
    """Fisher's exact test of direction x compartment, Feces against each other compartment.

    The crude pooled rows add the three conditions' tables, so they mix condition with
    compartment when the conditions sample the compartments unevenly. The stratified
    rows test a common odds ratio across the conditions instead.
    """
    rows = []
    nan = float("nan")
    for name in (ALL_CONTRASTS, ORIENTED):
        profiles, _ = population(corpus, name)
        target = profiles[profiles["condition"].isin(CONDITIONS.values())]
        by_condition: dict[tuple[str, str], list[list[list[int]]]] = {}
        strata = [(short, [full]) for short, full in CONDITIONS.items()]
        strata.append(("Pooled", list(CONDITIONS.values())))
        for short, members in strata:
            in_condition = target[target["condition"].isin(members)]
            feces = in_condition[in_condition["compartment"].eq(FECES)]
            for other in COMPARTMENTS[1:]:
                rest = in_condition[in_condition["compartment"].eq(other)]
                for unit in ("experiment", "study"):
                    table = [
                        direction_counts(paired_deltas(feces), unit),
                        direction_counts(paired_deltas(rest), unit),
                    ]
                    if short != "Pooled":
                        by_condition.setdefault((other, unit), []).append(table)
                    usable = sum(table[1]) > 0 and sum(table[0]) > 0
                    odds, p = stats.fisher_exact(table) if usable else (nan, nan)
                    conditional, low, high = nan, nan, nan
                    if usable:
                        estimate = stats.contingency.odds_ratio(table)
                        interval = estimate.confidence_interval()
                        conditional = float(estimate.statistic)
                        low, high = float(interval.low), float(interval.high)
                    shared = set(feces["study_id"]) & set(rest["study_id"])
                    note = "tied pairs excluded; counts treat each unit as independent"
                    if short == "Pooled":
                        note += "; crude pooling mixes condition with compartment"
                    if shared:
                        note += "; some studies sampled both compartments"
                    rows.append(
                        interaction_row(name, short, unit, other, table)
                        | {
                            "odds_ratio": float(odds),
                            "conditional_or": conditional,
                            "conditional_or_ci95_low": low,
                            "conditional_or_ci95_high": high,
                            "fisher_p_two_sided": float(p),
                            "studies_in_both_compartments": len(shared),
                            "note": note,
                        }
                    )
        for (other, unit), tables in by_condition.items():
            exact = stratified_exact(tables)
            summed = [
                [sum(table[0][0] for table in tables), sum(table[0][1] for table in tables)],
                [sum(table[1][0] for table in tables), sum(table[1][1] for table in tables)],
            ]
            rows.append(
                interaction_row(name, POOLED_STRATIFIED, unit, other, summed)
                | {
                    "odds_ratio": exact["mantel_haenszel_or"],
                    "conditional_or": exact["conditional_or"],
                    "conditional_or_ci95_low": exact["low"],
                    "conditional_or_ci95_high": exact["high"],
                    "fisher_p_two_sided": exact["p"],
                    "studies_in_both_compartments": nan,
                    "note": (
                        "exact conditional test of a common odds ratio across PD, HIV and CRC; "
                        "odds_ratio is the Mantel-Haenszel estimate; tied pairs excluded"
                    ),
                }
            )
    return pd.DataFrame(rows)


def interaction_row(
    population_name: str, short: str, unit: str, other: str, table: Sequence[Sequence[int]]
) -> dict[str, Any]:
    nan = float("nan")
    condition = {
        "Pooled": "Pooled (PD + HIV + CRC)",
        POOLED_STRATIFIED: "Pooled (PD + HIV + CRC), stratified by condition",
    }.get(short, CONDITIONS.get(short, short))
    return {
        "contrast_population": population_name,
        "condition": condition,
        "short_condition": short,
        "unit": unit,
        "compartment_a": FECES,
        "compartment_b": other,
        "a_down": table[0][0],
        "a_up": table[0][1],
        "b_down": table[1][0],
        "b_up": table[1][1],
        "a_share_down": table[0][0] / sum(table[0]) if sum(table[0]) else nan,
        "b_share_down": table[1][0] / sum(table[1]) if sum(table[1]) else nan,
    }


# --------------------------------------------------------------------------- tables


def pairwise_table(
    sets: Mapping[str, set[str]],
    summary: Mapping[str, Any],
    universes: Mapping[str, int],
    nulls: Mapping[str, Mapping[str, Any]],
) -> pd.DataFrame:
    rows = []
    for a, b in PAIRS:
        name = f"{a} vs. {b}"
        row: dict[str, Any] = {
            "endpoint": "pairwise_jaccard_distance",
            "comparison": name,
            "condition_a": CONDITIONS[a],
            "condition_b": CONDITIONS[b],
            "n_a": len(sets[a]),
            "n_b": len(sets[b]),
            "intersection": summary["intersections"][name],
            "union": summary["unions"][name],
            "jaccard_similarity": 1 - summary["jaccard_distance"][name],
            "jaccard_distance": summary["jaccard_distance"][name],
        }
        for key, label in (
            ("U1", "u1_matrix"),
            ("U2", "u2_feces_cell"),
            ("U3", "u3_feces_dec_cell"),
        ):
            exact = hypergeom_pair(
                universes[key], len(sets[a]), len(sets[b]), summary["intersections"][name]
            )
            prefix = label.split("_")[0]
            row[f"{label}_N"] = universes[key]
            row[f"{prefix}_expected_intersection"] = exact["expected_intersection"]
            row[f"{prefix}_expected_jaccard_distance"] = exact["expected_jaccard_distance"]
            row[f"{prefix}_fold_enrichment"] = exact["fold_enrichment"]
            row[f"{prefix}_odds_ratio"] = exact["odds_ratio"]
            row[f"{prefix}_hypergeom_p"] = exact["p"]
        for key, prefix in NULL_COLUMNS.items():
            pair = nulls[key]["pairs"][name]
            row[f"{prefix}_expected_intersection"] = pair["expected_intersection"]
            row[f"{prefix}_expected_jaccard_distance"] = pair["expected_jaccard_distance"]
            row[f"{prefix}_p"] = pair["p"]
        rows.append(row)
    mean_row: dict[str, Any] = {
        "endpoint": "mean_pairwise_jaccard_distance",
        "comparison": "Mean across all 3 pairs",
        "jaccard_similarity": 1 - summary["mean_pairwise_jaccard_distance"],
        "jaccard_distance": summary["mean_pairwise_jaccard_distance"],
    }
    core_row: dict[str, Any] = {
        "endpoint": "three_way_core",
        "comparison": "3-way core (PD ∩ HIV ∩ CRC)",
        "n_a": len(sets["PD"]),
        "n_b": len(sets["HIV"]),
        "n_c": len(sets["CRC"]),
        "intersection": summary["three_way_core"],
        "union": summary["three_way_union"],
    }
    sizes = [len(sets[short]) for short in SHORT]
    for key, label in (("U1", "u1_matrix"), ("U2", "u2_feces_cell"), ("U3", "u3_feces_dec_cell")):
        prefix = label.split("_")[0]
        exact = three_way_exact(universes[key], sizes, summary["three_way_core"])
        mc = nulls[key]
        core_row[f"{label}_N"] = universes[key]
        core_row[f"{prefix}_expected_intersection"] = exact["expected_3way_intersection"]
        core_row[f"{prefix}_fold_enrichment"] = exact["fold_3way_enrichment"]
        core_row[f"{prefix}_hypergeom_p"] = exact["p_3way_exact"]
        mean_row[f"{label}_N"] = universes[key]
        mean_row[f"{prefix}_expected_jaccard_distance"] = mc["mean_pairwise_jaccard_distance"][
            "expected"
        ]
        mean_row[f"{prefix}_z"] = mc["mean_pairwise_jaccard_distance"]["z"]
        mean_row[f"{prefix}_mc_p"] = mc["mean_pairwise_jaccard_distance"]["p"]
    for key, prefix in NULL_COLUMNS.items():
        mean_row[f"{prefix}_expected_jaccard_distance"] = nulls[key][
            "mean_pairwise_jaccard_distance"
        ]["expected"]
        mean_row[f"{prefix}_p"] = nulls[key]["mean_pairwise_jaccard_distance"]["p"]
        core_row[f"{prefix}_expected_intersection"] = nulls[key]["three_way_core"]["expected"]
        core_row[f"{prefix}_p"] = nulls[key]["three_way_core"]["p"]
    rows.extend([mean_row, core_row])
    return pd.DataFrame(rows)


NULL_COLUMNS = {
    "U4": "u4_shared_trait",
    "U5": "u5_emp_weighted",
    "within": "within_cond_random_subset",
    "U6": "u6_signature_perm",
    "U6_size_matched": "u6_size_matched",
}
U6_SIZE_TOLERANCE = 0.1


def sensitivity_table(
    definitions: Sequence[tuple[str, str, Mapping[str, set[str]], int, int]],
) -> pd.DataFrame:
    rows = []
    for slug, label, sets, matrix_n, empirical_n in definitions:
        summary = set_summary(sets)
        for a, b in PAIRS:
            name = f"{a} vs. {b}"
            exact_m = hypergeom_pair(
                matrix_n, len(sets[a]), len(sets[b]), summary["intersections"][name]
            )
            exact_e = hypergeom_pair(
                empirical_n, len(sets[a]), len(sets[b]), summary["intersections"][name]
            )
            rows.append(
                {
                    "definition_slug": slug,
                    "definition_label": label,
                    "comparison": name,
                    "n_a": len(sets[a]),
                    "n_b": len(sets[b]),
                    "intersection": summary["intersections"][name],
                    "union": summary["unions"][name],
                    "jaccard_similarity": 1 - summary["jaccard_distance"][name],
                    "jaccard_distance": summary["jaccard_distance"][name],
                    "matrix_universe_N": matrix_n,
                    "matrix_expected_intersection": exact_m["expected_intersection"],
                    "matrix_expected_jaccard_distance": exact_m["expected_jaccard_distance"],
                    "matrix_fold_enrichment": exact_m["fold_enrichment"],
                    "matrix_hypergeom_p": exact_m["p"],
                    "empirical_universe_N": empirical_n,
                    "empirical_expected_intersection": exact_e["expected_intersection"],
                    "empirical_expected_jaccard_distance": exact_e["expected_jaccard_distance"],
                    "empirical_fold_enrichment": exact_e["fold_enrichment"],
                    "empirical_hypergeom_p": exact_e["p"],
                }
            )
        sizes = [len(sets[short]) for short in SHORT]
        exact_m = three_way_exact(matrix_n, sizes, summary["three_way_core"])
        exact_e = three_way_exact(empirical_n, sizes, summary["three_way_core"])
        rows.append(
            {
                "definition_slug": slug,
                "definition_label": label,
                "comparison": "Mean across all 3 pairs (3-way core)",
                "n_a": sizes[0],
                "n_b": sizes[1],
                "n_c": sizes[2],
                "intersection": summary["three_way_core"],
                "union": summary["three_way_union"],
                "jaccard_similarity": 1 - summary["mean_pairwise_jaccard_distance"],
                "jaccard_distance": summary["mean_pairwise_jaccard_distance"],
                "matrix_universe_N": matrix_n,
                "matrix_expected_intersection": exact_m["expected_3way_intersection"],
                "matrix_fold_enrichment": exact_m["fold_3way_enrichment"],
                "matrix_hypergeom_p": exact_m["p_3way_exact"],
                "empirical_universe_N": empirical_n,
                "empirical_expected_intersection": exact_e["expected_3way_intersection"],
                "empirical_fold_enrichment": exact_e["fold_3way_enrichment"],
                "empirical_hypergeom_p": exact_e["p_3way_exact"],
            }
        )
    return pd.DataFrame(rows)


def net_decreased(taxa: pd.DataFrame, universe: set[str], unit: str) -> dict[str, set[str]]:
    """Taxa with more decreased than increased ``unit``s (study_id or profile_id) per condition."""
    result = {}
    for short, full in CONDITIONS.items():
        feces = reported(taxa, full, None)
        feces = feces[feces["taxon_id"].isin(universe)]
        counts = feces.groupby(["taxon_id", "dir"])[unit].nunique().unstack(fill_value=0)
        for column in ("dec", "inc"):
            if column not in counts:
                counts[column] = 0
        result[short] = set(counts.index[counts["dec"] > counts["inc"]])
    return result


def direct_trait_subjects(edges: Path, trait: str) -> set[str]:
    """NCBITaxon subjects with a KG edge ``predicate -> object`` for a matrix trait column.

    Matrix counts include every annotated taxon at or below a node, so a genus can
    bear a trait only through its species or strains; this edge list separates
    direct annotation from inheritance (#565).
    """
    predicate, object_id = trait.split("|", 1)
    subjects = set()
    with edges.open(encoding="utf-8") as handle:
        for line in handle:
            if object_id not in line:
                continue
            fields = line.rstrip("\n").split("\t")
            if (
                fields[2] == predicate
                and fields[3] == object_id
                and fields[1].startswith("NCBITaxon:")
            ):
                subjects.add(fields[1])
    if not subjects:
        raise ValueError(f"no {predicate} -> {object_id} edges in {edges}")
    return subjects


def strictly_net_decreased(row: Mapping[str, Any], unit: str) -> bool:
    """Reported decreased in more studies or profiles (``unit``) than increased, in each of
    PD, HIV and CRC; a tie is not net decreased (#632)."""
    return all(row[f"{short}_dec_{unit}"] > row[f"{short}_inc_{unit}"] for short in SHORT)


def catalog_table(
    taxa: pd.DataFrame, matrix: pd.DataFrame, sets: Mapping[str, set[str]], direct: set[str]
) -> pd.DataFrame:
    union = set.union(*sets.values())
    feces = taxa[taxa["compartment"].eq(FECES) & taxa["taxon_id"].isin(union)]
    names = matrix["taxon_name"].to_dict()
    ranks = matrix["taxonomic_rank"].to_dict()
    bearing = matrix[TRAIT_ID].to_dict()
    rows = []
    for taxon in sorted(union):
        row: dict[str, Any] = {
            "taxon_id": taxon,
            "taxon_name": names[taxon],
            "taxonomic_rank": ranks[taxon],
            "matrix_trait_bearing_taxa_at_or_below": int(bearing[taxon]),
            "n_conditions_feces_dec": sum(taxon in sets[short] for short in SHORT),
            "is_3way_core": all(taxon in sets[short] for short in SHORT),
        }
        for short in SHORT:
            row[f"in_{short}_feces_dec"] = taxon in sets[short]
        for short, full in CONDITIONS.items():
            mine = feces[feces["taxon_id"].eq(taxon) & feces["condition"].eq(full)]
            for direction in ("dec", "inc"):
                chosen = mine[mine["dir"].eq(direction)]
                row[f"{short}_{direction}_studies"] = chosen["study_id"].nunique()
                row[f"{short}_{direction}_profiles"] = chosen["profile_id"].nunique()
        row["has_direct_trait_edge"] = taxon in direct
        row["is_direct_or_genus_species"] = (
            row["has_direct_trait_edge"] or row["taxonomic_rank"] in GENUS_SPECIES
        )
        row["is_strictly_net_decreased_all3"] = strictly_net_decreased(row, "studies")
        row["is_net_decreased_all3_by_profiles"] = strictly_net_decreased(row, "profiles")
        if not row["is_3way_core"]:
            row["core_class"] = "not in 3-way core"
        elif row["taxonomic_rank"] in GENUS_SPECIES:
            row["core_class"] = "3-way core, genus or species rank"
        else:
            row["core_class"] = "3-way core, family-to-phylum rank"
        rows.append(row)
    frame = pd.DataFrame(rows)
    frame["total_dec_studies"] = frame[[f"{short}_dec_studies" for short in SHORT]].sum(axis=1)
    frame = frame.sort_values(
        ["n_conditions_feces_dec", "total_dec_studies", "taxon_name"],
        ascending=[False, False, True],
    )
    return frame.drop(columns="total_dec_studies").reset_index(drop=True)


def profile_pair_table(taxa: pd.DataFrame, cellulose: set[str]) -> pd.DataFrame:
    rows = []
    feces = reported(taxa, None, "dec")
    feces = feces[feces["taxon_id"].isin(cellulose)]
    by_profile = {
        short: feces[feces["condition"].eq(full)]
        .groupby("profile_id")["taxon_id"]
        .agg(set)
        .tolist()
        for short, full in CONDITIONS.items()
    }
    by_study = {
        short: feces[feces["condition"].eq(full)].groupby("study_id")["taxon_id"].agg(set).tolist()
        for short, full in CONDITIONS.items()
    }
    for a, b in PAIRS:
        row: dict[str, Any] = {"comparison": f"{a} vs. {b}"}
        for unit, plural, groups in (
            ("profile", "profiles", by_profile),
            ("study", "studies", by_study),
        ):
            distances = [jaccard_distance(x, y) for x in groups[a] for y in groups[b]]
            overlap = [bool(x & y) for x in groups[a] for y in groups[b]]
            row[f"n_{plural}_a"] = len(groups[a])
            row[f"n_{plural}_b"] = len(groups[b])
            row[f"n_{unit}_pairs"] = len(distances)
            row[f"mean_{unit}_pair_jaccard_distance"] = float(np.mean(distances))
            row[f"mean_{unit}_pair_jaccard_similarity"] = 1 - float(np.mean(distances))
            row[f"frac_{unit}_pairs_with_overlap"] = float(np.mean(overlap))
        rows.append(row)
    return pd.DataFrame(rows)


def benchmark_table(taxa: pd.DataFrame, matrix: pd.DataFrame, labels: pd.DataFrame) -> pd.DataFrame:
    """Rank every trait with >= 10 feces-decreased taxa per condition by mean pairwise d_J."""
    traits = [column for column in matrix.columns if column not in ("taxon_name", "taxonomic_rank")]
    reported_taxa = {
        short: set(reported(taxa, full, "dec")["taxon_id"]) for short, full in CONDITIONS.items()
    }
    union = sorted(set.union(*reported_taxa.values()))
    bearing = matrix.loc[union, traits].to_numpy() > 0
    member = {short: np.isin(union, sorted(values)) for short, values in reported_taxa.items()}
    size = {short: (bearing & member[short][:, None]).sum(axis=0) for short in SHORT}
    keep = np.minimum.reduce([size[short] for short in SHORT]) >= 10
    frame = pd.DataFrame({"trait_id": np.asarray(traits)[keep]})
    missing = sorted(set(frame["trait_id"]) - set(labels.index))
    if missing:
        raise ValueError(f"traits without a taxa_trait_matrix_columns.tsv row: {missing}")
    frame["trait_label"] = labels.loc[frame["trait_id"], "object_name"].to_numpy()
    if frame["trait_label"].eq("").any():
        raise ValueError("a benchmark trait has an empty object_name")
    frame["trait_predicate"] = labels.loc[frame["trait_id"], "predicate_name"].to_numpy()
    frame["is_growth_medium"] = frame["trait_id"].str.contains("mediadive.medium", regex=False)
    frame["matrix_taxon_count"] = (matrix.loc[:, frame["trait_id"]] > 0).sum(axis=0).to_numpy()
    for short in SHORT:
        frame[f"{short}_dec_taxa"] = size[short][keep]
    all_three = member["PD"] & member["HIV"] & member["CRC"]
    any_three = member["PD"] | member["HIV"] | member["CRC"]
    frame["three_way_intersection"] = (bearing & all_three[:, None]).sum(axis=0)[keep]
    frame["three_way_union"] = (bearing & any_three[:, None]).sum(axis=0)[keep]
    frame["three_way_jaccard_similarity"] = (
        frame["three_way_intersection"] / frame["three_way_union"]
    )
    distances = []
    for a, b in PAIRS:
        inter = (bearing & (member[a] & member[b])[:, None]).sum(axis=0)[keep]
        union_ab = (bearing & (member[a] | member[b])[:, None]).sum(axis=0)[keep]
        frame[f"d_J_{a}_vs_{b}"] = 1 - inter / union_ab
        distances.append(frame[f"d_J_{a}_vs_{b}"])
    frame["mean_pairwise_jaccard_distance"] = pd.concat(distances, axis=1).mean(axis=1)
    frame["mean_pairwise_jaccard_similarity"] = 1 - frame["mean_pairwise_jaccard_distance"]
    frame = frame.sort_values(["mean_pairwise_jaccard_distance", "trait_id"]).reset_index(drop=True)
    frame["rank_by_convergence"] = np.arange(1, len(frame) + 1)
    biological = ~frame["is_growth_medium"]
    frame["rank_excluding_growth_media"] = pd.Series(pd.NA, index=frame.index, dtype="Int64")
    frame.loc[biological, "rank_excluding_growth_media"] = np.arange(1, int(biological.sum()) + 1)
    frame["percentile_excluding_growth_media"] = frame["rank_excluding_growth_media"] / int(
        biological.sum()
    )
    return frame


def four_pair_tables(
    taxa: pd.DataFrame, matrix: pd.DataFrame, studies: pd.DataFrame, dump: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Taxon drivers and study evidence for the four highlighted condition-trait pairs."""
    pairs = (
        ("1_PD_Feces_Cellulose_dec", "PD", TRAIT_ID, "dec"),
        ("2_HIV_Feces_Cellulose_dec", "HIV", TRAIT_ID, "dec"),
        ("3_CRC_Feces_Cellulose_dec", "CRC", TRAIT_ID, "dec"),
        ("4_PD_Feces_ModeratelyHalophilic_inc", "PD", HALOPHILIC_ID, "inc"),
    )
    names = {"dec": "decreased", "inc": "increased"}
    taxon_names = matrix["taxon_name"].to_dict()
    ranks = matrix["taxonomic_rank"].to_dict()
    info = studies.set_index("study_id")
    driver_frames: list[pd.DataFrame] = []
    evidence_rows: list[dict[str, Any]] = []
    summary: dict[str, Any] = {}
    for pair_id, short, trait, target in pairs:
        opposite = "inc" if target == "dec" else "dec"
        bearing = trait_taxa(matrix, trait)
        weights = matrix[trait].to_dict()
        feces = reported(taxa, CONDITIONS[short], None)
        feces = feces[feces["taxon_id"].isin(bearing)]
        hits = feces[feces["dir"].eq(target)]
        rows = []
        for taxon, group in hits.groupby("taxon_id"):
            other = feces[feces["taxon_id"].eq(taxon) & feces["dir"].eq(opposite)]
            study_ids = sorted(group["study_id"].unique())
            rows.append(
                {
                    "pair_id": pair_id,
                    "condition": CONDITIONS[short],
                    "trait_id": trait,
                    "target_direction": names[target],
                    "taxon_id": taxon,
                    "taxon_name": taxon_names[taxon],
                    "taxonomic_rank": ranks[taxon],
                    "kg_rolled_up_trait_count": int(weights[taxon]),
                    "target_dir_studies": len(study_ids),
                    "target_dir_profiles": group["profile_id"].nunique(),
                    "opposite_dir_studies": other["study_id"].nunique(),
                    "opposite_dir_profiles": other["profile_id"].nunique(),
                    "supporting_study_ids": "; ".join(study_ids),
                    "supporting_pmids": "; ".join(info.loc[study_ids, "pmid"]),
                    "supporting_dois": "; ".join(info.loc[study_ids, "doi"]),
                }
            )
        drivers = pd.DataFrame(rows).sort_values(
            ["target_dir_studies", "target_dir_profiles", "taxon_name"],
            ascending=[False, False, True],
        )
        driver_frames.append(drivers)
        for study_id, group in hits.groupby("study_id"):
            other = feces[feces["study_id"].eq(study_id) & feces["dir"].eq(opposite)]
            experiments = set(feces.loc[feces["study_id"].eq(study_id), "experiment_id"])
            rows_dump = dump[dump["experiment_key"].isin(experiments)]
            row = info.loc[study_id]
            evidence_rows.append(
                {
                    "pair_id": pair_id,
                    "condition": CONDITIONS[short],
                    "trait_id": trait,
                    "target_direction": names[target],
                    "study_id": study_id,
                    "pmid": row["pmid"],
                    "doi": row["doi"],
                    "year": row["publication_year"],
                    "authors": first_value(rows_dump, "Authors list"),
                    "journal": first_value(rows_dump, "Journal"),
                    "title": row["title"],
                    "sequencing_types": joined_values(rows_dump, "Sequencing type"),
                    "locations": joined_values(rows_dump, "Location of subjects"),
                    "n_target_taxa": group["taxon_id"].nunique(),
                    "n_opposite_taxa": other["taxon_id"].nunique(),
                    "target_taxa_list": taxon_list(group, taxon_names, ranks, weights),
                    "opposite_taxa_list": taxon_list(other, taxon_names, ranks, weights),
                }
            )
        summary[pair_id] = {
            "condition": CONDITIONS[short],
            "trait_id": trait,
            "target_direction": names[target],
            "body_site": FECES,
            "n_supporting_studies_in_target_dir": hits["study_id"].nunique(),
            "n_unique_target_taxa": hits["taxon_id"].nunique(),
            "top_12_taxa": drivers.head(12).to_dict("records"),
        }
    return pd.concat(driver_frames, ignore_index=True), pd.DataFrame(evidence_rows), summary


RANK_CODES = {
    "species": "S",
    "genus": "G",
    "family": "F",
    "order": "O",
    "class": "C",
    "phylum": "P",
}


def taxon_list(
    frame: pd.DataFrame,
    names: Mapping[Any, Any],
    ranks: Mapping[Any, Any],
    weights: Mapping[Any, Any],
) -> str:
    """``Name (rank code, w=trait-bearing taxa at or below)`` for each distinct taxon."""
    parts = [
        f"{names[taxon]} ({RANK_CODES.get(str(ranks[taxon]), '?')}, w={int(weights[taxon])})"
        for taxon in frame["taxon_id"].drop_duplicates()
    ]
    return "; ".join(sorted(parts))


def first_value(frame: pd.DataFrame, column: str) -> str:
    values = [value for value in frame[column] if value]
    return values[0] if values else ""


def joined_values(frame: pd.DataFrame, column: str) -> str:
    values = sorted(
        {part.strip() for value in frame[column] for part in value.split(",") if part.strip()}
    )
    return ", ".join(values)


def dump_vintage(path: Path) -> str:
    """The release stamp on the first line of full_dump.csv, e.g. 2026-04-24_00:41_UTC."""
    with path.open(encoding="utf-8") as handle:
        first = handle.readline().strip()
    match = re.match(r"^#\s*BugSigDB\s+([^,\s]+)", first)
    if match is None:
        raise ValueError(f"{path} does not start with a '# BugSigDB <vintage>' line")
    return match.group(1)


def load_dump(path: Path) -> pd.DataFrame:
    dump = pd.read_csv(path, dtype=str, keep_default_na=False, comment=None, skiprows=1)
    parts = dump["BSDB ID"].str.extract(r"^bsdb:([^/]+)/(\d+)/(\d+)$")
    dump["experiment_key"] = "bsdb:" + parts[0] + "-" + parts[1]
    return dump


# --------------------------------------------------------------------------- summary rows


def fmt_p(value: float) -> str:
    if math.isnan(value):
        return "NA"
    return f"{value:.3g}"


def tier(value: float) -> str:
    if math.isnan(value):
        return "not tested"
    for bound, label in ((1e-3, "p < 0.001"), (1e-2, "p < 0.01"), (5e-2, "p < 0.05")):
        if value < bound:
            return label
    return "not significant (p >= 0.05)"


def definition_range(scenarios: pd.DataFrame, population_name: str, short: str, unit: str) -> str:
    """How far a feces direction p moves across the condition-definition scenarios."""
    chosen = scenarios[
        scenarios["contrast_population"].eq(population_name)
        & scenarios["short_condition"].eq(short)
        & scenarios["unit"].eq(unit)
    ]["binom_p_two_sided"]
    return f"{fmt_p(float(chosen.min()))} to {fmt_p(float(chosen.max()))}"


def direction_summary_rows(
    compartments: pd.DataFrame, scenarios: pd.DataFrame
) -> list[dict[str, Any]]:
    rows = []
    feces = compartments[
        compartments["compartment"].eq(FECES)
        & compartments["contrast_population"].isin([ALL_CONTRASTS, ORIENTED])
    ]
    for _, row in feces.iterrows():
        for unit in ("experiment", "study"):
            if unit == "experiment":
                down, up, tied = (
                    row["paired_exp_down_dec_gt_inc"],
                    row["paired_exp_up_inc_gt_dec"],
                    row["paired_exp_tied"],
                )
                share, p = row["share_down_untied"], row["binom_p_two_sided"]
                extra = f"Wilcoxon 2-sided p = {fmt_p(row['wilcoxon_p_two_sided'])}"
            else:
                down, up, tied = row["study_down"], row["study_up"], row["study_tied"]
                share, p = row["study_share_down_untied"], row["study_binom_p_two_sided"]
                study_p = fmt_p(row["study_wilcoxon_p_two_sided"])
                extra = f"Wilcoxon 2-sided p on study means = {study_p}"
            population_name = row["contrast_population"]
            spread = definition_range(scenarios, population_name, row["short_condition"], unit)
            units = "experiments" if unit == "experiment" else "studies"
            mean_delta = row["mean_paired_delta_inc_minus_dec"]
            effect = (
                f"mean delta (inc - dec annotation share) = {mean_delta:.5f}"
                if unit == "experiment"
                else f"mean of study means = {row['study_mean_delta_inc_minus_dec']:.5f}"
            )
            rows.append(
                {
                    "result_category": "1_within_condition_direction",
                    "test_id": f"feces_{row['short_condition']}_{unit}_direction_{population_name}",
                    "contrast_population": population_name,
                    "unit_of_observation": (
                        "paired experiment" if unit == "experiment" else "study (mean delta)"
                    ),
                    "scope_or_comparison": (
                        f"{row['condition']}, Feces ({row['n_studies_dec']} dec / "
                        f"{row['n_studies_inc']} inc studies; {row['n_profiles_dec']} dec / "
                        f"{row['n_profiles_inc']} inc profiles)"
                    ),
                    "observed_statistic": f"{down} down vs {up} up ({tied} tied); {share:.1%} down",
                    "null_or_reference_expectation": (
                        "50% of untied units down (exact two-sided sign test)"
                    ),
                    "effect_size_metric": effect,
                    "p_value": p,
                    "significance_tier": tier(p),
                    "summary_interpretation": (
                        "Taxa annotated with cellulose catabolic process formed a larger "
                        "annotation share in the decreased than in the increased signature in "
                        f"{share:.1%} of untied {units}; {extra}. Across the condition "
                        "definitions of condition_definition_sensitivity.tsv this p runs from "
                        f"{spread}. Describes curated signature annotations, not measured "
                        "cellulose degradation."
                    ),
                }
            )
    return rows


def convergence_summary_rows(pairwise: pd.DataFrame) -> list[dict[str, Any]]:
    rows = []
    for _, row in pairwise.iterrows():
        if row["endpoint"] == "pairwise_jaccard_distance":
            rows.append(
                {
                    "result_category": "2_cross_condition_overlap",
                    "test_id": f"jaccard_{row['comparison'].replace(' vs. ', '_vs_')}",
                    "contrast_population": ALL_CONTRASTS,
                    "unit_of_observation": "taxon label",
                    "scope_or_comparison": (
                        f"{row['comparison']} (feces decreased, cellulose-annotated taxa)"
                    ),
                    "observed_statistic": (
                        f"{int(row['intersection'])}/{int(row['union'])} shared; "
                        f"d_J = {row['jaccard_distance']:.4f}"
                    ),
                    "null_or_reference_expectation": (
                        f"U1 expected d_J {row['u1_expected_jaccard_distance']:.4f}; "
                        f"U5 {row['u5_emp_weighted_expected_jaccard_distance']:.4f}; "
                        f"U6 {row['u6_signature_perm_expected_jaccard_distance']:.4f}"
                    ),
                    "effect_size_metric": f"fold over U1 = {row['u1_fold_enrichment']:.2f}",
                    "p_value": row["u6_signature_perm_p"],
                    "significance_tier": tier(row["u6_signature_perm_p"]),
                    "summary_interpretation": (
                        f"U1 exact p = {fmt_p(row['u1_hypergeom_p'])}; "
                        f"U5 p = {fmt_p(row['u5_emp_weighted_p'])}; "
                        "U6 (whole-signature permutation) p = "
                        f"{fmt_p(row['u6_signature_perm_p'])}. Overlap of taxonomic labels "
                        "ever reported decreased, not consistent loss."
                    ),
                }
            )
        else:
            endpoint = row["endpoint"]
            observed = (
                f"mean d_J = {row['jaccard_distance']:.4f}"
                if endpoint == "mean_pairwise_jaccard_distance"
                else f"{int(row['intersection'])}/{int(row['union'])} taxa in all three"
            )
            rows.append(
                {
                    "result_category": "2_cross_condition_overlap",
                    "test_id": endpoint,
                    "contrast_population": ALL_CONTRASTS,
                    "unit_of_observation": "taxon label",
                    "scope_or_comparison": (
                        "PD, HIV and CRC feces decreased, cellulose-annotated taxa"
                    ),
                    "observed_statistic": observed,
                    "null_or_reference_expectation": (
                        f"U4 {row['u4_shared_trait_p']:.3g}; U5 {row['u5_emp_weighted_p']:.3g}; "
                        f"U6 {row['u6_signature_perm_p']:.3g} (p-values)"
                    ),
                    "effect_size_metric": (
                        f"U6 expected {row['u6_signature_perm_expected_jaccard_distance']:.4f}"
                        if endpoint == "mean_pairwise_jaccard_distance"
                        else "U6 expected core "
                        f"{row['u6_signature_perm_expected_intersection']:.2f}"
                    ),
                    "p_value": row["u6_signature_perm_p"],
                    "significance_tier": tier(row["u6_signature_perm_p"]),
                    "summary_interpretation": (
                        "U6 keeps each study's co-reported taxa together; U1-U3 treat taxa "
                        "as independent draws, which overstates significance."
                    ),
                }
            )
    return rows


def interaction_summary_rows(interaction: pd.DataFrame) -> list[dict[str, Any]]:
    rows = []
    chosen = interaction[
        interaction["compartment_b"].eq(NON_FECAL_GI)
        & interaction["short_condition"].isin(["CRC", POOLED_STRATIFIED])
    ]
    for _, row in chosen.iterrows():
        stratified = row["short_condition"] == POOLED_STRATIFIED
        slug = "pooled_stratified" if stratified else row["short_condition"]
        population_name = row["contrast_population"]
        # Studies that sampled both compartments; a stratified row sums its conditions'.
        same = interaction[
            interaction["contrast_population"].eq(population_name)
            & interaction["unit"].eq(row["unit"])
            & interaction["compartment_b"].eq(row["compartment_b"])
            & interaction["short_condition"].isin(SHORT)
        ]
        shared = int(
            same["studies_in_both_compartments"].sum()
            if stratified
            else row["studies_in_both_compartments"]
        )
        estimate = (
            f"common odds ratio (conditional MLE) = {row['conditional_or']:.2f} "
            f"(95% CI {row['conditional_or_ci95_low']:.2f}-{row['conditional_or_ci95_high']:.2f})"
        )
        rows.append(
            {
                "result_category": "3_compartment_interaction",
                "test_id": f"fisher_feces_vs_gi_{slug}_{row['unit']}_{population_name}",
                "contrast_population": population_name,
                "unit_of_observation": "paired experiment"
                if row["unit"] == "experiment"
                else "study (mean delta)",
                "scope_or_comparison": f"{row['condition']}: Feces vs {NON_FECAL_GI}",
                "observed_statistic": (
                    f"Feces {row['a_down']} down : {row['a_up']} up vs "
                    f"GI {row['b_down']} down : {row['b_up']} up"
                ),
                "null_or_reference_expectation": (
                    "common odds ratio 1 across conditions (exact conditional test)"
                    if stratified
                    else "same down share in both compartments (Fisher's exact test)"
                ),
                "effect_size_metric": estimate,
                "p_value": row["fisher_p_two_sided"],
                "significance_tier": tier(row["fisher_p_two_sided"]),
                "summary_interpretation": interaction_interpretation(row, shared),
            }
        )
    return rows


def interaction_interpretation(row: pd.Series | Mapping[str, Any], shared: int = 0) -> str:
    """What one Feces-vs-GI row shows, given its p, interval, unit and population (#631).

    ``shared`` counts the studies that sampled both compartments, so a study-level row's
    two compartments are not fully independent samples.
    """
    p = row["fisher_p_two_sided"]
    low, high = row["conditional_or_ci95_low"], row["conditional_or_ci95_high"]
    if math.isnan(p):
        finding = "Not testable with these counts"
    elif p >= ALPHA:
        finding = (
            "No difference in direction between stool and GI mucosa/lumen samples "
            "detected at p < 0.05"
        )
    elif low <= 1 <= high:
        finding = (
            "Direction differs between stool and GI mucosa/lumen samples at p < 0.05, "
            "but the 95% interval of the odds ratio includes 1"
        )
    else:
        finding = "Direction differs between stool and GI mucosa/lumen samples"
    caveats = []
    if row["unit"] == "experiment":
        caveats.append("experiments of one study are treated as independent")
    elif shared:
        caveats.append(
            f"{shared} {'study' if shared == 1 else 'studies'} sampled both compartments"
        )
    if row["contrast_population"] == ALL_CONTRASTS:
        caveats.append("most CRC mucosa/lumen experiments do not compare patients with controls")
    caveats.append(
        "sampling differs between compartments, so this is not a test of anatomical specificity"
    )
    return f"{finding}; {'; '.join(caveats)}."


def benjamini_hochberg(p_values: pd.Series) -> pd.Series:
    """BH-adjusted p-values over the non-missing entries; NaN stays NaN."""
    values = p_values.astype(float)
    present = values.dropna()
    order = present.sort_values()
    ranked = order * len(order) / np.arange(1, len(order) + 1)
    adjusted = np.minimum.accumulate(ranked.to_numpy()[::-1])[::-1].clip(max=1.0)
    result = pd.Series(np.nan, index=values.index)
    result.loc[order.index] = adjusted
    return result


# --------------------------------------------------------------------------- output


def clean(value: Any) -> Any:
    """JSON-safe values with floats rounded to six significant digits.

    An infinite value is written as the string "inf" or "-inf", so that it is not mistaken
    for a missing one, which is null (#631).
    """
    if isinstance(value, dict):
        return {str(key): clean(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(item) for item in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, (float, np.floating)):
        number = float(value)
        if math.isnan(number):
            return None
        if math.isinf(number):
            return "inf" if number > 0 else "-inf"
        return float(f"{number:.6g}")
    if value is pd.NA:
        return None
    return value


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False, float_format=FLOAT_FORMAT, lineterminator="\n")


def write_json(payload: Mapping[str, Any], path: Path) -> None:
    path.write_text(
        json.dumps(clean(payload), indent=1, ensure_ascii=False) + "\n", encoding="utf-8"
    )


# --------------------------------------------------------------------------- main


def headline(corpus: Corpus, cellulose: set[str]) -> dict[str, Any]:
    """Label-dependent headline numbers, for comparing label sources."""
    _, taxa = population(corpus, ALL_CONTRASTS)
    sets = condition_sets(taxa, cellulose)
    summary = set_summary(sets)
    profiles, _ = population(corpus, ALL_CONTRASTS)
    rows = {}
    for short, full in CONDITIONS.items():
        feces = profiles[profiles["condition"].eq(full) & profiles["compartment"].eq(FECES)]
        tests = sign_tests(paired_deltas(feces))
        rows[short] = {
            "feces_studies_dec": feces.loc[feces["dir"].eq("dec"), "study_id"].nunique(),
            "feces_studies_inc": feces.loc[feces["dir"].eq("inc"), "study_id"].nunique(),
            "paired_down": tests["down"],
            "paired_up": tests["up"],
            "paired_tied": tests["tied"],
            "binom_p_two_sided": tests["binom_p_two_sided"],
            "wilcoxon_p_two_sided": tests["wilcoxon_p_two_sided"],
        }
    return {"sets": summary, "feces_direction": rows}


def label_corrections(paths: Mapping[str, Path]) -> pd.DataFrame:
    """Experiments in the three conditions whose study-modal label differs from their own."""
    metadata = read_tsv(paths["profile_metadata"], index_col=0)
    experiments = read_tsv(paths["experiments"])
    keys = ["study_id", "experiment_id"]
    pairs = (
        metadata[keys]
        .drop_duplicates(keys)
        .join(study_modal_labels(paths["full_dump"]), on="study_id")
    )
    merged = pairs.merge(
        experiments[[*keys, "condition", "body_site"]],
        on=keys,
        suffixes=("_study_modal", "_experiment"),
    )
    differs = merged["condition_study_modal"].ne(merged["condition_experiment"]) | merged[
        "body_site_study_modal"
    ].ne(merged["body_site_experiment"])
    touches = merged["condition_study_modal"].isin(CONDITIONS.values()) | merged[
        "condition_experiment"
    ].isin(CONDITIONS.values())
    return merged[differs & touches].sort_values(keys).reset_index(drop=True)


def label_inventory(
    paths: Mapping[str, Path], corpus: Corpus, dump: pd.DataFrame, condition_labels: str
) -> pd.DataFrame:
    """Every human condition label naming one of the three diseases, with its status (#561).

    Counts are of retained human experiments (every one has a retained profile).
    """
    experiments = read_tsv(paths["experiments"])
    retained = corpus.profiles[["study_id", "experiment_id"]].drop_duplicates()
    human = experiments.merge(retained, on=["study_id", "experiment_id"], validate="one_to_one")
    named = human[human["condition"].str.contains(NEAR_LABEL_PATTERN)]
    ontology_ids = dump.groupby("Condition")["EFO ID"].agg(lambda ids: ", ".join(sorted(set(ids))))
    counted = "counted in this run" if condition_labels == "synonyms" else "sensitivity only"
    rows = []
    for label, group in named.groupby("condition"):
        label = str(label)
        if label in CONDITIONS.values():
            status = "exact label (primary analysis)"
        elif label in SYNONYM_TO_CONDITION:
            status = f"synonym of {SYNONYM_TO_CONDITION[label]} ({counted})"
        elif "," in label:
            status = "multi-condition label (excluded)"
        else:
            status = "different condition (excluded)"
        rows.append(
            {
                "condition_label": label,
                "ontology_ids": ontology_ids.get(label, ""),
                "status": status,
                "studies": group["study_id"].nunique(),
                "experiments": len(group),
                "feces_experiments": int(group["body_site"].eq(FECES).sum()),
            }
        )
    frame = pd.DataFrame(rows)
    return frame.sort_values(["status", "experiments"], ascending=[True, False]).reset_index(
        drop=True
    )


def u6_units(
    profiles: pd.DataFrame, taxa: pd.DataFrame, cellulose: set[str]
) -> tuple[BoolArray, IntArray, dict[str, int], pd.DataFrame]:
    """U6 units: one per study and condition with a feces decreased profile."""
    feces_taxa = taxa[taxa["compartment"].eq(FECES)]
    dec_feces_profiles = profiles[profiles["compartment"].eq(FECES) & profiles["dir"].eq("dec")]
    units = (
        dec_feces_profiles[["study_id", "condition"]]
        .drop_duplicates()
        .sort_values(["study_id", "condition"])
    )
    u3_taxa = sorted(set(feces_taxa.loc[feces_taxa["dir"].eq("dec"), "taxon_id"]) & cellulose)
    unit_taxa = (
        feces_taxa[feces_taxa["dir"].eq("dec") & feces_taxa["taxon_id"].isin(cellulose)]
        .groupby(["study_id", "condition"])["taxon_id"]
        .agg(set)
    )
    unit_matrix = np.zeros((len(units), len(u3_taxa)), dtype=bool)
    position = {taxon: index for index, taxon in enumerate(u3_taxa)}
    for row, unit in enumerate(units.itertuples(index=False, name=None)):
        for taxon in unit_taxa.get(unit, set()):
            unit_matrix[row, position[taxon]] = True
    unit_study = pd.factorize(units["study_id"], sort=True)[0].astype(np.int64)
    unit_counts = {
        short: int(units["condition"].eq(full).sum()) for short, full in CONDITIONS.items()
    }
    return unit_matrix, unit_study, unit_counts, units


def scenario_corpora(corpus: Corpus) -> list[tuple[str, str, Corpus]]:
    """Condition and contrast definitions that change which experiments count (#561, #566).

    ``corpus`` must use the exact condition labels; each scenario is built from it.
    """
    if corpus.condition_labels != "exact":
        raise ValueError("scenarios are defined relative to the exact condition labels")
    merged = with_synonyms(corpus)
    return [
        ("exact_labels", "Exact condition labels (primary analysis)", corpus),
        ("synonyms_merged", "Synonym labels counted under their disease", merged),
        (
            "hiv_minus_vs_plus_counted",
            "Exact labels; the HIV- vs HIV+ experiment counted as disease vs healthy",
            reclassified(corpus, HIV_RECLASSIFIED),
        ),
        (
            "synonyms_and_hiv_minus_vs_plus",
            "Synonyms counted and the HIV- vs HIV+ experiment counted",
            reclassified(merged, HIV_RECLASSIFIED),
        ),
        (
            "synonyms_hiv_minus_vs_plus_and_art",
            "Synonyms, the HIV- vs HIV+ experiment and three HIV+ on ART vs uninfected counted",
            reclassified(merged, [*HIV_RECLASSIFIED, *HIV_ON_ART]),
        ),
        (
            "without_site_miscurations",
            "Exact labels without bsdb:24451087 (rectal mucosa curated as Feces)",
            without_studies(corpus, SITE_MISCURATIONS),
        ),
        (
            "crc_within_patient_excluded",
            "Exact labels; the CRC tumour-vs-adjacent-tissue experiment not counted as disease "
            "vs healthy",
            reclassified(corpus, CRC_WITHIN_PATIENT, "within_patient_tissue"),
        ),
    ]


def scenario_tables(
    corpus: Corpus, cellulose: set[str], permutations: int, seed: int
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Feces direction, overlap and the Feces-vs-GI interaction under each scenario."""
    direction_rows: list[dict[str, Any]] = []
    overlap_rows: list[dict[str, Any]] = []
    interaction_frames: list[pd.DataFrame] = []
    for slug, label, scenario in scenario_corpora(corpus):
        interaction = interaction_table(scenario)
        interaction = interaction[
            interaction["compartment_b"].eq(NON_FECAL_GI)
            & interaction["short_condition"].isin([*SHORT, POOLED_STRATIFIED])
        ]
        interaction_frames.append(interaction.assign(scenario=slug, scenario_label=label))
        for name in (ALL_CONTRASTS, ORIENTED):
            profiles, _ = population(scenario, name)
            feces = profiles[profiles["compartment"].eq(FECES)]
            strata = [(short, [full]) for short, full in CONDITIONS.items()]
            strata.append(("Pooled", list(CONDITIONS.values())))
            for short, members in strata:
                deltas = paired_deltas(feces[feces["condition"].isin(members)])
                for unit, values in (("experiment", deltas), ("study", study_votes(deltas))):
                    tests = sign_tests(values)
                    direction_rows.append(
                        {
                            "scenario": slug,
                            "scenario_label": label,
                            "contrast_population": name,
                            "short_condition": short,
                            "unit": unit,
                            "down": tests["down"],
                            "up": tests["up"],
                            "tied": tests["tied"],
                            "share_down_untied": tests["share_down_untied"],
                            "binom_p_two_sided": tests["binom_p_two_sided"],
                            "wilcoxon_p_two_sided": tests["wilcoxon_p_two_sided"],
                        }
                    )
        profiles, taxa = population(scenario, ALL_CONTRASTS)
        sets = condition_sets(taxa, cellulose)
        summary = set_summary(sets)
        feces_taxa = taxa[taxa["compartment"].eq(FECES)]
        u3 = len(set(feces_taxa.loc[feces_taxa["dir"].eq("dec"), "taxon_id"]) & cellulose)
        exact = three_way_exact(
            u3, [len(sets[short]) for short in SHORT], summary["three_way_core"]
        )
        unit_matrix, unit_study, unit_counts, _ = u6_units(profiles, taxa, cellulose)
        u6 = signature_null(
            np.random.default_rng(seed), unit_matrix, unit_study, unit_counts, permutations
        ).summary(summary)
        overlap_rows.append(
            {
                "scenario": slug,
                "scenario_label": label,
                **{f"n_{short}": len(sets[short]) for short in SHORT},
                **{
                    f"intersection_{a}_{b}": summary["intersections"][f"{a} vs. {b}"]
                    for a, b in PAIRS
                },
                "three_way_core": summary["three_way_core"],
                "three_way_union": summary["three_way_union"],
                "mean_pairwise_jaccard_distance": summary["mean_pairwise_jaccard_distance"],
                "u3_N": u3,
                "u3_three_way_exact_p": exact["p_3way_exact"],
                **{f"u6_units_{short}": unit_counts[short] for short in SHORT},
                "u6_expected_mean_jaccard_distance": u6["mean_pairwise_jaccard_distance"][
                    "expected"
                ],
                "u6_mean_jaccard_distance_p": u6["mean_pairwise_jaccard_distance"]["p"],
                "u6_expected_core": u6["three_way_core"]["expected"],
                "u6_core_p": u6["three_way_core"]["p"],
            }
        )
    interactions = pd.concat(interaction_frames, ignore_index=True)
    leading = ["scenario", "scenario_label"]
    interactions = interactions[[*leading, *(c for c in interactions if c not in leading)]]
    return pd.DataFrame(direction_rows), pd.DataFrame(overlap_rows), interactions


def input_hashes(paths: Mapping[str, Path]) -> dict[str, dict[str, Any]]:
    return {
        key: {"path": INPUTS[key], "sha256": file_sha256(path), "bytes": path.stat().st_size}
        for key, path in paths.items()
    }


def check_input_drift(
    hashes: Mapping[str, Mapping[str, Any]], pinned_path: Path, allow: bool, repin: bool
) -> tuple[list[str], dict[str, Any]]:
    """The inputs that differ from the pins, and the pins this run records.

    The pins are the provenance's ``pinned_inputs`` (or ``inputs`` in a provenance
    written before pins were kept apart). Only ``repin`` replaces them. ``allow`` runs
    on drifted inputs but carries the pins forward unchanged, so a later run without
    it still compares against them. An input the pins do not list yet is added.
    """
    if repin:
        return [], dict(hashes)
    if not pinned_path.exists():
        raise FileNotFoundError(
            f"no pinned provenance at {pinned_path}; pass --repin to pin this run's inputs"
        )
    recorded = json.loads(pinned_path.read_text(encoding="utf-8"))
    pinned = recorded.get("pinned_inputs", recorded.get("inputs", {}))
    drifted = [
        key
        for key, value in hashes.items()
        if key in pinned and pinned[key]["sha256"] != value["sha256"]
    ]
    if drifted and not allow:
        details = "; ".join(
            f"{key}: {pinned[key]['sha256'][:12]} pinned, {hashes[key]['sha256'][:12]} now"
            for key in drifted
        )
        raise ValueError(
            f"inputs differ from {pinned_path}: {details}. Rebuild them at the pinned "
            "versions, pass --allow-input-drift for a run that keeps the pins, or --repin "
            "to adopt these inputs."
        )
    return drifted, {key: pinned.get(key, value) for key, value in hashes.items()}


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    paths = {key: args.input_root / relative for key, relative in INPUTS.items()}
    missing = [str(path) for path in paths.values() if not path.exists()]
    if missing:
        raise FileNotFoundError(f"missing inputs: {missing}")
    hashes = input_hashes(paths)
    drifted, pinned_inputs = check_input_drift(
        hashes, args.pinned_provenance, args.allow_input_drift, args.repin
    )
    rng_seeds = {
        name: args.seed + offset
        for offset, name in enumerate(("U1", "U2", "U3", "U4", "U5", "within", "U6"))
    }

    corpus = load_corpus(paths, args.label_source, args.condition_labels)
    # Scenarios and their labels are defined against the exact condition labels.
    exact_corpus = (
        corpus
        if args.condition_labels == "exact"
        else load_corpus(paths, args.label_source, "exact")
    )
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    matrix = corpus.matrix
    cellulose = trait_taxa(matrix, TRAIT_ID)
    genus_species = trait_taxa(matrix, TRAIT_ID, GENUS_SPECIES)
    all_matrix_taxa = set(matrix.index)

    profiles, taxa = population(corpus, ALL_CONTRASTS)
    feces_taxa = taxa[taxa["compartment"].eq(FECES)]
    universes = {
        "U1": len(cellulose),
        "U2": len(set(feces_taxa["taxon_id"]) & cellulose),
        "U3": len(set(feces_taxa.loc[feces_taxa["dir"].eq("dec"), "taxon_id"]) & cellulose),
    }
    sets = condition_sets(taxa, cellulose)
    summary = set_summary(sets)
    sizes = {short: len(values) for short, values in sets.items()}
    core_and_union = f"core {summary['three_way_core']}, union {summary['three_way_union']}"
    print(f"[{args.label_source}] sets {sizes}, {core_and_union}")

    # Baseline: every matrix taxon reported decreased in feces, any trait.
    base_sets = condition_sets(taxa, all_matrix_taxa)
    base_summary = set_summary(base_sets)
    base_union = sorted(set.union(*base_sets.values()))
    base_masks = {short: np.isin(base_union, sorted(values)) for short, values in base_sets.items()}

    weights = u5_weights(feces_taxa, cellulose)

    unit_matrix, unit_study, unit_counts, units = u6_units(profiles, taxa, cellulose)

    n = args.permutations
    nulls: dict[str, Any] = {}
    nulls["U1"] = uniform_null(
        np.random.default_rng(rng_seeds["U1"]), universes["U1"], sizes, n
    ).summary(summary)
    nulls["U2"] = uniform_null(
        np.random.default_rng(rng_seeds["U2"]), universes["U2"], sizes, n
    ).summary(summary)
    nulls["U3"] = uniform_null(
        np.random.default_rng(rng_seeds["U3"]), universes["U3"], sizes, n
    ).summary(summary)
    nulls["U4"] = shared_trait_null(
        np.random.default_rng(rng_seeds["U4"]), base_masks, summary["three_way_union"], n
    ).summary(summary)
    nulls["U5"] = weighted_null(np.random.default_rng(rng_seeds["U5"]), weights, sizes, n).summary(
        summary
    )
    nulls["within"] = subset_null(
        np.random.default_rng(rng_seeds["within"]), base_masks, sizes, n
    ).summary(summary)
    u6_draws = signature_null(
        np.random.default_rng(rng_seeds["U6"]), unit_matrix, unit_study, unit_counts, n
    )
    nulls["U6"] = u6_draws.summary(summary)
    matched, accepted = u6_draws.size_matched(summary, U6_SIZE_TOLERANCE)
    # At least 100 matched draws at the default 10^5 permutations, scaled down for
    # shorter runs (#567).
    needed = min(100, max(10, n // 1000))
    if accepted < needed:
        raise ValueError(f"only {accepted} U6 draws match the observed set sizes (need {needed})")
    nulls["U6_size_matched"] = {
        **matched.summary(summary),
        "tolerance": f"each set size within max(2, {U6_SIZE_TOLERANCE:.0%} of observed)",
    }

    pairwise = pairwise_table(sets, summary, universes, nulls)

    # U4 exact tests: cellulose union against the baseline 3-way and >= 2-way sets.
    def u4_exact(observed: int, base_count: int) -> dict[str, float]:
        size = summary["three_way_union"]
        total = len(base_union)
        expected = size * base_count / total
        table = [
            [observed, size - observed],
            [base_count - observed, total - size - base_count + observed],
        ]
        return {
            "n_universe": total,
            "baseline_count": base_count,
            "expected_intersection": expected,
            "fold_enrichment": observed / expected,
            "odds_ratio": float(stats.fisher_exact(table)[0]),
            "p_value_one_sided": float(stats.hypergeom.sf(observed - 1, total, base_count, size)),
        }

    three_way = {
        "U4_exact_3way_vs_baseline": u4_exact(
            summary["three_way_core"], base_summary["three_way_core"]
        ),
        "U4_exact_ge2way_vs_baseline": u4_exact(
            summary["two_or_more_way"], base_summary["two_or_more_way"]
        ),
    }
    for key in ("U1", "U2", "U3"):
        three_way[key] = {
            "n_universe": universes[key],
            **three_way_exact(universes[key], list(sizes.values()), summary["three_way_core"]),
            "mc_mean_pairwise_jaccard_distance": nulls[key]["mean_pairwise_jaccard_distance"],
        }

    # Sensitivity definitions.
    z_positive = taxa[taxa["profile_id"].map(profiles["z"]) > 0]
    cellulose_any_site = set(taxa["taxon_id"]) & cellulose
    oriented_profiles, oriented_taxa = population(corpus, ORIENTED)
    strong_profiles, strong_taxa = population(corpus, STRONG)
    del oriented_profiles, strong_profiles
    feces_gs = set(feces_taxa["taxon_id"]) & genus_species
    definitions = [
        (
            "1_feces_dec_all_ranks",
            "Feces dec, all cellulose-annotated ranks (primary)",
            sets,
            universes["U1"],
            universes["U2"],
        ),
        (
            "2_feces_dec_z_gt_0",
            "Feces dec, profiles with cellulose z > 0 only",
            condition_sets(z_positive, cellulose),
            universes["U1"],
            universes["U2"],
        ),
        (
            "3_feces_dec_genus_species",
            "Feces dec, genus + species ranks only",
            condition_sets(taxa, genus_species),
            len(genus_species),
            len(feces_gs),
        ),
        (
            "4_feces_net_decreased",
            "Feces net-decreased (dec profiles > inc profiles)",
            net_decreased(taxa, cellulose, "profile_id"),
            universes["U1"],
            universes["U2"],
        ),
        (
            "4b_feces_net_decreased_by_studies",
            "Feces net-decreased (dec studies > inc studies)",
            net_decreased(taxa, cellulose, "study_id"),
            universes["U1"],
            universes["U2"],
        ),
        (
            "5_all_body_sites_dec",
            "All human body sites dec, cellulose-annotated",
            condition_sets(taxa, cellulose, None),
            universes["U1"],
            len(cellulose_any_site),
        ),
        (
            "6_baseline_all_feces_dec_taxa",
            "Baseline: all Feces dec matrix taxa (any trait)",
            base_sets,
            len(all_matrix_taxa),
            len(set(feces_taxa["taxon_id"])),
        ),
        (
            "7_feces_dec_disease_vs_healthy_oriented",
            "Feces dec, disease-vs-healthy contrasts only (oriented)",
            condition_sets(oriented_taxa, cellulose),
            universes["U1"],
            universes["U2"],
        ),
        (
            "8_feces_dec_disease_vs_healthy_strong",
            "Feces dec, disease-vs-healthy contrasts with strong comparison evidence",
            condition_sets(strong_taxa, cellulose),
            universes["U1"],
            universes["U2"],
        ),
        (
            "9_feces_dec_without_site_miscurations",
            "Feces dec, without bsdb:24451087 (rectal mucosa curated as Feces)",
            condition_sets(taxa[~taxa["study_id"].isin(SITE_MISCURATIONS)], cellulose),
            universes["U1"],
            universes["U2"],
        ),
        (
            "10_feces_dec_synonyms_merged"
            if args.condition_labels == "exact"
            else "10_feces_dec_exact_labels",
            "Feces dec, synonym condition labels counted under their disease"
            if args.condition_labels == "exact"
            else "Feces dec, exact condition labels only",
            condition_sets(
                population(
                    with_synonyms(exact_corpus)
                    if args.condition_labels == "exact"
                    else exact_corpus,
                    ALL_CONTRASTS,
                )[1],
                cellulose,
            ),
            universes["U1"],
            universes["U2"],
        ),
    ]
    sensitivity = sensitivity_table(definitions)

    compartments = compartment_table(corpus, cellulose)
    interaction = interaction_table(corpus)
    direct = direct_trait_subjects(paths["kg_edges"], TRAIT_ID)
    catalog = catalog_table(taxa, matrix, sets, direct)
    bearing_ranks = matrix.loc[sorted(cellulose), "taxonomic_rank"]
    direct_by_rank = {
        str(rank): {
            "trait_bearing_matrix_taxa": len(group),
            "with_direct_edge": int(group.index.isin(direct).sum()),
        }
        for rank, group in bearing_ranks.groupby(bearing_ranks)
    }
    scenario_direction, scenario_overlap, scenario_interaction = scenario_tables(
        exact_corpus, cellulose, n, rng_seeds["U6"]
    )
    dump = load_dump(paths["full_dump"])
    inventory = label_inventory(paths, corpus, dump, args.condition_labels)
    profile_pairs = profile_pair_table(taxa, cellulose)
    benchmark = benchmark_table(taxa, matrix, corpus.trait_labels)
    studies = read_tsv(paths["studies"])
    drivers, evidence, four_pairs = four_pair_tables(taxa, matrix, studies, dump)

    summary_rows = (
        direction_summary_rows(compartments, scenario_direction)
        + convergence_summary_rows(pairwise)
        + interaction_summary_rows(interaction)
    )
    significant = pd.DataFrame(summary_rows)
    significant.insert(
        list(significant.columns).index("significance_tier") + 1,
        "bh_q_within_table",
        benjamini_hochberg(significant["p_value"]),
    )

    # Label-source comparison, without the Monte Carlo nulls.
    other_source = "study_modal" if args.label_source == "experiment" else "experiment"
    comparison = {
        args.label_source: headline(corpus, cellulose),
        other_source: headline(load_corpus(paths, other_source, args.condition_labels), cellulose),
    }
    corrections = label_corrections(paths)

    target_sites = (
        profiles[profiles["condition"].isin(CONDITIONS.values())][
            ["condition", "body_site", "compartment"]
        ]
        .drop_duplicates()
        .sort_values(["condition", "compartment", "body_site"])
    )
    multi_condition = profiles[
        profiles["compartment"].eq(FECES)
        & ~profiles["condition"].isin(CONDITIONS.values())
        & profiles["condition"].str.contains(
            "|".join(map(re.escape, CONDITIONS.values())), regex=True
        )
    ]

    benchmark_row = benchmark[benchmark["trait_id"].eq(TRAIT_ID)].iloc[0]
    biological = benchmark[~benchmark["is_growth_medium"]]
    trait_label = corpus.trait_labels.at[TRAIT_ID, "object_name"]
    result = {
        "trait_id": TRAIT_ID,
        "trait_label": trait_label,
        "trait_column_label": corpus.trait_labels.at[TRAIT_ID, "column_label"],
        "label_source": args.label_source,
        "condition_labels": args.condition_labels,
        "contrast_classifier_version": CLASSIFIER_VERSION,
        "bugsigdb_vintage": dump_vintage(paths["full_dump"]),
        "conditions": list(CONDITIONS.values()),
        "condition_synonyms": {
            CONDITIONS[short]: list(labels) for short, labels in CONDITION_SYNONYMS.items()
        },
        "site_miscurations": SITE_MISCURATIONS,
        "hiv_contrasts_counted_in_sensitivities": {
            f"{experiment}": reason
            for (_, experiment), reason in {**HIV_RECLASSIFIED, **HIV_ON_ART}.items()
        },
        "definitions": {
            "profile": (
                "one study/experiment/direction aggregate of curated signature taxa (human hosts)"
            ),
            "feces": (
                "experiment body site exactly 'Feces'; multi-site labels such as "
                "'Feces,Blood' or 'Oral cavity,Feces' are grouped with Other non-GI"
            ),
            "condition": (
                "exact BugSigDB condition label of the experiment; condition_synonyms are "
                "counted only in the synonyms_merged sensitivity or with "
                "--condition-labels synonyms"
            ),
            "set": "matrix taxa annotated with the trait in at least one feces decreased profile",
            "annotation": (
                "a matrix taxon bears the trait when it or a descendant has the KG edge; "
                "has_direct_trait_edge marks the taxon's own edge"
            ),
            "annotation_share": (
                "trait count / sum of all retained non-MediaDive trait counts in the profile"
            ),
            "z": (
                "per-trait z-score of log1p(1000 * annotation share) over all 5,829 retained "
                "profiles"
            ),
            "paired_delta": "annotation share in the increased profile minus the decreased profile",
            "down": "the decreased signature has the larger annotation share (paired delta < 0)",
            "study_collapsed": (
                "one vote per study: the sign of the mean paired delta over its experiments; "
                "study_majority_* columns use the majority of experiment signs instead"
            ),
            "jaccard_distance": "1 - |intersection| / |union|",
            "monte_carlo_p": "(1 + #null at least as extreme) / (1 + permutations)",
            "bh_q_within_table": (
                "Benjamini-Hochberg over every p-value in significant_results_summary.tsv"
            ),
        },
        "universes": {
            "U1_all_cellulose_taxa_in_matrix": universes["U1"],
            "U2_all_cellulose_taxa_in_retained_feces_profiles": universes["U2"],
            "U3_all_cellulose_taxa_in_feces_decreased_profiles": universes["U3"],
            "U5_draws_from": "U3",
            "all_feces_decreased_taxa_PD_HIV_CRC_union": len(base_union),
            "U6_study_condition_units": len(units),
            "U6_studies": int(units["study_id"].nunique()),
            "U6_units_per_condition": unit_counts,
        },
        "observed": summary,
        "baseline_all_feces_decreased_taxa": base_summary,
        "pairwise_jaccard_significance": pairwise.to_dict("records"),
        "three_way_details": three_way,
        "null_models": nulls,
        "sensitivity_definitions": sensitivity.to_dict("records"),
        "compartment_stats": compartments.to_dict("records"),
        "compartment_direction_interaction": interaction.to_dict("records"),
        "profile_and_study_pair_jaccard": profile_pairs.to_dict("records"),
        "trait_benchmark_total_traits_tested": len(benchmark),
        "trait_benchmark_non_growth_medium_traits_tested": len(biological),
        "cellulose_rank_in_trait_benchmark": int(benchmark_row["rank_by_convergence"]),
        "cellulose_rank_excluding_growth_media": int(benchmark_row["rank_excluding_growth_media"]),
        "cellulose_percentile_excluding_growth_media": float(
            benchmark_row["percentile_excluding_growth_media"]
        ),
        "benchmark_median_mean_pairwise_jaccard_distance_excluding_growth_media": float(
            biological["mean_pairwise_jaccard_distance"].median()
        ),
        "benchmark_top_10_excluding_growth_media": biological.head(10)[
            [
                "trait_id",
                "trait_label",
                "trait_predicate",
                "mean_pairwise_jaccard_distance",
                "three_way_intersection",
                "three_way_union",
            ]
        ].to_dict("records"),
        "catalog_core_classes": catalog["core_class"].value_counts().to_dict(),
        "core_taxa_with_direct_trait_edge": int(
            catalog.loc[catalog["is_3way_core"], "has_direct_trait_edge"].sum()
        ),
        "direct_trait_edges_by_rank": direct_by_rank,
        "kg_direct_trait_subjects": len(direct),
        "condition_definition_overlap": scenario_overlap.to_dict("records"),
        "four_highlighted_pairs": {
            key: {k: v for k, v in value.items() if k != "top_12_taxa"}
            for key, value in four_pairs.items()
        },
        "label_source_comparison": comparison,
        "label_source_corrections_in_three_conditions": len(corrections),
        "target_condition_body_sites": target_sites.to_dict("records"),
        "excluded_multi_condition_feces_experiments": (
            multi_condition[["study_id", "experiment_id", "condition"]]
            .drop_duplicates()
            .to_dict("records")
        ),
        "condition_label_inventory": inventory.to_dict("records"),
        "monte_carlo": {"permutations": n, "seeds": rng_seeds},
    }

    outputs = {
        "pairwise_jaccard_significance.tsv": pairwise,
        "sensitivity_definitions_jaccard.tsv": sensitivity,
        "compartment_jaccard_and_direction_stats.tsv": compartments,
        "compartment_direction_interaction.tsv": interaction,
        "cellulose_decreased_taxa_catalog.tsv": catalog,
        "profile_and_study_pair_jaccard.tsv": profile_pairs,
        "trait_convergence_benchmark.tsv": benchmark,
        "significant_results_summary.tsv": significant,
        "four_highlighted_pairs_taxon_drivers.tsv": drivers,
        "four_highlighted_pairs_study_evidence.tsv": evidence,
        "label_source_corrections.tsv": corrections,
        "condition_definition_sensitivity.tsv": scenario_direction,
        "condition_definition_overlap.tsv": scenario_overlap,
        "condition_definition_interaction.tsv": scenario_interaction,
        "condition_label_inventory.tsv": inventory,
    }
    # Inputs are read many times above; refuse to write if any changed meanwhile.
    final_hashes = input_hashes(paths)
    changed = [key for key in hashes if final_hashes[key]["sha256"] != hashes[key]["sha256"]]
    if changed:
        raise RuntimeError(f"inputs changed during the run: {changed}; rerun once they are stable")
    for name, frame in outputs.items():
        write_tsv(frame, out / name)
    write_json(result, out / "cellulose_decr_summary.json")
    write_json(four_pairs, out / "four_highlighted_pairs_summary.json")

    modules = sorted(Path(__file__).resolve().parent.glob("*.py"))
    provenance = {
        "generated_at_utc": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "code_commit": args.code_commit,
        "script": "src/mtp/cellulose_decr/analysis.py",
        "bugsigdb_vintage": dump_vintage(paths["full_dump"]),
        "contrast_classifier_version": CLASSIFIER_VERSION,
        "label_source": args.label_source,
        "condition_labels": args.condition_labels,
        "permutations": n,
        "seed": args.seed,
        "seeds": rng_seeds,
        "scripts": {f"src/mtp/cellulose_decr/{path.name}": file_sha256(path) for path in modules},
        "inputs": hashes,
        "pinned_inputs": pinned_inputs,
        "inputs_drifted_from_pinned_provenance": drifted,
        "repinned": args.repin,
        "environment": {
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
        },
        "outputs": sorted(
            [*outputs, "cellulose_decr_summary.json", "four_highlighted_pairs_summary.json"]
        ),
    }
    write_json(provenance, out / "cellulose_decr_provenance.json")
    print(f"wrote {len(outputs) + 3} files to {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
