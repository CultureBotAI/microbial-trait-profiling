"""Filter and orient canonical derived tables to the 548-study contrast population.

When analysis and visualization scripts read the canonical unfiltered tables in
``data/derived/`` (such as ``study_trait_profiles_wide.tsv`` and ``bugsigdb_taxa.tsv``),
this helper restricts them to the canonical 548-study paired contrast population
(``disease_vs_healthy_tier1_2_subject_controls``) and reorients ``REVERSED`` contrasts
so that ``inc`` / ``increased`` always represents increased in disease versus control
(``direction_basis = "disease_vs_comparison"``).

When scripts are invoked by workflow arms or unit tests with pre-filtered or synthetic
input paths outside ``data/derived/``, ``is_canonical_derived_path`` returns ``False``
so those callers operate on their provided inputs directly.
"""

from __future__ import annotations

import functools
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from mtp import contrast_orientation as orientation

REPO_ROOT = Path(__file__).resolve().parents[2]
CANONICAL_DERIVED_DIR = (REPO_ROOT / "data/derived").resolve()
DEFAULT_RAW_DUMP_PATH = REPO_ROOT / "data/raw/full_dump.csv"
DEFAULT_PROFILES_PATH = REPO_ROOT / "data/derived/study_trait_profiles_wide.tsv"
DEFAULT_CONTRAST_POPULATION = "disease_vs_healthy_tier1_2_subject_controls"

ExperimentKey = tuple[str, str]


@dataclass(frozen=True)
class PopulationSelection:
    """Selected and reversed experiment keys for a contrast population."""

    contrast_population: str
    selected_keys: frozenset[ExperimentKey]
    reversed_keys: frozenset[ExperimentKey]
    human_corpus_keys: frozenset[ExperimentKey]


def is_canonical_derived_path(path: Path | str, filename: str) -> bool:
    """Return True iff ``path`` resolves to ``data/derived/<filename>`` in this repo."""
    try:
        resolved = Path(path).resolve()
    except OSError:
        return False
    return resolved == (CANONICAL_DERIVED_DIR / filename).resolve()


@functools.lru_cache(maxsize=8)
def get_population_selection(
    raw_dump_path: Path = DEFAULT_RAW_DUMP_PATH,
    profiles_path: Path = DEFAULT_PROFILES_PATH,
    contrast_population: str = DEFAULT_CONTRAST_POPULATION,
) -> PopulationSelection:
    """Compute the paired human 16S/WMS contrast selection and reversed experiment keys."""
    profiles = pd.read_csv(
        profiles_path,
        sep="\t",
        usecols=["study_id", "experiment_id", "inc_dec"],
        dtype=str,
    )
    directions_by_exp: dict[ExperimentKey, set[str]] = {}
    for study_id, experiment_id, inc_dec in profiles.itertuples(index=False, name=None):
        directions_by_exp.setdefault((str(study_id), str(experiment_id)), set()).add(str(inc_dec))
    human_corpus_keys = set(directions_by_exp.keys())
    paired_keys = {
        key for key, directions in directions_by_exp.items() if directions == {"inc", "dec"}
    }

    if contrast_population == "all_paired":
        return PopulationSelection(
            contrast_population=contrast_population,
            selected_keys=frozenset(paired_keys),
            reversed_keys=frozenset(),
            human_corpus_keys=frozenset(human_corpus_keys),
        )

    records = orientation.read_experiment_records(Path(raw_dump_path))
    classes = {key: orientation.classify(record) for key, record in records.items()}
    selected_keys: set[ExperimentKey] = set()
    reversed_keys: set[ExperimentKey] = set()
    for key in sorted(paired_keys):
        found = classes.get((str(key[0]), str(key[1])))
        if found is None or found.contrast_class != orientation.DISEASE_VS_HEALTHY:
            continue
        if (
            contrast_population != "disease_vs_healthy_oriented"
            and found.contrast_type != orientation.CASE_VS_CONTROL_SUBJECT
        ):
            continue
        if (
            contrast_population == "disease_vs_healthy_tier1_2_subject_controls"
            and found.control_tier
            not in {
                orientation.TIER_1_HEALTHY_NORMAL,
                orientation.TIER_2_CONDITION_NEGATIVE_SCREENED,
            }
        ):
            continue
        if found.orientation not in {orientation.NATIVE, orientation.REVERSED}:
            continue
        selected_keys.add((str(key[0]), str(key[1])))
        if found.orientation == orientation.REVERSED:
            reversed_keys.add((str(key[0]), str(key[1])))

    return PopulationSelection(
        contrast_population=contrast_population,
        selected_keys=frozenset(selected_keys),
        reversed_keys=frozenset(reversed_keys),
        human_corpus_keys=frozenset(human_corpus_keys),
    )


def orient_and_filter_profiles_df(
    profiles: pd.DataFrame,
    *,
    filter_to_selection: bool = True,
    raw_dump_path: Path = DEFAULT_RAW_DUMP_PATH,
    profiles_path: Path = DEFAULT_PROFILES_PATH,
    contrast_population: str = DEFAULT_CONTRAST_POPULATION,
) -> pd.DataFrame:
    """Filter a ``study_trait_profiles_wide`` DataFrame to 548 studies and orient ``inc_dec``."""
    selection = get_population_selection(
        raw_dump_path=raw_dump_path,
        profiles_path=profiles_path,
        contrast_population=contrast_population,
    )
    keys = list(
        zip(
            profiles["study_id"].astype(str),
            profiles["experiment_id"].astype(str),
            strict=False,
        )
    )
    if filter_to_selection:
        keep_mask = [key in selection.selected_keys for key in keys]
        filtered = profiles.loc[keep_mask].copy()
        keys = [key for key, keep in zip(keys, keep_mask, strict=False) if keep]
    else:
        filtered = profiles.copy()

    if selection.reversed_keys:
        flip_map = {"inc": "dec", "dec": "inc"}
        new_inc_dec = [
            flip_map.get(str(inc_dec), str(inc_dec))
            if key in selection.reversed_keys
            else str(inc_dec)
            for key, inc_dec in zip(keys, filtered["inc_dec"], strict=False)
        ]
        filtered["inc_dec"] = new_inc_dec
    return filtered


def orient_and_filter_taxa_df(
    taxa: pd.DataFrame,
    *,
    filter_to_selection: bool = True,
    preserve_non_corpus_rows: bool = False,
    raw_dump_path: Path = DEFAULT_RAW_DUMP_PATH,
    profiles_path: Path = DEFAULT_PROFILES_PATH,
    contrast_population: str = DEFAULT_CONTRAST_POPULATION,
) -> pd.DataFrame:
    """Filter a ``bugsigdb_taxa`` DataFrame to the 548-study cohort and orient ``direction``.

    When ``preserve_non_corpus_rows=True``, non-human and non-16S/WMS rows are kept in the
    table so downstream callers that audit excluded non-human/non-16S/WMS occurrences
    (via ``mtp.human_host_corpus.join_human_host_profiles``) continue to report them, while
    all human-host 16S/WMS rows are restricted to the 548-study selection and oriented.
    """
    selection = get_population_selection(
        raw_dump_path=raw_dump_path,
        profiles_path=profiles_path,
        contrast_population=contrast_population,
    )
    keys = list(
        zip(
            taxa["study_id"].astype(str),
            taxa["experiment_id"].astype(str),
            strict=False,
        )
    )
    if filter_to_selection:
        if preserve_non_corpus_rows:
            keep_mask = [
                (key in selection.selected_keys) if (key in selection.human_corpus_keys) else True
                for key in keys
            ]
        else:
            keep_mask = [key in selection.selected_keys for key in keys]
        filtered = taxa.loc[keep_mask].copy()
        keys = [key for key, keep in zip(keys, keep_mask, strict=False) if keep]
    else:
        filtered = taxa.copy()

    if selection.reversed_keys and "direction" in filtered.columns:
        flip_map = {
            "increased": "decreased",
            "decreased": "increased",
            "inc": "dec",
            "dec": "inc",
        }
        new_direction = [
            flip_map.get(str(direction), str(direction))
            if key in selection.reversed_keys
            else str(direction)
            for key, direction in zip(keys, filtered["direction"], strict=False)
        ]
        filtered["direction"] = new_direction
    return filtered
