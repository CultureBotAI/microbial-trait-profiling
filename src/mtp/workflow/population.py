"""Project an ingested corpus onto declared contrast populations (#658).

The raw-derived source collection stays immutable. Selection precedes trait
annotation and never filters on annotation totals. Direction relabelling uses
the existing classifier API; group metadata and taxon occurrences stay intact.
"""

from __future__ import annotations

import copy
import json
import shutil
from collections import Counter
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import yaml

from mtp import contrast_orientation as orientation
from mtp.datamodel import Study
from mtp.paths import file_sha256
from mtp.trait_guilds import HUMAN_HOST_SPECIES, Profile

SCHEMA_VERSION = "contrast-selection-v1"
ALL_PAIRED = "all_paired"
POPULATIONS = (
    ALL_PAIRED,
    "disease_vs_healthy_oriented",
    "case_vs_control_subject",
    "disease_vs_healthy_tier1_2_subject_controls",
)
NATIVE_DIRECTION_BASIS = "group_1_vs_group_0"
ORIENTED_DIRECTION_BASIS = "disease_vs_comparison"
ExperimentKey = tuple[str, str]


def _records(directory: Path) -> dict[str, dict[str, Any]]:
    if directory.is_symlink() or not directory.is_dir():
        raise ValueError("Population Study collection must be a regular directory")
    members = sorted(directory.iterdir())
    if not members:
        raise ValueError("Population Study collection must not be empty")
    records: dict[str, dict[str, Any]] = {}
    identifiers: set[str] = set()
    experiments: set[ExperimentKey] = set()
    for path in members:
        if path.is_symlink() or not path.is_file() or path.suffix != ".yaml":
            raise ValueError(f"Unexpected member of population Study collection: {path}")
        value = yaml.safe_load(path.read_text(encoding="utf-8"))
        study = Study.model_validate(value)
        if study.id in identifiers:
            raise ValueError(f"Duplicate population study ID: {study.id}")
        if not study.experiments:
            raise ValueError(f"Population source study has no experiments: {study.id}")
        identifiers.add(study.id)
        for experiment in study.experiments:
            key = (study.id, experiment.id)
            if key in experiments:
                raise ValueError(f"Duplicate population experiment ID: {key}")
            experiments.add(key)
        records[path.name] = value
    return records


def _ids(keys: list[ExperimentKey]) -> list[dict[str, str]]:
    return [{"study_id": study, "experiment_id": experiment} for study, experiment in sorted(keys)]


def _inventory(directory: Path) -> list[dict[str, str]]:
    return [
        {"name": path.name, "sha256": file_sha256(path)} for path in sorted(directory.iterdir())
    ]


def _counts(records: Mapping[str, dict[str, Any]]) -> dict[str, int]:
    experiments = [experiment for study in records.values() for experiment in study["experiments"]]
    signatures = [signature for experiment in experiments for signature in experiment["signatures"]]
    return {
        "studies": len(records),
        "experiments": len(experiments),
        "signatures": len(signatures),
        "taxon_occurrences": sum(len(signature["taxa"]) for signature in signatures),
    }


def _reason(
    population: str,
    human: bool,
    paired: bool,
    found: orientation.Classification | None,
) -> str:
    if population == ALL_PAIRED:
        return "identity"
    if not human:
        return "non_human"
    if not paired:
        return "unpaired"
    if found is None:
        raise ValueError("Population experiment has no classification")
    if found.contrast_class != orientation.DISEASE_VS_HEALTHY:
        return "contrast_class"
    if population != "disease_vs_healthy_oriented":
        if found.contrast_type != orientation.CASE_VS_CONTROL_SUBJECT:
            return "contrast_type"
        if (
            population == "disease_vs_healthy_tier1_2_subject_controls"
            and found.control_tier
            not in {
                orientation.TIER_1_HEALTHY_NORMAL,
                orientation.TIER_2_CONDITION_NEGATIVE_SCREENED,
            }
        ):
            return "control_tier"
    if found.orientation not in {orientation.NATIVE, orientation.REVERSED}:
        raise ValueError("Selected population experiment has an invalid orientation")
    return "selected"


def _direction_map(
    study_id: str, experiment: Mapping[str, Any], found: orientation.Classification
) -> dict[str, str]:
    # These probes contain direction identity only. There are no trait counts at
    # this stage; oriented_profiles neither reads nor filters their empty totals.
    probes = [
        Profile(
            study_id,
            str(experiment["id"]),
            direction,
            str(experiment.get("condition") or ""),
            str(experiment.get("body_site") or ""),
            {},
            0,
        )
        for direction in ("inc", "dec")
    ]
    relabelled = orientation.oriented_profiles(probes, {(study_id, str(experiment["id"])): found})
    if len(relabelled) != len(probes):
        raise ValueError("Classifier unexpectedly dropped a selected paired experiment")
    full_direction = {"inc": "increased", "dec": "decreased"}
    return {
        full_direction[before.inc_dec]: full_direction[after.inc_dec]
        for before, after in zip(probes, relabelled, strict=True)
    }


def _plan(
    source: Mapping[str, dict[str, Any]], raw_dump: Path, population: str
) -> tuple[dict[str, dict[str, Any]], dict[str, object]]:
    if population not in POPULATIONS:
        raise ValueError(f"Unknown contrast population: {population!r}")
    # Identity mode intentionally adds no classifier-specific CSV requirements.
    records = orientation.read_experiment_records(raw_dump) if population != ALL_PAIRED else {}
    classes = {key: orientation.classify(record) for key, record in records.items()}
    selected: dict[str, dict[str, Any]] = {}
    source_keys: list[ExperimentKey] = []
    selected_keys: list[ExperimentKey] = []
    reversed_keys: list[ExperimentKey] = []
    decisions: list[dict[str, object]] = []
    for name, study in source.items():
        kept_experiments = []
        for experiment in study["experiments"]:
            key = (str(study["id"]), str(experiment["id"]))
            source_keys.append(key)
            found = classes.get(key)
            if population != ALL_PAIRED and found is None:
                raise ValueError(f"Population experiment has no raw export record: {key}")
            directions = {str(signature["direction"]) for signature in experiment["signatures"]}
            paired = directions == {"increased", "decreased"}
            human = str(experiment.get("host_species") or "").strip() == HUMAN_HOST_SPECIES
            reason = _reason(population, human, paired, found)
            keep = reason in {"identity", "selected"}
            decisions.append(
                {
                    "study_id": key[0],
                    "experiment_id": key[1],
                    "host_species": experiment.get("host_species"),
                    "paired": paired,
                    "source_directions": sorted(directions),
                    "contrast_class": found.contrast_class if found else None,
                    "contrast_type": found.contrast_type if found else None,
                    "control_tier": found.control_tier if found else None,
                    "orientation": found.orientation if found else None,
                    "evidence": found.evidence if found else None,
                    "flags": found.flags if found else None,
                    "selected": keep,
                    "exclusion_reason": None if keep else reason,
                }
            )
            if not keep:
                continue
            selected_keys.append(key)
            projected = copy.deepcopy(experiment)
            if population != ALL_PAIRED:
                if found is None:
                    raise ValueError(f"Selected population experiment lacks classification: {key}")
                mapping = _direction_map(key[0], experiment, found)
                for signature in projected["signatures"]:
                    signature["direction"] = mapping[signature["direction"]]
                if found.orientation == orientation.REVERSED:
                    reversed_keys.append(key)
            kept_experiments.append(projected)
        if kept_experiments:
            projected_study = copy.deepcopy(study)
            projected_study["experiments"] = kept_experiments
            selected[name] = projected_study
    if not selected:
        raise ValueError(f"Contrast population {population!r} selects no experiments")
    audit: dict[str, object] = {
        "schema_version": SCHEMA_VERSION,
        "population": population,
        "input_direction_basis": NATIVE_DIRECTION_BASIS,
        "output_direction_basis": (
            NATIVE_DIRECTION_BASIS if population == ALL_PAIRED else ORIENTED_DIRECTION_BASIS
        ),
        "classifier": {
            "applied": population != ALL_PAIRED,
            "version": orientation.CLASSIFIER_VERSION,
            "source_sha256": file_sha256(Path(orientation.__file__)),
        },
        "raw_dump_sha256": file_sha256(raw_dump),
        "source_experiments": _ids(source_keys),
        "selected_experiments": _ids(selected_keys),
        "reversed_experiments": _ids(reversed_keys),
        "policy": {
            "human_only": population != ALL_PAIRED,
            "paired_only": population != ALL_PAIRED,
            "trait_total_filter": "none",
            "strong_evidence_only": False,
            "treated_disease_flag_filter": False,
            "preserve_group_metadata": True,
            "preserve_taxon_occurrences": True,
        },
        "counts": {
            "source": _counts(source),
            "selected": _counts(selected),
            "reversed_experiments": len(reversed_keys),
            "excluded_experiments": len(source_keys) - len(selected_keys),
            "exclusions_by_reason": dict(
                sorted(
                    Counter(
                        str(row["exclusion_reason"]) for row in decisions if not row["selected"]
                    ).items()
                )
            ),
        },
        "experiments": sorted(
            decisions, key=lambda row: (str(row["study_id"]), str(row["experiment_id"]))
        ),
    }
    return selected, audit


def _preflight(source: Path, raw: Path, selected: Path, audit: Path) -> None:
    if not raw.is_file():
        raise ValueError("Population raw dump must be a file")
    if audit.exists() or audit.is_symlink():
        raise ValueError("Population audit must be a new file")
    if selected.is_symlink() or (
        selected.exists() and (not selected.is_dir() or any(selected.iterdir()))
    ):
        raise ValueError("Selected population collection must be an empty regular directory")
    source_path, selected_path, audit_path = source.resolve(), selected.resolve(), audit.resolve()
    if (
        source_path == selected_path
        or source_path in selected_path.parents
        or selected_path in source_path.parents
        or source_path in audit_path.parents
        or selected_path in audit_path.parents
        or selected_path == audit_path
        or raw.resolve() == selected_path
    ):
        raise ValueError("Population outputs must not overlap their inputs or each other")


def select_contrasts(
    inputs: dict[str, Path], outputs: dict[str, Path], params: dict[str, object]
) -> None:
    """Write a source-preserving population projection and its auditable decision record."""
    if set(inputs) != {"source_studies", "raw_dump"}:
        raise ValueError("Population selection requires source_studies and raw_dump")
    if set(outputs) != {"studies", "contrast_selection"}:
        raise ValueError("Population selection requires studies and contrast_selection outputs")
    if set(params) != {"contrast_population"} or not isinstance(params["contrast_population"], str):
        raise ValueError("Population selection requires exactly one contrast_population string")
    source, raw = inputs["source_studies"], inputs["raw_dump"]
    target, audit_path = outputs["studies"], outputs["contrast_selection"]
    _preflight(source, raw, target, audit_path)
    selected, audit = _plan(_records(source), raw, params["contrast_population"])
    audit["source_files"] = _inventory(source)
    target.mkdir(parents=True, exist_ok=True)
    for name, record in selected.items():
        if params["contrast_population"] == ALL_PAIRED:
            shutil.copyfile(source / name, target / name)
        else:
            (target / name).write_text(
                yaml.safe_dump(record, sort_keys=False, allow_unicode=True), encoding="utf-8"
            )
    audit["selected_files"] = _inventory(target)
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    audit_path.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def validate_population(
    studies: Path,
    selection_path: Path,
    source_studies: Path,
    raw_dump: Path,
    population: str,
) -> dict[str, object]:
    """Reject any unrecorded selection, metadata/taxon change, or direction transformation."""
    expected, audit = _plan(_records(source_studies), raw_dump, population)
    actual = _records(studies)
    if json.dumps(actual, sort_keys=True, allow_nan=False) != json.dumps(
        expected, sort_keys=True, allow_nan=False
    ):
        raise ValueError("Population source reconciliation: selected Study projection differs")
    audit["source_files"] = _inventory(source_studies)
    audit["selected_files"] = _inventory(studies)
    if population == ALL_PAIRED and audit["source_files"] != audit["selected_files"]:
        raise ValueError("Population identity projection must preserve source YAML bytes")
    recorded_audit = json.loads(selection_path.read_text(encoding="utf-8"))
    if json.dumps(recorded_audit, sort_keys=True, allow_nan=False) != json.dumps(
        audit, sort_keys=True, allow_nan=False
    ):
        raise ValueError("Population selection audit differs from its source-derived contract")
    return {
        "population": population,
        "direction_basis": audit["output_direction_basis"],
        "counts": audit["counts"],
        "classifier": audit["classifier"],
    }
