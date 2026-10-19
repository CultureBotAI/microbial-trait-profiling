"""Reconcile declared raw, Study, and flattened sources independently (#651).

The checks read inputs directly; they do not rerun ingest or export writers.
Native leaf/rank parsing, metadata normalization, and logged exclusions are
explicit contracts, while repeated taxon observations retain multiplicity.
"""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping
from itertools import zip_longest
from pathlib import Path
from typing import Any

import yaml
from pydantic import BaseModel

from mtp.datamodel import Study
from mtp.datamodel.core_products_pydantic import (
    ExperimentOntologyTermRow,
    ExperimentTableRow,
    SignatureTableRow,
    StudyTableRow,
    TaxonOccurrenceRow,
)
from mtp.ingest import METAPHLAN_RANK_MAP

_EXPERIMENT_TEXT = {
    "experiment_name": "Experiment",
    "location_of_subjects": "Location of subjects",
    "host_species": "Host species",
    "body_site": "Body site",
    "condition": "Condition",
    "group_0_name": "Group 0 name",
    "group_1_name": "Group 1 name",
    "group_1_definition": "Group 1 definition",
    "antibiotics_exclusion": "Antibiotics exclusion",
    "sequencing_type": "Sequencing type",
    "variable_region_16s": "16S variable region",
    "sequencing_platform": "Sequencing platform",
    "data_transformation": "Data transformation",
    "statistical_test": "Statistical test",
    "mht_correction": "MHT correction",
    "matched_on": "Matched on",
    "confounders_controlled_for": "Confounders controlled for",
    "pielou": "Pielou",
    "shannon": "Shannon",
    "chao1": "Chao1",
    "simpson": "Simpson",
    "inverse_simpson": "Inverse Simpson",
    "richness": "Richness",
}
_STUDY_TEXT = {
    "title": "Title",
    "study_design": "Study design",
    "authors_list": "Authors list",
    "journal": "Journal",
    "keywords": "Keywords",
}
_TABLE_MODELS: dict[str, type[BaseModel]] = {
    "study_table": StudyTableRow,
    "experiments": ExperimentTableRow,
    "experiment_ontology_terms": ExperimentOntologyTermRow,
    "signatures": SignatureTableRow,
    "reported_taxa": TaxonOccurrenceRow,
}


def _clean(value: str) -> str | None:
    stripped = value.strip()
    return None if stripped in {"", "NA"} else stripped


def _integer(value: str) -> int | None:
    return int(value.strip()) if value.strip().isdigit() else None


def _floating(value: str) -> float | None:
    try:
        return float(value.strip())
    except ValueError:
        return None


def _skip(events: list[dict[str, str]], entity: str, reason: str, **identity: str) -> None:
    events.append({"level": "skip", "entity": entity, "reason": reason, **identity})


def _taxon_name(name: str) -> None:
    if not name.strip() or name.strip().startswith("taxon_"):
        raise ValueError(f"source reconciliation: blank or placeholder taxon name {name!r}")


def _source_taxa(
    row: Mapping[str, str], signature: str, events: list[dict[str, str]]
) -> list[dict[str, str]]:
    ids = row.get("NCBI Taxonomy IDs", "")
    if ids.strip() in {"", "NA"}:
        return []
    chains = [value.strip() for value in ids.split(";") if value.strip()]
    names = row.get("MetaPhlAn taxon names", "") or row.get("Taxon", "")
    name_chains = [value.strip() for value in names.split(",") if value.strip()]
    taxa: list[dict[str, str]] = []
    for index, chain in enumerate(chains):
        lineage = [value.strip() for value in chain.split("|") if value.strip()]
        if not lineage:
            continue
        leaf = lineage[-1]
        if leaf.isdigit():
            taxon_id = f"NCBITaxon:{leaf}"
        elif leaf.startswith("NCBITaxon:"):
            taxon_id = leaf
        else:
            _skip(
                events,
                "taxon",
                f"unparsable leaf taxon ID: {leaf!r}",
                id=leaf,
                signature=signature,
            )
            continue
        if index >= len(name_chains):
            raise ValueError(f"source reconciliation: missing taxon name for {signature}/{leaf}")
        name_parts = [value.strip() for value in name_chains[index].split("|") if value.strip()]
        name = name_parts[-1] if name_parts else ""
        rank = None
        for prefix, candidate in METAPHLAN_RANK_MAP.items():
            if name.startswith(prefix):
                name, rank = name[len(prefix) :], candidate
                break
        _taxon_name(name)
        if rank is None:
            words = len(name.split())
            rank = "strain" if words >= 3 else "species" if words == 2 else "genus"
        taxa.append({"id": taxon_id, "taxon_name": name, "taxonomic_rank": rank})
    return taxa


def _experiment_metadata(row: Mapping[str, str], identifier: str) -> dict[str, Any]:
    result: dict[str, Any] = {"id": identifier}
    result.update({name: _clean(row.get(column, "")) for name, column in _EXPERIMENT_TEXT.items()})
    for name, column in (
        ("body_site_ontology_id", "UBERON ID"),
        ("condition_ontology_id", "EFO ID"),
    ):
        result[name] = [
            value.strip() for value in row.get(column, "").split(",") if _clean(value)
        ] or None
    for group in (0, 1):
        result[f"group_{group}_sample_size"] = _integer(row.get(f"Group {group} sample size", ""))
    result["significance_threshold"] = _floating(row.get("Significance threshold", ""))
    result["lda_score_above"] = _floating(row.get("LDA Score above", ""))
    return result


def _study_metadata(row: Mapping[str, str], identifier: str) -> dict[str, Any]:
    doi = (
        row.get("DOI", "").strip().removeprefix("https://doi.org/").removeprefix("http://doi.org/")
    )
    url = row.get("URL", "").strip()
    return {
        "id": identifier,
        "source_record_id": row.get("Study", "").strip(),
        "pmid": _integer(row.get("PMID", "")),
        "publication_year": _integer(row.get("Year", "")),
        "doi": doi if doi.startswith("10.") else None,
        "url": url if url.startswith("http") else None,
        **{name: _clean(row.get(column, "")) for name, column in _STUDY_TEXT.items()},
    }


def _expected_studies(raw: Path) -> tuple[dict[str, dict[str, Any]], list[dict[str, str]], int]:
    events: list[dict[str, str]] = []
    grouped: dict[str, dict[str, list[dict[str, str]]]] = defaultdict(lambda: defaultdict(list))
    first_study_rows: dict[str, dict[str, str]] = {}
    with raw.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(line for line in handle if not line.startswith("#")))
    for row in rows:
        raw_id = row.get("BSDB ID", "").strip()
        if not raw_id:
            _skip(events, "row", "empty BSDB ID")
            continue
        parts = raw_id.removeprefix("bsdb:").split("/")
        study_id = f"bsdb:{parts[0]}"
        experiment_id = f"bsdb:{parts[0]}-{parts[1]}" if len(parts) >= 2 else raw_id
        grouped[study_id][experiment_id].append(row)
        first_study_rows.setdefault(study_id, row)
    expected: dict[str, dict[str, Any]] = {}
    for study_id, experiments in sorted(grouped.items()):
        experiment_records = []
        for experiment_id, source_rows in sorted(experiments.items()):
            signatures = []
            for row in source_rows:
                signature_id = row["BSDB ID"].replace("/", "-")
                direction = row.get("Abundance in Group 1", "").strip().lower()
                if direction not in {"increased", "decreased"}:
                    _skip(
                        events, "signature", f"unparsable direction: {direction!r}", id=signature_id
                    )
                    continue
                taxa = _source_taxa(row, signature_id, events)
                if not taxa:
                    _skip(events, "signature", "no parseable taxa", id=signature_id)
                    continue
                signatures.append(
                    {
                        "id": signature_id,
                        "direction": direction,
                        "signature_source": _clean(row.get("Source", "")),
                        "signature_description": _clean(row.get("Description", "")),
                        "taxa": taxa,
                    }
                )
            if not signatures:
                _skip(events, "experiment", "no valid signatures", id=experiment_id)
                continue
            experiment_records.append(
                {**_experiment_metadata(source_rows[0], experiment_id), "signatures": signatures}
            )
        if not experiment_records:
            _skip(events, "study", "no valid experiments/signatures", id=study_id)
            continue
        record = {
            **_study_metadata(first_study_rows[study_id], study_id),
            "experiments": experiment_records,
        }
        expected[study_id] = Study.model_validate(record).model_dump(mode="json")
    return expected, events, len(rows)


def _study_records(studies: Path) -> list[tuple[Path, Study]]:
    paths = sorted(studies.glob("*.yaml"))
    if not paths or set(studies.iterdir()) != set(paths):
        raise ValueError("Study collection must contain only its nonempty YAML membership")
    records: list[tuple[Path, Study]] = []
    seen: set[str] = set()
    for path in paths:
        study = Study.model_validate(yaml.safe_load(path.read_text(encoding="utf-8")))
        if study.id in seen:
            raise ValueError(f"Duplicate study ID in collection: {study.id}")
        seen.add(study.id)
        for experiment in study.experiments or []:
            for signature in experiment.signatures or []:
                for taxon in signature.taxa or []:
                    _taxon_name(taxon.taxon_name or "")
        records.append((path, study))
    return records


def validate_study_sources(raw: Path, studies: Path, ingest_log: Path) -> dict[str, object]:
    """Require exact raw-derived records and independently justified skip events."""
    expected, expected_events, raw_rows = _expected_studies(raw)
    actual = {study.id: study.model_dump(mode="json") for _, study in _study_records(studies)}
    if set(actual) != set(expected):
        raise ValueError("source reconciliation: Study membership differs from raw rows")
    for identifier in expected:
        if actual[identifier] != expected[identifier]:
            raise ValueError(
                f"source reconciliation: {identifier} differs from normalized raw rows"
            )
    events = json.loads(ingest_log.read_text(encoding="utf-8"))
    if not isinstance(events, list) or any(
        not isinstance(event, dict)
        or any(
            not isinstance(key, str) or not isinstance(value, str) for key, value in event.items()
        )
        for event in events
    ):
        raise ValueError("source reconciliation: ingest log must contain string event mappings")
    observed_skips = Counter(
        tuple(sorted(event.items())) for event in events if event.get("level") == "skip"
    )
    if observed_skips != Counter(tuple(sorted(event.items())) for event in expected_events):
        raise ValueError("source reconciliation: ingest skip events do not match raw exclusions")
    return {
        "studies": len(actual),
        "raw_rows": raw_rows,
        "reconciled_skip_events": len(expected_events),
    }


def _text(value: object) -> str:
    return "" if value is None else str(value)


def _export_rows(
    studies: Path,
) -> tuple[dict[str, list[dict[str, object]]], dict[tuple[str, str, str], Counter[str]]]:
    tables: dict[str, list[dict[str, object]]] = {name: [] for name in _TABLE_MODELS}
    profiles: dict[tuple[str, str, str], Counter[str]] = defaultdict(Counter)
    for path, study in _study_records(studies):
        study_signatures = study_taxa = 0
        for experiment in study.experiments or []:
            experiment_signatures = experiment_taxa = 0
            identifiers = {"study_id": study.id, "experiment_id": experiment.id}
            for signature in experiment.signatures or []:
                taxa = signature.taxa or []
                experiment_signatures += 1
                experiment_taxa += len(taxa)
                child_ids = {
                    **identifiers,
                    "signature_id": signature.id,
                    "direction": signature.direction,
                }
                tables["signatures"].append({**child_ids, "taxon_count": len(taxa)})
                direction = "inc" if signature.direction == "increased" else "dec"
                for taxon in taxa:
                    tables["reported_taxa"].append(
                        {
                            **child_ids,
                            "taxon_id": taxon.id,
                            "taxon_name": taxon.taxon_name,
                            "taxonomic_rank": taxon.taxonomic_rank,
                        }
                    )
                    profiles[(study.id, experiment.id, direction)][
                        (taxon.taxon_name or "").strip()
                    ] += 1
            study_signatures += experiment_signatures
            study_taxa += experiment_taxa
            tables["experiments"].append(
                {
                    **experiment.model_dump(mode="json"),
                    **identifiers,
                    "signature_count": experiment_signatures,
                    "taxon_count": experiment_taxa,
                }
            )
            for slot in ("condition_ontology_id", "body_site_ontology_id"):
                for term in getattr(experiment, slot) or []:
                    if term.strip():
                        tables["experiment_ontology_terms"].append(
                            {**identifiers, "slot": slot, "term_id": term.strip()}
                        )
        tables["study_table"].append(
            {
                **study.model_dump(mode="json"),
                "study_id": study.id,
                "experiment_count": len(study.experiments or []),
                "signature_count": study_signatures,
                "taxon_count": study_taxa,
                "input_file": path.name,
            }
        )
    return tables, profiles


def _compare_rows(path: Path, header: list[str], expected: Iterable[list[str]]) -> None:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        if next(reader, []) != header:
            raise ValueError(
                f"source reconciliation: {path.name} differs from Study source records"
            )
        for number, (actual, source) in enumerate(zip_longest(reader, expected), start=2):
            if actual != source:
                raise ValueError(
                    f"source reconciliation: {path.name}:{number} differs from Study source records"
                )


def validate_export_sources(studies: Path, paths: Mapping[str, Path]) -> dict[str, object]:
    """Check all seven sibling exports once, including metadata and name profiles."""
    required = {*_TABLE_MODELS, "taxon_profiles_wide", "taxon_profiles_long"}
    if not required.issubset(paths):
        raise ValueError("source reconciliation: export validation requires every sibling product")
    tables, profiles = _export_rows(studies)
    for name, model in _TABLE_MODELS.items():
        header = list(model.model_fields)
        expected = ([_text(row.get(column)) for column in header] for row in tables[name])
        _compare_rows(paths[name], header, expected)
    keys = sorted(profiles)
    features = sorted({feature for counts in profiles.values() for feature in counts})
    _compare_rows(
        paths["taxon_profiles_wide"],
        ["study_id", "experiment_id", "inc_dec", *features],
        ([*key, *[str(profiles[key][feature]) for feature in features]] for key in keys),
    )
    _compare_rows(
        paths["taxon_profiles_long"],
        ["taxon", *["_".join(key) for key in keys]],
        ([feature, *[str(profiles[key][feature]) for key in keys]] for feature in features),
    )
    return {
        "source_reconciled_exports": len(required),
        "source_taxon_occurrences": len(tables["reported_taxa"]),
    }
