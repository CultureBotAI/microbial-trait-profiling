"""BugSigDB CSV to per-study YAML ingestion."""

from __future__ import annotations

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import yaml

from mtp.datamodel import Study
from mtp.paths import DERIVED_DATA_DIR, RAW_DATA_DIR, STUDY_DATA_DIR

REQUIRED_EXPORTS = ("studies.csv", "experiments.csv", "signatures.csv")


def assert_required_exports(raw_dir: Path) -> None:
    """Fail fast when the expected BugSigDB CSV exports are missing."""
    missing = [name for name in REQUIRED_EXPORTS if not (raw_dir / name).exists()]
    if missing:
        missing_display = ", ".join(sorted(missing))
        raise FileNotFoundError(f"Missing required BugSigDB exports: {missing_display}")


def _read_csv(path: Path) -> list[dict[str, str]]:
    """Read a CSV, skipping comment lines that start with ``#``."""
    with path.open(encoding="utf-8", newline="") as handle:
        lines = (line for line in handle if not line.startswith("#"))
        return list(csv.DictReader(lines))


def _parse_study_id(bsdb_id: str) -> str:
    """Extract the numeric study id from a BSDB ID like ``bsdb:11/1/1``.

    The first number after ``bsdb:`` is the study number.
    """
    raw = bsdb_id.removeprefix("bsdb:").split("/")[0]
    return f"bsdb:{raw}"


def _parse_experiment_id(bsdb_id: str) -> str:
    """Extract study/experiment id from ``bsdb:11/1/1`` → ``bsdb:11-1``."""
    parts = bsdb_id.removeprefix("bsdb:").split("/")
    return f"bsdb:{parts[0]}-{parts[1]}" if len(parts) >= 2 else bsdb_id


def _parse_signature_id(bsdb_id: str) -> str:
    """Convert ``bsdb:11/1/1`` → ``bsdb:11-1-1`` (CURIE-safe delimiters)."""
    return bsdb_id.replace("/", "-")


def _clean(value: str) -> str | None:
    """Strip whitespace and convert 'NA' or empty strings to None."""
    v = value.strip()
    return None if v in ("", "NA") else v


def _clean_list(value: str) -> list[str] | None:
    """Split a comma-separated value into a list, filtering blanks and NA."""
    items = [v.strip() for v in value.split(",") if v.strip() not in ("", "NA")]
    return items if items else None


def _clean_int(value: str) -> int | None:
    """Parse an integer from a CSV value, returning None for non-numeric."""
    v = value.strip()
    return int(v) if v.isdigit() else None


def _clean_float(value: str) -> float | None:
    """Parse a float from a CSV value, returning None for non-numeric."""
    v = value.strip()
    if v in ("", "NA"):
        return None
    try:
        return float(v)
    except ValueError:
        return None


def _taxon_name_to_rank(name: str) -> str:
    """Heuristic for taxonomic rank from the taxon name string.

    BugSigDB signatures use names at varying ranks. This is a rough
    heuristic — proper resolution belongs in issue #29.
    """
    parts = name.strip().split()
    if len(parts) >= 3:
        return "strain"
    if len(parts) == 2:
        return "species"
    return "genus"


def ingest(
    raw_dir: Path | None = None,
    study_dir: Path | None = None,
    derived_dir: Path | None = None,
    *,
    full_dump_path: Path | None = None,
    log_path: Path | None = None,
) -> list[Path]:
    """Parse raw CSVs and emit per-study YAML files.

    Each YAML file is a bare Study validated through the Pydantic model
    before writing. Writes an ingest log to derived_dir/ingest_log.json
    documenting every skipped row, signature, experiment, and study.
    Returns the list of written file paths.
    """
    raw_dir = raw_dir or RAW_DATA_DIR
    study_dir = study_dir or STUDY_DATA_DIR
    derived_dir = derived_dir or DERIVED_DATA_DIR
    study_dir.mkdir(parents=True, exist_ok=True)
    derived_dir.mkdir(parents=True, exist_ok=True)

    log: list[dict[str, str]] = []

    full_dump = full_dump_path if full_dump_path is not None else raw_dir / "full_dump.csv"
    if full_dump_path is not None and not full_dump.is_file():
        raise FileNotFoundError(f"BugSigDB full dump not found: {full_dump}")
    if full_dump.exists():
        rows = _read_csv(full_dump)
    else:
        assert_required_exports(raw_dir)
        rows = _build_rows_from_separate_csvs(raw_dir)

    # Group rows by study id
    studies: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        bsdb_id = row.get("BSDB ID", "").strip()
        if not bsdb_id:
            log.append({"level": "skip", "entity": "row", "reason": "empty BSDB ID"})
            continue
        sid = _parse_study_id(bsdb_id)
        studies[sid].append(row)

    written: list[Path] = []
    for study_id, study_rows in sorted(studies.items()):
        record = _build_study_record(study_id, study_rows, log)
        if not record["experiments"]:
            log.append(
                {
                    "level": "skip",
                    "entity": "study",
                    "id": study_id,
                    "reason": "no valid experiments/signatures",
                }
            )
            print(f"skipped {study_id}: no valid experiments/signatures", file=sys.stderr)
            continue

        # Validate through Pydantic before writing
        Study.model_validate(record)

        filename = study_id.replace(":", "_").replace("/", "_") + ".yaml"
        dest = study_dir / filename
        dest.write_text(yaml.dump(record, default_flow_style=False, sort_keys=False))
        written.append(dest)

    # Write ingest log
    log_path = log_path if log_path is not None else derived_dir / "ingest_log.json"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(json.dumps(log, indent=2) + "\n")

    # Print summary and check for systemic problems
    from collections import Counter

    skips = [e for e in log if e["level"] == "skip"]
    infos = [e for e in log if e["level"] == "info"]
    skip_counts = Counter(e["reason"] for e in skips)
    info_counts = Counter(e["reason"] for e in infos)
    print(
        f"wrote {len(written)} studies, {len(skips)} skips, {len(infos)} info entries:",
        file=sys.stderr,
    )
    for reason, count in skip_counts.most_common():
        print(f"  {count:>5}  {reason}", file=sys.stderr)

    # Fail if ingest produced any placeholder-related info entries.
    # Placeholders are treated as evidence of a systemic parsing problem,
    # such as wrong delimiters or missing lookups, rather than an edge case.
    placeholder_reasons = Counter(
        {reason: count for reason, count in info_counts.items() if "placeholder" in reason}
    )
    placeholder_count = sum(placeholder_reasons.values())
    if placeholder_count > 0:
        most_common = placeholder_reasons.most_common(1)[0][0]
        msg = (
            f"Ingest produced {placeholder_count} placeholder taxon names. "
            f"Most common reason: {most_common}. "
            f"This indicates a systemic data parsing problem, not an edge case."
        )
        raise ValueError(msg)

    return written


def _build_rows_from_separate_csvs(raw_dir: Path) -> list[dict[str, str]]:
    """Merge the three separate CSV exports into a flat row list.

    This mimics what BugSigDBExports does in dump_release.R but in Python.
    For now we only use full_dump.csv — this is a stub for #27 raw CSV support.
    """
    msg = "Separate CSV merging not yet implemented. Please use full_dump.csv from BugSigDBExports."
    raise NotImplementedError(msg)


def _build_study_record(
    study_id: str,
    rows: list[dict[str, str]],
    log: list[dict[str, str]],
) -> dict[str, object]:
    """Build a nested study record dict from flat rows sharing a study id."""
    first = rows[0]

    # Group by experiment
    experiments: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        eid = _parse_experiment_id(row.get("BSDB ID", ""))
        experiments[eid].append(row)

    experiment_records: list[dict[str, object]] = []
    for exp_id, exp_rows in sorted(experiments.items()):
        exp_first = exp_rows[0]
        sig_records: list[dict[str, object]] = []
        for row in exp_rows:
            sig_id = _parse_signature_id(row.get("BSDB ID", ""))
            taxa_str = row.get("NCBI Taxonomy IDs", "")
            names_str = row.get("MetaPhlAn taxon names", "") or row.get("Taxon", "") or ""
            direction_raw = (row.get("Abundance in Group 1", "") or "").strip().lower()
            if direction_raw not in ("increased", "decreased"):
                log.append(
                    {
                        "level": "skip",
                        "entity": "signature",
                        "id": sig_id,
                        "reason": f"unparsable direction: {direction_raw!r}",
                    }
                )
                continue

            taxa = _parse_taxa(taxa_str, names_str, sig_id, log)
            if not taxa:
                log.append(
                    {
                        "level": "skip",
                        "entity": "signature",
                        "id": sig_id,
                        "reason": "no parseable taxa",
                    }
                )
                continue

            sig_records.append(
                {
                    "id": sig_id,
                    "direction": direction_raw,
                    "signature_source": _clean(row.get("Source", "")),
                    "signature_description": _clean(row.get("Description", "")),
                    "taxa": taxa,
                }
            )

        if not sig_records:
            log.append(
                {
                    "level": "skip",
                    "entity": "experiment",
                    "id": exp_id,
                    "reason": "no valid signatures",
                }
            )
            continue

        experiment_records.append(
            {
                "id": exp_id,
                "experiment_name": _clean(exp_first.get("Experiment", "")),
                "location_of_subjects": _clean(exp_first.get("Location of subjects", "")),
                "host_species": _clean(exp_first.get("Host species", "")),
                "body_site": _clean(exp_first.get("Body site", "")),
                "body_site_ontology_id": _clean_list(exp_first.get("UBERON ID", "")),
                "condition": _clean(exp_first.get("Condition", "")),
                "condition_ontology_id": _clean_list(exp_first.get("EFO ID", "")),
                "group_0_name": _clean(exp_first.get("Group 0 name", "")),
                "group_1_name": _clean(exp_first.get("Group 1 name", "")),
                "group_1_definition": _clean(exp_first.get("Group 1 definition", "")),
                "group_0_sample_size": _clean_int(exp_first.get("Group 0 sample size", "")),
                "group_1_sample_size": _clean_int(exp_first.get("Group 1 sample size", "")),
                "antibiotics_exclusion": _clean(exp_first.get("Antibiotics exclusion", "")),
                "sequencing_type": _clean(exp_first.get("Sequencing type", "")),
                "variable_region_16s": _clean(exp_first.get("16S variable region", "")),
                "sequencing_platform": _clean(exp_first.get("Sequencing platform", "")),
                "data_transformation": _clean(exp_first.get("Data transformation", "")),
                "statistical_test": _clean(exp_first.get("Statistical test", "")),
                "significance_threshold": _clean_float(exp_first.get("Significance threshold", "")),
                "mht_correction": _clean(exp_first.get("MHT correction", "")),
                "lda_score_above": _clean_float(exp_first.get("LDA Score above", "")),
                "matched_on": _clean(exp_first.get("Matched on", "")),
                "confounders_controlled_for": _clean(
                    exp_first.get("Confounders controlled for", "")
                ),
                "pielou": _clean(exp_first.get("Pielou", "")),
                "shannon": _clean(exp_first.get("Shannon", "")),
                "chao1": _clean(exp_first.get("Chao1", "")),
                "simpson": _clean(exp_first.get("Simpson", "")),
                "inverse_simpson": _clean(exp_first.get("Inverse Simpson", "")),
                "richness": _clean(exp_first.get("Richness", "")),
                "signatures": sig_records,
            }
        )

    year_str = first.get("Year", "").strip()
    pub_year: int | None = None
    if year_str.isdigit():
        pub_year = int(year_str)
    elif year_str and year_str not in ("", "NA"):
        log.append(
            {
                "level": "info",
                "entity": "study",
                "id": study_id,
                "reason": f"non-numeric publication year dropped: {year_str!r}",
            }
        )

    pmid_raw = first.get("PMID", "").strip()
    pmid: int | None = None
    if pmid_raw.isdigit():
        pmid = int(pmid_raw)
    elif pmid_raw and pmid_raw not in ("", "NA"):
        log.append(
            {
                "level": "info",
                "entity": "study",
                "id": study_id,
                "reason": f"non-numeric PMID dropped: {pmid_raw!r}",
            }
        )

    doi_raw = first.get("DOI", "").strip()
    doi_clean = doi_raw.removeprefix("https://doi.org/").removeprefix("http://doi.org/")
    doi: str | None = None
    if doi_clean.startswith("10."):
        doi = doi_clean
        if doi_clean != doi_raw and doi_raw not in ("", "NA"):
            log.append(
                {
                    "level": "info",
                    "entity": "study",
                    "id": study_id,
                    "reason": f"DOI normalized from URL: {doi_raw!r} → {doi_clean!r}",
                }
            )
    elif doi_raw and doi_raw not in ("", "NA"):
        log.append(
            {
                "level": "info",
                "entity": "study",
                "id": study_id,
                "reason": f"unrecognized DOI format dropped: {doi_raw!r}",
            }
        )

    url_raw = first.get("URL", "").strip()
    url: str | None = None
    if url_raw.startswith("http"):
        url = url_raw
    elif url_raw and url_raw not in ("", "NA"):
        log.append(
            {
                "level": "info",
                "entity": "study",
                "id": study_id,
                "reason": f"non-HTTP URL dropped: {url_raw!r}",
            }
        )

    return {
        "id": study_id,
        "source_record_id": first.get("Study", "").strip(),
        "pmid": pmid,
        "title": _clean(first.get("Title", "")),
        "publication_year": pub_year,
        "doi": doi,
        "url": url,
        "study_design": _clean(first.get("Study design", "")),
        "authors_list": _clean(first.get("Authors list", "")),
        "journal": _clean(first.get("Journal", "")),
        "keywords": _clean(first.get("Keywords", "")),
        "experiments": experiment_records,
    }


def _parse_taxa(
    tax_ids: str,
    names: str,
    sig_id: str,
    log: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Parse BugSigDB taxon ID and name strings into Taxon dicts.

    BugSigDB encodes taxa as:
    - Semicolons separate distinct taxa within one signature
    - Pipes separate lineage levels within one taxon (kingdom→...→leaf)
    - The last pipe-separated element is the most specific (leaf) taxon
    - MetaPhlAn names use the same delimiters with rank prefixes (k__, p__, etc.)
    """
    if not tax_ids or tax_ids.strip() in ("", "NA"):
        return []

    # Semicolons separate distinct taxa
    id_chains = [c.strip() for c in tax_ids.split(";") if c.strip()]
    # BugSigDB uses different delimiters: semicolons for IDs, commas for names
    name_chains = [c.strip() for c in names.split(",") if c.strip()] if names else []

    taxa: list[dict[str, str]] = []
    for i, id_chain in enumerate(id_chains):
        # Pipes separate lineage levels; take the last (most specific)
        lineage_ids = [x.strip() for x in id_chain.split("|") if x.strip()]
        if not lineage_ids:
            continue

        leaf_id = lineage_ids[-1]
        if leaf_id.isdigit():
            curie = f"NCBITaxon:{leaf_id}"
        elif leaf_id.startswith("NCBITaxon:"):
            curie = leaf_id
        else:
            log.append(
                {
                    "level": "skip",
                    "entity": "taxon",
                    "id": leaf_id,
                    "signature": sig_id,
                    "reason": f"unparsable leaf taxon ID: {leaf_id!r}",
                }
            )
            continue

        if len(lineage_ids) > 1:
            log.append(
                {
                    "level": "info",
                    "entity": "taxon",
                    "id": curie,
                    "signature": sig_id,
                    "reason": f"lineage with {len(lineage_ids)} levels, used leaf",
                }
            )

        # Extract name from MetaPhlAn lineage if available
        name: str
        rank: str | None = None
        if i < len(name_chains):
            name_parts = [x.strip() for x in name_chains[i].split("|") if x.strip()]
            if name_parts:
                leaf_name = name_parts[-1]
                # Parse MetaPhlAn rank prefix (k__, p__, g__, s__, etc.)
                name, rank = _parse_metaphlan_name(leaf_name)
            else:
                name = f"taxon_{leaf_id}"
                log.append(
                    {
                        "level": "info",
                        "entity": "taxon",
                        "id": curie,
                        "signature": sig_id,
                        "reason": "empty name chain, using placeholder",
                    }
                )
        else:
            name = f"taxon_{leaf_id}"
            log.append(
                {
                    "level": "info",
                    "entity": "taxon",
                    "id": curie,
                    "signature": sig_id,
                    "reason": f"no name chain at index {i}, using placeholder",
                }
            )

        if rank is None:
            rank = _taxon_name_to_rank(name)

        taxa.append(
            {
                "id": curie,
                "taxon_name": name,
                "taxonomic_rank": rank,
            }
        )

    return taxa


METAPHLAN_RANK_MAP: dict[str, str] = {
    "k__": "superkingdom",
    "p__": "phylum",
    "c__": "class",
    "o__": "order",
    "f__": "family",
    "g__": "genus",
    "s__": "species",
    "t__": "strain",
}


def _parse_metaphlan_name(raw: str) -> tuple[str, str | None]:
    """Parse a MetaPhlAn name like ``g__Enterococcus`` into (name, rank)."""
    for prefix, rank in METAPHLAN_RANK_MAP.items():
        if raw.startswith(prefix):
            return raw[len(prefix) :], rank
    return raw, None
