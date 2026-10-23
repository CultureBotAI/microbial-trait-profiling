"""Results the deck and the Tufte plates quote from tables other than the Table 1 sign tests.

The cover, Methods Plate II and Plate V typed Arm 2's test and call counts, the oxygen-guild
McNemar test and the post-selection rejection rate of the screened sign test (#579). A rerun
whose analysis code changed any of them left the text stale. Reading them here keeps the
text in step with the tables:

- ``trait_enrichment_results.tsv`` and ``trait_depletion_results.tsv`` (``just analyze`` and
  ``just analyze-depletion``) head their rows with ``# Key: value`` lines that give the tests
  and calls;
- ``guild_aerobe_anaerobe_correlation.tsv`` (``just guild-maps``) gives the McNemar test;
- ``calibration_criteria.tsv`` (``just calibration-report``) gives the rejection rate.
"""

from __future__ import annotations

import csv
import json
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from mtp.taxa_trait_matrix import STRAIN_ROLLUPS

ORA_ENRICHMENT = Path("data/derived/trait_enrichment_results.tsv")
ORA_DEPLETION = Path("data/derived/trait_depletion_results.tsv")
GUILD_CORRELATION = Path("data/derived/guild_maps/guild_aerobe_anaerobe_correlation.tsv")
CALIBRATION_CRITERIA = Path("data/derived/calibration/calibration_criteria.tsv")

# Shared column and basis definitions for taxa-coverage tables (#858, #864, #956, #982).
COVERAGE_TOTAL_TAXA_COL = "total_taxa"
COVERAGE_TRAIT_BEARING_TAXA_COL = "trait_bearing_taxa"
COVERAGE_ZERO_TRAIT_TAXA_COL = "zero_trait_taxa"
COVERAGE_TRAIT_BEARING_PCT_COL = "trait_bearing_pct"
COVERAGE_ZERO_TRAIT_PCT_COL = "zero_trait_pct"
COVERAGE_TOTAL_OCC_COL = "total_occurrences"
COVERAGE_TRAIT_BEARING_OCC_COL = "trait_bearing_occurrences"
COVERAGE_ZERO_TRAIT_OCC_COL = "zero_trait_occurrences"
COVERAGE_TRAIT_BEARING_OCC_PCT_COL = "trait_bearing_occurrences_pct"
COVERAGE_ZERO_TRAIT_OCC_PCT_COL = "zero_trait_occurrences_pct"
COVERAGE_ANNOTATION_BASIS = "any annotation [incl. MediaDive]"


@dataclass(frozen=True)
class CoverageRunInfo:
    """Provenance and profile/study counts for a taxa-coverage run (#858, #956, #982)."""

    strain_rollup: str
    population: Mapping[str, Any] | None
    selected_profiles: int | None
    retained_profiles: int
    retained_studies: int


def format_coverage_population(population: Mapping[str, Any] | str | None) -> str:
    """Format the rank-table population label from provenance.json (#982)."""
    if isinstance(population, str) and population.strip():
        return population.strip()
    if not isinstance(population, Mapping):
        return "population unrecorded in rank table"
    host_raw = str(population.get("host_species") or "").strip()
    host_label = "human" if host_raw.lower() == "homo sapiens" else host_raw
    seq_types = population.get("sequencing_types")
    if isinstance(seq_types, Sequence) and not isinstance(seq_types, (str, bytes)):
        seq_label = "/".join(str(s) for s in seq_types if str(s).strip())
    else:
        seq_label = ""
    parts = [p for p in ("all", host_label, seq_label, "signatures in rank table") if p]
    return " ".join(parts) if (host_label or seq_label) else "population unrecorded in rank table"


def read_coverage_run_info(
    coverage_dir: Path,
    *,
    profile_summary_rows: Sequence[Mapping[str, str]] | None = None,
    study_meta_rows: Sequence[Mapping[str, str]] | None = None,
    require_counts: bool = False,
) -> CoverageRunInfo:
    """Read and cross-check taxa-coverage provenance and profile counts (#858, #956, #982)."""
    prov_path = coverage_dir / "provenance.json"
    strain_rollup = "unrecorded"
    population: Mapping[str, Any] | None = None
    selected_profiles: int | None = None
    if prov_path.exists():
        data = json.loads(prov_path.read_text(encoding="utf-8"))
        raw_scope = data.get("strain_rollup")
        if raw_scope is None and isinstance(data.get("parameters"), Mapping):
            raw_scope = data["parameters"].get("strain_rollup")
        if raw_scope is not None and raw_scope != "":
            if not isinstance(raw_scope, str) or raw_scope not in STRAIN_ROLLUPS:
                raise ValueError(
                    f"Unknown strain_rollup in {prov_path}: {raw_scope!r} "
                    f"(expected one of {STRAIN_ROLLUPS})"
                )
            strain_rollup = raw_scope
        raw_pop = data.get("population")
        if isinstance(raw_pop, Mapping):
            population = raw_pop
        raw_counts = data.get("counts")
        if isinstance(raw_counts, Mapping) and "selected_profiles" in raw_counts:
            selected_profiles = int(raw_counts["selected_profiles"])

    if profile_summary_rows is None:
        summary_path = coverage_dir / "profile_trait_mapping_summary.tsv"
        if summary_path.exists():
            with summary_path.open(newline="", encoding="utf-8") as handle:
                lines = [line for line in handle if not line.startswith("#")]
                profile_summary_rows = list(csv.DictReader(lines, delimiter="\t"))

    summary_profiles = len(profile_summary_rows) if profile_summary_rows else None
    summary_studies = (
        len({r["study_id"] for r in profile_summary_rows if "study_id" in r})
        if profile_summary_rows
        else None
    )
    meta_profiles = len(study_meta_rows) if study_meta_rows else None
    meta_studies = (
        len({r["study_id"] for r in study_meta_rows if "study_id" in r})
        if study_meta_rows
        else None
    )

    if (
        selected_profiles is not None
        and summary_profiles is not None
        and selected_profiles != summary_profiles
    ):
        raise ValueError(
            f"taxa-coverage provenance.json selected_profiles ({selected_profiles}) "
            f"disagrees with profile_trait_mapping_summary.tsv ({summary_profiles} rows)"
        )
    if summary_profiles is not None and meta_profiles is not None:
        studies_disagree = (
            bool(summary_studies) and bool(meta_studies) and summary_studies != meta_studies
        )
        if summary_profiles != meta_profiles or studies_disagree:
            raise ValueError(
                f"taxa-coverage profile_trait_mapping_summary.tsv ({summary_profiles} profiles "
                f"across {summary_studies} studies) disagrees with "
                f"study_trait_profile_study_metadata.tsv ({meta_profiles} profiles across "
                f"{meta_studies} studies)"
            )
    if (
        selected_profiles is not None
        and meta_profiles is not None
        and selected_profiles != meta_profiles
    ):
        raise ValueError(
            f"taxa-coverage provenance.json selected_profiles ({selected_profiles}) "
            f"disagrees with study_trait_profile_study_metadata.tsv ({meta_profiles} rows)"
        )

    if summary_profiles is not None and summary_studies is not None:
        retained_profiles = summary_profiles
        retained_studies = summary_studies or (meta_studies or 0)
    elif selected_profiles is not None:
        retained_profiles = selected_profiles
        retained_studies = meta_studies or 0
    elif meta_profiles is not None and meta_studies is not None:
        retained_profiles = meta_profiles
        retained_studies = meta_studies
    elif require_counts:
        raise ValueError(
            "Cannot determine coverage profile count: provide profile_trait_mapping_summary.tsv "
            "or provenance.json"
        )
    else:
        retained_profiles = 0
        retained_studies = 0

    return CoverageRunInfo(
        strain_rollup=strain_rollup,
        population=population,
        selected_profiles=selected_profiles,
        retained_profiles=retained_profiles,
        retained_studies=retained_studies,
    )


@dataclass(frozen=True)
class OraCalls:
    """One ORA table's tests and calls, with the thresholds that define a call."""

    tests: int
    calls: int
    p_threshold: str
    fold: str

    @classmethod
    def from_header(cls, lines: Iterable[str]) -> OraCalls:
        """Read the ``# Key: value`` lines heading the table; the rows are not read. The
        calls are quoted as BH q, so a table whose calls are not FDR-controlled is refused."""
        header: dict[str, str] = {}
        for line in lines:
            if not line.startswith("#"):
                break
            key, _, value = line[1:].partition(":")
            header[key.strip()] = value.strip()
        if header.get("FDR") != "True":
            raise ValueError(
                f"the ORA table's calls are not FDR-controlled (FDR: {header.get('FDR')})"
            )
        try:
            return cls(
                tests=int(header["Total-Tests-Performed"]),
                calls=int(header["Significant-Hits"]),
                p_threshold=header["P-Threshold"],
                fold=header.get("Min-Fold") or header["Max-Fold"],
            )
        except (KeyError, ValueError) as error:
            raise ValueError(f"an ORA table header without its counts ({error})") from error


@dataclass(frozen=True)
class McNemar:
    """The paired oxygen-guild test: studies whose aerobe and anaerobe shares moved in
    opposite directions, each way, and the exact McNemar test on them."""

    aerobe_up_anaerobe_down: int
    anaerobe_up_aerobe_down: int
    odds_ratio: float
    p_value: float

    @property
    def discordant(self) -> int:
        return self.aerobe_up_anaerobe_down + self.anaerobe_up_aerobe_down

    @classmethod
    def from_rows(cls, rows: Iterable[Mapping[str, str]]) -> McNemar:
        discordant: dict[str, int] = {}
        test: Mapping[str, str] | None = None
        for row in rows:
            if row["statistic"] == "discordant_studies":
                discordant[row["comparison"]] = int(row["value"])
            elif row["statistic"] == "mcnemar_odds_ratio":
                test = row
        if test is None or not test["value"] or not test["p_value"]:
            raise ValueError("the guild table has no McNemar odds ratio and p-value")
        try:
            return cls(
                aerobe_up_anaerobe_down=discordant["aerobe_up_anaerobe_down"],
                anaerobe_up_aerobe_down=discordant["anaerobe_up_aerobe_down"],
                odds_ratio=float(test["value"]),
                p_value=float(test["p_value"]),
            )
        except KeyError as error:
            raise ValueError(f"the guild table has no {error} discordant count") from error


def screened_rejection_rate(rows: Iterable[Mapping[str, str]]) -> float:
    """The pooled per-test rejection rate of the screened fraction sign test under the
    complete null: what screening the traits on the same data does to its size. A run with
    fewer replicates than the protocol requires, such as the documented five-replicate quick
    pass, certifies nothing, so it is refused (#704)."""
    for row in rows:
        if (row["criterion"], row["procedure"]) == (
            "per_test_level",
            "screened_fraction_sign_test",
        ):
            recorded, required = int(row["replicates_recorded"]), int(row["replicates_required"])
            if recorded < required:
                raise ValueError(
                    f"the calibration ran {recorded} of the {required} replicates its protocol "
                    "requires; rerun just calibration-report at the declared budget"
                )
            return float(row["estimate"])
    raise ValueError("the calibration criteria have no per-test level for the screened sign test")
