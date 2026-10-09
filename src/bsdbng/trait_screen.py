"""Is there a trait outside Table 1 that shifts consistently? (#255, #265)

The directional-consistency arm has until now tested exactly the six traits of
manuscript Table 1. Those six were chosen as the intersection of the composition
and enrichment arms' significant sets on this same corpus, so the third arm has
been a third statistic computed on a feature set the other two arms had already
agreed on. It could not discover a seventh trait, and its agreement with the
other arms was not confirmation, because the features were picked for having
already agreed.

This module runs the same test over the whole retained trait space instead. The
statistic does not change: one vote per study per trait, the sign of the paired
increased-minus-decreased contrast, ties excluded, a two-sided exact binomial
test against a half-and-half null with an exact Clopper-Pearson interval. The
study contrasts and the sign tally come from ``bsdbng.table1_trait_maps``, so
this is literally the same rule rather than a second copy of it.

What does change is the family. It is now every retained non-MediaDive column of
the pinned matrix, declared in ``methods/trait_screen_manifest.json`` before the
first run, and corrected with Benjamini-Yekutieli rather than Benjamini-Hochberg
because the candidates are strongly and positively dependent. Because they are,
the nominal family size overstates how many independent chances to be wrong the
data contain, so the significant set is reported twice: once at the nominal
family size and once at the effective one (#265). Neither is the answer on its
own; their agreement or disagreement is.

The six Table 1 traits appear here as a flagged subset of the full table. Their
separately published six-trait q-values are produced by ``table1_trait_maps`` and
are untouched by this module; revising any manuscript number is #205.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections.abc import Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from statistics import fmean, median

from scipy.linalg import eigvalsh
from scipy.stats import binomtest, spearmanr

from bsdbng.paths import DERIVED_DATA_DIR, REPO_ROOT
from bsdbng.table1_trait_maps import study_trait_contrasts, tally_direction_signs
from bsdbng.trait_guilds import (
    DEFAULT_EXPERIMENTS_PATH,
    DEFAULT_FULL_DUMP_PATH,
    DEFAULT_PROFILE_PATH,
    TABLE1_TRAIT_LABELS,
    TABLE1_TRAITS,
    CorpusStats,
    Profile,
    column_term_id,
    load_profiles,
    read_bugsigdb_vintage,
    retained_trait_columns,
    write_stats_report,
)

VERSION = "trait-direction-screen-v1"
MANIFEST_PATH = REPO_ROOT / "methods" / "trait_screen_manifest.json"
DEFAULT_OUTPUT_DIR = DERIVED_DATA_DIR / "trait_screen"
DEFAULT_COLUMN_LABELS_PATH = DERIVED_DATA_DIR / "taxa_trait_matrix_columns.tsv"

CONFIDENCE_LEVEL = 0.95
Q_THRESHOLDS: tuple[float, ...] = (0.05, 0.001)

# The estimator whose effective family size drives the sensitivity correction.
# Galwey is reported beside it; Cheverud/Nyholt is reported and excluded.
EFFECTIVE_ESTIMATOR = "li_ji"

TESTED = "tested"
ALL_TIED = "all_tied"
NO_ELIGIBLE_STUDIES = "no_eligible_studies"

MIN_DEPENDENCE_UNITS = 3
MIN_DEPENDENCE_TRAITS = 3

CONTRAST_BASIS = "study_contrasts"
PROFILE_BASIS = "profiles"


@dataclass(frozen=True)
class TraitScreenRow:
    """One candidate trait's directional vote over studies.

    ``share_up`` and the interval, the p-value and both q-values are ``None``
    whenever the test was not defined. A trait every study ties on has no share
    of 0.0 and no p-value of 1.0; it has no result, and the table says so.
    """

    trait: str
    term_id: str
    trait_label: str
    object_name: str
    object_category: str
    in_table1: bool
    table1_slug: str
    studies_up: int
    studies_down: int
    studies_tied: int
    share_up: float | None
    ci_low: float | None
    ci_high: float | None
    p_value: float | None
    bh_q_nominal: float | None
    by_q_nominal: float | None
    by_q_effective: float | None
    test_status: str

    @property
    def studies_scored(self) -> int:
        return self.studies_up + self.studies_down


@dataclass(frozen=True)
class Dependence:
    """How much independent information a set of candidate traits carries.

    ``basis`` names what the correlation was computed over, because the answer
    differs: the study contrasts are what this arm actually tests, while the
    profiles are the basis quoted in #263 and #265 and are retained so those
    figures can be reproduced.
    """

    basis: str
    units: int
    traits: int
    constant_traits_excluded: int
    pc1_share: float
    eigenvalues: tuple[float, ...]
    li_ji: float
    galwey: float
    cheverud_nyholt: float
    pairs: int
    median_rho: float
    mean_rho: float
    share_rho_above_half: float
    share_rho_above_nine_tenths: float

    @property
    def effective(self) -> float:
        return self.li_ji if EFFECTIVE_ESTIMATOR == "li_ji" else self.galwey


@dataclass(frozen=True)
class ScreenResult:
    """The whole screen: per-trait rows, the family it was corrected in, and why."""

    rows: list[TraitScreenRow]
    studies: int
    family_size: int
    effective_family_size: float | None
    effective_family_size_used: int | None
    contrast_dependence: Dependence | None
    profile_dependence: Dependence | None


def read_header(profile_path: Path) -> list[str]:
    """The profile matrix's header row, without reading the matrix itself."""
    with profile_path.open(encoding="utf-8", newline="") as handle:
        line = handle.readline()
    if not line:
        msg = f"{profile_path} is empty"
        raise ValueError(msg)
    return line.rstrip("\n").split("\t")


def screen_candidates(profile_path: Path = DEFAULT_PROFILE_PATH) -> tuple[str, ...]:
    """Every retained trait column of the matrix, in header order.

    The candidate family is read from the corpus rather than written down, so a
    corpus change moves the family instead of silently mismatching it.
    """
    _, columns, _ = retained_trait_columns(read_header(profile_path))
    return tuple(columns)


def candidate_slots(columns: Sequence[str]) -> dict[str, tuple[str, ...]]:
    """One single-column trait slot per candidate, for ``load_profiles``.

    Keyed on the full column name. ``_select_trait_columns`` raises if a key
    resolves to more than one column; on the pinned matrix no retained column's
    term suffix collides with another's, so this passes. If a later corpus breaks
    that, the raise is the correct outcome and must not be suppressed: two
    columns summed into one slot would be double-counted against a denominator
    that counts them once.
    """
    return {column: (column,) for column in columns}


def harmonic(count: int) -> float:
    """The Benjamini-Yekutieli dependence factor, sum of 1/i up to ``count``."""
    return math.fsum(1.0 / index for index in range(1, count + 1))


def step_up(p_values: Sequence[float], *, method: str, family_size: int) -> list[float]:
    """Benjamini-Hochberg or Benjamini-Yekutieli at an explicit family size.

    ``scipy.stats.false_discovery_control`` always corrects at the number of
    p-values handed to it. #265 asks for the same procedure run as though the
    family held its effective number of tests, so the multiplier has to be a
    parameter. With ``family_size == len(p_values)`` this returns exactly what
    scipy returns, which ``tests/test_trait_screen.py`` asserts against scipy on
    random input for both methods.

    An adjusted p-value is never reported below its own raw p-value. That floor
    only binds when ``family_size < len(p_values)``, that is, on the effective-m
    sensitivity run, where the substitution would otherwise deflate a large-rank
    p-value below itself.
    """
    if method not in ("bh", "by"):
        msg = f"unknown correction method {method!r}"
        raise ValueError(msg)
    if family_size < 1:
        msg = "family_size must be at least 1"
        raise ValueError(msg)
    if not p_values:
        return []
    factor = family_size * (harmonic(family_size) if method == "by" else 1.0)
    order = sorted(range(len(p_values)), key=lambda index: p_values[index])
    adjusted = [0.0] * len(p_values)
    running = math.inf
    for rank, index in reversed(list(enumerate(order, start=1))):
        running = min(running, factor * p_values[index] / rank)
        adjusted[index] = min(1.0, max(running, p_values[index]))
    return adjusted


def effective_family_size(eigenvalues: Sequence[float]) -> tuple[float, float, float]:
    """Li & Ji, Galwey and Cheverud/Nyholt effective numbers of tests.

    Cheverud/Nyholt is computed so that it can be named and excluded rather than
    silently omitted (#265). Li & Ji (2005) report that its eigenvalue-variance
    estimate is overly large, giving excessively conservative results; Salyakina
    et al. (2005) found its nominal error rate off in both directions depending
    on the dependence structure. It is excluded as unreliable, not as lenient
    (#325, #357).

    All three return the trait count for an identity correlation matrix and 1.0
    for a perfectly correlated one.
    """
    count = len(eigenvalues)
    if count < 2:
        msg = "an effective family size needs at least two candidates"
        raise ValueError(msg)
    li_ji = math.fsum(
        (1.0 if abs(value) >= 1.0 else 0.0) + (abs(value) - math.floor(abs(value)))
        for value in eigenvalues
    )
    positive = [value for value in eigenvalues if value > 0.0]
    galwey = math.fsum(math.sqrt(value) for value in positive) ** 2 / math.fsum(positive)
    mean = math.fsum(eigenvalues) / count
    variance = math.fsum((value - mean) ** 2 for value in eigenvalues) / (count - 1)
    cheverud = 1.0 + (count - 1) * (1.0 - variance / count)
    return li_ji, galwey, cheverud


def dependence(matrix: Sequence[Sequence[float]], basis: str) -> Dependence | None:
    """Spearman dependence among candidate traits, and what it implies for the family.

    A candidate that is constant across the units is dropped from the correlation
    and counted, because Spearman is undefined against a constant series and
    would otherwise put a nan through the whole spectrum. It stays in the tested
    family; only the dependence estimate loses it.

    Returns ``None`` when there is too little to estimate from — too few units or
    too few varying candidates. The dependence is a diagnostic about the family,
    not the quantity under test, so a corpus too small to support it leaves the
    effective family size undefined rather than failing the screen. Every output
    then says so instead of carrying a fabricated number.
    """
    if not matrix or len(matrix) < MIN_DEPENDENCE_UNITS:
        return None
    width = len(matrix[0])
    keep = [index for index in range(width) if len({row[index] for row in matrix}) > 1]
    excluded = width - len(keep)
    if len(keep) < MIN_DEPENDENCE_TRAITS:
        return None
    reduced = [[row[index] for index in keep] for row in matrix]
    correlation, _ = spearmanr(reduced)
    flat = [float(correlation[i][j]) for i in range(len(keep)) for j in range(i + 1, len(keep))]
    eigenvalues = sorted((float(value) for value in eigvalsh(correlation)), reverse=True)
    li_ji, galwey, cheverud = effective_family_size(eigenvalues)
    total = math.fsum(eigenvalues)
    return Dependence(
        basis=basis,
        units=len(matrix),
        traits=len(keep),
        constant_traits_excluded=excluded,
        pc1_share=eigenvalues[0] / total,
        eigenvalues=tuple(eigenvalues),
        li_ji=li_ji,
        galwey=galwey,
        cheverud_nyholt=cheverud,
        pairs=len(flat),
        median_rho=float(median(flat)),
        mean_rho=float(fmean(flat)),
        share_rho_above_half=sum(value > 0.5 for value in flat) / len(flat),
        share_rho_above_nine_tenths=sum(value > 0.9 for value in flat) / len(flat),
    )


def table1_slug(column: str) -> str | None:
    """The Table 1 trait this matrix column carries, or ``None`` if it is not one.

    ``TABLE1_TRAITS`` is keyed by slug and valued by term, while the screen's
    candidates are raw matrix columns, so the six have to be matched through the
    same full-column-then-term-suffix rule the loader uses. Two of the six sit
    behind a METPO qualifier (``METPO:2000103|GO:0030245``), which is exactly the
    case a naive membership test gets wrong.
    """
    by_term = {term: slug for slug, terms in TABLE1_TRAITS.items() for term in terms}
    return by_term.get(column) or by_term.get(column_term_id(column))


@dataclass(frozen=True)
class ColumnLabel:
    """What a matrix column stands for, read from the trait matrix's own column index."""

    label: str
    object_name: str
    object_category: str


def load_column_labels(path: Path = DEFAULT_COLUMN_LABELS_PATH) -> dict[str, ColumnLabel]:
    """Human-readable names for matrix columns, from ``taxa_trait_matrix_columns.tsv``.

    A screen over hundreds of bare CURIEs is unreadable, and the pipeline already
    emits the labels beside the matrix it built, so they are read rather than
    re-derived. The file is a regenerable pipeline output; when it is absent the
    screen still runs and every candidate falls back to its column name, which
    the provenance record states.
    """
    if not path.exists():
        return {}
    labels: dict[str, ColumnLabel] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            labels[row["column_id"]] = ColumnLabel(
                label=row.get("column_label", "") or row["column_id"],
                object_name=row.get("object_name", ""),
                object_category=row.get("object_category", ""),
            )
    return labels


def run_screen(
    profiles: list[Profile],
    columns: Sequence[str],
    *,
    dependence_profiles: list[Profile] | None = None,
    labels: dict[str, ColumnLabel] | None = None,
) -> ScreenResult:
    """Vote every study once per candidate trait, then correct over the whole family.

    ``dependence_profiles`` is the unpaired human corpus used only for the
    secondary, profile-basis dependence estimate quoted in #263 and #265. The
    tested quantity and the correction are unaffected by it.
    """
    contrasts = study_trait_contrasts(profiles, columns)
    tallies = tally_direction_signs(contrasts, columns)

    tests = {}
    for trait in columns:
        up, down, _ = tallies[trait]
        if up + down:
            tests[trait] = binomtest(up, up + down, p=0.5)

    contrast_dependence = dependence(
        [[row.get(trait, 0.0) for trait in columns] for row in contrasts.values()],
        CONTRAST_BASIS,
    )
    profile_dependence = (
        dependence(
            [[profile.guild_counts[trait] for trait in columns] for profile in dependence_profiles],
            PROFILE_BASIS,
        )
        if dependence_profiles
        else None
    )

    family_size = len(columns)
    effective = contrast_dependence.effective if contrast_dependence is not None else None
    effective_used = max(1, min(family_size, round(effective))) if effective is not None else None

    # Untestable hypotheses are held at p = 1 inside every correction and
    # reported blank, the convention PR #181 set for the ratio family.
    held = [float(tests[trait].pvalue) if trait in tests else 1.0 for trait in columns]
    bh_nominal = step_up(held, method="bh", family_size=family_size)
    by_nominal = step_up(held, method="by", family_size=family_size)
    by_effective: list[float | None] = (
        list(step_up(held, method="by", family_size=effective_used))
        if effective_used is not None
        else [None] * len(columns)
    )

    named = labels or {}
    rows: list[TraitScreenRow] = []
    for index, trait in enumerate(columns):
        up, down, tied = tallies[trait]
        test = tests.get(trait)
        if test is None:
            status = ALL_TIED if tied else NO_ELIGIBLE_STUDIES
            interval_low = interval_high = share = None
            p_value = bh_value = by_value = by_effective_value = None
        else:
            status = TESTED
            interval = test.proportion_ci(confidence_level=CONFIDENCE_LEVEL)
            share = up / (up + down)
            interval_low, interval_high = float(interval.low), float(interval.high)
            p_value = float(test.pvalue)
            bh_value = bh_nominal[index]
            by_value = by_nominal[index]
            by_effective_value = by_effective[index]
        slug = table1_slug(trait)
        naming = named.get(trait, ColumnLabel(trait, "", ""))
        rows.append(
            TraitScreenRow(
                trait=trait,
                term_id=column_term_id(trait),
                trait_label=naming.label,
                object_name=naming.object_name,
                object_category=naming.object_category,
                in_table1=slug is not None,
                table1_slug=slug or "",
                studies_up=up,
                studies_down=down,
                studies_tied=tied,
                share_up=share,
                ci_low=interval_low,
                ci_high=interval_high,
                p_value=p_value,
                bh_q_nominal=bh_value,
                by_q_nominal=by_value,
                by_q_effective=by_effective_value,
                test_status=status,
            )
        )

    return ScreenResult(
        rows=rows,
        studies=len(contrasts),
        family_size=family_size,
        effective_family_size=effective,
        effective_family_size_used=effective_used,
        contrast_dependence=contrast_dependence,
        profile_dependence=profile_dependence,
    )


def significant(rows: Sequence[TraitScreenRow], field: str, threshold: float) -> set[str]:
    """The traits a q column puts at or below a threshold. Undefined is not significant."""
    chosen = set()
    for row in rows:
        value = getattr(row, field)
        if value is not None and value <= threshold:
            chosen.add(row.trait)
    return chosen


def _cell(value: float | int | str | bool | None, digits: str = ".6g") -> str:
    """One atomic value per cell; an undefined statistic is blank, never zero."""
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return format(value, digits)
    return str(value)


def _sha256(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def _family_rows(result: ScreenResult) -> list[dict[str, str]]:
    """The family record: its size, its dependence, and what that does to the set."""
    rows: list[dict[str, str]] = []

    def add(quantity: str, basis: str, value: float | int | str | None, note: str) -> None:
        rows.append({"quantity": quantity, "basis": basis, "value": _cell(value), "note": note})

    add("family_size_nominal", "", result.family_size, "every retained non-MediaDive trait")
    add(
        "family_size_effective",
        CONTRAST_BASIS,
        result.effective_family_size,
        f"{EFFECTIVE_ESTIMATOR}; rounded to {result.effective_family_size_used} for correction"
        if result.effective_family_size is not None
        else "undefined: the corpus is too small to estimate dependence from",
    )
    add(
        "dependence_status",
        CONTRAST_BASIS,
        "estimated" if result.contrast_dependence is not None else "unavailable",
        "an unavailable estimate leaves by_q_effective blank rather than guessing",
    )
    add("studies", "", result.studies, "studies contributing at least one contrast")

    bases = [
        record
        for record in (result.contrast_dependence, result.profile_dependence)
        if record is not None
    ]
    for record in bases:
        primary = "tested quantity" if record.basis == CONTRAST_BASIS else "reproduces #265"
        add("units", record.basis, record.units, primary)
        add("traits_in_dependence", record.basis, record.traits, "candidates with variation")
        add(
            "constant_traits_excluded",
            record.basis,
            record.constant_traits_excluded,
            "kept in the tested family, dropped from the correlation only",
        )
        add("effective_li_ji", record.basis, record.li_ji, "Li & Ji (2005)")
        add("effective_galwey", record.basis, record.galwey, "Galwey (2009)")
        add(
            "effective_cheverud_nyholt",
            record.basis,
            record.cheverud_nyholt,
            "reported and EXCLUDED: unreliable estimator (Li & Ji 2005; Salyakina et al. 2005)",
        )
        add("pc1_variance_share", record.basis, record.pc1_share, "first principal component")
        add("trait_pairs", record.basis, record.pairs, "unordered candidate pairs")
        add("median_spearman_rho", record.basis, record.median_rho, "between candidates")
        add("mean_spearman_rho", record.basis, record.mean_rho, "between candidates")
        add("share_rho_above_0.5", record.basis, record.share_rho_above_half, "")
        add("share_rho_above_0.9", record.basis, record.share_rho_above_nine_tenths, "")

    for threshold in Q_THRESHOLDS:
        nominal = significant(result.rows, "by_q_nominal", threshold)
        effective = significant(result.rows, "by_q_effective", threshold)
        shared = nominal & effective
        union = nominal | effective
        add(
            f"significant_by_q_{threshold}_nominal",
            "",
            len(nominal),
            f"BY at m = {result.family_size}",
        )
        add(
            f"significant_by_q_{threshold}_effective",
            "",
            len(effective),
            f"BY at m = {result.effective_family_size_used}",
        )
        add(
            f"membership_jaccard_q_{threshold}",
            "",
            len(shared) / len(union) if union else None,
            "1.0 means family size does not change which traits are called",
        )
        add(
            f"membership_symmetric_difference_q_{threshold}",
            "",
            len(union) - len(shared),
            "traits called under one family size and not the other",
        )
    return rows


def write_screen_outputs(
    result: ScreenResult,
    output_dir: Path,
    *,
    profile_path: Path,
    experiments_path: Path,
    full_dump_path: Path,
    corpus_stats: CorpusStats | None = None,
    code_commit: str | None = None,
) -> list[Path]:
    """Write the per-trait screen, the family record, the spectrum, hashes and a report."""
    output_dir.mkdir(parents=True, exist_ok=True)
    screen_path = output_dir / "trait_direction_screen.tsv"
    family_path = output_dir / "trait_screen_family.tsv"
    spectrum_path = output_dir / "trait_screen_correlation.tsv"
    provenance_path = output_dir / "trait_screen_provenance.json"
    report_path = output_dir / "trait_screen_report.md"
    vintage = read_bugsigdb_vintage(full_dump_path)

    fields = [
        "trait",
        "term_id",
        "trait_label",
        "object_name",
        "object_category",
        "in_table1",
        "table1_slug",
        "studies_up",
        "studies_down",
        "studies_tied",
        "studies_scored",
        "share_up",
        "ci_low",
        "ci_high",
        "binomial_p",
        "bh_q_nominal",
        "by_q_nominal",
        "by_q_effective",
        "test_status",
        "significant_by_q_0.05_nominal",
        "significant_by_q_0.05_effective",
        "significant_by_q_0.001_nominal",
        "significant_by_q_0.001_effective",
        "family_id",
        "family_size_nominal",
        "family_size_effective",
        "bugsigdb_vintage",
    ]
    with screen_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        for row in result.rows:
            writer.writerow(
                [
                    row.trait,
                    row.term_id,
                    row.trait_label,
                    row.object_name,
                    row.object_category,
                    _cell(row.in_table1),
                    row.table1_slug,
                    row.studies_up,
                    row.studies_down,
                    row.studies_tied,
                    row.studies_scored,
                    _cell(row.share_up),
                    _cell(row.ci_low),
                    _cell(row.ci_high),
                    _cell(row.p_value),
                    _cell(row.bh_q_nominal),
                    _cell(row.by_q_nominal),
                    _cell(row.by_q_effective),
                    row.test_status,
                    _cell(row.by_q_nominal is not None and row.by_q_nominal <= 0.05),
                    _cell(row.by_q_effective is not None and row.by_q_effective <= 0.05),
                    _cell(row.by_q_nominal is not None and row.by_q_nominal <= 0.001),
                    _cell(row.by_q_effective is not None and row.by_q_effective <= 0.001),
                    VERSION,
                    result.family_size,
                    result.effective_family_size_used,
                    vintage or "",
                ]
            )

    with family_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["quantity", "basis", "value", "note"])
        for entry in _family_rows(result):
            writer.writerow([entry["quantity"], entry["basis"], entry["value"], entry["note"]])

    with spectrum_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["basis", "component", "eigenvalue", "variance_share", "cumulative_share"])
        records = [
            record
            for record in (result.contrast_dependence, result.profile_dependence)
            if record is not None
        ]
        for record in records:
            total = math.fsum(record.eigenvalues)
            cumulative = 0.0
            for position, value in enumerate(record.eigenvalues, start=1):
                cumulative += value
                writer.writerow(
                    [
                        record.basis,
                        position,
                        _cell(value),
                        _cell(value / total),
                        _cell(cumulative / total),
                    ]
                )

    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    provenance = {
        "analysis_version": VERSION,
        "manifest": manifest["protocol"],
        "manifest_sha256": _sha256(MANIFEST_PATH),
        "comparison_spec": manifest["comparison_spec"],
        "profile_path": str(profile_path),
        "profile_sha256": _sha256(profile_path),
        "experiments_path": str(experiments_path),
        "experiments_sha256": _sha256(experiments_path),
        "analysis_source_sha256": _sha256(Path(__file__)),
        "shared_loader_source_sha256": _sha256(Path(__file__).with_name("trait_guilds.py")),
        "shared_contrast_source_sha256": _sha256(Path(__file__).with_name("table1_trait_maps.py")),
        "code_commit": code_commit,
        "bugsigdb_vintage": vintage,
        "input_filter_counts": asdict(corpus_stats) if corpus_stats else None,
        "family_size_nominal": result.family_size,
        "family_size_effective": result.effective_family_size,
        "family_size_effective_used": result.effective_family_size_used,
        "effective_estimator": EFFECTIVE_ESTIMATOR,
        "effective_estimator_excluded": "cheverud_nyholt",
        "correction_primary": "Benjamini-Yekutieli",
        "correction_also_reported": "Benjamini-Hochberg",
        "studies": result.studies,
        "candidates_with_a_resolved_label": sum(
            1 for row in result.rows if row.trait_label != row.trait
        ),
        "host": "Homo sapiens",
        "pairing": "study_id, experiment_id; both inc and dec required",
        "study_weighting": "one mean fraction per direction per study, then one vote per study",
        "null": "P(study contrast > 0 | non-tied) = 0.5",
        "confidence_level": CONFIDENCE_LEVEL,
        "thresholds_reported": list(Q_THRESHOLDS),
        "untestable_family_policy": "p=1 internally for the correction; p and q left blank in TSV",
        "table1_traits": TABLE1_TRAITS,
        "table1_note": (
            "flagged as a subset of this family; their separately published six-trait "
            "BH q-values are produced by table1_trait_maps and are not changed here (#205)"
        ),
    }
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")

    report_path.write_text(_report(result, vintage), encoding="utf-8")
    return [screen_path, family_path, spectrum_path, provenance_path, report_path]


def cheverud_note(record: Dependence) -> str:
    """Why Cheverud/Nyholt is excluded, claiming only what this record shows."""
    note = (
        "it is reported here only so that it can be named and excluded as an "
        "unreliable estimator, which Li & Ji (2005) found overly large and "
        "excessively conservative"
    )
    if record.cheverud_nyholt > max(record.li_ji, record.galwey):
        return (
            f"{note}. Here too it is the largest of the three, so it is not excluded "
            "for being lenient: leaving it out makes the effective-size run less strict."
        )
    return f"{note}."


def _dependence_paragraphs(result: ScreenResult) -> list[str]:
    """What the dependence actually measured, on each basis, without overstating it.

    The two bases disagree sharply on this corpus and the difference is the
    substantive point, so the text reports both rather than asserting one
    conclusion about "the" dependence.
    """
    contrast = result.contrast_dependence
    if contrast is None:
        return [
            "The dependence between candidates could not be estimated on this corpus, so "
            "no effective family size is reported and the effective-m column is blank "
            "throughout. The nominal correction is unaffected.",
        ]
    lines = [
        f"Between the quantities actually tested, the study-level contrasts, the "
        f"candidates are only weakly dependent: median Spearman "
        f"{contrast.median_rho:+.3f} over {contrast.pairs} pairs, with "
        f"{contrast.share_rho_above_half:.1%} above +0.5 and a first component holding "
        f"{contrast.pc1_share:.1%} of the variance. The effective number of tests is "
        f"{contrast.li_ji:.1f} by Li & Ji and {contrast.galwey:.1f} by Galwey against a "
        f"nominal {result.family_size}. Cheverud/Nyholt gives "
        f"{contrast.cheverud_nyholt:.1f}; {cheverud_note(contrast)}",
    ]
    profiles = result.profile_dependence
    if profiles is not None:
        lines.extend(
            [
                "",
                f"Between profiles, the basis quoted in #263 and #265, the same candidates "
                f"look far more dependent: median Spearman {profiles.median_rho:+.3f}, "
                f"{profiles.share_rho_above_half:.1%} above +0.5, first component "
                f"{profiles.pc1_share:.1%}, effective tests {profiles.li_ji:.1f} by Li & Ji "
                f"and {profiles.galwey:.1f} by Galwey. The gap between the two bases is "
                "itself a result: much of the between-trait correlation across profiles is "
                "shared annotation depth, and taking the paired within-study contrast "
                "removes a large part of it. The correction here uses the contrast basis "
                "because that is what is being tested; the profile figures are retained so "
                "the numbers in those issues can be checked against this code.",
            ]
        )
    return lines


def _report(result: ScreenResult, vintage: str | None) -> str:
    """A short human read of the screen: the six in context, and the strongest movers."""
    nominal = significant(result.rows, "by_q_nominal", 0.05)
    effective = significant(result.rows, "by_q_effective", 0.05)
    tested = [row for row in result.rows if row.test_status == TESTED]
    lines = [
        "# Directional consistency over the full retained trait space",
        "",
        f"Version: {VERSION}; {result.family_size} candidate traits, "
        f"{result.studies} studies, Benjamini-Yekutieli adjustment. "
        f"BugSigDB vintage: {vintage or 'unrecorded'}.",
        "",
        "Every retained non-MediaDive trait of the pinned matrix is tested, not just the "
        "six of Table 1. The statistic is unchanged: one vote per study per trait, the "
        "sign of the paired increased-minus-decreased contrast, ties excluded, two-sided "
        "exact binomial against 0.5. Only the family it is corrected in has changed, so "
        "the six traits' q-values here are not the six-trait q-values published elsewhere.",
        "",
        *_dependence_paragraphs(result),
        "",
        f"The reported result is BY q <= 0.05 at the nominal family size of "
        f"{result.family_size}: {len(nominal)} traits"
        + (
            f". The effective family size gives {len(effective)}, a sensitivity bound with "
            "no FDR guarantee (#265). The nominal set is contained in the effective one, so "
            "the difference is candidates gained, not candidates swapped. The count depends "
            "on the assumed dependence and is not stable; do not quote either number alone."
            if result.effective_family_size_used is not None
            else ". No effective family size was available, so no second count is given."
        ),
        "",
        "## The six Table 1 traits in the widened family",
        "",
        "| Trait | Up / down / tied studies | Share up (95% CI) | BY q (m = "
        f"{result.family_size}) | BY q (m = "
        f"{result.effective_family_size_used if result.effective_family_size_used else 'NA'}) "
        "| Status |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for row in result.rows:
        if not row.in_table1:
            continue
        manuscript = TABLE1_TRAIT_LABELS[row.table1_slug]
        lines.append(
            f"| {manuscript} | {row.studies_up}/{row.studies_down}/{row.studies_tied} | "
            f"{_cell(row.share_up, '.4g')} ({_cell(row.ci_low, '.4g')}, "
            f"{_cell(row.ci_high, '.4g')}) | {_cell(row.by_q_nominal, '.4g')} | "
            f"{_cell(row.by_q_effective, '.4g')} | {row.test_status} |"
        )
    strongest = sorted(
        (row for row in tested if row.p_value is not None),
        key=lambda row: (row.p_value if row.p_value is not None else 1.0, -row.studies_scored),
    )[:15]
    lines.extend(
        [
            "",
            "## The fifteen strongest candidates",
            "",
            "A display order, not a selection: every candidate is in the table, and no "
            "threshold was applied to build this list. Ordering is by p-value rather than "
            "by distance from an even split, because support varies by orders of magnitude "
            "across candidates and a lopsided share on a handful of studies is an "
            "extreme-value artefact rather than a finding. Scored is the number of non-tied "
            "studies, which is the sample size of that trait's test.",
            "",
            "| Trait | In Table 1 | Up / down / tied | Scored | Share up | BY q | BH q |",
            "| --- | --- | --- | --- | --- | --- | --- |",
        ]
    )
    for row in strongest:
        lines.append(
            f"| {row.trait_label} | {_cell(row.in_table1)} | "
            f"{row.studies_up}/{row.studies_down}/{row.studies_tied} | {row.studies_scored} | "
            f"{_cell(row.share_up, '.4g')} | {_cell(row.by_q_nominal, '.4g')} | "
            f"{_cell(row.bh_q_nominal, '.4g')} |"
        )
    lines.extend(
        [
            "",
            "This screen selects nothing and combines nothing. It measures relative "
            "annotation representation between increased and decreased signatures, not "
            "absolute trait abundance, and agreement with another arm on this shared input "
            "is not independent replication. Calibration is #182 and #263; the "
            "assertion-level sampling unit is #131.",
            "",
            "Complete statistics: [TSV](trait_direction_screen.tsv). "
            "Family size and dependence: [TSV](trait_screen_family.tsv). "
            "Eigenvalue spectrum: [TSV](trait_screen_correlation.tsv). "
            "Declared candidate manifest: "
            "[JSON](../../../methods/trait_screen_manifest.json). "
            "Source hashes and analysis rules: [JSON](trait_screen_provenance.json).",
        ]
    )
    return "\n".join(lines) + "\n"


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile-input", type=Path, default=DEFAULT_PROFILE_PATH)
    parser.add_argument("--experiments-input", type=Path, default=DEFAULT_EXPERIMENTS_PATH)
    parser.add_argument("--full-dump-input", type=Path, default=DEFAULT_FULL_DUMP_PATH)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--column-labels",
        type=Path,
        default=DEFAULT_COLUMN_LABELS_PATH,
        help="Trait matrix column index; candidates fall back to their CURIE without it",
    )
    parser.add_argument("--code-commit", help="Source revision, supplied by the Justfile recipe")
    args = parser.parse_args(argv)

    columns = screen_candidates(args.profile_input)
    slots = candidate_slots(columns)
    profiles, stats = load_profiles(
        args.profile_input, args.experiments_input, slots, reject_missing_counts=True
    )
    write_stats_report(stats, read_bugsigdb_vintage(args.full_dump_input))
    unpaired, _ = load_profiles(
        args.profile_input, args.experiments_input, slots, paired_only=False
    )
    result = run_screen(
        profiles,
        columns,
        dependence_profiles=unpaired,
        labels=load_column_labels(args.column_labels),
    )
    paths = write_screen_outputs(
        result,
        args.output_dir,
        profile_path=args.profile_input,
        experiments_path=args.experiments_input,
        full_dump_path=args.full_dump_input,
        corpus_stats=stats,
        code_commit=args.code_commit,
    )
    for path in paths:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
