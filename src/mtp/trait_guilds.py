"""Shared loading, filtering, and guild aggregation for trait guild figures.

Both the disease-by-body-site guild heatmap and the gram-stain Sankey start from
the same corpus definition, so the input filters agreed in issues #121, #122, and
#129 are applied here once rather than in each figure script.
"""

from __future__ import annotations

import csv
import sys
from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

from mtp.paths import DERIVED_DATA_DIR, RAW_DATA_DIR

DEFAULT_PROFILE_PATH = DERIVED_DATA_DIR / "study_trait_profiles_wide.tsv"
DEFAULT_EXPERIMENTS_PATH = DERIVED_DATA_DIR / "experiments.tsv"
DEFAULT_FULL_DUMP_PATH = RAW_DATA_DIR / "full_dump.csv"

PROFILE_ID_COLUMNS = ("study_id", "experiment_id", "inc_dec")
PROFILE_METADATA_COLUMNS = frozenset(
    {"study_id", "experiment_id", "inc_dec", "taxa_count", "taxon_count", "species_count"}
)
HUMAN_HOST_SPECIES = "Homo sapiens"
MEDIADIVE_MARKER = "mediadive.medium"

# Oxygen guilds, fixed a priori and identical to both analysis arms. Labels were
# verified against data/kg-microbe/merged-kg_nodes.tsv.
OXYGEN_GUILDS: dict[str, tuple[str, ...]] = {
    "aerobe": ("METPO:1000602", "METPO:1000604", "METPO:1000608"),
    "anaerobe": ("METPO:1000603", "METPO:1000607"),
}

# Gram stain classes. METPO:1000699 is absent from every script on the analysis
# branches; both labels were read from the KG-Microbe nodes file.
GRAM_GUILDS: dict[str, tuple[str, ...]] = {
    "gram positive": ("METPO:1000698",),
    "gram negative": ("METPO:1000699",),
}

# The six traits of manuscript Table 1, each its own single-member group rather
# than pooled into a guild. Slugs, display names, and colours are taken verbatim
# from Jed's six_shared_trait_selection.tsv on origin/clustered-heatmap-generation
# so the same trait is the same colour in both arms' manuscript figures.
# Two of the six are carried in the profile matrix as METPO-qualified columns
# (METPO:2000103|GO:0030245 and METPO:2000103|pathways:arsenate_reduction), which
# is why term lookup below also matches the segment after the final pipe.
TABLE1_TRAITS: dict[str, tuple[str, ...]] = {
    "facultatively_aerobic": ("METPO:1000608",),
    "cellulose_catabolic_process": ("GO:0030245",),
    "coccus_shaped": ("METPO:1000668",),
    "sphere_shaped": ("METPO:1000683",),
    "arsenate_reduction": ("pathways:arsenate_reduction",),
    "moderately_halophilic": ("METPO:1000623",),
}

TABLE1_TRAIT_LABELS: dict[str, str] = {
    "facultatively_aerobic": "facultatively aerobic (METPO:1000608)",
    "cellulose_catabolic_process": "cellulose catabolic process (GO:0030245)",
    "coccus_shaped": "coccus shaped (METPO:1000668)",
    "sphere_shaped": "sphere shaped (METPO:1000683)",
    "arsenate_reduction": "arsenate reduction (pathways:arsenate_reduction)",
    "moderately_halophilic": "moderately halophilic (METPO:1000623)",
}

TABLE1_TRAIT_COLORS: dict[str, str] = {
    "facultatively_aerobic": "#0072b2",
    "cellulose_catabolic_process": "#d55e00",
    "coccus_shaped": "#009e73",
    "sphere_shaped": "#cc79a7",
    "arsenate_reduction": "#e69f00",
    "moderately_halophilic": "#6a3d9a",
}

DIRECTIONS = ("inc", "dec")


@dataclass(frozen=True)
class Profile:
    """One study/experiment/direction row with its guild sums and row total."""

    study_id: str
    experiment_id: str
    inc_dec: str
    condition: str
    body_site: str
    guild_counts: dict[str, int]
    total: int

    def fraction(self, guild: str) -> float:
        """Return this guild's share of the profile's total trait signal."""
        if self.total <= 0:
            return 0.0
        return self.guild_counts.get(guild, 0) / self.total


@dataclass(frozen=True)
class CorpusStats:
    """Row counts dropped at each filter stage, for reporting and audit."""

    profile_rows: int
    experiments_seen: int
    unmatched_experiments: int
    dropped_non_human: int
    dropped_unpaired: int
    mediadive_columns_excluded: int
    trait_columns_used: int
    kept_experiments: int
    kept_profiles: int


class ExperimentMetadata:
    """Per-experiment condition, body site, and host species from experiments.tsv."""

    def __init__(self, rows: dict[tuple[str, str], dict[str, str]]) -> None:
        self._rows = rows

    @classmethod
    def load(cls, path: Path) -> ExperimentMetadata:
        required = {"study_id", "experiment_id", "condition", "body_site", "host_species"}
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = required.difference(reader.fieldnames or [])
            if missing:
                msg = (
                    f"{path} is missing {sorted(missing)}. Regenerate it with "
                    "`just study-tsvs` on a checkout that includes issue #142."
                )
                raise ValueError(msg)
            rows = {(r["study_id"], r["experiment_id"]): r for r in reader}
        return cls(rows)

    def get(self, study_id: str, experiment_id: str) -> dict[str, str] | None:
        return self._rows.get((study_id, experiment_id))


def read_bugsigdb_vintage(full_dump_path: Path) -> str | None:
    """Return the export timestamp recorded in the dump's leading comment line."""
    if not full_dump_path.exists():
        return None
    with full_dump_path.open(encoding="utf-8", errors="replace") as handle:
        first = handle.readline().strip()
    if not first.startswith("#"):
        return None
    body = first.lstrip("#").strip()
    head = body.split(",", 1)[0].strip()
    # The header reads "# BugSigDB <timestamp>, License: ..."; keep the timestamp.
    prefix = "BugSigDB "
    if head.startswith(prefix):
        head = head[len(prefix) :].strip()
    return head or None


def column_term_id(column: str) -> str:
    """Return the term a profile column stands for.

    Trait columns are either a bare CURIE (``METPO:1000608``) or a METPO-qualified
    pair (``METPO:2000103|GO:0030245``). The qualifier says how the term was
    attached, not which term it is, so lookups key on the segment after the final
    pipe. Only ``mediadive.medium:372`` shares a suffix with another column across
    the whole header (2,610 columns on the pinned v1.3.1 matrix), and MediaDive
    columns are dropped before any lookup runs.
    """
    return column.rsplit("|", 1)[-1]


def retained_trait_columns(header: Sequence[str]) -> tuple[list[int], list[str], int]:
    """The trait columns a profile matrix retains, in header order.

    One rule, in one place: drop the three identifier columns, drop every
    MediaDive growth-medium column (#121), keep the rest. What survives is both
    the denominator universe of every fraction in this package and the candidate
    space of the widened trait screen (#255), so the two cannot drift apart.

    Returns the retained positions, their column names, and how many MediaDive
    columns were excluded.
    """
    positions: list[int] = []
    columns: list[str] = []
    mediadive_excluded = 0
    for index, column in enumerate(header):
        if column in PROFILE_METADATA_COLUMNS:
            continue
        if MEDIADIVE_MARKER in column.lower():
            mediadive_excluded += 1
            continue
        positions.append(index)
        columns.append(column)
    return positions, columns, mediadive_excluded


def _select_trait_columns(
    header: list[str], guilds: dict[str, tuple[str, ...]]
) -> tuple[list[int], dict[str, list[int]], int]:
    """Map trait columns to positions, dropping MediaDive growth-medium columns."""
    wanted: dict[str, set[str]] = {g: set(ids) for g, ids in guilds.items()}
    guild_positions: dict[str, list[int]] = {guild: [] for guild in guilds}
    matched_columns: defaultdict[tuple[str, str], list[str]] = defaultdict(list)
    total_positions, total_columns, mediadive_excluded = retained_trait_columns(header)

    for index, column in zip(total_positions, total_columns, strict=True):
        term = column_term_id(column)
        for guild, ids in wanted.items():
            if column in ids or term in ids:
                guild_positions[guild].append(index)
                matched_columns[(guild, column if column in ids else term)].append(column)

    for guild, positions in guild_positions.items():
        if not positions:
            msg = f"No trait column in the profile matrix matches guild {guild!r}"
            raise ValueError(msg)

    # Orthogonal check: a term that resolved to more than one column would be
    # counted twice in the guild sum and once in the denominator.
    ambiguous = {key: cols for key, cols in matched_columns.items() if len(cols) > 1}
    if ambiguous:
        detail = "; ".join(
            f"{guild}/{term} -> {sorted(cols)}" for (guild, term), cols in sorted(ambiguous.items())
        )
        msg = f"Trait terms resolved to more than one profile column: {detail}"
        raise ValueError(msg)

    return total_positions, guild_positions, mediadive_excluded


def load_profiles(
    profile_path: Path = DEFAULT_PROFILE_PATH,
    experiments_path: Path = DEFAULT_EXPERIMENTS_PATH,
    guilds: dict[str, tuple[str, ...]] | None = None,
    *,
    human_only: bool = True,
    paired_only: bool = True,
    reject_missing_counts: bool = False,
) -> tuple[list[Profile], CorpusStats]:
    """Load profile rows joined to experiment metadata, with corpus filters applied."""
    guild_map = OXYGEN_GUILDS if guilds is None else guilds
    metadata = ExperimentMetadata.load(experiments_path)

    with profile_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if header is None:
            msg = f"{profile_path} is empty"
            raise ValueError(msg)
        if header[:3] != list(PROFILE_ID_COLUMNS):
            msg = f"{profile_path} must start with {PROFILE_ID_COLUMNS}"
            raise ValueError(msg)

        total_positions, guild_positions, mediadive_excluded = _select_trait_columns(
            header, guild_map
        )

        profile_rows = 0
        unmatched = 0
        dropped_non_human = 0
        candidates: list[Profile] = []

        for row in reader:
            if not row:
                continue
            profile_rows += 1
            study_id, experiment_id, inc_dec = row[0], row[1], row[2]
            if inc_dec not in DIRECTIONS:
                continue

            meta = metadata.get(study_id, experiment_id)
            if meta is None:
                unmatched += 1
                continue
            if human_only and meta["host_species"].strip() != HUMAN_HOST_SPECIES:
                dropped_non_human += 1
                continue

            if reject_missing_counts:
                for i in total_positions:
                    if i >= len(row) or not row[i].strip():
                        msg = (
                            f"Missing trait count for {study_id}, {experiment_id}, {inc_dec}: "
                            f"{header[i]}; missing is not zero"
                        )
                        raise ValueError(msg)
            values = [int(row[i]) if row[i].strip() else 0 for i in total_positions]
            guild_counts = {
                guild: sum(int(row[i]) if row[i].strip() else 0 for i in positions)
                for guild, positions in guild_positions.items()
            }
            candidates.append(
                Profile(
                    study_id=study_id,
                    experiment_id=experiment_id,
                    inc_dec=inc_dec,
                    condition=meta["condition"].strip(),
                    body_site=meta["body_site"].strip(),
                    guild_counts=guild_counts,
                    total=sum(values),
                )
            )

    kept, dropped_unpaired = _apply_pair_filter(candidates, paired_only=paired_only)
    experiments_seen = len({(p.study_id, p.experiment_id) for p in candidates})
    stats = CorpusStats(
        profile_rows=profile_rows,
        experiments_seen=experiments_seen,
        unmatched_experiments=unmatched,
        dropped_non_human=dropped_non_human,
        dropped_unpaired=dropped_unpaired,
        mediadive_columns_excluded=mediadive_excluded,
        trait_columns_used=len(total_positions),
        kept_experiments=len({(p.study_id, p.experiment_id) for p in kept}),
        kept_profiles=len(kept),
    )
    return kept, stats


def _apply_pair_filter(profiles: list[Profile], *, paired_only: bool) -> tuple[list[Profile], int]:
    """Keep only experiments carrying both an inc and a dec profile (#122)."""
    if not paired_only:
        return profiles, 0
    by_experiment: defaultdict[tuple[str, str], set[str]] = defaultdict(set)
    for profile in profiles:
        by_experiment[(profile.study_id, profile.experiment_id)].add(profile.inc_dec)
    matched = {key for key, seen in by_experiment.items() if seen >= set(DIRECTIONS)}
    kept = [p for p in profiles if (p.study_id, p.experiment_id) in matched]
    return kept, len(profiles) - len(kept)


def format_stats_report(stats: CorpusStats, vintage: str | None) -> str:
    """Render the corpus filter cascade so no reduction is silent."""
    lines = [
        f"BugSigDB vintage:            {vintage or 'unknown'}",
        f"Profile rows read:           {stats.profile_rows}",
        f"Experiments after join:      {stats.experiments_seen}",
        f"Dropped, no metadata match:  {stats.unmatched_experiments}",
        f"Dropped, non-human host:     {stats.dropped_non_human}",
        f"Dropped, unpaired:           {stats.dropped_unpaired}",
        f"MediaDive columns excluded:  {stats.mediadive_columns_excluded}",
        f"Trait columns used:          {stats.trait_columns_used}",
        f"Experiments kept:            {stats.kept_experiments}",
        f"Profiles kept:               {stats.kept_profiles}",
    ]
    return "\n".join(lines)


def write_stats_report(
    stats: CorpusStats, vintage: str | None, stream: TextIO = sys.stdout
) -> None:
    """Print the corpus filter cascade so no reduction is silent."""
    print(format_stats_report(stats, vintage), file=stream)
