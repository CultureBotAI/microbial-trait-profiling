"""Shared loading, filtering, and guild aggregation for trait guild figures.

Both the disease-by-body-site guild heatmap and the gram-stain Sankey start from
the same corpus definition, so the input filters agreed in issues #121, #122, and
#129 are applied here once rather than in each figure script.
"""

from __future__ import annotations

import csv
import sys
from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

from mtp.count_input import (
    COMPLETE_COUNT_MODE,
    LEGACY_SPARSE_ZERO_MODE,
    CountInputMode,
    open_input_text,
    profile_count_reader,
)
from mtp.paths import DERIVED_DATA_DIR, RAW_DATA_DIR
from mtp.trait_key_policy import display_trait_label, read_recorded_member_ids

DEFAULT_PROFILE_PATH = DERIVED_DATA_DIR / "study_trait_profiles_wide.tsv"
DEFAULT_EXPERIMENTS_PATH = DERIVED_DATA_DIR / "experiments.tsv"
DEFAULT_FULL_DUMP_PATH = RAW_DATA_DIR / "full_dump.csv"

PROFILE_ID_COLUMNS = ("study_id", "experiment_id", "inc_dec")
PROFILE_METADATA_COLUMNS = frozenset(
    {"study_id", "experiment_id", "inc_dec", "taxa_count", "taxon_count", "species_count"}
)
HUMAN_HOST_SPECIES = "Homo sapiens"
MEDIADIVE_MARKER = "mediadive.medium"
# How count_input.parse_nonnegative_integer starts a blank-count error.
_READER_MISSING_COUNT_PREFIX = "Missing count in "

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

# The six traits named in manuscript Table 1, each its own single-member group rather
# than pooled into a guild. Used in the matched-corpus Table 1-trait analog
# sensitivity analysis (not a reproduction of published Table 1 baseline values).
# Slugs, display names, and colours are taken verbatim from Jed's
# six_shared_trait_selection.tsv on origin/clustered-heatmap-generation so the
# same trait is the same colour in both arms' manuscript figures.
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


def table1_trait_label(
    trait_name: str,
    column_metadata_path: Path | None = None,
    member_ids_by_trait: Mapping[str, Sequence[str] | str] | None = None,
) -> str:
    """Return the display label for a Table 1 trait slug from recorded member IDs."""
    base_label = TABLE1_TRAIT_LABELS[trait_name]
    terms = TABLE1_TRAITS.get(trait_name, ())
    if not terms:
        return base_label
    term_id = terms[0]
    members = (
        dict(member_ids_by_trait)
        if member_ids_by_trait is not None
        else read_recorded_member_ids(column_metadata_path)
    )
    return display_trait_label(term_id, base_label, members.get(term_id))


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
        """Return this guild's share of the profile's total trait signal.

        For paired contrasts across experiments, callers must filter zero-total
        pairs via ``drop_zero_total_pairs`` to satisfy the shared comparison
        basis of methods/three_arm_statistics.md (#228).
        """
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


@dataclass(frozen=True)
class PairExclusion:
    """Experiments kept for a paired contrast, and what was removed."""

    profiles: list[Profile]
    zero_total_experiments: int
    studies_before: int
    studies_after: int

    @property
    def studies_lost(self) -> int:
        return self.studies_before - self.studies_after


def drop_zero_total_pairs(profiles: Sequence[Profile]) -> PairExclusion:
    """Remove whole experiments whose retained-annotation total is zero on a side.

    ``Profile.fraction`` returns 0.0 when the total is zero, so such a profile
    enters a contrast as a measured share of nothing. It is not: the share is
    undefined. Scoring it as zero makes the contrast exactly zero, which the
    sign test then counts as a tie, so all-tied studies and zero-annotation
    studies collapse into one number - the collapse the shared comparison basis
    of #178 forbids (#294).

    The whole experiment goes, not the offending side, so the corpus stays
    paired. On the pinned corpus this removes 55 of 2,257 experiments and takes
    the study count from 932 to 922.
    """
    by_experiment: defaultdict[tuple[str, str], list[Profile]] = defaultdict(list)
    for profile in profiles:
        by_experiment[(profile.study_id, profile.experiment_id)].append(profile)
    kept: list[Profile] = []
    zero_total = 0
    for sides in by_experiment.values():
        if any(side.total <= 0 for side in sides):
            zero_total += 1
        else:
            kept.extend(sides)
    return PairExclusion(
        profiles=kept,
        zero_total_experiments=zero_total,
        studies_before=len({p.study_id for p in profiles}),
        studies_after=len({p.study_id for p in kept}),
    )


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
    the whole header (2,653 columns on the all-host v1.3.1 matrix; 2,610 on the
    canonical human table, the committed ``.gz`` (#656)), and MediaDive columns are
    dropped before any lookup runs.
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


@dataclass(frozen=True)
class TraitColumnSelection:
    """Where a guild's declared terms were found in the profile matrix header."""

    total_positions: list[int]
    guild_positions: dict[str, list[int]]
    mediadive_excluded: int
    #: guild -> ((declared term id, raw matrix column), ...), in header order. The
    #: raw column is the only key a result can be joined back to the matrix on.
    matched_columns: dict[str, tuple[tuple[str, str], ...]]


def _select_trait_columns(
    header: list[str], guilds: dict[str, tuple[str, ...]]
) -> TraitColumnSelection:
    """Map trait columns to positions, dropping MediaDive growth-medium columns."""
    wanted: dict[str, set[str]] = {g: set(ids) for g, ids in guilds.items()}
    guild_positions: dict[str, list[int]] = {guild: [] for guild in guilds}
    matched: defaultdict[str, list[tuple[str, str]]] = defaultdict(list)
    total_positions, total_columns, mediadive_excluded = retained_trait_columns(header)

    for index, column in zip(total_positions, total_columns, strict=True):
        term = column_term_id(column)
        for guild, ids in wanted.items():
            if column in ids or term in ids:
                guild_positions[guild].append(index)
                matched[guild].append((column if column in ids else term, column))

    for guild, positions in guild_positions.items():
        if not positions:
            msg = f"No trait column in the profile matrix matches guild {guild!r}"
            raise ValueError(msg)

    # Orthogonal check: a term that resolved to more than one column would be
    # counted twice in the guild sum and once in the denominator.
    per_term: defaultdict[tuple[str, str], list[str]] = defaultdict(list)
    for guild, pairs in matched.items():
        for term, column in pairs:
            per_term[(guild, term)].append(column)
    ambiguous = {key: cols for key, cols in per_term.items() if len(cols) > 1}
    if ambiguous:
        detail = "; ".join(
            f"{guild}/{term} -> {sorted(cols)}" for (guild, term), cols in sorted(ambiguous.items())
        )
        msg = f"Trait terms resolved to more than one profile column: {detail}"
        raise ValueError(msg)

    return TraitColumnSelection(
        total_positions=total_positions,
        guild_positions=guild_positions,
        mediadive_excluded=mediadive_excluded,
        matched_columns={guild: tuple(pairs) for guild, pairs in sorted(matched.items())},
    )


def _read_profile_header(profile_path: Path) -> list[str]:
    """Return the profile matrix header, rejecting a file that is not one."""
    with open_input_text(profile_path) as handle:
        header = next(csv.reader(handle, delimiter="\t"), None)
    if header is None:
        msg = f"{profile_path} is empty"
        raise ValueError(msg)
    if header[:3] != list(PROFILE_ID_COLUMNS):
        msg = f"{profile_path} must start with {PROFILE_ID_COLUMNS}"
        raise ValueError(msg)
    return header


def resolve_trait_columns(
    profile_path: Path = DEFAULT_PROFILE_PATH,
    guilds: dict[str, tuple[str, ...]] | None = None,
) -> dict[str, tuple[tuple[str, str], ...]]:
    """Return the (declared term, raw matrix column) pairs each guild resolved to.

    A consumer that needs to join a result back to the matrix needs the raw
    column name, not the declared CURIE. Two of the six Table 1 traits are
    carried as METPO-qualified composites (``METPO:2000103|GO:0030245`` and
    ``METPO:2000103|pathways:arsenate_reduction``), so the declared term matches
    the segment after the final pipe rather than the whole column, and a result
    exported under the bare term cannot be joined to the matrix at all (#294).
    """
    header = _read_profile_header(profile_path)
    selection = _select_trait_columns(header, OXYGEN_GUILDS if guilds is None else guilds)
    return selection.matched_columns


def load_profiles(
    profile_path: Path = DEFAULT_PROFILE_PATH,
    experiments_path: Path = DEFAULT_EXPERIMENTS_PATH,
    guilds: dict[str, tuple[str, ...]] | None = None,
    *,
    human_only: bool = True,
    paired_only: bool = True,
    count_input_mode: CountInputMode = COMPLETE_COUNT_MODE,
    reject_missing_counts: bool | None = None,
) -> tuple[list[Profile], CorpusStats]:
    """Load profile rows joined to experiment metadata, with corpus filters applied."""
    guild_map = OXYGEN_GUILDS if guilds is None else guilds
    metadata = ExperimentMetadata.load(experiments_path)
    effective_mode: CountInputMode = (
        count_input_mode
        if reject_missing_counts is None
        else (COMPLETE_COUNT_MODE if reject_missing_counts else LEGACY_SPARSE_ZERO_MODE)
    )

    count_reader = profile_count_reader(profile_path, mode=effective_mode)
    header = list(count_reader.header)
    selection = _select_trait_columns(header, guild_map)
    total_positions = selection.total_positions
    guild_positions = selection.guild_positions
    mediadive_excluded = selection.mediadive_excluded

    profile_rows = 0
    unmatched = 0
    dropped_non_human = 0
    candidates: list[Profile] = []
    count_offset = len(PROFILE_ID_COLUMNS)

    try:
        for parsed_row in count_reader:
            profile_rows += 1
            study_id, experiment_id, inc_dec = parsed_row.key
            if inc_dec not in DIRECTIONS:
                msg = (
                    f"Unknown profile direction in {profile_path}:{parsed_row.line_number}: "
                    f"{inc_dec!r}; expected one of {DIRECTIONS}"
                )
                raise ValueError(msg)

            meta = metadata.get(study_id, experiment_id)
            if meta is None:
                unmatched += 1
                continue
            if human_only and meta["host_species"].strip() != HUMAN_HOST_SPECIES:
                dropped_non_human += 1
                continue

            values = [parsed_row.counts[i - count_offset] for i in total_positions]
            guild_counts = {
                guild: sum(parsed_row.counts[i - count_offset] for i in positions)
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
    except ValueError as exc:
        # Relabel only the strict reader's missing-count error, matched as a
        # prefix: a path that contains the phrase must not relabel any other
        # error (#1300). The reader's file:line:column and its statement that
        # missing is not an observed zero are kept as they are.
        message = str(exc)
        if message.startswith(_READER_MISSING_COUNT_PREFIX):
            detail = message.removeprefix(_READER_MISSING_COUNT_PREFIX)
            raise ValueError(f"Missing trait count in {detail}") from exc
        raise

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


def format_pair_exclusion(exclusion: PairExclusion) -> str:
    """Render the zero-total pair exclusion in one line, counts and rule together.

    Every arm that applies the shared comparison basis reports the same reduction,
    so the wording lives here rather than being retyped per caller (#294).
    """
    return (
        f"Excluded {exclusion.zero_total_experiments} experiment pairs with a zero "
        f"retained-annotation total on a side; studies {exclusion.studies_before} "
        f"-> {exclusion.studies_after} ({exclusion.studies_lost} lost). The share of "
        "a zero total is undefined, not zero, so the whole pair is dropped rather "
        "than entering the contrast as a tie."
    )
