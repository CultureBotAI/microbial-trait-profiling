"""Which corpus a run actually read, recorded so a number can be traced to it (#222).

"The pinned v1.3.1 corpus" names only the BugSigDB dump. The matrix every arm
consumes is a *join* of that dump with a KG-Microbe build, and the KG side is
pinned by Drive file id and checksum but carries no release tag. Two runs can
therefore both be "on the pinned corpus" and read different matrices, which is
what happened between the 2,290 and 2,333 MediaDive column counts: the retained
trait set was 317 in both, and the difference lay entirely in growth-medium
columns that every arm discards before computing a total.

So that a reported count can be traced back to the file it came from, this module
records the identity of both sides - the matrix itself and the KG-Microbe inputs
behind it - in one record that every recipe writes into its provenance JSON.

What it deliberately does not do is decide which corpus is correct. Whether the
shared corpus is the committed ``.gz`` or a regeneration is a project-wide choice
(#115), and until it is settled the honest thing is to say precisely which file a
given number came from.
"""

from __future__ import annotations

import csv
import json
from collections.abc import Mapping
from dataclasses import asdict, dataclass, field
from pathlib import Path

from mtp.paths import KG_MICROBE_DIR, file_sha256
from mtp.trait_guilds import MEDIADIVE_MARKER, PROFILE_ID_COLUMNS, retained_trait_columns

DEFAULT_KG_PROVENANCE_PATH = KG_MICROBE_DIR / "provenance.json"
#: Pinned KG-Microbe build folder used by scripts/download_kg_microbe.py (#336).
PINNED_KG_MICROBE_RELEASE = "20260120"
_PINNED_20260120_FILES: dict[str, tuple[str, str]] = {
    "merged-kg_nodes.tsv": (
        "1lDVjE67E82LMUdNtqWGvyczrHvNu3BdY",  # pragma: allowlist secret
        "dd76b4656f9341db25a7c4afa96cd58e",  # pragma: allowlist secret
    ),
    "merged-kg_edges.tsv": (
        "1iIVd-FEleeRIv7CPS0l1F732p3J9JO5d",  # pragma: allowlist secret
        "555ecf5ba5d89288cdc44e368375374f",  # pragma: allowlist secret
    ),
}
#: What a release tag reads as when the sidecar does not carry one. Recorded
#: rather than omitted, because "not known" and "not applicable" differ.
UNKNOWN_RELEASE = "unknown: the KG-Microbe download records no release tag"

__all__ = [
    "DEFAULT_KG_PROVENANCE_PATH",
    "PINNED_KG_MICROBE_RELEASE",
    "UNKNOWN_RELEASE",
    "CorpusIdentity",
    "KgInput",
    "MatrixIdentity",
    "corpus_identity",
    "measure_matrix",
    "read_kg_inputs",
    "summary_lines",
    "verify_provenance_source_hashes",
    "write_identity",
]


@dataclass(frozen=True)
class MatrixIdentity:
    """The trait matrix a run read, by content rather than by path."""

    path: str
    sha256: str
    bytes: int
    header_columns: int
    identifier_columns: int
    mediadive_columns: int
    retained_trait_columns: int
    profile_rows: int

    def check(self) -> None:
        """The three column groups must account for the whole header."""
        total = self.identifier_columns + self.mediadive_columns + self.retained_trait_columns
        if total != self.header_columns:
            msg = (
                f"{self.path}: {self.identifier_columns} identifier + "
                f"{self.mediadive_columns} MediaDive + {self.retained_trait_columns} retained "
                f"= {total}, not the {self.header_columns} columns in the header"
            )
            raise ValueError(msg)
        if self.identifier_columns != len(PROFILE_ID_COLUMNS):
            msg = (
                f"{self.path}: expected {len(PROFILE_ID_COLUMNS)} identifier columns "
                f"{PROFILE_ID_COLUMNS}, got {self.identifier_columns}"
            )
            raise ValueError(msg)


@dataclass(frozen=True)
class KgInput:
    """One KG-Microbe file as the download sidecar describes it."""

    name: str
    declared_md5: str | None
    declared_size: int | None
    fetched_at: str | None
    drive_file_id: str | None
    present: bool


@dataclass(frozen=True)
class CorpusIdentity:
    """Both sides of the join, plus the BugSigDB pin."""

    matrix: MatrixIdentity
    bugsigdb_vintage: str | None
    kg_release_tag: str
    kg_provenance_path: str
    kg_inputs: list[KgInput] = field(default_factory=list)

    def as_record(self) -> dict[str, object]:
        return {
            "matrix": asdict(self.matrix),
            "bugsigdb_vintage": self.bugsigdb_vintage,
            "kg_microbe": {
                "release_tag": self.kg_release_tag,
                "provenance_path": self.kg_provenance_path,
                "provenance_source": (
                    "declared by download sidecar (not measured at matrix build time)"
                ),
                "inputs": [asdict(item) for item in self.kg_inputs],
            },
            "note": (
                "A matrix is a join of the BugSigDB dump with a KG-Microbe build. "
                "Both sides are recorded because the vintage alone does not identify it."
            ),
        }


def measure_matrix(profile_path: Path) -> MatrixIdentity:
    """Read ``profile_path`` once for its shape, and hash it."""
    rows = 0
    with profile_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        for row in reader:
            if row:
                rows += 1
    unique_header = set(header)
    identifier = sum(1 for column in unique_header if column in PROFILE_ID_COLUMNS)
    _, retained_cols, _ = retained_trait_columns(header)
    retained = len(set(retained_cols))
    mediadive = sum(
        1
        for column in unique_header
        if column not in PROFILE_ID_COLUMNS and MEDIADIVE_MARKER in column.lower()
    )
    identity = MatrixIdentity(
        path=str(profile_path),
        sha256=file_sha256(profile_path),
        bytes=profile_path.stat().st_size,
        header_columns=len(header),
        identifier_columns=identifier,
        mediadive_columns=mediadive,
        retained_trait_columns=retained,
        profile_rows=rows,
    )
    identity.check()
    return identity


def _matches_pinned_20260120(files: dict[str, object]) -> bool:
    """True when the sidecar's file IDs or MD5s match the pinned 20260120 KG build."""
    if not _PINNED_20260120_FILES.keys() <= files.keys():
        return False
    for name, (expected_id, expected_md5) in _PINNED_20260120_FILES.items():
        entry = files.get(name)
        if not isinstance(entry, dict):
            return False
        if entry.get("file_id") != expected_id and entry.get("md5") != expected_md5:
            return False
    return True


def read_kg_inputs(provenance_path: Path) -> tuple[str, list[KgInput]]:
    """The KG-Microbe files the sidecar describes, and its release tag if any."""
    if not provenance_path.is_file():
        return f"{UNKNOWN_RELEASE} ({provenance_path} is absent)", []
    try:
        recorded = json.loads(provenance_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return f"{UNKNOWN_RELEASE} ({provenance_path} is not valid JSON)", []
    if not isinstance(recorded, dict):
        return f"{UNKNOWN_RELEASE} ({provenance_path} is not an object)", []

    files = recorded.get("files")
    tag = recorded.get("release") or recorded.get("release_tag") or recorded.get("version")
    if isinstance(tag, str) and tag:
        release = tag
    elif isinstance(files, dict) and _matches_pinned_20260120(files):
        release = PINNED_KG_MICROBE_RELEASE
    else:
        release = UNKNOWN_RELEASE

    inputs: list[KgInput] = []
    if isinstance(files, dict):
        for name, entry in sorted(files.items()):
            details = entry if isinstance(entry, dict) else {}
            inputs.append(
                KgInput(
                    name=name,
                    declared_md5=_as_str(details.get("md5")),
                    declared_size=_as_int(details.get("size")),
                    fetched_at=_as_str(details.get("fetched_at")),
                    drive_file_id=_as_str(details.get("file_id")),
                    present=(provenance_path.parent / name).is_file(),
                )
            )
    return release, inputs


def _as_str(value: object) -> str | None:
    return value if isinstance(value, str) else None


def _as_int(value: object) -> int | None:
    return value if isinstance(value, int) else None


def corpus_identity(
    profile_path: Path,
    *,
    bugsigdb_vintage: str | None,
    kg_provenance_path: Path = DEFAULT_KG_PROVENANCE_PATH,
) -> CorpusIdentity:
    """The full identity record a recipe writes into its provenance JSON."""
    release, inputs = read_kg_inputs(kg_provenance_path)
    return CorpusIdentity(
        matrix=measure_matrix(profile_path),
        bugsigdb_vintage=bugsigdb_vintage,
        kg_release_tag=release,
        kg_provenance_path=str(kg_provenance_path),
        kg_inputs=inputs,
    )


def write_identity(
    identity: CorpusIdentity,
    output_dir: Path,
    *,
    filename: str = "corpus_identity.json",
    code_commit: str | None = None,
) -> Path:
    """Write the identity record beside a recipe's other outputs.

    Recipes that already build a provenance JSON embed the record instead; this
    is for the figure recipes, which wrote no provenance at all, so a figure's
    numbers could not be traced to the matrix behind them.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / filename
    record = {**identity.as_record(), "code_commit": code_commit}
    path.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    return path


def summary_lines(identity: CorpusIdentity) -> list[str]:
    """The identity as a recipe prints it, so a run states its own corpus."""
    matrix = identity.matrix
    lines = [
        f"Matrix:            {matrix.path}",
        f"  sha256           {matrix.sha256}",
        f"  columns          {matrix.header_columns} "
        f"({matrix.identifier_columns} identifier, {matrix.mediadive_columns} MediaDive, "
        f"{matrix.retained_trait_columns} retained)",
        f"  profile rows     {matrix.profile_rows}",
        f"BugSigDB vintage:  {identity.bugsigdb_vintage or 'unknown'}",
        f"KG-Microbe:        {identity.kg_release_tag} (sidecar: {identity.kg_provenance_path})",
    ]
    lines.extend(
        f"  {item.name:<22} declared md5 (sidecar) {item.declared_md5 or '(not recorded)'}"
        f"{'' if item.present else '  [file not present]'}"
        for item in identity.kg_inputs
    )
    return lines


def verify_provenance_source_hashes(
    record: Mapping[str, object],
    repo_root: Path | None = None,
) -> dict[str, dict[str, str]]:
    """Return any recorded source-file SHA-256 entries that do not match disk (#499).

    Inspects ``record["source_hashes"]`` (mapping relative paths to SHA-256 digests)
    as well as top-level ``justfile_sha256`` if present. Returns a mapping of
    relative path -> ``{"recorded": ..., "actual": ...}`` for every mismatched or
    missing file (empty dict when all recorded source hashes are fresh).
    """
    root = repo_root if repo_root is not None else Path(__file__).resolve().parents[2]
    candidates: dict[str, str] = {}
    source_hashes = record.get("source_hashes")
    if isinstance(source_hashes, Mapping):
        for rel_path, digest in source_hashes.items():
            if isinstance(rel_path, str) and isinstance(digest, str):
                candidates[rel_path] = digest
    justfile_digest = record.get("justfile_sha256")
    if isinstance(justfile_digest, str) and "Justfile" not in candidates:
        candidates["Justfile"] = justfile_digest

    mismatches: dict[str, dict[str, str]] = {}
    for rel_path, recorded_sha in candidates.items():
        target = root / rel_path
        actual_sha = file_sha256(target) if target.is_file() else "(missing)"
        if actual_sha != recorded_sha:
            mismatches[rel_path] = {
                "recorded": recorded_sha,
                "actual": actual_sha,
            }
    return mismatches
