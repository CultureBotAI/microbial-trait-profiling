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
from dataclasses import asdict, dataclass, field
from pathlib import Path

from mtp.paths import KG_MICROBE_DIR, file_sha256
from mtp.trait_guilds import MEDIADIVE_MARKER, PROFILE_ID_COLUMNS

DEFAULT_KG_PROVENANCE_PATH = KG_MICROBE_DIR / "provenance.json"
#: What a release tag reads as when the sidecar does not carry one. Recorded
#: rather than omitted, because "not known" and "not applicable" differ.
UNKNOWN_RELEASE = "unknown: the KG-Microbe download records no release tag"

__all__ = [
    "DEFAULT_KG_PROVENANCE_PATH",
    "UNKNOWN_RELEASE",
    "CorpusIdentity",
    "KgInput",
    "MatrixIdentity",
    "corpus_identity",
    "measure_matrix",
    "read_kg_inputs",
    "summary_lines",
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
                "inputs": [asdict(item) for item in self.kg_inputs],
            },
            "note": (
                "A matrix is a join of the BugSigDB dump with a KG-Microbe build. "
                "Both sides are recorded because the vintage alone does not identify it."
            ),
        }


def measure_matrix(profile_path: Path) -> MatrixIdentity:
    """Read ``profile_path`` once for its shape, and hash it."""
    identifier = mediadive = retained = 0
    rows = 0
    with profile_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        for column in header:
            if column in PROFILE_ID_COLUMNS:
                identifier += 1
            elif MEDIADIVE_MARKER in column.lower():
                mediadive += 1
            else:
                retained += 1
        for _ in reader:
            rows += 1
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

    tag = recorded.get("release") or recorded.get("release_tag") or recorded.get("version")
    release = tag if isinstance(tag, str) and tag else UNKNOWN_RELEASE

    files = recorded.get("files")
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
        f"KG-Microbe:        {identity.kg_release_tag}",
    ]
    lines.extend(
        f"  {item.name:<22} md5 {item.declared_md5 or '(not recorded)'}"
        f"{'' if item.present else '  [file not present]'}"
        for item in identity.kg_inputs
    )
    return lines
