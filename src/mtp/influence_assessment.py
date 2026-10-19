"""Adapt #190 execution evidence without manufacturing a scientific verdict (#671).

Native manifests describe a diagnostic run, not an assessment of every procedure
in the integrated table. Their completion and review handoff therefore remain
separate from explicitly supplied, procedure-specific assessment verdicts.
"""

from __future__ import annotations

import csv
import gzip
import json
import re
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

from mtp.paths import file_sha256

Detail = tuple[str, str, str, str]
ASSESSMENT_VERSION = "taxon-influence-assessment-v1"
NATIVE_VERSION = "190-native-v2"
INPUT_NAMES = ("study_trait_profiles_wide.tsv", "experiments.tsv")


@dataclass(frozen=True)
class InfluenceAssessment:
    """Procedure standing and separately scoped execution/provenance detail."""

    verdicts: dict[str, str]
    default: str
    detail: tuple[Detail, ...]


def _object(value: object, name: str) -> dict[str, object]:
    if not isinstance(value, dict) or any(not isinstance(key, str) for key in value):
        raise ValueError(f"{name} must be an object")
    return dict(value)


def _text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a nonempty string")
    return value


def _integer(value: object, name: str, minimum: int = 0) -> int:
    if type(value) is not int or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return value


def _hash(value: object, name: str) -> str:
    text = _text(value, name)
    if re.fullmatch(r"[0-9a-f]{64}", text) is None:
        raise ValueError(f"{name} must be a SHA-256 digest")
    return text


def _hashes(value: object, name: str) -> dict[str, str]:
    hashes = _object(value, name)
    if not hashes:
        raise ValueError(f"{name} must not be empty")
    return {key: _hash(item, f"{name}.{key}") for key, item in hashes.items()}


def _read(path: Path) -> dict[str, object]:
    return _object(json.loads(path.read_text(encoding="utf-8")), path.name)


def _rows(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else Path.open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t", strict=True)
        header = reader.fieldnames
        if not header or len(header) != len(set(header)) or any(not key for key in header):
            raise ValueError(f"invalid table header: {path.name}")
        rows: list[dict[str, str]] = []
        for row in reader:
            if any(
                not isinstance(key, str) or not isinstance(value, str) for key, value in row.items()
            ):
                raise ValueError(f"malformed table row: {path.name}")
            rows.append(dict(row))
    return rows


def _artifact(root: Path, name: str, expected: str) -> Path:
    path = (root / name).resolve()
    if Path(name).is_absolute() or not path.is_relative_to(root.resolve()):
        raise ValueError(f"artifact escapes the evidence directory: {name}")
    if not path.is_file() or file_sha256(path) != expected:
        raise ValueError(f"missing or changed evidence artifact: {name}")
    return path


def _detail(kind: str, scope: str, value: Mapping[str, object]) -> list[Detail]:
    # JSON values preserve nested protocol fields and booleans without guessing
    # that one of them is a procedure-level verdict.
    return [
        (kind, scope, key, item if isinstance(item, str) else json.dumps(item, sort_keys=True))
        for key, item in sorted(value.items())
    ]


def _inputs(content: Mapping[str, object], expected: Mapping[str, str] | None) -> dict[str, str]:
    hashes = _hashes(content.get("input_sha256"), "input_sha256")
    for name in INPUT_NAMES:
        if name not in hashes:
            raise ValueError(f"input_sha256 lacks {name}")
    for name, digest in (expected or {}).items():
        if hashes.get(name) != digest:
            raise ValueError(f"influence evidence belongs to different input bytes: {name}")
    return hashes


def _explicit(
    content: Mapping[str, object], issue: str, expected: Mapping[str, str] | None
) -> InfluenceAssessment:
    hashes = _inputs(content, expected)
    raw = _object(content.get("procedure_verdicts"), "procedure_verdicts")
    if not raw:
        raise ValueError("procedure_verdicts must not be empty")
    detail = _detail("input", "assessment", hashes)
    verdicts = {}
    for procedure, value in sorted(raw.items()):
        _text(procedure, "procedure")
        entry = _object(value, f"procedure_verdicts.{procedure}")
        verdict = _text(entry.get("verdict"), "verdict")
        _text(entry.get("basis"), "basis")
        if entry.get("applies_to") != "all_rows_of_procedure":
            raise ValueError(
                "procedure verdict must explicitly apply to all_rows_of_procedure; "
                "a scenario-specific finding cannot stand for an entire procedure"
            )
        verdicts[procedure] = f"{verdict} {issue}"
        detail.extend(_detail("procedure", procedure, entry))
    if "native_manifest_sha256" in content:
        digest = _hash(content["native_manifest_sha256"], "native_manifest_sha256")
        detail.append(("source", "assessment", "declared_native_manifest_sha256", digest))
    return InfluenceAssessment(
        verdicts, f"no procedure-level influence verdict supplied {issue}", tuple(detail)
    )


def _scope_detail(path: Path, results: Mapping[str, object]) -> list[Detail]:
    rows = _rows(path)
    identities = [
        (row.get("scenario", ""), row.get("scope", ""), row.get("background", "")) for row in rows
    ]
    if not identities or any(not all(identity) for identity in identities):
        raise ValueError("native_ora_families.tsv lacks scenario/scope/background identities")
    if len(set(identities)) != len(identities):
        raise ValueError("native_ora_families.tsv contains duplicate family identities")
    scenarios = {identity[0] for identity in identities}
    scopes = {identity[1] for identity in identities}
    backgrounds = {identity[2] for identity in identities}
    if (
        "baseline" not in scenarios
        or len(scenarios) != _integer(results.get("scenarios"), "results.scenarios", 1)
        or len(scopes) != _integer(results.get("scopes"), "results.scopes", 1)
        or len(rows) != _integer(results.get("ora_families"), "results.ora_families", 1)
        or backgrounds != {"fixed_baseline", "recomputed_remaining"}
        or len(rows) != len(scenarios) * len(scopes) * len(backgrounds)
    ):
        raise ValueError("native family scope/scenario inventory disagrees with the manifest")
    detail: list[Detail] = []
    for row, identity in zip(rows, identities, strict=True):
        scope = json.dumps(identity)
        detail.extend(_detail("native_ora_family", scope, row))
    return detail


def _review(
    path: Path, native_path: Path, native: Mapping[str, object]
) -> tuple[list[Detail], Path]:
    review = _read(path)
    manifest_hash = _hash(review.get("native_manifest_sha256"), "native_manifest_sha256")
    if manifest_hash != file_sha256(native_path):
        raise ValueError("publication manifest belongs to a different native manifest")
    _hash(review.get("publisher_sha256"), "publisher_sha256")
    reconciliation = _object(review.get("source_reconciliation"), "source_reconciliation")
    if _integer(reconciliation.get("mismatches"), "source_reconciliation.mismatches") != 0:
        raise ValueError("publication source reconciliation reports mismatches")
    _hashes(reconciliation.get("sha256"), "source_reconciliation.sha256")
    for name in ("source_contribution_rows", "profile_reported_taxon_rows", "profile_count_rows"):
        _integer(reconciliation.get(name), f"source_reconciliation.{name}", 1)
    _hashes(review.get("interpretation_source_sha256"), "interpretation_source_sha256")
    tables = _object(review.get("tables"), "tables")
    required = {
        "native_table1_stats.tsv.gz",
        "native_composition_comparison.tsv.gz",
        "native_ora_fold_comparison.tsv.gz",
        "native_ora_families.tsv",
        "scenarios.tsv",
        "interpretation_handoff.tsv",
    }
    if not required.issubset(tables):
        raise ValueError("publication manifest lacks required review tables")
    for name, item in tables.items():
        entry = _object(item, f"tables.{name}")
        row_count = _integer(entry.get("rows"), f"tables.{name}.rows", 1)
        table = _artifact(path.parent, name, _hash(entry.get("sha256"), f"tables.{name}.sha256"))
        if len(_rows(table)) != row_count:
            raise ValueError(f"publication table row count disagrees with manifest: {name}")
    results = _object(native.get("results"), "results")
    output_hashes = _hashes(results.get("output_sha256"), "results.output_sha256")
    family = _artifact(
        path.parent, "native_ora_families.tsv", output_hashes["native_ora_families.tsv"]
    )
    input_hashes = _hashes(results.get("input_sha256"), "results.input_sha256")
    _artifact(path.parent, "scenarios.tsv", input_hashes["scenarios.tsv"])
    detail = _detail("publication", "review", review)
    detail += [
        ("source", "publication", "path", str(path)),
        ("source", "publication", "sha256", file_sha256(path)),
    ]
    for row in _rows(path.parent / "interpretation_handoff.tsv"):
        finding = row.get("finding_id")
        if not finding:
            raise ValueError("interpretation handoff lacks finding_id")
        detail.extend(_detail("interpretation_flag", finding, row))
    return detail, family


def _native(
    path: Path,
    content: Mapping[str, object],
    issue: str,
    expected: Mapping[str, str] | None,
    review_path: Path | None,
) -> InfluenceAssessment:
    protocol = _object(content.get("protocol"), "protocol")
    if protocol.get("version") != NATIVE_VERSION:
        raise ValueError("incompatible native influence protocol version")
    for method in ("composition", "ora"):
        if not _object(protocol.get(method), f"protocol.{method}"):
            raise ValueError(f"protocol.{method} must not be empty")
    completed = content.get("completed")
    if not isinstance(completed, bool):
        raise ValueError("completed must be a boolean")
    hashes = _inputs(content, expected)
    for key in ("protocol_sha256", "contribution_manifest_sha256"):
        _hash(content.get(key), key)
    _hashes(content.get("code_sha256"), "code_sha256")
    commit = _text(content.get("code_commit"), "code_commit")
    if re.fullmatch(r"[0-9a-f]{40}", commit) is None:
        raise ValueError("code_commit must be a full Git revision")
    if not isinstance(content.get("code_dirty"), bool):
        raise ValueError("code_dirty must be a boolean")
    detail = _detail(
        "execution",
        "native_influence",
        {key: value for key, value in content.items() if key not in {"protocol", "results"}},
    )
    detail += _detail("input", "native_influence", hashes)
    detail += _detail("protocol", "native_influence", protocol)
    detail += [
        ("source", "native_manifest", "path", str(path)),
        ("source", "native_manifest", "sha256", file_sha256(path)),
    ]
    if not completed:
        return InfluenceAssessment(
            {}, f"native influence execution incomplete {issue}", tuple(detail)
        )
    results = _object(content.get("results"), "results")
    for key in (
        "scenarios",
        "scopes",
        "composition_rows",
        "composition_comparisons",
        "ora_comparisons",
        "ora_families",
    ):
        _integer(results.get(key), f"results.{key}", 1)
    detail += _detail("execution_results", "native_influence", results)
    if review_path is not None:
        review_detail, family_path = _review(review_path, path, content)
        detail += review_detail
        detail.append(("verification", "native_influence", "artifacts", "compact review tables"))
    else:
        for group in ("input_sha256", "output_sha256"):
            declared = _hashes(results.get(group), f"results.{group}")
            for name, expected_hash in declared.items():
                _artifact(path.parent, name, expected_hash)
        output_hashes = _hashes(results.get("output_sha256"), "results.output_sha256")
        family_path = _artifact(
            path.parent, "native_ora_families.tsv", output_hashes["native_ora_families.tsv"]
        )
        detail.append(
            ("verification", "native_influence", "artifacts", "declared native artifacts")
        )
        detail.append(("publication", "review", "status", "not supplied"))
    detail += _scope_detail(family_path, results)
    return InfluenceAssessment(
        {},
        f"native influence results reported; no procedure-level verdict supplied {issue}",
        tuple(detail),
    )


def read_influence_assessment(
    path: Path | None,
    issue: str = "#190",
    *,
    expected_input_sha256: Mapping[str, str] | None = None,
) -> InfluenceAssessment:
    """Read native/review evidence or an explicit, versioned procedure assessment.

    Missing input says nothing about whether execution occurred elsewhere. A
    malformed or incompatible source never becomes a positive assessment.
    """
    if path is None or not path.exists():
        return InfluenceAssessment({}, f"pending {issue}: influence evidence not supplied", ())
    source: tuple[Detail, ...] = (
        ("source", "supplied_influence", "path", str(path)),
        ("source", "influence_adapter", "sha256", file_sha256(Path(__file__))),
    )
    try:
        content = _read(path)
        source += (("source", "supplied_influence", "sha256", file_sha256(path)),)
        if "schema_version" in content:
            if content["schema_version"] != ASSESSMENT_VERSION:
                raise ValueError("incompatible influence assessment schema_version")
            result = _explicit(content, issue, expected_input_sha256)
        elif "native_manifest_sha256" in content and "tables" in content:
            native_path = path.parent / "run_manifest.json"
            result = _native(native_path, _read(native_path), issue, expected_input_sha256, path)
        else:
            candidate = path.parent / "review" / "publication_manifest.json"
            review = candidate if candidate.exists() else None
            result = _native(path, content, issue, expected_input_sha256, review)
        return InfluenceAssessment(result.verdicts, result.default, source + result.detail)
    except (OSError, UnicodeError, ValueError, KeyError, csv.Error, EOFError) as error:
        return InfluenceAssessment({}, f"unreadable {issue}: {error}", source)
