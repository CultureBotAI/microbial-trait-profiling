"""Compare core-bound strata with the preserved direct-path outputs (#727, #1511).

Inputs: historical and regenerated stratification directories. Output: JSON to
stdout with exact file comparisons and complete/disjoint partition accounting.
Only output-directory references and validated additive statistic/scope/basis
labels in the three cross-stratum summary reports may differ.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
from pathlib import Path

from mtp.paths import file_sha256
from mtp.stratification import BEST_STAT_NOTE, PRIMARY_STRATA

ALLOWED_BASIS_VALUES = {"FDR q", "nominal p", "unknown", "unavailable"}
ALLOWED_BEST_BASIS_VALUES = {"FDR q", "nominal p", "unavailable"}
ALLOWED_STAT_TYPES = {"fdr_q", "p_value", "unavailable"}
ALLOWED_SCOPE_VALUES = {"global", "profile", "none", "unknown", "unavailable"}
ADDITIVE_SUMMARY_SUFFIXES: dict[str, set[str]] = {
    "best_stat_type": ALLOWED_STAT_TYPES,
    "best_stat_scope": ALLOWED_SCOPE_VALUES,
    "best_stat_basis": ALLOWED_BEST_BASIS_VALUES,
    "significance_basis": ALLOWED_BASIS_VALUES,
    "fdr_scope": ALLOWED_SCOPE_VALUES,
}
ADDITIVE_DECILE_COLUMNS: dict[str, set[str]] = {
    "threshold_statistic": ALLOWED_BASIS_VALUES,
    "hit_statistic": ALLOWED_BASIS_VALUES,
    "fdr_scope": ALLOWED_SCOPE_VALUES,
}
LEGACY_HIT_DESCRIPTION = (
    "Sig Traits: any test passing the analyzer's declared p/q and fold thresholds."
)
VALID_HIT_DESCRIPTIONS = {
    LEGACY_HIT_DESCRIPTION,
    "FDR hits: traits with any test passing the declared FDR q and fold thresholds.",
    "Nominal hits: traits with any test passing the declared nominal p and fold thresholds.",
    "Hits (statistic not declared): traits with any test passing the analyzer's thresholds.",
}
LEGACY_DECILE_HEADER = (
    f"{'Decile':<8} {'Min K':<8} {'Max K':<8} {'Total Traits':<14} {'Sig Traits':<12} {'% Sig':<8}"
)
VALID_DECILE_HEADERS = {
    LEGACY_DECILE_HEADER,
    f"{'Decile':<8} {'Min K':<8} {'Max K':<8} {'Total Traits':<14} {'FDR hits':<24} {'% Hits':<8}",
    f"{'Decile':<8} {'Min K':<8} {'Max K':<8} {'Total Traits':<14} "
    f"{'Nominal hits':<24} {'% Hits':<8}",
    f"{'Decile':<8} {'Min K':<8} {'Max K':<8} {'Total Traits':<14} "
    f"{'Hits (statistic not declared)':<29} {'% Hits':<8}",
}
SUPPORT_STATE_NOTE = (
    "Support-state names retain 'significant'; interpret each stratum using "
    "its hit statistic and FDR scope above."
)


def _split_comments_and_csv(path: Path) -> tuple[dict[str, str], list[str], list[list[str]]] | None:
    comments: dict[str, str] = {}
    body_lines: list[str] = []
    in_comments = True
    for line in path.read_text(encoding="utf-8").splitlines():
        if in_comments and line.startswith("#"):
            key, sep, value = line[1:].strip().partition(": ")
            if not sep or not key or key in comments:
                return None
            comments[key] = value
        else:
            in_comments = False
            body_lines.append(line)
    reader = csv.reader(io.StringIO("\n".join(body_lines)), delimiter="\t", strict=True)
    header = next(reader, [])
    if not header or len(header) != len(set(header)) or any(not col for col in header):
        return None
    rows: list[list[str]] = []
    for row in reader:
        if len(row) != len(header):
            return None
        rows.append(row)
    return comments, header, rows


def _valid_probability(value: str) -> bool:
    try:
        parsed = float(value)
    except ValueError:
        return False
    return math.isfinite(parsed) and 0.0 < parsed <= 1.0


def _summary_tsv_scientifically_equal(old: Path, new: Path) -> bool:
    old_parsed = _split_comments_and_csv(old)
    new_parsed = _split_comments_and_csv(new)
    if old_parsed is None or new_parsed is None:
        return False
    old_comments, old_header, old_rows = old_parsed
    new_comments, new_header, new_rows = new_parsed
    if "Evaluated-Strata" not in old_comments or old_comments.get(
        "Evaluated-Strata"
    ) != new_comments.get("Evaluated-Strata"):
        return False
    strata = [s.strip() for s in old_comments["Evaluated-Strata"].split(",") if s.strip()]
    for key, old_val in old_comments.items():
        if new_comments.get(key) != old_val:
            return False
    for key, new_val in new_comments.items():
        if key in old_comments:
            continue
        if key == "Best-Statistic":
            if new_val != BEST_STAT_NOTE:
                return False
            continue
        matched = False
        for stratum in strata:
            prefix = f"{stratum}-"
            if key.startswith(prefix):
                field = key[len(prefix) :]
                if field in {"Threshold-Statistic", "Hit-Statistic"}:
                    matched = new_val in ALLOWED_BASIS_VALUES
                elif field == "FDR-Scope":
                    matched = new_val in ALLOWED_SCOPE_VALUES
                elif field == "Stat-Threshold":
                    matched = _valid_probability(new_val)
                break
        if not matched:
            return False

    old_set = set(old_header)
    if [col for col in new_header if col in old_set] != old_header:
        return False
    allowed_added: dict[str, set[str]] = {
        f"{stratum}_{suffix}": allowed
        for stratum in strata
        for suffix, allowed in ADDITIVE_SUMMARY_SUFFIXES.items()
    }
    added_cols = [col for col in new_header if col not in old_set]
    if any(col not in allowed_added for col in added_cols):
        return False
    new_idx = {col: idx for idx, col in enumerate(new_header)}
    for row in new_rows:
        for col in added_cols:
            if row[new_idx[col]] not in allowed_added[col]:
                return False
        for stratum in strata:
            type_col = f"{stratum}_best_stat_type"
            basis_col = f"{stratum}_best_stat_basis"
            scope_col = f"{stratum}_best_stat_scope"
            if type_col in new_idx:
                stat_type = row[new_idx[type_col]]
                if basis_col in new_idx:
                    expected_basis = {
                        "fdr_q": "FDR q",
                        "p_value": "nominal p",
                        "unavailable": "unavailable",
                    }.get(stat_type)
                    if row[new_idx[basis_col]] != expected_basis:
                        return False
                if scope_col in new_idx:
                    stat_scope = row[new_idx[scope_col]]
                    if stat_type == "fdr_q" and stat_scope not in {"global", "profile", "unknown"}:
                        return False
                    if stat_type == "p_value" and stat_scope != "none":
                        return False
                    if stat_type == "unavailable" and stat_scope != "unavailable":
                        return False
    projected_new = [[row[new_idx[col]] for col in old_header] for row in new_rows]
    return projected_new == old_rows


def _deciles_tsv_scientifically_equal(old: Path, new: Path) -> bool:
    old_parsed = _split_comments_and_csv(old)
    new_parsed = _split_comments_and_csv(new)
    if old_parsed is None or new_parsed is None:
        return False
    old_comments, old_header, old_rows = old_parsed
    new_comments, new_header, new_rows = new_parsed
    if old_comments or new_comments:
        return False
    old_set = set(old_header)
    if [col for col in new_header if col in old_set] != old_header:
        return False
    added_cols = [col for col in new_header if col not in old_set]
    if any(col not in ADDITIVE_DECILE_COLUMNS for col in added_cols):
        return False
    new_idx = {col: idx for idx, col in enumerate(new_header)}
    for row in new_rows:
        for col in added_cols:
            if row[new_idx[col]] not in ADDITIVE_DECILE_COLUMNS[col]:
                return False
        if (
            "threshold_statistic" in new_idx
            and "hit_statistic" in new_idx
            and row[new_idx["threshold_statistic"]] != row[new_idx["hit_statistic"]]
        ):
            return False
    projected_new = [[row[new_idx[col]] for col in old_header] for row in new_rows]
    return projected_new == old_rows


def _normalize_summary_txt(text: str) -> list[str] | None:
    cleaned = text.replace(f"\n{BEST_STAT_NOTE}\n", "")
    normalized: list[str] = []
    in_deciles = False
    for line in cleaned.splitlines():
        if line.startswith("  Hit statistic: "):
            payload = line[len("  Hit statistic: ") :]
            hit_stat, sep, scope_part = payload.partition(" | FDR scope: ")
            if (
                not sep
                or hit_stat not in ALLOWED_BASIS_VALUES
                or scope_part not in ALLOWED_SCOPE_VALUES
            ):
                return None
            continue
        if line in VALID_HIT_DESCRIPTIONS:
            normalized.append(LEGACY_HIT_DESCRIPTION)
            continue
        if line in VALID_DECILE_HEADERS:
            normalized.append(LEGACY_DECILE_HEADER)
            in_deciles = True
            continue
        if in_deciles:
            if not line.strip():
                in_deciles = False
                normalized.append(line)
                continue
            parts = line.split()
            if len(parts) == 7 and parts[6] == "%":
                pct_str = parts[5]
            elif len(parts) == 6 and parts[5].endswith("%"):
                pct_str = parts[5][:-1]
            else:
                return None
            try:
                normalized.append(
                    f"{int(parts[0]):<8} {int(parts[1]):<8} {int(parts[2]):<8} "
                    f"{int(parts[3]):<14} {int(parts[4]):<12} {float(pct_str):<8.1f}%"
                )
            except ValueError:
                return None
            continue
        if line == SUPPORT_STATE_NOTE:
            continue
        normalized.append(line)
    return normalized


def _summary_txt_scientifically_equal(old: Path, new: Path, before: Path, after: Path) -> bool:
    old_text = old.read_text(encoding="utf-8").replace(str(before), "<OUTPUT>")
    new_text = new.read_text(encoding="utf-8").replace(str(after), "<OUTPUT>")
    if old_text == new_text:
        return True
    old_norm = _normalize_summary_txt(old_text)
    new_norm = _normalize_summary_txt(new_text)
    return old_norm is not None and old_norm == new_norm


def compare(before: Path, after: Path) -> dict[str, object]:
    """Fail on changed scientific products or incomplete/overlapping partitions."""
    reserved = {"core_consumption.json", "core_consumption.coverage.tsv"}
    old_files = {path.name for path in before.iterdir() if path.is_file()}
    new_files = {path.name for path in after.iterdir() if path.is_file()} - reserved
    if len(old_files) != 24 or old_files != new_files:
        raise ValueError("Expected the same 24 stratification products in both directories")
    comparisons = []
    for name in sorted(old_files):
        old, new = before / name, after / name
        old_hash, new_hash = file_sha256(old), file_sha256(new)
        equal = old_hash == new_hash
        if not equal:
            if name == "stratification_sensitivity_summary.txt":
                equal = _summary_txt_scientifically_equal(old, new, before, after)
            elif name == "stratification_sensitivity_summary.tsv":
                equal = _summary_tsv_scientifically_equal(old, new)
            elif name == "annotation_density_deciles.tsv":
                equal = _deciles_tsv_scientifically_equal(old, new)
        if not equal:
            raise ValueError(f"Scientific output differs: {name}")
        comparisons.append(
            {
                "file": name,
                "historical_sha256": old_hash,
                "current_sha256": new_hash,
                "byte_identical": old_hash == new_hash,
                "scientifically_identical": equal,
            }
        )
    memberships: dict[str, set[tuple[str, str, str]]] = {}
    for stratum in PRIMARY_STRATA:
        with (after / f"study_trait_profiles_{stratum}.tsv").open() as handle:
            rows = csv.DictReader(handle, delimiter="\t")
            keys = [(row["study_id"], row["experiment_id"], row["inc_dec"]) for row in rows]
        if len(keys) != len(set(keys)):
            raise ValueError(f"Duplicate profile keys in {stratum}")
        memberships[stratum] = set(keys)
    feces, other = memberships["feces"], memberships["feces_excluded"]
    if feces & other or feces | other != memberships["pooled"]:
        raise ValueError("Feces partitions must be disjoint and cover the pooled profiles")
    return {
        "issue": 727,
        "historical": str(before),
        "current": str(after),
        "passed": True,
        "partition_profiles": {key: len(value) for key, value in memberships.items()},
        "complete_disjoint_partitions": True,
        "products": comparisons,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("historical", type=Path)
    parser.add_argument("current", type=Path)
    args = parser.parse_args()
    print(json.dumps(compare(args.historical.resolve(), args.current.resolve()), indent=2))


if __name__ == "__main__":
    main()
