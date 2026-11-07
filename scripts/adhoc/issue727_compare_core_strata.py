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
from typing import NamedTuple

from mtp.paths import file_sha256
from mtp.stratification import (
    BEST_STAT_NOTE,
    PRIMARY_STRATA,
)

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
EXPECTED_HIT_DESCRIPTIONS = {
    "FDR q": "FDR hits: traits with any test passing the declared FDR q and fold thresholds.",
    "nominal p": (
        "Nominal hits: traits with any test passing the declared nominal p and fold thresholds."
    ),
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
EXPECTED_DECILE_HEADERS = {
    "FDR q": (
        f"{'Decile':<8} {'Min K':<8} {'Max K':<8} {'Total Traits':<14} "
        f"{'FDR hits':<24} {'% Hits':<8}"
    ),
    "nominal p": (
        f"{'Decile':<8} {'Min K':<8} {'Max K':<8} {'Total Traits':<14} "
        f"{'Nominal hits':<24} {'% Hits':<8}"
    ),
}
SUPPORT_STATE_NOTE = (
    "Support-state names retain 'significant'; interpret each stratum using "
    "its hit statistic and FDR scope above."
)


class _StratumTruth(NamedTuple):
    stat_type: str | None
    fdr_scope: str
    stat_threshold: str | None
    hit_statistic: str
    trait_types: dict[str, str]
    trait_sig_basis: dict[str, str]


def _parse_optional_finite(raw: str) -> float | None:
    stripped = raw.strip()
    if not stripped:
        return None
    try:
        parsed = float(stripped)
    except ValueError:
        return None
    return parsed if math.isfinite(parsed) else None


def _independent_row_stat_type(
    q_val: float | None, p_val: float | None, stat_type: str | None
) -> str:
    if q_val is not None:
        if stat_type == "p_value":
            if p_val is None or q_val != p_val:
                raise ValueError("nominal Stat-Type contradicts best_qvalue: expected copied p")
            return "p_value"
        if stat_type is None and p_val is not None and q_val == p_val:
            raise ValueError(
                "Missing Stat-Type header with best_qvalue == best_pvalue: "
                "cannot distinguish fdr_q from copied nominal p"
            )
        return "fdr_q"
    if p_val is not None:
        return "p_value"
    return "unavailable"


def _expected_labels_from_tested(directory: Path) -> dict[str, _StratumTruth]:
    truth: dict[str, _StratumTruth] = {}
    for stratum in PRIMARY_STRATA:
        tested_path = directory / f"trait_enrichment_traits_tested_{stratum}.tsv"
        if not tested_path.exists():
            continue
        parsed = _split_comments_and_csv(tested_path)
        if parsed is None:
            raise ValueError(f"Malformed tested traits table: {tested_path}")
        comments, header, rows = parsed
        raw_stat_type = comments.get("Stat-Type")
        if raw_stat_type is not None and raw_stat_type not in {"fdr_q", "p_value"}:
            raise ValueError(f"Invalid Stat-Type in {tested_path}: {raw_stat_type!r}")
        raw_fdr = comments.get("FDR")
        if raw_fdr is not None and raw_fdr not in {"bh", "none"}:
            raise ValueError(f"Invalid FDR in {tested_path}: {raw_fdr!r}")
        raw_scope = comments.get("FDR-Scope")
        if raw_scope is not None and raw_scope not in {"global", "profile", "none"}:
            raise ValueError(f"Invalid FDR-Scope in {tested_path}: {raw_scope!r}")
        stat_threshold = comments.get("Stat-Threshold") or comments.get("P-Threshold")
        if stat_threshold is not None and not _valid_probability(stat_threshold):
            raise ValueError(
                f"Stat-Threshold out of range (0, 1] in {tested_path}: {stat_threshold!r}"
            )
        if raw_stat_type == "p_value" or raw_fdr == "none":
            if raw_stat_type == "fdr_q" or raw_fdr == "bh":
                raise ValueError(
                    f"Contradictory Stat-Type {raw_stat_type!r} and FDR {raw_fdr!r} "
                    f"in {tested_path}"
                )
            stat_type: str | None = "p_value"
            fdr_scope = "none"
        elif raw_stat_type == "fdr_q" or raw_fdr == "bh":
            if raw_scope == "none":
                raise ValueError(
                    f"Contradictory Stat-Type 'fdr_q' and FDR-Scope 'none' in {tested_path}"
                )
            stat_type = "fdr_q"
            fdr_scope = raw_scope if raw_scope in {"global", "profile"} else "unknown"
        else:
            stat_type = None
            fdr_scope = raw_scope if raw_scope in {"global", "profile", "none"} else "unknown"
        hit_statistic = {"p_value": "nominal p", "fdr_q": "FDR q"}.get(
            stat_type or "unknown", "unknown"
        )
        col_idx = {col: idx for idx, col in enumerate(header)}
        trait_types: dict[str, str] = {}
        trait_sig_basis: dict[str, str] = {}
        for row in rows:
            trait_name = row[col_idx["trait"]]
            q_val = (
                _parse_optional_finite(row[col_idx["best_qvalue"]])
                if "best_qvalue" in col_idx
                else None
            )
            p_val = (
                _parse_optional_finite(row[col_idx["best_pvalue"]])
                if "best_pvalue" in col_idx
                else None
            )
            trait_types[trait_name] = _independent_row_stat_type(q_val, p_val, stat_type)
            trait_sig_basis[trait_name] = (
                hit_statistic if int(row[col_idx["n_profiles_tested"]]) > 0 else "unavailable"
            )
        truth[stratum] = _StratumTruth(
            stat_type=stat_type,
            fdr_scope=fdr_scope,
            stat_threshold=stat_threshold,
            hit_statistic=hit_statistic,
            trait_types=trait_types,
            trait_sig_basis=trait_sig_basis,
        )
    return truth


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


def _summary_tsv_scientifically_equal(
    old: Path, new: Path, truth: dict[str, _StratumTruth] | None = None
) -> bool:
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
                stratum_truth = truth.get(stratum) if truth is not None else None
                if field in {"Threshold-Statistic", "Hit-Statistic"}:
                    matched = new_val in ALLOWED_BASIS_VALUES and (
                        stratum_truth is None
                        or stratum_truth.hit_statistic == "unknown"
                        or new_val == stratum_truth.hit_statistic
                    )
                elif field == "FDR-Scope":
                    matched = new_val in ALLOWED_SCOPE_VALUES and (
                        stratum_truth is None
                        or stratum_truth.fdr_scope == "unknown"
                        or new_val == stratum_truth.fdr_scope
                    )
                elif field == "Stat-Threshold":
                    matched = _valid_probability(new_val) and (
                        stratum_truth is None
                        or stratum_truth.stat_threshold is None
                        or float(new_val) == float(stratum_truth.stat_threshold)
                    )
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
            stratum_truth = truth.get(stratum) if truth is not None else None
            if stratum_truth is not None and "trait" in new_idx:
                trait_name = row[new_idx["trait"]]
                if (
                    type_col in new_idx
                    and stratum_truth.stat_type is not None
                    and row[new_idx[type_col]]
                    != stratum_truth.trait_types.get(trait_name, "unavailable")
                ):
                    return False
                if (
                    scope_col in new_idx
                    and type_col in new_idx
                    and row[new_idx[type_col]] == "fdr_q"
                    and stratum_truth.fdr_scope != "unknown"
                    and row[new_idx[scope_col]] != stratum_truth.fdr_scope
                ):
                    return False
                sig_col = f"{stratum}_significance_basis"
                if (
                    sig_col in new_idx
                    and stratum_truth.hit_statistic != "unknown"
                    and row[new_idx[sig_col]]
                    != stratum_truth.trait_sig_basis.get(trait_name, "unavailable")
                ):
                    return False
                fdr_col = f"{stratum}_fdr_scope"
                if (
                    fdr_col in new_idx
                    and stratum_truth.fdr_scope != "unknown"
                    and row[new_idx[fdr_col]] != stratum_truth.fdr_scope
                ):
                    return False
    projected_new = [[row[new_idx[col]] for col in old_header] for row in new_rows]
    return projected_new == old_rows


def _deciles_tsv_scientifically_equal(
    old: Path, new: Path, truth: dict[str, _StratumTruth] | None = None
) -> bool:
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
    pooled_truth = truth.get("pooled") if truth is not None else None
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
        if pooled_truth is not None:
            if (
                "hit_statistic" in new_idx
                and pooled_truth.hit_statistic != "unknown"
                and row[new_idx["hit_statistic"]] != pooled_truth.hit_statistic
            ):
                return False
            if (
                "fdr_scope" in new_idx
                and pooled_truth.fdr_scope != "unknown"
                and row[new_idx["fdr_scope"]] != pooled_truth.fdr_scope
            ):
                return False
    projected_new = [[row[new_idx[col]] for col in old_header] for row in new_rows]
    return projected_new == old_rows


def _normalize_summary_txt(
    text: str, truth: dict[str, _StratumTruth] | None = None
) -> list[str] | None:
    cleaned = text.replace(f"\n{BEST_STAT_NOTE}\n", "")
    normalized: list[str] = []
    in_deciles = False
    current_stratum: str | None = None
    pooled_truth = truth.get("pooled") if truth is not None else None
    for line in cleaned.splitlines():
        if line.startswith("Stratum:"):
            current_stratum = line[len("Stratum:") :].split("|", 1)[0].strip()
        if line.startswith("  Hit statistic: "):
            payload = line[len("  Hit statistic: ") :]
            hit_stat, sep, scope_part = payload.partition(" | FDR scope: ")
            if (
                not sep
                or hit_stat not in ALLOWED_BASIS_VALUES
                or scope_part not in ALLOWED_SCOPE_VALUES
            ):
                return None
            if truth is not None and current_stratum is not None and current_stratum in truth:
                stratum_truth = truth[current_stratum]
                if (
                    stratum_truth.hit_statistic != "unknown"
                    and hit_stat != stratum_truth.hit_statistic
                ):
                    return None
                if stratum_truth.fdr_scope != "unknown" and scope_part != stratum_truth.fdr_scope:
                    return None
            continue
        if line in VALID_HIT_DESCRIPTIONS:
            if (
                pooled_truth is not None
                and pooled_truth.hit_statistic in EXPECTED_HIT_DESCRIPTIONS
                and line != EXPECTED_HIT_DESCRIPTIONS[pooled_truth.hit_statistic]
            ):
                return None
            normalized.append(LEGACY_HIT_DESCRIPTION)
            continue
        if line in VALID_DECILE_HEADERS:
            if (
                pooled_truth is not None
                and pooled_truth.hit_statistic in EXPECTED_DECILE_HEADERS
                and line != EXPECTED_DECILE_HEADERS[pooled_truth.hit_statistic]
            ):
                return None
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


def _summary_txt_scientifically_equal(
    old: Path,
    new: Path,
    before: Path,
    after: Path,
    truth: dict[str, _StratumTruth] | None = None,
) -> bool:
    old_text = old.read_text(encoding="utf-8").replace(str(before), "<OUTPUT>")
    new_text = new.read_text(encoding="utf-8").replace(str(after), "<OUTPUT>")
    if old_text == new_text:
        return True
    old_norm = _normalize_summary_txt(old_text)
    new_norm = _normalize_summary_txt(new_text, truth=truth)
    return old_norm is not None and new_norm is not None and old_norm == new_norm


def compare(before: Path, after: Path) -> dict[str, object]:
    """Fail on changed scientific products or incomplete/overlapping partitions."""
    reserved = {"core_consumption.json", "core_consumption.coverage.tsv"}
    old_files = {path.name for path in before.iterdir() if path.is_file()}
    new_files = {path.name for path in after.iterdir() if path.is_file()} - reserved
    if len(old_files) != 24 or old_files != new_files:
        raise ValueError("Expected the same 24 stratification products in both directories")
    truth = _expected_labels_from_tested(after)
    comparisons = []
    for name in sorted(old_files):
        old, new = before / name, after / name
        old_hash, new_hash = file_sha256(old), file_sha256(new)
        equal = old_hash == new_hash
        if not equal:
            if name == "stratification_sensitivity_summary.txt":
                equal = _summary_txt_scientifically_equal(old, new, before, after, truth)
            elif name == "stratification_sensitivity_summary.tsv":
                equal = _summary_tsv_scientifically_equal(old, new, truth)
            elif name == "annotation_density_deciles.tsv":
                equal = _deciles_tsv_scientifically_equal(old, new, truth)
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
