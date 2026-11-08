"""Noel baseline identity, typed missingness and consumed-byte checks (#1029)."""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
from pathlib import Path
from typing import Any

import pytest
from scipy.stats import false_discovery_control

from mtp import influence_extension as extension


def families() -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    baseline, native = [], []
    for scope in ("pooled", "feces"):
        pvalues = [0.001, 0.01, 0.02, None, 0.1, 0.8]
        filled = [1.0 if p is None else p for p in pvalues]
        bh = false_discovery_control(filled, method="bh")
        by = false_discovery_control(filled, method="by")
        for index, p in enumerate(pvalues):
            row = {
                "scope": scope,
                "trait_key": f"METPO:fixture{index}",
                "family_size": 6,
                "p": "" if p is None else p,
                "q_bh_six": "" if p is None else float(bh[index]),
                "q_by_six": "" if p is None else float(by[index]),
            }
            baseline.append(row)
            native.append(
                {
                    "scenario": "baseline",
                    "scope": scope,
                    "representation": "shared_fraction",
                    "trait_key": row["trait_key"],
                    "denominator_trait": "",
                    "unit": "fraction_difference",
                    "correction": "BH6",
                    "before_p": row["p"],
                    "before_q": row["q_bh_six"],
                }
            )
    return baseline, native


def write_native(path: Path, rows: list[dict[str, Any]]) -> str:
    text = io.StringIO(newline="")
    writer = csv.DictWriter(text, fieldnames=list(rows[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
    raw = gzip.compress(text.getvalue().encode(), mtime=0)
    path.write_bytes(raw)
    return hashlib.sha256(raw).hexdigest()


def check(
    tmp_path: Path, baseline: list[dict[str, Any]], native: list[dict[str, Any]]
) -> dict[str, Any]:
    path = tmp_path / "effects.tsv.gz"
    digest = write_native(path, native)
    result: dict[str, Any] = extension.check_noel_baseline(path, digest, baseline)
    return result


def test_checks_every_scope_family_and_labels_by_as_derived(tmp_path: Path) -> None:
    baseline, native = families()
    # Another representation and another scenario are not members of this family.
    native += [
        dict(native[0], representation="native_noel_zero_fraction", before_p="bad"),
        dict(native[0], scenario="removal", before_q="bad"),
    ]
    result = check(tmp_path, baseline, native)
    assert result["status"] == "matched"
    assert result["rows"] == 12
    assert result["families"] == 2
    assert result["by_source"] == "derived_from_noel_p_values"
    assert result["relative_tolerance"] == 1e-10
    assert result["absolute_tolerance"] == 0.0
    assert len(result["comparisons"]) == 12
    assert all(row["status"] == "matched" for row in result["comparisons"])
    missing = [r for r in result["comparisons"] if r["trait_key"] == "METPO:fixture3"]
    assert all(r["noel_p"] is None and r["noel_q_by_six_derived"] is None for r in missing)


@pytest.mark.parametrize("field", ["p", "q_bh_six", "q_by_six"])
def test_refuses_sidecar_mismatch_in_each_compared_quantity(tmp_path: Path, field: str) -> None:
    baseline, native = families()
    baseline[-2][field] = 0.51  # Feces, not merely the pooled family.
    with pytest.raises(ValueError, match=field):
        check(tmp_path, baseline, native)


@pytest.mark.parametrize("field", ["before_p", "before_q"])
def test_refuses_changed_native_values(tmp_path: Path, field: str) -> None:
    baseline, native = families()
    native[0][field] = 0.5
    with pytest.raises(ValueError, match="mismatch"):
        check(tmp_path, baseline, native)


@pytest.mark.parametrize(
    "damage", ["duplicate", "missing", "extra_scope", "wrong_family", "wrong_unit", "denominator"]
)
def test_refuses_incomplete_or_ambiguous_native_family(tmp_path: Path, damage: str) -> None:
    baseline, native = families()
    if damage == "duplicate":
        native.append(dict(native[0]))
    elif damage == "missing":
        native.pop(3)  # Even an untestable slot remains a family member.
    elif damage == "extra_scope":
        native.append(dict(native[0], scope="unexpected"))
    elif damage == "wrong_family":
        native[0]["correction"] = "BH5"
    elif damage == "wrong_unit":
        native[0]["unit"] = "log2_ratio_contrast"
    else:
        native[0]["denominator_trait"] = "METPO:other"
    with pytest.raises(ValueError, match=r"Noel"):
        check(tmp_path, baseline, native)


@pytest.mark.parametrize("value", ["nan", "inf", "-0.1", "1.1", "bad", ""])
def test_refuses_bad_p_values_and_inconsistent_missingness(tmp_path: Path, value: str) -> None:
    baseline, native = families()
    native[0]["before_p"] = value
    with pytest.raises(ValueError, match=r"Noel"):
        check(tmp_path, baseline, native)


@pytest.mark.parametrize(("native_p", "sidecar_p"), [(0.0, 1e-15), (1e-100, 1e-20)])
def test_absolute_tolerance_does_not_hide_tiny_p_mismatches(
    tmp_path: Path, native_p: float, sidecar_p: float
) -> None:
    baseline, native = families()
    native[0]["before_p"] = native_p
    # Keep both q columns consistent with the changed native family. Only p is
    # wrong in the sidecar, so a broad absolute tolerance cannot hide behind an
    # unrelated BY mismatch after accepting that p difference.
    pooled = native[:6]
    filled = [1.0 if row["before_p"] == "" else row["before_p"] for row in pooled]
    bh = false_discovery_control(filled, method="bh")
    by = false_discovery_control(filled, method="by")
    for index, row in enumerate(pooled):
        if row["before_p"] != "":
            row["before_q"] = float(bh[index])
            baseline[index]["q_bh_six"] = float(bh[index])
            baseline[index]["q_by_six"] = float(by[index])
    baseline[0]["p"] = sidecar_p
    with pytest.raises(ValueError, match="pooled/METPO:fixture0/p"):
        check(tmp_path, baseline, native)


def test_refuses_by_values_computed_after_dropping_untestable_slot(tmp_path: Path) -> None:
    baseline, native = families()
    for scope in ("pooled", "feces"):
        tested = [row for row in baseline if row["scope"] == scope and row["p"] != ""]
        assert len(tested) == 5
        five_slot_by = false_discovery_control([row["p"] for row in tested], method="by")
        for row, qvalue in zip(tested, five_slot_by, strict=True):
            row["q_by_six"] = float(qvalue)
    assert len(baseline) == len(native) == 12
    with pytest.raises(ValueError, match="q_by_six"):
        check(tmp_path, baseline, native)


def test_tolerates_only_declared_relative_roundoff(tmp_path: Path) -> None:
    baseline, native = families()
    baseline[0]["p"] *= 1 + 1e-11
    assert check(tmp_path, baseline, native)["status"] == "matched"
    baseline[0]["p"] *= 1 + 1e-8
    with pytest.raises(ValueError, match="mismatch"):
        check(tmp_path, baseline, native)


def test_consumed_gzip_bytes_match_the_verified_digest(tmp_path: Path) -> None:
    baseline, native = families()
    path = tmp_path / "effects.tsv.gz"
    digest = write_native(path, native)
    native[0]["before_p"] = 0.002
    write_native(path, native)
    with pytest.raises(ValueError, match=r"hash|digest"):
        extension.check_noel_baseline(path, digest, baseline)


def test_consumed_gzip_hash_cannot_be_replaced_by_a_fresh_disk_hash(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    baseline, native = families()
    path = tmp_path / "effects.tsv.gz"
    digest = write_native(path, native)
    original_bytes = path.read_bytes()
    different = bytearray(original_bytes)
    different[4] ^= 1  # A gzip mtime byte: parsed TSV rows and total size stay identical.
    parser_bytes = bytes(different)
    assert len(parser_bytes) == len(original_bytes)
    assert gzip.decompress(parser_bytes) == gzip.decompress(original_bytes)
    assert hashlib.sha256(parser_bytes).hexdigest() != digest
    assert extension.check_noel_baseline(path, digest, baseline)["status"] == "matched"
    original_open = Path.open
    binary_reads: list[Path] = []

    def parser_then_disk(current: Path, mode: str = "r", *args: Any, **kwargs: Any) -> Any:
        if current == path and mode == "rb":
            binary_reads.append(current)
            if len(binary_reads) == 1:
                return io.BytesIO(parser_bytes)
        return original_open(current, mode, *args, **kwargs)

    monkeypatch.setattr(Path, "open", parser_then_disk)
    with pytest.raises(ValueError, match=r"hash|digest"):
        extension.check_noel_baseline(path, digest, baseline)
    assert binary_reads == [path]
    with original_open(path, "rb") as persisted:
        assert persisted.read() == original_bytes


def test_native_effect_rows_are_streamed_without_read_bytes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    baseline, native = families()
    path = tmp_path / "effects.tsv.gz"
    digest = write_native(path, native)

    def refuse_read_bytes(_path: Path) -> bytes:
        raise AssertionError("Do not materialize the whole native effects file")

    monkeypatch.setattr(Path, "read_bytes", refuse_read_bytes)
    result = extension.check_noel_baseline(path, digest, baseline)
    assert result["rows"] == 12
