"""Reject incomplete or scientifically changed ORA refreshes (#720)."""

from __future__ import annotations

import csv
import io
import json
import subprocess
import sys
from pathlib import Path

import pytest

from mtp.paths import file_sha256


@pytest.fixture
def refresh_manifest(tmp_path: Path, request: pytest.FixtureRequest) -> Path:
    """Use actual producer output, including mixed metadata/table line endings."""
    settings: dict[str, str] = getattr(request, "param", {})
    zero = settings.get("zero")
    draws = 0 if zero in {"draws", "one_draw"} else 2
    dec_draws = 2 if zero == "one_draw" else draws
    annotated = 0 if zero else 2
    extra_trait = settings.get("extra_trait", "")
    extra_name = "mediadive.medium:1" if extra_trait == "medium" else "TraitB"
    extra_header = f"\t{extra_name}" if extra_trait else ""
    extra_inc = ("\t0" if extra_trait == "zero" else "\t1") if extra_trait else ""
    extra_dec = "\t0" if extra_trait else ""
    decreased = f"s1\te1\tdec\t{dec_draws}\t0{extra_dec}\n"
    if settings.get("unpaired") == "only":
        decreased = ""
    profiles = tmp_path / "profiles.tsv"
    matrix = tmp_path / "matrix.tsv"
    profiles.write_text(
        f"study_id\texperiment_id\tinc_dec\ttaxa_count\tTraitA{extra_header}\n"
        f"s1\te1\tinc\t{draws}\t{annotated}{extra_inc}\n"
        f"{decreased}",
        encoding="utf-8",
    )
    positive = 0 if zero == "traits" else 1
    matrix.write_text(
        f"taxon_id\tTraitA{extra_header}\n1\t{positive}{extra_inc}\n"
        f"2\t{positive}{extra_dec}\n3\t0{extra_dec}\n4\t0{extra_dec}\n",
        encoding="utf-8",
    )
    current = tmp_path / "current.tsv"
    subprocess.run(
        [
            sys.executable,
            "scripts/trait_enrichment.py",
            "--input",
            str(profiles),
            "--taxa-matrix",
            str(matrix),
            "--output",
            str(tmp_path / "report.txt"),
            "--all-results-output",
            str(current),
            "--p-threshold",
            settings.get("threshold", "1.0"),
            "--min-fold",
            settings.get("min_fold", "0.0"),
            "--max-fold",
            settings.get("max_fold", "0.67"),
            "--tail",
            settings.get("tail", "upper"),
            "--fdr-scope",
            settings.get("scope", "global"),
            "--no-fdr" if settings.get("fdr") == "False" else "--fdr",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    (tmp_path / "historical.tsv").write_bytes(
        b"".join(
            line
            for line in current.read_bytes().splitlines(keepends=True)
            if not line.startswith((b"# Profiles-SHA256:", b"# Taxa-Matrix-SHA256:"))
        )
    )
    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "pairs": [
                    {
                        "label": "fixture-upper",
                        "historical": "historical.tsv",
                        "current": "current.tsv",
                        "profiles": "profiles.tsv",
                        "taxa_matrix": "matrix.tsv",
                    }
                ]
            }
        ),
        encoding="utf-8",
    )
    return manifest


def run_audit(manifest: Path, output: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            sys.executable,
            "scripts/adhoc/issue720_audit_ora_refresh.py",
            "--manifest",
            str(manifest),
            "--output",
            str(output),
        ],
        capture_output=True,
        text=True,
        check=False,
    )


def test_refresh_audit_accepts_only_new_fingerprints(refresh_manifest: Path) -> None:
    output = refresh_manifest.parent / "audit.json"
    historical = refresh_manifest.parent / "historical.tsv"
    before = file_sha256(historical)
    result = run_audit(refresh_manifest, output)
    assert result.returncode == 0, result.stderr
    report = json.loads(output.read_text(encoding="utf-8"))
    assert report["passed"]
    pair = report["pairs"][0]
    assert pair["table_body_equal"]
    assert pair["current"]["rows"] == 2
    assert pair["current"]["significant_rows"] == 2
    assert set(pair["metadata_changes"]) == {"Profiles-SHA256", "Taxa-Matrix-SHA256"}
    assert file_sha256(historical) == before
    assert pair["current"]["body_sha256"] == pair["historical"]["body_sha256"]


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("body", "Table body bytes differ"),
        ("truncated", "Total-Tests-Performed: declared 2, observed 1"),
        ("selected_count", "Significant-Hits: declared 3, observed 2"),
        ("metadata", "Unexpected metadata change: Background-Universe-Definition"),
        ("profiles", "Profiles-SHA256: current export does not match supplied input"),
        ("matrix", "Taxa-Matrix-SHA256: current export does not match supplied input"),
        ("missing_fingerprint", "Profiles-SHA256: current export does not match supplied input"),
    ],
)
def test_refresh_audit_rejects_changed_or_incomplete_results(
    refresh_manifest: Path, mutation: str, message: str
) -> None:
    current = refresh_manifest.parent / "current.tsv"
    contents = current.read_bytes()
    if mutation == "body":
        mutate_table(current, "p_value", "0.25")
    elif mutation == "truncated":
        current.write_bytes(b"".join(contents.splitlines(keepends=True)[:-1]))
    elif mutation == "selected_count":
        current.write_bytes(contents.replace(b"# Significant-Hits: 2", b"# Significant-Hits: 3"))
    elif mutation == "metadata":
        current.write_bytes(contents.replace(b"mixed ranks", b"different ranks"))
    elif mutation in {"profiles", "matrix"}:
        path = refresh_manifest.parent / f"{mutation}.tsv"
        path.write_bytes(path.read_bytes() + b"\n")
    elif mutation == "missing_fingerprint":
        current.write_bytes(
            b"".join(
                line
                for line in contents.splitlines(keepends=True)
                if not line.startswith(b"# Profiles-SHA256:")
            )
        )
    output = refresh_manifest.parent / "audit.json"
    result = run_audit(refresh_manifest, output)
    assert result.returncode == 1, result.stderr
    report = json.loads(output.read_text(encoding="utf-8"))
    assert not report["passed"]
    assert any(message in error for error in report["pairs"][0]["errors"])


@pytest.mark.parametrize("hard_link", [False, True])
def test_refresh_audit_cannot_overwrite_historical_input(
    refresh_manifest: Path, hard_link: bool
) -> None:
    historical = refresh_manifest.parent / "historical.tsv"
    output = historical
    if hard_link:
        output = refresh_manifest.parent / "audit.json"
        output.hardlink_to(historical)
    before = file_sha256(historical)
    result = run_audit(refresh_manifest, output)
    assert result.returncode != 0
    assert "must not overwrite any input" in result.stderr
    assert file_sha256(historical) == before


def mutate_table(path: Path, column: str, value: str | None) -> None:
    lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
    start = next(index for index, line in enumerate(lines) if not line.startswith("#"))
    rows = list(csv.reader(lines[start:], delimiter="\t"))
    index = rows[0].index(column)
    if value is None:
        for row in rows:
            del row[index]
    else:
        rows[1][index] = value
    output = io.StringIO(newline="")
    csv.writer(output, delimiter="\t").writerows(rows)
    path.write_text("".join(lines[:start]) + output.getvalue(), encoding="utf-8", newline="")


def mutate_both_tables(manifest: Path, column: str, value: str | None) -> None:
    """Keep both bodies identical so comparison cannot conceal invalid results (#724)."""
    for filename in ("historical.tsv", "current.tsv"):
        mutate_table(manifest.parent / filename, column, value)


@pytest.mark.parametrize(
    ("column", "value", "message"),
    [
        ("p_value", "not-a-number", "invalid p_value"),
        ("p_value", "nan", "invalid p_value"),
        ("fdr_q", "2", "invalid fdr_q"),
        ("fold", "inf", "invalid fold"),
        ("expected", "-1", "invalid expected"),
        ("expected", "0", "testable rows require positive expected"),
        ("count", "-1", "invalid count"),
        ("count", "0.5", "invalid count"),
        ("count", "3", "count exceeds background_k"),
        ("testable_profiles", "3", "testable_profiles differs from metadata"),
        ("background_m", "5", "background_m or testable_profiles differs from metadata"),
        ("background_k", "0", "background_k must be positive"),
        ("background_k", "5", "background_k must be positive"),
        ("study_id", " ", "empty study_id"),
        ("inc_dec", "unknown", "invalid inc_dec"),
        ("correction_family", "unknown", "correction_family differs from FDR-Scope"),
        ("eligibility_status", "untestable", "invalid eligibility_status"),
        ("fold", "-1", "invalid fold"),
    ],
)
def test_refresh_audit_rejects_identical_malformed_rows(
    refresh_manifest: Path, column: str, value: str, message: str
) -> None:
    mutate_both_tables(refresh_manifest, column, value)
    output = refresh_manifest.parent / "audit.json"
    result = run_audit(refresh_manifest, output)
    assert result.returncode == 1, result.stdout
    report = json.loads(output.read_text(encoding="utf-8"))
    assert not report["passed"]
    assert any(message in error for error in report["pairs"][0]["errors"])


@pytest.mark.parametrize(
    "column", ["background_m", "background_k", "correction_family", "eligibility_status"]
)
def test_refresh_audit_rejects_identical_incomplete_headers(
    refresh_manifest: Path, column: str
) -> None:
    mutate_both_tables(refresh_manifest, column, None)
    output = refresh_manifest.parent / "audit.json"
    result = run_audit(refresh_manifest, output)
    assert result.returncode == 1, result.stdout
    report = json.loads(output.read_text(encoding="utf-8"))
    assert not report["passed"]
    assert "missing or duplicate complete-result columns" in report["pairs"][0]["errors"][0]


@pytest.mark.parametrize(
    ("old", "new", "message"),
    [
        (
            "P-Threshold: 1.0",
            "P-Threshold: 0.0",
            "is_significant disagrees with declared thresholds",
        ),
        ("P-Threshold: 1.0", "P-Threshold: nan", "invalid P-Threshold"),
        ("FDR: True", "FDR: False", "fdr_q differs from p_value with FDR disabled"),
        ("FDR-Scope: global", "FDR-Scope: unknown", "missing or invalid FDR-Scope"),
        ("Total-Profiles-Input: 2", "Total-Profiles-Input: 1", "Testable-Profiles exceeds"),
    ],
)
def test_refresh_audit_rejects_identical_invalid_settings(
    refresh_manifest: Path, old: str, new: str, message: str
) -> None:
    for filename in ("historical.tsv", "current.tsv"):
        path = refresh_manifest.parent / filename
        content = path.read_bytes()
        assert old.encode() in content
        path.write_bytes(content.replace(old.encode(), new.encode()))
    output = refresh_manifest.parent / "audit.json"
    result = run_audit(refresh_manifest, output)
    assert result.returncode == 1, result.stdout
    report = json.loads(output.read_text(encoding="utf-8"))
    assert any(message in error for error in report["pairs"][0]["errors"])


@pytest.mark.parametrize(
    ("refresh_manifest", "status", "rows"),
    [
        ({"scope": "profile"}, "discoveries", 2),
        ({"tail": "lower"}, "discoveries", 2),
        ({"tail": "lower", "scope": "profile"}, "discoveries", 2),
        ({"fdr": "False"}, "discoveries", 2),
        ({"tail": "lower", "scope": "profile", "fdr": "False"}, "discoveries", 2),
        ({"threshold": "0.0"}, "zero_discoveries", 2),
        ({"zero": "traits"}, "zero_eligible_tests", 0),
        ({"zero": "draws"}, "zero_eligible_tests", 0),
        ({"zero": "one_draw"}, "discoveries", 1),
        ({"unpaired": "only"}, "zero_eligible_tests", 0),
        ({"extra_trait": "yes"}, "discoveries", 4),
        ({"extra_trait": "medium"}, "discoveries", 2),
        ({"extra_trait": "zero"}, "discoveries", 2),
        pytest.param(
            {"tail": "lower", "max_fold": "inf"}, "discoveries", 2, id="lower-infinite-cutoff"
        ),
        pytest.param({"min_fold": "inf"}, "zero_discoveries", 2, id="upper-infinite-cutoff"),
    ],
    indirect=["refresh_manifest"],
)
def test_refresh_audit_preserves_supported_producer_modes_and_empty_results(
    refresh_manifest: Path, status: str, rows: int
) -> None:
    output = refresh_manifest.parent / "audit.json"
    result = run_audit(refresh_manifest, output)
    assert result.returncode == 0, result.stderr + result.stdout
    report = json.loads(output.read_text(encoding="utf-8"))
    assert report["passed"]
    current = report["pairs"][0]["current"]
    assert current["rows"] == rows
    assert current["metadata"]["Result-Status"] == status


@pytest.mark.parametrize(
    ("column", "value", "message"),
    [
        ("inc_dec", "dec", "duplicate profile/trait"),
        ("study_id", "unknown-study", "profile absent from eligible inputs"),
        ("trait", "UnknownTrait", "trait absent from eligible inputs"),
        ("background_k", "3", "background_k differs from input matrix"),
    ],
)
def test_refresh_audit_binds_identical_rows_to_input_population(
    refresh_manifest: Path, column: str, value: str, message: str
) -> None:
    """Reproduce population loss concealed by equal table bodies (#725)."""
    mutate_both_tables(refresh_manifest, column, value)
    output = refresh_manifest.parent / "audit.json"
    result = run_audit(refresh_manifest, output)
    assert result.returncode == 1, result.stdout
    report = json.loads(output.read_text(encoding="utf-8"))
    assert any(message in error for error in report["pairs"][0]["errors"])


@pytest.mark.parametrize("refresh_manifest", [{"extra_trait": "yes"}], indirect=True)
def test_refresh_audit_rejects_identical_missing_trait_family(refresh_manifest: Path) -> None:
    """Adjust declared counts too; an absent entire trait must still fail (#725)."""
    for filename in ("historical.tsv", "current.tsv"):
        path = refresh_manifest.parent / filename
        content = path.read_bytes()
        content = b"".join(
            line for line in content.splitlines(keepends=True) if b"\tTraitB\t" not in line
        )
        content = content.replace(b"Total-Tests-Performed: 4", b"Total-Tests-Performed: 2")
        content = content.replace(b"Significant-Hits: 4", b"Significant-Hits: 2")
        path.write_bytes(content)
    output = refresh_manifest.parent / "audit.json"
    result = run_audit(refresh_manifest, output)
    assert result.returncode == 1, result.stdout
    report = json.loads(output.read_text(encoding="utf-8"))
    assert any(
        "incomplete eligible profile/trait family" in error
        for error in report["pairs"][0]["errors"]
    )


@pytest.mark.parametrize(
    ("field", "replacement"),
    [("Total-Profiles-Input", "3"), ("Background-Mode", "LegacyProfileCounts")],
)
def test_refresh_audit_checks_population_metadata(
    refresh_manifest: Path, field: str, replacement: str
) -> None:
    for filename in ("historical.tsv", "current.tsv"):
        path = refresh_manifest.parent / filename
        lines = path.read_bytes().splitlines(keepends=True)
        prefix = f"# {field}: ".encode()
        path.write_bytes(
            b"".join(
                prefix + replacement.encode() + b"\n" if line.startswith(prefix) else line
                for line in lines
            )
        )
    output = refresh_manifest.parent / "audit.json"
    result = run_audit(refresh_manifest, output)
    assert result.returncode == 1, result.stdout
    report = json.loads(output.read_text(encoding="utf-8"))
    assert any(field in error for error in report["pairs"][0]["errors"])


@pytest.mark.parametrize(
    ("filename", "column", "value", "fingerprint", "message"),
    [
        ("profiles.tsv", "inc_dec", "dec", "Profiles-SHA256", "duplicate input profile"),
        ("matrix.tsv", "taxon_id", "2", "Taxa-Matrix-SHA256", "duplicate taxon_id"),
    ],
)
def test_refresh_audit_rejects_duplicate_input_identities(
    refresh_manifest: Path, filename: str, column: str, value: str, fingerprint: str, message: str
) -> None:
    path = refresh_manifest.parent / filename
    mutate_table(path, column, value)
    current = refresh_manifest.parent / "current.tsv"
    prefix = f"# {fingerprint}: ".encode()
    current.write_bytes(
        b"".join(
            prefix + file_sha256(path).encode() + b"\n" if line.startswith(prefix) else line
            for line in current.read_bytes().splitlines(keepends=True)
        )
    )
    output = refresh_manifest.parent / "audit.json"
    result = run_audit(refresh_manifest, output)
    assert result.returncode == 1, result.stdout
    report = json.loads(output.read_text(encoding="utf-8"))
    assert any(message in error for error in report["pairs"][0]["errors"])
