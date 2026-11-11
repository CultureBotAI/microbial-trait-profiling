"""Explicit ORA background contracts and unchanged default outputs (#746)."""

import csv
import hashlib
import importlib.util
from pathlib import Path
from types import ModuleType

import pytest

ROOT = Path(__file__).resolve().parents[1]


def load_script(root: Path = ROOT) -> ModuleType:
    spec = importlib.util.spec_from_file_location("ora_746", root / "scripts/trait_enrichment.py")
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_case(module: ModuleType, root: Path, extra: list[str]) -> tuple[bytes, str]:
    profiles = root / "profiles.tsv"
    profiles.write_text(
        "study_id\texperiment_id\tinc_dec\ttaxa_count\tA\tB\n"
        "s\te\tinc\t1\t1\t0\ns\te\tdec\t1\t0\t1\n"
    )
    matrix = root / "matrix.tsv"
    matrix.write_text("taxon_id\tA\tB\nt1\t1\t0\nt2\t0\t1\nt3\t1\t1\nt4\t0\t1\n")
    outputs = [root / name for name in ("report.txt", "hits.tsv", "all.tsv")]
    module.main(
        [
            "--input",
            str(profiles),
            "--taxa-matrix",
            str(matrix),
            "--output",
            str(outputs[0]),
            "--tsv-output",
            str(outputs[1]),
            "--all-results-output",
            str(outputs[2]),
            *extra,
        ]
    )
    return b"\0".join(path.read_bytes() for path in outputs), outputs[2].read_text()


def test_default_background_bytes(tmp_path: Path) -> None:
    module = load_script()
    payload, text = run_case(module, tmp_path, [])
    expected_sha = "".join(
        (
            "ce54c",
            "89e87",
            "6a040",
            "af051",
            "80f39",
            "0fa9f",
            "94f13",
            "20087",
            "a7cdd",
            "23e38",
            "91b4b",
            "50802",
            "2ccf",
        )
    )
    assert hashlib.sha256(payload).hexdigest() == expected_sha
    assert "# Background-Mode: TaxaMatrix (matrix.tsv)" in text
    explicit, _ = run_case(module, tmp_path, ["--background-scope", "all_matrix"])
    assert payload == explicit


@pytest.mark.parametrize("contents", ["t1\nt2\n", "taxon_id\tnote\nt1\ta\nt2\tb\nt1\tc\n"])
def test_profile_corpus_background(tmp_path: Path, contents: str) -> None:
    corpus = tmp_path / "corpus.tsv"
    corpus.write_text(contents)
    payload, text = run_case(
        load_script(),
        tmp_path,
        ["--background-scope", "profile_corpus", "--corpus-taxa-file", str(corpus)],
    )
    rows = list(
        csv.DictReader(
            [line for line in text.splitlines() if not line.startswith("#")], delimiter="\t"
        )
    )
    assert len(rows) == 4
    assert {row["background_m"] for row in rows} == {"2"}
    assert {row["background_k"] for row in rows} == {"1"}
    assert "# Background-Mode: ProfileCorpusTaxaMatrix (matrix.tsv; corpus=corpus.tsv)" in text
    assert "# Background-Scope: profile_corpus" in text
    assert hashlib.sha256(corpus.read_bytes()).hexdigest() in text
    assert b"background_scope=profile_corpus" in payload


@pytest.mark.parametrize("contents", [None, "", "taxon_id\n", "unknown\n", "t1\nunknown\n"])
def test_profile_corpus_fails_closed(tmp_path: Path, contents: str | None) -> None:
    extra = ["--background-scope", "profile_corpus"]
    if contents is not None:
        corpus = tmp_path / "corpus.tsv"
        corpus.write_text(contents)
        extra += ["--corpus-taxa-file", str(corpus)]
    with pytest.raises(ValueError, match=r"corpus|taxon"):
        run_case(load_script(), tmp_path, extra)
    assert not (tmp_path / "all.tsv").exists()


def test_background_methods_contract() -> None:
    text = (ROOT / "methods/ora_enrichment.md").read_text()
    for required in (
        "--background-scope all_matrix",
        "--background-scope profile_corpus",
        "--corpus-taxa-file",
        "16S",
        "WMS",
        "across hosts",
        "hypergeometric",
        "K",
        "M",
        "#746",
        "not in vivo",
    ):
        assert required in text


def test_missing_corpus_file_is_value_error(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="corpus taxon file is missing"):
        load_script().read_corpus_taxa(tmp_path / "absent.tsv")


@pytest.mark.parametrize("scope", [[], ["--background-scope", "all_matrix"]])
def test_all_matrix_rejects_corpus_file(tmp_path: Path, scope: list[str]) -> None:
    corpus = tmp_path / "corpus.tsv"
    corpus.write_text("t1\nt2\n")
    error = ""
    try:
        run_case(load_script(), tmp_path, [*scope, "--corpus-taxa-file", str(corpus)])
    except ValueError as exc:
        error = str(exc)
    assert error == "--corpus-taxa-file requires --background-scope profile_corpus"
    for name in ("report.txt", "hits.tsv", "all.tsv"):
        assert not (tmp_path / name).exists()
