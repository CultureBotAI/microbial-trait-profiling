"""Explicit ORA background contracts and unchanged default outputs (#746)."""

import csv
import hashlib
import importlib.util
from pathlib import Path
from types import ModuleType
from unittest.mock import patch

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


def support_case(
    module: ModuleType, root: Path, *, scope: str, draws: int, count: int, tail: str
) -> str:
    profiles = root / "profiles.tsv"
    profiles.write_text(
        "study_id\texperiment_id\tinc_dec\ttaxa_count\tA\n"
        f"s\te\tinc\t{draws}\t{count}\ns\te\tdec\t0\t0\n"
    )
    matrix = root / "matrix.tsv"
    matrix.write_text("taxon_id\tA\nt1\t1\nt2\t0\nt3\t0\nt4\t1\nt5\t0\n")
    corpus = root / "corpus.tsv"
    corpus.write_text("t1\nt2\nt3\n")
    extra = (
        ["--background-scope", scope, "--corpus-taxa-file", str(corpus)]
        if scope == "profile_corpus"
        else []
    )
    module.main(
        [
            "--input",
            str(profiles),
            "--taxa-matrix",
            str(matrix),
            "--output",
            str(root / "report.txt"),
            "--tsv-output",
            str(root / "hits.tsv"),
            "--all-results-output",
            str(root / "all.tsv"),
            "--tail",
            tail,
            *extra,
        ]
    )
    return (root / "all.tsv").read_text()


@pytest.mark.parametrize(
    ("scope", "draws", "count"),
    [
        ("profile_corpus", 3, 0),
        ("all_matrix", 5, 0),
        ("all_matrix", 5, 1),
    ],
)
@pytest.mark.parametrize("tail", ["upper", "lower"])
def test_full_support_refused_before_outputs(
    tmp_path: Path, scope: str, draws: int, count: int, tail: str
) -> None:
    error = ""
    try:
        support_case(load_script(), tmp_path, scope=scope, draws=draws, count=count, tail=tail)
    except ValueError as exc:
        error = str(exc)
    assert "hypergeometric support" in error, "impossible draw must be refused, including k=0"
    for name in ("report.txt", "hits.tsv", "all.tsv"):
        assert not (tmp_path / name).exists(), "refusal must precede every output"


def test_full_support_boundary_is_inclusive(tmp_path: Path) -> None:
    error, text = "", ""
    try:
        text = support_case(
            load_script(), tmp_path, scope="profile_corpus", draws=3, count=1, tail="lower"
        )
    except ValueError as exc:
        error = str(exc)
    assert not error, f"a valid draw at the lower support boundary must pass: {error}"
    rows = list(
        csv.DictReader(
            filter(lambda line: not line.startswith("#"), text.splitlines()), delimiter="\t"
        )
    )
    assert len(rows) == 1
    assert float(rows[0]["p_value"]) == 1.0
    assert float(rows[0]["fdr_q"]) == 1.0
    assert float(rows[0]["fold"]) == 1.0


@pytest.mark.parametrize(("tail", "expected"), [("upper", (0.5, 1.0)), ("lower", (1.0, 0.5))])
def test_corpus_probabilities_and_family(
    tmp_path: Path, tail: str, expected: tuple[float, float]
) -> None:
    corpus = tmp_path / "corpus.tsv"
    corpus.write_text("t1\nt2\n")
    _, text = run_case(
        load_script(),
        tmp_path,
        ["--background-scope", "profile_corpus", "--corpus-taxa-file", str(corpus), "--tail", tail],
    )
    rows = list(
        csv.DictReader(
            filter(lambda line: not line.startswith("#"), text.splitlines()), delimiter="\t"
        )
    )
    assert len(rows) == 4
    for row in rows:
        positive = int(row["count"]) == 1
        assert float(row["p_value"]) == pytest.approx(expected[0 if positive else 1])
        assert float(row["fdr_q"]) == 1.0, "BH family includes zero-count cells"
        assert float(row["expected"]) == 0.5
        assert float(row["fold"]) == (2.0 if positive else 0.0)


def test_corpus_hash_binds_parsed_bytes(tmp_path: Path) -> None:
    corpus = tmp_path / "corpus.tsv"
    original = b"t1\nt2\n"
    replacement = b"t1\nt3\n"
    corpus.write_bytes(original)
    read_text, read_bytes = Path.read_text, Path.read_bytes
    reads: list[str] = []

    def changing_text(path: Path, *_args: object, **_kwargs: object) -> str:
        # Keep the hook narrow: only the corpus changes between parse and later hash.
        if path == corpus:
            value = original.decode() if not reads else replacement.decode()
            reads.append("text")
            corpus.write_bytes(replacement)
            return value
        return read_text(path)

    def changing_bytes(path: Path) -> bytes:
        if path == corpus:
            value = original if not reads else replacement
            reads.append("bytes")
            corpus.write_bytes(replacement)
            return value
        return read_bytes(path)

    with (
        patch.object(Path, "read_text", changing_text),
        patch.object(Path, "read_bytes", changing_bytes),
    ):
        _, text = run_case(
            load_script(),
            tmp_path,
            ["--background-scope", "profile_corpus", "--corpus-taxa-file", str(corpus)],
        )
    assert reads == ["bytes"], "corpus must be parsed and hashed from one captured buffer"
    digest = hashlib.sha256(original).hexdigest()
    for name in ("hits.tsv", "all.tsv"):
        assert f"# Corpus-Taxa-SHA256: {digest}" in (tmp_path / name).read_text()
    rows = list(
        csv.DictReader(
            filter(lambda line: not line.startswith("#"), text.splitlines()), delimiter="\t"
        )
    )
    assert {row["background_k"] for row in rows} == {"1"}, (
        "selected rows must come from the hashed buffer"
    )


def test_zero_corpus_K_is_outside_family(tmp_path: Path) -> None:
    profiles, matrix, corpus = (
        tmp_path / name for name in ("profiles.tsv", "matrix.tsv", "corpus.tsv")
    )
    profiles.write_text(
        "study_id\texperiment_id\tinc_dec\ttaxa_count\tA\tC\ns\te\tinc\t1\t1\t0\ns\te\tdec\t0\t0\t0\n"
    )
    matrix.write_text("taxon_id\tA\tC\nt1\t1\t0\nt2\t0\t0\nt3\t0\t0\nt4\t0\t1\n")
    corpus.write_text("t1\nt2\nt3\n")
    output = tmp_path / "all.tsv"
    load_script().main(
        [
            "--input",
            str(profiles),
            "--taxa-matrix",
            str(matrix),
            "--output",
            str(tmp_path / "report.txt"),
            "--all-results-output",
            str(output),
            "--background-scope",
            "profile_corpus",
            "--corpus-taxa-file",
            str(corpus),
        ]
    )
    rows = list(
        csv.DictReader(
            filter(lambda line: not line.startswith("#"), output.read_text().splitlines()),
            delimiter="\t",
        )
    )
    assert len(rows) == 1, "K=0 corpus traits and N=0 profiles are outside the family"
    assert rows[0]["trait"] == "A"
    assert float(rows[0]["p_value"]) == pytest.approx(1 / 3)
    assert float(rows[0]["fdr_q"]) == pytest.approx(1 / 3)
