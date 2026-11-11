"""Verify #746 scratch copies of source/tests/methods; print verdicts, 2026-11-10."""

import argparse
import contextlib
import importlib.util
import io
from pathlib import Path
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = Path("scripts/trait_enrichment.py")
TEST = Path("tests/test_trait_enrichment.py")
DOC = Path("methods/ora_enrichment.md")


def check(root: Path) -> None:
    spec = importlib.util.spec_from_file_location("checks_746", root / TEST)
    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot load tests")
    checks = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(checks)
    case = root / "case"
    case.mkdir()
    checks.test_default_background_bytes(case)
    checks.test_profile_corpus_background(case, "t1\nt2\n")
    checks.test_background_methods_contract()
    for index, scope in enumerate(([], ["--background-scope", "all_matrix"])):
        rejection_case = root / f"rejection_{index}"
        rejection_case.mkdir()
        checks.test_all_matrix_rejects_corpus_file(rejection_case, scope)
    corpus = case / "corpus.tsv"
    corpus.write_text("t1\nt2\nunknown\n")
    try:
        checks.run_case(
            checks.load_script(),
            case,
            ["--background-scope", "profile_corpus", "--corpus-taxa-file", str(corpus)],
        )
    except ValueError:
        pass
    else:
        raise AssertionError("Unknown corpus IDs accepted")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--canary", action="store_true", help="Run control and new flag mutant only"
    )
    parser.add_argument("--scratch-dir", type=Path)
    args = parser.parse_args()
    if args.scratch_dir is not None:
        args.scratch_dir.mkdir(parents=True, exist_ok=True)
    originals = {path: (ROOT / path).read_bytes() for path in (SCRIPT, TEST, DOC)}
    mutants = [
        (
            SCRIPT,
            'if args.background_scope == "all_matrix" and args.corpus_taxa_file is not None:',
            "if False:",
        ),
        (
            SCRIPT,
            "background_taxa_universe = corpus_rows",
            "background_taxa_universe = matrix_taxa_total",
        ),
        (SCRIPT, "if taxon not in corpus_taxa:", "if False:"),
        (SCRIPT, "if unknown:", "if False:"),
        (SCRIPT, "ProfileCorpusTaxaMatrix (", "TaxaMatrix ("),
        (SCRIPT, "# Background-Scope:", "# Scope:"),
        (SCRIPT, "{file_sha256(args.corpus_taxa_file)}", "unbound"),
        (DOC, "--background-scope profile_corpus", "--background-scope other"),
    ]
    if args.canary:
        mutants = mutants[:1]
    try:
        with TemporaryDirectory(dir=args.scratch_dir) as directory:
            for index in range(len(mutants) + 1):
                root = Path(directory) / str(index)
                for path, content in originals.items():
                    target = root / path
                    target.parent.mkdir(parents=True, exist_ok=True)
                    target.write_bytes(content)
                if index:
                    path, before, after = mutants[index - 1]
                    target = root / path
                    source = target.read_text()
                    if before not in source:
                        raise RuntimeError(f"Missing anchor {index}")
                    target.write_text(source.replace(before, after))
                try:
                    with contextlib.redirect_stdout(io.StringIO()):
                        check(root)
                except AssertionError:
                    if not index:
                        raise
                    print(f"mutant {index}: killed by AssertionError")
                else:
                    if index:
                        raise RuntimeError(f"mutant {index} survived")
    finally:
        for path, content in originals.items():
            if (ROOT / path).read_bytes() != content:
                raise RuntimeError(f"Changed worktree: {path}")
    print(f"{len(mutants)}/{len(mutants)} killed; worktree targets byte-identical")


if __name__ == "__main__":
    main()
