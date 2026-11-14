"""Verify #746/#2308-#2310 scratch mutants; print assertion kills, 2026-11-11.

Inputs: this checkout's source, tests and methods; outputs: temporary bounded
copies under --scratch-dir and stdout verdicts. Live worktree bytes stay intact.
"""

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
    for background_scope, draws, count in (
        ("profile_corpus", 3, 0),
        ("all_matrix", 5, 0),
        ("all_matrix", 5, 1),
    ):
        for tail in ("upper", "lower"):
            support_case = root / f"support-{background_scope}-{draws}-{count}-{tail}"
            support_case.mkdir()
            checks.test_full_support_refused_before_outputs(
                support_case, background_scope, draws, count, tail
            )
    boundary_case = root / "inclusive-boundary"
    boundary_case.mkdir()
    checks.test_full_support_boundary_is_inclusive(boundary_case)
    for tail, expected in (("upper", (0.5, 1.0)), ("lower", (1.0, 0.5))):
        probability_case = root / f"probabilities-{tail}"
        probability_case.mkdir()
        checks.test_corpus_probabilities_and_family(probability_case, tail, expected)
    binding_case = root / "corpus-byte-binding"
    binding_case.mkdir()
    checks.test_corpus_hash_binds_parsed_bytes(binding_case)
    zero_family_case = root / "zero-corpus-family"
    zero_family_case.mkdir()
    checks.test_zero_corpus_K_is_outside_family(zero_family_case)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--canary", action="store_true", help="Run control and new flag mutant only"
    )
    parser.add_argument("--scratch-dir", type=Path)
    parser.add_argument(
        "--baseline-script",
        type=Path,
        help="Optional read-only historical script for six-mode default byte parity",
    )
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
        (SCRIPT, "{corpus_sha256}", "unbound"),
        (DOC, "--background-scope profile_corpus", "--background-scope other"),
        (SCRIPT, "if not lower <= count <= upper:", "if count > upper:"),
        (
            SCRIPT,
            "if not lower <= count <= upper:",
            "if count > 0 and not lower <= count <= upper:",
        ),
        (
            SCRIPT,
            "if not lower <= count <= upper:",
            "if not (lower <= count <= upper and (lower == 0 or count != lower)):",
        ),
        (
            SCRIPT,
            "row_counts, background_taxa_universe, trait_totals, n_draws",
            "row_counts, matrix_taxa_total, trait_totals, n_draws",
        ),
        (SCRIPT, "p_values = hypergeom_cdf_many(", "p_values = hypergeom_sf_many("),
        (SCRIPT, "{corpus_sha256}", "{file_sha256(args.corpus_taxa_file)}"),
        (
            SCRIPT,
            "candidate_traits = [t for t in trait_columns if background_counts[t] > 0]",
            "candidate_traits = [t for t in trait_columns if background_counts[t] >= 0]",
        ),
    ]
    if args.canary:
        mutants = mutants[:1]
    try:
        with TemporaryDirectory(dir=args.scratch_dir) as directory:
            if args.baseline_script is not None:
                baseline_spec = importlib.util.spec_from_file_location(
                    "baseline_746", args.baseline_script
                )
                tests_spec = importlib.util.spec_from_file_location(
                    "parity_checks_746", ROOT / TEST
                )
                if (
                    baseline_spec is None
                    or baseline_spec.loader is None
                    or tests_spec is None
                    or tests_spec.loader is None
                ):
                    raise RuntimeError("Cannot load baseline parity inputs")
                baseline = importlib.util.module_from_spec(baseline_spec)
                tests = importlib.util.module_from_spec(tests_spec)
                baseline_spec.loader.exec_module(baseline)
                tests_spec.loader.exec_module(tests)
                current = tests.load_script()
                parity_root = Path(directory) / "parity"
                parity_root.mkdir()
                for tail in ("upper", "lower"):
                    for mode in ("none", "global", "profile"):
                        extra = ["--tail", tail] + (
                            ["--no-fdr"] if mode == "none" else ["--fdr", "--fdr-scope", mode]
                        )
                        with contextlib.redirect_stdout(io.StringIO()):
                            before, _ = tests.run_case(baseline, parity_root, extra)
                            after, _ = tests.run_case(current, parity_root, extra)
                        if before != after:
                            raise AssertionError(
                                f"Default output byte parity failed: {tail}/{mode}"
                            )
                        print(f"default byte parity {tail}/{mode}: passed ({len(after)} bytes)")
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
