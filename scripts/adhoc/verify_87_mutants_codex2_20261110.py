"""Verify #87 decision-record mutants in scratch copies; print verdicts, 2026-11-10."""

from __future__ import annotations

import importlib.util
from pathlib import Path
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[2]
TEST = Path("tests/test_bugsigdb_access_docs_87.py")
METHODS = Path("docs/bugsigdb-access-methods.md")
BENCHMARK = Path("docs/bugsigdb-access-benchmark.md")


def main() -> None:
    spec = importlib.util.spec_from_file_location("checks_87", ROOT / TEST)
    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot load documentation checks")
    checks = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(checks)
    originals = {path: (ROOT / path).read_bytes() for path in (TEST, METHODS, BENCHMARK)}
    mutants = [
        (METHODS, "| Freshness cadence |", "| Speed |"),
        (METHODS, "`BugSigDBExports/*.gmt`", "`other sets`"),
        (METHODS, "not a service-level guarantee", "a service-level guarantee"),
        (METHODS, "`Confounders controlled for`", "`Covariates`"),
        (METHODS, "absent from v1.2.1", "present in v1.2.1"),
        (METHODS, "2024-02-06", "2026-04-24"),
        (METHODS, "oldid=259555", "oldid=0"),
        (BENCHMARK, "not a benchmark of v1.2.1", "a benchmark of v1.2.1"),
    ]
    try:
        with TemporaryDirectory() as directory:
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
                        raise RuntimeError(f"Missing mutation anchor {index}")
                    target.write_text(source.replace(before, after, 1))
                try:
                    checks.check_docs(root)
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
                raise RuntimeError(f"Worktree changed: {path}")
    print("8/8 killed; worktree sources byte-identical")


if __name__ == "__main__":
    main()
