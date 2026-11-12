"""Verify #525 in scratch copies; inputs: generator/READMEs; output: verdicts, 2026-11-09."""

from __future__ import annotations

import importlib.util
from pathlib import Path
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[2]
TEST = Path("tests/test_enrichment_followup_readme_525.py")


def main() -> None:
    spec = importlib.util.spec_from_file_location("readme_checks", ROOT / TEST)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load module spec for {ROOT / TEST}")
    checks = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(checks)
    script = checks.SCRIPT
    reports = [checks.REPORTS / name for name in checks.NAMES]
    originals = {path: (ROOT / path).read_bytes() for path in [TEST, script, *reports]}
    mutants = [
        (script, "Historical source figure:", "Source figure:"),
        (script, ".relative_to(REPO_ROOT).as_posix()", ".relative_to(REPO_ROOT)"),
        (script, "not interchangeable", "interchangeable"),
        (reports[0], "Historical source figure:", "Source figure:"),
        (reports[1], "prefixed artifacts", "artifacts"),
        (
            reports[2],
            "data/derived/enrichment_infographic.png",
            "data\\derived\\enrichment_infographic.png",
        ),
        (reports[3], "separate output directory", "same output directory"),
        (reports[3], "CultureBotAI/microbial-trait-profiling", "CultureBotAI/bsdbng/blob/legacy"),
    ]
    try:
        with TemporaryDirectory() as directory:
            base = Path(directory)
            for index in range(len(mutants) + 1):
                root = base / str(index)
                for path, content in originals.items():
                    target = root / path
                    target.parent.mkdir(parents=True, exist_ok=True)
                    target.write_bytes(content)
                if index:
                    path, before, after = mutants[index - 1]
                    target = root / path
                    source = target.read_text()
                    if before not in source:
                        raise RuntimeError(f"Expected substring not found in {target}: {before!r}")
                    target.write_text(source.replace(before, after, 1))
                    try:
                        checks.check_reports(root, root)
                    except AssertionError:
                        print(f"mutant {index}: killed by AssertionError")
                    else:
                        raise AssertionError(f"mutant {index} survived")
                else:
                    checks.check_reports(root, root)
    finally:
        for path, content in originals.items():
            if (ROOT / path).read_bytes() != content:
                raise RuntimeError(f"Worktree file modified during verification: {path}")
    print("8/8 killed; worktree targets byte-identical")


if __name__ == "__main__":
    main()
