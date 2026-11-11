"""Verify #1740 scratch-copy mutants; inputs: tests/methods; output: verdicts, 2026-11-09."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

import pytest

ROOT = Path(__file__).resolve().parents[2]
TEST = Path("tests/test_cellulose_numeric_snapshot_capture.py")
DOC = Path("methods/cellulose_decr.md")


def check(root: Path) -> None:
    spec = importlib.util.spec_from_file_location("pin_mutant", root / TEST)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load module spec for {root / TEST}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    try:
        spec.loader.exec_module(module)
        with pytest.MonkeyPatch.context() as patch:
            module.test_shallow_ci_disclosure(patch)
    finally:
        sys.modules.pop(spec.name, None)


def main() -> None:
    originals = {path: (ROOT / path).read_bytes() for path in (TEST, DOC)}
    mutants = [
        (DOC, "only full local clones", "all clones"),
        (DOC, "fetch-depth: 1", "fetch-depth: 0"),
        (TEST, "Real-Git object pin: only full local clones check it; ", "Pin: "),
        (TEST, "CI uses shallow checkouts (fetch-depth: 1) and skips it", "Skipped"),
        (TEST, "Check the real-Git pin in only full local clones.", "Check the pin."),
        (
            TEST,
            'pytest.skip(\n                "Real-Git',
            'repository.rev_parse("0672ee7d")\n            pytest.skip(\n                "Real-Git',
        ),
    ]
    try:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            for path, content in originals.items():
                (root / path).parent.mkdir(parents=True, exist_ok=True)
                (root / path).write_bytes(content)
            check(root)
            for index, (path, before, after) in enumerate(mutants, 1):
                # A fresh directory prevents Python bytecode reuse between mutants.
                scratch = root / str(index)
                for relative, content in originals.items():
                    target = scratch / relative
                    target.parent.mkdir(parents=True, exist_ok=True)
                    target.write_bytes(content)
                target = scratch / path
                source = target.read_text()
                if before not in source:
                    raise RuntimeError(f"Expected substring not found in {target}: {before!r}")
                target.write_text(source.replace(before, after, 1))
                try:
                    check(scratch)
                except AssertionError:
                    print(f"mutant {index}: killed by AssertionError")
                else:
                    raise AssertionError(f"mutant {index} survived")
    finally:
        for path, content in originals.items():
            if (ROOT / path).read_bytes() != content:
                raise RuntimeError(f"Worktree file modified during verification: {path}")
    print("6/6 killed; worktree targets byte-identical")


if __name__ == "__main__":
    main()
