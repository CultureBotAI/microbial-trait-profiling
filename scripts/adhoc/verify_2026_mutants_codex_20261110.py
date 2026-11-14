"""Prove the six #2026 documentation guards using bounded copies (2026-11-10).

Inputs: the declared docs, canonical modules and regression test; --scratch-dir.
Outputs: retained control/mutant trees, input hashes and verdicts.json. Calls the
actual pytest test function directly; no live file is edited. --canary runs one
mutant through the same path as the full six-mutant batch.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from pathlib import Path
from tempfile import mkdtemp
from types import ModuleType

ROOT = Path(__file__).resolve().parents[2]
TEST = Path("tests/test_contributor_docs_mtp_paths_2026.py")
DOCUMENTS = (
    ".claude/skills/review-open-issues/SKILL.md",
    "methods/README.md",
    "methods/integrated_evidence.md",
    "methods/table1_log_ratio_sensitivity.md",
    "methods/table1_paired_log_ratios.md",
    "docs/rosalind-six-trait-interpretation.md",
)


def load_tests(root: Path) -> ModuleType:
    """Load the copied test without leaving interpreter caches in the source tree."""
    sys.dont_write_bytecode = True
    spec = importlib.util.spec_from_file_location("checks2026", root / TEST)
    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot load the #2026 regression tests")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if root != module.ROOT:
        raise RuntimeError("The test is not bound to the requested scratch root")
    if tuple(module.CANONICAL_REFERENCES) != DOCUMENTS:
        raise RuntimeError("The regression test must cover exactly the six declared documents")
    return module


def stage(root: Path, originals: dict[Path, bytes]) -> ModuleType:
    """Copy only the bounded input inventory, then load that copy's real test."""
    root.mkdir()
    for path, content in originals.items():
        target = root / path
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open("xb") as stream:
            stream.write(content)
    return load_tests(root)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scratch-dir", type=Path, required=True)
    parser.add_argument("--canary", action="store_true")
    args = parser.parse_args()
    scratch = args.scratch_dir.resolve()
    if scratch == ROOT or scratch.is_relative_to(ROOT):
        raise ValueError("Scratch copies must be outside the live worktree")
    originals = {path: (ROOT / path).read_bytes() for path in (TEST, *map(Path, DOCUMENTS))}
    tests = load_tests(ROOT)
    for references in tests.CANONICAL_REFERENCES.values():
        for canonical in references:
            path = Path(canonical)
            if not canonical.startswith("src/mtp/") or not path.suffix == ".py":
                raise RuntimeError(f"Unexpected canonical input: {canonical}")
            if not (ROOT / path).resolve().is_relative_to(ROOT):
                raise RuntimeError(f"Canonical input escapes the worktree: {canonical}")
            originals[path] = (ROOT / path).read_bytes()
    scratch.mkdir(parents=True, exist_ok=True)
    run = Path(mkdtemp(prefix="verify2026_", dir=scratch))
    hashes = {str(path): hashlib.sha256(content).hexdigest() for path, content in originals.items()}
    (run / "inputs.json").write_text(json.dumps(hashes, indent=2) + "\n", encoding="utf-8")
    verdicts: list[dict[str, str]] = []
    documents = DOCUMENTS[:1] if args.canary else DOCUMENTS
    try:
        control = stage(run / "control", originals)
        for document in DOCUMENTS:
            control.test_active_producer_references_use_canonical_mtp(document)
        verdicts.append(
            {
                "case": "control",
                "verdict": "passed",
                "detail": "All six real test-function calls passed against the unmodified copies.",
            }
        )
        for index, document in enumerate(documents, start=1):
            root = run / f"mutant-{index}"
            copied = stage(root, originals)
            source = originals[Path(document)].decode("utf-8")
            for canonical, count in copied.CANONICAL_REFERENCES[document].items():
                if source.count(canonical) != count:
                    raise RuntimeError(f"Mutation anchor count differs: {document}: {canonical}")
                source = source.replace(canonical, canonical.replace("src/mtp/", "src/bsdbng/"))
            (root / document).write_text(source, encoding="utf-8")
            for other in DOCUMENTS:
                if other != document:
                    copied.test_active_producer_references_use_canonical_mtp(other)
            try:
                copied.test_active_producer_references_use_canonical_mtp(document)
            except AssertionError as exc:
                detail = str(exc).strip()
                if not detail.startswith(f"canonical producer reference count: {document}:"):
                    raise RuntimeError(
                        f"Unexpected assertion site for {document}: {detail}"
                    ) from exc
                verdicts.append(
                    {
                        "case": document,
                        "verdict": "killed by AssertionError",
                        "detail": detail,
                    }
                )
                print(f"KILLED {document}: {detail}")
            else:
                raise RuntimeError(f"Mutant survived: {document}")
    finally:
        for path, content in originals.items():
            if (ROOT / path).read_bytes() != content:
                raise RuntimeError(f"Live input changed during verification: {path}")
    (run / "verdicts.json").write_text(json.dumps(verdicts, indent=2) + "\n", encoding="utf-8")
    print(f"{len(documents)}/{len(documents)} killed; live inputs unchanged; evidence: {run}")


if __name__ == "__main__":
    main()
