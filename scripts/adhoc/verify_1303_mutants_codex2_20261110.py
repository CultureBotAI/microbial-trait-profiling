"""Verify unavailable-revision contracts in scratch copies (#1303, 2026-11-10).

Inputs: repository helper/tests, --scratch-dir and optional --canary.
Outputs: persistent per-case source copies and verdicts.json; originals unchanged.
"""

import argparse
import importlib.util
import json
from pathlib import Path
from tempfile import mkdtemp

import pytest

ROOT = Path(__file__).resolve().parents[2]
HELPER = Path("scripts/adhoc/issue1295_1300_probe_loader_mutants.py")
TEST = Path("tests/test_issue1295_1300_probe_loader_mutants_1303.py")


def check(root: Path) -> None:
    spec = importlib.util.spec_from_file_location("checks1303", root / TEST)
    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot load tests")
    checks = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(checks)
    # The test loader's default ROOT is bound from its isolated file location.
    for index, missing in enumerate((("head",), ("base",), ("head", "base"), ())):
        case = root / f"summary_{index}"
        case.mkdir()
        with pytest.MonkeyPatch.context() as patch:
            checks.test_summary_preserves_available_cases(case, patch, missing)
    for mode in ("missing", "git_other", "permission"):
        with pytest.MonkeyPatch.context() as patch:
            checks.test_revision_errors_are_narrow(patch, mode)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scratch-dir", type=Path, required=True)
    parser.add_argument("--canary", action="store_true")
    args = parser.parse_args()
    args.scratch_dir.mkdir(parents=True, exist_ok=True)
    run = Path(mkdtemp(prefix="verify1303_", dir=args.scratch_dir))
    originals = {path: (ROOT / path).read_bytes() for path in (HELPER, TEST)}
    mutants = [
        ('"status": "unavailable"', '"status": "verified"'),
        ('"verified": False', '"verified": True'),
        ('return int(not all(row["verified"] for row in summary))', "return 0"),
        ("for case in args.cases or CASES]", "for case in (args.cases or CASES)[:1]]"),
        ('"selected": 0', '"selected": 1'),
        ("except GitCommandError as exc:", "except PermissionError as exc:"),
        ("if exc.status == 1:", "if True:"),
        ('"revision": str(exc)', '"revision": "wrong"'),
    ]
    if args.canary:
        mutants = mutants[:1]
    verdicts = []
    for index in range(len(mutants) + 1):
        root = run / str(index)
        for path, content in originals.items():
            target = root / path
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(content)
        if index:
            before, after = mutants[index - 1]
            source = (root / HELPER).read_text()
            if source.count(before) != 1:
                raise RuntimeError(f"Nonunique mutant anchor {index}")
            (root / HELPER).write_text(source.replace(before, after))
        try:
            check(root)
        except AssertionError:
            if not index:
                raise
            verdict = "killed: AssertionError"
        else:
            if index:
                raise RuntimeError(f"Surviving mutant {index}")
            verdict = "control passed"
        verdicts.append({"case": index, "verdict": verdict})
        print(index, verdict)
    for path, content in originals.items():
        if (ROOT / path).read_bytes() != content:
            raise RuntimeError(f"Worktree changed: {path}")
    (run / "verdicts.json").write_text(json.dumps(verdicts, indent=2) + "\n")
    print(f"{len(mutants)}/{len(mutants)} killed; evidence: {run}")


if __name__ == "__main__":
    main()
