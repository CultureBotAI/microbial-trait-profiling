"""Verify #816 hash-binding mutants on scratch copies; print verdicts, 2026-11-10."""

from __future__ import annotations

import contextlib
import importlib.util
import io
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
from types import ModuleType

import pytest

ROOT = Path(__file__).resolve().parents[2]
TARGET = Path("src/mtp/integrated_evidence.py")
TEST = Path("tests/test_integrated_evidence.py")


def load(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {path.name}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def check(root: Path) -> None:
    names = ("mtp.integrated_evidence", "checks_816")
    saved = {name: sys.modules.get(name) for name in names}
    try:
        load(names[0], root / TARGET)
        checks = load(names[1], root / TEST)
        case = root / "case"
        case.mkdir()
        for shape in ("nested", "top_level"):
            for missing in ("profile_sha256", "experiments_sha256"):
                for change in ("absent", "renamed"):
                    checks.test_status_missing_input_keys_816(case, shape, missing, change)
            for bad in ("", "wrong", None):
                checks.test_status_rejects_invalid_hashes_816(case, shape, bad)
        invalid_inputs: tuple[object, ...] = ({}, None, [])
        for inputs in invalid_inputs:
            checks.test_present_invalid_or_empty_inputs_never_fall_back_816(case, inputs)
        checks.test_real_sensitivity_producer_binding_816(case)
        for binding in ("incomplete_inputs", "inputs_not_declared"):
            with pytest.MonkeyPatch.context() as patch:
                checks.test_strict_cli_rejects_unbound_status_816(case, patch, binding)
    finally:
        for name, previous in saved.items():
            if previous is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = previous


def main() -> None:
    originals = {path: (ROOT / path).read_bytes() for path in (TARGET, TEST)}
    mutants = [
        (
            "missing = sorted(name for name in expected_input_sha256 if name not in claimed)",
            "missing = []",
        ),
        ('source_meta["binding"] = "incomplete_inputs"', 'source_meta["binding"] = "matched"'),
        ('source_meta["missing_input_keys"] = missing', 'source_meta["missing_input_keys"] = []'),
        ("if not required or not isinstance(actual, str) or actual != required:", "if False:"),
        (
            'claimed_inputs = content.get("inputs") if "inputs" in content else content',
            'claimed_inputs = content.get("inputs")',
        ),
        (
            'claimed_inputs = content.get("inputs") if "inputs" in content else content',
            'claimed_inputs = content.get("inputs") or content',
        ),
        (
            'for name in ("profile_sha256", "experiments_sha256")\n',
            'for name in ("profile_sha256", "experiments_sha256", "full_dump_sha256")\n',
        ),
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
                    before, after = mutants[index - 1]
                    target = root / TARGET
                    source = target.read_text()
                    if before not in source:
                        raise RuntimeError(f"Missing mutation anchor {index}")
                    target.write_text(source.replace(before, after, 1))
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
                raise RuntimeError(f"Worktree changed: {path}")
    print("7/7 killed; worktree sources byte-identical")


if __name__ == "__main__":
    main()
