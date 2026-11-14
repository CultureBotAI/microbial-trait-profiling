"""Require every real infographic mutation anchor to match current source (#2336)."""

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
VERIFIERS = (
    ("verify_2035_2047_mutants_agy_20261111.py", "rel_path"),
    ("verify_1922_mutants_agy_20261108.py", "target_rel"),
)


def _anchor_cases() -> list[tuple[str, Path, str]]:
    """Read the verifiers' specifications without invoking their mutation runners."""
    cases: list[tuple[str, Path, str]] = []
    for filename, target_attribute in VERIFIERS:
        path = ROOT / "scripts/adhoc" / filename
        spec = importlib.util.spec_from_file_location(f"anchor2336_{path.stem}", path)
        assert spec is not None
        assert spec.loader is not None
        verifier = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = verifier
        spec.loader.exec_module(verifier)
        mutants = verifier.MUTANTS
        assert mutants, f"{filename} must define mutation specifications"
        for mutant in mutants:
            name = mutant.name
            target = getattr(mutant, target_attribute)
            anchor = mutant.old
            assert isinstance(name, str)
            assert isinstance(target, Path)
            assert isinstance(anchor, str)
            assert anchor, f"{filename}:{name} must define a nonempty source anchor"
            cases.append((f"{path.stem}:{name}", ROOT / target, anchor))
    return cases


ANCHOR_CASES = _anchor_cases()


@pytest.mark.parametrize(
    ("case_name", "target_path", "anchor"),
    ANCHOR_CASES,
    ids=[case[0] for case in ANCHOR_CASES],
)
def test_real_mutation_anchor_matches_current_source_once(
    case_name: str, target_path: Path, anchor: str
) -> None:
    source = target_path.read_text(encoding="utf-8")
    count = source.count(anchor)
    assert count == 1, (
        f"{case_name}: expected exactly one anchor in {target_path.relative_to(ROOT)}, "
        f"found {count}"
    )
