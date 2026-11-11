"""Check #270 scratch-copy mutants against cohort tests; print verdicts, 2026-11-10."""

from __future__ import annotations

import contextlib
import importlib.util
import io
from pathlib import Path
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = Path("scripts/create_six_shared_trait_cytoscape_heatmaps.py")
TEST = Path("tests/test_six_shared_trait_cytoscape_heatmaps.py")


def check(root: Path) -> None:
    spec = importlib.util.spec_from_file_location("checks_270", root / TEST)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load module spec for {root / TEST}")
    checks = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(checks)
    for cohort in ("oral", "skin"):
        for index, effect in enumerate((1e-7, 0.0, float("nan"), float("inf"))):
            case = root / f"{cohort}-{index}"
            case.mkdir()
            source = checks.oral_contrast_source.__wrapped__(case)
            checks.test_cohort_row_filter_and_audit_270(case, source, cohort, effect)
    for count in (0, 1):
        case = root / f"sparse-{count}"
        case.mkdir()
        source = checks.oral_contrast_source.__wrapped__(case)
        checks.test_skin_sparse_cohort_has_audit_270(case, source, count)
    checks.test_skin_exact_selection_and_source_preservation_270(source)


def main() -> None:
    originals = {path: (ROOT / path).read_bytes() for path in (SCRIPT, TEST)}
    mutants = [
        (
            'cells["displayed_in_heatmap"] & np.isfinite(effect)',
            'cells["selected_after_bh"] & np.isfinite(effect)',
        ),
        ("effect.ne(0)", "effect.round(2).ne(0)"),
        ("& np.isfinite(effect)", "& True"),
        (
            'rows["finite_nonzero_displayed_traits"] > 0',
            'rows["finite_nonzero_displayed_traits"] >= 0',
        ),
        ("vmax=1.0,", "vmax=2.0,"),
        (
            'family_size = int((source["view_key"] == '
            '"individual_conditions_direction_contrast").sum())',
            "family_size = len(cells)",
        ),
        (
            'SKIN_BODY_SITE_GROUPS = frozenset({"Skin"})',
            'SKIN_BODY_SITE_GROUPS = frozenset({"Wrong site"})',
        ),
        ("term.strip().casefold() in diseases", "False"),
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
                    target = root / SCRIPT
                    source = target.read_text()
                    # Limit shared spellings to the cohort writer, not unrelated views.
                    if index <= 6:
                        prefix, source = source.split("def write_cohort_trait_heatmaps(", 1)
                        prefix += "def write_cohort_trait_heatmaps("
                    else:
                        prefix = ""
                    if before not in source:
                        raise RuntimeError(f"Expected substring not found in {target}: {before!r}")
                    target.write_text(prefix + source.replace(before, after, 1))
                try:
                    with contextlib.redirect_stdout(io.StringIO()):
                        check(root)
                except AssertionError:
                    if index <= 0:
                        raise
                    print(f"mutant {index}: killed by AssertionError")
                else:
                    if index != 0:
                        raise AssertionError(f"mutant {index} survived")
    finally:
        for path, content in originals.items():
            if (ROOT / path).read_bytes() != content:
                raise RuntimeError(f"Worktree file modified during verification: {path}")
    print("8/8 killed; worktree sources byte-identical")


if __name__ == "__main__":
    main()
