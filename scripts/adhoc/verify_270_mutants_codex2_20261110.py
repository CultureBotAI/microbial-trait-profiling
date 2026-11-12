"""Check #270/#2312 mutants using two-file scratch copies; 2026-11-11.

Inputs: cohort heatmap producer and its shipped test file. Outputs: printed
verdicts and optional exclusive JSON evidence. Temporary files use TMPDIR.
The unchanged legacy check must miss both #2312 mutants; only the new signed,
single-trait cases may kill them. No real data or production files are modified.
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import importlib.util
import io
import json
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = Path("scripts/create_six_shared_trait_cytoscape_heatmaps.py")
TEST = Path("tests/test_six_shared_trait_cytoscape_heatmaps.py")
sys.dont_write_bytecode = True


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


LEGACY_MUTANTS = (
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
)
BOUNDARY_MUTANTS = (
    ("M9_positive_effects_only", "effect.ne(0)", "effect.gt(0)"),
    (
        "M10_requires_two_nonzero_traits",
        'rows["finite_nonzero_displayed_traits"] > 0',
        'rows["finite_nonzero_displayed_traits"] > 1',
    ),
)


def _stage(
    root: Path,
    originals: dict[Path, bytes],
    mutation: tuple[str, str] | None = None,
    *,
    cohort_only: bool = True,
) -> None:
    for path, content in originals.items():
        target = root / path
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(content)
    if mutation is None:
        return
    before, after = mutation
    target = root / SCRIPT
    source = target.read_text(encoding="utf-8")
    # Preserve the original verifier's restriction to the cohort writer.
    if cohort_only:
        prefix, source = source.split("def write_cohort_trait_heatmaps(", 1)
        prefix += "def write_cohort_trait_heatmaps("
    else:
        prefix = ""
    if before not in source:
        raise RuntimeError(f"Expected substring not found in {target}: {before!r}")
    target.write_text(prefix + source.replace(before, after, 1), encoding="utf-8")


def _boundary_check(root: Path, cohort: str, effect: float) -> None:
    spec = importlib.util.spec_from_file_location("checks_2312", root / TEST)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load module spec for {root / TEST}")
    checks = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(checks)
    case = root / "signed-single-trait"
    case.mkdir()
    source = checks.oral_contrast_source.__wrapped__(case)
    checks.test_single_signed_displayed_trait_retains_cohort_row_2312(case, source, cohort, effect)


def _boundary_cases(
    directory: Path,
    originals: dict[Path, bytes],
    mutant: tuple[str, str, str] | None,
) -> list[dict[str, object]]:
    verdicts: list[dict[str, object]] = []
    name = mutant[0] if mutant is not None else "control"
    for cohort in ("oral", "skin"):
        for effect in (-1e-7, 1e-7):
            sign = "negative" if effect < 0 else "positive"
            root = directory / f"{cohort}-{sign}"
            _stage(root, originals, None if mutant is None else (mutant[1], mutant[2]))
            expected_kill = mutant is not None and (
                name == "M10_requires_two_nonzero_traits" or effect < 0
            )
            detail = "Single displayed signed trait retained with exact values and source audit"
            exception = ""
            try:
                with contextlib.redirect_stdout(io.StringIO()):
                    _boundary_check(root, cohort, effect)
            except AssertionError as error:
                if not expected_kill:
                    raise RuntimeError(
                        f"Unexpected boundary failure: {name}/{cohort}/{sign}"
                    ) from error
                detail = str(error) or "New single-trait retention assertion rejected the mutation"
                exception = "AssertionError"
            else:
                if expected_kill:
                    raise RuntimeError(f"Boundary mutant survived: {name}/{cohort}/{sign}")
            verdict: dict[str, object] = {
                "name": name,
                "cohort": cohort,
                "effect": effect,
                "test": "test_single_signed_displayed_trait_retains_cohort_row_2312",
                "status": "killed" if expected_kill else "passed",
                "exception": exception,
                "detail": detail,
            }
            verdicts.append(verdict)
            print(f"{name}/{cohort}/{sign}: {verdict['status']} ({detail})")
    return verdicts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--canary",
        action="store_true",
        help="Run controls plus the first new mutant's legacy proof and all signed cases.",
    )
    parser.add_argument("--output", type=Path, help="Write exclusive JSON verdicts to this path.")
    arguments = parser.parse_args()
    originals = {path: (ROOT / path).read_bytes() for path in (SCRIPT, TEST)}
    verdicts: list[dict[str, object]] = []
    selected_legacy = () if arguments.canary else LEGACY_MUTANTS
    selected_boundary = BOUNDARY_MUTANTS[:1] if arguments.canary else BOUNDARY_MUTANTS
    try:
        with TemporaryDirectory(prefix="verify-270-") as directory:
            scratch = Path(directory)
            control = scratch / "legacy-control"
            _stage(control, originals)
            with contextlib.redirect_stdout(io.StringIO()):
                check(control)
            verdicts.append(
                {
                    "name": "legacy-control",
                    "status": "passed",
                    "detail": "The unchanged legacy check passed all original cohort cases",
                }
            )
            verdicts.extend(_boundary_cases(scratch / "boundary-control", originals, None))
            for index, mutation in enumerate(selected_legacy, 1):
                root = scratch / f"legacy-mutant-{index}"
                _stage(root, originals, mutation, cohort_only=index <= 6)
                try:
                    with contextlib.redirect_stdout(io.StringIO()):
                        check(root)
                except AssertionError as error:
                    detail = str(error) or "Unchanged legacy check rejected the original mutation"
                    verdicts.append(
                        {
                            "name": f"M{index}",
                            "status": "killed",
                            "exception": "AssertionError",
                            "detail": detail,
                        }
                    )
                    print(f"mutant {index}: killed by AssertionError ({detail})")
                else:
                    raise RuntimeError(f"Original mutant {index} survived")
            for mutant in selected_boundary:
                root = scratch / f"{mutant[0]}-legacy"
                _stage(root, originals, (mutant[1], mutant[2]))
                try:
                    with contextlib.redirect_stdout(io.StringIO()):
                        check(root)
                except AssertionError as error:
                    raise RuntimeError(
                        f"{mutant[0]} unexpectedly failed the legacy check"
                    ) from error
                verdicts.append(
                    {
                        "name": mutant[0],
                        "test": "unchanged legacy check",
                        "status": "survived_legacy",
                        "detail": "Legacy cases do not pin the signed single-trait boundary",
                    }
                )
                print(f"{mutant[0]}: survived unchanged legacy check")
                verdicts.extend(_boundary_cases(scratch / mutant[0], originals, mutant))
    finally:
        for path, content in originals.items():
            if (ROOT / path).read_bytes() != content:
                raise RuntimeError(f"Worktree file modified during verification: {path}")
    report = {
        "verdicts": verdicts,
        "legacy_mutants_killed": len(selected_legacy),
        "boundary_mutants_killed": len(selected_boundary),
        "source_byte_identical": True,
        "source_sha256": {
            path.as_posix(): hashlib.sha256(content).hexdigest()
            for path, content in originals.items()
        },
    }
    if arguments.output is not None:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        with arguments.output.open("x", encoding="utf-8") as handle:
            json.dump(report, handle, indent=2)
            handle.write("\n")
    print(
        f"{len(selected_legacy)} original and {len(selected_boundary)} boundary mutants killed; "
        "worktree sources byte-identical"
    )


if __name__ == "__main__":
    main()
