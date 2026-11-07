"""Warn about, or explicitly remove, the 23 historical output files (#1705).

Cleanup is a dry run unless --delete is supplied. No directory, wildcard,
symlink or current figure is a deletion target. This module renders nothing.
"""

from __future__ import annotations

import argparse
import stat
import sys
from collections.abc import Sequence
from pathlib import Path

RETIRED_OUTPUTS = (
    *(
        Path("figures_tufte") / f"{stem}.{suffix}"
        for stem in (
            "fig1_tufte_convergence_slope_and_overlap",
            "fig2_tufte_core_guild_dot_table",
            "fig3_tufte_why_only_stool",
            "fig4_tufte_trait_landscape_and_four_pairs",
            "fig5_tufte_three_statistical_arms",
            "fig6_tufte_combining_results_integrated_evidence",
        )
        for suffix in ("svg", "png", "pdf")
    ),
    Path("figures_tufte/fig7_tufte_convergence_network_and_venn.svg"),
    Path("figures_tufte/fig8_literature_support_four_pairs.svg"),
    Path("figures_tufte/fig9_set_theory_17_core_vs_pairwise.svg"),
    Path("figures_tufte/cellulose_decr_tufte_figures.pdf"),
    Path("cellulose_decr_slide_deck.html"),
)


def existing_outputs(analysis_dir: Path) -> list[Path]:
    """Inspect only the literal historical names, including dangling symlinks."""
    return [
        analysis_dir / relative
        for relative in RETIRED_OUTPUTS
        if (analysis_dir / relative).exists() or (analysis_dir / relative).is_symlink()
    ]


def warn(analysis_dir: Path) -> None:
    """Leave existing files in place, warning that old copies may lack the banner."""
    paths = existing_outputs(analysis_dir)
    if not paths:
        return
    print(
        "WARNING: retired historical cellulose outputs remain; older copies may be "
        "unmarked and their numbers are unverified (#1572). No files were removed.",
        file=sys.stderr,
    )
    for path in paths:
        print(f"  {path}", file=sys.stderr)
    print(
        "Preview cleanup with just cellulose-decr-cleanup-retired; "
        "add --delete to remove exactly the listed historical files. "
        "Use --analysis-dir to select a non-default output directory.",
        file=sys.stderr,
    )


def cleanup(analysis_dir: Path, *, delete: bool = False) -> list[Path]:
    """Preview or delete exact regular-file targets with non-symlink parents.

    Validate every candidate before reporting or removing any, so an unsafe target
    refuses both modes (#1718). The caller must not concurrently replace files.
    Empty directories and all unlisted files are preserved.
    """
    paths = existing_outputs(analysis_dir)
    for path in paths:
        if not stat.S_ISREG(path.lstat().st_mode):
            raise ValueError(f"Cleanup requires a regular non-symlink file: {path}")
        for parent in (path.parent, *path.parent.parents):
            if parent.is_symlink():
                raise ValueError(f"Cleanup refuses a symlink parent: {parent}")
    for path in paths:
        print(f"{'DELETE' if delete else 'WOULD DELETE'} {path}")
        if delete:
            path.unlink()
    return paths


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("warn", "cleanup"))
    parser.add_argument("--analysis-dir", type=Path, default=Path("data/analysis/cellulose_decr"))
    parser.add_argument(
        "--delete", action="store_true", help="Remove only the 23 named retired files"
    )
    args = parser.parse_args(argv)
    if args.action == "warn":
        if args.delete:
            parser.error("--delete is only valid with cleanup")
        warn(args.analysis_dir)
    else:
        cleanup(args.analysis_dir, delete=args.delete)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
