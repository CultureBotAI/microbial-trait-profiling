"""Warn about old unmarked outputs; cleanup is explicit and exact (#1705)."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

from mtp.cellulose_decr import retired_outputs

STEMS = (
    "fig1_tufte_convergence_slope_and_overlap",
    "fig2_tufte_core_guild_dot_table",
    "fig3_tufte_why_only_stool",
    "fig4_tufte_trait_landscape_and_four_pairs",
    "fig5_tufte_three_statistical_arms",
    "fig6_tufte_combining_results_integrated_evidence",
)
EXPECTED = {
    Path("figures_tufte") / f"{stem}.{suffix}" for stem in STEMS for suffix in ("png", "svg", "pdf")
} | {
    Path("figures_tufte/fig7_tufte_convergence_network_and_venn.svg"),
    Path("figures_tufte/fig8_literature_support_four_pairs.svg"),
    Path("figures_tufte/fig9_set_theory_17_core_vs_pairwise.svg"),
    Path("figures_tufte/cellulose_decr_tufte_figures.pdf"),
    Path("cellulose_decr_slide_deck.html"),
}
ROOT = Path(__file__).resolve().parents[1]


def seed(analysis: Path) -> set[Path]:
    """Mix the complete retired set with deliberately adjacent current products."""
    current = {
        Path("figures/fig1.png"),
        Path("figures/cellulose_decr_figures.pdf"),
        Path("figures_tufte/current.svg"),
        Path("figures_tufte/FIG1_TUFTE_CONVERGENCE_SLOPE_AND_OVERLAP.jpg"),
        Path("other/cellulose_decr_slide_deck.html"),
    }
    for relative in EXPECTED | current:
        path = analysis / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f"keep exact bytes {relative}\n")
    return current


def test_cleanup_has_exactly_twenty_three_named_paths() -> None:
    assert set(retired_outputs.RETIRED_OUTPUTS) == EXPECTED
    assert len(retired_outputs.RETIRED_OUTPUTS) == 23


@pytest.mark.parametrize("relative", sorted(EXPECTED), ids=str)
def test_warn_names_each_stale_product_without_changing_it(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], relative: Path
) -> None:
    path = tmp_path / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"old unmarked output")
    assert retired_outputs.main(["warn", "--analysis-dir", str(tmp_path)]) == 0
    error = capsys.readouterr().err
    assert "WARNING" in error
    assert str(relative) in error
    assert "just cellulose-decr-cleanup-retired" in error
    assert "--delete" in error
    assert path.read_bytes() == b"old unmarked output"


def test_no_warning_or_directories_for_missing_retired_outputs(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    missing = tmp_path / "absent"
    assert retired_outputs.main(["warn", "--analysis-dir", str(missing)]) == 0
    assert not capsys.readouterr().err
    assert not missing.exists()


@pytest.mark.parametrize("delete", [False, True])
def test_real_cleanup_recipe_is_dry_by_default_and_preserves_current_files(
    tmp_path: Path, delete: bool
) -> None:
    analysis = tmp_path / "analysis"
    current = seed(analysis)
    expected_bytes = {
        relative: (analysis / relative).read_bytes() for relative in EXPECTED | current
    }
    # The shim only selects this test interpreter; the actual recipe and module execute.
    shim = tmp_path / "bin" / "uv"
    shim.parent.mkdir()
    shim.write_text('#!/bin/sh\nshift 2\nexec "$TEST_PYTHON" "$@"\n')
    shim.chmod(0o700)
    result = subprocess.run(
        [
            "just",
            "--justfile",
            str(ROOT / "Justfile"),
            "cellulose-decr-cleanup-retired",
            "--analysis-dir",
            str(analysis),
            *(["--delete"] if delete else []),
        ],
        cwd=ROOT,
        env={
            **os.environ,
            "PATH": f"{shim.parent}{os.pathsep}{os.environ['PATH']}",
            "TEST_PYTHON": sys.executable,
        },
        check=False,
        capture_output=True,
        text=True,
        timeout=15,
    )
    assert result.returncode == 0, result.stderr
    for relative, expected in expected_bytes.items():
        path = analysis / relative
        if delete and relative in EXPECTED:
            assert not path.exists(), relative
        else:
            assert path.read_bytes() == expected, relative
    prefix = "DELETE " if delete else "WOULD DELETE "
    listed = {
        Path(line.removeprefix(prefix)).relative_to(analysis)
        for line in result.stdout.splitlines()
        if line.startswith(prefix)
    }
    assert listed == EXPECTED
    assert (analysis / "figures_tufte").is_dir()


@pytest.mark.parametrize(
    "kind", ["file-symlink", "dangling-symlink", "directory", "parent-symlink", "nested-parent"]
)
def test_delete_refuses_unsafe_paths_before_removing_any_file(tmp_path: Path, kind: str) -> None:
    analysis = tmp_path / "analysis"
    seed(analysis)
    target = analysis / "cellulose_decr_slide_deck.html"
    target.unlink()
    outside = tmp_path / "outside"
    outside.write_text("external data")
    if kind == "file-symlink":
        target.symlink_to(outside)
    elif kind == "dangling-symlink":
        target.symlink_to(tmp_path / "absent-target")
    elif kind == "directory":
        target.mkdir()
    elif kind == "nested-parent":
        target.write_text("historical deck")
        (analysis / "figures_tufte").rename(tmp_path / "real-figures")
        (analysis / "figures_tufte").symlink_to(tmp_path / "real-figures", target_is_directory=True)
    else:
        analysis.rename(tmp_path / "real-analysis")
        analysis.symlink_to(tmp_path / "real-analysis", target_is_directory=True)
    with pytest.raises(ValueError, match=r"regular|symlink"):
        retired_outputs.main(["cleanup", "--analysis-dir", str(analysis), "--delete"])
    assert outside.read_text() == "external data"
    for relative in EXPECTED - {Path("cellulose_decr_slide_deck.html")}:
        assert (analysis / relative).is_file()


def test_default_recipe_warns_before_the_current_renderer(tmp_path: Path) -> None:
    # Copy the Justfile so its default relative output tree is this isolated fixture.
    (tmp_path / "Justfile").write_bytes((ROOT / "Justfile").read_bytes())
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True, timeout=10)
    analysis = tmp_path / "data/analysis/cellulose_decr"
    seed(analysis)
    shim = tmp_path / "bin" / "uv"
    shim.parent.mkdir()
    shim.write_text(
        "#!/bin/sh\nshift 2\n"
        'case "$*" in\n'
        '  *mtp.cellulose_decr.retired_outputs*) exec "$TEST_PYTHON" "$@" ;;\n'
        "  *) echo CURRENT_RENDERER >&2 ;;\n"
        "esac\n"
    )
    shim.chmod(0o700)
    # git rev-parse in the recipe requires a real HEAD, but no fixture commit is needed.
    git_shim = shim.parent / "git"
    git_shim.write_text('#!/bin/sh\ncase "$1" in rev-parse) echo fixture-head ;; esac\n')
    git_shim.chmod(0o700)
    result = subprocess.run(
        ["just", "--justfile", str(tmp_path / "Justfile"), "cellulose-decr-figures"],
        cwd=tmp_path,
        env={
            **os.environ,
            "PATH": f"{shim.parent}{os.pathsep}{os.environ['PATH']}",
            "TEST_PYTHON": sys.executable,
        },
        check=False,
        capture_output=True,
        text=True,
        timeout=15,
    )
    assert result.returncode == 0, result.stderr
    assert result.stderr.index("WARNING") < result.stderr.index("CURRENT_RENDERER")
    assert "just cellulose-decr-cleanup-retired" in result.stderr
    assert all((analysis / relative).exists() for relative in EXPECTED)
