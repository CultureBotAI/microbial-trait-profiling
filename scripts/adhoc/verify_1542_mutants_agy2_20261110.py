#!/usr/bin/env python3
"""Verify #1542 tests kill subtitle-overflow and header-collision mutants on scratch copies.
Inputs:
  - Repository checkout containing src/mtp/cellulose_decr/tufte_plates.py
    and tests/test_cellulose_decr_figures.py.
Outputs:
  - Exit 0 when all 6 mutants fail with AssertionError and canonical files remain untouched.
Issue: #1542
Date: 2026-11-10
"""

from __future__ import annotations

import hashlib
import os
import shutil
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
WATCHED_FILES = (
    "src/mtp/cellulose_decr/tufte_plates.py",
    "tests/test_cellulose_decr_figures.py",
)


@dataclass(frozen=True)
class MutantSpec:
    mutant_id: str
    old_text: str
    new_text: str


MUTANTS: tuple[MutantSpec, ...] = (
    MutantSpec(
        mutant_id="M01_disable_subtitle_wrapping",
        old_text="textwrap.fill(line, width=136, break_long_words=False, break_on_hyphens=False)",
        new_text="line",
    ),
    MutantSpec(
        mutant_id="M02_plate2_unwrapped_single_line_overflow",
        old_text=(
            '"(93 background taxa) strengthens the overlap:\\n'
            '62.5% of lost fiber-digesting genera "'
        ),
        new_text=(
            '"(93 background taxa) strengthens the overlap: 62.5% of lost fiber-digesting genera "'
        ),
    ),
    MutantSpec(
        mutant_id="M03_plate4_unwrapped_single_line_overflow",
        old_text=(
            '"pathways (cellulose, xylan, chitin) rank #1,\\n'
            '#2, and #3 among metabolic pathways for "'
        ),
        new_text=(
            '"pathways (cellulose, xylan, chitin) rank #1, '
            '#2, and #3 among metabolic pathways for "'
        ),
    ),
    MutantSpec(
        mutant_id="M04_subtitle_collides_with_title",
        old_text="subtitle_y = 0.898 if is_multiline_subtitle else 0.885",
        new_text="subtitle_y = 0.935 if is_multiline_subtitle else 0.885",
    ),
    MutantSpec(
        mutant_id="M05_subtitle_collides_with_header_rule",
        old_text=(
            "effective_header_rule_y = min(header_rule_y, 0.846) "
            "if is_multiline_subtitle else header_rule_y"
        ),
        new_text="effective_header_rule_y = 0.875",
    ),
    MutantSpec(
        mutant_id="M06_header_rule_collides_with_panel_top",
        old_text=(
            "effective_header_rule_y = min(header_rule_y, 0.846) "
            "if is_multiline_subtitle else header_rule_y"
        ),
        new_text="effective_header_rule_y = 0.835",
    ),
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _prepare_scratch_tree(scratch_root: Path) -> None:
    for top_entry in (".git", "uv.lock", "pyproject.toml", "data", "methods"):
        src_entry = REPO_ROOT / top_entry
        if src_entry.exists():
            (scratch_root / top_entry).symlink_to(src_entry, target_is_directory=src_entry.is_dir())
    shutil.copytree(REPO_ROOT / "src", scratch_root / "src")
    for rel_file in WATCHED_FILES:
        dst = scratch_root / rel_file
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(REPO_ROOT / rel_file, dst)


def main() -> int:
    before_hashes = {rel: _sha256(REPO_ROOT / rel) for rel in WATCHED_FILES}
    env = dict(os.environ)
    env["GIT_CONFIG_GLOBAL"] = "/dev/null"
    env["GIT_CONFIG_NOSYSTEM"] = "1"
    env["MPLBACKEND"] = "Agg"

    for spec in MUTANTS:
        with tempfile.TemporaryDirectory(prefix=f"mut_{spec.mutant_id}_") as tmp_dir:
            scratch = Path(tmp_dir)
            _prepare_scratch_tree(scratch)
            target_file = scratch / "src/mtp/cellulose_decr/tufte_plates.py"
            source_text = target_file.read_text(encoding="utf-8")
            if spec.old_text not in source_text:
                raise RuntimeError(f"{spec.mutant_id}: target snippet not found")
            mutated = source_text.replace(spec.old_text, spec.new_text, 1)
            if spec.mutant_id in (
                "M02_plate2_unwrapped_single_line_overflow",
                "M03_plate4_unwrapped_single_line_overflow",
            ):
                mutated = mutated.replace(
                    "textwrap.fill(line, width=136, break_long_words=False, "
                    "break_on_hyphens=False)",
                    "line",
                    1,
                )
            target_file.write_text(mutated, encoding="utf-8")
            run_env = dict(env)
            run_env["PYTHONPATH"] = str(scratch / "src")
            out_path = scratch / "pytest_out.txt"
            argv = [
                sys.executable,
                "-m",
                "pytest",
                str(scratch / "tests/test_cellulose_decr_figures.py"),
                "-k",
                "test_tufte_plate_headers_and_subtitles_fit_inside_fixed_packet_canvas",
                "-q",
            ]
            old_cwd = Path.cwd()
            try:
                os.chdir(scratch)
                with out_path.open("wb") as out_f:
                    pid = os.posix_spawn(
                        sys.executable,
                        argv,
                        run_env,
                        file_actions=[
                            (os.POSIX_SPAWN_DUP2, out_f.fileno(), 1),
                            (os.POSIX_SPAWN_DUP2, out_f.fileno(), 2),
                        ],
                    )
                    _, status = os.waitpid(pid, 0)
            finally:
                os.chdir(old_cwd)
            returncode = os.waitstatus_to_exitcode(status)
            combined = out_path.read_text(encoding="utf-8", errors="replace")
            if returncode == 0 or "AssertionError" not in combined:
                print(
                    f"SURVIVED or wrong failure for {spec.mutant_id} (rc={returncode}):\n"
                    f"{combined}",
                    file=sys.stderr,
                )
                return 1
            print(f"KILLED {spec.mutant_id} (rc={returncode}, AssertionError confirmed)")

    after_hashes = {rel: _sha256(REPO_ROOT / rel) for rel in WATCHED_FILES}
    if before_hashes != after_hashes:
        raise RuntimeError("Canonical repository files changed during mutant verification!")
    print(f"All {len(MUTANTS)} mutants killed; canonical worktree files verified unchanged.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
