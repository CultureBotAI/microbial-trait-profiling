"""The complete plate build validates and retains its inputs before writing (#2379)."""

from __future__ import annotations

import csv
import hashlib
import shutil
from pathlib import Path

import pytest

from mtp.ncasm_tufte_plates import build, plates_common, run_tables, story_plots
from mtp.ncasm_tufte_plates.data_status import InputFile

FIXTURE_RUN = Path(__file__).parent / "fixtures/ncasm_tufte_plates/run_tier12_off_511c891"


@pytest.fixture(autouse=True)
def _small_plates(monkeypatch: pytest.MonkeyPatch) -> None:
    """Keep the real build path, with lower raster resolution for these output tests."""
    monkeypatch.setattr(plates_common, "PLATE_DPI", 40)
    monkeypatch.setattr(build, "PLATE_DPI", 40)
    monkeypatch.setattr(story_plots, "SLIDE_DPI", 40)


@pytest.mark.parametrize("existing", [False, True], ids=["absent-output", "existing-output"])
def test_label_refusal_leaves_the_whole_build_output_untouched(
    tmp_path: Path, existing: bool
) -> None:
    run_root = tmp_path / "run"
    shutil.copytree(FIXTURE_RUN, run_root)
    screen = run_root / run_tables.DIRECTION_SCREEN
    with screen.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows[:40]:
        row["in_table1"] = "true"
    with screen.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    output = tmp_path / "plates"
    if existing:
        output.mkdir()
        (output / "keep.txt").write_bytes(b"existing output must remain untouched\n")
    before = {path.name: path.read_bytes() for path in output.iterdir()} if existing else {}

    with pytest.raises(run_tables.RunTableError, match="too little of plate 21 for the data"):
        build.build_plates(output, audit_random_marks=True, run_root=run_root)

    assert output.exists() == existing
    after = {path.name: path.read_bytes() for path in output.iterdir()} if output.exists() else {}
    assert after == before


def test_build_draws_from_the_bytes_validated_before_drawing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Changing a table after preflight cannot change the validated plot or its hash."""
    run_root = tmp_path / "run"
    shutil.copytree(FIXTURE_RUN, run_root)
    screen = run_root / run_tables.DIRECTION_SCREEN
    original_hash = hashlib.sha256(screen.read_bytes()).hexdigest()
    reads: list[Path] = []
    read = run_tables._read

    def replace_after_read(root: Path, relative: Path) -> tuple[bytes, InputFile]:
        result = read(root, relative)
        reads.append(relative)
        if relative == run_tables.DIRECTION_SCREEN:
            (root / relative).write_bytes(b"replaced after preflight\n")
        return result

    monkeypatch.setattr(run_tables, "_read", replace_after_read)
    output = tmp_path / "plates"

    manifest = build.build_plates(output, audit_random_marks=True, run_root=run_root)

    assert reads == [
        run_tables.COMPARISON_RUN,
        run_tables.ARM1_POOLED_STATS,
        run_tables.DIRECTION_SCREEN,
    ]
    inputs = manifest["plates"][20]["data_status"]["inputs"]
    hashes = {item["path"]: item["sha256"] for item in inputs}
    assert hashes[run_tables.DIRECTION_SCREEN.as_posix()] == original_hash
    assert (output / build.MANIFEST).is_file()
    assert (output / build.COMBINED_PDF).stat().st_size > 0
