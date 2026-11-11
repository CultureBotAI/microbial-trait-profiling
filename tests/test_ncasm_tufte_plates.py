"""The NCASM plate builders ported from slide-wizard: outputs, repeatability, hygiene (#2264).

None of these tests checks a number on a plate against data. The plates read no data
files; #775-#779 record what is wrong with their content.
"""

from __future__ import annotations

import json
import subprocess
import sys
import xml.etree.ElementTree as ET
from collections.abc import Iterator
from pathlib import Path

import matplotlib as mpl
import pytest

from mtp.ncasm_tufte_plates import build, plates_common, plates_part1, plates_part2, story_plots

STORED_TEXT = (
    Path(__file__).parent / "fixtures/ncasm_tufte_plates/stored_plate_text_2026-10-22.json"
)
TEST_DPI = 40
Built = tuple[Path, list[build.ManifestEntry]]


def _svg_text(svg: Path) -> list[str]:
    root = ET.parse(svg).getroot()
    return ["".join(node.itertext()) for node in root.iter("{http://www.w3.org/2000/svg}text")]


@pytest.fixture(scope="module")
def built(tmp_path_factory: pytest.TempPathFactory) -> Iterator[Built]:
    """All 30 plates, drawn once at 40 dpi: the layout is the same and the build is quick.

    The resolution stays patched for the whole module, so the repeatability test draws
    at the same size."""
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(plates_common, "PLATE_DPI", TEST_DPI)
        patch.setattr(build, "PLATE_DPI", TEST_DPI)
        patch.setattr(story_plots, "SLIDE_DPI", TEST_DPI)
        out_dir = tmp_path_factory.mktemp("ncasm_plates")
        yield out_dir, build.build_plates(out_dir)


def test_build_writes_thirty_plates_a_manifest_and_one_combined_pdf(built: Built) -> None:
    out_dir, manifest = built

    assert [entry["plate_num"] for entry in manifest] == list(range(1, 31))
    assert [entry["roman"] for entry in manifest] == plates_common.ROMAN
    assert len({entry["stem"] for entry in manifest}) == 30
    for entry in manifest:
        assert entry["stem"].startswith(f"plate_{entry['plate_num']:02d}_")
        assert entry["speaker_notes"].strip()
        for key in ("png_path", "svg_path"):
            path = out_dir / entry[key]
            assert not Path(entry[key]).is_absolute()
            assert path.is_file()
            assert path.stat().st_size > 0
    # Plates I-XX also have a vector PDF each; the data-story plots keep theirs in a subfolder.
    assert [entry["pdf_path"] is not None for entry in manifest] == [True] * 20 + [False] * 10
    for entry in manifest[:20]:
        pdf_path = entry["pdf_path"]
        assert pdf_path is not None
        assert (out_dir / pdf_path).stat().st_size > 0

    assert sorted(path.name for path in out_dir.glob("*.pdf") if "NCASM" in path.name) == [
        build.COMBINED_PDF
    ]
    assert (out_dir / build.COMBINED_PDF).stat().st_size > 0
    story_dir = out_dir / build.STORY_DIR
    for suffix in ("png", "svg", "pdf"):
        assert len(list(story_dir.glob(f"plot_*.{suffix}"))) == 10

    on_disk = json.loads((out_dir / build.MANIFEST).read_text(encoding="utf-8"))
    assert on_disk == manifest


def test_data_story_plates_are_copies_of_the_data_story_plots(built: Built) -> None:
    out_dir, manifest = built
    for entry in manifest[20:]:
        plot_stem = entry["stem"].split("_", 2)[2]
        for suffix in ("png", "svg"):
            plate = (out_dir / f"{entry['stem']}.{suffix}").read_bytes()
            plot = (out_dir / build.STORY_DIR / f"{plot_stem}.{suffix}").read_bytes()
            assert plate == plot


def test_plates_one_to_twenty_draw_the_text_of_the_stored_renders(built: Built) -> None:
    """The port's text equals the 2026-10-22 renders' text, node for node.

    This is a record of what the plates said when they were made. A deliberate change to
    a plate's wording or numbers has to change the stored reference in the same commit."""
    out_dir, manifest = built
    stored = json.loads(STORED_TEXT.read_text(encoding="utf-8"))

    assert [entry["stem"] for entry in manifest[:20]] == list(stored)
    for entry in manifest[:20]:
        assert _svg_text(out_dir / entry["svg_path"]) == stored[entry["stem"]]["text"]


def test_plate_canvas_is_the_width_the_stored_renders_were_drawn_at(built: Built) -> None:
    out_dir, manifest = built
    assert (plates_common.PLATE_WIDTH_IN, plates_common.PLATE_HEIGHT_IN) == (13.33, 7.5)
    for entry in (manifest[0], manifest[20]):
        root = ET.parse(out_dir / entry["svg_path"]).getroot()
        assert (root.get("width"), root.get("height")) == ("959.76pt", "540pt")


def test_drawing_a_plate_twice_writes_the_same_bytes(built: Built, tmp_path: Path) -> None:
    out_dir, manifest = built
    with mpl.style.context("default"):
        plates_part1.build_plate_01(tmp_path)
        plates_part2.build_plate_14(tmp_path)
        story = story_plots.PLOT_BUILDERS[0](tmp_path / build.STORY_DIR)

    for entry in (manifest[0], manifest[13]):
        for suffix in ("png", "svg", "pdf"):
            name = f"{entry['stem']}.{suffix}"
            assert (tmp_path / name).read_bytes() == (out_dir / name).read_bytes(), name
    for suffix in ("png", "svg", "pdf"):
        name = f"{story['stem']}.{suffix}"
        again = (tmp_path / build.STORY_DIR / name).read_bytes()
        assert again == (out_dir / build.STORY_DIR / name).read_bytes(), name


def test_building_leaves_matplotlib_settings_as_they_were(built: Built, tmp_path: Path) -> None:
    del built  # drawn under the same patch; this test only needs the low resolution
    before = dict(mpl.rcParams)
    with mpl.style.context("default"):
        plates_part1.build_plate_01(tmp_path)
    assert dict(mpl.rcParams) == before


def test_main_refuses_a_folder_that_already_holds_files(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    (tmp_path / "plate_01_someone_elses.png").write_bytes(b"not ours")

    assert build.main(["--output-dir", str(tmp_path)]) == 1

    assert "--overwrite" in capsys.readouterr().out
    assert [path.name for path in tmp_path.iterdir()] == ["plate_01_someone_elses.png"]


def test_default_output_folder_is_not_the_decks_plate_folder() -> None:
    assert Path("data/derived/ncasm_tufte_plates") == build.DEFAULT_OUTPUT_DIR
    assert build.parse_args([]).output_dir != Path("data/derived/tufte_plates")


def test_importing_the_builders_writes_nothing_and_leaves_the_import_path_alone(
    tmp_path: Path,
) -> None:
    """The originals created two folders and edited ``sys.path`` on import."""
    package = Path(build.__file__).parent
    code = "\n".join(
        [
            "import json, sys",
            "from pathlib import Path",
            "import matplotlib",
            "path_before = list(sys.path)",
            "rc_before = dict(matplotlib.rcParams)",
            "import mtp.ncasm_tufte_plates.build",
            "import mtp.ncasm_tufte_plates.__main__",
            "print(json.dumps({",
            "    'path_unchanged': sys.path == path_before,",
            "    'rc_unchanged': dict(matplotlib.rcParams) == rc_before,",
            "    'cwd': sorted(p.name for p in Path.cwd().iterdir()),",
            "}))",
        ]
    )
    listing_before = sorted(path.name for path in package.iterdir() if path.name != "__pycache__")

    result = subprocess.run(
        [sys.executable, "-c", code], cwd=tmp_path, check=True, capture_output=True, text=True
    )

    assert json.loads(result.stdout) == {"path_unchanged": True, "rc_unchanged": True, "cwd": []}
    listing_after = sorted(path.name for path in package.iterdir() if path.name != "__pycache__")
    assert listing_after == listing_before
