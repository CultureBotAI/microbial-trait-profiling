"""The NCASM plate builders ported from slide-wizard: outputs, repeatability, hygiene (#2264).

None of these tests checks a number on a plate against data. The plates read no data
files; #775-#779 record what is wrong with their content.
"""

from __future__ import annotations

import ast
import json
import os
import subprocess
import sys
import xml.etree.ElementTree as ET
from collections.abc import Iterator
from pathlib import Path

import matplotlib as mpl
import matplotlib.image as mpimg
import numpy as np
import pytest

from mtp.ncasm_tufte_plates import (
    build,
    data_status,
    plates_common,
    plates_part1,
    plates_part2,
    story_plots,
)

STORED_TEXT = (
    Path(__file__).parent / "fixtures/ncasm_tufte_plates/stored_plate_text_2026-10-22.json"
)
TEST_DPI = 40
Built = tuple[Path, list[build.ManifestEntry]]
PACKAGE = Path(build.__file__).parent
CODE_COMMIT = "0123456789abcdef0123456789abcdef01234567+dirty"


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
        manifest = build.build_plates(out_dir, audit_random_marks=True, code_commit=CODE_COMMIT)
        yield out_dir, manifest["plates"]


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
    assert on_disk == {
        "schema_version": 1,
        "generator": "python -m mtp.ncasm_tufte_plates",
        "code_commit": CODE_COMMIT,
        "combined_pdf": build.COMBINED_PDF,
        "plates": manifest,
    }


def test_manifest_says_where_each_plates_numbers_and_marks_come_from(built: Built) -> None:
    """No plate reads a data file yet: #776 for the numbers, #775 for the random points."""
    _out_dir, manifest = built
    for entry in manifest[:20]:
        assert entry["data_status"] == {
            "numbers": "typed_literals",
            "marks": "typed_literals",
            "inputs": [],
            "issues": [776, 777, 778, 779],
        }
    for entry in manifest[20:]:
        assert entry["data_status"] == {
            "numbers": "typed_literals",
            "marks": "random_generator",
            "inputs": [],
            "issues": [775, 776, 777, 778, 779],
        }


def _banner_pixels(png: Path) -> int:
    """Pixels of the banner's fill colour in the strip above the plate header."""
    image = mpimg.imread(png)[..., :3]
    strip = image[: max(1, round(image.shape[0] * 0.025))]
    accent = np.array([0x9E, 0x2A, 0x2B]) / 255.0
    return int(np.count_nonzero(np.abs(strip - accent).max(axis=2) < 0.02))


def test_plates_drawn_from_random_points_carry_the_banner_and_the_others_do_not(
    built: Built,
) -> None:
    out_dir, manifest = built
    assert "random draws, not data (#775)" in data_status.RANDOM_MARKS_BANNER
    for entry in manifest:
        random_marks = entry["data_status"]["marks"] == "random_generator"
        banner = _banner_pixels(out_dir / entry["png_path"])
        assert (banner > 50) == random_marks, (entry["stem"], banner)


def test_build_stops_before_drawing_when_random_points_are_not_asked_for(tmp_path: Path) -> None:
    out_dir = tmp_path / "plates"

    with pytest.raises(build.RandomMarksError, match="--audit-random-marks") as refusal:
        build.build_plates(out_dir)

    assert "plates 21, 22, 23, 24, 25, 26, 27, 28, 29, 30 " in str(refusal.value)
    assert not out_dir.exists()


def test_main_reports_the_refusal_and_exits_nonzero(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    out_dir = tmp_path / "plates"

    assert build.main(["--output-dir", str(out_dir)]) == 1

    assert "--audit-random-marks" in capsys.readouterr().out
    assert not out_dir.exists()


def _is_random_generator(node: ast.AST) -> bool:
    if isinstance(node, ast.Attribute):
        if node.attr in {"default_rng", "RandomState", "Generator"}:
            return True
        return (
            node.attr == "random"
            and isinstance(node.value, ast.Name)
            and node.value.id in {"np", "numpy"}
        )
    if isinstance(node, ast.Import):
        return any(alias.name.split(".")[0] == "random" for alias in node.names)
    if isinstance(node, ast.ImportFrom):
        return node.module in {"random", "numpy.random"}
    return False


def _random_generator_calls(tree: ast.AST) -> list[int]:
    """Line numbers where the code reaches for a random generator."""
    return [
        node.lineno
        for node in ast.walk(tree)
        if isinstance(node, ast.Attribute | ast.Import | ast.ImportFrom)
        and _is_random_generator(node)
    ]


def test_only_the_gated_data_story_plots_reach_for_a_random_generator() -> None:
    """A builder may draw random points only if the build gates it and stamps it (#775).

    Every function in the package that touches a random generator must be one of
    ``RANDOM_MARK_PLOTS`` and must save with ``random_marks=True``; every other builder
    must not touch one. A new random draw anywhere else fails here."""
    gated = {plot.__name__ for plot in story_plots.RANDOM_MARK_PLOTS}
    assert gated <= {plot.__name__ for plot in story_plots.PLOT_BUILDERS}

    using_random: set[str] = set()
    stamped: set[str] = set()
    for module in sorted(PACKAGE.glob("*.py")):
        tree = ast.parse(module.read_text(encoding="utf-8"))
        functions = [node for node in tree.body if isinstance(node, ast.FunctionDef)]
        inside = {line for function in functions for line in _random_generator_calls(function)}
        assert set(_random_generator_calls(tree)) <= inside, f"{module.name}: module-level draw"
        for function in functions:
            if _random_generator_calls(function):
                assert module.name == "story_plots.py", f"{module.name}:{function.name}"
                using_random.add(function.name)
            for node in ast.walk(function):
                if (
                    isinstance(node, ast.Call)
                    and isinstance(node.func, ast.Name)
                    and node.func.id == "save_and_verify_story_plot"
                ):
                    flags = {keyword.arg: keyword.value for keyword in node.keywords}
                    flag = flags["random_marks"]
                    assert isinstance(flag, ast.Constant)
                    if flag.value is True:
                        stamped.add(function.name)

    assert using_random == gated
    assert stamped == gated


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

    assert build.main(["--output-dir", str(tmp_path), "--audit-random-marks"]) == 1

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


def test_recipe_records_this_checkouts_commit_and_forwards_its_arguments(tmp_path: Path) -> None:
    capture = tmp_path / "arguments.txt"
    shim = tmp_path / "uv"
    shim.write_text('#!/bin/sh\nprintf "%s\\n" "$@" > "$NCASM_PLATES_CAPTURE"\n')
    shim.chmod(0o700)
    justfile = Path(__file__).resolve().parents[1] / "Justfile"
    environment = {
        **os.environ,
        "PATH": str(tmp_path) + os.pathsep + os.environ["PATH"],
        "NCASM_PLATES_CAPTURE": str(capture),
    }

    result = subprocess.run(
        [
            "just",
            "--justfile",
            str(justfile),
            "ncasm-tufte-plates",
            "--audit-random-marks",
            "--output-dir",
            "some folder",
        ],
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    arguments = capture.read_text().splitlines()
    assert arguments[:5] == ["run", "python", "-m", "mtp.ncasm_tufte_plates", "--code-commit"]
    assert len(arguments[5].removesuffix("+dirty")) == 40
    assert arguments[6:] == ["--audit-random-marks", "--output-dir", "some folder"]

    refused = subprocess.run(
        ["just", "--justfile", str(justfile), "ncasm-tufte-plates", "--code-commit=abc"],
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )
    assert refused.returncode != 0
    assert "drop --code-commit=abc" in refused.stderr
