"""The NCASM plate builders ported from slide-wizard: outputs, repeatability, hygiene (#2264).

Plates 1-20 and 22-30 read no data files, and none of these tests checks their numbers
against data; #775-#779 record what is wrong with their content. Plate 21 is drawn from
run tables, and the tests at the end check every number it prints against an independent
reading of those tables.
"""

from __future__ import annotations

import ast
import csv
import hashlib
import json
import os
import re
import shutil
import statistics
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
    run_tables,
    story_from_tables,
    story_plots,
)

STORED_TEXT = (
    Path(__file__).parent / "fixtures/ncasm_tufte_plates/stored_plate_text_2026-10-22.json"
)
FIXTURE_RUN = Path(__file__).parent / "fixtures/ncasm_tufte_plates/run_tier12_off_511c891"
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
        manifest = build.build_plates(
            out_dir, audit_random_marks=True, run_root=FIXTURE_RUN, code_commit=CODE_COMMIT
        )
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
        "run": {
            "comparison_id": "comparison-disease_vs_healthy_tier1_2_subject_controls-3d40e23048d0",
            "contrast_population": "disease_vs_healthy_tier1_2_subject_controls",
            "direction_basis": "disease_vs_comparison",
            "core_run_id": "run-a2d30f9ca1dd45aeb98085f39badc4bf",
            "git_commit": "511c891229b659e207cfe0f4281a97c3fc039a80",
            "git_dirty": False,
        },
        "combined_pdf": build.COMBINED_PDF,
        "plates": manifest,
    }


def test_manifest_says_where_each_plates_numbers_and_marks_come_from(built: Built) -> None:
    """Plate 21 reads run tables; the rest do not: #776 for numbers, #775 for random points."""
    _out_dir, manifest = built
    for entry in manifest[:20]:
        assert entry["data_status"] == {
            "numbers": "typed_literals",
            "marks": "typed_literals",
            "inputs": [],
            "issues": [776, 777, 778, 779],
        }
    from_tables: data_status.DataStatus = {
        "numbers": "run_tables",
        "marks": "run_tables",
        "inputs": [
            {
                "path": relative.as_posix(),
                "sha256": hashlib.sha256((FIXTURE_RUN / relative).read_bytes()).hexdigest(),
            }
            for relative in (
                run_tables.COMPARISON_RUN,
                run_tables.ARM1_POOLED_STATS,
                run_tables.DIRECTION_SCREEN,
            )
        ],
        "issues": [],
    }
    assert manifest[20]["data_status"] == from_tables
    for entry in manifest[21:]:
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

    with pytest.raises(build.PlateBuildError, match="--audit-random-marks") as refusal:
        build.build_plates(out_dir, run_root=FIXTURE_RUN)

    assert str(refusal.value).startswith("plates 22, 23, 24, 25, 26, 27, 28, 29, 30: ")
    assert not out_dir.exists()


def test_build_stops_before_drawing_when_the_table_plate_has_no_run(tmp_path: Path) -> None:
    out_dir = tmp_path / "plates"

    with pytest.raises(build.PlateBuildError, match="--run") as refusal:
        build.build_plates(out_dir, audit_random_marks=True)

    assert str(refusal.value).startswith("plate 21: drawn from run tables")
    assert not out_dir.exists()


def test_main_reports_the_refusal_and_exits_nonzero(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    out_dir = tmp_path / "plates"

    assert build.main(["--output-dir", str(out_dir)]) == 1

    said = capsys.readouterr().out
    assert "--audit-random-marks" in said
    assert "--run" in said
    assert not out_dir.exists()

    missing = ["--output-dir", str(out_dir), "--audit-random-marks", "--run", str(tmp_path)]
    assert build.main(missing) == 1
    assert "comparison/comparison_run.json is missing" in capsys.readouterr().out


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
    assert gated == {plot.random.__name__ for plot in build.STORY_PLOTS if plot.random}
    assert [plot.number for plot in build.STORY_PLOTS] == list(range(1, 11))
    assert all((plot.random is None) != (plot.from_run is None) for plot in build.STORY_PLOTS)

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


# --- Plate 21: every number against an independent reading of the run tables (#776) ---


def _table(relative: Path, root: Path = FIXTURE_RUN) -> list[dict[str, str]]:
    with (root / relative).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _signed(value: float) -> str:
    return f"{value:+.2f}".replace("-", "\u2212")


def _expected_plate_21_text(root: Path = FIXTURE_RUN) -> set[str]:
    """Every digit-bearing string plate 21 should print, computed here from the TSVs.

    This reads the tables with ``csv`` and formats the values itself; it shares no code
    with ``run_tables`` or ``story_from_tables``."""
    arm1 = {row["trait"]: row for row in _table(run_tables.ARM1_POOLED_STATS, root)}
    screen = _table(run_tables.DIRECTION_SCREEN, root)
    both = [row for row in screen if row["trait"] in arm1]
    table1 = [row for row in both if row["in_table1"] == "true"]
    by_hits = [row for row in screen if float(row["by_q_nominal"]) <= 0.05]
    table1_hits = [row for row in table1 if float(row["by_q_nominal"]) <= 0.05]
    comparisons = {row["pair_count"] for row in arm1.values()}
    studies = {row["study_count"] for row in arm1.values()}
    assert len(comparisons) == len(studies) == 1
    n_comparisons, n_studies = int(comparisons.pop()), int(studies.pop())
    r = statistics.correlation(
        [float(arm1[row["trait"]]["mean_increased_minus_decreased"]) for row in both],
        [float(row["share_up"]) for row in both],
    )
    commit = json.loads((root / run_tables.COMPARISON_RUN).read_text())["git_commit"]
    words = {6: "Six"}
    labels = set()
    for row in table1:
        name = row["object_name"].replace("_", " ")
        shift = float(arm1[row["trait"]]["mean_increased_minus_decreased"])
        labels.add(
            f"{name[0].upper()}{name[1:]}  {_signed(shift)} \u00b7 "
            f"{100 * float(row['share_up']):.1f}% up "
            f"({row['studies_up']}:{row['studies_down']}) \u00b7 "
            f"BY q = {float(row['by_q_nominal']):.3f}"
        )
    return labels | {
        "DATA STORY PLOT 01 OF 10  \u00b7  TUFTE JOINT SCATTER + MARGINAL HISTOGRAMS "
        "(CF. TUFTE_023, TUFTE_101)",
        f"Where the {words[len(table1)]} Table 1 Traits Sit Among {len(both)} Traits Tested "
        "by Both Arms",
        f"One dot per trait: mean paired shift across {n_comparisons:,} comparisons "
        "(\u0394z) against the share of studies in which the trait rises.",
        f"Source: Arm 1 pooled shift ({n_comparisons:,} comparisons, {n_studies} studies) and "
        f"Arm 3 study sign screen ({len(screen)} traits); run {commit[:7]}.",
        "Tufte Data Story (1/10)",
        f"A. {len(both)} Traits: Mean Paired Shift (\u0394z) vs. Share of Studies Voting Up",
        f"Gray dots = {len(both) - len(table1)} other traits; dark rings = the {len(by_hits)} "
        "traits passing BY q \u2264 0.05 in the sign screen",
        f"Arm 1 Mean Paired Shift \u0394z, Increased \u2212 Decreased "
        f"({n_comparisons:,} Comparisons)",
        "Arm 3 Share of Studies Voting Up (%)",
        f"{len(both)} / {len(screen)}",
        f"{len(screen) - len(both)} screen traits lack an Arm 1 row",
        f"r = {r:.2f}",
        f"{len(by_hits)} / {len(screen)}",
        "Pass BY q \u2264 0.05 (sign test)",
        f"{len(table1_hits)} of the {len(table1)} Table 1 traits pass",
    }


def _is_axis_tick(text: str) -> bool:
    return re.fullmatch(r"[+\u2212]?\d+(\.\d+)?%?", text) is not None


def _digit_bearing_text(svg: Path) -> set[str]:
    return {
        text
        for text in _svg_text(svg)
        if any(character.isdigit() for character in text) and not _is_axis_tick(text)
    }


def test_plate_21_prints_only_numbers_read_from_the_run_tables(built: Built) -> None:
    """Every string with a digit on the plate equals one computed here from the tables.

    Axis tick labels are the only digits left out. A number typed into the builder, or a
    table value formatted from the wrong column, fails here."""
    out_dir, manifest = built
    printed = _digit_bearing_text(out_dir / manifest[20]["svg_path"])

    assert printed == _expected_plate_21_text()
    # One value spelled out, so the comparison above cannot pass on two matching mistakes.
    assert (
        "Cellulose catabolic process  \u22120.25 \u00b7 40.6% up (181:265) \u00b7 BY q = 0.041"
        in printed
    )


def test_plate_21_draws_one_dot_per_trait_at_the_tabulated_values() -> None:
    cloud = run_tables.load_trait_cloud(FIXTURE_RUN)
    arm1 = {row["trait"]: row for row in _table(run_tables.ARM1_POOLED_STATS)}
    screen = {row["trait"]: row for row in _table(run_tables.DIRECTION_SCREEN)}

    assert [point.trait for point in cloud.points] == list(arm1)
    for point in cloud.points:
        assert point.mean_shift == float(arm1[point.trait]["mean_increased_minus_decreased"])
        assert point.share_up == float(screen[point.trait]["share_up"])
        assert point.studies_up == int(screen[point.trait]["studies_up"])
        assert point.studies_down == int(screen[point.trait]["studies_down"])
        assert point.in_table1 == (screen[point.trait]["in_table1"] == "true")
    assert cloud.comparisons == 935
    assert cloud.arm1_studies == 548
    assert (len(cloud.points), cloud.screen_traits, len(cloud.table1)) == (311, 317, 6)


def _edited_run(tmp_path: Path, relative: Path, edit: dict[str, dict[str, str]]) -> Path:
    """A copy of the fixture run with some cells of one table replaced, keyed by trait."""
    root = tmp_path / "run"
    shutil.copytree(FIXTURE_RUN, root)
    rows = _table(relative, root)
    for row in rows:
        row.update(edit.get(row["trait"], {}))
    with (root / relative).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return root


def test_plate_21_follows_the_tables_when_they_change(tmp_path: Path) -> None:
    """Change one table cell and the plate, its note and its input hash change with it."""
    cellulose = "METPO:2000103|GO:0030245"
    root = _edited_run(
        tmp_path,
        run_tables.DIRECTION_SCREEN,
        {cellulose: {"studies_up": "200", "share_up": "0.4301", "by_q_nominal": "0.2"}},
    )
    with pytest.MonkeyPatch.context() as patch, mpl.style.context("default"):
        patch.setattr(story_plots, "SLIDE_DPI", TEST_DPI)
        drawn = story_from_tables.plot_01_trait_cloud(tmp_path / "out", root)

    printed = _digit_bearing_text(Path(drawn.info["svg_path"]))
    assert printed == _expected_plate_21_text(root)
    assert (
        "Cellulose catabolic process  \u22120.25 \u00b7 43.0% up (200:265) \u00b7 BY q = 0.200"
        in printed
    )
    assert "3 / 317" in printed
    assert "0 of the 6 Table 1 traits pass" in printed
    assert "Of the 6 Table 1 traits, 0 pass (none)" in drawn.speaker_notes
    original = hashlib.sha256((FIXTURE_RUN / run_tables.DIRECTION_SCREEN).read_bytes())
    hashes = {entry["path"]: entry["sha256"] for entry in drawn.inputs}
    assert hashes[run_tables.DIRECTION_SCREEN.as_posix()] != original.hexdigest()


@pytest.mark.parametrize(
    ("relative", "edit", "message"),
    [
        (
            run_tables.ARM1_POOLED_STATS,
            {"METPO:1000608": {"stratum_kind": "body_site"}},
            "stratum_kind is not constant",
        ),
        (
            run_tables.ARM1_POOLED_STATS,
            {"METPO:1000608": {"mean_increased_minus_decreased": "nan"}},
            "is not finite",
        ),
        (
            run_tables.ARM1_POOLED_STATS,
            {"METPO:1000608": {"trait": "METPO:0000000"}},
            "lacks pooled-shift traits",
        ),
        (
            run_tables.DIRECTION_SCREEN,
            {"METPO:1000608": {"test_status": "untestable"}},
            "test_status is not constant",
        ),
        (
            run_tables.DIRECTION_SCREEN,
            {"METPO:1000608": {"studies_up": "12.5"}},
            "is not a count",
        ),
    ],
)
def test_loader_refuses_tables_that_are_not_what_the_plate_assumes(
    tmp_path: Path, relative: Path, edit: dict[str, dict[str, str]], message: str
) -> None:
    root = _edited_run(tmp_path, relative, edit)

    with pytest.raises(run_tables.RunTableError, match=message):
        run_tables.load_trait_cloud(root)


def test_loader_refuses_an_unfinished_comparison_and_a_missing_table(tmp_path: Path) -> None:
    root = tmp_path / "run"
    shutil.copytree(FIXTURE_RUN, root)
    record = json.loads((root / run_tables.COMPARISON_RUN).read_text())
    (root / run_tables.COMPARISON_RUN).write_text(json.dumps(record | {"status": "failed"}))
    with pytest.raises(run_tables.RunTableError, match="not completed"):
        run_tables.load_trait_cloud(root)

    (root / run_tables.DIRECTION_SCREEN).unlink()
    with pytest.raises(run_tables.RunTableError, match=r"trait_direction_screen\.tsv is missing"):
        run_tables.load_trait_cloud(root)
