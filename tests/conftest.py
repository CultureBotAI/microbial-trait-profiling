"""Suite-wide guards."""

import importlib.util
import json
import os
import sys
from collections.abc import Callable, Iterator
from pathlib import Path

import pytest

from bsdbng.paths import REPO_ROOT

# The cellulose_decr figure tests import matplotlib, from the viz extra. Without it they
# are left out rather than failing at collection (#631).
collect_ignore = (
    []
    if importlib.util.find_spec("matplotlib")
    else [
        "test_cellulose_decr_figures.py",
        "test_cellulose_nominal_labels.py",
        "test_cellulose_review_layout.py",
        "test_cellulose_review_wording.py",
        "test_cellulose_slide_figures.py",
        "test_cellulose_statement_consumers.py",
    ]
)

# Every real pipeline input and output: raw downloads and their provenance,
# study YAML, derived tables, and the KG-Microbe snapshot.
GUARDED_DATA_DIR = REPO_ROOT / "data"
_PROTECTED_OUTPUT_DIRS = (GUARDED_DATA_DIR, REPO_ROOT / "trait-visualization")


def _snapshot(root: Path) -> dict[str, tuple[int, int]]:
    if not root.is_dir():
        return {}
    return {
        str(path.relative_to(root)): (path.stat().st_mtime_ns, path.stat().st_size)
        for path in root.rglob("*")
        if path.is_file()
    }


@pytest.fixture(autouse=True, scope="session")
def _tests_refuse_protected_output_writes() -> Iterator[None]:
    """Reject audited in-process writes, regardless of how an output path was stored.

    Defaults, composite constants and dynamically constructed paths all reach these
    filesystem operations. Child processes, native code and descriptor-relative writes
    can bypass Python's audit information; the independent snapshot remains a backstop.
    Hooks cannot be removed, so deactivate this one even when session teardown fails.
    """
    # Protect each root entry as well as its storage when the root itself is a link.
    roots = tuple(
        spelling
        for root in _PROTECTED_OUTPUT_DIRS
        for spelling in (Path(os.path.abspath(root)), root.resolve())
    )
    active = True
    write_flags = os.O_WRONLY | os.O_RDWR | os.O_CREAT | os.O_TRUNC | os.O_APPEND
    mutations = {
        "os.mkdir",
        "os.remove",
        "os.rmdir",
        "os.rename",
        "os.link",
        "os.symlink",
        "os.chmod",
        "os.chown",
        "os.truncate",
        "os.utime",
    }

    def refuse_write(event: str, arguments: tuple[object, ...]) -> None:
        if not active:
            return
        if event == "open":
            flags = arguments[2]
            if not isinstance(flags, int) or not flags & write_flags:
                return
            paths = arguments[:1]
        elif event in mutations:
            # Renaming removes the source; links only modify the destination.
            paths = arguments[:2] if event == "os.rename" else arguments[:1]
            if event in {"os.link", "os.symlink"}:
                paths = arguments[1:2]
        else:
            return
        for value in paths:
            if not isinstance(value, str | bytes | os.PathLike):
                continue
            path = Path(os.fsdecode(value))
            absolute = Path(os.path.abspath(path))
            # Unlink/rename affect the final entry, while writable opens affect its target.
            resolved = (
                path.parent.resolve() / path.name
                if event in {"os.remove", "os.rename"}
                else path.resolve()
            )
            if any(
                candidate.is_relative_to(root)
                for root in roots
                for candidate in (absolute, resolved)
            ):
                raise AssertionError(f"tests must not write protected output: {path} ({event})")

    sys.addaudithook(refuse_write)
    try:
        yield
    finally:
        active = False


@pytest.fixture(autouse=True, scope="session")
def _tests_never_write_real_data(_tests_refuse_protected_output_writes: None) -> Iterator[None]:
    """Detect changes under data/ and trait-visualization/ (#417, #421, #1026).

    A test that falls back to a default path writes over real pipeline data:
    ingest() once wrote its fixture log over data/derived/ingest_log.json, and
    ingest() and download_exports() also default to data/studies and data/raw.
    """
    before = {root: _snapshot(root) for root in _PROTECTED_OUTPUT_DIRS}
    yield
    after = {root: _snapshot(root) for root in _PROTECTED_OUTPUT_DIRS}
    changed = sorted(
        str(root / name)
        for root in _PROTECTED_OUTPUT_DIRS
        for name in before[root].keys() | after[root].keys()
        if before[root].get(name) != after[root].get(name)
    )
    assert not changed, (
        f"tests changed real files under protected output roots: {', '.join(changed)}"
    )


@pytest.fixture
def rendered_title_clearance(monkeypatch: pytest.MonkeyPatch) -> Callable[[Path], float]:
    """Measure title-to-header ink clearance in layout pixels after real rendering.

    Diagnostic colours change no geometry. The render's actual scale converts PNG
    pixels to layout pixels; no production margin constants enter the measurement.
    """
    import numpy as np
    import plotly.graph_objects as go
    from PIL import Image

    render = go.Figure.write_image
    scales: dict[Path, int] = {}

    def diagnostic_render(figure: go.Figure, path: str, *, scale: int) -> None:
        figure.update_layout(title_font_color="#ff00ff", title_subtitle_font_color="#ff00ff")
        if figure.layout.annotations:
            for annotation in figure.layout.annotations:
                annotation.font.color = "#00ffff"
        else:
            figure.update_yaxes(tickfont_color="#00ffff")
        scales[Path(path)] = scale
        render(figure, path, scale=scale)

    monkeypatch.setattr(go.Figure, "write_image", diagnostic_render)

    def measure(path: Path) -> float:
        with Image.open(path) as image:
            pixels = np.asarray(image.convert("RGB"))
        red, green, blue = (pixels[:, :, channel] for channel in range(3))
        title_rows = np.where((red > 200) & (blue > 200) & (green < 80))[0]
        header_rows = np.where((green > 200) & (blue > 200) & (red < 80))[0]
        assert title_rows.size, "No detectable title/subtitle glyphs"
        assert header_rows.size, "No detectable header/row glyphs"
        clearance = (int(header_rows.min()) - int(title_rows.max())) / scales[path]
        path.with_suffix(".geometry.json").write_text(
            json.dumps(
                {
                    "title_bottom_px": int(title_rows.max()),
                    "header_top_px": int(header_rows.min()),
                    "render_scale": scales[path],
                    "clearance_layout_px": clearance,
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        return clearance

    return measure
