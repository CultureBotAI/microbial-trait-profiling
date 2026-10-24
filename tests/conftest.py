"""Suite-wide guards."""

import importlib.util
import json
from collections.abc import Callable, Iterator
from pathlib import Path

import pytest

from bsdbng.paths import REPO_ROOT

# The cellulose_decr figure tests import matplotlib, from the viz extra. Without it they
# are left out rather than failing at collection (#631).
collect_ignore = (
    [] if importlib.util.find_spec("matplotlib") else ["test_cellulose_decr_figures.py"]
)

# Every real pipeline input and output: raw downloads and their provenance,
# study YAML, derived tables, and the KG-Microbe snapshot.
GUARDED_DATA_DIR = REPO_ROOT / "data"


def _snapshot(root: Path) -> dict[str, tuple[int, int]]:
    if not root.is_dir():
        return {}
    return {
        str(path.relative_to(root)): (path.stat().st_mtime_ns, path.stat().st_size)
        for path in root.rglob("*")
        if path.is_file()
    }


@pytest.fixture(autouse=True, scope="session")
def _tests_never_write_real_data() -> Iterator[None]:
    """Fail the run if any test creates, modifies or deletes a file under data/ (#417, #421).

    A test that falls back to a default path writes over real pipeline data:
    ingest() once wrote its fixture log over data/derived/ingest_log.json, and
    ingest() and download_exports() also default to data/studies and data/raw.
    """
    before = _snapshot(GUARDED_DATA_DIR)
    yield
    after = _snapshot(GUARDED_DATA_DIR)
    changed = sorted(
        name for name in before.keys() | after.keys() if before.get(name) != after.get(name)
    )
    assert not changed, f"tests changed real files under {GUARDED_DATA_DIR}: {', '.join(changed)}"


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
        figure.update_layout(title_font_color="#ff00ff")
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
