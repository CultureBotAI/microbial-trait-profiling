"""Keep plotting defaults independent across investigation scripts (#723)."""

from __future__ import annotations

import importlib.util
import sys
from io import BytesIO
from pathlib import Path
from typing import Literal

import matplotlib as mpl
import pytest
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

MODULE_PATH = (
    Path(__file__).resolve().parents[1] / "scripts" / "regenerate_trait_profile_investigation.py"
)
SPEC = importlib.util.spec_from_file_location("trait_profile_investigation", MODULE_PATH)
assert SPEC is not None
assert SPEC.loader is not None
investigation = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(investigation)
PlotSetting = Literal["figure.dpi", "savefig.dpi", "font.size"]


def test_sequential_scripts_start_with_independent_plotting_defaults(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A figure after a settings-changing script matches an isolated default figure."""
    observed: list[tuple[object, object, object, tuple[int, int]]] = []
    calls: list[tuple[str, str, Path, list[str]]] = []

    def run_path(path_name: str, *, run_name: str) -> dict[str, object]:
        calls.append((path_name, run_name, Path.cwd(), sys.argv[:]))
        if Path(path_name).name == "change_plot_settings.py":
            mpl.rcParams.update({"figure.dpi": 160, "savefig.dpi": 200, "font.size": 19})
        else:
            figure = Figure(figsize=(2, 1))
            canvas = FigureCanvasAgg(figure)
            figure.savefig(BytesIO(), format="png")
            observed.append(
                (
                    figure.dpi,
                    mpl.rcParams["savefig.dpi"],
                    mpl.rcParams["font.size"],
                    canvas.get_width_height(),
                )
            )
        return {}

    monkeypatch.setattr(investigation.runpy, "run_path", run_path)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", ["caller", "--existing-argument"])
    with mpl.rc_context():
        mpl.rcdefaults()
        default_figure = Figure(figsize=(2, 1))
        default_canvas = FigureCanvasAgg(default_figure)
        expected = (
            default_figure.dpi,
            mpl.rcParams["savefig.dpi"],
            mpl.rcParams["font.size"],
            default_canvas.get_width_height(),
        )
        caller_settings: dict[PlotSetting, int] = {
            "figure.dpi": 143,
            "savefig.dpi": 211,
            "font.size": 14,
        }
        for key, value in caller_settings.items():
            mpl.rcParams[key] = value

        investigation.run_script("change_plot_settings.py")
        investigation.run_script("probe_plot_settings.py", "--example", "value")

        assert observed == [expected]
        assert {key: mpl.rcParams[key] for key in caller_settings} == caller_settings

    assert Path.cwd() == tmp_path
    assert sys.argv == ["caller", "--existing-argument"]
    assert calls == [
        (
            str(MODULE_PATH.parent / "change_plot_settings.py"),
            "__main__",
            MODULE_PATH.parent.parent,
            [str(MODULE_PATH.parent / "change_plot_settings.py")],
        ),
        (
            str(MODULE_PATH.parent / "probe_plot_settings.py"),
            "__main__",
            MODULE_PATH.parent.parent,
            [str(MODULE_PATH.parent / "probe_plot_settings.py"), "--example", "value"],
        ),
    ]


@pytest.mark.parametrize("error_type", [ValueError, SystemExit])
def test_failed_script_restores_plotting_and_process_state(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, error_type: type[BaseException]
) -> None:
    """A failed script cannot contaminate a subsequent caller or investigation."""

    def run_path(path_name: str, *, run_name: str) -> dict[str, object]:
        assert Path(path_name).name == "failing_plot.py"
        assert run_name == "__main__"
        mpl.rcParams.update({"figure.dpi": 160, "savefig.dpi": 200, "font.size": 19})
        raise error_type("script failed")

    monkeypatch.setattr(investigation.runpy, "run_path", run_path)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", ["caller"])
    with mpl.rc_context():
        caller_settings: dict[PlotSetting, int] = {
            "figure.dpi": 143,
            "savefig.dpi": 211,
            "font.size": 14,
        }
        for key, value in caller_settings.items():
            mpl.rcParams[key] = value

        with pytest.raises(error_type, match="script failed"):
            investigation.run_script("failing_plot.py")

        assert {key: mpl.rcParams[key] for key in caller_settings} == caller_settings
    assert Path.cwd() == tmp_path
    assert sys.argv == ["caller"]
