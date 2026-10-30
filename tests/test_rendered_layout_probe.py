"""Native subtitle measurement and trustworthy layout-probe outcomes (#958)."""

import importlib.util
import json
import sys
from collections.abc import Callable
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Literal

import plotly.graph_objects as go
import pytest

SCRIPT = Path(__file__).parents[1] / "scripts/adhoc/issue839_probe_guild_layout_mutants.py"
SPEC = importlib.util.spec_from_file_location("layout_probe_958", SCRIPT)
assert SPEC is not None
assert SPEC.loader is not None
probe = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(probe)


def test_native_subtitle_ink_is_included_in_clearance(
    tmp_path: Path, rendered_title_clearance: Callable[[Path], float]
) -> None:
    gaps = []
    for name, subtitle in (("main", ""), ("native", "Native subtitle first line<br>Second line")):
        figure = go.Figure(
            go.Heatmap(
                z=[[index] for index in range(8)],
                x=["Trait"],
                y=[f"Condition {index}" for index in range(8)],
                showscale=False,
            )
        )
        figure.update_layout(
            width=700,
            height=500,
            margin={"l": 130, "r": 40, "t": 170, "b": 80},
            title={
                "text": "Main title",
                "y": 0.97,
                "yref": "container",
                "yanchor": "top",
                "font": {"size": 24},
                "subtitle": {"text": subtitle, "font": {"size": 20, "color": "#222222"}},
            },
            annotations=[
                {
                    "text": "Column header",
                    "x": 0.5,
                    "y": 1.05,
                    "xref": "paper",
                    "yref": "paper",
                    "showarrow": False,
                }
            ],
        )
        path = tmp_path / f"{name}.png"
        figure.write_image(str(path), scale=2)
        gaps.append(rendered_title_clearance(path))
    assert gaps[0] - gaps[1] >= 30, f"Native subtitle ink was omitted from clearance: {gaps}"


def fake_pytest(*, failure: str, captured: list[list[str]]) -> Callable[..., pytest.ExitCode]:
    """Drive the real probe with actual pytest CallInfo exceptions, without rendering."""
    calls = 0

    def run(arguments: list[str], **kwargs: Any) -> pytest.ExitCode:
        nonlocal calls
        captured.append(arguments)
        mutated = calls > 0
        site_header_only = calls == 3
        calls += 1
        nodes = [
            f"tests/test_trait_guilds.py::test_the_plotting_area_starts_below_the_last_subtitle_line[{value}]"
            for value in (False, True)
        ] + [
            f"tests/test_table1_trait_maps.py::test_the_site_names_sit_below_the_last_subtitle_line[{value}]"
            for value in (False, True)
        ]
        for index, node in enumerate(nodes):
            if mutated and index == 3 and failure == "missing":
                continue

            def execute(index: int = index) -> None:
                if mutated and not (site_header_only and index == 0):
                    if index == 3:
                        if failure == "crash":
                            raise RuntimeError("Synthetic renderer crash")
                        if failure == "unrelated_assertion":
                            raise AssertionError("Missing input data")
                        if failure == "survived_case":
                            return
                    raise AssertionError("Title/header clearance below 30 layout px")

            phase: Literal["setup", "call"] = (
                "setup" if mutated and index == 3 and failure == "setup" else "call"
            )
            info = pytest.CallInfo.from_call(execute, phase)
            for plugin in kwargs.get("plugins", []):
                plugin.pytest_runtest_makereport(SimpleNamespace(nodeid=node), info)
                if mutated and index == 3 and failure == "teardown":

                    def teardown() -> None:
                        raise RuntimeError("Synthetic cleanup error")

                    plugin.pytest_runtest_makereport(
                        SimpleNamespace(nodeid=node),
                        pytest.CallInfo.from_call(teardown, "teardown"),
                    )
                if mutated and index == 3 and failure == "collection":
                    plugin.pytest_collectreport(
                        SimpleNamespace(failed=True, longreprtext="Synthetic collection error")
                    )
        return pytest.ExitCode.TESTS_FAILED if mutated else pytest.ExitCode.OK

    return run


@pytest.mark.parametrize(
    "failure",
    [
        "geometry",
        "crash",
        "unrelated_assertion",
        "setup",
        "teardown",
        "missing",
        "collection",
        "survived_case",
    ],
)
def test_probe_distinguishes_assertion_rejection_from_renderer_crash(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, failure: str
) -> None:
    captured: list[list[str]] = []
    monkeypatch.setattr(probe.pytest, "main", fake_pytest(failure=failure, captured=captured))
    monkeypatch.setattr(sys, "argv", [str(SCRIPT), "--output-dir", str(tmp_path)])
    status = probe.main()
    assert status == (0 if failure == "geometry" else 1), (
        "Only the intended call-phase geometry assertion may count as a rejected mutant"
    )
    verdicts = json.loads((tmp_path / "verdicts.json").read_text())
    for row in verdicts:
        assert len(row["tests"]) == 4, "Every parameterized layout case needs its own verdict"
    assert len(captured) == 5


def test_probe_disables_pytest_cache(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    captured: list[list[str]] = []
    monkeypatch.setattr(probe.pytest, "main", fake_pytest(failure="geometry", captured=captured))
    monkeypatch.setattr(sys, "argv", [str(SCRIPT), "--output-dir", str(tmp_path)])
    assert probe.main() == 0
    for arguments in captured:
        assert any(
            arguments[index : index + 2] == ["-p", "no:cacheprovider"]
            for index in range(len(arguments) - 1)
        ), "Deliberate mutant failures must not update the normal pytest cache"
