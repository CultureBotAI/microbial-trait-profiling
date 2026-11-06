"""Figure 2's ORA streaming and pre-change figure contract (#1505)."""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterator
from pathlib import Path
from types import TracebackType
from typing import Any, TextIO

import pandas as pd
import plotly.graph_objects as go
import pytest

from mtp.deck_visualizations import generate_figure2_calibrated_ora_infographic

FIXTURES = Path(__file__).parent / "fixtures" / "ora_figure_1505"


def _inputs(tmp_path: Path, combined: bool) -> dict[str, Path]:
    contents = {
        "matrix": "taxon_id\tMETPO:1000608\tMETPO:1000668\nNCBI:1\t1\t0\nNCBI:2\t0\t1\n",
        "experiments": "body_site\nFeces\nFeces\nBlood\n",
        "enrichment": (
            "# FDR: True\n# Stat-Type: fdr_q\n\n"
            "trait\tinc_dec\tfdr_q\n"
            "coccus shaped\tinc\t0.01\n"
            "coccus shaped\tinc\t0.02\n"
            "facultatively aerobic\tinc\t0.04\n"
            "not significant\tinc\t0.051\n"
            + ("cellulose catabolism\tdec\t0.03\n" if combined else "")
        ),
        "depletion": (
            "# Stat-Type: fdr_q\n# FDR: True\n"
            "trait\tfdr_q\ncellulose catabolism\t0.03\nnot significant\t0.8\n"
        ),
    }
    paths = {}
    for name, content in contents.items():
        path = tmp_path / f"{name}.tsv"
        path.write_text(content, encoding="utf-8")
        paths[name] = path
    return paths


def _capture_figure(monkeypatch: pytest.MonkeyPatch, paths: dict[str, Path], output: Path) -> bytes:
    captured: list[bytes] = []

    def capture(figure: go.Figure, filename: str, **kwargs: Any) -> None:
        assert filename == str(output)
        assert kwargs == {"scale": 2}
        captured.append(figure.to_json().encode("utf-8"))

    monkeypatch.setattr(go.Figure, "write_image", capture)
    assert (
        generate_figure2_calibrated_ora_infographic(
            paths["matrix"],
            paths["experiments"],
            output,
            paths["enrichment"],
            paths["depletion"],
        )
        == output
    )
    assert len(captured) == 1
    return captured[0]


@pytest.mark.parametrize("combined", [False, True], ids=["separate", "combined"])
def test_figure2_plotly_bytes_match_prechange_snapshot(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, combined: bool
) -> None:
    """Pin every Plotly figure byte before its expensive renderer is invoked."""
    paths = _inputs(tmp_path, combined)
    actual = _capture_figure(monkeypatch, paths, tmp_path / "figure.png")
    name = "combined.json" if combined else "separate.json"
    actual_path = tmp_path / name
    actual_path.write_bytes(actual)
    expected_path = FIXTURES / name
    assert expected_path.exists(), f"Pre-change capture is ready at {actual_path}"
    assert actual == expected_path.read_bytes().removesuffix(b"\n")


def test_figure2_does_not_parse_unused_depletion_body(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """An unused table's adjusted header matters; its body is never parsed."""
    paths = _inputs(tmp_path, combined=True)
    paths["depletion"].write_text(
        '# Stat-Type: fdr_q\n# FDR: True\ntrait\tfdr_q\n"unterminated\t0.01\n',
        encoding="utf-8",
    )
    actual = _capture_figure(monkeypatch, paths, tmp_path / "figure.png")
    assert actual == (FIXTURES / "combined.json").read_bytes().removesuffix(b"\n")


class _ReadOnceFile:
    """Reject rereads and unbounded full-text allocations at the source handle."""

    def __init__(self, handle: TextIO) -> None:
        self.handle = handle
        self.characters_read = 0

    def __enter__(self) -> _ReadOnceFile:
        return self

    def __exit__(
        self,
        _exc_type: type[BaseException] | None,
        _exc_value: BaseException | None,
        _traceback: TracebackType | None,
    ) -> None:
        self.close()

    @property
    def mode(self) -> str:
        return self.handle.mode

    @property
    def closed(self) -> bool:
        return self.handle.closed

    def __iter__(self) -> Iterator[str]:
        return self

    def __next__(self) -> str:
        line = self.readline()
        if not line:
            raise StopIteration
        return line

    def readable(self) -> bool:
        return True

    def read(self, size: int = -1) -> str:
        assert size >= 0, "ORA input must not be allocated as one full text string"
        value = self.handle.read(size)
        self.characters_read += len(value)
        return value

    def readline(self, size: int = -1) -> str:
        value = self.handle.readline(size)
        self.characters_read += len(value)
        return value

    def seek(self, offset: int, whence: int = 0) -> int:
        raise AssertionError(f"ORA input must not be reread (seek {offset}, {whence})")

    def close(self) -> None:
        self.handle.close()


@pytest.mark.parametrize("combined", [False, True], ids=["separate", "combined"])
@pytest.mark.parametrize("shared_source", [None, "enrichment", "depletion"])
def test_figure2_streams_each_ora_source_once(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    combined: bool,
    shared_source: str | None,
) -> None:
    paths = _inputs(tmp_path, combined)
    if shared_source is not None:
        paths["enrichment"] = paths["depletion"] = paths[shared_source]
    sources = {paths["enrichment"], paths["depletion"]}
    lengths = {path: len(path.read_text(encoding="utf-8")) for path in sources}
    opens: Counter[Path] = Counter()
    parse_calls: Counter[Path] = Counter()
    readers: dict[Path, _ReadOnceFile] = {}
    original_open = Path.open
    original_read_csv = pd.read_csv

    def open_source(path: Path, *args: Any, **kwargs: Any) -> Any:
        handle = original_open(path, *args, **kwargs)
        if path in sources:
            opens[path] += 1
            assert opens[path] == 1, "Validated ORA source must not be reopened"
            reader = _ReadOnceFile(handle)
            readers[path] = reader
            return reader
        return handle

    def read_csv(source: Any, *args: Any, **kwargs: Any) -> Any:
        for path, reader in readers.items():
            if getattr(source, "_source", None) is reader:
                parse_calls[path] += 1
        return original_read_csv(source, *args, **kwargs)

    monkeypatch.setattr(Path, "open", open_source)
    monkeypatch.setattr(pd, "read_csv", read_csv)
    _capture_figure(monkeypatch, paths, tmp_path / "figure.png")
    assert opens == Counter(dict.fromkeys(sources, 1))
    for path, reader in readers.items():
        unused = combined and shared_source is None and path == paths["depletion"]
        if unused:
            assert reader.characters_read == len("# Stat-Type: fdr_q\n# FDR: True\ntrait\tfdr_q\n")
            assert parse_calls[path] == 0
        else:
            assert reader.characters_read == lengths[path]
            assert parse_calls[path] == 1
        assert reader.closed


@pytest.mark.parametrize("nominal_source", ["enrichment", "depletion"])
def test_figure2_refuses_nominal_even_when_depletion_would_be_unused(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, nominal_source: str
) -> None:
    paths = _inputs(tmp_path, combined=True)
    source = paths[nominal_source]
    source.write_text(
        source.read_text(encoding="utf-8").replace("# FDR: True", "# FDR: False"),
        encoding="utf-8",
    )

    def unexpected_render(*_args: Any, **_kwargs: Any) -> None:
        pytest.fail("Nominal ORA must fail before rendering")

    monkeypatch.setattr(go.Figure, "write_image", unexpected_render)
    with pytest.raises(ValueError, match="Nominal ORA"):
        generate_figure2_calibrated_ora_infographic(
            paths["matrix"],
            paths["experiments"],
            tmp_path / "figure.png",
            paths["enrichment"],
            paths["depletion"],
        )
