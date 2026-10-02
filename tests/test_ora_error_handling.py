"""Regression checks for failures previously hidden by broad exception handlers."""

import importlib.util
import sys
from pathlib import Path
from types import ModuleType, SimpleNamespace

import pytest


def _load_script(name: str) -> ModuleType:
    path = Path(__file__).resolve().parents[1] / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


ENRICHMENT = _load_script("trait_enrichment")
DEPLETION = _load_script("trait_depletion")
HEATMAP = _load_script("create_trait_enrichment_heatmap")
CYTOSCAPE = _load_script("generate_cytoscape_diagram")
SANKEY = _load_script("generate_six_trait_sankey")


@pytest.mark.parametrize(
    ("module", "tail", "expected"),
    [(ENRICHMENT, "sf", [1.0, 5 / 6, 1 / 6]), (DEPLETION, "cdf", [1 / 6, 5 / 6, 1.0])],
)
def test_hypergeometric_reference_values(
    module: ModuleType, tail: str, expected: list[float]
) -> None:
    # Drawing two of four taxa, two annotated: counts 0, 1, 2 have mass 1/6, 4/6, 1/6.
    scalar = getattr(module, f"hypergeom_{tail}")
    vector = getattr(module, f"hypergeom_{tail}_many")
    assert [scalar(count, 4, 2, 2) for count in [0, 1, 2]] == pytest.approx(expected)
    assert vector([0, 1, 2], 4, [2, 2, 2], 2) == pytest.approx(expected)
    assert vector([], 4, [], 2) == []


@pytest.mark.parametrize(("module", "tail"), [(ENRICHMENT, "sf"), (DEPLETION, "cdf")])
def test_hypergeometric_rejects_invalid_inputs(module: ModuleType, tail: str) -> None:
    scalar = getattr(module, f"hypergeom_{tail}")
    vector = getattr(module, f"hypergeom_{tail}_many")
    with pytest.raises(ValueError, match="Invalid hypergeometric"):
        scalar(1, 4, 5, 2)
    with pytest.raises(ValueError, match="Invalid hypergeometric"):
        vector([1], 4, [5], 2)
    with pytest.raises(ValueError, match="same length"):
        vector([0, 1], 4, [2], 2)
    with pytest.raises(ValueError, match="same length"):
        vector([], 4, [2], 2)


@pytest.mark.parametrize(("module", "tail"), [(ENRICHMENT, "sf"), (DEPLETION, "cdf")])
def test_hypergeometric_computation_failure_is_visible(
    module: ModuleType, tail: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    def fail(*_args: object) -> float:
        raise FloatingPointError("calculation failed")

    monkeypatch.setattr(module, "_scipy_hypergeom", SimpleNamespace(**{tail: fail}))
    with pytest.raises(FloatingPointError, match="calculation failed"):
        getattr(module, f"hypergeom_{tail}")(1, 4, 2, 2)
    with pytest.raises(FloatingPointError, match="calculation failed"):
        getattr(module, f"hypergeom_{tail}_many")([1], 4, [2], 2)


@pytest.mark.parametrize("module", [CYTOSCAPE, SANKEY])
def test_missing_optional_label_files_keep_canonical_labels(
    module: ModuleType, tmp_path: Path
) -> None:
    labels = module.load_trait_labels(tmp_path / "absent-meta.tsv", tmp_path / "absent-map.tsv")
    assert labels["METPO:1000668"] == "coccus shaped"


@pytest.mark.parametrize("module", [CYTOSCAPE, SANKEY])
@pytest.mark.parametrize("source", ["metadata", "mapping"])
def test_malformed_present_label_files_fail(
    module: ModuleType, source: str, tmp_path: Path
) -> None:
    metadata = tmp_path / "metadata.tsv"
    mapping = tmp_path / "mapping.tsv"
    malformed = metadata if source == "metadata" else mapping
    malformed.write_text("wrong_column\ninvalid\n", encoding="utf-8")
    with pytest.raises(ValueError, match=r"missing column_id|at least two columns"):
        module.load_trait_labels(metadata, mapping)


@pytest.mark.parametrize("module", [CYTOSCAPE, SANKEY])
def test_valid_label_files_preserve_labels(module: ModuleType, tmp_path: Path) -> None:
    metadata = tmp_path / "metadata.tsv"
    mapping = tmp_path / "mapping.tsv"
    metadata.write_text(
        "column_id\tcolumn_label\ncustom:X\tcapable of: sentinel\n", encoding="utf-8"
    )
    mapping.write_text("custom:Y\tother\n", encoding="utf-8")
    labels = module.load_trait_labels(metadata, mapping)
    assert labels["custom:X"] == "sentinel"
    assert labels["custom:Y"] == "other"


def test_heatmap_rejects_malformed_mapping(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    mapping = tmp_path / "mapping.tsv"
    mapping.write_text("unexpected\tcolumns\ncustom:X\tsentinel\n", encoding="utf-8")
    monkeypatch.setattr(HEATMAP, "METPO_MAPPING_PATH", mapping)
    with pytest.raises(ValueError, match="requires METPO_ID and Description"):
        HEATMAP.load_metpo_label_lookup()


def test_heatmap_without_optional_files_retains_raw_trait(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    hit = HEATMAP.ParsedHit("study_1", "exp1", "inc", "custom:X", 0.01, 2.0)
    matrix, metadata = HEATMAP.build_matrix([hit], 0, 0, 0, 0, 0.0)
    assert matrix.columns.tolist() == ["custom:X"]
    assert matrix.loc["study_1|exp1|inc", "custom:X"] == pytest.approx(2.0)
    assert metadata.loc["study_1|exp1|inc", "study_id"] == "study_1"
