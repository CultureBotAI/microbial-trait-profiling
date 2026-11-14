"""Legacy consumers share canonical and explicitly historical facultative labels (#1234)."""

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import pandas as pd
import pytest

from mtp.trait_key_policy import FACULTATIVE_RAW_LABEL, FACULTATIVE_UNION_LABEL

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = (
    "create_six_shared_trait_sankeys",
    "create_pd_cellulose_analysis",
    "create_pd_ecological_reporter_combined",
    "create_pd_significant_condition_taxon_trait_network",
    "investigate_facultative_aerobe_expansion",
)


def load_script(name: str, root: Path = ROOT) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, root / "scripts" / f"{name}.py")
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    sys.path.insert(0, str(root / "scripts"))
    try:
        spec.loader.exec_module(module)
    finally:
        sys.path.pop(0)
    return module


def check_module(module: ModuleType, tmp_path: Path) -> None:
    resolver = getattr(module, "resolve_facultative_display_label", None)
    assert callable(resolver)
    assert resolver() == FACULTATIVE_UNION_LABEL
    assert resolver("phenotype") == FACULTATIVE_UNION_LABEL
    assert resolver("off") == FACULTATIVE_RAW_LABEL
    with pytest.raises(ValueError, match="policy"):
        resolver("unknown")
    source = tmp_path / "labels.tsv"
    frame = pd.DataFrame(
        {
            "trait_id": ["METPO:1000608"] * 2 + ["other"],
            "trait_label": [FACULTATIVE_RAW_LABEL, FACULTATIVE_UNION_LABEL, "aerobic"],
            "value": [1.0, 2.0, 3.0],
        }
    )
    frame.to_csv(source, sep="\t", index=False)
    loaded = module.read_tsv(source)
    assert loaded["trait_label"].tolist() == [FACULTATIVE_UNION_LABEL] * 2 + ["aerobic"]
    assert loaded["trait_id"].tolist() == frame["trait_id"].tolist()
    assert loaded["value"].tolist() == frame["value"].tolist()
    historical = module.read_tsv(source, policy="off")
    assert historical["trait_label"].tolist() == [FACULTATIVE_RAW_LABEL] * 2 + ["aerobic"]
    if hasattr(module, "TRAIT_COLORS"):
        assert (
            module.TRAIT_COLORS[FACULTATIVE_RAW_LABEL]
            == module.TRAIT_COLORS[FACULTATIVE_UNION_LABEL]
        )
    if hasattr(module, "SHARED_TRAITS"):
        trait = next(item for item in module.SHARED_TRAITS if item.trait_id == "METPO:1000608")
        assert trait.label == FACULTATIVE_UNION_LABEL
    if hasattr(module, "SIX_TRAITS"):
        assert FACULTATIVE_UNION_LABEL in module.SIX_TRAITS
        assert FACULTATIVE_RAW_LABEL not in module.SIX_TRAITS
        assert FACULTATIVE_UNION_LABEL in module.INCREASED_SIDE_TRAITS


@pytest.mark.parametrize("name", SCRIPTS)
def test_legacy_facultative_policy(name: str, tmp_path: Path) -> None:
    check_module(load_script(name), tmp_path)
