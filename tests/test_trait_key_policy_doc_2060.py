"""Pin the five-positional-argument recipe example and mutation phase checks (#2060)."""

import importlib.util
import re
import shlex
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
COMMAND = (
    "just table1-taxon-influence-from <kg> data/derived "
    "data/derived/table1_taxon_influence phenotype raw_object_v1"
)


def documented_command() -> str:
    text = (ROOT / "methods/trait_key_policy.md").read_text()
    paragraph = text.split("Standalone contribution,", 1)[1].split("\n\n", 1)[0]
    assert "`trait_key_policy` as its fifth positional parameter" in " ".join(paragraph.split())
    assert "pass all five arguments when overriding it" in " ".join(paragraph.split())
    examples: list[str] = [
        match.group(1)
        for match in re.finditer(r"`(just table1-taxon-influence-from [^`]+)`", paragraph)
    ]
    assert examples == [COMMAND]
    return examples[0]


def test_fifth_parameter_documentation() -> None:
    assert documented_command() == COMMAND


def test_example_matches_recipe_parameters() -> None:
    headers = [
        line
        for line in (ROOT / "Justfile").read_text().splitlines()
        if line.startswith("table1-taxon-influence-from ")
    ]
    assert len(headers) == 1
    assert headers[0].endswith(":")
    parameters = shlex.split(headers[0][:-1])[1:]
    assert parameters == [
        "kg_dir",
        "inputs=data/derived",
        "output=data/derived/table1_taxon_influence",
        "strain_rollup=phenotype",
        "trait_key_policy=facultative_union_v1",
    ]
    args = shlex.split(documented_command())[2:]
    assert len(args) == 5
    assert args[0] == "<kg>"
    assert args[1:4] == [parameter.split("=", 1)[1] for parameter in parameters[1:4]]
    assert args[4] == "raw_object_v1"


def test_verifier_rejects_setup_error(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    path = ROOT / "scripts/adhoc/verify_2060_mutants_codex2_20261111.py"
    spec = importlib.util.spec_from_file_location("verifier_2060_setup_test", path)
    assert spec is not None
    assert spec.loader is not None
    verifier = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = verifier
    try:
        spec.loader.exec_module(verifier)
        (tmp_path / "tests").mkdir()
        test_rel = "tests/setup_probe_2060.py"
        (tmp_path / test_rel).write_text(
            "import pytest\n"
            "@pytest.fixture(autouse=True)\n"
            "def broken_setup():\n"
            "    raise ValueError('AssertionError in setup')\n"
            "def test_probe():\n"
            "    pass\n"
        )
        monkeypatch.setattr(verifier, "TEST", test_rel)
        with monkeypatch.context() as context:
            context.chdir(tmp_path)
            context.setattr(sys, "path", list(sys.path))
            context.setenv("MPLBACKEND", "Agg")
            context.setenv("PYTHONPATH", str(tmp_path / "src"))
            verifier.worker(str(tmp_path), "test_probe")
        with pytest.raises(RuntimeError, match="Setup/collection/teardown error or skip"):
            verifier.validate_result(tmp_path, "test_probe", mutant=True)
    finally:
        sys.modules.pop(spec.name, None)
