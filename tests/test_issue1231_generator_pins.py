"""Tests that issue1231_generate_historical_union_core_fixture loads unsplit pins (#1258)."""

from __future__ import annotations

import ast
import importlib.util
import json
import re
import tarfile
from pathlib import Path
from types import ModuleType

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
GENERATOR_SCRIPT = (
    REPO_ROOT / "scripts" / "adhoc" / "issue1231_generate_historical_union_core_fixture.py"
)
UNION_FIXTURE_ARCHIVE = (
    REPO_ROOT / "tests" / "fixtures" / "trait_key_policy" / "historical_union_core.tar.gz"
)


def _load_generator_module() -> ModuleType:
    spec = importlib.util.spec_from_file_location("issue1231_generator", GENERATOR_SCRIPT)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _find_chunk_joined_hex_literals(source_text: str) -> list[str]:
    """Return any AST Constant/Call nodes that join short hex chunks to form a digest."""
    tree = ast.parse(source_text)
    offenders: list[str] = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        if not (
            isinstance(node.func, ast.Attribute)
            and node.func.attr == "join"
            and isinstance(node.func.value, ast.Constant)
            and node.func.value.value == ""
            and len(node.args) == 1
            and isinstance(node.args[0], (ast.Tuple, ast.List))
        ):
            continue
        elements = node.args[0].elts
        if elements and all(
            isinstance(elt, ast.Constant)
            and isinstance(elt.value, str)
            and bool(re.fullmatch(r"[0-9a-fA-F]{2,16}", elt.value))
            for elt in elements
        ):
            joined = "".join(
                elt.value
                for elt in elements
                if isinstance(elt, ast.Constant) and isinstance(elt.value, str)
            )
            offenders.append(joined)
    return offenders


def test_generator_does_not_split_hex_pins_into_chunks() -> None:
    """#1258: generator script must not use ''.join(...) over hex chunks to evade detect-secrets."""
    source_text = GENERATOR_SCRIPT.read_text(encoding="utf-8")
    offenders = _find_chunk_joined_hex_literals(source_text)
    assert offenders == [], (
        f"Found chunk-split hex literals in {GENERATOR_SCRIPT.name}: {offenders}"
    )

    module = _load_generator_module()
    assert getattr(module, "EXPECTED_MAIN_REFERENCE", None) == "44897c6"


def test_load_audited_producer_pins_reads_and_validates_committed_fixture() -> None:
    """#1258: load_audited_producer_pins derives commit, sha256, and size from provenance."""
    module = _load_generator_module()
    assert hasattr(module, "load_audited_producer_pins"), (
        "Generator must expose load_audited_producer_pins(fixture_path)"
    )

    commit, matrix_sha256, matrix_size = module.load_audited_producer_pins(UNION_FIXTURE_ARCHIVE)
    assert commit == module.EXPECTED_MAIN_COMMIT
    assert matrix_sha256 == module.EXPECTED_MAIN_MATRIX_SHA256
    assert matrix_size == module.EXPECTED_MAIN_MATRIX_SIZE

    assert re.fullmatch(r"[0-9a-f]{40}", commit)
    assert commit.startswith("44897c6")
    assert re.fullmatch(r"[0-9a-f]{64}", matrix_sha256)
    assert matrix_size == 28555

    with tarfile.open(UNION_FIXTURE_ARCHIVE, mode="r:gz") as tar:
        gen_member = tar.extractfile("generation.json")
        run_member = tar.extractfile("core/run.json")
        assert gen_member is not None
        assert run_member is not None
        generation = json.loads(gen_member.read().decode("utf-8"))
        run_record = json.loads(run_member.read().decode("utf-8"))

    assert commit == generation["producer_commit"] == run_record["git_commit"]
    assert matrix_sha256 == generation["producer_matrix_inventory"]["sha256"]
    assert matrix_size == generation["producer_matrix_inventory"]["size_bytes"]


def test_load_audited_producer_pins_rejects_inconsistent_or_malformed_provenance(
    tmp_path: Path,
) -> None:
    """#1258: load_audited_producer_pins rejects tampered commit, sha256, path, or run.json."""
    module = _load_generator_module()
    bundle_dir = tmp_path / "bundle"
    module.extract_fixture_bundle(UNION_FIXTURE_ARCHIVE, bundle_dir)

    # 1. Directory bundle works identically when unmodified
    assert module.load_audited_producer_pins(bundle_dir) == (
        module.EXPECTED_MAIN_COMMIT,
        module.EXPECTED_MAIN_MATRIX_SHA256,
        module.EXPECTED_MAIN_MATRIX_SIZE,
    )

    gen_path = bundle_dir / "generation.json"
    run_path = bundle_dir / "core" / "run.json"
    orig_gen = json.loads(gen_path.read_text(encoding="utf-8"))
    orig_run = json.loads(run_path.read_text(encoding="utf-8"))

    # 2. Wrong commit prefix (even when core/run.json has the same 40-hex commit)
    bad_prefix_gen = dict(orig_gen, producer_commit="0" * 40)
    bad_prefix_run = dict(orig_run, git_commit="0" * 40)
    gen_path.write_text(json.dumps(bad_prefix_gen), encoding="utf-8")
    run_path.write_text(json.dumps(bad_prefix_run), encoding="utf-8")
    with pytest.raises(ValueError, match="starting with 44897c6"):
        module.load_audited_producer_pins(bundle_dir)

    # 3. Commit mismatch between generation.json and core/run.json
    gen_path.write_text(json.dumps(orig_gen), encoding="utf-8")
    bad_run_commit = dict(orig_run, git_commit="44897c6" + "0" * 33)
    run_path.write_text(json.dumps(bad_run_commit), encoding="utf-8")
    with pytest.raises(ValueError, match="git_commit differs"):
        module.load_audited_producer_pins(bundle_dir)
    run_path.write_text(json.dumps(orig_run), encoding="utf-8")

    # 4. Wrong matrix path in generation.json
    bad_path_inv = dict(orig_gen["producer_matrix_inventory"], path="src/mtp/other.py")
    gen_path.write_text(
        json.dumps(dict(orig_gen, producer_matrix_inventory=bad_path_inv)), encoding="utf-8"
    )
    with pytest.raises(ValueError, match=r"must target src/mtp/taxa_trait_matrix\.py"):
        module.load_audited_producer_pins(bundle_dir)

    # 5. Matrix sha256 mismatch against core/run.json code_inventory
    bad_sha_inv = dict(orig_gen["producer_matrix_inventory"], sha256="a" * 64)
    gen_path.write_text(
        json.dumps(dict(orig_gen, producer_matrix_inventory=bad_sha_inv)), encoding="utf-8"
    )
    with pytest.raises(ValueError, match="code_inventory disagrees"):
        module.load_audited_producer_pins(bundle_dir)

    # 6. Matrix size_bytes mismatch against core/run.json code_inventory
    bad_size_inv = dict(orig_gen["producer_matrix_inventory"], size_bytes=99999)
    gen_path.write_text(
        json.dumps(dict(orig_gen, producer_matrix_inventory=bad_size_inv)), encoding="utf-8"
    )
    with pytest.raises(ValueError, match="code_inventory disagrees"):
        module.load_audited_producer_pins(bundle_dir)
