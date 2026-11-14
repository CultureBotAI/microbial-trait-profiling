"""PD policy provenance binds the metadata bytes consumed, without claiming core verification."""

from __future__ import annotations

import hashlib
import inspect
import json
import sys
from pathlib import Path
from types import ModuleType

import pandas as pd
import pytest

import test_contribution_trait_policy as fixtures
from mtp import contribution_inputs
from mtp.trait_key_policy import TraitKeyPolicy
from test_contribution_trait_policy import pd_contributors_mod as pd_contributors_mod


def _prepare_pd_main(
    root: Path,
    monkeypatch: pytest.MonkeyPatch,
    module: ModuleType,
    mode: str,
) -> tuple[Path, Path, Path, TraitKeyPolicy, str]:
    policy: TraitKeyPolicy = "facultative_union_v1" if mode == "default" else "raw_object_v1"
    with_manifest = mode in {"manifest", "both", "explicit_with_sources"}
    with_columns = mode in {"columns", "both", "explicit_with_sources"}
    inputs = root / "inputs"
    protocol = fixtures._write_synthetic_contribution_inputs(
        inputs, policy, with_run_manifest=with_manifest, with_column_metadata=with_columns
    )
    output = root / "output"
    output.mkdir()
    metadata = pd.DataFrame(
        [
            {
                "study_id": "S",
                "experiment_id": "E",
                "inc_dec": direction,
                "profile_id": f"S_E_{direction}",
                "condition": module.PD_CONDITION,
                "body_site_group": module.FECES,
            }
            for direction in ("inc", "dec")
        ]
    )
    traits = pd.DataFrame(
        [{"trait_key": fixtures.TRAIT_608, "trait_slug": "facultatively_aerobic"}]
    )
    cells = pd.DataFrame(columns=["rank", "taxon", "trait_slug"])
    selected = pd.DataFrame(
        columns=[
            "scope",
            "rank",
            "taxon",
            "name",
            "display_name",
            "pd_feces_profiles",
            "other_feces_profiles",
            "pd_other_feces_counts",
            "pd_feces_studies",
            "other_feces_studies",
            "pd_specificity_fraction",
            "pd_vs_other_feces_log2_ratio",
            "pd_feces_total_profiles",
            "other_feces_total_profiles",
            "pd_feces_prevalence",
            "other_feces_prevalence",
            "pd_vs_other_feces_log2_prevalence_ratio",
            "pd_minus_other_feces_prevalence_difference",
            "pd_other_feces_display",
        ]
    )
    for name, value in {
        "REPO_ROOT": root,
        "KG_DIR": inputs,
        "DERIVED_DIR": inputs,
        "KG_PROTOCOL_PATH": protocol,
        "OUTPUT_DIR": output,
    }.items():
        monkeypatch.setattr(module, name, value)
    monkeypatch.setattr(module, "load_profile_metadata", lambda: metadata.copy())
    monkeypatch.setattr(module, "require_current_trait_taxa_evidence", lambda *_a: None)
    monkeypatch.setattr(module, "read_tsv", lambda *_a: pd.DataFrame())
    monkeypatch.setattr(module, "trait_directions", lambda *_a: pd.DataFrame())
    monkeypatch.setattr(module, "check_trait_intervals", lambda *_a: None)
    monkeypatch.setattr(module, "trait_frame", lambda *_a: traits.copy())
    monkeypatch.setattr(
        module,
        "build_contributor_tables",
        lambda _c, cov, _t, _n: (cells.copy(), cov.copy(), selected.copy()),
    )
    monkeypatch.setattr(module, "augment_selected_taxa", lambda frame, _m: frame)
    monkeypatch.setattr(module, "build_direction_differences", lambda *_a: pd.DataFrame())
    monkeypatch.setattr(module, "select_large_difference_taxa", lambda *_a: selected.copy())
    monkeypatch.setattr(
        module,
        "build_trait_cocarriage_tables",
        lambda *_a: (pd.DataFrame(), pd.DataFrame(), pd.DataFrame()),
    )
    monkeypatch.setattr(module, "plot_contributor_heatmap", lambda *_a: None)
    monkeypatch.setattr(module, "plot_trait_cocarriage_network", lambda *_a: None)
    monkeypatch.setattr(module, "write_note", lambda *_a: None)
    arguments = ["pd"]
    if mode.startswith("explicit"):
        arguments += ["--trait-key-policy", policy]
    monkeypatch.setattr(sys, "argv", arguments)
    basis = (
        "explicit"
        if mode.startswith("explicit")
        else "core_manifest"
        if with_manifest
        else "column_metadata"
        if with_columns
        else "default"
    )
    return inputs, protocol, output, policy, basis


@pytest.mark.parametrize(
    "mode", ["explicit", "default", "columns", "manifest", "both", "explicit_with_sources"]
)
def test_pd_manifest_records_policy_resolution_sources(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    pd_contributors_mod: ModuleType,
    mode: str,
) -> None:
    inputs, _protocol, output, policy, basis = _prepare_pd_main(
        tmp_path, monkeypatch, pd_contributors_mod, mode
    )
    expected_hashes = {
        role: hashlib.sha256(path.read_bytes()).hexdigest()
        for role, path in (
            ("core_manifest", inputs / "run.json"),
            ("taxa_trait_matrix_columns.tsv", inputs / "taxa_trait_matrix_columns.tsv"),
        )
        if path.is_file()
    }

    pd_contributors_mod.main()

    manifest = json.loads((output / f"{pd_contributors_mod.STEM}_manifest.json").read_text())
    assert manifest["trait_key_policy"] == policy
    assert manifest.get("trait_key_policy_resolution") == {
        "policy": policy,
        "basis": basis,
        "source_sha256": expected_hashes,
    }
    assert manifest["ledger_profiles"] == 2


def test_pd_main_uses_one_captured_policy_resolution(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, pd_contributors_mod: ModuleType
) -> None:
    inputs, _protocol, output, policy, basis = _prepare_pd_main(
        tmp_path, monkeypatch, pd_contributors_mod, "both"
    )
    paths = {
        "core_manifest": inputs / "run.json",
        "taxa_trait_matrix_columns.tsv": inputs / "taxa_trait_matrix_columns.tsv",
    }
    first = {path: path.read_bytes() for path in paths.values()}
    reads = dict.fromkeys(first, 0)
    original = Path.read_bytes

    def changing_read(path: Path) -> bytes:
        if path not in first:
            return original(path)
        reads[path] += 1
        value = original(path)
        changed = value.replace(b"raw_object_v1", b"facultative_union_v1")
        path.write_bytes(changed if changed != value else value + b"\n")
        return value

    monkeypatch.setattr(Path, "read_bytes", changing_read)
    try:
        pd_contributors_mod.main()
    except ValueError as exc:
        raise AssertionError(
            "PD main must reuse the captured policy, not re-read metadata"
        ) from exc

    manifest = json.loads((output / f"{pd_contributors_mod.STEM}_manifest.json").read_text())
    assert manifest.get("trait_key_policy_resolution") == {
        "policy": policy,
        "basis": basis,
        "source_sha256": {
            role: hashlib.sha256(first[path]).hexdigest() for role, path in paths.items()
        },
    }
    assert reads == dict.fromkeys(first, 1)
    assert all(original(path) != value for path, value in first.items())


def test_pd_loader_refuses_policy_that_contradicts_captured_resolution(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, pd_contributors_mod: ModuleType
) -> None:
    inputs, protocol, _output, _policy, _basis = _prepare_pd_main(
        tmp_path, monkeypatch, pd_contributors_mod, "both"
    )
    resolver = getattr(contribution_inputs, "resolve_inputs_trait_key_policy_record", None)
    assert callable(resolver), "The captured policy resolver must be available"
    assert (
        "policy_resolution"
        in inspect.signature(pd_contributors_mod.load_contribution_ledger).parameters
    )
    resolution = resolver(inputs)
    message = ""
    try:
        pd_contributors_mod.load_contribution_ledger(
            inputs,
            inputs,
            protocol,
            "phenotype",
            trait_key_policy="facultative_union_v1",
            policy_resolution=resolution,
        )
    except ValueError as exc:
        message = str(exc)
    assert "contradicts" in message


def test_pd_policy_help_and_loader_do_not_claim_core_verification(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    pd_contributors_mod: ModuleType,
) -> None:
    monkeypatch.setattr(sys, "argv", ["pd", "--help"])
    with pytest.raises(SystemExit) as exc:
        pd_contributors_mod.main()
    assert exc.value.code == 0
    help_text = " ".join(capsys.readouterr().out.split()).lower()
    assert "policy metadata" in help_text
    assert "caller" in help_text
    assert "verify the core" in help_text
    assert "bind" in help_text
    assert "inputs" in help_text
    doc = " ".join((pd_contributors_mod.load_contribution_ledger.__doc__ or "").split()).lower()
    assert "policy metadata" in doc
    assert "caller" in doc
    assert "verify the core" in doc
    assert "bind" in doc
