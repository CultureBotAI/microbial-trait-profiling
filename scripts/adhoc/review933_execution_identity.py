"""Probe execution-identity changes during PR #876 analysis (2026-10-23).

Inputs: existing synthetic workflow fixture and presence/depth implementation.
Outputs: JSON evidence of source/policy hash changes and real consumer verification.
Issue #953; refs #933/#922. Only a controlled in-memory hash provider changes; no source or
policy file is edited, and no production inputs are used.
"""

from __future__ import annotations

import argparse
import json
import runpy
from pathlib import Path
from typing import Any

import pytest

import mtp.presence_depth_sensitivity as sensitivity
from mtp.paths import REPO_ROOT
from mtp.workflow.consumers import verify_consumer_record


def probe(core: Path, output: Path, key: str) -> dict[str, Any]:
    """Change the observed hash during computation and check the published binding."""
    hashes = sensitivity.implementation_sources()
    original = hashes[key]
    changed = "f" * 64 if original != "f" * 64 else "e" * 64
    native = sensitivity.build_screen
    observed = []

    def inventory() -> dict[str, str]:
        observed.append(hashes[key])
        return dict(hashes)

    def analysis(*args: Any, **kwargs: Any) -> list[sensitivity.Result]:
        result = native(*args, **kwargs)
        hashes[key] = changed
        return result

    result: dict[str, Any] = {"changed_source": key, "original_hash": original}
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(sensitivity, "implementation_sources", inventory)
        patch.setattr(sensitivity, "build_screen", analysis)
        try:
            binding = sensitivity.run_sensitivity(core, output)
            verify_consumer_record(binding)
            record = json.loads((output / "provenance.json").read_text())
            result.update(
                rejected=False,
                consumer_verified=True,
                recorded_changed_hash=record["implementation_sources_sha256"][key] == changed,
            )
        except ValueError as error:
            result.update(rejected=True, consumer_verified=False, error=str(error))
    result.update(inventory_observations=observed, output_exists=output.exists())
    return result


def main() -> int:
    """Run both probes against a fresh synthetic core and retain their observations."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    inputs = output / "inputs"
    inputs.mkdir()
    tests = runpy.run_path(str(REPO_ROOT / "tests/test_presence_depth_sensitivity.py"))
    core = tests["_core"](inputs, True)
    results = [
        probe(core, output / name, key)
        for name, key in (
            ("code", "src/mtp/presence_depth_sensitivity.py"),
            ("policy", "methods/multiple_testing_policy_v1.yaml"),
        )
    ]
    rendered = json.dumps(results, indent=2, allow_nan=False) + "\n"
    (output / "evidence.json").write_text(rendered, encoding="utf-8")
    print(rendered)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
