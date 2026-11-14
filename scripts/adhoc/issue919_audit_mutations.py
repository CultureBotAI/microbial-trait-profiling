"""Check #919 survivor mutations against audit tests (2026-10-23).

Inputs: current label-support source and synthetic tests. Outputs: pytest logs
and JSON mutation outcomes in a new directory outside this checkout. Mutates
only a scratch module copy; run each mutation in a fresh interpreter (#1184).
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import sys
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from types import ModuleType

import pytest
from mutation_probe_isolation import isolated_module, verify_checkout_identity

ROOT = Path(__file__).resolve().parents[2]
TEST = ROOT / "tests/test_label_support_audit.py"


class ImportOrigin:
    """Require both module and direct function imports to use the scratch copy."""

    def __init__(self, module: ModuleType) -> None:
        self.module = module
        self.verified = False

    def pytest_collection_modifyitems(self, items: list[pytest.Item]) -> None:
        for item in items:
            test_module = getattr(item, "module", None)
            if (
                test_module is None
                or Path(test_module.__file__).resolve() != TEST.resolve()
                or test_module.label_support_audit is not self.module
                or test_module.run_label_support_audit is not self.module.run_label_support_audit
            ):
                raise AssertionError("Collected test does not use the isolated audit module")
        self.verified = bool(items)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--mutation", choices=("all_input_taxa", "expanded_focal_family"), required=True
    )
    args = parser.parse_args()
    verify_checkout_identity(ROOT)
    if "mtp.label_support_audit" in sys.modules:
        raise ValueError("Run each mutation in a fresh interpreter before importing the audit")
    source = ROOT / "src/mtp/label_support_audit.py"
    original = source.read_bytes()
    text = original.decode("utf-8")
    bound = 'bound_taxa = set(occurrences.loc[pair_keys.isin(paired), "taxon_id"])'
    focal = "FOCAL_TRAITS = tuple(AUDIT_TRAITS)\n"
    if text.count(bound) != 1 or text.count(focal) != 1 or text.count("FOCAL_TAXA = {") != 1:
        raise ValueError("Mutation targets changed; inspect the current implementation")
    mutants = {
        "all_input_taxa": (
            text.replace(bound, 'bound_taxa = set(occurrences["taxon_id"])'),
            "test_all_support_rules_share_zero_total_filtered_pairs",
        ),
        "expanded_focal_family": (
            text.replace(focal, "").replace("FOCAL_TAXA = {", focal + "FOCAL_TAXA = {"),
            "test_sidecar_records_sources_reconciliations_and_study_means",
        ),
    }
    name = args.mutation
    mutated, test = mutants[name]
    output = io.StringIO()
    with isolated_module(
        ROOT, "mtp.label_support_audit", mutated.encode("utf-8"), args.output_dir, (TEST,)
    ) as module:
        origin = ImportOrigin(module)
        with redirect_stdout(output), redirect_stderr(output):
            code = pytest.main(
                [
                    "-q",
                    "-p",
                    "no:cacheprovider",
                    "--basetemp",
                    str(args.output_dir.resolve() / "fixtures"),
                    f"{TEST}::{test}",
                ],
                plugins=[origin],
            )
        if not origin.verified:
            raise AssertionError("Pytest did not verify the isolated audit import")
        loaded_origin = str(Path(str(module.__file__)).resolve())
    if source.read_bytes() != original:
        raise AssertionError("Canonical audit source changed during the isolated probe")
    (args.output_dir / f"{name}.log").write_text(output.getvalue(), encoding="utf-8")
    record = {
        "original_source_sha256": hashlib.sha256(original).hexdigest(),
        "restored_source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "canonical_unchanged": True,
        "loaded_module_origin": loaded_origin,
        "mutation": name,
        "pytest_exit": int(code),
        "killed": code == pytest.ExitCode.TESTS_FAILED,
    }
    (args.output_dir / f"{name}.json").write_text(
        json.dumps(record, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(record, indent=2))
    return 0 if record["killed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
