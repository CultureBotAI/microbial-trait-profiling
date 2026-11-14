"""Pin the docs/workflow.md site-stats provenance schema 1.1/1.2 contract (#1028)."""

from __future__ import annotations

import importlib.util
import json
import re
import sys
from pathlib import Path
from typing import Any

import pytest

ROOT = Path(__file__).resolve().parents[1]
WORKFLOW_DOC = ROOT / "docs/workflow.md"
PROVENANCE_JSON = (
    ROOT / "trait-visualization/trait-profile-investigation/stratified_stats_provenance.json"
)
VERIFIER_SCRIPT = ROOT / "scripts/adhoc/verify_1028_mutants_agy2_20261111.py"
SECTION_HEADING = "### Regenerate historical body-site composition tables (#824)"
LEGACY_REASON = "not refreshed by pooled/site statistics mode"
EXPECTED_SCHEMA_1_2_SENTENCE = (
    "Schema version 1.2 distinguishes `existing_stratum_table_not_refreshed` from "
    "`not_produced_by_site_stats_mode` (the summaries and figure)."
)
EXPECTED_SCHEMA_1_1_SENTENCE = (
    "The committed "
    "`trait-visualization/trait-profile-investigation/stratified_stats_provenance.json` "
    "written by the final rerun (#1004) is schema version 1.1, with all eight of its "
    "`not_regenerated` entries using the single legacy reason "
    '`"not refreshed by pooled/site statistics mode"`; the next '
    "`just trait-profile-site-stats` run rewrites it compatibly under schema version 1.2, "
    "and nothing reads it."
)


def _site_stats_section(markdown: str) -> str:
    assert SECTION_HEADING in markdown
    after_heading = markdown.split(SECTION_HEADING, 1)[1]
    section = after_heading.split("\n## ", 1)[0]
    return re.sub(r"\s+", " ", section).strip()


def _load_verifier() -> Any:
    spec = importlib.util.spec_from_file_location("verify_1028_mutants", VERIFIER_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {VERIFIER_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_workflow_documents_committed_schema_1_1_and_schema_1_2_site_stats_provenance() -> None:
    """docs/workflow.md states the committed 1.1 record (#1004) and 1.2 rewrite (#1028)."""
    section = _site_stats_section(WORKFLOW_DOC.read_text(encoding="utf-8"))
    record = json.loads(PROVENANCE_JSON.read_text(encoding="utf-8"))

    assert record["schema_version"] == "1.1"
    not_regenerated = record["not_regenerated"]
    assert isinstance(not_regenerated, list)
    assert len(not_regenerated) == 8
    assert [entry["reason"] for entry in not_regenerated] == [LEGACY_REASON] * 8

    assert EXPECTED_SCHEMA_1_2_SENTENCE in section
    assert EXPECTED_SCHEMA_1_1_SENTENCE in section
    assert f"{EXPECTED_SCHEMA_1_2_SENTENCE} {EXPECTED_SCHEMA_1_1_SENTENCE}" in section


def test_verifier_rejects_setup_error_and_call_value_error(tmp_path: Path) -> None:
    """Mutant verifier rejects setup errors and call-phase non-AssertionErrors."""
    verifier = _load_verifier()
    doc_text = WORKFLOW_DOC.read_text(encoding="utf-8")
    target = verifier.NODE

    setup_error_test = (
        "import pytest\n"
        "@pytest.fixture(autouse=True)\n"
        "def broken_setup() -> None:\n"
        '    raise ValueError("AssertionError in setup")\n'
        f"def {target}() -> None:\n"
        "    assert True\n"
    )
    with pytest.raises(RuntimeError, match="Setup/collection/teardown error or skip"):
        verifier.run_case(
            tmp_path,
            "setup_error_probe",
            target,
            doc_text,
            True,
            test_override=setup_error_test,
        )

    call_value_error_test = (
        f'def {target}() -> None:\n    raise ValueError("AssertionError in call phase")\n'
    )
    with pytest.raises(RuntimeError, match="Not a target call-phase AssertionError"):
        verifier.run_case(
            tmp_path,
            "call_value_error_probe",
            target,
            doc_text,
            True,
            test_override=call_value_error_test,
        )
