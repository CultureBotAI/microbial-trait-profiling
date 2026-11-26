"""Pin the docs/workflow.md site-stats provenance contract (#1028, #2432).

Since #857 a full render deletes stratified_stats_provenance.json, because it rewrites the
tables the record lists; `just trait-profile-site-stats` writes a fresh schema 1.2 record for
its own run. The committed schema 1.1 record from the final rerun (#1004) was retired in the
#2397 regeneration (#2432). The deletion itself is tested in
tests/test_trait_profile_infographic.py::test_other_modes_invalidate_old_site_provenance.
"""

from __future__ import annotations

import importlib.util
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
EXPECTED_SCHEMA_1_2_SENTENCE = (
    "Schema version 1.2 distinguishes `existing_stratum_table_not_refreshed` from "
    "`not_produced_by_site_stats_mode` (the summaries and figure)."
)
EXPECTED_RETIRED_SENTENCE = (
    "No record is committed: a full render (`just trait-profile-investigation`) deletes "
    "`trait-visualization/trait-profile-investigation/stratified_stats_provenance.json`, "
    "because it rewrites the tables the record lists, so the schema 1.1 record from the "
    "final rerun (#1004) was retired in the #2397 regeneration (#2432); each "
    "`just trait-profile-site-stats` run writes a fresh record under schema version 1.2, "
    "and nothing reads it."
)
RETIRED_COMMITTED_CLAIM = "is schema version 1.1, with all eight of its"


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


def test_workflow_documents_retired_record_and_schema_1_2_site_stats_provenance() -> None:
    """No record is committed after a full render; docs/workflow.md says so (#1028, #2432)."""
    section = _site_stats_section(WORKFLOW_DOC.read_text(encoding="utf-8"))

    assert not PROVENANCE_JSON.exists()
    assert EXPECTED_SCHEMA_1_2_SENTENCE in section
    assert EXPECTED_RETIRED_SENTENCE in section
    assert f"{EXPECTED_SCHEMA_1_2_SENTENCE} {EXPECTED_RETIRED_SENTENCE}" in section
    assert RETIRED_COMMITTED_CLAIM not in section


def test_verifier_rejects_setup_error_and_call_value_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Mutant verifier rejects setup errors and call-phase non-AssertionErrors.

    The dated #1028 verifier copies the committed record into each scratch case. That record
    is retired (#2432), so a placeholder stands in for it under a temporary root; these two
    probes replace the test file and never read the record.
    """
    verifier = _load_verifier()
    doc_text = WORKFLOW_DOC.read_text(encoding="utf-8")
    placeholder_root = tmp_path / "root"
    (placeholder_root / verifier.PROVENANCE).parent.mkdir(parents=True)
    (placeholder_root / verifier.PROVENANCE).write_text("{}\n", encoding="utf-8")
    monkeypatch.setattr(verifier, "ROOT", placeholder_root)
    runs = tmp_path / "runs"
    runs.mkdir()
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
            runs,
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
            runs,
            "call_value_error_probe",
            target,
            doc_text,
            True,
            test_override=call_value_error_test,
        )
