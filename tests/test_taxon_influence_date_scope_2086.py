"""Regression tests for #2086: date-scope #188/#189 references in taxon influence docs."""

from __future__ import annotations

from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
RESULTS_MD = REPO_ROOT / "methods" / "table1_taxon_influence_results.md"
METHODS_MD = REPO_ROOT / "methods" / "table1_taxon_influence.md"
NATIVE_MD = REPO_ROOT / "methods" / "table1_native_influence.md"
AUDIT_MD = REPO_ROOT / "methods" / "biological_interpretation_audit.md"

CHECKPOINT_RELATIVE_LINK = "biological_interpretation_audit.md#status-checkpoint-2026-10-30"
CHECKPOINT_DOC_LABEL = "methods/biological_interpretation_audit.md#status-checkpoint-2026-10-30"


def _diagnostic_paragraph(text: str, lead_sentence: str) -> str:
    paragraphs = [p.strip() for p in text.split("\n\n") if p.strip()]
    matches = [p for p in paragraphs if p.startswith(lead_sentence)]
    assert len(matches) == 1, f"Expected unique paragraph starting with {lead_sentence!r}"
    return matches[0]


def _assert_status_checkpoint_closures(
    audit_text: str, expected_closures: tuple[tuple[str, str], ...]
) -> None:
    assert "## Status checkpoint (2026-10-30)" in audit_text, (
        "Missing status checkpoint section in biological_interpretation_audit.md"
    )
    checkpoint = audit_text.split("## Status checkpoint (2026-10-30)", 1)[1].split("\n## ", 1)[0]
    for issue_num, expected_date in expected_closures:
        target = f"https://github.com/CultureBotAI/microbial-trait-profiling/issues/{issue_num}"
        rows = [line for line in checkpoint.splitlines() if f"]({target})" in line]
        assert len(rows) == 1, f"Expected unique status checkpoint row for #{issue_num}"
        assert "| Closed |" in rows[0], f"Expected #{issue_num} to be Closed in status checkpoint"
        assert expected_date in rows[0], (
            f"Expected #{issue_num} closure date {expected_date} in status checkpoint"
        )


def _assert_date_scoped_taxon_influence_docs(
    results_text: str, methods_text: str, audit_text: str
) -> None:
    _assert_status_checkpoint_closures(audit_text, (("188", "2026-10-11"), ("189", "2026-10-08")))

    results_diag = _diagnostic_paragraph(
        results_text, "These are provisional annotation influence diagnostics."
    )
    methods_diag = _diagnostic_paragraph(
        methods_text, "This is a provisional annotation influence diagnostic."
    )

    for label, full_text, diag in (
        ("table1_taxon_influence_results.md", results_text, results_diag),
        ("table1_taxon_influence.md", methods_text, methods_diag),
    ):
        flat_lower_full = " ".join(full_text.split()).lower()
        for forbidden in (
            "remain unreviewed",
            "remain subject to #188",
            "feature semantics remain subject to #189",
            "the #188/#189 prerequisites",
            "prerequisites and review of resolved source semantics remain required",
        ):
            assert forbidden not in flat_lower_full, (
                f"Timeless present-tense phrase {forbidden!r} must not appear in {label}"
            )

        flat_diag = " ".join(diag.split())
        lower_diag = flat_diag.lower()
        assert "historical `190-v2` run" in flat_diag, (
            f"Expected date-scoping to the historical 190-v2 run in {label}"
        )
        assert f"({CHECKPOINT_RELATIVE_LINK})" in flat_diag, (
            f"Expected relative link to {CHECKPOINT_RELATIVE_LINK} in {label}"
        )
        assert CHECKPOINT_DOC_LABEL in flat_diag, (
            f"Expected explicit citation of {CHECKPOINT_DOC_LABEL} in {label}"
        )
        assert "#188" in flat_diag, f"Expected #188 reference in {label}"
        assert "2026-10-11" in flat_diag, f"Expected #188 closure date (2026-10-11 UTC) in {label}"
        assert "#189" in flat_diag, f"Expected #189 reference in {label}"
        assert "2026-10-08" in flat_diag, f"Expected #189 closure date (2026-10-08 UTC) in {label}"
        paired_188 = (
            "#188 (Closed 2026-10-11 UTC)"
            if label == "table1_taxon_influence_results.md"
            else "#188 Closed 2026-10-11 UTC"
        )
        paired_189 = (
            "#189 (Closed 2026-10-08 UTC)"
            if label == "table1_taxon_influence_results.md"
            else "#189 Closed 2026-10-08 UTC"
        )
        assert paired_188 in flat_diag, f"Expected paired {paired_188!r} in {label}"
        assert paired_189 in flat_diag, f"Expected paired {paired_189!r} in {label}"
        assert "does not retroactively recertify the historical `190-v2` outputs" in flat_diag, (
            f"Expected non-recertification statement for historical 190-v2 outputs in {label}"
        )
        assert "establish biological validation" in lower_diag, (
            f"Expected explicit biological-validation non-equivalence caveat in {label}"
        )

    methods_handoff = _diagnostic_paragraph(
        methods_text, "For #178, join on the **full qualified trait key**"
    )
    flat_handoff = " ".join(methods_handoff.split())
    lower_handoff = flat_handoff.lower()
    assert "as date-scoped above for #188 and #189" in lower_handoff, (
        "Expected back-reference to the earlier #188/#189 date-scoped paragraph in "
        "table1_taxon_influence.md"
    )
    assert "historical `190-v2` run" in flat_handoff, (
        "Expected closing handoff paragraph in table1_taxon_influence.md to name "
        "the historical 190-v2 run"
    )
    assert "(table1_native_influence.md)" in flat_handoff, (
        "Expected link to table1_native_influence.md in table1_taxon_influence.md"
    )
    assert "historical `190-native-v1` snapshot" in flat_handoff, (
        "Expected closing handoff paragraph in table1_taxon_influence.md to name "
        "the historical 190-native-v1 snapshot"
    )
    assert "does not retroactively recertify either historical output set" in lower_handoff, (
        "Expected closing non-recertification statement for both historical output sets in "
        "table1_taxon_influence.md"
    )
    for required_clause in (
        "matched-corpus analog scope (#210; not a reproduction of the published table 1 baseline)",
        "native-owner review",
        (
            "output-specific verification of contrast orientation and feature semantics for any "
            "`190-v2` or native outputs adopted in a manuscript"
        ),
        "review of resolved source semantics",
        "remain required before manuscript adoption or assigning mechanisms",
    ):
        assert required_clause in lower_handoff, (
            f"Expected remaining requirement {required_clause!r} in table1_taxon_influence.md"
        )


def test_taxon_influence_docs_date_scope_188_and_189_to_historical_190_v2_run() -> None:
    """#2086: both taxon influence methods docs date-scope #188/#189 to historical 190-v2."""
    _assert_date_scoped_taxon_influence_docs(
        RESULTS_MD.read_text(encoding="utf-8"),
        METHODS_MD.read_text(encoding="utf-8"),
        AUDIT_MD.read_text(encoding="utf-8"),
    )


@pytest.mark.parametrize(
    ("case_id", "target_doc", "old", "new"),
    [
        pytest.param(
            "results-timeless-remain-unreviewed",
            "results",
            "were unreviewed at execution time",
            "remain unreviewed",
            id="results-timeless-remain-unreviewed",
        ),
        pytest.param(
            "results-missing-historical-190-v2",
            "results",
            "In the historical `190-v2` run, contrast orientation (#188)",
            "In this run, contrast orientation (#188)",
            id="results-missing-historical-190-v2",
        ),
        pytest.param(
            "results-missing-checkpoint-link",
            "results",
            f"({CHECKPOINT_RELATIVE_LINK})",
            "(biological_interpretation_audit.md)",
            id="results-missing-checkpoint-link",
        ),
        pytest.param(
            "results-missing-non-recertification",
            "results",
            (
                "closing them in the tracker does not retroactively recertify the "
                "historical `190-v2` outputs or establish biological validation."
            ),
            "closing them in the tracker completes the historical review.",
            id="results-missing-non-recertification",
        ),
        pytest.param(
            "results-wrong-188-date",
            "results",
            "2026-10-11",
            "2026-10-12",
            id="results-wrong-188-date",
        ),
        pytest.param(
            "methods-timeless-remain-subject",
            "methods",
            "were scoped to #188 and feature\nsemantics to #189",
            "remain subject to #188; feature semantics remain subject to #189",
            id="methods-timeless-remain-subject",
        ),
        pytest.param(
            "methods-missing-historical-190-v2",
            "methods",
            "In the historical\n`190-v2` run, contrast eligibility",
            "In the current\nrun, contrast eligibility",
            id="methods-missing-historical-190-v2",
        ),
        pytest.param(
            "methods-missing-checkpoint-link",
            "methods",
            (
                "UTC; #189 Closed 2026-10-08 UTC; see the status checkpoint in\n"
                f"[`methods/{CHECKPOINT_RELATIVE_LINK}`]({CHECKPOINT_RELATIVE_LINK})"
            ),
            (
                "UTC; #189 Closed 2026-10-08 UTC; see the status checkpoint in\n"
                f"[`methods/{CHECKPOINT_RELATIVE_LINK}`](biological_interpretation_audit.md)"
            ),
            id="methods-missing-checkpoint-link",
        ),
        pytest.param(
            "methods-missing-non-recertification",
            "methods",
            (
                "closing them in the tracker does not retroactively recertify the "
                "historical `190-v2`\noutputs or establish biological validation.\n"
                "Increased/decreased"
            ),
            "closing them in the tracker updates the historical run.\nIncreased/decreased",
            id="methods-missing-non-recertification",
        ),
        pytest.param(
            "methods-wrong-189-date",
            "methods",
            "#189 Closed 2026-10-08 UTC",
            "#189 Closed 2026-10-09 UTC",
            id="methods-wrong-189-date",
        ),
        pytest.param(
            "methods-swapped-188-189-dates",
            "methods",
            "(#188 Closed 2026-10-11\nUTC; #189 Closed 2026-10-08 UTC",
            "(#188 Closed 2026-10-08\nUTC; #189 Closed 2026-10-11 UTC",
            id="methods-swapped-188-189-dates",
        ),
        pytest.param(
            "methods-handoff-timeless-prerequisites",
            "methods",
            (
                "As date-scoped above for #188 and #189 in the\n"
                "historical `190-v2` run (and in "
                "[the native rerun extension](table1_native_influence.md)\n"
                "for the historical `190-native-v1` snapshot), closing those tracker issues does\n"
                "not retroactively recertify either historical output set. Matched-corpus analog\n"
                "scope (#210; not a reproduction of the published Table 1 baseline), native-owner\n"
                "review, output-specific verification of contrast orientation and feature "
                "semantics\n"
                "for any `190-v2` or native outputs adopted in a manuscript, and review of "
                "resolved\n"
                "source semantics remain required before manuscript adoption or assigning\n"
                "mechanisms."
            ),
            (
                "Matched-corpus analog scope (#210; not a reproduction of the published Table 1 "
                "baseline), native-owner review, the #188/#189 prerequisites and review of "
                "resolved source semantics remain required before manuscript adoption or assigning "
                "mechanisms."
            ),
            id="methods-handoff-timeless-prerequisites",
        ),
        pytest.param(
            "methods-handoff-missing-backref-to-188-189",
            "methods",
            "As date-scoped above for #188 and #189 in the\n",
            "In the\n",
            id="methods-handoff-missing-backref-to-188-189",
        ),
        pytest.param(
            "methods-handoff-missing-190-v2",
            "methods",
            "historical `190-v2` run (",
            "earlier run (",
            id="methods-handoff-missing-190-v2",
        ),
        pytest.param(
            "methods-handoff-missing-190-native-v1",
            "methods",
            "for the historical `190-native-v1` snapshot)",
            "for the native extension)",
            id="methods-handoff-missing-190-native-v1",
        ),
        pytest.param(
            "methods-handoff-missing-native-link",
            "methods",
            "[the native rerun extension](table1_native_influence.md)",
            "the native rerun extension",
            id="methods-handoff-missing-native-link",
        ),
        pytest.param(
            "methods-handoff-missing-non-recertification",
            "methods",
            (
                "closing those tracker issues does\n"
                "not retroactively recertify either historical output set."
            ),
            "those tracker issues are closed.",
            id="methods-handoff-missing-non-recertification",
        ),
        pytest.param(
            "methods-handoff-deleted-analog-scope",
            "methods",
            (
                "Matched-corpus analog\n"
                "scope (#210; not a reproduction of the published Table 1 baseline), "
            ),
            "",
            id="methods-handoff-deleted-analog-scope",
        ),
        pytest.param(
            "methods-handoff-deleted-native-owner-review",
            "methods",
            "native-owner\nreview, ",
            "",
            id="methods-handoff-deleted-native-owner-review",
        ),
        pytest.param(
            "methods-handoff-deleted-orientation-and-semantics-requirement",
            "methods",
            (
                "output-specific verification of contrast orientation and feature semantics\n"
                "for any `190-v2` or native outputs adopted in a manuscript, "
            ),
            "",
            id="methods-handoff-deleted-orientation-and-semantics-requirement",
        ),
        pytest.param(
            "methods-handoff-deleted-orientation-only",
            "methods",
            "output-specific verification of contrast orientation and feature semantics",
            "output-specific verification of feature semantics",
            id="methods-handoff-deleted-orientation-only",
        ),
        pytest.param(
            "methods-handoff-deleted-resolved-source-semantics",
            "methods",
            ", and review of resolved\nsource semantics remain required",
            " remain required",
            id="methods-handoff-deleted-resolved-source-semantics",
        ),
        pytest.param(
            "audit-missing-checkpoint-heading",
            "audit",
            "## Status checkpoint (2026-10-30)",
            "## Historical status table",
            id="audit-missing-checkpoint-heading",
        ),
        pytest.param(
            "audit-188-not-closed",
            "audit",
            (
                "| [#188](https://github.com/CultureBotAI/microbial-trait-profiling/issues/188)"
                ", contrast orientation | Closed | 2026-10-11 |"
            ),
            (
                "| [#188](https://github.com/CultureBotAI/microbial-trait-profiling/issues/188)"
                ", contrast orientation | Open | 2026-10-11 |"
            ),
            id="audit-188-not-closed",
        ),
        pytest.param(
            "audit-189-wrong-date",
            "audit",
            (
                "| [#189](https://github.com/CultureBotAI/microbial-trait-profiling/issues/189)"
                ", trait semantics | Closed | 2026-10-08 |"
            ),
            (
                "| [#189](https://github.com/CultureBotAI/microbial-trait-profiling/issues/189)"
                ", trait semantics | Closed | 2026-10-10 |"
            ),
            id="audit-189-wrong-date",
        ),
        pytest.param(
            "results-missing-diagnostic-lead",
            "results",
            "These are provisional annotation influence diagnostics.",
            "Historical annotation influence diagnostics.",
            id="results-missing-diagnostic-lead",
        ),
    ],
)
def test_taxon_influence_date_scope_guard_rejects_unscoped_or_inconsistent_text(
    case_id: str, target_doc: str, old: str, new: str
) -> None:
    """#2086: guard rejects timeless present tense, missing links/dates, or missing caveats."""
    results_text = RESULTS_MD.read_text(encoding="utf-8")
    methods_text = METHODS_MD.read_text(encoding="utf-8")
    audit_text = AUDIT_MD.read_text(encoding="utf-8")
    if target_doc == "results":
        assert results_text.count(old) == 1, f"Expected unique anchor in results for {case_id}"
        results_text = results_text.replace(old, new, 1)
    elif target_doc == "methods":
        assert methods_text.count(old) == 1, f"Expected unique anchor in methods for {case_id}"
        methods_text = methods_text.replace(old, new, 1)
    else:
        assert audit_text.count(old) == 1, f"Expected unique anchor in audit for {case_id}"
        audit_text = audit_text.replace(old, new, 1)

    with pytest.raises(AssertionError):
        _assert_date_scoped_taxon_influence_docs(results_text, methods_text, audit_text)


def _assert_date_scoped_native_influence_doc(native_text: str, audit_text: str) -> None:
    _assert_status_checkpoint_closures(
        audit_text, (("188", "2026-10-11"), ("189", "2026-10-08"), ("182", "2026-10-12"))
    )

    lower_full = native_text.lower()
    for forbidden in (
        "#188 orientation and #189 semantics review remain required",
        "this work does not resolve those issues or #182/#184",
        "remain unreviewed",
    ):
        assert forbidden not in lower_full, (
            f"Stale unresolved wording {forbidden!r} must not appear in table1_native_influence.md"
        )

    handoff = _diagnostic_paragraph(
        native_text, "Supply these alongside the original ledger/shared effects to #178."
    )
    flat_handoff = " ".join(handoff.split())
    lower_handoff = flat_handoff.lower()
    assert "historical `190-native-v1` committed snapshot" in flat_handoff, (
        "Expected date-scoping to the historical 190-native-v1 snapshot in "
        "table1_native_influence.md"
    )
    assert f"({CHECKPOINT_RELATIVE_LINK})" in flat_handoff, (
        f"Expected relative link to {CHECKPOINT_RELATIVE_LINK} in table1_native_influence.md"
    )
    assert CHECKPOINT_DOC_LABEL in flat_handoff, (
        f"Expected explicit citation of {CHECKPOINT_DOC_LABEL} in table1_native_influence.md"
    )
    for paired_closure in (
        "#188 Closed 2026-10-11 UTC",
        "#189 Closed 2026-10-08 UTC",
        "closing #182 on 2026-10-12 UTC",
    ):
        assert paired_closure in flat_handoff, (
            f"Expected paired tracker closure {paired_closure!r} in table1_native_influence.md"
        )
    assert (
        "does not retroactively recertify the historical `190-native-v1` outputs" in flat_handoff
    ), "Expected non-recertification statement for historical 190-native-v1 outputs"
    assert "establish biological validation" in lower_handoff, (
        "Expected explicit biological-validation non-equivalence caveat"
    )
    for version in ("190-native-v1", "190-native-v2", "190-native-v3"):
        assert f"`{version}`" in flat_handoff, (
            f"Expected the manuscript-adoption requirements to name {version}"
        )
    for required_clause in (
        "native model owner review",
        (
            "output-specific verification of contrast orientation and trait semantics for any "
            "native results (`190-native-v1`, `190-native-v2` or `190-native-v3`) adopted in a "
            "manuscript"
        ),
        "source-count validity (#184)",
        "full-pipeline calibration of these native diagnostics",
        "remain required for manuscript adoption",
        "this handoff does not complete those reviews or resolve #184",
    ):
        assert required_clause in lower_handoff, (
            f"Expected remaining manuscript-adoption clause {required_clause!r} in "
            "table1_native_influence.md"
        )


def test_native_influence_doc_date_scopes_188_and_189_and_keeps_remaining_requirements() -> None:
    """#2086: table1_native_influence.md date-scopes #188/#189/#182 and preserves open work."""
    _assert_date_scoped_native_influence_doc(
        NATIVE_MD.read_text(encoding="utf-8"),
        AUDIT_MD.read_text(encoding="utf-8"),
    )


@pytest.mark.parametrize(
    ("case_id", "target_doc", "old", "new"),
    [
        pytest.param(
            "native-timeless-188-189-remain-required",
            "native",
            (
                "In the historical\n`190-native-v1` committed snapshot "
                "(`methods/table1-native-results/`)"
            ),
            (
                "#188 orientation and #189 semantics review remain required for manuscript "
                "adoption; this work does not resolve those issues or #182/#184. "
                "In the historical\n`190-native-v1` committed snapshot "
                "(`methods/table1-native-results/`)"
            ),
            id="native-timeless-188-189-remain-required",
        ),
        pytest.param(
            "native-missing-historical-190-native-v1",
            "native",
            "historical\n`190-native-v1` committed snapshot",
            "current native snapshot",
            id="native-missing-historical-190-native-v1",
        ),
        pytest.param(
            "native-missing-checkpoint-link",
            "native",
            f"({CHECKPOINT_RELATIVE_LINK})",
            "(biological_interpretation_audit.md)",
            id="native-missing-checkpoint-link",
        ),
        pytest.param(
            "native-missing-non-recertification",
            "native",
            (
                "retroactively recertify the historical `190-native-v1` outputs or establish\n"
                "biological validation."
            ),
            "require further historical review.",
            id="native-missing-non-recertification",
        ),
        pytest.param(
            "native-wrong-188-date",
            "native",
            "2026-10-11",
            "2026-10-14",
            id="native-wrong-188-date",
        ),
        pytest.param(
            "native-wrong-189-date",
            "native",
            "2026-10-08",
            "2026-10-09",
            id="native-wrong-189-date",
        ),
        pytest.param(
            "native-wrong-182-date",
            "native",
            "2026-10-12",
            "2026-10-15",
            id="native-wrong-182-date",
        ),
        pytest.param(
            "native-requirement-omits-v2",
            "native",
            "\n`190-native-v2` or ",
            "\n",
            id="native-requirement-omits-v2",
        ),
        pytest.param(
            "native-dropped-owner-review-requirement",
            "native",
            "Native model owner review,",
            "Prior handoff review,",
            id="native-dropped-owner-review-requirement",
        ),
        pytest.param(
            "native-swapped-188-189-dates",
            "native",
            "#188 Closed 2026-10-11 UTC; #189 Closed\n2026-10-08 UTC",
            "#188 Closed 2026-10-08 UTC; #189 Closed\n2026-10-11 UTC",
            id="native-swapped-188-189-dates",
        ),
        pytest.param(
            "native-deleted-owner-review-requirement",
            "native",
            "Native model owner review, output-specific",
            "Output-specific",
            id="native-deleted-owner-review-requirement",
        ),
        pytest.param(
            "native-deleted-orientation-semantics-requirement",
            "native",
            (
                "output-specific verification of\n"
                "contrast orientation and trait semantics for any native results "
                "(`190-native-v1`,\n`190-native-v2` or `190-native-v3`) adopted in a "
                "manuscript, "
            ),
            "",
            id="native-deleted-orientation-semantics-requirement",
        ),
        pytest.param(
            "native-deleted-orientation-semantics-version-moved",
            "native",
            (
                "output-specific verification of\n"
                "contrast orientation and trait semantics for any native results "
            ),
            "reporting of any native outputs ",
            id="native-deleted-orientation-semantics-version-moved",
        ),
        pytest.param(
            "native-deleted-source-count-validity-requirement",
            "native",
            "source-count validity\n(#184), ",
            "",
            id="native-deleted-source-count-validity-requirement",
        ),
        pytest.param(
            "native-deleted-calibration-requirement",
            "native",
            ", and full-pipeline calibration of these native diagnostics",
            "",
            id="native-deleted-calibration-requirement",
        ),
        pytest.param(
            "native-deleted-does-not-resolve-184",
            "native",
            "; this handoff does not complete those reviews or\nresolve #184.",
            ".",
            id="native-deleted-does-not-resolve-184",
        ),
        pytest.param(
            "audit-182-not-closed",
            "audit",
            (
                "| [#182](https://github.com/CultureBotAI/microbial-trait-profiling/issues/182)"
                ", calibration | Closed | 2026-10-12 |"
            ),
            (
                "| [#182](https://github.com/CultureBotAI/microbial-trait-profiling/issues/182)"
                ", calibration | Open | 2026-10-12 |"
            ),
            id="audit-182-not-closed",
        ),
    ],
)
def test_native_influence_date_scope_guard_rejects_unscoped_or_inconsistent_text(
    case_id: str, target_doc: str, old: str, new: str
) -> None:
    """#2086: guard rejects stale #188/#189/#182 phrasing or missing requirements."""
    native_text = NATIVE_MD.read_text(encoding="utf-8")
    audit_text = AUDIT_MD.read_text(encoding="utf-8")
    if target_doc == "native":
        assert native_text.count(old) == 1, f"Expected unique anchor in native for {case_id}"
        native_text = native_text.replace(old, new, 1)
    else:
        assert audit_text.count(old) == 1, f"Expected unique anchor in audit for {case_id}"
        audit_text = audit_text.replace(old, new, 1)

    with pytest.raises(AssertionError):
        _assert_date_scoped_native_influence_doc(native_text, audit_text)
