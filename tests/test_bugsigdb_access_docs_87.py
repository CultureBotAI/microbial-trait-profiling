"""Keep the #87 access decision complete and its benchmark scope explicit."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METHODS = Path("docs/bugsigdb-access-methods.md")
BENCHMARK = Path("docs/bugsigdb-access-benchmark.md")


def check_docs(root: Path) -> None:
    methods = (root / METHODS).read_text()
    benchmark = (root / BENCHMARK).read_text()
    assert "## Decision matrix (#87)" in methods
    matrix = methods.split("## Decision matrix (#87)", 1)[1].split("## ", 1)[0]
    assert (
        "| Access path | Fidelity | Freshness cadence | Field completeness | Endpoint stability |"
        in matrix
    )
    for name in (
        "Help:Export",
        "BugSigDBExports/full_dump.csv",
        "BugSigDBExports/*.gmt",
        "10.5281/zenodo.10627578",
    ):
        rows = [line for line in matrix.splitlines() if line.startswith("|") and name in line]
        assert len(rows) == 1, name
        assert len(rows[0].split("|")) == 7, name
        assert all(cell.strip() for cell in rows[0].split("|")[1:-1])
    assert "Special:Ask" in matrix
    assert "hourly" in matrix
    assert "not a service-level guarantee" in matrix
    assert "immutable" in matrix
    assert "devel" in matrix
    assert "## Field-by-field completeness" in methods
    fields = methods.split("## Field-by-field completeness", 1)[1].split("## ", 1)[0]
    for field in (
        "Study page name",
        "Study design",
        "PMID",
        "DOI",
        "URL",
        "Authors list",
        "Title",
        "Journal",
        "Year",
        "Abstract",
        "Keywords",
        "State",
        "Reviewer",
        "Curation method",
        "Curation protocol",
        "submittedBy",
        "submissionDate",
        "Review state",
        "Quality control tag",
        "Experiment page name",
        "Study",
        "Location of subjects",
        "Host species",
        "Body site",
        "UBERON ID",
        "Condition",
        "EFO ID",
        "Group 0 name",
        "Group 1 name",
        "Group 1 definition",
        "Group 0 sample size",
        "Group 1 sample size",
        "Antibiotics exclusion",
        "Sequencing type",
        "16S variable region",
        "Sequencing platform",
        "Data transformation",
        "Statistical test",
        "Significance threshold",
        "MHT correction",
        "LDA Score above",
        "Matched on",
        "Confounders controlled for",
        "Pielou",
        "Shannon",
        "Chao1",
        "Simpson",
        "Inverse Simpson",
        "Richness",
        "Signature page name",
        "Experiment",
        "Source",
        "Curated date",
        "Curator",
        "Revision editor",
        "Description",
        "Abundance in Group 1",
        "MetaPhlAn taxon names",
        "NCBI Taxonomy IDs",
        "BSDB ID",
    ):
        assert f"`{field}`" in fields, field
    assert "absent from v1.2.1" in fields
    assert "not separate GMT columns" in fields
    assert "v1.2.1" in methods
    assert "2024-02-06" in methods
    assert "2024-02-04_18:58_UTC" in methods
    assert "## Use-case verdicts" in methods
    for verdict in (
        "Canonical ingestion",
        "Live curation audit",
        "Signature-set analysis",
        "Historical reproduction",
    ):
        assert verdict in methods
    assert "untested" in methods
    assert "bugsigdbr" in methods
    assert "#62" in methods
    for evidence in (
        "oldid=259555",
        "5815c2638d6b56134ac8db19edc2bb3c4b17f1ca",
        "2523676c1bb667e524dad5a8b76bc0a97bbb4bbf",
        "api/records/10627578",
    ):
        assert evidence in methods
    assert "bugsigdb-access-methods.md#decision-matrix-87" in benchmark
    assert "not a benchmark of v1.2.1" in benchmark
    assert "10.5281/zenodo.10627578" in benchmark
    assert "single-run" in benchmark


def test_access_decision_record() -> None:
    check_docs(ROOT)
