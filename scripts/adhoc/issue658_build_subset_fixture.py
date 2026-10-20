"""Add complete contrast metadata to the small raw-to-core fixture.

Input: tests/fixtures/workflow/full_dump.csv. Output: contrast_population.csv
beside it. Synthetic native/reversed/excluded cases serve #658 / #660; 2026-10-19.
"""

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FIXTURES = ROOT / "tests/fixtures/workflow"
METADATA = (
    "Study design",
    "PMID",
    "DOI",
    "EFO ID",
    "UBERON ID",
    "Group 1 definition",
    "Group 0 sample size",
    "Group 1 sample size",
    "Description",
)


def main() -> None:
    with (FIXTURES / "full_dump.csv").open() as stream:
        reader = csv.DictReader(line for line in stream if not line.startswith("#"))
        header = [*(reader.fieldnames or []), *METADATA]
        rows = list(reader)
    for row in rows:
        row.update(dict.fromkeys(METADATA, ""))
        row.update(
            {
                "Study design": "case-control",
                "Condition": "Asthma",
                "Group 0 name": "Healthy controls",
                "Group 1 name": "Asthma patients",
            }
        )
        if row["BSDB ID"].startswith("bsdb:1/2/"):
            row["Group 0 name"], row["Group 1 name"] = "Asthma patients", "Healthy controls"
        if row["BSDB ID"].startswith("bsdb:3/"):
            row["Group 0 name"], row["Group 1 name"] = "Placebo group", "Drug treatment group"
        row["Group 1 definition"] = row["Group 1 name"]
    with (FIXTURES / "contrast_population.csv").open("w", newline="") as stream:
        stream.write("# BugSigDB contrast-selection fixture-658, License: CC0 synthetic data\n")
        writer = csv.DictWriter(stream, fieldnames=header)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
