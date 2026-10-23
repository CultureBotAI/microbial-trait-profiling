"""Write an explicit #190 influence assessment for comparison-stage tests (#844, 2026-10-22).

The native influence protocol only runs on the production KG bytes, so a synthetic test
core cannot produce Noel's native manifest. This writes the other shape integrated
evidence reads, a versioned procedure assessment, bound to the supplied inputs.

Inputs: the count profiles and experiments table the assessment is bound to.
Outputs: one taxon-influence-assessment-v1 JSON carrying those files' SHA-256 digests.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from mtp.influence_assessment import ASSESSMENT_VERSION
from mtp.paths import file_sha256


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profiles", type=Path, required=True)
    parser.add_argument("--experiments", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    assessment = {
        "schema_version": ASSESSMENT_VERSION,
        "input_sha256": {
            "study_trait_profiles_wide.tsv": file_sha256(args.profiles),
            "experiments.tsv": file_sha256(args.experiments),
        },
        "procedure_verdicts": {
            "fraction_sign_test": {
                "verdict": "fixture_assessment",
                "basis": "synthetic comparison-stage test fixture (#844)",
                "applies_to": "all_rows_of_procedure",
            }
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(assessment, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
