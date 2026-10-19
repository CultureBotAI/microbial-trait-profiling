"""Audit #664 before/after stratification artifacts (2026-10-19).

Inputs: old and corrected stratification output directories. Output: JSON on
stdout with tested-universe accounting, density rates, and an example denominator.
This preserves the calculation behind the issue/PR production validation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from mtp.enrichment_results import parse_metadata_lines


def read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", comment="#")


def file_hash(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def audit(before: Path, after: Path) -> dict[str, object]:
    """Require complete accounting and return reviewable numerical evidence."""
    example = "METPO:1000487"
    strata: dict[str, object] = {}
    for stratum in ("pooled", "feces", "feces_excluded"):
        filename = f"trait_enrichment_traits_overall_ranked_{stratum}.tsv"
        old_ranked = read_table(before / filename).set_index("trait")
        new_ranked = read_table(after / filename).set_index("trait")
        tested = read_table(after / f"trait_enrichment_traits_tested_{stratum}.tsv")
        metadata_lines = []
        with (after / f"trait_enrichment_all_results_{stratum}.tsv").open() as handle:
            for line in handle:
                if not line.startswith("#"):
                    break
                metadata_lines.append(line)
        meta = parse_metadata_lines(metadata_lines)
        results_filename = f"trait_enrichment_all_results_{stratum}.tsv"
        before_hash = file_hash(before / results_filename)
        after_hash = file_hash(after / results_filename)
        if before_hash != after_hash:
            raise ValueError(f"{stratum}: analyzer results differ; inspect before comparison")
        tests = int(meta["Total-Tests-Performed"])
        hits = int(meta["Significant-Hits"])
        if tested["n_profiles_tested"].sum() != tests or tested["hit_rows"].sum() != hits:
            raise ValueError(
                f"{stratum}: complete summary does not account for analyzer tests/hits"
            )
        eligible = tested.loc[tested["n_profiles_tested"] > 0]
        entry: dict[str, object] = {
            "tests": tests,
            "discoveries": hits,
            "eligible_traits": len(eligible),
            "untested_input_traits": len(tested) - len(eligible),
            "zero_discovery_traits": int(eligible["hit_rows"].eq(0).sum()),
            "one_discovery_traits": int(eligible["hit_rows"].eq(1).sum()),
            "old_ranked_traits": len(old_ranked),
            "new_ranked_traits": len(new_ranked),
            "unchanged_all_results_sha256": after_hash,
        }
        if example in old_ranked.index and example in new_ranked.index:
            entry["example_trait"] = example
            entry["old_prevalence"] = old_ranked.loc[example, "prevalence_summary"]
            entry["corrected_prevalence"] = new_ranked.loc[example, "prevalence_summary"]
            entry["old_trait_total"] = int(old_ranked.loc[example, "trait_total"])
            entry["reference_k"] = int(new_ranked.loc[example, "trait_total"])
        strata[stratum] = entry
    deciles = read_table(after / "annotation_density_deciles.tsv")
    pooled = read_table(after / "trait_enrichment_traits_tested_pooled.tsv")
    if deciles["total_traits"].sum() != pooled["n_profiles_tested"].gt(0).sum():
        raise ValueError("Density diagnostic omitted eligible traits")
    return {
        "before": str(before.resolve()),
        "after": str(after.resolve()),
        "strata": strata,
        "density_deciles": deciles.to_dict(orient="records"),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("before", type=Path)
    parser.add_argument("after", type=Path)
    args = parser.parse_args()
    print(json.dumps(audit(args.before, args.after), indent=2))


if __name__ == "__main__":
    main()
