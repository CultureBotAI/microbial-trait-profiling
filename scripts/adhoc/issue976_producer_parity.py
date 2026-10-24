"""Compare #976 shared interfaces with an immutable producer revision (2026-10-23).

Inputs: --baseline-ref and deterministic synthetic counts (no production data).
Outputs: isolated baseline source copies, fixture TSVs and parity.json. Every
native statistic/selection value and old BH/BY arithmetic is compared exactly.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path
from types import ModuleType

import numpy as np
import pandas as pd
from git import Repo

from mtp import native_composition
from mtp.enrichment_stats import benjamini_hochberg
from mtp.paths import REPO_ROOT, file_sha256
from mtp.trait_screen import step_up, step_up_array


def load_module(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    if module.__file__ is None or Path(module.__file__).resolve() != path.resolve():
        raise RuntimeError("Unexpected producer import identity")
    return module


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-ref", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    repository = Repo(REPO_ROOT)
    commit = repository.commit(args.baseline_ref)
    copied = {}
    for name, relative in {
        "producer": "scripts/create_trait_profile_infographic.py",
        "correction": "src/mtp/trait_screen.py",
    }.items():
        path = output / f"baseline_{name}.py"
        path.write_bytes((commit.tree / relative).data_stream.read())
        copied[name] = path
    old = load_module("issue976_original_producer", copied["producer"])
    previous = load_module("issue976_original_correction", copied["correction"])
    current_path = REPO_ROOT / "scripts/create_trait_profile_infographic.py"
    current = load_module("issue976_current_producer", current_path)
    watched = [current_path, Path(native_composition.__file__), Path(__file__)]
    before = {str(path): file_sha256(path) for path in watched}

    features = ["METPO:one", "qualified|GO:two", "ubiquitous", "rare", "mediadive.medium:X"]
    records = []
    metadata = []
    for pair in range(12):
        study = "shared-a" if pair < 3 else "shared-b" if pair < 5 else f"s{pair}"
        experiment = f"e{pair}"
        metadata.append(
            {
                "study_id": study,
                "experiment_id": experiment,
                "condition": "fixture",
                "body_site": "Feces",
                "host_species": "Homo sapiens",
            }
        )
        for direction in ("inc", "dec"):
            index = 2 * pair + int(direction == "dec")
            records.append(
                {
                    "study_id": study,
                    "experiment_id": experiment,
                    "inc_dec": direction,
                    features[0]: 1 + pair % 4 if index % 5 else 0,
                    features[1]: 1 + pair % 3 if index % 4 else 0,
                    "ubiquitous": 50 + index if index else 0,
                    "rare": int(0 < index < 6),
                    "mediadive.medium:X": 1000,
                }
            )
    profiles, experiments, labels = (
        output / name for name in ("profiles.tsv", "experiments.tsv", "labels.tsv")
    )
    pd.DataFrame(records).to_csv(profiles, sep="\t", index=False)
    pd.DataFrame(metadata).to_csv(experiments, sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "column_id": feature,
                "column_label": feature,
                "object_name": feature,
                "object_category": "biolink:PhenotypicQuality",
            }
            for feature in features
        ]
    ).to_csv(labels, sep="\t", index=False)
    outputs = []
    for producer in (old, current):
        counts = producer.load_profile_matrix(profiles)
        pairs = producer.core_pair_metadata(profiles, experiments)
        names = producer.load_trait_labels(labels)
        selected = producer.select_stratum_traits(counts, names)
        standardized = producer.build_zscore_matrix(counts, selected.index)
        paired = producer.build_pair_matrix_from_profiles(standardized, pairs)
        stats = producer.compute_stratum_trait_stats(
            counts, pairs, names, stratum_kind="pooled", stratum="pooled"
        )
        outputs.append(
            {
                "counts": counts,
                "pairs": pairs,
                "labels": names,
                "selected": selected,
                "standardized": standardized,
                "paired": paired,
                "stats": stats,
            }
        )
    for name in outputs[0]:
        pd.testing.assert_frame_equal(outputs[0][name], outputs[1][name], check_exact=True)
    if outputs[0]["stats"].to_csv(index=False) != outputs[1]["stats"].to_csv(index=False):
        raise AssertionError("Native exported statistics changed")
    if "ubiquitous" in outputs[1]["selected"].index or "rare" in outputs[1]["selected"].index:
        raise AssertionError("Synthetic selection does not exercise excluded denominator traits")
    if (outputs[1]["counts"].sum(axis=1) == 0).sum() != 1:
        raise AssertionError("Synthetic population does not exercise a retained zero-total profile")
    generator = np.random.default_rng(976)
    comparisons = 0
    for size in (0, 1, 2, 17, 317, 4096):
        values = generator.random(size)
        if size >= 6:
            values[:6] = [0.0, np.nextafter(0.0, 1.0), 1.0, 0.001, 0.001, 1e-300]
        for family_size in sorted({1, max(1, size // 2), max(1, size), max(1, size * 2)}):
            for method in ("bh", "by"):
                expected = previous.step_up(values.tolist(), method=method, family_size=family_size)
                actual = step_up_array(values, method=method, family_size=family_size)
                np.testing.assert_array_equal(actual, expected)
                np.testing.assert_array_equal(
                    step_up(values.tolist(), method=method, family_size=family_size), expected
                )
                if size:
                    np.testing.assert_array_equal(
                        benjamini_hochberg(values, method=method, family_size=family_size), expected
                    )
                comparisons += 1
    after = {str(path): file_sha256(path) for path in watched}
    if before != after:
        raise AssertionError("Canonical producer/helper bytes changed during parity probe")
    record = {
        "baseline_commit": commit.hexsha,
        "baseline_sha256": {name: file_sha256(path) for name, path in copied.items()},
        "canonical_before": before,
        "canonical_after": after,
        "verified": True,
        "native_frames_exact": list(outputs[0]),
        "native_export_exact": True,
        "profiles": len(outputs[1]["counts"]),
        "pairs": len(outputs[1]["pairs"]),
        "studies": int(outputs[1]["pairs"].study_id.nunique()),
        "selected_features": list(outputs[1]["selected"].index),
        "step_up_comparisons": comparisons,
    }
    (output / "parity.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps({"verified": True, "step_up_comparisons": comparisons}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
