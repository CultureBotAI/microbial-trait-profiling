"""Plot reported group sample sizes for the experiments behind the oral heatmap."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator, ScalarFormatter

from create_condition_body_site_guild_maps import FACULTATIVE_CONFIG, v3_body_site_display_label
from create_six_shared_trait_cytoscape_heatmaps import select_oral_contrast_cells

ROOT = Path(__file__).resolve().parents[1]
VIS = ROOT / "trait-visualization"
SIX = VIS / "trait-profile-investigation" / "six-shared-traits"
OUTPUT = SIX / "oral" / "sample-sizes"
KEYS = ["study_id", "experiment_id"]
SIZE_COLUMNS = ["Group 0 sample size", "Group 1 sample size"]
MISSING = {"", "na", "nan", "none", "unknown"}


def read_raw_dump(path: Path) -> pd.DataFrame:
    """Skip only the export preamble, preserving URL fragments and quoted newlines."""
    with path.open(encoding="utf-8", newline="") as handle:
        while True:
            position = handle.tell()
            line = handle.readline()
            if not line.startswith("#"):
                handle.seek(position)
                break
        return pd.read_csv(handle, dtype=str, keep_default_na=False)


def oral_experiments(metadata: pd.DataFrame, cells: pd.DataFrame) -> pd.DataFrame:
    """Join the exact saved oral condition/site rows to their contributing profiles."""
    oral = select_oral_contrast_cells(cells)
    row_keys = set(zip(oral["condition"], oral["body_site_group"], strict=True))
    profiles = metadata.loc[metadata["inc_dec"].isin(["inc", "dec"])].copy()
    profiles["heatmap_body_site_group"] = [
        v3_body_site_display_label(row.body_site, row.body_site_group, FACULTATIVE_CONFIG)[0]
        for row in profiles.itertuples()
    ]
    profiles = profiles.loc[
        [(row.condition, row.heatmap_body_site_group) in row_keys for row in profiles.itertuples()]
    ]
    found = set(zip(profiles["condition"], profiles["heatmap_body_site_group"], strict=True))
    if found != row_keys:
        raise ValueError(f"Oral heatmap rows missing from profile metadata: {row_keys - found}")
    columns = [*KEYS, "condition", "body_site", "heatmap_body_site_group", "host_species"]
    if profiles[KEYS].eq("").any().any():
        raise ValueError("Missing study or experiment ID in oral profile metadata")
    distinct = profiles[columns].drop_duplicates()
    if distinct.duplicated(KEYS).any():
        raise ValueError("Conflicting oral profile metadata within an experiment")
    directions = profiles.groupby(KEYS)["inc_dec"].agg(lambda x: ";".join(sorted(set(x))))
    result = distinct.merge(directions.rename("profile_directions"), on=KEYS, validate="1:1")
    return result.sort_values(KEYS).reset_index(drop=True)


def parse_size(value: str) -> tuple[int | None, str]:
    """Require a positive integer count; never map missing or invalid sizes to zero."""
    if value.strip().lower() in MISSING:
        return None, "missing"
    try:
        number = float(value)
    except ValueError:
        return None, "invalid"
    if not math.isfinite(number) or number <= 0 or not number.is_integer():
        return None, "invalid"
    return int(number), "valid"


def sample_size_audit(experiments: pd.DataFrame, raw: pd.DataFrame) -> pd.DataFrame:
    """Deduplicate signatures, rejecting ambiguous group counts instead of summing them."""
    raw = raw.copy()
    parts = raw["BSDB ID"].str.extract(r"^(bsdb:[^/]+)/([^/]+)/[^/]+$")
    if parts.isna().any().any():
        bad = raw.loc[parts.isna().any(axis=1), "BSDB ID"].unique().tolist()
        raise ValueError(f"Malformed BSDB ID in source dump: {bad[:20]}")
    raw["study_id"] = parts[0]
    raw["experiment_id"] = parts[0] + "-" + parts[1]
    groups = dict(iter(raw.groupby(KEYS)))
    records = []
    for record in experiments.to_dict("records"):
        group = groups.get((record["study_id"], record["experiment_id"]))
        record["source_signature_rows"] = 0 if group is None else len(group)
        reasons = []
        for i, column in enumerate(SIZE_COLUMNS):
            observed = [] if group is None else sorted(set(group[column].str.strip()))
            record[f"group_{i}_source_values"] = json.dumps(observed)
            parsed = [parse_size(value) for value in observed]
            counts = {value for value, status in parsed if status == "valid"}
            statuses = {status for _, status in parsed}
            if group is None:
                reason = "unmatched"
            elif len(counts) > 1:
                reason = "conflicting"
            elif "invalid" in statuses:
                reason = "invalid"
            elif "missing" in statuses or not counts:
                reason = "missing"
            else:
                reason = "valid"
            record[f"group_{i}_n"] = next(iter(counts)) if reason == "valid" else None
            record[f"group_{i}_status"] = reason
            if reason != "valid":
                reasons.append(f"group_{i}:{reason}")
        record["size_status"] = ";".join(reasons) if reasons else "included"
        record["total_n"] = record["group_0_n"] + record["group_1_n"] if not reasons else None
        for name in ["Group 0 name", "Group 1 name", "Study design"]:
            record[name] = "" if group is None else "; ".join(sorted(set(group[name])))
        records.append(record)
    audit = pd.DataFrame(records)
    for column in ["group_0_n", "group_1_n", "total_n"]:
        audit[column] = audit[column].astype("Int64")
    return audit


def summarize(audit: pd.DataFrame) -> dict[str, int | float]:
    """Describe observed experiment counts and the location of the external reference."""
    valid = audit.loc[audit["size_status"].eq("included")]
    sizes = valid["total_n"].astype(int)
    if sizes.empty:
        raise ValueError("No oral experiments have two valid, unambiguous group sizes")
    return {
        "oral_experiments": len(audit),
        "oral_studies": audit["study_id"].nunique(),
        "experiments_with_sizes": len(valid),
        "studies_with_sizes": valid["study_id"].nunique(),
        "experiments_excluded": len(audit) - len(valid),
        "minimum_n": int(sizes.min()),
        "q25_n": float(sizes.quantile(0.25)),
        "median_n": float(sizes.median()),
        "q75_n": float(sizes.quantile(0.75)),
        "maximum_n": int(sizes.max()),
        "below_200": int(sizes.lt(200).sum()),
        "equal_200": int(sizes.eq(200).sum()),
        "above_200": int(sizes.gt(200).sum()),
        "percent_below_200": float(100 * sizes.lt(200).mean()),
        "above_500": int(sizes.gt(500).sum()),
    }


def plot_histogram(audit: pd.DataFrame, summary: dict, output: Path) -> None:
    """Show the common range and the complete tail without dropping large experiments."""
    sizes = audit.loc[audit["size_status"].eq("included"), "total_n"].astype(int)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.6), gridspec_kw={"width_ratios": [1.3, 1]})
    fig.subplots_adjust(left=0.075, right=0.975, top=0.70, bottom=0.28, wspace=0.24)
    fig.suptitle(
        "How large are the oral BugSigDB experiments?",
        x=0.075,
        ha="left",
        y=0.96,
        fontsize=19,
        fontweight="bold",
        color="#172b4d",
    )
    fig.text(
        0.075,
        0.86,
        f"{summary['experiments_with_sizes']:,} experiments with reported group sizes "
        f"from {summary['studies_with_sizes']:,} studies  |  "
        f"Median n = {summary['median_n']:g}\n"
        f"200 exceeds {summary['percent_below_200']:.1f}% of reported experiment totals",
        fontsize=11,
        linespacing=1.7,
        color="#334155",
        va="top",
    )
    axes[0].hist(
        sizes[sizes <= 500],
        bins=np.arange(0, 501, 25),
        color="#218a8a",
        edgecolor="white",
        linewidth=0.7,
    )
    axes[0].set_xlim(0, 500)
    axes[0].set_title(
        f"Detail: n ≤ 500 ({summary['above_500']} larger experiments at right)",
        fontsize=10,
        loc="left",
        pad=12,
    )
    low, high = min(2, int(sizes.min())), max(250, int(sizes.max())) * 1.01
    axes[1].hist(
        sizes, bins=np.geomspace(low, high, 24), color="#218a8a", edgecolor="white", linewidth=0.7
    )
    axes[1].set_xscale("log")
    axes[1].set_xlim(low, high)
    ticks = [n for n in [2, 10, 50, 200, 1000, 5000, 20000, 100000] if low <= n <= high]
    axes[1].set_xticks(ticks)
    axes[1].xaxis.set_major_formatter(ScalarFormatter())
    axes[1].set_title("All reported totals (logarithmic x-axis)", fontsize=10, loc="left", pad=12)
    for ax in axes:
        ax.axvline(200, color="#cb2535", linewidth=2.3, zorder=5)
        ax.set_xlabel("Reported samples per experiment: group 0 + group 1", fontsize=10)
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_axisbelow(True)
        ax.grid(axis="y", color="#e2e8f0", linewidth=0.6)
    axes[0].set_ylabel("Number of experiments")
    fig.text(
        0.075,
        0.15,
        "Red line: reference study with 200 samples (100 disease + 100 healthy)",
        color="#b51f2e",
        fontsize=11,
        fontweight="bold",
    )
    fig.text(
        0.075,
        0.058,
        f"Oral heatmap cohort: {summary['oral_experiments']} experiments; "
        f"{summary['experiments_excluded']} lack usable group-size totals. "
        "Each experiment counted once.\n"
        "Reported group-size sums may reuse participants within or across experiments. "
        "Reference is not added to the histogram.",
        fontsize=9,
        color="#475569",
        linespacing=1.6,
    )
    with plt.rc_context({"svg.fonttype": "none", "svg.hashsalt": "oral-experiment-sizes"}):
        for suffix in ["png", "pdf", "svg"]:
            fig.savefig(
                output / f"oral_experiment_sample_sizes.{suffix}", dpi=200, facecolor="white"
            )
    plt.close(fig)


def main() -> None:
    """Generate an auditable oral experiment sample-size figure from local inputs."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dump", type=Path, default=ROOT / "data/raw/full_dump.csv")
    parser.add_argument(
        "--metadata", type=Path, default=VIS / "study_trait_profile_study_metadata.tsv"
    )
    parser.add_argument(
        "--contrasts",
        type=Path,
        default=SIX / "six_shared_trait_condition_body_site_direction_contrasts.tsv",
    )
    parser.add_argument("--output-dir", type=Path, default=OUTPUT)
    args = parser.parse_args()
    metadata = pd.read_csv(args.metadata, sep="\t", dtype=str, keep_default_na=False)
    cells = pd.read_csv(args.contrasts, sep="\t")
    raw = read_raw_dump(args.raw_dump)
    audit = sample_size_audit(oral_experiments(metadata, cells), raw)
    summary = summarize(audit)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    audit.to_csv(args.output_dir / "oral_experiment_sample_size_audit.tsv", sep="\t", index=False)
    pd.DataFrame([summary]).to_csv(
        args.output_dir / "oral_experiment_sample_size_summary.tsv", sep="\t", index=False
    )
    # Store provenance beside the audit, without embedding machine-specific absolute paths.
    inputs = []
    for role, path in [
        ("raw_dump", args.raw_dump),
        ("profiles", args.metadata),
        ("contrasts", args.contrasts),
    ]:
        inputs.append(
            {
                "role": role,
                "filename": path.name,
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
        )
    pd.DataFrame(inputs).to_csv(args.output_dir / "inputs.tsv", sep="\t", index=False)
    plot_histogram(audit, summary, args.output_dir)
    print(json.dumps(summary, indent=2))
    print(audit["size_status"].value_counts().to_string())


if __name__ == "__main__":
    main()
