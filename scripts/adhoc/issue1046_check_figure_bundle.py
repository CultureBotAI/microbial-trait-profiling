"""Audit exported figure hashes/coverage and optional final-size reports (#1046, 2026-10-28).

Inputs: completed figure bundle and optional existing design-check reports.
Outputs: a new evidence directory; sources are read-only. No scientific replay.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
from collections import Counter
from pathlib import Path
from typing import Any

import yaml
from PIL import Image

from mtp.workflow.output_roots import require_durable_output_path

VERSION = "publication-figures-v1"
FIGURES = {"forest", "quadrant", "three-test"}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def identity(raw: bytes) -> dict[str, Any]:
    return {"sha256": hashlib.sha256(raw).hexdigest(), "size_bytes": len(raw)}


def read_bound(
    path: Path,
    snapshots: dict[Path, dict[str, Any]],
    targets: dict[Path, Path],
    expected: Any = None,
) -> bytes:
    """Hash exactly the bytes passed to parsers and retain the resolved identity."""
    require(
        path.is_file() and not path.is_symlink(), "Audit inputs must be regular nonsymlink files"
    )
    target = path.resolve(strict=True)
    raw = path.read_bytes()
    actual = identity(raw)
    require(path.resolve(strict=True) == target, "Audit input target changed while reading")
    if expected is not None:
        require(actual == expected, "Audit input hash or size identity differs")
    if path in snapshots:
        require(
            actual == snapshots[path] and target == targets[path], "Audit input identity changed"
        )
    snapshots[path], targets[path] = actual, target
    return raw


def json_object(raw: bytes) -> dict[str, Any]:
    def unique(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            require(key not in result, "Duplicate audit manifest field")
            result[key] = value
        return result

    result = json.loads(raw, object_pairs_hook=unique)
    require(isinstance(result, dict), "Audit manifest must be an object")
    parsed: dict[str, Any] = result
    return parsed


def inventory(directory: Path, expected: set[str]) -> None:
    require(
        {path.name for path in directory.iterdir()} == expected, "Figure bundle inventory differs"
    )


def pages_and_exports(record: dict[str, Any], name: str) -> tuple[list[dict[str, Any]], set[str]]:
    require(
        record.get("version") == VERSION and record.get("figure") == name, "Figure identity differs"
    )
    commit = record.get("code_commit")
    require(
        isinstance(commit, str)
        and re.fullmatch(r"[0-9a-f]{40}", commit) is not None
        and type(record.get("code_dirty")) is bool,
        "Figure producer identity is malformed",
    )
    pages, rendering = record.get("pages"), record.get("rendering")
    if (
        not isinstance(pages, list)
        or not pages
        or not all(isinstance(page, dict) for page in pages)
    ):
        raise ValueError("Figure page inventory must be nonempty")
    require(
        all(type(page.get("page")) is int for page in pages)
        and [page["page"] for page in pages] == list(range(1, len(pages) + 1)),
        "Figure page inventory must be contiguous and unique",
    )
    require(
        isinstance(rendering, dict)
        and type(rendering.get("page_count")) is int
        and rendering["page_count"] == len(pages)
        and rendering.get("vector_formats") == ["pdf", "svg"]
        and rendering.get("raster_formats") == ["png"]
        and type(rendering.get("dpi")) is int
        and rendering["dpi"] >= 300,
        "Figure rendering/page inventory differs",
    )
    for page in pages:
        require(
            isinstance(page.get("row_ids"), list)
            and all(isinstance(row, str) and bool(row) for row in page["row_ids"])
            and type(page.get("mark_count")) is int
            and page["mark_count"] >= 0
            and isinstance(page.get("details"), dict),
            "Figure page coverage metadata is malformed",
        )
    exports = {f"{name}.pdf"} | {
        f"page-{page['page']:03d}.{suffix}" for page in pages for suffix in ("pdf", "svg", "png")
    }
    require(
        isinstance(record.get("outputs"), dict) and set(record["outputs"]) == exports,
        "Figure export inventory differs from its pages",
    )
    for item in record["outputs"].values():
        require(
            isinstance(item, dict)
            and set(item) == {"sha256", "size_bytes"}
            and isinstance(item["sha256"], str)
            and re.fullmatch(r"[0-9a-f]{64}", item["sha256"]) is not None
            and type(item["size_bytes"]) is int
            and item["size_bytes"] >= 0,
            "Malformed figure export identity",
        )
    return pages, exports


def protect_table_sources(output: Path, aliases: dict[Path, Path]) -> None:
    """Keep #1206 audit evidence outside source bundles and aliases (#1208)."""
    parents = {alias.parent.resolve(strict=True) for alias in aliases} | {
        target.parent for target in aliases.values()
    }
    require(
        all(
            not output.is_relative_to(parent) and not parent.is_relative_to(output)
            for parent in parents
        ),
        "Audit output must be disjoint from parent trait evidence source directories",
    )


def expected_table_traits(
    record: dict[str, Any],
    snapshots: dict[Path, dict[str, Any]],
    targets: dict[Path, Path],
    output: Path,
) -> set[str]:
    """Read the bound parent family, independently of the page coverage declarations."""
    sources = record.get("input_files")
    if not isinstance(sources, list):
        raise ValueError("Table requires bound parent trait evidence")
    evidence = [
        item
        for item in sources
        if isinstance(item, dict)
        and isinstance(item.get("path"), str)
        and Path(item["path"]).name == "per_trait_evidence.tsv"
    ]
    require(bool(evidence), "Table requires bound parent trait evidence")
    aliases = {Path(item["path"]): Path(item["path"]).resolve(strict=True) for item in evidence}
    require(len(set(aliases.values())) == 1, "Ambiguous distinct parent trait evidence tables")
    protect_table_sources(output, aliases)
    contents: list[bytes] = []
    for item in evidence:
        require(
            set(item) == {"path", "sha256", "size_bytes"},
            "Malformed bound parent trait evidence identity",
        )
        contents.append(
            read_bound(
                aliases[Path(item["path"])],
                snapshots,
                targets,
                {key: item[key] for key in ("sha256", "size_bytes")},
            )
        )
    for alias, target in aliases.items():
        require(alias.resolve(strict=True) == target, "Parent trait evidence target changed")
        require(
            alias not in targets or targets[alias] == target, "Parent trait evidence target changed"
        )
        targets[alias] = target
    require(all(raw == contents[0] for raw in contents), "Conflicting parent trait evidence")
    rows = list(csv.reader(io.StringIO(contents[0].decode("utf-8")), delimiter="\t"))
    require(
        bool(rows) and len(rows[0]) == len(set(rows[0])) and "trait" in rows[0],
        "Malformed parent trait evidence columns",
    )
    column = rows[0].index("trait")
    require(
        all(len(row) == len(rows[0]) and bool(row[column]) for row in rows[1:]),
        "Malformed parent trait evidence rows",
    )
    traits = [row[column] for row in rows[1:]]
    require(bool(traits) and len(traits) == len(set(traits)), "Invalid parent trait family")
    return set(traits)


def table_coverage(pages: list[dict[str, Any]], traits: set[str]) -> None:
    ids = [row for page in pages for row in page["row_ids"]]
    require(
        set(ids) == traits and len(ids) == len(traits), "Table does not cover the full trait family"
    )
    for page in pages:
        details = page["details"]
        require(
            type(details.get("retained_family_count")) is int
            and details["retained_family_count"] == len(traits)
            and type(details.get("page_row_count")) is int
            and details["page_row_count"] == len(page["row_ids"])
            and page["mark_count"] == len(page["row_ids"]),
            "Table page coverage differs from the bound trait family",
        )


def quadrant_coverage(pages: list[dict[str, Any]]) -> None:
    for page in pages:
        populations = page["details"].get("populations")
        require(
            isinstance(populations, dict) and set(populations) == {"all", "stool"},
            "Quadrant population coverage is incomplete",
        )
        counts: dict[str, int] = {}
        for population, detail in populations.items():
            require(
                isinstance(detail, dict)
                and type(detail.get("n_comparisons")) is int
                and detail["n_comparisons"] >= 0,
                "Quadrant population comparison count is invalid",
            )
            counts[population] = detail["n_comparisons"]
        observed = Counter(row.partition(":")[0] for row in page["row_ids"])
        require(
            page["mark_count"] == sum(counts.values())
            and len(page["row_ids"]) == sum(counts.values())
            and all(observed[population] == count for population, count in counts.items())
            and set(observed).issubset(counts),
            "Quadrant marks and identities do not cover every population comparison",
        )


def audit(bundle: Path, output: Path, design_reports: Path | None) -> dict[str, Any]:
    if output.exists() or output.is_symlink():
        raise FileExistsError("Choose a fresh figure-check evidence directory")
    output = require_durable_output_path(output, allowed=False)
    if output.is_relative_to(bundle.resolve()) or bundle.resolve().is_relative_to(output):
        raise ValueError("Check reports must be disjoint from the figure bundle")
    if design_reports is not None:
        reports_root = design_reports.resolve(strict=True)
        require(
            not output.is_relative_to(reports_root) and not reports_root.is_relative_to(output),
            "Check output must be disjoint from existing design reports",
        )
    snapshots: dict[Path, dict[str, Any]] = {}
    targets: dict[Path, Path] = {}
    memberships: dict[Path, set[str]] = {}
    root_path = bundle / "figures.json"
    root_record = json_object(read_bound(root_path, snapshots, targets))
    require(
        root_record.get("version") == VERSION and root_record.get("status") == "completed",
        "Incomplete or unsupported figure bundle",
    )
    entries = root_record.get("figures")
    if not isinstance(entries, dict) or not entries or not set(entries).issubset(FIGURES):
        raise ValueError("Figure inventory must be a nonempty supported selection")
    memberships[bundle] = {"figures.json", *entries}
    inventory(bundle, memberships[bundle])
    figures: dict[str, Any] = {}
    report_identities: dict[str, dict[str, Any]] = {}
    for name, entry in entries.items():
        require(
            isinstance(entry, dict)
            and set(entry) == {"path", "sha256", "size_bytes"}
            and entry["path"] == f"{name}/provenance.json",
            "Figure provenance path differs",
        )
        directory = bundle / name
        require(
            directory.is_dir() and not directory.is_symlink(),
            "Figure directory must not be a symlink",
        )
        record_path = bundle / entry["path"]
        record = json_object(
            read_bound(
                record_path,
                snapshots,
                targets,
                {key: entry[key] for key in ("sha256", "size_bytes")},
            )
        )
        pages, exports = pages_and_exports(record, name)
        memberships[directory] = {"provenance.json", *exports}
        inventory(directory, memberships[directory])
        for filename, expected in record["outputs"].items():
            read_bound(directory / filename, snapshots, targets, expected)
        all_ids = [row for page in pages for row in page["row_ids"]]
        if len(all_ids) != len(set(all_ids)):
            raise ValueError("Figure page coverage has duplicate identities")
        if name == "three-test":
            traits = expected_table_traits(record, snapshots, targets, output)
            table_coverage(pages, traits)
        elif name == "quadrant":
            quadrant_coverage(pages)
        verdicts: Counter[str] = Counter()
        for page in pages:
            prefix = record_path.parent / f"page-{page['page']:03d}"
            png = read_bound(prefix.with_suffix(".png"), snapshots, targets)
            with Image.open(io.BytesIO(png)) as image:
                if min(image.info["dpi"]) < record["rendering"]["dpi"] - 0.02:
                    raise ValueError("PNG resolution differs from provenance")
            if design_reports is not None:
                report = design_reports / f"{name}-page-{page['page']:03d}.yaml"
                checked = yaml.safe_load(read_bound(report, snapshots, targets))
                require(
                    isinstance(checked, dict)
                    and checked.get("verdict") in {"pass", "pass_with_warnings", "fail"},
                    "Malformed design-check report",
                )
                actual = {item["format"]: item["sha256"] for item in checked["inputs"]}
                expected = {
                    suffix: snapshots[prefix.with_suffix(f".{suffix}")]["sha256"]
                    for suffix in ("pdf", "png")
                }
                if len(checked["inputs"]) != 2 or actual != expected:
                    raise ValueError("Design-check report does not bind the exported page")
                verdicts[checked["verdict"]] += 1
                report_identities[report.name] = snapshots[report]
        figures[name] = {
            "rows": len(all_ids),
            "marks": sum(page["mark_count"] for page in record["pages"]),
            "pages": len(record["pages"]),
            "design_check_verdicts": dict(verdicts),
            "provenance": snapshots[record_path],
            "producer_commit": record["code_commit"],
            "producer_dirty": record["code_dirty"],
        }
    for directory, expected_members in memberships.items():
        inventory(directory, expected_members)
    for path in snapshots:
        read_bound(path, snapshots, targets)
    for alias, target in targets.items():
        require(alias.resolve(strict=True) == target, "Audit input target changed")
    result = {
        "figures": figures,
        "bundle_manifest": snapshots[root_path],
        "design_reports": report_identities,
    }
    output.mkdir(parents=True)
    with (output / "audit.json").open("x", encoding="utf-8") as handle:
        handle.write(json.dumps(result, indent=2) + "\n")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--design-reports", type=Path)
    args = parser.parse_args()
    record = audit(args.bundle, args.output, args.design_reports)
    print(json.dumps(record, indent=2))
    return int(
        any(item["design_check_verdicts"].get("fail", 0) for item in record["figures"].values())
    )


if __name__ == "__main__":
    raise SystemExit(main())
