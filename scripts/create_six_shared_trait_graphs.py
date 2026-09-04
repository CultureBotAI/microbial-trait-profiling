"""Create Cytoscape.js graph versions of the six shared trait Sankey outputs."""
# ruff: noqa: E501

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = REPO_ROOT / "trait-visualization" / "trait-profile-investigation" / "six-shared-traits"

HTML_FILENAME = "six_shared_traits_profile_graphs.html"
DATA_FILENAME = "six_shared_traits_profile_graph_data.json"
MANIFEST_FILENAME = "six_shared_trait_graph_manifest.tsv"
CYTOSCAPE_CDN = "https://cdn.jsdelivr.net/npm/cytoscape/dist/cytoscape.min.js"


@dataclass(frozen=True)
class GraphView:
    """One Cytoscape graph view derived from a six-shared-traits Sankey TSV."""

    key: str
    label: str
    focus_direction: str
    condition_mode: str
    display_edges_file: str


GRAPH_VIEWS = (
    GraphView(
        key="increased_grouped",
        label="Increased, grouped conditions",
        focus_direction="increased",
        condition_mode="grouped condition labels",
        display_edges_file="six_shared_traits_increased_profile_sankey_display_edges.tsv",
    ),
    GraphView(
        key="decreased_grouped",
        label="Decreased, grouped conditions",
        focus_direction="decreased",
        condition_mode="grouped condition labels",
        display_edges_file="six_shared_traits_decreased_profile_sankey_display_edges.tsv",
    ),
    GraphView(
        key="increased_individual_conditions",
        label="Increased, individual conditions",
        focus_direction="increased",
        condition_mode="raw individual condition labels",
        display_edges_file=(
            "six_shared_traits_increased_individual_conditions_profile_sankey_display_edges.tsv"
        ),
    ),
    GraphView(
        key="decreased_individual_conditions",
        label="Decreased, individual conditions",
        focus_direction="decreased",
        condition_mode="raw individual condition labels",
        display_edges_file=(
            "six_shared_traits_decreased_individual_conditions_profile_sankey_display_edges.tsv"
        ),
    ),
)


def read_display_edges(path: Path) -> pd.DataFrame:
    """Read one display-edge TSV and coerce graph-relevant numeric fields."""
    frame = pd.read_csv(path, sep="\t")
    numeric_columns = [
        "trait_rank",
        "focus_enriched_profile_count",
        "focus_enriched_study_count",
        "background_no_enriched_profile_count",
        "focus_enriched_mean_guild_z",
        "focus_enriched_minus_comparator_guild_z",
        "rank_score",
    ]
    for column in numeric_columns:
        frame[column] = pd.to_numeric(frame[column], errors="coerce").fillna(0.0)
    return frame


def value_for_color(row: pd.Series) -> float:
    """Return the nonnegative signal value used for edge color intensity."""
    contrast = row.get("focus_enriched_minus_comparator_guild_z")
    mean = row.get("focus_enriched_mean_guild_z")
    value = contrast if pd.notna(contrast) else mean
    if pd.isna(value):
        return 0.0
    return max(0.0, float(value))


def short_label(value: object, *, max_length: int = 34) -> str:
    """Return a compact label for graph nodes."""
    text = str(value or "").strip() or "unknown"
    if len(text) <= max_length:
        return text
    return f"{text[: max_length - 3].rstrip()}..."


def stable_id(view_key: str, kind: str, *parts: object) -> str:
    """Return a compact, deterministic Cytoscape element id."""
    raw = "\x1f".join(str(part) for part in parts)
    digest = hashlib.sha1(raw.encode("utf-8")).hexdigest()[:14]
    return f"{view_key}-{kind}-{digest}"


def joined_unique_count(values: pd.Series) -> int:
    """Count unique semicolon-delimited identifiers from summarized edge rows."""
    unique_values: set[str] = set()
    for value in values.dropna():
        for token in str(value).split(";"):
            token = token.strip()
            if token:
                unique_values.add(token)
    return len(unique_values)


def node_stats(frame: pd.DataFrame, column: str) -> pd.DataFrame:
    """Summarize graph node support across displayed trait-specific links."""
    records: list[dict[str, object]] = []
    for label, group in frame.groupby(column, sort=False):
        records.append(
            {
                column: label,
                "links": len(group),
                "traits": group["trait_id"].nunique(),
                "n": int(group["focus_enriched_profile_count"].sum()),
                "s": joined_unique_count(group["focus_enriched_study_ids"]),
                "rank_score": float(group["rank_score"].max()),
            }
        )
    return pd.DataFrame.from_records(records).set_index(column)


def hex_to_rgb(color: str) -> tuple[int, int, int]:
    """Convert #rrggbb into RGB integers."""
    color = color.strip().lstrip("#")
    return tuple(int(color[index : index + 2], 16) for index in (0, 2, 4))


def rgb_to_hex(rgb: tuple[float, float, float]) -> str:
    """Convert RGB floats to #rrggbb."""
    return "#" + "".join(f"{max(0, min(255, round(channel))):02x}" for channel in rgb)


def lighten_color(color: str, fraction: float = 0.18) -> str:
    """Return a light tint from the same hue family."""
    rgb = hex_to_rgb(color)
    mixed = tuple(255 - ((255 - channel) * fraction) for channel in rgb)
    return rgb_to_hex(mixed)


def color_from_signal(base_color: str, signal: float, shared_signal_max: float) -> str:
    """Interpolate from a light tint to the base trait color."""
    light = hex_to_rgb(lighten_color(base_color))
    base = hex_to_rgb(base_color)
    weight = 0.0 if shared_signal_max <= 0 else min(1.0, signal / shared_signal_max)
    mixed = tuple((light[i] * (1.0 - weight)) + (base[i] * weight) for i in range(3))
    return rgb_to_hex(mixed)


def parallel_offsets(count: int, *, spacing: float = 14.0, max_offset: float = 42.0) -> list[float]:
    """Return top-to-bottom offsets for body-site trait bundle nodes."""
    if count <= 1:
        return [0.0]
    center = (count - 1) / 2
    return [
        min(max_offset, max(-max_offset, (center - position) * spacing))
        for position in range(count)
    ]


def weighted_signal(group: pd.DataFrame) -> float:
    """Return count-weighted signal for a body-site trait trunk."""
    values = group["color_value"].astype(float)
    weights = group["focus_enriched_profile_count"].astype(float)
    if float(weights.sum()) <= 0:
        return float(values.mean()) if not values.empty else 0.0
    return float((values * weights).sum() / weights.sum())


def position_labels(
    labels: list[str],
    *,
    step: float,
    shared_span: float,
    top: float = 80.0,
) -> dict[str, float]:
    """Return centered y positions for one graph column."""
    if not labels:
        return {}
    local_span = max(0.0, (len(labels) - 1) * step)
    start = top + ((shared_span - local_span) / 2)
    return {label: start + (index * step) for index, label in enumerate(labels)}


def ordered_index(stats: pd.DataFrame) -> list[str]:
    """Return node labels ordered by trait breadth and support."""
    return (
        stats.sort_values(
            ["traits", "n", "links", "rank_score"],
            ascending=[False, False, False, False],
        )
        .index.astype(str)
        .tolist()
    )


def trait_legend_records(frames: list[pd.DataFrame]) -> list[dict[str, object]]:
    """Return one legend row per trait in display order."""
    all_edges = pd.concat(frames, ignore_index=True)
    legend = (
        all_edges[["trait_id", "trait_label", "trait_slug", "trait_color", "trait_rank"]]
        .drop_duplicates("trait_id")
        .sort_values("trait_rank")
    )
    return [
        {
            "traitId": row.trait_id,
            "label": row.trait_label,
            "slug": row.trait_slug,
            "color": row.trait_color,
            "rank": int(row.trait_rank),
            "lightColor": lighten_color(row.trait_color),
        }
        for row in legend.itertuples(index=False)
    ]


def graph_ranges(source_counts: pd.Series, trunk_counts: pd.Series) -> dict[str, float]:
    """Return min/max ranges for Cytoscape mapData styling."""
    source_positive = source_counts[source_counts > 0]
    trunk_positive = trunk_counts[trunk_counts > 0]
    source_min = float(source_positive.min()) if not source_positive.empty else 1.0
    source_max = float(source_positive.max()) if not source_positive.empty else source_min
    trunk_min = float(trunk_positive.min()) if not trunk_positive.empty else 1.0
    trunk_max = float(trunk_positive.max()) if not trunk_positive.empty else trunk_min
    if source_min == source_max:
        source_max = source_min + 1.0
    if trunk_min == trunk_max:
        trunk_max = trunk_min + 1.0
    return {
        "sourceNMin": source_min,
        "sourceNMax": source_max,
        "trunkNMin": trunk_min,
        "trunkNMax": trunk_max,
    }


def build_view_payload(
    view: GraphView,
    frame: pd.DataFrame,
    *,
    shared_signal_max: float,
) -> tuple[dict[str, object], dict[str, object]]:
    """Build Cytoscape elements and manifest metadata for one graph view."""
    frame = frame.copy()
    frame["color_value"] = frame.apply(value_for_color, axis=1)
    condition_stats = node_stats(frame, "condition")
    body_site_stats = node_stats(frame, "body_site_group")
    condition_order = ordered_index(condition_stats)
    body_site_order = ordered_index(body_site_stats)

    condition_step = 58.0
    body_step = 68.0
    shared_span = max(
        max(0.0, (len(condition_order) - 1) * condition_step),
        max(0.0, (len(body_site_order) - 1) * body_step),
    )
    condition_y = position_labels(
        condition_order,
        step=condition_step,
        shared_span=shared_span,
    )
    body_site_y = position_labels(
        body_site_order,
        step=body_step,
        shared_span=shared_span,
    )

    bundle_records: list[dict[str, object]] = []
    for (body_site, trait_id), group in frame.groupby(["body_site_group", "trait_id"]):
        base_color = str(group["trait_color"].iloc[0])
        signal = weighted_signal(group)
        bundle_records.append(
            {
                "body_site_group": body_site,
                "trait_id": trait_id,
                "trait_label": str(group["trait_label"].iloc[0]),
                "trait_color": base_color,
                "trait_rank": int(group["trait_rank"].iloc[0]),
                "n": int(group["focus_enriched_profile_count"].sum()),
                "bg": int(group["background_no_enriched_profile_count"].sum()),
                "s": joined_unique_count(group["focus_enriched_study_ids"]),
                "links": len(group),
                "signal": signal,
                "color": color_from_signal(base_color, signal, shared_signal_max),
                "rank_score": float(group["rank_score"].max()),
            }
        )
    bundles = pd.DataFrame.from_records(bundle_records)
    bundles["bundle_rank"] = 1
    bundles["bundle_offset"] = 0.0
    for _body_site, group in bundles.groupby("body_site_group", sort=False):
        ordered = (
            group.sort_values(
                ["n", "rank_score", "trait_rank"],
                ascending=[False, False, True],
            )
            .index.astype(int)
            .tolist()
        )
        for rank, (index, offset) in enumerate(
            zip(ordered, parallel_offsets(len(ordered)), strict=True),
            start=1,
        ):
            bundles.loc[index, "bundle_rank"] = rank
            bundles.loc[index, "bundle_offset"] = offset

    condition_ids = {label: stable_id(view.key, "condition", label) for label in condition_order}
    body_site_ids = {label: stable_id(view.key, "body", label) for label in body_site_order}
    bundle_ids = {
        (row.body_site_group, row.trait_id): stable_id(
            view.key,
            "bundle",
            row.body_site_group,
            row.trait_id,
        )
        for row in bundles.itertuples(index=False)
    }

    elements: list[dict[str, object]] = []
    for label in condition_order:
        stats = condition_stats.loc[label]
        y = condition_y[label]
        elements.append(
            {
                "data": {
                    "id": condition_ids[label],
                    "nodeType": "condition",
                    "label": label,
                    "fullLabel": label,
                    "displayLabel": (
                        f"{short_label(label)}\n"
                        f"links={int(stats.links)}, traits={int(stats.traits)}, "
                        f"n={int(stats.n):,}"
                    ),
                    "links": int(stats.links),
                    "traits": int(stats.traits),
                    "n": int(stats.n),
                    "s": int(stats.s),
                    "presetX": 120,
                    "presetY": y,
                },
                "position": {"x": 120, "y": y},
            }
        )
    for row in bundles.itertuples(index=False):
        y = body_site_y[row.body_site_group] + float(row.bundle_offset)
        x = 760
        elements.append(
            {
                "data": {
                    "id": bundle_ids[(row.body_site_group, row.trait_id)],
                    "nodeType": "bundle",
                    "label": f"{row.trait_label} -> {row.body_site_group}",
                    "fullLabel": f"{row.trait_label} into {row.body_site_group}",
                    "displayLabel": "",
                    "bodySite": row.body_site_group,
                    "traitId": row.trait_id,
                    "traitLabel": row.trait_label,
                    "color": row.color,
                    "baseColor": row.trait_color,
                    "n": int(row.n),
                    "bg": int(row.bg),
                    "s": int(row.s),
                    "links": int(row.links),
                    "signal": round(float(row.signal), 4),
                    "presetX": x,
                    "presetY": y,
                },
                "position": {"x": x, "y": y},
            }
        )
    for label in body_site_order:
        stats = body_site_stats.loc[label]
        y = body_site_y[label]
        elements.append(
            {
                "data": {
                    "id": body_site_ids[label],
                    "nodeType": "body_site",
                    "label": label,
                    "fullLabel": label,
                    "displayLabel": (
                        f"{short_label(label)}\n"
                        f"links={int(stats.links)}, traits={int(stats.traits)}, "
                        f"n={int(stats.n):,}"
                    ),
                    "links": int(stats.links),
                    "traits": int(stats.traits),
                    "n": int(stats.n),
                    "s": int(stats.s),
                    "presetX": 1010,
                    "presetY": y,
                },
                "position": {"x": 1010, "y": y},
            }
        )

    for index, row in enumerate(frame.itertuples(index=False), start=1):
        base_color = str(row.trait_color)
        signal = float(row.color_value)
        source = condition_ids[str(row.condition)]
        target = bundle_ids[(row.body_site_group, row.trait_id)]
        elements.append(
            {
                "data": {
                    "id": stable_id(
                        view.key,
                        "source-edge",
                        index,
                        row.condition,
                        row.body_site_group,
                        row.trait_id,
                    ),
                    "edgeType": "source",
                    "source": source,
                    "target": target,
                    "label": f"{row.trait_label}: n={int(row.focus_enriched_profile_count)}",
                    "condition": row.condition,
                    "bodySite": row.body_site_group,
                    "traitId": row.trait_id,
                    "traitLabel": row.trait_label,
                    "baseColor": base_color,
                    "color": color_from_signal(base_color, signal, shared_signal_max),
                    "n": int(row.focus_enriched_profile_count),
                    "bg": int(row.background_no_enriched_profile_count),
                    "s": int(row.focus_enriched_study_count),
                    "signal": round(signal, 4),
                    "contrast": round(
                        float(row.focus_enriched_minus_comparator_guild_z),
                        4,
                    ),
                    "rankScore": round(float(row.rank_score), 4),
                }
            }
        )

    for row in bundles.itertuples(index=False):
        elements.append(
            {
                "data": {
                    "id": stable_id(
                        view.key,
                        "trunk-edge",
                        row.body_site_group,
                        row.trait_id,
                    ),
                    "edgeType": "trunk",
                    "source": bundle_ids[(row.body_site_group, row.trait_id)],
                    "target": body_site_ids[row.body_site_group],
                    "label": f"{row.trait_label}: summed n={int(row.n)}",
                    "bodySite": row.body_site_group,
                    "traitId": row.trait_id,
                    "traitLabel": row.trait_label,
                    "baseColor": row.trait_color,
                    "color": row.color,
                    "n": int(row.n),
                    "bg": int(row.bg),
                    "s": int(row.s),
                    "signal": round(float(row.signal), 4),
                    "sourceLinks": int(row.links),
                }
            }
        )

    ranges = graph_ranges(
        frame["focus_enriched_profile_count"],
        bundles["n"],
    )
    ranges["signalMax"] = round(shared_signal_max, 4)

    payload = {
        "key": view.key,
        "label": view.label,
        "focusDirection": view.focus_direction,
        "conditionMode": view.condition_mode,
        "sourceFile": view.display_edges_file,
        "elements": elements,
        "ranges": ranges,
        "stats": {
            "conditions": len(condition_order),
            "bodySites": len(body_site_order),
            "traitBodyBundles": len(bundles),
            "sourceEdges": len(frame),
            "trunkEdges": len(bundles),
            "totalNodes": len([element for element in elements if "position" in element]),
            "totalEdges": len([element for element in elements if "position" not in element]),
            "displayedTraits": frame["trait_id"].nunique(),
        },
    }
    manifest = {
        "view_key": view.key,
        "view_label": view.label,
        "focus_direction": view.focus_direction,
        "condition_mode": view.condition_mode,
        "source_display_edges_path": (OUTPUT_DIR / view.display_edges_file)
        .relative_to(REPO_ROOT)
        .as_posix(),
        "condition_node_count": len(condition_order),
        "body_site_node_count": len(body_site_order),
        "trait_body_bundle_count": len(bundles),
        "source_edge_count": len(frame),
        "trunk_edge_count": len(bundles),
        "html_path": (OUTPUT_DIR / HTML_FILENAME).relative_to(REPO_ROOT).as_posix(),
        "data_path": (OUTPUT_DIR / DATA_FILENAME).relative_to(REPO_ROOT).as_posix(),
    }
    return payload, manifest


def build_graph_payload(output_dir: Path = OUTPUT_DIR) -> tuple[dict[str, object], pd.DataFrame]:
    """Build the complete multi-view graph payload and manifest table."""
    frames = [read_display_edges(output_dir / view.display_edges_file) for view in GRAPH_VIEWS]
    all_signals = pd.concat(
        [frame.apply(value_for_color, axis=1) for frame in frames],
        ignore_index=True,
    )
    shared_signal_max = max(0.4, float(all_signals.max() if not all_signals.empty else 0.0))

    views: dict[str, object] = {}
    manifest_records: list[dict[str, object]] = []
    for view, frame in zip(GRAPH_VIEWS, frames, strict=True):
        payload, manifest = build_view_payload(
            view,
            frame,
            shared_signal_max=shared_signal_max,
        )
        views[view.key] = payload
        manifest_records.append(manifest)

    payload = {
        "title": "Six shared traits Cytoscape graph explorer",
        "description": (
            "Graph versions of the all-six-trait grouped and individual-condition "
            "Sankey maps. Condition links fan into body-site trait bundles; each "
            "body-site node has at most six incoming trait trunks."
        ),
        "cytoscapePackage": CYTOSCAPE_CDN,
        "sharedSignalMax": round(shared_signal_max, 4),
        "traits": trait_legend_records(frames),
        "viewOrder": [view.key for view in GRAPH_VIEWS],
        "views": views,
    }
    return payload, pd.DataFrame.from_records(manifest_records)


def html_template() -> str:
    """Return the static Cytoscape graph explorer template."""
    return r"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>Six shared traits Cytoscape graph explorer</title>
  <script src="__CYTOSCAPE_CDN__"></script>
  <style>
    :root {
      --ink: #0f172a;
      --muted: #64748b;
      --line: #cbd5e1;
      --panel: #f8fafc;
      --accent: #334155;
    }
    * { box-sizing: border-box; }
    body {
      margin: 0;
      font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      color: var(--ink);
      background: #ffffff;
    }
    .topbar {
      display: grid;
      grid-template-columns: minmax(320px, 1fr) auto;
      gap: 18px;
      align-items: end;
      padding: 16px 18px 12px;
      border-bottom: 1px solid var(--line);
      background: #ffffff;
    }
    h1 {
      margin: 0 0 4px;
      font-size: 20px;
      line-height: 1.2;
      letter-spacing: 0;
    }
    .subtitle {
      margin: 0;
      max-width: 1120px;
      color: var(--muted);
      font-size: 13px;
      line-height: 1.35;
    }
    .controls {
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
      justify-content: flex-end;
      align-items: center;
    }
    select, input[type="search"], button {
      height: 34px;
      border: 1px solid #94a3b8;
      background: #ffffff;
      color: var(--ink);
      font: inherit;
      font-size: 13px;
      border-radius: 6px;
    }
    select, input[type="search"] { padding: 0 10px; }
    input[type="search"] { width: 220px; }
    button {
      padding: 0 11px;
      cursor: pointer;
    }
    button:hover { background: var(--panel); }
    .workspace {
      display: grid;
      grid-template-columns: minmax(680px, 1fr) 320px;
      min-height: calc(100vh - 86px);
    }
    #cy {
      width: 100%;
      height: calc(100vh - 86px);
      min-height: 720px;
      background: #ffffff;
    }
    aside {
      border-left: 1px solid var(--line);
      background: var(--panel);
      padding: 14px;
      overflow: auto;
      max-height: calc(100vh - 86px);
    }
    .section {
      border-bottom: 1px solid #e2e8f0;
      padding: 0 0 13px;
      margin-bottom: 14px;
    }
    .section:last-child { border-bottom: 0; }
    .section h2 {
      margin: 0 0 8px;
      font-size: 13px;
      line-height: 1.2;
      color: var(--accent);
      letter-spacing: 0;
    }
    .stats {
      display: grid;
      grid-template-columns: 1fr 1fr;
      gap: 6px;
    }
    .stat {
      padding: 7px 8px;
      border: 1px solid #e2e8f0;
      border-radius: 6px;
      background: #ffffff;
    }
    .stat b {
      display: block;
      font-size: 15px;
      line-height: 1.1;
    }
    .stat span {
      color: var(--muted);
      font-size: 11px;
    }
    .trait-list {
      display: grid;
      gap: 5px;
    }
    .trait-filter {
      display: grid;
      grid-template-columns: 18px 14px 1fr;
      gap: 7px;
      align-items: center;
      font-size: 12px;
    }
    .swatch {
      width: 14px;
      height: 14px;
      border-radius: 3px;
      border: 1px solid rgba(15, 23, 42, 0.25);
    }
    .range-row {
      display: grid;
      grid-template-columns: 1fr auto;
      gap: 8px;
      align-items: center;
      color: var(--muted);
      font-size: 12px;
    }
    input[type="range"] { width: 100%; }
    .legend-item {
      display: grid;
      grid-template-columns: 72px 1fr;
      gap: 9px;
      align-items: center;
      color: var(--muted);
      font-size: 12px;
      margin: 5px 0;
    }
    .line-sample {
      height: 0;
      border-top: solid #334155;
      border-radius: 99px;
    }
    .details {
      min-height: 140px;
      color: #1e293b;
      font-size: 12px;
      line-height: 1.45;
    }
    .details dl {
      display: grid;
      grid-template-columns: 90px 1fr;
      gap: 5px 8px;
      margin: 0;
    }
    .details dt {
      color: var(--muted);
    }
    .details dd {
      margin: 0;
      overflow-wrap: anywhere;
    }
    .note {
      color: var(--muted);
      font-size: 11px;
      line-height: 1.35;
      margin: 0;
    }
    .error {
      margin: 20px;
      padding: 12px;
      border: 1px solid #fca5a5;
      background: #fef2f2;
      color: #7f1d1d;
      border-radius: 6px;
      display: none;
    }
    @media (max-width: 980px) {
      .topbar { grid-template-columns: 1fr; }
      .controls { justify-content: flex-start; }
      .workspace { grid-template-columns: 1fr; }
      #cy { height: 70vh; min-height: 560px; }
      aside {
        max-height: none;
        border-left: 0;
        border-top: 1px solid var(--line);
      }
    }
  </style>
</head>
<body>
  <div id="load-error" class="error">
    Cytoscape.js did not load. This explorer uses the Cytoscape.js package from jsDelivr; open it with network access or vendor the package locally.
  </div>
  <header class="topbar">
    <div>
      <h1>Six shared traits Cytoscape graph explorer</h1>
      <p class="subtitle">
        Condition nodes connect into trait/body-site bundle nodes, then each bundle connects to its body-site node as one summed trunk. Edge width uses <code>n=</code> positive-z focus-profile support; color shade uses the shared positive z-score scale.
      </p>
    </div>
    <div class="controls">
      <select id="view-select" aria-label="Graph view"></select>
      <input id="search" type="search" placeholder="Search labels" aria-label="Search labels" />
      <button id="layout-layered" type="button">Layered</button>
      <button id="layout-force" type="button">Force</button>
      <button id="fit" type="button">Fit</button>
      <button id="download" type="button">PNG</button>
    </div>
  </header>
  <main class="workspace">
    <section id="cy" aria-label="Cytoscape graph"></section>
    <aside>
      <section class="section">
        <h2>View</h2>
        <div id="stats" class="stats"></div>
      </section>
      <section class="section">
        <h2>Traits</h2>
        <div id="trait-filters" class="trait-list"></div>
      </section>
      <section class="section">
        <h2>Minimum Source Edge n</h2>
        <div class="range-row">
          <input id="min-n" type="range" min="0" max="1" value="0" step="1" />
          <span id="min-n-value">0</span>
        </div>
      </section>
      <section class="section">
        <h2>Legend</h2>
        <div class="legend-item"><span class="line-sample" style="border-top-width: 2px;"></span><span>small source edge</span></div>
        <div class="legend-item"><span class="line-sample" style="border-top-width: 7px;"></span><span>large source edge</span></div>
        <div class="legend-item"><span class="line-sample" style="border-top-width: 12px;"></span><span>body-site trunk, summed per trait color</span></div>
        <p class="note">Bundle nodes are the small colored diamonds before each body-site node. Each body-site can receive up to six trunks, one for each trait color scale.</p>
      </section>
      <section class="section">
        <h2>Selected Element</h2>
        <div id="details" class="details">Select a node or edge to inspect counts and labels.</div>
      </section>
      <section class="section">
        <h2>Data</h2>
        <p class="note" id="data-note"></p>
      </section>
    </aside>
  </main>
  <script>
    const GRAPH_DATA = __GRAPH_DATA__;
    let cy = null;

    function escapeHtml(value) {
      return String(value ?? "")
        .replaceAll("&", "&amp;")
        .replaceAll("<", "&lt;")
        .replaceAll(">", "&gt;")
        .replaceAll('"', "&quot;")
        .replaceAll("'", "&#039;");
    }

    function currentView() {
      return GRAPH_DATA.views[document.getElementById("view-select").value];
    }

    function edgeWidthRange(view, edgeType) {
      if (edgeType === "trunk") {
        return [view.ranges.trunkNMin, view.ranges.trunkNMax, 2.2, 12.0];
      }
      return [view.ranges.sourceNMin, view.ranges.sourceNMax, 1.1, 7.0];
    }

    function graphStyle(view) {
      const sourceRange = edgeWidthRange(view, "source");
      const trunkRange = edgeWidthRange(view, "trunk");
      return [
        {
          selector: "node",
          style: {
            "label": "data(displayLabel)",
            "font-family": "Inter, system-ui, sans-serif",
            "font-size": 9,
            "line-height": 1.05,
            "text-wrap": "wrap",
            "text-max-width": 112,
            "text-valign": "center",
            "text-halign": "center",
            "color": "#0f172a",
            "background-color": "#ffffff",
            "border-color": "#475569",
            "border-width": 1.4,
            "width": "label",
            "height": "label",
            "padding": 9,
            "shape": "roundrectangle"
          }
        },
        {
          selector: 'node[nodeType="condition"]',
          style: {
            "background-color": "#f8fafc",
            "text-margin-x": -8
          }
        },
        {
          selector: 'node[nodeType="body_site"]',
          style: {
            "background-color": "#ffffff",
            "text-margin-x": 8
          }
        },
        {
          selector: 'node[nodeType="bundle"]',
          style: {
            "label": "",
            "shape": "diamond",
            "background-color": "data(color)",
            "border-color": "data(baseColor)",
            "border-width": 1,
            "width": `mapData(n, ${trunkRange[0]}, ${trunkRange[1]}, 11, 30)`,
            "height": `mapData(n, ${trunkRange[0]}, ${trunkRange[1]}, 11, 30)`,
            "padding": 0
          }
        },
        {
          selector: "edge",
          style: {
            "line-color": "data(color)",
            "target-arrow-color": "data(color)",
            "target-arrow-shape": "triangle",
            "arrow-scale": 0.55,
            "line-cap": "round",
            "curve-style": "bezier",
            "font-size": 8,
            "color": "#334155",
            "text-background-color": "#ffffff",
            "text-background-opacity": 0.86,
            "text-background-padding": 2
          }
        },
        {
          selector: 'edge[edgeType="source"]',
          style: {
            "width": `mapData(n, ${sourceRange[0]}, ${sourceRange[1]}, ${sourceRange[2]}, ${sourceRange[3]})`,
            "opacity": 0.48,
            "z-index": 1
          }
        },
        {
          selector: 'edge[edgeType="trunk"]',
          style: {
            "width": `mapData(n, ${trunkRange[0]}, ${trunkRange[1]}, ${trunkRange[2]}, ${trunkRange[3]})`,
            "opacity": 0.9,
            "curve-style": "straight",
            "z-index": 3
          }
        },
        {
          selector: ":selected",
          style: {
            "border-width": 3,
            "border-color": "#111827",
            "line-color": "#111827",
            "target-arrow-color": "#111827",
            "z-index": 20
          }
        },
        {
          selector: ".matched",
          style: {
            "border-width": 3,
            "border-color": "#111827",
            "background-blacken": -0.08
          }
        },
        {
          selector: ".hidden",
          style: {
            "display": "none"
          }
        }
      ];
    }

    function populateViewSelect() {
      const select = document.getElementById("view-select");
      select.innerHTML = GRAPH_DATA.viewOrder
        .map((key) => `<option value="${key}">${escapeHtml(GRAPH_DATA.views[key].label)}</option>`)
        .join("");
    }

    function populateTraitFilters() {
      const container = document.getElementById("trait-filters");
      container.innerHTML = GRAPH_DATA.traits.map((trait) => `
        <label class="trait-filter">
          <input type="checkbox" value="${escapeHtml(trait.traitId)}" checked />
          <span class="swatch" style="background:${escapeHtml(trait.color)}"></span>
          <span>${escapeHtml(trait.label)}</span>
        </label>
      `).join("");
      container.querySelectorAll("input").forEach((input) => {
        input.addEventListener("change", applyFilters);
      });
    }

    function updateStats() {
      const view = currentView();
      const visibleSourceEdges = cy
        ? cy.edges('[edgeType="source"]').filter((edge) => !edge.hasClass("hidden")).length
        : view.stats.sourceEdges;
      const visibleTrunks = cy
        ? cy.edges('[edgeType="trunk"]').filter((edge) => !edge.hasClass("hidden")).length
        : view.stats.trunkEdges;
      document.getElementById("stats").innerHTML = `
        <div class="stat"><b>${view.stats.conditions}</b><span>condition nodes</span></div>
        <div class="stat"><b>${view.stats.bodySites}</b><span>body-site nodes</span></div>
        <div class="stat"><b>${visibleSourceEdges}</b><span>visible source edges</span></div>
        <div class="stat"><b>${visibleTrunks}</b><span>visible trunks</span></div>
      `;
      document.getElementById("data-note").innerHTML = `
        Source TSV: <code>${escapeHtml(view.sourceFile)}</code><br />
        Focus direction: ${escapeHtml(view.focusDirection)}<br />
        Conditions: ${escapeHtml(view.conditionMode)}<br />
        Shared shade scale max: ${escapeHtml(GRAPH_DATA.sharedSignalMax)}
      `;
    }

    function selectedTraitIds() {
      return new Set(
        Array.from(document.querySelectorAll('#trait-filters input:checked')).map(
          (input) => input.value
        )
      );
    }

    function applyFilters() {
      if (!cy) return;
      const traits = selectedTraitIds();
      const minN = Number(document.getElementById("min-n").value);
      document.getElementById("min-n-value").textContent = String(minN);

      cy.batch(() => {
        cy.elements().removeClass("hidden");
        cy.edges().forEach((edge) => {
          const traitVisible = traits.has(edge.data("traitId"));
          const sourceTooSmall = edge.data("edgeType") === "source" && Number(edge.data("n")) < minN;
          if (!traitVisible || sourceTooSmall) {
            edge.addClass("hidden");
          }
        });
        cy.nodes('[nodeType="bundle"]').forEach((node) => {
          const hasVisibleEdge = node
            .connectedEdges()
            .filter((edge) => !edge.hasClass("hidden")).length > 0;
          if (!traits.has(node.data("traitId")) || !hasVisibleEdge) {
            node.addClass("hidden");
          }
        });
        cy.nodes('[nodeType="condition"]').forEach((node) => {
          const hasVisibleSource = node
            .connectedEdges('[edgeType="source"]')
            .filter((edge) => !edge.hasClass("hidden")).length > 0;
          if (!hasVisibleSource) node.addClass("hidden");
        });
        cy.nodes('[nodeType="body_site"]').forEach((node) => {
          const hasVisibleTrunk = node
            .connectedEdges('[edgeType="trunk"]')
            .filter((edge) => !edge.hasClass("hidden")).length > 0;
          if (!hasVisibleTrunk) node.addClass("hidden");
        });
      });
      applySearch();
      updateStats();
    }

    function applySearch() {
      if (!cy) return;
      const term = document.getElementById("search").value.trim().toLowerCase();
      cy.elements().removeClass("matched");
      if (term.length < 2) return;
      cy.nodes().forEach((node) => {
        const text = `${node.data("fullLabel") || ""} ${node.data("traitLabel") || ""}`.toLowerCase();
        if (text.includes(term)) node.addClass("matched");
      });
    }

    function resetLayeredLayout() {
      const view = currentView();
      cy.nodes().forEach((node) => {
        node.position({ x: node.data("presetX"), y: node.data("presetY") });
      });
      cy.fit(undefined, 42);
      updateStats(view);
    }

    function runForceLayout() {
      cy.layout({
        name: "cose",
        animate: true,
        fit: true,
        padding: 54,
        nodeRepulsion: 180000,
        idealEdgeLength: 96,
        edgeElasticity: 70,
        gravity: 0.18,
        numIter: 1200
      }).run();
    }

    function renderDetails(element) {
      const data = element.data();
      const isNode = element.isNode();
      const type = isNode ? data.nodeType : data.edgeType;
      const rows = [
        ["type", type],
        ["label", data.fullLabel || data.label],
      ];
      if (data.traitLabel) rows.push(["trait", data.traitLabel]);
      if (data.condition) rows.push(["condition", data.condition]);
      if (data.bodySite) rows.push(["body site", data.bodySite]);
      if (data.n !== undefined) rows.push(["n", Number(data.n).toLocaleString()]);
      if (data.s !== undefined) rows.push(["s", Number(data.s).toLocaleString()]);
      if (data.bg !== undefined) rows.push(["bg", Number(data.bg).toLocaleString()]);
      if (data.links !== undefined) rows.push(["links", Number(data.links).toLocaleString()]);
      if (data.sourceLinks !== undefined) rows.push(["source links", Number(data.sourceLinks).toLocaleString()]);
      if (data.signal !== undefined) rows.push(["z shade", data.signal]);
      if (data.contrast !== undefined) rows.push(["contrast", data.contrast]);

      document.getElementById("details").innerHTML = `
        <dl>
          ${rows.map(([key, value]) => `<dt>${escapeHtml(key)}</dt><dd>${escapeHtml(value)}</dd>`).join("")}
        </dl>
      `;
    }

    function renderGraph() {
      const view = currentView();
      document.getElementById("min-n").max = String(Math.ceil(view.ranges.sourceNMax));
      document.getElementById("min-n").value = "0";
      document.getElementById("min-n-value").textContent = "0";
      if (cy) cy.destroy();
      cy = cytoscape({
        container: document.getElementById("cy"),
        elements: view.elements,
        style: graphStyle(view),
        layout: { name: "preset", fit: true, padding: 42 },
        wheelSensitivity: 0.25,
        minZoom: 0.18,
        maxZoom: 2.5
      });
      cy.on("select", "node, edge", (event) => renderDetails(event.target));
      cy.on("unselect", "node, edge", () => {
        if (cy.$(":selected").length === 0) {
          document.getElementById("details").textContent = "Select a node or edge to inspect counts and labels.";
        }
      });
      applyFilters();
      cy.fit(undefined, 42);
    }

    function downloadPng() {
      if (!cy) return;
      const view = currentView();
      const link = document.createElement("a");
      link.download = `${view.key}_cytoscape_graph.png`;
      link.href = cy.png({ full: true, scale: 2, bg: "#ffffff" });
      link.click();
    }

    function init() {
      if (typeof cytoscape === "undefined") {
        document.getElementById("load-error").style.display = "block";
        return;
      }
      populateViewSelect();
      populateTraitFilters();
      document.getElementById("view-select").addEventListener("change", renderGraph);
      document.getElementById("search").addEventListener("input", applySearch);
      document.getElementById("min-n").addEventListener("input", applyFilters);
      document.getElementById("layout-layered").addEventListener("click", resetLayeredLayout);
      document.getElementById("layout-force").addEventListener("click", runForceLayout);
      document.getElementById("fit").addEventListener("click", () => cy && cy.fit(undefined, 42));
      document.getElementById("download").addEventListener("click", downloadPng);
      renderGraph();
    }

    window.addEventListener("DOMContentLoaded", init);
  </script>
</body>
</html>
"""


def append_readme_section(output_dir: Path, manifest: pd.DataFrame) -> None:
    """Append or replace the graph explorer section in the generated README."""
    readme_path = output_dir / "README.md"
    if not readme_path.exists():
        return
    content = readme_path.read_text(encoding="utf-8")
    marker = "## Cytoscape Graph Explorer"
    content = content.split(marker, maxsplit=1)[0].rstrip()
    lines = [
        "",
        "",
        marker,
        "",
        f"- Interactive graph explorer: `{HTML_FILENAME}`",
        f"- Inline graph data audit copy: `{DATA_FILENAME}`",
        f"- Graph manifest: `{MANIFEST_FILENAME}`",
        "",
        "The graph explorer uses Cytoscape.js from jsDelivr and includes four "
        "views: increased/decreased grouped conditions and increased/decreased "
        "individual conditions. Each source condition link connects to a "
        "trait/body-site bundle node, then each bundle connects to the body-site "
        "node as one summed trait-colored trunk.",
        "",
        "| view | condition nodes | body-site nodes | trait/body-site bundles | source edges | trunks |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in manifest.itertuples(index=False):
        lines.append(
            f"| {row.view_label} | {row.condition_node_count} | "
            f"{row.body_site_node_count} | {row.trait_body_bundle_count} | "
            f"{row.source_edge_count} | {row.trunk_edge_count} |"
        )
    lines.append("")
    readme_path.write_text(content + "\n".join(lines), encoding="utf-8")


def write_six_shared_trait_graphs(
    output_dir: Path = OUTPUT_DIR,
    repo_root: Path = REPO_ROOT,
) -> pd.DataFrame:
    """Write graph explorer HTML, inline data JSON, and graph manifest TSV."""
    del repo_root
    output_dir.mkdir(parents=True, exist_ok=True)
    payload, manifest = build_graph_payload(output_dir)
    json_text = json.dumps(payload, indent=2)
    safe_inline_json = json_text.replace("</", "<\\/")

    (output_dir / DATA_FILENAME).write_text(json_text + "\n", encoding="utf-8")
    manifest.to_csv(output_dir / MANIFEST_FILENAME, sep="\t", index=False)
    html = (
        html_template()
        .replace("__CYTOSCAPE_CDN__", CYTOSCAPE_CDN)
        .replace("__GRAPH_DATA__", safe_inline_json)
    )
    (output_dir / HTML_FILENAME).write_text(html, encoding="utf-8")
    append_readme_section(output_dir, manifest)
    return manifest


def main() -> None:
    """Write Cytoscape graph outputs from existing six-shared-traits TSVs."""
    manifest = write_six_shared_trait_graphs()
    html_path = OUTPUT_DIR / HTML_FILENAME
    print(f"Wrote Cytoscape graph explorer: {html_path}")
    print(f"Wrote graph views: {len(manifest)}")


if __name__ == "__main__":
    main()
