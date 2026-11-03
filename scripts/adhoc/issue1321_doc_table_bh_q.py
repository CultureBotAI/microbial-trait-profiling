"""Relabel the documented six-trait density table with BH q-values (#1321); 2026-11-02.

Purpose: the section 5 table of docs/density_adjustment_audit.md records the raw and
adjusted Wilcoxon p-values and statuses of an earlier six-trait run, labelled from
p <= 0.05. This recomputes each row's BH q-values over the six documented traits with
mtp.density_adjustment.family_bh_q_values and its status with adjustment_status, to
show whether any documented label changes under the q-based rule.

Inputs: docs/density_adjustment_audit.md (the committed document; no data files).
Outputs: a TSV on stdout, or at --output, with one row per documented trait: trait,
raw and adjusted p, raw and adjusted BH q, the status read from p and from q.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

from mtp.density_adjustment import adjustment_status, family_bh_q_values

DOC = Path(__file__).resolve().parents[2] / "docs" / "density_adjustment_audit.md"
_SCIENTIFIC = re.compile(r"^([+-]?\d+(?:\.\d+)?)\s*\\times\s*10\^\{(-?\d+)\}$")
_ESCAPED_PIPE = "\x00"


def _number(cell: str) -> float:
    """Parse a documented value such as **$2.97 \\times 10^{-12}$** or $+356.65$."""
    text = cell.replace("*", "").replace("$", "").strip()
    match = _SCIENTIFIC.match(text)
    if match:
        return float(match.group(1)) * 10.0 ** int(match.group(2))
    return float(text)


def documented_rows(doc: Path) -> list[dict[str, object]]:
    """The six-trait rows of the section 5 table, in document order."""
    rows: list[dict[str, object]] = []
    for line in doc.read_text(encoding="utf-8").splitlines():
        if not line.startswith("| `METPO:"):
            continue
        cells = [
            cell.strip().replace(_ESCAPED_PIPE, "|")
            for cell in line.replace("\\|", _ESCAPED_PIPE).strip("|").split("|")
        ]
        rows.append(
            {
                "trait_id": cells[0].strip("`"),
                "trait_label": cells[1],
                "raw_mean_diff": _number(cells[2]),
                "raw_wilcoxon_p": _number(cells[3]),
                "adj_mean_diff": _number(cells[4]),
                "adj_wilcoxon_p": _number(cells[5]),
                "documented_status": cells[8].replace("*", "").strip(),
            }
        )
    if len(rows) != 6:
        raise ValueError(f"expected the six documented Table 1 rows, found {len(rows)}")
    return rows


def relabel(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    """Add BH q-values over the documented family and both statuses to each row."""
    raw_p = [float(str(row["raw_wilcoxon_p"])) for row in rows]
    adj_p = [float(str(row["adj_wilcoxon_p"])) for row in rows]
    raw_q = family_bh_q_values(raw_p, family_size=len(rows))
    adj_q = family_bh_q_values(adj_p, family_size=len(rows))
    out = []
    for row, rp, ap, rq, aq in zip(rows, raw_p, adj_p, raw_q, adj_q, strict=True):
        raw_mean = float(str(row["raw_mean_diff"]))
        adj_mean = float(str(row["adj_mean_diff"]))
        out.append(
            {
                "trait_id": row["trait_id"],
                "trait_label": row["trait_label"],
                "raw_wilcoxon_p": f"{rp:.3g}",
                "adj_wilcoxon_p": f"{ap:.3g}",
                "raw_wilcoxon_bh_q": f"{rq:.3g}",
                "adj_wilcoxon_bh_q": f"{aq:.3g}",
                "documented_status": row["documented_status"],
                "status_from_p": adjustment_status(raw_mean, rp, adj_mean, ap),
                "status_from_q": adjustment_status(raw_mean, rq, adj_mean, aq),
            }
        )
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--doc", type=Path, default=DOC)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    table = relabel(documented_rows(args.doc))
    handle = sys.stdout if args.output is None else args.output.open("w", encoding="utf-8")
    try:
        writer = csv.DictWriter(handle, fieldnames=list(table[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(table)
    finally:
        if args.output is not None:
            handle.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
