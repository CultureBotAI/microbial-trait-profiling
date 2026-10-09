#!/usr/bin/env python3
"""Forwarder to trait_enrichment.py configured for depletion (--tail lower)."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import trait_enrichment

hypergeom_cdf = trait_enrichment.hypergeom_cdf
hypergeom_cdf_many = trait_enrichment.hypergeom_cdf_many
hypergeom_sf = trait_enrichment.hypergeom_sf
hypergeom_sf_many = trait_enrichment.hypergeom_sf_many


def main(argv: list[str] | None = None) -> None:
    raw_args = sys.argv[1:] if argv is None else argv
    forwarded = ["--tail", "lower", *raw_args] if "--tail" not in raw_args else list(raw_args)
    trait_enrichment.main(forwarded)


if __name__ == "__main__":
    main()
