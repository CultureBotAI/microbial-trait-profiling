#!/usr/bin/env python3
"""Generate data/derived/taxa_trait_matrix.tsv."""

from __future__ import annotations

import sys

from mtp.taxa_trait_matrix import main

if __name__ == "__main__":
    sys.exit(main())
