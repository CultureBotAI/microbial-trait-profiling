#!/usr/bin/env python3
"""Generate data/derived/study_trait_profiles_wide.tsv."""

from __future__ import annotations

import sys

from mtp.study_trait_profiles import main

if __name__ == "__main__":
    sys.exit(main())
