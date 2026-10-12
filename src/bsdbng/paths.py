"""Project path helpers."""

from __future__ import annotations

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
RAW_DATA_DIR = REPO_ROOT / "data" / "raw"
STUDY_DATA_DIR = REPO_ROOT / "data" / "studies"
DERIVED_DATA_DIR = REPO_ROOT / "data" / "derived"
KG_MICROBE_DIR = REPO_ROOT / "data" / "kg-microbe"
SCHEMA_PATH = REPO_ROOT / "schema" / "bsdbng.yaml"
