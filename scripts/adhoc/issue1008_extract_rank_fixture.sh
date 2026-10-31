#!/usr/bin/env bash
# Purpose: retain one verbatim pinned KG rank-class record and its native header.
# Input: pinned merged-kg_edges.tsv; output: tests/fixtures/coverage/taxonomy_rank_class.tsv.
# Issue #1008; 2026-10-23. Reads only through the first matching record; no producer run.
set -euo pipefail
repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
source_edges="${1:?Pass the pinned merged-kg_edges.tsv path}"
fixture_path="$repo_root/tests/fixtures/coverage/taxonomy_rank_class.tsv"
temporary_fixture="$fixture_path.partial"
trap 'rm -f "$temporary_fixture"' EXIT
head -n 1 "$source_edges" > "$temporary_fixture"
rg --max-count 1 --fixed-strings 'urn:uuid:3e65a2cc-1904-4020-9ca7-1b766bf3b79e' \
    "$source_edges" >> "$temporary_fixture"
mv "$temporary_fixture" "$fixture_path"
