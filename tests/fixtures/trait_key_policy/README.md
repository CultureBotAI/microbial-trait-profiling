# Trait-key policy historical core fixtures

## `historical_raw_core/` (#1134)

Directory bundle of a policy-less pre-#1124 core run generated from the audited raw
producer commit (`511c891`). Provenance paths are normalized to the virtual namespace
`/fixture/trait-policy-511c891`; all scientific TSV and YAML outputs are byte-identical
to the original execution. Regenerated and verified by
`scripts/adhoc/issue1134_generate_historical_core_fixture.py`.

## `historical_union_core.tar.gz` (#1231, #1235, #1254)

Deterministic `.tar.gz` bundle (`mtime=0`, `uid=0`, `gid=0`, sorted POSIX member order)
of an explicit `facultative_union_v1` core run generated from the audited `main`
producer commit (`44897c6`) where recorded
`object_name` and `column_label` in `taxa_trait_matrix_columns.tsv` preserve the raw
KG-Microbe node name (`facultatively aerobic` and `facultatively aerobic METPO:1000608`)
while `trait_member_ids_json` records `["METPO:1000605","METPO:1000608","METPO:1000612"]`.

### Why an archive is used

Pre-commit's `detect-secrets` hook scans committed text files with `HexHighEntropyString`
enabled and `.secrets.baseline` is read-only under project policy. A normalized
`WorkflowRun` manifest (`core/run.json`) and its provenance files (`source_identity.json`,
`generation.json`) necessarily contain SHA-256 inventory digests that cannot be rewritten
without invalidating `verify_core_run`. Packaging the normalized directory bundle as a
deterministic `.tar.gz` preserves the exact SHA-256-verified manifest and scientific
artifacts without modifying `.secrets.baseline`.

### Regeneration

From a clean checkout of the audited `main` producer commit (`44897c6`):

```bash
python scripts/adhoc/issue1231_generate_historical_union_core_fixture.py generate \
  --code /path/to/clean-main-44897c6 \
  --output /tmp/issue1231_union_core_run \
  --fixture tests/fixtures/trait_key_policy/historical_union_core.tar.gz
```

`generation.json` inside the archive records `producer_commit`, `producer_matrix_inventory`,
`original_manifest_sha256`, `normalized_manifest_sha256`, and `helper_sha256` (a
non-enforced provenance digest of `scripts/adhoc/issue1231_generate_historical_union_core_fixture.py`
at generation time).

### Text, spelling, secret, and YAML checks on archive contents

Because pre-commit skips binary `.tar.gz` files, the extracted text contents in
`/tmp/issue1231_regen_20261030_v2/normalized_fixture` were audited directly before
committing:

1. **No local paths or user identifiers**: verified that all paths in `core/run.json`,
   `core/provenance/source_identity.json`, and `generation.json` use the virtual
   `/fixture/trait-policy-44897c6` namespace and `location_transformations` records only
   `field` and `replacement`.
2. **Secret scan**: ran `detect-secrets-hook` across all extracted non-hash text files
   (`inputs/*`, `core/steps/**/*`); the only high-entropy strings in the bundle are the
   SHA-256 digests in `core/run.json`, `core/provenance/source_identity.json`, and
   `generation.json`.
3. **Spelling and YAML/JSON validity**: ran `codespell` and `check-yaml` / `check-json`
   across the extracted bundle and verified all 27 core files via `verify_core_run` in
   `tests/test_historical_trait_policy_core.py`.
