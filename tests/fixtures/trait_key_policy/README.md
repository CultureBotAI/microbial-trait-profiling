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
`WorkflowRun` manifest (`core/run.json`) and its provenance files
(`core/steps/identify_sources/provenance/source_identity.json`, `generation.json`)
necessarily contain SHA-256 inventory digests that cannot be rewritten without
invalidating `verify_core_run`. Packaging the normalized directory bundle as a
deterministic `.tar.gz` preserves the exact SHA-256-verified manifest and scientific
artifacts without modifying `.secrets.baseline`.

### Regeneration

From the repository with a clean checkout of the audited `main` producer commit
(`44897c6`) and a fresh durable output directory `$ORIGINAL_OUTPUT`:

```sh
PYTHONDONTWRITEBYTECODE=1 uv run python \
  scripts/adhoc/issue1231_generate_historical_union_core_fixture.py generate \
  --code "$MAIN_44897C6_CODE" \
  --output "$ORIGINAL_OUTPUT" \
  --fixture tests/fixtures/trait_key_policy/historical_union_core.tar.gz
```

`generation.json` inside the archive records `producer_commit`, `producer_matrix_inventory`,
`original_manifest_sha256`, `normalized_manifest_sha256`, and `helper_sha256` (a
non-enforced provenance digest of
`scripts/adhoc/issue1231_generate_historical_union_core_fixture.py` at generation time;
`tests/test_historical_trait_policy_core.py` checks its type and 64-character hex length
per #1146).

### Text, spelling, secret, and YAML checks on archive contents

Because pre-commit skips binary `.tar.gz` files, extract the archive to a fresh scratch
directory `$EXTRACTED_FIXTURE` and run the text, secret, spelling, and schema checks
reproducibly:

```sh
mkdir -p "$EXTRACTED_FIXTURE"
tar -xzf tests/fixtures/trait_key_policy/historical_union_core.tar.gz -C "$EXTRACTED_FIXTURE"
```

1. **No local paths or user identifiers**: verify that all paths in `core/run.json`,
   `core/steps/identify_sources/provenance/source_identity.json`, and `generation.json`
   use the virtual `/fixture/trait-policy-44897c6` namespace and `location_transformations`
   records only `field` and `replacement`:
   ```sh
   ! rg '/Users/|/home/|/private/|/var/folders/' "$EXTRACTED_FIXTURE"
   ```
2. **Secret scan**: run `detect-secrets-hook` across all extracted non-hash text files
   (the only high-entropy strings in the bundle are the SHA-256 digests in
   `core/run.json`, `core/steps/identify_sources/provenance/source_identity.json`, and
   `generation.json`):
   ```sh
   uv run detect-secrets-hook --baseline .secrets.baseline \
     "$EXTRACTED_FIXTURE"/inputs/* \
     "$EXTRACTED_FIXTURE"/core/steps/ingest_bugsigdb/studies/*.yaml \
     "$EXTRACTED_FIXTURE"/core/steps/taxa_matrix/derived/*.tsv \
     "$EXTRACTED_FIXTURE"/core/steps/count_profiles/derived/*.tsv \
     "$EXTRACTED_FIXTURE"/core/steps/ora_profiles/derived/*.tsv
   ```
3. **Spelling and core verification**: run `codespell` across `$EXTRACTED_FIXTURE` and
   verify all 27 core files via `verify_core_run`:
   ```sh
   uv run codespell "$EXTRACTED_FIXTURE"
   uv run pytest -q tests/test_historical_trait_policy_core.py
   ```
