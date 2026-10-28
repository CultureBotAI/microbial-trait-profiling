# Portable historical core fixture (#1134)

This is a **location-normalized test copy**, not a byte-identical execution
record. The synthetic core was generated and fully verified using the clean
`511c891` producer whose matrix implementation is independently pinned in
`../historical_implementations.json`. It records `strain_rollup=phenotype`, no
trait-key policy, and the historical seven-column feature metadata.

The KG has three distinct contributors: a direct `METPO:1000608` annotation on
taxon 101, a direct BactoTraits `METPO:1000612` annotation on taxon 102, and a
BacDive strain `METPO:1000605` annotation linked to taxon 103. Genus 100 therefore
has one contributor in each of the three separate raw columns. These are small
synthetic inputs, not research results.

The generation helper copies the run and inputs, then replaces only provenance
locations with the portable `/fixture/trait-policy-511c891` namespace. It
updates the source-identity artifact's hash, size and inventory after changing
its consumed-file paths. It also normalizes the recorded interpreter path,
which the workflow fingerprint explicitly excludes. No scientific TSV/YAML or
step publication is changed. `generation.json` records every location
transformation, both source-identity hashes, the original and copied manifest
hashes, the exact helper hash, and the verified producer identity. The original
generation remains unchanged in durable development storage.

Regenerate from the repository with a clean historical checkout and two fresh
durable output directories:

```sh
PYTHONDONTWRITEBYTECODE=1 uv run python \
  scripts/adhoc/issue1134_generate_historical_core_fixture.py generate \
  --code "$HISTORICAL_CODE" --output "$ORIGINAL_OUTPUT" --fixture "$FIXTURE_OUTPUT"
```

The helper imports the selected historical code with bytecode writing disabled,
checks its Git identity and audited matrix hash, generates the original, runs
its real verifier, and fully verifies the normalized copy with the same code.
Copy the resulting `core/`, `inputs/` and `generation.json` here after verification.
The CI tests additionally relocate this committed copy and call the current
`verify_core_run` with a strict test-only location resolver. Nothing is mocked;
no original-path fallback is permitted. Inventory corruption and a matrix count
changed with a matching replacement inventory must both fail verification.

For the negative control, export main commit `409a9c6` into a separate source
directory and run:

```sh
PYTHONDONTWRITEBYTECODE=1 uv run python \
  scripts/adhoc/issue1134_generate_historical_core_fixture.py verify \
  --code "$CONTROL_CODE" --fixture "$FIXTURE_OUTPUT" --output "$CONTROL_REPORT" \
  --expected-error 'Matrix feature vocabulary/order differs from consumed KG evidence'
```

That verifier rejects the historical raw matrix because it applies the implicit
union rule. The current verifier accepts the same fixture using the recorded
historical producer identity. All helper and test commands used for this change
were run through the coordination heavy-run wrapper with a durable scratch root.

The historical YAML serializer emits unindented block lists. A dedicated
pre-commit YAML hook checks these files with that indentation rule and all other
repository YAML rules unchanged. Reformatting would invalidate their recorded
artifact hashes; the fixture bytes are preserved instead.
