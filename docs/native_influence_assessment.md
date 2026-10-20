# Native influence standing in integrated evidence

Issue [#671](https://github.com/CultureBotAI/microbial-trait-profiling/issues/671)
separates completion of the native influence diagnostic from a scientific
assessment of any procedure. `completed=true` never means clear, calibrated,
robust, or causally validated.

`just integrated-evidence --influence-status PATH` accepts:

- A `190-native-v2` `native_manifest.json`. If its sibling
  `review/publication_manifest.json` exists, the adapter checks that review
  handoff and its compact tables. Otherwise it checks the native manifest's
  declared input and output artifacts.
- A compact review `publication_manifest.json`, beside its `run_manifest.json`
  and the six published review tables. This supports reviewing the published
  bundle without the large source tables.
- An explicit `taxon-influence-assessment-v1` assessment as described below.

The adapter checks input profile and experiment hashes against the integrated
run. It checks native protocol compatibility, completion field type, source
hash formats, review-to-native identity, review table bytes and row counts,
and the complete scenario/scope/background family inventory. The review's
reported source reconciliation must have zero mismatches. The adapter does not
repeat that producer's reconciliation against unavailable full source tables;
it preserves that report, its source hashes, and the verification scope.

Issue #707 additionally requires the complete native producer artifact inventory,
including each selected scenario's input counts and ORA profile output. Native
composition and comparison table row counts must match the native result counts.
The compact bundle's complete composition and ORA comparison counts must also
match those native counts; its intentionally selected six-trait statistics table
is not required to have the full composition table's row count. A bundle missing
analysis products cannot acquire completed-result standing from a family inventory
alone.

The assessment sidecar keeps execution completion, code revision and dirty
state, protocol and runtime metadata, source hashes,  diagnostic scenario and
scope detail, and any scoped interpretation flags. A lineage-specific flag
stays attached to its feature and scenario; it is not copied into every row's
procedure standing.

| Input state | Standing on a row without an explicit procedure verdict |
| --- | --- |
| No input file | `pending #190: influence evidence not supplied` |
| Valid native manifest with `completed=false` | `native influence execution incomplete #190` |
| Completed and compatible native/review evidence | `native influence results reported; no procedure-level verdict supplied #190` |
| Invalid JSON, incompatible contract, changed artifact, inconsistent inventory or mismatched input corpus | `unreadable #190: ...` |
| Explicit assessment that omits this procedure | `no procedure-level influence verdict supplied #190` |

`pending` means no applicable assessment was supplied or mapped. It does not
establish that the analysis never ran elsewhere. The generated data dictionary
uses this meaning too.

An explicit procedure assessment must identify the input hashes, the verdict's
basis, and that it applies to every integrated row of that procedure. For example,
the structural contract is:

```json
{
  "schema_version": "taxon-influence-assessment-v1",
  "input_sha256": {
    "study_trait_profiles_wide.tsv": "<actual 64-character SHA-256>",
    "experiments.tsv": "<actual 64-character SHA-256>"
  },
  "procedure_verdicts": {
    "fraction_sign_test": {
      "verdict": "<the assessment's actual verdict>",
      "applies_to": "all_rows_of_procedure",
      "basis": "<the assessment's procedure-wide evidential basis>"
    }
  }
}
```

The example is a contract template, not an assessment of these results. A
scenario-specific verdict cannot be promoted through this interface. A global
`taxon_influence_status` string is refused because it has no procedure scope.
An optional `native_manifest_sha256` is retained as a declared evidence pointer;
it is not independently authenticated when the assessment is supplied alone.
The adapter does not invent an assessment from native numerical results.

## Bounded production verification

On 2026-10-19 the integrated table was regenerated using the completed Noel run
from #660, without rerunning native influence. The source native manifest hash
was `e44f93deb12e2a78993d77234b227fb8f648957c99f097fc49d67a502d898f4d`;
the review manifest hash was
`3eff59a0c00b0c92ac76f54d6398976e812f58df479124347a28e0150e1554b7`.
The regenerated sidecar retained 25 scenarios, four scopes, 200 ORA family
inventory rows and both existing interpretation flags. The procedure-verdict
mapping remained empty, while standing correctly reported completed native
results with no procedure-level verdict.

Reproduce the bounded check with the completed run's paths:

```bash
just integrated-evidence \
  --profile-input "$NOEL_RUN/inputs/study_trait_profiles_wide.tsv" \
  --experiments-input "$NOEL_RUN/inputs/experiments.tsv" \
  --full-dump-input "$NOEL_RUN/inputs/full_dump.csv" \
  --influence-status "$NOEL_RUN/table1_taxon_influence/native_manifest.json" \
  --output-dir "$NEW_EVIDENCE_DIR"
```

Regression coverage checks absent, incomplete, completed, explicit scoped and
malformed evidence, cross-corpus mismatches, changed native/review artifacts,
path traversal, scoped handoff export and generated wording.
