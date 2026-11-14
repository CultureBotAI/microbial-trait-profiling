# Study-split assessment of exploratory condition triples (#854)

Protocol `study-split-selection-v2` in `study_split_selection_protocol.json`
implements the owner's phenotype-primary decision of 2026-10-25 (#1076). The
original declaration is preserved byte-for-byte as
`study_split_selection_protocol_v1.json`; v2 names its path, version and SHA-256,
and the reader verifies those bytes before accepting the amendment. All within-scope
statistical rules remain unchanged. New families and run manifests carry v2;
historical v1 outputs require their recorded implementation for replay.

This is a post-results amendment. The historical corpus had already been inspected,
and both off and phenotype selection-scan results were known. The exact disclosure
is `amendments[0].known_when_amended` in
`cellulose_decr_selection_scan_protocol_v2.json`: it distinguishes the recorded
off evidence from the phenotype result used outside the repository's recorded
evidence. The amendment does not assert that study-split results were unseen.
This remains an internal exploratory assessment, not external replication or
retrospective preregistration.

`phenotype` is primary (`primary_phenotype`); `off` is the declared sensitivity
(`declared_sensitivity`). Both must be reported side by side without selecting the
more favorable result. The core's verified record supplies the execution scope.
The explicit `absent_means: off` declaration applies only to historical core
records without a strain-rollup parameter. Calls to the selection and held-out
helpers that omit their scope use the current primary, `phenotype`; explicit
`off` remains supported. Unknown scopes are refused even for empty selections.

For each seeded split, all experiments and signatures of a study stay together.
The training half supplies the Noel directional trait screen over all body sites
and the exact-Feces condition inventory. Among its ten highest-ranked decreasing
traits, the scan API chooses
one study-disjoint condition triple per trait by smallest observed mean pairwise
Jaccard distance. Exact ties use the sorted condition labels. This is a declared
selection rule for the split analysis, not a reconstruction of the undocumented
historical selection. Cellulose and PD/HIV/colorectal cancer are never forced in.

Only the held-out half supplies the observed overlap and U6 whole-study-signature
null pool. Both overlap statistics enter one declared family per split, including
selected hypotheses that cannot be tested. BY is primary and BH is a comparison.
The retained core vocabulary defines the candidate screen size; the selected
family size is the number of selected triples times two. These are different
quantities and both are recorded. Changing the strain scope changes the core
vocabulary, rather than a hard-coded family size.
This is `mt-supp-study-split-u6`: the two overlap statistics for every frozen
selected triple within one split, including untestable members. Read its actual
membership and size from that split's correction record in `split_results.json`.
The BY-primary and companion BH rejection flags use q <=0.05 and q <=0.001;
no family is formed by pooling repeated splits.

The shared v1.2 family record copies the verified strain scope, ordered feature-filter
identity and source profile representation. Each member names its exact three
condition labels, Feces stratum and split ID as explicit scientific coordinates.
Study-split raw hypotheses retain the foundation's compatible
`legacy_free_text_v1` reason vocabulary; the selected statistical rules and
untestable reasons are unchanged. Corrections record the v1.2 engine and include
the reason-vocabulary field in their serialized identity. Replay historical
outputs with their recorded implementation and locked environment; unchanged
numerical rules do not imply byte-identical records across contract versions.

A held-out triple with shared studies is untestable because U6 allocates distinct
studies to its condition blocks. If any null draw has undefined mean distance,
the distance test is untestable: undefined draws must not silently count as
non-extreme. The core-size statistic can remain testable. No replacement is
chosen using held-out information. The U6 model still assumes exchangeability of
study signatures and does not establish causal, physiological or clinical effects.

Repeated splits overlap and their p-values are dependent. Every split is retained,
including empty selections. Each split records whether candidates were selected,
no trait passed the training rule, fewer than three conditions were eligible, or
no testable study-disjoint triple remained. Stability reports contain descriptive selection,
testability and rejection counts only. There is no averaged p-value, minimum
p-value, best-split inference, across-split FDR guarantee or independent
replication claim. Trait-bearing taxa refer to database annotations of reported
taxa; overlap does not measure pathway activity in the patient.

Use `just selection-study-splits <comparison-run.json> <new-output-directory>` only
after the lead has verified the Tier 1+2 core, arms and comparison. Optional
split/draw/seed overrides are recorded as off-protocol, for fixtures or explicitly
declared sensitivity runs. Add `--verify` with the original parameter overrides
to replay the full analysis from verified sources and compare every output.
The native arm calls and outputs remain unchanged.

After producing and verifying one v2 run for each scope, build the joint report
with the module CLI (the existing Justfile recipe continues to run one scope):

```sh
uv run python -m mtp.cellulose_decr.study_split_run \
  --phenotype-run evidence/phenotype/study-splits \
  --off-run evidence/off/study-splits \
  --output-dir evidence/study-split-scopes
```

The reader checks the exact saved output inventory and hashes the same buffers
that it parses. It requires the current v2 declaration, equal split parameters,
implementation and environment, the same population and direction, trait-key
policy, raw corpus and KG hashes, core producer identity and core parameters
other than strain scope. Each scope must name a distinct core run, and the study
partitions must match. The selected traits and triples may differ by scope.
Per-split family scope and direction must agree with their run. The saved
stability TSV must equal the descriptive counts reconstructed from the split
results. The reader performs no selection or null simulation.

`study_split_scopes.tsv` contains both saved stability tables with scope and role
columns, phenotype first. `study_split_scopes.json` records both run-manifest
hashes, their input and output hashes, common comparison contract and split counts;
empty-selection runs remain represented even when their stability table has no
rows. `README.md` states the roles and inferential limits. Output directories must
be new or empty and must not overlap either source run. The joint `--verify` mode
reconstructs and compares this report without writing. It differs from a single
run's `--verify`, which replays the null calculations against verified upstream
sources.

The joint report checks recorded evidence, not live upstream bundles. Its hashes
detect changed bytes relative to a saved manifest; its cross-checks detect
inconsistent recorded scope and stability. They do not authenticate a coordinated
rewrite of manifests and their dependent outputs. Use full per-run replay for
live source-to-result verification, and preserve externally recorded manifest
hashes for an independent historical identity.

The sidecar binds every Git-tracked canonical `src/mtp` Python file through the
shared multiplicity code inventory, the consumed protocols and methods, and the
dependency lock. Untracked scratch files and editor locks do not change that
inventory. Unresolved tracked Python conflicts, a missing tracked module, or an
inventory that excludes the shared `mtp.multiplicity.stage` module are rejected.
The runner also explicitly hashes both study-split modules and the other paths
in its `IMPLEMENTATIONS` inventory. It records the repository commit and dirty
state as historical execution context, preserved during replay. Output publication
or a commit changing only unbound documentation may change that context without
changing the scientific identity checked by replay. This methods file is explicitly
bound, so edits to it do change the recorded identity.
This includes transitive helpers for pairing and sign tallies, even when a helper
change happens to leave a particular fixture's numbers unchanged. Publication
writes the exact UTF-8 bytes that are hashed, independently of host newline rules.

Input hashes come from the verified upstream inventories, including the comparison
manifest digest captured before verification. They are checked before reading the
analysis inputs, when constructing their source records, and before publication.
The runner requires a recorded digest for every consumed source; a later disk
snapshot cannot replace the upstream identity. Code and protocol identity are
captured before upstream verification and checked again before analysis and
publication. The comparison-bound invocation does not consume a multiplicity
sidecar manifest; it shares the verified-source and correction implementation.

These entry and final checks do not attest Python bytecode imported before run
entry, so use a fresh CLI process. They also do not bind every non-streamed read
to its consumed bytes against transient edit-and-restore races; that shared
limitation is tracked in #983. The recorded results remain an internal exploratory
assessment under the statistical limitations above.

The #883 implementation is merged into main. Update this author branch by
merging `origin/main`, never by rebasing. The #974/#977/#979 fixes and tests include
merge-carried changes
that a rebase can silently drop. Preserve their ancestry and explicit regression
files, then repeat the required checks after integration.
