# Study-split assessment of exploratory condition triples (#854)

The protocol in `study_split_selection_protocol.json` is declared before this
sidecar is run on real data. The historical corpus has already been inspected.
This is an internal exploratory assessment, not external replication or a claim
that the historical choices were preregistered.

The primary scope is `off`, matching the historical lead selection and the talk's
comparison. `phenotype` is a declared sensitivity. Other scopes are rejected;
both accepted scopes must be reported separately without selecting the more
favorable result. The core's verified record supplies the scope.

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

The v1.1 family record copies the verified strain scope, ordered feature-filter
identity and source profile representation. Each member names its exact three
condition labels, Feces stratum and split ID as explicit scientific coordinates.

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

The sidecar binds every Git-tracked canonical `src/mtp` Python file through the
shared multiplicity code inventory, the consumed protocols and methods, and the
dependency lock. Untracked scratch files and editor locks do not change that
inventory. Unresolved tracked Python conflicts, a missing tracked module, or an
inventory that excludes the executing stage are rejected. It records the repository
commit and dirty state as historical execution context, preserved during replay. Output publication
or a documentation-only commit may change that context without changing the
scientific identity checked by replay.
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

Update this author branch by merging the final #883 head and later `origin/main`,
never by rebasing. The #974/#977/#979 fixes and tests include merge-carried changes
that a rebase can silently drop. Preserve their ancestry and explicit regression
files, then repeat the required checks after integration.
