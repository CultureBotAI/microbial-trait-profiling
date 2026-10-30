# Pinned-parent study Wilcoxon sensitivity audit

Issue #1073 changes the secondary study-level signed-rank test from the historical
`study-mean-delta-z-wilcoxon-v1` contract to the v2 contract described in
[three_test_table.md](three_test_table.md). This audit quantifies that change from
completed parent outputs. It reads the persisted study means and does not
recompute Δz, run an owner arm, or rewrite the parent bundle.

## Invocation and inputs

Run the actual recipe from the intended producer checkout, supplying the
independently obtained SHA-256 of `three_test_run.json`:

```sh
just three-test-wilcoxon-audit \
  --parent /durable/parent/three_test_run.json \
  --parent-sha256 FULL_PARENT_SHA256 \
  --scope phenotype \
  --output /durable/new-wilcoxon-audit \
  --code-commit FULL_PRODUCER_COMMIT
```

`--parent` also accepts a parent directory. Scope is required and must be `off`
or `phenotype`; it must match the pinned manifest. Output must be a fresh,
durable directory disjoint from every parent file alias and its target. No
scratch directory is needed. A supplied producer commit must match the executing
checkout. The record also reports whether that checkout is dirty; qualifying
scientific runs use a clean committed producer.

The four inputs are `three_test_run.json`, `per_trait_evidence.tsv`,
`composition_pair_delta_z.tsv`, and `composition_study_delta_z.tsv`. The caller
pins the manifest; its output inventory binds the tables. Validation checks the
declared family and scope, pair/study identity sets, and the old
Wilcoxon statistic, p-values, BH/BY values and support counts before computing
the new secondary test. All retained traits remain in the correction family;
untestable traits contribute internal p=1 and have blank reported inference.

## What is verified

The verification scope is
`pinned_parent_outputs_and_legacy_reproduction_only`. Parent run IDs, historical
producer commits, dirty status, population and direction are **pinned
declarations**, not a fresh verification of their historical execution. The
audit does not reopen the upstream comparison, core, arm, KG, or BugSigDB
sources, inspect the historical Git tree, or claim to reproduce every parent
analysis. Reproduction of the old Wilcoxon columns establishes the numerical
baseline for this bounded sensitivity calculation.

The new audit has its own producer record: executing commit, dirty status,
Python/numerical-library versions, lockfile and implementation inventory. That
inventory includes the statistical methods document and this audit protocol.
Input paths are serialized relative to the pinned parent root as filenames,
SHA-256 digests and sizes; absolute machine-local paths are not public inputs.
Input byte hashes, resolved targets and producer identity are checked again
before exclusive directory publication. A validation failure creates no output
ancestors; publication never replaces an existing destination.

This audit requires nonempty parent pair and study ledgers. Empty-population
audits are outside its input scope; this restriction makes no claim about the
upstream pipeline's handling of empty vectors or populations.

## Comparing old and new inference

The TSV retains archived, reproduced-v1 and new-v2 statistics separately. Old
columns must reproduce within relative tolerance `1e-12` and absolute tolerance
`1e-15`. Change counts compare v2 to the **reproduced v1** values, so accepted
last-bit differences in archived serialization are not attributed to the new
method.

For each of p, BH q and BY q, the summary reports exact floating-point inequality
and material inequality beyond the same declared tolerances, among traits
tested under both contracts. It also reports the largest absolute difference
and the comparable trait count. Status transitions are reported separately;
they are not silently counted as numeric differences between finite p-values
and missing values. BH/BY call changes use `q <= 0.05`, counting an untestable
trait as not rejected, and therefore include relevant status transitions.

The new method counts distinguish exact untied, exhaustive tied, asymptotic,
and untestable outcomes according to the v2 contract. These results remain
secondary inference on alignment between reported taxon shifts and trait
annotations; they are not measurements of in vivo physiology.

The fresh bundle contains only `wilcoxon_audit.tsv` and
`wilcoxon_audit_run.json`. The JSON contains the parent pin, relative input
inventory, numerical audit summary, own producer identity and the TSV hash.
