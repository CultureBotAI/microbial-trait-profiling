# Run-bound per-trait evidence (#1044)

This additive table consumes one verified core, its three completed arms and their
comparison. It preserves distinct estimands, populations and correction families.
`phenotype` is primary; `off` is sensitivity. The historical F inputs are development
and reconciliation evidence. Tracked result tables wait for the selected new main
basis P. No generated table changes any native arm call.

## Declared composition quantities

The denominator is every retained non-MediaDive annotation column. With annotation
fraction f, define L = log1p(1000 f), then standardize each trait over the paired
input profiles using population SD (ddof=0). The experiment contrast is increased
minus decreased z. The study contrast is the mean of its experiment contrasts.
These are annotation shares and standardized annotation contrasts, not measured
physiology or patient-level abundance.

Persist both experiment-weighted and study-weighted means. Rank each signed mean
in ascending and descending order and its magnitude in descending order; exact
rank ties share the minimum rank. No primary ranking has been selected.

Reuse `JedCompositionAdapter.run_inference` for the full retained-family,
two-sided study-mean t-test and its Student-t 95% interval. Report BY and BH over
that whole family. The adapter's historical zero-total convention is retained:
zero retained totals give zero fractions. Per-trait profile SD <=1e-12 gives z=0.
At least two study observations and study sample SD >1e-12 are required for an
estimable t-test. Otherwise retain descriptive effects when defined, leave reported
t/p/q/interval blank and record the reason; retain p=1 inside the correction.
Do not replace the selected native pooled-composition BH/BY family with the full
family. Keep its original results in separately named columns.

The secondary Wilcoxon test uses the study-level mean delta-z: two-sided signed
rank, `zero_method="wilcox"`, `method="approx"`, continuity correction disabled.
Zero differences are removed from ranking. All-zero or absent contrasts are
untestable. Correct across the full retained catalog with BY and companion BH,
holding untestable traits at internal p=1. It is secondary and never a headline.
It does not reuse the six-trait raw-count density-adjustment Wilcoxon results.

## ORA declaration: ora-study-enrichment-direction-v1

Declared before computation; accepted through the coordination review on 2026-10-25.
Use only the complete native upper-tail export's `is_significant` flag: global
separate-tail BH q <=0.001 and fold >=1.5. Count increased and decreased enriched
lists separately, retaining explicit zeros after family completeness validation.
Lower-tail depletion counts (native q <=0.001 and fold <=0.5) are separate columns;
enrichment among decreased lists is not depletion.

For each trait, use experiment pairs for which both directions have N>0 and the
trait has background K>0. Let H(e,d,t) be the binary native enrichment flag. For
study s, D(s,t) = sum_e[H(e,inc,t)-H(e,dec,t)] / eligible_pair_count(s,t).
Give each study one equally weighted vote. Determine its sign from the integer
numerator; exact zero is a tie, with no floating-point tie tolerance. Report
positive, negative and tied studies, eligible studies/pairs, excluded pairs,
mean study D, positive-share and its exact Clopper-Pearson 95% interval.

Apply a two-sided exact binomial sign test to non-tied studies under probability
0.5. BY is primary and BH is a companion, over every trait in the bound retained
catalog separately per rollup/population. All-tied, no-eligible and no-background
traits stay in the family at internal p=1, with blank reported p/q and explicit
status. Missing rows must fail completeness checks; they are never imputed as
non-hits. Positive/negative/mixed describes enrichment preference only. There is
no depletion-direction test and this column never decides a headline call.

Every row records:
`uncalibrated (methods/calibration.md excludes native ORA)`.

The test runs on thresholded calls whose global BH threshold depends on all traits
and lists, so hits are not independent across traits; it describes enrichment
preference only. BY does not establish the validity of the underlying sign-test
or taxon-sampling model. Repeated experiments are reduced within study; shared
cohorts and cross-study dependence remain possible.

The native urn counts distinct matched taxa against all rows of the taxon-by-trait
matrix, including unannotated rows. It is not host-filtered or subset-filtered.
Mixed ranks and nested taxa can inherit the same labels and are not independent
biological draws. Record M, K and native eligible-list support rather than hard-code
numbers. The ORA harness's complete-null flag is not biological calibration;
closing #551 aligned integration and calibration with production, without proving
mixed-signal FDR control or pooled trait significance.

## Native votes and label support

Reuse Noel's raw-fraction study votes, exact binomial p and nominal full-family BY.
Ties use its native tolerance 1e-12. Noel excludes a complete experiment pair if
one retained annotation total is zero; preserve this eligibility difference from
historical transformed composition. Effective-family BY is not the nominal result.

Label-support imports retain the existing audit's focal/oxygen coverage. A carrier
is not a source record. Report direct-at-reported-taxon and inherited support
only for supported, reconciled audit traits; mark unsupported traits explicitly.
Direct and inherited reported-taxon counts can overlap. Unique carrier IDs are
deduplicated across reported ancestors and descendants. The importer checks the
audit's bound inputs, complete focal matrix cells, signature occurrences and
carrier-count coherence; it does not independently reconstruct the KG annotations.
Do not equate unique source names or carrier taxa with original record counts.
Madin source resolution requires its separately pinned table. Missing original
record identifiers remain unresolved.

`METPO:1000603` means anaerobic; `METPO:1000607` means obligately anaerobic.
The anaerobe guild combines both. Preserve qualified trait IDs; never substitute
one term for the other based on a display label. The historical six-trait delta-z
script uses neither oxygen term; the retained expanded summary uses 1000603,
whereas the obligately-anaerobic plate uses 1000607.

## Execution and provenance

Use materialized completed `off` or `phenotype` run bundles with their recorded
paths available. The exploratory `full` rollup is not yet supported by the shared
source-verification contract.
Verification retains upstream execution identities and hashes, population,
direction basis, catalog and rollup scope. A newer sidecar may consume supported
historical runs; it does not rewrite their producer identity. A supplied
multiplicity sidecar must pass its own verification and bind the same comparison.
Record the new producer Git commit, dirty state, implementation hash, dependency
lock hash, input hashes/run IDs and declared methods. Native values remain
separately named; each new test has its own family. Never compare p-values across
tests or treat agreement on shared data as independent replication.

All generated outputs and temporary scratch must use durable directories outside
system temporary roots and remain disjoint from source bundles. Historical
cellulose t=-5.24 and full-family BY q=0.00047 are reproduction targets for the
off Tier1+2 unit, not values to tune the implementation toward. Report a mismatch
with its source and cause before producing downstream claims.

The temporary-root guard follows the shared workflow policy: it checks fixed
macOS temporary roots, the process temporary directory and `TMPDIR`, both as
written and after resolving symlinks. Unlisted roots can escape detection (#1051),
including `/tmp` on Linux when the process temporary directory is redirected
elsewhere, and `/var/tmp` when neither dynamic root identifies it. This limitation does
not authorize temporary storage for these runs. Choose the explicit durable
scratch directory outside the current process temporary directory; the sidecar
then routes internal temporary files into that validated scratch directory.

## Commands and output contract

```bash
just three-test-table --comparison-run "$COMPARISON" --output "$OUTPUT" \
  --scratch-dir "$SCRATCH" --code-commit "$(git rev-parse HEAD)"
just three-test-verify --comparison-run "$COMPARISON" --output "$OUTPUT" \
  --scratch-dir "$SCRATCH"
```

Optional `--multiplicity-run` and `--label-support` arguments bind matching,
verified auxiliary evidence. Supply the same auxiliary inputs during replay.
An absent carrier audit is explicitly unavailable, not a zero carrier count.
The carrier population is reported taxa observed in the bound count-profile
population; it must not be silently described as the eligible human study-test
population. Carrier and source-record counts are distinct.

| Output | Contract |
| --- | --- |
| `per_trait_evidence.tsv` | One row per qualified retained trait, in catalog order; explicit family, support, status and scope fields. `study_t_*` and `study_wilcoxon_*` are full-family tests; `noel_*` are unchanged native votes; `native_selected_*` retain selected-family composition; `ora_*` separate native calls and study preference; `carrier_*` declare bounded support. |
| `composition_pair_delta_z.tsv` | All retained traits by `(study_id, experiment_id)`; increased minus decreased, using the same transformed profiles as the full-family observed t-test. |
| `composition_study_delta_z.tsv` | All retained traits by `study_id`; equal mean of that study's eligible experiment contrasts. |
| `three_test_run.json` | Completed provenance: producer commit/dirty state/code inventory, dependency lock and numerical library versions, input hashes, original run IDs/commits, analysis bases, ordered trait family and output hashes. Contains no machine-local input/output paths. |

Verification checks the exact output inventory and regenerates results from the
verified source bundles. It compares the complete record, including numerical
environment and implementation identity. Editing both a table and its local
output checksum therefore does not make an altered result pass replay.

`just test-three-test-table` covers synthetic contracts. The separate
`test-three-test-off-reproduction` recipe requires a completed, verified historical
off sidecar and its original comparison; it checks the rounded cellulose target
and the 317-trait/548-study/935-pair scope. This data-dependent check does not run
any owner arm or reconstruct missing source data.

```bash
just test-three-test-off-reproduction "$COMPARISON" "$OUTPUT" "$SCRATCH" \
  511c891 "$MULTIPLICITY" "$LABEL_SUPPORT"
```

The final two arguments are optional; use empty strings when those auxiliary
inputs were not supplied to the original sidecar. Direct pytest invocation must
also supply a nonempty hexadecimal `MTP_THREE_TEST_BASIS` before replay.
