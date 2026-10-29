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

`study_t_primary_correction` records BY. For every trait, record the number and
fraction of paired input profiles with a positive annotation count. The native
selected family excludes traits present in fewer than six profiles or more than
95% of profiles, among other eligibility rules. Full-family results retain those
traits. Because fractions sum to one, even a constant positive annotation count
can produce a contrast when the retained totals differ; this closure effect is
especially relevant to near-ubiquitous traits.

For overlapping native-selected traits, check profile support, experiment/study
counts, both means and study SD against the full-family calculation. Where both
study t-tests are estimable, also require agreement of raw t, p and the study
95% interval; the differently sized families' q values are not compared.
Use the native whitespace-stripped human-host rule and require canonical `inc`
and `dec` directions, rejecting aliases that the native population would omit.
Keep the declared full-family SD threshold of 1e-12. A native study test excluded
by its broader `np.isclose(SD, 0)` rule is explicitly identified in
`native_selected_overlap_status`; any unexplained population or estimate mismatch
fails before publication. The native profile standardization also treats exact
zero SD differently from the adapter's 1e-12 threshold, so resulting mismatches
are rejected rather than silently reconciled.

The secondary Wilcoxon test was declared at the study level before computation:
each study contributes one mean delta-z, regardless of its number of experiments.
The current `study-mean-delta-z-wilcoxon-v2` contract uses two-sided signed rank
with `zero_method="wilcox"` and no continuity correction. In a private test vector,
exclude study means with absolute value <=1e-12, using the same zero tolerance as
Noel's study votes. Do not round remaining values or modify stored experiment or
study contrasts, descriptive means, ranks, t-tests or confidence intervals.

Method selection uses the number of remaining studies and exact equality of
absolute magnitudes. With no tied magnitudes, use SciPy's exact signed-rank
distribution through 50 active studies. With ties and at most 13 active studies,
use `PermutationMethod(n_resamples=np.inf)` to enumerate every sign assignment;
the untied signed-rank distribution is not exact for ties. Otherwise use the
asymptotic normal approximation with tie correction and no continuity correction.
These explicit thresholds make the policy independent of SciPy's automatic
method selection. Exact sign enumeration is deterministic and uses no random
sample. These p-values retain the signed-rank test's symmetry and independence
assumptions; exact computation alone does not establish biological calibration.
The 50/13 cutoffs follow the documented SciPy 1.17 policy, applied explicitly
after tolerance filtering: [SciPy 1.17.0 Wilcoxon notes](https://docs.scipy.org/doc/scipy-1.17.0/reference/generated/scipy.stats.wilcoxon.html#notes).
The implementation and tests use locked SciPy 1.17.1; its version is also recorded
in each run's numerical environment.

Record the selected method, zero tolerance, both cutoffs and active/zero study
counts in each trait row. Nonempty vectors with no active studies have reason
`all_tied`; absent vectors have reason `insufficient_support`; both have method
`not_tested`. Leave their reported statistic/p/q blank and retain internal p=1.
Correct across the full retained catalog with BY and companion BH. The test is
secondary and never a headline; it does not reuse the six-trait raw-count
density-adjustment Wilcoxon results. Study means give one observation per study,
so additional experiments do not add signed-rank observations.
`study_wilcoxon_statistic` is the two-sided minimum of the
positive- and negative-rank sums, min(W+, W−); it carries no sign.

The legacy `study-mean-delta-z-wilcoxon-v1` contract used exact-zero removal and a
normal approximation for every nonempty, nonzero sample. Historical v1 evidence
must retain that identity. Compare v1 and v2 on the same stored study vectors and
full retained families; do not rewrite historical p/q fields or infer independent
replication from the two calculations.

## ORA declaration: ora-study-enrichment-direction-v1

Declared before computation; accepted through the coordination review on 2026-10-25.
Use only the complete native upper-tail export's `is_significant` flag: global
separate-tail BH q <=0.001 and fold >=1.5. Count increased and decreased enriched
lists separately, retaining explicit zeros after family completeness validation.
For K=0, enrichment and depletion counts are blank because no native test exists;
these are not observed zero-hit counts.
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
non-hits. The preference label is the raw, unadjusted majority of non-tied study
votes: positive, negative, mixed for equal opposing votes, or none when every
eligible study ties. With no eligible studies it is unavailable. There is
no depletion-direction test and this column never decides a headline call.

Every row records:
`uncalibrated (methods/calibration.md excludes native ORA)`.

The test runs on thresholded calls whose global BH threshold depends on all traits
and lists, so hits are not independent across traits; it describes enrichment
preference only. BY does not establish the validity of the underlying sign-test
or taxon-sampling model. Repeated experiments are reduced within study; shared
cohorts and cross-study dependence remain possible.

Thresholded calls also depend on list size N: the same trait-bearing share can
produce a hit in a larger list and no hit in a smaller list. The directional
summary does not adjust for this difference in detection power. Columns
`ora_direction_pairs_n_inc_gt_dec`, `ora_direction_pairs_n_inc_lt_dec` and
`ora_direction_pairs_n_inc_eq_dec` count the eligible pairs with larger increased,
larger decreased and equal-sized lists, respectively. These scope-level counts
are repeated for each K>0 trait; their sum is its eligible-pair count. For K=0
all three are zero because the trait has no eligible directional pairs.

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

Label-support imports retain the existing audit's 14-trait `AUDIT_TRAITS`
coverage (`focal_and_selected_oxygen_traits_only`) by default, and extend to
any declared `--trait-list` recorded in the audit's `provenance.json`
(`trait_list`, #1132) with `carrier_catalog_coverage` naming
`focal_and_selected_oxygen_traits_plus_<name>`. A carrier
is not a source record. Report direct-at-reported-taxon and inherited support
only for supported, reconciled audit traits; mark unsupported traits explicitly.
Direct and inherited reported-taxon counts can overlap. Unique carrier IDs are
deduplicated across reported ancestors and descendants. The importer checks the
audit's bound inputs, complete focal matrix cells, signature occurrences,
carrier-count coherence, and per-cell reconciliation between
`unresolved_madin_carriers` in `taxon_trait_label_support.tsv` and
`madin_resolution` in `carrier_source_evidence.tsv`; it does not independently
reconstruct the KG annotations. Do not equate unique source names or carrier
taxa with original record counts. Madin source resolution requires its
separately pinned table, and `carrier_madin_status` is reported per trait
(#1151): `not_supplied` when no label-support sidecar is supplied,
`not_in_focal_audit` for traits outside the audited catalog, `not_applicable`
for audited traits with no bound `infores:madin_etal` carrier evidence rows,
`not_supplied_source_support_unknown` when bound `infores:madin_etal` carriers
exist without `--madin-csv`, `unresolved` when `--madin-csv` was verified but
at least one bound `infores:madin_etal` carrier row remained
`madin_resolution=unknown`, and `pinned_verified` when every bound
`infores:madin_etal` carrier row resolved (`madin_resolution=resolved`; this
per-trait `pinned_verified` reuses the audit-wide `provenance.json`
`madin_status` name with a narrower meaning requiring all bound
`infores:madin_etal` carriers of that trait to resolve).
Missing original record identifiers remain unresolved.

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
Its role is verified co-presence: its inventory and replay are checked, but its
adjusted values are not imported into, or cross-checked against, this table.
`native_selected_study_q_by` is reconstructed from the native selected-family
raw p values. The manifest records `multiplicity_binding` as
`verified_co_presence_only` or `not_supplied`. Optional-sidecar replay retains
the multiplicity stage's existing producer-code and policy identity requirements.
Record the new producer Git commit, dirty state, implementation hash, dependency
lock hash, input hashes/run IDs and declared methods. Native values remain
separately named; each new test has its own family. Never compare p-values across
tests or treat agreement on shared data as independent replication.

Each generation or complete table replay verifies the core and each arm once.
The comparison and optional multiplicity replay reuse those verified records;
all existing execution, binding, policy, inventory and output checks still apply.
The input inventory explicitly includes declared file core artifacts consumed
during verification, including external KG nodes/edges and reported taxa, along
with the native and optional evidence files. File hashes are rechecked before
publication to detect changes during adaptation. A declared file remains a
required input if it disappears or changes filesystem type; it is never silently
removed from the inventory. The same declared storage kind controls source-root
protection before verification.

All generated outputs and temporary scratch must use durable directories outside
system temporary roots and remain disjoint from source bundles. Historical
cellulose t=-5.24 and full-family BY q=0.00047 are reproduction targets for the
off Tier1+2 unit, not values to tune the implementation toward. Report a mismatch
with its source and cause before producing downstream claims.

The temporary-root guard calls `mtp.workflow.output_roots.temporary_root`: it checks fixed
macOS temporary roots, the process temporary directory and `TMPDIR`, both as
written and after resolving symlinks. Unlisted roots can escape detection (#1051),
including `/tmp` on Linux when the process temporary directory is redirected
elsewhere, and `/var/tmp` when neither dynamic root identifies it. This limitation does
not authorize temporary storage for these runs. Choose the explicit durable
scratch directory outside the current process temporary directory; the sidecar
then routes internal temporary files into that validated scratch directory.
An explicit `--scratch-dir` underneath inherited `TMPDIR` is refused before any
verification or write, even when that temporary directory has a durable-looking name.

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

Publication stages all tables and the completed manifest in a fresh durable
sibling directory, rechecks input and producer identities, then publishes the
directory by an exclusive atomic rename. A failed late check leaves no final
output directory; a competing output destination is never overwritten. The
dirty-state flag covers tracked changes, while source-file hashes continue to
bind the implementation. Untracked logs or scratch files do not change that flag.
Exclusive publication uses Darwin `renamex_np` or Linux `renameat2`. An
unsupported platform, C library or filesystem operation fails closed; it never
falls back to a rename that could replace an existing destination.

Before expensive source verification, a private probe in the nearest existing
output ancestor checks that exclusive rename both rejects an existing destination and
can publish to an absent one. The probe follows checks that keep output and scratch
separate from declared source directories, external core files and every symlink
hop. These declaration-only checks also precede the first replay verification.
The probe never creates missing output parents or targets the requested output.
Output parents are created only when staging the verified result for publication. Final
publication still uses exclusive rename. Staging is created with an exclusive
plain `mkdir`, so the filesystem applies `0777 & ~umask` without a worker or a
query that changes the process's mask. Its hidden directory can be readable under
that mask, but only the final destination names a published bundle. A failure
before publication removes this run's staging; a creation collision or a path
created by another process after publication is preserved.
Optional count columns use nullable
integers: populated counts serialize as integers, and unavailable counts stay blank.

An uncatchable termination such as SIGKILL can leave hidden
`.<name>-staging-*` directories in the output parent or private `.<name>-probe-*`
directories in its nearest existing ancestor. They are
unpublished intermediates, even if they contain a manifest. After confirming that
their producer has stopped, they may be removed; do not remove an active staging
directory. Only the final requested output directory is a published bundle.

Verification checks the exact output inventory and regenerates results from the
verified source bundles. It compares the complete record, including numerical
environment and implementation identity. Editing both a table and its local
output checksum therefore does not make an altered result pass replay.

`just test-three-test-table` covers synthetic contracts. The separate
`test-three-test-off-reproduction` recipe requires a completed, verified historical
off sidecar and its original comparison; it checks the rounded cellulose target
and the 317-trait/548-study/935-pair scope. This data-dependent check does not run
any owner arm or reconstruct missing source data.

`just probe-three-test-guards baseline --output-dir OUTPUT --scratch-dir SCRATCH`
runs the carrier/ORA positive control. The same recipe accepts named guard
mutations listed by `--help`, loading isolated source copies and requiring the
intended regression failures. Each invocation needs fresh durable directories;
its logs and canonical before/after hashes show that source files were unchanged.

```bash
just test-three-test-off-reproduction "$COMPARISON" "$OUTPUT" "$SCRATCH" \
  511c891 "$MULTIPLICITY" "$LABEL_SUPPORT"
```

The final two arguments are optional; use empty strings when those auxiliary
inputs were not supplied to the original sidecar. Direct pytest invocation must
also supply a nonempty hexadecimal `MTP_THREE_TEST_BASIS` before replay.
