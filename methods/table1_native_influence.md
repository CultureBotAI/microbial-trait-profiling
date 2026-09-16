# Native Table 1 lineage influence

This completes the statistical rerun extension to [#190's contribution audit](table1_taxon_influence.md).
The source ledger is unchanged. The
[native protocol](table1_native_influence_protocol.json) applies the same rule-declared
exclusions as the shared-fraction rerun, with four exact scopes, 317 retained columns
and 4,670 human paired profiles from 959 studies. Nothing is selected by the new effects,
directions or p-values. See the generated
[native summary and review tables](table1-native-results/README.md) for results.

## Reproduction

Run the complete pipeline with `just table1-taxon-influence <pinned-KG-directory>`.
It now performs the source reconciliation, shared effects and native extension.
For an already completed, unchanged contribution audit, use
`just table1-native-influence <derived-input-directory> <audit-output-directory>`.
Both paths use the same native implementation. The second recipe checks the
four input SHA-256 hashes against the completed source audit, reconstructs every
remaining denominator column from the sparse removals, and verifies the original
and remaining six-trait counts and full totals against the scenario ledger.
Changed inputs, duplicate removal entries, negative counts, missing profile
records or inconsistent totals stop the run. No replacement KG is downloaded.

`native_manifest.json` records the source manifest hash, protocol, runtime versions,
code commit and source-file hashes, removal-ledger input hashes and result hashes.
It starts with `completed: false`; a failed/interrupted run must not be integrated.
The historical contribution manifest's `not_rerun` statements describe that earlier
stage, not the native extension. Only a completed native manifest establishes that
the native files exist and match the recorded inputs.

`just table1-native-summary <audit-output-directory>` verifies all native input and
output artifact hashes and publishes the compact review bundle in
`methods/table1-native-results/`. This is also the final step of the native recipe.
The bundle contains every six-trait comparison, native Table 1 statistics, full
native family sizes and run manifests. Large full-profile data remain in the
generated-data directory under the repository's existing convention.
The bundle also includes the complete `scenarios.tsv` inventory, with lineage
names and selection/support status for candidates outside the declared panel.

## Composition Contract

The numerical kernels follow the harmonized composition implementation at
`b931b3a` (`create_trait_profile_infographic.py`): full non-MediaDive row totals,
zero-total fractions equal to zero, profile support at least six and at most
95% of profiles, positive sample variance, and no top-N trait limit. All columns
are exported, including those not selected. This uses the exact #190 scopes,
not the differently grouped historical manuscript strata.

Within each scope and scenario, transform with `log1p(1000 * fraction)` and
refit the population mean/SD. Centering cancels in paired differences. The native
two-sided one-sample t test weights experiment contrasts equally. Fewer than
three pairs or contrast SD within `1e-8` of zero is untestable. BH uses every
finite native p-value among selected traits in that scope/scenario. Untestable
selected traits retain the native q=1 convention, missing p and explicit status;
unselected traits have no native p/q. q=1 does not make an untestable result null
evidence. Tests do not use the six-trait-only family.

Native mean intervals target experiment-weighted effects; exact binomial intervals
describe the positive share of non-tied experiments. Experiments from one study
are not necessarily independent, so native uncertainty has the existing #133
limitation. Separate study-weighted comparisons average experiments within study
and then studies equally, using the original audit's conditional mean and
directional-share intervals. Their sign-test p-values are **not native t-test
p-values**; their BH6 families are separate for each scope/scenario/reference and
before/after/matched/change population.

Two reference policies are exported: refitted z-scores and contrasts divided by
the fixed baseline SD. The latter isolates change on a common scale; a change in
refitted SD is not itself removed biological signal. If the baseline SD is zero,
the fixed-scale comparison is unavailable. Refitted comparisons retain the native
zero-column-SD convention (zero contrasts) symmetrically before and after, so a
constant baseline does not create artificial pair gains. Native test columns always refer to
the native refitted model, even on the fixed-reference comparison row.

## ORA Contract

The numerical kernels follow Ilyas's enrichment/depletion implementation at
`3ae7ac8`, before the separate study-dependence modeling changes. For each declared
scope, the population is the sum of counts over **all** retained columns and both
directions of its human paired profiles. Each nonzero-total profile is a draw.
The fixed-baseline population and the recomputed remaining population are separate
analyses, not candidates for a shared ranking. This is a matched-corpus native
model comparison, not reproduction of outputs built with other host/site filters.

Enrichment uses the inclusive upper hypergeometric tail on positive cells only.
Depletion uses the inclusive lower tail on every retained trait, including zero
cells, in nonzero-total profiles. Global BH is computed before extracting the six
traits, separately by scope/scenario/background/tail. Reported membership uses
q<=0.001 or q<=0.05 plus the existing fold limits (>=1.5 or <=0.67). Native ORA
does not give equal weights to studies; it tests individual count profiles.

Zero-total profiles remain visible with no tests. A zero-background trait remains
in the depletion family with p=q=1, but has undefined fold, explicit status and no
significance membership. No pseudocount or artificial infinite fold is introduced.
Each profile also has an exact central 95% **null count interval**, not a confidence
interval for its fold. Dependent propagated count units mean this is a conditional
count-model diagnostic, not validated false-discovery control (#182/#184).

Study-weighted paired fold differences are supplied as an additional descriptive
representation with conditional mean/directional-share intervals and separate BH6
sign-test families. Both profile totals and the expected trait count must be
positive for those pairs. Before/after own populations, common eligible pairs,
lost/gained pairs and eligibility components remain separate. Loss of significant
ORA cells alone is not an explanation of the signal.

## Files And Handoff

| File | Contents |
| --- | --- |
| `native_composition.tsv.gz` | All 317 traits, native selection, effects, t-test p/q, intervals, directional support and native family sizes for every scope/scenario. |
| `native_composition_comparison.tsv` | All six traits; native before/after fields prefixed `native_`; separate study-weighted comparisons and two reference policies. |
| `<scenario>/native_ora_profiles.tsv.gz` | All six traits and both directions; observed/expected counts, folds, both tails' p/q/status/membership and null count intervals, separately by scope/background. Join baseline and scenario rows by study/experiment/direction/trait/scope/background. |
| `native_ora_families.tsv` | Complete native family sizes, active profile counts and population totals. |
| `native_ora_fold_comparison.tsv` | Six-trait study-weighted fold comparisons and native significant-cell counts by direction at both thresholds. |
| `native_manifest.json`, `native_summary.md` | Reproducibility record and measured interpretation. |

Supply these alongside the original ledger/shared effects to #178. Join on the
full qualified trait key, exact scope, scenario and reference/background; keep
units, weights, missingness and correction families. Opposite and untestable
results belong in the integration, not just the five originally upward traits.
For #179, retain annotation-level interpretation, unresolved higher-rank source
allocation, taxonomy coverage and corpus restrictions. Native model owner review,
#188 orientation and #189 semantics review remain required for manuscript adoption;
this work does not resolve those issues or #182/#184. Taxonomic removal is not a
causal intervention or proof of taxonomy-independent mechanisms.

## Public Checksum Review

The secret scanner's reviewed baseline contains exact public checksum values, not
credentials: the BugSigDB export revision and dump checksum in `download.py`, the
two pinned KG file checksums, and the two native source-code revisions. Windows
and POSIX filename spellings are recorded because detect-secrets 1.5.0 uses native
path separators. Changed values and other lines remain scanned. This replaces the
earlier whole-line checksum exemptions and addresses the relevant #180 check gate.
