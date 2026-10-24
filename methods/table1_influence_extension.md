# Additional lineage exclusions and study omissions

Issue [#771](https://github.com/CultureBotAI/microbial-trait-profiling/issues/771).
The executable rule is
[`table1_influence_extension_protocol.json`](table1_influence_extension_protocol.json).
This rule was written after inspecting historical results and the omissions of
the original panel. It defines an exploratory sensitivity analysis. It is not a
preregistered discovery test or independent replication.
The reviewed rule is pinned by its SHA-256 in the implementation before this
sidecar's production execution. Changing the rule while retaining its version
is rejected; any future change needs a reviewed declaration and new pin.

Run on a verified core after the lead finishes the rerun:

```bash
just table1-influence-extension --core-run RUN/core/run.json \
  --noel-run RUN/arms/noel/arm_run.json --output-dir RUN/sidecars/influence-extension
```

The command requires a fresh output directory and uses the core's actual strain
rollup scope, population and direction basis through the public
`run_strain_rollup` and `load_graph` interfaces. It reconciles
every retained profile and matrix count against the source-taxonomy graph before
calculating an effect. It records core consumption, input/producer/output hashes,
the protocol and the complete retained-feature denominator. The lead performs
production execution; implementation tests use synthetic data only.

For each of the six focal traits, each increased/decreased direction, each
family/genus rank and each declared scope, the extension selects the top three
resolved lineages with at least three supporting studies. Ranking uses their
propagated annotation-count contribution, with taxon ID breaking ties. The
ranking population matches the declared `contributing_universe`, before
dropping pairs with a zero annotation total. The
scopes are pooled, exact Feces, exact Saliva and exact Parkinson's disease/Feces.
The candidate table retains unresolved and below-threshold support. Each scope
uses the union of its selections; no p-value or effect enters the ranking.
Every scope also includes a fixed joint exclusion of Lachnospiraceae
(`NCBITaxon:186803`) and Oscillospiraceae (`NCBITaxon:216572`), informed by the
historical cellulose observation. The joint scenario remains visible even
when a scope has no support. Its six-trait effects use the union of source
taxa, so a source with multiple ancestral paths is removed once.

Each exclusion subtracts that lineage's source contributions from every
reported ancestor, across the full denominator. Direct higher-rank assertions
whose membership in the excluded lineage is unresolved remain and are flagged.
Pairs with a zero retained-annotation total on either side are excluded as whole
pairs. Eligibility is recomputed after exclusion. The output separates changes
on common experiment support from changes in which pairs remain eligible.

`lineage_removed_annotation_mass.tsv` reports each scenario, direction and retained
trait on the eligible pre-exclusion profile pairs used by the effect calculation.
It records total retained annotation mass before and after removal and the removed
mass for each retained trait and full denominator. These are annotation counts summed
over reported profile occurrences, not distinct taxa, cell counts or abundance.
Denominator totals repeat for each trait row; within a scenario and direction,
sum the per-trait removed counts to reconcile that total, not the repeated totals.
The effects table repeats the removed mass summed across directions and explicitly
labels `no_removed_annotation_mass`, including a declared lineage absent from the
source graph. A supported removal with unchanged fraction effects can therefore
be distinguished from an exclusion that removed no annotations.

The effect reduction is `1 - effect_after / effect_before`, exported as both
`effect_reduction_fraction` and `effect_reduction_percent` (100 times that ratio).
The effects are the full-data study-mean increased-minus-decreased shared fractions
in the same scope, with eligibility recomputed after removal. Separate `matched_`
fields use only experiments eligible both before and after exclusion; these isolate
the denominator/removal change from eligibility changes. A zero or missing
before-effect, or missing after-effect, has an explicit status and undefined ratio.
The ratio is not clipped: negative values indicate a larger shift in the original
direction, values between zero and one indicate attenuation, and values above one
indicate reversal. It is a descriptive change, not a confidence level or a causal
fraction attributable to a lineage. Historical percentages must be reproduced on
the same scope, cohort and estimator before citation; pooled raw count pairs alone
do not determine this equal-study-weight effect ratio.

The study-omission table removes all experiments from one eligible study at a
time. Each trait's contrast is the increased-minus-decreased annotation fraction,
averaged equally over experiments within study. Every omission recomputes all
six exact two-sided study sign tests, followed by BH-six and BY-six over the
same complete six-member family. Untestable members occupy their family slot
internally at p=1 but retain undefined exported p/q values. Other traits' original
p-values are never held fixed. Tied study means have absolute value at most
1e-12 and do not enter the sign-test denominator.

The sidecar's `influence_extension_summary.md` includes the verified same-core Noel report
when `--noel-run` is supplied, followed by the full additional lineage list and
the six-fraction LOSO ranges. The source report is not edited. Omitting the Noel
manifest leaves an explicit missing-report notice in that combined report.
The declared
`190-native-v2` panel, its native outputs and all three arm plans remain unchanged.
`lineage_fraction_effects.tsv` compares the additional exclusions against the
baseline on total and matched support, and its effects are also included in
the report. Each LOSO row includes the full-data shared-fraction baseline,
BH-six and BY-six. The full-data shared-fraction baseline and BH-six use the
`190-native-v2` panel's `shared_fraction` family in `effects.tsv.gz`
(`pooled_fraction_effects.tsv` for pooled rows),
recomputed on this same core and scope; all-trait native
composition q-values are a different family. `fraction_loso_summary.tsv` gives
descriptive ranges across the dependent study omissions. These conditional
six-test adjusted values do not control selection of the six traits or searching
across lineages, scopes and omissions. They supplement the full-family
multiple-testing sidecar; they do not replace it.

Effect-reduction ratios are undefined when the absolute before-effect is at most
the shared `1e-12` tie tolerance, including floating-point cancellation residue.
The effect itself is retained; this rule prevents an effectively zero denominator
from becoming a huge percentage. Values outside that tolerance remain unclipped.

Each fraction-effect row records JSON `lineage_accounting`: the graph presence,
removed annotation mass and removed focal-trait mass for every declared lineage.
These independent exclusions use the same eligible pre-exclusion profiles as the
joint mass audit. Individual masses can overlap through shared descendants and
must not be summed to recover the joint union, which remains deduplicated.
`absent_from_graph` is distinct from `present_in_graph` with zero removed mass
in a scope. Presence uses recorded nodes, annotation subjects and taxonomy edges;
ancestor traversal's inclusion of the requested ID alone is not evidence of presence.
The combined report displays the individual accounting alongside the union totals.

`just probe-influence-extension-guards baseline OUTPUT` checks the relevant synthetic
controls. Repeat with `loso_bh_by`, `total_matched_percent`, `fraction_as_percent`,
`annotation_focal_mass`, `exact_zero_only`, `post_protocol_check`,
`assumed_graph_presence` and `summed_union_mass`, each with a fresh output directory.
The retained helper loads one isolated weakened module per process, requires the
intended test assertion rather than a setup error, and records canonical source
hashes before and after. These probes verify guards, not biological calibration.

`just test-influence-extension` covers both unit and executed-core integration
tests. Core consumption explicitly binds the reported-taxon table and both KG
inputs in addition to profiles, matrix, metadata and source identity.

All effects concern alignment of reported taxon shifts with organism-level
trait annotations. They do not measure in vivo physiology, metabolite flux,
patient-level abundance or a causal disease mechanism.
