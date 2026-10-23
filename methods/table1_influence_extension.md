# Additional lineage exclusions and study omissions

Issue [#771](https://github.com/CultureBotAI/microbial-trait-profiling/issues/771).
The executable rule is
[`table1_influence_extension_protocol.json`](table1_influence_extension_protocol.json).
This rule was written after inspecting historical results and the omissions of
the original panel. It defines an exploratory sensitivity analysis. It is not a
preregistered discovery test or independent replication.

Run on a verified core after the lead finishes the rerun:

```bash
just table1-influence-extension --core-run RUN/core/run.json \
  --noel-run RUN/arms/noel/arm_run.json --output-dir RUN/sidecars/influence-extension
```

The command requires a fresh output directory and uses the core's actual strain
rollup scope, population and direction basis. Its scope integration depends on
PR #828's public `run_strain_rollup` and `load_graph` interfaces. It reconciles
every retained profile and matrix count against the source-taxonomy graph before
calculating an effect. It records core consumption, input/producer/output hashes,
the protocol and the complete retained-feature denominator. The lead performs
production execution; implementation tests use synthetic data only.

For each of the six focal traits, each increased/decreased direction, each
family/genus rank and each declared scope, the extension selects the top three
resolved lineages with at least three supporting studies. Ranking uses their
propagated annotation-count contribution, with taxon ID breaking ties. The
scopes are pooled, exact Feces, exact Saliva and exact Parkinson's disease/Feces.
The candidate table retains unresolved and below-threshold support. Each scope
uses the union of its selections; no p-value or effect enters the ranking.

Each exclusion subtracts that lineage's source contributions from every
reported ancestor, across the full denominator. Direct higher-rank assertions
whose membership in the excluded lineage is unresolved remain and are flagged.
Pairs with a zero retained-annotation total on either side are excluded as whole
pairs. Eligibility is recomputed after exclusion. The output separates changes
on common experiment support from changes in which pairs remain eligible.

The study-omission table removes all experiments from one eligible study at a
time. Each trait's contrast is the increased-minus-decreased annotation fraction,
averaged equally over experiments within study. Every omission recomputes all
six exact two-sided study sign tests, followed by BH-six and BY-six over the
same complete six-member family. Untestable members occupy their family slot
internally at p=1 but retain undefined exported p/q values. Other traits' original
p-values are never held fixed. Tied study means have absolute value at most
1e-12 and do not enter the sign-test denominator.

The sidecar's `native_summary.md` includes the verified same-core Noel baseline
when `--noel-run` is supplied, followed by the full additional lineage list and
the six-fraction LOSO ranges. The source report is not edited. Omitting the Noel
manifest leaves an explicit missing-baseline notice in that combined report.
The declared
`190-native-v2` panel, its native outputs and all three arm plans remain unchanged.
`lineage_fraction_effects.tsv` compares the additional exclusions against the
baseline on total and matched support; `fraction_loso_summary.tsv` gives
descriptive ranges across the dependent study omissions. These conditional
six-test adjusted values do not control selection of the six traits or searching
across lineages, scopes and omissions. They supplement the full-family
multiple-testing sidecar; they do not replace it.

All effects concern alignment of reported taxon shifts with organism-level
trait annotations. They do not measure in vivo physiology, metabolite flux,
patient-level abundance or a causal disease mechanism.
