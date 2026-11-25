# Native lineage influence completion

> [!IMPORTANT]
> **Matched-Corpus Table 1-Trait Analog (Issue #210)**
> The native composition and ORA reruns / taxon contribution audits reported here are a matched-corpus Table 1-trait analog sensitivity analysis, not influence on the published manuscript Table 1. They use the pinned human paired KG-Microbe profile corpus and are not a reproduction of the published Table 1 baseline (they do not reproduce published Table 1 baseline values such as facultatively aerobic $q = 2.8 \times 10^{-10}$).

These reruns extend the reconciled 190-v2 audit using protocol 190-native-v2. The 25 rule-declared exclusion scenarios, four scopes, raw trait keys and human paired corpus are unchanged. This is matched-corpus sensitivity analysis, not a reproduction of published results with other filters. See run_manifest.json for hashes, the implementation revision, model policies and result files.

## Composition

Native effects are experiment-weighted increased-minus-decreased z-score contrasts. Native p-values test the mean; native BH families contain all supported variable traits, not just Table 1. Zero-total profiles remain as zero fractions. Study-weighted comparisons, sign tests and conditional intervals are separately labeled in native_composition_comparison.tsv.gz. Refitted scales and fixed baseline SD comparisons must not be conflated.

| Trait | Pooled native effect | Native p | Native q | Largest matched change, fixed baseline SD | Exclusion |
| --- | --- | --- | --- | --- | --- |
| facultatively_aerobic | 0.201152 | 5.37744e-13 | 1.81757e-10 | -0.0343541 | exclude_family_NCBITaxon:543 |
| cellulose_catabolic_process | -0.163088 | 6.959e-09 | 1.17607e-06 | -0.0115341 | exclude_family_NCBITaxon:31979 |
| coccus_shaped | 0.160342 | 2.68484e-08 | 3.02492e-06 | -0.0540449 | exclude_family_NCBITaxon:1300 |
| sphere_shaped | 0.146733 | 2.80024e-07 | 2.36621e-05 | -0.0236589 | exclude_family_NCBITaxon:1300 |
| arsenate_reduction | 0.149326 | 7.40527e-07 | 4.51358e-05 | -0.174593 | exclude_genus_NCBITaxon:561 |
| moderately_halophilic | 0.140286 | 8.01228e-07 | 4.51358e-05 | -0.0481342 | exclude_family_NCBITaxon:81852 |

The displayed maximum is descriptive selection within the rule-declared panel, not an explanation or causal attribution. The TSV retains all opposite, null and untestable results, eligibility components, intervals and family membership.

## Baseline-versus-exclusion native q-values

This table directly checks whether any predeclared genus/family exclusion causes a Table 1-named trait that is significant at pooled baseline to lose native q<=0.05 membership. These baseline q-values are calculated on the matched-corpus Table 1-trait analog, not a reproduction of the published Table 1 baseline. Q-values come from the complete supported-trait family, not a six-trait recalculation. `NA` means the exclusion result was not tested or not selected.

| Trait | Baseline q | Highest exclusion q | Highest-q exclusion | Exclusions losing q<=0.05 | Losing exclusions |
| --- | ---: | ---: | --- | ---: | --- |
| facultatively_aerobic | 1.81757e-10 | 4.16632e-06 | exclude_family_NCBITaxon:543 | 0 | none |
| cellulose_catabolic_process | 1.17607e-06 | 5.76818e-06 | exclude_family_NCBITaxon:1300 | 0 | none |
| coccus_shaped | 3.02492e-06 | 0.0250569 | exclude_family_NCBITaxon:1300 | 0 | none |
| sphere_shaped | 2.36621e-05 | 0.0154112 | exclude_family_NCBITaxon:1300 | 0 | none |
| arsenate_reduction | 4.51358e-05 | 0.806779 | exclude_genus_NCBITaxon:561 | 2 | exclude_family_NCBITaxon:543<br>exclude_genus_NCBITaxon:561 |
| moderately_halophilic | 4.51358e-05 | 0.0171181 | exclude_family_NCBITaxon:81852 | 0 | none |

## ORA

Native enrichment and depletion were rerun against fixed-baseline and recomputed scope backgrounds separately. The global BH families use every eligible retained trait/profile cell. Enrichment tests positive cells; depletion includes zeros in nonzero-total profiles. All six traits remain in the export, not only significant hits. Fold limits are >=1.5 for enrichment and <=0.67 for depletion.

The profile table includes exact 95% central hypergeometric null count intervals: these are NOT confidence intervals for fold effects. The fold-comparison TSV instead gives study-weighted paired fold differences, conditional mean intervals, directional support and matched versus eligibility changes. It also reports native hit counts at q<=0.001 and q<=0.05 separately by curated direction.

| Trait | Background | Largest absolute matched fold change | Exclusion | Lost pairs |
| --- | --- | --- | --- | --- |
| facultatively_aerobic | fixed_baseline | -0.0665669 | exclude_family_NCBITaxon:1300 | 24 |
| facultatively_aerobic | recomputed_remaining | -0.0332522 | exclude_family_NCBITaxon:1300 | 24 |
| cellulose_catabolic_process | fixed_baseline | -0.103828 | exclude_family_NCBITaxon:31979 | 20 |
| cellulose_catabolic_process | recomputed_remaining | -0.146493 | exclude_family_NCBITaxon:31979 | 20 |
| coccus_shaped | fixed_baseline | -0.0977837 | exclude_family_NCBITaxon:1300 | 24 |
| coccus_shaped | recomputed_remaining | 0.0787379 | exclude_family_NCBITaxon:90964 | 18 |
| sphere_shaped | fixed_baseline | 0.0823701 | exclude_family_NCBITaxon:33958 | 47 |
| sphere_shaped | recomputed_remaining | 0.0803242 | exclude_family_NCBITaxon:33958 | 47 |
| arsenate_reduction | fixed_baseline | -0.279762 | exclude_genus_NCBITaxon:561 | 6 |
| arsenate_reduction | recomputed_remaining | -0.303311 | exclude_genus_NCBITaxon:561 | 6 |
| moderately_halophilic | fixed_baseline | -0.0554214 | exclude_family_NCBITaxon:81852 | 7 |
| moderately_halophilic | recomputed_remaining | -0.0513827 | exclude_family_NCBITaxon:81852 | 7 |

## Integration Limits

For #178, use qualified trait keys, exact scope, scenario, background/reference, effect units and correction-family fields. Keep native and shared statistics separate; do not combine ORA backgrounds in a common ranking. ORA profile-count tests are not independent corroboration of study-level effects. The native composition t intervals assume independent experiments; study-weighted summaries are supplied separately, without silently changing the native model (#133).

For #179, these are annotation influence diagnostics. Higher-rank direct assertions that cannot be allocated remain retained and flagged in the source ledger. Missing source identities/ranks are not repaired by statistical reruns. Disease/control orientation (#188), feature semantics (#189), source-count validity (#184) and full-pipeline calibration (#182) remain independent review requirements. No causal, oxygen-exposure, SCFA or taxonomy-independent mechanism follows from these results. Native kernel/ORA owner review and integration sign-off are still required before manuscript adoption.

## Cross-scale interpretation flags (#213)

The study-weighted raw-fraction mean and the native transformed analysis are different estimands. The raw-fraction directional q-value is a third statistic: it tests the share of non-tied study directions, not whether the mean-effect interval excludes zero. The machine-readable `interpretation_handoff.tsv` keeps those fields separate for #178.

| Trait | Study-weighted raw-fraction finding | Native transformed finding | Interpretation |
| --- | --- | --- | --- |
| Arsenate reduction | Excluding Escherichia changes the matched mean by -0.0061 pp [-0.0088, -0.0034]; the remaining mean is -0.0005 pp [-0.0016, 0.0006]. Escherichia supplies 58.7% of increased-arm counts. | q changes from 4.51358e-05 to 0.806779 after the genus exclusion. | Lineage-sensitive on this corpus; do not describe a lineage-independent mechanism or treat removal as causal attribution. |
| Sphere shaped | Baseline mean 0.0191 pp [-0.0137, 0.0519] crosses zero; the separate directional sign-test q is 0.00101414. | Baseline native q is 2.36621e-05. | Cross-scale non-reproduction of the mean interval, not a refutation or independent replication of the native result. |

## Review Bundle

The `.tsv.gz` comparison files contain all six traits, scenarios and scopes, including null and untestable results. `native_table1_stats.tsv.gz` is the six-trait subset of the native composition table; its q-values were computed over the FULL native families, whose sizes remain attached. The separate ORA family table includes both backgrounds and both tail-family sizes. `interpretation_handoff.tsv` carries the two measured #213 flags without collapsing their scales or tests. `scenarios.tsv` supplies lineage names, support and selection status for the full genus/family inventory, including lineages outside the rule-declared panel.

The full 317-trait composition table, profile-level ORA results and source ledger are generated under `data/derived/table1_taxon_influence/`, following the repository's untracked generated-data convention. Regenerate the default input chain with `just table1-taxon-influence`, or use `just table1-taxon-influence-from <pinned-KG-directory>` for already materialized inputs. The recipe writes a candidate review bundle under the untracked run directory; this committed directory is an explicitly reviewed snapshot. This compact review bundle does not discard any six-trait influence comparison. `run_manifest.json` records hashes for the full data; `publication_manifest.json` identifies these review copies. All full-run file hashes were checked before publication.

Publication also re-sums the complete source-contribution ledger through the profile/reported-taxon links and requires exact equality with every exported profile/trait count. The publication manifest records those file hashes and row counts.
