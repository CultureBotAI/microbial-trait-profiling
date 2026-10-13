# Native lineage influence completion

These reruns extend the reconciled 190-v2 audit using protocol 190-native-v2. The 24 rule-declared exclusion scenarios, four scopes, raw trait keys and human paired corpus are unchanged. This is matched-corpus sensitivity analysis, not a reproduction of published results with other filters. See run_manifest.json for hashes, the implementation revision, model policies and result files. This historical manifest predates the explicit dirty-tree field: its commit identifies the revision but does not prove the checkout was clean. The current publisher requires a recorded clean state.

## Composition

Native effects are experiment-weighted increased-minus-decreased z-score contrasts. Native p-values test the mean; native BH families contain all supported variable traits, not just Table 1. Zero-total profiles remain as zero fractions. Study-weighted comparisons, sign tests and conditional intervals are separately labeled in native_composition_comparison.tsv.gz. Refitted scales and fixed baseline SD comparisons must not be conflated.

| Trait | Pooled native effect | Native p | Native q | Largest matched change, fixed baseline SD | Exclusion |
| --- | --- | --- | --- | --- | --- | --- |
| facultatively_aerobic | 0.18631 | 5.35241e-12 | 1.68601e-09 | -0.0544385 | exclude_family_NCBITaxon:543 |
| cellulose_catabolic_process | -0.164876 | 2.10484e-09 | 3.31512e-07 | 0.00949518 | exclude_family_NCBITaxon:543 |
| coccus_shaped | 0.156906 | 2.35116e-08 | 2.46872e-06 | -0.0657366 | exclude_family_NCBITaxon:1300 |
| sphere_shaped | 0.139652 | 5.18261e-07 | 3.33761e-05 | -0.0305592 | exclude_family_NCBITaxon:1300 |
| arsenate_reduction | 0.147797 | 5.2978e-07 | 3.33761e-05 | -0.197758 | exclude_family_NCBITaxon:543 |
| moderately_halophilic | 0.131709 | 1.7926e-06 | 9.41115e-05 | -0.0530599 | exclude_family_NCBITaxon:81852 |

The displayed maximum is descriptive selection within the rule-declared panel, not an explanation or causal attribution. The TSV retains all opposite, null and untestable results, eligibility components, intervals and family membership.

## Baseline-versus-exclusion native q-values

This table directly checks whether any predeclared genus/family exclusion causes a Table 1 trait that is significant at pooled baseline to lose native q<=0.05 membership. Q-values come from the complete supported-trait family, not a six-trait recalculation. `NA` means the exclusion result was not tested or not selected.

| Trait | Baseline q | Highest exclusion q | Highest-q exclusion | Exclusions losing q<=0.05 | Losing exclusions |
| --- | ---: | ---: | --- | ---: | --- |
| facultatively_aerobic | 1.68601e-09 | 3.43643e-05 | exclude_family_NCBITaxon:543 | 0 | none |
| cellulose_catabolic_process | 3.31512e-07 | 3.39141e-06 | exclude_family_NCBITaxon:1300 | 0 | none |
| coccus_shaped | 2.46872e-06 | 0.0185388 | exclude_family_NCBITaxon:1300 | 0 | none |
| sphere_shaped | 3.33761e-05 | 0.0164878 | exclude_family_NCBITaxon:1300 | 0 | none |
| arsenate_reduction | 3.33761e-05 | 0.673778 | exclude_genus_NCBITaxon:561 | 2 | exclude_family_NCBITaxon:543<br>exclude_genus_NCBITaxon:561 |
| moderately_halophilic | 9.41115e-05 | 0.0234308 | exclude_family_NCBITaxon:81852 | 0 | none |

## ORA

Native enrichment and depletion were rerun against fixed-baseline and recomputed scope backgrounds separately. The global BH families use every eligible retained trait/profile cell. Enrichment tests positive cells; depletion includes zeros in nonzero-total profiles. All six traits remain in the export, not only significant hits. Fold limits are >=1.5 for enrichment and <=0.67 for depletion.

The profile table includes exact 95% central hypergeometric null count intervals: these are NOT confidence intervals for fold effects. The fold-comparison TSV instead gives study-weighted paired fold differences, conditional mean intervals, directional support and matched versus eligibility changes. It also reports native hit counts at q<=0.001 and q<=0.05 separately by curated direction.

| Trait | Background | Largest absolute matched fold change | Exclusion | Lost pairs |
| --- | --- | --- | --- | --- |
| facultatively_aerobic | fixed_baseline | -0.0693003 | exclude_family_NCBITaxon:543 | 38 |
| facultatively_aerobic | recomputed_remaining | -0.0329404 | exclude_genus_NCBITaxon:561 | 4 |
| cellulose_catabolic_process | fixed_baseline | -0.102376 | exclude_family_NCBITaxon:31979 | 20 |
| cellulose_catabolic_process | recomputed_remaining | -0.147721 | exclude_family_NCBITaxon:31979 | 20 |
| coccus_shaped | fixed_baseline | -0.100947 | exclude_family_NCBITaxon:1300 | 26 |
| coccus_shaped | recomputed_remaining | 0.0930277 | exclude_genus_NCBITaxon:1279 | 22 |
| sphere_shaped | fixed_baseline | 0.0836635 | exclude_family_NCBITaxon:33958 | 48 |
| sphere_shaped | recomputed_remaining | 0.0821474 | exclude_family_NCBITaxon:33958 | 48 |
| arsenate_reduction | fixed_baseline | -0.268529 | exclude_genus_NCBITaxon:561 | 4 |
| arsenate_reduction | recomputed_remaining | -0.297373 | exclude_genus_NCBITaxon:561 | 4 |
| moderately_halophilic | fixed_baseline | -0.0511097 | exclude_family_NCBITaxon:81852 | 8 |
| moderately_halophilic | recomputed_remaining | -0.0474935 | exclude_family_NCBITaxon:81852 | 8 |

## Integration Limits

For #178, use qualified trait keys, exact scope, scenario, background/reference, effect units and correction-family fields. Keep native and shared statistics separate; do not combine ORA backgrounds in a common ranking. ORA profile-count tests are not independent corroboration of study-level effects. The native composition t intervals assume independent experiments; study-weighted summaries are supplied separately, without silently changing the native model (#133).

For #179, these are annotation influence diagnostics. Higher-rank direct assertions that cannot be allocated remain retained and flagged in the source ledger. Missing source identities/ranks are not repaired by statistical reruns. Disease/control orientation (#188), feature semantics (#189), source-count validity (#184) and full-pipeline calibration (#182) remain independent review requirements. No causal, oxygen-exposure, SCFA or taxonomy-independent mechanism follows from these results. Native kernel/ORA owner review and integration sign-off are still required before manuscript adoption.

## Review Bundle

The `.tsv.gz` comparison files contain all six traits, scenarios and scopes, including null and untestable results. `native_table1_stats.tsv.gz` is the six-trait subset of the native composition table; its q-values were computed over the FULL native families, whose sizes remain attached. The separate ORA family table includes both backgrounds and both tail-family sizes. `scenarios.tsv` supplies lineage names, support and selection status for the full genus/family inventory, including lineages outside the rule-declared panel.

The full 317-trait composition table, profile-level ORA results and source ledger are generated under `data/derived/table1_taxon_influence/`, following the repository's untracked generated-data convention. Regenerate the default input chain with `just table1-taxon-influence`, or use `just table1-taxon-influence-from <pinned-KG-directory>` for already materialized inputs. The recipe writes a candidate review bundle under the untracked run directory; this committed directory is an explicitly reviewed snapshot. This compact review bundle does not discard any six-trait influence comparison. `run_manifest.json` records hashes for the full data; `publication_manifest.json` identifies these review copies. All full-run file hashes were checked before publication.

Publication also re-sums the complete source-contribution ledger through the profile/reported-taxon links and requires exact equality with every exported profile/trait count. The publication manifest records those file hashes and row counts.
