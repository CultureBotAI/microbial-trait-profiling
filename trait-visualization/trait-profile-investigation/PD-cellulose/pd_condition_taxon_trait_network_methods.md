# FDR-tiered condition-taxon-trait evidence landscape

## Statistical universe

- Human fecal conditions with at least `3` studies containing both increased and decreased signatures: `66`.
- Exact genus/species taxa required reporting in at least `3` paired studies and at least one of the six KG-Microbe trait annotations.
- The resulting condition-taxon test family contained `288` eligible comparisons (`288` with finite nominal p), and the condition-trait test family contained `396` eligible condition-trait comparisons (`392` with finite nominal p; `6` PD trait comparisons shown in Panel C).
- Taxon direction was tested with an exact McNemar/binomial test on discordant paired studies.
- Trait contrasts were tested across paired study-level increased-minus-decreased z scores with a two-sided Wilcoxon signed-rank test.
- Reported `p_value` columns and Panel C `nominal p` labels are nominal (uncorrected) across `2688` finite nominal p-values in the condition-taxon table (`288` support-eligible with finite nominal p in the `288`-row BH family) and `392` finite nominal p-values in the condition-trait table (`392` support-eligible with finite nominal p in the `396`-row BH family); `q_value` is Benjamini-Hochberg adjusted globally within the complete condition-taxon family and within the complete condition-trait family (holding untestable support rows at `p = 1.0`).
- In the edge export, `displayed_p_value_count` counts finite p-values in exported edges separately by inferential edge class; `n_tests` is the number of finite nominal p-values across every row of the source table, including rows outside support eligibility that are not in the BH family, and `bh_family_size` preserves its support-eligible padded BH family. Annotation edges have no inferential counts.
- FDR-significant results require `q <= 0.05`.
- Results with `0.05 < q <= 0.20` are shown as exploratory, not significant.

## Results shown

- Eligible conditions represented in the evidence landscape: `66`.
- FDR-selected condition-taxon results: `0` involving `0` genera and `0` species.
- FDR-selected condition-trait results: `0`.
- Exploratory condition-taxon results: `1`.
- Exploratory condition-trait results: `6`.
- Conditions in the significant-only node/edge tables: `0`.
- Dashed KG-Microbe annotation links: `0`.
- Network nodes/edges: `0` / `0`.

## Interpretation boundary

Filled points are the only results called FDR-significant. Hollow amber points expose potentially useful follow-up hypotheses but are not discoveries. Taxon tests measure directional asymmetry of reporting in curated BugSigDB signatures; they are not abundance-effect meta-analyses because harmonized study effect sizes and standard errors are unavailable. KG-Microbe annotations do not establish activity, expression, strain-level co-occurrence, or causality.

[All condition-taxon tests](pd_condition_taxon_trait_network_condition_taxon_tests.tsv) | [All condition-trait tests](pd_condition_taxon_trait_network_condition_trait_tests.tsv) | [Condition eligibility](pd_condition_taxon_trait_network_condition_eligibility.tsv) | [Selected taxa](pd_condition_taxon_trait_network_taxon_selection.tsv) | [Edges](pd_condition_taxon_trait_network_edges.tsv) | [Nodes](pd_condition_taxon_trait_network_nodes.tsv) | [Vector figure](pd_condition_taxon_trait_network.svg)
