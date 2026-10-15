# Prompt 1 — Parkinson's disease: moderately halophilic and arsenate reduction

Shared by Jed Dongjin Kim-Ozaeta on 2026-10-14. Wording unchanged (quotation marks normalized to straight quotes); see [README](README.md) for status.

Figure: [`pd_condition_taxon_trait_network.png`](https://github.com/CultureBotAI/microbial-trait-profiling/blob/205257c5ec887bdf145e2c66e12f17ee005fafa0/trait-visualization/trait-profile-investigation/PD-cellulose/pd_condition_taxon_trait_network.png) (Panel B), on branch `clustered-heatmap-generation` at `205257c`.

---

Conduct a critical, publication-quality literature investigation of the following two Parkinson's disease signals in Panel B of:

trait-visualization/trait-profile-investigation/PD-cellulose/pd_condition_taxon_trait_network.png

Local results:

1. Moderately halophilic (METPO:1000623)
   - Mean increased-minus-decreased study contrast: +0.5296
   - Median contrast: +0.5440
   - Paired studies: 45
   - Positive/negative/zero studies: 33/8/4
   - Raw p = 0.0001428
   - Globally BH-adjusted q = 0.04798
   - This is the only FDR-significant condition–trait result.

2. Arsenate reduction (METPO:2000103 | pathways:arsenate_reduction)
   - Mean contrast: +0.5061
   - Median contrast: +0.1886
   - Paired studies: 45
   - Positive/negative/zero studies: 24/6/15
   - Raw p = 0.001593
   - Globally BH-adjusted q = 0.1641
   - This is exploratory and not FDR-significant.

First inspect the local methods, complete results, and study-level contrasts:

- pd_condition_taxon_trait_network_methods.md
- pd_condition_taxon_trait_network_condition_trait_tests.tsv
- pd_condition_taxon_trait_network_study_trait_contrasts.tsv
- pd_ecological_reporter_figures.md
- relevant taxon–trait and KG-Microbe audit tables in the same directory

Research questions:

1. Which exact genera, species, and strains generate each trait signal? Determine whether the two signals are driven by the same taxa.
2. Trace every relevant KG-Microbe/METPO annotation to its original source. Assess whether the evidence is strain-specific, species-level, experimentally demonstrated, predicted, or inherited through taxonomy.
3. For "moderately halophilic," distinguish true moderate halophily from salt tolerance, halotolerance, osmotic-stress resistance, and growth under ordinary physiological salinity. Evaluate whether the label is biologically credible for the implicated gut taxa.
4. For "arsenate reduction," distinguish respiratory arsenate reduction (for example arr genes) from arsenate detoxification (for example ars systems), nonspecific reductase activity, and computational pathway prediction.
5. Search human PD microbiome studies, meta-analyses, shotgun metagenomics, metatranscriptomics, metabolomics, isolate studies, and relevant animal or mechanistic experiments for support, contradiction, or null findings.
6. Investigate plausible biological contexts—including intestinal redox conditions, osmolarity and salt exposure, diet, constipation and transit time, inflammation, environmental arsenic exposure, medication, and cross-feeding—but do not present plausibility as demonstrated causality.
7. Determine whether the literature independently supports these traits or merely reports the same taxa used to construct the BugSigDB signal. Explicitly flag circular evidence.
8. Examine heterogeneity by geography, sequencing method, disease stage, constipation, levodopa and other PD medications, diet, and study design.

Interpretation constraint: these are traits annotated to taxa appearing in curated differential-abundance signatures. They are not direct measurements of pathway abundance, gene expression, arsenate reduction, salt adaptation, metabolite flux, or causal activity.

Deliver:

- A concise executive verdict for each signal: supported, partially supported, unsupported, or presently unresolvable.
- A table of contributing taxa and the provenance and specificity of each trait annotation.
- A literature evidence table with claim, population or model, assay, direction, directness, limitations, PMID, DOI, and stable link.
- Separate sections for direct human evidence, functional evidence, organism-level trait evidence, mechanisms, contradictory evidence, and confounders.
- Exact databases, search strings, search date, and inclusion criteria.
- A comparison explaining whether the two PD signals reflect one ecological pattern or two independent phenomena.
- Recommendations for analyses or experiments that could validate the signals.
- A cautious, manuscript-ready interpretation paragraph.

Prioritize primary research and systematic reviews. Verify every citation and do not invent references. Clearly distinguish evidence from inference.
