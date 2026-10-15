# Prompt 3 — Colorectal cancer and cellulose catabolic process

Shared by Jed Dongjin Kim-Ozaeta on 2026-10-14. Saved verbatim; see [README](README.md) for status.

Figure: [`pd_condition_taxon_trait_network.png`](https://github.com/CultureBotAI/microbial-trait-profiling/blob/205257c5ec887bdf145e2c66e12f17ee005fafa0/trait-visualization/trait-profile-investigation/PD-cellulose/pd_condition_taxon_trait_network.png) (Panel B), on branch `clustered-heatmap-generation` at `205257c`.

---

Perform a publication-quality literature investigation of the colorectal cancer–cellulose catabolic process point in Panel B of:

trait-visualization/trait-profile-investigation/PD-cellulose/pd_condition_taxon_trait_network.png

Local result:

- Trait: cellulose catabolic process
- Identifier: METPO:2000103 | GO:0030245
- Mean increased-minus-decreased study contrast: −0.8206
- Median contrast: −0.8155
- Paired studies: 24
- Positive/negative/zero studies: 3/17/4
- Raw p = 0.001944
- Globally BH-adjusted q = 0.1641
- Status: exploratory, not FDR-significant

First inspect the local methods, complete result table, study-level contrasts, and taxon–trait provenance:

- pd_condition_taxon_trait_network_methods.md
- pd_condition_taxon_trait_network_condition_trait_tests.tsv
- pd_condition_taxon_trait_network_study_trait_contrasts.tsv
- relevant KG-Microbe and taxon–trait audit tables in the same directory

Interpret the negative contrast precisely: taxa annotated with cellulose catabolism were more strongly represented among taxa reported as decreased than among taxa reported as increased in colorectal cancer. This is not a direct measurement of cellulose degradation, gene expression, metabolite production, or causal protection.

Research questions:

1. Identify the taxa driving the signal and verify the experimental or computational basis of their cellulose-catabolism annotations.
2. Distinguish true cellulose degradation from utilization of soluble fiber, resistant starch, hemicellulose, other plant polysaccharides, or downstream cross-feeding products.
3. Determine whether independent colorectal-cancer studies consistently report depletion of these organisms or functions.
4. Review shotgun metagenomic, metatranscriptomic, metabolomic, culture-based, and dietary evidence concerning cellulose/fiber degradation, CAZymes, butyrate and other fermentation products, mucosal integrity, inflammation, and tumor-associated metabolism.
5. Separate potentially protective fiber-degrading ecology from simplistic claims that cellulose degradation itself prevents cancer.
6. Examine differences between fecal and mucosal samples, adenoma and carcinoma, tumor stage and location, treatment-naive and treated patients, screening and symptomatic cohorts, and longitudinal versus cross-sectional designs.
7. Evaluate diet, bowel preparation, stool transit, occult bleeding, antibiotics, geography, age, obesity, and sequencing methodology as possible confounders.
8. Distinguish source studies already included in the BugSigDB-derived signal from independent validation, and identify contradictory or null results.
9. Assess whether a few recurrent taxa dominate the result and whether species- or strain-level variation undermines genus-level interpretation.

Deliver:

- An executive verdict: supported, partially supported, unsupported, or unresolvable.
- A contributing-taxon and annotation-provenance table.
- A literature evidence matrix containing cohort, sample type, cancer stage, assay, finding, direction, directness, limitations, PMID, DOI, and stable link.
- Separate sections for taxonomic observations, direct functional measurements, metabolites and mechanisms, conflicting evidence, and confounders.
- Exact databases, search strings, search date, and selection criteria.
- Recommendations for validation using CAZyme profiles, pathway-resolved metagenomics, expression data, metabolites, diet, and clinical covariates.
- A cautious manuscript-ready paragraph that treats q = 0.164 as hypothesis-generating rather than significant.

Prioritize primary studies and systematic reviews. Verify every citation, avoid causal overstatement, and explicitly distinguish measured microbial function from function inferred through taxonomic annotation.
