# Prompt 2 — COVID-19 and cellulose catabolic process

Shared by Jed Dongjin Kim-Ozaeta on 2026-10-14. Wording unchanged (quotation marks normalized to straight quotes); see [README](README.md) for status.

Figure: [`pd_condition_taxon_trait_network.png`](https://github.com/CultureBotAI/microbial-trait-profiling/blob/205257c5ec887bdf145e2c66e12f17ee005fafa0/trait-visualization/trait-profile-investigation/PD-cellulose/pd_condition_taxon_trait_network.png) (Panel B), on branch `clustered-heatmap-generation` at `205257c`.

---

Perform a deep, critical literature investigation of the COVID-19–cellulose catabolic process point in Panel B of:

trait-visualization/trait-profile-investigation/PD-cellulose/pd_condition_taxon_trait_network.png

Local result:

- Trait: cellulose catabolic process
- Identifier: METPO:2000103 | GO:0030245
- Mean increased-minus-decreased study contrast: −1.3070
- Median contrast: −1.2985
- Paired studies: 11
- Positive/negative/zero studies: 0/10/1
- Raw p = 0.001953
- Globally BH-adjusted q = 0.1641
- Status: exploratory, not FDR-significant

Begin by inspecting:

- pd_condition_taxon_trait_network_methods.md
- pd_condition_taxon_trait_network_condition_trait_tests.tsv
- pd_condition_taxon_trait_network_study_trait_contrasts.tsv
- relevant taxon–trait and KG-Microbe audit tables in the same directory

A negative contrast means the cellulose-catabolism annotation was more strongly represented among taxa reported as decreased than among taxa reported as increased in COVID-19 studies. It does not directly demonstrate reduced cellulose degradation or pathway activity.

Research questions:

1. Which taxa drive this signal, how consistently were they decreased, and what is the original evidence that each taxon can catabolize cellulose?
2. Does the annotation represent true cellulose degradation, growth on cellulose, carbohydrate-active enzymes, broader plant-polysaccharide utilization, or an inferred ontology assignment?
3. What do human COVID-19 fecal microbiome studies report about cellulolytic or fiber-degrading organisms, CAZymes, functional pathways, fermentation products, short-chain fatty acids, and ecological cross-feeding?
4. Do shotgun metagenomics, metatranscriptomics, metabolomics, culture studies, or longitudinal recovery studies directly support reduced cellulose-related function?
5. Separate acute infection, severe or hospitalized COVID-19, post-acute disease/long COVID, and recovery cohorts.
6. Critically evaluate antibiotics, antivirals, hospitalization, altered diet, reduced food intake, geography, age, disease severity, sampling time, stool consistency, and pre-existing metabolic disease as alternative explanations.
7. Distinguish papers already represented in the contributing BugSigDB study set from genuinely independent replication.
8. Identify contradictory and null evidence, publication bias, taxonomic-resolution problems, and cases where functional inference was made solely from 16S data.

Deliver:

- An executive verdict on whether the literature supports the direction of the signal.
- A contributing-taxon table with direction, recurrence, annotation provenance, and confidence.
- An evidence table with cohort, sample size, disease phase, method, reported result, directness, limitations, PMID, DOI, and link.
- Separate syntheses of taxonomic evidence, direct functional evidence, metabolites, mechanisms, confounding, and contradictory findings.
- Exact databases, search queries, search date, and inclusion criteria.
- Recommended validation analyses, such as CAZyme profiling, pathway-resolved metagenomics, transcriptomics, metabolomics, dietary adjustment, or longitudinal testing.
- A cautious manuscript-ready paragraph that explicitly labels q = 0.164 as exploratory.

Use primary sources and high-quality systematic reviews, verify every citation, and distinguish direct measurements from taxon-based functional predictions.
