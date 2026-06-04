# Gut Anaerobe / Fermentation Depletion Investigation

This investigation uses the existing study trait profile clustered heatmap matrix, study metadata, taxa guild heatmap metadata, and derived taxa/trait matrices.

## Scope

- Guild traits: anaerobic (METPO:1000603), obligately anaerobic (METPO:1000607), Fermentation (METPO:1002005), cellulose catabolic process (METPO:2000103|GO:0030245), xylan catabolic process (METPO:2000103|GO:0045493)
- Target taxa: Faecalibacterium prausnitzii, Anaerobutyricum hallii, Lachnospira eligens, Agathobacter rectalis, Roseburia intestinalis, Ruminococcus bromii, Roseburia inulinivorans
- Gut context: body-site labels matching feces/faeces/stool/colon/caecum/cecum/gut/intestin/digestive, plus matching grouped body-site labels.

## Main Signal

The selected guild traits are consistently higher in decreased profiles than in increased profiles, which supports the gut anaerobe/fermentation depletion hypothesis.

| trait | decreased - increased z |
| --- | ---: |
| cellulose catabolic process | 0.132 |
| xylan catabolic process | 0.103 |
| anaerobic | 0.077 |
| obligately anaerobic | 0.074 |
| Fermentation | 0.059 |

## Strongest Depletion-Biased Taxa

| taxon | decreased / increased profiles | depletion bias |
| --- | ---: | ---: |
| Agathobacter rectalis | 92 / 33 | 0.472 |
| Lachnospira eligens | 101 / 39 | 0.443 |
| Roseburia intestinalis | 76 / 31 | 0.421 |
| Roseburia inulinivorans | 62 / 29 | 0.363 |
| Ruminococcus bromii | 61 / 38 | 0.232 |
| Faecalibacterium prausnitzii | 177 / 111 | 0.229 |
| Anaerobutyricum hallii | 86 / 63 | 0.154 |

## Recurrent Target-Taxon Depletions By Condition

| taxon | condition | body site | decreased / increased profiles | decreased studies |
| --- | --- | --- | ---: | ---: |
| Faecalibacterium prausnitzii | Parkinson's disease | Feces | 16 / 0 | 11 |
| Roseburia intestinalis | Parkinson's disease | Feces | 12 / 0 | 7 |
| Faecalibacterium prausnitzii | COVID-19 | Feces | 11 / 0 | 6 |
| Lachnospira eligens | Colorectal cancer | Feces | 9 / 0 | 5 |
| Ruminococcus bromii | COVID-19 | Feces | 9 / 0 | 5 |
| Agathobacter rectalis | COVID-19 | Feces | 7 / 0 | 5 |
| Anaerobutyricum hallii | Colorectal cancer | Feces | 5 / 0 | 5 |
| Roseburia intestinalis | Colorectal cancer | Feces | 5 / 0 | 5 |
| Agathobacter rectalis | Colorectal cancer | Feces | 4 / 0 | 3 |
| Agathobacter rectalis | Parkinson's disease | Feces | 3 / 0 | 3 |

## Gut Contexts With High Decreased-Profile Guild Signal

| condition | profiles | studies | decreased mean guild z | decreased - increased z |
| --- | ---: | ---: | ---: | ---: |
| Graft versus host disease | 13 | 4 | 1.199 | 1.141 |
| Acute graft vs. host disease | 11 | 2 | 1.063 | 0.915 |
| Colorectal carcinoma | 15 | 3 | 0.961 | 1.286 |
| Autoimmune disease | 8 | 2 | 0.922 | 0.732 |
| Age-related macular degeneration | 10 | 3 | 0.913 | 0.729 |
| Inflammatory bowel disease | 42 | 5 | 0.751 | 0.639 |
| Lewy body dementia | 23 | 2 | 0.723 | 0.427 |
| Attention deficit hyperactivity disorder | 6 | 3 | 0.658 | 0.719 |

## Body-Site Check

| body_site_group | profiles | decreased mean guild z | decreased - increased z |
| --- | ---: | ---: | ---: |
| Colon | 139 | 0.362 | 0.331 |
| Feces | 7595 | 0.228 | 0.138 |
| Caecum | 165 | 0.155 | 0.052 |
| Other | 391 | -0.119 | -0.039 |

## Top Decreased Profiles To Inspect

| profile | condition | body site | guild z | title |
| --- | --- | --- | ---: | --- |
| bsdb:38424099\|bsdb:38424099-6\|dec | Diet | Feces | 2.602 | Effect of gut microbiome modulation on muscle function and cognition: the PROMOTe rando... |
| bsdb:36721179\|bsdb:36721179-2\|dec | Response to transplant | Feces | 2.340 | Gut microbiota-derived metabolites mediate the neuroprotective effect of melatonin in c... |
| bsdb:41270737\|bsdb:41270737-7\|dec | Treatment outcome measurement | Feces | 2.319 | Interindividual variability in gut microbiome mediates the efficacy of resistant starch... |
| bsdb:41378921\|bsdb:41378921-1\|dec | Restless legs syndrome | Feces | 2.253 | Analysis of gut microbiota in Restless Legs Syndrome: searching for a metagenomic signa... |
| bsdb:23274889\|bsdb:23274889-1\|dec | Autoimmune disease | Feces | 2.224 | Fecal microbiota composition differs between children with β-cell autoimmunity and thos... |
| bsdb:40441865\|bsdb:40441865-3\|dec | Acute graft vs. host disease | Feces | 2.203 | Commensal <i>Bacteroides</i> T6SS alleviate GI-aGVHD via mediating gut microbiota compo... |
| bsdb:40308807\|bsdb:40308807-3\|dec | Chronic hepatitis B virus infection | Feces | 2.203 | Machine learning algorithms reveal gut microbiota signatures associated with chronic he... |
| bsdb:35985814\|bsdb:35985814-3\|dec | Clinical treatment | Feces | 2.203 | Inflammation in Children with CKD Linked to Gut Dysbiosis and Metabolite Imbalance |
| bsdb:34852831\|bsdb:34852831-1\|dec | Response to transplant | Feces | 2.203 | Effect of fecal microbiota transplantation in patients with slow transit constipation a... |
| bsdb:33071933\|bsdb:33071933-3\|dec | Parkinson's disease | Feces | 2.203 | Changes in Gastrointestinal Microbiome Composition in PD: A Pivotal Role of Covariates |

## Suggested Follow-Up

1. Treat the target taxa as a depletion guild and test enrichment by condition/body site rather than analyzing each taxon in isolation.
2. Prioritize conditions where decreased profiles have high guild scores and target taxa are repeatedly decreased across multiple studies.
3. Re-run the same summaries after adding explicit SCFA/butyrate production annotations, because the current selected heatmap traits capture fermentation and fiber degradation but do not fully encode butyrate production.

## Files

- `guild_trait_direction_shift.tsv`
- `gut_body_site_guild_summary.tsv`
- `gut_condition_guild_summary.tsv`
- `candidate_taxa_direction_summary.tsv`
- `candidate_taxa_trait_presence.tsv`
- `candidate_taxa_condition_direction.tsv`
- `candidate_taxa_profile_occurrences.tsv`
- `top_decreased_guild_profiles.tsv`
- `figure_methods.md`
- `candidate_taxa_depletion_bias.png`
- `condition_decreased_guild_signal.png`
- `guild_trait_direction_shift.png`
- `anatomical_sites_diagram.png`
- `target_taxa_condition_breakdown_bar.png`
