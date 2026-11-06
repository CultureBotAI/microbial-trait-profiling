# Taxa Associated With the PD Lollipop Traits

Source: `pd_six_trait_direction_lollipop.png` and its companion summary. The recomputed summary matches the saved source table.

Scope: 126 increased and 116 decreased PD fecal profiles. Occurrences are deduplicated by taxon and profile; names and traits come from the local KG-Microbe taxa-trait matrix.

All six lollipop traits are traced whatever their sign. Each trait's observed direction is the profile direction in which its mean z is higher, from the sign of `increased_minus_decreased_mean_z` in `pd_six_trait_direction_summary.tsv`; the tables carry it in `trait_direction`. Observed direction (increased minus decreased mean z): moderately halophilic increased (+0.509); arsenate reduction increased (+0.390); facultatively aerobic increased (+0.348); coccus shaped decreased (-0.076); sphere shaped decreased (-0.128); cellulose catabolic process decreased (-0.711). A lollipop value is a descriptive difference in mean z, not evidence of statistical significance or of gene expression.

These are observed annotation carriers, not estimated causal contributions to the z-score difference. Profile z>0 describes relative trait enrichment. The matrix counts trait-bearing taxa at or below each reported taxon (src/bsdbng/taxa_trait_matrix.py, compute_trait_counts). Different descendant taxa can supply different annotations, including for species with strain-level descendants, so one reported label carrying several trait annotations does not establish that one organism has those traits together. Counts across taxonomic ranks must not be added. Raw KG edges are not present locally, so original assertion and strain provenance were not independently verified.

## moderately halophilic

Observed direction: increased (+0.509); its mean z is higher in increased than in decreased PD fecal profiles. Annotated taxon labels observed: 79 in increased profiles, 53 in decreased profiles; 100 across both directions.

Top genus/species labels by occurrence in positive-z increased profiles, the trait's observed direction:

| Taxon | Rank | Increased profiles | Decreased profiles | Increased positive-z profiles | Decreased positive-z profiles |
| --- | --- | ---: | ---: | ---: | ---: |
| Lactobacillus | genus | 15 | 6 | 11 | 5 |
| Escherichia | genus | 9 | 2 | 9 | 2 |
| Alistipes | genus | 8 | 0 | 8 | 0 |
| Streptococcus | genus | 10 | 5 | 7 | 2 |
| Enterococcus | genus | 7 | 1 | 7 | 1 |
| Klebsiella | genus | 6 | 1 | 6 | 1 |
| Serratia | genus | 6 | 0 | 6 | 0 |
| Enterobacter | genus | 4 | 0 | 4 | 0 |
| Klebsiella pneumoniae | species | 4 | 0 | 4 | 0 |
| Acinetobacter | genus | 2 | 0 | 2 | 0 |
| Bacillus | genus | 2 | 0 | 2 | 0 |
| Citrobacter | genus | 2 | 0 | 2 | 0 |

## arsenate reduction

Observed direction: increased (+0.390); its mean z is higher in increased than in decreased PD fecal profiles. Annotated taxon labels observed: 23 in increased profiles, 13 in decreased profiles; 26 across both directions.

Top genus/species labels by occurrence in positive-z increased profiles, the trait's observed direction:

| Taxon | Rank | Increased profiles | Decreased profiles | Increased positive-z profiles | Decreased positive-z profiles |
| --- | --- | ---: | ---: | ---: | ---: |
| Escherichia | genus | 9 | 2 | 9 | 2 |
| Citrobacter | genus | 2 | 0 | 2 | 0 |
| Escherichia coli | species | 1 | 2 | 1 | 2 |

## facultatively aerobic

Observed direction: increased (+0.348); its mean z is higher in increased than in decreased PD fecal profiles. Annotated taxon labels observed: 105 in increased profiles, 89 in decreased profiles; 149 across both directions.

Top genus/species labels by occurrence in positive-z increased profiles, the trait's observed direction:

| Taxon | Rank | Increased profiles | Decreased profiles | Increased positive-z profiles | Decreased positive-z profiles |
| --- | --- | ---: | ---: | ---: | ---: |
| Bifidobacterium | genus | 24 | 2 | 18 | 1 |
| Lactobacillus | genus | 15 | 6 | 14 | 5 |
| Streptococcus | genus | 10 | 5 | 10 | 5 |
| Desulfovibrio | genus | 11 | 3 | 9 | 0 |
| Escherichia | genus | 9 | 2 | 9 | 2 |
| Enterococcus | genus | 7 | 1 | 7 | 1 |
| Klebsiella | genus | 6 | 1 | 6 | 0 |
| Serratia | genus | 6 | 0 | 6 | 0 |
| Bifidobacterium longum | species | 6 | 0 | 4 | 0 |
| Enterobacter | genus | 4 | 0 | 4 | 0 |
| Klebsiella pneumoniae | species | 4 | 0 | 4 | 0 |
| Burkholderiales bacterium 1_1_47 | species | 4 | 0 | 3 | 0 |

## coccus shaped

Observed direction: decreased (-0.076); its mean z is higher in decreased than in increased PD fecal profiles. Annotated taxon labels observed: 108 in increased profiles, 94 in decreased profiles; 149 across both directions.

Top genus/species labels by occurrence in positive-z decreased profiles, the trait's observed direction:

| Taxon | Rank | Increased profiles | Decreased profiles | Increased positive-z profiles | Decreased positive-z profiles |
| --- | --- | ---: | ---: | ---: | ---: |
| Blautia | genus | 2 | 16 | 0 | 14 |
| Ruminococcus | genus | 2 | 9 | 1 | 7 |
| Coprococcus | genus | 1 | 11 | 1 | 6 |
| Butyricicoccus | genus | 1 | 6 | 1 | 5 |
| Streptococcus | genus | 10 | 5 | 10 | 5 |
| Blautia hansenii | species | 0 | 3 | 0 | 3 |
| Actinomyces | genus | 1 | 3 | 1 | 2 |
| Brevibacterium | genus | 0 | 3 | 0 | 2 |
| Dialister | genus | 0 | 2 | 0 | 2 |
| Gemella | genus | 0 | 2 | 0 | 2 |
| Granulicatella | genus | 1 | 2 | 1 | 2 |
| Ruminococcus gauvreauii | species | 0 | 2 | 0 | 2 |

## sphere shaped

Observed direction: decreased (-0.128); its mean z is higher in decreased than in increased PD fecal profiles. Annotated taxon labels observed: 53 in increased profiles, 49 in decreased profiles; 72 across both directions.

Top genus/species labels by occurrence in positive-z decreased profiles, the trait's observed direction:

| Taxon | Rank | Increased profiles | Decreased profiles | Increased positive-z profiles | Decreased positive-z profiles |
| --- | --- | ---: | ---: | ---: | ---: |
| Ruminococcus | genus | 2 | 9 | 2 | 7 |
| Butyricicoccus | genus | 1 | 6 | 1 | 6 |
| Streptococcus | genus | 10 | 5 | 10 | 5 |
| Brevibacterium | genus | 0 | 3 | 0 | 3 |
| Dialister | genus | 0 | 2 | 0 | 2 |
| Ruminococcus gauvreauii | species | 0 | 2 | 0 | 2 |
| Enterococcus | genus | 7 | 1 | 5 | 1 |
| Howardella | genus | 0 | 1 | 0 | 1 |
| Lactococcus | genus | 2 | 2 | 2 | 0 |
| Veillonella | genus | 5 | 2 | 5 | 0 |
| Phocaeicola | genus | 0 | 1 | 0 | 0 |

## cellulose catabolic process

Observed direction: decreased (-0.711); its mean z is higher in decreased than in increased PD fecal profiles. Annotated taxon labels observed: 49 in increased profiles, 38 in decreased profiles; 58 across both directions.

Top genus/species labels by occurrence in positive-z decreased profiles, the trait's observed direction:

| Taxon | Rank | Increased profiles | Decreased profiles | Increased positive-z profiles | Decreased positive-z profiles |
| --- | --- | ---: | ---: | ---: | ---: |
| Roseburia | genus | 4 | 33 | 4 | 33 |
| Blautia | genus | 2 | 16 | 2 | 16 |
| Clostridium | genus | 3 | 10 | 2 | 10 |
| Lachnospira | genus | 2 | 9 | 2 | 9 |
| Ruminococcus | genus | 2 | 9 | 2 | 9 |
| Bacteroides | genus | 4 | 6 | 3 | 6 |
| Roseburia faecis | species | 0 | 6 | 0 | 6 |
| Butyrivibrio | genus | 2 | 3 | 2 | 3 |
| Anaerostipes | genus | 2 | 2 | 0 | 2 |
| Bacteroides cellulosilyticus | species | 1 | 2 | 1 | 2 |
| Ruminococcus bromii | species | 2 | 2 | 1 | 2 |
| Anaerostipes hadrus | species | 0 | 1 | 0 | 1 |

## Complete Evidence

- [All six trait carrier lists](pd_six_trait_taxon_carriers.tsv)
- [Taxon/profile/study evidence](pd_six_trait_taxon_profile_evidence.tsv)
- [PD fecal profile set and trait directions used](pd_six_trait_taxon_evidence_provenance.json)
