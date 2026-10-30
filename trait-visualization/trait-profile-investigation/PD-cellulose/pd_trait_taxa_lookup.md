# Taxa Associated With the PD Lollipop Traits

Source: `pd_six_trait_direction_lollipop.png` and its companion summary. The recomputed summary matches the saved source table.

Scope: 70 increased and 70 decreased PD fecal profiles. Occurrences are deduplicated by taxon and profile; names and traits come from the local KG-Microbe taxa-trait matrix.

All six lollipop traits are traced whatever their sign. Each trait's observed direction is the profile direction in which its mean z is higher, from the sign of `increased_minus_decreased_mean_z` in `pd_six_trait_direction_summary.tsv`; the tables carry it in `trait_direction`. Observed direction (increased minus decreased mean z): moderately halophilic increased (+0.495); arsenate reduction increased (+0.437); facultatively aerobic increased (+0.187); coccus shaped increased (+0.047); sphere shaped decreased (-0.026); cellulose catabolic process decreased (-0.728). A lollipop value is a descriptive difference in mean z, not evidence of statistical significance or of gene expression.

These are observed annotation carriers, not estimated causal contributions to the z-score difference. Profile z>0 describes relative trait enrichment. The matrix counts trait-bearing taxa at or below each reported taxon (src/bsdbng/taxa_trait_matrix.py, compute_trait_counts). Different descendant taxa can supply different annotations, including for species with strain-level descendants, so one reported label carrying several trait annotations does not establish that one organism has those traits together. Counts across taxonomic ranks must not be added. Raw KG edges are not present locally, so original assertion and strain provenance were not independently verified.

## moderately halophilic

Observed direction: increased (+0.495); its mean z is higher in increased than in decreased PD fecal profiles. Annotated taxon labels observed: 65 in increased profiles, 41 in decreased profiles; 85 across both directions.

Top genus/species labels by occurrence in positive-z increased profiles, the trait's observed direction:

| Taxon | Rank | Increased profiles | Decreased profiles | Increased positive-z profiles | Decreased positive-z profiles |
| --- | --- | ---: | ---: | ---: | ---: |
| Lactobacillus | genus | 9 | 4 | 5 | 4 |
| Streptococcus | genus | 7 | 3 | 4 | 2 |
| Alistipes | genus | 5 | 0 | 4 | 0 |
| Enterococcus | genus | 4 | 0 | 4 | 0 |
| Escherichia | genus | 4 | 1 | 4 | 1 |
| Klebsiella pneumoniae | species | 3 | 0 | 3 | 0 |
| Acinetobacter | genus | 2 | 0 | 2 | 0 |
| Corynebacterium | genus | 2 | 0 | 2 | 0 |
| Klebsiella | genus | 2 | 1 | 2 | 1 |
| Citrobacter | genus | 1 | 0 | 1 | 0 |
| Enterobacter | genus | 1 | 0 | 1 | 0 |
| Escherichia coli | species | 1 | 0 | 1 | 0 |

## arsenate reduction

Observed direction: increased (+0.437); its mean z is higher in increased than in decreased PD fecal profiles. Annotated taxon labels observed: 21 in increased profiles, 10 in decreased profiles; 24 across both directions.

Top genus/species labels by occurrence in positive-z increased profiles, the trait's observed direction:

| Taxon | Rank | Increased profiles | Decreased profiles | Increased positive-z profiles | Decreased positive-z profiles |
| --- | --- | ---: | ---: | ---: | ---: |
| Escherichia | genus | 4 | 1 | 4 | 1 |
| Citrobacter | genus | 1 | 0 | 1 | 0 |
| Escherichia coli | species | 1 | 0 | 1 | 0 |

## facultatively aerobic

Observed direction: increased (+0.187); its mean z is higher in increased than in decreased PD fecal profiles. Annotated taxon labels observed: 82 in increased profiles, 69 in decreased profiles; 122 across both directions.

Top genus/species labels by occurrence in positive-z increased profiles, the trait's observed direction:

| Taxon | Rank | Increased profiles | Decreased profiles | Increased positive-z profiles | Decreased positive-z profiles |
| --- | --- | ---: | ---: | ---: | ---: |
| Bifidobacterium | genus | 17 | 0 | 17 | 0 |
| Lactobacillus | genus | 9 | 4 | 9 | 3 |
| Streptococcus | genus | 7 | 3 | 7 | 3 |
| Desulfovibrio | genus | 6 | 0 | 5 | 0 |
| Bifidobacterium longum | species | 4 | 0 | 4 | 0 |
| Enterococcus | genus | 4 | 0 | 4 | 0 |
| Escherichia | genus | 4 | 1 | 4 | 1 |
| Klebsiella pneumoniae | species | 3 | 0 | 3 | 0 |
| Burkholderiales bacterium 1_1_47 | species | 4 | 0 | 2 | 0 |
| Bacteroides | genus | 3 | 6 | 2 | 2 |
| Corynebacterium | genus | 2 | 0 | 2 | 0 |
| Klebsiella | genus | 2 | 1 | 2 | 0 |

## coccus shaped

Observed direction: increased (+0.047); its mean z is higher in increased than in decreased PD fecal profiles. Annotated taxon labels observed: 92 in increased profiles, 83 in decreased profiles; 137 across both directions.

Top genus/species labels by occurrence in positive-z increased profiles, the trait's observed direction:

| Taxon | Rank | Increased profiles | Decreased profiles | Increased positive-z profiles | Decreased positive-z profiles |
| --- | --- | ---: | ---: | ---: | ---: |
| Streptococcus | genus | 7 | 3 | 7 | 3 |
| Blautia obeum | species | 5 | 0 | 5 | 0 |
| Butyricimonas | genus | 5 | 2 | 4 | 2 |
| Enterococcus | genus | 4 | 0 | 4 | 0 |
| Escherichia | genus | 4 | 1 | 4 | 1 |
| Veillonella | genus | 4 | 1 | 4 | 0 |
| Parasutterella | genus | 3 | 0 | 2 | 0 |
| Peptoniphilus | genus | 3 | 0 | 2 | 0 |
| Acidaminococcus | genus | 2 | 0 | 2 | 0 |
| Acinetobacter | genus | 2 | 0 | 2 | 0 |
| Anaerostipes | genus | 2 | 1 | 2 | 0 |
| Megasphaera | genus | 2 | 0 | 2 | 0 |

## sphere shaped

Observed direction: decreased (-0.026); its mean z is higher in decreased than in increased PD fecal profiles. Annotated taxon labels observed: 45 in increased profiles, 41 in decreased profiles; 63 across both directions.

Top genus/species labels by occurrence in positive-z decreased profiles, the trait's observed direction:

| Taxon | Rank | Increased profiles | Decreased profiles | Increased positive-z profiles | Decreased positive-z profiles |
| --- | --- | ---: | ---: | ---: | ---: |
| Ruminococcus | genus | 1 | 5 | 1 | 4 |
| Streptococcus | genus | 7 | 3 | 7 | 3 |
| Butyricicoccus | genus | 1 | 3 | 1 | 2 |
| Dialister | genus | 0 | 2 | 0 | 2 |
| Brevibacterium | genus | 0 | 1 | 0 | 1 |
| Phocaeicola | genus | 0 | 1 | 0 | 0 |
| Veillonella | genus | 4 | 1 | 4 | 0 |

## cellulose catabolic process

Observed direction: decreased (-0.728); its mean z is higher in decreased than in increased PD fecal profiles. Annotated taxon labels observed: 43 in increased profiles, 33 in decreased profiles; 53 across both directions.

Top genus/species labels by occurrence in positive-z decreased profiles, the trait's observed direction:

| Taxon | Rank | Increased profiles | Decreased profiles | Increased positive-z profiles | Decreased positive-z profiles |
| --- | --- | ---: | ---: | ---: | ---: |
| Roseburia | genus | 1 | 29 | 1 | 29 |
| Blautia | genus | 1 | 12 | 1 | 11 |
| Clostridium | genus | 2 | 9 | 1 | 9 |
| Bacteroides | genus | 3 | 6 | 1 | 6 |
| Roseburia faecis | species | 0 | 5 | 0 | 5 |
| Ruminococcus | genus | 1 | 5 | 1 | 5 |
| Lachnospira | genus | 1 | 4 | 1 | 4 |
| Butyrivibrio | genus | 1 | 3 | 1 | 3 |
| Bacteroides cellulosilyticus | species | 1 | 2 | 1 | 2 |
| Anaerostipes | genus | 2 | 1 | 0 | 1 |
| Anaerostipes hadrus | species | 0 | 1 | 0 | 1 |
| Lachnoclostridium | genus | 1 | 1 | 1 | 1 |

## Complete Evidence

- [All six trait carrier lists](pd_six_trait_taxon_carriers.tsv)
- [Taxon/profile/study evidence](pd_six_trait_taxon_profile_evidence.tsv)
- [PD fecal profile set and trait directions used](pd_six_trait_taxon_evidence_provenance.json)
