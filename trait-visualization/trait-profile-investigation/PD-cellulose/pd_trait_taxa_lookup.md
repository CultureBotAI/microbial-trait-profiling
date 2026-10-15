# Taxa Associated With the PD Lollipop Traits

Source: `pd_six_trait_direction_lollipop.png` and its companion summary. The recomputed summary matches the saved source table.

Scope: 130 increased and 123 decreased PD fecal profiles. Occurrences are deduplicated by taxon and profile; names and traits come from the local KG-Microbe taxa-trait matrix.

The four positive-side traits are moderately halophilic, facultatively aerobic, arsenate reduction, and coccus shaped. Coccus shaped is nearly neutral (+0.017 mean-z difference). A positive lollipop value is not evidence of statistical significance or gene expression. Cellulose breakdown is the strongest negative-side trait.

These are observed annotation carriers, not estimated causal contributions to the z-score difference. Profile z>0 describes relative trait enrichment. The matrix counts trait-bearing taxa at or below each reported taxon (src/bsdbng/taxa_trait_matrix.py, compute_trait_counts). Different descendant taxa can supply different annotations, including for species with strain-level descendants. An all-four matrix intersection therefore does not establish co-occurrence of four traits in one organism. Counts across taxonomic ranks must not be added. Raw KG edges are not present locally, so original assertion and strain provenance were not independently verified.

## facultatively aerobic

Mean-z difference (increased - decreased): +0.452. 110 annotated taxon labels observed in increased profiles; 157 across both directions.

Top genus/species labels by occurrence in positive-z increased profiles:

| Taxon | Rank | Increased profiles | Decreased profiles | Increased positive-z profiles | Increased positive-z studies |
| --- | --- | ---: | ---: | ---: | ---: |
| Bifidobacterium | genus | 21 | 2 | 19 | 10 |
| Lactobacillus | genus | 15 | 6 | 14 | 11 |
| Streptococcus | genus | 11 | 5 | 11 | 10 |
| Desulfovibrio | genus | 11 | 3 | 9 | 5 |
| Escherichia | genus | 9 | 2 | 9 | 5 |
| Enterococcus | genus | 7 | 2 | 7 | 7 |
| Klebsiella | genus | 6 | 1 | 6 | 4 |
| Serratia | genus | 6 | 0 | 6 | 3 |
| Streptococcus mutans | species | 5 | 0 | 5 | 2 |
| Bifidobacterium longum | species | 5 | 0 | 4 | 3 |
| Actinomyces sp. oral taxon 448 | species | 4 | 0 | 4 | 1 |
| Enterobacter | genus | 4 | 0 | 4 | 3 |

## moderately halophilic

Mean-z difference (increased - decreased): +0.537. 80 annotated taxon labels observed in increased profiles; 103 across both directions.

Top genus/species labels by occurrence in positive-z increased profiles:

| Taxon | Rank | Increased profiles | Decreased profiles | Increased positive-z profiles | Increased positive-z studies |
| --- | --- | ---: | ---: | ---: | ---: |
| Lactobacillus | genus | 15 | 6 | 11 | 10 |
| Escherichia | genus | 9 | 2 | 9 | 5 |
| Streptococcus | genus | 11 | 5 | 8 | 8 |
| Alistipes | genus | 8 | 0 | 8 | 8 |
| Enterococcus | genus | 7 | 2 | 7 | 7 |
| Klebsiella | genus | 6 | 1 | 6 | 4 |
| Serratia | genus | 6 | 0 | 6 | 3 |
| Streptococcus mutans | species | 5 | 0 | 5 | 2 |
| Enterobacter | genus | 4 | 0 | 4 | 3 |
| Klebsiella pneumoniae | species | 4 | 0 | 4 | 3 |
| Acinetobacter | genus | 2 | 0 | 2 | 2 |
| Bacillus | genus | 2 | 0 | 2 | 1 |

## coccus shaped

Mean-z difference (increased - decreased): +0.017. 114 annotated taxon labels observed in increased profiles; 163 across both directions.

Top genus/species labels by occurrence in positive-z increased profiles:

| Taxon | Rank | Increased profiles | Decreased profiles | Increased positive-z profiles | Increased positive-z studies |
| --- | --- | ---: | ---: | ---: | ---: |
| Streptococcus | genus | 11 | 5 | 11 | 10 |
| Peptoniphilus | genus | 9 | 0 | 9 | 1 |
| Enterococcus | genus | 7 | 2 | 5 | 5 |
| Blautia obeum | species | 5 | 0 | 5 | 2 |
| Streptococcus mutans | species | 5 | 0 | 5 | 2 |
| Escherichia | genus | 9 | 2 | 4 | 3 |
| Veillonella | genus | 5 | 2 | 4 | 3 |
| Streptococcus salivarius | species | 4 | 0 | 4 | 4 |
| Streptococcus vestibularis | species | 3 | 0 | 3 | 1 |
| Turicibacter | genus | 3 | 1 | 3 | 1 |
| Mediterraneibacter torques | species | 3 | 1 | 2 | 2 |
| Methanobrevibacter smithii | species | 3 | 0 | 2 | 2 |

## sphere shaped

Mean-z difference (increased - decreased): -0.083. 54 annotated taxon labels observed in increased profiles; 74 across both directions.

Top genus/species labels by occurrence in positive-z increased profiles:

| Taxon | Rank | Increased profiles | Decreased profiles | Increased positive-z profiles | Increased positive-z studies |
| --- | --- | ---: | ---: | ---: | ---: |
| Streptococcus | genus | 11 | 5 | 11 | 10 |
| Peptoniphilus | genus | 9 | 0 | 9 | 1 |
| Enterococcus | genus | 7 | 2 | 6 | 6 |
| Veillonella | genus | 5 | 2 | 5 | 4 |
| Acidaminococcus | genus | 2 | 0 | 2 | 2 |
| Acinetobacter | genus | 2 | 0 | 2 | 2 |
| Lactococcus | genus | 2 | 3 | 2 | 2 |
| Megasphaera | genus | 2 | 0 | 2 | 2 |
| Ruminococcus | genus | 2 | 9 | 2 | 2 |
| Butyricicoccus | genus | 1 | 6 | 1 | 1 |
| Streptococcus anginosus | species | 1 | 0 | 1 | 1 |
| Victivallis | genus | 1 | 0 | 1 | 1 |

## cellulose catabolic process

Mean-z difference (increased - decreased): -0.609. 50 annotated taxon labels observed in increased profiles; 60 across both directions.

Top genus/species labels by occurrence in positive-z increased profiles:

| Taxon | Rank | Increased profiles | Decreased profiles | Increased positive-z profiles | Increased positive-z studies |
| --- | --- | ---: | ---: | ---: | ---: |
| Blautia obeum | species | 5 | 0 | 5 | 2 |
| Roseburia | genus | 4 | 31 | 4 | 4 |
| Bacteroides | genus | 3 | 6 | 3 | 2 |
| Clostridium | genus | 3 | 10 | 2 | 2 |
| Blautia | genus | 2 | 15 | 2 | 2 |
| Butyrivibrio | genus | 2 | 3 | 2 | 2 |
| Lachnospira | genus | 2 | 9 | 2 | 2 |
| Ruminococcus | genus | 2 | 9 | 2 | 2 |
| Mediterraneibacter torques | species | 3 | 1 | 1 | 1 |
| Bacillus | genus | 2 | 0 | 1 | 1 |
| Ruminococcus bromii | species | 2 | 2 | 1 | 1 |
| Bacteroides cellulosilyticus | species | 1 | 2 | 1 | 1 |

## arsenate reduction

Mean-z difference (increased - decreased): +0.391. 23 annotated taxon labels observed in increased profiles; 27 across both directions.

Top genus/species labels by occurrence in positive-z increased profiles:

| Taxon | Rank | Increased profiles | Decreased profiles | Increased positive-z profiles | Increased positive-z studies |
| --- | --- | ---: | ---: | ---: | ---: |
| Escherichia | genus | 9 | 2 | 9 | 5 |
| Citrobacter | genus | 2 | 0 | 2 | 2 |
| Escherichia coli | species | 1 | 2 | 1 | 1 |

## Taxa Mapped to All Four Positive-Side Annotations

Annotation intersection is separate from all four traits having z>0 in a profile.

| Taxon | Rank | Increased profiles | Decreased profiles | Increased profiles with all four z>0 |
| --- | --- | ---: | ---: | ---: |
| Enterobacteriaceae | family | 18 | 2 | 5 |
| Pseudomonadota | phylum | 15 | 1 | 4 |
| Escherichia | genus | 9 | 2 | 4 |
| Clostridia | class | 8 | 6 | 5 |
| Actinomycetota | phylum | 8 | 0 | 2 |
| Eubacteriales | order | 7 | 7 | 1 |
| Bacillota | phylum | 5 | 14 | 5 |
| Bacilli | class | 4 | 3 | 3 |
| Enterobacterales | order | 3 | 0 | 0 |
| Gammaproteobacteria | class | 3 | 1 | 1 |
| Citrobacter | genus | 2 | 0 | 1 |
| Burkholderiales | order | 2 | 0 | 1 |
| Betaproteobacteria | class | 1 | 0 | 0 |
| Actinomycetes | class | 1 | 0 | 1 |
| Alphaproteobacteria | class | 1 | 0 | 1 |
| Microbacteriaceae | family | 1 | 0 | 1 |
| Pseudomonadales | order | 1 | 0 | 1 |
| Bacillaceae | family | 1 | 1 | 0 |
| Escherichia coli | species | 1 | 2 | 1 |
| Micrococcaceae | family | 0 | 2 | 0 |
| Bacillales | order | 0 | 2 | 0 |
| Comamonadaceae | family | 0 | 2 | 0 |
| Staphylococcaceae | family | 0 | 1 | 0 |

## Complete Evidence

- [All six trait carrier lists](pd_six_trait_taxon_carriers.tsv)
- [Taxon/profile/study evidence](pd_six_trait_taxon_profile_evidence.tsv)
- [All-four annotation intersection](pd_four_positive_trait_shared_taxa.tsv)
