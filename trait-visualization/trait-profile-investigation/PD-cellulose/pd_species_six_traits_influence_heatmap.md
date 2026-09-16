# PD Species-by-Six-Trait Influence Heatmap

![Heatmap](pd_species_six_traits_influence_heatmap.png)

Includes 90 species and all six lollipop traits: 60 occur in the local PD profiles and 46 have species-level DOI-cited PD evidence. The union preserves locally observed species while adding cited species that were absent from the local profiles. Each estimated trait cell is specific to the species-trait combination: the analysis removes the species' complete non-MediaDive annotation vector from profiles where it was reported, recomputes the six log-fraction scores using the unchanged human-corpus reference, forms paired increased-minus-decreased experiment contrasts, averages experiments within study, and weights studies equally.

`inc` and `dec` require at least 3 supporting studies and a 95% Student t interval entirely above or below zero. `mixed` retains supported estimates whose interval crosses zero; hatched `low support` cells do not receive a forced direction. Gray means the trait matrix has no annotation at or below that reported species label. Pale yellow marks cited species not observed in the local PD profiles, so no local influence estimate is fabricated; cross-hatching marks an exact cited species that is also absent from the local trait matrix. The first right-hand annotation retains raw increased/decreased profile occurrence counts. The second reports published PD direction and reference keys, deliberately separate from the computed trait influence.

The curated evidence table records current NCBI-mapped names alongside source-reported names, direction, evidence scope, and notes. Broader genus/family findings remain in that table as context but are not projected onto species rows. Conflicting genus-level patterns for *Prevotella*, *Streptococcus*, and *Ruminococcus* are explicitly retained rather than over-collapsed. DOI metadata and links are normalized in the source table.

This is a reported-species influence diagnostic, not a causal attribution or strain-level measurement. Removing a reported species changes both trait numerators and the nonlinear full denominator. Parent/child propagation and unresolved source assertions remain the limitations documented in #190; contrast orientation and trait semantics remain bounded by #188 and #189.

## DOI sources

- `R21` — Romano et al. Meta-analysis of the Parkinson's disease gut microbiome suggests alterations linked to intestinal inflammation. [10.1038/s41531-021-00156-z](https://doi.org/10.1038/s41531-021-00156-z)
- `W22` — Wallen et al. Metagenomics of Parkinson's disease implicates the gut microbiome in multiple disease mechanisms. [10.1038/s41467-022-34667-x](https://doi.org/10.1038/s41467-022-34667-x)
- `N24` — Nishiwaki et al. Meta-analysis of shotgun sequencing of gut microbiota in Parkinson's disease. [10.1038/s41531-024-00724-z](https://doi.org/10.1038/s41531-024-00724-z)
- `R25` — Romano et al. Machine learning-based meta-analysis reveals gut microbiome alterations associated with Parkinson's disease. [10.1038/s41467-025-56829-3](https://doi.org/10.1038/s41467-025-56829-3)
- `M26` — Menozzi et al. Microbiome signature of Parkinson's disease in healthy and genetically at-risk individuals. [10.1038/s41591-026-04318-5](https://doi.org/10.1038/s41591-026-04318-5)

[Cell statistics](pd_species_six_traits_influence_heatmap_cells.tsv) | [Taxon row annotations](pd_species_six_traits_influence_heatmap_taxa.tsv) | [Trait summary](pd_species_six_traits_influence_heatmap_traits.tsv) | [Manifest](pd_species_six_traits_influence_heatmap_manifest.json) | [Literature evidence](pd_literature_taxon_evidence.tsv) | [DOI sources](pd_literature_sources.tsv) | [Vector figure](pd_species_six_traits_influence_heatmap.svg)
