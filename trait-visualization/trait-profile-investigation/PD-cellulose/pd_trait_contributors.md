# PD Trait Contributor Dot Heatmap

![Contributor heatmap](pd_trait_contributors_species_family_dot_heatmap.png)

The refreshed heatmap keeps the species/family contributor view and adds a PD specificity strip beside each taxon. The strip reports exact reported taxon-label profiles in PD feces versus other feces, so it checks whether a contributor is specific to the PD fecal setting rather than simply common across the fecal corpus.

Direction 5 is added as a trait co-carriage network:

![Trait co-carriage network](pd_trait_contributors_trait_cocarriage_network.png)

Top PD-skewed reported taxon labels:

| rank | taxon | PD/other feces profiles | log2 ratio |
| --- | --- | ---: | ---: |
| family | Erwiniaceae (1903409) | 1/0 | 1.58 |
| family | Listeriaceae (186820) | 1/1 | 0.00 |
| species | Streptococcus mutans (1309) | 5/5 | 0.00 |
| species | Roseburia faecis (301302) | 6/18 | -1.51 |
| species | Listeria monocytogenes (1639) | 0/1 | -1.58 |

Largest trait co-carriage pairs among selected contributors:

| trait pair | shared taxa | examples |
| --- | ---: | --- |
| facultatively aerobic + coccus shaped | 28 | Bacillaceae (186817); Micrococcaceae (1268); Staphylococcaceae (90964); Enterobacteriaceae (543); Streptococcaceae (1300) |
| facultatively aerobic + moderately halophilic | 28 | Bacillaceae (186817); Micrococcaceae (1268); Staphylococcaceae (90964); Enterobacteriaceae (543); Streptococcaceae (1300) |
| coccus shaped + sphere shaped | 21 | Bacillaceae (186817); Micrococcaceae (1268); Staphylococcaceae (90964); Enterobacteriaceae (543); Streptococcaceae (1300) |
| coccus shaped + moderately halophilic | 19 | Bacillaceae (186817); Micrococcaceae (1268); Staphylococcaceae (90964); Enterobacteriaceae (543); Streptococcaceae (1300) |
| facultatively aerobic + sphere shaped | 16 | Bacillaceae (186817); Micrococcaceae (1268); Staphylococcaceae (90964); Enterobacteriaceae (543); Streptococcaceae (1300) |

Audit files:

- `pd_trait_contributors_cells.tsv`: long-form dot heatmap cells.
- `pd_trait_contributors_selected_taxa.tsv`: selected taxa with specificity and co-carriage counts.
- `pd_trait_contributors_taxon_specificity.tsv`: PD-vs-other feces specificity table.
- `pd_trait_contributors_trait_cocarriage_taxa.tsv`: taxon-by-trait presence table.
- `pd_trait_contributors_trait_cocarriage_edges.tsv`: pairwise trait sharing table.
- `pd_trait_contributors_species_family_dot_heatmap.svg`: vector contributor heatmap.
- `pd_trait_contributors_trait_cocarriage_network.svg`: vector co-carriage network.
