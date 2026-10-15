# Autism and Autoimmune Disease Follow-Up: Aerobe v2

This follow-up starts from aerobe v2 trait enrichment in increased profiles and inspects the BugSigDB taxa reported in the autism/autism spectrum disorder and autoimmune disease condition labels.

## Plot Guild Specificity

| plot | guild-specific? | basis |
| --- | --- | --- |
| `autism_autoimmune_taxa_signed_heatmap.png` | no | taxon direction-count heatmap; selected traits are annotation context |

## Included Condition Labels

| condition group | exact condition label | profiles | studies | increased / decreased profiles | increased - decreased trait z | body-site groups |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| Autism / autism spectrum disorder | Autism | 14 | 4 | 8 / 6 | 0.749 | Feces |
| Autism / autism spectrum disorder | Autism spectrum disorder | 13 | 3 | 8 / 5 | -0.031 | Feces |
| Autoimmune disease | Autoimmune disease | 8 | 2 | 4 / 4 | 0.419 | Feces |

## Top Taxa

The tables below count unique taxon/profile occurrences after deduplicating the same taxon within the same study/experiment/direction profile.

### Autism / autism spectrum disorder

Top increased taxa:

| taxon | rank | profiles | studies | selected traits |
| --- | --- | ---: | ---: | --- |
| Parabacteroides | genus | 4 | 3 |  |
| Clostridium | genus | 3 | 3 | microaerophilic |
| Pseudomonadota | phylum | 3 | 2 | facultatively aerobic; aerobic; microaerophilic |
| Coprococcus | genus | 3 | 2 |  |
| Bifidobacteriaceae | family | 3 | 1 | facultatively aerobic; aerobic; microaerophilic |
| Enterococcus | genus | 3 | 1 | facultatively aerobic; aerobic; microaerophilic |
| Microbacteriaceae | family | 3 | 1 | facultatively aerobic; aerobic; microaerophilic |
| Alistipes | genus | 3 | 1 |  |
| Bacteroides | genus | 2 | 2 | facultatively aerobic; aerobic; microaerophilic |
| Bacteroidota | phylum | 2 | 2 | facultatively aerobic; aerobic; microaerophilic |
| Corynebacteriaceae | family | 2 | 2 | facultatively aerobic; aerobic; microaerophilic |
| Erysipelotrichaceae | family | 2 | 2 | facultatively aerobic; aerobic; microaerophilic |

Top decreased taxa:

| taxon | rank | profiles | studies | selected traits |
| --- | --- | ---: | ---: | --- |
| Bacillota | phylum | 2 | 2 | facultatively aerobic; aerobic; microaerophilic |
| Lactobacillus | genus | 2 | 2 | facultatively aerobic; aerobic; microaerophilic |
| Streptococcus | genus | 2 | 2 | facultatively aerobic; aerobic; microaerophilic |
| Lachnospiraceae | family | 2 | 2 | facultatively aerobic; microaerophilic |
| Dialister | genus | 2 | 2 | aerobic |
| Prevotella | genus | 2 | 2 |  |
| Ruminococcus | genus | 2 | 2 |  |
| Veillonella | genus | 2 | 2 |  |
| Actinomyces | genus | 1 | 1 | facultatively aerobic; aerobic; microaerophilic |
| Actinomycetota | phylum | 1 | 1 | facultatively aerobic; aerobic; microaerophilic |
| Aerococcus | genus | 1 | 1 | facultatively aerobic; aerobic; microaerophilic |
| Bacteroidales | order | 1 | 1 | facultatively aerobic; aerobic; microaerophilic |

### Autoimmune disease

Top increased taxa:

| taxon | rank | profiles | studies | selected traits |
| --- | --- | ---: | ---: | --- |
| Bacteroides | genus | 3 | 2 | facultatively aerobic; aerobic; microaerophilic |
| Bacteroides ovatus | species | 3 | 2 |  |
| Gammaproteobacteria | class | 3 | 1 | facultatively aerobic; aerobic; microaerophilic |
| Pseudomonadota | phylum | 3 | 1 | facultatively aerobic; aerobic; microaerophilic |
| Veillonellaceae | family | 3 | 1 | aerobic |
| Bacilli | class | 2 | 1 | facultatively aerobic; aerobic; microaerophilic |
| Bacillota | phylum | 2 | 1 | facultatively aerobic; aerobic; microaerophilic |
| Bacteroidaceae | family | 2 | 1 | facultatively aerobic; aerobic; microaerophilic |
| Bacteroidales | order | 2 | 1 | facultatively aerobic; aerobic; microaerophilic |
| Bacteroidia | class | 2 | 1 | facultatively aerobic; aerobic; microaerophilic |
| Bacteroidota | phylum | 2 | 1 | facultatively aerobic; aerobic; microaerophilic |
| Clostridia | class | 2 | 1 | facultatively aerobic; aerobic; microaerophilic |

Top decreased taxa:

| taxon | rank | profiles | studies | selected traits |
| --- | --- | ---: | ---: | --- |
| Bacteroides fragilis | species | 3 | 1 |  |
| Eubacteriaceae | family | 3 | 1 |  |
| Bacteroides thetaiotaomicron | species | 2 | 2 |  |
| Haemophilus parainfluenzae | species | 2 | 1 | facultatively aerobic; microaerophilic |
| Lachnospiraceae | family | 2 | 1 | facultatively aerobic; microaerophilic |
| Oscillospiraceae | family | 2 | 1 | facultatively aerobic; microaerophilic |
| Porphyromonadaceae | family | 2 | 1 | aerobic |
| Eubacterium | genus | 2 | 1 |  |
| Faecalibacterium sp. DJF_VR20 | species | 2 | 1 |  |
| Lachnospira eligens | species | 2 | 1 |  |
| Mediterraneibacter gnavus | species | 2 | 1 |  |
| human intestinal firmicute CO19 | species | 2 | 1 |  |

## Red/Blue Pattern Labels

These lists make the profile-similarity interpretation auditable. Shared red labels have positive signed counts in both exact `Autism` and `Autoimmune disease`, supporting the increased-profile similarity readout. Shared blue labels would have negative signed counts in both exact labels. The ASD blue-contrast list marks taxa where the broader `Autism spectrum disorder` label is decreased-skewed while at least one exact comparison label is neutral or increased-skewed.

### Shared red labels: exact Autism and Autoimmune disease

| heatmap label | Autism | ASD | Autoimmune disease | selected traits | in heatmap |
| --- | ---: | ---: | ---: | --- | --- |
| Bacteroides genus | 1 | 1 | 2 | aerobic; facultatively aerobic; microaerophilic | yes |
| Pseudomonadota phylum | 3 | 0 | 3 | aerobic; facultatively aerobic; microaerophilic | no |
| Enterobacteriaceae family | 1 | 0 | 2 | aerobic; facultatively aerobic; microaerophilic | no |
| Bacteroidota phylum | 1 | 0 | 1 | aerobic; facultatively aerobic; microaerophilic | no |
| Bacteroidia class | 1 | 0 | 1 | aerobic; facultatively aerobic; microaerophilic | no |
| Fusobacteriaceae family | 1 | 0 | 1 |  | no |

### Shared blue labels: exact Autism and Autoimmune disease

| heatmap label | Autism | ASD | Autoimmune disease | selected traits | in heatmap |
| --- | ---: | ---: | ---: | --- | --- |
| Lachnospiraceae family | -1 | -1 | -1 | facultatively aerobic; microaerophilic | no |

### ASD blue-contrast labels

| heatmap label | Autism | ASD | Autoimmune disease | selected traits | in heatmap |
| --- | ---: | ---: | ---: | --- | --- |
| Prevotella genus | 0 | -2 | 0 |  | yes |
| Lactobacillus genus | 0 | -1 | 1 | aerobic; facultatively aerobic; microaerophilic | yes |
| Corynebacterium genus | 1 | -1 | 0 | aerobic; facultatively aerobic; microaerophilic | yes |
| Mediterraneibacter gnavus species | 0 | -1 | -1 |  | yes |
| Ruminococcus genus | -1 | -1 | 0 |  | yes |
| Lachnospira eligens species | 0 | -1 | -2 |  | yes |
| Aerococcus genus | 0 | -1 | 0 | aerobic; facultatively aerobic; microaerophilic | yes |
| Prevotellaceae family | 0 | -1 | -1 | aerobic | no |
| Lactococcus genus | 0 | -1 | 0 | aerobic; facultatively aerobic; microaerophilic | no |
| Listeria genus | 0 | -1 | 0 | aerobic; facultatively aerobic; microaerophilic | no |
| Methylobacterium genus | 0 | -1 | 0 | aerobic; facultatively aerobic; microaerophilic | no |
| Microbacterium genus | 0 | -1 | 0 | aerobic; facultatively aerobic; microaerophilic | no |
| Pseudomonas genus | 0 | -1 | 0 | aerobic; facultatively aerobic; microaerophilic | no |
| Staphylococcus genus | 0 | -1 | 0 | aerobic; facultatively aerobic; microaerophilic | no |
| Gemella morbillorum species | 0 | -1 | 0 | facultatively aerobic; microaerophilic | no |

## Inferred Guild Evidence

This section uses a data-derived trait-enrichment view: for each exact condition label, it identifies v2 traits enriched in that condition's increased profiles, then lists taxa observed in those same profiles that also belong to the KG-Microbe-wide candidate universe and carry the enriched trait annotations. The universe contains 1,992 taxa with at least one core oxygen/respiration criterion annotation, not the fixed focal taxa list.

| condition | increased profiles | universe taxa / observed taxa | enriched traits | supporting taxa |
| --- | ---: | ---: | --- | --- |
| Autism | 8 | 24 / 36 | aerobic (+1.05); facultatively aerobic (+0.78); microaerophilic (+0.42) | Pseudomonadota (3 traits, 3 profiles); Enterococcus (3 traits, 3 profiles); Bifidobacteriaceae (3 traits, 3 profiles); Microbacteriaceae (3 traits, 3 profiles); Corynebacteriaceae (3 traits, 2 profiles); Bacteroidota (3 traits, 2 profiles); Bifidobacterium (3 traits, 2 profiles); Erysipelotrichaceae (3 traits, 1 profile); Lactobacillus (3 traits, 1 profile); Corynebacterium (3 traits, 1 profile); Actinomycetes (3 traits, 1 profile); Bacillaceae (3 traits, 1 profile) |
| Autoimmune disease | 4 | 45 / 84 | facultatively aerobic (+0.60); aerobic (+0.45); microaerophilic (+0.21) | Bacteroides (3 traits, 3 profiles); Pseudomonadota (3 traits, 3 profiles); Gammaproteobacteria (3 traits, 3 profiles); Bacillota (3 traits, 2 profiles); Streptococcaceae (3 traits, 2 profiles); Streptococcus (3 traits, 2 profiles); Bacteroidales (3 traits, 2 profiles); Clostridia (3 traits, 2 profiles); Eubacteriales (3 traits, 2 profiles); Lactobacillales (3 traits, 2 profiles); Bacteroidia (3 traits, 2 profiles); Enterobacteriaceae (3 traits, 2 profiles) |

## Interpretation

This v2 follow-up uses only the narrowed aerobe v2 selected traits: facultatively aerobic, aerobic, microaerophilic. It supports a label-dependent signal check rather than a broad causal claim.

Included condition labels with positive increased-minus-decreased v2 trait-enrichment shifts are `Autism` (0.749), `Autoimmune disease` (0.419).

At the taxa-list level, overlap should still be interpreted cautiously. The signed heatmap counts direction balance by taxon/profile occurrence, while the guild annotation columns only indicate whether those taxa carry one or more of the v2 oxygen-tolerance traits.

`autism_autoimmune_taxa_signed_pattern_labels.tsv` extracts the heatmap labels that directly support the red/blue readout above so the visual pattern can be audited without manually reading the plot.

The most defensible interpretation is that the exact autism and autoimmune disease labels both show taxa annotated as facultative or aerotolerant enriched among increased profiles. This analysis does not show that autism and autoimmune disease are causally linked, driven by the same taxa, or even consistently connected across broader ASD metadata. The autoimmune subset is small, so the taxa-level results are candidates for follow-up.

## Heatmap

`autism_autoimmune_taxa_signed_heatmap.png` shows the top 45 genus/species taxa selected from 210 observed taxa. Cell values are `increased unique taxon-profile count - decreased unique taxon-profile count`; red is increased-skewed and blue is decreased-skewed. Colored x-axis labels and the marker strip above the heatmap identify taxa included in the red/blue pattern-label extraction.

## Files

- `autism_autoimmune_profile_summary.tsv`
- `autism_autoimmune_taxa_direction_summary.tsv`
- `autism_autoimmune_taxa_direction_summary_grouped.tsv`
- `autism_autoimmune_taxa_signed_balance.tsv`
- `autism_autoimmune_taxa_signed_heatmap_matrix.tsv`
- `autism_autoimmune_taxa_signed_pattern_labels.tsv`
- `autism_autoimmune_taxa_signed_heatmap.png`
- `autism_autoimmune_inferred_guild_summary.tsv`
- `autism_autoimmune_inferred_guild_trait_evidence.tsv`
- `autism_autoimmune_inferred_guild_taxa_evidence.tsv`
- `kg_microbe_candidate_taxon_universe.tsv`
