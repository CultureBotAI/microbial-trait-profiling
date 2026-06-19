# Figure Methods and Interpretations

These notes describe every PNG in this investigation folder. Interpretations are descriptive summaries of the labels and patterns visible in the generated figures; they do not re-estimate differential abundance beyond the BugSigDB profile directions used by the scripts.

## Shared Inputs

Study-profile z-scores were read from `trait-visualization/study_trait_profile_clustered_heatmap_matrix.tsv` and joined to study metadata from `trait-visualization/study_trait_profile_study_metadata.tsv`. The facultative-aerobe v2 guild score was the row-wise mean of facultatively aerobic, aerobic, and microaerophilic traits. Direction labels (`increased` and `decreased`) were taken from BugSigDB's curated `Abundance in Group 1` signature field.

Focal taxa were defined a priori by NCBI Taxonomy ID: Parvimonas micra, Klebsiella pneumoniae, Streptococcus anginosus, Escherichia coli, Fusobacterium nucleatum, Enterocloster bolteae, Eggerthella lenta, Prevotella intermedia, Porphyromonas gingivalis, and Veillonella dispar.

## `candidate_taxa_increased_profile_tendency.png`

Methods: The script intersected the focal taxon IDs with `trait-visualization/taxa_guild_discovery_taxa_metadata.tsv`, then calculated `(increased profile count - decreased profile count) / total profile count` for each focal taxon. No condition or body-site filter was applied, so this is an overall direction-count view rather than a v2 guild-score plot.

Interpretation: All focal taxa lean toward increased profiles. The top labels are Parvimonas micra, Klebsiella pneumoniae, Streptococcus anginosus, Escherichia coli, and Fusobacterium nucleatum, followed by other oral/gut aerotolerant candidates. The theme is a repeated increased-profile tendency among potential pathobionts and oxygen-tolerant taxa, not a single-condition effect.

## `condition_increased_guild_signal.png`

Methods: Study profiles from all body sites were grouped by BugSigDB condition label. Groups with fewer than five profiles were excluded, broad metadata labels were removed, and top conditions were selected by positive increased-minus-decreased v2 guild z-score and positive increased-profile mean guild z-score. Bar length is the condition-level contrast, with paired standardized mean difference and paired t-test labels when matched study/experiment profiles were available.

Interpretation: The strongest labels include lung transplantation, age-related macular degeneration, multiple myeloma, tuberculosis, sarcopenia, graft-versus-host disease, acute graft-versus-host disease, and autism. These themes point toward immune disturbance, transplant/oncology contexts, infection, aging-related labels, and neurodevelopmental labels. Because the figure uses all body sites, the labels should be followed by body-site checks before making gut-specific claims.

## `guild_trait_direction_shift.png`

Methods: For each v2 guild trait, the script calculated mean z-score in all increased profiles and mean z-score in all decreased profiles. Bars show `increased mean z - decreased mean z`; positive bars indicate traits enriched in increased profiles.

Interpretation: Facultatively aerobic, aerobic, and microaerophilic traits all shift upward in increased profiles. Facultatively aerobic is the largest component, with aerobic and microaerophilic providing a consistent oxygen-tolerance pattern. The v2 restriction removes gram-negative, motility, and nitrate-reduction terms, so this figure is a cleaner oxygen-tolerance readout than the original broader guild.

## `anatomical_sites_diagram.png`

Methods: Study profiles were grouped by `body_site_group`, requiring at least 10 profiles per group. The generic `Other` group was excluded from the diagram. Bubble and bar length encode profile count, and color encodes the increased-profile facultative-aerobe v2 guild z-score at each body-site group.

Interpretation: Nasopharynx, saliva, mouth, and subgingival dental plaque are visually prominent, matching oxygen-exposed oral and airway sites. Feces and colon have large profile counts but weaker increased-profile v2 signal than the oral/airway labels. This figure reinforces that oral and airway ecology can drive the v2 guild signal and should be separated from gut-only interpretations.

## `target_taxa_condition_breakdown_bar.png`

Methods: Focal taxon occurrences were loaded from `data/derived/taxa.tsv`, deduplicated by taxon and study profile, joined to scored study profiles, and counted by condition for increased-profile occurrences. A condition can contribute multiple counts if multiple focal taxa appear in the same profile.

Interpretation: Colorectal cancer dominates the target-taxon occurrence labels, especially through Fusobacterium nucleatum, Parvimonas micra, Prevotella intermedia, and Escherichia coli. Other recurrent themes include COVID-19, type II diabetes mellitus, Parkinson's disease, smoking behavior, rheumatoid arthritis, and oral plaque or saliva conditions. The figure shows where the fixed focal taxa recur, not where the v2 guild score itself is strongest.

## `autism_autoimmune_taxa_signed_heatmap.png`

Methods: The follow-up script selected Autism, Autism spectrum disorder, and Autoimmune disease profiles, deduplicated taxon/profile occurrences, and built a signed condition-by-taxon matrix where each cell is increased-profile count minus decreased-profile count. Rows are exact condition labels and columns are the top genus/species taxa selected by signed count and profile support.

Interpretation: Autism spectrum disorder is the densest row and has many negative cells, showing that many selected taxa are more often decreased than increased for that exact label. Autism and Autoimmune disease have fewer profiles and sparser signals. The figure's main theme is label heterogeneity: related condition labels can have different direction patterns, so interpretation should preserve exact labels when possible.
