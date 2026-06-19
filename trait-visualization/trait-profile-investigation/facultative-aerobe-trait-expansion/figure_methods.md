# Figure Methods and Interpretations

These notes describe every PNG in this investigation folder. Interpretations are descriptive summaries of the labels and patterns visible in the generated figures; they do not re-estimate differential abundance beyond the BugSigDB profile directions used by the scripts.

## Shared Inputs

Study-profile z-scores were read from `trait-visualization/study_trait_profile_clustered_heatmap_matrix.tsv` and joined to study metadata from `trait-visualization/study_trait_profile_study_metadata.tsv`. The facultative-aerobe guild score was the row-wise mean of facultatively aerobic, aerobic, gram negative, motile, microaerophilic, and nitrate_reduction traits. Direction labels (`increased` and `decreased`) were taken from BugSigDB's curated `Abundance in Group 1` signature field.

Focal taxa were defined a priori by NCBI Taxonomy ID: Parvimonas micra, Klebsiella pneumoniae, Streptococcus anginosus, Escherichia coli, Fusobacterium nucleatum, Enterocloster bolteae, Eggerthella lenta, Prevotella intermedia, Porphyromonas gingivalis, and Veillonella dispar.

## `candidate_taxa_increased_profile_tendency.png`

Methods: The script intersected the focal taxon IDs with `trait-visualization/taxa_guild_discovery_taxa_metadata.tsv`, then calculated `(increased profile count - decreased profile count) / total profile count` for each focal taxon. No condition or body-site filter was applied, so this is an overall BugSigDB direction-count view rather than a guild-score plot.

Interpretation: Every focal taxon has a positive increased-profile tendency. The strongest labels include Parvimonas micra, Klebsiella pneumoniae, Streptococcus anginosus, Escherichia coli, and Fusobacterium nucleatum, a mixed set of gut, oral, and opportunistic/pathobiont taxa. The theme is not a single disease label; it is a cross-context tendency for aerotolerant or facultative taxa to appear in increased signatures.

## `condition_increased_guild_signal.png`

Methods: Study profiles from all body sites were grouped by BugSigDB condition label. Groups with fewer than five profiles were excluded, broad metadata labels were removed, and the top eligible conditions were selected by positive increased-minus-decreased guild z-score and positive increased-profile mean guild z-score. Bar length is the condition-level increased-minus-decreased guild contrast, with paired standardized mean difference and paired t-test labels when matched study/experiment profiles were available.

Interpretation: The strongest condition labels include lung transplantation, female infertility, multiple myeloma, tuberculosis, autism, graft-versus-host disease, colorectal carcinoma, and vulvovaginitis. These labels cluster around immune disruption, transplantation, oncology, infection, and reproductive or mucosal contexts. Several labels can mix body sites, so the figure is best read as a condition-screening view rather than evidence that all signals are gut-specific.

## `guild_trait_direction_shift.png`

Methods: For each selected guild trait, the script calculated the mean z-score in all increased profiles and the mean z-score in all decreased profiles. Bars show `increased mean z - decreased mean z`; positive bars indicate traits enriched in increased profiles.

Interpretation: All six selected traits shift upward in increased profiles. The largest shifts are for facultatively aerobic and aerobic traits, with gram negative, microaerophilic, motile, and nitrate_reduction providing supporting but smaller signals. The figure shows that the broad guild behaves coherently, though it combines oxygen tolerance with envelope, motility, and respiration-associated features.

## `anatomical_sites_diagram.png`

Methods: Study profiles were grouped by `body_site_group`, requiring at least 10 profiles per group. The generic `Other` group was excluded from the diagram. Bubble and bar length encode profile count, and color encodes the increased-profile facultative-aerobe guild z-score at each body-site group.

Interpretation: Oxygen-exposed sites dominate the visual theme: nasopharynx, saliva, subgingival dental plaque, and mouth carry prominent increased-profile signal. Feces and colon are still shown because they are high-count gut contexts, but their signal is weaker than the oral or airway-associated labels. This supports splitting oral/airway patterns from gut interpretations.

## `target_taxa_condition_breakdown_bar.png`

Methods: Focal taxon occurrences were loaded from `data/derived/taxa.tsv`, deduplicated by taxon and study profile, joined to scored study profiles, and counted by condition for increased-profile occurrences. A condition can contribute multiple counts if multiple focal taxa appear in the same profile.

Interpretation: Colorectal cancer is the dominant condition theme for focal increased taxa, especially for Fusobacterium nucleatum, Parvimonas micra, Prevotella intermedia, and Escherichia coli. Other recurrent labels include COVID-19, type II diabetes mellitus, Parkinson's disease, smoking behavior, rheumatoid arthritis, and oral plaque contexts. The figure therefore mixes gut cancer/metabolic/neuro labels with oral inflammatory labels.

## `autism_autoimmune_taxa_signed_heatmap.png`

Methods: The follow-up script selected Autism, Autism spectrum disorder, and Autoimmune disease profiles, deduplicated taxon/profile occurrences, and built a signed condition-by-taxon matrix where each cell is increased-profile count minus decreased-profile count. Rows are exact condition labels and columns are the top genus/species taxa selected by signed count and profile support.

Interpretation: The Autism spectrum disorder row carries the densest signal and shows many negative cells, indicating more decreased than increased profile counts for many selected taxa. Autism and Autoimmune disease have fewer profiles and sparser cells, with isolated positive and negative pockets. The main theme is label-specific heterogeneity: the exact autism spectrum label differs visibly from the broader Autism and Autoimmune disease labels, so these rows should not be collapsed without checking direction balance.
