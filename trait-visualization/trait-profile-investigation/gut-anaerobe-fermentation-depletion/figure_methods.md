# Figure Methods and Interpretations

These notes describe every PNG in this investigation folder. Interpretations are descriptive summaries of the labels and patterns visible in the generated figures; they do not re-estimate differential abundance beyond the BugSigDB profile directions used by the scripts.

## Shared Inputs

Study-profile z-scores were read from `trait-visualization/study_trait_profile_clustered_heatmap_matrix.tsv` and joined to study metadata from `trait-visualization/study_trait_profile_study_metadata.tsv`. The gut anaerobe/fermentation trait-enrichment score was the row-wise mean of anaerobic, obligately anaerobic, Fermentation, cellulose catabolic process, and xylan catabolic process traits. Gut-context profiles were defined by body-site or grouped body-site labels matching `feces`, `faeces`, `stool`, `colon`, `caecum`, `cecum`, `gut`, `intestin`, or `digestive`.

Focal taxa were defined a priori by NCBI Taxonomy ID: Faecalibacterium prausnitzii, Anaerobutyricum hallii, Lachnospira eligens, Agathobacter rectalis, Roseburia intestinalis, Ruminococcus bromii, and Roseburia inulinivorans.

## `candidate_taxa_depletion_bias.png`

Methods: The script intersected the focal taxon IDs with `trait-visualization/taxa_guild_discovery_taxa_metadata.tsv`, then calculated `(decreased profile count - increased profile count) / total profile count` for each focal taxon. No condition or body-site filter was applied, so this is an overall direction-count view rather than a condition-specific guild-score plot.

Interpretation: All focal gut anaerobe candidates are depletion-biased. Agathobacter rectalis, Lachnospira eligens, Roseburia intestinalis, and Roseburia inulinivorans have the largest depletion bias, while Faecalibacterium prausnitzii and Anaerobutyricum hallii still lean decreased despite larger mixed profile counts. The figure's theme is a broad depletion tendency among gut anaerobe and fiber-fermentation taxa.

## `condition_decreased_guild_signal.png`

Methods: Gut-context study profiles were grouped by BugSigDB condition label. Groups with fewer than five profiles were excluded, `Clinical treatment` was removed from the displayed summary, and top conditions were selected by positive decreased-minus-increased guild z-score and positive decreased-profile mean guild z-score. Bar length is the condition-level contrast, with paired standardized mean difference and paired t-test labels when matched study/experiment profiles were available.

Interpretation: The leading labels include colorectal carcinoma, graft-versus-host disease, dysmenorrhea, acute graft-versus-host disease, spinal cord injury, diarrhea, prediabetes syndrome, and response to allogeneic hematopoietic stem cell transplant. These themes cluster around gut malignancy, transplant or mucosal injury, diarrheal/inflammatory states, and metabolic labels. The highest bars identify contexts where decreased profiles are especially enriched for anaerobe/fermentation traits.

## `guild_trait_direction_shift.png`

Methods: For each selected guild trait, the script calculated the mean z-score in all decreased profiles and the mean z-score in all increased profiles. Bars show `decreased mean z - increased mean z`; positive bars indicate traits enriched in decreased profiles.

Interpretation: All five selected traits are enriched in decreased profiles. Cellulose and xylan catabolic processes have the largest shifts, followed by anaerobic, obligately anaerobic, and Fermentation. The theme is not only oxygen intolerance; the original trait set also emphasizes fiber-degradation and fermentation capabilities that are characteristic of gut commensal anaerobes.

## `anatomical_sites_diagram.png`

Methods: Gut-context profiles were grouped by `body_site_group`, requiring at least 10 profiles per group. The generic `Other` group was excluded from the diagram. Bubble and bar length encode profile count, and color encodes the decreased-profile gut anaerobe/fermentation guild z-score at each body-site group.

Interpretation: Colon, feces, and caecum are the main anatomical labels, which is consistent with the gut-context filter and the biology of anaerobe/fiber-fermentation traits. Colon shows the strongest trait enrichment in decreased profiles, with feces providing the largest profile-count context. The figure supports interpreting this trait set primarily as a lower-gut depletion signal.

## `target_taxa_condition_breakdown_bar.png`

Methods: Focal taxon occurrences were loaded from `data/derived/taxa.tsv`, deduplicated by taxon and study profile, joined to scored study profiles, and counted by condition for decreased-profile occurrences in gut-context profiles. A condition can contribute multiple counts if multiple focal taxa appear in the same profile.

Interpretation: Parkinson's disease, COVID-19, and colorectal cancer dominate the recurrent focal-taxon depletion labels, with additional gut or inflammatory labels contributing smaller counts. The target-taxon view is consistent with the trait-enrichment plots but answers a different question: where the named anaerobe candidates are repeatedly decreased, not where the composite trait-enrichment score is highest.
