# PD-cellulose draft for *Microbiome trait profiles across health and disease*

## Purpose and interpretation boundary

This document proposes where the Parkinson's disease (PD)-cellulose analysis should enter the current manuscript and supplies manuscript-ready Results, Discussion, Methods, and Limitations text. It also provides a separate interpretation and takeaway for every diagram currently in this directory.

The central claim should remain deliberately narrow:

> PD fecal signatures show a reproducible descriptive opposition between cellulose catabolic process, which is stronger in decreased profiles, and a three-trait group comprising moderately halophilic, facultatively aerobic, and arsenate reduction, which is stronger in increased profiles. The pattern is direction-aware and distributed across taxa. It is not a validated diagnostic signature, a measurement of cellulose degradation, or evidence that a single organism causes the contrast.

The strongest globally adjusted result is PD-moderately-halophilic, not PD-cellulose. In the global screen of 336 eligible condition-trait tests, moderately halophilic is the only FDR-significant PD trait (`q=0.048`). The PD-cellulose contrast is directionally coherent and has a nominal paired-study `p=0.0124`, but it does not remain significant after the global correction (`q=0.597`). That distinction should be visible wherever the cellulose result is discussed.

Profile counts in this analysis are curated BugSigDB differential-abundance signature occurrences. They are not patients, sample prevalence, gene expression, metabolic flux, or direct cellulose-degradation measurements.

## Recommended location in the paper

### Primary Results location

Insert a new Results subsection in **Combined interpretation by condition and body-site**, immediately after the autism-by-autoimmune case study and its Figure 6, and immediately before **Discussion**. In the current v1 PDF, this is at the end of page 17.

This is the best fit because the PD analysis:

- re-projects the six cross-framing convergent traits onto one condition and one body site;
- depends on direction-aware increased-minus-decreased contrasts from the composition framing;
- tests whether the resulting condition-level pattern is specific when compared with other fecal conditions; and
- uses taxon-to-trait links as annotation provenance rather than reverting to a taxon-first disease signature.

Suggested subsection title:

> **Parkinson's disease reveals a direction-aware ecological-reporter pattern**

The preferred main figure is [`pd_ecological_reporter_combined.png`](pd_ecological_reporter_combined.png). It can become new **Figure 7**. The current Methods figures numbered 7A and 7B would then become Figures 8A and 8B. If renumbering is undesirable, retain this figure as a supplementary figure and keep the Results subsection in the main text.

Use [`pd_condition_taxon_trait_network.png`](pd_condition_taxon_trait_network.png) as the inferential companion, preferably in the supplement. It supplies the global FDR context that prevents the descriptive ecological-reporter figure from being read as a biomarker claim.

### Secondary Discussion location

Add the proposed Discussion paragraph under **Biological synthesis** on current PDF page 18. The cleanest location is after the sentence ending with the autism-by-autoimmune case study and before the recommendation that trait-level meta-analyses run both framings.

### Methods location

Add a short **PD fecal case-study analysis** subsection after **Integration procedure** on current PDF page 23 and before Table 2. This keeps the condition-selection rule, equal-study weighting, composite definition, and global correction family auditable.

### Limitations location

Add the proposed limitations paragraph to **Limitations** on current PDF page 19, after the paragraph about KG-Microbe and METPO annotation completeness.

### Optional abstract update

If the PD analysis is promoted to a main Results figure, the current abstract sentence ending in "the autism by autoimmune contrast" could be expanded to:

> The combined view reshapes condition-by-body-site trait modules and clarifies case studies including the autism-by-autoimmune contrast and a direction-aware PD fecal pattern in which cellulose catabolic process is stronger in decreased signatures while moderately halophilic, facultatively aerobic, and arsenate reduction are stronger in increased signatures.

Do not call the pattern PD-specific, diagnostic, or FDR-confirmed as a whole in the abstract.

## Evidence hierarchy for the PD interpretation

| Evidence level | Result | Interpretation permitted |
| --- | --- | --- |
| Global FDR-supported | PD-moderately-halophilic: mean paired-study contrast `+0.530`, raw `p=0.000143`, global `q=0.048`, `n=45` paired studies | Moderately halophilic is stronger in PD increased profiles in the globally corrected six-trait condition screen. |
| Exploratory after global correction | PD-arsenate reduction: `+0.506`, `q=0.164`; PD-*Roseburia* genus direction: `q=0.192` | Follow-up hypotheses only; neither is an FDR-significant discovery. |
| Directionally coherent but not globally significant | PD-facultatively-aerobic: `+0.463`, `q=0.597`; PD-cellulose: `-0.651`, `q=0.597` | These traits help describe the PD pattern, but should not be called globally significant. |
| Descriptive study-level composite | Three-trait increased-side composite `+0.499` (95% descriptive interval `0.221` to `0.778`); ecological-reporter composite `+1.150` (`0.535` to `1.765`) | Pattern summary only. The composite is data-derived and is not a validated score or an independently corrected test. |
| Nominal PD-versus-other-feces comparison | Decreased cellulose-positive rate `54.5%` in PD versus `44.8%` in other feces; odds ratio `1.48`, Fisher `p=0.0391`; z-score Mann-Whitney `p=0.0232` | PD decreased profiles contain more cellulose signal than the same direction in other fecal conditions, before study clustering and global multiplicity control. |
| Study-clustered interaction | Cellulose z-score PD-by-decreased coefficient `+0.400`, `p=0.086`; positive-z interaction `+0.189`, `p=0.053` | Suggestive extra PD shift in the decreased direction, but no conventional `p<0.05` interaction. |
| Annotation provenance | Recurrent cellulose carriers include Lachnospiraceae, *Roseburia*, *Blautia*, and Clostridiaceae; no condition-taxon test passes global `q<=0.05` | Taxa help trace where annotations enter the profile score. They are not inferential drivers or causal organisms. |

## Manuscript-ready Results subsection

### Parkinson's disease reveals a direction-aware ecological-reporter pattern

We used the six cross-framing convergent traits as fixed axes for a fecal Parkinson's disease case study. At the profile level, the PD subset contained 130 increased and 123 decreased profiles, contributed by 50 and 46 studies, respectively. Cellulose catabolic process was present at positive z-score in 35 increased profiles but 67 decreased profiles, and its mean z-score was `-0.015` in increased profiles compared with `0.594` in decreased profiles. The resulting profile-weighted increased-minus-decreased contrast was `-0.609`, making cellulose the clearest decreased-profile-skewed member of the six-trait set.

The paired, equal-study-weighted analysis gave the same directional conclusion without allowing studies containing more profiles to dominate. Across 45 studies represented in both directions, the mean cellulose contrast was `-0.651` (descriptive 95% t interval `-1.124` to `-0.178`). Moderately halophilic (`+0.530`), facultatively aerobic (`+0.463`), and arsenate reduction (`+0.506`) pointed in the opposite direction, whereas coccus shaped (`+0.137`) and sphere shaped (`+0.133`) were smaller and their descriptive intervals crossed zero. We summarized the opposing pattern as an exploratory ecological-reporter composite: the mean contrast of moderately halophilic, facultatively aerobic, and arsenate reduction minus the cellulose contrast. The PD composite was `+1.150` (descriptive 95% interval `0.535` to `1.765`). Because its component traits were selected from the PD pattern, this composite is a visualization aid rather than an independently validated score.

The cross-condition comparison did not support uniqueness to PD. Among the 13 fecal condition labels retained because at least one of the six individual traits had a descriptive study-level interval excluding zero, PD occupied the quadrant combining decreased-profile cellulose skew with an increased-side three-trait composite, but several other conditions occupied the same region and some had larger composite values. The global evidence screen was still more restrictive. Across 336 eligible condition-trait tests, only PD-moderately-halophilic passed Benjamini-Hochberg correction (`q=0.048`). PD-arsenate reduction was exploratory (`q=0.164`), while PD-cellulose and PD-facultatively-aerobic each had `q=0.597`. Thus, the multi-trait display describes a coherent PD ecological-reporter pattern, but the pattern as a whole is neither PD-specific nor globally FDR-confirmed.

Within direction, PD decreased fecal profiles had a higher cellulose-positive rate than other decreased fecal profiles (`54.5%` versus `44.8%`; odds ratio `1.48`, Fisher `p=0.0391`), and their raw z-score distributions also differed nominally (Mann-Whitney `p=0.0232`). The converse contrast appeared among increased profiles, where the cellulose-positive rate was `26.9%` in PD and `36.1%` in other fecal profiles. A feces-only regression with study-clustered standard errors estimated an additional PD-associated shift in the decreased direction, but the interaction did not cross the conventional threshold for either the continuous cellulose z-score (`+0.400`, `p=0.086`) or the positive-z indicator (`+0.189`, `p=0.053`). These analyses support concentration of the PD cellulose signal in decreased profiles, but not a definitive PD-by-direction interaction.

The taxonomic provenance was distributed rather than dominated by one resolved species. Frequently recurring cellulose-annotated labels in PD decreased profiles included Lachnospiraceae (39 decreased versus 10 increased profiles), *Roseburia* (31 versus 4), *Blautia* (15 versus 2), and Clostridiaceae (13 versus 1). Other lineages showed different or balanced patterns, including Pseudomonadota (1 decreased versus 15 increased) and Oscillospiraceae (13 versus 15), indicating that cellulose annotation alone does not impose one direction. No condition-taxon comparison passed global FDR correction; the strongest taxon result, the *Roseburia* genus, remained exploratory (`q=0.192`). Species-removal influence diagnostics likewise did not assign any of the 30 retained species-trait cells a supported direction because every 95% interval crossed zero. Together, these results favor a community-level, annotation-mediated interpretation over attribution to a single microbial species.

**Results takeaway.** PD fecal signatures contain a direction-aware four-trait opposition: cellulose catabolic process is stronger in decreased profiles, while moderately halophilic, facultatively aerobic, and arsenate reduction are stronger in increased profiles. The signal is useful as an ecological-reporter hypothesis and as an example of information recovered by paired trait analysis. Only the moderately-halophilic component is globally FDR-supported in the current cross-condition screen, so cellulose should remain a descriptive component of the pattern rather than a standalone PD biomarker.

## Suggested main-figure caption

**Figure 7. Parkinson's disease ecological-reporter hypothesis, cross-condition specificity, and taxon-to-trait provenance** (`clustered-heatmap-generation` at `205257c`). **A**, Equal-study-weighted increased-minus-decreased contrasts for the six convergent traits across 13 human fecal condition labels having a descriptive 95% interval excluding zero for at least one individual trait. The ecological-reporter composite is the mean contrast of moderately halophilic, facultatively aerobic, and arsenate reduction minus the cellulose catabolic process contrast; the composite does not qualify a condition for inclusion. Dark borders mark descriptive individual-trait intervals excluding zero. **B**, The same conditions positioned by reversed cellulose contrast on the x-axis and the three-trait increased-side composite on the y-axis. PD lies in the upper-right quadrant but is not unique within it. Point area represents paired-study count. **C**, Exact genus and species labels reported in PD fecal BugSigDB signatures crossed with local KG-Microbe annotations for the six traits. Circle size represents distinct-study support and color represents the log2 ratio of increased to decreased study support; an `x` indicates no local annotation. Panel C is annotation provenance, not a significance or activity map. All contrasts summarize curated differential-abundance signatures rather than patient-level measurements.

## Draft-paper figure description: condition-taxon-trait evidence landscape

**Figure X. Globally adjusted human fecal condition-taxon and condition-trait evidence landscape.** The analysis includes 56 conditions with at least three studies containing both increased and decreased BugSigDB signatures. Benjamini-Hochberg correction was applied separately across the complete condition-taxon and condition-trait test families. Filled points are FDR-significant (`q<=0.05`), hollow amber points are exploratory (`0.05<q<=0.20`) and are not claimed as discoveries, and grey points show the remaining tested background. **A**, Results for all 264 eligible condition-taxon comparisons. The x-axis is the directional study fraction, with negative values indicating more decreased-only than increased-only studies; the y-axis is `-log10(global q)`. Circles denote genera, triangles denote species, and point area represents the number of reporting studies. No condition-taxon comparison passes the FDR threshold. The strongest result, decreased reporting of the *Roseburia* genus in PD, remains exploratory (directional fraction `-0.333`, `q=0.192`). **B**, Results for all 336 eligible condition-trait comparisons. The x-axis is the mean paired study-level increased-minus-decreased z-score contrast, the y-axis is `-log10(global q)`, and point area represents the number of paired studies. PD-moderately-halophilic is the only FDR-significant result (mean contrast `+0.530`, raw `p=0.000143`, `q=0.048`, `n=45`). PD-arsenate reduction (`+0.506`, `q=0.164`), COVID-19-cellulose (`-1.307`, `q=0.164`), and colorectal-cancer-cellulose (`-0.821`, `q=0.164`) are exploratory. **C**, All six PD trait contrasts, with unadjusted `p`, globally adjusted `q`, and paired-study count printed beside each estimate. Cellulose catabolic process is stronger in decreased PD profiles (contrast `-0.651`, raw `p=0.0124`) but is not significant after global correction (`q=0.597`); arsenate reduction is exploratory, and moderately halophilic is the only FDR-supported PD trait. The figure therefore supports a direction-aware multi-trait PD pattern while limiting the formal discovery claim to moderately halophilic; taxon points describe directional reporting in curated signatures and do not establish abundance effects, microbial activity, or causality.

## Manuscript-ready Discussion insertion

The PD case study shows why condition-level trait interpretation should preserve direction and distinguish pattern description from multiplicity-controlled discovery. PD fecal signatures combined a decreased-profile cellulose component with increased-profile moderately-halophilic, facultatively-aerobic, and arsenate-reduction components. This opposition was stable across profile-weighted and equal-study-weighted summaries, but it did not define a unique PD region among fecal conditions and only the moderately-halophilic component survived the global condition-trait FDR correction. The appropriate biological interpretation is therefore an ecological-reporter hypothesis: taxon composition and its annotated trait potential may jointly reflect a changed gut environment associated with PD, including possible differences in diet, intestinal transit, medication, inflammation, or oxygen and redox conditions, without identifying which of those pressures generated the pattern.

The distributed taxonomic provenance reinforces this interpretation. Cellulose signal was associated with several recurrent, partly overlapping labels rather than one supported species-level influence, and no condition-taxon result survived global correction. The trait contrast may therefore be more portable than any individual taxon association, but it should not be equated with measured cellulose degradation or metabolic activity. More generally, the PD example demonstrates the value of the paper's combined reporting convention: a direction-aware composition analysis can expose a biologically interpretable multi-trait configuration, while the global inferential layer states which individual components are robust enough to claim as discoveries.

## Manuscript-ready Methods insertion

### PD fecal case-study analysis

The PD case study was restricted to human fecal BugSigDB profiles and to the six traits in the cross-framing convergent core. A profile-level descriptive summary compared mean composition z-scores and positive-z rates between increased and decreased PD signatures. For study-weighted analyses, increased and decreased profile z-scores were first averaged within each study and trait, and the paired contrast was defined as mean increased z minus mean decreased z. Each study therefore contributed one contrast irrespective of its number of experiments or profiles. Descriptive 95% Student t intervals were calculated across study contrasts.

The exploratory three-trait increased-side composite was the mean of the moderately-halophilic, facultatively-aerobic, and arsenate-reduction contrasts. The ecological-reporter composite was that three-trait mean minus the cellulose catabolic process contrast. These components were selected because they formed the interval-excluding-zero positive and negative sides of the PD pattern; neither composite was used to qualify conditions for the cross-condition display or treated as a validated diagnostic score. Conditions entered that display only when they had at least three paired studies and at least one of the six individual traits had a descriptive 95% study-level interval excluding zero.

For global condition-trait inference, each eligible human fecal condition required at least three studies containing both increased and decreased signatures. Each condition-trait distribution of paired study contrasts was tested against zero with a two-sided Wilcoxon signed-rank test. Benjamini-Hochberg correction was applied across the full family of 336 eligible condition-trait tests. Exact genus and species identifiers were kept separate in the condition-taxon analysis; taxa required at least three reporting studies and were tested with an exact McNemar/binomial test on discordant paired studies, followed by a separate global Benjamini-Hochberg correction across 264 eligible condition-taxon tests.

The PD-versus-other-feces sensitivity analysis compared positive-z rates with two-sided Fisher exact tests and raw z-score distributions with two-sided Mann-Whitney tests within direction. A feces-only ordinary least-squares model included indicators for PD, decreased direction, and their interaction, with standard errors clustered by study. The interaction tested whether decreased profiles had an additional PD-associated cellulose shift beyond the main PD and direction effects. Taxon-carrier, removal-influence, contributor, and co-carriage analyses were treated as annotation-provenance or sensitivity analyses rather than independent inferential tests.

## Manuscript-ready Limitations insertion

The PD case study adds further interpretation limits. Its units are curated differential-abundance signatures rather than patients, and repeated profiles within studies are not independent unless the analysis explicitly collapses or clusters them by study. KG-Microbe annotations indicate reported trait potential, not gene expression, substrate use, or metabolite flux; consequently, a cellulose catabolic process score is not a direct measurement of cellulose degradation. The ecological-reporter composite was constructed from the observed PD pattern and compared descriptively across conditions, so it requires external, prospectively specified validation before any diagnostic or specificity claim. Medication, diet, intestinal transit, geography, sequencing method, and other cohort differences can also generate or modify the observed pattern.

## Diagram-by-diagram interpretations and takeaways

### 1. [`pd_cellulose_direction_specificity.png`](pd_cellulose_direction_specificity.png)

**What it shows.** The left panel compares the two PD directions. There are 35 cellulose-positive increased PD fecal profiles from 22 studies and 67 cellulose-positive decreased profiles from 32 studies. The cellulose focus-minus-opposite-direction z contrast is larger for the decreased edge (`1.58`) than for the increased edge (`0.75`). The right panel compares PD with all other fecal profiles within direction: `26.9%` versus `36.1%` are cellulose-positive among increased profiles, while `54.5%` versus `44.8%` are positive among decreased profiles.

**Interpretation.** The cellulose signal is asymmetric. It is concentrated in PD decreased profiles rather than elevated across both PD directions. The decreased-direction PD-versus-other-feces comparison is nominally significant (odds ratio `1.48`, Fisher `p=0.0391`), but this panel does not account for repeated profiles within studies or global multiple testing.

**Takeaway.** Use this as the clearest descriptive overview of the PD-cellulose observation. Phrase the result as "decreased-profile concentration," not "PD-specific enrichment" or "biomarker." Recommended role: supplementary descriptive figure.

### 2. [`pd_cellulose_cytoscape_highlight.png`](pd_cellulose_cytoscape_highlight.png)

**What it shows.** The orange PD-to-cellulose-to-feces path is highlighted against faded links for the other five shared traits. The decreased view has the thicker cellulose path and the larger edge contrast (`n=67`, `s=32`, `delta z=1.58`) compared with the increased view (`n=35`, `s=22`, `delta z=0.75`).

**Interpretation.** This is a network rendering of the same direction-specific counts and edge scores shown elsewhere. It contributes visual orientation, not an independent statistical result.

**Takeaway.** The network makes the path from condition to trait to body site intuitive, but the thicker decreased edge must not be interpreted as causal or statistically significant by itself. Recommended role: graphical overview or supplement; omit if space is limited because it is redundant with the direction-specific bar figure.

### 3. [`pd_cellulose_enrichment_test.png`](pd_cellulose_enrichment_test.png)

**What it shows.** The feces-only group means are `0.09` for other increased, `-0.01` for PD increased, `0.30` for other decreased, and `0.59` for PD decreased profiles. In the study-clustered model, the PD effect among increased profiles is `-0.10` (`p=0.475`), the decreased effect among other fecal profiles is `+0.21` (`p=0.00011`), and the extra PD effect in decreased profiles is `+0.40` (`p=0.086`).

**Interpretation.** Some cellulose elevation in decreased profiles is a broader fecal-signature phenomenon, not uniquely a PD phenomenon. PD appears to add to that decreased-direction effect, but the PD-by-decreased interaction interval includes zero. The binary positive-z model tells the same story more narrowly (`+0.189`, `p=0.053`).

**Takeaway.** This is the key guardrail against overclaiming specificity. The PD direction interaction is suggestive, not established. Recommended role: statistical supplementary figure accompanying the descriptive direction panel.

### 4. [`pd_cellulose_taxon_driver_bar.png`](pd_cellulose_taxon_driver_bar.png)

**What it shows.** The diagram counts unique PD fecal profiles containing reported taxa that carry a local cellulose annotation. Lachnospiraceae (`39` decreased, `10` increased), *Roseburia* (`31`, `4`), *Blautia* (`15`, `2`), Clostridiaceae (`13`, `1`), *Clostridium* (`10`, `3`), *Ruminococcus* (`9`, `2`), and *Lachnospira* (`9`, `2`) favor decreased profiles. Pseudomonadota (`1`, `15`) and Actinomycetota (`0`, `8`) favor increased profiles, while Oscillospiraceae (`13`, `15`) is nearly balanced.

**Interpretation.** The aggregate cellulose contrast is assembled from multiple, directionally heterogeneous taxonomic labels. The strong decreased recurrence of Lachnospiraceae and *Roseburia* is consistent with the bulk cellulose pattern, but the opposing lineages show that a cellulose annotation does not dictate one disease direction. Exact labels span species, genus, family, order, class, and phylum and therefore cannot be added as disjoint contributions.

**Takeaway.** Refer to these rows as recurrent reported cellulose carriers, not causal "drivers." *Roseburia faecis* (`6` decreased, `0` increased profiles) is a follow-up candidate, not a validated species biomarker. Recommended role: supplementary provenance figure.

### 5. [`pd_cellulose_taxon_six_trait_heatmap.png`](pd_cellulose_taxon_six_trait_heatmap.png)

**What it shows.** The rows are the recurrent cellulose-carrier labels from the preceding bar plot, and the columns are all six convergent traits. Orange cells denote a local KG-Microbe annotation; printed values are annotation counts. The heatmap shows that some cellulose carriers are annotated only narrowly, such as *Roseburia* and *Roseburia faecis*, whereas broad labels such as Bacillota, Pseudomonadota, Clostridia, Bacilli, and Eubacteriales carry several of the six annotations.

**Interpretation.** The PD-cellulose observation is embedded in overlapping trait potential. Some reported taxa can influence multiple trait scores, so the six PD contrasts should not be interpreted as six taxonomically independent signals. Very large annotation counts for broad ranks reflect graph coverage and descendant propagation rather than abundance or effect strength.

**Takeaway.** This figure is an annotation-provenance map. It explains possible trait co-carriage but does not show trait expression, within-sample activity, or statistical association. Recommended role: supplement next to the carrier bar plot.

### 6. [`pd_six_trait_direction_lollipop.png`](pd_six_trait_direction_lollipop.png)

**What it shows.** In the profile-weighted PD summary, moderately halophilic (`+0.54`), facultatively aerobic (`+0.45`), and arsenate reduction (`+0.39`) favor increased profiles; cellulose catabolic process (`-0.61`) strongly favors decreased profiles. Coccus shaped is near zero (`+0.02`) and sphere shaped is modestly negative (`-0.08`).

**Interpretation.** The diagram establishes the basic opposition that motivated the ecological-reporter analysis. It is descriptive and weights profiles rather than studies. This matters: in the paired equal-study-weighted analysis, sphere shaped changes sign to a small positive estimate (`+0.133`), and cellulose is `-0.651` in the 45-study condition screen or `-0.698` in the separate removal-influence input. These values answer related but different questions and should not be mixed.

**Takeaway.** Cellulose is the largest and most stable negative-direction trait across estimands, while the morphology traits are small and unstable. For the manuscript, prefer the study-weighted values in the ecological-reporter figure and use this lollipop only as a descriptive orientation panel.

### 7. [`pd_condition_taxon_trait_network.png`](pd_condition_taxon_trait_network.png)

**What it shows.** Panel A displays all 264 eligible condition-taxon tests; Panel B displays all 336 condition-trait tests; Panel C isolates all six PD trait tests. Filled points pass global `q<=0.05`, hollow amber points are exploratory at `0.05<q<=0.20`, and grey points show the remaining tested background.

**Interpretation.** The figure sharply narrows the defensible inferential claims. No condition-taxon result is FDR-significant. PD-*Roseburia* is exploratory (`q=0.192`). Among condition-trait tests, PD-moderately-halophilic is the only FDR-significant result (`+0.530`, `q=0.048`). PD-arsenate reduction is exploratory (`q=0.164`), whereas PD-cellulose (`-0.651`, `q=0.597`), facultatively aerobic (`+0.463`, `q=0.597`), coccus shaped (`q=0.889`), and sphere shaped (`q=0.889`) are not globally significant.

**Takeaway.** This is the principal inferential figure and should be cited whenever the multi-trait PD pattern is described. It supports one FDR discovery and several descriptive or exploratory components, not a globally significant PD-cellulose biomarker. Recommended role: high-priority supplement or main statistical companion.

### 8. [`pd_ecological_reporter_combined.png`](pd_ecological_reporter_combined.png)

**What it shows.** Panel A compares equal-study-weighted six-trait fingerprints and the ecological-reporter composite across 13 qualifying fecal condition labels. Panel B places those conditions on two summary axes: reversed cellulose contrast and the three-trait increased-side composite. Panel C maps reported PD genera and species to KG-Microbe annotations and BugSigDB direction support.

**Interpretation.** PD has the proposed opposing pattern: moderately halophilic `+0.53`, facultatively aerobic `+0.46`, arsenate reduction `+0.51`, cellulose `-0.65`, and composite `+1.15`, with `n=45` paired studies. Coccus and sphere morphology estimates are small. PD sits in the upper-right quadrant in Panel B, but colorectal cancer, COVID-19, stroke, graft-versus-host disease, ulcerative colitis, and other labels share parts of that region. The pattern therefore agrees with the ecological-reporter hypothesis but is not specific to PD. Panel C shows where trait annotations enter the score; it does not show activity or statistical mediation.

**Takeaway.** This is the best single main-text figure because it joins effect direction, cross-condition specificity, and taxonomic provenance. Its claim must stay at the level of a data-derived ecological-reporter hypothesis, with the FDR landscape cited as the inferential boundary.

### 9. [`pd_species_six_traits_influence_heatmap.png`](pd_species_six_traits_influence_heatmap.png)

**What it shows.** Five reported species pass the strict display filters: *Streptococcus salivarius*, *Methanobrevibacter smithii*, *Escherichia coli*, *Roseburia faecis*, and *Faecalibacterium prausnitzii*. Each cell is the full-data paired trait contrast minus the contrast after removing that species' complete non-MediaDive annotation vector. Nine of 30 cells exceed the display threshold of absolute influence greater than `0.01`.

**Interpretation.** All 30 audited species-trait intervals cross zero and are classified as mixed, so none of the aggregate PD trait directions can be assigned reliably to one displayed species. The largest visible direct cell is the *E. coli*-arsenate reduction influence (`+0.1012`), but its interval crosses zero. *R. faecis* has a negative cellulose influence (`-0.0302`) and *F. prausnitzii* has a negative facultatively-aerobic influence (`-0.1033`), again without supported direction. Four of the five BugSigDB recurrence directions agree with the targeted literature direction; *E. coli* differs, using a broader *Escherichia-Shigella* literature proxy.

**Takeaway.** The aggregate pattern is not attributable to a single displayed species under this sensitivity analysis. Cell values are non-additive influence diagnostics, not causal effects. Recommended role: mechanistic/provenance supplement.

### 10. [`pd_exploratory_taxa_six_traits_influence_heatmap.png`](pd_exploratory_taxa_six_traits_influence_heatmap.png)

**What it shows.** The expanded view lowers the support threshold to two studies and includes 30 of 51 eligible taxa: 8 species, 14 genera, and 8 families. It contains 180 audited cells, of which 71 exceed the absolute `0.01` display threshold. The audit classifications are 167 mixed, 6 decreased, 1 increased, and 6 low-support cells.

**Interpretation.** Rank expansion reveals candidate patterns hidden by the conservative species-only view, including strong-looking family or genus cells. It also increases ambiguity because a genus or family can aggregate descendants with different directions, and selection prioritizes literature-matched or well-supported rows. Even interval-excluding-zero cells in this relaxed view arise from an exploratory, selected analysis rather than a globally corrected taxon-by-trait screen.

**Takeaway.** Use this figure to generate follow-up taxa and rank-level hypotheses, not to name confirmed drivers. The strict species figure is the better source for conservative attribution statements. Recommended role: exploratory supplement only.

### 11. [`pd_trait_contributors_species_family_dot_heatmap.png`](pd_trait_contributors_species_family_dot_heatmap.png)

**What it shows.** The plot retains 11 of 53 candidate taxa - 3 species and 8 families - having at least one absolute increased-versus-decreased full-cell contributor-share difference of 5 percentage points. Red and blue circle areas show increased and decreased contributor fractions; the grey halo preserves unresolved or ambiguous mass. The right strip reports PD-minus-other-feces occurrence prevalence, normalized by 253 PD and 3,419 other fecal profiles.

**Interpretation.** Several taxa contribute different shares by direction, but the unresolved fraction is substantial for some cells, especially species-level sphere shaped (`80.3%` unresolved in increased and `70.6%` in decreased) and decreased arsenate reduction (`65.7%`). Lachnospiraceae has the largest displayed PD-versus-other-feces prevalence difference (`49/253` versus `307/3,419`, `+10.4` percentage points). These are selected propagated annotation counts, not independent observations, and the right strip is descriptive rather than study-clustered.

**Takeaway.** The figure supports distributed, rank-dependent contributor structure and makes missing resolution visible. It should not be used to test individual taxa or to sum parent and descendant labels. Recommended role: detailed provenance supplement.

### 12. [`pd_trait_contributors_trait_cocarriage_network.png`](pd_trait_contributors_trait_cocarriage_network.png)

**What it shows.** Only 3 of the 15 possible trait pairs change their selected-taxon co-carriage count between directions. Facultatively aerobic plus coccus shaped changes from 27 decreased to 28 increased shared taxa; facultatively aerobic plus moderately halophilic also changes from 27 to 28; coccus shaped plus moderately halophilic changes from 18 to 19. In all three cases, the one additional increased-direction family is Bifidobacteriaceae.

**Interpretation.** Co-carriage among the selected contributors is largely stable across direction. The displayed differences are one-taxon count changes, not evidence for broad network rewiring. The counts reflect shared annotations among a selected label set and do not establish that traits are expressed together within samples.

**Takeaway.** The network is more informative as a negative result: most pairwise co-carriage structure does not change, and the three visible changes are minimal. Recommended role: optional exploratory supplement; omit from the main paper unless the stability result is discussed explicitly.

## Overall analytical interpretation

The PD-cellulose work is strongest when presented as a layered case study rather than as a single-significance claim:

1. **Direction is the robust descriptive feature.** Cellulose repeatedly falls on the decreased-profile side, while moderately halophilic, facultatively aerobic, and arsenate reduction fall on the increased-profile side.
2. **Specificity is limited.** Other fecal conditions share the same two-axis region, and the formal PD-by-direction cellulose interaction is suggestive but not conventionally significant.
3. **Global inference is narrower than the visual pattern.** Only PD-moderately-halophilic passes the complete condition-trait FDR screen. PD-cellulose remains a coherent but globally non-significant component.
4. **Taxonomic provenance is distributed.** Multiple reported taxa and ranks carry the relevant annotations; no taxon survives the global condition-taxon screen, and the strict species-removal analysis supports no single species-level direction.
5. **The biological language should remain ecological.** The figures are consistent with a community trait pattern that reports a changed gut environment. They do not identify the host pressure, prove microbial activity, or establish causality.

The paper-level message is therefore not "PD causes loss of cellulose degradation." It is:

> A direction-aware trait synthesis exposes an opposing, community-level PD fecal pattern that taxon-only or direction-agnostic summaries would obscure, while global correction and influence diagnostics sharply limit the result to an ecological-reporter hypothesis rather than a biomarker or mechanism.
