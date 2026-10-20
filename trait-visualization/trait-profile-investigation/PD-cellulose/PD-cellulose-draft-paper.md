# PD-cellulose draft for *Microbiome trait profiles across health and disease*

## Purpose and interpretation boundary

This document proposes where the Parkinson's disease (PD)-cellulose analysis should enter the current manuscript and supplies manuscript-ready Results, Discussion, Methods, and Limitations text. It also provides a separate interpretation and takeaway for every diagram currently in this directory.

The central claim should remain deliberately narrow:

> PD fecal signatures show a reproducible descriptive opposition between cellulose catabolic process, which is stronger in decreased profiles, and a three-trait group comprising moderately halophilic, facultatively aerobic, and arsenate reduction, which is stronger in increased profiles. The pattern is direction-aware and distributed across taxa. It is not a validated diagnostic signature, a measurement of cellulose degradation, or evidence that a single organism causes the contrast.

The strongest globally adjusted result is PD-moderately-halophilic, not PD-cellulose. In the global screen of 408 eligible condition-trait tests, no PD trait is FDR-significant; moderately halophilic has the smallest global q and is exploratory (`q=0.0795`). The PD-cellulose contrast is directionally coherent and has a nominal paired-study `p=0.0070`, but it is not significant after the global correction (`q=0.477`). That distinction should be visible wherever the cellulose result is discussed.

Profile counts in this analysis are curated BugSigDB differential-abundance signature occurrences. They are not patients, sample prevalence, gene expression, metabolic flux, or direct cellulose-degradation measurements.

**Label sources.** All numbers here use each profile's own experiment condition and body site (#555), including the ecological-reporter tables behind Figure 7 (`pd_ecological_reporter_combined_*.tsv`), which `just pd-ecological-reporter-combined` now builds from the profile inputs (#585), and the species-removal influence, exploratory influence and contributor tables, which were regenerated with per-experiment labels and all six traits, each with its observed direction in the lollipop (#584).

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
| Global FDR-supported | None: no condition-trait or condition-taxon test passes global `q<=0.05` | No PD trait or taxon can be called globally significant. |
| Exploratory after global correction | PD-moderately-halophilic: mean paired-study contrast `+0.497`, raw `p=0.000195`, global `q=0.0795`, `n=47` paired studies; PD-arsenate reduction: `+0.520`, `q=0.0893`; PD-*Roseburia* genus direction: `q=0.116` | Follow-up hypotheses only; none is an FDR-significant discovery. |
| Directionally coherent but not globally significant | PD-facultatively-aerobic: `+0.360`, `q=0.860`; PD-cellulose: `-0.677`, `q=0.477` | These traits help describe the PD pattern, but should not be called globally significant. |
| Descriptive study-level composite | Three-trait increased-side composite `+0.459` (95% descriptive interval `0.194` to `0.724`); ecological-reporter composite `+1.136` (`0.538` to `1.734`); 47 paired studies | Pattern summary only. The composite is data-derived and is not a validated score or an independently corrected test. |
| Nominal PD-versus-other-feces comparison | Decreased cellulose-positive rate `58.0%` in PD versus `45.9%` in other feces; odds ratio `1.63`, Fisher `p=0.0131`; z-score Mann-Whitney `p=0.0072` | PD decreased profiles contain more cellulose signal than the same direction in other fecal conditions, before study clustering and global multiplicity control. |
| Study-clustered interaction | Cellulose z-score PD-by-decreased coefficient `+0.487`, `p=0.042`; positive-z interaction `+0.223`, `p=0.026` | Suggestive extra PD shift in the decreased direction at nominal `p<0.05`, without correction for the other tests in this case study; not an established interaction. |
| Annotation provenance | Recurrently reported cellulose-annotated taxa include Lachnospiraceae, *Roseburia*, *Blautia*, and Clostridiaceae; no condition-taxon test passes global `q<=0.05` | Taxa help trace where annotations enter the profile score. They are not inferential drivers or causal organisms. |

## Manuscript-ready Results subsection

### Parkinson's disease reveals a direction-aware ecological-reporter pattern

We used the six cross-framing convergent traits as fixed axes for a fecal Parkinson's disease case study. At the profile level, the PD subset contained 128 increased and 119 decreased profiles, contributed by 53 and 49 studies, respectively. Cellulose catabolic process was present at positive z-score in 34 increased profiles but 69 decreased profiles, and its mean z-score was `-0.033` in increased profiles compared with `0.668` in decreased profiles. The resulting profile-weighted increased-minus-decreased contrast was `-0.701`, making cellulose the clearest decreased-profile-skewed member of the six-trait set.

The paired, equal-study-weighted analysis gave the same directional conclusion without allowing studies containing more profiles to dominate. With each profile's own experiment labels, the global condition screen covers 47 PD studies represented in both directions. There, the mean cellulose contrast was `-0.677`; moderately halophilic (`+0.497`), facultatively aerobic (`+0.360`), and arsenate reduction (`+0.520`) pointed in the opposite direction; and coccus shaped (`+0.007`) and sphere shaped (`+0.090`) were near zero. The ecological-reporter tables behind Figure 7 hold the same study contrasts and add descriptive 95% t intervals: cellulose `-1.139` to `-0.214`, moderately halophilic `0.263` to `0.732`, facultatively aerobic `0.024` to `0.695`, and arsenate reduction `0.092` to `0.948`, while the intervals for coccus shaped and sphere shaped cross zero. We summarized the opposing pattern as an exploratory ecological-reporter composite: the mean contrast of moderately halophilic, facultatively aerobic, and arsenate reduction minus the cellulose contrast. In the reporter tables the PD composite was `+1.136` (descriptive 95% interval `0.538` to `1.734`). Because its component traits were selected from the PD pattern, this composite is a visualization aid rather than an independently validated score.

The cross-condition comparison did not support uniqueness to PD. Among the 16 fecal condition labels retained because at least one of the six individual traits had a descriptive study-level interval excluding zero (reporter tables), PD occupied the quadrant combining decreased-profile cellulose skew with an increased-side three-trait composite, but 12 other conditions occupied the same region; 9 of them had a larger ecological-reporter composite, and 10 a larger three-trait composite. The global evidence screen was still more restrictive. Across 408 eligible condition-trait tests, no result passed Benjamini-Hochberg correction. PD-moderately-halophilic (`q=0.0795`) and PD-arsenate reduction (`q=0.0893`) were exploratory, while PD-cellulose had `q=0.477` and PD-facultatively-aerobic `q=0.860`. Thus, the multi-trait display describes a coherent PD ecological-reporter pattern, but the pattern as a whole is neither PD-specific nor globally FDR-confirmed.

Within direction, PD decreased fecal profiles had a higher cellulose-positive rate than other decreased fecal profiles (`58.0%` versus `45.9%`; odds ratio `1.63`, Fisher `p=0.0131`), and their raw z-score distributions also differed nominally (Mann-Whitney `p=0.0072`). The converse contrast appeared among increased profiles, where the cellulose-positive rate was `26.6%` in PD and `36.8%` in other fecal profiles. A feces-only regression with study-clustered standard errors estimated an additional PD-associated shift in the decreased direction, with nominal `p<0.05` for both the continuous cellulose z-score (`+0.487`, `p=0.042`) and the positive-z indicator (`+0.223`, `p=0.026`); these p-values are not corrected for the other tests in this case study. These analyses support concentration of the PD cellulose signal in decreased profiles, but not a definitive PD-by-direction interaction.

The taxonomic provenance was distributed rather than dominated by one resolved species. Frequently recurring cellulose-annotated labels in PD decreased profiles included Lachnospiraceae (40 decreased versus 9 increased profiles), *Roseburia* (33 versus 4), *Blautia* (16 versus 2), and Clostridiaceae (13 versus 2). Other lineages showed different or balanced patterns, including Pseudomonadota (1 decreased versus 16 increased) and Oscillospiraceae (12 versus 15), indicating that cellulose annotation alone does not impose one direction. No condition-taxon comparison passed global FDR correction; the strongest taxon result, the *Roseburia* genus, remained exploratory (`q=0.116`). Species-removal influence diagnostics likewise did not assign any of the 30 retained species-trait cells a supported direction because every 95% interval crossed zero. Together, these results favor a community-level, annotation-mediated interpretation over attribution to a single microbial species.

**Results takeaway.** PD fecal signatures contain a direction-aware four-trait opposition: cellulose catabolic process is stronger in decreased profiles, while moderately halophilic, facultatively aerobic, and arsenate reduction are stronger in increased profiles. The signal is useful as an ecological-reporter hypothesis and as an example of information recovered by paired trait analysis. No component is globally FDR-supported in the current cross-condition screen (moderately halophilic and arsenate reduction are exploratory), so cellulose should remain a descriptive component of the pattern rather than a standalone PD biomarker.

## Suggested main-figure caption

**Figure 7. Parkinson's disease ecological-reporter pattern across fecal conditions, and taxon-to-trait provenance** (`just pd-ecological-reporter-combined`). **A**, Equal-study-weighted increased-minus-decreased contrasts for the six convergent traits across 16 human fecal condition labels having a descriptive 95% interval excluding zero for at least one individual trait. The ecological-reporter composite is the mean contrast of moderately halophilic, facultatively aerobic, and arsenate reduction minus the cellulose catabolic process contrast; the composite does not qualify a condition for inclusion. Dark borders mark descriptive individual-trait intervals excluding zero. **B**, The same conditions positioned by reversed cellulose contrast on the x-axis and the three-trait increased-side composite on the y-axis. PD lies in the upper-right quadrant but is not unique within it. Point area represents paired-study count. **C**, Exact genus and species labels reported in PD fecal BugSigDB signatures crossed with local KG-Microbe annotations for the six traits. Circle size is the number of increased plus decreased study reports (a study reporting the taxon in both directions counts twice) and color represents the log2 ratio of increased to decreased study support; an `x` indicates no local annotation. Panel C is annotation provenance, not a significance or activity map. All contrasts summarize curated differential-abundance signatures rather than patient-level measurements. All three panels read the ecological-reporter tables.

## Draft-paper figure description: condition-taxon-trait evidence landscape

**Figure X. Globally adjusted human fecal condition-taxon and condition-trait evidence landscape.** The analysis includes 68 conditions with at least three studies containing both increased and decreased BugSigDB signatures. Benjamini-Hochberg correction was applied separately across the complete condition-taxon and condition-trait test families. Filled points are FDR-significant (`q<=0.05`), hollow amber points are exploratory (`0.05<q<=0.20`) and are not claimed as discoveries, and grey points show the remaining tested background. **A**, Results for all 287 eligible condition-taxon comparisons. The x-axis is the directional study fraction, with negative values indicating more decreased-only than increased-only studies; the y-axis is `-log10(global q)`. Circles denote genera, triangles denote species, and point area represents the number of reporting studies. No condition-taxon comparison passes the FDR threshold. The strongest result, decreased reporting of the *Roseburia* genus in PD, remains exploratory (directional fraction `-0.340`, `q=0.116`). **B**, Results for all 408 eligible condition-trait comparisons. The x-axis is the mean paired study-level increased-minus-decreased z-score contrast, the y-axis is `-log10(global q)`, and point area represents the number of paired studies. No condition-trait comparison passes the FDR threshold. PD-moderately-halophilic (mean contrast `+0.497`, raw `p=0.000195`, `q=0.0795`, `n=47`), PD-arsenate reduction (`+0.520`, `q=0.0893`), COVID-19-cellulose (`-1.298`, `q=0.0893`), colorectal-cancer-cellulose (`-0.851`, `q=0.0893`), and colorectal-cancer-moderately-halophilic (`+0.448`, `q=0.0942`) are exploratory. **C**, All six PD trait contrasts, with unadjusted `p`, globally adjusted `q`, and paired-study count printed beside each estimate. Cellulose catabolic process is stronger in decreased PD profiles (contrast `-0.677`, raw `p=0.0070`) but is not significant after global correction (`q=0.477`); moderately halophilic and arsenate reduction are exploratory, and no PD trait is FDR-supported. The figure therefore supports a direction-aware multi-trait PD pattern without a formal discovery claim for any single component; taxon points describe directional reporting in curated signatures and do not establish abundance effects, microbial activity, or causality.

## Manuscript-ready Discussion insertion

The PD case study shows why condition-level trait interpretation should preserve direction and distinguish pattern description from multiplicity-controlled discovery. PD fecal signatures combined a decreased-profile cellulose component with increased-profile moderately-halophilic, facultatively-aerobic, and arsenate-reduction components. This opposition was stable across profile-weighted and equal-study-weighted summaries, but it did not define a unique PD region among fecal conditions and no component survived the global condition-trait FDR correction (moderately halophilic and arsenate reduction remained exploratory). The appropriate biological interpretation is therefore an ecological-reporter hypothesis: taxon composition and its annotated trait potential may jointly reflect a changed gut environment associated with PD, including possible differences in diet, intestinal transit, medication, inflammation, or oxygen and redox conditions, without identifying which of those pressures generated the pattern.

The distributed taxonomic provenance reinforces this interpretation. Cellulose signal was associated with several recurrent, partly overlapping labels rather than one supported species-level influence, and no condition-taxon result survived global correction. The trait contrast may therefore be more portable than any individual taxon association, but it should not be equated with measured cellulose degradation or metabolic activity. More generally, the PD example demonstrates the value of the paper's combined reporting convention: a direction-aware composition analysis can expose a biologically interpretable multi-trait configuration, while the global inferential layer states which individual components are robust enough to claim as discoveries.

## Manuscript-ready Methods insertion

### PD fecal case-study analysis

The PD case study was restricted to human fecal BugSigDB profiles and to the six traits in the cross-framing convergent core. A profile-level descriptive summary compared mean composition z-scores and positive-z rates between increased and decreased PD signatures. For study-weighted analyses, increased and decreased profile z-scores were first averaged within each study and trait, and the paired contrast was defined as mean increased z minus mean decreased z. Each study therefore contributed one contrast irrespective of its number of experiments or profiles. Descriptive 95% Student t intervals were calculated across study contrasts.

The exploratory three-trait increased-side composite was the mean of the moderately-halophilic, facultatively-aerobic, and arsenate-reduction contrasts. The ecological-reporter composite was that three-trait mean minus the cellulose catabolic process contrast. These components were selected because they formed the interval-excluding-zero positive and negative sides of the PD pattern; neither composite was used to qualify conditions for the cross-condition display or treated as a validated diagnostic score. Conditions entered that display only when they had at least three paired studies and at least one of the six individual traits had a descriptive 95% study-level interval excluding zero.

For global condition-trait inference, each eligible human fecal condition required at least three studies containing both increased and decreased signatures. Each condition-trait distribution of paired study contrasts was tested against zero with a two-sided Wilcoxon signed-rank test. Benjamini-Hochberg correction was applied across the full family of 408 eligible condition-trait tests. Exact genus and species identifiers were kept separate in the condition-taxon analysis; taxa required at least three reporting studies and were tested with an exact McNemar/binomial test on discordant paired studies, followed by a separate global Benjamini-Hochberg correction across 287 eligible condition-taxon tests.

The PD-versus-other-feces sensitivity analysis compared positive-z rates with two-sided Fisher exact tests and raw z-score distributions with two-sided Mann-Whitney tests within direction. A feces-only ordinary least-squares model included indicators for PD, decreased direction, and their interaction, with standard errors clustered by study. The interaction tested whether decreased profiles had an additional PD-associated cellulose shift beyond the main PD and direction effects. Annotated-taxon, removal-influence, contributor, and co-carriage analyses were treated as annotation-provenance or sensitivity analyses rather than independent inferential tests.

## Manuscript-ready Limitations insertion

The PD case study adds further interpretation limits. Its units are curated differential-abundance signatures rather than patients, and repeated profiles within studies are not independent unless the analysis explicitly collapses or clusters them by study. KG-Microbe annotations indicate reported trait potential, not gene expression, substrate use, or metabolite flux; consequently, a cellulose catabolic process score is not a direct measurement of cellulose degradation. The ecological-reporter composite was constructed from the observed PD pattern and compared descriptively across conditions, so it requires external, prospectively specified validation before any diagnostic or specificity claim. Medication, diet, intestinal transit, geography, sequencing method, and other cohort differences can also generate or modify the observed pattern.

## Diagram-by-diagram interpretations and takeaways

### 1. [`pd_cellulose_direction_specificity.png`](pd_cellulose_direction_specificity.png)

**What it shows.** The left panel compares the two PD directions. There are 34 cellulose-positive increased PD fecal profiles from 22 studies and 69 cellulose-positive decreased profiles from 34 studies. The cellulose focus-minus-opposite-direction z contrast is larger for the decreased edge (`1.61`) than for the increased edge (`0.66`). The right panel compares PD with all other fecal profiles within direction: `26.6%` versus `36.8%` are cellulose-positive among increased profiles, while `58.0%` versus `45.9%` are positive among decreased profiles.

**Interpretation.** The cellulose signal is asymmetric. It is concentrated in PD decreased profiles rather than elevated across both PD directions. The decreased-direction PD-versus-other-feces comparison has nominal Fisher `p<0.05` (odds ratio `1.63`, `p=0.0131`), but this panel does not account for repeated profiles within studies or global multiple testing.

**Takeaway.** Use this as the clearest descriptive overview of the PD-cellulose observation. Phrase the result as "decreased-profile concentration," not "PD-specific enrichment" or "biomarker." Recommended role: supplementary descriptive figure.

### 2. [`pd_cellulose_cytoscape_highlight.png`](pd_cellulose_cytoscape_highlight.png)

**What it shows.** The orange PD-to-cellulose-to-feces path is highlighted against faded links for the other five shared traits. The decreased view has the thicker cellulose path and the larger edge contrast (`n=69`, `s=34`, `delta z=1.61`) compared with the increased view (`n=34`, `s=22`, `delta z=0.66`).

**Interpretation.** This is a network rendering of the same per-direction counts and edge scores shown elsewhere. It contributes visual orientation, not an independent statistical result.

**Takeaway.** The network makes the path from condition to trait to body site intuitive, but the thicker decreased edge must not be interpreted as causal or statistically significant by itself. Recommended role: graphical overview or supplement; omit if space is limited because it is redundant with the by-direction bar figure.

### 3. [`pd_cellulose_enrichment_test.png`](pd_cellulose_enrichment_test.png)

**What it shows.** The feces-only group means are `0.10` for other increased, `-0.03` for PD increased, `0.31` for other decreased, and `0.67` for PD decreased profiles. In the study-clustered model, the PD effect among increased profiles is `-0.13` (`p=0.327`), the decreased effect among other fecal profiles is `+0.21` (`p=0.000098`), and the extra PD effect in decreased profiles is `+0.49` (`p=0.042`).

**Interpretation.** Some cellulose elevation in decreased profiles is a broader fecal-signature phenomenon, not uniquely a PD phenomenon. PD appears to add to that decreased-direction effect; the PD-by-decreased interaction interval (`0.02` to `0.96`) excludes zero at a nominal level, without correction for the other tests in this case study. The binary positive-z model tells the same story (`+0.223`, `p=0.026`).

**Takeaway.** This is the key guardrail against overclaiming specificity. The PD direction interaction is suggestive, not established. Recommended role: statistical supplementary figure accompanying the descriptive direction panel.

### 4. [`pd_cellulose_taxon_driver_bar.png`](pd_cellulose_taxon_driver_bar.png)

**What it shows.** The diagram counts unique PD fecal profiles containing reported taxa that carry a local cellulose annotation. Lachnospiraceae (`40` decreased, `9` increased), *Roseburia* (`33`, `4`), *Blautia* (`16`, `2`), Clostridiaceae (`13`, `2`), *Clostridium* (`10`, `3`), *Ruminococcus* (`9`, `2`), and *Lachnospira* (`9`, `2`) favor decreased profiles. Pseudomonadota (`1`, `16`) and Actinomycetota (`0`, `8`) favor increased profiles, while Oscillospiraceae (`12`, `15`) is nearly balanced.

**Interpretation.** The aggregate cellulose contrast is assembled from multiple, directionally heterogeneous taxonomic labels. The strong decreased recurrence of Lachnospiraceae and *Roseburia* is consistent with the bulk cellulose pattern, but the opposing lineages show that a cellulose annotation does not dictate one disease direction. Exact labels span species, genus, family, order, class, and phylum and therefore cannot be added as disjoint contributions.

**Takeaway.** Refer to these rows as recurrently reported taxa annotated with cellulose catabolic process, not as causal "drivers" or as organisms shown to degrade cellulose. *Roseburia faecis* (`6` decreased, `0` increased profiles) is a follow-up candidate, not a validated species biomarker. Recommended role: supplementary provenance figure.

### 5. [`pd_cellulose_taxon_six_trait_heatmap.png`](pd_cellulose_taxon_six_trait_heatmap.png)

**What it shows.** The rows are the recurrently reported cellulose-annotated labels from the preceding bar plot, and the columns are all six convergent traits. Orange cells denote a local KG-Microbe annotation; printed values are annotation counts. The heatmap shows that some cellulose-annotated taxa are annotated only narrowly, such as *Roseburia* and *Roseburia faecis*, whereas broad labels such as Bacillota, Pseudomonadota, Clostridia, Bacilli, and Eubacteriales carry several of the six annotations.

**Interpretation.** The PD-cellulose observation is embedded in overlapping trait potential. Some reported taxa can influence multiple trait scores, so the six PD contrasts should not be interpreted as six taxonomically independent signals. Very large annotation counts for broad ranks reflect graph coverage and descendant propagation rather than abundance or effect strength.

**Takeaway.** This figure is an annotation-provenance map. It explains possible trait co-carriage but does not show trait expression, within-sample activity, or statistical association. Recommended role: supplement next to the carrier bar plot.

### 6. [`pd_six_trait_direction_lollipop.png`](pd_six_trait_direction_lollipop.png)

**What it shows.** In the profile-weighted PD summary, moderately halophilic (`+0.49`), facultatively aerobic (`+0.35`), and arsenate reduction (`+0.42`) favor increased profiles; cellulose catabolic process (`-0.70`) strongly favors decreased profiles. Coccus shaped (`-0.08`) and sphere shaped (`-0.13`) are modestly negative.

**Interpretation.** The diagram establishes the basic opposition that motivated the ecological-reporter analysis. It is descriptive and weights profiles rather than studies. This matters: in the paired equal-study-weighted condition screen, which uses per-experiment labels, coccus shaped and sphere shaped change sign to small positive estimates (`+0.007` and `+0.090`), and cellulose is `-0.677` across 47 paired studies; the ecological-reporter tables behind Figure 7 hold the same study contrasts. The species-removal influence analysis, also on per-experiment labels but with its own fixed human reference, gives a cellulose baseline of `-0.691`. These values come from different analyses and answer related but different questions, so they should not be mixed.

**Takeaway.** Cellulose is the largest and most stable negative-direction trait across estimands, while the morphology traits are small and unstable. For the manuscript, prefer study-weighted values (the condition screen and the ecological-reporter figure) and use this lollipop only as a descriptive orientation panel.

### 7. [`pd_condition_taxon_trait_network.png`](pd_condition_taxon_trait_network.png)

**What it shows.** Panel A displays all 287 eligible condition-taxon tests; Panel B displays all 408 condition-trait tests; Panel C isolates all six PD trait tests. Filled points pass global `q<=0.05`, hollow amber points are exploratory at `0.05<q<=0.20`, and grey points show the remaining tested background.

**Interpretation.** The figure sharply narrows the defensible inferential claims. No condition-taxon result is FDR-significant. PD-*Roseburia* is exploratory (`q=0.116`). No condition-trait result is FDR-significant either: PD-moderately-halophilic (`+0.497`, `q=0.0795`) and PD-arsenate reduction (`+0.520`, `q=0.0893`) are exploratory, whereas PD-cellulose (`-0.677`, `q=0.477`), facultatively aerobic (`+0.360`, `q=0.860`), coccus shaped (`q=1.000`), and sphere shaped (`q=1.000`) are not globally significant.

**Takeaway.** This is the principal inferential figure and should be cited whenever the multi-trait PD pattern is described. It supports no FDR discovery, two exploratory components and several descriptive ones, not a globally significant PD-cellulose biomarker. Recommended role: high-priority supplement or main statistical companion.

### 8. [`pd_ecological_reporter_combined.png`](pd_ecological_reporter_combined.png)

**What it shows.** Panel A compares equal-study-weighted six-trait contrasts and the ecological-reporter composite across 16 qualifying fecal condition labels. Panel B places those conditions on two summary axes: reversed cellulose contrast and the three-trait increased-side composite. Panel C maps reported PD genera and species to KG-Microbe annotations and BugSigDB direction support.

**Interpretation.** PD has the proposed opposing pattern: moderately halophilic `+0.50`, facultatively aerobic `+0.36`, arsenate reduction `+0.52`, cellulose `-0.68`, and composite `+1.14`, with `n=47` paired studies. Coccus and sphere morphology estimates are small. PD sits in the upper-right quadrant in Panel B, but colorectal cancer, COVID-19, stroke, graft-versus-host disease, chronic kidney disease, and seven other labels share that region. The pattern therefore agrees with the ecological-reporter hypothesis but is not specific to PD. Panel C shows where trait annotations enter the score; it does not show activity or statistical mediation.

**Takeaway.** This is the best single main-text figure because it joins effect direction, the cross-condition comparison, and taxonomic provenance. Its claim must stay at the level of a data-derived ecological-reporter hypothesis, with the FDR landscape cited as the inferential boundary.

### 9. [`pd_species_six_traits_influence_heatmap.png`](pd_species_six_traits_influence_heatmap.png)

**What it shows.** Five reported species pass the strict display filters: *Bifidobacterium longum*, *Methanobrevibacter smithii*, *Escherichia coli*, *Roseburia faecis*, and *Faecalibacterium prausnitzii*. Each cell is the full-data paired trait contrast minus the contrast after removing that species' complete non-MediaDive annotation vector. Eight of 30 cells exceed the display threshold of absolute influence greater than `0.01`. Each column carries the trait's observed direction in the six-trait lollipop: moderately halophilic, arsenate reduction and facultatively aerobic are higher in increased profiles; coccus shaped, sphere shaped and cellulose catabolic process are higher in decreased profiles.

**Interpretation.** All 30 audited species-trait intervals cross zero and are classified as mixed, so none of the aggregate PD trait directions can be assigned reliably to one displayed species. The largest visible direct cells are the *F. prausnitzii*-facultatively aerobic influence (`-0.0986`) and the *E. coli*-arsenate reduction influence (`+0.0966`), but both intervals cross zero. *R. faecis* has a negative cellulose influence (`-0.0288`), again without supported direction. Four of the five BugSigDB recurrence directions agree with the targeted literature direction; *E. coli* differs, using a broader *Escherichia-Shigella* literature proxy.

**Takeaway.** The aggregate pattern is not attributable to a single displayed species under this sensitivity analysis. Cell values are non-additive influence diagnostics, not causal effects. Recommended role: mechanistic/provenance supplement.

### 10. [`pd_exploratory_taxa_six_traits_influence_heatmap.png`](pd_exploratory_taxa_six_traits_influence_heatmap.png)

**What it shows.** The expanded view lowers the support threshold to two studies and includes 30 of 51 eligible taxa: 10 species, 12 genera, and 8 families. It contains 180 audited cells, of which 68 exceed the absolute `0.01` display threshold. The audit classifications are 161 mixed, 6 decreased, 1 increased, and 12 low-support cells.

**Interpretation.** Rank expansion reveals candidate patterns hidden by the conservative species-only view, including strong-looking family or genus cells. It also increases ambiguity because a genus or family can aggregate descendants with different directions, and selection prioritizes literature-matched or well-supported rows. Even interval-excluding-zero cells in this relaxed view arise from an exploratory, selected analysis rather than a globally corrected taxon-by-trait screen.

**Takeaway.** Use this figure to generate follow-up taxa and rank-level hypotheses, not to name confirmed drivers. The strict species figure is the better source for conservative attribution statements. Recommended role: exploratory supplement only.

### 11. [`pd_trait_contributors_species_family_dot_heatmap.png`](pd_trait_contributors_species_family_dot_heatmap.png)

**What it shows.** The contributor counts cover the 198 PD fecal profiles from human experiments with both directions (99 increased, 99 decreased, 45 studies). The plot retains 11 of 53 candidate taxa - 3 species and 8 families - having at least one absolute increased-versus-decreased full-cell contributor-share difference of 5 percentage points. Red and blue circle areas show increased and decreased contributor fractions; the grey halo preserves unresolved or ambiguous mass. Trait labels carry each trait's observed direction in the lollipop. The right strip reports PD-minus-other-feces occurrence prevalence, normalized by 247 PD and 3,357 other fecal profiles.

**Interpretation.** Several taxa contribute different shares by direction, but the unresolved fraction is substantial for some cells, especially species-level sphere shaped (`80.3%` unresolved in increased and `70.6%` in decreased) and decreased arsenate reduction (`65.7%`). Lachnospiraceae has the largest displayed PD-versus-other-feces prevalence difference (`49/247` versus `307/3,357`, `+10.7` percentage points). These are selected propagated annotation counts, not independent observations, and the right strip is descriptive rather than study-clustered.

**Takeaway.** The figure supports distributed, rank-dependent contributor structure and makes missing resolution visible. It should not be used to test individual taxa or to sum parent and descendant labels. Recommended role: detailed provenance supplement.

### 12. [`pd_trait_contributors_trait_cocarriage_network.png`](pd_trait_contributors_trait_cocarriage_network.png)

**What it shows.** Only 3 of the 15 possible trait pairs change their selected-taxon co-carriage count between directions. Facultatively aerobic plus coccus shaped changes from 27 decreased to 28 increased shared taxa; facultatively aerobic plus moderately halophilic also changes from 27 to 28; coccus shaped plus moderately halophilic changes from 18 to 19. In all three cases, the one additional increased-direction family is Bifidobacteriaceae.

**Interpretation.** Co-carriage among the selected contributors is largely stable across direction. The displayed differences are one-taxon count changes, not evidence for broad network rewiring. The counts reflect shared annotations among a selected label set and do not establish that traits are expressed together within samples.

**Takeaway.** The network is more informative as a negative result: most pairwise co-carriage structure does not change, and the three visible changes are minimal. Recommended role: optional exploratory supplement; omit from the main paper unless the stability result is discussed explicitly.

## Overall analytical interpretation

The PD-cellulose work is strongest when presented as a layered case study rather than as a single-significance claim:

1. **Direction is the robust descriptive feature.** Cellulose repeatedly falls on the decreased-profile side, while moderately halophilic, facultatively aerobic, and arsenate reduction fall on the increased-profile side.
2. **Specificity is limited.** Other fecal conditions share the same two-axis region, and the formal PD-by-direction cellulose interaction reaches only nominal, uncorrected `p<0.05` (`p=0.042` for the z-score, `p=0.026` for the positive-z indicator), so it is suggestive rather than established.
3. **Global inference is narrower than the visual pattern.** No PD trait passes the complete condition-trait FDR screen; PD-moderately-halophilic and PD-arsenate reduction are exploratory (`q=0.0795` and `q=0.0893`). PD-cellulose remains a coherent but globally non-significant component.
4. **Taxonomic provenance is distributed.** Multiple reported taxa and ranks carry the relevant annotations; no taxon survives the global condition-taxon screen, and the strict species-removal analysis supports no single species-level direction.
5. **The biological language should remain ecological.** The figures are consistent with a community trait pattern that reports a changed gut environment. They do not identify the host pressure, prove microbial activity, or establish causality.

The paper-level message is therefore not "PD causes loss of cellulose degradation." It is:

> A direction-aware trait synthesis exposes an opposing, community-level PD fecal pattern that taxon-only or direction-agnostic summaries would obscure, while global correction and influence diagnostics sharply limit the result to an ecological-reporter hypothesis rather than a biomarker or mechanism.
