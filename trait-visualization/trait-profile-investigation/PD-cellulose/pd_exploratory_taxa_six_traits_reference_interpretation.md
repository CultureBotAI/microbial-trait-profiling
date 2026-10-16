# PD Taxa Reference Interpretation

This table joins the displayed species influence rows to DOI-backed PD taxa literature, BugSigDB PD occurrence direction, KG-Microbe six-trait annotations, and the strongest current influence value.

Rows follow the heatmap order: descending BugSigDB PD increased/decreased log2 profile-count ratio, with count-based and alphabetical tie-breakers.

## Agreement Summary

| Taxon | BugSigDB inc/dec | BugSigDB direction | Literature value | Match scope | Hamming distance | Agreement | Strongest six-trait influence |
| --- | ---: | --- | --- | --- | ---: | --- | --- |
| Faecalibacterium prausnitzii | 0/13 | decreased | decreased | Faecalibacterium; genus/species; exact | 0 | agrees | facultatively aerobic (-0.10; mixed) |
| Roseburia faecis | 0/6 | decreased | decreased | Roseburia; genus; genus_or_parent | 0 | agrees | cellulose catabolic process (-0.03; mixed) |
| Streptococcus salivarius | 4/0 | increased | increased | Streptococcus; genus; genus_or_parent | 0 | agrees | facultatively aerobic (+0.01; mixed) |
| Escherichia coli | 1/2 | decreased | increased | Escherichia-Shigella; genus group; exact | 1 | differs | arsenate reduction (+0.10; mixed) |
| Methanobrevibacter smithii | 3/0 | increased | increased | Methanobrevibacter; genus; exact | 0 | agrees | coccus shaped (+0.03; mixed) |
| Streptococcus mutans | 5/0 | increased | increased | Streptococcus; genus; genus_or_parent | 0 | agrees | coccus shaped (+0.05; mixed) |
| Blautia obeum | 5/0 | increased | decreased | Blautia; genus; genus_or_parent | 1 | differs | cellulose catabolic process (+0.12; mixed) |
| Ruminococcus gauvreauii | 0/2 | decreased | decreased | Ruminococcus; genus; genus_or_parent | 0 | agrees | sphere shaped (-0.03; low support) |
| Roseburia | 4/31 | decreased | decreased | Roseburia; genus; exact | 0 | agrees | cellulose catabolic process (-0.15; mixed) |
| Faecalibacterium | 4/32 | decreased | decreased | Faecalibacterium; genus/species; exact | 0 | agrees | facultatively aerobic (-0.09; dec) |
| Lactobacillus | 15/6 | increased | increased | Lactobacillus; genus; exact | 0 | agrees | facultatively aerobic (+0.06; mixed) |
| Streptococcus | 11/5 | increased | increased | Streptococcus; genus; exact | 0 | agrees | coccus shaped (+0.10; mixed) |
| Blautia | 2/15 | decreased | decreased | Blautia; genus; exact | 0 | agrees | coccus shaped (-0.07; mixed) |
| Bifidobacterium | 21/2 | increased | increased | Bifidobacterium; genus; exact | 0 | agrees | moderately halophilic (-0.01; mixed) |
| Coprococcus | 1/10 | decreased | decreased | Coprococcus; genus; exact | 0 | agrees | coccus shaped (-0.01; mixed) |
| Escherichia | 9/2 | increased | increased | Escherichia-Shigella; genus group; exact | 0 | agrees | arsenate reduction (+0.02; mixed) |
| Ruminococcus | 2/9 | decreased | decreased | Ruminococcus; genus; exact | 0 | agrees | cellulose catabolic process (-0.14; mixed) |
| Enterococcus | 7/2 | increased | increased | Enterococcus; genus; exact | 0 | agrees | facultatively aerobic (-0.01; mixed) |
| Bacteroides | 3/6 | decreased | decreased | Bacteroides fragilis group; species group; synonym | 0 | agrees | cellulose catabolic process (+0.01; mixed) |
| Clostridium | 3/10 | decreased | not listed | not listed; not listed; none | NA | not in literature reference | coccus shaped (+0.02; mixed) |
| Butyrivibrio | 2/3 | decreased | not listed | not listed; not listed; none | NA | not in literature reference | cellulose catabolic process (-0.04; mixed) |
| Butyricicoccus | 1/6 | decreased | not listed | not listed; not listed; none | NA | not in literature reference | sphere shaped (-0.04; mixed) |
| Lachnospiraceae | 10/39 | decreased | decreased | Lachnospiraceae; family; exact | 0 | agrees | cellulose catabolic process (-0.13; dec) |
| Enterobacteriaceae | 18/2 | increased | increased | Enterobacteriaceae; family; exact | 0 | agrees | arsenate reduction (+0.20; mixed) |
| Lactobacillaceae | 12/5 | increased | increased | Lactobacillus; genus; exact | 0 | agrees | arsenate reduction (-0.03; mixed) |
| Bifidobacteriaceae | 16/0 | increased | increased | Bifidobacterium; genus; exact | 0 | agrees | coccus shaped (+0.01; mixed) |
| Oscillospiraceae | 15/13 | increased | not listed | not listed; not listed; none | NA | not in literature reference | cellulose catabolic process (+0.07; mixed) |
| Rikenellaceae | 7/0 | increased | not listed | not listed; not listed; none | NA | not in literature reference | moderately halophilic (+0.04; mixed) |
| Veillonellaceae | 3/3 | balanced | not listed | not listed; not listed; none | NA | not in literature reference | sphere shaped (-0.02; mixed) |
| Peptostreptococcaceae | 5/3 | increased | not listed | not listed; not listed; none | NA | not in literature reference | sphere shaped (-0.02; mixed) |

## Methods: Literature Review and Direction Comparison

The literature direction table is a targeted narrative review updated 2026-09-22; it is not a systematic review or meta-analysis. PubMed, PubMed Central, and publisher records were searched using `Parkinson disease`, `gut microbiome`, each taxon name or synonym, and direction terms (`increased`, `decreased`, `enriched`, or `depleted`). Existing DOI-backed PD microbiome references were checked backward for additional direction evidence.

Eligible evidence comprised human PD-versus-control fecal or intestinal microbiome studies and cross-cohort human meta-analyses that reported an explicit taxon direction. Animal-only, intervention-only, non-gut, and severity-only findings were excluded from increased/decreased calls. Association-only evidence is retained as `associated` and is not treated as directional.

Literature direction is assigned at the most specific supported taxonomic level. An exact species result is preferred; otherwise a genus or family result is used as a flagged parent-level proxy. `variable` is used when eligible reports conflict at the matched level. The reference table records the matched term, reference rank, match type, evidence statement, and DOI list so broader-rank comparisons remain visible.

BugSigDB direction is a recurrence summary, not a pooled abundance effect: `increased` when increased-profile count is greater than decreased-profile count, `decreased` for the reverse, and `balanced` for a tie. `agrees` and `differs` are assigned only against directional literature calls. For those calls, literature and BugSigDB direction are each encoded as one binary categorical position and compared with Hamming distance: 0 means the directions are similar and 1 means they differ. Balanced, variable, associated, and unlisted cases have no defined Hamming distance and remain separate. This is a descriptive direction check, not a statistical test of agreement, and it does not resolve cohort, medication, geography, sequencing, or taxonomic-version heterogeneity.

## Taxon Details

### Faecalibacterium prausnitzii
- Literature match: Faecalibacterium (decreased; genus/species; exact); source(s): [10.1016/j.parkreldis.2016.08.019](https://doi.org/10.1016/j.parkreldis.2016.08.019); [10.1007/s11427-016-9001-4](https://doi.org/10.1007/s11427-016-9001-4); [10.1002/mds.28119](https://doi.org/10.1002/mds.28119)
- BugSigDB PD direction: decreased (0/13 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: facultatively aerobic
- Interpretation: BugSigDB PD occurrence agrees with literature (decreased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: facultatively aerobic. Strongest current influence: facultatively aerobic (-0.10 fixed-reference z, mixed).

### Roseburia faecis
- Literature match: Roseburia (decreased; genus; genus_or_parent); source(s): [10.1002/mds.26307](https://doi.org/10.1002/mds.26307); [10.1002/mds.28119](https://doi.org/10.1002/mds.28119); [10.1038/s41467-022-34667-x](https://doi.org/10.1038/s41467-022-34667-x)
- BugSigDB PD direction: decreased (0/6 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: cellulose catabolic process
- Interpretation: BugSigDB PD occurrence agrees with literature (decreased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: cellulose catabolic process. Strongest current influence: cellulose catabolic process (-0.03 fixed-reference z, mixed).

### Streptococcus salivarius
- Literature match: Streptococcus (increased; genus; genus_or_parent); source(s): [10.1007/s11427-016-9001-4](https://doi.org/10.1007/s11427-016-9001-4); [10.1038/s41467-022-34667-x](https://doi.org/10.1038/s41467-022-34667-x)
- BugSigDB PD direction: increased (4/0 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: coccus shaped; facultatively aerobic
- Interpretation: BugSigDB PD occurrence agrees with literature (increased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: coccus shaped; facultatively aerobic. Strongest current influence: facultatively aerobic (+0.01 fixed-reference z, mixed).

### Escherichia coli
- Literature match: Escherichia-Shigella (increased; genus group; exact); source(s): [10.1007/s11427-016-9001-4](https://doi.org/10.1007/s11427-016-9001-4)
- BugSigDB PD direction: decreased (1/2 increased/decreased profiles); Hamming distance: 1 (different); agreement: differs
- KG-Microbe direct six-trait annotations: arsenate reduction; coccus shaped; facultatively aerobic; moderately halophilic
- Interpretation: BugSigDB PD occurrence is decreased, but the literature reference is increased; review cohort, taxonomy, and direction labels before interpreting. KG-Microbe direct six-trait annotations: arsenate reduction; coccus shaped; facultatively aerobic; moderately halophilic. Strongest current influence: arsenate reduction (+0.10 fixed-reference z, mixed).

### Methanobrevibacter smithii
- Literature match: Methanobrevibacter (increased; genus; exact); source(s): [10.1038/s41467-022-34667-x](https://doi.org/10.1038/s41467-022-34667-x)
- BugSigDB PD direction: increased (3/0 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: coccus shaped
- Interpretation: BugSigDB PD occurrence agrees with literature (increased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: coccus shaped. Strongest current influence: coccus shaped (+0.03 fixed-reference z, mixed).

### Streptococcus mutans
- Literature match: Streptococcus (increased; genus; genus_or_parent); source(s): [10.1007/s11427-016-9001-4](https://doi.org/10.1007/s11427-016-9001-4); [10.1038/s41467-022-34667-x](https://doi.org/10.1038/s41467-022-34667-x)
- BugSigDB PD direction: increased (5/0 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: coccus shaped; facultatively aerobic; moderately halophilic
- Interpretation: BugSigDB PD occurrence agrees with literature (increased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: coccus shaped; facultatively aerobic; moderately halophilic. Strongest current influence: coccus shaped (+0.05 fixed-reference z, mixed).

### Blautia obeum
- Literature match: Blautia (decreased; genus; genus_or_parent); source(s): [10.1002/mds.26307](https://doi.org/10.1002/mds.26307); [10.1007/s11427-016-9001-4](https://doi.org/10.1007/s11427-016-9001-4)
- BugSigDB PD direction: increased (5/0 increased/decreased profiles); Hamming distance: 1 (different); agreement: differs
- KG-Microbe direct six-trait annotations: cellulose catabolic process; coccus shaped
- Interpretation: BugSigDB PD occurrence is increased, but the literature reference is decreased; review cohort, taxonomy, and direction labels before interpreting. KG-Microbe direct six-trait annotations: cellulose catabolic process; coccus shaped. Strongest current influence: cellulose catabolic process (+0.12 fixed-reference z, mixed).

### Ruminococcus gauvreauii
- Literature match: Ruminococcus (decreased; genus; genus_or_parent); source(s): [10.1007/s11427-016-9001-4](https://doi.org/10.1007/s11427-016-9001-4)
- BugSigDB PD direction: decreased (0/2 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: coccus shaped; sphere shaped
- Interpretation: BugSigDB PD occurrence agrees with literature (decreased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: coccus shaped; sphere shaped. Strongest current influence: sphere shaped (-0.03 fixed-reference z, low support).

### Roseburia
- Literature match: Roseburia (decreased; genus; exact); source(s): [10.1002/mds.26307](https://doi.org/10.1002/mds.26307); [10.1002/mds.28119](https://doi.org/10.1002/mds.28119); [10.1038/s41467-022-34667-x](https://doi.org/10.1038/s41467-022-34667-x)
- BugSigDB PD direction: decreased (4/31 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: cellulose catabolic process
- Interpretation: BugSigDB PD occurrence agrees with literature (decreased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: cellulose catabolic process. Strongest current influence: cellulose catabolic process (-0.15 fixed-reference z, mixed).

### Faecalibacterium
- Literature match: Faecalibacterium (decreased; genus/species; exact); source(s): [10.1016/j.parkreldis.2016.08.019](https://doi.org/10.1016/j.parkreldis.2016.08.019); [10.1007/s11427-016-9001-4](https://doi.org/10.1007/s11427-016-9001-4); [10.1002/mds.28119](https://doi.org/10.1002/mds.28119)
- BugSigDB PD direction: decreased (4/32 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: facultatively aerobic
- Interpretation: BugSigDB PD occurrence agrees with literature (decreased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: facultatively aerobic. Strongest current influence: facultatively aerobic (-0.09 fixed-reference z, dec).

### Lactobacillus
- Literature match: Lactobacillus (increased; genus; exact); source(s): [10.1371/journal.pone.0142164](https://doi.org/10.1371/journal.pone.0142164); [10.1002/mds.26942](https://doi.org/10.1002/mds.26942)
- BugSigDB PD direction: increased (15/6 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: facultatively aerobic; moderately halophilic
- Interpretation: BugSigDB PD occurrence agrees with literature (increased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: facultatively aerobic; moderately halophilic. Strongest current influence: facultatively aerobic (+0.06 fixed-reference z, mixed).

### Streptococcus
- Literature match: Streptococcus (increased; genus; exact); source(s): [10.1007/s11427-016-9001-4](https://doi.org/10.1007/s11427-016-9001-4); [10.1038/s41467-022-34667-x](https://doi.org/10.1038/s41467-022-34667-x)
- BugSigDB PD direction: increased (11/5 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: coccus shaped; facultatively aerobic; moderately halophilic; sphere shaped
- Interpretation: BugSigDB PD occurrence agrees with literature (increased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: coccus shaped; facultatively aerobic; moderately halophilic; sphere shaped. Strongest current influence: coccus shaped (+0.10 fixed-reference z, mixed).

### Blautia
- Literature match: Blautia (decreased; genus; exact); source(s): [10.1002/mds.26307](https://doi.org/10.1002/mds.26307); [10.1007/s11427-016-9001-4](https://doi.org/10.1007/s11427-016-9001-4)
- BugSigDB PD direction: decreased (2/15 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: cellulose catabolic process; coccus shaped
- Interpretation: BugSigDB PD occurrence agrees with literature (decreased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: cellulose catabolic process; coccus shaped. Strongest current influence: coccus shaped (-0.07 fixed-reference z, mixed).

### Bifidobacterium
- Literature match: Bifidobacterium (increased; genus; exact); source(s): [10.1002/mds.26942](https://doi.org/10.1002/mds.26942); [10.1016/j.parkreldis.2016.08.019](https://doi.org/10.1016/j.parkreldis.2016.08.019)
- BugSigDB PD direction: increased (21/2 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: facultatively aerobic
- Interpretation: BugSigDB PD occurrence agrees with literature (increased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: facultatively aerobic. Strongest current influence: moderately halophilic (-0.01 fixed-reference z, mixed).

### Coprococcus
- Literature match: Coprococcus (decreased; genus; exact); source(s): [10.1002/mds.26307](https://doi.org/10.1002/mds.26307)
- BugSigDB PD direction: decreased (1/10 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: coccus shaped
- Interpretation: BugSigDB PD occurrence agrees with literature (decreased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: coccus shaped. Strongest current influence: coccus shaped (-0.01 fixed-reference z, mixed).

### Escherichia
- Literature match: Escherichia-Shigella (increased; genus group; exact); source(s): [10.1007/s11427-016-9001-4](https://doi.org/10.1007/s11427-016-9001-4)
- BugSigDB PD direction: increased (9/2 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: arsenate reduction; coccus shaped; facultatively aerobic; moderately halophilic
- Interpretation: BugSigDB PD occurrence agrees with literature (increased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: arsenate reduction; coccus shaped; facultatively aerobic; moderately halophilic. Strongest current influence: arsenate reduction (+0.02 fixed-reference z, mixed).

### Ruminococcus
- Literature match: Ruminococcus (decreased; genus; exact); source(s): [10.1007/s11427-016-9001-4](https://doi.org/10.1007/s11427-016-9001-4)
- BugSigDB PD direction: decreased (2/9 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: cellulose catabolic process; coccus shaped; sphere shaped
- Interpretation: BugSigDB PD occurrence agrees with literature (decreased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: cellulose catabolic process; coccus shaped; sphere shaped. Strongest current influence: cellulose catabolic process (-0.14 fixed-reference z, mixed).

### Enterococcus
- Literature match: Enterococcus (increased; genus; exact); source(s): [10.1007/s11427-016-9001-4](https://doi.org/10.1007/s11427-016-9001-4)
- BugSigDB PD direction: increased (7/2 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: coccus shaped; facultatively aerobic; moderately halophilic; sphere shaped
- Interpretation: BugSigDB PD occurrence agrees with literature (increased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: coccus shaped; facultatively aerobic; moderately halophilic; sphere shaped. Strongest current influence: facultatively aerobic (-0.01 fixed-reference z, mixed).

### Bacteroides
- Literature match: Bacteroides fragilis group (decreased; species group; synonym); source(s): [10.1371/journal.pone.0142164](https://doi.org/10.1371/journal.pone.0142164)
- BugSigDB PD direction: decreased (3/6 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: cellulose catabolic process; facultatively aerobic
- Interpretation: BugSigDB PD occurrence agrees with literature (decreased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: cellulose catabolic process; facultatively aerobic. Strongest current influence: cellulose catabolic process (+0.01 fixed-reference z, mixed).

### Clostridium
- Literature match: not listed (not listed; not listed; none); source(s): none
- BugSigDB PD direction: decreased (3/10 increased/decreased profiles); Hamming distance: NA (not applicable); agreement: not in literature reference
- KG-Microbe direct six-trait annotations: cellulose catabolic process
- Interpretation: No curated PD literature match was found; treat the row as BugSigDB/KG-Microbe only. KG-Microbe direct six-trait annotations: cellulose catabolic process. Strongest current influence: coccus shaped (+0.02 fixed-reference z, mixed).

### Butyrivibrio
- Literature match: not listed (not listed; not listed; none); source(s): none
- BugSigDB PD direction: decreased (2/3 increased/decreased profiles); Hamming distance: NA (not applicable); agreement: not in literature reference
- KG-Microbe direct six-trait annotations: cellulose catabolic process
- Interpretation: No curated PD literature match was found; treat the row as BugSigDB/KG-Microbe only. KG-Microbe direct six-trait annotations: cellulose catabolic process. Strongest current influence: cellulose catabolic process (-0.04 fixed-reference z, mixed).

### Butyricicoccus
- Literature match: not listed (not listed; not listed; none); source(s): none
- BugSigDB PD direction: decreased (1/6 increased/decreased profiles); Hamming distance: NA (not applicable); agreement: not in literature reference
- KG-Microbe direct six-trait annotations: coccus shaped; sphere shaped
- Interpretation: No curated PD literature match was found; treat the row as BugSigDB/KG-Microbe only. KG-Microbe direct six-trait annotations: coccus shaped; sphere shaped. Strongest current influence: sphere shaped (-0.04 fixed-reference z, mixed).

### Lachnospiraceae
- Literature match: Lachnospiraceae (decreased; family; exact); source(s): [10.1002/mds.26307](https://doi.org/10.1002/mds.26307); [10.1002/mds.26942](https://doi.org/10.1002/mds.26942)
- BugSigDB PD direction: decreased (10/39 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: cellulose catabolic process; coccus shaped; facultatively aerobic; sphere shaped
- Interpretation: BugSigDB PD occurrence agrees with literature (decreased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: cellulose catabolic process; coccus shaped; facultatively aerobic; sphere shaped. Strongest current influence: cellulose catabolic process (-0.13 fixed-reference z, dec).

### Enterobacteriaceae
- Literature match: Enterobacteriaceae (increased; family; exact); source(s): [10.1016/j.parkreldis.2016.08.019](https://doi.org/10.1016/j.parkreldis.2016.08.019); [10.1002/mds.26069](https://doi.org/10.1002/mds.26069)
- BugSigDB PD direction: increased (18/2 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: arsenate reduction; coccus shaped; facultatively aerobic; moderately halophilic; sphere shaped
- Interpretation: BugSigDB PD occurrence agrees with literature (increased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: arsenate reduction; coccus shaped; facultatively aerobic; moderately halophilic; sphere shaped. Strongest current influence: arsenate reduction (+0.20 fixed-reference z, mixed).

### Lactobacillaceae
- Literature match: Lactobacillus (increased; genus; exact); source(s): [10.1371/journal.pone.0142164](https://doi.org/10.1371/journal.pone.0142164); [10.1002/mds.26942](https://doi.org/10.1002/mds.26942)
- BugSigDB PD direction: increased (12/5 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: coccus shaped; facultatively aerobic; moderately halophilic; sphere shaped
- Interpretation: BugSigDB PD occurrence agrees with literature (increased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: coccus shaped; facultatively aerobic; moderately halophilic; sphere shaped. Strongest current influence: arsenate reduction (-0.03 fixed-reference z, mixed).

### Bifidobacteriaceae
- Literature match: Bifidobacterium (increased; genus; exact); source(s): [10.1002/mds.26942](https://doi.org/10.1002/mds.26942); [10.1016/j.parkreldis.2016.08.019](https://doi.org/10.1016/j.parkreldis.2016.08.019)
- BugSigDB PD direction: increased (16/0 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: coccus shaped; facultatively aerobic; moderately halophilic
- Interpretation: BugSigDB PD occurrence agrees with literature (increased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: coccus shaped; facultatively aerobic; moderately halophilic. Strongest current influence: coccus shaped (+0.01 fixed-reference z, mixed).

### Oscillospiraceae
- Literature match: not listed (not listed; not listed; none); source(s): none
- BugSigDB PD direction: increased (15/13 increased/decreased profiles); Hamming distance: NA (not applicable); agreement: not in literature reference
- KG-Microbe direct six-trait annotations: cellulose catabolic process; coccus shaped; facultatively aerobic; sphere shaped
- Interpretation: No curated PD literature match was found; treat the row as BugSigDB/KG-Microbe only. KG-Microbe direct six-trait annotations: cellulose catabolic process; coccus shaped; facultatively aerobic; sphere shaped. Strongest current influence: cellulose catabolic process (+0.07 fixed-reference z, mixed).

### Rikenellaceae
- Literature match: not listed (not listed; not listed; none); source(s): none
- BugSigDB PD direction: increased (7/0 increased/decreased profiles); Hamming distance: NA (not applicable); agreement: not in literature reference
- KG-Microbe direct six-trait annotations: moderately halophilic
- Interpretation: No curated PD literature match was found; treat the row as BugSigDB/KG-Microbe only. KG-Microbe direct six-trait annotations: moderately halophilic. Strongest current influence: moderately halophilic (+0.04 fixed-reference z, mixed).

### Veillonellaceae
- Literature match: not listed (not listed; not listed; none); source(s): none
- BugSigDB PD direction: balanced (3/3 increased/decreased profiles); Hamming distance: NA (not applicable); agreement: not in literature reference
- KG-Microbe direct six-trait annotations: coccus shaped; sphere shaped
- Interpretation: No curated PD literature match was found; treat the row as BugSigDB/KG-Microbe only. KG-Microbe direct six-trait annotations: coccus shaped; sphere shaped. Strongest current influence: sphere shaped (-0.02 fixed-reference z, mixed).

### Peptostreptococcaceae
- Literature match: not listed (not listed; not listed; none); source(s): none
- BugSigDB PD direction: increased (5/3 increased/decreased profiles); Hamming distance: NA (not applicable); agreement: not in literature reference
- KG-Microbe direct six-trait annotations: coccus shaped; facultatively aerobic; sphere shaped
- Interpretation: No curated PD literature match was found; treat the row as BugSigDB/KG-Microbe only. KG-Microbe direct six-trait annotations: coccus shaped; facultatively aerobic; sphere shaped. Strongest current influence: sphere shaped (-0.02 fixed-reference z, mixed).
