# PD Taxa Reference Interpretation

This table joins the displayed species influence rows to DOI-backed PD taxa literature, BugSigDB PD occurrence direction, KG-Microbe six-trait annotations, and the strongest current influence value.

Rows follow the heatmap order: descending BugSigDB PD increased/decreased log2 profile-count ratio, with count-based and alphabetical tie-breakers. The arrow after the strongest-influence trait gives its observed direction in the PD six-trait lollipop: up where its mean z is higher in increased profiles, down where it is higher in decreased profiles.

## Agreement Summary

| Taxon | BugSigDB inc/dec | BugSigDB direction | Literature value | Match scope | Hamming distance | Agreement | Strongest six-trait influence |
| --- | ---: | --- | --- | --- | ---: | --- | --- |
| Faecalibacterium prausnitzii | 0/9 | decreased | decreased | Faecalibacterium; genus/species; exact | 0 | agrees | facultatively aerobic ↑ (-0.13; mixed) |
| Roseburia faecis | 0/5 | decreased | decreased | Roseburia; genus; genus_or_parent | 0 | agrees | cellulose catabolic process ↓ (-0.04; mixed) |
| Methanobrevibacter smithii | 3/0 | increased | increased | Methanobrevibacter; genus; exact | 0 | agrees | coccus shaped ↑ (+0.04; mixed) |
| Blautia obeum | 5/0 | increased | decreased | Blautia; genus; genus_or_parent | 1 | differs | cellulose catabolic process ↓ (+0.15; mixed) |
| Escherichia coli | 1/0 | increased | increased | Escherichia-Shigella; genus group; exact | 0 | agrees | arsenate reduction ↑ (+0.14; low support) |
| Streptococcus mutans | 1/0 | increased | increased | Streptococcus; genus; genus_or_parent | 0 | agrees | coccus shaped ↑ (+0.04; low support) |
| Streptococcus salivarius | 1/0 | increased | increased | Streptococcus; genus; genus_or_parent | 0 | agrees | facultatively aerobic ↑ (+0.01; low support) |
| Roseburia | 1/29 | decreased | decreased | Roseburia; genus; exact | 0 | agrees | cellulose catabolic process ↓ (-0.16; mixed) |
| Faecalibacterium | 0/24 | decreased | decreased | Faecalibacterium; genus/species; exact | 0 | agrees | facultatively aerobic ↑ (-0.13; dec) |
| Blautia | 1/12 | decreased | decreased | Blautia; genus; exact | 0 | agrees | sphere shaped ↓ (+0.03; mixed) |
| Bifidobacterium | 17/0 | increased | increased | Bifidobacterium; genus; exact | 0 | agrees | moderately halophilic ↑ (-0.02; mixed) |
| Lactobacillus | 9/4 | increased | increased | Lactobacillus; genus; exact | 0 | agrees | sphere shaped ↓ (-0.01; dec) |
| Streptococcus | 7/3 | increased | increased | Streptococcus; genus; exact | 0 | agrees | coccus shaped ↑ (+0.13; mixed) |
| Bacteroides | 3/6 | decreased | decreased | Bacteroides fragilis group; species group; synonym | 0 | agrees | cellulose catabolic process ↓ (+0.02; mixed) |
| Ruminococcus | 1/5 | decreased | decreased | Ruminococcus; genus; exact | 0 | agrees | sphere shaped ↓ (-0.10; mixed) |
| Alistipes | 5/0 | increased | increased | Alistipes; genus/species; exact | 0 | agrees | moderately halophilic ↑ (+0.01; mixed) |
| Enterococcus | 4/0 | increased | increased | Enterococcus; genus; exact | 0 | agrees | moderately halophilic ↑ (+0.05; mixed) |
| Clostridium | 2/9 | decreased | not listed | not listed; not listed; none | NA | not in literature reference | coccus shaped ↑ (+0.03; mixed) |
| Veillonella | 4/1 | increased | not listed | not listed; not listed; none | NA | not in literature reference | sphere shaped ↓ (-0.01; mixed) |
| Lachnospiraceae | 6/30 | decreased | decreased | Lachnospiraceae; family; exact | 0 | agrees | cellulose catabolic process ↓ (-0.11; mixed) |
| Bifidobacteriaceae | 12/0 | increased | increased | Bifidobacterium; genus; exact | 0 | agrees | coccus shaped ↑ (+0.03; mixed) |
| Prevotellaceae | 0/8 | decreased | decreased | Prevotella; genus/family; exact | 0 | agrees | facultatively aerobic ↑ (+0.02; mixed) |
| Enterobacteriaceae | 11/0 | increased | increased | Enterobacteriaceae; family; exact | 0 | agrees | arsenate reduction ↑ (+0.16; inc) |
| Lactobacillaceae | 7/3 | increased | increased | Lactobacillus; genus; exact | 0 | agrees | arsenate reduction ↑ (-0.03; mixed) |
| Oscillospiraceae | 10/7 | increased | not listed | not listed; not listed; none | NA | not in literature reference | cellulose catabolic process ↓ (+0.08; mixed) |
| Clostridiaceae | 1/8 | decreased | not listed | not listed; not listed; none | NA | not in literature reference | sphere shaped ↓ (+0.01; mixed) |
| Veillonellaceae | 3/3 | balanced | not listed | not listed; not listed; none | NA | not in literature reference | sphere shaped ↓ (-0.05; mixed) |
| Rikenellaceae | 4/0 | increased | not listed | not listed; not listed; none | NA | not in literature reference | moderately halophilic ↑ (+0.03; mixed) |
| Pasteurellaceae | 0/4 | decreased | not listed | not listed; not listed; none | NA | not in literature reference | facultatively aerobic ↑ (-0.08; mixed) |
| Peptostreptococcaceae | 2/2 | balanced | not listed | not listed; not listed; none | NA | not in literature reference | sphere shaped ↓ (-0.03; mixed) |

## Methods: Literature Review and Direction Comparison

The literature direction table is a targeted narrative review updated 2026-09-22; it is not a systematic review or meta-analysis. PubMed, PubMed Central, and publisher records were searched using `Parkinson disease`, `gut microbiome`, each taxon name or synonym, and direction terms (`increased`, `decreased`, `enriched`, or `depleted`). Existing DOI-backed PD microbiome references were checked backward for additional direction evidence.

Eligible evidence comprised human PD-versus-control fecal or intestinal microbiome studies and cross-cohort human meta-analyses that reported an explicit taxon direction. Animal-only, intervention-only, non-gut, and severity-only findings were excluded from increased/decreased calls. Association-only evidence is retained as `associated` and is not treated as directional.

Literature direction is assigned at the most specific supported taxonomic level. An exact species result is preferred; otherwise a genus or family result is used as a flagged parent-level proxy. `variable` is used when eligible reports conflict at the matched level. The reference table records the matched term, reference rank, match type, evidence statement, and DOI list so broader-rank comparisons remain visible.

BugSigDB direction is a recurrence summary, not a pooled abundance effect: `increased` when increased-profile count is greater than decreased-profile count, `decreased` for the reverse, and `balanced` for a tie. `agrees` and `differs` are assigned only against directional literature calls. For those calls, literature and BugSigDB direction are each encoded as one binary categorical position and compared with Hamming distance: 0 means the directions are similar and 1 means they differ. Balanced, variable, associated, and unlisted cases have no defined Hamming distance and remain separate. This is a descriptive direction check, not a statistical test of agreement, and it does not resolve cohort, medication, geography, sequencing, or taxonomic-version heterogeneity.

## Taxon Details

### Faecalibacterium prausnitzii
- Literature match: Faecalibacterium (decreased; genus/species; exact); source(s): [10.1016/j.parkreldis.2016.08.019](https://doi.org/10.1016/j.parkreldis.2016.08.019); [10.1007/s11427-016-9001-4](https://doi.org/10.1007/s11427-016-9001-4); [10.1002/mds.28119](https://doi.org/10.1002/mds.28119)
- BugSigDB PD direction: decreased (0/9 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: facultatively aerobic
- Interpretation: BugSigDB PD occurrence agrees with literature (decreased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: facultatively aerobic. Strongest current influence: facultatively aerobic (lollipop direction increased; -0.13 fixed-reference z, mixed).

### Roseburia faecis
- Literature match: Roseburia (decreased; genus; genus_or_parent); source(s): [10.1002/mds.26307](https://doi.org/10.1002/mds.26307); [10.1002/mds.28119](https://doi.org/10.1002/mds.28119); [10.1038/s41467-022-34667-x](https://doi.org/10.1038/s41467-022-34667-x)
- BugSigDB PD direction: decreased (0/5 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: cellulose catabolic process
- Interpretation: BugSigDB PD occurrence agrees with literature (decreased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: cellulose catabolic process. Strongest current influence: cellulose catabolic process (lollipop direction decreased; -0.04 fixed-reference z, mixed).

### Methanobrevibacter smithii
- Literature match: Methanobrevibacter (increased; genus; exact); source(s): [10.1038/s41467-022-34667-x](https://doi.org/10.1038/s41467-022-34667-x)
- BugSigDB PD direction: increased (3/0 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: coccus shaped
- Interpretation: BugSigDB PD occurrence agrees with literature (increased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: coccus shaped. Strongest current influence: coccus shaped (lollipop direction increased; +0.04 fixed-reference z, mixed).

### Blautia obeum
- Literature match: Blautia (decreased; genus; genus_or_parent); source(s): [10.1002/mds.26307](https://doi.org/10.1002/mds.26307); [10.1007/s11427-016-9001-4](https://doi.org/10.1007/s11427-016-9001-4)
- BugSigDB PD direction: increased (5/0 increased/decreased profiles); Hamming distance: 1 (different); agreement: differs
- KG-Microbe direct six-trait annotations: cellulose catabolic process; coccus shaped
- Interpretation: BugSigDB PD occurrence is increased, but the literature reference is decreased; review cohort, taxonomy, and direction labels before interpreting. KG-Microbe direct six-trait annotations: cellulose catabolic process; coccus shaped. Strongest current influence: cellulose catabolic process (lollipop direction decreased; +0.15 fixed-reference z, mixed).

### Escherichia coli
- Literature match: Escherichia-Shigella (increased; genus group; exact); source(s): [10.1007/s11427-016-9001-4](https://doi.org/10.1007/s11427-016-9001-4)
- BugSigDB PD direction: increased (1/0 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: arsenate reduction; coccus shaped; facultatively aerobic; moderately halophilic
- Interpretation: BugSigDB PD occurrence agrees with literature (increased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: arsenate reduction; coccus shaped; facultatively aerobic; moderately halophilic. Strongest current influence: arsenate reduction (lollipop direction increased; +0.14 fixed-reference z, low support).

### Streptococcus mutans
- Literature match: Streptococcus (increased; genus; genus_or_parent); source(s): [10.1007/s11427-016-9001-4](https://doi.org/10.1007/s11427-016-9001-4); [10.1038/s41467-022-34667-x](https://doi.org/10.1038/s41467-022-34667-x)
- BugSigDB PD direction: increased (1/0 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: coccus shaped; facultatively aerobic; moderately halophilic
- Interpretation: BugSigDB PD occurrence agrees with literature (increased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: coccus shaped; facultatively aerobic; moderately halophilic. Strongest current influence: coccus shaped (lollipop direction increased; +0.04 fixed-reference z, low support).

### Streptococcus salivarius
- Literature match: Streptococcus (increased; genus; genus_or_parent); source(s): [10.1007/s11427-016-9001-4](https://doi.org/10.1007/s11427-016-9001-4); [10.1038/s41467-022-34667-x](https://doi.org/10.1038/s41467-022-34667-x)
- BugSigDB PD direction: increased (1/0 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: coccus shaped; facultatively aerobic
- Interpretation: BugSigDB PD occurrence agrees with literature (increased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: coccus shaped; facultatively aerobic. Strongest current influence: facultatively aerobic (lollipop direction increased; +0.01 fixed-reference z, low support).

### Roseburia
- Literature match: Roseburia (decreased; genus; exact); source(s): [10.1002/mds.26307](https://doi.org/10.1002/mds.26307); [10.1002/mds.28119](https://doi.org/10.1002/mds.28119); [10.1038/s41467-022-34667-x](https://doi.org/10.1038/s41467-022-34667-x)
- BugSigDB PD direction: decreased (1/29 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: cellulose catabolic process
- Interpretation: BugSigDB PD occurrence agrees with literature (decreased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: cellulose catabolic process. Strongest current influence: cellulose catabolic process (lollipop direction decreased; -0.16 fixed-reference z, mixed).

### Faecalibacterium
- Literature match: Faecalibacterium (decreased; genus/species; exact); source(s): [10.1016/j.parkreldis.2016.08.019](https://doi.org/10.1016/j.parkreldis.2016.08.019); [10.1007/s11427-016-9001-4](https://doi.org/10.1007/s11427-016-9001-4); [10.1002/mds.28119](https://doi.org/10.1002/mds.28119)
- BugSigDB PD direction: decreased (0/24 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: facultatively aerobic
- Interpretation: BugSigDB PD occurrence agrees with literature (decreased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: facultatively aerobic. Strongest current influence: facultatively aerobic (lollipop direction increased; -0.13 fixed-reference z, dec).

### Blautia
- Literature match: Blautia (decreased; genus; exact); source(s): [10.1002/mds.26307](https://doi.org/10.1002/mds.26307); [10.1007/s11427-016-9001-4](https://doi.org/10.1007/s11427-016-9001-4)
- BugSigDB PD direction: decreased (1/12 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: cellulose catabolic process; coccus shaped
- Interpretation: BugSigDB PD occurrence agrees with literature (decreased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: cellulose catabolic process; coccus shaped. Strongest current influence: sphere shaped (lollipop direction decreased; +0.03 fixed-reference z, mixed).

### Bifidobacterium
- Literature match: Bifidobacterium (increased; genus; exact); source(s): [10.1002/mds.26942](https://doi.org/10.1002/mds.26942); [10.1016/j.parkreldis.2016.08.019](https://doi.org/10.1016/j.parkreldis.2016.08.019)
- BugSigDB PD direction: increased (17/0 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: facultatively aerobic
- Interpretation: BugSigDB PD occurrence agrees with literature (increased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: facultatively aerobic. Strongest current influence: moderately halophilic (lollipop direction increased; -0.02 fixed-reference z, mixed).

### Lactobacillus
- Literature match: Lactobacillus (increased; genus; exact); source(s): [10.1371/journal.pone.0142164](https://doi.org/10.1371/journal.pone.0142164); [10.1002/mds.26942](https://doi.org/10.1002/mds.26942)
- BugSigDB PD direction: increased (9/4 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: facultatively aerobic; moderately halophilic
- Interpretation: BugSigDB PD occurrence agrees with literature (increased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: facultatively aerobic; moderately halophilic. Strongest current influence: sphere shaped (lollipop direction decreased; -0.01 fixed-reference z, dec).

### Streptococcus
- Literature match: Streptococcus (increased; genus; exact); source(s): [10.1007/s11427-016-9001-4](https://doi.org/10.1007/s11427-016-9001-4); [10.1038/s41467-022-34667-x](https://doi.org/10.1038/s41467-022-34667-x)
- BugSigDB PD direction: increased (7/3 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: coccus shaped; facultatively aerobic; moderately halophilic; sphere shaped
- Interpretation: BugSigDB PD occurrence agrees with literature (increased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: coccus shaped; facultatively aerobic; moderately halophilic; sphere shaped. Strongest current influence: coccus shaped (lollipop direction increased; +0.13 fixed-reference z, mixed).

### Bacteroides
- Literature match: Bacteroides fragilis group (decreased; species group; synonym); source(s): [10.1371/journal.pone.0142164](https://doi.org/10.1371/journal.pone.0142164)
- BugSigDB PD direction: decreased (3/6 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: cellulose catabolic process; facultatively aerobic
- Interpretation: BugSigDB PD occurrence agrees with literature (decreased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: cellulose catabolic process; facultatively aerobic. Strongest current influence: cellulose catabolic process (lollipop direction decreased; +0.02 fixed-reference z, mixed).

### Ruminococcus
- Literature match: Ruminococcus (decreased; genus; exact); source(s): [10.1007/s11427-016-9001-4](https://doi.org/10.1007/s11427-016-9001-4)
- BugSigDB PD direction: decreased (1/5 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: cellulose catabolic process; coccus shaped; sphere shaped
- Interpretation: BugSigDB PD occurrence agrees with literature (decreased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: cellulose catabolic process; coccus shaped; sphere shaped. Strongest current influence: sphere shaped (lollipop direction decreased; -0.10 fixed-reference z, mixed).

### Alistipes
- Literature match: Alistipes (increased; genus/species; exact); source(s): [10.1093/brain/awaa201](https://doi.org/10.1093/brain/awaa201); [10.1038/s41467-022-34667-x](https://doi.org/10.1038/s41467-022-34667-x)
- BugSigDB PD direction: increased (5/0 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: moderately halophilic
- Interpretation: BugSigDB PD occurrence agrees with literature (increased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: moderately halophilic. Strongest current influence: moderately halophilic (lollipop direction increased; +0.01 fixed-reference z, mixed).

### Enterococcus
- Literature match: Enterococcus (increased; genus; exact); source(s): [10.1007/s11427-016-9001-4](https://doi.org/10.1007/s11427-016-9001-4)
- BugSigDB PD direction: increased (4/0 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: coccus shaped; facultatively aerobic; moderately halophilic; sphere shaped
- Interpretation: BugSigDB PD occurrence agrees with literature (increased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: coccus shaped; facultatively aerobic; moderately halophilic; sphere shaped. Strongest current influence: moderately halophilic (lollipop direction increased; +0.05 fixed-reference z, mixed).

### Clostridium
- Literature match: not listed (not listed; not listed; none); source(s): none
- BugSigDB PD direction: decreased (2/9 increased/decreased profiles); Hamming distance: NA (not applicable); agreement: not in literature reference
- KG-Microbe direct six-trait annotations: cellulose catabolic process
- Interpretation: No curated PD literature match was found; treat the row as BugSigDB/KG-Microbe only. KG-Microbe direct six-trait annotations: cellulose catabolic process. Strongest current influence: coccus shaped (lollipop direction increased; +0.03 fixed-reference z, mixed).

### Veillonella
- Literature match: not listed (not listed; not listed; none); source(s): none
- BugSigDB PD direction: increased (4/1 increased/decreased profiles); Hamming distance: NA (not applicable); agreement: not in literature reference
- KG-Microbe direct six-trait annotations: coccus shaped; sphere shaped
- Interpretation: No curated PD literature match was found; treat the row as BugSigDB/KG-Microbe only. KG-Microbe direct six-trait annotations: coccus shaped; sphere shaped. Strongest current influence: sphere shaped (lollipop direction decreased; -0.01 fixed-reference z, mixed).

### Lachnospiraceae
- Literature match: Lachnospiraceae (decreased; family; exact); source(s): [10.1002/mds.26307](https://doi.org/10.1002/mds.26307); [10.1002/mds.26942](https://doi.org/10.1002/mds.26942)
- BugSigDB PD direction: decreased (6/30 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: cellulose catabolic process; coccus shaped; facultatively aerobic; sphere shaped
- Interpretation: BugSigDB PD occurrence agrees with literature (decreased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: cellulose catabolic process; coccus shaped; facultatively aerobic; sphere shaped. Strongest current influence: cellulose catabolic process (lollipop direction decreased; -0.11 fixed-reference z, mixed).

### Bifidobacteriaceae
- Literature match: Bifidobacterium (increased; genus; exact); source(s): [10.1002/mds.26942](https://doi.org/10.1002/mds.26942); [10.1016/j.parkreldis.2016.08.019](https://doi.org/10.1016/j.parkreldis.2016.08.019)
- BugSigDB PD direction: increased (12/0 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: coccus shaped; facultatively aerobic; moderately halophilic
- Interpretation: BugSigDB PD occurrence agrees with literature (increased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: coccus shaped; facultatively aerobic; moderately halophilic. Strongest current influence: coccus shaped (lollipop direction increased; +0.03 fixed-reference z, mixed).

### Prevotellaceae
- Literature match: Prevotella (decreased; genus/family; exact); source(s): [10.1002/mds.26069](https://doi.org/10.1002/mds.26069); [10.1186/s13073-017-0428-y](https://doi.org/10.1186/s13073-017-0428-y)
- BugSigDB PD direction: decreased (0/8 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: coccus shaped; sphere shaped
- Interpretation: BugSigDB PD occurrence agrees with literature (decreased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: coccus shaped; sphere shaped. Strongest current influence: facultatively aerobic (lollipop direction increased; +0.02 fixed-reference z, mixed).

### Enterobacteriaceae
- Literature match: Enterobacteriaceae (increased; family; exact); source(s): [10.1016/j.parkreldis.2016.08.019](https://doi.org/10.1016/j.parkreldis.2016.08.019); [10.1002/mds.26069](https://doi.org/10.1002/mds.26069)
- BugSigDB PD direction: increased (11/0 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: arsenate reduction; coccus shaped; facultatively aerobic; moderately halophilic; sphere shaped
- Interpretation: BugSigDB PD occurrence agrees with literature (increased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: arsenate reduction; coccus shaped; facultatively aerobic; moderately halophilic; sphere shaped. Strongest current influence: arsenate reduction (lollipop direction increased; +0.16 fixed-reference z, inc).

### Lactobacillaceae
- Literature match: Lactobacillus (increased; genus; exact); source(s): [10.1371/journal.pone.0142164](https://doi.org/10.1371/journal.pone.0142164); [10.1002/mds.26942](https://doi.org/10.1002/mds.26942)
- BugSigDB PD direction: increased (7/3 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: coccus shaped; facultatively aerobic; moderately halophilic; sphere shaped
- Interpretation: BugSigDB PD occurrence agrees with literature (increased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: coccus shaped; facultatively aerobic; moderately halophilic; sphere shaped. Strongest current influence: arsenate reduction (lollipop direction increased; -0.03 fixed-reference z, mixed).

### Oscillospiraceae
- Literature match: not listed (not listed; not listed; none); source(s): none
- BugSigDB PD direction: increased (10/7 increased/decreased profiles); Hamming distance: NA (not applicable); agreement: not in literature reference
- KG-Microbe direct six-trait annotations: cellulose catabolic process; coccus shaped; facultatively aerobic; sphere shaped
- Interpretation: No curated PD literature match was found; treat the row as BugSigDB/KG-Microbe only. KG-Microbe direct six-trait annotations: cellulose catabolic process; coccus shaped; facultatively aerobic; sphere shaped. Strongest current influence: cellulose catabolic process (lollipop direction decreased; +0.08 fixed-reference z, mixed).

### Clostridiaceae
- Literature match: not listed (not listed; not listed; none); source(s): none
- BugSigDB PD direction: decreased (1/8 increased/decreased profiles); Hamming distance: NA (not applicable); agreement: not in literature reference
- KG-Microbe direct six-trait annotations: cellulose catabolic process
- Interpretation: No curated PD literature match was found; treat the row as BugSigDB/KG-Microbe only. KG-Microbe direct six-trait annotations: cellulose catabolic process. Strongest current influence: sphere shaped (lollipop direction decreased; +0.01 fixed-reference z, mixed).

### Veillonellaceae
- Literature match: not listed (not listed; not listed; none); source(s): none
- BugSigDB PD direction: balanced (3/3 increased/decreased profiles); Hamming distance: NA (not applicable); agreement: not in literature reference
- KG-Microbe direct six-trait annotations: coccus shaped; sphere shaped
- Interpretation: No curated PD literature match was found; treat the row as BugSigDB/KG-Microbe only. KG-Microbe direct six-trait annotations: coccus shaped; sphere shaped. Strongest current influence: sphere shaped (lollipop direction decreased; -0.05 fixed-reference z, mixed).

### Rikenellaceae
- Literature match: not listed (not listed; not listed; none); source(s): none
- BugSigDB PD direction: increased (4/0 increased/decreased profiles); Hamming distance: NA (not applicable); agreement: not in literature reference
- KG-Microbe direct six-trait annotations: moderately halophilic
- Interpretation: No curated PD literature match was found; treat the row as BugSigDB/KG-Microbe only. KG-Microbe direct six-trait annotations: moderately halophilic. Strongest current influence: moderately halophilic (lollipop direction increased; +0.03 fixed-reference z, mixed).

### Pasteurellaceae
- Literature match: not listed (not listed; not listed; none); source(s): none
- BugSigDB PD direction: decreased (0/4 increased/decreased profiles); Hamming distance: NA (not applicable); agreement: not in literature reference
- KG-Microbe direct six-trait annotations: coccus shaped; facultatively aerobic; moderately halophilic
- Interpretation: No curated PD literature match was found; treat the row as BugSigDB/KG-Microbe only. KG-Microbe direct six-trait annotations: coccus shaped; facultatively aerobic; moderately halophilic. Strongest current influence: facultatively aerobic (lollipop direction increased; -0.08 fixed-reference z, mixed).

### Peptostreptococcaceae
- Literature match: not listed (not listed; not listed; none); source(s): none
- BugSigDB PD direction: balanced (2/2 increased/decreased profiles); Hamming distance: NA (not applicable); agreement: not in literature reference
- KG-Microbe direct six-trait annotations: coccus shaped; facultatively aerobic; sphere shaped
- Interpretation: No curated PD literature match was found; treat the row as BugSigDB/KG-Microbe only. KG-Microbe direct six-trait annotations: coccus shaped; facultatively aerobic; sphere shaped. Strongest current influence: sphere shaped (lollipop direction decreased; -0.03 fixed-reference z, mixed).
