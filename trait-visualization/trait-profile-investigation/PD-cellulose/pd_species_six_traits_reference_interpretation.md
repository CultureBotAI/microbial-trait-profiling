# PD Taxa Reference Interpretation

This table joins the displayed species influence rows to DOI-backed PD taxa literature, BugSigDB PD occurrence direction, KG-Microbe six-trait annotations, and the strongest current influence value.

Rows follow the heatmap order: descending BugSigDB PD increased/decreased log2 profile-count ratio, with count-based and alphabetical tie-breakers.

## Agreement Summary

| Taxon | BugSigDB inc/dec | BugSigDB direction | Literature value | Match scope | Hamming distance | Agreement | Strongest six-trait influence |
| --- | ---: | --- | --- | --- | ---: | --- | --- |
| Streptococcus salivarius | 4/0 | increased | increased | Streptococcus; genus; genus_or_parent | 0 | agrees | facultatively aerobic (+0.01; mixed) |
| Methanobrevibacter smithii | 3/0 | increased | increased | Methanobrevibacter; genus; exact | 0 | agrees | coccus shaped (+0.03; mixed) |
| Escherichia coli | 1/2 | decreased | increased | Escherichia-Shigella; genus group; exact | 1 | differs | arsenate reduction (+0.10; mixed) |
| Roseburia faecis | 0/6 | decreased | decreased | Roseburia; genus; genus_or_parent | 0 | agrees | cellulose catabolic process (-0.03; mixed) |
| Faecalibacterium prausnitzii | 0/13 | decreased | decreased | Faecalibacterium; genus/species; exact | 0 | agrees | facultatively aerobic (-0.10; mixed) |

## Methods: Literature Review and Direction Comparison

The literature direction table is a targeted narrative review updated 2026-09-22; it is not a systematic review or meta-analysis. PubMed, PubMed Central, and publisher records were searched using `Parkinson disease`, `gut microbiome`, each taxon name or synonym, and direction terms (`increased`, `decreased`, `enriched`, or `depleted`). Existing DOI-backed PD microbiome references were checked backward for additional direction evidence.

Eligible evidence comprised human PD-versus-control fecal or intestinal microbiome studies and cross-cohort human meta-analyses that reported an explicit taxon direction. Animal-only, intervention-only, non-gut, and severity-only findings were excluded from increased/decreased calls. Association-only evidence is retained as `associated` and is not treated as directional.

Literature direction is assigned at the most specific supported taxonomic level. An exact species result is preferred; otherwise a genus or family result is used as a flagged parent-level proxy. `variable` is used when eligible reports conflict at the matched level. The reference table records the matched term, reference rank, match type, evidence statement, and DOI list so broader-rank comparisons remain visible.

BugSigDB direction is a recurrence summary, not a pooled abundance effect: `increased` when increased-profile count is greater than decreased-profile count, `decreased` for the reverse, and `balanced` for a tie. `agrees` and `differs` are assigned only against directional literature calls. For those calls, literature and BugSigDB direction are each encoded as one binary categorical position and compared with Hamming distance: 0 means the directions are similar and 1 means they differ. Balanced, variable, associated, and unlisted cases have no defined Hamming distance and remain separate. This is a descriptive direction check, not a statistical test of agreement, and it does not resolve cohort, medication, geography, sequencing, or taxonomic-version heterogeneity.

## Taxon Details

### Streptococcus salivarius
- Literature match: Streptococcus (increased; genus; genus_or_parent); source(s): [10.1007/s11427-016-9001-4](https://doi.org/10.1007/s11427-016-9001-4); [10.1038/s41467-022-34667-x](https://doi.org/10.1038/s41467-022-34667-x)
- BugSigDB PD direction: increased (4/0 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: coccus shaped; facultatively aerobic
- Interpretation: BugSigDB PD occurrence agrees with literature (increased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: coccus shaped; facultatively aerobic. Strongest current influence: facultatively aerobic (+0.01 fixed-reference z, mixed).

### Methanobrevibacter smithii
- Literature match: Methanobrevibacter (increased; genus; exact); source(s): [10.1038/s41467-022-34667-x](https://doi.org/10.1038/s41467-022-34667-x)
- BugSigDB PD direction: increased (3/0 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: coccus shaped
- Interpretation: BugSigDB PD occurrence agrees with literature (increased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: coccus shaped. Strongest current influence: coccus shaped (+0.03 fixed-reference z, mixed).

### Escherichia coli
- Literature match: Escherichia-Shigella (increased; genus group; exact); source(s): [10.1007/s11427-016-9001-4](https://doi.org/10.1007/s11427-016-9001-4)
- BugSigDB PD direction: decreased (1/2 increased/decreased profiles); Hamming distance: 1 (different); agreement: differs
- KG-Microbe direct six-trait annotations: arsenate reduction; coccus shaped; facultatively aerobic; moderately halophilic
- Interpretation: BugSigDB PD occurrence is decreased, but the literature reference is increased; review cohort, taxonomy, and direction labels before interpreting. KG-Microbe direct six-trait annotations: arsenate reduction; coccus shaped; facultatively aerobic; moderately halophilic. Strongest current influence: arsenate reduction (+0.10 fixed-reference z, mixed).

### Roseburia faecis
- Literature match: Roseburia (decreased; genus; genus_or_parent); source(s): [10.1002/mds.26307](https://doi.org/10.1002/mds.26307); [10.1002/mds.28119](https://doi.org/10.1002/mds.28119); [10.1038/s41467-022-34667-x](https://doi.org/10.1038/s41467-022-34667-x)
- BugSigDB PD direction: decreased (0/6 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: cellulose catabolic process
- Interpretation: BugSigDB PD occurrence agrees with literature (decreased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: cellulose catabolic process. Strongest current influence: cellulose catabolic process (-0.03 fixed-reference z, mixed).

### Faecalibacterium prausnitzii
- Literature match: Faecalibacterium (decreased; genus/species; exact); source(s): [10.1016/j.parkreldis.2016.08.019](https://doi.org/10.1016/j.parkreldis.2016.08.019); [10.1007/s11427-016-9001-4](https://doi.org/10.1007/s11427-016-9001-4); [10.1002/mds.28119](https://doi.org/10.1002/mds.28119)
- BugSigDB PD direction: decreased (0/13 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: facultatively aerobic
- Interpretation: BugSigDB PD occurrence agrees with literature (decreased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: facultatively aerobic. Strongest current influence: facultatively aerobic (-0.10 fixed-reference z, mixed).
