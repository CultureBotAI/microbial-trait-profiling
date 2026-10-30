# PD Taxa Reference Interpretation

This table joins the displayed species influence rows to DOI-backed PD taxa literature, BugSigDB PD occurrence direction, KG-Microbe six-trait annotations, and the strongest current influence value.

Rows follow the heatmap order: descending BugSigDB PD increased/decreased log2 profile-count ratio, with count-based and alphabetical tie-breakers. The arrow after the strongest-influence trait gives its observed direction in the PD six-trait lollipop: up where its mean z is higher in increased profiles, down where it is higher in decreased profiles.

## Agreement Summary

| Taxon | BugSigDB inc/dec | BugSigDB direction | Literature value | Match scope | Hamming distance | Agreement | Strongest six-trait influence |
| --- | ---: | --- | --- | --- | ---: | --- | --- |
| Methanobrevibacter smithii | 3/0 | increased | increased | Methanobrevibacter; genus; exact | 0 | agrees | coccus shaped ↑ (+0.04; mixed) |
| Roseburia faecis | 0/5 | decreased | decreased | Roseburia; genus; genus_or_parent | 0 | agrees | cellulose catabolic process ↓ (-0.04; mixed) |
| Faecalibacterium prausnitzii | 0/9 | decreased | decreased | Faecalibacterium; genus/species; exact | 0 | agrees | facultatively aerobic ↑ (-0.13; mixed) |

## Methods: Literature Review and Direction Comparison

The literature direction table is a targeted narrative review updated 2026-09-22; it is not a systematic review or meta-analysis. PubMed, PubMed Central, and publisher records were searched using `Parkinson disease`, `gut microbiome`, each taxon name or synonym, and direction terms (`increased`, `decreased`, `enriched`, or `depleted`). Existing DOI-backed PD microbiome references were checked backward for additional direction evidence.

Eligible evidence comprised human PD-versus-control fecal or intestinal microbiome studies and cross-cohort human meta-analyses that reported an explicit taxon direction. Animal-only, intervention-only, non-gut, and severity-only findings were excluded from increased/decreased calls. Association-only evidence is retained as `associated` and is not treated as directional.

Literature direction is assigned at the most specific supported taxonomic level. An exact species result is preferred; otherwise a genus or family result is used as a flagged parent-level proxy. `variable` is used when eligible reports conflict at the matched level. The reference table records the matched term, reference rank, match type, evidence statement, and DOI list so broader-rank comparisons remain visible.

BugSigDB direction is a recurrence summary, not a pooled abundance effect: `increased` when increased-profile count is greater than decreased-profile count, `decreased` for the reverse, and `balanced` for a tie. `agrees` and `differs` are assigned only against directional literature calls. For those calls, literature and BugSigDB direction are each encoded as one binary categorical position and compared with Hamming distance: 0 means the directions are similar and 1 means they differ. Balanced, variable, associated, and unlisted cases have no defined Hamming distance and remain separate. This is a descriptive direction check, not a statistical test of agreement, and it does not resolve cohort, medication, geography, sequencing, or taxonomic-version heterogeneity.

## Taxon Details

### Methanobrevibacter smithii
- Literature match: Methanobrevibacter (increased; genus; exact); source(s): [10.1038/s41467-022-34667-x](https://doi.org/10.1038/s41467-022-34667-x)
- BugSigDB PD direction: increased (3/0 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: coccus shaped
- Interpretation: BugSigDB PD occurrence agrees with literature (increased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: coccus shaped. Strongest current influence: coccus shaped (lollipop direction increased; +0.04 fixed-reference z, mixed).

### Roseburia faecis
- Literature match: Roseburia (decreased; genus; genus_or_parent); source(s): [10.1002/mds.26307](https://doi.org/10.1002/mds.26307); [10.1002/mds.28119](https://doi.org/10.1002/mds.28119); [10.1038/s41467-022-34667-x](https://doi.org/10.1038/s41467-022-34667-x)
- BugSigDB PD direction: decreased (0/5 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: cellulose catabolic process
- Interpretation: BugSigDB PD occurrence agrees with literature (decreased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: cellulose catabolic process. Strongest current influence: cellulose catabolic process (lollipop direction decreased; -0.04 fixed-reference z, mixed).

### Faecalibacterium prausnitzii
- Literature match: Faecalibacterium (decreased; genus/species; exact); source(s): [10.1016/j.parkreldis.2016.08.019](https://doi.org/10.1016/j.parkreldis.2016.08.019); [10.1007/s11427-016-9001-4](https://doi.org/10.1007/s11427-016-9001-4); [10.1002/mds.28119](https://doi.org/10.1002/mds.28119)
- BugSigDB PD direction: decreased (0/9 increased/decreased profiles); Hamming distance: 0 (similar); agreement: agrees
- KG-Microbe direct six-trait annotations: facultatively aerobic
- Interpretation: BugSigDB PD occurrence agrees with literature (decreased); treat this as a higher-priority taxon for trait-level follow-up. KG-Microbe direct six-trait annotations: facultatively aerobic. Strongest current influence: facultatively aerobic (lollipop direction increased; -0.13 fixed-reference z, mixed).
