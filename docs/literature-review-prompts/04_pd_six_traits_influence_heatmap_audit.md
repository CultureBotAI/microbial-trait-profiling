# Prompt 4 — PD exploratory taxa × six traits influence heatmap: evidence audit

Shared by Jed Dongjin Kim-Ozaeta on 2026-10-14. Saved verbatim; see [README](README.md) for status.

Figure: [`pd_exploratory_taxa_six_traits_influence_heatmap.png`](https://github.com/CultureBotAI/microbial-trait-profiling/blob/205257c5ec887bdf145e2c66e12f17ee005fafa0/trait-visualization/trait-profile-investigation/PD-cellulose/pd_exploratory_taxa_six_traits_influence_heatmap.png), on branch `clustered-heatmap-generation` at `205257c`.

---

Conduct a deep, systematic literature review of the data represented in:

trait-visualization/trait-profile-investigation/PD-cellulose/pd_exploratory_taxa_six_traits_influence_heatmap.png

Also inspect any source tables, metadata, README files, or citation files in the same directory that were used to generate the heatmap.

## Objective

Evaluate every taxon–trait relationship shown in the heatmap in the context of Parkinson's disease (PD):

1. Independently verify relationships that already have references.
2. For relationships marked "no reference," "no evidence," missing, or otherwise unsupported, conduct a comprehensive search for relevant published evidence.
3. Clearly distinguish "no literature found" from evidence showing that a relationship is absent or contradicted.

First, determine exactly what each heatmap cell represents. Do not assume whether it denotes:

- a taxon possessing a microbial trait;
- a trait affecting the taxon's abundance or activity;
- a taxon or trait influencing PD;
- an association observed in PD cohorts; or
- a proposed mechanistic relationship.

Use the figure legend and accompanying files to resolve this before assessing the evidence.

## Literature-search requirements

- Search PubMed/MEDLINE, Europe PMC, Google Scholar, Crossref, and relevant microbiology resources.
- Search each accepted taxon name together with historical names, synonyms, reclassifications, strain names, and relevant higher or lower taxonomic levels.
- Combine those terms with each of the six trait names and their synonyms.
- Where relevant, also combine them with terms such as:
  "Parkinson's disease," "parkinsonism," "gut microbiome," "cellulose," "cellulolytic," "cellulose degradation," "fermentation," "short-chain fatty acids," "intestinal barrier," "neuroinflammation," and "alpha-synuclein."
- Search primary experimental studies first. Use authoritative genome, culture, or taxonomic resources when primary phenotypic literature is unavailable.
- Use reviews only to locate primary studies; do not treat a review as the sole support when the underlying study can be identified.
- Search literature available through the current date.

## Evidence-assessment rules

For every heatmap relationship:

- Verify that the cited paper exists and that its DOI, PMID, title, authors, and year are correct.
- Confirm that the paper actually supports the claim attributed to it.
- Record whether the evidence is direct or indirect.
- Record the studied organism's exact taxonomic level and strain.
- Do not automatically generalize strain-level findings to an entire species, genus, or family.
- Treat evidence from a taxonomic parent, relative, metagenome-assembled genome, predicted pathway, or computational annotation as indirect.
- Distinguish demonstrated phenotype from genomic potential, inferred function, ecological correlation, and mechanistic speculation.
- For PD evidence, distinguish human observational studies, intervention studies, animal models, in vitro experiments, and general biological plausibility.
- Note the relevant body site, cohort, experimental conditions, substrate, and direction of association.
- Identify conflicting findings and plausible reasons for them.
- Do not invent citations or infer support solely from article titles or abstracts.
- Label inaccessible or insufficiently described evidence as "unable to verify."
- If no relevant literature is found after the full search, report "no evidence located," document the searches performed, and do not interpret this as evidence of absence.

## Required output

### 1. Evidence-audit table

Create one row per taxon–trait relationship with these columns:

- taxon shown in heatmap
- current accepted taxon name
- historical names/synonyms searched
- taxonomic rank
- trait
- precise claim encoded by the heatmap
- heatmap's current evidence status
- existing reference(s)
- existing reference verified? (yes/partial/no/unable)
- newly identified reference(s)
- evidence type
- directness (direct/indirect)
- evidence strength (strong/moderate/weak/none/contradictory)
- organism or strain actually studied
- experimental or clinical context
- PD relevance (direct/indirect/not established)
- conclusion
- caveats
- DOI
- PMID
- stable URL

### 2. Citation-verification section

For every existing reference, state whether it:

- directly supports the heatmap claim;
- supports only a narrower or broader claim;
- is miscited or taxonomically outdated;
- contradicts the claim; or
- cannot be verified.

### 3. Previously unsupported relationships

Provide a dedicated table for cells originally marked as having no reference. Classify each as:

- new direct evidence found;
- new indirect evidence found;
- conflicting evidence found; or
- no evidence located after systematic searching.

Include the exact search strings and databases used for relationships where no evidence was located.

### 4. Narrative synthesis

Summarize:

- which taxon–trait relationships are best supported;
- which require correction or qualification;
- which remain unsupported;
- whether the six traits provide biologically plausible mechanisms connecting the taxa to PD;
- which conclusions are associations versus experimentally demonstrated mechanisms; and
- the most important knowledge gaps and recommended follow-up studies.

### 5. Proposed heatmap revisions

Provide a concise list of recommended changes to the underlying data and heatmap, including:

- references to add or replace;
- taxonomic names to update;
- claims that should be narrowed;
- cells that should be labeled "indirect evidence";
- cells that remain "no evidence located"; and
- any relationships contradicted by the literature.

Return both a readable Markdown report and a machine-readable TSV/CSV evidence table. Use inline citations and provide complete references with DOI, PMID, and stable links. Preserve uncertainty and avoid overstating causal relevance to Parkinson's disease.
