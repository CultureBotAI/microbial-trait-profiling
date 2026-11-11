# ORA background scope (#746)

The default `--background-scope all_matrix` retains the established urn:
M is every row of `taxa_trait_matrix.tsv`; K is the number of those rows positive
for the tested trait. The matrix contains BugSigDB-reported taxa mapped across hosts
and sequencing types, including mixed taxonomic ranks and inherited annotations.
It is not the universe of all organisms in KG-Microbe. The default keeps existing
report and TSV bytes unchanged, including `TaxaMatrix (filename)` metadata.

Profile eligibility is a different choice. Human 16S/WMS experiments define the
analysis profiles; that filter does not automatically filter the matrix background.
For an explicitly matched reference, use `--background-scope profile_corpus` with
`--corpus-taxa-file`: a newline list of IDs (optionally headed `taxon_id`) or a TSV
with a `taxon_id` column. Supply the distinct matrix taxon IDs represented in the
eligible profile corpus, not organism names or all KG nodes. This script does not
infer host or sequencing eligibility from the file's name.

In corpus mode, only matching matrix rows contribute to both M and K. Duplicate
IDs in the corpus list are collapsed; an absent, empty, or out-of-matrix taxon
list fails closed. Existing checks still reject profile counts incompatible with
the selected urn. The aggregated profile input cannot prove membership of each
reported taxon, so preparing a compatible corpus list remains the caller's duty.

Changing M or K changes the reference for every hypergeometric p-value and can
change the K>0 testing family and BH-adjusted q-values. Numerical values need not
all differ (for example, identical reference sets give identical results). Declare
the scope explicitly when comparing runs; never interpret a changed background as
a new physiological observation. Profile draws N and observed trait counts k are
unchanged. These are organism-level trait annotations, not in vivo activity.

Corpus-mode reports identify `background_scope=profile_corpus` and
`ProfileCorpusTaxaMatrix (filename; corpus=filename)`. Both TSV exports additionally
record `Background-Scope` and `Corpus-Taxa-SHA256`, alongside the existing full
matrix/profile fingerprints and M/K counts. In legacy default outputs the unchanged
`Background-Mode: TaxaMatrix (...)` identifies all_matrix; no new metadata lines
are inserted, preserving committed output compatibility.
