# Offline raw-to-core fixture for issue #629

These synthetic CSV/TSV inputs exercise the production builders without downloads.
The increased profile `bsdb:1 / bsdb:1-1 / inc` reports genus G twice and its
species A once. G has three distinct descendant taxa annotated with
`METPO:1000608`, including A. Duplicate source edges for A do not add another
taxon count. Deduplicating reported G,G,A to G,A gives a count-profile value of
4, an ORA-profile value of 2, and an ORA denominator N of 2. A contributes to
both G and A, so these are overlapping annotations rather than four organisms.

Other profiles cover a culture-medium-only annotation, a present-in-KG taxon
with no retained annotations, a taxon absent from the KG, an incomplete
increased/decreased pair, a family-level report, and a mouse-host experiment.
The raw core retains MediaDive columns and incomplete pairs. Count profiles
include all hosts by default; ORA profiles include humans by default.

The matrix builder emits all reported taxa as rows. Consequently, the KG-absent
taxon remains matrix-matched, has zero annotation counts, and contributes to
ORA N. A separate supplied-matrix test removes this row to test actual matrix
lookup failure. Matrix matching and KG membership are different coverage facts.
