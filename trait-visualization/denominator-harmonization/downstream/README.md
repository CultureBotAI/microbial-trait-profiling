# Issue #174: downstream replay

The baseline and harmonized runs use the exact matrix cited in #174, the
same auxiliary files, and the original `a2b2df9` renderer defaults. Outputs
are under `trait-visualization/trait-profile-investigation/`; baseline TSVs
are retained in `before/`. Baseline images are reproducible in
`.tmp/issue174/downstream-before/` and are not duplicated in Git.

## Main consequences

- The six pooled q <= .001 traits and their directions remain unchanged.
  This does not imply unchanged profile-level z-scores or downstream edges.
- The combined high-support Sankey still displays 144 trait/direction edges,
  72 per direction. Six edges are replaced: the increased-direction
  `Respiratory / airway disease | Feces` pair leaves and
  `Colorectal neoplasia | Lower GI tract` enters, across all six traits.
  Among 138 retained edges, 42 change positive-z focus-profile support.
  The original top-12 selection per direction and support rules are unchanged.
- All six view-specific BH-selected heatmap sets are identical. Pattern
  heatmaps retain 15/6/4/5 displayed cells for increased grouped, decreased
  grouped, increased individual, and decreased individual views respectively.
  Signed-direction views show 413 -> 414 grouped and 597 -> 599 individual
  supported cells. A filled signed-direction cell is not a significance claim.
- The PD cellulose raw taxon-driver table (60 taxa) and taxon/trait-presence
  table (108 rows) are unchanged. Positive-z profile counts do change.
- Condition modules retain the same 80 modules and 60 selected traits.
  All 60 selection variances/scores change; transformed heatmap values change.
- The follow-up retains the same 45 x 34 binary organism/trait heatmap.
  All four reported guild-fraction summaries change. Ten of 481 condition
  depth-outlier flags change, and five of ten profile-neighbor memberships
  are replaced. Depth and cosine changes reflect removal of MediaDive
  columns from those descriptive vectors, not a new inferential model.

## Heatmap correction families

Each row is its own original Fisher/BH view family, not the composition
trait family or ORA global family. Family sizes can change because the
original view-building/support rules depend on the corrected z-scores.
BH is applied before the final p/q display or selection cutoff. There is no
new BH correction across this comparison table or across selected cells.

| View | Tested cells before / after | BH-selected cells before / after |
| --- | ---: | ---: |
| Increased grouped | 188 / 190 | 5 / 5 |
| Decreased grouped | 185 / 186 | 4 / 4 |
| Increased individual | 117 / 117 | 2 / 2 |
| Decreased individual | 105 / 105 | 4 / 4 |
| Grouped direction contrast | 747 / 747 | 8 / 8 |
| Individual direction contrast | 2,379 / 2,388 | 4 / 4 |

The selected-cell identities, not just counts, are unchanged. Pattern
heatmaps retain exploratory raw p < .01; reported selected sets use
q <= .05 plus the original support filter. The full p/q values, support,
display status, signed effects where applicable, and membership are in
[heatmap_cells_before_after.tsv](heatmap_cells_before_after.tsv).

## PD-cellulose

The original PD comparison contains 130 increased profiles from 50 studies
and 123 decreased profiles from 46 studies. Other fecal profiles number
1,766 increased and 1,653 decreased. These populations are held fixed and
are not restricted to the core matched-pair subset.

| Quantity | Before | After |
| --- | ---: | ---: |
| PD increased positive-z profiles | 35 / 130 | 36 / 130 |
| PD decreased positive-z profiles | 67 / 123 | 68 / 123 |
| Increased PD vs other feces, Fisher p | .036666 | .057582 |
| Decreased PD vs other feces, Fisher p | .039135 | .030928 |
| Continuous-z PD-by-direction interaction effect | .400398 | .409777 |
| Continuous-z interaction p | .085963 | .080037 |
| Positive-z interaction effect | .188516 | .185929 |
| Positive-z interaction p | .053321 | .055280 |

These PD p-values are unadjusted exploratory statistics. The clustered
interaction models keep 3,672 profiles and 666 study clusters. Both
interaction p-values remain above .05. Avoid claiming that the PD-specific
direction interaction is established. In particular, do not repeat the old
increased-profile raw p < .05 statement with the corrected denominator.
Figure wording/mechanistic interpretation remains subject to #179 review.

## Reproduction and safeguards

```bash
just sync
just composition-denominator-audit
just composition-input-replay /path/to/data /path/to/kg
just composition-downstream /path/to/data
just check
```

The raw replay reconstructs all retained counts exactly from the pinned
BugSigDBExports v1.3.1 and KG files. The downstream recipe uses the frozen
matrix, not `/path/to/data/derived/study_trait_profiles_wide.tsv`.
Auxiliary inputs are hashed in `before_provenance.json` and
`after_provenance.json`; common inputs must match between stages.
The recipe also checks them against `inputs/downstream_auxiliary_sha256.json`
in the parent report directory. The auxiliary `experiments.csv` is the
August 26 snapshot recorded in the original download provenance, not the
v1.3.1 full dump; its MD5 matches that snapshot's ETag exactly. This fixed
auxiliary vocabulary is used identically in both runs, not downloaded live.
The PD baseline reproduces its committed numerical summary exactly. The
follow-up is a controlled rerun on the issue's frozen 5,829-profile matrix,
not a claim to reproduce every historical follow-up artifact: Thalassemia
has no rows in this matrix's joined follow-up metadata, leaving two of the
three requested target labels represented and four guild-summary rows.
The base commit plus source hashes identify the rerun code; the comparison
provenance records shared helpers, versions, and the final cache check.

Sankey/PD loaders reject cached totals incompatible with the full retained
matrix. Guild-map loaders check row identity, totals, and cached transformed
values against the raw matrix. Zero-retained-signal profiles remain zero
fractions; their standardized scores need not be zero. Neither missing
values nor an obsolete cache are silently treated as current biological data.

`artifact_comparison.tsv` compares 138 tables; 130 differ exactly and eight
are identical. Exact differences include changed paths, row order, and
ordering within joined label lists, not only scientific changes.
`changed_columns.tsv` reports numeric changes with rtol=1e-9, atol=1e-12
on matched keys; tables without unique inferred keys retain full snapshots
but no misleading positional comparison. `heatmap_family_comparison.tsv`
preserves view-specific family sizes and cutoffs.

Adopt these corrected artifacts in #172/#173 integration work. Keep the
separate interpretation, study-dependence, taxon-influence, and layout
changes separate from this controlled denominator comparison.
