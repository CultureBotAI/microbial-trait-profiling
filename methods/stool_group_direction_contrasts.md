# Declared stool condition-group contrasts (#1093)

Trait-list mode adds a stool-group table alongside the existing full-corpus
contrast tables. The default six-trait figures and tables keep their existing
path and definitions. Run the declared mode with explicit inputs:

```sh
just trait-list-direction-contrasts \
  methods/trait_lists/headline_ten_phenotype.txt \
  data/derived/study_trait_profiles_wide.tsv.gz \
  trait-visualization/study_trait_profile_study_metadata.tsv \
  "$TRAIT_COLUMNS" "$ROW4_OUTPUT" --stool-group-count 6
```

`TRAIT_COLUMNS` names the column-label table for the supplied profiles;
`ROW4_OUTPUT` names a durable output directory. There is no automatic `.gz`
fallback: pandas reads the explicitly supplied gzip and provenance hashes its
compressed bytes. The trait file is parsed in declared order and its SHA-256 is
recorded. Its historical filename now contains nine traits; it does not select
traits based on this analysis. `METPO:1000603` denotes annotated **anaerobic**,
not obligately anaerobic.

## Population and cell quantity

Use the canonical human 16S/WMS profile metadata and its matching full profile
table. For the pre-union development inputs, this is 5,629 profiles from 1,047
studies with 344 retained non-MediaDive traits. F's full-phenotype core has 7,111
profiles across all hosts and is a different population. Supplied metadata keys
select the consumed profiles; the script does not turn arbitrary metadata into
the canonical human population. Cached totals must reconcile to all retained
input traits, including those absent from the declared display list.

Stool is **literal metadata `body_site == "Feces"`**. Multi-site strings such as
`Feces,Blood`, whitespace variants, and grouped anatomy labels do not qualify.
The trait-list path preserves the raw field for this comparison and refuses
metadata without it. Neither `body_site_group` nor the grouped-mode Feces label
substitutes for the raw field.

For each profile, divide annotation counts by the total across all retained
non-MediaDive traits and transform each displayed trait with
`log1p(1000 * fraction)`. Standardize each trait over the complete consumed
profile corpus with population standard deviation (`ddof=0`), **before** applying
the stool filter. Positive z therefore means above that corpus's average, not
above the stool-only average. Zero-total profiles retain the native handling.

For each grouped condition and trait, the selected table reports
`100 * (positive_increased / all_increased - positive_decreased / all_decreased)`
in percentage points. Fisher's two-sided test uses those positive/nonpositive
profile counts and is reported nominally. The native BH values remain available
for the complete stool-group cell family; they do not select groups. Profile
lists are the inference units here, not studies or patients. These are contrasts
in trait annotations attached to reported taxa, not measurements of physiology.

## Declared group rule and outputs

`--stool-group-count` is a positive integer with default 6. A group is eligible
when at least one declared-trait cell passes the existing support filter:
at least three positive-z profiles in either direction and at least two distinct
studies. Eligible groups are ordered by distinct stool-study count descending,
then condition-group name ascending. The first requested number is retained;
if fewer are eligible, all eligible groups are retained. Neither p nor q enters
this ranking. Every available cell of a chosen group is retained with its own
support flag; unsupported cells are not silently promoted. As in the native
contrast tables, a trait with no positive-z profiles in either direction has no
edge/cell, rather than an invented zero-effect cell.

Additional outputs are:

- `trait_list_stool_increased_grouped_edges.tsv` and its decreased counterpart:
  edges rebuilt from the exact-stool profiles, with native study/profile support.
- `trait_list_stool_all_group_direction_contrasts.tsv`: the entire stool family.
- `trait_list_stool_group_direction_contrasts.tsv`: cells of the selected groups,
  with `increased_minus_decreased_positive_percentage_points` added.
- `trait_list_stool_group_selection.tsv`: every represented group's stool-study
  count, eligibility, eligible rank and selection status.

`trait_list_provenance.json` records the exact stool rule, standardization
population, requested group count, support thresholds, ranking, chosen groups
with counts, trait-list hash, input hashes and additional output hashes. The
original full-corpus contrast and edge tables remain available. Inputs with
separate `METPO:1000605` or `METPO:1000612` columns are labelled `pre-union`;
absence of those columns is labelled `undetermined`, not asserted to prove a
particular canonicalization policy.

## Pre-union development reference; regeneration belongs to #1133

The tracked gzip predates #1124's facultative union. Reading it through newer
code does not merge its three distinct facultative columns. Issue #1133 owns its
regeneration at P; #1093 does not refresh it. Final row 4 results must use the
canonical table rebuilt in the P worktree, with fresh provenance.

On the tracked canonical pre-union inputs, exact stool retains 3,493 profiles in
664 studies. The historical six study counts, in group order, were
94, 88, 61, 52, 46 and 45 under grouped-mode Feces. Literal Feces gives
93, 88, 60, 52, 46 and 45. The differences are the Parkinson's experiment
`bsdb:35645785-1` (`Feces,Blood`) and the HIV experiment `bsdb:32024712-1`
(`Oral cavity,Feces`), each contributing one increased and one decreased profile.
This is an intended population-definition change; the values are not tuned to
the historical reference and are not P results.

## Protecting density outputs

Declared-list density runs refuse a destination inside any `density_adjustment/`
path, including case variants and symlink aliases, or a directory already holding
`table1_density_adjusted_effects.tsv`. That filename is reserved regardless of
its current contents or file type, so partial or malformed earlier results are
also protected. Use a fresh destination for a declared-list rerun. The audit API
checks the destination before reading profile inputs or creating
output directories; the module CLI also checks before loading the declared trait
list. The default six-trait density command is unchanged.
