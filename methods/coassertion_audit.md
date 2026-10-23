# Phenotype co-assertions across strain scopes (#846)

The KG audit records co-asserted annotation classes under `off`, `phenotype`,
and `full`, using the matrix's shared scope and negation policies. `full` is a
virtual annotation census only. It does not launch the full-scope pipeline or
authorize a full-scope inferential result.

```bash
just audit-coassertions --edges KG/merged-kg_edges.tsv \
  --nodes KG/merged-kg_nodes.tsv --output-dir RESULTS/coassertion-audit
```

Run this only on the exact KG used by the final verified cores. The generated
`README.md` is the count-bearing methods supplement; its table is constructed
from `coassertion_summary.tsv`, never copied from historical issue numbers.
The lead independently verified the pinned-KG counts below during review #948.
The final rerun still requires its own hashed audit products; a methods table
does not replace the run-specific provenance record.

| Pair | Left class | Right class |
| --- | --- | --- |
| Aerobic and anaerobic | `METPO:1000602` | `METPO:1000603` |
| Obligately aerobic and anaerobic | `METPO:1000606` | `METPO:1000603` |
| Mesophilic and psychrophilic | `METPO:1000615` | `METPO:1000614` |
| Facultatively anaerobic and obligately aerobic | `METPO:1000605` | `METPO:1000606` |
| Aerobe and anaerobe guilds | `trait_guilds.OXYGEN_GUILDS["aerobe"]` | `trait_guilds.OXYGEN_GUILDS["anaerobe"]` |

These are literal class or guild memberships, without added ontology ancestry.
The four exact-class pairs are separate from the broader canonical guild pair.
The audit records the actual KG labels for every class and fails if their
metadata are missing. A missing taxon name remains explicitly missing.

The exact-class pair has the ID `aerobic_class_and_anaerobic_class`. It differs
from the label-support sidecar's `aerobic_and_anaerobic`, which compares
{`METPO:1000602`, `METPO:1000606`} with
{`METPO:1000603`, `METPO:1000607`, `METPO:1000611`}. These are different annotation
sets and their counts must not be treated as the same measure.

## Pinned-KG review counts

The lead's [independent review and recount](https://github.com/CultureBotAI/microbial-trait-profiling/issues/948)
used KG edges with SHA-256
`d13812a436bd0b1c51cfe20c97076912a131f72b2b8e8ca3bb2db15095ab9ebb`,
also recorded in `methods/table1-native-results/contribution_manifest.json`.
The counts below are distinct NCBITaxon carrier nodes before ancestor closure.

| Pair | off | phenotype | virtual full |
| --- | ---: | ---: | ---: |
| Aerobic class and anaerobic class | 141 | 222 | 222 |
| Obligately aerobic and anaerobic | 9 | 44 | 44 |
| Mesophilic and psychrophilic | 5 | 210 | 210 |
| Facultatively anaerobic and obligately aerobic | 0 | 74 | 74 |
| Aerobe guild and anaerobe guild | 420 | 613 | 613 |

No individual strain subject co-asserted both sides of any audited pair in that
snapshot. Every additional co-asserting taxon under `phenotype` therefore required
annotations from different original subjects. This is annotation union evidence,
not evidence of incompatible physiology in an individual organism. These are
lead-reviewed counts from the pinned KG, not a new analysis run by this PR.

`coassertion_summary.tsv` counts distinct NCBITaxon carrier nodes before
taxonomic ancestor closure. Each node contributes once to a pair, independent
of the number of source records or mapped strains. `coasserted_taxon_nodes.tsv`
lists each node, its exact left/right annotations, original direct support,
and supporting strain IDs. It distinguishes a direct taxon subject carrying
both sides, a single strain carrying both, and a union requiring different
subjects. The first two flags may overlap and must not be added. A strain
mapped to two taxa can contribute one carrier to each; duplicate links do not
increase counts. `coasserted_strain_nodes.tsv` retains unlinked strains, which
add no taxon carrier. `coassertion_evidence.tsv` retains original edge IDs,
subjects, predicate/object IDs and individual primary source fields.
Unlinked co-asserting strains retain their source evidence with
`origin=unmapped_strain` and an empty target TaxID; no taxon is invented.

These inventories do not count ancestor-derived matrix overlap: a genus row
can carry both classes through different descendants even when no carrier
node co-asserts them. Nor do they count BugSigDB signatures, studies or
patients. They complement the reported-taxon census in the #773 label-support
audit. Input and implementation hashes, parser repair counts, class/pair
definitions and exact output hashes are recorded in `provenance.json`.
Its `code_dirty_at_execution` flag records whole-repository Git status, including
modified data or documentation; the implementation hashes identify the code bytes.

## Biological meaning

Co-assertion is not sufficient evidence of mutually exclusive physiological
states in one organism. Broad aerobic and anaerobic annotations can describe
compatible capabilities. Different strains, culture conditions, source
definitions, or errors can also account for overlap. A shared NCBITaxon node
does not establish a shared experimental subject. Even same-strain annotations
need their underlying measurement contexts before calling them contradictory.
The current union rule records available annotations; it does not estimate the
fraction of strains or the abundance of a phenotype in a patient.

## Proposed oxygen-guild sensitivity; no resolution selected

The proposed primary sensitivity excludes the **carrier nodes** that co-assert
both canonical oxygen guilds under that run's scope, before taxonomic closure.
Remove each excluded carrier's contribution from every retained feature, then
reconstruct the affected reported-taxon vectors and profiles. Keep the same
reported taxa, oriented experiments and initial study population. Removing a
whole carrier from the retained denominator is explicit; deleting only an
oxygen label would define a different estimand.

Run the existing native oxygen-guild estimators with and without this exclusion
on verified `off` and `phenotype` cores separately. Record counts of removed
carriers, affected reported taxa by rank, profiles and studies. A sensitivity
profile with zero retained denominator must be unavailable, not an imputed
zero. Report attrition and compare both estimates on the intersection of
eligible pairs, retaining all paired differences before study aggregation.
Keep pooled and body-site results distinct and use the native family rules;
report it as a declared sensitivity on shared evidence, not independent
replication or a new primary claim.

Excluding every *reported taxon* whose matrix row carries both guilds would
also remove higher-rank taxa with different descendants. That is a separate,
more aggressive sensitivity and must be named separately if the owner requests
it. Strain-majority or fraction-based replacements would additionally depend
on uneven strain sampling and missingness; they are proposals requiring an
owner decision, not automatic corrections.

The audit neither executes these sensitivities nor changes the matrix, feature
catalog, arm plans or primary calls. Issue #846 remains open until the lead
publishes the run-bound audit and oxygen sensitivity, or the owner records a
resolution decision. The reviewed counts above do not settle that sensitivity.
