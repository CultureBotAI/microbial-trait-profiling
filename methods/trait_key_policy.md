# Versioned trait-key policy (#1129)

New core runs bind `trait_key_policy` in the LinkML `WorkflowRun`, the resolved
workflow parameters, and both matrix and feature-catalog producer steps. The
parameter participates in the execution fingerprint. Scientific verification
checks the same policy against the consumed KG, every matrix cell, metadata,
catalog and profile contracts. Changing the field alone cannot relabel a run.

`raw_object_v1` preserves each raw object CURIE. `facultative_union_v1` maps
`METPO:1000605`, `METPO:1000608` and `METPO:1000612` to `METPO:1000608` before
per-taxon set deduplication and ancestor rollup. Existing predicate handling is
unchanged: negated METPO-object assertions and relation-value keys retain their
predicates; positive METPO-object assertions remain bare class keys. A direct and
strain-derived annotation
of the same trait at the same taxon contributes once. Strain-rollup scope is an
independent parameter. Neither policy converts annotations to measured physiology.

The union is a local analytical mapping, not an upstream ontology merger.
The [pinned METPO template](https://github.com/berkeleybop/metpo/blob/c6f49473730488c3aa71c3da2ef8dd4931c809a7/src/templates/metpo_sheet.tsv)
maps the source vocabularies as follows:

| Raw object | Source | Source label |
| --- | --- | --- |
| `METPO:1000605` | BacDive | facultative anaerobe |
| `METPO:1000605` | Madin et al. | facultative |
| `METPO:1000608` | BacDive | facultative aerobe |
| `METPO:1000608` | Madin et al. | facultative |
| `METPO:1000612` | BactoTraits | Ox_facultative_aerobe_anaerobe |

The shared interpretation is an annotation for growth with or without oxygen.
Original KG assertions and literal node labels remain unchanged. This mapping
does not re-audit each underlying strain experiment or establish oxygen use in vivo.

Matrix column metadata and the coverage module's feature catalog carry
`trait_key_policy`, `trait_member_ids_json` and `trait_source_vocabularies_json`.
The JSON arrays describe the declared mapping; they do not imply that every member
or source supplied evidence to every cell. Unmapped features have singleton member
arrays and no invented source-vocabulary mappings. Existing profile-coverage
counts retain their definitions and remain linked to this catalog by the core run.

Historical F manifests lack an explicit policy. The verifier recognizes only
audited matrix-implementation SHA-256 and byte-size pairs from their recorded code
inventory. The raw producer bytes cover the development interval
`f8a57e2..0e49e87`: from `f8a57e2` (2026-10-22, inclusive), including F
(`511c891`), until just before `0e49e87`. These select `raw_object_v1`.
At `0e49e87`, #1124 introduces the distinct implicit union producer, which
selects `facultative_union_v1`. Earlier cores whose matrix implementation
predates `f8a57e2` are outside this audited range and are refused as unknown
historical implementations. Both recognized producers retain their historical
metadata layouts. Unknown producer bytes, partial declarations, and
contradictory step bindings fail closed. The compatibility resolver neither edits
old records nor executes historical code, and it skips no inventory or scientific
checks. Adding another historical producer requires an explicit compatibility review.

Current runs emit the extended metadata layout even when explicitly selecting
`raw_object_v1`. Reproducing a historical policy is distinct from reproducing the
historical output serialization. Historical counts in methods notes remain labelled
as F references; current union results must be regenerated before use for basis P.
