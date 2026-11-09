# Exact-end validation of identifiers and digests

The five authoritative LinkML schemas (`mtp`, `bsdbng`, `workflow`,
`core_products` and `multiple_testing`) require patterned values to reach the
actual end of the string (#1316). Python's `$` also matches before a final LF;
a bare `^...$` therefore accepted a digest or identifier followed by that newline.
Such a value could validate and still differ from a computed digest or a taxon key.

Patterns now end with `(?![\s\S])$`. The negative lookahead requires that no
character remain. Leading anchors and pattern bodies are preserved. Input is
never stripped or normalized: a final LF fails validation. The permissive bodies
of Study identifiers and DOIs remain permissive; this correction is not a general
whitespace prohibition or a new DOI syntax rule.

All five Python models are generated from these schemas by `just gen-pydantic`.
Their custom validators compile the source pattern with Python `re` and apply it
to scalars or each item of a list. LinkML's JSON Schema generator carries the same
pattern into its output. The regression tests exercise generated models, schema
validation and LinkML-generated JSON Schema, including valid controls. The pattern
inventory tests cover every declared and generated pattern and reject bare-dollar
regressions. `just check-generated` verifies that the models match their sources.

This validation correction changes no fields, defaults, scientific quantities or
canonical serialization of valid records. Schema format versions remain unchanged.
Source and model digests change, so new execution provenance changes as expected.
Historical records and archives are not rewritten to conceal invalid values.
The committed-record audit compares old and new validation outcomes, records
inherited failures separately, and identifies any newly rejected values. Only a
confirmed trailing newline in an identifier or digest warrants a data correction
under this task; other newly failing records require review.
