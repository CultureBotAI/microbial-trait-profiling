# Input population checks for the ORA refresh audit

Date: 2026-10-19. Follow-up issue:
[#725](https://github.com/CultureBotAI/microbial-trait-profiling/issues/725).

The earlier audit checked each result row and compared complete table bodies,
but equality could conceal the same population error in both files. Five
actual-producer fixture mutations reproduced false passes: duplicate tests,
unknown profiles, unknown traits, incorrect background trait counts, and a
missing entire trait family whose declared counts had also been reduced.

The audit now derives the expected family from the explicitly supplied inputs:

- Experiment pairing is determined before excluding profiles with zero
  `taxa_count`. Only paired profiles with a positive draw count are testable.
- Retained traits come from the profile header, excluding MediaDive columns.
  A trait is eligible only when it has positive counts in at least one row of
  the supplied taxa matrix; a missing matrix column has background count zero.
- Matrix taxon identifiers must be nonempty and unique. Its observed number of
  taxa defines `M`; positive annotation rows define each trait's `K`.
- Declared paired-profile counts, testable-profile counts and `M` must match
  these inputs. Every result row must have an eligible profile and trait, use
  the matrix's `K`, and occur exactly once in the complete profile-by-trait
  family. This also rejects a completely omitted eligible trait.

Validation streams the tables. It retains input identity indexes and one bit
per eligible profile for each trait, rather than storing every result-row key
or loading complete result tables into data frames. Zero-discovery output and
zero-eligible families remain valid. Input identity and draw-count validation
is deliberately stricter than the producer's treatment of ignored unpaired
rows: malformed identifiers or draw counts in those rows are rejected.

These checks establish input-defined population membership and background
counts, in addition to the existing byte equality, metadata, fingerprint,
numeric-domain and significance-rule checks. They do not independently
recalculate profile trait counts, effect estimates, hypergeometric probabilities
or BH-adjusted values. Those result values must remain byte-identical to the
historical export. The audit does not establish physiological activity or
causality, and its input fingerprints are not a substitute for an execution
record from the producer.

The five-pair production audit passed. Its JSON is byte-identical to the
committed `2026-10-19-ora-provenance-refresh-audit.json` snapshot (`cmp` exit 0).
The check used the existing manifest and historical/refreshed analysis files;
its separate output is `issue725-audit.json` under the original issue #720
output root, with log `logs/issue725-audit.log`. The analysis files were not
modified. The linked implementation PR records focused, full-check and CI
results for the code change.
