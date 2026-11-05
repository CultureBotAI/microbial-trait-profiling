# Multiple-testing policy

This is the project-wide rule for every p-value and q-value the repository reports. It
applies to every arm, sidecar, script, figure, deck and document. It sits above two
narrower documents, which keep their detail:
- [`multiple_testing.md`](multiple_testing.md), the contract for the additive
  full-family sidecar (policy v1);
- the arm notes: [`trait_direction_screen.md`](trait_direction_screen.md),
  [`three_arm_statistics.md`](three_arm_statistics.md),
  [`three_test_table.md`](three_test_table.md) and
  [`trait_coupling_sensitivity.md`](trait_coupling_sensitivity.md);
- the density audit, [`../docs/density_adjustment_audit.md`](../docs/density_adjustment_audit.md);
- the sidecar guide, [`../docs/multiple_testing_sidecars.md`](../docs/multiple_testing_sidecars.md).

Where one of them disagrees with this policy, it is corrected in a follow-up change, which
is tracked as an issue. If a correction would change a reported number, the note keeps the
current practice, and the change is listed under [Choices not adopted](#choices-not-adopted).

[`multiple_testing_families.yaml`](multiple_testing_families.yaml) is the registry of
every correction family in the code. `tests/test_multiple_testing_registry.py` keeps
it in step with the code (see [Enforcement](#enforcement)).

This policy changes no computation and no reported number. It records the current
procedures, names the deliberate differences between families, and sets the reporting
rules. Procedure changes that would move reported numbers are listed under
[Choices not adopted](#choices-not-adopted) and need an owner decision.

## Rules

**R1. Every corrected value names its family.** A q-value is reported with four things:
- its family: an identifier and a one-line definition;
- the family's size in the run that produced it;
- its procedure, BH or BY;
- its threshold.

In a table, these go in columns or in a provenance file beside the table. Examples are
`family_id`, `family_size` and `correction_method` in the Table 1 shift table, and the
density provenance. In prose, a figure or a deck, they go in the text or the caption.

**R2. Family sizes come from the run, never from a constant.** A family's size is the
number of hypotheses the run declared, counted from the executed inventory. An untestable
member is handled in one of two ways, and the output says which:
- it stays in the family, held at p = 1 internally, and is counted;
- or it is left out by a declared selection rule recorded with the output (R6).

Family sizes printed in documents and plates describe the run they were taken from.
Quote them with that run, never as the current size.

**R3. An uncorrected p-value is labelled nominal.**
- Its column name, a field beside it, or its caption or legend says "nominal".
- Any display that shows uncorrected p-values states how many tests it shows.
- A threshold applied to a nominal p is a nominal threshold, and is labelled as one.
- Calibration and audit p-values are diagnostics of a procedure, not inferential
  results, and are named as diagnostics.

**R4. Procedures.**
- **BY is the primary procedure for a screen**, that is, a family holding every candidate
  in a declared catalog. Examples: all retained traits, all ordered ratios, all
  profile × trait × tail cells. BY controls the false discovery rate under arbitrary
  dependence. Trait annotations are strongly dependent through shared taxa and curation
  density, and that dependence has not been shown to be of the kind under which BH is
  valid.
- **BH is reported beside BY** where both are produced.
- **BH is primary only where the registry records that choice and its reason**, in the
  entry's `deliberate_difference`. Today this covers:
  - pre-specified focal families, such as the six Table 1 traits and the coupling focal
    families;
  - the density reported-trait family (#1321);
  - the native, stratified and site ORA families, which are directional views;
  - heatmap and map cell families;
  - selected native composition and influence families (R6);
  - the cellulose summary and the PD script families, which are exploratory;
  - sensitivity analyses that mirror a native BH family.

  Moving these families to BY may change reported numbers. See
  [Choices not adopted](#choices-not-adopted).
- **A family with primary `none`** is either a diagnostic (calibration, audit or
  migration) or a sensitivity bound (R9).

**R5. One hypothesis may sit in several families.** The same raw p-value may be
corrected in more than one family. Examples:
- the full screen and the six-trait family;
- a selected family and its full catalog;
- a focal family and its full family.

Each q column then names its family under R1. A q from one family is never relabelled as,
or compared as if it were, a q from another. The registry records which family is
primary, and why the others exist (`deliberate_difference`).

**R6. Selected families declare their selection.** Some families are built after a
support, variance or significance filter (`size_source.kind: selected_rows` in the
registry). Each one must record its selection rule, and the size of the unfiltered
inventory, with its output. The registry gives the selection rule and the recipe for that
count for every selected family. Producers that do not yet write the count are tracked in
#1450. Correcting within a selected family does not repair the selection, and
the text says so. Significance, effect-size and presentation filters never shrink a
family after testing.

**R7. Thresholds.**
- Each family's threshold is the one its registry entry declares.
- q ≤ 0.05 (inclusive) is the default for a new family.
- Existing families keep their declared thresholds. For example, native ORA calls a hit
  at q ≤ 0.001 together with a fold cutoff. Some exploratory and calibration families
  use other levels.
- A display that thresholds on q states the family and the threshold.

**R8. Documents.**
- A significance statement in `methods/` or `docs/` quotes q and its family.
- A statement read from p alone says it is nominal. This includes historical analyses
  that were read from uncorrected p.
- A q-value recomputed for a document says it was recomputed, and from what.
- No causal or severity word goes on a statistic that does not test it. "Driven by" is
  not attached to a principal component's share; "severe" is not attached to a VIF
  without a declared cutoff.

**R9. Effective family sizes are sensitivity bounds.** Effective test counts (Li & Ji,
Galwey) are never a primary family size. They are reported as sensitivity bounds that
carry no FDR guarantee, as [`trait_direction_screen.md`](trait_direction_screen.md) rules
(#265).

## The family registry

[`multiple_testing_families.yaml`](multiple_testing_families.yaml) has one entry per
correction family. Each entry is anchored to code by path and qualified function name,
not by line number. An entry gives:

| Field | Content |
|---|---|
| `family_id` | A new `mt-*` identifier, declared by the registry. Outputs keep their own labels. |
| `area`, `hypotheses`, `family_definition` | What is tested, and which tests are corrected together. |
| `size_source` | `executed_inventory`, or `selected_rows` with its selection rule and count rule (R2, R6). |
| `primary_procedure`, `comparison_procedure`, `threshold` | R4, R7. |
| `outputs`, `consumers` | Where the values are written and read. |
| `policy_rule` | The rules that apply. |
| `deliberate_difference` | Why one hypothesis is also in another family (R5), and why a BH-primary family is BH-primary (R4). |
| `procedure_role`, `inference_role` | What the primary designation and the family's results are used for. |
| `conformance_status`, `open_issue` | Whether the family conforms, and the issue that tracks it if not. |
| `inventory_ids`, `joint_status_outputs` | The inventory records behind the entry, and outputs that combine several families' calls. |

Observed family sizes are not stored in the registry, because they belong to an
execution (R2). The registry was built from a quote-checked inventory of every correction
and its observed family sizes. That inventory is not published; the registry entries are
the published record.

To add or change a correction:
1. add or update its registry entry in the same pull request;
2. say which rule it follows;
3. record a `deliberate_difference` if the hypothesis is already corrected in another
   family.

## Enforcement

`tests/test_multiple_testing_registry.py` parses every Python file under `src/` and
`scripts/`, including ignored files. It resolves import and assignment aliases and local
wrapper functions. It fails when:
- a correction call has no registry entry for its enclosing function;
- a registry anchor does not exist, or an infrastructure anchor gives no reason;
- the literal procedures in an anchored function do not match its registered families,
  counted per procedure, so a single BY call changed to BH is caught;
- an entry is malformed: a required field is missing or empty, a procedure is outside
  {BY, BH, none}, a `family_id` is duplicated, an anchor or consumer does not resolve in
  the source, or an `open_issue` value is not an integer. The test runs offline, so it does
  not check that an issue exists.

Discovery reads files lazily and keeps compact per-file summaries. Bounded alias
expansion, together with cycle and size limits, keeps it fast and its memory small.

The test checks registration, anchors and literal procedures. It recognises a correction
by the name of the function it calls, so it does not detect any of these (#1443):
- hand-rolled arithmetic that is not named like a correction, such as a step-up loop or
  `p * m` Bonferroni;
- a correction function re-exported under another name;
- `functools.partial` and `getattr` dispatch;
- calls between scripts that import each other by bare module name;
- methods called on objects.

These need source review. The test does not show that a correction is statistically valid,
or that a family's membership matches its definition.

## Choices not adopted

Each of the following would change reported q-values, and so a pinned analysis rerun.
None is adopted here; each is recorded for the next analysis version:
- **D1.** Make BY primary, with BH as the comparison, for the BH-primary families that do
  not already export BY. These are the six-trait Table 1 family, the density family, the
  coupling, heatmap and map cell families, the PD script families, the native influence
  families and the cellulose summary. (Where BY is already exported, as for presence/depth
  focal families, designating it primary changes no number.)
- **D2.** Report the joint-tail BY family as the ORA q, instead of the native per-tail
  BH families.
- **D3.** Hold untestable members at p = 1 in the families that correct over tested rows
  only: native influence composition, the PD script families, the coupling ORA profile and
  the cellulose summary. The core native composition families already hold them at p = 1.
