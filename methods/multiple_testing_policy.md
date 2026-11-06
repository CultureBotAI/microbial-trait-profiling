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
every correction family in the code. `tests/test_multiple_testing_registry.py` checks
it against the correction calls it can recognise in the code. The test has known limits,
listed under [Enforcement](#enforcement).

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
  - heatmap cell families;
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
| `launcher` | A family flag: the family's anchor starts its computation without itself making the correction call. |
| `derived_procedures` | BY computed from a literal BH result in the supported form (see Enforcement), with its reason. |
| `size_source.unfiltered_count_rule` | For a selected family, how to count the unfiltered inventory (R6). |

The registry also pins `literal_procedure_counts`: for each discovered function, how many
literal BY and BH calls it makes.

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
wrapper functions. Among other checks, it fails when:
- a correction call has no registry entry for its enclosing function;
- a registry anchor does not exist, or an infrastructure anchor gives no reason;
- the literal procedures in a discovered function do not match the per-function counts
  pinned in `literal_procedure_counts`. A change that alters a count is caught; a swap
  between two calls that keeps the counts is not;
- a family's declared primary or comparison procedure does not appear among its own
  anchors' procedures, or a site where every call is literal lacks a declared procedure;
- a `derived_procedures` claim lacks supported arithmetic in its anchor. The supported
  form is a literal BH result multiplied by the same vector length's harmonic sum and
  clipped to [0, 1] to derive BY. A reason alone does not establish the derivation.
  The arithmetic is checked only where an anchor claims it in `derived_procedures`;
  an undeclared arithmetic derivation does not cause a registration failure;
- a launcher is declared where none is needed, or a selected family's selection rule merely
  repeats its definition or lacks an unfiltered-count rule;
- an entry is malformed. That is: a required field is missing or empty; a procedure is
  outside {BY, BH, none}; a `family_id` or a BH-primary reason is duplicated; an
  `outputs` item or a `policy_rule` value is invalid; an anchor or consumer does not resolve
  in the source; or an `open_issue` value is not an integer. The test runs offline, so it
  does not check that an issue exists.

Discovery reads files lazily and keeps compact per-file summaries. Alias expansion,
source size and summary size have explicit limits. The derivation recognizer runs only
in functions whose own resolved calls include both `benjamini_hochberg` and
`numpy.clip`. A separate per-function `Limits.derivation_steps` budget defaults to
100,000 work units (#1581). It charges AST visits, expression/statement analysis,
loop fixed-point passes, declaration-name processing, and the sizes of copied,
joined or invalidated fact environments before doing that work. Exceeding the budget
refuses discovery with the source path, anchor and limit; it does not silently skip
the anchor or credit a partial witness. This deterministic work bound is not a
wall-time or RSS guarantee.
`global` and `nonlocal` declarations conservatively invalidate the declared names
(#1609), including read-only vectors: declaration-aware resolution of enclosing
bindings is outside the recognizer. An undeclared free-vector read may be accepted;
this is a documented false-negative boundary, not proof that a declared read is
invalid Python or invalid arithmetic. Removing the invalidation would also admit
declared callable bindings that the recognizer cannot verify.
Derivation facts are reduced while each file's AST is resident; the summaries retain
neither the AST nor source text. The check establishes a supported calculation in an
anchor, not the ancestry of each family's output. At `run_composition_calibration`,
the observed family's derived BY is accepted on the strength of the null-replicate
calculation; the observed arithmetic is not checked. A calculation in a shared anchor
can therefore satisfy several families' claims without establishing each output's
derivation.

The test checks registration, anchors and literal procedures. It recognises a correction
by the name of the function it calls, so it does not detect any of these (#1443):
- hand-rolled arithmetic that is not named like a correction, such as a step-up loop or
  `p * m` Bonferroni;
- a correction function re-exported under another name;
- `functools.partial` and `getattr` dispatch;
- calls between scripts that import each other by bare module name;
- methods called on objects. Conversely, any `obj.correct()` or `obj.adjust()` call is
  reported as a correction even when it is not one.

Within a registered function, it also does not detect:
- a procedure chosen from a variable or a loop;
- BY derived arithmetically from BH where no anchor claims it in
  `derived_procedures`; the arithmetic check validates claims rather than discovering
  all derived procedures;
- arithmetic outside the supported local BH-to-BY form. Opaque helper results,
  mutation through calls (including closures), and alternative forms such as
  `np.minimum`, a list-comprehension harmonic sum or a `pd.Series` wrapper are not
  interpreted. Derivations inside unsupported `try` or `match` blocks, or stored only
  through a subscript or attribute rather than a local name, are not credited.
  Escaping `break` or `continue` paths inside those unsupported blocks are not propagated;
- other binding forms, including starred targets, `def`/`class` rebinding (including
  decorated definitions), `global`
  writes from another function and context managers that suppress a `raise`, are not
  modelled;
- whether a witnessed calculation reaches an output or is subsequently used. The
  check rejects constant-false branches, discarded expression statements and calculations
  after explicit block terminators, but does not prove general reachability or evaluate
  nonconstant conditions or loop exit/else path feasibility;
- a swap between two calls that keeps the per-function counts, for example exchanging the
  BH and BY calls that fill the screen's reported columns;
- a new family added as loop data that makes no new call;
- a family deleted from the registry, or never registered, when its function also hosts
  another registered family;
- a primary and comparison procedure exchanged, when the function makes both calls.

These need source review. The test does not show that a correction is statistically valid,
or that a family's membership matches its definition.

## Choices not adopted

Each of the following would change reported q-values, and so a pinned analysis rerun.
None is adopted here; each is recorded for the next analysis version:
- **D1.** Make BY primary, with BH as the comparison, for the BH-primary families that do
  not already export BY. These are every BH-primary entry in the registry whose
  `comparison_procedure` is empty, other than the ORA families (D2), for example:
  - the six-trait Table 1 family and the density family;
  - the coupling and heatmap cell families;
  - the native pair-t and study-t composition families;
  - the ratio-sensitivity S0, S1 and S4 fraction families, and the population fraction
    family;
  - the source-lineage influence fraction families;
  - the PD script families and the cellulose summary.

  Where BY is already exported, as for `mt-core-three-test-native-selected`, designating it
  primary changes no number.
- **D2.** Report the joint-tail BY family as the ORA q, instead of the native per-tail
  BH families.
- **D3.** Hold untestable members at p = 1 in the families that correct over tested rows
  only: native influence composition, the PD script families, the coupling ORA profile and
  the cellulose summary. The core native composition families already hold them at p = 1.
