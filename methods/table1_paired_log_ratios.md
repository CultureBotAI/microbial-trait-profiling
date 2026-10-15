# Paired log-ratio direction consistency for Table 1-trait pairs (matched-corpus analog)

Run `just table1-compositional-stats`. The implementation is
`src/bsdbng/trait_composition.py`. It writes `trait_log_ratio_stats.tsv`,
`trait_log_ratio_studies.tsv`, `trait_log_ratio_provenance.json`, and a readable
`trait_log_ratio_report.md` to `data/derived/table1_composition/`.

The shared loader retains human samples and paired increased/decreased profiles,
resolves qualified trait identifiers, and excludes MediaDive columns. This
analysis uses the six traits named in manuscript Table 1 and all 15 unordered
pairs as a matched-corpus analog sensitivity analysis (not a reproduction of
published Table 1 baseline values). It does not assume any one of those traits is
biologically stable.

## Hypothesis and sampling unit

For numerator A and reference B, each eligible experiment contributes

`d = log2(A_inc / B_inc) - log2(A_dec / B_dec)`.

Because `(A / total) / (B / total) = A / B`, unrelated denominator components and
common per-profile scaling cancel. A positive value means the A/B annotation
ratio is higher in the increased signature. It does not identify which member
changed absolutely. Log-ratios of annotation counts are not sequencing-count
models or measurements of patient-level trait abundance.

Experiment contrasts are averaged within each study before scoring their sign;
every eligible study receives one vote regardless of its number of conditions
or experiments. A two-sided exact binomial test evaluates
`P(study mean d > 0 | eligible, non-tied) = 0.5`. Effects within 1e-12 log2 units
of zero count as ties. The reported mean and median log2 changes include all
eligible studies, including ties. The exact 95% binomial interval describes the
fraction of non-tied studies with a positive contrast, not uncertainty in the
mean effect. Independence between studies remains an assumption.

## Zeros, coverage and multiplicity

An experiment is eligible for A/B only if both annotations are positive on both
sides. No pseudocount is added. An ineligible experiment is excluded as a whole
pair; remaining eligible experiments determine its study's mean. This conditions
inference on annotation coverage and can exclude different populations for
different ratios. The study table records every study, eligible/excluded pair
counts and its mean when defined; the summary also identifies partial coverage.

All-tied comparisons and comparisons without any positive complete pairs retain
separate statuses and blank undefined estimates/p/q-values. Neither state is
encoded as a 100% decrease. Internally these untestable hypotheses receive p=1
for adjustment so they remain in the declared 15-comparison family.

Benjamini-Yekutieli is used because the pair tests share components; its guarantee
allows arbitrary dependence given valid input p-values and a fixed family
([SciPy documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.false_discovery_control.html)).
The six traits were selected using overlapping data, so these remain exploratory
follow-up comparisons: BY does not remove that selection, annotation bias, or
dependence between studies. There is no combined FDR guarantee with the original
six-trait fraction tests or other analysis arms. No automatic significance
threshold is imposed by this report.

Each invocation defines its own 15-pair family. Running several body-site strata
does not establish FDR control across the union of those separate runs.

## Reproduction and interpretation

The provenance file records input hashes, the supplied Git revision, separate
hashes of the analysis module and of the shared loader `trait_guilds.py` (both
taken from the working files, so they reflect uncommitted changes), the BugSigDB
timestamp, input filter counts, trait IDs, study weighting, zero policy, and
testing family. The recipe appends `+dirty` to the recorded Git revision when the
tracked tree has uncommitted changes, because the revision alone would otherwise
describe code that was not the code that ran. Original
fraction-based outputs are retained for their original hypotheses. Comparing a
ratio result with them requires both the reference and its coverage exclusions;
a changed conclusion cannot automatically be attributed to normalization alone.

The [pinned-corpus validation table](table1_log_ratio_validation.md) records all
15 ratio results and the original fraction directions from the same input files
as a matched-corpus Table 1-trait analog sensitivity analysis.
Matched-population sensitivity analysis is tracked in
[#183](https://github.com/CultureBotAI/microbial-trait-profiling/issues/183).

For an exact body-site stratum, use a separate output directory:

```bash
just table1-compositional-stats --body-site Feces --output-dir data/derived/table1_composition_feces
```

No top-N condition selection, display masking or Other pooling defines that
inference corpus. `--body-site` matches the source label exactly, which also
excludes composite labels: of 221 distinct labels over 4,435 experiments, 59 are
comma-joined composites covering 141 experiments, and 18 of those name Feces
alongside another site (37 experiments) against 2,778 experiments labelled
exactly `Feces`. Those composite experiments are absent from every single-site
stratum rather than assigned to one, so a stratum is a clean subset, not a
partition of the corpus. This recipe rejects missing trait columns and blank/truncated
counts in retained trait columns instead of treating them as observed zeros. The analysis
rejects duplicate profile keys, incomplete pairs, mismatched pair metadata and
negative counts. The output is data and Markdown, so no figure renderer is needed.

References: [Morton et al. (2019), reference frames](https://doi.org/10.1038/s41467-019-10656-5)
and [Lovell et al. (2015), relative-data associations](https://doi.org/10.1371/journal.pcbi.1004075).
