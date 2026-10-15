# Superseded branch-label outputs

These four tables are preserved byte-for-byte from commit
`925b91c79fae4c052fb9abe3608296ccaf668d00` on
`clustered-heatmap-generation`. They record the earlier grouped-figure branch
labels and their underlying calculations; they are not outputs of the current
clustering method.

The reconciliation of issues #198, #220, and #221 replaces that method's
zero-filled, input-order-dependent distances with canonical ordering,
pairwise-complete Spearman distances, and a separate block for rows without a
defined neighbor. The current tree is a descriptive ordering, without branch
significance labels. Trait-pattern tests on branches selected using the same
trait patterns do not provide independent evidence for those patterns.

Use the current `six_shared_trait_signed_direction_cluster_tree_nodes.tsv`,
`six_shared_trait_signed_direction_cluster_tree_leaves.tsv`, and
`six_shared_trait_signed_direction_cluster_comparison.tsv` in the parent
six-trait directory for current clustering diagnostics. The archived tables
must not be combined with the regenerated tree or cited as current results.
