# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [3.0.0] - 2026 October 28

#### Added
- Read typed node, branch and root-stem annotations with `parse_meta=True`. Unknown keys are retained in `extra_meta`. Repeated keys on the same node or branch raise `TreeError`.
- `MISSING` represents an explicitly missing annotation value, distinct from unset None. It retains its identity through copying and multiprocessing.
- `Tree.rooted` records True, False or unspecified rootedness. `Tree.root_suffix` holds root labels and optional stem data separately from the topology.
- `read_tree_text()` returns complete original tree descriptions. `skip_trees()` discards burn-in without constructing Trees.
- `Tree.unroot()` explicitly declares a tree unrooted. `collapse_bifurcating_root()` remains a structural operation, with `deroot()` as its deprecated alias.

#### Changed
- Requires Python 3.10 or newer.
- Consolidated duplicate tests and corrected an invalid quantile-error expectation.
- Nexus output uses the Tree's rootedness and `print_rootedness=True` to emit known R/U declarations. False suppresses the declaration. Printing does not reroot the tree.
- Printing options are passed directly to `newick()` or `nexus()`. Omitted precision or `precision=None` retains full stored float precision. Use a numeric precision for significant-digit output.
- With `print_meta=True`, None selections include populated annotations, empty selections omit them, and explicit selections use complete keys. TreeSet forwards the same controls. Branch/stem annotations require a printed length.
- Readers, counting and burn-in skipping share nested-comment handling. Raw wrappers retain comments and suffixes. TREE/UTREE declarations and translation tables are handled consistently across blocks.
- Rooting, construction, pruning, extraction, grafting and SPR now maintain rootedness and root-suffix data according to the documented editing rules.
- Copies isolate retained mutable annotations. `copyattr=False` also omits node metadata. Extracting at the root copies its suffix, while proper subtrees keep the former incoming branch separate.
- Merged branches retain summed lengths and fresh non-length fields. Subdivision defaults to `copy_attrs="none"`. Explicit `"lower"` and `"both"` retain independent annotation copies.
- Grafting and `TreeSummary.update()` transfer donor records. Continue with the recipient afterward. Retained input representatives must remain unchanged while a summary uses them. See the API ownership rules.
- Length edits retain annotations without rescaling their statistical meaning. Interval write-back replaces only the exact produced keys and preserves unrelated metadata. Low-level length-slot writes require explicit cache invalidation.
- Unset `Nodestruct.state` and unestimated `Branchstruct.length_sd` now default to None, distinct from an empty state or zero SD.
- Faster parsing, reading, counting and burn-in skipping, including annotation-heavy files with metadata parsing off.
- Shortened the README and moved detailed upgrade examples to recipes. The API reference includes a 3.0 migration table.

#### Removed
- `PrintSpec`, `configure_basic_printing()`, `configure_sumtree_printing()`, `Tree.set_print_spec()` and `Tree.get_print_spec()`. Pass options to each output call.
- Reader/factory `strip_comments` and reader/Tree `below_root`. Use keyword-only `parse_meta`, raw reading, root-node metadata and `root_suffix` as appropriate.
- Public `.ci`, writer `ci_labels` and the `"ci"` output-group shortcut. Store/select complete keys such as `height_95%_CI` in `extra_meta`. Literal `ci` remains ordinary metadata. Metadata/declaration switches are keyword-only.
- `Nodestruct.nleaves` and its constructor argument. Descendant counts remain in summary records and topology-derived results.
- The unused `Branchstruct.merge(check_compat=...)` argument.

#### Fixed
- Source distributions include the public documentation and the fixtures needed to run the tests.
- Online means and standard deviations retain small differences between nearly equal observations, including when parallel summaries are merged.
- All-compatible consensus now considers the highest-frequency below-cutoff bipartition instead of skipping it.
- Malformed completed roots raise `TreeError`, and reusable parsers recover after errors. Duplicate translated taxon names are rejected.
- Pruning clears the old root suffix when the root is removed, while retaining the remaining Tree's rootedness and node data.
- Root-frequency selection no longer raises a TypeError when the selected root moves. Minimum-variance and root-frequency selection establish rootedness even when the chosen root is unchanged.
- Summary outputs have independent branches and mutable annotations. Accumulating bipartitions leaves input branches and their statistics unchanged.
- Selected invalid `extra_meta` containers raise an informative `TreeError`. Suppressed metadata is not validated.
- Copies retain nested mutable annotations and nondefault root-suffix fields, including labels selected through a custom field.
- Leaf SPR retains moved node metadata, including when it was the only annotated node.
- Grafting retains donor metadata when the recipient initially has none and gives the graftpoint default node data.
- Renaming leaves or internal nodes moves their metadata to the final identifier, including renaming during grafting.
- Merging root counts refreshes cached root frequencies and sorted summaries.
- Branch-length write-back refreshes derived length caches after assigning the root edges. Metadata setters retain valid topology caches.

See the [API migration table](docs/api.md#upgrading-to-30) and [recipes](docs/recipes.md) for replacement calls and ownership examples.

---

## [2.3.1] - 2026 September 22

#### Fixed
- `TreeSummary` and `CAHeightEstimator` now preserve the exact observed node height when all observations are effectively identical. The median and credible-interval endpoints therefore equal the observed height instead of the midpoint of its log bucket.

---

## [2.3.0] - 2026 April 4

#### Added
- Multi-site parsimony: `parsimony_possible_states()`, `parsimony_assign_fits()`, and `parsimony_count_changes()` now accept tuple-valued states in addition to single strings. When `.state` is set to a tuple of strings (e.g., `("mink", "A", "G")`), the Hartigan algorithm runs independently at each site. The parsimony score and change counts are aggregated across all sites. Single-string states continue to work exactly as before.
- `fitpref` in `parsimony_assign_fits()` and `parsimony_count_changes()` now accepts a tuple of per-site preferences in addition to a single string (which is broadcast to all sites). Individual entries may be `None` for no preference at that site.

---

## [2.2.2] - 2026 April 2

#### Changed
- `Tree.nodedict` is now always a plain `dict` (never `None`), but may be empty on a fresh tree. Code that assumed every node had an entry immediately after tree construction should use `tree.set_node_attribute(...)` or call `tree.ensure_nodedict()` first.
- `Tree.newick()` / `Tree.nexus()` now skip node meta-comments for nodes that have no `Nodestruct` entry, rather than assuming all nodes carry metadata.
- `Tree.cladedict()` now returns `{clade: SummaryNodestruct}` instead of `{clade: Nodestruct}`. This only affects callers that inspect the returned struct type directly.
- `Tree.set_blens_from_heights()` now raises `TreeError` when `nodedict` is empty, with a clearer error message about missing height metadata.

#### Added
- `Tree.ensure_nodedict()`: populates `Nodestruct` entries for every node. Use this when an algorithm needs to read or write attributes on all nodes (e.g., parsimony workflows).
- `SummaryNodestruct`: new class for per-clade summary accumulators used by `TreeSummary` and `CAHeightEstimator`. This separates cross-tree summary state from per-tree node metadata (`Nodestruct`).

#### Fixed
- `Tree.parsimony_possible_states()` now calls `ensure_nodedict()` internally, so parsimony methods no longer fail on trees with a sparse `nodedict`. The method still raises `TreeError` if no node has a non-empty `.state`.
- Tree-editing methods (`graft`, `split_off_children`, `add_node_on_branch`, `add_leaf`, `remove_leaf`, `collapse_bifurcating_root`, `reroot`) now keep `nodedict` consistent when adding or removing nodes.
- `cluster_n()` no longer writes temporary scratch fields onto `Branchstruct` during clustering.
- `Tree.get_node_attribute()` now returns the supplied default when a node has no `Nodestruct` entry, instead of raising an error.
- `Tree.set_node_attribute()` now creates a `Nodestruct` on demand for the target node.

---

## [2.2.1] - 2026 April 1

#### Fixed
- `Nodestruct` now includes the slot fields required by parsimony methods (`state`, `primary_set`, `secondary_set`, `optimal_set`, `fit`, `wasambig`), fixing failures when using `parsimony_possible_states`, `parsimony_assign_fits`, and `parsimony_count_changes` with slot-based node objects.

#### Changed
- Test suite reorganized and modernized: old tests migrated to pytest style, test files consolidated into a single combined test module, and test data moved into dedicated subdirectories under `tests/data/`.

---

## [2.2.0] - 2026 April 1

#### Added
- `Tree.randtree` now simulates proper Yule or coalescent tree topologies, producing more realistic tree structure. Default is to create an ultrametric tree with exactly the user-specified tree height. Setting `rate_sd > 0` adds lognormal noise to branch lengths, breaking ultrametricity.

#### Fixed
- NEXUS format output now consistently uses spaces for indentation (previously mixed spaces and tabs).

---

## [2.0.0] - 2026 March 17

This is a major release with breaking API changes. See the migration guide in `README.md` under "Upgrading from 1.x to 2.x" for full details.

#### Added
- Explicit summary-tree pipeline via `SummaryTreeBuilder` and `TreePostProcessor`, with optional `CAHeightEstimator` for common-ancestor heights and credible intervals.
- Credible interval support via `QuantileAccumulator` and the `trackci` / `ci_probs` parameters in `TreeSummary`.
- `PrintSpec` class plus helpers `configure_basic_printing(...)` and `configure_sumtree_printing(...)` for consistent, configurable Newick/NEXUS output.

#### Changed
- `Tree.graft(...)` now grafts onto a specified edge (parent → child) and uses `connect_length` / `connect_branchstruct` for the connecting branch.
- `Tree.spr(...)` updated to match the new grafting model and by default preserves total tree length.
- `Tree.subtree(...)` / `Tree.prune_subtree(...)` now return `(subtree, basal_branch)`; leaf subtrees return the leaf name as `str`.
- `Tree.reroot(...)` signature updated to use the same distance/fraction placement convention as other editing methods.
- `Tree.newick(...)` / `Tree.nexus(...)` now treat `None` as "use PrintSpec defaults".

#### Renamed (terminology: `depth` --> `height`)
All names relating to node heights — the distance from the tips to an internal node — have been renamed from `depth`-based to `height`-based, aligning with the convention in most phylogenetic software. **This is a breaking change; the old names are removed.**

| v1.x name | v2.x name |
|---|---|
| `TreeSummary(..., trackdepth=True, ...)` | `TreeSummary(..., trackheight=True, ...)` |
| `Tree.nodedepth(node)` | `Tree.nodeheight(node)` |
| `Tree.set_blens_from_depths()` | `Tree.set_blens_from_heights()` |
| `TreePostProcessor.set_mean_node_depths(tree)` | `TreePostProcessor.set_mean_node_heights(tree)` |
| `CADepthEstimator` | `CAHeightEstimator` |
| `build_sumtree(..., blen="meandepth", ...)` | `build_sumtree(..., blen="cladeheight", ...)` |
| `build_sumtree(..., blen="cadepth", ...)` | `build_sumtree(..., blen="caheight", ...)` |
| Node attributes `depth`, `depth_sd`, `depth_median` | `height`, `height_sd`, `height_median` |

#### Deprecated
- `Tree.deroot()` → use `Tree.collapse_bifurcating_root()`.

#### Fixed
- Multiple refactorings for improved performance (speed and memory) and bug fixes across tree editing, summary-tree construction, and output formatting.
