# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## Unreleased

#### Fixed
- An annotation following an extra root-closing parenthesis now raises `TreeError` instead of leaking an empty-stack `IndexError`.
- Translated Newick input now rejects duplicate decoded taxon names just as untranslated input does, including distinct translation IDs naming the same leaf.
- Shared Newick parsing now requires one completed topology, releases reusable parse state after errors, and preserves structural root labels/stem lengths in `root_suffix` without adding graph edges.
- Pruning now clears the old root suffix when surviving-child promotion removes that root, while preserving the remaining Tree's rootedness and surviving node data through leaf, subtree and wrapper pruning.
- Root-frequency selection now passes both bipartition sides to the current edge-lookup API, repairing a TypeError when the selected root must move. Its existing numerical positioning algorithm is unchanged.
- Successful minimum-variance and root-frequency selection now establish rootedness even when the chosen root already matches and no root movement occurs. Root geometry and unfinished root-frequency calculations are unchanged.
- Summary builders and branch annotation write-back now materialize independent output branches and mutable annotations, leaving quantile accumulators and scratch statistics in the summary. Repeated outputs can be edited independently.
- Bipartition accumulation now owns fresh records instead of adopting input branches, leaving sample lengths, annotations and existing accumulator state unchanged. Later observations reuse the accumulator.
- Selected invalid `extra_meta` containers now raise an owner-specific `TreeError`; metadata-disabled or empty-selection output bypasses annotation validation.
- Tree copies retain nondefault root-suffix fields, including root labels selected through a configured native field.
- Struct and Tree copies now isolate retained nested mutable annotation data. `copyattr=False` omits node metadata as well as branch attributes.
- Single-leaf SPR now preserves an independent copy of the moved node metadata, including when it was the tree's only annotated node, without changing leaf extraction or branch-length calculations.
- Grafting now retains donor node metadata even when the recipient initially has none, and initializes the inserted graftpoint with default node data.
- Leaf and internal-node renaming now moves node metadata to the final identifier, including duplicate-name resolution and donor renaming during grafting.
- Merging `TreeSummary` root counts now refreshes cached root frequencies and sorted root summaries, including newly observed root positions.
- `TreePostProcessor.set_biplen_on_existing_tree()` now invalidates length-derived caches after its final root-edge assignments, so warmed distance queries reflect the new lengths.

#### Removed
- Removed reader/factory `strip_comments` arguments and reader/Tree `below_root` storage. Use keyword-only `parse_meta`, explicit raw descriptions and structured root-node/root-suffix fields.
- Removed writer `ci_labels` and the `"ci"` output-group shortcut. Select complete exact annotation keys; literal `ci` is ordinary metadata. Trailing metadata/declaration switches are keyword-only, so removed positional calls fail instead of rebinding silently.
- Removed public `.ci` storage from ordinary and summary structs. Intervals now use complete keys such as `height_95%_CI` or `length_95%_CI` in lazy `extra_meta` dictionaries, alongside unknown annotations.
- Removed `Nodestruct.nleaves` and its constructor argument. Descendant counts remain in `SummaryNodestruct` and topology-derived calculations; summary construction no longer leaks bookkeeping into exported tree metadata.
- Breaking API change: removed `PrintSpec`, `configure_basic_printing()`, `configure_sumtree_printing()`, `Tree.set_print_spec()`, and `Tree.get_print_spec()`. Specify printing options directly in `Tree.newick()` and `Tree.nexus()`; the tree no longer stores formatting defaults.

#### Changed
- Shared parsing checks completion once and reports failed token/stack operations locally, avoiding redundant bounds checks on every successful child transition while retaining malformed-root rejection.
- Verified supported BEAST/MrBayes input ownership and selected semantic round trips, including explicit translation/label settings, typed values and Nexus rootedness versus intentionally undeclared Newick output.
- Ordinary branch-length parsing keeps direct structural lookup; only actual root stems use the suffix helper, avoiding a helper call and temporary stack slice for every represented edge.
- File readers now share buffered nested-comment framing, reusable annotation-aware parsers and lexical count/skip paths. `read_tree_text()` and legacy raw wrappers return complete descriptions; Nexus TREE/UTREE and translation context remain stable across blocks.
- `Tree.from_string(..., parse_meta=True)` imports typed positional node/branch/stem annotations into native fields and lazy `extra_meta`. Both metadata modes retain one leading R/U declaration, skip ordinary comments, and preserve structural root labels/stems. Same-owner repeated keys and unsupported metadata positions/directives raise `TreeError`.
- Added a reusable typed key/value metacomment decoder with exact-key duplicate rejection, recursive tuple collections, explicit missing values and doubled-quote strings. It has no topology or owner interpretation.
- Active API, README and recipes now document implemented ownership/root/edit contracts and explicit-copy workflows, correct the Nexus-default and summary-method descriptions, and keep parser/reader migrations distinct from completed behavior.
- Documented and verified direct length assignment, caller scaling loops and height-to-length conversion: retain stored annotations/rootedness/suffixes without statistical rescaling, and refresh length-derived caches. Raw slot writes require explicit invalidation.
- Documented and verified annotation attachment on existing reversed edges: objects and values are retained without transformation, while directional/rooted-clade statistics remain the caller's scientific responsibility. Merge/subdivision rules apply separately.
- Documented and verified complete SPR preservation of rootedness and surviving node/edge attachment, suffix retention only while the old root survives, and fresh merged/split/connecting fields. Qualified its existing root-adjacent connection-length adjustment.
- Documented and tested graft preservation of recipient rootedness at actual graft entry, independently of donor state, suffix information and connection length, including leaf grafts.
- Documented and verified explicit graft connection precedence: donor suffix data never contributes implicitly; a deliberately prepared suffix Branchstruct can be transferred through the existing connection argument.
- Documented consuming graft ownership for donor contents and supplied connecting Branchstructs, with verified explicit-copy usage when callers retain independently usable originals.
- Documented first-representative reference retention with `store_trees=True`, including the caller's no-edit lifetime through summary transfers and explicit-copy editing workflow. Read-only/cache work remains allowed.
- Documented consuming `TreeSummary.update()` ownership across branch, clade, root and topology records. Donor scalar counts remain available for progress; independent statistical reuse after transfer is unsupported.
- Constructors/builders now establish rootedness from their semantics: bipartition/NJ/consensus/MBC outputs are unrooted; clade/random-ancestral/UPGMA/MCC/HIPSTR outputs are rooted. Generic Tree/star/branch-table construction remains unspecified. Height estimation alone does not select a root.
- Extraction at the source root now independently copies its root suffix; proper subtrees have no suffix and keep the copied former incoming branch separate. Copy/extraction reuse the same suffix retention rules.
- Tree results from `subtree()` and `prune_subtree()` inherit source rootedness at extraction time. Leaf extraction retains the existing string-plus-basal-branch interface.
- Successful public `reroot()` establishes rootedness, including a current-root no-op; computational root movement and `match_nodes()` use shared mechanics without selecting an interpretation. Actual root replacement clears the old suffix while node data stays attached.
- `remove_branch()` now gives reattached branches fresh non-length fields and removes the deleted internal node's metadata, while preserving its existing length redistribution and unrelated data.
- Branch subdivision now defaults to `copy_attrs="none"`, including automatic graft/reroot subdivisions. Explicit `"lower"` and `"both"` retain independent annotation copies without rescaling statistics; an explicit upper Branchstruct still overrides the upper segment.
- Serial branch joins now consistently retain only summed structural lengths, including non-root parent removal in `remove_leaf()`. `Branchstruct.merge()` returns fresh non-length defaults; its unused `check_compat` argument has been removed.
- Computational interval write-back updates only the exact produced non-None height/length keys, preserving unrelated target metadata and copying retained mutable result values. Generic whole-field setters still replace `extra_meta`.
- Nexus output now reads `Tree.rooted` and uses keyword-only `print_rootedness=True` to include known state; False suppresses it. The removed `rooted` value argument is not a direct keyword rename. TreeSet output supports mixed states.
- Omitted precision and `precision=None` now emit full stored float precision in all four writers; explicit numeric precision still requests significant digits. Annotation floats retain a decimal point or exponent when reduced precision would look integral. sumt continues to request seven digits explicitly.
- Unset `Nodestruct.state` and `Branchstruct.length_sd` now default to `None`, distinguishing unlabelled/unestimated data from an actual empty string or estimated zero. Parsimony retains its existing no-state behavior.
- `Tree.set_branch_attribute()` now invalidates only length-derived caches for structural length writes and retains valid Tree caches for stored metadata updates, avoiding repeated topology recomputation during annotation.
- Printing is configured at the point of export to simplify the workflow and remove dependence on earlier configuration calls. `Tree.newick()` and `Tree.nexus()` have fixed defaults: lengths and labels enabled, `labelfield="label"`, `precision=None`, and `print_meta=False`. Omitted precision and `None` retain full stored float precision; the other scalar settings retain their explicit-None default behavior. Each output call is independent.
- With `print_meta=True`, omitted or `None` attribute selections now include all populated annotations, excluding internal accumulators and temporary branch-height caches. Empty sequences suppress that kind of metadata; explicit selections include only the named fields. Use complete interval keys in explicit selections; missing requested keys are skipped.
- Metadata strings are quoted, booleans use `true`/`false`, and collections use braces with deterministic ordering for sets. Unset None values are omitted; explicit MISSING, empty strings/collections, zero and `False` are retained.

#### Added
- Regression coverage for explicit no-copy branch subdivision and root-ID-only renaming, preserving suffix/state/node attachment without treating a rename as root replacement.
- `Tree.unroot()` deliberately establishes an unrooted interpretation, including two-leaf and unchanged-root cases. Structural collapse clears a displaced rooted interpretation to None; the deprecated `deroot()` alias keeps structural semantics.
- Programmatic `root_suffix` label/stem output under the existing explicit controls. Selected stem metadata requires an emitted stem length; an absent length never becomes an invented zero. Root-node annotations remain separate.
- Exact metadata-key resolution and selected owner/key value validation, typed float/MISSING/empty-value output, quoted keys and independently balanced brackets in annotation strings. Emitted branch annotations require emitted structural lengths.
- `Tree.rooted` and optional `Tree.root_suffix` data fields and copy support, preserving absent versus explicit-zero stem lengths. Constructor/edit transitions follow the contracts above; input interpretation remains part of the parser migration.
- Shared public `MISSING` annotation marker, preserving identity through copying and multiprocessing serialization. `None` continues to mean unset; required explicit node heights reject either missing state before changing lengths.
- `TreeSet.newick()` and `TreeSet.nexus()` accept and forward explicit metadata selections and `print_meta`.
- Regression tests for explicit and automatic printing options, metadata and CI output, independent calls, and TreeSet forwarding.
- Complete credible-interval computation and Nexus metadata export examples in the README, API reference, and recipes. The README example now specifies `ci_probs` as well as `trackci=True` to actually request intervals.

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
